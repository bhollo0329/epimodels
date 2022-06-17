#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
from numpy.random import binomial, hypergeometric
import xarray as xr
import xsimlab as xs
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

# from ..setup import epi
from episimlab.models.epi_model import EpiModel
from episimlab.foi import BaseFOI
from episimlab.compt_model import ComptModel
from episimlab.utils import (
    get_var_dims, group_dict_by_var, discrete_time_approx as dta,
    IntPerDay, get_rng, any_negative, visualize_compt_graph
)
from episimlab.partition import Partition, TravelPatRepeatDaily, ContactsFromCSV
from episimlab.setup.sto import SetupStochasticFromToggle
from episimlab.setup.seed import SeedGenerator
from episimlab.setup.greek import (
    gamma, sigma, rho, mu
)
import logging

logging.basicConfig(level=logging.DEBUG)



@xs.process
class SetupVaccineDoses:
    """Initialize vaccine doses"""
    DIMS = ('age', 'risk')
    max_daily_doses = xs.variable(dims=DIMS, global_name='max_daily_doses', intent='out')
    doses_delivered = xs.variable(global_name='doses_delivered', intent='out')
    eligible_pop = xs.variable(global_name='eligible_pop', intent='out')
    state = xs.global_ref('state', intent='in')
    _coords = xs.group_dict('coords')

    @property
    def dims(self):
        return self.DIMS

    @property
    def coords(self):
        return {k: v for k, v in group_dict_by_var(self._coords).items()
                if k in self.dims}

    @property
    def S(self):
        return self.state.loc[dict(compt='S')].sum('vertex')

    @property
    def eligible_pop(self):
        return self.state.loc[dict(compt=['V', 'E', 'Ev', 'Pa', 'Py', 'Ia', 'Iy', 'Ih', 'R', 'D'])].sum(
            ['compt', 'vertex'])

    def initialize(self):
        self.max_daily_doses = xr.DataArray(
            [[0, 0],
             [0, 0],
             [10, 20],
             [30, 40],
             [50, 100],
             [50, 100],
             [50, 100]],
            dims=self.dims, coords=self.coords)

    def run_step(self):

        assert self.eligible_pop.shape == self.max_daily_doses.shape

        flat_doses = self.max_daily_doses.values.flatten()
        flat_pop = self.eligible_pop.values.flatten()
        actual_doses = np.array([min(i) for i in zip(flat_doses, flat_pop)])
        doses_delivered = []
        for i in zip(self.S.values.flatten(), flat_pop, actual_doses):
            try:
                doses_delivered.append(np.random.hypergeometric(i[0], i[1], i[2]))
            except ValueError:
                doses_delivered.append(0)

        self.doses_delivered = xr.DataArray(
            np.array(doses_delivered).reshape(self.max_daily_doses.shape),
            dims=self.dims, coords=self.coords
        )


@xs.process
class RateS2E(BaseFOI):
    """FOI that provides a `rate_S2E`"""
    TAGS = ('model::ElevenComptV1', 'FOI')
    PHI_DIMS = ('age0', 'age1', 'risk0', 'risk1', 'vertex0', 'vertex1',)
    I_COMPT_LABELS = ('Ia', 'Iy', 'Pa', 'Py')
    S_COMPT_LABELS = ('S')
    rate_S2E = xs.variable(intent='out', groups=['edge_weight'])

    def run_step(self):
        self.rate_S2E = self.foi.sum('compt')


@xs.process
class RateS2V:
    """Vaccination dosage model"""
    rate_S2V = xs.variable(global_name='rate_S2V', groups=['edge_weight'], intent='out')
    doses_delivered = xs.global_ref('doses_delivered', intent='in')
    eff_vaccine = xs.variable(global_name='eff_vaccine', intent='in')

    def run_step(self):
        self.rate_S2V = xr.apply_ufunc(binomial, self.doses_delivered, self.eff_vaccine)


@xs.process
class RateV2Ev(BaseFOI):
    """FOI that provides a `rate_V2Ev`"""
    TAGS = ('model::ElevenComptV1', 'FOI')
    I_COMPT_LABELS = ('Ia', 'Iy', 'Pa', 'Py')
    S_COMPT_LABELS = ('V')

    # reference phi, beta from global environment
    phi = xs.global_ref('phi', intent='in')
    beta = xs.global_ref('reduced_beta', intent='in')
    omega = xs.global_ref('omega', intent='in')
    rate_V2Ev = xs.variable(intent='out', groups=['edge_weight'])

    def run_step(self):
        self.rate_V2Ev = self.foi.sum('compt')


@xs.process
class BetaReduction:
    beta = xs.global_ref('beta', intent='in')
    reduced_beta = xs.variable(global_name='reduced_beta', intent='out')
    beta_reduction = xs.variable(global_name='beta_reduction', intent='in')

    def run_step(self):
        self.reduced_beta = self.beta * self.beta_reduction


@xs.process
class RateEv2P:
    """Provide a `rate_Ev2P`"""
    rate_Ev2P = xs.variable(global_name='rate_Ev2P', intent='out')
    sigma = xs.variable(global_name='sigma', intent='in')
    state = xs.global_ref('state', intent='in')
    int_per_day = xs.global_ref('int_per_day', intent='in')

    def run_step(self):
        self.rate_Ev2P = dta(self.sigma, self.int_per_day) * self.state.loc[dict(compt='Ev')]


@xs.process
class RateE2P:
    """Provide a `rate_E2P`"""
    rate_E2P = xs.variable(global_name='rate_E2P', intent='out')
    sigma = xs.global_ref('sigma', intent='in')
    state = xs.global_ref('state', intent='in')
    int_per_day = xs.global_ref('int_per_day', intent='in')

    def run_step(self):
        self.rate_E2P = dta(self.sigma, self.int_per_day) * self.state.loc[dict(compt='E')]
        # DEBUG
        assert not any_negative(self.rate_E2P, raise_err=True)


@xs.process
class RateE2Py:
    """Provide a `rate_E2Py`"""
    rate_E2Py = xs.variable(global_name='rate_E2Py', groups=['edge_weight'], intent='out')
    tau = xs.variable(global_name='tau', intent='in')
    rate_E2P = xs.global_ref('rate_E2P', intent='in')

    def run_step(self):
        self.rate_E2Py = self.tau * self.rate_E2P


@xs.process
class RateEv2Py:
    """Provide a `rate_Ev2Py"""
    rate_Ev2Py = xs.variable(global_name='rate_Ev2Py', groups=['edge_weight'], intent='out')
    tau_v = xs.variable(global_name='tau_v', intent='in')
    rate_Ev2P = xs.global_ref('rate_Ev2P', intent='in')

    def run_step(self):
        self.rate_Ev2Py = self.tau_v * self.rate_Ev2P


@xs.process
class RateE2Pa:
    """Provide a `rate_E2Pa`"""
    rate_E2Pa = xs.variable(global_name='rate_E2Pa', groups=['edge_weight'], intent='out')
    tau = xs.variable(global_name='tau', intent='in')
    rate_E2P = xs.global_ref('rate_E2P', intent='in')

    def run_step(self):
        self.rate_E2Pa = (1 - self.tau) * self.rate_E2P
        # DEBUG
        assert not any_negative(self.rate_E2Pa, raise_err=True)


@xs.process
class RateEv2Pa:
    """Provide a `rate_Ev2Pa"""
    rate_Ev2Pa = xs.variable(global_name='rate_Ev2Pa', groups=['edge_weight'], intent='out')
    tau_v = xs.variable(global_name='tau_v', intent='in')
    rate_Ev2P = xs.global_ref('rate_Ev2P', intent='in')

    def run_step(self):
        self.rate_Ev2Pa = (1 - self.tau_v) * self.rate_Ev2P


@xs.process
class RatePy2Iy:
    """Provide a `rate_Py2Iy`"""
    rate_Py2Iy = xs.variable(global_name='rate_Py2Iy', groups=['edge_weight'], intent='out')
    rho_Iy = xs.variable(global_name='rho_Iy', intent='in')
    state = xs.global_ref('state', intent='in')

    def run_step(self):
        self.rate_Py2Iy = self.rho_Iy * self.state.loc[dict(compt='Py')]


@xs.process
class RatePa2Ia:
    """Provide a `rate_Pa2Ia`"""
    rate_Pa2Ia = xs.variable(global_name='rate_Pa2Ia', groups=['edge_weight'], intent='out')
    rho_Ia = xs.variable(global_name='rho_Ia', intent='in')
    state = xs.global_ref('state', intent='in')

    def run_step(self):
        self.rate_Pa2Ia = self.rho_Ia * self.state.loc[dict(compt='Pa')]


@xs.process
class RateIy2Ih:
    """Provide a `rate_Iy2Ih`"""
    rate_Iy2Ih = xs.variable(global_name='rate_Iy2Ih', groups=['edge_weight'], intent='out')
    eta = xs.variable(global_name='eta', intent='in')
    pi = xs.global_ref('pi', intent='in')
    state = xs.global_ref('state', intent='in')
    int_per_day = xs.global_ref('int_per_day', intent='in')

    def run_step(self):
        self.rate_Iy2Ih = self.pi * dta(self.eta, self.int_per_day) * self.state.loc[dict(compt='Iy')]


@xs.process
class RateIh2D:
    """Provide a `rate_Ih2D`"""
    rate_Ih2D = xs.variable(global_name='rate_Ih2D', groups=['edge_weight'], intent='out')
    mu = xs.variable(global_name='mu', intent='in')
    nu = xs.global_ref('nu', intent='in')
    state = xs.global_ref('state', intent='in')
    int_per_day = xs.global_ref('int_per_day', intent='in')

    def run_step(self):
        self.rate_Ih2D = self.nu * dta(self.mu, self.int_per_day) * self.state.loc[dict(compt='Ih')]


@xs.process
class RateIh2R:
    """Provide a `rate_Ih2R`"""
    rate_Ih2R = xs.variable(global_name='rate_Ih2R', groups=['edge_weight'], intent='out')
    gamma_Ih = xs.variable(global_name='gamma_Ih', intent='in')
    nu = xs.global_ref('nu', intent='in')
    state = xs.global_ref('state', intent='in')
    int_per_day = xs.global_ref('int_per_day', intent='in')

    def run_step(self):
        self.rate_Ih2R = (1 - self.nu) * dta(self.gamma_Ih, self.int_per_day) * self.state.loc[dict(compt='Ih')]


@xs.process
class RateIy2R:
    """Provide a `rate_Iy2R`"""
    rate_Iy2R = xs.variable(global_name='rate_Iy2R', groups=['edge_weight'], intent='out')
    gamma_Iy = xs.variable(global_name='gamma_Iy', intent='in')
    pi = xs.global_ref('pi', intent='in')
    state = xs.global_ref('state', intent='in')

    def run_step(self):
        self.rate_Iy2R = (
                self.gamma_Iy *
                self.state.loc[dict(compt='Iy')] *
                (1 - self.pi))


@xs.process
class RateIa2R:
    """Provide a `rate_Ia2R`"""
    rate_Ia2R = xs.variable(global_name='rate_Ia2R', groups=['edge_weight'], intent='out')
    gamma_Ia = xs.variable(global_name='gamma_Ia', intent='in')
    state = xs.global_ref('state', intent='in')

    def run_step(self):
        self.rate_Ia2R = self.gamma_Ia * self.state.loc[dict(compt='Ia')]


@xs.process
class SetupComptGraph:
    """Generate an 11-node compartment graph"""
    compt_graph = xs.global_ref('compt_graph', intent='out')

    def get_compt_graph(self) -> nx.DiGraph:
        g = nx.DiGraph()
        g.add_nodes_from([
            ('S', {"color": "red"}),
            ('V', {"color": "black"}),
            ('E', {"color": "black"}),
            ('Ev', {"color": "black"}),
            ('Pa', {"color": "orange"}),
            ('Py', {"color": "blue"}),
            ('Ia', {"color": "green"}),
            ('Iy', {"color": "purple"}),
            ('Ih', {"color": "yellow"}),
            ('R', {"color": "green"}),
            ('D', {"color": "blue"}),
        ])
        g.add_edges_from([
            ('S', 'E', {"priority": 0}),
            ('S', 'V', {"priority": 0}),
            ('V', 'Ev', {"priority": 1}),
            ('Ev', 'Pa', {"priority": 2}),
            ('Ev', 'Py', {"priority": 2}),
            ('E', 'Pa', {"priority": 3}),
            ('E', 'Py', {"priority": 3}),
            ('Pa', 'Ia', {"priority": 4}),
            ('Py', 'Iy', {"priority": 5}),
            ('Ia', 'R', {"priority": 6}),
            ('Iy', 'R', {"priority": 7}),
            ('Iy', 'Ih', {"priority": 7}),
            ('Ih', 'R', {"priority": 8}),
            ('Ih', 'D', {"priority": 8}),
        ])
        return g

    def vis(self, path=None):
        return visualize_compt_graph(self.compt_graph, path=path)

    def initialize(self):
        self.compt_graph = self.get_compt_graph()
        self.vis()
