from episimlab.models.example_sir import SetupPhi
from episimlab.models.partition_v1 import *
from episimlab.compt_model import ComptModel
from episimlab.utils import IntPerDay
from episimlab.setup.sto import SetupStochasticFromToggle
from episimlab.setup.seed import SeedGenerator
import networkx as nx
import pandas as pd
import numpy as np
import scipy
from scipy import stats
from decimal import Decimal
import math
import datetime
from copy import deepcopy
import multiprocessing as mp

######################################################################################
#  BUILD EPISIMLAB MODEL                                                             #
######################################################################################

@xs.process
class SetupComptGraphNoVis:
    """A single process in the model. Defines the directed graph `compt_graph`
    that defines the compartments and allowed transitions between them.
    """
    compt_graph = xs.global_ref('compt_graph', intent='out')

    def get_compt_graph(self) -> nx.DiGraph:
        g = nx.DiGraph()
        g.add_nodes_from([
            ('S', {"color": "red"}),
            ('E', {"color": "black"}),
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
            ('E', 'Pa', {"priority": 1}),
            ('E', 'Py', {"priority": 1}),
            ('Pa', 'Ia', {"priority": 2}),
            ('Py', 'Iy', {"priority": 3}),
            ('Ia', 'R', {"priority": 4}),
            ('Iy', 'R', {"priority": 5}),
            ('Iy', 'Ih', {"priority": 5}),
            ('Ih', 'R', {"priority": 6}),
            ('Ih', 'D', {"priority": 6}),
        ])
        return g

    def initialize(self):
        self.compt_graph = self.get_compt_graph()


@xs.process
class SetupCoordsOneVertex:
    """Initialize state coordinates. Imports compartment coordinates from the
    compartment graph.
    """
    compt = xs.index(dims=('compt'), global_name='compt_coords', groups=['coords'])
    age = xs.index(dims=('age'), global_name='age_coords', groups=['coords'])
    risk = xs.index(dims=('risk'), global_name='risk_coords', groups=['coords'])
    vertex = xs.index(dims=('vertex'), global_name='vertex_coords', groups=['coords'])
    compt_graph = xs.global_ref('compt_graph', intent='in')

    def initialize(self):
        self.compt = self.compt_graph.nodes
        self.age = ['0-4', '5-17', '18-49', '50-64', '65+']
        self.risk = ['low', 'high']
        self.vertex = ['Austin']


@xs.process
class SetupPolymodPhi(SetupPhi):
    """Use the setup class from the example but replace the random contacts with polymod contacts."""

    RANDOM_PHI_DATA = np.array([
        [2.160940833918119, 0.5973413405271149, 0.3822025191217617, 0.3523966597811896, 0.1897560907154107],
        [2.164117384279739, 8.146970087503425, 2.431391745980527, 1.885100325362032, 0.8929089276963574],
        [4.11917591841824, 5.406032123327019, 10.236283859611031, 6.706928609713956, 2.39037020719936],
        [0.8085660224830034, 0.7366175482973526, 1.70068845968991, 3.0596349022713474, 1.2034596887951288],
        [0.2808456029008987, 0.2261465786834747, 0.2102197287500719, 0.5003303235100943, 1.2322626069423352]
    ])

    phi = xs.global_ref('phi', intent='out')
    _coords = xs.group_dict('coords')


@xs.process
class SetupKappa:
    kappa = xs.variable(global_name='kappa', intent='out')

    def initialize(self):
        self.kappa = 1

    @xs.runtime(args=['step'])
    def run_step(self, step):
        if step <= 30:
            self.kappa = 1
        elif (step > 30 and step < 90):
            self.kappa = 1 - 0.95
        elif (step > 120 and step < 180):
            self.kappa = 1 - 0.95
        elif (step > 240 and step < 300):
            self.kappa = 1 - 0.95
        elif (step > 360):
            self.kappa = 1 - 0.95
        else:
            self.kappa = 1


@xs.process
class SetupBeta:
    beta = xs.variable(global_name='beta', intent='out')
    kappa = xs.global_ref('kappa', intent='in')

    def initialize(self):
        self.beta = 0.035

    def run_step(self):
        """Change beta at certain time points to simulate changing use of NPIs
        """
        self.beta = 0.035 * self.kappa


class SEPIRNoVis(EpiModel):
    """Nine-compartment SEIR model from Episimlab V1"""
    TAGS = ('SEIR', 'compartments::9', 'contact-partitioning')
    DATA_DIR = './tests/data'
    PROCESSES = {
        # Core processes
        'compt_model': ComptModel,
        'setup_sto': SetupStochasticFromToggle,
        'setup_seed': SeedGenerator,
        'setup_compt_graph': SetupComptGraphNoVis,
        'int_per_day': IntPerDay,
        'setup_coords': SetupCoordsOneVertex,
        'setup_state': SetupState,
        'setup_phi': SetupPolymodPhi,

        # Calculate greeks used by edge weight processes
        'setup_omega': SetupOmega,
        'setup_pi': SetupPiDefault,
        'setup_nu': SetupNuDefault,
        'setup_mu': mu.SetupStaticMuIh2D,
        'setup_gamma_Ih': gamma.SetupGammaIh,
        'setup_gamma_Ia': gamma.SetupGammaIa,
        'setup_gamma_Iy': gamma.SetupGammaIy,
        'setup_sigma': sigma.SetupStaticSigmaFromExposedPara,
        'setup_rho_Ia': rho.SetupRhoIa,
        'setup_rho_Iy': rho.SetupRhoIy,

        # Used for RateE2Pa and RateE2Py
        'rate_E2P': RateE2P,

        # All the expected edge weights
        'rate_S2E': RateS2E,
        'rate_E2Pa': RateE2Pa,
        'rate_E2Py': RateE2Py,
        'rate_Pa2Ia': RatePa2Ia,
        'rate_Py2Iy': RatePy2Iy,
        'rate_Ia2R': RateIa2R,
        'rate_Iy2R': RateIy2R,
        'rate_Iy2Ih': RateIy2Ih,
        'rate_Ih2R': RateIh2R,
        'rate_Ih2D': RateIh2D,
    }

    RUNNER_DEFAULTS = dict(
        clocks={
            'step': pd.date_range(start='3/11/2020', end='4/1/2020', freq='24H')
        },
        input_vars={
            'setup_sto__sto_toggle': 0,
            'setup_seed__seed_entropy': 12345,
            'rate_S2E__beta': 0.35,
            'rate_Iy2Ih__eta': 0.169492,
            'rate_E2Py__tau': 0.57,
            'rate_E2Pa__tau': 0.57,
            'setup_rho_Ia__tri_Pa2Ia': 2.3,
            'setup_rho_Iy__tri_Py2Iy': 2.3,
            'setup_sigma__tri_exposed_para': [1.9, 2.9, 3.9],
            'setup_gamma_Ih__tri_Ih2R': [9.4, 10.7, 12.8],
            'setup_gamma_Ia__tri_Iy2R_para': [3.0, 4.0, 5.0],
            'setup_mu__tri_Ih2D': [5.2, 8.1, 10.1],
        },
        output_vars={
            'compt_model__state': 'step'
        }
    )

    def plot(self, show=True):
        plot = self.out_ds['compt_model__state'].sum(['age', 'risk', 'vertex']).plot.line(x='step', aspect=2, size=9)
        if show:
            plt.show()

######################################################################################
#  SIMULATE DATA WITH EPISIMLAB                                                      #
######################################################################################

def simulate_epidemic_hospitalizations(beta):

    model_1 = SEPIRNoVis()
    input_vars = {
        'rate_S2E__beta': beta
    }
    model_result = model_1.run(input_vars=input_vars)
    hosp_data = deepcopy(model_result['compt_model__state'].sel(dict(compt='Ih')).sum(['age', 'risk', 'vertex']).to_numpy())

    return hosp_data

######################################################################################
#  BUILD MCMC FRAMEWORK                                                              #
######################################################################################

def uniform_prior(prediction):
    return 1


def hosp_error_loglik(prediction, data, sigma=1):
    """ pr(data, time=t | prediction, time = t) """

    loglik = 0
    for pred, data in zip(prediction, data):
        loglik += scipy.stats.norm.logpdf(x=data, loc=pred, scale=sigma)

    return loglik


def logit(x):
    return math.log(x / (1 - x))


def ilogit(x):
    return 1 / (1 + math.exp(-x))


def logit_sampling_distribution(x, step_size):
    """ for sampling parameters bounded by (0, 1) """

    assert x <= 1
    assert x >= 0

    logit_x = logit(x)
    new_logit_x = np.random.normal(logit_x, scale=step_size, size=1)
    new_x = ilogit(new_logit_x)

    return new_x


def posterior(prediction, data):
    prior = uniform_prior(prediction)
    likelihood = hosp_error_loglik(prediction, data)
    posterior = prior * likelihood

    return posterior


def acceptance(t0, t1):
    """t0 and t1 are log posterior probabilities"""

    acceptance_prob = Decimal(t1 - t0).exp()
    if np.random.uniform(0, 1) < acceptance_prob:
        return True
    else:
        return False


def MCMC(n_iterations, step_size, parameter_name, parameter_start, data):
    # instantiate a new model object to run
    fit_model = SEPIRNoVis()

    # update clocks assuming 1 day per data point
    fit_model.RUNNER_DEFAULTS['clocks'] = {
        'step': pd.date_range(
            start=datetime.datetime(2020, 3, 11),
            end=datetime.datetime(2020, 3, 11) + datetime.timedelta(days=(len(data) - 1)),
            freq='24H')
    }

    par_chain = np.zeros(n_iterations)
    post_chain = np.zeros(n_iterations)

    # first run:
    fit_result = fit_model.run(input_vars={parameter_name: parameter_start})
    fit_array = fit_result['compt_model__state'].sel(dict(compt='Ih')).sum(['age', 'risk', 'vertex']).to_numpy()

    try:
        assert len(fit_array) == len(data)
    except AssertionError:
        print('model output length {}'.format(len(fit_array)))
        print('data length {}'.format(len(data)))
        return

    post_chain[0] = posterior(fit_array, data)
    par_chain[0] = parameter_start

    n_accept = 0

    for i in range(1, n_iterations):

        new_param = logit_sampling_distribution(par_chain[i - 1], step_size=step_size)
        new_result = fit_model.run(input_vars={parameter_name: new_param})
        new_array = new_result['compt_model__state'].sel(dict(compt='Ih')).sum(['age', 'risk', 'vertex']).to_numpy()
        new_posterior = posterior(new_array, data)

        if acceptance(new_posterior, post_chain[i - 1]):
            post_chain[i] = new_posterior
            par_chain[i] = new_param
            n_accept += 1

        else:
            post_chain[i] = post_chain[i - 1]
            par_chain[i] = par_chain[i - 1]

    print('The final acceptance rate after {} iterations is {}'.format(n_iterations, n_accept / n_iterations))

    chains = pd.DataFrame.from_dict(
        {
            'posterior': post_chain,
            parameter_name: par_chain
        }
    )

    return chains

def main(savepath, n_chains=10):

    # MCMC inputs
    sim_data = simulate_epidemic_hospitalizations(beta=0.5)
    n_iter = 10000
    step_size = 0.005
    parameter_name = 'beta'
    parameter_start = np.random.uniform(0, 1, size=1)

    # run chains in parallel on single node
    task_list = [(n_iter, step_size, parameter_name, parameter_start, sim_data) for n in n_chains]
    pool = mp.Pool(n_chains)
    results = [pool.apply_async(MCMC, t) for t in task_list]
    pool.close()
    outputs = []
    chain_number = 0

    # combine results, add indicator column for chain number
    for r in results:
        out = r.get()
        out['chain'] = chain_number
        chain_number += 1
        outputs.append(out)

    # save outputs as a single file
    fit_data = pd.concat(outputs)
    fit_data.to_csv(savepath)

if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Example Bayesian MCMC with episimlab.')
    parser.add_argument('-o', '--outpath', help='Path to save output MCMC chain')

    opts = parser.parse_args()

    main(savepath=opts.outpath)


