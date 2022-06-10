######################################################################################
#  BUILD EPISIMLAB MODEL
######################################################################################

#######################################################################################
#  Import the system of equations
#       The system of equations for the SEPIAR model is defined in module partion_v1.
#       Rather than re-define those equations here, we import their definitions.
#       The SetupComptGraph class defines the order of the equations, and contains
#       a method that draws the graph for every run. We don't want to create an image
#       with every run of a large simulation, so we extend this class to nullify
#       the function that makes the graph.
#######################################################################################

from .vaccine_seven_ages import *

@xs.process
class SetupComptGraphNoVis(SetupComptGraph):
    """A single process in the model. Defines the directed graph `compt_graph`
    that defines the compartments and allowed transitions between them.
    """

    def vis(self):
        pass

#######################################################################################
#  Import the compartment model engine.
#       This module contains the low-level function that manage running the model
#       and organizing the model outputs
#######################################################################################

from episimlab.compt_model import ComptModel

#######################################################################################
#  Import additional modules to help with setup
#   - module IntPerDay contains methods to discretize rates
#   - module SeedGenerator contains methods to manage random seeds and ensure that
#     stochastic simulations are unique
#   - module SetupStochasticFromToggle processes the boolean to determine if random
#     states will be used
#######################################################################################

from episimlab.utils import IntPerDay
from episimlab.setup.sto import SetupStochasticFromToggle
from episimlab.setup.seed import SeedGenerator
from episimlab.setup.state import SetupStateWithRiskFromCSV

#######################################################################################
#  Setup partitioned contact matrix (based on travel patterns and baseline contacts)
#       This is handled by classes in the module episimlab.partition, without
#       modification. If adding axes other than age and risk to the contact matrix, you
#       need to extend class Partition and modify property phi_dims, as well as modify
#       the force of infection calculation (episimlab.foi.BaseFOI) to include these
#       dimensions.
#######################################################################################

from episimlab.partition.contacts import ContactsFromNetCDF
from episimlab.partition.partition import Partition
from episimlab.partition.behavior_change import BehaviorChange  # contains school calendar info (no travel to schools when school is closed)
import pandas as pd

@xs.process
class Intervention(Partition):
    # we can add in an intervetion process here if we want
    pass

######################################################################################
#  Define the dimensions of the array to hold model results
#     - number and names of age groups
#     - number and names of risk groups
#     - number and names of vertices (e.g. ZCTAs)
#     - number and names of compartments
#  Some of these attributes, like number and names of age groups and vertices,
#  are defined by external data structures containing stratified contact information
######################################################################################

@xs.process
class SetupCoords:
    """Initialize state coordinates. Imports compartment coordinates from the
    compartment graph.
    """
    #travel_pat = xs.global_ref('travel_pat', intent='in')
    census_fp = xs.variable(intent='in')
    compt = xs.index(dims=('compt'), global_name='compt_coords', groups=['coords'])
    age = xs.index(dims=('age'), global_name='age_coords', groups=['coords'])
    risk = xs.index(dims=('risk'), global_name='risk_coords', groups=['coords'])
    vertex = xs.index(dims=('vertex'), global_name='vertex_coords', groups=['coords'])
    compt_graph = xs.global_ref('compt_graph', intent='in')

    def read_census_csv(self):
        df = pd.read_csv(
            self.census_fp, dtype={'GEOID': str}
        ).drop('Unnamed: 0', axis=1).rename(columns={'GEOID': 'vertex', 'age_bin': 'age'})
        return df

    def initialize(self):
        df = self.read_census_csv()
        self.compt = self.compt_graph.nodes
        self.age = df['age'].unique()
        self.risk = df['risk'].unique()
        self.vertex = df['vertex'].unique()  # intentionally excludes schools


#######################################################################################
#  Combine all the pieces into a model. This must be the final step, after all other
#  components are defined.
#######################################################################################

class SEPIRSevenAgesNoVis(EpiModel):
    """Nine-compartment SEIR model from Episimlab V2"""
    TAGS = ('SEIR', 'compartments::9', 'contact-partitioning-seven-ages')
    DATA_DIR = './tests/data'
    PROCESSES = {
        # Core processes
        'setup_compt_graph': SetupComptGraphNoVis,
        'compt_model': ComptModel,
        'int_per_day': IntPerDay,
        'setup_coords': SetupCoords,
        'setup_state': SetupStateWithRiskFromCSV,
        'setup_sto': SetupStochasticFromToggle,
        'setup_seed': SeedGenerator,

        # Contact partitioning
        'setup_travel': BehaviorChange,
        'setup_contacts': ContactsFromCSV,
        'partition': Partition,

        # Calculate greeks used by edge weight processes
        'beta_reduction': BetaReduction,
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

        # calculate vaccine doses
        'setup_doses': SetupVaccineDoses,

        # used for RateE2Pa and RateE2Py
        'rate_E2P': RateE2P,
        'rate_Ev2P': RateEv2P,

        # all the expected edge weights
        'rate_S2E': RateS2E,
        'rate_S2V': RateS2V,
        'rate_V2Ev': RateV2Ev,
        'rate_E2Pa': RateE2Pa,
        'rate_Ev2Pa': RateEv2Pa,
        'rate_E2Py': RateE2Py,
        'rate_Ev2Py': RateEv2Py,
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
            'beta_reduction': 1.0,  # no beta reduction
            'rate_S2V__eff_vaccine': 0.8,
            'rate_Iy2Ih__eta': 0.169492,
            'rate_E2Py__tau': 0.57,
            'rate_E2Pa__tau': 0.57,
            'rate_Ev2Py__tau_v': 0.055,
            'rate_Ev2Pa__tau_v': 0.055,
            'setup_rho_Ia__tri_Pa2Ia': 2.3,
            'setup_rho_Iy__tri_Py2Iy': 2.3,
            'setup_sigma__tri_exposed_para': [1.9, 2.9, 3.9],
            'setup_gamma_Ih__tri_Ih2R': [9.4, 10.7, 12.8],
            'setup_gamma_Ia__tri_Iy2R_para': [3.0, 4.0, 5.0],
            'setup_mu__tri_Ih2D': [5.2, 8.1, 10.1],
            # these file paths must be added at runtime
            'travel_pat_fp': None,
            'contacts_fp': None,
            'census_fp': None,
            'hosp_catchment_fp': None
        },
        output_vars={
            'compt_model__state': 'step',
            'compt_model__tm': 'step'
        }
    )

    def plot(self, show=True):
        plot = self.out_ds['compt_model__state'].sum(['age', 'risk', 'vertex']).plot.line(x='step', aspect=2, size=9)
        if show:
            plt.show()
