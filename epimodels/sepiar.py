from episimlab.models.example_sir import SetupPhi
from episimlab.models.partition_v1 import *
from episimlab.compt_model import ComptModel
from episimlab.utils import IntPerDay
from episimlab.setup.sto import SetupStochasticFromToggle
from episimlab.setup.seed import SeedGenerator
import networkx as nx
import pandas as pd
import numpy as np

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
