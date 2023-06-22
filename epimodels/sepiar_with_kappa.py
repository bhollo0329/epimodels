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

from episimlab.models.partition_v1 import *

@xs.process
class SetupComptGraphNoVis(SetupComptGraph):
    """A single process in the model. Defines the directed graph `compt_graph`
    that defines the compartments and allowed transitions between them.
    """

    def vis(self):
        pass

#######################################################################################
#  Define the parameter vectors with seven age groups and one risk group
#  original ages: 0-4, 5-17, 18-49, 50-64, 65+
#  new ages: 0-4,
#            5-10 (copy 5-17)
#            11-13 (copy 5-17)
#            14-17, (copy 5-17)
#            18-49, 50-64, 65+
#  parameters gamma, sigma, mu and rho are handled appropriately in the import
#  from partition_v1 (these do not have age and risk-specific structure)
#######################################################################################


@xs.process
class SetupKappa:
    """Follows same file handling pattern as TravelPatFromCSV"""

    DIMS = ['vertex']
    kappa_fp = xs.variable(global_name='kappa_fp', intent='in')
    kappa = xs.variable(dims=DIMS, global_name='kappa', intent='out')
    _coords = xs.group_dict('coords')

    @property
    def dims(self):
        return self.DIMS

    @property
    def coords(self):
        return {k: v for k, v in group_dict_by_var(self._coords).items()
                if k in self.dims}

    def initialize(self):
        self.run_step(None, None)

    def get_date_mask_kappa(self, date: pd.Series, step_start, step_end) -> pd.Series:
        """Given timestamps `step_start` and `step_end`, returns a mask
        for the travel dataframe. Special handling for NaT and cases where
        `step_start` equals `step_end`.
        """
        isnull = (pd.isnull(step_start), pd.isnull(step_end))
        assert not all(isnull), f"both of `step_start` and `step_end` are null (NaT)"

        if isnull[0]:
            mask = (date == step_end)
        elif isnull[1]:
            mask = (date == step_start)
        elif step_start == step_end:
            mask = (date == step_start)
        else:
            assert step_start <= step_end
            mask = (date >= step_start) & (date < step_end)
        return mask

    @xs.runtime(args=('step_end', 'step_start'))
    def run_step(self, step_start, step_end):

        df = self.get_kappa_df()

        # Both step_start and step_end will be None for initialize
        if step_start is None and step_end is None:
            df = df[df['date'] == df['date'].min()]
        else:
            df = df[self.get_date_mask_kappa(df['date'], step_start, step_end)]

        # Validation
        if df.empty:
            raise ValueError(f'No travel data between {step_start} and {step_end}')
        logging.info(f'The date in Partition.get_travel_df is {df["date"].unique()}')

        self.kappa = self.get_kappa_da(df)

    def get_kappa_df(self):

        kappa_df = pd.read_csv(self.kappa_fp)
        assert 'vertex' in kappa_df.columns
        assert 'kappa' in kappa_df.columns
        kappa_df['date'] = pd.to_datetime(kappa_df['date'])

        return kappa_df

    def get_kappa_da(self, df):

        df = df[['vertex', 'kappa']].set_index('vertex')
        ds = xr.Dataset.from_dataframe(df)
        ds = ds.rename({'vertex': 'vertex0'})
        da = ds['kappa']

        return da

## TO DO: TESTING, MAKE SURE THAT self.beta = original_beta * self.kappa
"""
@xs.process
class SetupBeta:
    DIMS = ['vertex']
    beta_0 = xs.variable(global_name='beta_0', intent='in')
    beta = xs.variable(global_name='beta', intent='out')
    kappa = xs.global_ref('kappa', intent='in')
    _coords = xs.group_dict('coords')

    @property
    def dims(self):
        return self.DIMS

    @property
    def coords(self):
        return {k: v for k, v in group_dict_by_var(self._coords).items()
                if k in self.dims}

    def initialize(self):
        self.beta = xr.DataArray(data=self.beta_0, dims=self.dims, coords=self.coords)

    def run_step(self):
        """Change beta at certain time points to simulate changing use of NPIs
        """
        beta = self.beta_0 * self.kappa
        self.beta = xr.DataArray(data=beta, dims=self.omega_dims, coords=self.omega_coords)
'''

@xs.process
class SetupOmega:
    """Set value of omega"""
    omega = xs.global_ref('omega', intent='out')
    _coords = xs.group_dict('coords')

    @property
    def coords(self):
        return group_dict_by_var(self._coords)

    @property
    def omega_dims(self):
        return get_var_dims(RateS2E, 'omega')

    @property
    def omega_coords(self):
        return {dim: self.coords[dim.rstrip('01')] for dim in self.omega_dims}

    def initialize(self):
        da = xr.DataArray(data=0., dims=self.omega_dims, coords=self.omega_coords)

        # One could set specific values per compartment here
        da.loc[dict(compt='Ia')] = 0.666666667
        da.loc[dict(compt='Iy')] = 1.
        da.loc[dict(compt='Pa')] = [0.91117513, 0.91117513, 0.91117513, 0.91117513, 0.92460653, 0.95798887, 0.98451149]
        da.loc[dict(compt='Py')] = [1.36676269, 1.36676269, 1.36676269, 1.36676269, 1.3869098, 1.43698331, 1.47676724]

        self.omega = da


@xs.process
class SetupNuDefault:
    """Provide a default value for nu"""
    DIMS = ['age']
    nu = xs.variable(dims=DIMS, global_name='nu', intent='out')
    _coords = xs.group_dict('coords')

    @property
    def dims(self):
        return self.DIMS

    @property
    def coords(self):
        return {k: v for k, v in group_dict_by_var(self._coords).items()
                if k in self.dims}

    def initialize(self):
        self.nu = xr.DataArray(
            [0.02878229, 0.09120554, 0.09120554, 0.09120554, 0.02241002, 0.07886779, 0.17651128],
            dims=self.dims, coords=self.coords)


@xs.process
class SetupPiDefault:
    """Provide a default value for pi"""
    DIMS = ('risk', 'age')
    pi = xs.variable(dims=DIMS, global_name='pi', intent='out')
    _coords = xs.group_dict('coords')

    @property
    def dims(self):
        return self.DIMS

    @property
    def coords(self):
        return {k: v for k, v in group_dict_by_var(self._coords).items()
                if k in self.dims}

    def initialize(self):
        self.pi = xr.DataArray(np.array([
            [5.92915812e-04, 4.55900959e-04, 4.55900959e-04, 4.55900959e-04, 2.78247788e-02, 5.95202276e-02, 7.03344654e-02],
            [5.91898663e-03, 4.55299354e-03, 4.55299354e-03, 4.55299354e-03, 2.57483139e-01, 5.07631836e-01, 5.84245731e-01]]),
            dims=self.dims, coords=self.coords)

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
from episimlab.partition.travel_pat import TravelPatFromCSV  # contains school calendar info (no travel to schools when school is closed)
import pandas as pd

@xs.process
class Intervention(Partition):
    # we can add in an intervention process here if we want
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
class SetupCoordsOneVertex:
    """Initialize state coordinates. Imports compartment coordinates from the
    compartment graph.
    """
    travel_pat = xs.global_ref('travel_pat', intent='in')
    compt = xs.index(dims=('compt'), global_name='compt_coords', groups=['coords'])
    age = xs.index(dims=('age'), global_name='age_coords', groups=['coords'])
    risk = xs.index(dims=('risk'), global_name='risk_coords', groups=['coords'])
    vertex = xs.index(dims=('vertex'), global_name='vertex_coords', groups=['coords'])
    compt_graph = xs.global_ref('compt_graph', intent='in')

    def initialize(self):
        self.compt = self.compt_graph.nodes
        self.age = self.travel_pat.coords['age0'].values
        self.risk = ['low', 'high']
        self.vertex = self.travel_pat.coords['vertex0'].values  # intentionally excludes schools


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
        'compt_model': ComptModel,
        'setup_sto': SetupStochasticFromToggle,
        'setup_seed': SeedGenerator,
        'setup_compt_graph': SetupComptGraphNoVis,
        'int_per_day': IntPerDay,
        'setup_coords': SetupCoordsOneVertex,
        'setup_state': SetupStateWithRiskFromCSV,

        # Contact partitioning
        'setup_travel': TravelPatFromCSV,
        'setup_contacts': ContactsFromCSV,
        'partition': Partition,

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
        'setup_kappa': SetupKappa,
        'setup_beta': SetupBeta,

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
            'setup_beta__beta_0': 0.35,
            'rate_Iy2Ih__eta': 0.169492,
            'rate_E2Py__tau': 0.57,
            'rate_E2Pa__tau': 0.57,
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
            'kappa_fp': None
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
