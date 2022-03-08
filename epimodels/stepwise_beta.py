from epimodels.sepiar import *

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


class SEPIRNoVisStepwise(SEPIRNoVis):
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
            'rate_S2E__beta': 0.035,
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
