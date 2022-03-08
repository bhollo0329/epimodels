from epimodels.stepwise_beta import SEPIRNoVisStepwise
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
#  SIMULATE DATA WITH EPISIMLAB                                                      #
######################################################################################

def simulate_epidemic_hospitalizations(beta):

    sim = SEPIRNoVisStepwise()
    input_vars = {
        'setup_seed__seed_entropy': 12345
    }
    sim_result = sim.run(input_vars=input_vars)

    return sim

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
    fit_model = SEPIRNoVisStepwise()

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
    parameter_start = np.random.uniform(0, 1, size=1).item()

    # run chains in parallel on single node
    task_list = [(n_iter, step_size, parameter_name, parameter_start, sim_data) for _ in range(n_chains)]
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


