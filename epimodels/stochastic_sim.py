from epimodels.sepiar import SEPIRNoVis
from epimodels.sepiar_seven_ages import SEPIRSevenAgesNoVis
import xarray as xr
import numpy as np
import multiprocessing as mp

######################################################################################
#  SIMULATE DATA WITH EPISIMLAB                                                      #
######################################################################################


def simulate_epidemic(seed_entropy):

    sim = SEPIRNoVis()
    input_vars = {
        'setup_seed__seed_entropy': seed_entropy
    }
    sim_result = sim.run(input_vars=input_vars)

    return sim_result

def simulate_granular_epidemic(seed_entropy):

    sim = SEPIRSevenAgesNoVis()
    input_vars = {
        'setup_seed__seed_entropy': seed_entropy
    }
    sim_result = sim.run(input_vars=input_vars)

    return sim_result


def entropy_generator(n_sims):

    seeds = np.random.SeedSequence(entropy=(8888))
    entropy_list = seeds.generate_state(n_sims)

    return entropy_list


def main(savepath, simfunc, n_parallel=2, n_sims=10):

    sim_entropy = entropy_generator(n_sims)

    pool = mp.Pool(n_parallel)
    results = [pool.apply_async(simfunc, (s, )) for s in sim_entropy]
    pool.close()
    outputs = []
    sim_number = 0

    # combine results, add indicator column for chain number
    for r in results:
        out = r.get()
        out = out.assign_coords({'index': sim_number})
        sim_number += 1
        outputs.append(out)

    # save outputs as a single file
    all_sims = xr.concat(outputs, dim='index')
    all_sims.to_zarr(store=savepath)

if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Example simulation.')
    parser.add_argument('-f', '--simfunc', help='Name of model function to simulate, one of "simulate_epidemic" or "simulate_granular_epidemic"')
    parser.add_argument('-o', '--outpath', help='Path to save simulation output')

    opts = parser.parse_args()

    main(savepath=opts.outpath)
