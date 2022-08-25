""" Command line arguments:
[1] number of walkers (I think should be twice the number of MPI processes)
[2] minimum kappa (float)
"""

import os
import os.path
import sys
from sys import argv
import subprocess
from time import time

import numpy as np
import h5py

import emcee
from schwimmbad import MPIPool

NWALKERS = int(argv[1])
KAPPA_MIN = float(argv[2])
WRITE_PERIOD = 10 # after this many likelihood evaluations the chisq is written to file

OUT_BASE = '/scratch/gpfs/lthiele/BCM_MCMC_TNG_chains_kappamin%.3f'%KAPPA_MIN
os.makedirs(OUT_BASE, exist_ok=True)

RANK = int(os.environ['SLURM_PROCID'])
WORLD_SIZE = int(os.environ['SLURM_NTASKS'])

# load the target, shape [[Dark, Hydro], subsample, kappa-bin]
ops = [np.load('/home/lthiele/pdf_baryon/measure_pdfs/ops_for_mcmc_%s_zs1.0341.npz'%s)['ops'] \
       for s in ['Dark', 'Hydro']]
kappa = np.load('/home/lthiele/pdf_baryon/measure_pdfs/ops_for_mcmc_Dark_zs1.0341.npz')['kappa']
MIN_IDX = np.argmin(np.fabs(kappa-KAPPA_MIN))
kappa = kappa[MIN_IDX:]
x_all = 2.0 * (ops[1] - ops[0]) / (ops[1] + ops[0]) # shape [subsample, kappa-bin]
# implement kappa cut
x_all = x_all[:, MIN_IDX:]
x_avg = np.mean(x_all, axis=0)
x_cov = np.cov(x_all, rowvar=False) / x_all.shape[0]
x_covinv = np.linalg.inv(x_cov)

# the parameters
param_names = ['log_M_c', 'log_M_1_z0_cen', 'eta', 'beta',
               'theta_inn', 'theta_out', 'log_M_inn', 'log_M_r', ]

# fiducial values, from bispectrum fit
theta_0 = [11.672356390589902, 12.998276528979659, 0.7950988565203113, 3.7422487294724927,
           0.11559464687549903, 2.2228694416908383, 9.013004910858838, 16.0, ]

# prior ranges, we'll also use these to initialize the chains as we don't really have much information
theta_priors = [(9.0, 13.0), (9.0, 15.0), (0.1, 3.0), (0.5, 5.0),
                (0.05, 0.5), (0.5, 3.0), (8.0, 15.0), (13.0, 18.0), ]

def log_prior (theta) :
    for t, p in zip(theta, theta_priors) :
        if not p[0] < t < p[1] :
            return -np.inf
    return 0.0


NUM_STEPS = 0 # changes
txt_fname = '%s/info_%d.dat'%(OUT_BASE, RANK)
chisq_arr = np.full(WRITE_PERIOD, -1.0)
time_arr = np.full(WRITE_PERIOD, -1.0)

def log_likelihood (theta) :

    global NUM_STEPS
    global chisq_arr
    global time_arr

    fname = '%s/sample_%d_%d'%(OUT_BASE, RANK, NUM_STEPS)

    binary = './run_hmpdf'

    args = [binary, fname, *map(lambda x: '%g'%x, theta)]

    start = time()
    returned = subprocess.run(args)
    end = time()
    runtime = end - start # in seconds

    if returned.returncode != 0 :
        raise RuntimeError

    theory_dmo, theory_hydro = np.loadtxt(fname, unpack=True)
    x_theory = 2.0 * (theory_hydro-theory_dmo) / (theory_hydro+theory_dmo)
    x_theory = x_theory[MIN_IDX:]

    delta_x = x_theory - x_avg
    chisq = np.einsum('i,ij,j->', delta_x, x_covinv, delta_x)

    idx_in_arrs = NUM_STEPS % WRITE_PERIOD
    chisq_arr[idx_in_arrs] = chisq
    time_arr[idx_in_arrs] = runtime

    NUM_STEPS += 1
    if NUM_STEPS % WRITE_PERIOD == 0 :
        with open(txt_fname, 'a') as f :
            np.savetxt(f, np.stack([chisq_arr, time_arr], -1), header='chisq time[sec]')

    return -0.5 * chisq


def log_probability (theta) :
    
    lp = log_prior(theta)
    if not np.isfinite(lp) :
        return lp
    ll = log_likelihood(theta)
    return lp + ll


def get_initial () :
    rng = np.random.default_rng()
    return np.array([[rng.uniform(*p) for p in theta_priors] for _ in range(NWALKERS)])

with MPIPool() as pool :
    
    if not pool.is_master() :
        pool.wait()
        sys.exit(0)

    initial = get_initial()
    ndim = initial.shape[1]

    chain_fname = '%s/chain.hdf5'%OUT_BASE

    h5_exists = os.path.isfile(chain_fname)

    if not h5_exists :
        with h5py.File(chain_fname, 'w') as f :
            header = f.create_group('Header')
            header.attrs.create('parameters', param_names)
            header.attrs.create('kappa_min', KAPPA_MIN)
            header.create_dataset('kappa', data=kappa)
            header.create_dataset('target_pdf', data=x_avg)
            header.create_dataset('cov_pdf', data=x_cov)

    backend = emcee.backends.HDFBackend(chain_fname)
    if not h5_exists :
        backend.reset(NWALKERS, ndim)

    sampler = emcee.EnsembleSampler(NWALKERS, ndim, log_probability,
                                    pool=pool, backend=backend)

    sampler.run_mcmc(None if h5_exists else initial, 100000)
