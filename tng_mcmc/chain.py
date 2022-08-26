""" Command line arguments:
[1] number of walkers (I think should be twice the number of MPI processes)
[2] minimum kappa (float)
[3] baryon mode (one of BCM, TOT_CONC, BAR_CONC)
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
BARYON_MODE = argv[3]

# do not change the order!
baryon_modes = ['BCM', 'TOT_CONC', 'BAR_CONC', ]
BARYON_MODE_IDX = baryon_modes.index(BARYON_MODE)

OUT_BASE = '/scratch/gpfs/lthiele/%s_MCMC_TNG_chains_kappamin%.3f'%(BARYON_MODE, KAPPA_MIN)
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
if BARYON_MODE == 'BCM' :
    param_names = ['log_M_c', 'log_M_1_z0_cen', 'eta', 'beta',
                   'theta_inn', 'theta_out', 'log_M_inn', 'log_M_r', ]
elif BARYON_MODE in ['TOT_CONC', 'BAR_CONC', ] :
    param_names = ['A', 'B', 'C',
                   'd', 'e', 'f', ]
else :
    raise NotImplementedError

if BARYON_MODE == 'BCM' :
    # fiducial values, from bispectrum fit
    theta_0 = [11.672356390589902, 12.998276528979659, 0.7950988565203113, 3.7422487294724927,
               0.11559464687549903, 2.2228694416908383, 9.013004910858838, 16.0, ]
# the following come from profile fits
elif BARYON_MODE == 'TOT_CONC' :
    theta_0 = [1.375023001e+01, -8.688706e-02, -8.4451454e-01,
               2.705265e-02, 5.2513069e-01, 1.12049399e+00, ]
elif BARYON_MODE == 'BAR_CONC' :
    # note theta_0[4] is outside of prior here, not a problem...
    theta_0 = [4.57478936e+00, -1.67755068e+00, -2.84129704e+00,
               4.90993037e-04, -1.16223295e+01, -6.83565016e+00, ]
    

# prior ranges, we'll also use these to initialize the chains as we don't really have much information
if BARYON_MODE == 'BCM' :
    theta_priors = [(9.0, 13.0), (9.0, 15.0), (0.1, 3.0), (0.5, 5.0),
                    (0.05, 0.5), (0.5, 3.0), (8.0, 15.0), (13.0, 18.0), ]
elif BARYON_MODE == 'TOT_CONC' :
    theta_priors = [(8.0, 16.0), (-0.2, 0.0), (-3.0, 0.0),
                    (0.0, 0.1), (0.0, 1.0), (-10.0, 10.0), ]
elif BARYON_MODE == 'BAR_CONC' :
    theta_priors = [(0.1, 16.0), (-5.0, 0.0), (-3.0, 0.0),
                    (0.0, 0.1), (-20.0, 0.0), (-10.0, 10.0), ]

def log_prior (theta) :
    for t, p in zip(theta, theta_priors) :
        if not p[0] < t < p[1] :
            return -np.inf
    return 0.0


NUM_STEPS = 0 # changes
samples_fname = '%s/samples_%.hdf5'%(OUT_BASE, RANK)

# figure out if previous samples have been written already
if not os.path.isfile(samples_fname) :
    NUM_EXISTING = 0
else :
    with h5py.File(samples_fname, 'r') as f :
        dsets = list(f.keys())
    dsets = list(filter(lambda s: s.startswith('sample_'), dsets))
    NUM_EXISTING = len(dsets)

def log_likelihood (theta) :

    global NUM_STEPS
    global chisq_arr
    global time_arr

    fname = '%s/sample_%d_%d'%(OUT_BASE, RANK, NUM_STEPS)

    binary = './run_hmpdf'

    args = [binary, fname, str(BARYON_MODE_IDX), *map(lambda x: '%g'%x, theta)]

    start = time()
    returned = subprocess.run(args)
    end = time()
    runtime = end - start # in seconds

    if returned.returncode != 0 :
        raise RuntimeError

    theory_dmo, theory_hydro = np.loadtxt(fname, unpack=True)

    os.remove(fname) # don't crash the file system

    x_theory = 2.0 * (theory_hydro-theory_dmo) / (theory_hydro+theory_dmo)
    x_theory = x_theory[MIN_IDX:]

    delta_x = x_theory - x_avg
    chisq = np.einsum('i,ij,j->', delta_x, x_covinv, delta_x)

    # for the record
    with h5py.File(samples_fname, 'a') as f :
        dset = f.create_dataset('sample_%d'%(NUM_STEPS+NUM_EXISTING),
                                data=np.stack([theory_dmo, theory_hydro], -1))
        dset.attrs.create('chisq', chisq)
        dset.attrs.create('runtime', runtime)
        dset.attrs.create('theta', theta)

    NUM_STEPS += 1

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
