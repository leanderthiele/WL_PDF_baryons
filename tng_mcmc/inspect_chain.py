# Command line arguments:
# [1] kappa_min

from sys import argv
from glob import glob

import numpy as np
from matplotlib import pyplot as plt
import emcee
import corner
import h5py

KAPPA_MIN = float(argv[1])
OUT_BASE = '/scratch/gpfs/lthiele/BCM_MCMC_TNG_chains_kappamin%.3f'%KAPPA_MIN

chain_fname = '%s/chain.hdf5'%OUT_BASE

with h5py.File(chain_fname, 'r') as f :
    param_names = f['Header'].attrs['parameters']
    target_kappa = f['Header/kappa'][...]
    target_res = f['Header/target_pdf'][...]
    print('Accepted:')
    print(f['mcmc/accepted'][...])

# get the posterior visually
reader = emcee.backends.HDFBackend(chain_fname, read_only=True)
flat_samples = reader.get_chain(flat=True)
corner_fig = corner.corner(samples, labels=param_names, plot_datapoints=True)
corner.savefig('posterior_kappamin%.3f.pdf'%KAPPA_MIN)

# get the best fit
nprocs = len(glob('%s/info_*.dat'%OUT_BASE))
min_chisq = np.inf
best_proc = None
best_idx = None
for ii in range(1, procs) : # process 0 is not contibuting!
    chisq = np.loadtxt('info_%d.dat'%ii, usecols=(0,))
    if np.min(chisq) < min_chisq :
        best_proc = ii
        best_idx = np.argmin(chisq)
        min_chisq = np.min(chisq)
theory_dmo_bf, theory_hydro_bf = np.loadtxt('sample_%d_%d'%(best_proc, best_idx))
theory_res_bf = 2.0 * (theory_hydro_bf-theory_dmo_bf) / (theory_hydro_bf+theory_dmo_bf)
theory_kappa_edges = np.linspace(0.0, 0.2, num=33)
theory_kappa = 0.5 * (theory_kappa_edges[1:] + theory_kappa_edges[:-1])
fig_obs, ax_obs = plt.subplots()
ax_obs.plot(target_kappa, target_res, label='target')
ax_obs.plot(theory_kappa, theory_res_bf, label='theory best fit')
ax_obs.set_xlabel('kappa')
ax_obs.set_ylabel('2(hydro-DMO)/(hydro+DMO)')
fig_obs.savefig('bestfit_kappamin%.3f.pdf'%KAPPA_MIN, bbox_inches='tight')
