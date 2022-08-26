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

# get the theory we don't fit to
ops = [np.load('/home/lthiele/pdf_baryon/measure_pdfs/ops_for_mcmc_%s_zs1.0341.npz'%s)['ops'] \
       for s in ['Dark', 'Hydro']]
kappa = np.load('/home/lthiele/pdf_baryon/measure_pdfs/ops_for_mcmc_Dark_zs1.0341.npz')['kappa']
MIN_IDX = np.argmin(np.fabs(kappa-KAPPA_MIN))
kappa = kappa[MIN_IDX:]
x_all = 2.0 * (ops[1] - ops[0]) / (ops[1] + ops[0]) # shape [subsample, kappa-bin]
x_avg = np.mean(x_all, axis=0)
target_res_not_fit = x_avg[:len(x_avg)-len(target_res)]
target_kappa_not_fit = kappa[:len(x_avg)-len(target_res)]

# get the posterior visually
reader = emcee.backends.HDFBackend(chain_fname, read_only=True)
flat_samples = reader.get_chain(flat=True)
corner_fig = corner.corner(flat_samples, labels=param_names, plot_datapoints=True)
corner_fig.savefig('posterior_kappamin%.3f.pdf'%KAPPA_MIN)

# get the best fit
nprocs = len(glob('%s/info_*.dat'%OUT_BASE))
print('Found output from %d processes'%nprocs)
min_chisq = np.inf
best_proc = None
best_idx = None
for ii in range(1, nprocs) : # process 0 is not contibuting!
    chisq = np.loadtxt('%s/info_%d.dat'%(OUT_BASE, ii), usecols=(0,))
    if np.min(chisq) < min_chisq :
        best_proc = ii
        best_idx = np.argmin(chisq)
        min_chisq = np.min(chisq)
print('Best fit chisquared = %g'%min_chisq)
theory_dmo_bf, theory_hydro_bf = np.loadtxt('%s/sample_%d_%d'%(OUT_BASE, best_proc, best_idx), unpack=True)
theory_res_bf = 2.0 * (theory_hydro_bf-theory_dmo_bf) / (theory_hydro_bf+theory_dmo_bf)
theory_kappa_edges = np.linspace(0.0, 0.2, num=33)
theory_kappa = 0.5 * (theory_kappa_edges[1:] + theory_kappa_edges[:-1])
fig_obs, ax_obs = plt.subplots()
l = ax_obs.plot(target_kappa, target_res, label='target')
ax_obs.plot(target_kappa_not_fit, target_res_not_fit, color=plt.getp(l[0], 'color'), linestyle='dashed')
ax_obs.plot(theory_kappa, theory_res_bf, label='theory best fit')
ax_obs.set_xlabel('kappa')
ax_obs.set_ylabel('2(hydro-DMO)/(hydro+DMO)')
ax_obs.legend(loc='upper left')
fig_obs.savefig('bestfit_kappamin%.3f.pdf'%KAPPA_MIN, bbox_inches='tight')
