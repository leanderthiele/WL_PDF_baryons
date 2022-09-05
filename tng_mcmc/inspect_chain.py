""" Command line arguments:
    [1] baryon mode
    [2] sim (TNG, BAHAMAS_fid, BAHAMAS_lo, BAHAMAS_hi)
    [3] kappa_min
"""

from sys import argv
from glob import glob

import numpy as np
from matplotlib import pyplot as plt
import emcee
import corner
import h5py

# burn in
DISCARD = 100

try :
    BARYON_MODE = argv[1].upper()
    SIM = argv[2]
    KAPPA_MIN = float(argv[3])
except Exception :
    print(__doc__)
    raise
OUT_BASE = '/scratch/gpfs/lthiele/%s_MCMC_%s_chains_kappamin%.3f'%(BARYON_MODE, SIM, KAPPA_MIN)

chain_fname = '%s/chain.hdf5'%OUT_BASE

with h5py.File(chain_fname, 'r') as f :
    param_names = f['Header'].attrs['parameters']
    target_kappa = f['Header/kappa'][...]
    target_res = f['Header/target_pdf'][...]
    accepted = f['mcmc/accepted'][...]
    print('Accepted:')
    print(accepted)
    Naccepted = np.sum(accepted)

# get the theory we don't fit to
if SIM == 'TNG' :
    ops = [np.load('/home/lthiele/pdf_baryon/measure_pdfs/ops_for_mcmc_%s_zs1.0341.npz'%s)['ops'] \
           for s in ['Dark', 'Hydro']]
    kappa = np.load('/home/lthiele/pdf_baryon/measure_pdfs/ops_for_mcmc_Dark_zs1.0341.npz')['kappa']
else :
    ops = [np.load('/home/lthiele/pdf_baryon/BAHAMAS/ops_for_mcmc_%s_zs1.0000.npz'%s)['ops'] \
           for s in ['DMO', 'hydro_'+SIM.split('_')[-1]]]
    kappa = np.load('/home/lthiele/pdf_baryon/BAHAMAS/ops_for_mcmc_DMO_zs1.0000.npz')['kappa']
x_all = 2.0 * (ops[1] - ops[0]) / (ops[1] + ops[0]) # shape [subsample, kappa-bin]
x_avg = np.mean(x_all, axis=0)
target_res_not_fit = x_avg[:len(x_avg)-len(target_res)+1]
target_kappa_not_fit = kappa[:len(x_avg)-len(target_res)+1]

# get the posterior visually
reader = emcee.backends.HDFBackend(chain_fname, read_only=True)
flat_samples = reader.get_chain(discard=DISCARD, flat=True)
flat_all_samples = reader.get_chain(flat=True)
Nsamples = flat_all_samples.shape[0]
print('Nsamples: %d, acceptance rate: %g'%(Nsamples, Naccepted/Nsamples))
corner_fig = corner.corner(flat_samples, labels=param_names, plot_datapoints=True)
corner_fig.savefig('posterior_%s_%s_kappamin%.3f.pdf'%(BARYON_MODE, SIM, KAPPA_MIN))

# get the best fit
nprocs = len(glob('%s/samples_*.hdf5'%OUT_BASE))
print('Found output from %d processes'%nprocs)
min_chisq = np.inf
best_proc = None
best_idx = None
for ii in range(1, nprocs) : # process 0 is not contibuting!
    with h5py.File('%s/samples_%d.hdf5'%(OUT_BASE, ii), 'r') as f :
        chisq = np.array([f[d].attrs['chisq'] for d in sorted(f.keys(), key=lambda s: s.split('_')[1])])
    if np.min(chisq) < min_chisq :
        best_proc = ii
        best_idx = np.argmin(chisq)
        min_chisq = np.min(chisq)
print('Best fit chisquared = %g'%min_chisq)
with h5py.File('%s/samples_%d.hdf5'%(OUT_BASE, best_proc), 'r') as f :
    d = f['sample_%d'%best_idx]
    theory_dmo_bf, theory_hydro_bf = d[...]
    theta_bf = d.attrs['theta']
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
ax_obs.set_xlim(0,0.2)
fig_obs.savefig('bestfit_%s_%s_kappamin%.3f.pdf'%(BARYON_MODE, SIM, KAPPA_MIN), bbox_inches='tight')
