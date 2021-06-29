import numpy as np
from matplotlib import pyplot as plt
from os import system
from sys import argv

from scipy.interpolate import interp1d

RECOMPUTE = bool(int(argv[1])) if len(argv)>1 else True

zs_indices = [12, 22, 34]
zs = np.loadtxt('kappaTNG_PDF/zs.dat')[zs_indices]

sim_DMO   = np.load('kappaTNG_PDF/DMO/PDF_smooth1arcmin.npy')[zs_indices,...]
sim_hydro = np.load('kappaTNG_PDF/hydro/PDF_smooth1arcmin.npy')[zs_indices,...]

if True : # enforce <kappa> = 0 exactly
    for ii in range(len(zs_indices)) :
        kavg_DMO = np.sum(sim_DMO[ii,:,0]*sim_DMO[ii,:,1]) / np.sum(sim_DMO[ii,:,1])
        kavg_hydro = np.sum(sim_hydro[ii,:,0]*sim_hydro[ii,:,1]) / np.sum(sim_hydro[ii,:,1])

        kprime_DMO = sim_DMO[ii,:,0] - kavg_DMO
        kprime_hydro = sim_hydro[ii,:,0] - kavg_hydro

        i_DMO = interp1d(kprime_DMO, sim_DMO[ii,:,1],
                         kind='cubic', bounds_error=False, fill_value=0)
        i_hydro = interp1d(kprime_hydro, sim_hydro[ii,:,1],
                           kind='cubic', bounds_error=False, fill_value=0)

        sim_DMO[ii,:,1] = i_DMO(sim_DMO[ii,:,0])
        sim_hydro[ii,:,1] = i_hydro(sim_hydro[ii,:,0])

theory_DMO   = np.empty((sim_DMO.shape[0], sim_DMO.shape[1]))
theory_hydro = np.empty((sim_DMO.shape[0], sim_DMO.shape[1]))

for ii, zsource in enumerate(zs) :
    edg_file = 'hmpdf_ws/edges_zs%.4f.dat'%zsource

    tmp_file_DMO   = 'hmpdf_ws/op_DMO_zs%.4f.dat'%zsource
    tmp_file_hydro = 'hmpdf_ws/op_hydro_zs%.4f.dat'%zsource

    if RECOMPUTE :
        for hydro in [0, 1] :
            print('---------------' +
                  ('HYDRO' if hydro else 'DMO')
                  + '---------------')
            system('./run_hmpdf %.4f %s %s %d'%(zsource, edg_file,
                                                tmp_file_hydro if hydro else tmp_file_DMO,
                                                hydro))

    theory_hydro[ii, :] = np.loadtxt(tmp_file_hydro)
    theory_DMO[ii, :]   = np.loadtxt(tmp_file_DMO)

fig_res, ax_res = plt.subplots()

fig_log, ax_log = plt.subplots(ncols=len(zs_indices))
ax_log = ax_log.flatten()

for ii, zsource in enumerate(zs) :
    lsim = ax_res.plot(sim_DMO[ii, :, 0],
                       (sim_hydro[ii, :, 1]-sim_DMO[ii, :, 1])/sim_DMO[ii, :, 1],
                       linestyle='solid', label='$z_s = %.2f$'%zsource)
    ax_res.plot(sim_DMO[ii, :, 0], (theory_hydro[ii, :]-theory_DMO[ii, :])/theory_DMO[ii, :],
                linestyle='dashed', color=plt.getp(lsim[0], 'color'))

    for hydro, pair in enumerate(zip([sim_DMO, sim_hydro], [theory_DMO, theory_hydro])) :
        sim, theory = pair
        lsim = ax_log[ii].semilogy(sim[ii, :, 0],
                                   sim[ii, :, 1] / np.sum(sim[ii, :, 1]),
                                   linestyle='solid',
                                   label='hydro' if hydro else 'DMO')
        ax_log[ii].semilogy(sim[ii, :, 0], theory[ii, :],
                            linestyle='dashed', color=plt.getp(lsim[0], 'color'))
        ax_log[ii].set_title('$z_s = %.2f$'%zsource)
    ax_log[ii].set_xlim(0.0, 0.5)
    ax_log[ii].set_ylim(1e-8, 1e-1)

ax_res.set_xlim(0.0, 0.5)
ax_res.set_ylim(-0.2, 0.2)
ax_res.set_xlabel('convergence $\kappa$')
ax_res.set_ylabel('$p_{\sf hydro} / p_{\sf DMO} - 1$')
ax_res.legend()
fig_res.savefig('hydro_vs_DMO.pdf', bbox_inches='tight')
fig_log.savefig('smoothing_calib.pdf', bbox_inches='tight')
