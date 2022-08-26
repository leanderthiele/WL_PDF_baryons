import numpy as np
from matplotlib import pyplot as plt
from os import system
from sys import argv

from scipy.interpolate import interp1d

SIM = 'BAHAMAShires_fid'
OBS = 'PK' # only relevant for Arico20 BCM

RECOMPUTE = bool(int(argv[1])) if len(argv)>1 else True

# TODO
if SIM == 'TNG' :
    zs_indices = [12, 22, 34]
elif SIM.startswith('BAHAMAShires') :
    zs_indices = [0, ]
elif SIM.startswith('BAHAMAS_') :
    zs_indices = [0, 1, 2]
else :
    raise RuntimeError

zs = np.atleast_1d(np.loadtxt('./%s_PDF/zs.dat'%SIM))[zs_indices]

# don't change
# NOTE: anything before underscore is assumed to uniquely identify the DMO version,
#       anything after underscore identifies the hydro model
sims = ['TNG',
        'BAHAMAS_fid', 'BAHAMAS_lo', 'BAHAMAS_hi',
        'BAHAMAShires_fid', 'BAHAMAShires_lo', 'BAHAMAShires_hi',
       ]
obss = ['PK', 'QK', 'PK_QK', ]

#sim_DMO   = np.load('kappaTNG_PDF/DMO/PDF_smooth1arcmin.npy')[zs_indices,...]
#sim_hydro = np.load('kappaTNG_PDF/hydro/PDF_smooth1arcmin.npy')[zs_indices,...]
# TODO
if SIM == 'TNG' :
    sim_BCM = [np.loadtxt('%s_PDF/BCM/op_BCM_zs%.4f.dat'%(SIM, zs[ii]))
               for ii in range(len(zs_indices))]
else :
    # for BAHAMAS, we don't have a BCM realization available I believe
    sim_BCM = None

sim_DMO   = [np.loadtxt('%s_PDF/Dark/op_Dark_zs%.4f.dat'%(SIM, zs[ii]))
             for ii in range(len(zs_indices))]
sim_hydro = [np.loadtxt('%s_PDF/Hydro/op_Hydro_zs%.4f.dat'%(SIM, zs[ii]))
             for ii in range(len(zs_indices))]

if True : # enforce <kappa> = 0 exactly
    for ii in range(len(zs_indices)) :
        edges = np.loadtxt('./hmpdf_ws/edges_zs%.4f.dat'%zs[ii])
        centers = 0.5 * (edges[1:] + edges[:-1])

        kavg_DMO = np.sum(centers*sim_DMO[ii]) / np.sum(sim_DMO[ii])
        kavg_hydro = np.sum(centers*sim_hydro[ii]) / np.sum(sim_hydro[ii])
        if sim_BCM is not None :
            kavg_BCM = np.sum(centers*sim_BCM[ii]) / np.sum(sim_BCM[ii])

        kprime_DMO = centers - kavg_DMO
        kprime_hydro = centers - kavg_hydro
        if sim_BCM is not None :
            kprime_BCM = centers - kavg_BCM

        i_DMO = interp1d(kprime_DMO, sim_DMO[ii],
                         kind='cubic', bounds_error=False, fill_value=0)
        i_hydro = interp1d(kprime_hydro, sim_hydro[ii],
                           kind='cubic', bounds_error=False, fill_value=0)
        if sim_BCM is not None :
            i_BCM = interp1d(kprime_BCM, sim_BCM[ii],
                         kind='cubic', bounds_error=False, fill_value=0)

        sim_DMO[ii] = i_DMO(centers)
        sim_hydro[ii] = i_hydro(centers)
        if sim_BCM is not None :
            sim_BCM[ii] = i_BCM(centers)

        sim_DMO[ii] /= np.sum(sim_DMO[ii])
        sim_hydro[ii] /= np.sum(sim_hydro[ii])
        if sim_BCM is not None :
            sim_BCM[ii] /= np.sum(sim_BCM[ii])

        # FIXME
#        plt.semilogy(sim_hydro[ii,:,1])
#        plt.semilogy(sim_BCM[ii])
#        plt.show()

theory_DMO   = np.empty((len(sim_DMO), len(sim_DMO[0])))
theory_hydro = np.empty((len(sim_DMO), len(sim_DMO[0])))

for ii, zsource in enumerate(zs) :
    edg_file = 'hmpdf_ws/edges_zs%.4f.dat'%zsource

    tmp_file_DMO   = 'hmpdf_ws/op_%s_DMO_zs%.4f.dat'%(SIM.split('_')[0], zsource)
    tmp_file_hydro = 'hmpdf_ws/op_%s_%s_hydro_zs%.4f.dat'%(SIM, OBS, zsource)

    if RECOMPUTE :
        for hydro in [0, 1] :
            print('---------------' +
                  ('HYDRO' if hydro else 'DMO')
                  + '---------------')
            mcuts_file = 'hmpdf_ws/mass_cuts_DMO.dat' # DMO is the 'ficticious' mass we integrate over
            cmd = './run_hmpdf %.4f %s %s %s %d %d %d'%(zsource, edg_file, mcuts_file,
                                                        tmp_file_hydro if hydro else tmp_file_DMO,
                                                        hydro, sims.index(SIM), obss.index(OBS))
            print(cmd)
            system(cmd)

    theory_hydro[ii, :] = np.loadtxt(tmp_file_hydro)
    theory_DMO[ii, :]   = np.loadtxt(tmp_file_DMO)

fig_res, ax_res = plt.subplots()

fig_log, ax_log = plt.subplots(ncols=len(zs_indices))
if len(zs_indices) == 1 :
    ax_log = np.array([ax_log, ])
ax_log = ax_log.flatten()

for ii, zsource in enumerate(zs) :
    edges = np.loadtxt('./hmpdf_ws/edges_zs%.4f.dat'%zsource)
    centers = 0.5 * (edges[1:] + edges[:-1])

    lsim = ax_res.plot(centers,
                       (sim_hydro[ii]-sim_DMO[ii])/(0.5*(sim_DMO[ii]+sim_hydro[ii])),
                       linestyle='solid', label='$z_s = %.2f$'%zsource)
    if sim_BCM is not None :
        ax_res.plot(centers,
                    (sim_BCM[ii]-sim_DMO[ii])/(0.5*(sim_DMO[ii]+sim_BCM[ii])),
                    linestyle='dotted', color=plt.getp(lsim[0], 'color'))
    ax_res.plot(centers,
                (theory_hydro[ii]-theory_DMO[ii])/(0.5*(theory_DMO[ii]+theory_hydro[ii])),
                linestyle='dashed', color=plt.getp(lsim[0], 'color'))

    for hydro, (sim, theory) in enumerate(zip([sim_DMO, sim_hydro, sim_BCM], [theory_DMO, theory_hydro, None])) :
        if sim is not None :
            lsim = ax_log[ii].semilogy(centers,
                                       sim[ii] / np.sum(sim[ii]),
                                       linestyle='solid',
                                       label=['DMO', 'hydro', 'BCM (at sim level)'][hydro])
        if theory is not None :
            ax_log[ii].semilogy(centers, theory[ii, :],
                                linestyle='dashed', color=plt.getp(lsim[0], 'color'),
                                label='theory (%s)'%['DMO', 'hydro'][hydro])
        ax_log[ii].set_title('$z_s = %.2f$'%zsource)
    ax_log[ii].set_xlim(0.0, 0.2)
    ax_log[ii].set_ylim(1e-6, 1e-1)

ax_log[0].legend()
ax_res.set_xlim(0.0, 0.2)
ax_res.set_ylim(-0.3, 0.3)
ax_res.set_xlabel('convergence $\kappa$')
ax_res.set_ylabel('$2 (p_{\sf hydro} - p_{\sf DMO}) / (p_{\sf hydro} + p_{\sf DMO})$')
ax_res.legend()
fig_res.savefig('hydro_vs_DMO_%s_%s.pdf'%(SIM, OBS), bbox_inches='tight')
fig_log.savefig('smoothing_calib_%s_%s.pdf'%(SIM, OBS), bbox_inches='tight')
