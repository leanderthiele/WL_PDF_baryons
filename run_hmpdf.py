import numpy as np
from matplotlib import pyplot as plt
from sys import argv
from os import system
from os.path import isfile

hydro = int(argv[1])

zs_indices = [12, 22, 34]
zs = np.loadtxt('kappaTNG_PDF/zs.dat')[zs_indices]

sim_DMO = np.load('kappaTNG_PDF/DMO/PDF_smooth1arcmin.npy')[zs_indices,...]

theory = np.empty((sim_DMO.shape[0], sim_DMO.shape[1]))

def edges(centres) :
    out = np.empty(len(centres)+1)
    out[1:-1] = 0.5 * (centres[1:]+centres[:-1])
    out[0]  = centres[0] - 0.5*(centres[1]-centres[0])
    out[-1] = centres[-1] + 0.5*(centres[-1]-centres[-2])
    return out

for ii, zsource in enumerate(zs) :
    edg_file = 'hmpdf_ws/edges_zs%.4f.dat'%zsource
    if not isfile(edg_file) :
        binedges = edges(sim_DMO[ii, :, 0])
        np.savetxt(edg_file, binedges)

    tmp_file = 'hmpdf_ws/op_zs%.4f.dat'%zsource

    system('./run_hmpdf %.4f %s %s %d'%(zsource, edg_file, tmp_file, hydro))

    theory[ii, :] = np.loadtxt(tmp_file)

fig, ax = plt.subplots()

for ii, zsource in enumerate(zs) :
    lsim = ax.semilogy(sim_DMO[ii, :, 0], sim_DMO[ii, :, 1]/np.sum(sim_DMO[ii, :, 1]),
                       linestyle='solid', label='$z_s = %.2f$'%zsource)
    ax.semilogy(sim_DMO[ii, :, 0], theory[ii, :],
                linestyle='dashed', color=plt.getp(lsim[0], 'color'))

ax.legend(loc='upper right')
fig.savefig('comp_hmpdf_illustris.pdf', bbox_inches='tight')
