/* Command line arguments:
 * [1] output file ... will be a text file with the first column DMO and the second column hydro
 * Now the Arico20 parameters (all floating point numbers):
 * [2] M_c
 * [3] M_1_z0_cen
 * [4] eta
 * [5] beta
 * [6] theta_inn
 * [7] theta_out
 * [8] M_inn
 * [9] M_r
 * The masses are in log10(Msun/h) [to make sampling easier]
 * and all parameters are assumed constant across redshift.
 */

#define ARICO20

#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "hmpdf.h"
#include "utils.h"

#include "init_hmpdf.h"

const double zs = 1.034; // this is where we run

const int Nbins = 32;
const double kappa_min = 0.0,
             kappa_max = 0.2;

int main (int argc, char **argv)
{
    if (argc != 1/*executable*/ + 1/*outfile*/ + hmpdf_Arico20_Nparams)
        return -1;
    
    char **c = argv;

    char *outfile = *(c++);
    printf("%s\n", outfile);

    double params_Arico[hmpdf_Arico20_Nparams];
    for (int ii=0; ii<hmpdf_Arico20_Nparams; ++ii)
        params_Arico[ii] = atof(*(c++));

    // some masses are special cases, convert them to h-less units
    const double h = 0.6711; const double logh = log10(h);
    params_Arico[hmpdf_Arico20_M_c] = pow(10.0, params_Arico[hmpdf_Arico20_M_c] - logh);
    params_Arico[hmpdf_Arico20_M_1_z0_cen] -= logh;
    params_Arico[hmpdf_Arico20_M_inn] = pow(10.0, params_Arico[hmpdf_Arico20_M_inn] - logh);
    params_Arico[hmpdf_Arico20_M_r] = pow(10.0, params_Arico[hmpdf_Arico20_M_r] - logh);

    const int Nz_Arico = 1;
    const double z_Arico[] = { 0.5 };

    double op[2][Nbins];

    // NOTE if we are only varying parameters pertaining to the hydro computation, it's not
    //      necessary to recompute the DMO version every time.
    //      However, the cost is marginal compared to hydro and it's easier to introduce additional
    //      nuisance parameters later on if we wish.

    for (int hydro=0; hydro<2; ++hydro)
    {
        hmpdf_obj *d = hmpdf_new(); if (!d) return 1;
        int status;
        // TODO maybe make pixel side or something in the k-filter a hyperparameter too?
        status = init_hmpdf(d, zs, /*for_cov=*/0, (hydro) ? params_Arico : NULL, Nz_Arico, z_Arico);
        if (status) return 2;

        double binedges[Nbins+1];
        linspace(Nbins+1, kappa_min, kappa_max, binedges);

        status = hmpdf_get_op(d, Nbins, binedges, op[hydro], /*incl_2h=*/1, /*noisy=*/0);
        if (status) return 3;

        status = hmpdf_delete(d);
        if (status) return 4;
    }

    savetxt(outfile, Nbins, 2, op[0], op[1]);
}
