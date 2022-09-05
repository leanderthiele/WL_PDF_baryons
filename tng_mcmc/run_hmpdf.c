/* Command line arguments:
 * [1] output file ... will be a text file with the first column DMO and the second column hydro
 * [2] baryonic correction mode ... integer, 0...BCM, 1...TOT_CONC, 2...BAR_CONC
 * [3] simulation ... integer, 0...TNG, 1...BAHAMAS
 * If [2]==BCM :
 *     Now the Arico20 parameters (all floating point numbers):
 *     [ 3] M_c
 *     [ 4] M_1_z0_cen 
 *     [ 5] eta
 *     [ 6] beta
 *     [ 7] theta_inn
 *     [ 8] theta_out
 *     [ 9] M_inn
 *     [10] M_r
 * If [2]==TOT_CONC / BAR_CONC
 *     Now the concentration model parameters (there are 6)
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
    char **c = argv+1;

    char *outfile = *(c++);
    enum BARYON_MODES baryon_mode = atoi(*(c++));
    enum SIMS sim = atoi(*(c++));

    if (argc != 1/*executable*/+1/*outfile*/+1/*baryon mode*/+1/*simulation*/
                +((baryon_mode==BCM) ? hmpdf_Arico20_Nparams : 6)/*theta*/)
        return -1;

    double theta[32]; // fill this with the parameter vector, definitely smaller than 32

    for (int ii=0; ii<((baryon_mode==BCM) ? hmpdf_Arico20_Nparams : 6); ++ii)
        theta[ii] = atof(*(c++));

    if (baryon_mode==BCM)
    {
        // some masses are special cases, convert them to h-less units
        const double h = 0.6711; const double logh = log10(h);
        theta[hmpdf_Arico20_M_c] = pow(10.0, theta[hmpdf_Arico20_M_c] - logh);
        theta[hmpdf_Arico20_M_1_z0_cen] -= logh;
        theta[hmpdf_Arico20_M_inn] = pow(10.0, theta[hmpdf_Arico20_M_inn] - logh);
        theta[hmpdf_Arico20_M_r] = pow(10.0, theta[hmpdf_Arico20_M_r] - logh);
    }

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
        // I think we have concluded this is a bad idea because it would destroy the DMO fit
        status = init_hmpdf(d, zs, /*for_cov=*/0, /*baryon_mode=*/(hydro) ? baryon_mode : DMO, sim,
                           (hydro) ? theta : NULL, Nz_Arico, z_Arico);
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
