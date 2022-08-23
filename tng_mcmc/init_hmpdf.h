#ifndef INIT_HMPDF_H
#define INIT_HMPDF_H

#include <math.h>
#include <omp.h>

#include "hmpdf.h"

#define ARICO20

double
k_filter(double k, double z, void *p)
{
    (void)p;
    double R = 0.025 * pow(1.0 + z, 1.0);
    return pow(1.0 + pow(k*R, 1.0), -1.5);
}


static double
conc_DM[] = { 5.71, -0.087, -0.47,
              7.85, -0.081, -0.71,
              1.29729068e+01, -7.05810865e-02, -1.07969465e+00, 
              5.92621692e-04, 4.23256450e-01, -7.01235269e+00 };

int
init_hmpdf (hmpdf_obj *d, double zs, int for_cov, double *params_Arico, int Nz_Arico, double *z_Arico)
// If params_Arico == NULL, the DMO version is run
{
    int status;

    status = hmpdf_init(d, "illustris_cosmo.ini",
                        hmpdf_kappa, zs,
                        hmpdf_rout_scale, 2.5,
                        hmpdf_gaussian_fwhm, 2.0 * sqrt(M_LN2),
                        hmpdf_Arico20_params, params_Arico,
                        hmpdf_Arico20_Nz, Nz_Arico,
                        hmpdf_Arico20_z, z_Arico,
                        hmpdf_N_threads, (int)(omp_get_max_threads()),
                        hmpdf_pixel_side, 0.29,
                        hmpdf_Duffy08_conc_params, conc_DM,
                        hmpdf_custom_k_filter, &k_filter,
                        hmpdf_signal_max, 0.6,
                        hmpdf_N_theta, 300,
                        /* for covariance matrix stability */
                        hmpdf_N_signal, (for_cov) ? 1024L : 4096L,
                        hmpdf_N_phi, 10000,
                        hmpdf_phi_pwr, 5.0,
                        hmpdf_verbosity, 1);

    if (status)
        return status;

    return 0;
}

#endif
