#ifndef INIT_HMPDF_H
#define INIT_HMPDF_H

#include <math.h>
#include <omp.h>

#include "hmpdf.h"

#define ARICO20

// the different ways we approximate baryonic physics corrections
// DO NOT CHANGE THE ORDER
enum BARYON_MODES {
                    BCM, // baryonic correction model
                    TOT_CONC, // total concentration model changed
                    BAR_CONC, // two NFW components, only the baryonic one is free to vary
                    DMO, // no baryonic correction
                  };

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

// refers to DM component in the hydro sim, we first try to keep this one constant
static double
conc_hydro_DM[] = { 5.71, -0.087, -0.47,
                    7.85, -0.081, -0.71,
                    1.29506089e+01, -8.866189e-02, -8.6722558e-02,
                    2.486489e-02, 5.510048e-01, 1.11423926e+00 };

int
init_hmpdf (hmpdf_obj *d, double zs, int for_cov, enum BARYON_MODES baryon_mode,
            const double *params, int Nz_Arico, const double *z_Arico)
{
    int status;

    // use these for various concentration models which we may fit for
    double conc_params[12];
    if (baryon_mode==TOT_CONC || baryon_mode==BAR_CONC)
    {
        double *p = conc_params;
        for (int ii=0; ii<6; ++ii)
            *(p++) = conc_DM[ii]; // these are for the mass definitions that are barely used
        for (int ii=0; ii<6; ++ii)
            *(p++) = params[ii]; // these are the important concentration parameters
    }

    status = hmpdf_init(d, "illustris_cosmo.ini",
                        hmpdf_kappa, zs,
                        hmpdf_rout_scale, 2.5,
                        hmpdf_gaussian_fwhm, 2.0 * sqrt(M_LN2),

                        /* begin baryon stuff */
                        hmpdf_Arico20_params, (baryon_mode==BCM) ? params : NULL,
                        hmpdf_Arico20_Nz, Nz_Arico,
                        hmpdf_Arico20_z, z_Arico,
                        hmpdf_Duffy08_conc_params, (baryon_mode==TOT_CONC) ? conc_params : conc_DM,
                        hmpdf_DM_conc_params, (baryon_mode==BAR_CONC) ? conc_hydro_DM : NULL,
                        hmpdf_bar_conc_params, (baryon_mode==BAR_CONC) ? conc_params : NULL,
                        /* end baryon stuff */

                        hmpdf_N_threads, (int)(omp_get_max_threads()),
                        hmpdf_N_z, 40, /*this should still be accurate enough and fits well with our compute layout*/
                        hmpdf_pixel_side, 0.29, // TODO make this a hyperparameter?
                                                // rather not because otherwise we may compromise DMO fit
                        hmpdf_custom_k_filter, &k_filter,
                        hmpdf_signal_max, 0.6,
                        hmpdf_N_theta, 1024, // TODO will this be too slow?
                        /* for covariance matrix stability */
                        hmpdf_N_signal, (for_cov) ? 512L : 4096L,
                        hmpdf_N_phi, 10000,
                        hmpdf_phi_pwr, 5.0,
                        hmpdf_phi_jitter, 0.05,
                        hmpdf_verbosity, 0);

    if (status)
        return status;

    return 0;
}

#endif
