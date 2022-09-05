// Command line arguments :
//      [1] source redshift   (float)
//      [2] file w/ edges     (char *)
//      [3] file w/ mass cuts (char *)
//      [4] outfile name      (char *)
//      [5] hydro             (int)
//      [6] simulation        (int) -- 0=TNG,
//                                     1=BAHAMAS_fid, 2=BAHAMAS_loAGN, 3=BAHAMAS_hiAGN
//                                     4=BAHAMAShires_fid, 5=BAHAMAShires_loAGN, 6=BAHAMAShires_hiAGN
//      [7] BCM data          (int) -- 0=PK, 1=QK, 2=PK+QK

// If ARICO20 is defined, use the updated 8-parameter BCM
#define ARICO20

// If SPLIT is defined, use separate NFW profiles for DM and baryons
// #define SPLIT

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>

#include "hmpdf.h"
#include "utils.h"

#include "params_arico.h"

double
k_filter(double k, double z, void *p)
{
    (void)p;
    double R = 0.025 * pow(1.0 + z, 1.0);
    return pow(1.0 + pow(k*R, 1.0), -1.5);
}

static double
hydro_hmf_corr_params[] = { 0.0150287, -0.43318792, 0.01964109, -0.30671152, 0.06731122 };

double
hydro_hmf_corr(double z, double M, void *p)
{
    (void)p;
    M = log(M/2e12);
    return pow(1.0+z, hydro_hmf_corr_params[0])
           * (1.0 + exp(-gsl_pow_2(M/log(10.0)-1.0))
                    * exp(hydro_hmf_corr_params[1] * M)
                    * (hydro_hmf_corr_params[2]
                       + hydro_hmf_corr_params[3] * M
                       + hydro_hmf_corr_params[4] * gsl_pow_2(M)));
}

static double
hydro_conc_resc_params[] = { 0.31457869, 0.05458279, 1.49058383, -0.14363158, 0.04816107, -0.0039402 };

double
hydro_conc_resc(double z, double M, void *p)
{
    (void)p;

    if (M > 1e15)
        return 1.0;
    else if (M < 1e12)
        return 1.2;

    M = log(M/2e12);
    return 1.0 + exp( - gsl_pow_2((M-log(10.0))/log(100.0)) )
                 * pow(1.0 + z, hydro_conc_resc_params[0])
                 * (hydro_conc_resc_params[1]
                        * pow(1.0 + z, hydro_conc_resc_params[2])
                    + hydro_conc_resc_params[3] * M
                    + hydro_conc_resc_params[4] * gsl_pow_2(M)
                    + hydro_conc_resc_params[5] * gsl_pow_3(M) );
}

static double
mass_resc_params[] = { -0.14246705, -0.7154184,  -0.52396824,  1.17713498 };

double
mass_resc(double z, double M, void *p)
{
    (void)p;
    M = log(M/2e12);
    double out = 1.0 + mass_resc_params[0]
                       * pow(1.0+z, mass_resc_params[1])
                       * exp( mass_resc_params[2]
                             * gsl_pow_2(M - mass_resc_params[3]) );

    return out;
}

// define the concentration model

static double
conc_Duffy08[] = { 5.71, -0.087, -0.47,
                   7.85, -0.081, -0.71,
                  10.14, -0.081, -1.01, 0.0, 0.0, 0.0 };

// FIXME update
// concentration in hydro simulation (total matter)
static double
conc_hydro_tot[] = { 5.71, -0.087, -0.47,
                     7.85, -0.081, -0.71,
                     1.375023001e+01, -8.688706e-02, -8.4451454e-01,
                     2.705265e-02, 5.2513069e-01, 1.12049399e+00 };

// FIXME update
// concentration in DMO simulation
static double
conc_DM[] = { 5.71, -0.087, -0.47,
              7.85, -0.081, -0.71,
              1.29729068e+01, -7.05810865e-02, -1.07969465e+00, 
              5.92621692e-04, 4.23256450e-01, -7.01235269e+00 };

// this only refers to the DM component
static double
conc_hydro_DM[] = { 5.71, -0.087, -0.47,
                    7.85, -0.081, -0.71,
                    1.29506089e+01, -8.866189e-02, -8.6722558e-02,
                    2.486489e-02, 5.510048e-01, 1.11423926e+00 };

// this only refers to the baryonic component -- only need this in terms of the internal mass
// FIXME need to write code for this!
static double
conc_hydro_bar[] = { 5.71, -0.087, -0.47,
                     7.85, -0.081, -0.71,
                     4.57478936e+00, -1.67755068e+00, -2.84129704e+00,
                     4.90993037e-04, -1.16223295e+01, -6.83565016e+00 };
                     

// define the mass cuts function
static double
mass_cuts(double z, void *p)
{
    double out;
    interp1d *i = (interp1d *)p;

    if (interp1d_eval(i, z, &out))
        printf("interp1d_eval failed in mass_cuts\n");

    // did interpolation in log
    return exp(out);
}

// define the bias rescaling function
static double
bias_resc(double z, double M, void *p)
{
    (void)p;
    return 0.1;
}

int
main(int argc, char **argv)
{
#ifdef SPLIT
    printf("Using the 2-component NFW profiles, no BCM!\n");
#else
    printf("Using the BCM model\n");
#endif

    double zs = atof(argv[1]);
    char *edg = argv[2];
    char *mcut = argv[3];
    char *out = argv[4];
    int hydro = atoi(argv[5]);
    enum SIM sim = atoi(argv[6]);
    enum OBS obs = atoi(argv[7]);

    // load the binedges
    int Nbins;
    double **x0 = loadtxt(edg, &Nbins, 1);
    Nbins -= 1;
    double *binedges = x0[0];

    // interpolate the mass cuts
    int Nz;
    double **x1 = loadtxt(mcut, &Nz, 2);
    double *mcuts_z = x1[0];
    double *mcuts_M = x1[1];
    // better to interpolate in log I believe
    for (int ii=0; ii<Nz; ii++)
        mcuts_M[ii] = log(mcuts_M[ii]);
    // NOTE that mass cuts is not called in parallel so no need for separate accelerators
    interp1d *mcuts_interpolator;
    if (new_interp1d(Nz, mcuts_z, mcuts_M, mcuts_M[0], mcuts_M[Nz-1],
                     interp_linear, NULL, &mcuts_interpolator))
        printf("new_interp1d failed\n");

#ifndef SPLIT
    // set the Arico BCM
    double params_Arico[80];
    int Nz_Arico;
    double z_Arico[10];
    populate_params(&Nz_Arico, z_Arico, params_Arico, sim, obs);
#endif // SPLIT

    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return 1;

    char ini_file[64];
    if (sim == TNG)
        sprintf(ini_file, "./illustris_cosmo.ini");
    else
        sprintf(ini_file, "./bahamas_cosmo.ini");

    // need to use the function here to be able to have preprocessor directives inside call
    if (hmpdf_init_fct(d, ini_file,
                   hmpdf_kappa, zs,

                   hmpdf_rout_scale, 2.5,

                   // TODO we don't have the Gaussian filter applied to the maps included here!!!
                   // It is W(\theta) \propto \exp(-\theta / \theta_G) where \theta_G=1 arcmin
                   // We need the FWHM of this filter
                   // FIXME for some reason this suppresses the PDF *way* too much compared to the
                   //       sims! (even when I take out the k-filter)
                   // TODO if I manually decrease the filter size by \pi, the PDFs look fairly reasonable
                   //      with very good agreement at large zs and visibly worse for zs~0.5
                   hmpdf_gaussian_fwhm, (sim==TNG) ? 2.0 * sqrt(M_LN2) * 1.0 : 1e-3/*some small value*/,

//                   hmpdf_warn_is_err, 0,

#ifndef SPLIT
                   hmpdf_Arico20_Nz, Nz_Arico,
                   hmpdf_Arico20_z, z_Arico,
                   hmpdf_Arico20_params, (hydro) ? params_Arico : NULL,
#endif // SPLIT

                   hmpdf_N_threads, 4,
                   hmpdf_pixel_side, (sim==TNG) ? 0.29 // == 5 deg / 1024 pixels
                                     : (sim>=BAHAMAShires_fid && sim<= BAHAMAShires_hiAGN) ? 0.5*0.1666666
                                     : 0.1*0.882, // hacky BAHAMAS lo-res

#ifndef SPLIT
                   hmpdf_Duffy08_conc_params, (sim==TNG) ? conc_DM : conc_Duffy08,
#else
                   // in hydro case, we use the concentration for the total profile for mass conversions only
                   hmpdf_Duffy08_conc_params, (hydro) ? conc_hydro_tot : conc_DM,
                   hmpdf_DM_conc_params, (hydro) ? conc_hydro_DM : NULL,
                   hmpdf_bar_conc_params, (hydro) ? conc_hydro_bar : NULL,
                   hmpdf_mass_resc, (hydro) ? &mass_resc : NULL,
#endif // SPLIT
//                   hmpdf_massfunc_corr, (hydro) ? &hydro_hmf_corr : NULL,
//                 NOTE that the massfunc_corr seems to give very similar results to the mass rescaling,
//                      which is encouraging
//                   hmpdf_conc_resc, (hydro) ? &hydro_conc_resc : NULL,
//                   hmpdf_mass_cuts, &mass_cuts,
//                   hmpdf_mass_cuts_params, (void *)mcuts_interpolator,
//                   hmpdf_bias_resc, (hydro) ? &bias_resc : NULL,
                   hmpdf_custom_k_filter, (sim==TNG) ? &k_filter : NULL,
                   hmpdf_signal_max, 3.0*binedges[Nbins],
                   hmpdf_N_signal, 4096L,
//                   hmpdf_N_M, 100,
//                   hmpdf_N_z, 100,

                   // N_theta=500 is not completely converged (1000 is, but quite a bit more expensive )
                   // This is mostly an issue for zs=0.5,
                   // but even there qualitatively the result is prety converged.
                   hmpdf_N_theta, 500,
                   hmpdf_verbosity, 0,
                   hmpdf_end_configs))
        return 2;

    double op[Nbins];
    if (hmpdf_get_op(d, Nbins, binedges, op, 1, 0))
        return 3;

    savetxt(out, Nbins, 1, op);

    if (hmpdf_delete(d))
        return 4;

    return 0;
}
