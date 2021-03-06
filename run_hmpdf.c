// Command line arguments :
//      [1] source redshift   (float)
//      [2] file w/ edges     (char *)
//      [3] file w/ mass cuts (char *)
//      [4] outfile name      (char *)
//      [5] hydro             (int)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>

#include "hmpdf.h"
#include "utils.h"

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
    return 1.0 + mass_resc_params[0]
                 * pow(1.0+z, mass_resc_params[1])
                 * exp( mass_resc_params[2]
                        * gsl_pow_2(M - mass_resc_params[3]) );
}

// define the concentration model
static double
conc_DM[] = { 5.71, -0.087, -0.47,
              7.85, -0.081, -0.71,
              1.29729068e+01, -7.05810865e-02, -1.07969465e+00,  5.92621692e-04, 4.23256450e-01, -7.01235269e+00 };

static double
conc_hydro[] = { 5.71, -0.087, -0.47,
                 7.85, -0.081, -0.71,
                 12.79550921,
                 -0.08862842, -0.84718926,  0.02557855,  0.59983028,  1.09660615 };

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

int
main(int argc, char **argv)
{
    double zs = atof(argv[1]);
    char *edg = argv[2];
    char *mcut = argv[3];
    char *out = argv[4];
    int hydro = atoi(argv[5]);

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

    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return 1;

    if (hmpdf_init(d, "./illustris_cosmo.ini",
                   hmpdf_kappa, zs,

                   hmpdf_N_threads, 4,
                   hmpdf_pixel_side, 0.29,
                   hmpdf_Duffy08_conc_params, (hydro) ? conc_hydro : conc_DM,
//                   hmpdf_massfunc_corr, (hydro) ? &hydro_hmf_corr : NULL,
//                   hmpdf_conc_resc, (hydro) ? &hydro_conc_resc : NULL,
//                   NOTE it is possible that the implementation of conc_resc in hmpdf
//                        is buggy
                   hmpdf_mass_resc, (hydro) ? &mass_resc : NULL,
//                   hmpdf_mass_cuts, &mass_cuts,
//                   hmpdf_mass_cuts_params, (void *)mcuts_interpolator,
                   hmpdf_custom_k_filter, &k_filter,
                   hmpdf_signal_max, 3.0*binedges[Nbins],
                   hmpdf_N_signal, 4096,
                   hmpdf_N_M, 100,
                   hmpdf_N_z, 100,
                   hmpdf_N_theta, 300,
                   hmpdf_verbosity, 0))
        return 2;

    double op[Nbins];
    if (hmpdf_get_op(d, Nbins, binedges, op, 1, 0))
        return 3;

    savetxt(out, Nbins, 1, op);

    if (hmpdf_delete(d))
        return 4;

    return 0;
}
