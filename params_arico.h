#ifndef PARAMS_ARICO_H
#define PARAMS_ARICO_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "hdf5/serial/hdf5.h"

#include "hmpdf.h"

// the simulation, don't change order
enum SIM { TNG,
           BAHAMAS_fid, BAHAMAS_loAGN, BAHAMAS_hiAGN,
           BAHAMAShires_fid, BAHAMAShires_loAGN, BAHAMAShires_hiAGN };

// the observation fit to, don't change order
enum OBS { PK, QK, PK_QK };

#ifdef ARICO20
int read_attr (hid_t grp, const char *name, double *out)
{
    int status;

    hid_t attr = H5Aopen(grp, name, H5P_DEFAULT);
    if (attr<0) return 1;
    
    status = H5Aread(attr, H5T_NATIVE_DOUBLE, out);
    if (status<0) return 1;

    status = H5Aclose(attr);
    if (status<0) return 1;

    return 0;
}

int read_attrs_from_h5 (const char *grp_name, double *params)
{
    hid_t file = H5Fopen("hydro_simulations_fit_power_and_bispectra.hdf5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file<0) return 1;

    hid_t group = H5Gopen(file, grp_name, H5P_DEFAULT);
    if (group<0) return 2; 

    if (read_attr(group, "eta", params+hmpdf_Arico20_eta)) return 3;
    if (read_attr(group, "M_c", params+hmpdf_Arico20_M_c)) return 4;
    if (read_attr(group, "beta", params+hmpdf_Arico20_beta)) return 5;
    if (read_attr(group, "M1_z0_cen", params+hmpdf_Arico20_M_1_z0_cen)) return 6;
    if (read_attr(group, "M_r", params+hmpdf_Arico20_M_r)) return 7;
    if (read_attr(group, "M_inn", params+hmpdf_Arico20_M_inn)) return 8;
    if (read_attr(group, "theta_inn", params+hmpdf_Arico20_theta_inn)) return 9;
    if (read_attr(group, "theta_out", params+hmpdf_Arico20_theta_out)) return 10;

    int status;

    status = H5Gclose(group);
    if (status<0) return 11;

    status = H5Fclose(file);
    if (status<0) return 12;

    return 0;
}
#endif // ARICO20

void populate_params (int *Nz, double *z, double *params, enum SIM sim, enum OBS obs)
{
    double h;
    if (sim == TNG) h = 0.6774;
    else h = 0.7;

#ifdef ARICO20
    char sim_str[32];
    switch (sim)
    {
        case TNG :           sprintf(sim_str, "tng"); break;
        case BAHAMAShires_fid :
        case BAHAMAS_fid :   sprintf(sim_str, "bahamas"); break;
        case BAHAMAShires_loAGN :
        case BAHAMAS_loAGN : sprintf(sim_str, "bahamas_lowAGN"); break;
        case BAHAMAShires_hiAGN :
        case BAHAMAS_hiAGN : sprintf(sim_str, "bahamas_hiAGN"); break;
    }

    char obs_str[8];
    switch (obs)
    {
        case PK :    sprintf(obs_str, "pk"); break;
        case QK :    sprintf(obs_str, "qk"); break;
        case PK_QK : sprintf(obs_str, "pk_qk"); break;
    }

    *Nz = 3;
    z[0] = 0.0; z[1] = 1.0; z[2] = 2.0;

    for (int ii=0; ii<*Nz; ++ii)
    {
        char grp_name[256];
        sprintf(grp_name, "%s/%s/%.1f/pars", obs_str, sim_str, z[ii]);

        int status = read_attrs_from_h5(grp_name, params);

        if (status)
        {
            fprintf(stderr, "read_attrs_from_h5 failed with error %d\n", status);
            exit(-1);
        }

        // do a few transformations (mostly transform Msun/h -> Msun as these are internal code units)
        params[hmpdf_Arico20_M_c] /= h;
        params[hmpdf_Arico20_M_1_z0_cen] = log10(params[hmpdf_Arico20_M_1_z0_cen]/h);
        params[hmpdf_Arico20_M_r] /= h;
        params[hmpdf_Arico20_M_inn] /= h;

        // advance position in the array to the next redshift
        params += hmpdf_Arico20_Nparams;
    }
#else
    #error
    *Nz = 1;
    if (sim == TNG)
    {
        params[hmpdf_Arico20_eta] = 0.14; /* vary: 0.1 (mm) -- 0.10 (m) -- 0.14 (fid) -- 0.18 (p) -- 0.5 (pp) */
        params[hmpdf_Arico20_M_c] = 2.3e13/h; /* vary: 0.1 (mm) -- 1.15 (m) -- 2.3 (fid) -- 4.6 (p) -- 10 (p1) -- 20 (p2) -- 40 (p3) -- 600.0 (pp) */
        params[hmpdf_Arico20_beta] = 4.09; /* vary: 0.2 (mm) -- 2.5 (m) -- 4.09 (fid) -- 5.5 (p) -- 9.0 (pp) */
        params[hmpdf_Arico20_M_1_z0_cen] = log10(2.2e10/h); /* vary: 0.1 (mm) -- 1.4 (m) -- 2.2 (fid) -- 3.0 (p) -- 100.0 (pp) */
    }
#endif
}

#endif
