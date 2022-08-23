
#include <stdlib.h>

#define ARICO20
#include "hmpdf.h"
#include "utils.h"

#include "init_hmpdf.h"
#include "params_arico.h"

int main (int argc, char **argv)
// argv[1] ... 0 for DMO, 1 for hydro
// argv[2] ... source redshift
{
    int hydro = atoi(argv[1]);
    double zs = atof(argv[2]);

    int Nbins = 32;
    double *binedges = malloc((Nbins+1)*sizeof(double));
    linspace(Nbins+1, 0.0, 0.2, binedges);

    double params_Arico[200];
    int Nz_Arico;
    double z_Arico[10];

    if (hydro)
        populate_params(&Nz_Arico, z_Arico, params_Arico, TNG, QK);

    hmpdf_obj *d = hmpdf_new();
    if (!d)
        return 1;

    int status;

    status = init_hmpdf(d, zs, /*for_cov=*/1, (hydro) ? params_Arico : NULL, Nz_Arico, z_Arico); 
    if (status)
        return 2;

    double op[Nbins];
    status = hmpdf_get_op(d, Nbins, binedges, op, 1, 0); 
    if (status)
        return 3;

    char buffer[512];
    sprintf(buffer, "op_%s_zs%.4f.bin", (hydro) ? "hydro" : "DMO", zs);
    tofile(buffer, Nbins, 1, op);

    double *cov = malloc(Nbins * Nbins * sizeof(double));
    status = hmpdf_get_cov(d, Nbins, binedges, cov, 0);
    if (status)
        return 4;

    sprintf(buffer, "cov_%s_zs%.4f.bin", (hydro) ? "hydro" : "DMO", zs);
    tofile(buffer, Nbins*Nbins, 1, cov);

    // get some diagnostics for debugging
    int Nphi;
    double *phi, *phiweights, *corr_diagn;
    status = hmpdf_get_cov_diagnostics(d, &Nphi, &phi, &phiweights, &corr_diagn);
    if (status)
        return 5;

    sprintf(buffer, "cov_diagn_%s_zs%.4f.bin", (hydro) ? "hydro" : "DMO", zs);
    tofile(buffer, Nphi, 3, phi, phiweights, corr_diagn);

    status = hmpdf_delete(d);
    if (status)
        return 6;
    free(binedges); free(cov); free(phi); free(phiweights); free(corr_diagn);
}
