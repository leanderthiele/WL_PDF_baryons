
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

    int Nbins;
    double **x0 = loadtxt("edges.dat", &Nbins, 1);
    Nbins -= 1;
    double *binedges = x0[0];

    double params_Arico[200];
    int Nz_Arico;
    double z_Arico[10];

    if (hydro)
        populate_params(&Nz_Arico, z_Arico, params_Arico, TNG, QK);

    hmpdf_obj *d = hmpdf_new();
    if (!d)
        return 1;

    int status;

    status = init_hmpdf(d, zs, (hydro) ? params_Arico : NULL, Nz_Arico, z_Arico); 
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

    sprintf(buffer, "cov_%s_zs%.4f.bin", (hydro) ? "hydro" : "DMO", zs);
    tofile(buffer, Nbins*Nbins, 1, cov);

    status = hmpdf_delete(d);
    if (status)
        return 4;
    free(cov);
}
