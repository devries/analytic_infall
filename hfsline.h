/* Solvable 4 parameter hfs line model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_main
   1: v_lsr
   2: sigma
   3: Tex
*/

void hfsline_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax, int ncomp, double *comp_voff, double *comp_relint);
void hfsline_free();
double *hfsline_getfit();
double hfsline_evaluate(double *params);
