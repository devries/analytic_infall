/* Solvable 5 parameter hill model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: v_in
   3: sigma
   4: Tpeak
*/

double get_min_lsr();
double get_max_lsr();
void hill5_multicomponent_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax, int ncomp, double *comp_voff, double *comp_relint);
void hill5_multicomponent_free();
double *hill5_multicomponent_getfit();
double hill5_multicomponent_evaluate(double *params);
