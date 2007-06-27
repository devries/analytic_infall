/* Solvable 5 parameter hill model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: v_in
   3: sigma
   4: Tpeak
*/

void hill5_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax);
void hill5_free();
double *hill5_getfit();
double hill5_evaluate(double *params);
