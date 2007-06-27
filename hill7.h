/* Solvable 7 parameter hill model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_core
   1: tau_env
   2: v_lsr
   3: v_in -- infall in the envelope
   4: sigma
   5: To
   6: Tpeak
*/

void hill7_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax);
void hill7_free();
double *hill7_getfit();
double hill7_evaluate(double *params);
