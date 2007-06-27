/* Solvable 6 parameter hill model with infall in the core, calculates fit 
   and returns rms residual.

   parameters:
   0: tau_core
   1: tau_env
   2: v_lsr
   3: v_in -- infall in the core
   4: sigma
   6: Tpeak
*/

void hill6_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax);
void hill6_free();
double *hill6_getfit();
double hill6_evaluate(double *params);
