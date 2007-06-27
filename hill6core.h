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

void hill6core_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax);
void hill6core_free();
double *hill6core_getfit();
double hill6core_evaluate(double *params);
