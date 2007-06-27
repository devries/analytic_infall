/* Solvable 5 parameter two layer model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: v_in
   3: sigma
   4: Tr
*/

void twolayer5_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax);
void twolayer5_free();
double *twolayer5_getfit();
double twolayer5_evaluate(double *params);
