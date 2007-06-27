/* Solvable 6 parameter two layer model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: v_in
   3: sigma
   4: Tr
   5: Tf
*/

void twolayer6_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax);
void twolayer6_free();
double *twolayer6_getfit();
double twolayer6_evaluate(double *params);
