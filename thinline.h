/* Solvable 4 parameter thin line model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: sigma
   3: Tex
*/

void thinline_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax);
void thinline_free();
double *thinline_getfit();
double thinline_evaluate(double *params);
