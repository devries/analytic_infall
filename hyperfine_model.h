/* Solvable hyperfine model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: sigma
   3: Tex
*/
typedef struct _hyperfine_struct {
  double *velocity_array;
  double *temperature_array;
  double *hyperfine_array;
  double frequency;
  int nchan;
  double vrange[2];
  double lsrrange[2];
  int n_components; /* Number of hyperfine components */
  double *comp_voff_array; /* hyperfine component velocities */
  double *comp_relint_array; /* Hyperfine component relative intensities */
} hyperfine_struct;

double get_min_lsr(hyperfine_struct *st);
double get_max_lsr(hyperfine_struct *st);
void hyperfine_init(hyperfine_struct *st, int channels, double *varray, double *tarray, double nu, double vmin, double vmax, int ncomp, double *comp_voff, double *comp_relint);
void hyperfine_free(hyperfine_struct *st);
double *hyperfine_getfit(hyperfine_struct *st);
double hyperfine_evaluate(hyperfine_struct *st, double *params);
