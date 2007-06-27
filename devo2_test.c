#include "devo2.h"

double parabola(double *vec) {
  double res = 3.0+(vec[0]-1.0)*(vec[0]-1.0)+(vec[1]-2.0)*(vec[1]-2.0);
  return res;
}

int main() {
  devo2_struct st;
  double min[2];
  double max[2];

  min[0] = 0.0;
  max[0] = 2.0;

  min[1] = 0.0;
  max[1] = 5.0;

  devo2_init(&st,2,min,max,15,0.4,0.8,parabola);

  do {
    printf("%lf: %lf,%lf\n",st.best_score, st.best_vector[0],st.best_vector[1]);
    fflush(stdout);
    devo2_step(&st);
  } while(st.best_score>3.000001);

  printf("\n");
  printf("Converged to %lf after %d generations.\n", st.best_score, st.gen_count);
  printf("Parameter: %lf,%lf\n", st.best_vector[0],st.best_vector[1]);

  devo2_free(&st);

  exit(0);
}

  
