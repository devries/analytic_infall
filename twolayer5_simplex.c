#include <stdio.h>
#include <stdlib.h>
#include "twolayer5.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

extern int read_data(FILE *fpin, int ncol, int maxlen, double *retdata[]);

double twolayer5_gsl(const gsl_vector *x, void *junk) {
  return twolayer5_evaluate(x->data);
}

int main(int argc, char *argv[]) {
  FILE *fpin;
  FILE *fpout;
  double *input_data[2];
  int nchan;
  gsl_vector *x;
  gsl_vector *step_size;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  double size;
  const double eps = 1.0e-6;
  double *model_spectrum;
  int i;
  int status;

  if(argc!=11) {
    fprintf(stderr, "Usage: %s <inputfilename> <frequency> <vmin> <vmax> <tau> <vlsr> <vin> <sigma> <tr> <outputfile>\n", argv[0]);
    exit(1);
  }

  fpin = fopen(argv[1],"r");
  nchan = read_data(fpin,2,132,input_data);
  fclose(fpin);

  x = gsl_vector_alloc(5);
  gsl_vector_set(x,0,atof(argv[5]));
  gsl_vector_set(x,1,atof(argv[6]));
  gsl_vector_set(x,2,atof(argv[7]));
  gsl_vector_set(x,3,atof(argv[8]));
  gsl_vector_set(x,4,atof(argv[9]));
  step_size = gsl_vector_alloc(5);
  gsl_vector_set(step_size,0,0.1);
  gsl_vector_set(step_size,1,0.005);
  gsl_vector_set(step_size,2,0.005);
  gsl_vector_set(step_size,3,0.005);
  gsl_vector_set(step_size,4,0.05);

  minex_func.f = twolayer5_gsl;
  minex_func.n = 5;
  minex_func.params = NULL;

  twolayer5_init(nchan,input_data[0],input_data[1],atof(argv[2]),atof(argv[3]),atof(argv[4]));

  s = gsl_multimin_fminimizer_alloc(T,5);
  gsl_multimin_fminimizer_set(s,&minex_func,x,step_size);

  i=0;
  do {
    i++;
    status=gsl_multimin_fminimizer_iterate(s);

    if(status) break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size,eps);

    if(status == GSL_SUCCESS) {
      printf("converged in %d iterations\n", i);
    }

  } while(status==GSL_CONTINUE && i<40000);

  if(status!=GSL_SUCCESS) {
    printf("Failed to converge after %d iterations\n", i);
    exit(1);
  }

  fpout = fopen(argv[10],"w");
  fprintf(fpout,"# Tau:   %g\n",gsl_vector_get(s->x,0));
  fprintf(fpout,"# Vlsr:  %g\n",gsl_vector_get(s->x,1));
  fprintf(fpout,"# Vin:   %g\n",gsl_vector_get(s->x,2));
  fprintf(fpout,"# sigma: %g\n",gsl_vector_get(s->x,3));
  fprintf(fpout,"# Tr:    %g\n",gsl_vector_get(s->x,4));
  fprintf(fpout,"# Attained Chisq: %g\n",twolayer5_gsl(s->x,NULL));

  model_spectrum = twolayer5_getfit();

  for(i=0;i<nchan;i++) {
    fprintf(fpout,"%g\t%g\t%g\n", input_data[0][i],input_data[1][i],model_spectrum[i]);
  }
  fclose(fpout);

  gsl_vector_free(x);
  gsl_vector_free(step_size);
  gsl_multimin_fminimizer_free(s);
  twolayer5_free();
  free(input_data[0]);
  free(input_data[1]);

  exit(0);
}
