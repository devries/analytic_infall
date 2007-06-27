#include <stdio.h>
#include <stdlib.h>
#include "devo2.h"
#include "hill5.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#define MAXITER 40000

extern int read_data(FILE *fpin, int ncol, int maxlen, double *retdata[]);

double hill5_gsl(const gsl_vector *x, void *junk) {
  return hill5_evaluate(x->data);
}

int main(int argc, char *argv[]) {
  FILE *fpin;
  FILE *fpout;
  int popingen;
  int genpercheck;
  int checkperconv;
  double *input_data[2];
  int nchan;
  double min[5];
  double max[5];
  double devo_fit[5];
  devo2_struct dstruct;
  int times_attained;
  double min_attained;
  double *model_spectrum;
  int i;
  gsl_vector *x;
  gsl_vector *step_size;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  double size;
  const double eps = 1.0e-6;
  int status;
  double vmax;
  double vmin;

  if(argc!=9) {
    fprintf(stderr, "Usage: %s <inputfilename> <frequency> <vmin> <vmax> <popingeneration> <generationspercheck> <checksperconv> <outputfile>\n", argv[0]);
    exit(1);
  }

  popingen = atoi(argv[5]);
  genpercheck = atoi(argv[6]);
  checkperconv = atoi(argv[7]);
  vmin=atof(argv[3]);
  vmax=atof(argv[4]);

  fpin = fopen(argv[1],"r");
  nchan = read_data(fpin,2,132,input_data);
  fclose(fpin);

  /* tau range */
  min[0] = 0.1;
  max[0] = 15.0;

  /* vlsr range */
  min[1] = vmin+2.0*(vmax-vmin)/6.0;
  max[1] = vmax-2.0*(vmax-vmin)/6.0;

  /* vin range */
  min[2] = 0.01;
  max[2] = (vmax-vmin)/3.0;

  /* sigma range */
  min[3] = 0.05;
  max[3] = (vmax-vmin)/3.0;

  /* tpeak range */
  min[4] = 2.75;
  max[4] = 40.0;

  hill5_init(nchan,input_data[0],input_data[1],atof(argv[2]),vmin,vmax);

  devo2_init(&dstruct,5,min,max,popingen,0.2,0.8,hill5_evaluate);

  times_attained=1;
  min_attained=dstruct.best_score;

  printf("Initial Result: %lf\n", dstruct.best_score);
  while(times_attained<checkperconv) {
    for(i=0;i<genpercheck;i++) {
      devo2_step(&dstruct);
    }

    printf("Attained Chi Squared: %lf\n",dstruct.best_score);
    printf("tau: %lg, Vlsr: %lg, Vin: %lg,\nsigma: %lg, tpeak: %lg\n", dstruct.best_vector[0], dstruct.best_vector[1], dstruct.best_vector[2], dstruct.best_vector[3], dstruct.best_vector[4]);


    if((min_attained-dstruct.best_score)/min_attained<1.0e-4) {
      min_attained=dstruct.best_score;
      times_attained++;
    }
    else {
      min_attained=dstruct.best_score;
      times_attained=1;
    }
  }

  printf("Switching to simplex method.\n");
  x = gsl_vector_alloc(5);
  gsl_vector_set(x,0,dstruct.best_vector[0]);
  gsl_vector_set(x,1,dstruct.best_vector[1]);
  gsl_vector_set(x,2,dstruct.best_vector[2]);
  gsl_vector_set(x,3,dstruct.best_vector[3]);
  gsl_vector_set(x,4,dstruct.best_vector[4]);
  
  devo2_free(&dstruct);

  step_size = gsl_vector_alloc(5);
  gsl_vector_set(step_size,0,0.1);
  gsl_vector_set(step_size,1,0.005);
  gsl_vector_set(step_size,2,0.005);
  gsl_vector_set(step_size,3,0.005);
  gsl_vector_set(step_size,4,0.05);

  minex_func.f = hill5_gsl;
  minex_func.n = 5;
  minex_func.params = NULL;

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

  } while(status==GSL_CONTINUE && i<MAXITER);

  if(status!=GSL_SUCCESS) {
    printf("Failed to converge after %d iterations\n", i);
    exit(1);
  }

  fpout = fopen(argv[8],"w");
  fprintf(fpout,"# Tau:   %g\n",gsl_vector_get(s->x,0));
  fprintf(fpout,"# Vlsr:  %g\n",gsl_vector_get(s->x,1));
  fprintf(fpout,"# Vin:   %g\n",gsl_vector_get(s->x,2));
  fprintf(fpout,"# sigma: %g\n",gsl_vector_get(s->x,3));
  fprintf(fpout,"# Tpeak: %g\n",gsl_vector_get(s->x,4));
  fprintf(fpout,"# Attained Chisq: %g\n",hill5_gsl(s->x,NULL));

  model_spectrum = hill5_getfit();

  for(i=0;i<nchan;i++) {
    fprintf(fpout,"%g\t%g\t%g\n", input_data[0][i],input_data[1][i],model_spectrum[i]);
  }
  fclose(fpout);

  gsl_vector_free(x);
  gsl_vector_free(step_size);
  gsl_multimin_fminimizer_free(s);
  hill5_free();
  free(input_data[0]);
  free(input_data[1]);

  exit(0);
}
