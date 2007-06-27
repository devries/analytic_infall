#include <stdio.h>
#include <stdlib.h>
#include "devo2.h"
#include "hill7.h"

extern int read_data(FILE *fpin, int ncol, int maxlen, double *retdata[]);

int main(int argc, char *argv[]) {
  FILE *fpin;
  FILE *fpout;
  int popingen;
  int genpercheck;
  int checkperconv;
  double *input_data[2];
  int nchan;
  double min[7];
  double max[7];
  devo2_struct dstruct;
  int times_attained;
  double min_attained;
  double *model_spectrum;
  int i;
  double vmin;
  double vmax;

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

  /* tau_core range */
  min[0] = 0.1;
  max[0] = 15.0;

  /* tau_env range */
  min[1] = 0.1;
  max[1] = 15.0;

  /* vlsr range */
  min[2] = vmin+2.0*(vmax-vmin)/6.0;
  max[2] = vmax-2.0*(vmax-vmin)/6.0;

  /* vin range */
  min[3] = 0.01;
  max[3] = (vmax-vmin)/3.0;

  /* sigma range */
  min[4] = 0.05;
  max[4] = (vmax-vmin)/3.0;

  /* to range */
  min[5] = 2.75;
  max[5] = 10.0;

  /* tpeak range */
  min[6] = 2.75;
  max[6] = 20.0;

  hill7_init(nchan,input_data[0],input_data[1],atof(argv[2]),vmin,vmax);

  devo2_init(&dstruct,7,min,max,popingen,0.1,0.8,hill7_evaluate);

  times_attained=1;
  min_attained=dstruct.best_score;

  printf("Initial Result: %lf\n", dstruct.best_score);
  while(times_attained<checkperconv) {
    for(i=0;i<genpercheck;i++) {
      devo2_step(&dstruct);
    }

    printf("Attained: %lf\n",dstruct.best_score);
    printf("tau_core: %lg, tau_env: %lg, Vlsr: %lg,\nVin: %lg, sigma: %lg, to: %lg,\ntpeak: %lg\n", dstruct.best_vector[0], dstruct.best_vector[1], dstruct.best_vector[2], dstruct.best_vector[3], dstruct.best_vector[4],dstruct.best_vector[5],dstruct.best_vector[6]);

    if((min_attained-dstruct.best_score)/min_attained<1.0e-4) {
      min_attained=dstruct.best_score;
      times_attained++;
    }
    else {
      min_attained=dstruct.best_score;
      times_attained=1;
    }
  }

  fpout = fopen(argv[8],"w");
  fprintf(fpout,"# Tau_c: %g\n",dstruct.best_vector[0]);
  fprintf(fpout,"# Tau_e: %g\n",dstruct.best_vector[1]);
  fprintf(fpout,"# Vlsr:  %g\n",dstruct.best_vector[2]);
  fprintf(fpout,"# Vin:   %g\n",dstruct.best_vector[3]);
  fprintf(fpout,"# sigma: %g\n",dstruct.best_vector[4]);
  fprintf(fpout,"# To:    %g\n",dstruct.best_vector[5]);
  fprintf(fpout,"# Tpeak: %g\n",dstruct.best_vector[6]);
  fprintf(fpout,"# Attained Chisq: %g\n",dstruct.best_score);

  hill7_evaluate(dstruct.best_vector);
  model_spectrum = hill7_getfit();

  for(i=0;i<nchan;i++) {
    fprintf(fpout,"%g\t%g\t%g\n", input_data[0][i],input_data[1][i],model_spectrum[i]);
  }
  fclose(fpout);

  devo2_free(&dstruct);
  hill7_free();
  free(input_data[0]);
  free(input_data[1]);

  exit(0);
}
