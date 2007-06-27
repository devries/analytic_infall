#include <stdio.h>
#include <stdlib.h>
#include "devo2.h"
#include "twolayer5.h"

extern int read_data(FILE *fpin, int ncol, int maxlen, double *retdata[]);

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
  devo2_struct dstruct;
  int times_attained;
  double min_attained;
  double *model_spectrum;
  int i;

  if(argc!=9) {
    fprintf(stderr, "Usage: %s <inputfilename> <frequency> <vmin> <vmax> <popingeneration> <generationspercheck> <checksperconv> <outputfile>\n", argv[0]);
    exit(1);
  }

  popingen = atoi(argv[5]);
  genpercheck = atoi(argv[6]);
  checkperconv = atoi(argv[7]);

  fpin = fopen(argv[1],"r");
  nchan = read_data(fpin,2,132,input_data);
  fclose(fpin);

  /* tau range */
  min[0] = 0.1;
  max[0] = 15.0;

  /* vlsr range */
  min[1] = atof(argv[3]);
  max[1] = atof(argv[4]);

  /* vin range */
  min[2] = 0.01;
  max[2] = (max[1]-min[1])/3.0;

  /* sigma range */
  min[3] = 0.05;
  max[3] = (max[1]-min[1])/3.0;

  /* tr range */
  min[4] = 2.75;
  max[4] = 30.0;

  twolayer5_init(nchan,input_data[0],input_data[1],atof(argv[2]),min[1],max[1]);

  devo2_init(&dstruct,5,min,max,popingen,0.2,0.8,twolayer5_evaluate);

  times_attained=1;
  min_attained=dstruct.best_score;

  printf("Initial Result: %lf\n", dstruct.best_score);
  while(times_attained<checkperconv) {
    for(i=0;i<genpercheck;i++) {
      devo2_step(&dstruct);
    }

    printf("Attained: %lf\n",dstruct.best_score);
    printf("tau: %lg, Vlsr: %lg, Vin: %lg,\nsigma: %lg, tr: %lg\n", dstruct.best_vector[0], dstruct.best_vector[1], dstruct.best_vector[2], dstruct.best_vector[3], dstruct.best_vector[4]);

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
  fprintf(fpout,"# Tau:   %g\n",dstruct.best_vector[0]);
  fprintf(fpout,"# Vlsr:  %g\n",dstruct.best_vector[1]);
  fprintf(fpout,"# Vin:   %g\n",dstruct.best_vector[2]);
  fprintf(fpout,"# sigma: %g\n",dstruct.best_vector[3]);
  fprintf(fpout,"# Tr:    %g\n",dstruct.best_vector[4]);
  fprintf(fpout,"# Attained Chisq: %g\n",dstruct.best_score);

  twolayer5_evaluate(dstruct.best_vector);
  model_spectrum = twolayer5_getfit();

  for(i=0;i<nchan;i++) {
    fprintf(fpout,"%g\t%g\t%g\n", input_data[0][i],input_data[1][i],model_spectrum[i]);
  }
  fclose(fpout);

  devo2_free(&dstruct);
  twolayer5_free();
  free(input_data[0]);
  free(input_data[1]);

  exit(0);
}
