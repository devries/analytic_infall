#include <stdio.h>
#include <stdlib.h>
#include "hill5.h"

extern int read_data(FILE *fpin, int ncol, int maxlen, double *retdata[]);

int main(int argc, char *argv[]) {
  FILE *fpin;
  FILE *fpout;
  double *input_data[2];
  int nchan;
  double params[5];
  double *model_spectrum;
  int i;

  if(argc!=9) {
    fprintf(stderr, "Usage: %s <inputfilename> <frequency> <tau> <vlsr> <vin> <sigma> <tpeak> <outputfile>\n",argv[0]);
    exit(1);
  }

  fpin = fopen(argv[1],"r");
  nchan = read_data(fpin,2,132,input_data);
  fclose(fpin);
  params[0] = atof(argv[3]);
  params[1] = atof(argv[4]);
  params[2] = atof(argv[5]);
  params[3] = atof(argv[6]);
  params[4] = atof(argv[7]);

  hill5_init(nchan,input_data[0],input_data[1],atof(argv[2]),input_data[0][0],input_data[0][nchan-1]);

  hill5_evaluate(params);

  fpout = fopen(argv[8],"w");
  fprintf(fpout,"# Tau:   %g\n",params[0]);
  fprintf(fpout,"# Vlsr:  %g\n",params[1]);
  fprintf(fpout,"# Vin:   %g\n",params[2]);
  fprintf(fpout,"# sigma: %g\n",params[3]);
  fprintf(fpout,"# Tpeak: %g\n",params[4]);
  model_spectrum = hill5_getfit();

  for(i=0;i<nchan;i++) {
    fprintf(fpout,"%g\t%g\n",input_data[0][i],model_spectrum[i]);
  }
  fclose(fpout);

  hill5_free();
  free(input_data[0]);
  free(input_data[1]);

  exit(0);
}
