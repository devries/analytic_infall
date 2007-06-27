#include <stdio.h>
#include <stdlib.h>

extern int read_data(FILE *fpin, int ncol, int maxlen, double *retdata[]);

int main() {
  FILE *fpin;
  double **bucket;
  int len;
  int i;

  fpin=fopen("read_data_test.dat","r");
  bucket = malloc(2*sizeof(double*));
  len=read_data(fpin,2,132,bucket);
  fclose(fpin);

  for(i=0;i<len;i++) {
    printf("%lf*%lf=%lf\n",bucket[0][i],bucket[1][i],bucket[0][i]*bucket[1][i]);
  }

  exit(0);
}
