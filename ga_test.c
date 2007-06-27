#include "grow_array.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
  ga_array g;
  int i;
  double *regarr;

  ga_init(&g);

  for(i=0;i<75;i++) {
    ga_set(&g,i,(double)i);
  }

  for(i=2;i<74;i+=2) {
    printf("%d = %lf\n", i, ga_get(&g,i));
  }

  regarr=ga_regularize(&g,75);

  printf("Element 25: %lf\n", regarr[25]);

  ga_free(&g);

  exit(0);
}

