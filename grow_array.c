#include "grow_array.h"

static double *ga_location(ga_array *ga_array_ptr, size_t index) {
  size_t skips = index/GARRAYBLOCK;
  size_t subindex = index%GARRAYBLOCK;
  size_t i;
  ga_array *subarray;

  subarray = ga_array_ptr;

  for(i=0;i<skips;i++) {
    if(subarray->next==NULL) {
      subarray->next = malloc(sizeof(ga_array));
      if(subarray->next==NULL) return NULL;
      ga_init(subarray->next);
    }
      subarray=subarray->next;
  }

  return &(subarray->data[subindex]);
}

double ga_get(ga_array *array_ptr, size_t index) {
  double *loc;

  loc = ga_location(array_ptr,index);

  return *loc;
}

int ga_set(ga_array *array_ptr, size_t index, double data) {
  double *loc;

  loc = ga_location(array_ptr,index);

  if(loc==NULL) return -1;

  *loc = data;

  return 0;
}

double *ga_regularize(ga_array *array_ptr, size_t length) {
  double *result;
  size_t i;

  result = malloc(length*sizeof(double));
  if(result==NULL) return NULL;

  for(i=0;i<length;i++) {
    result[i] = ga_get(array_ptr,i);
  }

  return result;
}

void ga_free(ga_array *array_ptr) {
  ga_array *thisone;
  ga_array *nextone;

  thisone = array_ptr->next;

  while(thisone!=NULL) {
    nextone=thisone->next;
    free(thisone);
    thisone=nextone;
  }
}
