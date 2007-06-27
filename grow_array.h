#ifndef _GROW_ARRAY_H
#define _GROW_ARRAY_H
#include <stdlib.h>

#define GARRAYBLOCK 32

typedef struct _ga_array {
  double data[GARRAYBLOCK];
  struct _ga_array *next;
} ga_array;

#define ga_init(array_ptr) {(array_ptr)->next = NULL;}

extern double ga_get(ga_array *array_ptr, size_t index);

extern int ga_set(ga_array *array_ptr, size_t index, double data);

extern double *ga_regularize(ga_array *array_ptr, size_t length);

extern void ga_free(ga_array *array_ptr);

#endif /* _GROW_ARRAY_H */
