#include "devo2.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

void devo2_init(devo2_struct *st, int num_param, double *min_param, double *max_param, int num_vecs, double mutation_prob, double diff_factor, double (*evaluate)(double *param_vector)) {
  int i, j;

  /* Check num_vecs Argument */
  if(num_vecs<4) {
    fprintf(stderr, "Too few parameter vectors to evolve.\n");
    exit(1);
  }

  /* Assign some basic values */
  st->num_param=num_param;
  st->num_vecs=num_vecs;
  st->mutation_prob = mutation_prob;
  st->diff_factor = diff_factor;
  st->gen_count=1;
  st->minimization_function = evaluate;
  st->best_vector=NULL;
  st->best_score=DBL_MAX;

  /* Allocate memory and assign pointers */
  st->param_matrix1=(double*)malloc(num_param*num_vecs*sizeof(double));
  st->param_matrix2=(double*)malloc(num_param*num_vecs*sizeof(double));
  st->scores = (double*)malloc(num_vecs*sizeof(double));
  st->trial = (double*)malloc(num_param*sizeof(double));
  st->prev_gen=st->param_matrix1;
  st->next_gen=st->param_matrix2;

  /* Check to see if memory allocation succeeded */
  if(st->param_matrix1==NULL || st->param_matrix2==NULL || st->scores==NULL || st->trial==NULL) {
    fprintf(stderr, "ERROR in devo2_init: Unable to allocate sufficient memory\n");
    exit(1);
  }

  /* Initialize Random Number Generator */
  gsl_rng_env_setup();
  st->rng = gsl_rng_alloc(gsl_rng_default);

  /* Initialize parent generation */
  for(i=0;i<num_vecs;i++) {
    for(j=0;j<num_param;j++) {
      st->prev_gen[i*num_param+j]=gsl_rng_uniform(st->rng)*(max_param[j]-min_param[j])+min_param[j];
    }
    st->scores[i]=(*st->minimization_function)(&(st->prev_gen[i*num_param]));
    if(st->scores[i]<st->best_score) {
      st->best_score = st->scores[i];
      st->best_vector = &(st->prev_gen[i*num_param]);
    }
  }

  return;
}

void devo2_step(devo2_struct *st) {
  /* Begin Differential Evolution */
  int i, j, k;
  int a, b, c;
  double score;
  double *temp;

  for(i=0;i<st->num_vecs;i++) {
    /* Select other parent vectors for mutation */
    do a=(int)gsl_rng_uniform_int(st->rng,(long int)st->num_vecs); while(a==i);
    do b=(int)gsl_rng_uniform_int(st->rng,(long int)st->num_vecs); while(b==i || b==a);
    do c=(int)gsl_rng_uniform_int(st->rng,(long int)st->num_vecs); while(c==i || c==a || c==b);

    /* Randomly pick parameter to adjust */
    j=(int)gsl_rng_uniform_int(st->rng,(long int)st->num_param);

    /* Create a Child Parameter Vector */
    for(k=1;k<=st->num_param;k++) {
      if(gsl_rng_uniform(st->rng) < st->mutation_prob || k==st->num_param) {
	/* Mutate parameter */
	st->trial[j]=st->prev_gen[c*st->num_param+j]+st->diff_factor*(st->prev_gen[a*st->num_param+j]-st->prev_gen[b*st->num_param+j]);
      }
      else {
	/* Copy parameter */
	st->trial[j]=st->prev_gen[i*st->num_param+j];
      }
      j=(j+1)%st->num_param;	/* Started on random param, must use mod */
    }

    /* Evaluate Child, and retain either child or parent or exit */
    score=(*st->minimization_function)(st->trial);

    /* Evaluate whether to pass parent or child to next generation */
    if(score<=st->scores[i]) {
      for(j=0;j<st->num_param;j++) {
	st->next_gen[i*st->num_param+j]=st->trial[j];
      }
      st->scores[i]=score;
      /* Check if it is best score */
      if(score<st->best_score) {
	st->best_score = score;
	st->best_vector = &(st->next_gen[i*st->num_param]);
      }
    }
    else {
      for(j=0;j<st->num_param;j++) {
	st->next_gen[i*st->num_param+j]=st->prev_gen[i*st->num_param+j];
      }
    }
  }

  /* Swap Generations */
  temp=st->next_gen;
  st->next_gen=st->prev_gen;
  st->prev_gen=temp;
  st->gen_count++;
  return;
}

void devo2_free(devo2_struct *st) {
  free(st->param_matrix1);
  free(st->param_matrix2);
  free(st->scores);
  free(st->trial);
  gsl_rng_free(st->rng);
}
