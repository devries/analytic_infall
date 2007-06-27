/***********************************************************************
	Differential Evolution Genetic Algorythm

This routine minimises a vector of parameters using differential
evolution as described by Kenneth Price and Rainer Storn in Dr. Dobb's
Journal #264, 1997. pps 18-24,78. 

The input parameters are as follows:

int num_param		Number of parameters in vector.

double *min_param	Minimum value of parameters in search space. 

double *max_param	Maximum value of parameters in search space.

int num_vecs		Number of vectors in each generation. Must be
			greater than 4. Suggested: 5*num_param or more.

double mutation_prob	Probability of mutation. At least one
			parameter in the vector will mutate, however
			the others will mutate with this probability.
			Suggested: 0.1 to start, want as large as
			possible.

double diff_factor	Factor by which differential parameter
			distances are modified before adding to
			parent. Suggest: 0.4 to 1.2

double request_min	Requested minimum. If the minimization
			achieves this value or lower the program
			terminates. 

int max_generations	Maximum allowed generations before
			termination.

long *ranseed		Random number seed (using Numerical Recipes
			ran3.c) Should be negative to initialize the
			random number generator.

double *attained_min	Minimization attained by routine.

double *attained_param	Vector of parameters attained in minimization.

double (*evaluate)(double *param_vector)	Pointer to function which
					evaluates parameter space and
					returns value to minimize. 

Routine Return Values:

int			Number of generations performed. If negative
			then the total number of generations was
			exhausted without sufficiently minimizing the 
			function. In this case, the min_param and
			max_param vectors are replaced with the
			minimum and maximum spread of the final
			generation of parameters, and the best
			parameter vector attained as well as its score
			is returned in attained_param and
			attained_min.

************************************************************************/
#ifndef _D_EVO_H
#define _D_EVO_H
#include <gsl/gsl_rng.h>

typedef struct _devo2_struct {
  int num_param;		/* Number of parameters */
  int num_vecs;			/* Number of vectors in a generation */
  double mutation_prob;
  double diff_factor;
  double *param_matrix1;	/* Memory for one generation */
  double *param_matrix2;	/* Memory for another generation */
  double *prev_gen;		/* pointer to previous generation */
  double *next_gen;		/* pointer to next generation */
  double *scores;		/* scores of previous generation */
  int gen_count;		/* number of generatione */
  double (*minimization_function)(double *param_vector); /* min func */
  gsl_rng *rng;    		/* Random number generator */
  double *trial;		/* Some space to put a trial vector */
  double best_score;		/* Score of best vector */
  double *best_vector;		/* Pointer to best vector */
} devo2_struct;

extern void devo2_init(devo2_struct *st, int num_param, double *min_param, double *max_param, int num_vecs, double mutation_prob, double diff_factor, double (*evaluate)(double *param_vector));

extern void devo2_step(devo2_struct *st);

extern void devo2_free(devo2_struct *st);

#endif /*!_D_EVO_H*/
