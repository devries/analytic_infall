#include "hyperfine_model.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define H 6.6260755e-27
#define K 1.380658e-16
#define TBG 2.73
#define NSIG 2.0

static double jfunc(double t, double nu) {
  double to;

  if(nu<1.0e-6) return t;

  to = H*nu/K;
  return to/(exp(to/t)-1.0);
}

double get_min_lsr(hyperfine_struct *st) {
    return st->lsrrange[0];
}

double get_max_lsr(hyperfine_struct *st) {
    return st->lsrrange[1];
}

void hyperfine_init(hyperfine_struct *st, int channels, double *varray, double *tarray, double nu, double vmin, double vmax, int ncomp, double *comp_voff, double *comp_relint) {
  int i;
  double vcompmin;
  double vcompmax;
  double vcompint;

  st->hyperfine_array = malloc(channels*sizeof(double));
  if(st->hyperfine_array==NULL) {
    fprintf(stderr, "hyperfine_init: Out of memory.\n");
    exit(1);
  }

  st->velocity_array=varray;
  st->temperature_array=tarray;
  st->nchan = channels;
  st->frequency = nu;
  st->vrange[0]=vmin;
  st->vrange[1]=vmax;

  vcompmin=comp_voff[0];
  vcompmax=comp_voff[0];
  for(i=1;i<ncomp;i++) {
    if(comp_voff[i]<vcompmin) {
      vcompmin=comp_voff[i];
    }
    else if(comp_voff[i]>vcompmax) {
      vcompmax=comp_voff[i];
    }
  }
  vcompint = vcompmax-vcompmin;
  /* CHANGE THESE VALUES */
  /*
  st->lsrrange[0]=vmin-vcompmin+2.0*(vmax-vmin-vcompint)/6.0;
  st->lsrrange[1]=vmax-vcompmax-2.0*(vmax-vmin-vcompint)/6.0;
  */
  st->lsrrange[0]=vmin;
  st->lsrrange[1]=vmax;

  st->n_components = ncomp;
  st->comp_voff_array = comp_voff;
  st->comp_relint_array = comp_relint;
}

void hyperfine_free(hyperfine_struct *st) {
  free(st->hyperfine_array);
}

double *hyperfine_getfit(hyperfine_struct *st) {
  return st->hyperfine_array;
}

double hyperfine_model(hyperfine_struct *st, double tau,double v_lsr, double sigma, double tex) {
  int i, j;
  double tauc; /* Optical depth at a particular channel */
  double resrms;
  double resid;

  resrms=0.0;
  for(i=0;i<st->nchan;i++) {
    tauc=0.0;
    for(j=0;j<st->n_components;j++) {

      tauc += tau*st->comp_relint_array[j]*exp(-pow((st->velocity_array[i]-v_lsr-st->comp_voff_array[j])/sigma,2.0)/2.0);
    }

    st->hyperfine_array[i]=(jfunc(tex,st->frequency)-jfunc(TBG,st->frequency))*(1.0-exp(-tauc));
    resid=st->temperature_array[i]-st->hyperfine_array[i];
    if(st->velocity_array[i]>st->vrange[0] && st->velocity_array[i]<st->vrange[1]) {
      resrms+=resid*resid;
    }
  }

  return resrms;
}

/* Solvable 4 parameter hyperfine model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: sigma
   3: Tex
*/

double hyperfine_evaluate(hyperfine_struct *st, double *params) {
  double tau = params[0];
  double v_lsr = params[1];
  double sigma = params[2];
  double tex = params[3];
  double resrms;

  if(sigma<0.0) {
    return DBL_MAX;
  }
  if(v_lsr<st->lsrrange[0] || v_lsr>st->lsrrange[1]) {
    return DBL_MAX;
  }
  if(tau<0.0) {
    return DBL_MAX;
  }

  resrms = hyperfine_model(st,tau,v_lsr,sigma,tex);

  return resrms;
}
