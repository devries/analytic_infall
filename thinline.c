#include "thinline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define H 6.6260755e-27
#define K 1.380658e-16
#define TBG 2.73
#define NSIG 2.0

static double *velocity_array;
static double *temperature_array;
static double *thinline_array;
static double frequency;
static int nchan;
static double vrange[2];
static double lsrrange[2];

static double jfunc(double t, double nu) {
  double to;

  if(nu<1.0e-6) return t;

  to = H*nu/K;
  return to/(exp(to/t)-1.0);
}

void thinline_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax) {
  thinline_array = malloc(channels*sizeof(double));
  if(thinline_array==NULL) {
    fprintf(stderr, "hill5_init: Out of memory.\n");
    exit(1);
  }

  velocity_array=varray;
  temperature_array=tarray;
  nchan = channels;
  frequency = nu;
  vrange[0]=vmin;
  vrange[1]=vmax;
  lsrrange[0]=vmin+2.0*(vmax-vmin)/6.0;
  lsrrange[1]=vmax-2.0*(vmax-vmin)/6.0;
}

void thinline_free() {
  free(thinline_array);
}

double *thinline_getfit() {
  return thinline_array;
}

double thinline_model(double tau, double v_lsr, double sigma, double tex) {
  int i;
  double resrms;
  double resid;
  double tauc;      /* Optical Depth at Channel */

  resrms=0.0;
  for(i=0;i<nchan;i++) {
    tauc = tau*exp(-pow((velocity_array[i]-v_lsr)/sigma,2.0)/2.0);

    thinline_array[i]=(jfunc(tex,frequency)-jfunc(TBG,frequency))*(1.0-exp(-tauc));
    resid=temperature_array[i]-thinline_array[i];
    if(velocity_array[i]>vrange[0] && velocity_array[i]<vrange[1]) {
      resrms+=resid*resid;
    }
  }

  return resrms;
}

/* Solvable 4 parameter thin line model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: sigma
   3: Tex
*/

double thinline_evaluate(double *params) {
  double tau = params[0];
  double v_lsr = params[1];
  double sigma = params[2];
  double tex = params[3];
  double resrms;

  if(sigma<0.0) {
    return DBL_MAX;
  }
  if(v_lsr<lsrrange[0] || v_lsr>lsrrange[1]) {
    return DBL_MAX;
  }
  if(tau<0.0) {
    return DBL_MAX;
  }

  resrms = thinline_model(tau,v_lsr,sigma,tex);

  return resrms;
}
