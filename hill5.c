#include "hill5.h"
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
static double *hill_array;
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

void hill5_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax) {
  hill_array = malloc(channels*sizeof(double));
  if(hill_array==NULL) {
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

void hill5_free() {
  free(hill_array);
}

double *hill5_getfit() {
  return hill_array;
}

double hill5_model(double tau,double v_lsr, double v_in, double sigma, double tpeak) {
  double tauf;
  double taur;
  int i;
  double vr;
  double vf;
  double resrms;
  double resid;
  double subf;
  double subr;

  vf = v_lsr+v_in;
  vr = v_lsr-v_in;
  resrms=0.0;
  for(i=0;i<nchan;i++) {
    tauf = tau*exp(-pow((velocity_array[i]-vf)/sigma,2.0)/2.0);
    taur = tau*exp(-pow((velocity_array[i]-vr)/sigma,2.0)/2.0);

    if(tauf>1.0e-4) {
      subf = (1.0-exp(-tauf))/tauf;
    }
    else {
      subf = 1.0;
    }
    if(taur>1.0e-4) {
      subr = (1.0-exp(-taur))/taur;
    }
    else {
      subr = 1.0;
    }

    hill_array[i]=(jfunc(tpeak,frequency)-jfunc(TBG,frequency))*(subf-exp(-tauf)*subr);
    resid=temperature_array[i]-hill_array[i];
    if(velocity_array[i]>vrange[0] && velocity_array[i]<vrange[1]) {
      resrms+=resid*resid;
    }
  }

  return resrms;
}

/* Solvable 5 parameter hill model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: v_in
   3: sigma
   4: Tpeak
*/

double hill5_evaluate(double *params) {
  double tau = params[0];
  double v_lsr = params[1];
  double v_in = params[2];
  double sigma = params[3];
  double tpeak = params[4];
  double resrms;

  if(sigma<0.0) {
    return DBL_MAX;
  }
  if(v_in<0.0) {
    return DBL_MAX;
  }
  if(v_lsr<lsrrange[0] || v_lsr>lsrrange[1]) {
    return DBL_MAX;
  }
  if(v_in>NSIG*sigma) {
    return DBL_MAX;
  }
  if(tau<0.0) {
    return DBL_MAX;
  }

  resrms = hill5_model(tau,v_lsr,v_in,sigma,tpeak);

  return resrms;
}
