#include "hill6.h"
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

void hill6_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax) {
  hill_array = malloc(channels*sizeof(double));
  if(hill_array==NULL) {
    fprintf(stderr, "hill6_init: Out of memory.\n");
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

void hill6_free() {
  free(hill_array);
}

double *hill6_getfit() {
  return hill_array;
}

double hill6_model(double tau_core, double tau_env, double v_lsr, double v_in, double sigma, double tpeak) {
  double taufr;
  double tauF;
  /* double tauR; */
  int i;
  /* double vr; */
  double vf;
  double resrms;
  double resid;
  double sub;

  vf = v_lsr+v_in;
  /* vr = v_lsr-v_in; */
  resrms=0.0;
  for(i=0;i<nchan;i++) {
    taufr = tau_core*exp(-pow((velocity_array[i]-v_lsr)/sigma,2.0)/2.0);
    tauF = tau_env*exp(-pow((velocity_array[i]-vf)/sigma,2.0)/2.0);
    /* tauR = tau_env*exp(-pow((velocity_array[i]-vr)/sigma,2.0)/2.0); */


    if(taufr>1.0e-4) {
      sub = (1.0-exp(-taufr))/taufr;
    }
    else {
      sub = 1.0;
    }

    hill_array[i]=(jfunc(tpeak,frequency)-jfunc(TBG,frequency))*sub*(1.0-exp(-taufr))*exp(-tauF);
    resid=temperature_array[i]-hill_array[i];
    if(velocity_array[i]>vrange[0] && velocity_array[i]<vrange[1]) {
      resrms+=resid*resid;
    }
  }

  return resrms;
}

/* Solvable 6 parameter hill model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_core
   1: tau_env
   2: v_lsr
   3: v_in -- infall in the envelope
   4: sigma
   5: Tpeak
*/

double hill6_evaluate(double *params) {
  double tau_core = params[0];
  double tau_env = params[1];
  double v_lsr = params[2];
  double v_in = params[3];
  double sigma = params[4];
  double tpeak = params[5];
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
  if(tau_core<0.0) {
    return DBL_MAX;
  }
  if(tau_env<0.0) {
    return DBL_MAX;
  }

  resrms = hill6_model(tau_core,tau_env,v_lsr,v_in,sigma,tpeak);

  return resrms;
}
