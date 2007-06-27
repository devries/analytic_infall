#include "hill6core.h"
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

void hill6core_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax) {
  hill_array = malloc(channels*sizeof(double));
  if(hill_array==NULL) {
    fprintf(stderr, "hill6core_init: Out of memory.\n");
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

void hill6core_free() {
  free(hill_array);
}

double *hill6core_getfit() {
  return hill_array;
}

double hill6core_model(double tau_core, double tau_env, double v_lsr, double v_in, double sigma, double tpeak) {
  double tauFR;
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
    tauFR = tau_env*exp(-pow((velocity_array[i]-v_lsr)/sigma,2.0)/2.0);
    tauf = tau_core*exp(-pow((velocity_array[i]-vf)/sigma,2.0)/2.0);
    taur = tau_core*exp(-pow((velocity_array[i]-vr)/sigma,2.0)/2.0);


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

    hill_array[i]=(jfunc(tpeak,frequency)-jfunc(TBG,frequency))*(subf-exp(-tauf)*subr)*exp(-tauFR);
    resid=temperature_array[i]-hill_array[i];
    if(velocity_array[i]>vrange[0] && velocity_array[i]<vrange[1]) {
      resrms+=resid*resid;
    }
  }

  return resrms;
}

/* Solvable 6 parameter hill model with infall in the core, calculates fit
   and returns sumsquare residual.

   parameters:
   0: tau_core
   1: tau_env
   2: v_lsr
   3: v_in -- infall in the core
   4: sigma
   5: Tpeak
*/

double hill6core_evaluate(double *params) {
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

  resrms = hill6core_model(tau_core,tau_env,v_lsr,v_in,sigma,tpeak);

  return resrms;
}
