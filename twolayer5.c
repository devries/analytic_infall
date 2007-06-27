#include "twolayer5.h"
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
static double *twolayer_array;
static double frequency;
static int nchan;
static double vrange[2];

static double jfunc(double t, double nu) {
  double to;

  if(nu<1.0e-6) return t;

  to = H*nu/K;
  return to/(exp(to/t)-1.0);
}

void twolayer5_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax) {
  twolayer_array = malloc(channels*sizeof(double));
  if(twolayer_array==NULL) {
    fprintf(stderr, "twolayer5_init: Out of memory.\n");
    exit(1);
  }

  velocity_array=varray;
  temperature_array=tarray;
  nchan = channels;
  frequency = nu;
  vrange[0]=vmin;
  vrange[1]=vmax;
}

void twolayer5_free() {
  free(twolayer_array);
}

double *twolayer5_getfit() {
  return twolayer_array;
}

static double twolayer5_model(double tau, double v_lsr, double v_in, double sigma, double tr) {
  double tauf;
  double taur;
  int i;
  double vr;
  double vf;
  double resrms;
  double resid;

  vf = v_lsr+v_in;
  vr = v_lsr-v_in;
  resrms=0.0;
  for(i=0;i<nchan;i++) {
    tauf = tau*exp(-pow((velocity_array[i]-vf)/sigma,2.0)/2.0);
    taur = tau*exp(-pow((velocity_array[i]-vr)/sigma,2.0)/2.0);

    twolayer_array[i]=(jfunc(tr,frequency)-jfunc(TBG,frequency))*(1.0-exp(-taur))*exp(-tauf);
    resid=temperature_array[i]-twolayer_array[i];
    if(velocity_array[i]>vrange[0] && velocity_array[i]<vrange[1]) {
      resrms+=resid*resid;
    }
  }

  return resrms;
}

/* Solvable 5 parameter two layer model, calculates fit and returns rms
   residual.

   parameters:
   0: tau_0
   1: v_lsr
   2: v_in
   3: sigma
   4: Tr
*/

double twolayer5_evaluate(double *params) {
  double tau = params[0];
  double v_lsr = params[1];
  double v_in = params[2];
  double sigma = params[3];
  double tr = params[4];
  double resrms;

  if(sigma<0.0) {
    return DBL_MAX;
  }
  if(v_in<0.0) {
    return DBL_MAX;
  }
  if(v_lsr<vrange[0] || v_lsr>vrange[1]) {
    return DBL_MAX;
  }
  if(v_in>NSIG*sigma) {
    return DBL_MAX;
  }

  resrms = twolayer5_model(tau,v_lsr,v_in,sigma,tr);

  return resrms;
}

