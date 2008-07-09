#include "hfsline.h"
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
static double *hfsline_array;
static double frequency;
static int nchan;
static double vrange[2];
static double lsrrange[2];
static int n_components; /* Number of hyperfine components */
static double *comp_voff_array; /* hyperfine component velocities */
static double *comp_relint_array; /* relative component intensity */
static double thin_tex; /* Thin line excitation temperature */

static double jfunc(double t, double nu) {
  double to;

  if(nu<1.0e-6) return t;

  to = H*nu/K;
  return to/(exp(to/t)-1.0);
}

void hfsline_init(int channels, double *varray, double *tarray, double nu, double vmin, double vmax, int ncomp, double *comp_voff, double *comp_relint, double default_tex) {
  int i;
  double vcompmin;
  double vcompmax;
  double vcompint;

  hfsline_array = malloc(channels*sizeof(double));
  if(hfsline_array==NULL) {
    fprintf(stderr, "hill5_init: Out of memory.\n");
    exit(1);
  }

  velocity_array=varray;
  temperature_array=tarray;
  nchan = channels;
  frequency = nu;
  vrange[0]=vmin;
  vrange[1]=vmax;

  vcompmin=comp_voff[0];
  vcompmax=comp_voff[0];
  for(i=1;i<ncomp;i++) {
    if(comp_voff[i]<vcompmin) {
      vcompmin = comp_voff[i];
    }
    else if(comp_voff[i]>vcompmax) {
      vcompmax = comp_voff[i];
    }
  }
  vcompint = vcompmax-vcompmin;

  lsrrange[0]=vmin-vcompmin+(vmax-vmin-vcompint);
  lsrrange[1]=vmax-vcompmax-(vmax-vmin-vcompint);

  n_components = ncomp;
  comp_voff_array = comp_voff;
  comp_relint_array = comp_relint;

  thin_tex = default_tex;
}

void hfsline_free() {
  free(hfsline_array);
}

double *hfsline_getfit() {
  return hfsline_array;
}

double hfsline_model(double tau, double v_lsr, double sigma, double tex) {
  int i;
  int j;
  double resrms;
  double resid;
  double tauc;      /* Optical Depth at Channel */

  resrms=0.0;
  for(i=0;i<nchan;i++) {
    tauc = 0.0;
    for(j=0;j<n_components;j++) {
      tauc += tau*comp_relint_array[j]*exp(-pow((velocity_array[i]-v_lsr-comp_voff_array[j])/sigma,2.0)/2.0);
    }

    hfsline_array[i]=(jfunc(tex,frequency)-jfunc(TBG,frequency))*(1.0-exp(-tauc));
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

double hfsline_evaluate(double *params) {
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

  resrms = hfsline_model(tau,v_lsr,sigma,tex);

  return resrms;
}

double hfsthinline_evaluate(double *params) {
  double tau = params[0];
  double v_lsr = params[1];
  double sigma = params[2];
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

  /* Use static excitation of thin line */
  resrms = hfsline_model(tau,v_lsr,sigma,thin_tex);

  return resrms;
}
