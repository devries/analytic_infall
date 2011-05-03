#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include "devo2.h"
#include "hyperfine_model.h"
#include "thinline.h"
#include "fitsio.h"
#include <string.h>

#define MAXITER 40000

static hyperfine_struct hyperfine_st;
static hyperfine_struct thin_st;
static double tex_thin; /* Thin line excitation temperature */

double hyperfine_gsl(const gsl_vector *x, void *junk) {
  return hyperfine_evaluate(&hyperfine_st,x->data);
}

double thinline_gsl(const gsl_vector *x, void *junk) {
  double params[4];
  params[0]=x->data[0];
  params[1]=x->data[1];
  params[2]=x->data[2];
  params[3]=tex_thin; /* give the actual number */

  return hyperfine_evaluate(&thin_st,params);
}

double hyperfine_devo(double *params) {
  return hyperfine_evaluate(&hyperfine_st,params);
}

void checkfits(int stat) {
  if(stat) {
    fits_report_error(stderr,stat);
    exit(stat);
  }
  return;
}

int main(int argc, char *argv[]) {
  fitsfile *fitsin;
  fitsfile *fitsout;
  fitsfile *tauout;
  fitsfile *vlsrout;
  fitsfile *sigmaout;
  fitsfile *texout;
  fitsfile *chisqout;
  char fitsoutname[80];
  char tauoutname[80];
  char vlsroutname[80];
  char sigmaoutname[80];
  char texoutname[80];
  char chisqoutname[80];
  char historyline[80];
  int status, nfound, anynull;
  int status_thin, status_hyperfine;
  int naxis;
  long *naxes;
  char card[80];
  static char *includes[] = {"*"};
  static char *excludes[] = {"NAXIS","BITPIX","NAXIS#","BSCALE","BZERO","BLANK","SIMPLE","EXTEND"};
  static int nexcludes = 8;
  double nullval=-1.0e+20;
  int popingen;
  int genpercheck;
  int checkperconv;
  int nchan;
  double min[4];
  double max[4];
  double devo_fit[4];
  double *velocity_spectrum;
  double *input_spectrum;
  devo2_struct dstruct;
  double frequency;
  int times_attained;
  double min_attained;
  double *model_spectrum;
  double chisq;
  long i, j, k;			/* counting variables */
  long crpix;
  double cdelt, crval;
  long fpixel[3];
  long lpixel[3];
  long inc[] = {1L,1L,1L};
  int anynul;
  gsl_vector *x_hyperfine;
  gsl_vector *x_thin;
  gsl_vector *step_size_hyperfine;
  gsl_vector *step_size_thin;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s_hyperfine = NULL;
  gsl_multimin_fminimizer *s_thin = NULL;
  gsl_multimin_function minex_func_hyperfine;
  gsl_multimin_function minex_func_thin;
  double size;
  const double eps = 1.0e-6;
  double vmax;
  double vmin;
  double minpk;
  int notnoise;
  double zero=0.0;
  int ncomponents;
  double *component_voffs;
  double *component_relints;
  char line[132];
  char *token;
  FILE *fpin;

  if(argc!=14) {
    fprintf(stderr,"Usage: %s <infile> <hyperfinefile> <frequency> <thin_tex> <vmin> <vmax> <popingen> <genpercheck> <checkperconv> <x1> <y1> <minpeak> <outbase>\n",argv[0]);
    exit(1);
  }

  status = 0;
  checkfits(fits_open_file(&fitsin,argv[1],READONLY,&status));

  checkfits(fits_get_img_dim(fitsin,&naxis,&status));

  if(naxis!=3) {
    fprintf(stderr,"Expecting 3-D input map.");
    exit(1);
  }

  naxes=(long*)malloc(naxis*sizeof(long));
  if(naxes==NULL) {
    fprintf(stderr, "Unable to allocate space to hold axis SIZES. How unusual.\n");
    exit(1);
  }

  checkfits(fits_get_img_size(fitsin,naxis,naxes,&status));

  sprintf(fitsoutname,"!%s-fits.fits",argv[13]);
  sprintf(tauoutname,"!%s-tau.fits",argv[13]);
  sprintf(vlsroutname,"!%s-vlsr.fits",argv[13]);
  sprintf(sigmaoutname,"!%s-sigma.fits",argv[13]);
  sprintf(texoutname,"!%s-tex.fits",argv[13]);
  sprintf(chisqoutname,"!%s-chisq.fits",argv[13]);

  checkfits(fits_create_file(&fitsout,fitsoutname,&status));
  checkfits(fits_create_file(&tauout,tauoutname,&status));
  checkfits(fits_create_file(&vlsrout,vlsroutname,&status));
  checkfits(fits_create_file(&sigmaout,sigmaoutname,&status));
  checkfits(fits_create_file(&texout,texoutname,&status));
  checkfits(fits_create_file(&chisqout,chisqoutname,&status));

  checkfits(fits_create_img(fitsout,FLOAT_IMG,3,naxes,&status));
  checkfits(fits_create_img(tauout,FLOAT_IMG,2,naxes,&status));
  checkfits(fits_create_img(vlsrout,FLOAT_IMG,2,naxes,&status));
  checkfits(fits_create_img(sigmaout,FLOAT_IMG,2,naxes,&status));
  checkfits(fits_create_img(texout,FLOAT_IMG,2,naxes,&status));
  checkfits(fits_create_img(chisqout,FLOAT_IMG,2,naxes,&status));

  /* Move to beginning of input header */
  checkfits(fits_read_record(fitsin,0,card,&status));

  while(fits_find_nextkey(fitsin,includes,1,excludes,nexcludes,card,&status)==0) {
    checkfits(fits_write_record(fitsout,card,&status));
    checkfits(fits_write_record(tauout,card,&status));
    checkfits(fits_write_record(vlsrout,card,&status));
    checkfits(fits_write_record(sigmaout,card,&status));
    checkfits(fits_write_record(texout,card,&status));
    checkfits(fits_write_record(chisqout,card,&status));
  }
  if(status==KEY_NO_EXIST) {
    status=0;
  }
  else {
    checkfits(status);
  }

  /* Delete axis 3 headers for 2-D images */
  checkfits(fits_delete_key(tauout,"CTYPE3",&status));
  checkfits(fits_delete_key(vlsrout,"CTYPE3",&status));
  checkfits(fits_delete_key(sigmaout,"CTYPE3",&status));
  checkfits(fits_delete_key(texout,"CTYPE3",&status));
  checkfits(fits_delete_key(chisqout,"CTYPE3",&status));
  checkfits(fits_delete_key(tauout,"CRVAL3",&status));
  checkfits(fits_delete_key(vlsrout,"CRVAL3",&status));
  checkfits(fits_delete_key(sigmaout,"CRVAL3",&status));
  checkfits(fits_delete_key(texout,"CRVAL3",&status));
  checkfits(fits_delete_key(chisqout,"CRVAL3",&status));
  checkfits(fits_delete_key(tauout,"CDELT3",&status));
  checkfits(fits_delete_key(vlsrout,"CDELT3",&status));
  checkfits(fits_delete_key(sigmaout,"CDELT3",&status));
  checkfits(fits_delete_key(texout,"CDELT3",&status));
  checkfits(fits_delete_key(chisqout,"CDELT3",&status));
  checkfits(fits_delete_key(tauout,"CRPIX3",&status));
  checkfits(fits_delete_key(vlsrout,"CRPIX3",&status));
  checkfits(fits_delete_key(sigmaout,"CRPIX3",&status));
  checkfits(fits_delete_key(texout,"CRPIX3",&status));
  checkfits(fits_delete_key(chisqout,"CRPIX3",&status));

  /* Adjust the units in each output image */
  checkfits(fits_update_key(tauout,TSTRING,"BUNIT","tau","",&status));
  checkfits(fits_update_key(vlsrout,TSTRING,"BUNIT","km/s","vlsr",&status));
  checkfits(fits_update_key(sigmaout,TSTRING,"BUNIT","km/s","sigma",&status));
  checkfits(fits_update_key(texout,TSTRING,"BUNIT","K","Tex",&status));
  checkfits(fits_update_key(chisqout,TSTRING,"BUNIT","rms","chisq",&status));

  /* add history cards to indicate what you did */
  sprintf(historyline,"--- Hyperfine Map Fit ---");
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Frequency: %s",argv[3]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Thin Line Tex: %s",argv[4]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Line Interval: (%s,%s)",argv[5],argv[6]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Popingen: %s",argv[7]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Genpercheck: %s",argv[8]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Checkperconv: %s",argv[9]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Devo position: (%s,%s)",argv[10],argv[11]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  sprintf(historyline,"  Minimum peak channel: %s",argv[12]);
  checkfits(fits_write_history(fitsout,historyline,&status));
  checkfits(fits_write_history(tauout,historyline,&status));
  checkfits(fits_write_history(vlsrout,historyline,&status));
  checkfits(fits_write_history(sigmaout,historyline,&status));
  checkfits(fits_write_history(texout,historyline,&status));
  checkfits(fits_write_history(chisqout,historyline,&status));

  /* construct velocity spectrum */
  velocity_spectrum=(double*)malloc(naxes[2]*sizeof(double));
  input_spectrum = (double*)malloc(naxes[2]*sizeof(double));

  checkfits(fits_read_key(fitsin,TLONG,"CRPIX3",&crpix,NULL,&status));
  checkfits(fits_read_key(fitsin,TDOUBLE,"CDELT3",&cdelt,NULL,&status));
  checkfits(fits_read_key(fitsin,TDOUBLE,"CRVAL3",&crval,NULL,&status));

  /* convert to km/s */
  cdelt/=1.0e+3;
  crval/=1.0e+3;

  for(i=0;i<naxes[2];i++) {
    velocity_spectrum[i] = (double)(i+1-crpix)*cdelt+crval;
  }

  /* Read in initial fit spectrum */
  fpixel[0]=atol(argv[10]);
  fpixel[1]=atol(argv[11]);
  fpixel[2]=1L;
  lpixel[0]=fpixel[0];
  lpixel[1]=fpixel[1];
  lpixel[2]=naxes[2];

  checkfits(fits_read_subset(fitsin,TDOUBLE,fpixel,lpixel,inc,&nullval,input_spectrum,&anynul,&status));

  if(anynul) {
    fprintf(stderr,"Undefined values in input file.\n");
  }

  popingen = atoi(argv[7]);
  genpercheck = atoi(argv[8]);
  checkperconv = atoi(argv[9]);
  vmin=atof(argv[5]);
  vmax=atof(argv[6]);
  minpk = atof(argv[12]);

  /* Read in hyperfine component file */
  fpin = fopen(argv[2],"r");
  if(fgets(line,132,fpin)==NULL) {
      fprintf(stderr, "Improperly formatted hypefine data file (empty).\n");
      exit(1);
  }
  token = strtok(line, " \f\n\r\t\v");
  ncomponents = atoi(token);
  component_voffs = malloc(ncomponents*sizeof(double));
  component_relints = malloc(ncomponents*sizeof(double));
  if(component_voffs == NULL || component_relints == NULL) {
      fprintf(stderr, "ERROR: Out of Memory");
      exit(1);
  }
  for(i=0;i<ncomponents;i++) {
      if(fgets(line,132,fpin)==NULL) {
          fprintf(stderr, "Hyperfine data file too short.\n");
          exit(1);
      }
      token = strtok(line, " \f\n\r\t\v");
      component_voffs[i] = atof(token);
      token = strtok(NULL, " \f\n\r\t\v");
      component_relints[i] = atof(token);
  }
  fclose(fpin);

  frequency = atof(argv[3]);
  tex_thin = atof(argv[4]);
  hyperfine_init(&hyperfine_st,naxes[2],velocity_spectrum,input_spectrum,frequency,vmin,vmax,ncomponents,component_voffs,component_relints);
  
  /* tau range */
  min[0] = 0.1;
  max[0] = 30.0;

  /* vlsr range is calculated in the hyperfine routine */
  min[1] = get_min_lsr(&hyperfine_st);
  max[1] = get_max_lsr(&hyperfine_st);

  /* sigma range */
  min[2] = 0.05;
  max[2] = (vmax-vmin)/3.0;

  /* tex range */
  min[3] = 2.75;
  max[3] = 40.0;


  devo2_init(&dstruct,4,min,max,popingen,0.2,0.3,hyperfine_devo);

  times_attained=1;
  min_attained=dstruct.best_score;

  printf("Initial Chi Squared: %lf\n", dstruct.best_score);
  while(times_attained<checkperconv) {
    for(i=0;i<genpercheck;i++) {
      devo2_step(&dstruct);
    }

    printf("Attained Chi Squared: %lf\n",dstruct.best_score);
    printf("tau: %lg, Vlsr: %lg, sigma: %lg, tex: %lg\n", dstruct.best_vector[0], dstruct.best_vector[1], dstruct.best_vector[2], dstruct.best_vector[3]);


    if((min_attained-dstruct.best_score)/min_attained<1.0e-4) {
      min_attained=dstruct.best_score;
      times_attained++;
    }
    else {
      min_attained=dstruct.best_score;
      times_attained=1;
    }
  }

  printf("Initial results obtained. Fitting Map.\n");

  devo_fit[0]=dstruct.best_vector[0];
  devo_fit[1]=dstruct.best_vector[1];
  devo_fit[2]=dstruct.best_vector[2];
  devo_fit[3]=dstruct.best_vector[3];

  devo2_free(&dstruct);
  hyperfine_free(&hyperfine_st);

  x_hyperfine = gsl_vector_alloc(4);
  step_size_hyperfine = gsl_vector_alloc(4);
  minex_func_hyperfine.f=hyperfine_gsl;
  minex_func_hyperfine.n=4;
  minex_func_hyperfine.params=NULL;

  x_thin = gsl_vector_alloc(3);
  step_size_thin = gsl_vector_alloc(3);
  minex_func_thin.f=thinline_gsl;
  minex_func_thin.n=3;
  minex_func_thin.params=NULL;

  for(i=0;i<naxes[0];i++) {
    for(j=0;j<naxes[1];j++) {
      fpixel[0]=i+1;
      fpixel[1]=j+1;
      lpixel[0]=fpixel[0];
      lpixel[1]=fpixel[1];

      checkfits(fits_read_subset(fitsin,TDOUBLE,fpixel,lpixel,inc,&nullval,input_spectrum,&anynul,&status));

      notnoise=0;
      for(k=0;k<naxes[2];k++) {
        if(velocity_spectrum[k]>vmin && velocity_spectrum[k]<vmax) {
          if(input_spectrum[k]>=minpk) {
            notnoise=1;
            break;
          }
        }
      }

      if(notnoise) {
        hyperfine_init(&thin_st,naxes[2],velocity_spectrum,input_spectrum,frequency,vmin,vmax,ncomponents,component_voffs,component_relints);
        gsl_vector_set(x_thin,0,0.1);
        gsl_vector_set(x_thin,1,devo_fit[1]);
        gsl_vector_set(x_thin,2,devo_fit[2]);
        gsl_vector_set(step_size_thin,0,0.1);
        gsl_vector_set(step_size_thin,1,0.005);
        gsl_vector_set(step_size_thin,2,0.005);

        s_thin = gsl_multimin_fminimizer_alloc(T,3);
        gsl_multimin_fminimizer_set(s_thin,&minex_func_thin,x_thin,step_size_thin);

        k=0;
        do {
          k++;
          status_thin=gsl_multimin_fminimizer_iterate(s_thin);
          if(status_thin) break;
          size = gsl_multimin_fminimizer_size(s_thin);
          status_thin = gsl_multimin_test_size(size,eps);

          if(status_thin == GSL_SUCCESS) {
            printf("(%d,%d) --- Thin line converged in %d iterations\n", fpixel[0],fpixel[1],k);
          }
        } while(status_thin==GSL_CONTINUE && k<40000);

        hyperfine_init(&hyperfine_st,naxes[2],velocity_spectrum,input_spectrum,frequency,vmin,vmax,ncomponents,component_voffs,component_relints);
        
        gsl_vector_set(x_hyperfine,0,devo_fit[0]);
        gsl_vector_set(x_hyperfine,1,devo_fit[1]);
        gsl_vector_set(x_hyperfine,2,devo_fit[2]);
        gsl_vector_set(x_hyperfine,3,devo_fit[3]);
        gsl_vector_set(step_size_hyperfine,0,0.1);
        gsl_vector_set(step_size_hyperfine,1,0.005);
        gsl_vector_set(step_size_hyperfine,2,0.005);
        gsl_vector_set(step_size_hyperfine,3,0.05);

        s_hyperfine = gsl_multimin_fminimizer_alloc(T,4);
        gsl_multimin_fminimizer_set(s_hyperfine,&minex_func_hyperfine,x_hyperfine,step_size_hyperfine);

        k=0;
        do {
          k++;
          status_hyperfine=gsl_multimin_fminimizer_iterate(s_hyperfine);
          if(status_hyperfine) break;
          size = gsl_multimin_fminimizer_size(s_hyperfine);
          status_hyperfine = gsl_multimin_test_size(size,eps);

          if(status_hyperfine == GSL_SUCCESS) {
            printf("(%d,%d) --- Thick Hyperfine converged in %d iterations\n", fpixel[0],fpixel[1],k);
          }
        } while(status_hyperfine==GSL_CONTINUE && k<40000);

        if(status_thin!=GSL_SUCCESS && status_hyperfine!=GSL_SUCCESS) {
          printf("(%d,%d) --- Thin and thick hyperfine failed to converge after %d iterations\n", fpixel[0],fpixel[1],k);
          status=0;
          checkfits(fits_write_pixnull(tauout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
          checkfits(fits_write_pixnull(vlsrout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
          checkfits(fits_write_pixnull(sigmaout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
          checkfits(fits_write_pixnull(texout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
          checkfits(fits_write_pixnull(chisqout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
        }
        else if(status_hyperfine!=GSL_SUCCESS) {
          status = 0;
          checkfits(fits_write_pix(tauout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_thin->x,0),&status));
          checkfits(fits_write_pix(vlsrout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_thin->x,1),&status));
          checkfits(fits_write_pix(sigmaout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_thin->x,2),&status));
          checkfits(fits_write_pix(texout,TDOUBLE,fpixel,1,&tex_thin,&status));
          chisq=thinline_gsl(s_thin->x,NULL);
          checkfits(fits_write_pix(chisqout,TDOUBLE,fpixel,1,&chisq,&status));
          model_spectrum=hyperfine_getfit(&thin_st);
        
          checkfits(fits_write_subset(fitsout,TDOUBLE,fpixel,lpixel,model_spectrum,&status));
        }
        else {
          if(thinline_gsl(s_thin->x,NULL)<hyperfine_gsl(s_hyperfine->x,NULL)) {
            status = 0;
            checkfits(fits_write_pix(tauout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_thin->x,0),&status));
            checkfits(fits_write_pix(vlsrout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_thin->x,1),&status));
            checkfits(fits_write_pix(sigmaout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_thin->x,2),&status));
            checkfits(fits_write_pix(texout,TDOUBLE,fpixel,1,&tex_thin,&status));
            chisq=thinline_gsl(s_thin->x,NULL);
            checkfits(fits_write_pix(chisqout,TDOUBLE,fpixel,1,&chisq,&status));
            model_spectrum=hyperfine_getfit(&thin_st);
          
            checkfits(fits_write_subset(fitsout,TDOUBLE,fpixel,lpixel,model_spectrum,&status));
          }
          else {
            status = 0;
            checkfits(fits_write_pix(tauout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_hyperfine->x,0),&status));
            checkfits(fits_write_pix(vlsrout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_hyperfine->x,1),&status));
            checkfits(fits_write_pix(sigmaout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_hyperfine->x,2),&status));
            checkfits(fits_write_pix(texout,TDOUBLE,fpixel,1,gsl_vector_ptr(s_hyperfine->x,3),&status));
            chisq=hyperfine_gsl(s_hyperfine->x,NULL);
            checkfits(fits_write_pix(chisqout,TDOUBLE,fpixel,1,&chisq,&status));
            model_spectrum=hyperfine_getfit(&hyperfine_st);
          
            checkfits(fits_write_subset(fitsout,TDOUBLE,fpixel,lpixel,model_spectrum,&status));
          }
        }

        gsl_multimin_fminimizer_free(s_thin);
        gsl_multimin_fminimizer_free(s_hyperfine);
        hyperfine_free(&thin_st);
        hyperfine_free(&hyperfine_st);
      }
      else {
        printf("(%d,%d) --- No Signal\n",fpixel[0],fpixel[1]);
        status=0;
        checkfits(fits_write_pixnull(tauout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
        checkfits(fits_write_pixnull(vlsrout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
        checkfits(fits_write_pixnull(sigmaout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
        checkfits(fits_write_pixnull(texout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
        checkfits(fits_write_pixnull(chisqout,TDOUBLE,fpixel,1,&nullval,&nullval,&status));
      }
    }
  }

  gsl_vector_free(x_hyperfine);
  gsl_vector_free(x_thin);
  gsl_vector_free(step_size_hyperfine);
  gsl_vector_free(step_size_thin);

  free(velocity_spectrum);
  free(input_spectrum);
  free(naxes);
  
  checkfits(fits_close_file(fitsin,&status));
  checkfits(fits_close_file(fitsout,&status));
  checkfits(fits_close_file(tauout,&status));
  checkfits(fits_close_file(vlsrout,&status));
  checkfits(fits_close_file(sigmaout,&status));
  checkfits(fits_close_file(texout,&status));
  checkfits(fits_close_file(chisqout,&status));

  exit(0);
}

