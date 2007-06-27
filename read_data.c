#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "grow_array.h"

/*****************************************************************
 ** Reads text file fpin, returns number of data lines, sets retdata to a
 ** two dimensional array of ncol columns and nlines of data. Each line of
 ** the data file must be less than maxlen long. retdata is an initialized
 ** array of pointers to doubles of length ncol. The characters '#' and '!'
 ** as the first character of a line indicate a comment.
 *****************************************************************/

int read_data(FILE *fpin, int ncol, int maxlen, double *retdata[]) {
  int i;
  int nlines=0;
  char *line;
  ga_array *arrays;
  char *token;

  line = malloc(maxlen*sizeof(char));
  arrays = malloc(ncol*sizeof(ga_array));

  if(line==NULL || arrays == NULL) {
    fprintf(stderr, "ERROR in read_data: Out of memory.\n");
    exit(1);
  }

  while(fgets(line,maxlen,fpin)) {
    i=(int)strlen(line);
    if(line[i-1]!='\n') {
      fprintf(stderr, "ERROR in read_data: Data file has lines that are too long.\n");
      exit(1);
    }

    if(line[0]!='#' && line[0]!='!') {
      for(i=0;i<ncol;i++) {
	if(i==0) {
	  token = strtok(line," \f\n\r\t\v");
	  if(token==NULL) {
	    nlines--;
	    break; /* Blank line */
	  }	    
	}
	else {
	  token = strtok(NULL," \f\n\r\t\v");
	  
	  if(token==NULL) {
	    fprintf(stderr, "ERROR in read_data: Data file is not properly formatted.\n");
	    exit(1);
	  }
	}

	if(ga_set(&arrays[i],nlines,atof(token))<0) {
	  fprintf(stderr, "ERROR in read_data: Out of memory.\n");
	  exit(1);
	}
      }
      nlines++;
    }
  }

  for(i=0;i<ncol;i++) {
    retdata[i]=ga_regularize(&arrays[i],nlines);
    ga_free(&arrays[i]);
  }

  free(arrays);
  free(line);

  return nlines;
}
