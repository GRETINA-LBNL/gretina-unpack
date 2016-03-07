#ifndef __RDGECALFILE_H
#define __RDGECALFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>

#include "Unpack.h"

/*************************************************************************/

int RdGeCalFile(const char *fn, GRETINAVariables* gVar) {
 
  /* Declarations */
  FILE *fp;
  int i1, nn;
  float f1, f2;
  char *st, str[128];

  /* Open file */
  fp = fopen(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  }
  printf("\"%s\" open....", fn);

  /* Read values */
  nn = 0;
  st = fgets(str, 64, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* Empty line, do nothing */
    } else {
      sscanf(str, "%i %f %f", &i1, &f1, &f2);
      if (i1>=0 && i1<MAXCHANNELS) {
        gVar->ehiGeOffset[i1]=f1;
        gVar->ehiGeGain[i1]=f2;
	// printf("det %3.3i: offset=%9.3f, gain=%9.6f\n", 
	//         i1, offset[i1], gain[i1]);
        fflush(stdout);
      }
      nn++;
    }
    
    /* Attempt to read next line */
    st = fgets(str, 64, fp);
  }
  printf("Read %i gain calibration coefficients.\n", nn);

  /* Done! */
  fclose(fp);
  return (0);
}

/*************************************************************************/

int RdDNLGeCalFile(const char *fn, float *dnl0, float *dnl1, float *dnl2,
		   float *dnl3, float *dnl4, float *dnl5, float *dnl6,
		   float *dnl7) {
  
  /* Declarations */
  FILE *fp;
  int i1, i2, i3, nn;
  float f1, f2, f3, f4, f5, f6=0, f7=0, f8=0;
  char *st, str[256];

  /* Open file */
  fp = fopen(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  }
  printf("\"%s\" open.\n", fn);
  
  /* Read values */
  nn = 0;
  st = fgets(str, 128, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* empty line, do nothing */
    } else {
      sscanf(str, "%i %i %i %f %f %f %f %f", 
	     &i1, &i2, &i3, &f1, &f2, &f3, &f4, &f5);
      if (i3>=0 && i3<MAXCHANNELS) {
	dnl0[i3] = f1;
	dnl1[i3] = f2;
	dnl2[i3] = f3;
	dnl3[i3] = f4;
	dnl4[i3] = f5;
	dnl5[i3] = f6;
	dnl6[i3] = f7;
	dnl7[i3] = f8;
      }
      nn++;
    }

    /* Attempt to read next line */
    st = fgets(str, 128, fp);
  }
  printf("Read %i DNL correction calibration coefficients.\n", nn);

  /* Done! */
  fclose(fp);
  return (0);
}

/*************************************************************************/

int RdDNLLookUpCalFile(const char *fn, GRETINAVariables* gVar) {
  
  /* Declarations */
  FILE *fp;
  int i1, i2, nn;
  float f[40];

  char *st, str[2048];

  /* Open file */
  fp = fopen(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  }
  printf("\"%s\" open.\n", fn);
  
  /* Read values */
  nn = 0;
  st = fgets(str, 2048, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* empty line, do nothing */
    } else {
      sscanf(str, "%i %i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
	     &i1, &i2, 
	     &f[0], &f[1], &f[2], &f[3], &f[4],
	     &f[5], &f[6], &f[7], &f[8], &f[9],
	     &f[10], &f[11], &f[12], &f[13], &f[14],
	     &f[15], &f[16], &f[17], &f[18], &f[19],
	     &f[20], &f[21], &f[22], &f[23], &f[24],
	     &f[25], &f[26], &f[27], &f[28], &f[29],
	     &f[30], &f[31], &f[32], &f[33], &f[34],
	     &f[35], &f[36], &f[37], &f[38], &f[39]);
      if (i1>0 && i1<=MAXCRYSTALS) {
	for (int m=0; m<40; m++) {
	  int elecID = (int)(((i1-1)*40) + m);
	  gVar->dnlLU[elecID][i2] = f[m];
	  // cout << elecID << " " << i2 << " " << f[m] << endl;
	}
      }
      nn++;
    }

    /* Attempt to read next line */
    st = fgets(str, 2048, fp);
  }
  printf("Read %i DNL lookup values.\n", nn);

  /* Done! */
  fclose(fp);
  return (0);
}

/*************************************************************************/

int RdGeDINOCalFile(int xtal, const char *fn, GRETINAVariables *gVar) {
  /* Read in a dino calibration file...these are 
     written and read one crystal at a time right now...*/
  
  FILE *fp;
  int i1, i2, nn;
  float f1, f2;
  char *st, str[128];
  
  /* Open file */
  fp = fopen(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\"\n", fn);
    exit(1);
  }
  printf(" %i ", xtal);
  
  /* Read values */
  nn = 0;
  st = fgets(str, 100, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* empty line, do nothing */
    } else {
      sscanf(str, "%i %i %f %f", &i1, &i2, &f1, &f2);
      if (i1>=0) {
	gVar->dinoFactor[xtal-1][i1][i2] = (f1)/100;
	gVar->dinoFactor[xtal-1][i2][i1] = (f2)/100;
	// printf("segment pair %i_%i: %9.5f,  %i_%i: %9.5f\n", 
	//         i1, i2, dinoFactor[xtal-1][i1][i2], 
	//         i2, i1, dinoFactor[xtal-1][i2][i1]);
      }
      nn++;
    }

    /* Attempt to read next line */
    st = fgets(str, 64, fp);
  }
 
  /* Done */
  fclose(fp);
  
  /* Now pack the dino calibration into a matrix to actually use. */
  
/*   TArrayD temp(37, 36); */
/*   for (int i=0; i<35; i++) { */
/*     temp[i] = 1; */
/*     for (int j=0; j<35; j++) { */
/*       temp[(i+(j*36))] = dinoFactor[xtal-1][j][i]; */
/*     } */
/*   } */
  
/*   matrixB[xtal-1].SetMatrixArray(temp.GetArray()); */

/*   dinoMatrix[xtal-1] = matrixB[xtal-1]*Transpose(matrixB[xtal-1]); */
/*   dinoMatrix[xtal-1].Invert(); */
/*   dinoMatrix[xtal-1] *= Transpose(matrixB[xtal-1]); */


  return (0);
}

/*************************************************************************/

int RdGeBaseLineFile(const char *fn) {
  
  /* Declarations */
  FILE *fp;
  int i1, nn;
  float f1, f2, f3, f4, f5, f6, f7, f8;
  char *st, str[256];

  /* Open file */
  fp = fopen(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  }
  printf("\"%s\" open.\n", fn);
  
  /* Read values */
  nn = 0;
  st = fgets(str, 128, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* empty line, do nothing */
    } else {
      sscanf(str, "%i %f %f %f %f %f %f %f %f", 
	     &i1, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8);
      if (i1>=0 && i1<MAXCHANNELS) {
	gWf->baseline[i1] = f3;
	gWf->tau[i1] = f5;
      }
      nn++;
    }

    /* Attempt to read next line */
    st = fgets(str, 128, fp);
  }
  printf("Read %i channel baseline values.\n", nn);

  /* Done! */
  fclose(fp);
  return (0);
}

#endif
