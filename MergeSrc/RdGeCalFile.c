/* $Id: RdGeCalFile.c,v 1.1 2011/04/20 00:00:21 crawford Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>

int
RdGeCalFile(char *fn, float *offset, float *gain)
{
  /* process  cal files */

  /* declarations */

  FILE           *fp;
  int             i1, nn;
  float           f1, f2;
  char           *st, str[128];

  /* open file */

  fp = fopen(fn, "r");
  if (fp == NULL)
  {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  };
  printf("\"%s\" open.\n", fn);

  /* read values */

  nn = 0;
  st = fgets(str, 64, fp);
  while (st != NULL)
  {
    if (str[0] == 35)
    {
      /* '#' comment line, do nothing */

    } else if (str[0] == 59)
    {
      /* ';' comment line, do nothing */

    } else if (str[0] == 10)
    {
      /* empty line, do nothing */

    } else
    {
      sscanf(str, "%i %f %f", &i1, &f1, &f2);
      if (i1>=0 /*&& i1<=320*/)
        {
        offset[i1]=f1;
        gain[i1]=f2;
	// printf("det %3.3i: offset=%9.3f, gain=%9.6f\n", 
	//         i1, offset[i1], gain[i1]);
        fflush(stdout);
        };
      nn++;
    };
    
    /* attempt to read next line */
    
    st = fgets(str, 64, fp);
  };
  printf("Read %i gain calibration coefficients.\n", nn);

  /* done */

  fclose(fp);
  return (0);

};


int
RdDNLGeCalFile(char *fn, float *dnl0, float *dnl1, float *dnl2,
	       float *dnl3, float *dnl4, float *dnl5, float *dnl6,
	       float *dnl7)
{
  /* process  cal files */
  
  /* declarations */

  FILE           *fp;
  int             i1=0, i2=0, i3=0, nn;
  float           f1=0, f2=0, f3=0, f4=0, f5=0, f6=0, f7=0, f8=0;
  char           *st, str[256];

  /* open file */

  fp = fopen(fn, "r");
  if (fp == NULL)
  {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  };
  printf("\"%s\" open.\n", fn);

  /* read values */

  nn = 0;
  st = fgets(str, 128, fp);
  while (st != NULL)
  {
    if (str[0] == 35)
    {
      /* '#' comment line, do nothing */

    } else if (str[0] == 59)
    {
      /* ';' comment line, do nothing */

    } else if (str[0] == 10)
    {
      /* empty line, do nothing */

    } else
    {
      sscanf(str, "%i %i %i %f %f %f %f %f %f %f %f", 
	     &i1, &i2, &i3, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8);
      //cout << i1 << " " << i2 << " " << i3 << " " << f1 << " " << f2 << " " << f3 << " " << f4 << " " << f5 << " " << f6 << " " << f7 << " " << f8 << endl;
      if (i3>=0) {
	dnl0[i3] = f1;
	dnl1[i3] = f2;
	dnl2[i3] = f3;
	dnl3[i3] = f4;
	dnl4[i3] = f5;
	dnl5[i3] = f6;
	dnl6[i3] = f7;
	dnl7[i3] = f8;
      };
      nn++;
    };
    
    /* attempt to read next line */
    
    st = fgets(str, 128, fp);
  };
  printf("Read %i DNL correction calibration coefficients.\n", nn);

  /* done */

  fclose(fp);
  return (0);

};

