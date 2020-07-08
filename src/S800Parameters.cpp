/*     File:  Parameters.cpp
       This file contains the bulk of Parameter stuff, 
       manipulations, etc.
       Author: Heather Crawford
       Date: October 2010                                */

#include "S800Parameters.h"

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>

#include "TString.h"

#define DEBUGS800 0


/************************************************************/
/* S800ImMap Class - Functions                              */
/************************************************************/

void S800ImMap::Initialize(S800Full *s800) {
  /* Parameters */
  m_top = s800;

  /* Variables */
  maxcoefficients = S800_II_TRACK_COEFFICIENTS;
  maxparameters = S800_II_TRACK_PARAMETERS;
  maxorder = 0.0;

  /* Internal use values */
  for (int i=0; i<S800_II_TRACK_PARAMETERS; i++) {
    maxcoefficient[i] = 0;
    for (int j=0; j<S800_II_TRACK_COEFFICIENTS; j++) {
      order[i][j] = 0;
      coefficient[i][j] = 0;
      for (int k=0; k<S800_II_TRACK_PARAMETERS; k++) {
	exponent[i][k][j] = 0;
      }
    }
  }

  // LoadInverseMap(imMapFilename);
}

void S800ImMap::LoadInverseMap(TString filename) {
  cout << "Loading S800 Im inverse map... ";
  FILE *file;
  file = fopen(filename.Data(), "r");
  if (file == NULL) {
    printf("\n Im inverse map file %s was not found!\n", filename.Data());
    return;
  }

  char line[400];

  unsigned short index;
  unsigned short icoefficient;
  unsigned short iorder;
  unsigned short iexp0, iexp1, iexp2, iexp3, iexp4, iexp5;
  int parameter = 0;

  char* retVal;

  while (!feof(file)) {
    retVal = fgets(line, 400, file);
    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, "COEFFICIENT", 11) == 0) {
      parameter = 0;
      while (!feof(file)) {
	retVal = fgets(line, 400, file);
	while (strncmp(line, "#", 1) == 0) {
	  sscanf(line, "%hu %hu %hu %hu %hu %hu %hu %hu %hu", &index, 
		 &icoefficient, &iorder, &iexp0, &iexp1, &iexp2, &iexp3, 
		 &iexp4, &iexp5);
	  if (iorder > maxorder) { maxorder = (double)iorder; }
	  if (index > maxcoefficients) {
	    cout << "Error: too many coefficients in map" << endl;
	    cout << " Please increase S800_TRACK_COEFFICIENTS and "<< endl
		 << " recompile. " << endl;
	    break;
	  }
	  
	  if (parameter >= S800_TRACK_PARAMETERS) {
	    cout << "Invalid parameter number, must be between 0 and " 
		 << S800_TRACK_PARAMETERS-1 << endl;
	  }
	  if (index >= S800_TRACK_COEFFICIENTS) {
	    cout << "Too many coefficients - maximum is " 
		 << S800_TRACK_COEFFICIENTS << endl;
	  }
	  if (index >= maxcoefficient[parameter]) { 
	    maxcoefficient[parameter] = index+1;
	  }
	  if (iorder > maxorder) { maxorder = iorder; }
	  order[parameter][index] = iorder;
	  exponent[parameter][0][index] = iexp0;
	  exponent[parameter][1][index] = iexp1;
	  exponent[parameter][2][index] = iexp2;
	  exponent[parameter][3][index] = iexp3;
	  exponent[parameter][4][index] = iexp4;
	  exponent[parameter][5][index] = iexp5;
	  coefficient[parameter][index] = icoefficient;

	  retVal = fgets(line, 400, file);
	}
	retVal = fgets(line, 400, file);
	parameter++;
	if (parameter > maxparameters) {
	  cout << "Error: too many parameters in map" << endl;
	  cout << " Please increase S800_TRACK_PARAMETERS and " << endl
	       << " recompile. " << endl;
	  break;
	}
      }
    }
  }
  fclose(file);  
  
  cout << "Done!" << endl;  
}

void S800ImMap::Reset() {
  /* Internal use values */
  for (int i=0; i<S800_II_TRACK_PARAMETERS; i++) {
    maxcoefficient[i] = 0;
    for (int j=0; j<S800_II_TRACK_COEFFICIENTS; j++) {
      order[i][j] = 0;
      coefficient[i][j] = 0;
      for (int k=0; k<S800_II_TRACK_PARAMETERS; k++) {
	exponent[i][k][j] = 0;
      }
    }
  }
}

double S800ImMap::Calculate(int calcorder, int parameter, double *input) {
  double cumul = 0;
  double multiplicator;

  for (int index=0; index<maxcoefficient[parameter]; index++) {
    if (calcorder < order[parameter][index]) { break; }
    
    multiplicator = 1;
    
    for (int nex=0; nex<S800_II_TRACK_PARAMETERS; nex++) {
      if (exponent[parameter][nex][index] != 0) {
	multiplicator *= pow(input[nex], exponent[parameter][nex][index]);
      }
    }

    cumul += multiplicator * coefficient[parameter][index];
  }

  return cumul;
}

/************************************************************/
/* S800Map Class - Functions                                */
/************************************************************/

void S800Map::Initialize(S800Full *s800) {
  /* Parameters */
  m_top = s800;

  /* Variables */
  maxcoefficients = S800_TRACK_COEFFICIENTS;
  maxparameters = S800_TRACK_PARAMETERS;
  maxorder = 0.0;

  /* Internal use values */
  for (int i=0; i<S800_TRACK_PARAMETERS; i++) {
    maxcoefficient[i] = 0;
    for (int j=0; j<S800_TRACK_COEFFICIENTS; j++) {
      order[i][j] = 0;
      coefficient[i][j] = 0;
      for (int k=0; k<S800_TRACK_PARAMETERS; k++) {
	exponent[i][k][j] = 0;
      }
    }
  }

  // LoadInverseMap(filename);

}

void S800Map::LoadInverseMap(TString filename) {
  cout << "Loading S800 inverse map " << filename.Data() << "...";

  FILE *file;
  if ( (file = fopen(filename.Data(), "r")) == NULL ) {
    printf("\n Inverse map file %s was not found!\n", filename.Data());
    return;
  }

  char line[400];

  unsigned short index;
  double icoefficient;
  unsigned short iorder;
  unsigned short iexp0, iexp1, iexp2, iexp3, iexp4, iexp5;
  int parameter = 0;
  char* retVal;

  retVal = fgets(line, 400, file);
  while (strncmp(line, "     I", 6) != 0) { 
    retVal = fgets(line, 400, file);
  }
  parameter = 0;
  
  while (!feof(file)) {
    retVal = fgets(line, 400, file);
    while (strncmp(line, "    ---", 7) != 0) {
      sscanf(line, "%hd %lf %hd %hd %hd %hd %hd %hd %hd", &index, 
	     &icoefficient, &iorder, &iexp0, &iexp1, &iexp2, &iexp3, 
	     &iexp4, &iexp5);
      if (iorder > maxorder) { maxorder = (double)iorder; }
      if (index > maxcoefficients) {
	cout << "Error: too many coefficients in map" << endl;
	cout << " Please increase S800_TRACK_COEFFICIENTS and "<< endl
	     << " recompile. " << endl;
	break;
      }

      if (parameter >= S800_TRACK_PARAMETERS) {
	cout << "Invalid parameter number, must be between 0 and " 
	     << S800_TRACK_PARAMETERS-1 << endl;
      }
      if (index >= S800_TRACK_COEFFICIENTS) {
	cout << "Too many coefficients - maximum is " 
	     << S800_TRACK_COEFFICIENTS << endl;
      }
      if (index >= maxcoefficient[parameter]) { 
	maxcoefficient[parameter] = index+1;
      }
      order[parameter][index] = iorder;
      exponent[parameter][0][index] = iexp0;
      exponent[parameter][1][index] = iexp1;
      exponent[parameter][2][index] = iexp2;
      exponent[parameter][3][index] = iexp3;
      exponent[parameter][4][index] = iexp4;
      exponent[parameter][5][index] = iexp5;
      coefficient[parameter][index] = icoefficient;
      retVal = fgets(line, 400, file);
    }
    retVal = fgets(line, 400, file);
    parameter++;
    if (parameter > maxparameters) {
      cout << "Error: too many parameters in map" << endl;
      cout << " Please increase S800_TRACK_PARAMETERS and " << endl
	   << " recompile. " << endl;
      break;
    }
  }
  fclose(file);  
 
  cout << "Done!" << endl;  
}

void S800Map::Reset() {
  /* Internal use values */
  for (int i=0; i<S800_TRACK_PARAMETERS; i++) {
    //maxcoefficient[i] = 0;
    for (int j=0; j<S800_TRACK_COEFFICIENTS; j++) {
      //order[i][j] = 0;
      //coefficient[i][j] = 0;
      for (int k=0; k<S800_TRACK_PARAMETERS; k++) {
	//exponent[i][k][j] = 0;
      }
    }
  }
}

double S800Map::Calculate(int calcorder, int parameter, double *input) {
  double cumul = 0;
  double multiplicator;

  for (int index=0; index<maxcoefficient[parameter]; index++) {
    if (calcorder < order[parameter][index]) { break; }
    multiplicator = 1;

    for (int nex=0; nex<S800_TRACK_PARAMETERS; nex++) {
      if (exponent[parameter][nex][index] != 0) {
	multiplicator *= pow(input[nex], exponent[parameter][nex][index]);
      }
    }

    cumul += multiplicator * coefficient[parameter][index];
  }

  return cumul;
}

/************************************************************/
/* S800FpCrdcPads Class - Functions                         */
/************************************************************/

UShort_t* S800FpCrdcPads::UnpackRawDataFastClassic(UShort_t *p) {
  /* Classic unpacking without : m_thresh checking & m_samples parameter */
  /* Purpose:
        Simple way for event parsing a la Daniel's original 
	(otherwise one has to set all m_thresh[] to 0, and 
	introducing an additional variable in UnpackRawDataFast
	would be time consuming [test in inner loop])  */
  /* For S800_FP_CRDC_RAW_PACKET only!  (not checked)
     Input data structure -> see class description in Parameter.h
     REMARK: all if-conditions seem to be needed, input data is 'messed' */
   UShort_t SubLen;
   UShort_t maxsample;
   UShort_t data, sample, begin, previous, chn, padChn, value;
   int thresh;

  m_minChn = S800_FP_CRDC_CHANNELS;
  m_maxChn = 0;
  m_width = 0;
  memset(m_dprojX, 0, sizeof(double)*S800_FP_CRDC_CHANNELS);

  m_padlist.clear();

  SubLen = *p; p += 2;
  SubLen -= 2;

  thresh = *p; ++p;
  SubLen--;

  /* First data defines 'begin' sample in time.
     REMARK: pedestal subtraction done in inner loop, 
             only accept values > pedestal   */
  if (SubLen) {
    data = *p; p++;
    SubLen--;
    begin = (data & 0x7FC0) >> 6;
    previous = begin;
    maxsample = begin + S800_FP_CRDC_MAX_WIDTH;
  }

  while (SubLen > 0) {
    sample = (data & 0x7FC0) >> 6;
    chn = (data & 0x003F);

    if ( (sample < previous) ||
	 (sample > maxsample) ) {
      while (SubLen) { /* Ignore this data */
	data = *p; p++; 
	SubLen--;
	if (data & 0x8000) { break; }
      } 
    } else {
      m_width = sample - begin;
      while (SubLen) { /* Process this data */
	data = *p; p++;
	SubLen--;
	if (data & 0x8000) {
	  previous = sample;
	  break;
	}

	padChn = ((data & 0x0C00) >> 10)*64 + chn;
	value = data & 0x03FF;
	
	if (padChn < S800_FP_CRDC_CHANNELS &&
	    value > thresh) {
	  if (m_dprojX[padChn] == 0) {
	    m_padlist.push_back(padChn);
	  }
	  m_dprojX[padChn] += value - ped[padChn];

	  if (padChn > m_maxChn) { m_maxChn = padChn; }
	  if (padChn < m_minChn) { m_minChn = padChn; }
	}
      }
    }
  }
  m_width += 1;
  m_maxChn += 1;

  return p;
}

UShort_t* S800FpCrdcPads::UnpackRawDataFast(UShort_t *p) {
  /* For S800_FP_CRDC_RAW_PACKET only (not checked ! )
     Input data structure-> see class description Parameters.h
     REMARK: All if conditions needed, if input is somehow 'messed' */
  UShort_t SubLen;
  UShort_t maxsample;
  UShort_t data, sample, begin, previous, chn, padChn, value;
  int thresh;

  m_minChn = S800_FP_CRDC_CHANNELS;
  m_maxChn = 0;
  m_width = 0;
  memset(m_dprojX, 0, sizeof(double)*S800_FP_CRDC_CHANNELS);
  memset(m_samples, 0, sizeof(int)*S800_FP_CRDC_CHANNELS);

  memset(m_samplemax, 0, sizeof(int)*S800_FP_CRDC_CHANNELS);
  memset(m_samplemin, 0x01, sizeof(int)*S800_FP_CRDC_CHANNELS);

  m_padlist.clear();

  SubLen = *p; p += 2;
  SubLen -= 2;
  
  thresh = *p; ++p;
  SubLen--;
  
  /* First data defines 'begin' sample in time. */
  /* REMARK: pedestal subtraction done in inner loop, 
     and we only accept values > pedestal */
  if (SubLen) {
    data = *p; p++; 
    SubLen--;
    begin = (data & 0x7FC0) >> 6;
    previous = begin;
    maxsample = begin+S800_FP_CRDC_MAX_WIDTH;
  }

  while (SubLen > 0) {
    sample = (data & 0x7FC0) >> 6;
    chn = (data & 0x003F);
    
    if ( (sample < previous) || (sample > maxsample) ) {
      while (SubLen) { /* Ignore this data */
	data = *p; p++;
	SubLen--;
	if (data & 0x8000) { break; }
      } 
    } else {
      m_width = sample - begin;
      while (SubLen) { /* Process this data */
	data = *p; p++;
	SubLen--;
	
	if (data & 0x8000) {
	  previous = sample;
	  break;
	}

	padChn = ((data & 0x0C00) >> 10)*64 + chn;
	value = data & 0x03FF;

	if ( padChn < S800_FP_CRDC_CHANNELS &&
	     value > thresh && 
	     value > m_thresh[padChn] + ped[padChn] ) {
	  if (m_dprojX[padChn] == 0) {
	    m_padlist.push_back(padChn);
	  }

	  m_dprojX[padChn] += value - ped[padChn];
	  m_samples[padChn]++;
	  
	  if (value > m_samplemax[padChn]) { m_samplemax[padChn] = value; }
	  if (value < m_samplemin[padChn]) { m_samplemin[padChn] = value; }
	  
	  if (padChn > m_maxChn) { m_maxChn = padChn; }
	  if (padChn < m_minChn) { m_minChn = padChn; }
	}
      }
    }
  }
  m_width += 1;
  m_maxChn += 1;

  return p;
}

void S800FpCrdcPads::CalcDelta(int center, S800FpCrdcCalc* calc) {
  /* Idea: center is result for x-position calculation in crdc, and
     we want to know the distance of all hit pads to this center. 
     If center is <= 0, then delta calculation is omitted.
     Deleting of data for bad pads will be done anyway. */
  for (m_padlist_iter = m_padlist.begin();
       m_padlist_iter != m_padlist.end();
       m_padlist_iter++) {
    if (calc->isbadpad(*m_padlist_iter)) {
      cal[*m_padlist_iter] = 0;
      samples[*m_padlist_iter] = 0;
      samplemax[*m_padlist_iter] = 0;
      samplemin[*m_padlist_iter] = 0;
    } else if (center > 0) {
      delta[*m_padlist_iter] = *m_padlist_iter - center;
    }
  }
}

void S800FpCrdcPads::BuildLookUp() {
  /* Build look-up table */
  for (int j = 0; j < 0x0400; j++) {
    for (int k = 0; k < S800_FP_CRDC_CHANNELS; k++) {
      m_caltable[k][j] = (int)(slope[k]*(double)j + offset[k]);
    }
  }
}


void S800FpCrdcPads::Initialize(S800Full *s800) {
  m_top = s800;
  
  for (int i=0; i<S800_FP_CRDC_CHANNELS; i++) {
    /* Parameters */
    cal[i] = -65000.;
    raw[i] = -65000.;
    samples[i] = -65000;
    samplemax[i] = -65000.;
    samplemin[i] = -65000.;
    delta[i] = -65000.;

    /* Variables */
    ped[i] = 0;
    minsamples[i] = 1;
    m_thresh[i] = 0;
    slope[i] = 1;
    offset[i] = 0;
  }

}

void S800FpCrdcPads::Reset() {
  for (int i=0; i<S800_FP_CRDC_CHANNELS; i++) {
    /* Parameters */
    cal[i] = -65000.;
    raw[i] = -65000.;
    samples[i] = -65000;
    samplemax[i] = -65000.;
    samplemin[i] = -65000.;
    delta[i] = -65000.;
  }
}

UShort_t* S800FpCrdcPads::UnpackRawData(UShort_t *p) {
  p = UnpackRawDataFast(p);

  /* This does the calibration. */
  for (m_padlist_iter = m_padlist.begin();
       m_padlist_iter != m_padlist.end();
       m_padlist_iter++) {
    double value;
    value = m_dprojX[*m_padlist_iter]/m_width;
    raw[*m_padlist_iter] = value; /* Bad pads, 'minsample' pads get raw, 
				     but no calibration */
    if (m_samplemax[*m_padlist_iter] > 0) {
      samplemax[*m_padlist_iter] = m_samplemax[*m_padlist_iter];
      samplemin[*m_padlist_iter] = m_samplemin[*m_padlist_iter];
    }

    if (m_samples[*m_padlist_iter]+1 > minsamples[*m_padlist_iter]) {
      cal[*m_padlist_iter] = 
	( (slope[*m_padlist_iter]*
	   (double)m_dprojX[*m_padlist_iter]/(double)m_width) +
	  (offset[*m_padlist_iter]) );
      samples[*m_padlist_iter] = m_samples[*m_padlist_iter];
    } else {
      /* Remove this pad from the list, so it can't be maxpad...
	 set it invalid or zero.  It won't be extrapolated as
	 a bad pad in x_gravity, because assumption is that 
	 pad < minsamples has not seen something. */
      cal[*m_padlist_iter] = -65000.;
      m_padlist_iter = m_padlist.erase(m_padlist_iter);
    }
  }
  return p;
}

UShort_t* S800FpCrdcPads::UnpackRawDataClassic(UShort_t *p) {
  p = UnpackRawDataFastClassic(p);
  memset(cal, 0, sizeof(double)*S800_FP_CRDC_CHANNELS);

  for (m_padlist_iter = m_padlist.begin();
       m_padlist_iter != m_padlist.end();
       m_padlist_iter++) {
    raw[*m_padlist_iter] = (double)m_dprojX[*m_padlist_iter]/(double)m_width;
    cal[*m_padlist_iter] = 
      ( (slope[*m_padlist_iter]*
	 (double)m_dprojX[*m_padlist_iter]/(double)m_width) +
	(offset[*m_padlist_iter]) );
  }

  return p;
}

void S800FpCrdcPads::AverageCal(int integration) {
  /***** NOT DEFINED ANYWHERE, OR USED *****/
}

/************************************************************/
/* S800FpCrdcCalc Class - Functions                         */
/************************************************************/

void S800FpCrdcCalc::FindPadMax() {
  /* Find pad at avalanche site maximum (2 non-zero adjacent pads
     required).  Same algorithm as in classic SeGA SpecTcl.  Bad
     pad status NOT checked for neightbours. 

     If we fail:  padmax becomes INVALID, padsum is still valid

     REMARK: pad 0 and 224 (most left/right) cannot be padmax */
  
  double my_maxpad = 0;
 
  int maxChn = m_pad->m_maxChn - 1; /* m_maxChn is NOT inclusive */
  padsum = 0;
  
  if (m_pad->m_padlist.empty()) { return; }

  /* The <m_padlist> is not sorted!  Includes bad pads */
  for (m_pad->m_padlist_iter = m_pad->m_padlist.begin(); 
       m_pad->m_padlist_iter != m_pad->m_padlist.end();
       m_pad->m_padlist_iter++) {
    if (!m_badpad[*m_pad->m_padlist_iter]) {
      padsum += m_pad->cal[*m_pad->m_padlist_iter];
      if (m_pad->cal[*m_pad->m_padlist_iter] > my_maxpad) {
	if (*m_pad->m_padlist_iter != m_pad->m_minChn && 
	    m_pad->cal[*m_pad->m_padlist_iter-1] != 1000 &&
	    m_pad->cal[*m_pad->m_padlist_iter-1] > 0 &&
	    *m_pad->m_padlist_iter != maxChn &&
	    m_pad->cal[*m_pad->m_padlist_iter+1] != -65000. && 
	    m_pad->cal[*m_pad->m_padlist_iter+1] > 0) {
	  my_maxpad = m_pad->cal[*m_pad->m_padlist_iter];
	  padmax = *m_pad->m_padlist_iter;
	}
      }
    }
  }

  /* Classic SpecTcl does padsum only if padmax is found */
  if (padmax == -65000.) {
    padsum = -65000.;
  }

  if (padmax != -65000.) {
    maxpad[(int)padmax] = my_maxpad;
    maxpad0[(int)padmax] = m_pad->cal[(int)padmax];
    if (!m_badpad[(int)padmax-1] && m_pad->cal[(int)padmax-1] != -65000.) {
      maxpad0[(int)padmax-1] = m_pad->cal[(int)padmax-1];
    }
    if (!m_badpad[(int)padmax+1] && m_pad->cal[(int)padmax+1] != -65000.) {
      maxpad0[(int)padmax+1] = m_pad->cal[(int)padmax+1];
    }
  }
}


void S800FpCrdcCalc::CalcSkew(int left, int right, double center) {
  /***** NOT DEFINED ANYWHERE, OR USED *****/
}

void S800FpCrdcCalc::CheckBadPads() {

  cout << "Checking CRDC bad pads... ";

  badpads = 0;
  for (int j=0; j<99; j++) {
    if (badpad[j] > 0) { badpads++; }
  }

  for (int j=0; j<S800_FP_CRDC_CHANNELS; j++) {
    m_badpad[j] = false;
  }
  
  if (badpads != 0) {
    cout << endl << " Bad pads (" << badpads << ") :";
  }
  for (int j=0; j<badpads; j++) {
    if ( ((int)badpad[j] >= 0) && ((int)badpad[j] < S800_FP_CRDC_CHANNELS) ) {
      m_badpad[(int)badpad[j]] = true;
      cout << " " << badpad[j];
    }
  }
  if (badpads != 0) {
    cout << endl;     
  }

  cout << "Done!" << endl;
}

void S800FpCrdcCalc::Initialize(S800Full *s800, S800FpCrdcPads *pad, S800FpCrdc *fpcrdc) {
  m_top = s800;

  m_pad = pad;
  m_parent = fpcrdc;

  /* Parameter values */
  padmax = -65000.;
  padsum = -65000.;
  padused = -65000.;
  
  for (int i=0; i<S800_FP_CRDC_CHANNELS; i++) {
    maxpad[i] = -65000.;
    maxpad0[i] = -65000.;
    padcalc[i] = -65000.;
  }
  
  saturationWidth = -65000.;
  x_gravity = -65000.;
  
  /* Variables */
  gravity_width = 15;
  fit_width = 9;
  method = 1;
  saturation = 1024;
  m_saturation = 1.00;

  badpads = 0;
  for (int i=0; i<99; i++) {
    badpad[i] = 0;
  }

}

void S800FpCrdcCalc::Reset() {
  /* Parameter values */
  padmax = -65000.;
  padsum = -65000.;
  padused = -65000.;
  
  for (int i=0; i<S800_FP_CRDC_CHANNELS; i++) {
    maxpad[i] = -65000.;
    maxpad0[i] = -65000.;
    padcalc[i] = -65000.;
  }
  
  saturationWidth = -65000.;
  x_gravity = -65000.;
}

void S800FpCrdcCalc::CalculateGravity() {
  /* New approach: 
      1) Check for saturation
      2) Correct padmax as CoM of all (neighbouring) saturated pads.  
      Here the saturation can be scaled with m_saturation.
      (e.g. all pads xx% saturated and more are taken into account)
      3) For CoM with padmax as center extrapolate bad pads */
  
  FindPadMax();
  if (padmax == -65000.) {
    m_pad->CalcDelta(-1, this);
    return;
  }

  /* Rerun padmax: figure out CoM of the saturated pads
     as the new padmax
     Remember: m_pad->cal[chn] only valid in [minChn, maxChn] */
  int l = (int)padmax, r = (int)padmax;

  while(1) {
    while( (l > m_pad->m_minChn) && 
	   (m_pad->cal[l] != -65000.) &&
	   (m_pad->cal[l] >= saturation*m_saturation) ) {
      l--;
    }
    /* 2nd chance for bad pads */
    if ( (m_badpad[l]) && 
	 (l > m_pad->m_minChn) &&
	 (m_pad->cal[l-1] != -65000.) &&
	 (m_pad->cal[l-1] >= saturation*m_saturation) ) {
      l--;
    } else {
      break; 
    }
  }

  while(1) {
    while( (r+1 < m_pad->m_maxChn) &&
	   (m_pad->cal[r] != -65000.) && 
	   (m_pad->cal[r] >= saturation*m_saturation) ) {
      r++;
    }
    /* 2nd chance for bad pads */
    if ( (m_badpad[r]) && 
	 (r+2 < m_pad->m_maxChn) &&
	 (m_pad->cal[r+1] != -65000.) &&
	 (m_pad->cal[r+1] >= saturation*m_saturation) ) {
      r++;
    } else {
      break;
    }
  }

  padmax = (l+r)/2.;
  saturationWidth = r-l;

  /* Calculate region for CoM calculation */
  int lowpad = (int)padmax - (int)gravity_width/2;
  int highpad = lowpad + (int)gravity_width;
  
  if (lowpad < m_pad->m_minChn) { lowpad = m_pad->m_minChn; }
  if (highpad >= m_pad->m_maxChn) { highpad = m_pad->m_maxChn; }
  
  /* Calculates pad center of gravity.
     - extrapolate -> bad pads will be extrapolated from their 
     not bad neighbours
     - saturation -> weight of saturated pads is limited to saturation */
  int sum = 0, mom = 0;
  padused = 0;
  
  for (int i=lowpad; i<highpad; i++) {
    double padamp;
    
    if (!m_badpad[i]) {
      if (m_pad->cal[i] != -65000.) {
	padamp = m_pad->cal[i];
      } else {
	padamp = 0;
      }
    } else {
      padamp = 0;
      if ( (i > m_pad->m_minChn) && (!m_badpad[i-1]) ) {
	if (m_pad->cal[i-1] != -65000.) {
	  padamp += (int)(0.5*m_pad->cal[i-1]);
	}
      }
      if ( (i < m_pad->m_maxChn-1) && (!m_badpad[i+1]) ) {
	if (m_pad->cal[i+1] != -65000.) {
	  padamp += (int)(0.5*m_pad->cal[i+1]);
	}
      }
    }
    
    if (padamp > saturation) {
      padamp = saturation;
    }

    sum += (int) padamp;
    mom += i*(int)padamp;
    
    padcalc[i] = padamp;
    if (padamp > 1) {
      padused += 1;
    }
  } /* End of for (int i=lowpad; i<highpad; i++) */
  
  x_gravity = (double)mom/sum;
  
  m_pad->CalcDelta((int)x_gravity, this);
}

void S800FpCrdcCalc::CalculateGravityClassic() {
  /* Classic gravity without bad pad interpolation, etc. */
  FindPadMax();

  if (padmax == -65000.) {
    return;
  }

  /* Calculate region for CoM calculation */
  int lowpad = (int)padmax - (int)gravity_width/2;
  int highpad = lowpad + (int)gravity_width;
  
  if (lowpad < m_pad->m_minChn) { lowpad = m_pad->m_minChn; }
  if (highpad >= m_pad->m_maxChn) { highpad = m_pad->m_maxChn; }
  
  /* Calculated pad center of gravity */
  int sum = 0, mom = 0;
  
  for (int i=lowpad; i<highpad; i++) {
    if ( (!m_badpad[i]) && (m_pad->cal[i] != -65000.) ) {
      sum += (int) m_pad->cal[i];
      mom += i*(int) m_pad->cal[i];
    }
  }

  x_gravity = (double)mom/sum;
}

/************************************************************/
/* S800FpCrdc Class - Functions                             */
/************************************************************/

void S800FpCrdc::Initialize(S800Full *s800) {
  m_top = s800;

  S800FpCrdc *m_crdc = this;
  
  pad.Initialize(m_top);
  calc.Initialize(m_top, &pad, m_crdc);

  /* Parameters */
  x = -65000.;
  y = -65000.;
  tac = -65000.;
  anode = -65000.;
  evtlength = -65000.;

  /* Variables */
  x_offset = -300; /* mm */
  x_slope = 2.54;  /* mm/pad */
  y_offset = 130;  /* mm */
  y_slope = -0.09; /* mm/ns */

}

void S800FpCrdc::Reset() {
  pad.Reset();
  calc.Reset();

  /* Parameters */
  x = -65000.;
  y = -65000.;
  tac = -65000.;
  anode = -65000.;
  evtlength = -65000.;
}

void S800FpCrdc::CalculatePosition() {
  /***** NOT DEFINED ANYWHERE, OR USED *****/  
}

UShort_t* S800FpCrdc::Unpack(UShort_t *p) {
  evtlength = *p;
  
  calc.method = 1;

  p += 3;
  if (calc.method == 1) { /* Classic treatment */
    p = pad.UnpackRawDataClassic(p);
    calc.CalculateGravityClassic();
  } else if (calc.method == 2) {
    p = pad.UnpackRawData(p);
    calc.CalculateGravity();
  }

  p += 2; /* SubLength (=4), SubTag (=S800_FP_CRDC_ANODE_PACKET) */
  anode = *p; p++;
  tac = *p; p++;
  
  if (calc.x_gravity != -65000.) {
    x = (calc.x_gravity * x_slope) + x_offset;
  }
  
  y = (tac * y_slope) + y_offset;
  y += (x * xy_slope) + xy_offset;
  
  return p;
}

/************************************************************/
/* S800Tppac Class - Functions                              */
/************************************************************/

void S800Tppac::CheckBadStrips() {
  badxstrips = 0; badystrips = 0;
  for (int j=0; j<20; j++) {
    if (badxstrip[j] >= 0) { badystrips++; }
    if (badystrip[j] >= 0) { badxstrips++; }
  }  

  for (int i=0; i<64; i++) {
    m_badxstrip[i] = false;
    m_badystrip[i] = false;
  }

  for (int i=0; i<badxstrips; i++) {
    m_badxstrip[(int)badxstrip[i]] = true;
  }
  for (int i=0; i<badystrips; i++) {
    m_badystrip[(int)badystrip[i]] = true;
  }
}

void S800Tppac::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  for (int i=0; i<64; i++) {
    xstrip[i] = -65000.;
    ystrip[i] = -65000.;
  }
  xsum = -65000.;
  ysum = -65000.;
  xmult = -65000.;
  ymult = -65000.;
  x = -65000.;
  y = -65000.;

  /* Variables */
  badxstrips = 0; /* number */
  badystrips = 0; /* number */
  width = 32; /* strips */
  xcenter = 40.64; /* mm */
  ycenter = 40.64; /* mm */
  for (int i=0; i<20; i++) {
    badxstrip[i] = 0; /* pad */
    badystrip[i] = 0; /* pad */
  }
}

void S800Tppac::Reset() {
  /* Parameters */
  for (int i=0; i<64; i++) {
    xstrip[i] = -65000.;
    ystrip[i] = -65000.;
  }
  xsum = -65000.;
  ysum = -65000.;
  xmult = -65000.;
  ymult = -65000.;
  x = -65000.;
  y = -65000.;
}

void S800Tppac::Calculate(int xleft, int xright, int yleft, int yright) {
  int maxpad = 0, padmax = 0;
  
  /* x-plane: xright is not > xleft if no xstrip has been read */
  if (xright > xleft) {
    for (int i=xleft; i<xright; i++) {
      if ( (!m_badxstrip[i]) && (xstrip[i] != -65000.) && (xstrip[i] > maxpad) ) {
	maxpad = (int) xstrip[i];
	padmax = i;
      }
    }
    
    int lowpad = padmax - (int)width/2;
    int highpad = lowpad + (int)width;
    
    if (lowpad < xleft) { lowpad = xleft; }
    if (highpad > xright) { highpad = xright; }
    if (highpad == 64) { highpad--; }

    int sum = 0, mom = 0;
    /* Reject bad pads stored in the badpad array */
    for (int i=lowpad; i<highpad; i++) {
      if ( (!m_badxstrip[i]) && (xstrip[i] != -65000.) ) {
	sum += (int)xstrip[i];
	mom += i*(int)xstrip[i];
      }
    }

    xsum = (double)sum;
    x = ((double)mom/xsum)*1.27 - xcenter; /* Pitch of strips is 1/20" */
  } /* End of if (xright > xleft) */

  /* y-plane */
  maxpad = 0;
  padmax = 0;
  
  if (yright > yleft) {
    for (int i=yleft; i<yright; i++) {
      if ( (!m_badystrip[i]) && (ystrip[i] != -65000.) && (ystrip[i] > maxpad) ) {
	maxpad = (int)ystrip[i];
	padmax = i;
      }
    }

    int lowpad = padmax - (int)width/2;
    int highpad = lowpad + (int)width;
    
    if (lowpad < yleft) { lowpad = yleft; }
    if (highpad > yright) { highpad = yright; }
    if (highpad == 64) { highpad--; }

    int sum = 0, mom = 0;
    /* Reject bad pads stored in the badpad array */
    for (int i=lowpad; i<highpad; i++) {
      if ( (!m_badystrip[i]) && (ystrip[i] != -65000.) ) {
	int yst = ystrip[i];
	sum += yst;
	mom += i*yst;
      }
    }

    ysum = (double)sum;
    y = ((double)mom/ysum)*1.27 - ycenter; /* Pitch of strips is 1/20" */
  } /* End of if (yright > yleft) */
}

/************************************************************/
/* S800Track Class - Functions                              */
/************************************************************/

void S800Track::Initialize(S800Full *s800) {
  m_top = s800;
 
  tppac1.Initialize(m_top);
  tppac2.Initialize(m_top);

  /* Parameters */
  xim = -65000.;
  yim = -65000.;
  aim = -65000.;
  bim = -65000.;
  dim = -65000.;
  xta = -65000.;
  yta = -65000.;
  ata = -65000.;
  bta = -65000.;

  order = 5; /* order */

  map.Initialize(m_top);
}

void S800Track::Reset() {
  tppac1.Reset();
  tppac2.Reset();

  /* Parameters */
  xim = -65000.;
  yim = -65000.;
  aim = -65000.;
  bim = -65000.;
  dim = -65000.;
  xta = -65000.;
  yta = -65000.;
  ata = -65000.;
  bta = -65000.;
  
  map.Reset();
}

UShort_t* S800Track::UnpackRawData(UShort_t *p) {
  int SubLen;
  int data = 0, sample, begin = 0, previous = 0, chn, det, value;
  
#if (S800_II_TRACK_CHANNELS > 4*64)
#error This number of S800_II_TRACK_CHANNELS cannot be set, sorry
#endif

  for (int i=0; i<4; i++) {
    memset((char*)m_projX[i], 0, sizeof(int)*64);
    m_minChn[i] = 64;
    m_maxChn[i] = 0;
  }

  m_width = 0;
  
  p += 2; /* PktLen and PktTag S800_II_TRACK_PACKET */
  
  SubLen = *p; p += 2; /* SubTag is S800_II_TRACK_RAW_PACKET */
  SubLen -= 2;
  
  m_thresh = *p; ++p;
  SubLen--;

  if (SubLen) {
    data = *p; p++; 
    SubLen--;
    begin = (data & 0x7FC0) >> 6;
    previous = begin;
  }

  while (SubLen > 0) {
    sample = (data & 0x7FC0) >> 6;
    chn = (data & 0x003F);
    
    if ( (sample<begin) || (sample < previous) || 
	 (sample > begin + S800_II_TRACK_MAX_WIDTH) ) {
      while (SubLen) { /* Ignore this data */
	data = *p; p++;
	SubLen--;
	if (data & 0x8000) { break; };
      } 
    } else {
      while(SubLen) { /* Process this data */
	data = *p; p++;
	SubLen--;
	if (data & 0x8000) {
	  previous = sample;
	  m_width = sample - begin;
	  break;
	}
	
	det = (data & 0x0C00) >> 10;
	value = data & 0x03FF;
	
	if (value > m_thresh) {
	  chn = (det%2 ? tppacmapy[chn] : tppacmapx[chn]);
	  m_projX[det][chn] += value;
	  if (chn > m_maxChn[det]) { m_maxChn[det] = chn; }
	  if (chn < m_minChn[det]) { m_minChn[det] = chn; }
	}
      }
    }
  }

  m_width += 1;
  m_maxChn[0] += 1;
  m_maxChn[1] += 1;
  m_maxChn[2] += 1;
  m_maxChn[3] += 1;

  int mult = 0;
  for (int chn = m_minChn[0]; chn < m_maxChn[0]; chn++) {
    tppac1.xstrip[chn] = m_projX[0][chn];
    if (tppac1.xstrip[chn] > 0) { mult++; }
  }
  tppac1.xmult = mult;

  mult = 0;
  for (int chn = m_minChn[1]; chn < m_maxChn[1]; chn++) {
    tppac1.ystrip[chn] = m_projX[1][chn];
    if (tppac1.ystrip[chn] > 0) { mult++; }
  }
  tppac1.ymult = mult;
    
  /* Calculate position */
  tppac1.Calculate(m_minChn[0], m_maxChn[0], m_minChn[1], m_maxChn[1]);
  
  mult = 0;
  for (int chn = m_minChn[2]; chn < m_maxChn[2]; chn++) {
    tppac2.xstrip[chn] = m_projX[2][chn];
    if (tppac2.xstrip[chn] > 0) { mult++; }
  }
  tppac2.xmult = mult;
	
  mult = 0;
  for (int chn = m_minChn[3]; chn < m_maxChn[3]; chn++) {
    tppac2.ystrip[chn] = m_projX[3][chn];
    if (tppac2.ystrip[chn] > 0) { mult++; }
  }
  tppac2.ymult = mult;
  
  /* Calculate position */
  tppac2.Calculate(m_minChn[2], m_maxChn[2], m_minChn[3], m_maxChn[3]);

  return p;
}

void S800Track::CalculateTracking() {
  /* Note: xim, aim, ... have to be valid in every case!!
     This method is called by S800Image::Calculate(), where 
     xim, aim, etc. are set. */
  double input[S800_II_TRACK_PARAMETERS];
  if (!map.WasLoaded()) { return; }
  input[0] = xim;
  input[1] = aim;
  input[2] = yim;
  input[3] = bim;
  input[5] = dim;
  if ((double)order > (double)map.maxorder) { order = map.maxorder; }
  xta = map.Calculate((int)order, 0, input);
  ata = map.Calculate((int)order, 1, input);
  yta = map.Calculate((int)order, 2, input);
  bta = map.Calculate((int)order, 3, input);
}

/************************************************************/
/* S800FpIonChamber Class - Functions                       */
/************************************************************/

void S800FpIonChamber::BuildLookUp() {
  for (int i=0; i<S800_FP_IC_CHANNELS; i++) {
    for (int chn = 0; chn < 0x0fff; chn++) {
      caltable_de[chn] = (UShort_t) (de_slope*chn + de_offset);
      caltable_xtilt[chn] = exp(de_xtilt * (chn));
      caltable_cal[i][chn] = (UShort_t)(slope[i]*chn + offset[i]);
    }
  }
}

void S800FpIonChamber::Initialize(S800Full *s800) {
  m_top = s800;

  m_activeChannels = S800_FP_IC_CHANNELS;
  
  /* Parameters */
  for (int i=0; i<S800_FP_IC_CHANNELS; i++) {
    raw[i] = -65000.;
    cal[i] = -65000.;
  }
  tac1 = -65000.;
  tac2 = -65000.;
  sum = -65000.;
  de = -65000.;

  /* Variables */
  de_x0tilt = -300.; /* crdcx */
  de_xtilt = 0.; /* 1/crdcx */
  de_ytilt = 0.; /* 1/crdcx */
  de_slope = 1.; /* MeV/ch */
  de_offset = 0.; /* MeV */

  for (int i=0; i<S800_FP_IC_CHANNELS; i++) {
    slope[i] = 1.; /* au */
    offset[i] = 0.; /* au */
    active[i] = 1.; /* bool */
  }
}

void S800FpIonChamber::Reset() {
  /* Parameters */
  for (int i=0; i<S800_FP_IC_CHANNELS; i++) {
    raw[i] = -65000.;
    cal[i] = -65000.;
  }
  tac1 = -65000.;
  tac2 = -65000.;
  sum = -65000.;
  de = -65000.;
}

UShort_t* S800FpIonChamber::Unpack(UShort_t *p) {
  int SubLen;
  int data;
  double tmpsum = 0;
  
  int chn, caldta;

  p += 2; /* PktLen, PktTag = S800_FP_IC_ENERGY_PACKET */
  
  SubLen = *p; p += 2; /* SubLen, SubTag */
  SubLen -= 2;
  
  if (SubLen) { /* At least one data entry to make sum valid */
    while (SubLen--) {
      data = *p; p++;
      
      chn = data >> 12;
      data &= 0x0fff;
      caldta = caltable_cal[chn][data];
      
      raw[chn] = data;
      cal[chn] = caldta;

      if (active[chn]) {
	tmpsum += caldta;
      }
    }
    
    sum = (double)tmpsum/m_activeChannels;
  }
  
#ifdef CS800_PARSE_IC_TIME_PACKET
  /* Skip this */
  p += *p;
#endif
  
  /* S800_FP_IC_TIME_PACKET not in data stream, omitted!! */
  /* Old style IC packet NOT SUPPORTED */

  return p;
}

void S800FpIonChamber::CalibrateIcSum() {
  if ( (sum != -65000.) && (m_top->fp.track.afp != -65000.) ) {
    de = sum;
    if (m_top->fp.crdc1.y != -65000.) {
      de += de * de_ytilt * m_top->fp.crdc1.y;
    }

    if ( (m_top->fp.crdc1.x != -65000.) &&
	 (m_top->fp.crdc1.x < de_x0tilt) ) {
      de *= exp(de_xtilt * (de_x0tilt - m_top->fp.crdc1.x));
    }

    de = de*de_slope + de_offset;
  } else {
    m_top->InvalidBitReg |= S800_INVALID_FP_IONCH;
  }
}

/************************************************************/
/* S800FpScintillator Class - Functions                     */
/************************************************************/

void S800FpScintillator::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  de_up = -65000.;
  de_down = -65000.;
  time_up = -65000.;
  time_down = -65000.;
  de = -65000.;
  time = -65000.;
  pos = -65000.;

  /* This class has no variables, no need for any
     variable initialization function. */
}

void S800FpScintillator::Reset() {
  /* Parameters */
  de_up = -65000.;
  de_down = -65000.;
  time_up = -65000.;
  time_down = -65000.;
  de = -65000.;
  time = -65000.;
  pos = -65000.;
}

void S800FpScintillator::Calculate() {
  /* I shouldn't calculated the time from only one PMT
     because both times are not calibrated.  Therefore 
     t1 approx (t1+t2)/2 is NOT valid. */
  if ( (time_up != -65000.) && (time_down != -65000.) ) { /* de_xx valid too */
    de = sqrt( (de_up*de_up) + (de_down*de_down) ) / 2;
    time = (time_up + time_down) / 2;
    
    if ( (de_down != 0) && (de_up != 0) ) {
      pos = log(de_up/de_down);
    } else {
      pos = -1001; /* de_up == de_down --> pos is 0 */
    }
  }
}

/************************************************************/
/* S800FpTrack Class - Functions                            */
/************************************************************/

void S800FpTrack::Initialize(S800Full *s800) {
  m_top = s800;

  map.Initialize(m_top);

  /* Parameters */
  xfp = -65000.;
  afp = -65000.;
  yfp = -65000.;
  bfp = -65000.;
  ata = -65000.;
  yta = -65000.;
  bta = -65000.;
  dta = -65000.;
  azita = -65000.;
  scatter = -65000.;
  energy = -65000.;
  ptot = -65000.;
  ppar = -65000.;
  ptra = -65000.;

  ata_cor = -65000.;
  bta_cor = -65000.;
  azita_cor = -65000.;
  scatter_cor = -65000.;

  /* Variables */
  anglea = 0; /* degree */
  angleb = 0; /* degree */

  anglea_cor = 0; /* degree */
  angleb_cor = 0; /* degree */
  bta_ytacor = 0; /* degree/yta */
  ata_dtacor = 0; /* degree/dta */

  brho = 3; /* Tm */
  mass = 1; /* amu */
  charge = 1; /* q */
  order = 5; /* order */
  zfp = 0; /* m */
  gecorr = 0; /* keV/dta */

}

void S800FpTrack::Reset() {
  /* Parameters */
  xfp = -65000.;
  afp = -65000.;
  yfp = -65000.;
  bfp = -65000.;
  ata = -65000.;
  yta = -65000.;
  bta = -65000.;
  dta = -65000.;
  azita = -65000.;
  scatter = -65000.;
  energy = -65000.;
  ptot = -65000.;
  ppar = -65000.;
  ptra = -65000.;

  ata_cor = -65000.;
  bta_cor = -65000.;
  azita_cor = -65000.;
  scatter_cor = -65000.;

  map.Reset();
}

void S800FpTrack::CalculateTracking() {
  double input[S800_TRACK_PARAMETERS];
  double pi = 3.14159265;
  double amu = 931.5016;
  double beta, betagamma0, gamma0, gamma, energy0, ptot0;
  if (!map.WasLoaded()) { return; }

  /* x (dispersive) axis is reversed between COSY, and the
     local coordinate system, where x>0 is down. */
  if (afp == -65000.) { return; }
  
  input[0] = -xfp;
  input[1] = -afp;
  input[2] = yfp;
  input[3] = bfp;
  if ((double)order > (double)map.maxorder) { order = map.maxorder; }

  ata = map.Calculate((int)order, 0, input);
  yta = map.Calculate((int)order, 1, input);
  bta = map.Calculate((int)order, 2, input);
  dta = map.Calculate((int)order, 3, input);
  ata += anglea; /* Add theta offset in radians */
  bta += angleb; /* Add phi offset in radians */
  
  ata_cor = ata + dta*ata_dtacor;
  bta_cor = bta + yta*bta_ytacor;
  
  ata_cor += anglea_cor/180.*pi;
  bta_cor += angleb_cor/180.*pi;
  
  double xsin = sin(ata);
  double ysin = sin(bta);
  
  if ( (xsin > 0) && (ysin > 0) ) {
    azita = atan(ysin/xsin);
  } else if ( (xsin < 0) && (ysin > 0) ) {
    azita = pi - atan(ysin/fabs(xsin));
  } else if ( (xsin < 0) && (ysin < 0) ) {
    azita = pi + atan(fabs(ysin)/fabs(xsin));
  } else if ( (xsin > 0) && (ysin < 0) ) {
    azita = 2*pi - atan(fabs(ysin)/xsin);
  } else {
    azita = 0.0;
  }
  
  scatter = asin(sqrt( (xsin*xsin) + (ysin*ysin) ))*1000; /* mrad */
  
  /* Same for corrected ata_cor and bta_cor */
  xsin = sin(ata_cor);
  ysin = sin(bta_cor);
  
  if ( (xsin > 0) && (ysin > 0) ) {
    azita_cor = atan(ysin/xsin);
  } else if ( (xsin < 0) && (ysin > 0) ) {
    azita_cor = pi - atan(ysin/fabs(xsin));
  } else if ( (xsin < 0) && (ysin < 0) ) {
    azita_cor = pi + atan(fabs(ysin)/fabs(xsin));
  } else if ( (xsin > 0) && (ysin < 0) ) {
    azita_cor = 2*pi - atan(fabs(ysin)/xsin);
  } else {
    azita_cor = 0.0;
  }
  scatter_cor = asin(sqrt( (xsin*xsin) + (ysin*ysin) ))*1000; /* mrad */

  // cout << brho << " " << charge << " " << mass << endl;

  betagamma0 = brho / 3.107 * charge / mass;
  gamma0 = sqrt( (betagamma0 * betagamma0) + 1);
  beta0 = betagamma0 / gamma0;
  energy0 = mass * amu * (gamma0 - 1); /* MeV */
  ptot0 = energy0 * sqrt(1 + 2 * mass * amu / energy0); /* MeV/c */
  energy = (1 + dta) * energy0; /* MeV */
  gamma = 1 + energy / mass / amu;
  beta = sqrt(1 - 1 / gamma / gamma);
  deltabeta = beta;// - beta0;
  ptot = energy * sqrt(1 + 2 * mass * amu / energy); /* MeV/c */
  ppar = ptot * cos(scatter / 1000);
  ptra = ptot * sin(scatter / 1000);
  ptot -= ptot0;
  ppar -= ptot0;
}

/************************************************************/
/* S800FpHodo Class - Functions                             */
/************************************************************/

void S800FpHodo::Initialize(S800Full *s800) {
  m_top = s800;
 
  /* Parameters */
  for (int i=0; i<32; i++) {
    raw[i] = -65000.;
    cal[i] = -65000.;
  }
  max = -65000.;
  sum = -65000.;
  mult = -65000.;

  /* Variables */
  for (int i=0; i<32; i++) {
    slope[i] = 1; /* au */
    offset[i] = 0; /* au */
    xCorr[i] = 0;
  }

}

void S800FpHodo::Reset() {
  /* Parameters */
  for (int i=0; i<32; i++) {
    raw[i] = -65000.;
    cal[i] = -65000.;
  }
  max = -65000.;
  sum = -65000.;
  mult = -65000.;
}

UShort_t* S800FpHodo::Unpack(UShort_t *p) {
  /* HodoTag, HodoEvLen, HodoID1, en, ...., en
     OR
     HodoTag, HodoEvLen, HodoID2, en, ...., en
     Raw and calibrated energies for each channel.
     max shows maximum energy detected in the 32 CsI(Na)
     sum is sum of calibrated energies of all detectors
     mult is multiplicity of detectors */
  UShort_t Len;
  UShort_t id;

  Len = *p; ++p;
  ++p; /* Skipping packet name */
  
  id = *p; ++p;
  Len -= 3;

  if (max == -65000.) { max = 0; }
  if (sum == -65000.) { sum = 0; }
  if (mult == -65000.) { mult = 0; }
  
  while (Len) {
    UShort_t data = *p; ++p;
    Len--;

    UShort_t en = data & 0x0fff;
    UShort_t chn = id*16 + (data >> 12);
    
    raw[chn] = en;
    mult++;
  }

  return p;  
}

void S800FpHodo::CorrectX() {
  for (Int_t i=0; i<32; i++) {
    if (raw[i] > 0) {
      cal[i] = raw[i] - xCorr[i]*m_top->fp.crdc1.x;
      cal[i] *= slope[i];
      cal[i] += offset[i];
      
      if (cal[i] > max) {
	max = cal[i];
      }
      
      if (cal[i] > 100) {
	sum += cal[i];
      }
    }
  }
}

/************************************************************/
/* S800FocalPlane Class - Functions                         */
/************************************************************/

void S800FocalPlane::Initialize(S800Full *s800) {
  m_top = s800;

  e1.Initialize(m_top);
  m_e[0] = &e1;
#ifdef S800_LINK_E2
  e2.Initialize(m_top);
  m_e[1] = &e2;
#endif
#ifdef S800_LINK_E3
  e3.Initialize(m_top);
  m_e[2] = &e3;
#endif

  ic.Initialize(m_top);
  crdc1.Initialize(m_top);
  crdc2.Initialize(m_top);
  track.Initialize(m_top);
  hodo.Initialize(m_top);
  
  /* Variables */
  gap = 1073.0; /* mm */

}

void S800FocalPlane::Reset() {
  e1.Reset();
#ifdef S800_LINK_E2
  e2.Reset();
#endif
#ifdef S800_LINK_E3
  e3.Reset();
#endif

  ic.Reset();
  crdc1.Reset();
  crdc2.Reset();
  track.Reset();
  hodo.Reset();
}

void S800FocalPlane::LoadValues() {
  if ( (crdc1.x != -65000.) && (crdc2.x != -65000.) ) {
    track.afp = atan((crdc2.x - crdc1.x) / gap);
    track.bfp = atan((crdc2.y - crdc1.y) / gap);
    /* xfp, yfp calculated in meters */  
    track.xfp = crdc1.x / 1000 + track.zfp*tan(track.afp);
    track.yfp = crdc1.y / 1000 + track.zfp*tan(track.bfp); 
  } else {
    /* afp, bfp, xfp and yfp are still invalid. */
  }
}

UShort_t* S800FocalPlane::ScintUnpack(UShort_t *p) {
  int PktLen = *p; p += 2;
  /* PktTag = S800_FP_SCINT_PACKET */
  PktLen -= 2;
  
  int value[12] = {-1, -1, -1, -1, -1, -1,
		   -1, -1, -1, -1, -1, -1};
  int PktRead[6]; /* There are up to 6 packets (e1-e3, up/down) */
  int PktCt = 0;

  while (PktLen) {
    int en = *p; ++p;
    int ti = *p; ++p;
    
    int chn = en >> 12;
    value[chn] = en & 0x0fff;
    value[chn+6] = ti & 0x0fff;

    PktLen -= 2;
    PktRead[PktCt] = chn;
    PktCt++;
  }

  if (m_top->trigger.s800 != -65000.) {
    while(PktCt--) {
      int chn = PktRead[PktCt];
#ifndef S800_LINK_E2
      if ((chn>>1) == 1) { continue; }
#endif
#ifndef S800_LINK_E3
      if ((chn>>1) == 2) { continue; }
#endif
      switch(chn % 2) { /* odd:down even:up */
      case 0: /* even (up) */
	m_e[chn >> 1]->de_up = value[chn];
	if (value[chn+6]) {
	  m_e[chn >> 1]->time_up = value[chn+6] - m_top->trigger.s800;
	}
	break;

      default: /* odd (down) */
	m_e[chn >> 1]->de_down = value[chn];
	if (value[chn+6]) {
	  m_e[chn >> 1]->time_down = value[chn+6] - m_top->trigger.s800;
	}
      }
    }

    e1.Calculate();
#ifdef S800_LINK_E2
    e2.Calculate();
#endif
#ifdef S800_LINK_E3
    e3.Calculate();
#endif
  } else {
    /* One could make the energy data valid...but for now, nothing */
  }

  return p;
}

/************************************************************/
/* S800ClassicPpac Class - Functions                        */
/************************************************************/

void S800ClassicPpac::Initialize(S800Full *s800) {
  m_top = s800;
 
  /* Parameters */
  up = -65000.;
  down = -65000.;
  left = -65000.;
  right = -65000.;

  x = -65000.;
  y = -65000.;

  /* Variables */
  z = 0; /* mm */
  u_offset = 0; 
  u_slope = 1;
  d_offset = 0;
  d_slope = 1;
  l_offset = 0;
  l_slope = 1;
  r_offset = 0;
  r_slope = 1;
  x_offset = 0; /* mm */
  x_slope = 1; /* mm/ch */
  y_offset = 0; /* mm */
  y_slope = 1; /* mm/ch */

}

void S800ClassicPpac::Reset() {
  /* Parameters */
  up = -65000.;
  down = -65000.;
  left = -65000.;
  right = -65000.;

  x = -65000.;
  y = -65000.;
}

UShort_t* S800ClassicPpac::Unpack(UShort_t *p, int ChnOffSet) {
  int ch;
  int read = 0;
  ch = ((*p) & 0xf000) >> 12;
  while( (ch >= ChnOffSet) && (ch <= ChnOffSet+3) && (read < 4) ) {
    if (ch == ChnOffSet+0) { up = ((*p++)&0x0fff)*u_slope + u_offset; }
    if (ch == ChnOffSet+1) { down = ((*p++)&0x0fff)*d_slope + d_offset; }
    if (ch == ChnOffSet+2) { left = ((*p++)&0x0fff)*l_slope + l_offset; }
    if (ch == ChnOffSet+3) { right = ((*p++)&0x0fff)*r_slope + r_offset; }

    ch = ((*p)&0xf000) >> 12;
    read++;
  }

  if ( (left != -65000.) && (right != -65000.) ) {
    if ((right + left) > 0) {
      x = ((right - left)/(right + left))*x_slope + x_offset;
    }
  }
  if ( (up != -65000.) && (down != -65000.) ) {
    if ((up + down) > 0) {
      y = ((up - down)/(up + down))*y_slope + y_offset;
    }
  }

  return p;
}

/************************************************************/
/* S800Image Class - Functions                              */
/************************************************************/

void S800Image::Initialize(S800Full *s800) {
  m_top = s800;

#ifdef S800_LINK_IMAGE_TRACK
  track.Initialize(m_top);
#endif

#ifdef S800_LINK_IMAGE_PPAC
  ppac1.Initialize(m_top);
  ppac2.Initialize(m_top);
#endif
  
  /* Parameters */
  x = -65000.; /* mm */
  theta = -65000.; /* mrad */
  y = -65000.; /* mm */
  phi = -65000.; /* mrad */
  dp = -65000.; /* mom */

  /* Variables */
  z = 120.6; /* mm */
  gap = 482.6; /* mm */
  dpphi = 0; /* slope */
  
  trigger = false;
}

void S800Image::Reset() {

#ifdef S800_LINK_IMAGE_TRACK
  track.Reset();
#endif

#ifdef S800_LINK_IMAGE_PPAC
  ppac1.Reset();
  ppac2.Reset();
#endif
  
  /* Parameters */
  x = -65000.; /* mm */
  theta = -65000.; /* mrad */
  y = -65000.; /* mm */
  phi = -65000.; /* mrad */
  dp = -65000.; /* mom */

  trigger = false;
}

void S800Image::Calculate() {
 
#ifdef S800_LINK_IMAGE_TRACK
  if ( (track.tppac1.x != -65000.) && (track.tppac2.x != -65000.) &&
       (track.tppac1.y != -65000.) && (track.tppac2.y != -65000.) ) {
    /* Calculate positions and angles */
    theta = atan((track.tppac2.x - track.tppac1.x)/gap)*1000;
    phi = atan((track.tppac2.y - track.tppac1.y)/gap)*1000;

    x = track.tppac1.x + z*tan(theta/1000);
    y = track.tppac1.y + z*tan(phi/1000);

    /* Load values for tracking calculation */
    track.xim = x;
    track.yim = y;
    track.aim = theta;
    track.bim = phi;
    /* FIX THIS - need to determine the best way to 
       calculate the energy deviation!! */
    track.dim = 0.0;
    /* Call tracking calculation */
    track.CalculateTracking();
  }
#endif
#ifdef S800_LINK_IMAGE_PPAC
  if ( (ppac1.x != -65000.) && (ppac2.x != -65000.) &&
       (ppac1.y != -65000.) && (ppac2.y != -65000.) ) {
    theta = atan((ppac2.x - ppac1.x)/gap)*1000;
    phi = atan((ppac2.y - ppac1.y)/gap)*1000;
    x = ppac1.x + z*tan(theta/1000);
    y = ppac1.y + z*tan(phi/1000);
  }
#endif
  return;
}

/************************************************************/
/* S800TimeStamp Class - Functions                          */
/************************************************************/

void S800TimeStamp::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  stampHigh = -65000.;
  stampMid = -65000.;
  stampLow = -65000.;
  diff = -65000.;

  m_prevStamp = 0;
}

void S800TimeStamp::Reset() {
  /* Parameters */
  stampHigh = -65000.;
  stampMid = -65000.;
  stampLow = -65000.;
  stampLowest = -65000.;
  timestamp = -65000.;
  diff = -65000.;
}

UShort_t* S800TimeStamp::Unpack(UShort_t *p) {
  UShort_t word1, word2, word3, word4;

  p += 2; /* PktLen, PktTag = S800_TIMESTAMP_PACKET 0x5803 */
  
  word1 = *p; ++p;
  word2 = *p; ++p;
  word3 = *p; ++p;
  word4 = *p; ++p;
  
  long long int stamp;
  stamp = (long long int)(word1 + 
			  65536.*word2 + 
			  65536.*65536.*word3 + 
			  65536.*65536.*65536.*word4);
  diff = stamp - m_prevStamp;
  m_prevStamp = (Double_t)stamp;
  
  stampLowest = word1;
  stampLow = word2;
  stampMid = word3;
  stampHigh = word4;

  timestamp = stamp*8;

  return p;  
}

/* S800 Event number */

void S800EventNumber::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  eventNum = -65000;
}

void S800EventNumber::Reset() {
  /* Parameters */
  eventNum = -65000;
}

UShort_t* S800EventNumber::Unpack(UShort_t *p) {
  UShort_t word1, word2, word3;

  p += 2; /* PktLen, PktTag = S800_TIMESTAMP_PACKET 0x5804 */
  
  word1 = *p; ++p;
  word2 = *p; ++p;
  word3 = *p; ++p;
  
  long long int number;
  number = (long long int)(word1 + 
			   65536.*word2 + 
			   65536.*65536.*word3);

  eventNum = number;

  return p;  
}

/************************************************************/
/* S800Trigger Class - Functions                            */
/************************************************************/

void S800Trigger::Initialize(S800Full *m_s800) {
  m_top = m_s800;

  /* Parameters */
  reg = -65000.;
  s800 = -65000.;
  external1 = -65000.;
  external2 = -65000.;
  secondary = -65000.;
#ifdef S800_LINK_TIMESTAMP
  timestamp.Initialize(m_top);
#endif
}

void S800Trigger::Reset() {
  /* Parameters */
  reg = -65000.;
  s800 = -65000.;
  external1 = -65000.;
  external2 = -65000.;
  secondary = -65000.;
#ifdef S800_LINK_TIMESTAMP
  timestamp.Reset();
#endif
}

UShort_t* S800Trigger::Unpack(UShort_t *p) {
  int PktLen = *p; p += 2; /* PktTag = S800_TRIGGER_PACKET */
  reg = *p; ++p;
  PktLen -= 3;
  
  int value[4] = {-1, -1, -1, -1};
  
  m_top->InvalidBitReg = S800_INVALID_NOFEVENTS;
  
  while (PktLen--) {
    int tmp = *p; ++p;
    value[(tmp>>12) - 8] = tmp & 0x0fff;
  }

  if (value[0] != -1) {
    s800 = value[0];
  } else {
    m_top->InvalidBitReg |= S800_INVALID_TRIGGER_S800;
  }

  if (value[1] != -1) {
    external1 = value[1];
  } else {
    m_top->InvalidBitReg |= S800_INVALID_TRIGGER_EXT1;
  }

  if (value[2] != -1) {
    external2 = value[2];
  } else {
    m_top->InvalidBitReg |= S800_INVALID_TRIGGER_EXT2;
  }

  if (value[3] != -1) {
    secondary = value[3];
  } else {
    m_top->InvalidBitReg |= S800_INVALID_TRIGGER_SEC;
  }

  return p;  
}

/************************************************************/
/* S800TimeOfFlight - Functions                             */
/************************************************************/

void S800TimeOfFlight::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  rf = -65000.;
  obj = -65000.;
  xfp = -65000.;
  xfp_obj = -65000.;
  rfe1 = -65000.;
  obje1 = -65000.;
  xfpe1 = -65000.;
#ifdef S800_LINK_E2
  obje2 = -65000.;
  xfpe2 = -65000.;
#endif

#ifdef S800_LINK_TOFTAC
  tac_obj = -65000.;
  tac_obje1 = -65000.;
  tac_xfp = -65000.;
  tac_xfpe1 = -65000.;
#endif

  mesypat = 0;
  mesye1up = -65000.;  mesye1down = -65000.;
  mesyobj = -65000.;  mesyxfp = -65000.;
  mesyrf = -65000.;
  mesycrdc1 = -65000.;  mesycrdc2 = -65000.;
  mesyref = -65000.;
  mesyobje1 = -65000.;  mesyxfpe1 = -65000.;
  mesyrfe1 = -65000.;
  mesyhitsum = -65000.;

  mesyrefchn = S800MESY_E1REF;
  mselobj = -2000.;
  mselxfp = 2500.;
  mselrf = 0.;
  mselref = 0.;

  for (Int_t i=0; i<S800MESY_MAXCHN; i++) { mesyhit[i] = -65000; mesymindiff[i] = -65000.; }
  me1up.clear(); me1down.clear();
  mobj.clear(); mxfp.clear(); mrf.clear();
  mcrdc1.clear(); mcrdc2.clear(); mhoth.clear();
  mref.clear();
  mobje1.clear(); mxfpe1.clear(); mrfe1.clear();

#ifdef S800_LINK_DIAMOND
  diaor = -65000.;
  dia1 = -65000.;
  dia2 = -65000.;
  dia3 = -65000.;
  dia4 = -65000.;
  dia1RF = -65000.;
  dia2RF = -65000.;
  dia3RF = -65000.;
  dia4RF = -65000.;
  diaRF = -65000.;
  dia1Cor = -65000.;
  dia2Cor = -65000.;
  dia3Cor = -65000.;
  dia3Cor = -65000.;
  diaCor = -65000.;
#endif

  /* Variables */
  obj_shift = 0; /* chn */
  rfCorrection = 0; /* rf/dta */
  rfe1Correction = 0; /* rf/afp */
  objCorrection = 0; /* obj/theta */
  diaCorrection = 0; /* dia/theta */
  diax = 0; /* dia/mm */
  dia2ndCor = 0; /* 1/theta */
  dia2ndCorOffset = 0; /* 1/theta */
  dia1offset = 0; /* time */
  dia2offset = 0; /* time */
  dia3offset = 0; /* time */
  dia4offset = 0; /* time */
  xfpCorrection = 0; /* xfp/theta */
  obje1Correction = 0; /* obje1/afp */
  xfpe1Correction = 0; /* xfpe1/afp */
  obje2Correction = 0; /* obje2/afp */
  xfpe2Correction = 0; /* xfpe2/afp */
  im2objCorrection = 0; /* tof/theta */
  
  tacobjCorrection = 0; /* obj/theta */
  tacobje1Correction = 0; /* obje1/afp */
  tacxfpCorrection = 0; /* xfp/theta */
  tacxfpe1Correction = 0; /* xfpe1/afp */  

}

void S800TimeOfFlight::Reset() {
  /* Parameters */
  rf = -65000.;
  obj = -65000.;
  xfp = -65000.;
  xfp_obj = -65000.;
  rfe1 = -65000.;
  obje1 = -65000.;
  xfpe1 = -65000.;
#ifdef S800_LINK_E2
  obje2 = -65000.;
  xfpe2 = -65000.;
#endif

#ifdef S800_LINK_TOFTAC
  tac_obj = -65000.;
  tac_obje1 = -65000.;
  tac_xfp = -65000.;
  tac_xfpe1 = -65000.;
#endif
  
  mesypat = 0;
  mesye1up = -65000.;  mesye1down = -65000.;
  mesyobj = -65000.;  mesyxfp = -65000.;
  mesyrf = -65000.;
  mesycrdc1 = -65000.;  mesycrdc2 = -65000.;
  mesyref = -65000.;
  mesyobje1 = -65000.;  mesyxfpe1 = -65000.;
  mesyrfe1 = -65000.;
  mesyhitsum = -65000.;

  mesyrefchn = S800MESY_E1REF;
  mselobj = -2000.;
  mselxfp = 2500.;
  mselrf = 0.;
  mselref = 0.;

  for (Int_t i=0; i<S800MESY_MAXCHN; i++) { mesyhit[i] = -65000; mesymindiff[i] = -65000.; }
  me1up.clear(); me1down.clear();
  mobj.clear(); mxfp.clear(); mrf.clear();
  mcrdc1.clear(); mcrdc2.clear(); mhoth.clear();
  mref.clear();
  mobje1.clear(); mxfpe1.clear(); mrfe1.clear();
  
#ifdef S800_LINK_DIAMOND
  diaor = -65000.;
  dia1 = -65000.;
  dia2 = -65000.;
  dia3 = -65000.;
  dia4 = -65000.;
  dia1RF = -65000.;
  dia2RF = -65000.;
  dia3RF = -65000.;
  dia4RF = -65000.;
  diaRF = -65000.;
  dia1Cor = -65000.;
  dia2Cor = -65000.;
  dia3Cor = -65000.;
  dia3Cor = -65000.;
  diaCor = -65000.;
#endif
}

UShort_t* S800TimeOfFlight::Unpack(UShort_t *p) {
  int PktLen = *p; p += 2; /* PktTag = S800_TOF_PACKET */
  PktLen -= 2;

  int value[16] = {-1, -1, -1, -1,
		   -1, -1, -1, -1,
		   -1, -1, -1, -1,
		   -1, -1, -1, -1};
  
  while(PktLen--) {
    int tmp = *p; ++p;
    value[tmp >> 12] = tmp & 0x0fff;
  }

  if (m_top->trigger.s800 != -65000.) {
    if (value[12] != -1) { /* Cyclotron RF */
      rf = value[12] - m_top->trigger.s800;
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_RF;
    }

    if (value[13] != -1) { /* Object scintillator */
      obj = value[13] - m_top->trigger.s800 + obj_shift;
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_OBJ;
    }

    if (value[14] != -1) { /* Extended FP scintillator */
      xfp = value[14] - m_top->trigger.s800;
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_XFP;
    }

#ifdef S800_LINK_DIAMOND
    if (value[15] != -1) { /* Diamond detector OR */
      diaor = value[15] - m_top->trigger.s800;
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_DIAOR;
    }

    /* 4 diamonds... */
    if (value[4] != -1) { /* Diamond detector 1 */
      dia1 = value[4] - m_top->trigger.s800 + dia1offset;
      tac_xfp = value[4];
      if (rf != -65000.) {
	dia1RF = dia1 - rf;
	diaRF = dia1RF;
      }
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_DIA1;
    }

    if (value[5] != -1) { /* Diamond detector 2 */
      dia2 = value[5] - m_top->trigger.s800 + dia2offset;
      tac_obj = value[5];
      if (rf != -65000.) {
	dia2RF = dia2 - rf;
	if (diaRF == -65000.) {
	  diaRF = dia2RF;
	} else {
	  diaRF = -65000.;
	}
      }
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_DIA2;
    }

    if (value[6] != -1) { /* Diamond detector 3 */
      dia3 = value[6] - m_top->trigger.s800 + dia3offset;
      if (rf != -65000.) {
	dia3RF = dia3 - rf;
	if (diaRF == -65000.) {
	  diaRF = dia3RF;
	} else {
	  diaRF = -65000.;
	}
      }
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_DIA3;
    }

    if (value[7] != -1) { /* Diamond detector 4 */
      dia4 = value[7] - m_top->trigger.s800 + dia4offset;
      if (rf != -65000.) {
	dia4RF = dia4 - rf;
	if (diaRF == -65000.) {
	  diaRF = dia4RF;
	} else {
	  diaRF = -65000.;
	}
      }
    } else {
      m_top->InvalidBitReg |= S800_INVALID_TOF_DIA4;
    }
#endif

  } /* End of if (s800.trigger.s800 != -65000.) */
  else {
    /* All data is invalid without trigger.s800 */
    m_top->InvalidBitReg |= 
      S800_INVALID_TOF_RF |
      S800_INVALID_TOF_OBJ |
      S800_INVALID_TOF_XFP |
      S800_INVALID_TOF_DIAOR |
      S800_INVALID_TOF_DIA1 |
      S800_INVALID_TOF_DIA2 |
      S800_INVALID_TOF_DIA3 | 
      S800_INVALID_TOF_DIA4;
  }

  return p;  
}

UShort_t* S800TimeOfFlight::UnpackMesy(UShort_t *p) {
  
  int PktLen = *p; p += 2; /* PktTag = S800_MESYTEC_TDC_PACKET */
  PktLen -= 2;

  int32_t traw[S800MESY_MAXCHN][S800MESY_MAXHIT];
  int32_t Tref = 0x8000;
  bool bad = false;

  mesypat = 0;
  mesyhitsum = 0;

  while (PktLen) {
    uint16_t chn, hit;
    chn = *p & 0xff; /* Channel */
    hit = (*p >> 8); /* Hit number 0 - 15 */
    ++p; PktLen--;

    mesypat = (int32_t) mesypat | (int32_t)(1<<chn);
    mesyhitsum++;

    if ( chn >= S800MESY_MAXCHN || hit >= S800MESY_MAXHIT ) {
      p++; PktLen--;
      printf("S800TimeOfFlight::UnpackMesy");
      printf(" channel %d OR hit %d out of range", chn, hit);
      printf(" (MAXCHN: %d, MAXHIT: %d)\n", S800MESY_MAXCHN, S800MESY_MAXHIT);
      continue;
    }

    if (mesyhit[chn] != -65000) { mesyhit[chn]++; }
    else { mesyhit[chn] = 1; }

    if (mesyhit[chn] != hit+1) {
      printf("S800TimeOfFlight::UnpackMesy -- channel %d", chn);
      printf(" hit %d but mesyhit[chn] is %d\n", hit, mesyhit[chn]);
      bad = true;
    }

    traw[chn][hit] = (int32_t)*p;
    ++p; PktLen--;
  
    if (chn == mesyrefchn) {
      if (hit==0) { Tref = traw[chn][hit]; }
      else if (TMath::Abs(traw[chn][hit]-mselref) < TMath::Abs(Tref-mselref)) { Tref = traw[chn][hit]; }
    }
  }

  if (mesyhit[S800MESY_E1REF] > 0) { mesyref = Tref; }

  if (!bad) {
    for (int chn=0; chn<S800MESY_MAXCHN; chn++) {
      if (mesyhit[chn] > 0) {
	for (int hit = 0; hit<mesyhit[chn]; hit++) {
	  /* Find minimum time difference between multiple hits */
	  if (hit == 1) {
	    mesymindiff[chn] = TMath::Abs(traw[chn][1] - traw[chn][0]);
	  }
	  if (hit > 1 && TMath::Abs(traw[chn][hit]-traw[chn][hit-1])<mesymindiff[chn]) {
	    mesymindiff[chn] = TMath::Abs(traw[chn][hit]-traw[chn][hit-1]);
	  }
	
	  switch(chn) {
	  case S800MESY_E1UP:
	    me1up.push_back(traw[chn][hit]-Tref);
	    mesye1up = me1up[0];
	    break;
	  case S800MESY_E1DOWN:
	    me1down.push_back(traw[chn][hit]-Tref);
	    mesye1down = me1down[0];
	    break;
	  case S800MESY_XFP:
	    mxfp.push_back(traw[chn][hit]-Tref);
	    if (hit==0) { mesyxfp = mxfp[0]; }
	    else if (TMath::Abs(mxfp[hit]-mselxfp) < TMath::Abs(mesyxfp-mselxfp)) { mesyxfp = mxfp[hit]; }
	    break;
	  case S800MESY_OBJ:
	    mobj.push_back(traw[chn][hit]-Tref);
	    if (hit==0) { mesyobj = mobj[0]; }
	    else if (TMath::Abs(mobj[hit]-mselobj) < TMath::Abs(mesyobj-mselobj)) { mesyobj = mobj[hit]; }
	    break;
	  case S800MESY_RF:
	    mrf.push_back(traw[chn][hit]-Tref);
	    if (hit==0) { mesyrf = mrf[0]; }
	    else if (TMath::Abs(mrf[hit]-mselrf) < TMath::Abs(mesyrf-mselrf)) { mesyrf = mrf[hit]; }
	    break;
	  case S800MESY_CRDC1:
	    mcrdc1.push_back(traw[chn][hit]-Tref);
	    mesycrdc1 = mcrdc1[0];
	    break;
	  case S800MESY_CRDC2:
	    mcrdc2.push_back(traw[chn][hit]-Tref);
	    mesycrdc2 = mcrdc2[0];
	    break;
	  case S800MESY_HODO:
	    mhoth.push_back(traw[chn][hit]-Tref);
	    break;
	  case S800MESY_E1REF:
	    mref.push_back(traw[chn][hit]);
	    break;
	  default:
	    break;
	  }
	}
      }
    }
  } /* if (!bad) */
  
  return p;
}


void S800TimeOfFlight::CalculateTOF() {
  /* Parameters which need correction with CRDC data */
  if (m_top->fp.track.afp != -65000.) { /* crdc1 and crdc2 are valid too */
#ifdef S800_LINK_DIAMOND
    if (diaor != -65000.) {
      diaor += ( (diaCorrection * m_top->fp.track.afp) + 
		 (diax * m_top->fp.crdc1.x) );
    }
    
    if (dia1 != -65000.) {
      dia1Cor = ( (dia1) + 
		  (diaCorrection * m_top->fp.track.afp) +
		  (diax * m_top->fp.crdc1.x) );
      dia1Cor += ( (dia1Cor - dia2ndCorOffset) *
		   (m_top->fp.track.afp * dia2ndCor) );
      diaCor = dia1Cor;
      diaor = dia1;
    }

    if (dia2 != -65000.) {
      dia2Cor = ( (dia2) + 
		  (diaCorrection * m_top->fp.track.afp) +
		  (diax * m_top->fp.crdc1.x) );
      dia2Cor += ( (dia2Cor - dia2ndCorOffset) *
		   (m_top->fp.track.afp * dia2ndCor) );
      if (diaCor == -65000.) {
	diaCor = dia2Cor;
      } else {
	diaCor = -65000.;
      }
      diaor = dia2;
    }

    if (dia3 != -65000.) {
      dia3Cor = ( (dia3) + 
		  (diaCorrection * m_top->fp.track.afp) +
		  (diax * m_top->fp.crdc1.x) );
      dia3Cor += ( (dia3Cor - dia2ndCorOffset) *
		   (m_top->fp.track.afp * dia2ndCor) );
      if (diaCor == -65000.) {
	diaCor = dia3Cor;
      } else {
	diaCor = -65000.;
      }
      diaor = dia3;
    }

    if (dia4 != -65000.) {
      dia4Cor = ( (dia4) + 
		  (diaCorrection * m_top->fp.track.afp) +
		  (diax * m_top->fp.crdc1.x) );
      dia4Cor += ( (dia4Cor - dia2ndCorOffset) *
		   (m_top->fp.track.afp * dia2ndCor) );
      if (diaCor == -65000.) {
	diaCor = dia4Cor;
      } else {
	diaCor = -65000.;
      }
      diaor = dia4;
    }
#endif

    if (rf != -65000.) {
      rfe1 = ( (rf) + (rfe1Correction * m_top->fp.track.afp) + 
	       (rfCorrection * m_top->fp.crdc1.x) );
    }

    if (obj != -65000.) {
      obje1 = ( (obj) + (obje1Correction * m_top->fp.track.afp) +
		(objCorrection * m_top->fp.crdc1.x) );
#ifdef S800_LINK_E2
      if (m_top->fp.e2.time != -65000.) {
	obje2 = ( (obj) - (m_top->fp.e2.time) +
		  (obje2Correction * m_top->fp.crdc1.x) );
      }
#endif
    }

    if (xfp != -65000.) {
      xfpe1 = ( (xfp) + (xfpe1Correction * m_top->fp.track.afp) +
		(xfpCorrection * m_top->fp.crdc1.x) );
#ifdef S800_LINK_E2
      if (m_top->fp.e2.time != -65000.) {
	xfpe2 = ( (xfp) - (m_top->fp.e2.time) +
		  (xfpe2Correction * m_top->fp.crdc1.calc.x_gravity) );
      }
#endif
    }

#ifdef S800_LINK_TOFTAC
    if (tac_obj != -65000.) {
      tac_obje1 = ( (tac_obj) + (tacobje1Correction * m_top->fp.track.afp) +
		    (tacobjCorrection * m_top->fp.crdc1.x) );
    }

    if (tac_xfp != -65000.) {
      tac_xfpe1 = ( (tac_xfp) + (tacxfpe1Correction * m_top->fp.track.afp) +
		    (tacxfpCorrection * m_top->fp.crdc1.x) );
    }
#endif

  } /* End of if (s800.fp.track.afp != -65000.) */

}

void S800TimeOfFlight::CalculateMesyTOF() {
  /* Parameters which need correction with CRDC data */
  if (m_top->fp.track.afp != -65000.) { /* crdc1 and crdc2 are valid too */

    if (mesyobj != -65000.) {
      mesyobje1 = ( (mesyobj) + (mobje1Correction*m_top->fp.track.afp) + 
		    (mobjCorrection*m_top->fp.crdc1.x) ); 
    }
    if (mesyxfp != -65000.) {
      mesyxfpe1 = ( (mesyxfp) + (mxfpe1Correction*m_top->fp.track.afp) + 
		    (mxfpCorrection*m_top->fp.crdc1.x) ); 
    }
    if (mesyrf != -65000.) {
      mesyrfe1 = ( (mesyrf) + (mrfe1Correction*m_top->fp.track.afp) + 
		    (mrfCorrection*m_top->fp.crdc1.x) ); 
    }
    for (int i=0; i<mobj.size(); i++) {
      mobje1.push_back( (mobj[i]) + (mobje1Correction*m_top->fp.track.afp) +
			(mobjCorrection*m_top->fp.crdc1.x) );
    }
    for (int i=0; i<mxfp.size(); i++) {
      mxfpe1.push_back( (mxfp[i]) + (mxfpe1Correction*m_top->fp.track.afp) +
			(mxfpCorrection*m_top->fp.crdc1.x) );
    }
    for (int i=0; i<mrf.size(); i++) {
      mrfe1.push_back( (mrf[i]) + (mrfe1Correction*m_top->fp.track.afp) +
			(mrfCorrection*m_top->fp.crdc1.x) );
    }
  }
}
	
   
 

/************************************************************/
/* S800Pin - Functions                                      */
/************************************************************/

void S800Pin::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  e = -65000.;
  t = -65000.;
  ecorr = -65000.;
  dE = -65000.; /* MeV */
  
  /* Variables */
  e_offset = 0; /* MeV */
  e_slope = 1; /* MeV/ch */
  
  pin_cor = 0; 
  
}

void S800Pin::Reset() {
  /* Parameters */
  e = -65000.;
  t = -65000.;
  ecorr = -65000.;
  dE = -65000.; /* MeV */
}

UShort_t* S800Pin::Unpack(UShort_t *p) {
  int PktLen = *p; 
  p+=2;
  PktLen -= 2;

  if (PktLen) {
    e = *p;
    ++p;
    dE = e*e_slope + e_offset; /* ToDo: Lookup table */
  }

  return p;
}

void S800Pin::PinCorrect() {
  if ( (e != -65000.) && (m_top->fp.crdc1.x != -65000.) ) {
    ecorr = ( (e) + (pin_cor * m_top->fp.crdc1.x) );
  }
}

/************************************************************/
/* S800LaBr - Functions                                      */
/************************************************************/

void S800LaBr::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  for (int index = 0; index<16; index++) {
    energy[index] = -65000.;
    time[index] = -65000.;
  }
}

void S800LaBr::Reset() {
  for (int index = 0; index<16; index++) {
    energy[index] = -65000.;
    time[index] = -65000.;
  }
}

/************************************************************/
/* S800Object - Functions                                   */
/************************************************************/

void S800Object::Initialize(S800Full *s800) {
  m_top = s800;
  
#ifdef S800_LINK_OBJECT_PIN
  pin.Initialize(m_top);
#endif
}

void S800Object::Reset() {
#ifdef S800_LINK_OBJECT_PIN
  pin.Reset();
#endif
}

/************************************************************/
/* S800Ppac - Functions                                     */
/************************************************************/

void S800Ppac::Initialize(S800Full *s800) {
  m_top = s800;

  /* Parameters */
  up = -65000.;
  down = -65000.;
  left = -65000.;
  right = -65000.;
  a = -65000.;
  time = -65000.;
  sumx = -65000.;
  sumy = -65000.;
  x = -65000.;
  y = -65000.;

  /* Variables */
  z = 0; /* mm */
  u_offset = 0; 
  u_slope = 1;
  d_offset = 0; 
  d_slope = 1;
  l_offset = 0; 
  l_slope = 1;
  r_offset = 0; 
  r_slope = 1;
  x_offset = 0; /* mm */
  x_slope = 1; /* mm/ch */
  y_offset = 0; /* mm */
  y_slope = 1; /* mm/ch */

}

void S800Ppac::Reset() {
  /* Parameters */
  up = -65000.;
  down = -65000.;
  left = -65000.;
  right = -65000.;
  a = -65000.;
  time = -65000.;
  sumx = -65000.;
  sumy = -65000.;
  x = -65000.;
  y = -65000.;
}

UShort_t* S800Ppac::Unpack(UShort_t *p) {
  up = (*p++) & 0x0fff;
  down = (*p++) & 0x0fff;
  left = (*p++) & 0x0fff;
  right = (*p++) & 0x0fff;
  time = (*p++) & 0x0fff;

  return p;
} 

void S800Ppac::Calibrate() {
  left = (left * l_slope) + l_offset;
  right = (right * r_slope) + r_offset;
  up = (up * u_slope) + u_offset;
  down = (down * d_slope) + d_offset;

  sumx = left + right;
  sumy = up + down;
}

void S800Ppac::CalculatePosition() {
  if (sumx > 0) {
    x = (((right - left) / sumx) * x_slope) + x_offset;
  }
  if (sumy > 0) {
    y = (((up - down) / sumy) * y_slope) + y_offset;
  }
}

/************************************************************/
/* S800Target - Functions                                   */
/************************************************************/

void S800Target::Initialize(S800Full *s800) {
  m_top = s800;

#ifdef S800_LINK_TARGET_PPAC
  ppac2.Initialize(m_top);
#endif
  
#ifdef S800_LINK_TARGET_PIN
  pin1.Initialize(m_top);
  pin2.Initialize(m_top);
#endif

#ifdef S800_LINK_TARGET
  E = -65000.;
#endif
  
}

void S800Target::Reset() {
  
#ifdef S800_LINK_TARGET_PPAC
  ppac2.Reset();
#endif
  
#ifdef S800_LINK_TARGET_PIN
  pin1.Reset();
  pin2.Reset();
#endif

#ifdef S800_LINK_TARGET
  E = -65000.;
#endif
  
}

void S800Target::Calculate() {
#ifdef S800_LINK_TARGET_PIN
  E = pin1.dE + pin2.dE;
#endif
}

/************************************************************/
/* S800 - Functions                                         */
/************************************************************/

S800Full::S800Full() {
  /* Variables */
  debugMode = 0; /* debug */
  evtreset = 0; /* bool */

  /* Parameters */
  invalid = -65000; /* Pattern */
  evtlength = -65000.; /* words */
  evtcount = -65000.; /* # */
}

void S800Full::Initialize() {

  /* Variables */
  debugMode = 0; /* debug */
  evtreset = 0; /* bool */

  /* Parameters */
  invalid = -65000; /* Pattern */
  evtlength = -65000.; /* words */
  evtcount = -65000.; /* # */
  
  evtNumber = 0;
  lastEvtNumber = 0;
  reportEvtIncrease = 1;
  
  fp.Initialize(this);
  target.Initialize(this);
  im.Initialize(this);
  obj.Initialize(this);
  labr.Initialize(this);
  trigger.Initialize(this);
  tof.Initialize(this);
#ifndef S800_LINK_TIMESTAMP
  ts.Initialize(this);
#endif
  evtnum.Initialize(this);

}

void S800Full::Reset() {
  /* Parameters */
  invalid = -65000; /* Pattern */
  evtlength = -65000.; /* words */
  evtcount = -65000.; /* # */

  evtNumber = 0;
  lastEvtNumber = 0;
  reportEvtIncrease = 1;
  
  fp.Reset();
  target.Reset();
  im.Reset();
  obj.Reset();
  labr.Reset();
  trigger.Reset();
  tof.Reset();
#ifndef S800_LINK_TIMESTAMP
  ts.Reset();
#endif
  evtnum.Reset();
}

void S800Full::DoInvalidBitReg() {
  invalid = InvalidBitReg; /* 32 bit register */
  
  /* TOF and TRIGGER treated in Unpacker() methods --
     invalid, if value not in event, or trigger.s800 is missing */
  
  /* Focal plane scintillators s800.fp.eX....
     invalid, if one or both PMT times are missing */
  if (fp.e1.time == -65000.) { invalid |= S800_INVALID_FP_SCINT_E1; }
#ifdef S800_LINK_E2
  if (fp.e2.time == -65000.) { invalid |= S800_INVALID_FP_SCINT_E2; }
#endif
#ifdef S800_LINK_E3
  if (fp.e3.time == -65000.) { invalid |= S800_INVALID_FP_SCINT_E3; }
#endif

  /* CRDC s800.fp.crdcX...
     A tac information for y is available in each event */
  if (fp.crdc1.x == -65000.) { invalid |= S800_INVALID_FP_CRDC1; }
  if (fp.crdc2.x == -65000.) { invalid |= S800_INVALID_FP_CRDC2; }
  if ( (fp.crdc1.x == -65000.) || (fp.crdc2.x == -65000.) ) { 
    invalid |= S800_INVALID_FP_CRDCS; 
  }

  /* Ion chamber s800.fp.ic...
     Invalid bit done in S800FpIonChamber::CalibrateIcSum() */
  
  /* Tracking ppac...
     Invalid if x or y position is missing */
#ifdef S800_LINK_IMAGE_TRACK
  if ( (im.trac.tppac1.x == -65000.) || (im.trac.tppac1.y == -65000.) ) {
    invalid |= S800_INVALID_IM_TPPAC1;
  }
  if ( (im.trac.tppac2.x == -65000.) || (im.trac.tppac2.y == -65000.) ) {
    invalid |= S800_INVALID_IM_TPPAC2;
  }
#endif

}

void S800Full::getAndProcessS800(FILE *inf, Int_t length) {

  //unsigned short *buf, int nw) {
  Reset();
  memset(s800data, 0, sizeof(s800data));
  
  Int_t siz = fread(s800data, 1, length, inf);
  if (siz != length) {
    cout << ALERTTEXT;
    printf("getAndProcessS800() failed in bytes read.\n");
    cout << RESET_COLOR;  fflush(stdout);
  }

  UShort_t *p = s800data;
  p++; p++;

  int words = *p; /* Word count is inclusive */

  //  int words = nw;
  
  while (words) {
    UShort_t PktLen = *p; p++;
    UShort_t PktTag = *p; p++;

    int PktWords = PktLen;
    
    if (PktTag == S800_PACKET) {
      if (DEBUGS800) { cout << "Got S800 packet...unpacking!" << endl; }
      /* Start parsing the S800 event. A fixed order
	 for different subpackets is assumed. */
      p++; /* S800 version tag */
      PktWords -= 3; /* PktTag, PktLen and S800 version */
      evtlength = PktLen;

      if (DEBUGS800) { printf("PktLen: %d\n", PktLen); }

      while (PktWords) {
	
	UShort_t SubPktLen = *p;
	UShort_t SubPktTag = *(p+1);

	if (DEBUGS800) { printf("SubPktTag = 0x%x\n", SubPktTag); }
	if (DEBUGS800) { printf("SubPktLen = %d\n", SubPktLen); cin.get(); }

	switch(SubPktTag) {
	case S800_TRIGGER_PACKET :
	  if (DEBUGS800) { cout << "S800 Trigger packet " << endl; }
	  p = trigger.Unpack(p);
	  break;
	  
	  //case S800_EVENT_NUMBER_PACKET :
	  //if (DEBUGS800) { cout << "S800 Event Number packet " << endl; }
	  //p = evtnum.Unpack(p);
	  //break;

	case S800_TOF_PACKET :
	  if (DEBUGS800) { cout << "S800 TOF packet " << endl; }
	  p = tof.Unpack(p);
	  break;
	  
	case S800_FP_SCINT_PACKET :
	  if (DEBUGS800) { cout << "S800 FP Scint packet " << endl; }
	  p = fp.ScintUnpack(p);
	  break;
	  
	case S800_FP_IC_PACKET :
	  if (DEBUGS800) { cout << "S800 IC packet " << endl; }
	  p = fp.ic.Unpack(p);
	  break;
	  
	case S800_FP_IC_TIME_PACKET :
	  if (DEBUGS800) { cout << "S800 FP IC time packet " << endl; }
	  break;
	  
	case S800_FP_CRDC_PACKET :
	  if (DEBUGS800) { cout << "S800 FP CRDC packet " << endl; }
	  if (*(p+2) == 0) {
	    p = fp.crdc1.Unpack(p);
	  } else {
	    p = fp.crdc2.Unpack(p);
	  }
	  break;
	  
	case S800_MESYTEC_TDC_PACKET :
	  if (DEBUGS800) { cout << "S800 FP MESYTOF packet " << endl; }
	  p = tof.UnpackMesy(p);
	  break;

	case S800_FP_HODO_PACKET : 
	  if (DEBUGS800) { cout << "S800 FP Hodoscope packet " << endl; }
#ifdef S800_PARSE_HODO
#ifdef S800_LINK_HODO
	  p = fp.hodo.Unpack(p);
#else 
	  p += SubPktLen;
#endif
#else
	  p += SubPktLen;
#endif
	  break;

	case S800_II_TRACK_PACKET :
	  if (DEBUGS800) { cout << "S800 II Track packet " << endl; }
#ifdef S800_PARSE_IMAGE_TRACK
#ifdef S800_LINK_IMAGE_TRACK
	  p = im.track.UnpackRawData(p);
#else
	  p += SubPktLen;
#endif
#else
	  p += SubPktLen;
#endif
	  break;

	case S800_OB_PIN_PACKET :
	  if (DEBUGS800) { cout << "S800 OB PIN packet " << endl; }
#ifdef S800_PARSE_OBJECT_PIN
#ifdef S800_LINK_OBJECT_PIN
	  p = obj.pin.Unpack(p);
#else
	  p += SubPktLen;
#endif
#else
	  p += SubPktLen;
#endif
	  break;
	  
	case S800_TIMESTAMP_PACKET :
	  if (DEBUGS800) { cout << "S800 Timestamp packet " << endl; }
#ifdef S800_PARSE_TIMESTAMP
#ifdef S800_LINK_TIMESTAMP
	  p = trigger.timestamp.Unpack(p);
#else
	  p = ts.Unpack(p);
#endif
#endif
	  break;
	  
	case S800_II_PPAC_PACKET :
	  if (DEBUGS800) { cout << "S800 II PPAC packet " << endl; }
#ifdef S800_PARSE_IMAGE_PPAC
#ifdef S800_LINK_IMAGE_PPAC
	  p += 2;
	  if (SubPktLen != 2) {
	    p = im.ppac1.Unpack(p, 0);
	    p = im.ppac2.Unpack(p, 4);		  
	  }
#else
	  p += SubPktLen;
#endif
#else
	  p += SubPktLen;
#endif
	  break;

	case S800_EVENTNUMBER_PACKET :
	  {
	    uint64_t uint64_evt;
	    p += 2;
	    
	    uint64_evt = (uint64_t) *p;
	    ++p;
	    uint64_evt += (uint64_t) (*p) << 16;
	    ++p;
	    uint64_evt += (uint64_t) (*p) << 32;
	    ++p;
	    evtNumber = (Double_t)uint64_evt;

	    if (lastEvtNumber && reportEvtIncrease!= 0) {
	      if ((lastEvtNumber+1) != uint64_evt) {
		printf("S800 event number increase error!\n");
		printf("  Last: %f, current: %lu\n", lastEvtNumber, uint64_evt);
	      }
	      lastEvtNumber = uint64_evt;
	    }
	  }
	  break;

	case S800_GALOTTE_PACKET : /* This is a MCP ? */
	  {
	    if (DEBUGS800) { cout << "S800 GALOTTE packet " << endl; }
	    p += 2; /* Length and tag. */
	    UShort_t len = SubPktLen;
	    len -= 2;
	    while (len) {
	      switch((*p) >> 12) {
	      case 0:
		tof.dia1 = (*p) & 0x0fff;
		break;
	      case 1:
		tof.dia2 = (*p) & 0x0fff;
		break;
	      case 2:
		tof.dia3 = (*p) & 0x0fff;
		break;
	      case 3:
		tof.dia4 = (*p) & 0x0fff;
		break;
	      default:
		break;
		//cout << "Suspicious channel in GALOTTE detector: "
		//    << ((*p) >> 12) << endl;
	      }
	      p++;
	      len--;
	    }
	  }
	
	  break;

	case S800_LABR_PACKET :
	  {
	    p += 2; /* Length and tag. */
	    UShort_t len = SubPktLen;
	    len -= 2;
	    while (len) {
	      labr.energy[(Int_t)(*p>>12)] = (*p) & 0x0fff;
	      p++;
	      labr.time[(Int_t)(*p>>12)] = (*p) & 0x0fff;
	      p++;
	    }
	  }
	  
	  break;

	case 0x5805 :
	  {
	    p += SubPktLen;
	  }
	  break;

	default:
	  cout << "Unexpected SubPktTag!! 0x" << std::hex << SubPktTag << std::dec << endl;
	  p += SubPktLen;
	  break;
	  
	} /* End of SubPktTag switch */

	PktWords -= SubPktLen;

      }

      words -= PktLen;

#ifdef S800_PARSE_TARGET_PPAC
#error S800_PARSE_TARGET_PPAC not implemented so far
#endif
      
#ifdef S800_PARSE_TARGET_PIN
#error S800_PARSE_TARGET_PIN not implemented so far
#endif

      /* Now some calculations...the order is important! */
      
      fp.LoadValues(); /* Does angles from CRDC and 
				position at the focal plane */

      fp.track.CalculateTracking(); /* Using a map... */
      im.Calculate();
      tof.CalculateTOF();
      tof.CalculateMesyTOF();
      fp.ic.CalibrateIcSum();
      fp.hodo.CorrectX();
      DoInvalidBitReg();

    } else {
      p += PktLen - 2;
      words -= PktLen;
    }

  } /* while(words) */
}

Double_t S800Full::GetNextValue(FILE *file) {
 
  /* This function assumes the line is a string descriptor
     and a value -- returns the value. */
  char line[100];
  char junk[300];
  double value;

  char* retVal;
  
  retVal = fgets(line, 100, file);
  sscanf(line, "%s %lf", junk, &value);
  return value;
}

void S800Full::InitializeS800Variables(TString inputFilename) {

  FILE *in;
  if ( (in = fopen(inputFilename.Data(), "r")) == NULL ) {
    perror("Input variable file could not be opened.");
    return;
  }

  cout << "Reading variable initialization file: "
       << inputFilename.Data() << endl;

  char line[300];
  char filename5[300];
  double value[20];
  int intvalue[5];

  char* retVal;

  while ( !feof(in) ) {
    retVal = fgets(line, 300, in);
    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, "S800", 4) == 0) {
      debugMode = GetNextValue(in);
      evtreset = GetNextValue(in);

      retVal = fgets(line, 300, in); /* Comment */
      fp.gap = GetNextValue(in);

      retVal = fgets(line, 300, in); /* Comment */
      for (int i=0; i<S800_FP_IC_CHANNELS; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%d %lf %lf %lf", &intvalue[0], &value[0], &value[1], &value[2]);
	fp.ic.slope[intvalue[0]] = value[0];
	fp.ic.offset[intvalue[0]] = value[1];
	fp.ic.active[intvalue[0]] = value[2];
      }
      fp.ic.de_x0tilt = GetNextValue(in);
      fp.ic.de_xtilt = GetNextValue(in);
      fp.ic.de_ytilt = GetNextValue(in);
      fp.ic.de_slope = GetNextValue(in);
      fp.ic.de_offset = GetNextValue(in);

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.x_offset = value[0];
      fp.crdc2.x_offset = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.x_slope = value[0];
      fp.crdc2.x_slope = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.y_offset = value[0];
      fp.crdc2.y_offset = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.y_slope = value[0];
      fp.crdc2.y_slope = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.xy_offset = value[0];
      fp.crdc2.xy_offset = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.xy_slope = value[0];
      fp.crdc2.xy_slope = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.calc.gravity_width = value[0];
      fp.crdc2.calc.gravity_width = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.calc.fit_width = value[0];
      fp.crdc2.calc.fit_width = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.calc.method = value[0];
      fp.crdc2.calc.method = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.calc.saturation = value[0];
      fp.crdc2.calc.saturation = value[1];

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      fp.crdc1.calc.badpads = value[0];
      fp.crdc2.calc.badpads = value[1];

      retVal = fgets(line, 300, in); /* Blank line */
      retVal = fgets(line, 300, in); /* Comment */

      for (int i=0; i<10; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[0],
	       &value[1], &value[2], &value[3], &value[4], &value[5], 
	       &value[6], &value[7], &value[8], &value[9]);
	for (int j=0; j<10; j++) {
	  if ( ((i*10) + j) < 99) {
	    fp.crdc1.calc.badpad[(i*10) + j] = value[j];
	  }
	}
      }

      retVal = fgets(line, 300, in); /* Blank line */
      retVal = fgets(line, 300, in); /* Comment */
      for (int i=0; i<10; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[0],
	       &value[1], &value[2], &value[3], &value[4], &value[5], 
	       &value[6], &value[7], &value[8], &value[9]);
	for (int j=0; j<10; j++) {
	  if ( ((i*10) + j) < 99) {
	    fp.crdc2.calc.badpad[(i*10) + j] = value[j];
	  }
	}
      }

      retVal = fgets(line, 300, in); /* Blank line */
      retVal = fgets(line, 300, in); /* Comment */
      for (int i=0; i<S800_FP_CRDC_CHANNELS; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%d %lf %lf %lf %lf %lf", &intvalue[0], &value[0], &value[1],
	       &value[2], &value[3], &value[4]);
	fp.crdc1.pad.ped[intvalue[0]] = value[0];
	fp.crdc1.pad.slope[intvalue[0]] = value[1];
	fp.crdc1.pad.offset[intvalue[0]] = value[2];
	fp.crdc1.pad.minsamples[intvalue[0]] = value[3];
	fp.crdc1.pad.m_thresh[intvalue[0]] = value[4];
      }
      
      retVal = fgets(line, 300, in); /* Comment */
      for (int i=0; i<S800_FP_CRDC_CHANNELS; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%d %lf %lf %lf %lf %lf", &intvalue[0], &value[0], &value[1],
	       &value[2], &value[3], &value[4]);
	fp.crdc2.pad.ped[intvalue[0]] = value[0];
	fp.crdc2.pad.slope[intvalue[0]] = value[1];
	fp.crdc2.pad.offset[intvalue[0]] = value[2];
	fp.crdc2.pad.minsamples[intvalue[0]] = value[3];
	fp.crdc2.pad.m_thresh[intvalue[0]] = value[4];
      }

      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%s", filename5);
      fp.track.map.mapFilename = filename5;
      cout << "fpMapFilename " << fp.track.map.mapFilename << endl;
      fp.track.anglea = GetNextValue(in);
      fp.track.angleb = GetNextValue(in);
      fp.track.brho = GetNextValue(in);
      fp.track.mass = GetNextValue(in);
      fp.track.charge = GetNextValue(in);
      fp.track.order = GetNextValue(in);
      fp.track.zfp = GetNextValue(in);
      fp.track.gecorr = GetNextValue(in);
      fp.track.anglea_cor = GetNextValue(in);
      fp.track.angleb_cor = GetNextValue(in);
      fp.track.ata_dtacor = GetNextValue(in);
      fp.track.bta_ytacor = GetNextValue(in);
      fp.track.map.maxcoefficients = GetNextValue(in);
      fp.track.map.maxparameters = GetNextValue(in);
      fp.track.map.maxorder = GetNextValue(in);

      retVal = fgets(line, 300, in); /* Comment */
      for (int i=0; i<32; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%d %lf %lf %lf", &intvalue[0], &value[0], &value[1], &value[2]);
	fp.hodo.xCorr[intvalue[0]] = value[0];
	fp.hodo.slope[intvalue[0]] = value[1];
	fp.hodo.offset[intvalue[0]] = value[2];
      }

      retVal = fgets(line, 300, in); /* Comment */
#ifdef S800_LINK_TARGET_PIN
      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf", &value[0], &value[1], &value[2]);
      target.pin1.e_offset[0] = value[0];
      target.pin1.e_slope[0] = value[1];
      target.pin1.pin_cor[0] = value[2];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf", &value[0], &value[1], &value[2]);
      target.pin2.e_offset[1] = value[0];
      target.pin2.e_slope[1] = value[1];
      target.pin2.pin_cor[1] = value[2];
#endif
#ifdef S800_LINK_TARGET_PPAC
      retVal = fgets(line, 300, in); /* Comment */
      target.ppac2.z = GetNextValue(in);
      target.ppac2.uoffset = GetNextValue(in);
      target.ppac2.uslope = GetNextValue(in);
      target.ppac2.doffset = GetNextValue(in);
      target.ppac2.dslope = GetNextValue(in);
      target.ppac2.loffset = GetNextValue(in);
      target.ppac2.lslope = GetNextValue(in);
      target.ppac2.roffset = GetNextValue(in);
      target.ppac2.rslope = GetNextValue(in);
      target.ppac2.xoffset = GetNextValue(in);
      target.ppac2.xslope = GetNextValue(in);
      target.ppac2.yoffset = GetNextValue(in);
      target.ppac2.yslope = GetNextValue(in);
#endif

      retVal = fgets(line, 300, in); /* Comment */
      im.z = GetNextValue(in);
      im.gap = GetNextValue(in);
      im.dpphi = GetNextValue(in);
#ifdef S800_LINK_IMAGE_TRACK
#ifdef S800_LINK_IMAGE_PPAC
#error Tppac OR ppac, not both...
#endif
      retVal = fgets(line, 300, in); /* Comment */
      im.track.order = GetNextValue(in);
      retVal = fgets(line, 300, in);
      sscanf(line, "%s", filename6);
      im.track.map.imMapFilename = filename6;
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.track.tppac1.badxstrips = value[0];
      im.track.tppac2.badxstrips = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[0], 
	     &value[1], &value[2], &value[3], &value[4], &value[5], &value[6], 
	     &value[7]; &value[8], &value[9]);
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[10], 
	     &value[11], &value[12], &value[13], &value[14], &value[15], 
	     &value[16], &value[17]; &value[18], &value[19]);
      for (int i=0; i<20; i++) {
	im.track.tppac1.badxstrip[i] = value[i];
      }
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[0], 
	     &value[1], &value[2], &value[3], &value[4], &value[5], &value[6], 
	     &value[7]; &value[8], &value[9]);
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[10], 
	     &value[11], &value[12], &value[13], &value[14], &value[15], 
	     &value[16], &value[17]; &value[18], &value[19]);
      for (int i=0; i<20; i++) {
	im.track.tppac2.badxstrip[i] = value[i];
      }	
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.track.tppac1.badystrips = value[0];
      im.track.tppac2.badystrips = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[0], 
	     &value[1], &value[2], &value[3], &value[4], &value[5], &value[6], 
	     &value[7]; &value[8], &value[9]);
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[10], 
	     &value[11], &value[12], &value[13], &value[14], &value[15], 
	     &value[16], &value[17]; &value[18], &value[19]);
      for (int i=0; i<20; i++) {
	im.track.tppac1.badystrip[i] = value[i];
      }
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[0], 
	     &value[1], &value[2], &value[3], &value[4], &value[5], &value[6], 
	     &value[7]; &value[8], &value[9]);
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &value[10], 
	     &value[11], &value[12], &value[13], &value[14], &value[15], 
	     &value[16], &value[17]; &value[18], &value[19]);
      for (int i=0; i<20; i++) {
	im.track.tppac2.badystrip[i] = value[i];
      }		      
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.track.tppac1.width = value[0];
      im.track.tppac2.width = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.track.tppac1.xcenter = value[0];
      im.track.tppac2.xcenter = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.track.tppac1.ycenter = value[0];
      im.track.tppac2.ycenter = value[1];
      im.track.map.maxcoefficients = GetNextValue(in);
      im.track.map.maxparameters = GetNextValue(in);
      im.track.map.maxorder = GetNextValue(in);
#endif
#ifdef S800_LINK_IMAGE_PPAC
      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.z = value[0];
      im.ppac2.z = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.u_offset = value[0];
      im.ppac2.u_offset = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.u_slope = value[0];
      im.ppac2.u_slope = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.d_offset = value[0];
      im.ppac2.d_offset = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.d_slope = value[0];
      im.ppac2.d_slope = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.l_offset = value[0];
      im.ppac2.l_offset = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.l_slope = value[0];
      im.ppac2.l_slope = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.r_offset = value[0];
      im.ppac2.r_offset = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.r_slope = value[0];
      im.ppac2.r_slope = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.x_offset = value[0];
      im.ppac2.x_offset = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.x_slope = value[0];
      im.ppac2.x_slope = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.y_offset = value[0];
      im.ppac2.y_offset = value[1];
      retVal = fgets(line, 300, in);
      sscanf(line, "%lf %lf", &value[0], &value[1]);
      im.ppac1.y_slope = value[0];
      im.ppac2.y_slope = value[1];
#endif
#ifdef S800_LINK_OBJECT_PIN
      retVal = fgets(line, 300, in); /* Comment */
      obj.pin.e_offset = GetNextValue(in);
      obj.pin.e_slope = GetNextValue(in);
      obj.pin.pin_cor = GetNextValue(in);
#endif

      retVal = fgets(line, 300, in); /* Comment */      
      tof.obj_shift = GetNextValue(in);
      tof.rfCorrection = GetNextValue(in);
      tof.rfe1Correction = GetNextValue(in);
      tof.objCorrection = GetNextValue(in);
      tof.obje1Correction = GetNextValue(in);
      tof.xfpCorrection = GetNextValue(in);
      tof.xfpe1Correction = GetNextValue(in); 
      tof.mrfCorrection = GetNextValue(in);
      tof.mrfe1Correction = GetNextValue(in);
      tof.mobjCorrection = GetNextValue(in);
      tof.mobje1Correction = GetNextValue(in);
      tof.mxfpCorrection = GetNextValue(in);
      tof.mxfpe1Correction = GetNextValue(in); 
      tof.tacobjCorrection = GetNextValue(in);
      tof.tacobje1Correction = GetNextValue(in);
      tof.tacxfpCorrection = GetNextValue(in);
      tof.tacxfpe1Correction = GetNextValue(in);
      tof.diaCorrection = GetNextValue(in);
      tof.diax = GetNextValue(in);
      tof.dia2ndCor = GetNextValue(in);
      tof.dia2ndCorOffset = GetNextValue(in);
      tof.dia1offset = GetNextValue(in);
      tof.dia2offset = GetNextValue(in);
      tof.dia3offset = GetNextValue(in);
      tof.dia4offset = GetNextValue(in);
      tof.obje2Correction = GetNextValue(in);
      tof.xfpe2Correction = GetNextValue(in);
      tof.im2objCorrection = GetNextValue(in);
    }
  }

  fclose(in);
}

/****************************************************/

void S800Full::UpdateS800RunVariables(TString filename) {
  FILE *input;
  if ( (input = fopen(filename.Data(), "r")) == NULL ) {
    cout << "No run variable file was found.  Moving on without it. " << endl;
    return;
  }
  
  cout << "Run variables set from : " << endl;
  cout << "    " << filename.Data() << endl;
  
  char line[300];
  char junk[300];
  Double_t value;
  Int_t chn;

  while (!feof(input)) {
    char* str = fgets(line, 300, input);
    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, "fp.crdc1.y_slope", 16) == 0) {
      sscanf(line, "%s %lf", junk, &value);
      fp.crdc1.y_slope = value;
    }
    if (strncmp(line, "fp.crdc2.y_slope", 16) == 0) {
      sscanf(line, "%s %lf", junk, &value);
      fp.crdc2.y_slope = value;
    }
    if (strncmp(line, "fp.ic.slope", 11) == 0) {
      sscanf(line, "%s %d %lf", junk, &chn, &value);
      fp.ic.slope[chn] = value;
    }
  }
  fclose(input);
}

Float_t S800Full::getDoppler(TVector3 xyz, Float_t beta, GRETINAVariables *gVar) {
 
  Float_t gamma = 1/TMath::Sqrt(1. - beta*beta);

  if (fp.track.yta < -1000) {
    xyz.SetY(xyz.Y() + fp.track.yta*1000.);
  }

  Float_t xSin = 0, ySin = 0;
  Float_t s800Theta, s800Phi;
 
  xSin = gVar->ataCor*TMath::DegToRad();
  ySin = gVar->btaCor*TMath::DegToRad();
  if (fp.track.ata > -1000.) { /* Valid ata parameter value */
    xSin += fp.track.ata;
    ySin += fp.track.bta;
  }
  xSin = TMath::Sin(xSin);  ySin = TMath::Sin(ySin);

  if (xSin > 0 && ySin > 0) {
    s800Phi = -TMath::ATan(ySin/xSin);
  } else if (xSin < 0 && ySin > 0) {
    s800Phi = -TMath::Pi() + TMath::ATan(ySin/TMath::Abs(xSin));
  } else if (xSin < 0 && ySin < 0) {
    s800Phi = -TMath::Pi() - TMath::ATan(TMath::Abs(ySin)/TMath::Abs(xSin));
  } else if (xSin > 0 && ySin < 0) {
    s800Phi = -TMath::TwoPi() + TMath::ATan(TMath::Abs(ySin)/xSin);
  } else {
    s800Phi = 0.0;
  }
  s800Theta = TMath::ASin(TMath::Sqrt(xSin*xSin + ySin*ySin));
  
  Double_t gTheta, gPhi;
  gTheta = xyz.Theta();
  if (xyz.Phi() < 0) { gPhi = xyz.Phi() + TMath::TwoPi(); }
  else { gPhi = xyz.Phi(); }
  
  Float_t xxSin = TMath::Sin(fp.track.ata);
  Float_t yySin = TMath::Sin(fp.track.bta);
  
  if (xxSin > 0 && yySin > 0) {
    fp.track.azita = TMath::ATan(yySin/xxSin);
  } else if (xxSin < 0 && yySin > 0) {
    fp.track.azita = (TMath::Pi() - TMath::ATan(yySin/TMath::Abs(xxSin)));
  } else if (xxSin < 0 && yySin < 0) {
    fp.track.azita = (TMath::Pi() + TMath::ATan(TMath::Abs(yySin)/TMath::Abs(xxSin)));
  } else if (xxSin > 0 && yySin < 0) {
    fp.track.azita = (TMath::TwoPi() - TMath::ATan(TMath::Abs(yySin)/xxSin));
  } else {
    fp.track.azita = 0.0;
  }
  fp.track.scatter = TMath::ASin(TMath::Sqrt(xxSin*xxSin + yySin*yySin))*1000;

  Float_t cosDop = (TMath::Sin(gTheta)*TMath::Sin(s800Theta)*
		    (TMath::Sin(gPhi)*TMath::Sin(s800Phi) + TMath::Cos(gPhi)*TMath::Cos(s800Phi)) +
		    TMath::Cos(gTheta)*TMath::Cos(s800Theta));
  Float_t doppler = gamma*(1-beta*cosDop);

  if (fp.track.dta > -1000.) {
    Float_t dp_p = gamma/(1 + gamma) * fp.track.dta;
    beta *= (1 + dp_p/(gamma*gamma));
    gamma = 1/TMath::Sqrt(1. - beta*beta);
  }
  doppler = gamma*(1-beta*cosDop);

  if (fp.track.azita == 0) { doppler = 0.; }
  
  return doppler;
}

/************************************************************/
/* S800 Scaler - Functions                                  */
/************************************************************/

S800Scaler::S800Scaler() {
  start = -65000;
  end = -65000;
  timestamp = -65000;
  globalTS = -65000;
  nScaler = -65000;

  for (UInt_t ui = 0; ui<MAXNSCALER; ui++) {
    value[ui] = -65000;
    valueOld[ui] = -65000;
    rate[ui] = -65000;
  }
}

void S800Scaler::Initialize() {
  start = -65000;
  end = -65000;
  timestamp = -65000;
  globalTS = -65000;
  nScaler = -65000;

  for (UInt_t ui = 0; ui<MAXNSCALER; ui++) {
    value[ui] = -65000;
    valueOld[ui] = -65000;
    rate[ui] = -65000;
  }
}


void S800Scaler::InitializeS800ScalerParameters(const char* inputFilename) {
  FILE *in;
  if ( (in = fopen(inputFilename, "r")) == NULL ) {
    perror("S800 scaler definition file could not be opened.");
  }
  
  cout << "Reading S800 scaler definition file: "
       << inputFilename << endl;
  
  char line[300];
  char name[500];
  int intvalue;
  int i=-1;  

  char* retVal;
  TString temp;

  retVal = fgets(line, 300, in);

  while ( retVal != NULL ) {

    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, "Scaler Mapping", 14) == 0) {
      retVal = fgets(line, 300, in);

      while ( retVal != NULL ) {
	sscanf(line, "%d %s", &intvalue, name);

	temp = name;
	if (intvalue == (i+1)) {
	  Id.push_back(temp);
	  i++;
	} else if (intvalue < i+1) {
	  cout << "Please fix scaler definition file, put entries in numerical order." << endl;
	} else if (intvalue > i+1) {
	  temp = "blank";
	  while (i<intvalue) { Id.push_back(temp); i++; }
	  temp = name;
	  Id.push_back(temp);
	  i++;
	}
	
	retVal = fgets(line, 300, in);
      }
    }
  }

  fclose(in);
}

void S800Scaler::Reset() {
  start = -65000;
  end = -65000;
  timestamp = -65000;
  globalTS = -65000;
  nScaler = -65000;

  for (UInt_t ui = 0; ui<MAXNSCALER; ui++) {
    // value[ui] = -65000;
    valueOld[ui] = -65000;
    rate[ui] = -65000;
  }
}


Bool_t S800Scaler::getAndProcessS800Aux(FILE *inf, Int_t length, long long int ts) {

  UInt_t ui=0;

  Bool_t scaler = 0;

  Reset();
  memset(auxData, 0, sizeof(auxData));
  
  Int_t siz = fread(auxData, 1, length, inf);
  if (siz != length) {
    cout << ALERTTEXT;
    printf("getAndProcessS800() failed in bytes read.\n");
    cout << RESET_COLOR;  fflush(stdout);
  }

  UShort_t *p = auxData;

  ULong64_t upper, lower, size, type;
  
  /* Read the header of the aux data. */
  lower = (ULong64_t)*p;  p++;
  upper = (ULong64_t)*p;  p++;
  size = (ULong64_t)(upper<<16 | lower);
  lower = (ULong64_t)*p;  p++;
  upper = (ULong64_t)*p;  p++;
  type = (ULong64_t)(upper<<16 | lower);

  if (type == PHYSICS_EVENT_COUNT) {

    scaler = 0;

  } else if (type == SCALER_EVENT) {
    
    /* First 8 bytes are start in little-endian form */
    lower = (ULong64_t)*p;  p++;
    upper = (ULong64_t)*p;  p++;
    start = (ULong64_t)(upper<<16 | lower);
    
    /* Next 8 bytes are the end */
    lower = (ULong64_t)*p;  p++;
    upper = (ULong64_t)*p;  p++;
    end = (ULong64_t)(upper<<16 | lower);
    
    /* Next 8 bytes are the timestamp */
    lower = (ULong64_t)*p;  p++;
    upper = (ULong64_t)*p;  p++;
    timestamp = (ULong64_t)(upper<<16 | lower);

    globalTS = ts;

    /* Next 8 bytes are the # of scalers */
    lower = (ULong64_t)*p;  p++;
    upper = (ULong64_t)*p;  p++;
    nScaler = (ULong64_t)(upper<<16 | lower);

    if (nScaler > MAXNSCALER) {
      cerr << "Too many scaler values found.  Something is weird. " << endl;
    }
    
    /* Now scaler values are stored every 8 bytes */
    for (ui = 0; ui<nScaler; ui++) {
      valueOld[ui] = value[ui];
      lower = (ULong64_t)*p;  p++;
      upper = (ULong64_t)*p;  p++;
      value[ui] = (ULong64_t)(upper<<16 | lower);
      rate[ui] = value[ui] - valueOld[ui];
    }

    if (0) {
      cout << "Scaler data: " << endl;
      cout << "   Start   -- " << start << endl;
      cout << "    End    -- " << end << endl;
      cout << "    TS     -- " << timestamp*8 << endl;
      cout << " Global TS -- " << globalTS << endl;
      cout << "  nScaler  -- " << nScaler << endl;
      
      cout << " Values: " << endl;
      for (ui = 0; ui<nScaler; ui++) {
	cout << "    " << ui << " = " << value[ui] << endl; 
      }
    }

    scaler = 1;

  } else {
    //cout << "Unknown auxiliary data type..." << endl;
    scaler = 0;
  }

  return scaler;
}

void S800Physics::Reset() {
  s800Ph.type = 0;
  s800Ph.crdc1_x = 0.0;
  s800Ph.crdc1_y = 0.0;
  s800Ph.crdc2_x = 0.0;
  s800Ph.crdc2_y = 0.0;
  s800Ph.ic_sum = 0.0;
  s800Ph.tof_xfp = 0.0;
  s800Ph.tof_obj = 0.0;
  s800Ph.rf = 0.0;
  s800Ph.trigger = 0;
  s800Ph.ic_de = 0.0;
  s800Ph.tof_xfpe1 = 0.0;
  s800Ph.tof_obje1 = 0.0;
  s800Ph.tof_rfe1 = 0.0;
  s800Ph.ata = 0.0;
  s800Ph.bta = 0.0;
  s800Ph.dta = 0.0;
  s800Ph.yta = 0.0;
}

void S800Full::getPhysics(FILE *inf) {
  size_t siz = fread(&phys.s800Ph, 1, sizeof(s800Phys), inf);
  //printf("GOT IT: %x\n", phys.s800Ph.type);

  fp.crdc1.x = phys.s800Ph.crdc1_x;
  fp.crdc1.y = phys.s800Ph.crdc1_y;
  fp.crdc2.x = phys.s800Ph.crdc2_x;
  fp.crdc2.y = phys.s800Ph.crdc2_y;
  fp.ic.sum = phys.s800Ph.ic_sum;
  tof.xfp = phys.s800Ph.tof_xfp;
  tof.obj = phys.s800Ph.tof_obj;
  fp.track.ata = phys.s800Ph.ata;
  fp.track.bta = phys.s800Ph.bta;
  Double_t xsin = sin(fp.track.ata);
  Double_t ysin = sin(fp.track.bta);
  fp.track.scatter = asin(sqrt( (xsin*xsin) + (ysin*ysin) ))*1000; /* mrad */
}
