#ifndef __BGSCONDITIONS_H
#define __BGSCONDITIONS_H

#include "BGSParamList.h"

/* BGS sorting conditions are all defined here... */

/* Energy range for recoils. */
#define r_mine  5000.
#define r_maxe 14000.

/* Conversion electron limits. */
#define pe_mine       100. // parent c.e. E lower limit (keV)
#define pe_maxe      2000. // parent c.e. E upper limit (keV)
#define pe_mint  0.000005  // parent c.e. min lifetime (sec) 
#define pe_maxt  2.000000  // parent c.e. max lifetime (sec)

#define de_mine       100. // daughter c.e. E lower limit (keV)
#define de_maxe      2000. // daughter c.e. E upper limit (keV)
#define de_mint  0.000002  // daughter c.e. min lifetime (sec)
#define de_maxt  0.500000  // daughter c.e. max lifetime (sec)

#define ge_mine       100. // granddaughter c.e. E lower limit (keV)
#define ge_maxe      2000. // granddaughter c.e. E upper limit (keV)
#define ge_mint  0.000002  // granddaughter c.e. min lifetime (sec)
#define ge_maxt  0.100000  // granddaughter c.e. max lifetime (sec)

/* Alpha limits */

#define pa_mine      5000. // parent alpha E lower limit (keV)
#define pa_maxe      9000. // parent alpha E upper limit (keV)
#define pa_mint  0.000010  // parent alpha min lifetime (sec) 
#define pa_maxt   300.000  // parent alpha max lifetime (sec)

#define da_mine      5000. // daughter alpha E lower limit (keV)
#define da_maxe      9000. // daughter alpha E upper limit (keV)
#define da_mint  0.000010  // daughter alpha min lifetime (sec)
#define da_maxt  9000.000  // daughter alpha max lifetime (sec)

#define ga_mine      5000. // granddaughter alpha E lower limit (keV)
#define ga_maxe      9000. // granddaughter alpha E upper limit (keV)
#define ga_mint  0.000010  // granddaughter alpha min lifetime (sec)
#define ga_maxt  0.000009  // granddaughter alpha max lifetime (sec)

/* These are for escapes */

#define px_mine       100. // parent alpha escape E lower limit (keV)
#define px_maxe      7500. // parent alpha escape E upper limit (keV)
#define px_mint   pa_mint  // parent alpha escape min lifetime (sec) 
#define px_maxt   pa_maxt  // parent alpha escape max lifetime (sec)

#define dx_mine       100. // daughter alpha escape E lower limit (keV)
#define dx_maxe      7500. // daughter alpha escape E upper limit (keV)
#define dx_mint   da_mint  // daughter alpha escape min lifetime (sec)
#define dx_maxt   da_maxt  // daughter alpha escape max lifetime (sec)

#define gx_mine       500. // granddaughter alpha escape E lower limit (keV)
#define gx_maxe      7500. // granddaughter alpha escape E upper limit (keV)
#define gx_mint   ga_mint  // granddaughter alpha escape min lifetime (sec)
#define gx_maxt   ga_maxt  // granddaughter alpha escape max lifetime (sec)

/* These are for fissions */

#define f_mine     100000  // minimum energy for fission
#define f_maxe     300000  // maximum energy for fission
#define f_mint   0.000010  // fission min lifetime (sec)
#define f_maxt   25.00000  // fission max lifetime (sec)

/* These are for the MWPC and associated TDCs */

#define mwa_min       50  // MWPC anode minimum channel
#define mwa_max      4000  // MWPC anode minimum channel

#define mwa_fp_tdc_min    50  // MWPC anode - DSSD TDC minimum channel number
#define mwa_fp_tdc_max  4000  // MWPC anode - DSSD TDC maximum channel number
#define mwc_fp_tdc_min     0  // MWPC cathode - DSSD TDC minimum channel number
#define mwc_fp_tdc_max  4096  // MWPC cathode - DSSD TDC maximum channel number

int GeTDC_min[nge] = { 500, 500, 500, 500, 500, 
		       500, 500, 500, 500, 500, 
		       500, 500};

int GeTDC_max[nge] = { 3500, 3500, 3500, 3500, 3500,
		       3500, 3500, 3500, 3500, 3500,
		       3500, 3500};

#define blen 10 // depth of correlation ring buffers. 

#define ADRC 0.994 // correction for calibration with external source

/* Macros for gating conditions go here... */

/* These are used for reconstructing events.  They should
   be equal to or wider than the energy gates used in the 
   correlation searches. */

/* Beta or conversion electron conditions */

#define recon_beta   1    // non-zero allows reconstructing beta/c.e. 
                          // > 1 allows reconstruct with same backstrip
#define beta_mine   50    // minimum total electron energy
#define beta_maxe 2000    // maximum total electron energy
#define beta_match 100    // front and back energies must agree within x (keV)

/* Alpha conditions */

#define recon_alph   27    // non-zero allows reconstructing alphas
                           // > 1 allows reconstruct with same backstrip
#define alph_mine   2000   // minimum total alpha energy
#define alph_maxe   15000  // maximum total alpha energy
#define alph_minore 400    // minimum origin energy for alphas
#define alph_mintee 4000   // minimum terminus energy for alphas
#define alph_minterme 3000 // minimum terminus alpha energy for reconstruction
#define alph_match  500    // front and back energies must agree within x (keV)

/* Fission conditions */

#define recon_fiss      0  // non-zero allows reconstructing fission
                           // > 1 allows reconstruct with same backstrip
#define fiss_mine   50000  // minimum total fission energy
#define fiss_maxe  400000  // maximum total fission energy
#define fiss_minore 40000  // minimum origin energy for fission
#define fiss_match  50000  // front and back energies must agree within x (keV)

/* Punchthru conditions */

#define pttdc_min  500  // punchthru minimum TDC channel
#define pttdc_max 3500  // punchthru maximum TDC channel

/* Upstream TDC conditions */

#define ustdc_min  500  // minimum channel for US TDC
#define ustdc_max 3500  // maximum channel for US TDC

/* FP chip top conditions */

#define fttdc_min    0
#define fttdc_max 4000

/* FP chip west conditions */

#define fwtdc_min    0
#define fwtdc_max 4000

/* FP chip east conditions */

#define fetdc_min    0
#define fetdc_max 4000

/* FP chip back conditions */

#define fbtdc_min    0
#define fbtdc_max 4000

#endif
