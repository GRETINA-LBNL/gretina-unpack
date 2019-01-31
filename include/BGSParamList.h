#ifndef __BGSPARAMLIST_H
#define __BGSPARAMLIST_H

/* Parameter list for C3 FP Detector */

/* From the CAEN V830 scaler */
#define iusec     1     /* bits 00-15 10 MHz pulse train */
#define iUSEC     2     /* bits 16-31 10 MHz pulse train */
#define ireqt     3     /* bits 00-15 of requested triggers */
#define iREQT     4     /* bits 16-31 of requested triggers */
#define iacct     5     /* bits 00-15 of accepted triggers */
#define iACCT     6     /* bits 16-31 of accepted triggers */
#define istmp     7     /* bits 00-15 of GRETINA timestamp */
#define iSTMP     8     /* bits 16-31 of GRETINA timestamp */

/* From the SIS3806 scaler */
#define iuswh     24    /* bits 00-15 of time since wheel */
#define iUSWH     25    /* bits 16-31 of time since wheel */
#define iusbp     26    /* bits 00-15 of time since pause */
#define iUSBP     27    /* bits 16-31 of time since pause */
#define ibit6     28    /* user bits from SIS3806 scaler */

/* From the CAEN V262 I/O module */
#define iV262     29    /* CAEN V262 I/O module status */
#define ierr1     30    /* readout error location */
#define ierr2     31    /* readout error type */

#define nch       3
#define nfp       32    /* # of signals in FP strip groups */
#define npt       32    /* # of punchthru signals */
#define nus       24    /* # of upstream signals */
#define nge       12    /* # of gamma detector signals */
#define nws       8     /* # of wing detector signals */

/* V785 ADC #0,1  N568 AMP #0,1 */
#define iftl      32    /* index of 0th lo E signal from front of top chip */
#define ifth      64    /* index of 0th hi E signal from front of top chip */

/* V785 ADC #2,3  N568 AMP #2,3 */
#define ifel      96    /* index of 0th lo E signal from front of west chip */
#define ifeh      128   /* index of 0th hi E signal from front of west chip */

/* V785 ADC #4,5  N568 AMP #4,5 */
#define ifwl      160   /* index of 0th lo E signal from front of east chip */
#define ifwh      192   /* index of 0th hi E signal from front of east chip */

/* V785 ADC #6,7  N568 AMP #6,7 */
#define ifbl      224   /* index of 0th lo E signal from back of chips */
#define ifbh      256   /* index of 0th hi E signal from back of chips */

/* V785 ADC #8  N568 AMP #8,9 */
#define iush      288   /* index of 0th hi E upstream signal */
#define iwdh      312   /* index of 0th hi E east wing signal */

/* V785 ADC #9  N568 AMP #10,11 */
#define iptl      320   /* index of 0th lo E punchthru signal */

/* V785 ADC #10 */
#define igee      352   /* index of 0th Ge energy (E) signal */
#define iree      380   /* index of Rutherford east signal */
#define irwe      381   /* index of Rutherford west signal */
#define imae      382   /* index of MWAC anode signal */

/* V785 TDC */
#define iget      400   /* index of 0th Ge TDC signal */

#define iftt      384   /* TDC for front of FP top chip */
#define ifwt      385   /* TDC for front of FP west chip */
#define ifet      386   /* TDC for front of FP east chip */
#define ifbt      387   /* TDC for back of FP detector chips */
#define iust      388   /* TDC for upstream */
#define imat      389   /* TDC for MWAC anode */
#define iptt      390   /* TDC for punchthru detector */

/* These are four-character mnemonics to identify the words in
   the data.  For the first group of 32, CAPITAL LETTERS denote the 
   most significant sixteen bits of the 32-bit scaler, lower case
   is for the least significant 16 bits.  Numbers (XXX:) denote 
   unused words.  Update these to reflect changes to the above
   macros. */

#define taglist								\
  "000:", "usec", "USEC", "reqt", "REQT", "acct", "ACCT", "stmp",	\
  "STMP", "009:", "010:", "011:", "012:", "013:", "014:", "015:",	\
  "016:", "017:", "018:", "019:", "020:", "021:", "022:", "023:",	\
  "uswh", "USWH", "usbp", "USBP", "bit6", "V262", "err1", "err2",	\
									\
  "tl00", "tl01", "tl02", "tl03", "tl04", "tl05", "tl06", "tl06",	\
  "tl08", "tl09", "tl10", "tl11", "tl12", "tl13", "tl14", "tl15",	\
  "tl16", "tl17", "tl18", "tl19", "tl20", "tl21", "tl22", "tl23",	\
  "tl24", "tl25", "tl26", "tl27", "tl28", "tl29", "tl30", "tl31",	\
									\
  "th00", "th01", "th02", "th03", "th04", "th05", "th06", "th06",	\
  "th08", "th09", "th10", "th11", "th12", "th13", "th14", "th15",	\
  "th16", "th17", "th18", "th19", "th20", "th21", "th22", "th23",	\
  "th24", "th25", "th26", "th27", "th28", "th29", "th30", "th31",	\
									\
  "wl00", "wl01", "wl02", "wl03", "wl04", "wl05", "wl06", "wl06",	\
  "wl08", "wl09", "wl10", "wl11", "wl12", "wl13", "wl14", "wl15",	\
  "wl16", "wl17", "wl18", "wl19", "wl20", "wl21", "wl22", "wl23",	\
  "wl24", "wl25", "wl26", "wl27", "wl28", "wl29", "wl30", "wl31",	\
									\
  "wh00", "wh01", "wh02", "wh03", "wh04", "wh05", "wh06", "wh06",	\
  "wh08", "wh09", "wh10", "wh11", "wh12", "wh13", "wh14", "wh15",	\
  "wh16", "wh17", "wh18", "wh19", "wh20", "wh21", "wh22", "wh23",	\
  "wh24", "wh25", "wh26", "wh27", "wh28", "wh29", "wh30", "wh31",	\
									\
  "el00", "el01", "el02", "el03", "el04", "el05", "el06", "el06",	\
  "el08", "el09", "el10", "el11", "el12", "el13", "el14", "el15",	\
  "el16", "el17", "el18", "el19", "el20", "el21", "el22", "el23",	\
  "el24", "el25", "el26", "el27", "el28", "el29", "el30", "el31",	\
  									\
  "eh00", "eh01", "eh02", "eh03", "eh04", "eh05", "eh06", "eh06",	\
  "eh08", "eh09", "eh10", "eh11", "eh12", "eh13", "eh14", "eh15",	\
  "eh16", "eh17", "eh18", "eh19", "eh20", "eh21", "eh22", "eh23",	\
  "eh24", "eh25", "eh26", "eh27", "eh28", "eh29", "eh30", "eh31",	\
  									\
  "bl00", "bl01", "bl02", "bl03", "bl04", "bl05", "bl06", "bl06",	\
  "bl08", "bl09", "bl10", "bl11", "bl12", "bl13", "bl14", "bl15",	\
  "bl16", "bl17", "bl18", "bl19", "bl20", "bl21", "bl22", "bl23",	\
  "bl24", "bl25", "bl26", "bl27", "bl28", "bl29", "bl30", "bl31",	\
									\
  "bh00", "bh01", "bh02", "bh03", "bh04", "bh05", "bh06", "bh06",	\
  "bh08", "bh09", "bh10", "bh11", "bh12", "bh13", "bh14", "bh15",	\
  "bh16", "bh17", "bh18", "bh19", "bh20", "bh21", "bh22", "bh23",	\
  "bh24", "bh25", "bh26", "bh27", "bh28", "bh29", "bh30", "bh31",	\
  									\
  "uh00", "uh01", "uh02", "uh03", "uh04", "uh05", "uh06", "uh06",	\
  "uh08", "uh09", "uh10", "uh11", "uh12", "uh13", "uh14", "uh15",	\
  "uh16", "uh17", "uh18", "uh19", "uh20", "uh21", "uh22", "uh23",	\
  "EW00", "EW01", "EW02", "EW03", "WW00", "WW01", "WW02", "WW03",	\
									\
  "pt00", "pt01", "pt02", "pt03", "pt04", "pt05", "pt06", "pt06",	\
  "pt08", "pt09", "pt10", "pt11", "pt12", "pt13", "pt14", "pt15",	\
  "pt16", "pt17", "pt18", "pt19", "pt20", "pt21", "pt22", "pt23",	\
  "pt24", "pt25", "pt26", "pt27", "pt28", "pt29", "pt30", "pt31",	\
									\
  "ge00", "ge01", "ge02", "ge03", "ge04", "ge05", "ge06", "ge07",	\
  "ge08", "ge09", "ge10", "ge11", "364:", "365:", "366:", "367:",	\
  "368:", "369:", "370:", "371:", "372:", "373:", "374:", "375:",	\
  "376:", "377:", "378:", "379:", "rute", "rutw", "mwae", "383:",	\
  									\
  "fptt", "fpwt", "fpet", "fpbt", "_ust", "mwat", "390:", "391:",	\
  "392:", "393:", "394:", "395:", "396:", "397:", "398:", "399:",	\
  "gt00", "gt01", "gt02", "gt03", "gt04", "gt05", "gt06", "gt07",	\
  "gt08", "gt09", "gt10", "gt11", "412:", "413:", "414:", "415:",	\
  									\
  "rloe", "rhie", "ploe", "phie", "dloe", "dhie", "floe", "fhie",	\
  "radt", "rfdt", "tf1l", "tf1h", "tf2l", "tf2h", "430:", "431:",	\
  "pl00", "pl01", "pl02", "pl03", "pl04", "pl05", "pl06", "pl07",	\
  "pl08", "pl09", "pl10", "pl11", "pl12", "pl13", "pl14", "pl15",	\
  "pl16", "pl17", "pl18", "pl19", "pl20", "pl21", "pl22", "pl23",	\
  "pl24", "pl25", "pl26", "pl27", "pl28", "pl29", "pl30", "pl31",	\
  "shut", "upst", "anti", "mwpc", "fisc", "mfpe", "mupe", "bsiz",	\
  "rell", "reul", "rwll", "rwul", "tsig", "psig", "pmul", "adrc"

#endif
