#ifndef __DEFINES_H
#define __DEFINES_H

#define DEBUG 0

#define EB_DIFF_TIME 1000

#define MAXQUADS 40
#define MAXCRYSTALS MAXQUADS*4
#define MAXCHANNELS MAXCRYSTALS*40

#define MAXDETPOS 30
#define MAXCRYSTALNUM 3

#define MAX_TRACE_LENGTH 3000

#define BUFSIZE 2048
#define MAX_SEGS 8              /* max. number of segments to take in events */
#define MAX_INTPTS (2*MAX_SEGS) /* max. number of interaction points */

#define MAX_SIM_GAMMAS 10

/********** GEB Header Types **********/

#define DECOMP 1
#define RAW    0x2
#define RAWHISTORY 0x2B
#define TRACK  3
#define BGS    4
#define S800   5
#define S800AUX 6
#define GRETSCALER 7
#define BANK29 8
#define S800PHYSICS 9
#define S800AUX_TS 10
#define G4SIM 11
#define CHICO 12
#define DGS 14
#define DGSTRIG 15
#define DFMA 16
#define PWALL 17
#define PWALLAUX 18
#define LENDA 21

/********** Superpulse information **********/

/*! Defines the tracelength to be used -- 108 samples, or 1.08us. */
#define AVG_TR_LENGTH 100
/*! Defines the assumed channels to have between waveform starts -- with 
    the definition of AVG_TR_LENGTH 108, this means two bins between
    channels in the .spn files. */
#define AVG_TR_STRIDE 110

#endif
