/***** CHICO Definitions *****/

#define PPAC_NUM 20
/*#define PPACMATCH 10*/
#define PPACMATCH 5
#define ANGNUM 49
/*#define CH2NS 0.025*/
#define CH2NS 0.100

#define ANODE_E_VSN 1
#define ANODE_T_VSN 2
#define CATHODE_T_VSN 3

#define DATA 0x00000000
#define COUNTMASK 0x3f00
#define COUNTSHIFT 8

#define QDCTYPESHIFT 24
#define QDCGEOSHIFT 27
#define QDCHEADER 2
#define QDCTYPEMASK 0x7000000
#define QDCTRAILER 4
#define QDCGEOMASK 0xf8000000
#define QDCCHANMASK 0x3f0000
#define QDCCHANSHIFT 16
#define QDCDATAMASK 0x0fff

#define TDCTYPEMASK 0xf8000000
#define TDCGEOMASK 0x0000001f
#define TDCCHANMASK 0x03f80000
#define TDCCHANSHIFT 19
#define TDCDATAMASK 0x0007ffff
#define TDCHEADER 0x40000000
#define TDCTRAILER 0x80000000

#define ANODE_REFCH 24
#define CATHODE_REFCH 80
#define RFCH 25





#define beamE 1088.0
//#define beamE 1028.5
#define beamA 136.0
#define targetA 198.0
#define P0 16602.2
//#define P0 16141.8
#define Mtot 334.0
#define BEAM 1
#define TARGET 2
#define LEFT 1
#define RIGHT 2



#define GTPosOffsetZ 5.334
#define GTPosOffsetY1 14.71
#define GTPosOffsetY2 -27.46

