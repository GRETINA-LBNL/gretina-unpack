/********************************************************************
* src/LendaDict.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error src/LendaDict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtableLendaDict();
extern void G__cpp_setup_inheritanceLendaDict();
extern void G__cpp_setup_typetableLendaDict();
extern void G__cpp_setup_memvarLendaDict();
extern void G__cpp_setup_globalLendaDict();
extern void G__cpp_setup_memfuncLendaDict();
extern void G__cpp_setup_funcLendaDict();
extern void G__set_cpp_environmentLendaDict();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "src/LENDA-DDAS.h"
#include "src/LENDA-Controls.h"
#include "src/ddasChannel.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__LendaDictLN_TClass;
extern G__linked_taginfo G__LendaDictLN_TBuffer;
extern G__linked_taginfo G__LendaDictLN_TMemberInspector;
extern G__linked_taginfo G__LendaDictLN_TObject;
extern G__linked_taginfo G__LendaDictLN_TNamed;
extern G__linked_taginfo G__LendaDictLN_vectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR;
extern G__linked_taginfo G__LendaDictLN_string;
extern G__linked_taginfo G__LendaDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_maplEstringcOintcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOintgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_maplEstringcOdoublecOlesslEstringgRcOallocatorlEpairlEconstsPstringcOdoublegRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_maplEintcOintcOlesslEintgRcOallocatorlEpairlEconstsPintcOintgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__LendaDictLN_maplEunsignedsPintcOpairlEdoublecOdoublegRcOlesslEunsignedsPintgRcOallocatorlEpairlEconstsPunsignedsPintcOpairlEdoublecOdoublegRsPgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEstringcOallocatorlEstringgRsPgR;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__LendaDictLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TElementActionTlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TElementPosActionTlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTSymlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTRow_constlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTRowlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTColumn_constlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTDiag_constlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTFlat_constlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTSub_constlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTSparseRow_constlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTSparseDiag_constlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTColumnlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTDiaglEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTFlatlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTSublEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTSparseRowlEdoublegR;
extern G__linked_taginfo G__LendaDictLN_TMatrixTSparseDiaglEdoublegR;
extern G__linked_taginfo G__LendaDictLN_ddasChannel;
extern G__linked_taginfo G__LendaDictLN_ddasEvent;
extern G__linked_taginfo G__LendaDictLN_vectorlEddasChannelmUcOallocatorlEddasChannelmUgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEddasChannelmUcOallocatorlEddasChannelmUgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlEddasChannelmUcOallocatorlEddasChannelmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_lendaChannel;
extern G__linked_taginfo G__LendaDictLN_vectorlEintcOallocatorlEintgRsPgR;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_lendaBar;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaChannelcOallocatorlElendaChannelgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaChannelcOallocatorlElendaChannelgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlElendaChannelcOallocatorlElendaChannelgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_lendaEvent;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaBarcOallocatorlElendaBargRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaBarcOallocatorlElendaBargRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlElendaBarcOallocatorlElendaBargRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_lendaFilter;
extern G__linked_taginfo G__LendaDictLN_lendaSettings;
extern G__linked_taginfo G__LendaDictLN_maplEstringcOboolcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOboolgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_maplEintcOstringcOlesslEintgRcOallocatorlEpairlEconstsPintcOstringgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_mapInfo;
extern G__linked_taginfo G__LendaDictLN_lendaPacker;
extern G__linked_taginfo G__LendaDictLN_maplEintcOmapInfocOlesslEintgRcOallocatorlEpairlEconstsPintcOmapInfogRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_maplEstringcOlendaBarcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOlendaBargRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_multimaplEintcOrefTimeContainercOlesslEintgRcOallocatorlEpairlEconstsPintcOrefTimeContainergRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaChannelmUcOallocatorlElendaChannelmUgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaChannelmUcOallocatorlElendaChannelmUgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlElendaChannelmUcOallocatorlElendaChannelmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRcOallocatorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRcOallocatorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRsPgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRcOallocatorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_lesslEvectorlEintcOallocatorlEintgRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_allocatorlEpairlEconstsPvectorlEintcOallocatorlEintgRsPgRcOdoublegRsPgR;
extern G__linked_taginfo G__LendaDictLN_maplEvectorlEintcOallocatorlEintgRsPgRcOdoublecOlesslEvectorlEintcOallocatorlEintgRsPgRsPgRcOallocatorlEpairlEconstsPvectorlEintcOallocatorlEintgRsPgRcOdoublegRsPgRsPgR;
extern G__linked_taginfo G__LendaDictLN_pairlEvectorlEintcOallocatorlEintgRsPgRcOdoublegR;
extern G__linked_taginfo G__LendaDictLN_maplEvectorlEintcOallocatorlEintgRsPgRcOdoublecOlesslEvectorlEintcOallocatorlEintgRsPgRsPgRcOallocatorlEpairlEconstsPvectorlEintcOallocatorlEintgRsPgRcOdoublegRsPgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_maplEvectorlEintcOallocatorlEintgRsPgRcOdoublecOlesslEvectorlEintcOallocatorlEintgRsPgRsPgRcOallocatorlEpairlEconstsPvectorlEintcOallocatorlEintgRsPgRcOdoublegRsPgRsPgRcLcLreverse_iterator;
extern G__linked_taginfo G__LendaDictLN_pairlEmaplEvectorlEintcOallocatorlEintgRsPgRcOdoublecOlesslEvectorlEintcOallocatorlEintgRsPgRsPgRcOallocatorlEpairlEconstsPvectorlEintcOallocatorlEintgRsPgRcOdoublegRsPgRsPgRcLcLiteratorcOboolgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaBarmUcOallocatorlElendaBarmUgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaBarmUcOallocatorlElendaBarmUgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlElendaBarmUcOallocatorlElendaBarmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaSettingsmUcOallocatorlElendaSettingsmUgRsPgR;
extern G__linked_taginfo G__LendaDictLN_vectorlElendaSettingsmUcOallocatorlElendaSettingsmUgRsPgRcLcLiterator;
extern G__linked_taginfo G__LendaDictLN_reverse_iteratorlEvectorlElendaSettingsmUcOallocatorlElendaSettingsmUgRsPgRcLcLiteratorgR;

/* STUB derived class for protected member access */
typedef vector<ddasChannel*,allocator<ddasChannel*> > G__vectorlEddasChannelmUcOallocatorlEddasChannelmUgRsPgR;
typedef vector<lendaChannel,allocator<lendaChannel> > G__vectorlElendaChannelcOallocatorlElendaChannelgRsPgR;
typedef vector<lendaBar,allocator<lendaBar> > G__vectorlElendaBarcOallocatorlElendaBargRsPgR;
typedef vector<lendaChannel*,allocator<lendaChannel*> > G__vectorlElendaChannelmUcOallocatorlElendaChannelmUgRsPgR;
typedef vector<vector<double,allocator<double> >,allocator<vector<double,allocator<double> > > > G__vectorlEvectorlEdoublecOallocatorlEdoublegRsPgRcOallocatorlEvectorlEdoublecOallocatorlEdoublegRsPgRsPgRsPgR;
typedef vector<vector<unsigned short,allocator<unsigned short> >,allocator<vector<unsigned short,allocator<unsigned short> > > > G__vectorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRcOallocatorlEvectorlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgRsPgRsPgR;
typedef map<vector<int,allocator<int> >,double,less<vector<int,allocator<int> > >,allocator<pair<const vector<int,allocator<int> >,double> > > G__maplEvectorlEintcOallocatorlEintgRsPgRcOdoublecOlesslEvectorlEintcOallocatorlEintgRsPgRsPgRcOallocatorlEpairlEconstsPvectorlEintcOallocatorlEintgRsPgRcOdoublegRsPgRsPgR;
typedef vector<lendaBar*,allocator<lendaBar*> > G__vectorlElendaBarmUcOallocatorlElendaBarmUgRsPgR;
typedef vector<lendaSettings*,allocator<lendaSettings*> > G__vectorlElendaSettingsmUcOallocatorlElendaSettingsmUgRsPgR;