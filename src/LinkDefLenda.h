#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ddasChannel+;

#pragma link C++ class ddasEvent+;
#pragma link C++ class std::vector<ddasChannel*>+;

#pragma link C++ class lendaChannel+;

#pragma link C++ class lendaBar+;
#pragma link C++ class std::vector<lendaChannel>+;
#pragma link C++ class std::vector<lendaChannel*>+;

#pragma link C++ class vector<vector<Double_t>>+;
#pragma link C++ class vector<vector<UShort_t>>+;
#pragma link C++ class map<vector<Int_t>, Double_t>+;
#pragma link C++ class lendaEvent+;
#pragma link C++ class std::vector<lendaBar>+;
#pragma link C++ class std::vector<lendaBar*>+;

#pragma link C++ class lendaFilter+;

#pragma link C++ class mapInfo+;
#pragma link C++ class lendaPacker+;

#pragma link C++ class lendaSettings+;
#pragma link C++ class vector<lendaSettings*>+;

#endif
