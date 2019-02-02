
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class rotationMatrix+;

#pragma link C++ struct globalHeader+;

#pragma link C++ struct mode3DataPacket+;
#pragma link C++ class g3ChannelEvent+;
#pragma link C++ class g3Waveform+;
#pragma link C++ class g3CrystalEvent+;
#pragma link C++ class g3OUT+;
#pragma link C++ class vector<g3ChannelEvent>+;
#pragma link C++ class vector<g3Waveform>+;
#pragma link C++ class vector<g3CrystalEvent>+;
#pragma link C++ struct mode3HistoryPacket+;
#pragma link C++ struct historyEvent+;
#pragma link C++ class g3HistoryEvent+;

#pragma link C++ class Bank88+;

#pragma link C++ struct ipOLD+;
#pragma link C++ struct ip+;
#pragma link C++ struct mode2Old+;
#pragma link C++ struct mode2ABCD1234+;
#pragma link C++ struct mode2ABCD5678+;
#pragma link C++ struct mode2ABCD6789+;
#pragma link C++ class vector<ip>+;
#pragma link C++ class g2IntPt+;
#pragma link C++ class g2CrystalEvent+;
#pragma link C++ class vector<g2IntPt>+;
#pragma link C++ class g2OUT+;
#pragma link C++ class vector<g2CrystalEvent>+;

#pragma link C++ struct trackedGamma+;
#pragma link C++ class  GTrackEvent+;
#pragma link C++ class  vector<trackedGamma>+;
#pragma link C++ class  g1GammaEvent+;
#pragma link C++ class  g1OUT+;
#pragma link C++ class  vector<g1GammaEvent>+;

#pragma link C++ struct g4Sim_emittedGamma+;
#pragma link C++ struct g4Sim_abcd1234+;
#pragma link C++ class g4SimOUT+;

#pragma link C++ class gHistos+;

#pragma link C++ class GRETINAVariables+;

#pragma link C++ class GRETINA+;

#pragma link C++ class controlVariables+;
#pragma link C++ class counterVariables+;

/* SuperPulse.h / SuperPulse.cpp */
#pragma link C++ class SuperPulse+;

/* INLCorrections.h / INLCorrections.cpp */
#pragma link C++ class INLCorrection+;

/* Histos.h / Histos.cpp */
#pragma link C++ class Histos+;

/* Track.h / Track.cpp */
#pragma link C++ struct PM+;

#pragma link C++ struct TrackClusterIP+;
#pragma link C++ struct TrackCluster+;
#pragma link C++ struct TrackShell+;
#pragma link C++ struct TrackCtrl+;
#pragma link C++ struct TrackStat+;
#pragma link C++ struct TrackPerm+;
#pragma link C++ struct Track+;


#endif
