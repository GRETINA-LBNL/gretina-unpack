{
  
  TString gretinaLib = "lib/libGRETINA.so";
  const char *searchP = "./";
  const char* foundlib;
  foundlib = gSystem->Which(searchP, gretinaLib, EAccessMode::kFileExists);
  if (foundlib) {
    gSystem->Load("lib/libGRETINA.so");
    cout << "Loaded: libGRETINA.so, ";
    gSystem->Load("lib/libGODDESS.so");
    cout << " and libGODDESS.so";
    cout << endl;
    // gSystem->Load("lib/libphosWall.so");
    // cout << " libphosWall.so, ";
    // gSystem->Load("lib/libLenda.so");
    //cout << " libLenda.so, ";
    //gSystem->Load("lib/libS800.so");
    //cout << "and libS800.so" << endl;

  } else {
    gretinaLib = "lib/libGRETINA.dll";
    foundlib = gSystem->Which(searchP, gretinaLib, EAccessMode::kFileExists);
    if (foundlib) {
      gSystem->Load("lib/libGRETINA.dll");
      cout << "Loaded: libGRETINA.dll, ";
      //gSystem->Load("lib/libS800.dll");
      //cout << "and libS800.dll" << endl;
    } else {
      gretinaLib = "lib/libGRETINA.dylib";
      foundlib = gSystem->Which(searchP, gretinaLib, EAccessMode::kFileExists);
      if (foundlib) {
	gSystem->Load("lib/libGRETINA.dylib");
	cout << "Loaded: libGRETINA.dylib, ";
	//gSystem->Load("lib/libS800.dylib");
	//cout << "and libS800.dylib" << endl;
      } else { 
	cout << "Could not load shared libraries!" << endl;
      }
    }
  }

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetNumberContours(99);
}
