#include "UnpackUtilities.h"

Int_t OpenInputFile(FILE** inf, controlVariables* ctrl, TString runNumber) {
  
  if (ctrl->fileType != "f" && ctrl->fileType != "f1" && ctrl->fileType != "f2")  {

    if (!ctrl->directory.EndsWith("/")) { ctrl->directory = ctrl->directory + "/"; }
      
    if (ctrl->fileType == "g") {
      ctrl->fileName = ctrl->directory + "Run" + runNumber + "/Global.dat";
    } else if (ctrl->fileType == "gr") {
      ctrl->fileName = ctrl->directory + "Run" + runNumber + "/GlobalRaw.dat";
    } else {
      cerr << "WHAT???" << endl;
      return 0;
    }
    
    if (ctrl->compressedFile) {
      
      if (ctrl->noHFC) {
	if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) {
	  ctrl->fileName = ctrl->fileName + ".gz";
	}
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  fclose(*inf);
	  ctrl->fileName = "zcat " + ctrl->fileName;	  
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
	
      } else if (!ctrl->noHFC) {
	
	if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) { 
	  ctrl->fileName = ctrl->fileName + ".gz";
	}
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  fclose(*inf);
	  ctrl->fileName = "./GEB_HFC -z -p " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
      }
      
    } else if (ctrl->compressedFileB) {
      
      if (ctrl->noHFC) {
	if (!ctrl->fileName.EndsWith(".bz2")) { 
	  ctrl->fileName = ctrl->fileName + ".bz2";
	}	 
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  fclose(*inf);
	  ctrl->fileName = "bzcat " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
	
      } else if (!ctrl->noHFC) {
	
	if (!ctrl->fileName.EndsWith(".bz2")) { 
	  ctrl->fileName = ctrl->fileName + ".bz2";
	}
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  ctrl->fileName = "./GEB_HFC -bz -p " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
      }
      
    } else if (ctrl->noHFC) {
      
      *inf = fopen(ctrl->fileName.Data(), "r");
      
    } else {
      
      *inf = fopen(ctrl->fileName.Data(), "r");
      if (*inf) {
	ctrl->fileName = "./GEB_HFC -p " + ctrl->fileName;
	*inf = popen(ctrl->fileName.Data(), "r");
      }
    }
    
  } else if (ctrl->fileType == "f" || ctrl->fileType == "f1" || ctrl->fileType == "f2") {
    
    if (ctrl->compressedFile) {

      if (ctrl->noHFC) {
	if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) {
	  ctrl->fileName = ctrl->fileName + ".gz";
	}
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  fclose(*inf);
	  ctrl->fileName = "zcat " + ctrl->fileName;
	  *inf = NULL;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
	
      } else if (!ctrl->noHFC) {
	
	if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) { 
	  ctrl->fileName = ctrl->fileName + ".gz";
	}
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  fclose(*inf);
	  ctrl->fileName = "./GEB_HFC -z -p " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
      }
      
      
    } else if (ctrl->compressedFileB) {
      if (ctrl->noHFC) {
	if (!ctrl->fileName.EndsWith(".bz2")) { 
	  ctrl->fileName = ctrl->fileName + ".bz2";
	}	 
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  fclose(*inf);
	  ctrl->fileName = "bzcat " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
      } else if (!ctrl->noHFC) {
	if (!ctrl->fileName.EndsWith(".bz2")) { 
	  ctrl->fileName = ctrl->fileName + ".bz2";
	}
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  ctrl->fileName = "./GEB_HFC -bz -p " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
      }
    } else if (ctrl->noHFC) {
      *inf = fopen(ctrl->fileName.Data(), "r");
    } else {
      *inf = fopen(ctrl->fileName.Data(), "r");
      if (*inf) {
	ctrl->fileName = "./GEB_HFC -p " + ctrl->fileName;
	*inf = popen(ctrl->fileName.Data(), "r");
      }      
    }
  }
  
  if (!*inf) {

    printf("Cannot open: %s \n", ctrl->fileName.Data());
    return(2);

  } else {
    
    printf("Opened: %s \n", ctrl->fileName.Data());  
    
    if (ctrl->fileType != "f" && ctrl->fileType != "f1" && ctrl->fileType != "f2") {
      ctrl->outfileName = (ctrl->directory + "Run" + runNumber + "/Run" + runNumber +
			   ctrl->outputSuffix + ".root");
    } else {
      
      if (ctrl->outfileName == "") {
	ctrl->outfileName = ("./ROOTFiles/test.root");
      } else { /* Do nothing, we have a filename. */ } 

    }
    
    return(0);
  }
  
  return(0);
}

int ProcessEvent(Float_t currTS, controlVariables* ctrl, counterVariables* cnt) {

  int badCrystal = 0;

  if (gret->g3Temp.size() > 0) { gret->analyzeMode3(ctrl); }

  /* Write the histograms/tree */
  if (currTS > 0 && badCrystal >= 0) {
    if (ctrl->withTREE) { 
      teb->Fill(); cnt->treeWrites++;
    }
  } 
  
  return (badCrystal);
}

void ResetEvent(controlVariables* ctrl, counterVariables* cnt) {
    
  /* Clear event structures. */
  if ((cnt->getEventBit(RAW)) || (cnt->getEventBit(DECOMP)) ||
      (cnt->getEventBit(BANK88)) || (cnt->getEventBit(RAWHISTORY)) ||
      (cnt->getEventBit(TRACK)) || (cnt->getEventBit(GRETSCALER))) {
    gret->Reset();
  }
  /* Reset temporary crystal event structures. */
  cnt->event = 0x0000;
}



