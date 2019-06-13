#ifndef __TREE_H
#define __TREE_H

/****************************************************/

void InitializeTree() {
  TTree::SetMaxTreeSize(100*Long64_t(2000000000));  
  teb = new TTree("teb", "Tree - event build data");
}

#ifdef WITH_S800
void InitializeTreeS800(controlVariables* ctrl) {

  /* Certain branches we'll always assume we need. */
  /* FP IC calibrated/calculated */
  teb->Branch("fp.ic.de", &(s800->fp.ic.de), "fp.ic.de/D");
  teb->Branch("fp.ic.sum", &(s800->fp.ic.sum), "fp.ic.sum/D");
  /* FP CRDC1 calculated */
  teb->Branch("fp.crdc1.x", &(s800->fp.crdc1.x), "fp.crdc1.x/D");
  teb->Branch("fp.crdc1.y", &(s800->fp.crdc1.y), "fp.crdc1.y/D");
  /* FP CRDC2 calculated */
  teb->Branch("fp.crdc2.x", &(s800->fp.crdc2.x), "fp.crdc2.x/D");
  teb->Branch("fp.crdc2.y", &(s800->fp.crdc2.y), "fp.crdc2.y/D");
  /* S800 Timestamp */
  teb->Branch("ts.timestamp" , &(s800->ts.timestamp), "ts.timestamp/L");
  teb->Branch("s800.evtnum", &(s800->evtnum.eventNum), "s800.evtnum/L");
  /* S800 Trigger */
  teb->Branch("trigger.s800", &(s800->trigger.s800), "trigger.s800/D");
  teb->Branch("trigger.reg", &(s800->trigger.reg), "trigger.reg/D");
  /* S800 Time of Flight */
  teb->Branch("tof.rf", &(s800->tof.rf), "tof.rf/D");
  teb->Branch("tof.obj", &(s800->tof.obj), "tof.obj/D");
  teb->Branch("tof.xfp", &(s800->tof.xfp), "tof.xfp/D");
  teb->Branch("tof.xfp_obj", &(s800->tof.xfp_obj), "tof.xfp_obj/D");
  teb->Branch("tof.rfe1", &(s800->tof.rfe1), "tof.rfe1/D");
  teb->Branch("tof.obje1", &(s800->tof.obje1), "tof.obje1/D");
  teb->Branch("tof.xfpe1", &(s800->tof.xfpe1), "tof.xfpe1/D");
  teb->Branch("tof.mesyrf", &(s800->tof.mesyrf), "tof.mesyrf/D");
  teb->Branch("tof.mesyobj", &(s800->tof.mesyobj), "tof.mesyobj/D");
  teb->Branch("tof.mesyxfp", &(s800->tof.mesyxfp), "tof.mesyxfp/D");
  teb->Branch("tof.mesyrfe1", &(s800->tof.mesyrfe1), "tof.mesyrfe1/D");
  teb->Branch("tof.mesyobje1", &(s800->tof.mesyobje1), "tof.mesyobje1/D");
  teb->Branch("tof.mesyxfpe1", &(s800->tof.mesyxfpe1), "tof.mesyxfpe1/D");
#ifdef S800_LINK_E2
  teb->Branch("tof.obje2", &(s800->tof.obje2), "tof.obje2/D");
  teb->Branch("tof.xfpe2", &(s800->tof.xfpe2), "tof.xfpe2/D");
#endif
#ifdef S800_LINK_TOFTAC
  teb->Branch("tof.tac_obj", &(s800->tof.tac_obj), "tof.tac_obj/D");
  teb->Branch("tof.tac_obje1", &(s800->tof.tac_obje1), "tof.tac_obje1/D");
  teb->Branch("tof.tac_xfp", &(s800->tof.tac_xfp), "tof.tac_xfp/D");
  teb->Branch("tof.tac_xfpe1", &(s800->tof.tac_xfpe1), "tof.tac_xfpe1/D");
#endif
#ifdef S800_LINK_DIAMOND
  teb->Branch("tof.diaor", &(s800->tof.diaor), "tof.diaor/D");
  teb->Branch("tof.dia1", &(s800->tof.dia1), "tof.dia1/D");
  teb->Branch("tof.dia2", &(s800->tof.dia2), "tof.dia2/D");
  teb->Branch("tof.dia3", &(s800->tof.dia3), "tof.dia3/D");
  teb->Branch("tof.dia4", &(s800->tof.dia4), "tof.dia4/D");
  teb->Branch("tof.dia1RF", &(s800->tof.dia1RF), "tof.dia1RF/D");
  teb->Branch("tof.dia2RF", &(s800->tof.dia2RF), "tof.dia2RF/D");
  teb->Branch("tof.dia3RF", &(s800->tof.dia3RF), "tof.dia3RF/D");
  teb->Branch("tof.dia4RF", &(s800->tof.dia4RF), "tof.dia4RF/D");
  teb->Branch("tof.diaRF", &(s800->tof.diaRF), "tof.diaRF/D");
  teb->Branch("tof.dia1Cor", &(s800->tof.dia1Cor), "tof.dia1Cor/D");
  teb->Branch("tof.dia2Cor", &(s800->tof.dia2Cor), "tof.dia2Cor/D");
  teb->Branch("tof.dia3Cor", &(s800->tof.dia3Cor), "tof.dia3Cor/D");
  teb->Branch("tof.dia4Cor", &(s800->tof.dia4Cor), "tof.dia4Cor/D");
  teb->Branch("tof.diaCor", &(s800->tof.diaCor), "tof.diaCor/D");
#endif
  
  /* S800 FP Scintillators */
  if (ctrl->E1_RAW) {
    teb->Branch("fp.e1.de_up", &(s800->fp.e1.de_up), "fp.e1.de_up/D");
    teb->Branch("fp.e1.de_down", &(s800->fp.e1.de_down), "fp.e1.de_down/D");
    teb->Branch("fp.e1.time_up", &(s800->fp.e1.time_up), "fp.e1.time_up/D");
    teb->Branch("fp.e1.time_down", &(s800->fp.e1.time_down), "fp.e1.time_down/D");
  }
  if (ctrl->E1_CAL) {
    teb->Branch("fp.e1.de", &(s800->fp.e1.de), "fp.e1.de/D");
    teb->Branch("fp.e1.time", &(s800->fp.e1.time), "fp.e1.time/D");
    teb->Branch("fp.e1.pos", &(s800->fp.e1.pos), "fp.e1.pos/D");
  }
  
#ifdef S800_LINK_E2
  if (ctrl->E2_RAW) {
    teb->Branch("fp.e2.de_up", &(s800->fp.e2.de_up), "fp.e2.de_up/D");
    teb->Branch("fp.e2.de_down", &(s800->fp.e2.de_down), "fp.e2.de_down/D");
    teb->Branch("fp.e2.time_up", &(s800->fp.e2.time_up), "fp.e2.time_up/D");
    teb->Branch("fp.e2.time_down", &(s800->fp.e2.time_down), "fp.e2.time_down/D");
  }
  if (ctrl->E2_CAL) {
    teb->Branch("fp.e2.de", &(s800->fp.e2.de), "fp.e2.de/D");
    teb->Branch("fp.e2.time", &(s800->fp.e2.time), "fp.e2.time/D");
    teb->Branch("fp.e2.pos", &(s800->fp.e2.pos), "fp.e2.pos/D");
  }
#endif

#ifdef S800_LINK_E3
  if (ctrl->E3_RAW) {
    teb->Branch("fp.e3.de_up", &(s800->fp.e3.de_up), "fp.e3.de_up/D");
    teb->Branch("fp.e3.de_down", &(s800->fp.e3.de_down), "fp.e3.de_down/D");
    teb->Branch("fp.e3.time_up", &(s800->fp.e3.time_up), "fp.e3.time_up/D");
    teb->Branch("fp.e3.time_down", &(s800->fp.e3.time_down), "fp.e3.time_down/D");
  }
  if (ctrl->E3_CAL) {
    teb->Branch("fp.e3.de", &(s800->fp.e3.de), "fp.e3.de/D");
    teb->Branch("fp.e3.time", &(s800->fp.e3.time), "fp.e3.time/D");
    teb->Branch("fp.e3.pos", &(s800->fp.e3.pos), "fp.e3.pos/D");
  }
#endif
  
  /* FP IC raw */
  if (ctrl->IC_RAW) {
    teb->Branch("fp.ic.raw", &(s800->fp.ic.raw), "fp.ic.raw[16]/D");
    teb->Branch("fp.ic.tac1", &(s800->fp.ic.tac1), "fp.ic.tac1/D");
    teb->Branch("fp.ic.tac2", &(s800->fp.ic.tac2), "fp.ic.tac2/D");
    teb->Branch("fp.ic.cal", &(s800->fp.ic.cal), "fp.ic.cal[16]/D");
  }

  /* FP CRDC1 raw pads */
  if (ctrl->CRDC1_RAW_PADS) {
    teb->Branch("fp.crdc1.pad.raw", &(s800->fp.crdc1.pad.raw[0]), "fp.crdc1.pad.raw[224]/D");
    teb->Branch("fp.crdc1.pad.cal", &(s800->fp.crdc1.pad.cal[0]), "fp.crdc1.pad.cal[224]/D");
    teb->Branch("fp.crdc1.pad.samples", &(s800->fp.crdc1.pad.samples[0]),
		"fp.crdc1.pad.samples[224]/I");
    teb->Branch("fp.crdc1.pad.samplemax", &(s800->fp.crdc1.pad.samplemax[0]),
		"fp.crdc1.pad.samplemax[224]/D");
    teb->Branch("fp.crdc1.pad.samplemin", &(s800->fp.crdc1.pad.samplemin[0]),
		"fp.crdc1.pad.samplemin[224]/D");
    teb->Branch("fp.crdc1.pad.delta", &(s800->fp.crdc1.pad.delta[0]),
		"fp.crdc1.pad.delta[224]/D");
    teb->Branch("fp.crdc1.anode", &(s800->fp.crdc1.anode), "fp.crdc1.anode/D");
    teb->Branch("fp.crdc1.evtlength", &(s800->fp.crdc1.evtlength), "fp.crdc1.evtlength/D");
  }
  
  if (ctrl->CRDC1_RAW_CALC) {
    teb->Branch("fp.crdc1.calc.padmax", &(s800->fp.crdc1.calc.padmax),
		"fp.crdc1.calc.padmax/D");
    teb->Branch("fp.crdc1.calc.padsum", &(s800->fp.crdc1.calc.padsum),
		"fp.crdc1.calc.padsum/D");
    teb->Branch("fp.crdc1.calc.padused", &(s800->fp.crdc1.calc.padused),
		"fp.crdc1.calc.padused/D");
    teb->Branch("fp.crdc1.calc.saturationWidth", &(s800->fp.crdc1.calc.saturationWidth),
		"fp.crdc1.calc.saturationWidth/D");
    teb->Branch("fp.crdc1.calc.x_gravity", &(s800->fp.crdc1.calc.x_gravity),
		"fp.crdc1.calc.x_gravity/D");
    teb->Branch("fp.crdc1.tac", &(s800->fp.crdc1.tac), "fp.crdc1.tac/D");
    teb->Branch("fp.crdc1.calc.maxpad", &(s800->fp.crdc1.calc.maxpad[0]),
		"fp.crdc1.calc.maxpad[224]/D");
    teb->Branch("fp.crdc1.calc.maxpad0", &(s800->fp.crdc1.calc.maxpad0[0]),
		"fp.crdc1.calc.maxpad0[224]/D");
    teb->Branch("fp.crdc1.calc.padcalc", &(s800->fp.crdc1.calc.padcalc[0]),
		"fp.crdc1.calc.padcalc[224]/D");
  }
  
  if (ctrl->CRDC2_RAW_PADS) {
    teb->Branch("fp.crdc2.pad.raw", &(s800->fp.crdc2.pad.raw[0]),
		"fp.crdc2.pad.raw[224]/D");
    teb->Branch("fp.crdc2.pad.cal", &(s800->fp.crdc2.pad.cal[0]),
		"fp.crdc2.pad.cal[224]/D");
    teb->Branch("fp.crdc2.pad.samples", &(s800->fp.crdc2.pad.samples[0]),
		"fp.crdc2.pad.samples[224]/I");
    teb->Branch("fp.crdc2.pad.samplemax", &(s800->fp.crdc2.pad.samplemax[0]),
		"fp.crdc2.pad.samplemax[224]/D");
    teb->Branch("fp.crdc2.pad.samplemin", &(s800->fp.crdc2.pad.samplemin[0]),
		"fp.crdc2.pad.samplemin[224]/D");
    teb->Branch("fp.crdc2.pad.delta", &(s800->fp.crdc2.pad.delta[0]),
		"fp.crdc2.pad.delta[224]/D");
    teb->Branch("fp.crdc2.anode", &(s800->fp.crdc2.anode), "fp.crdc2.anode/D");
    teb->Branch("fp.crdc2.evtlength", &(s800->fp.crdc2.evtlength), "fp.crdc2.evtlength/D");
  }

  if (ctrl->CRDC2_RAW_CALC) {
    teb->Branch("fp.crdc2.calc.padmax", &(s800->fp.crdc2.calc.padmax),
		"fp.crdc2.calc.padmax/D");
    teb->Branch("fp.crdc2.calc.padsum", &(s800->fp.crdc2.calc.padsum),
		"fp.crdc2.calc.padsum/D");
    teb->Branch("fp.crdc2.calc.padused", &(s800->fp.crdc2.calc.padused),
		"fp.crdc2.calc.padused/D");
    teb->Branch("fp.crdc2.calc.saturationWidth", &(s800->fp.crdc2.calc.saturationWidth),
	      "fp.crdc2.calc.saturationWidth/D");
    teb->Branch("fp.crdc2.tac", &(s800->fp.crdc2.tac), "fp.crdc2.tac/D");
    teb->Branch("fp.crdc2.calc.x_gravity", &(s800->fp.crdc2.calc.x_gravity),
		"fp.crdc2.calc.x_gravity/D");
    teb->Branch("fp.crdc2.calc.maxpad", &(s800->fp.crdc2.calc.maxpad[0]),
		"fp.crdc2.calc.maxpad[224]/D");
    teb->Branch("fp.crdc2.calc.maxpad0", &(s800->fp.crdc2.calc.maxpad0[0]),
		"fp.crdc2.calc.maxpad0[224]/D");
    teb->Branch("fp.crdc2.calc.padcalc", &(s800->fp.crdc2.calc.padcalc[0]),
		"fp.crdc2.calc.padcalc[224]/D");
  }
  
 /* S800 FP Track */
  //  if (ctrl->FP_TRACK_RAW) {
  if (1) { 
    teb->Branch("fp.track.xfp", &(s800->fp.track.xfp), "fp.track.xfp/D");
    teb->Branch("fp.track.afp", &(s800->fp.track.afp), "fp.track.afp/D");
    teb->Branch("fp.track.yfp", &(s800->fp.track.yfp), "fp.track.yfp/D");
    teb->Branch("fp.track.bfp", &(s800->fp.track.bfp), "fp.track.bfp/D");
    teb->Branch("fp.track.ata", &(s800->fp.track.ata), "fp.track.ata/D");
    teb->Branch("fp.track.yta", &(s800->fp.track.yta), "fp.track.yta/D");
    teb->Branch("fp.track.bta", &(s800->fp.track.bta), "fp.track.bta/D");
    teb->Branch("fp.track.dta", &(s800->fp.track.dta), "fp.track.dta/D");
    teb->Branch("fp.track.azita", &(s800->fp.track.azita), "fp.track.azita/D");
    teb->Branch("fp.track.scatter", &(s800->fp.track.scatter), "fp.track.scatter/D");
    teb->Branch("fp.track.energy", &(s800->fp.track.energy), "fp.track.energy/D");
    teb->Branch("fp.track.deltabeta", &(s800->fp.track.deltabeta), "fp.track.deltabeta/D");
    teb->Branch("fp.track.ptot", &(s800->fp.track.ptot), "fp.track.ptot/D");
    teb->Branch("fp.track.ppar", &(s800->fp.track.ppar), "fp.track.ppar/D");
    teb->Branch("fp.track.ptra", &(s800->fp.track.ptra), "fp.track.ptra/D");
  }

  if (ctrl->FP_TRACK_COR) {
    teb->Branch("fp.track.ata_cor", &(s800->fp.track.ata_cor), "fp.track.ata_cor/D");
    teb->Branch("fp.track.bta_cor", &(s800->fp.track.bta_cor), "fp.track.bta_cor/D");
    teb->Branch("fp.track.azita_cor", &(s800->fp.track.azita_cor), "fp.track.azita_cor/D");
    teb->Branch("fp.track.scatter_cor", &(s800->fp.track.scatter_cor), "fp.track.scatter_cor/D");
  }

  /* S800 Hodoscope */
  if (ctrl->HODO_RAW) {
    teb->Branch("fp.hodo.raw", &(s800->fp.hodo.raw[0]), "fp.hodo.raw[32]/D");
    teb->Branch("fp.hodo.max", &(s800->fp.hodo.max), "fp.hodo.max/D");
    teb->Branch("fp.hodo.sum", &(s800->fp.hodo.sum), "fp.hodo.sum/D");
  }

  if (ctrl->HODO_CAL) {
    teb->Branch("fp.hodo.cal", &(s800->fp.hodo.cal[0]), "fp.hodo.cal[32]/D");
  }

}
#endif /* WITH_S800 */

#ifdef WITH_CHICO
void InitializeTreeCHICO() {
  teb->Branch("particle.id", &(chico->particle.id), "particle.id/I");
  teb->Branch("particle.t", &(chico->particle.t), "particle.t/l");
  teb->Branch("particle.tof", &(chico->particle.tof), "particle.tof/D");
  teb->Branch("particle.mass", &(chico->particle.mass), "particle.mass/D");
  teb->Branch("particle.thetaL", &(chico->particle.thetaL), "particle.thetaL/I");
  teb->Branch("particle.thetaR", &(chico->particle.thetaR), "particle.thetaR/I");
  teb->Branch("particle.phiL", &(chico->particle.phiL), "particle.phiL/I");
  teb->Branch("particle.phiR", &(chico->particle.phiR), "particle.phiR/I");
  teb->Branch("particle.fThetaL", &(chico->particle.fThetaL), "particle.fThetaL/D");
  teb->Branch("particle.fThetaR", &(chico->particle.fThetaR), "particle.fThetaR/D");
  teb->Branch("particle.fPhiL", &(chico->particle.fPhiL), "particle.fPhiL/D");
  teb->Branch("particle.fPhiR", &(chico->particle.fPhiR), "particle.fPhiR/D");
  // teb->Branch("particle.eL", &(chico->particle.eL), "particle.eL/I");
  // teb->Branch("particle.eR", &(chico->particle.eR), "particle.eR/I");
  teb->Branch("particle.rf", &(chico->particle.rf), "particle.rf/D");
  teb->Branch("particle.idParticle", &(chico->idParticle), "particle.idParticle/I");
  teb->Branch("particle.dL", &(chico->particle.dL), "particle.dL/D");
  teb->Branch("particle.dR", &(chico->particle.dR), "particle.dR/D");
  teb->Branch("particle.massP", &(chico->particle.massP), "particle.massP/D");
  teb->Branch("particle.massT", &(chico->particle.massT), "particle.massT/D");
  teb->Branch("particle.betaP", &(chico->particle.betaP), "particle.betaP/D");
  teb->Branch("particle.betaT", &(chico->particle.betaT), "particle.betaT/D");
  teb->Branch("particle.Qvalue", &(chico->particle.Qvalue), "particle.Qvalue/D");
  teb->Branch("particle.pgCosT", &(chico->particle.pgCosT), "particle.pgCosT/D");
  teb->Branch("particle.pgCosP", &(chico->particle.pgCosP), "particle.pgCosP/D");
}
#endif /* WITH_CHICO */

#ifdef WITH_PWALL
void InitializeTreePhosWall() {
  teb->Branch("pwall.ts", &(phosWall->timestamp), "pwall.ts/L");
  teb->Branch("pwall.calc", "phosWallCalc", &(phosWall->calc));  
  teb->Branch("pwall.aux", "phosWallAux", &(phosWall->aux));  
}
#endif

#ifdef WITH_LENDA
void InitializeTreeLENDA() {
  teb->Branch("lenda", "lendaEvent", &lendaEv);
}
#endif

#ifdef WITH_GOD
void InitializeTreeGODDESS() {
  teb->Branch("god.ts", &(goddess->ts), "god.ts/L");
  teb->Branch("godEv", "goddessOut", &(goddess->eventOut));
  teb->Branch("god.eventRaw", "goddessEvent", &(goddess->rawGoddess));
}
#endif

void InitializeTreeMode1() {
  teb->Branch("g1", "g1OUT", &(gret->g1out));
}

void InitializeTreeMode2() {
  teb->Branch("g2", "g2OUT", &(gret->g2out));
}

void InitializeTreeMode3() {
  teb->Branch("g3", "g3OUT", &(gret->g3out));
}

void InitializeTreeSimulation() {
  teb->Branch("gSim", "g4SimOUT", &(gret->gSimOut));
}

void InitializeTreeBank88() {
  teb->Branch("b88", "Bank88", &(gret->b88));
}

void InitializeTreeHistory() {
  teb->Branch("g3H", "g3HistoryEvent", &(gret->g3H));
}

#endif
