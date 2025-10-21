#ifndef __TREE_H
#define __TREE_H

/****************************************************/

void InitializeTree() {
  TTree::SetMaxTreeSize(100*Long64_t(2000000000));  
  teb = new TTree("teb", "Tree - event build data");
}

void InitializeTreeMode3() {
  teb->Branch("g3", "g3OUT", &(gret->g3out));
}

void InitializeTreeBank88() {
  teb->Branch("b88", "Bank88", &(gret->b88));
}

void InitializeTreeHistory() {
  teb->Branch("g3H", "g3HistoryEvent", &(gret->g3H));
}

#endif
