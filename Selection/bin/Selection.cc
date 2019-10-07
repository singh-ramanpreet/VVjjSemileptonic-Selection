#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include <TClonesArray.h>

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

int main (int ac, char** av) {

  std::string inputFile = av[1];

  const   float MUON_MASS = 0.1056583745;
  //lepton cuts
  //const float LEP_PT_VETO_CUT = 20;
  //const float EL_PT_CUT = 35;
  //const float EL_ETA_CUT = 2.5;
  //const float MU_PT_CUT = 35;
  //const float MU_ETA_CUT = 2.4;

  //int run=0, ls=0, evt=0;
  double lep1_pt=0, lep1_eta=0, lep1_phi=0, lep1_m=0;//, lep1_iso=0;
  double lep2_pt=0, lep2_eta=0, lep2_phi=0, lep2_m=0;//, lep2_iso=0;
  //double dilep_m=0;

  //baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  
  TFile *f = TFile::Open(TString(inputFile),"read");
  TTree *t = (TTree *)f->Get("Events");

  //t->SetBranchAddress("Info", &info);            TBranch *infoBr = t->GetBranch("Info");
  t->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = t->GetBranch("Muon");

  std::cout << "nentries: " << t->GetEntries() << std::endl;

  int nSelected=0;

  for (uint i=0; i < t->GetEntries(); i++) {
    //run=0; ls=0; evt=0;
    lep1_pt=0; lep1_eta=0; lep1_phi=0; lep1_m=0; //lep1_iso=0;
    lep2_pt=0; lep2_eta=0; lep2_phi=0; lep2_m=0; //lep2_iso=0;
    //dilep_m=0;
    
    muonArr->Clear();
    muonBr->GetEntry(i);

    for (int j=0; j < muonArr->GetEntries(); j++) {
      const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[j]);

      if ( mu->pt < 35 ) continue;
      if ( abs(mu->eta) > 2.4 ) continue;

      if ( mu->pt > lep1_pt ) {

	lep2_pt=lep1_pt;
	lep2_eta=lep1_eta;
	lep2_phi=lep1_phi;
	lep2_m=lep1_m;
	//lep2_iso=lep1_iso;

	lep1_pt=mu->pt;
	lep1_eta=mu->eta;
	lep1_phi=mu->phi;
	lep1_m=MUON_MASS;
	//lep1_iso = mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), 0.0);

      }
      else if ( mu->pt > lep2_pt ) {

	lep2_pt=mu->pt;
	lep2_eta=mu->eta;
	lep2_phi=mu->phi;
	lep2_m=MUON_MASS;
	//lep2_iso = mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), 0.0);

      }
    }

    if (lep1_pt<0||lep2_pt<0) continue;

    TLorentzVector lep1(0,0,0,0);
    lep1.SetPtEtaPhiM( lep1_pt, lep1_eta, lep1_phi, lep1_m );

    TLorentzVector lep2(0,0,0,0);
    lep2.SetPtEtaPhiM( lep2_pt, lep2_eta, lep2_phi, lep2_m );

    TLorentzVector dilep = lep1+lep2;

    if (dilep.M() < 40 || dilep.M() > 200) continue;
    
    nSelected++;
    
  }

  std::cout << nSelected << std::endl;

  delete f;
  f=0; t=0;

  return 0;

}
