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

#include "WVJJAna/Selection/interface/WVJJData.hh"

float deltaR(float eta1, float phi1, float eta2, float phi2) {

  float dPhi = fabs(phi1-phi2);
  if (dPhi>6.283185308) dPhi -= 6.283185308;
  if (dPhi>3.141592654) dPhi = 3.141592654 - dPhi;

  float dEta = fabs(eta1-eta2);

  return sqrt( dPhi*dPhi + dEta*dEta );

}

int main (int ac, char** av) {

  std::string inputFile = av[1];

  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;
  const float W_MASS = 80.385;
  const float Z_MASS = 91.1876;

  //lepton cuts
  const float LEP_PT_VETO_CUT = 20;
  const float EL_PT_CUT = 35;
  const float EL_ETA_CUT = 2.5;
  const float MU_PT_CUT = 35;
  const float MU_ETA_CUT = 2.4;

  //ak8 jet cuts
  const float AK8_MIN_PT = 200;
  const float AK8_MAX_ETA = 2.4;
  const float AK8_MIN_SDM = 40;
  const float AK8_MAX_SDM = 150;

  //ak4 jet cuts
  const float AK4_PT_VETO_CUT = 20;
  const float AK4_ETA_CUT = 2.4;
  const float AK4_PT_CUT = 30;
  const float AK4_JJ_MIN_M = 40.0;
  const float AK4_JJ_MAX_M = 150.0;
  const float VBF_MJJ_CUT= 500;

  //cleaning cuts
  const float AK8_DR_CUT = 1.0;
  const float AK4_AK8_DR_CUT = 0.8;
  const float AK4_DR_CUT = 0.3;

  //2016 csv tag thresholds
  const float CSV_LOOSE_2016 = 0.5426;
  const float CSV_MEDIUM_2016 = 0.8484;
  const float CSV_TIGHT_2016 = 0.9535;

  //misc
  const int RUN_BOUND = 278820;
  
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  TClonesArray *vertexArr   = new TClonesArray("baconhep::TVertex");
  TClonesArray *muonArr     = new TClonesArray("baconhep::TMuon");
  TClonesArray *electronArr = new TClonesArray("baconhep::TElectron");
  TClonesArray *AK4Arr      = new TClonesArray("baconhep::TJet");
  TClonesArray *PuppiAK8Arr = new TClonesArray("baconhep::TJet");
  TClonesArray *PuppiAK8AddArr = new TClonesArray("baconhep::TAddJet");

  
  TFile *f = TFile::Open(TString(inputFile),"read");
  TTree *t = (TTree *)f->Get("Events");

  t->SetBranchAddress("Info", &info);            TBranch *infoBr = t->GetBranch("Info");
  t->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = t->GetBranch("Muon");
  t->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = t->GetBranch("Electron");
  t->SetBranchAddress("AK4CHS", &AK4Arr);        TBranch *AK4Br = t->GetBranch("AK4CHS");
  t->SetBranchAddress("AK8Puppi", &PuppiAK8Arr); TBranch *puppiAK8Br = t->GetBranch("AK8Puppi");
  t->SetBranchAddress("AddAK8Puppi", &PuppiAK8AddArr); TBranch *puppiAK8AddBr = t->GetBranch("AddAK8Puppi");
  
  std::cout << "nentries: " << t->GetEntries() << std::endl;

  TFile *of = new TFile("outfile.root","RECREATE");
  TTree *ot = new TTree("Events","Events");
  WVJJData* WVJJTree = new WVJJData(ot);

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> goodAK4Jets;

  //int nSelected=0;

  for (uint i=0; i < t->GetEntries(); i++) {
    infoBr->GetEntry(i);
    WVJJTree->clearVars();

    tightMuon.clear();
    tightEle.clear();

    WVJJTree->run = info->runNum;
    WVJJTree->evt = info->evtNum;
    WVJJTree->ls = info->lumiSec;
    
    // LEPTON SELECTION

    muonArr->Clear();
    muonBr->GetEntry(i);

    int nTightEle=0;
    int nTightMu=0;
    int nVetoLeps=0;

    for (int j=0; j < muonArr->GetEntries(); j++) {
      const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[j]);

      if ( abs(mu->eta) > MU_ETA_CUT ) continue;
      if ( mu->pt < LEP_PT_VETO_CUT ) continue;
      // loose muon ID
      nVetoLeps++;
      // tight muon ID
      if ( mu->pt < MU_PT_CUT ) continue;
      nTightMu++;
      tightMuon.push_back(TLorentzVector(0,0,0,0));
      tightMuon.back().SetPtEtaPhiM(mu->pt, mu->eta, mu->phi, MUON_MASS);

      if ( mu->pt > WVJJTree->lep1_pt ) {

	WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	WVJJTree->lep2_m = WVJJTree->lep1_m;
	WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	WVJJTree->lep2_q = WVJJTree->lep1_q;

	WVJJTree->lep1_pt = mu->pt;
	WVJJTree->lep1_eta = mu->eta;
	WVJJTree->lep1_phi = mu->phi;
	WVJJTree->lep1_m = MUON_MASS;
	WVJJTree->lep1_iso = (mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), 0.0))/(mu->pt);
	WVJJTree->lep1_q = mu->q;

      }
      else if ( mu->pt > WVJJTree->lep2_pt ) {

	WVJJTree->lep2_pt = mu->pt;
	WVJJTree->lep2_eta = mu->eta;
	WVJJTree->lep2_phi = mu->phi;
	WVJJTree->lep2_m = MUON_MASS;
	WVJJTree->lep2_iso = (mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), 0.0))/(mu->pt);
	WVJJTree->lep2_q = mu->q;
      }
    }

    //muon scale variations
    

    //if ( WVJJTree->lep2_pt > 0 && nVetoLeps > 2 ) continue;
    //else if ( WVJJTree->lep1_pt > 0 && nVetoLeps > 1 ) continue;
    
    //if no muons found, check for electrons
    if (WVJJTree->lep1_pt<0) {

      electronArr->Clear();
      electronBr->GetEntry(i);

      for (int j=0; j < electronArr->GetEntries(); j++) {
	const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[j]);

	if ( abs(ele->eta) > EL_ETA_CUT ) continue;
	if ( ele->pt < LEP_PT_VETO_CUT ) continue;
	//loose electron ID
	nVetoLeps++;
	//loose muon ID
	nTightEle++;
	tightEle.push_back(TLorentzVector(0,0,0,0));
	tightEle.back().SetPtEtaPhiM(ele->pt,ele->eta,ele->phi,ELE_MASS);

	if ( ele->pt < EL_PT_CUT ) continue;

	if ( ele->pt > WVJJTree->lep1_pt ) {

	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;

	  WVJJTree->lep1_pt = ele->pt;
	  WVJJTree->lep1_eta = ele->eta;
	  WVJJTree->lep1_phi = ele->phi;
	  WVJJTree->lep1_m = ELE_MASS;
	  //electron ISO
	  WVJJTree->lep1_q = ele->q;

	}
	else if ( ele->pt > WVJJTree->lep2_pt ) {

	  WVJJTree->lep2_pt = ele->pt;
	  WVJJTree->lep2_eta = ele->eta;
	  WVJJTree->lep2_phi = ele->phi;
	  WVJJTree->lep2_m = ELE_MASS;
	  //electron ISO
	  WVJJTree->lep2_q = ele->q;

	}
      }
      //electron scale variations
      WVJJTree->lep1_pt_scaleUp = 1.01 * WVJJTree->lep1_pt;
      WVJJTree->lep1_pt_scaleDn = 0.99 * WVJJTree->lep1_pt;
    }

    //if(!(WWTree->l_pt1>0)) continue;
    //if ((nTightMu+nTightEle)==0) continue; //no leptons with required ID
    //if((nLooseEle+nLooseMu)>2) continue;
    //if(nTightMu>0 && nLooseEle>0) continue;
    //if(nTightEle>0 && nLooseMu>0) continue;
    //if(nTightMu==1 && nLooseMu>1) continue;
    //if(nTightEle==1 && nLooseEle>1) continue;

    //check conditions
    if (WVJJTree->lep1_pt < 0) continue; // no lepton canddiates
    if (nVetoLeps>2) continue; // too many leptons
    if ((nTightMu==1||nTightEle==1) && nVetoLeps>1) continue;
    
    //if ( WVJJTree->lep2_pt > 0 && nVetoLeps > 2 ) continue;
    //else if ( WVJJTree->lep1_pt > 0 && nVetoLeps > 1 ) continue;

    if (WVJJTree->lep1_pt > 0 && WVJJTree->lep2_pt > 0) {
      
      TLorentzVector lep1(0,0,0,0);
      lep1.SetPtEtaPhiM( WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );
      
      TLorentzVector lep2(0,0,0,0);
      lep2.SetPtEtaPhiM( WVJJTree->lep2_pt, WVJJTree->lep2_eta, WVJJTree->lep2_phi, WVJJTree->lep2_m );
      
      TLorentzVector dilep = lep1+lep2;
       
      WVJJTree->dilep_m = dilep.M();
      WVJJTree->dilep_pt = dilep.Pt();
      WVJJTree->dilep_eta = dilep.Eta();
      WVJJTree->dilep_phi = dilep.Phi();

      //dilepton scale variations

    }

    //lepton ID/iso/trigger efficiencies

    // MET

    // AK8

    PuppiAK8Arr->Clear();
    puppiAK8Br->GetEntry(i);
    PuppiAK8AddArr->Clear();
    puppiAK8AddBr->GetEntry(i);

    float dmW = 3000.0;
    int nGoodFatJet=0;
    
    for (int j=0; j<PuppiAK8Arr->GetEntries(); j++) {
      const baconhep::TJet *ak8jet = (baconhep::TJet*)((*PuppiAK8Arr)[j]);
      const baconhep::TAddJet *ak8addjet = (baconhep::TAddJet*)((*PuppiAK8AddArr)[j]);

      if ( ak8jet->pt < AK8_MIN_PT || fabs(ak8jet->eta) > AK8_MAX_ETA ) continue;
      if ( ak8addjet->mass_sd0 < AK8_MIN_SDM || ak8addjet->mass_sd0 > AK8_MAX_SDM ) continue;
      if ( fabs(ak8addjet->mass_sd0 - W_MASS) > dmW ) continue;

      bool isClean=true;
      //lepton cleaning
      for ( std::size_t k=0; k<tightEle.size(); k++) {
	if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		   ak8jet->eta, ak8jet->phi) < AK8_DR_CUT)
	  isClean = false;
      }
      for ( std::size_t k=0; k<tightMuon.size(); k++) {
	if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		   ak8jet->eta, ak8jet->phi) < AK8_DR_CUT)
	  isClean = false;
      }
      if ( isClean == false ) continue;
      
      WVJJTree->bos_PuppiAK8_m_sd0 = ak8addjet->mass_sd0;
      WVJJTree->bos_PuppiAK8_m_sd0_corr = ak8addjet->mass_sd0; //add puppi weight
      WVJJTree->bos_PuppiAK8_tau2tau1 = ak8addjet->tau2/ak8addjet->tau1;
      WVJJTree->bos_PuppiAK8_pt_ungroomed = ak8jet->pt;
      WVJJTree->bos_PuppiAK8_eta_ungroomed = ak8jet->eta;
      WVJJTree->bos_PuppiAK8_phi_ungroomed = ak8jet->phi;
      WVJJTree->bos_PuppiAK8_e_ungroomed = 0; // ak8 jet energy

      dmW = fabs(ak8addjet->mass_sd0 - W_MASS);
      nGoodFatJet++;
    }

    if ( nGoodFatJet == 0 ) {
      
      AK4Arr->Clear();
      AK4Br->GetEntry(i);

      for (int j=0; j<AK4Arr->GetEntries(); j++) {
	const baconhep::TJet *ak4jet = (baconhep::TJet*)((*AK4Arr)[j]);

	//jet energy scale variations
	if ( ak4jet->pt <= AK4_PT_VETO_CUT ) continue;
	if ( fabs(ak4jet->eta) <= AK4_ETA_CUT ) continue;

	bool isClean=true;
	// object cleaning
	for ( std::size_t k=0; k<goodAK4Jets.size(); k++) {
	  if (deltaR(goodAK4Jets.at(k).Eta(), goodAK4Jets.at(k).Phi(),
		     ak4jet->eta, ak4jet->phi) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     ak4jet->eta,   ak4jet->phi) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     ak4jet->eta,   ak4jet->phi) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}

	if ( isClean == false ) continue;

	goodAK4Jets.push_back(TLorentzVector(0,0,0,0));
	goodAK4Jets.back().SetPtEtaPhiM(ak4jet->pt, ak4jet->eta, ak4jet->phi, ak4jet->mass);
	
      }
      


    }
    


    ot->Fill();
    
  }

  of->Write();
  of->Close();

  delete f;
  f=0; t=0;

  return 0;

}
