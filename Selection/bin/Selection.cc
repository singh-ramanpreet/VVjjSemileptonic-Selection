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
#include "TH2.h"
#include <TClonesArray.h>

#include "WVJJAna/Selection/interface/NanoAOD_MC.hh"
#include "WVJJAna/Selection/interface/NanoAOD_Weights.hh"
#include "WVJJAna/Selection/interface/Utils.hh"
#include "WVJJAna/Selection/interface/METzCalculator.h"

int main (int ac, char** av) {

  std::string inputFile = av[1];
  std::string outputFile = av[2];
  int isMC = atoi(av[3]);
  int era = atoi(av[4]);
  int nanoVersion = atoi(av[5]);

  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;
  const float W_MASS = 80.385;
  //const float Z_MASS = 91.1876;

  //lepton cuts
  const float LEP_PT_VETO_CUT = 10;
  const float EL_PT_CUT = 20;
  const float EL_ETA_CUT = 2.5;
  const float MU_PT_CUT = 20;
  const float MU_ETA_CUT = 2.4;

  //ak8 jet cuts
  const float AK8_MIN_PT = 200;
  const float AK8_MAX_ETA = 2.4;
  const float AK8_MIN_SDM = 40;
  const float AK8_MAX_SDM = 150;

  //ak4 jet cuts
  //const float AK4_PT_VETO_CUT = 20;
  const float AK4_ETA_CUT = 2.4;
  const float AK4_PT_CUT = 30;
  const float AK4_JJ_MIN_M = 40.0;
  const float AK4_JJ_MAX_M = 150.0;
  const float VBF_MJJ_CUT= 500;

  //cleaning cuts
  const float AK8_LEP_DR_CUT = 1.0;
  const float AK4_AK8_DR_CUT = 0.8;
  const float AK4_DR_CUT = 0.3;

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<int>goodJetIndex;

  //
  //
  //   INPUT/OUTPUT
  //
  //

  TFile *of = new TFile(outputFile.c_str(),"RECREATE");
  TTree *ot = new TTree("Events","Events");
  WVJJData* WVJJTree = new WVJJData(ot);  
  TH1F *totalEvents = new TH1F("TotalEvents","TotalEvents",2,-1,1);

  TFile *f=0;
  TTree *t=0, *r=0;

  std::ifstream ifs;
  ifs.open(inputFile.data());
  assert(ifs.is_open());
  std::string line;
  while (getline(ifs,line)) {
    std::stringstream ss(line);
    std::string filetoopen;
    ss >> filetoopen;
    
    f = TFile::Open(TString("root://cmseos.fnal.gov/")+TString(filetoopen),"read");
    //f = TFile::Open(TString("root://xrootd-cms.infn.it/")+TString(filetoopen),"read");
    t = (TTree *)f->Get("Events");
    r = (TTree *)f->Get("Runs");
    if (t==NULL) continue;
    
    std::cout << filetoopen << std::endl;
    
    NanoAOD_MC NanoReader = NanoAOD_MC(t);
    NanoAOD_Weights NanoWeightReader = NanoAOD_Weights(r);

    if (isMC==1) {

      for (uint i=0; i < r->GetEntries(); i++) {
        NanoWeightReader.GetEntry(i);
        totalEvents->SetBinContent(2,totalEvents->GetBinContent(2)+NanoWeightReader.genEventSumw_);
      }

    }

    for (uint i=0; i < t->GetEntries(); i++) {
      WVJJTree->clearVars();
      NanoReader.GetEntry(i);

      if (i%10000==0) std::cout <<"event " << i << std::endl;

      if (isMC==1) {
	WVJJTree->genWeight=NanoReader.Generator_weight;
      }
      
      if (era==2018) {

      	if ( NanoReader.HLT_IsoMu24 || NanoReader.HLT_IsoMu27 || NanoReader.HLT_IsoMu30 || NanoReader.HLT_Mu50 ) WVJJTree->trigger_1Mu = true;

	if ( NanoReader.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || NanoReader.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ||
	     NanoReader.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || NanoReader.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 ) WVJJTree->trigger_2Mu = true;

	if ( NanoReader.HLT_Ele27_WPTight_Gsf || NanoReader.HLT_Ele28_WPTight_Gsf || NanoReader.HLT_Ele32_WPTight_Gsf ||
	     NanoReader.HLT_Ele35_WPTight_Gsf || NanoReader.HLT_Ele38_WPTight_Gsf || NanoReader.HLT_Ele40_WPTight_Gsf ) WVJJTree->trigger_1El = true;
	
	if ( NanoReader.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || NanoReader.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ||
	     NanoReader.HLT_DiEle27_WPTightCaloOnly_L1DoubleEG || NanoReader.HLT_DoubleEle33_CaloIdL_MW ||
	     NanoReader.HLT_DoubleEle25_CaloIdL_MW || NanoReader.HLT_DoubleEle27_CaloIdL_MW ) WVJJTree->trigger_2El = true;
      	
      }
      else if (era==2017) {

      	if ( NanoReader.HLT_IsoMu24 || NanoReader.HLT_IsoMu27 || NanoReader.HLT_IsoMu30 || NanoReader.HLT_Mu50 ) WVJJTree->trigger_1Mu = true;


	if ( NanoReader.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || NanoReader.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ||
	     NanoReader.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || NanoReader.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 ) WVJJTree->trigger_2Mu = true;


	if ( NanoReader.HLT_Ele27_WPTight_Gsf || NanoReader.HLT_Ele32_WPTight_Gsf || NanoReader.HLT_Ele35_WPTight_Gsf ) WVJJTree->trigger_1El = true;

	if ( NanoReader.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || NanoReader.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ||
	     NanoReader.HLT_DiEle27_WPTightCaloOnly_L1DoubleEG || NanoReader.HLT_DoubleEle33_CaloIdL_MW || 
	     NanoReader.HLT_DoubleEle25_CaloIdL_MW ) WVJJTree->trigger_2El = true;
	
      }

      //std::cout << "passed trigger: ";
      //if (WVJJTree->trigger_1Mu) std::cout << "1 muon ";
      //if (WVJJTree->trigger_2Mu) std::cout << "2 muon ";
      //if (WVJJTree->trigger_1El) std::cout << "1 ele ";
      //if (WVJJTree->trigger_2El) std::cout << "2 ele ";
      //std::cout << std::endl;
      
      if ( ! ( WVJJTree->trigger_1Mu || WVJJTree->trigger_2Mu || WVJJTree->trigger_1El || WVJJTree->trigger_2El ) ) continue;



      tightMuon.clear();
      tightEle.clear();
      
      WVJJTree->run = NanoReader.run;
      WVJJTree->evt = NanoReader.event;
      WVJJTree->ls = NanoReader.luminosityBlock;
      
      WVJJTree->nPV = NanoReader.PV_npvsGood;
      WVJJTree->nPU_mean = NanoReader.Pileup_nPU;
      
      WVJJTree->puWeight = 1.0;//scaleFactor.GetPUWeight(info->nPUmean, 0);
      
      // LEPTON SELECTION
      
      int nTightEle=0;
      int nTightMu=0;
      int nVetoEle=0;
      int nVetoMu=0;
      
      for (uint j=0; j < NanoReader.nMuon; j++) {
	if ( abs(NanoReader.Muon_eta[j]) > MU_ETA_CUT ) continue;
	//using conservative uncertainty value of 3%
	if ( 1.03*NanoReader.Muon_pt[j] < LEP_PT_VETO_CUT ) continue;
	
	if (!NanoReader.Muon_looseId[j]) continue;
	if (NanoReader.Muon_pfRelIso04_all[j]>0.25) continue;
	nVetoMu++;

	//using conservative uncertainty value of 3%
	if ( 1.03*NanoReader.Muon_pt[j] < MU_PT_CUT ) continue;
	if (!NanoReader.Muon_tightId[j]) continue;
	if (NanoReader.Muon_pfRelIso04_all[j]>0.15) continue;

	nTightMu++;
	tightMuon.push_back(TLorentzVector(0,0,0,0));
	tightMuon.back().SetPtEtaPhiM(NanoReader.Muon_pt[j], NanoReader.Muon_eta[j], 
				      NanoReader.Muon_phi[j], MUON_MASS);
	
	if ( NanoReader.Muon_pt[j] > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = NanoReader.Muon_pt[j];
	  WVJJTree->lep1_eta = NanoReader.Muon_eta[j];
	  WVJJTree->lep1_phi = NanoReader.Muon_phi[j];
	  WVJJTree->lep1_m = MUON_MASS;
	  WVJJTree->lep1_iso = NanoReader.Muon_pfRelIso04_all[j];
	  WVJJTree->lep1_q = NanoReader.Muon_charge[j];
	  
	}
	else if ( NanoReader.Muon_pt[j] > WVJJTree->lep2_pt ) {
	  //if (WVJJTree->lep1_q*NanoReader.Muon_charge[j]>0) continue;

	  WVJJTree->lep2_pt = NanoReader.Muon_pt[j];
	  WVJJTree->lep2_eta = NanoReader.Muon_eta[j];
	  WVJJTree->lep2_phi = NanoReader.Muon_phi[j];
	  WVJJTree->lep2_m = MUON_MASS;
	  WVJJTree->lep2_iso = NanoReader.Muon_pfRelIso04_all[j];
	  WVJJTree->lep2_q = NanoReader.Muon_charge[j];
	  
	  }
      }
      
      for (uint j=0; j < NanoReader.nElectron; j++) {
	if ( abs(NanoReader.Electron_eta[j]) > EL_ETA_CUT ) continue;
	//using conservative uncertainty value of 3%
	if ( 1.03*NanoReader.Electron_pt[j] < LEP_PT_VETO_CUT ) continue;
	
	//cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
	if ( NanoReader.Electron_cutBased[j]<2 ) continue;
	nVetoEle++;

	//using conservative uncertainty value of 3%
	if ( 1.03*NanoReader.Electron_pt[j] < EL_PT_CUT ) continue;
	
	if (NanoReader.Electron_cutBased[j]<4) continue;
	nTightEle++;
	
	tightEle.push_back(TLorentzVector(0,0,0,0));
	tightEle.back().SetPtEtaPhiM(NanoReader.Electron_pt[j],NanoReader.Electron_eta[j],
				     NanoReader.Electron_phi[j],ELE_MASS);
	
	//don't try to select electrons unless we don't already
	//have muons
	if (WVJJTree->lep1_m == MUON_MASS) continue;
	
	if ( NanoReader.Electron_pt[j] > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = NanoReader.Electron_pt[j];
	  WVJJTree->lep1_eta = NanoReader.Electron_eta[j];
	  WVJJTree->lep1_phi = NanoReader.Electron_phi[j];
	  WVJJTree->lep1_m = ELE_MASS;
	  WVJJTree->lep1_iso = NanoReader.Electron_pfRelIso03_all[j];
	  WVJJTree->lep1_q = NanoReader.Electron_charge[j];
	  
	}
	else if ( NanoReader.Electron_pt[j] > WVJJTree->lep2_pt ) {

	  WVJJTree->lep2_pt = NanoReader.Electron_pt[j];
	  WVJJTree->lep2_eta = NanoReader.Electron_eta[j];
	  WVJJTree->lep2_phi = NanoReader.Electron_phi[j];
	  WVJJTree->lep2_m = ELE_MASS;
	  WVJJTree->lep2_iso = NanoReader.Electron_pfRelIso03_all[j];
	  WVJJTree->lep2_q = NanoReader.Electron_charge[j];
	  
	}
      }

      //check conditions
      if(!(WVJJTree->lep1_pt>0)) continue;
      if ((nTightMu+nTightEle)==0) continue; //no leptons with required ID
      if((nVetoEle+nVetoMu)>2) continue;
      if(nTightMu>0 && nVetoEle>0) continue;
      if(nTightEle>0 && nVetoMu>0) continue;
      if(nTightMu==1 && nVetoMu>1) continue;
      if(nTightEle==1 && nVetoEle>1) continue;

      //muon scale variations
      if (WVJJTree->lep1_m == MUON_MASS) {
	//https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceScaleResolRun2
	if (WVJJTree->lep1_eta<-2.1) {
	  WVJJTree->lep1_pt_scaleUp = 1.027 * WVJJTree->lep1_pt;
	  WVJJTree->lep1_pt_scaleUp = 0.973 * WVJJTree->lep1_pt;
	}
	else if (WVJJTree->lep1_eta<-1.2) {
	  WVJJTree->lep1_pt_scaleUp = 1.009 * WVJJTree->lep1_pt;
          WVJJTree->lep1_pt_scaleUp = 0.991 * WVJJTree->lep1_pt;
	}
	else if (WVJJTree->lep1_eta<1.2) {
	  WVJJTree->lep1_pt_scaleUp = 1.004 * WVJJTree->lep1_pt;
          WVJJTree->lep1_pt_scaleUp = 0.996 * WVJJTree->lep1_pt;
	}
	else if (WVJJTree->lep1_eta<2.1) {
	  WVJJTree->lep1_pt_scaleUp = 1.009 * WVJJTree->lep1_pt;
          WVJJTree->lep1_pt_scaleUp = 0.991 * WVJJTree->lep1_pt;
	}
	else {
	  WVJJTree->lep1_pt_scaleUp = 1.017 * WVJJTree->lep1_pt;
          WVJJTree->lep1_pt_scaleUp = 0.983 * WVJJTree->lep1_pt;
	}
	if (WVJJTree->lep2_pt>0) {
	  if (WVJJTree->lep2_eta<-2.1) {
	    WVJJTree->lep2_pt_scaleUp = 1.027 * WVJJTree->lep2_pt;
	    WVJJTree->lep2_pt_scaleUp = 0.973 * WVJJTree->lep2_pt;
	  }
	  else if (WVJJTree->lep2_eta<-1.2) {
	    WVJJTree->lep2_pt_scaleUp = 1.009 * WVJJTree->lep2_pt;
	    WVJJTree->lep2_pt_scaleUp = 0.991 * WVJJTree->lep2_pt;
	  }
	  else if (WVJJTree->lep2_eta<1.2) {
	    WVJJTree->lep2_pt_scaleUp = 1.004 * WVJJTree->lep2_pt;
	    WVJJTree->lep2_pt_scaleUp = 0.996 * WVJJTree->lep2_pt;
	  }
	  else if (WVJJTree->lep2_eta<2.1) {
	    WVJJTree->lep2_pt_scaleUp = 1.009 * WVJJTree->lep2_pt;
	    WVJJTree->lep2_pt_scaleUp = 0.991 * WVJJTree->lep2_pt;
	  }
	  else {
	    WVJJTree->lep2_pt_scaleUp = 1.017 * WVJJTree->lep2_pt;
	    WVJJTree->lep2_pt_scaleUp = 0.983 * WVJJTree->lep2_pt;
	  }
	}
      }
      
      //electron scale variations
      if (WVJJTree->lep1_m == ELE_MASS) {
	WVJJTree->lep1_pt_scaleUp = 1.01 * WVJJTree->lep1_pt;
	WVJJTree->lep1_pt_scaleDn = 0.99 * WVJJTree->lep1_pt;
	if (WVJJTree->lep2_pt>0) {
	  WVJJTree->lep2_pt_scaleUp = 1.01 * WVJJTree->lep2_pt;
	  WVJJTree->lep2_pt_scaleDn = 0.99 * WVJJTree->lep2_pt;
	}
      }
      
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

	lep1.SetPtEtaPhiM( WVJJTree->lep1_pt_scaleUp, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );
	lep2.SetPtEtaPhiM( WVJJTree->lep2_pt_scaleUp, WVJJTree->lep2_eta, WVJJTree->lep2_phi, WVJJTree->lep2_m );

	WVJJTree->dilep_m_scaleUp = (lep1+lep2).M();
	WVJJTree->dilep_pt_scaleUp = (lep1+lep2).Pt();

	lep1.SetPtEtaPhiM( WVJJTree->lep1_pt_scaleDn, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );
	lep2.SetPtEtaPhiM( WVJJTree->lep2_pt_scaleDn, WVJJTree->lep2_eta, WVJJTree->lep2_phi, WVJJTree->lep2_m );

	WVJJTree->dilep_m_scaleDn = (lep1+lep2).M();
	WVJJTree->dilep_pt_scaleDn = (lep1+lep2).Pt();
	
      }
      
      //lepton ID/iso/trigger efficiencies
      //WVJJTree->lep1_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
      //WVJJTree->lep1_idEffWeight = GetSFs_Lepton(WVJJTree->lep1_pt, WVJJTree->lep1_eta, hIDIsoEle);
      //WVJJTree->lep1_idEffWeight *= GetSFs_Lepton(WVJJTree->lep1_pt,WVJJTree->lep1_eta, hGSFCorrEle);
      //WVJJTree->lep1_idEffWeight *= GetSFs_Lepton(WVJJTree->lep1_pt,WVJJTree->lep1_eta, hTriggerEle);
      
      //if (WVJJTree->lep2_pt>0) {
	//WVJJTree->lep2_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep2_pt, WVJJTree->lep2_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
	//WVJJTree->lep2_idEffWeight = GetSFs_Lepton(WVJJTree->lep2_pt, WVJJTree->lep2_eta, hIDIsoEle);
	//WVJJTree->lep2_idEffWeight *= GetSFs_Lepton(WVJJTree->lep2_pt,WVJJTree->lep2_eta, hGSFCorrEle);
	//WVJJTree->lep2_idEffWeight *= GetSFs_Lepton(WVJJTree->lep2_pt,WVJJTree->lep2_eta, hTriggerEle);
      //}
      
      
      // MET
      
      if (NanoReader.MET_pt < 0) continue;
      WVJJTree->MET = NanoReader.MET_pt;
      WVJJTree->MET_phi = NanoReader.MET_phi;
      WVJJTree->MET_scaleUp = NanoReader.MET_pt_jesTotalUp;
      WVJJTree->MET_scaleDn = NanoReader.MET_pt_jesTotalDown;
      
      if (WVJJTree->lep2_pt<0) {
	
	TLorentzVector lep1(0,0,0,0);
	lep1.SetPtEtaPhiM( WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );

	TLorentzVector tempMet(0,0,0,0);
	tempMet.SetPxPyPzE(WVJJTree->MET*TMath::Cos(WVJJTree->MET_phi), 
			   WVJJTree->MET*TMath::Sin(WVJJTree->MET_phi), 
			   0.0, WVJJTree->MET);
	
	METzCalculator NeutrinoPz_type0;
	NeutrinoPz_type0.SetMET(tempMet);
	NeutrinoPz_type0.SetLepton(lep1);
	NeutrinoPz_type0.SetLeptonType(WVJJTree->lep1_m == ELE_MASS ? "el" : "mu");
	
	WVJJTree->neu_pz_type0 = NeutrinoPz_type0.Calculate();

	tempMet.SetPxPyPzE(WVJJTree->MET_scaleUp*TMath::Cos(WVJJTree->MET_phi), 
			   WVJJTree->MET_scaleUp*TMath::Sin(WVJJTree->MET_phi), 
			   0.0, WVJJTree->MET_scaleUp);
	
	NeutrinoPz_type0.SetMET(tempMet);
	WVJJTree->neu_pz_type0_scaleUp = NeutrinoPz_type0.Calculate();

	tempMet.SetPxPyPzE(WVJJTree->MET_scaleDn*TMath::Cos(WVJJTree->MET_phi), 
			   WVJJTree->MET_scaleDn*TMath::Sin(WVJJTree->MET_phi), 
			   0.0, WVJJTree->MET_scaleDn);	

	NeutrinoPz_type0.SetMET(tempMet);
	WVJJTree->neu_pz_type0_scaleDn = NeutrinoPz_type0.Calculate();

	tempMet.SetPxPyPzE(WVJJTree->MET*TMath::Cos(WVJJTree->MET_phi), 
			   WVJJTree->MET*TMath::Sin(WVJJTree->MET_phi), 
			   0.0, WVJJTree->MET);
	
	TLorentzVector neutrino(0,0,0,0);
	neutrino.SetPxPyPzE(tempMet.Px(), tempMet.Py(), WVJJTree->neu_pz_type0, 
			    sqrt(tempMet.Px()*tempMet.Px()+tempMet.Py()*tempMet.Py()+ WVJJTree->neu_pz_type0* WVJJTree->neu_pz_type0));
	
	TLorentzVector bosLep = lep1+neutrino;
	
	WVJJTree->dilep_m = bosLep.M();
	WVJJTree->dilep_pt = bosLep.Pt();
	WVJJTree->dilep_eta = bosLep.Eta();
	WVJJTree->dilep_phi = bosLep.Phi();
	
	//dilepton scale variations
	
      }
      
      // AK8
      
      float dmW = 3000.0;
      int nGoodFatJet=0;

      for (uint j=0; j<NanoReader.nFatJet; j++) {
	if ( ! (NanoReader.FatJet_pt[j]>AK8_MIN_PT || NanoReader.FatJet_pt_jesTotalUp[j]>AK8_MIN_PT || 
		NanoReader.FatJet_pt_jesTotalDown[j]>AK8_MIN_PT) ) continue;
	if ( fabs(NanoReader.FatJet_eta[j]) > AK8_MAX_ETA ) continue;
	
	if ( ! (NanoReader.FatJet_msoftdrop[j]>AK8_MIN_SDM || NanoReader.FatJet_msoftdrop_jesTotalUp[j]>AK8_MIN_SDM ||
		NanoReader.FatJet_msoftdrop_jesTotalDown[j]>AK8_MIN_SDM) ) continue;

	if ( ! (NanoReader.FatJet_msoftdrop[j]<AK8_MAX_SDM || NanoReader.FatJet_msoftdrop_jesTotalUp[j]<AK8_MAX_SDM ||
		NanoReader.FatJet_msoftdrop_jesTotalDown[j]<AK8_MAX_SDM) ) continue;
	
	bool isClean=true;
	//lepton cleaning
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     NanoReader.FatJet_eta[j], NanoReader.FatJet_phi[j]) < AK8_LEP_DR_CUT)
	    isClean = false;
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     NanoReader.FatJet_eta[j], NanoReader.FatJet_phi[j]) < AK8_LEP_DR_CUT)
	    isClean = false;
	}

	if ( isClean == false ) continue;

	if ( fabs(NanoReader.FatJet_msoftdrop[j] - W_MASS) > dmW ) continue;
	
	WVJJTree->bos_PuppiAK8_m_sd0 = NanoReader.FatJet_msoftdrop[j];
	WVJJTree->bos_PuppiAK8_m_sd0_corr = NanoReader.FatJet_msoftdrop[j];
	WVJJTree->bos_PuppiAK8_tau2tau1 = NanoReader.FatJet_tau2[j]/NanoReader.FatJet_tau1[j];
	WVJJTree->bos_PuppiAK8_pt = NanoReader.FatJet_pt[j];
	WVJJTree->bos_PuppiAK8_eta = NanoReader.FatJet_eta[j];
	WVJJTree->bos_PuppiAK8_phi = NanoReader.FatJet_phi[j];
	
	WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp = NanoReader.FatJet_msoftdrop_jesTotalUp[j];
	WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn = NanoReader.FatJet_msoftdrop_jesTotalDown[j];
	WVJJTree->bos_PuppiAK8_pt_scaleUp = NanoReader.FatJet_pt_jesTotalUp[j];
	WVJJTree->bos_PuppiAK8_pt_scaleDn = NanoReader.FatJet_pt_jesTotalDown[j];
	
	dmW = fabs(NanoReader.FatJet_msoftdrop[j] - W_MASS);
	nGoodFatJet++;
      }
      
      goodJetIndex.clear();
      for (uint j=0; j<NanoReader.nJet; j++) {
	
	//jet energy scale variations
	if ( NanoReader.Jet_pt[j] < AK4_PT_CUT && NanoReader.Jet_pt_jesTotalUp[j] < AK4_PT_CUT && 
	     NanoReader.Jet_pt_jesTotalDown[j] < AK4_PT_CUT) continue;
	//jet ID??
	
	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
	if (NanoReader.Jet_eta[j]<2.4 && NanoReader.Jet_pt[j]>30) {
	  if (NanoReader.Jet_btagDeepB[j] > 0.1241) WVJJTree->nBtag_loose++;
          if (NanoReader.Jet_btagDeepB[j] > 0.4184) WVJJTree->nBtag_medium++;
          if (NanoReader.Jet_btagDeepB[j] > 0.7527) WVJJTree->nBtag_tight++;
	}
	
	bool isClean=true;
	// object cleaning
	
	if (nGoodFatJet>0) {
	  if (deltaR(WVJJTree->bos_PuppiAK8_eta, WVJJTree->bos_PuppiAK8_phi,
		     NanoReader.Jet_eta[j], NanoReader.Jet_phi[j]) < AK4_AK8_DR_CUT) {
	    isClean = false;
	  }
	}
	
	for ( std::size_t k=0; k<goodJetIndex.size(); k++) {
	  if (deltaR(NanoReader.Jet_eta[goodJetIndex.at(k)], NanoReader.Jet_phi[goodJetIndex.at(k)],
		     NanoReader.Jet_eta[j], NanoReader.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     NanoReader.Jet_eta[j],   NanoReader.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     NanoReader.Jet_eta[j],   NanoReader.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	
	if ( isClean == false ) continue;

	if (NanoReader.Jet_pt[j]>30) WVJJTree->nJet30++;
	if (NanoReader.Jet_pt[j]>50) WVJJTree->nJet50++;

	goodJetIndex.push_back(j);
	
      }
      
      int nGoodDijet=0;
      
      uint sel1=1000, sel2=1000;
      if (nGoodFatJet==0) {
	TLorentzVector tmpV1, tmpV2;
	dmW=3000.0;
	for (uint j=0; j<goodJetIndex.size(); j++) {
	  if ( fabs( NanoReader.Jet_eta[goodJetIndex.at(j)] ) > AK4_ETA_CUT ) continue;
	  for(uint k=j+1; k<goodJetIndex.size(); k++) {
	    if ( fabs( NanoReader.Jet_eta[goodJetIndex.at(k)] ) > AK4_ETA_CUT ) continue;

	    TLorentzVector tmp1(0,0,0,0); 
	    tmp1.SetPtEtaPhiM( NanoReader.Jet_pt[goodJetIndex.at(j)], NanoReader.Jet_eta[goodJetIndex.at(j)],
			       NanoReader.Jet_phi[goodJetIndex.at(j)], NanoReader.Jet_mass[goodJetIndex.at(j)] );

	    TLorentzVector tmp2(0,0,0,0); 
	    tmp2.SetPtEtaPhiM( NanoReader.Jet_pt[goodJetIndex.at(k)], NanoReader.Jet_eta[goodJetIndex.at(k)],
			       NanoReader.Jet_phi[goodJetIndex.at(k)], NanoReader.Jet_mass[goodJetIndex.at(k)] );
			       
	    TLorentzVector tmpV=tmp1+tmp2;
	    
	    if (tmpV.M()<AK4_JJ_MIN_M || tmpV.M()>AK4_JJ_MAX_M) continue;
	    
	    if (fabs(tmpV.M()-W_MASS)>dmW) continue;
	    
	    WVJJTree->bos_j1_AK4_pt =  NanoReader.Jet_pt[goodJetIndex.at(j)];
	    WVJJTree->bos_j1_AK4_eta = NanoReader.Jet_eta[goodJetIndex.at(j)];
	    WVJJTree->bos_j1_AK4_phi = NanoReader.Jet_phi[goodJetIndex.at(j)];
	    WVJJTree->bos_j1_AK4_m =   NanoReader.Jet_mass[goodJetIndex.at(j)];
	    
	    WVJJTree->bos_j2_AK4_pt =  NanoReader.Jet_pt[goodJetIndex.at(k)];
	    WVJJTree->bos_j2_AK4_eta = NanoReader.Jet_eta[goodJetIndex.at(k)];
	    WVJJTree->bos_j2_AK4_phi = NanoReader.Jet_phi[goodJetIndex.at(k)];
	    WVJJTree->bos_j2_AK4_m =   NanoReader.Jet_mass[goodJetIndex.at(k)];
	    
	    WVJJTree->bos_AK4AK4_pt =  tmpV.Pt();
	    WVJJTree->bos_AK4AK4_eta = tmpV.Eta();
	    WVJJTree->bos_AK4AK4_phi = tmpV.Phi();
	    WVJJTree->bos_AK4AK4_m =   tmpV.M();
	    
	    sel1=j; sel2=k;
	    dmW=fabs(tmpV.M()-W_MASS);
	    nGoodDijet=1;
	    
	  }
	}
	
	if (nGoodDijet==0) continue;
	
	WVJJTree->bos_j1_AK4_pt_scaleUp = NanoReader.Jet_pt_jesTotalUp[goodJetIndex.at(sel1)];
	WVJJTree->bos_j1_AK4_pt_scaleDn = NanoReader.Jet_pt_jesTotalDown[goodJetIndex.at(sel1)];
	WVJJTree->bos_j1_AK4_m_scaleUp = NanoReader.Jet_mass_jesTotalUp[goodJetIndex.at(sel1)];
	WVJJTree->bos_j1_AK4_m_scaleDn = NanoReader.Jet_mass_jesTotalDown[goodJetIndex.at(sel1)];

	WVJJTree->bos_j2_AK4_pt_scaleUp = NanoReader.Jet_pt_jesTotalUp[goodJetIndex.at(sel2)];
	WVJJTree->bos_j2_AK4_pt_scaleDn = NanoReader.Jet_pt_jesTotalDown[goodJetIndex.at(sel2)];
	WVJJTree->bos_j2_AK4_m_scaleUp = NanoReader.Jet_mass_jesTotalUp[goodJetIndex.at(sel2)];
	WVJJTree->bos_j2_AK4_m_scaleDn = NanoReader.Jet_mass_jesTotalDown[goodJetIndex.at(sel2)];
	
	TLorentzVector tempBos1(0,0,0,0);
	TLorentzVector tempBos2(0,0,0,0);
	
	tempBos1.SetPtEtaPhiM(WVJJTree->bos_j1_AK4_pt_scaleUp, WVJJTree->bos_j1_AK4_eta,
			      WVJJTree->bos_j1_AK4_phi, WVJJTree->bos_j1_AK4_m_scaleUp);
	tempBos2.SetPtEtaPhiM(WVJJTree->bos_j2_AK4_pt_scaleUp, WVJJTree->bos_j2_AK4_eta,
			      WVJJTree->bos_j2_AK4_phi, WVJJTree->bos_j2_AK4_m_scaleUp);
	
	WVJJTree->bos_AK4AK4_pt_scaleUp = (tempBos1+tempBos2).Pt();
	WVJJTree->bos_AK4AK4_m_scaleUp = (tempBos1+tempBos2).M();
	
	tempBos1.SetPtEtaPhiM(WVJJTree->bos_j1_AK4_pt_scaleDn, WVJJTree->bos_j1_AK4_eta,
			      WVJJTree->bos_j1_AK4_phi, WVJJTree->bos_j1_AK4_m_scaleDn);
	tempBos2.SetPtEtaPhiM(WVJJTree->bos_j2_AK4_pt_scaleDn, WVJJTree->bos_j2_AK4_eta,
			      WVJJTree->bos_j2_AK4_phi, WVJJTree->bos_j2_AK4_m_scaleDn);
	
	WVJJTree->bos_AK4AK4_pt_scaleDn = (tempBos1+tempBos2).Pt();
	WVJJTree->bos_AK4AK4_m_scaleDn = (tempBos1+tempBos2).M();
	
      } //if (nGoodFatJet==0)
      
      //check we have a hadronic boson candidate
      if ( nGoodFatJet == 0 && nGoodDijet == 0 ) continue;
      
      float tmpMassMax = 0.0;
      int vbf1=-1, vbf2=-1;
      
      for (uint j=0; j<goodJetIndex.size(); j++) {
	if (j==sel1 || j==sel2) continue;
	for(uint k=j+1; k<goodJetIndex.size(); k++) {
	  if (k==sel1 || k==sel2) continue;

	  TLorentzVector tmp1(0,0,0,0);
	  tmp1.SetPtEtaPhiM( NanoReader.Jet_pt[goodJetIndex.at(j)], NanoReader.Jet_eta[goodJetIndex.at(j)],
			     NanoReader.Jet_phi[goodJetIndex.at(j)], NanoReader.Jet_mass[goodJetIndex.at(j)] );

	  TLorentzVector tmp2(0,0,0,0);
	  tmp2.SetPtEtaPhiM( NanoReader.Jet_pt[goodJetIndex.at(k)], NanoReader.Jet_eta[goodJetIndex.at(k)],
			     NanoReader.Jet_phi[goodJetIndex.at(k)], NanoReader.Jet_mass[goodJetIndex.at(k)] );

	  TLorentzVector tempVBF=tmp1+tmp2;

	  //require 2 jets be in opposite hemispheres
	  if ( NanoReader.Jet_eta[goodJetIndex.at(j)] * NanoReader.Jet_eta[goodJetIndex.at(k)] > 0 ) continue; 
	  if ( tempVBF.M() < VBF_MJJ_CUT ) continue;
	  if ( tempVBF.M() < tmpMassMax ) continue;
	  tmpMassMax = tempVBF.M();
	  vbf1=j; vbf2=k;
	}
      }
    
      if (vbf1==-1 && vbf2==-1) continue;
      
      WVJJTree->vbf1_AK4_pt = NanoReader.Jet_pt[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_eta = NanoReader.Jet_eta[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_phi = NanoReader.Jet_phi[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_m = NanoReader.Jet_mass[goodJetIndex.at(vbf1)];

      WVJJTree->vbf2_AK4_pt = NanoReader.Jet_pt[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_eta = NanoReader.Jet_eta[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_phi = NanoReader.Jet_phi[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_m = NanoReader.Jet_mass[goodJetIndex.at(vbf2)];

      TLorentzVector tempVBF1(0,0,0,0);
      TLorentzVector tempVBF2(0,0,0,0);

      tempVBF1.SetPtEtaPhiM(WVJJTree->vbf1_AK4_pt, WVJJTree->vbf1_AK4_eta,
			    WVJJTree->vbf1_AK4_phi, WVJJTree->vbf1_AK4_m);
      tempVBF2.SetPtEtaPhiM(WVJJTree->vbf2_AK4_pt, WVJJTree->vbf2_AK4_eta,
			    WVJJTree->vbf2_AK4_phi, WVJJTree->vbf2_AK4_m);

      WVJJTree->vbf_pt = (tempVBF1+tempVBF2).Pt();
      WVJJTree->vbf_eta = (tempVBF1+tempVBF2).Eta();
      WVJJTree->vbf_phi = (tempVBF1+tempVBF2).Phi();
      WVJJTree->vbf_m = (tempVBF1+tempVBF2).M();
      
      WVJJTree->vbf_deta = abs( WVJJTree->vbf2_AK4_eta - WVJJTree->vbf1_AK4_eta );
      
      WVJJTree->vbf1_AK4_pt_scaleUp = NanoReader.Jet_pt_jesTotalUp[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_pt_scaleDn = NanoReader.Jet_pt_jesTotalDown[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_m_scaleUp = NanoReader.Jet_mass_jesTotalUp[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_m_scaleDn = NanoReader.Jet_mass_jesTotalDown[goodJetIndex.at(vbf1)];

      WVJJTree->vbf2_AK4_pt_scaleUp = NanoReader.Jet_pt_jesTotalUp[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_pt_scaleDn = NanoReader.Jet_pt_jesTotalDown[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_m_scaleUp = NanoReader.Jet_mass_jesTotalUp[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_m_scaleDn = NanoReader.Jet_mass_jesTotalDown[goodJetIndex.at(vbf2)];

      tempVBF1.SetPtEtaPhiM(WVJJTree->vbf1_AK4_pt_scaleUp, WVJJTree->vbf1_AK4_eta,
			    WVJJTree->vbf1_AK4_phi, WVJJTree->vbf1_AK4_m_scaleUp);
      tempVBF2.SetPtEtaPhiM(WVJJTree->vbf2_AK4_pt_scaleUp, WVJJTree->vbf2_AK4_eta,
			    WVJJTree->vbf2_AK4_phi, WVJJTree->vbf2_AK4_m_scaleUp);
      
      WVJJTree->vbf_pt_scaleUp = (tempVBF1+tempVBF2).Pt();
      WVJJTree->vbf_m_scaleUp = (tempVBF1+tempVBF2).M();
      
      tempVBF1.SetPtEtaPhiM(WVJJTree->vbf1_AK4_pt_scaleDn, WVJJTree->vbf1_AK4_eta,
			    WVJJTree->vbf1_AK4_phi, WVJJTree->vbf1_AK4_m_scaleDn);
      tempVBF2.SetPtEtaPhiM(WVJJTree->vbf2_AK4_pt_scaleDn, WVJJTree->vbf2_AK4_eta,
			    WVJJTree->vbf2_AK4_phi, WVJJTree->vbf2_AK4_m_scaleDn);
      
      WVJJTree->vbf_pt_scaleDn = (tempVBF1+tempVBF2).Pt();
      WVJJTree->vbf_m_scaleDn = (tempVBF1+tempVBF2).M();
      
      TLorentzVector bosHad(0,0,0,0), bosHad_up(0,0,0,0), bosHad_dn(0,0,0,0);
      TLorentzVector bosLep(0,0,0,0), bosLep_up(0,0,0,0), bosLep_dn(0,0,0,0);
      
      //boosted event
      if (WVJJTree->bos_PuppiAK8_m_sd0_corr > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt, WVJJTree->bos_PuppiAK8_eta,
			    WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr);
	bosHad_up.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleUp, WVJJTree->bos_PuppiAK8_eta,
			       WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp);
	bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleDn, WVJJTree->bos_PuppiAK8_eta,
			       WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn);
      }
      //resolved event
      else if (WVJJTree->bos_AK4AK4_m > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt, WVJJTree->bos_AK4AK4_eta, 
			    WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m);
	bosHad_up.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleUp, WVJJTree->bos_AK4AK4_eta, 
			       WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleUp);
	bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleDn, WVJJTree->bos_AK4AK4_eta, 
			       WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleDn);
      }
      
      bosLep.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      //lepton variations not implemented
      bosLep_up.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      bosLep_dn.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      
      TLorentzVector diBos = bosHad+bosLep;
      TLorentzVector diBos_up = bosHad_up+bosLep_up;
      TLorentzVector diBos_dn = bosHad_dn+bosLep_dn;
      
      WVJJTree->dibos_m = diBos.M();
      WVJJTree->dibos_pt = diBos.Pt();
      WVJJTree->dibos_eta = diBos.Eta();
      WVJJTree->dibos_phi = diBos.Phi();
      
      WVJJTree->dibos_m_scaleUp = diBos_up.M();
      WVJJTree->dibos_m_scaleDn = diBos_dn.M();
      WVJJTree->dibos_pt_scaleUp = diBos_up.Pt();
      WVJJTree->dibos_pt_scaleDn = diBos_dn.Pt();
      
      if (WVJJTree->lep2_pt < 0) {
	WVJJTree->bosCent = std::min( std::min(bosHad.Eta(), bosLep.Eta()) - std::min(WVJJTree->vbf1_AK4_eta, WVJJTree->vbf2_AK4_eta), 
				      std::max(WVJJTree->vbf1_AK4_eta, WVJJTree->vbf2_AK4_eta) - std::max(bosHad.Eta(), bosLep.Eta()) );
      }
      
      WVJJTree->zeppLep = bosLep.Eta() - 0.5*(WVJJTree->vbf1_AK4_eta + WVJJTree->vbf2_AK4_eta);
      WVJJTree->zeppHad = bosHad.Eta() - 0.5*(WVJJTree->vbf1_AK4_eta + WVJJTree->vbf2_AK4_eta);

      if (isMC==1) {

	WVJJTree->nScaleWeight = NanoReader.nLHEScaleWeight;
	WVJJTree->nPdfWeight = NanoReader.nLHEPdfWeight;

	for (uint j=0; j<WVJJTree->nScaleWeight; j++) {
	  //LHE scale variation weights (w_var / w_nominal); [0] is MUR="0.5" MUF="0.5"; 
	  //[1] is MUR="0.5" MUF="1.0"; [2] is MUR="0.5" MUF="2.0"; [3] is MUR="1.0" MUF="0.5"; 
	  //[4] is MUR="1.0" MUF="2.0"; [5] is MUR="2.0" MUF="0.5"; [6] is MUR="2.0" MUF="1.0"; 
	  //[7] is MUR="2.0" MUF="2.0"
	  WVJJTree->scaleWeight[j]=NanoReader.LHEScaleWeight[j];
	}
	for (uint j=0; j<WVJJTree->nPdfWeight; j++) {
	  //LHE pdf variation weights (w_var / w_nominal) for LHA IDs 91400 - 91432
	  WVJJTree->pdfWeight[j]=NanoReader.LHEPdfWeight[j];
	}

      }

      WVJJTree->btagWeight = NanoReader.btagWeight_DeepCSVB;

      ot->Fill();
    }
    
    delete t; 
    delete r;
    delete f;
    t=0; 
    r=0;
    f=0;
  }

  of->Write();
  of->Close();
  
  return 0;

}
