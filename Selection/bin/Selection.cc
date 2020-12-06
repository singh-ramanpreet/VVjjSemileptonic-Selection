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

#include "WVJJAna/Selection/interface/NanoReader.hh"
#include "WVJJAna/Selection/interface/ScaleFactors.hh"
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

  ScaleFactors scaleFactor(era);

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

  Bool_t has_HLT_IsoMu22=false, has_HLT_IsoTkMu22=false, has_HLT_IsoMu24=false, has_HLT_IsoTkMu24=false,
    has_HLT_IsoMu27=false, has_HLT_IsoMu30=false, has_HLT_Mu50=false,
    has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ=false, has_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ=false,
    has_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ=false, has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8=false,
    has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8=false, has_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8=false,
    has_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8=false, has_HLT_Ele27_WPTight_Gsf=false, has_HLT_Ele28_WPTight_Gsf=false,
    has_HLT_Ele32_WPTight_Gsf=false, has_HLT_Ele35_WPLoose_Gsf=false, has_HLT_Ele35_WPTight_Gsf=false, has_HLT_Ele38_WPTight_Gsf=false,
    has_HLT_Ele40_WPTight_Gsf=false, has_HLT_Ele25_eta2p1_WPTight_Gsf=false, has_HLT_Ele27_eta2p1_WPLoose_Gsf=false,
    has_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ=false, has_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL=false,
    has_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG=false,
    has_HLT_DoubleEle33_CaloIdL_MW=false, has_HLT_DoubleEle25_CaloIdL_MW=false, has_HLT_DoubleEle27_CaloIdL_MW=false;

  std::ifstream ifs;
  ifs.open(inputFile.data());
  assert(ifs.is_open());
  std::string line;
  int lineCount=0;
  while (getline(ifs,line)) {
    std::stringstream ss(line);
    std::string filetoopen;
    ss >> filetoopen;
    
    lineCount+=1;

    f = TFile::Open(TString("root://cmseos.fnal.gov/")+TString(filetoopen),"read");
    //f = TFile::Open(TString("root://xrootd-cms.infn.it/")+TString(filetoopen),"read");
    t = (TTree *)f->Get("Events");
    r = (TTree *)f->Get("Runs");
    if (t==NULL) continue;
    
    //std::cout << filetoopen << std::endl;
    
    NanoReader nr = NanoReader(t);

    NanoReader NanoWeightReader = NanoReader(r);
    NanoWeightReader.GetEntry(0);
    float genEventSumw=0.0;
    if (isMC && nanoVersion == 6) {
      genEventSumw = *NanoWeightReader.genEventSumw_;
    }
    else if (isMC && nanoVersion == 7) {
      genEventSumw = *NanoWeightReader.genEventSumw;
    }
    totalEvents->SetBinContent(2,totalEvents->GetBinContent(2)+genEventSumw);

    //check if tree has these hlt branches
    if (lineCount == 1){
      has_HLT_IsoMu22 = t->GetBranchStatus("HLT_IsoMu22");
      has_HLT_IsoTkMu22 = t->GetBranchStatus("HLT_IsoTkMu22");
      has_HLT_IsoMu24 = t->GetBranchStatus("HLT_IsoMu24");
      has_HLT_IsoTkMu24 = t->GetBranchStatus("HLT_IsoTkMu24");
      has_HLT_IsoMu27 = t->GetBranchStatus("HLT_IsoMu27");
      has_HLT_IsoMu30 = t->GetBranchStatus("HLT_IsoMu30");
      has_HLT_Mu50 = t->GetBranchStatus("HLT_Mu50");

      has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = t->GetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
      has_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = t->GetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
      has_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = t->GetBranchStatus("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
      has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 = t->GetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
      has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 = t->GetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
      has_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 = t->GetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8");
      has_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 = t->GetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8");

      has_HLT_Ele27_WPTight_Gsf = t->GetBranchStatus("HLT_Ele27_WPTight_Gsf");
      has_HLT_Ele28_WPTight_Gsf = t->GetBranchStatus("HLT_Ele28_WPTight_Gsf");
      has_HLT_Ele32_WPTight_Gsf = t->GetBranchStatus("HLT_Ele32_WPTight_Gsf");
      has_HLT_Ele35_WPLoose_Gsf = t->GetBranchStatus("HLT_Ele35_WPLoose_Gsf");
      has_HLT_Ele35_WPTight_Gsf = t->GetBranchStatus("HLT_Ele35_WPTight_Gsf");
      has_HLT_Ele38_WPTight_Gsf = t->GetBranchStatus("HLT_Ele38_WPTight_Gsf");
      has_HLT_Ele40_WPTight_Gsf = t->GetBranchStatus("HLT_Ele40_WPTight_Gsf");
      has_HLT_Ele25_eta2p1_WPTight_Gsf = t->GetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf");
      has_HLT_Ele27_eta2p1_WPLoose_Gsf = t->GetBranchStatus("HLT_Ele27_eta2p1_WPLoose_Gsf");

      has_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = t->GetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
      has_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = t->GetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
      has_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG = t->GetBranchStatus("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG");
      has_HLT_DoubleEle33_CaloIdL_MW = t->GetBranchStatus("HLT_DoubleEle33_CaloIdL_MW");
      has_HLT_DoubleEle25_CaloIdL_MW = t->GetBranchStatus("HLT_DoubleEle25_CaloIdL_MW");
      has_HLT_DoubleEle27_CaloIdL_MW = t->GetBranchStatus("HLT_DoubleEle27_CaloIdL_MW");
    }
      
    for (uint i=0; i < t->GetEntries(); i++) {
    //for (uint i=0; i < 1000000; i++) {
      WVJJTree->clearVars();
      nr.GetEntry(i);

      if (i%100000==0) std::cout <<"file " << lineCount << ": event " << i << std::endl;

      if (isMC==1) {
        WVJJTree->genWeight=*nr.Generator_weight;
      }

      WVJJTree->trigger_1Mu = ((has_HLT_IsoMu22 ? *nr.HLT_IsoMu22 : 0) ||
                               (has_HLT_IsoTkMu22 ? *nr.HLT_IsoTkMu22 : 0) ||
                               (has_HLT_IsoMu24 ? *nr.HLT_IsoMu24 : 0) ||
                               (has_HLT_IsoTkMu24 ? *nr.HLT_IsoTkMu24 : 0) ||
                               (has_HLT_IsoMu27 ? *nr.HLT_IsoMu27 : 0) ||
                               (has_HLT_IsoMu30 ? *nr.HLT_IsoMu30 : 0) ||
                               (has_HLT_Mu50 ? *nr.HLT_Mu50 : 0));

      WVJJTree->trigger_2Mu = ((has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ ? *nr.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ : 0) ||
                               (has_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ ? *nr.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ : 0) ||
                               (has_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ ? *nr.HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ : 0) ||
                               (has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 ? *nr.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 : 0) ||
                               (has_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ? *nr.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 : 0) ||
                               (has_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 ? *nr.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 : 0) ||
                               (has_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 ? *nr.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 : 0));

      WVJJTree->trigger_1El = ((has_HLT_Ele27_WPTight_Gsf ? *nr.HLT_Ele27_WPTight_Gsf : 0) ||
                               (has_HLT_Ele28_WPTight_Gsf ? *nr.HLT_Ele28_WPTight_Gsf : 0) ||
                               (has_HLT_Ele32_WPTight_Gsf ? *nr.HLT_Ele32_WPTight_Gsf : 0) ||
                               (has_HLT_Ele35_WPLoose_Gsf ? *nr.HLT_Ele35_WPLoose_Gsf : 0) ||
                               (has_HLT_Ele35_WPTight_Gsf ? *nr.HLT_Ele35_WPTight_Gsf : 0) ||
                               (has_HLT_Ele38_WPTight_Gsf ? *nr.HLT_Ele38_WPTight_Gsf : 0) ||
                               (has_HLT_Ele40_WPTight_Gsf ? *nr.HLT_Ele40_WPTight_Gsf : 0) ||
                               (has_HLT_Ele25_eta2p1_WPTight_Gsf ? *nr.HLT_Ele25_eta2p1_WPTight_Gsf : 0) ||
                               (has_HLT_Ele27_eta2p1_WPLoose_Gsf ? *nr.HLT_Ele27_eta2p1_WPLoose_Gsf : 0));

      WVJJTree->trigger_2El = ((has_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ? *nr.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ : 0) ||
                               (has_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ? *nr.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL : 0) ||
                               (has_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG ? *nr.HLT_DiEle27_WPTightCaloOnly_L1DoubleEG : 0) ||
                               (has_HLT_DoubleEle33_CaloIdL_MW ? *nr.HLT_DoubleEle33_CaloIdL_MW : 0) ||
                               (has_HLT_DoubleEle25_CaloIdL_MW ? *nr.HLT_DoubleEle25_CaloIdL_MW : 0) ||
                               (has_HLT_DoubleEle27_CaloIdL_MW ? *nr.HLT_DoubleEle27_CaloIdL_MW : 0));

      //std::cout << "passed trigger: ";
      //if (WVJJTree->trigger_1Mu) std::cout << "1 muon ";
      //if (WVJJTree->trigger_2Mu) std::cout << "2 muon ";
      //if (WVJJTree->trigger_1El) std::cout << "1 ele ";
      //if (WVJJTree->trigger_2El) std::cout << "2 ele ";
      //std::cout << std::endl;
      
      if ( ! ( WVJJTree->trigger_1Mu || WVJJTree->trigger_2Mu || WVJJTree->trigger_1El || WVJJTree->trigger_2El ) ) continue;

      tightMuon.clear();
      tightEle.clear();

      WVJJTree->run = *nr.run;
      WVJJTree->evt = *nr.event;
      WVJJTree->ls = *nr.luminosityBlock;
      
      WVJJTree->nPV = *nr.PV_npvsGood;
      if (isMC) WVJJTree->nPU_mean = *nr.Pileup_nPU;
      
      WVJJTree->puWeight = scaleFactor.GetPUWeight(WVJJTree->nPU_mean, 0);

      // LEPTON SELECTION
      
      int nTightEle=0;
      int nTightMu=0;
      int nVetoEle=0;
      int nVetoMu=0;
      
      for (uint j=0; j < *nr.nMuon; j++) {
	if ( abs(nr.Muon_eta[j]) > MU_ETA_CUT ) continue;
	//using conservative uncertainty value of 3%
	if ( 1.03*nr.Muon_pt[j] < LEP_PT_VETO_CUT ) continue;
	
	if (!nr.Muon_looseId[j]) continue;
	if (nr.Muon_pfRelIso04_all[j]>0.25) continue;
	nVetoMu++;

	//using conservative uncertainty value of 3%
	if ( 1.03*nr.Muon_pt[j] < MU_PT_CUT ) continue;
	if (!nr.Muon_tightId[j]) continue;
	if (nr.Muon_pt[j] > 20 && abs(nr.Muon_dxy[j])>0.01) continue;
        if (nr.Muon_pt[j] < 20 && abs(nr.Muon_dxy[j])>0.02) continue;
        if (abs(nr.Muon_dz[j])>0.1) continue;
	if (nr.Muon_pfRelIso04_all[j]>0.15) continue;

	nTightMu++;
	tightMuon.push_back(TLorentzVector(0,0,0,0));
	tightMuon.back().SetPtEtaPhiM(nr.Muon_pt[j], nr.Muon_eta[j], 
				      nr.Muon_phi[j], MUON_MASS);
	
	if ( nr.Muon_pt[j] > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_dxy = WVJJTree->lep1_dxy;
          WVJJTree->lep2_dz = WVJJTree->lep1_dz;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = nr.Muon_pt[j];
	  WVJJTree->lep1_eta = nr.Muon_eta[j];
	  WVJJTree->lep1_phi = nr.Muon_phi[j];
	  WVJJTree->lep1_m = MUON_MASS;
	  WVJJTree->lep1_iso = nr.Muon_pfRelIso04_all[j];
	  WVJJTree->lep1_dxy = nr.Muon_dxy[j];
	  WVJJTree->lep1_dz = nr.Muon_dz[j];
	  WVJJTree->lep1_q = nr.Muon_charge[j];
	  
	}
	else if ( nr.Muon_pt[j] > WVJJTree->lep2_pt ) {
	  //if (WVJJTree->lep1_q*nr.Muon_charge[j]>0) continue;

	  WVJJTree->lep2_pt = nr.Muon_pt[j];
	  WVJJTree->lep2_eta = nr.Muon_eta[j];
	  WVJJTree->lep2_phi = nr.Muon_phi[j];
	  WVJJTree->lep2_m = MUON_MASS;
	  WVJJTree->lep2_iso = nr.Muon_pfRelIso04_all[j];
	  WVJJTree->lep2_dxy = nr.Muon_dxy[j];
	  WVJJTree->lep2_dz = nr.Muon_dz[j];
	  WVJJTree->lep2_q = nr.Muon_charge[j];
	  
	  }
      }
      
      for (uint j=0; j < *nr.nElectron; j++) {
	if ( abs(nr.Electron_eta[j]) > EL_ETA_CUT ) continue;
	//using conservative uncertainty value of 3%
	if ( 1.03*nr.Electron_pt[j] < LEP_PT_VETO_CUT ) continue;
	
	//cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
	if ( nr.Electron_cutBased[j]<2 ) continue;
	if ( abs(nr.Electron_eta[j])>1.5) {
          if (nr.Electron_sieie[j]>0.03) continue;
          if (nr.Electron_eInvMinusPInv[j]>0.014) continue;
          if (abs(nr.Electron_dxy[j])>0.1) continue;
          if (abs(nr.Electron_dz[j])>0.2) continue;
        }
        else {
          if (abs(nr.Electron_dxy[j])>0.05) continue;
          if (abs(nr.Electron_dz[j])>0.1) continue;
        }

	nVetoEle++;

	//using conservative uncertainty value of 3%
	if ( 1.03*nr.Electron_pt[j] < EL_PT_CUT ) continue;
	
	//if (nr.Electron_cutBased[j]<4) continue;
	if (!nr.Electron_mvaFall17V2Iso_WP90[j]) continue;
	nTightEle++;
	
	tightEle.push_back(TLorentzVector(0,0,0,0));
	tightEle.back().SetPtEtaPhiM(nr.Electron_pt[j],nr.Electron_eta[j],
				     nr.Electron_phi[j],ELE_MASS);
	
	//don't try to select electrons unless we don't already
	//have muons
	if (WVJJTree->lep1_m == MUON_MASS) continue;
	
	if ( nr.Electron_pt[j] > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_dxy = WVJJTree->lep1_dxy;
	  WVJJTree->lep2_dz = WVJJTree->lep1_dz;
	  WVJJTree->lep2_sip3d = WVJJTree->lep1_sip3d;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = nr.Electron_pt[j];
	  WVJJTree->lep1_eta = nr.Electron_eta[j];
	  WVJJTree->lep1_phi = nr.Electron_phi[j];
	  WVJJTree->lep1_m = ELE_MASS;
	  WVJJTree->lep1_iso = nr.Electron_pfRelIso03_all[j];
	  WVJJTree->lep1_dxy = nr.Electron_dxy[j];
	  WVJJTree->lep1_dz = nr.Electron_dz[j];
	  WVJJTree->lep1_sip3d = nr.Electron_sip3d[j];
	  WVJJTree->lep1_q = nr.Electron_charge[j];
	  
	}
	else if ( nr.Electron_pt[j] > WVJJTree->lep2_pt ) {

	  WVJJTree->lep2_pt = nr.Electron_pt[j];
	  WVJJTree->lep2_eta = nr.Electron_eta[j];
	  WVJJTree->lep2_phi = nr.Electron_phi[j];
	  WVJJTree->lep2_m = ELE_MASS;
	  WVJJTree->lep2_iso = nr.Electron_pfRelIso03_all[j];
	  WVJJTree->lep2_dxy = nr.Electron_dxy[j];
	  WVJJTree->lep2_dz = nr.Electron_dz[j];
	  WVJJTree->lep2_sip3d = nr.Electron_sip3d[j];
	  WVJJTree->lep2_q = nr.Electron_charge[j];
	  
	}
      }

      //check conditions

      bool passLepSel = true;
      if(!(WVJJTree->lep1_pt>0)) passLepSel=false;
      if ((nTightMu+nTightEle)==0) passLepSel=false; //no leptons with required ID
      if((nVetoEle+nVetoMu)>2) passLepSel=false;
      if(nTightMu>0 && nVetoEle>0) passLepSel=false;
      if(nTightEle>0 && nVetoMu>0) passLepSel=false;
      if(nTightMu==1 && nVetoMu>1) passLepSel=false;
      if(nTightEle==1 && nVetoEle>1) passLepSel=false;

      //if (0) {
      if (passLepSel==false && (nTightMu+nTightEle)==0 && ((nVetoEle>0) ^ (nVetoMu>0))) {
	if (nVetoEle>0) {

	  for (uint j=0; j < *nr.nElectron; j++) {
	    if ( abs(nr.Electron_eta[j]) > EL_ETA_CUT ) continue;
	    if ( 1.03*nr.Electron_pt[j] < EL_PT_CUT ) continue;
	    if ( nr.Electron_cutBased[j]<2 ) continue;

	    if ( nr.Electron_pt[j] > WVJJTree->lep1_pt ) {
	      WVJJTree->lep1_pt = nr.Electron_pt[j];
	      WVJJTree->lep1_eta = nr.Electron_eta[j];
	      WVJJTree->lep1_phi = nr.Electron_phi[j];
	      WVJJTree->lep1_m = ELE_MASS;
	      WVJJTree->lep1_iso = nr.Electron_pfRelIso03_all[j];
	      WVJJTree->lep1_dxy = nr.Electron_dxy[j];
	      WVJJTree->lep1_dz = nr.Electron_dz[j];
	      WVJJTree->lep1_sip3d = nr.Electron_sip3d[j];
	      WVJJTree->lep1_q = nr.Electron_charge[j];

	    }
	  }

	  tightEle.push_back(TLorentzVector(0,0,0,0));
	  tightEle.back().SetPtEtaPhiM(WVJJTree->lep1_pt, WVJJTree->lep1_eta, 
				       WVJJTree->lep1_phi, ELE_MASS);

	  WVJJTree->isAntiIso=true;

	  WVJJTree->lepFakeRate = scaleFactor.GetLeptonFakeWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, 11);

	}

	if (nVetoMu>0) {

	  for (uint j=0; j < *nr.nMuon; j++) {
	    if ( abs(nr.Muon_eta[j]) > MU_ETA_CUT ) continue;
	    if (nr.Muon_pfRelIso04_all[j]>0.25) continue;
	    if ( 1.03*nr.Muon_pt[j] < MU_PT_CUT ) continue;
	    if (!nr.Muon_tightId[j]) continue;

	    if ( nr.Muon_pt[j] > WVJJTree->lep1_pt ) {
	      WVJJTree->lep1_pt = nr.Muon_pt[j];
	      WVJJTree->lep1_eta = nr.Muon_eta[j];
	      WVJJTree->lep1_phi = nr.Muon_phi[j];
	      WVJJTree->lep1_m = MUON_MASS;
	      WVJJTree->lep1_iso = nr.Muon_pfRelIso04_all[j];
	      WVJJTree->lep1_dxy = nr.Muon_dxy[j];
	      WVJJTree->lep1_dz = nr.Muon_dz[j];
	      WVJJTree->lep1_q = nr.Muon_charge[j];
	    }
	  }

	  tightMuon.push_back(TLorentzVector(0,0,0,0));
	  tightMuon.back().SetPtEtaPhiM(WVJJTree->lep1_pt, WVJJTree->lep1_eta,
					WVJJTree->lep1_phi,MUON_MASS);

	  WVJJTree->isAntiIso=true;
	  WVJJTree->lepFakeRate = scaleFactor.GetLeptonFakeWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, 13);

	}
      }
      
      if (!passLepSel && !WVJJTree->isAntiIso) continue;


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
      WVJJTree->lep1_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
      if (WVJJTree->lep2_pt>0) {
	WVJJTree->lep2_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep2_pt, WVJJTree->lep2_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
      }
      
      
      // MET

      if (*nr.MET_pt < 0) continue;
      WVJJTree->MET = isMC ? *nr.MET_T1Smear_pt : *nr.MET_T1_pt;
      WVJJTree->MET_phi = isMC ? *nr.MET_T1Smear_phi : *nr.MET_T1_phi;
      if (era==2017) {
        WVJJTree->MET_2017 = *nr.METFixEE2017_pt;
      }

      WVJJTree->PuppiMET = *nr.PuppiMET_pt;
      WVJJTree->PuppiMET_phi = *nr.PuppiMET_phi;

      if (isMC) {
        //WVJJTree->MET_scaleUp = *nr.MET_pt_jesTotalUp;
        //WVJJTree->MET_scaleDn = *nr.MET_pt_jesTotalDown;
      } 


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

      for (uint j=0; j<*nr.nFatJet; j++) {

	if ( fabs(nr.FatJet_eta[j]) > AK8_MAX_ETA ) continue;

	if ( isMC ) {
	
	  if ( ! (nr.FatJet_pt_nom[j]>AK8_MIN_PT || nr.FatJet_pt_jesTotalUp[j]>AK8_MIN_PT || 
		  nr.FatJet_pt_jesTotalDown[j]>AK8_MIN_PT) ) continue;
	  
	  if ( ! (nr.FatJet_msoftdrop[j]>AK8_MIN_SDM || nr.FatJet_msoftdrop_jesTotalUp[j]>AK8_MIN_SDM ||
		  nr.FatJet_msoftdrop_jesTotalDown[j]>AK8_MIN_SDM) ) continue;
	  
	  if ( ! (nr.FatJet_msoftdrop[j]<AK8_MAX_SDM || nr.FatJet_msoftdrop_jesTotalUp[j]<AK8_MAX_SDM ||
		  nr.FatJet_msoftdrop_jesTotalDown[j]<AK8_MAX_SDM) ) continue;
	
	}
	
	else {

	  if ( nr.FatJet_pt_nom[j]<AK8_MIN_PT ) continue;
	  if ( nr.FatJet_msoftdrop[j]<AK8_MIN_SDM ) continue;
	  if ( nr.FatJet_msoftdrop[j]>AK8_MAX_SDM ) continue;

	}
	
	bool isClean=true;
	//lepton cleaning
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     nr.FatJet_eta[j], nr.FatJet_phi[j]) < AK8_LEP_DR_CUT)
	    isClean = false;
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     nr.FatJet_eta[j], nr.FatJet_phi[j]) < AK8_LEP_DR_CUT)
	    isClean = false;
	}

	if ( isClean == false ) continue;

	if ( fabs(nr.FatJet_msoftdrop[j] - W_MASS) > dmW ) continue;
	
	WVJJTree->bos_PuppiAK8_m_sd0 = nr.FatJet_msoftdrop[j];
	WVJJTree->bos_PuppiAK8_m_sd0_corr = nr.FatJet_msoftdrop[j];
	WVJJTree->bos_PuppiAK8_tau2tau1 = nr.FatJet_tau2[j]/nr.FatJet_tau1[j];
	WVJJTree->bos_PuppiAK8_pt = nr.FatJet_pt[j];
	WVJJTree->bos_PuppiAK8_eta = nr.FatJet_eta[j];
	WVJJTree->bos_PuppiAK8_phi = nr.FatJet_phi[j];

	if (isMC) {
	
	  WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp = nr.FatJet_msoftdrop_jesTotalUp[j];
	  WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn = nr.FatJet_msoftdrop_jesTotalDown[j];
	  WVJJTree->bos_PuppiAK8_pt_scaleUp = nr.FatJet_pt_jesTotalUp[j];
	  WVJJTree->bos_PuppiAK8_pt_scaleDn = nr.FatJet_pt_jesTotalDown[j];
	
	}

	dmW = fabs(nr.FatJet_msoftdrop[j] - W_MASS);
	nGoodFatJet++;
      }
      
      goodJetIndex.clear();
      for (uint j=0; j<*nr.nJet; j++) {
	
	//jet energy scale variations
	if ( isMC && ( nr.Jet_pt_nom[j] < AK4_PT_CUT && nr.Jet_pt_jesTotalUp[j] < AK4_PT_CUT && 
		       nr.Jet_pt_jesTotalDown[j] < AK4_PT_CUT ) ) continue;
	else if ( !isMC && nr.Jet_pt_nom[j] < AK4_PT_CUT ) continue;
	//jet ID??
	
	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
	if (fabs(nr.Jet_eta[j])<2.4 && nr.Jet_pt_nom[j]>30) {
	    
	 if (era==2018) {
	    if (nr.Jet_btagDeepB[j] > 0.1241) WVJJTree->nBtag_loose++;
	    if (nr.Jet_btagDeepB[j] > 0.4184) WVJJTree->nBtag_medium++;
	    if (nr.Jet_btagDeepB[j] > 0.7527) WVJJTree->nBtag_tight++;
          }
	  if (era==2017) {
	    if (nr.Jet_btagDeepB[j] > 0.1522) WVJJTree->nBtag_loose++;
	    if (nr.Jet_btagDeepB[j] > 0.4941) WVJJTree->nBtag_medium++;
	    if (nr.Jet_btagDeepB[j] > 0.8001) WVJJTree->nBtag_tight++;
          }
	  if (era==2016) {
	    if (nr.Jet_btagDeepB[j] > 0.2217) WVJJTree->nBtag_loose++;
	    if (nr.Jet_btagDeepB[j] > 0.6321) WVJJTree->nBtag_medium++;
	    if (nr.Jet_btagDeepB[j] > 0.8953) WVJJTree->nBtag_tight++;
          }
	  }

	bool isClean=true;
	// object cleaning
	
	if (nGoodFatJet>0) {
	  if (deltaR(WVJJTree->bos_PuppiAK8_eta, WVJJTree->bos_PuppiAK8_phi,
		     nr.Jet_eta[j], nr.Jet_phi[j]) < AK4_AK8_DR_CUT) {
	    isClean = false;
	  }
	}
	
	for ( std::size_t k=0; k<goodJetIndex.size(); k++) {
	  if (deltaR(nr.Jet_eta[goodJetIndex.at(k)], nr.Jet_phi[goodJetIndex.at(k)],
		     nr.Jet_eta[j], nr.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     nr.Jet_eta[j],   nr.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     nr.Jet_eta[j],   nr.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	
	if ( isClean == false ) continue;

	if (nr.Jet_pt_nom[j]>30) WVJJTree->nJet30++;
	if (nr.Jet_pt_nom[j]>50) WVJJTree->nJet50++;

	goodJetIndex.push_back(j);
	
      }
      
      int nGoodDijet=0;
      
      uint sel1=1000, sel2=1000;
      if (nGoodFatJet==0) {
	TLorentzVector tmpV1, tmpV2;
	dmW=3000.0;

	for (uint j=0; j<goodJetIndex.size(); j++) {
	  if ( fabs( nr.Jet_eta[goodJetIndex.at(j)] ) > AK4_ETA_CUT ) continue;
	  for(uint k=j+1; k<goodJetIndex.size(); k++) {
	    if ( fabs( nr.Jet_eta[goodJetIndex.at(k)] ) > AK4_ETA_CUT ) continue;

	    TLorentzVector tmp1(0,0,0,0); 
	    tmp1.SetPtEtaPhiM( nr.Jet_pt_nom[goodJetIndex.at(j)], nr.Jet_eta[goodJetIndex.at(j)],
			       nr.Jet_phi[goodJetIndex.at(j)], nr.Jet_mass[goodJetIndex.at(j)] );

	    TLorentzVector tmp2(0,0,0,0); 
	    tmp2.SetPtEtaPhiM( nr.Jet_pt_nom[goodJetIndex.at(k)], nr.Jet_eta[goodJetIndex.at(k)],
			       nr.Jet_phi[goodJetIndex.at(k)], nr.Jet_mass[goodJetIndex.at(k)] );
			       
	    TLorentzVector tmpV=tmp1+tmp2;
	    
	    if (tmpV.M()<AK4_JJ_MIN_M || tmpV.M()>AK4_JJ_MAX_M) continue;
	    
	    if (fabs(tmpV.M()-W_MASS)>dmW) continue;
	    
	    WVJJTree->bos_j1_AK4_pt =  nr.Jet_pt_nom[goodJetIndex.at(j)];
	    WVJJTree->bos_j1_AK4_eta = nr.Jet_eta[goodJetIndex.at(j)];
	    WVJJTree->bos_j1_AK4_phi = nr.Jet_phi[goodJetIndex.at(j)];
	    WVJJTree->bos_j1_AK4_m =   nr.Jet_mass[goodJetIndex.at(j)];
	    WVJJTree->bos_j1_AK4_qgid = nr.Jet_qgl[goodJetIndex.at(j)];
	    
	    WVJJTree->bos_j2_AK4_pt =  nr.Jet_pt_nom[goodJetIndex.at(k)];
	    WVJJTree->bos_j2_AK4_eta = nr.Jet_eta[goodJetIndex.at(k)];
	    WVJJTree->bos_j2_AK4_phi = nr.Jet_phi[goodJetIndex.at(k)];
	    WVJJTree->bos_j2_AK4_m =   nr.Jet_mass[goodJetIndex.at(k)];
	    WVJJTree->bos_j2_AK4_qgid = nr.Jet_qgl[goodJetIndex.at(k)];
	    
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

	if (isMC) {
	
	  WVJJTree->bos_j1_AK4_pt_scaleUp = nr.Jet_pt_jesTotalUp[goodJetIndex.at(sel1)];
	  WVJJTree->bos_j1_AK4_pt_scaleDn = nr.Jet_pt_jesTotalDown[goodJetIndex.at(sel1)];
	  WVJJTree->bos_j1_AK4_m_scaleUp = nr.Jet_mass_jesTotalUp[goodJetIndex.at(sel1)];
	  WVJJTree->bos_j1_AK4_m_scaleDn = nr.Jet_mass_jesTotalDown[goodJetIndex.at(sel1)];
	  
	  WVJJTree->bos_j2_AK4_pt_scaleUp = nr.Jet_pt_jesTotalUp[goodJetIndex.at(sel2)];
	  WVJJTree->bos_j2_AK4_pt_scaleDn = nr.Jet_pt_jesTotalDown[goodJetIndex.at(sel2)];
	  WVJJTree->bos_j2_AK4_m_scaleUp = nr.Jet_mass_jesTotalUp[goodJetIndex.at(sel2)];
	  WVJJTree->bos_j2_AK4_m_scaleDn = nr.Jet_mass_jesTotalDown[goodJetIndex.at(sel2)];
	
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
	
	}
	
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
	  tmp1.SetPtEtaPhiM( nr.Jet_pt_nom[goodJetIndex.at(j)], nr.Jet_eta[goodJetIndex.at(j)],
			     nr.Jet_phi[goodJetIndex.at(j)], nr.Jet_mass[goodJetIndex.at(j)] );

	  TLorentzVector tmp2(0,0,0,0);
	  tmp2.SetPtEtaPhiM( nr.Jet_pt_nom[goodJetIndex.at(k)], nr.Jet_eta[goodJetIndex.at(k)],
			     nr.Jet_phi[goodJetIndex.at(k)], nr.Jet_mass[goodJetIndex.at(k)] );

	  TLorentzVector tempVBF=tmp1+tmp2;

	  //require 2 jets be in opposite hemispheres
	  if ( nr.Jet_eta[goodJetIndex.at(j)] * nr.Jet_eta[goodJetIndex.at(k)] > 0 ) continue; 
	  if ( tempVBF.M() < VBF_MJJ_CUT ) continue;
	  if ( tempVBF.M() < tmpMassMax ) continue;
	  tmpMassMax = tempVBF.M();
	  vbf1=j; vbf2=k;
	}
      }
    
      if (vbf1==-1 && vbf2==-1) continue;
      
      WVJJTree->vbf1_AK4_pt = nr.Jet_pt_nom[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_eta = nr.Jet_eta[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_phi = nr.Jet_phi[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_m = nr.Jet_mass[goodJetIndex.at(vbf1)];
      WVJJTree->vbf1_AK4_qgid = nr.Jet_qgl[goodJetIndex.at(vbf1)];

      WVJJTree->vbf2_AK4_pt = nr.Jet_pt_nom[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_eta = nr.Jet_eta[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_phi = nr.Jet_phi[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_m = nr.Jet_mass[goodJetIndex.at(vbf2)];
      WVJJTree->vbf2_AK4_qgid = nr.Jet_qgl[goodJetIndex.at(vbf2)];

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

      if (isMC) {
      
      	WVJJTree->vbf1_AK4_pt_scaleUp = nr.Jet_pt_jesTotalUp[goodJetIndex.at(vbf1)];
      	WVJJTree->vbf1_AK4_pt_scaleDn = nr.Jet_pt_jesTotalDown[goodJetIndex.at(vbf1)];
      	WVJJTree->vbf1_AK4_m_scaleUp = nr.Jet_mass_jesTotalUp[goodJetIndex.at(vbf1)];
      	WVJJTree->vbf1_AK4_m_scaleDn = nr.Jet_mass_jesTotalDown[goodJetIndex.at(vbf1)];
      	
      	WVJJTree->vbf2_AK4_pt_scaleUp = nr.Jet_pt_jesTotalUp[goodJetIndex.at(vbf2)];
      	WVJJTree->vbf2_AK4_pt_scaleDn = nr.Jet_pt_jesTotalDown[goodJetIndex.at(vbf2)];
      	WVJJTree->vbf2_AK4_m_scaleUp = nr.Jet_mass_jesTotalUp[goodJetIndex.at(vbf2)];
      	WVJJTree->vbf2_AK4_m_scaleDn = nr.Jet_mass_jesTotalDown[goodJetIndex.at(vbf2)];
      	
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
      }

      TLorentzVector bosHad(0,0,0,0), bosHad_up(0,0,0,0), bosHad_dn(0,0,0,0);
      TLorentzVector bosLep(0,0,0,0), bosLep_up(0,0,0,0), bosLep_dn(0,0,0,0);
      
      //boosted event
      if (WVJJTree->bos_PuppiAK8_m_sd0_corr > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt, WVJJTree->bos_PuppiAK8_eta,
			    WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr);

	if (isMC) {
	  bosHad_up.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleUp, WVJJTree->bos_PuppiAK8_eta,
				 WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp);
	  bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleDn, WVJJTree->bos_PuppiAK8_eta,
				 WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn);
	}
      }
      //resolved event
      else if (WVJJTree->bos_AK4AK4_m > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt, WVJJTree->bos_AK4AK4_eta, 
			    WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m);

	if (isMC) {
	  bosHad_up.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleUp, WVJJTree->bos_AK4AK4_eta, 
				 WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleUp);
	  bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleDn, WVJJTree->bos_AK4AK4_eta, 
				 WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleDn);
	}
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
      
      	WVJJTree->nScaleWeight = *nr.nLHEScaleWeight;
      	WVJJTree->nPdfWeight = *nr.nLHEPdfWeight;
      	
      	for (uint j=0; j<WVJJTree->nScaleWeight; j++) {
      	  //LHE scale variation weights (w_var / w_nominal); [0] is MUR="0.5" MUF="0.5"; 
      	  //[1] is MUR="0.5" MUF="1.0"; [2] is MUR="0.5" MUF="2.0"; [3] is MUR="1.0" MUF="0.5"; 
      	  //[4] is MUR="1.0" MUF="2.0"; [5] is MUR="2.0" MUF="0.5"; [6] is MUR="2.0" MUF="1.0"; 
      	  //[7] is MUR="2.0" MUF="2.0"
      	  WVJJTree->scaleWeight[j]=nr.LHEScaleWeight[j];
      	}
      	for (uint j=0; j<WVJJTree->nPdfWeight; j++) {
      	  //LHE pdf variation weights (w_var / w_nominal) for LHA IDs 91400 - 91432
      	  WVJJTree->pdfWeight[j]=nr.LHEPdfWeight[j];
      	}
      	
      }
      
      if (isMC==1 && *nr.nLHEReweightingWeight!=0) {
      	WVJJTree->nAqgcWeight=*nr.nLHEReweightingWeight;
      
      	for (uint j=0; j<WVJJTree->nAqgcWeight; j++) {
      	  WVJJTree->aqgcWeight[j]=nr.LHEReweightingWeight[j];
      	}
      
      }
      
      if (isMC) {
	if (era==2016) {
	  WVJJTree->btagWeight = *nr.btagWeight_CMVA;
	}
	else {
        WVJJTree->btagWeight = *nr.btagWeight_DeepCSVB;
	}

	if (era!=2018) {
	  WVJJTree->L1PFWeight = *nr.L1PreFiringWeight_Nom;
	  WVJJTree->L1PFWeight_Up = *nr.L1PreFiringWeight_Up;
	  WVJJTree->L1PFWeight_Dn = *nr.L1PreFiringWeight_Dn;
	}
      }
      
      ot->Fill();
    }
    
    delete t;
    delete r;
    delete f;
    //t=0; 
    //r=0;
    //f=0;
  }

  of->Write();
  of->Close();
  
  return 0;

}
