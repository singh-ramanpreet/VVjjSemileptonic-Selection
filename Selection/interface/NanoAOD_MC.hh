//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 12 03:35:23 2020 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root
//////////////////////////////////////////////////////////

#ifndef NanoAOD_MC_h
#define NanoAOD_MC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class NanoAOD_MC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         btagWeight_CSVV2;
   Float_t         btagWeight_DeepCSVB;
   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Float_t         ChsMET_phi;
   Float_t         ChsMET_pt;
   Float_t         ChsMET_sumEt;
   UInt_t          nCorrT1METJet;
   Float_t         CorrT1METJet_area[29];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_eta[29];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_muonSubtrFactor[29];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_phi[29];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_rawPt[29];   //[nCorrT1METJet]
   UInt_t          nElectron;
   Float_t         Electron_deltaEtaSC[15];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[15];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[15];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[15];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[15];   //[nElectron]
   Float_t         Electron_dxy[15];   //[nElectron]
   Float_t         Electron_dxyErr[15];   //[nElectron]
   Float_t         Electron_dz[15];   //[nElectron]
   Float_t         Electron_dzErr[15];   //[nElectron]
   Float_t         Electron_eCorr[15];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[15];   //[nElectron]
   Float_t         Electron_energyErr[15];   //[nElectron]
   Float_t         Electron_eta[15];   //[nElectron]
   Float_t         Electron_hoe[15];   //[nElectron]
   Float_t         Electron_ip3d[15];   //[nElectron]
   Float_t         Electron_jetPtRelv2[15];   //[nElectron]
   Float_t         Electron_jetRelIso[15];   //[nElectron]
   Float_t         Electron_mass[15];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[15];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[15];   //[nElectron]
   Float_t         Electron_mvaFall17V1Iso[15];   //[nElectron]
   Float_t         Electron_mvaFall17V1noIso[15];   //[nElectron]
   Float_t         Electron_mvaFall17V2Iso[15];   //[nElectron]
   Float_t         Electron_mvaFall17V2noIso[15];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[15];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[15];   //[nElectron]
   Float_t         Electron_phi[15];   //[nElectron]
   Float_t         Electron_pt[15];   //[nElectron]
   Float_t         Electron_r9[15];   //[nElectron]
   Float_t         Electron_sieie[15];   //[nElectron]
   Float_t         Electron_sip3d[15];   //[nElectron]
   Float_t         Electron_mvaTTH[15];   //[nElectron]
   Int_t           Electron_charge[15];   //[nElectron]
   Int_t           Electron_cutBased[15];   //[nElectron]
   Int_t           Electron_cutBased_Fall17_V1[15];   //[nElectron]
   Int_t           Electron_jetIdx[15];   //[nElectron]
   Int_t           Electron_pdgId[15];   //[nElectron]
   Int_t           Electron_photonIdx[15];   //[nElectron]
   Int_t           Electron_tightCharge[15];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[15];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmapHEEP[15];   //[nElectron]
   Bool_t          Electron_convVeto[15];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[15];   //[nElectron]
   Bool_t          Electron_isPFcand[15];   //[nElectron]
   UChar_t         Electron_lostHits[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP80[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP90[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WPL[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP80[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP90[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WPL[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP80[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP90[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WPL[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP80[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP90[15];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WPL[15];   //[nElectron]
   UChar_t         Electron_seedGain[15];   //[nElectron]
   Bool_t          Flag_ecalBadCalibFilterV2;
   UInt_t          nFatJet;
   Float_t         FatJet_area[5];   //[nFatJet]
   Float_t         FatJet_btagCMVA[5];   //[nFatJet]
   Float_t         FatJet_btagCSVV2[5];   //[nFatJet]
   Float_t         FatJet_btagDDBvL[5];   //[nFatJet]
   Float_t         FatJet_btagDDBvL_noMD[5];   //[nFatJet]
   Float_t         FatJet_btagDDCvB[5];   //[nFatJet]
   Float_t         FatJet_btagDDCvB_noMD[5];   //[nFatJet]
   Float_t         FatJet_btagDDCvL[5];   //[nFatJet]
   Float_t         FatJet_btagDDCvL_noMD[5];   //[nFatJet]
   Float_t         FatJet_btagDeepB[5];   //[nFatJet]
   Float_t         FatJet_btagHbb[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_H4qvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_HbbvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_TvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_WvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHbbvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHccvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZbbvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_bbvsLight[5];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ccvsLight[5];   //[nFatJet]
   Float_t         FatJet_deepTag_H[5];   //[nFatJet]
   Float_t         FatJet_deepTag_QCD[5];   //[nFatJet]
   Float_t         FatJet_deepTag_QCDothers[5];   //[nFatJet]
   Float_t         FatJet_deepTag_TvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTag_WvsQCD[5];   //[nFatJet]
   Float_t         FatJet_deepTag_ZvsQCD[5];   //[nFatJet]
   Float_t         FatJet_eta[5];   //[nFatJet]
   Float_t         FatJet_mass[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop[5];   //[nFatJet]
   Float_t         FatJet_n2b1[5];   //[nFatJet]
   Float_t         FatJet_n3b1[5];   //[nFatJet]
   Float_t         FatJet_phi[5];   //[nFatJet]
   Float_t         FatJet_pt[5];   //[nFatJet]
   Float_t         FatJet_rawFactor[5];   //[nFatJet]
   Float_t         FatJet_tau1[5];   //[nFatJet]
   Float_t         FatJet_tau2[5];   //[nFatJet]
   Float_t         FatJet_tau3[5];   //[nFatJet]
   Float_t         FatJet_tau4[5];   //[nFatJet]
   Int_t           FatJet_jetId[5];   //[nFatJet]
   Int_t           FatJet_subJetIdx1[5];   //[nFatJet]
   Int_t           FatJet_subJetIdx2[5];   //[nFatJet]
   UInt_t          nFsrPhoton;
   Float_t         FsrPhoton_dROverEt2[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_eta[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_phi[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_pt[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_relIso03[3];   //[nFsrPhoton]
   Int_t           FsrPhoton_muonIdx[3];   //[nFsrPhoton]
   UInt_t          nGenJetAK8;
   Float_t         GenJetAK8_eta[15];   //[nGenJetAK8]
   Float_t         GenJetAK8_mass[15];   //[nGenJetAK8]
   Float_t         GenJetAK8_phi[15];   //[nGenJetAK8]
   Float_t         GenJetAK8_pt[15];   //[nGenJetAK8]
   UInt_t          nGenJet;
   Float_t         GenJet_eta[20];   //[nGenJet]
   Float_t         GenJet_mass[20];   //[nGenJet]
   Float_t         GenJet_phi[20];   //[nGenJet]
   Float_t         GenJet_pt[20];   //[nGenJet]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[130];   //[nGenPart]
   Float_t         GenPart_mass[130];   //[nGenPart]
   Float_t         GenPart_phi[130];   //[nGenPart]
   Float_t         GenPart_pt[130];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[130];   //[nGenPart]
   Int_t           GenPart_pdgId[130];   //[nGenPart]
   Int_t           GenPart_status[130];   //[nGenPart]
   Int_t           GenPart_statusFlags[130];   //[nGenPart]
   UInt_t          nSubGenJetAK8;
   Float_t         SubGenJetAK8_eta[13];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_mass[13];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_phi[13];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_pt[13];   //[nSubGenJetAK8]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   UInt_t          nGenVisTau;
   Float_t         GenVisTau_eta[3];   //[nGenVisTau]
   Float_t         GenVisTau_mass[3];   //[nGenVisTau]
   Float_t         GenVisTau_phi[3];   //[nGenVisTau]
   Float_t         GenVisTau_pt[3];   //[nGenVisTau]
   Int_t           GenVisTau_charge[3];   //[nGenVisTau]
   Int_t           GenVisTau_genPartIdxMother[3];   //[nGenVisTau]
   Int_t           GenVisTau_status[3];   //[nGenVisTau]
   Float_t         genWeight;
   Float_t         LHEWeight_originalXWGTUP;
   UInt_t          nLHEPdfWeight;
   Float_t         LHEPdfWeight[1];   //[nLHEPdfWeight]
   UInt_t          nLHEReweightingWeight;
   Float_t         LHEReweightingWeight[1];   //[nLHEReweightingWeight]
   UInt_t          nLHEScaleWeight;
   Float_t         LHEScaleWeight[1];   //[nLHEScaleWeight]
   UInt_t          nPSWeight;
   Float_t         PSWeight[4];   //[nPSWeight]
   UInt_t          nIsoTrack;
   Float_t         IsoTrack_dxy[5];   //[nIsoTrack]
   Float_t         IsoTrack_dz[5];   //[nIsoTrack]
   Float_t         IsoTrack_eta[5];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_all[5];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_chg[5];   //[nIsoTrack]
   Float_t         IsoTrack_phi[5];   //[nIsoTrack]
   Float_t         IsoTrack_pt[5];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_all[5];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_chg[5];   //[nIsoTrack]
   Int_t           IsoTrack_fromPV[5];   //[nIsoTrack]
   Int_t           IsoTrack_pdgId[5];   //[nIsoTrack]
   Bool_t          IsoTrack_isHighPurityTrack[5];   //[nIsoTrack]
   Bool_t          IsoTrack_isPFcand[5];   //[nIsoTrack]
   Bool_t          IsoTrack_isFromLostTrack[5];   //[nIsoTrack]
   UInt_t          nJet;
   Float_t         Jet_area[30];   //[nJet]
   Float_t         Jet_btagCMVA[30];   //[nJet]
   Float_t         Jet_btagCSVV2[30];   //[nJet]
   Float_t         Jet_btagDeepB[30];   //[nJet]
   Float_t         Jet_btagDeepC[30];   //[nJet]
   Float_t         Jet_btagDeepFlavB[30];   //[nJet]
   Float_t         Jet_btagDeepFlavC[30];   //[nJet]
   Float_t         Jet_chEmEF[30];   //[nJet]
   Float_t         Jet_chHEF[30];   //[nJet]
   Float_t         Jet_eta[30];   //[nJet]
   Float_t         Jet_jercCHF[30];   //[nJet]
   Float_t         Jet_jercCHPUF[30];   //[nJet]
   Float_t         Jet_mass[30];   //[nJet]
   Float_t         Jet_muEF[30];   //[nJet]
   Float_t         Jet_muonSubtrFactor[30];   //[nJet]
   Float_t         Jet_neEmEF[30];   //[nJet]
   Float_t         Jet_neHEF[30];   //[nJet]
   Float_t         Jet_phi[30];   //[nJet]
   Float_t         Jet_pt[30];   //[nJet]
   Float_t         Jet_qgl[30];   //[nJet]
   Float_t         Jet_rawFactor[30];   //[nJet]
   Float_t         Jet_bRegCorr[30];   //[nJet]
   Float_t         Jet_bRegRes[30];   //[nJet]
   Int_t           Jet_electronIdx1[30];   //[nJet]
   Int_t           Jet_electronIdx2[30];   //[nJet]
   Int_t           Jet_jetId[30];   //[nJet]
   Int_t           Jet_muonIdx1[30];   //[nJet]
   Int_t           Jet_muonIdx2[30];   //[nJet]
   Int_t           Jet_nConstituents[30];   //[nJet]
   Int_t           Jet_nElectrons[30];   //[nJet]
   Int_t           Jet_nMuons[30];   //[nJet]
   Int_t           Jet_puId[30];   //[nJet]
   Float_t         LHE_HT;
   Float_t         LHE_HTIncoming;
   Float_t         LHE_Vpt;
   UChar_t         LHE_Njets;
   UChar_t         LHE_Nb;
   UChar_t         LHE_Nc;
   UChar_t         LHE_Nuds;
   UChar_t         LHE_Nglu;
   UChar_t         LHE_NpNLO;
   UChar_t         LHE_NpLO;
   UInt_t          nLHEPart;
   Float_t         LHEPart_pt[6];   //[nLHEPart]
   Float_t         LHEPart_eta[6];   //[nLHEPart]
   Float_t         LHEPart_phi[6];   //[nLHEPart]
   Float_t         LHEPart_mass[6];   //[nLHEPart]
   Int_t           LHEPart_pdgId[6];   //[nLHEPart]
   Float_t         GenMET_phi;
   Float_t         GenMET_pt;
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_covXX;
   Float_t         MET_covXY;
   Float_t         MET_covYY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         MET_significance;
   Float_t         MET_sumEt;
   UInt_t          nMuon;
   Float_t         Muon_dxy[15];   //[nMuon]
   Float_t         Muon_dxyErr[15];   //[nMuon]
   Float_t         Muon_dz[15];   //[nMuon]
   Float_t         Muon_dzErr[15];   //[nMuon]
   Float_t         Muon_eta[15];   //[nMuon]
   Float_t         Muon_ip3d[15];   //[nMuon]
   Float_t         Muon_jetPtRelv2[15];   //[nMuon]
   Float_t         Muon_jetRelIso[15];   //[nMuon]
   Float_t         Muon_mass[15];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[15];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[15];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[15];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[15];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[15];   //[nMuon]
   Float_t         Muon_phi[15];   //[nMuon]
   Float_t         Muon_pt[15];   //[nMuon]
   Float_t         Muon_ptErr[15];   //[nMuon]
   Float_t         Muon_segmentComp[15];   //[nMuon]
   Float_t         Muon_sip3d[15];   //[nMuon]
   Float_t         Muon_softMva[15];   //[nMuon]
   Float_t         Muon_tkRelIso[15];   //[nMuon]
   Float_t         Muon_tunepRelPt[15];   //[nMuon]
   Float_t         Muon_mvaLowPt[15];   //[nMuon]
   Float_t         Muon_mvaTTH[15];   //[nMuon]
   Int_t           Muon_charge[15];   //[nMuon]
   Int_t           Muon_jetIdx[15];   //[nMuon]
   Int_t           Muon_nStations[15];   //[nMuon]
   Int_t           Muon_nTrackerLayers[15];   //[nMuon]
   Int_t           Muon_pdgId[15];   //[nMuon]
   Int_t           Muon_tightCharge[15];   //[nMuon]
   Int_t           Muon_fsrPhotonIdx[15];   //[nMuon]
   UChar_t         Muon_highPtId[15];   //[nMuon]
   Bool_t          Muon_inTimeMuon[15];   //[nMuon]
   Bool_t          Muon_isGlobal[15];   //[nMuon]
   Bool_t          Muon_isPFcand[15];   //[nMuon]
   Bool_t          Muon_isTracker[15];   //[nMuon]
   Bool_t          Muon_looseId[15];   //[nMuon]
   Bool_t          Muon_mediumId[15];   //[nMuon]
   Bool_t          Muon_mediumPromptId[15];   //[nMuon]
   UChar_t         Muon_miniIsoId[15];   //[nMuon]
   UChar_t         Muon_multiIsoId[15];   //[nMuon]
   UChar_t         Muon_mvaId[15];   //[nMuon]
   UChar_t         Muon_pfIsoId[15];   //[nMuon]
   Bool_t          Muon_softId[15];   //[nMuon]
   Bool_t          Muon_softMvaId[15];   //[nMuon]
   Bool_t          Muon_tightId[15];   //[nMuon]
   UChar_t         Muon_tkIsoId[15];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[15];   //[nMuon]
   UInt_t          nPhoton;
   Float_t         Photon_eCorr[10];   //[nPhoton]
   Float_t         Photon_energyErr[10];   //[nPhoton]
   Float_t         Photon_eta[10];   //[nPhoton]
   Float_t         Photon_hoe[10];   //[nPhoton]
   Float_t         Photon_mass[10];   //[nPhoton]
   Float_t         Photon_mvaID[10];   //[nPhoton]
   Float_t         Photon_mvaIDV1[10];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all[10];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg[10];   //[nPhoton]
   Float_t         Photon_phi[10];   //[nPhoton]
   Float_t         Photon_pt[10];   //[nPhoton]
   Float_t         Photon_r9[10];   //[nPhoton]
   Float_t         Photon_sieie[10];   //[nPhoton]
   Int_t           Photon_charge[10];   //[nPhoton]
   Int_t           Photon_cutBasedBitmap[10];   //[nPhoton]
   Int_t           Photon_cutBasedV1Bitmap[10];   //[nPhoton]
   Int_t           Photon_electronIdx[10];   //[nPhoton]
   Int_t           Photon_jetIdx[10];   //[nPhoton]
   Int_t           Photon_pdgId[10];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[10];   //[nPhoton]
   Bool_t          Photon_electronVeto[10];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[10];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[10];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[10];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[10];   //[nPhoton]
   Bool_t          Photon_pixelSeed[10];   //[nPhoton]
   UChar_t         Photon_seedGain[10];   //[nPhoton]
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nGenDressedLepton;
   Float_t         GenDressedLepton_eta[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_mass[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_phi[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_pt[4];   //[nGenDressedLepton]
   Int_t           GenDressedLepton_pdgId[4];   //[nGenDressedLepton]
   Bool_t          GenDressedLepton_hasTauAnc[4];   //[nGenDressedLepton]
   UInt_t          nSoftActivityJet;
   Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   UInt_t          nSubJet;
   Float_t         SubJet_btagCMVA[10];   //[nSubJet]
   Float_t         SubJet_btagCSVV2[10];   //[nSubJet]
   Float_t         SubJet_btagDeepB[10];   //[nSubJet]
   Float_t         SubJet_eta[10];   //[nSubJet]
   Float_t         SubJet_mass[10];   //[nSubJet]
   Float_t         SubJet_n2b1[10];   //[nSubJet]
   Float_t         SubJet_n3b1[10];   //[nSubJet]
   Float_t         SubJet_phi[10];   //[nSubJet]
   Float_t         SubJet_pt[10];   //[nSubJet]
   Float_t         SubJet_rawFactor[10];   //[nSubJet]
   Float_t         SubJet_tau1[10];   //[nSubJet]
   Float_t         SubJet_tau2[10];   //[nSubJet]
   Float_t         SubJet_tau3[10];   //[nSubJet]
   Float_t         SubJet_tau4[10];   //[nSubJet]
   UInt_t          nTau;
   Float_t         Tau_chargedIso[15];   //[nTau]
   Float_t         Tau_dxy[15];   //[nTau]
   Float_t         Tau_dz[15];   //[nTau]
   Float_t         Tau_eta[15];   //[nTau]
   Float_t         Tau_leadTkDeltaEta[15];   //[nTau]
   Float_t         Tau_leadTkDeltaPhi[15];   //[nTau]
   Float_t         Tau_leadTkPtOverTauPt[15];   //[nTau]
   Float_t         Tau_mass[15];   //[nTau]
   Float_t         Tau_neutralIso[15];   //[nTau]
   Float_t         Tau_phi[15];   //[nTau]
   Float_t         Tau_photonsOutsideSignalCone[15];   //[nTau]
   Float_t         Tau_pt[15];   //[nTau]
   Float_t         Tau_puCorr[15];   //[nTau]
   Float_t         Tau_rawAntiEle[15];   //[nTau]
   Float_t         Tau_rawAntiEle2018[15];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSe[15];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSjet[15];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSmu[15];   //[nTau]
   Float_t         Tau_rawIso[15];   //[nTau]
   Float_t         Tau_rawIsodR03[15];   //[nTau]
   Float_t         Tau_rawMVAnewDM2017v2[15];   //[nTau]
   Float_t         Tau_rawMVAoldDM[15];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v1[15];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v2[15];   //[nTau]
   Float_t         Tau_rawMVAoldDMdR032017v2[15];   //[nTau]
   Int_t           Tau_charge[15];   //[nTau]
   Int_t           Tau_decayMode[15];   //[nTau]
   Int_t           Tau_jetIdx[15];   //[nTau]
   Int_t           Tau_rawAntiEleCat[15];   //[nTau]
   Int_t           Tau_rawAntiEleCat2018[15];   //[nTau]
   UChar_t         Tau_idAntiEle[15];   //[nTau]
   UChar_t         Tau_idAntiEle2018[15];   //[nTau]
   UChar_t         Tau_idAntiMu[15];   //[nTau]
   Bool_t          Tau_idDecayMode[15];   //[nTau]
   Bool_t          Tau_idDecayModeNewDMs[15];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSe[15];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSjet[15];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSmu[15];   //[nTau]
   UChar_t         Tau_idMVAnewDM2017v2[15];   //[nTau]
   UChar_t         Tau_idMVAoldDM[15];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v1[15];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v2[15];   //[nTau]
   UChar_t         Tau_idMVAoldDMdR032017v2[15];   //[nTau]
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[36];   //[nTrigObj]
   Float_t         TrigObj_eta[36];   //[nTrigObj]
   Float_t         TrigObj_phi[36];   //[nTrigObj]
   Float_t         TrigObj_l1pt[36];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[36];   //[nTrigObj]
   Float_t         TrigObj_l2pt[36];   //[nTrigObj]
   Int_t           TrigObj_id[36];   //[nTrigObj]
   Int_t           TrigObj_l1iso[36];   //[nTrigObj]
   Int_t           TrigObj_l1charge[36];   //[nTrigObj]
   Int_t           TrigObj_filterBits[36];   //[nTrigObj]
   Int_t           genTtbarId;
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[11];   //[nSV]
   Float_t         SV_dlenSig[11];   //[nSV]
   Float_t         SV_dxy[11];   //[nSV]
   Float_t         SV_dxySig[11];   //[nSV]
   Float_t         SV_pAngle[11];   //[nSV]
   Int_t           Electron_genPartIdx[15];   //[nElectron]
   UChar_t         Electron_genPartFlav[15];   //[nElectron]
   Int_t           GenJetAK8_partonFlavour[15];   //[nGenJetAK8]
   UChar_t         GenJetAK8_hadronFlavour[15];   //[nGenJetAK8]
   Int_t           GenJet_partonFlavour[20];   //[nGenJet]
   UChar_t         GenJet_hadronFlavour[20];   //[nGenJet]
   Int_t           Jet_genJetIdx[30];   //[nJet]
   Int_t           Jet_hadronFlavour[30];   //[nJet]
   Int_t           Jet_partonFlavour[30];   //[nJet]
   Int_t           Muon_genPartIdx[15];   //[nMuon]
   UChar_t         Muon_genPartFlav[15];   //[nMuon]
   Int_t           Photon_genPartIdx[10];   //[nPhoton]
   UChar_t         Photon_genPartFlav[10];   //[nPhoton]
   Float_t         MET_fiducialGenPhi;
   Float_t         MET_fiducialGenPt;
   UChar_t         Electron_cleanmask[15];   //[nElectron]
   UChar_t         Jet_cleanmask[30];   //[nJet]
   UChar_t         Muon_cleanmask[15];   //[nMuon]
   UChar_t         Photon_cleanmask[10];   //[nPhoton]
   UChar_t         Tau_cleanmask[15];   //[nTau]
   Float_t         SV_chi2[11];   //[nSV]
   Float_t         SV_eta[11];   //[nSV]
   Float_t         SV_mass[11];   //[nSV]
   Float_t         SV_ndof[11];   //[nSV]
   Float_t         SV_phi[11];   //[nSV]
   Float_t         SV_pt[11];   //[nSV]
   Float_t         SV_x[11];   //[nSV]
   Float_t         SV_y[11];   //[nSV]
   Float_t         SV_z[11];   //[nSV]
   Int_t           Tau_genPartIdx[15];   //[nTau]
   UChar_t         Tau_genPartFlav[15];   //[nTau]
   Bool_t          L1simulation_step;
   Bool_t          HLTriggerFirstPath;
   Bool_t          HLT_DoubleEle25_CaloIdL_MW;
   Bool_t          HLT_DoubleEle27_CaloIdL_MW;
   Bool_t          HLT_DoubleEle33_CaloIdL_MW;
   Bool_t          HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Ele27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu37_Ele27_CaloIdL_MW;
   Bool_t          HLT_Mu37_TkMu27;
   Bool_t          HLT_DoubleMu4_3_Bs;
   Bool_t          HLT_DoubleMu4_3_Jpsi;
   Bool_t          HLT_DoubleMu4_JpsiTrk_Displaced;
   Bool_t          HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
   Bool_t          HLT_DoubleMu3_TkMu_DsTau3Mu;
   Bool_t          HLT_DoubleMu4_PsiPrimeTrk_Displaced;
   Bool_t          HLT_DoubleMu4_Mass3p8_DZ_PFHT350;
   Bool_t          HLT_Mu3_PFJet40;
   Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
   Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
   Bool_t          HLT_Mu7p5_Track2_Jpsi;
   Bool_t          HLT_Mu7p5_Track3p5_Jpsi;
   Bool_t          HLT_Mu7p5_Track7_Jpsi;
   Bool_t          HLT_Mu7p5_Track2_Upsilon;
   Bool_t          HLT_Mu7p5_Track3p5_Upsilon;
   Bool_t          HLT_Mu7p5_Track7_Upsilon;
   Bool_t          HLT_Mu3_L1SingleMu5orSingleMu7;
   Bool_t          HLT_DoublePhoton33_CaloIdL;
   Bool_t          HLT_DoublePhoton70;
   Bool_t          HLT_DoublePhoton85;
   Bool_t          HLT_Ele20_WPTight_Gsf;
   Bool_t          HLT_Ele15_WPLoose_Gsf;
   Bool_t          HLT_Ele17_WPLoose_Gsf;
   Bool_t          HLT_Ele20_WPLoose_Gsf;
   Bool_t          HLT_Ele20_eta2p1_WPLoose_Gsf;
   Bool_t          HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   Bool_t          HLT_Ele27_WPTight_Gsf;
   Bool_t          HLT_Ele28_WPTight_Gsf;
   Bool_t          HLT_Ele30_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf_L1EGMT;
   Bool_t          HLT_Ele38_WPTight_Gsf;
   Bool_t          HLT_Ele40_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t          HLT_IsoMu20;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_IsoMu30;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu25_TkMu0_Onia;
   Bool_t          HLT_Mu30_TkMu0_Psi;
   Bool_t          HLT_Mu30_TkMu0_Upsilon;
   Bool_t          HLT_Mu20_TkMu0_Phi;
   Bool_t          HLT_Mu25_TkMu0_Phi;
   Bool_t          HLT_Mu12;
   Bool_t          HLT_Mu15;
   Bool_t          HLT_Mu20;
   Bool_t          HLT_Mu27;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;
   Bool_t          HLT_OldMu100;
   Bool_t          HLT_TkMu100;
   Bool_t          HLT_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL;
   Bool_t          HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu8;
   Bool_t          HLT_Mu17;
   Bool_t          HLT_Mu19;
   Bool_t          HLT_Mu18_Mu9;
   Bool_t          HLT_Mu18_Mu9_DZ;
   Bool_t          HLT_Mu20_Mu10;
   Bool_t          HLT_Mu20_Mu10_DZ;
   Bool_t          HLT_Mu23_Mu12;
   Bool_t          HLT_Mu23_Mu12_DZ;
   Bool_t          HLTriggerFinalPath;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   Float_t         Jet_pt_raw[30];   //[nJet]
   Float_t         Jet_pt_nom[30];   //[nJet]
   Float_t         Jet_mass_raw[30];   //[nJet]
   Float_t         Jet_mass_nom[30];   //[nJet]
   Float_t         Jet_corr_JEC[30];   //[nJet]
   Float_t         Jet_corr_JER[30];   //[nJet]
   Float_t         MET_pt_nom;
   Float_t         MET_phi_nom;
   Float_t         MET_pt_jer;
   Float_t         MET_phi_jer;
   Float_t         Jet_pt_jerUp[30];   //[nJet]
   Float_t         Jet_mass_jerUp[30];   //[nJet]
   Float_t         MET_pt_jerUp;
   Float_t         MET_phi_jerUp;
   Float_t         Jet_pt_jesAbsoluteStatUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteStatUp[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteStatUp;
   Float_t         MET_phi_jesAbsoluteStatUp;
   Float_t         Jet_pt_jesAbsoluteScaleUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteScaleUp[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteScaleUp;
   Float_t         MET_phi_jesAbsoluteScaleUp;
   Float_t         Jet_pt_jesAbsoluteSampleUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteSampleUp[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteSampleUp;
   Float_t         MET_phi_jesAbsoluteSampleUp;
   Float_t         Jet_pt_jesAbsoluteFlavMapUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteFlavMapUp[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteFlavMapUp;
   Float_t         MET_phi_jesAbsoluteFlavMapUp;
   Float_t         Jet_pt_jesAbsoluteMPFBiasUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteMPFBiasUp[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteMPFBiasUp;
   Float_t         MET_phi_jesAbsoluteMPFBiasUp;
   Float_t         Jet_pt_jesFragmentationUp[30];   //[nJet]
   Float_t         Jet_mass_jesFragmentationUp[30];   //[nJet]
   Float_t         MET_pt_jesFragmentationUp;
   Float_t         MET_phi_jesFragmentationUp;
   Float_t         Jet_pt_jesSinglePionECALUp[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionECALUp[30];   //[nJet]
   Float_t         MET_pt_jesSinglePionECALUp;
   Float_t         MET_phi_jesSinglePionECALUp;
   Float_t         Jet_pt_jesSinglePionHCALUp[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionHCALUp[30];   //[nJet]
   Float_t         MET_pt_jesSinglePionHCALUp;
   Float_t         MET_phi_jesSinglePionHCALUp;
   Float_t         Jet_pt_jesFlavorQCDUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorQCDUp[30];   //[nJet]
   Float_t         MET_pt_jesFlavorQCDUp;
   Float_t         MET_phi_jesFlavorQCDUp;
   Float_t         Jet_pt_jesTimePtEtaUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimePtEtaUp[30];   //[nJet]
   Float_t         MET_pt_jesTimePtEtaUp;
   Float_t         MET_phi_jesTimePtEtaUp;
   Float_t         Jet_pt_jesRelativeJEREC1Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC1Up[30];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC1Up;
   Float_t         MET_phi_jesRelativeJEREC1Up;
   Float_t         Jet_pt_jesRelativeJEREC2Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC2Up[30];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC2Up;
   Float_t         MET_phi_jesRelativeJEREC2Up;
   Float_t         Jet_pt_jesRelativeJERHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJERHFUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativeJERHFUp;
   Float_t         MET_phi_jesRelativeJERHFUp;
   Float_t         Jet_pt_jesRelativePtBBUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtBBUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtBBUp;
   Float_t         MET_phi_jesRelativePtBBUp;
   Float_t         Jet_pt_jesRelativePtEC1Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC1Up[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC1Up;
   Float_t         MET_phi_jesRelativePtEC1Up;
   Float_t         Jet_pt_jesRelativePtEC2Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC2Up[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC2Up;
   Float_t         MET_phi_jesRelativePtEC2Up;
   Float_t         Jet_pt_jesRelativePtHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtHFUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtHFUp;
   Float_t         MET_phi_jesRelativePtHFUp;
   Float_t         Jet_pt_jesRelativeBalUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeBalUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativeBalUp;
   Float_t         MET_phi_jesRelativeBalUp;
   Float_t         Jet_pt_jesRelativeSampleUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeSampleUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativeSampleUp;
   Float_t         MET_phi_jesRelativeSampleUp;
   Float_t         Jet_pt_jesRelativeFSRUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeFSRUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativeFSRUp;
   Float_t         MET_phi_jesRelativeFSRUp;
   Float_t         Jet_pt_jesRelativeStatFSRUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatFSRUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativeStatFSRUp;
   Float_t         MET_phi_jesRelativeStatFSRUp;
   Float_t         Jet_pt_jesRelativeStatECUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatECUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativeStatECUp;
   Float_t         MET_phi_jesRelativeStatECUp;
   Float_t         Jet_pt_jesRelativeStatHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatHFUp[30];   //[nJet]
   Float_t         MET_pt_jesRelativeStatHFUp;
   Float_t         MET_phi_jesRelativeStatHFUp;
   Float_t         Jet_pt_jesPileUpDataMCUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpDataMCUp[30];   //[nJet]
   Float_t         MET_pt_jesPileUpDataMCUp;
   Float_t         MET_phi_jesPileUpDataMCUp;
   Float_t         Jet_pt_jesPileUpPtRefUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtRefUp[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtRefUp;
   Float_t         MET_phi_jesPileUpPtRefUp;
   Float_t         Jet_pt_jesPileUpPtBBUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtBBUp[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtBBUp;
   Float_t         MET_phi_jesPileUpPtBBUp;
   Float_t         Jet_pt_jesPileUpPtEC1Up[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC1Up[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC1Up;
   Float_t         MET_phi_jesPileUpPtEC1Up;
   Float_t         Jet_pt_jesPileUpPtEC2Up[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC2Up[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC2Up;
   Float_t         MET_phi_jesPileUpPtEC2Up;
   Float_t         Jet_pt_jesPileUpPtHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtHFUp[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtHFUp;
   Float_t         MET_phi_jesPileUpPtHFUp;
   Float_t         Jet_pt_jesPileUpMuZeroUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpMuZeroUp[30];   //[nJet]
   Float_t         MET_pt_jesPileUpMuZeroUp;
   Float_t         MET_phi_jesPileUpMuZeroUp;
   Float_t         Jet_pt_jesPileUpEnvelopeUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpEnvelopeUp[30];   //[nJet]
   Float_t         MET_pt_jesPileUpEnvelopeUp;
   Float_t         MET_phi_jesPileUpEnvelopeUp;
   Float_t         Jet_pt_jesSubTotalPileUpUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPileUpUp[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalPileUpUp;
   Float_t         MET_phi_jesSubTotalPileUpUp;
   Float_t         Jet_pt_jesSubTotalRelativeUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalRelativeUp[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalRelativeUp;
   Float_t         MET_phi_jesSubTotalRelativeUp;
   Float_t         Jet_pt_jesSubTotalPtUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPtUp[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalPtUp;
   Float_t         MET_phi_jesSubTotalPtUp;
   Float_t         Jet_pt_jesSubTotalScaleUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalScaleUp[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalScaleUp;
   Float_t         MET_phi_jesSubTotalScaleUp;
   Float_t         Jet_pt_jesSubTotalAbsoluteUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalAbsoluteUp[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalAbsoluteUp;
   Float_t         MET_phi_jesSubTotalAbsoluteUp;
   Float_t         Jet_pt_jesSubTotalMCUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalMCUp[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalMCUp;
   Float_t         MET_phi_jesSubTotalMCUp;
   Float_t         Jet_pt_jesTotalUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalUp[30];   //[nJet]
   Float_t         MET_pt_jesTotalUp;
   Float_t         MET_phi_jesTotalUp;
   Float_t         Jet_pt_jesTotalNoFlavorUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorUp[30];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorUp;
   Float_t         MET_phi_jesTotalNoFlavorUp;
   Float_t         Jet_pt_jesTotalNoTimeUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoTimeUp[30];   //[nJet]
   Float_t         MET_pt_jesTotalNoTimeUp;
   Float_t         MET_phi_jesTotalNoTimeUp;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeUp[30];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorNoTimeUp;
   Float_t         MET_phi_jesTotalNoFlavorNoTimeUp;
   Float_t         Jet_pt_jesFlavorZJetUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorZJetUp[30];   //[nJet]
   Float_t         MET_pt_jesFlavorZJetUp;
   Float_t         MET_phi_jesFlavorZJetUp;
   Float_t         Jet_pt_jesFlavorPhotonJetUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPhotonJetUp[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPhotonJetUp;
   Float_t         MET_phi_jesFlavorPhotonJetUp;
   Float_t         Jet_pt_jesFlavorPureGluonUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureGluonUp[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureGluonUp;
   Float_t         MET_phi_jesFlavorPureGluonUp;
   Float_t         Jet_pt_jesFlavorPureQuarkUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureQuarkUp[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureQuarkUp;
   Float_t         MET_phi_jesFlavorPureQuarkUp;
   Float_t         Jet_pt_jesFlavorPureCharmUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureCharmUp[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureCharmUp;
   Float_t         MET_phi_jesFlavorPureCharmUp;
   Float_t         Jet_pt_jesFlavorPureBottomUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureBottomUp[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureBottomUp;
   Float_t         MET_phi_jesFlavorPureBottomUp;
   Float_t         Jet_pt_jesTimeRunAUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunAUp[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunAUp;
   Float_t         MET_phi_jesTimeRunAUp;
   Float_t         Jet_pt_jesTimeRunBUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunBUp[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunBUp;
   Float_t         MET_phi_jesTimeRunBUp;
   Float_t         Jet_pt_jesTimeRunCUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunCUp[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunCUp;
   Float_t         MET_phi_jesTimeRunCUp;
   Float_t         Jet_pt_jesTimeRunDUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunDUp[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunDUp;
   Float_t         MET_phi_jesTimeRunDUp;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituUp[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupMPFInSituUp;
   Float_t         MET_phi_jesCorrelationGroupMPFInSituUp;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationUp[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupIntercalibrationUp;
   Float_t         MET_phi_jesCorrelationGroupIntercalibrationUp;
   Float_t         Jet_pt_jesCorrelationGroupbJESUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupbJESUp[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupbJESUp;
   Float_t         MET_phi_jesCorrelationGroupbJESUp;
   Float_t         Jet_pt_jesCorrelationGroupFlavorUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupFlavorUp[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupFlavorUp;
   Float_t         MET_phi_jesCorrelationGroupFlavorUp;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedUp[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_phi_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_pt_unclustEnUp;
   Float_t         MET_phi_unclustEnUp;
   Float_t         Jet_pt_jerDown[30];   //[nJet]
   Float_t         Jet_mass_jerDown[30];   //[nJet]
   Float_t         MET_pt_jerDown;
   Float_t         MET_phi_jerDown;
   Float_t         Jet_pt_jesAbsoluteStatDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteStatDown[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteStatDown;
   Float_t         MET_phi_jesAbsoluteStatDown;
   Float_t         Jet_pt_jesAbsoluteScaleDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteScaleDown[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteScaleDown;
   Float_t         MET_phi_jesAbsoluteScaleDown;
   Float_t         Jet_pt_jesAbsoluteSampleDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteSampleDown[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteSampleDown;
   Float_t         MET_phi_jesAbsoluteSampleDown;
   Float_t         Jet_pt_jesAbsoluteFlavMapDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteFlavMapDown[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteFlavMapDown;
   Float_t         MET_phi_jesAbsoluteFlavMapDown;
   Float_t         Jet_pt_jesAbsoluteMPFBiasDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteMPFBiasDown[30];   //[nJet]
   Float_t         MET_pt_jesAbsoluteMPFBiasDown;
   Float_t         MET_phi_jesAbsoluteMPFBiasDown;
   Float_t         Jet_pt_jesFragmentationDown[30];   //[nJet]
   Float_t         Jet_mass_jesFragmentationDown[30];   //[nJet]
   Float_t         MET_pt_jesFragmentationDown;
   Float_t         MET_phi_jesFragmentationDown;
   Float_t         Jet_pt_jesSinglePionECALDown[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionECALDown[30];   //[nJet]
   Float_t         MET_pt_jesSinglePionECALDown;
   Float_t         MET_phi_jesSinglePionECALDown;
   Float_t         Jet_pt_jesSinglePionHCALDown[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionHCALDown[30];   //[nJet]
   Float_t         MET_pt_jesSinglePionHCALDown;
   Float_t         MET_phi_jesSinglePionHCALDown;
   Float_t         Jet_pt_jesFlavorQCDDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorQCDDown[30];   //[nJet]
   Float_t         MET_pt_jesFlavorQCDDown;
   Float_t         MET_phi_jesFlavorQCDDown;
   Float_t         Jet_pt_jesTimePtEtaDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimePtEtaDown[30];   //[nJet]
   Float_t         MET_pt_jesTimePtEtaDown;
   Float_t         MET_phi_jesTimePtEtaDown;
   Float_t         Jet_pt_jesRelativeJEREC1Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC1Down[30];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC1Down;
   Float_t         MET_phi_jesRelativeJEREC1Down;
   Float_t         Jet_pt_jesRelativeJEREC2Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC2Down[30];   //[nJet]
   Float_t         MET_pt_jesRelativeJEREC2Down;
   Float_t         MET_phi_jesRelativeJEREC2Down;
   Float_t         Jet_pt_jesRelativeJERHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJERHFDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativeJERHFDown;
   Float_t         MET_phi_jesRelativeJERHFDown;
   Float_t         Jet_pt_jesRelativePtBBDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtBBDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtBBDown;
   Float_t         MET_phi_jesRelativePtBBDown;
   Float_t         Jet_pt_jesRelativePtEC1Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC1Down[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC1Down;
   Float_t         MET_phi_jesRelativePtEC1Down;
   Float_t         Jet_pt_jesRelativePtEC2Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC2Down[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtEC2Down;
   Float_t         MET_phi_jesRelativePtEC2Down;
   Float_t         Jet_pt_jesRelativePtHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtHFDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativePtHFDown;
   Float_t         MET_phi_jesRelativePtHFDown;
   Float_t         Jet_pt_jesRelativeBalDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeBalDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativeBalDown;
   Float_t         MET_phi_jesRelativeBalDown;
   Float_t         Jet_pt_jesRelativeSampleDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeSampleDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativeSampleDown;
   Float_t         MET_phi_jesRelativeSampleDown;
   Float_t         Jet_pt_jesRelativeFSRDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeFSRDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativeFSRDown;
   Float_t         MET_phi_jesRelativeFSRDown;
   Float_t         Jet_pt_jesRelativeStatFSRDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatFSRDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativeStatFSRDown;
   Float_t         MET_phi_jesRelativeStatFSRDown;
   Float_t         Jet_pt_jesRelativeStatECDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatECDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativeStatECDown;
   Float_t         MET_phi_jesRelativeStatECDown;
   Float_t         Jet_pt_jesRelativeStatHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatHFDown[30];   //[nJet]
   Float_t         MET_pt_jesRelativeStatHFDown;
   Float_t         MET_phi_jesRelativeStatHFDown;
   Float_t         Jet_pt_jesPileUpDataMCDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpDataMCDown[30];   //[nJet]
   Float_t         MET_pt_jesPileUpDataMCDown;
   Float_t         MET_phi_jesPileUpDataMCDown;
   Float_t         Jet_pt_jesPileUpPtRefDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtRefDown[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtRefDown;
   Float_t         MET_phi_jesPileUpPtRefDown;
   Float_t         Jet_pt_jesPileUpPtBBDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtBBDown[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtBBDown;
   Float_t         MET_phi_jesPileUpPtBBDown;
   Float_t         Jet_pt_jesPileUpPtEC1Down[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC1Down[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC1Down;
   Float_t         MET_phi_jesPileUpPtEC1Down;
   Float_t         Jet_pt_jesPileUpPtEC2Down[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC2Down[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtEC2Down;
   Float_t         MET_phi_jesPileUpPtEC2Down;
   Float_t         Jet_pt_jesPileUpPtHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtHFDown[30];   //[nJet]
   Float_t         MET_pt_jesPileUpPtHFDown;
   Float_t         MET_phi_jesPileUpPtHFDown;
   Float_t         Jet_pt_jesPileUpMuZeroDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpMuZeroDown[30];   //[nJet]
   Float_t         MET_pt_jesPileUpMuZeroDown;
   Float_t         MET_phi_jesPileUpMuZeroDown;
   Float_t         Jet_pt_jesPileUpEnvelopeDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpEnvelopeDown[30];   //[nJet]
   Float_t         MET_pt_jesPileUpEnvelopeDown;
   Float_t         MET_phi_jesPileUpEnvelopeDown;
   Float_t         Jet_pt_jesSubTotalPileUpDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPileUpDown[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalPileUpDown;
   Float_t         MET_phi_jesSubTotalPileUpDown;
   Float_t         Jet_pt_jesSubTotalRelativeDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalRelativeDown[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalRelativeDown;
   Float_t         MET_phi_jesSubTotalRelativeDown;
   Float_t         Jet_pt_jesSubTotalPtDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPtDown[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalPtDown;
   Float_t         MET_phi_jesSubTotalPtDown;
   Float_t         Jet_pt_jesSubTotalScaleDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalScaleDown[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalScaleDown;
   Float_t         MET_phi_jesSubTotalScaleDown;
   Float_t         Jet_pt_jesSubTotalAbsoluteDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalAbsoluteDown[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalAbsoluteDown;
   Float_t         MET_phi_jesSubTotalAbsoluteDown;
   Float_t         Jet_pt_jesSubTotalMCDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalMCDown[30];   //[nJet]
   Float_t         MET_pt_jesSubTotalMCDown;
   Float_t         MET_phi_jesSubTotalMCDown;
   Float_t         Jet_pt_jesTotalDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalDown[30];   //[nJet]
   Float_t         MET_pt_jesTotalDown;
   Float_t         MET_phi_jesTotalDown;
   Float_t         Jet_pt_jesTotalNoFlavorDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorDown[30];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorDown;
   Float_t         MET_phi_jesTotalNoFlavorDown;
   Float_t         Jet_pt_jesTotalNoTimeDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoTimeDown[30];   //[nJet]
   Float_t         MET_pt_jesTotalNoTimeDown;
   Float_t         MET_phi_jesTotalNoTimeDown;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeDown[30];   //[nJet]
   Float_t         MET_pt_jesTotalNoFlavorNoTimeDown;
   Float_t         MET_phi_jesTotalNoFlavorNoTimeDown;
   Float_t         Jet_pt_jesFlavorZJetDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorZJetDown[30];   //[nJet]
   Float_t         MET_pt_jesFlavorZJetDown;
   Float_t         MET_phi_jesFlavorZJetDown;
   Float_t         Jet_pt_jesFlavorPhotonJetDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPhotonJetDown[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPhotonJetDown;
   Float_t         MET_phi_jesFlavorPhotonJetDown;
   Float_t         Jet_pt_jesFlavorPureGluonDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureGluonDown[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureGluonDown;
   Float_t         MET_phi_jesFlavorPureGluonDown;
   Float_t         Jet_pt_jesFlavorPureQuarkDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureQuarkDown[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureQuarkDown;
   Float_t         MET_phi_jesFlavorPureQuarkDown;
   Float_t         Jet_pt_jesFlavorPureCharmDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureCharmDown[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureCharmDown;
   Float_t         MET_phi_jesFlavorPureCharmDown;
   Float_t         Jet_pt_jesFlavorPureBottomDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureBottomDown[30];   //[nJet]
   Float_t         MET_pt_jesFlavorPureBottomDown;
   Float_t         MET_phi_jesFlavorPureBottomDown;
   Float_t         Jet_pt_jesTimeRunADown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunADown[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunADown;
   Float_t         MET_phi_jesTimeRunADown;
   Float_t         Jet_pt_jesTimeRunBDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunBDown[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunBDown;
   Float_t         MET_phi_jesTimeRunBDown;
   Float_t         Jet_pt_jesTimeRunCDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunCDown[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunCDown;
   Float_t         MET_phi_jesTimeRunCDown;
   Float_t         Jet_pt_jesTimeRunDDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunDDown[30];   //[nJet]
   Float_t         MET_pt_jesTimeRunDDown;
   Float_t         MET_phi_jesTimeRunDDown;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituDown[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupMPFInSituDown;
   Float_t         MET_phi_jesCorrelationGroupMPFInSituDown;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationDown[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupIntercalibrationDown;
   Float_t         MET_phi_jesCorrelationGroupIntercalibrationDown;
   Float_t         Jet_pt_jesCorrelationGroupbJESDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupbJESDown[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupbJESDown;
   Float_t         MET_phi_jesCorrelationGroupbJESDown;
   Float_t         Jet_pt_jesCorrelationGroupFlavorDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupFlavorDown[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupFlavorDown;
   Float_t         MET_phi_jesCorrelationGroupFlavorDown;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedDown[30];   //[nJet]
   Float_t         MET_pt_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_phi_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_pt_unclustEnDown;
   Float_t         MET_phi_unclustEnDown;
   Float_t         FatJet_pt_raw[5];   //[nFatJet]
   Float_t         FatJet_pt_nom[5];   //[nFatJet]
   Float_t         FatJet_mass_raw[5];   //[nFatJet]
   Float_t         FatJet_mass_nom[5];   //[nFatJet]
   Float_t         FatJet_corr_JEC[5];   //[nFatJet]
   Float_t         FatJet_corr_JER[5];   //[nFatJet]
   Float_t         FatJet_corr_JMS[5];   //[nFatJet]
   Float_t         FatJet_corr_JMR[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_raw[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_nom[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_corr_JMR[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_corr_JMS[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_corr_PUPPI[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_nom[5];   //[nFatJet]
   Float_t         FatJet_pt_jerUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jerUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jmrUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jmsUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jerUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmrUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmsUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_jerUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_jmrUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_jmsUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteStatUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteStatUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteStatUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteScaleUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteScaleUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteScaleUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteSampleUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteSampleUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteSampleUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteFlavMapUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteFlavMapUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteFlavMapUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteMPFBiasUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteMPFBiasUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteMPFBiasUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFragmentationUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFragmentationUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFragmentationUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionECALUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionECALUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionECALUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionHCALUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionHCALUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionHCALUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorQCDUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorQCDUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorQCDUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimePtEtaUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimePtEtaUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimePtEtaUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC1Up[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC1Up[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC1Up[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC2Up[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC2Up[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC2Up[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJERHFUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJERHFUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJERHFUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtBBUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtBBUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtBBUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC1Up[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC1Up[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC1Up[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC2Up[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC2Up[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC2Up[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtHFUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtHFUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtHFUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeBalUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeBalUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeBalUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeSampleUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeSampleUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeSampleUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeFSRUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeFSRUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeFSRUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatFSRUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatFSRUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatFSRUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatECUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatECUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatECUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatHFUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatHFUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatHFUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpDataMCUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpDataMCUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpDataMCUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtRefUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtRefUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtRefUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtBBUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtBBUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtBBUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC1Up[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC1Up[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC1Up[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC2Up[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC2Up[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC2Up[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtHFUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtHFUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtHFUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpMuZeroUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpMuZeroUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpMuZeroUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpEnvelopeUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpEnvelopeUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpEnvelopeUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPileUpUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPileUpUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPileUpUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalRelativeUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalRelativeUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalRelativeUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPtUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPtUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPtUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalScaleUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalScaleUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalScaleUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalAbsoluteUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalAbsoluteUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalAbsoluteUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalMCUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalMCUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalMCUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoTimeUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoTimeUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoTimeUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorNoTimeUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorNoTimeUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorZJetUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorZJetUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorZJetUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPhotonJetUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPhotonJetUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPhotonJetUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureGluonUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureGluonUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureGluonUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureQuarkUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureQuarkUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureQuarkUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureCharmUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureCharmUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureCharmUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureBottomUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureBottomUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureBottomUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunAUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunAUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunAUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunBUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunBUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunBUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunCUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunCUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunCUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunDUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunDUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunDUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupMPFInSituUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupMPFInSituUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupIntercalibrationUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupIntercalibrationUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupbJESUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupbJESUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupbJESUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupFlavorUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupFlavorUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupFlavorUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupUncorrelatedUp[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupUncorrelatedUp[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp[5];   //[nFatJet]
   Float_t         FatJet_pt_jerDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jerDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jmrDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jmsDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jerDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmrDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jmsDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_jerDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_jmrDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_tau21DDT_jmsDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteStatDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteStatDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteStatDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteScaleDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteScaleDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteScaleDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteSampleDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteSampleDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteSampleDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteFlavMapDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteFlavMapDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteFlavMapDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesAbsoluteMPFBiasDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesAbsoluteMPFBiasDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesAbsoluteMPFBiasDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFragmentationDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFragmentationDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFragmentationDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionECALDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionECALDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionECALDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSinglePionHCALDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSinglePionHCALDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSinglePionHCALDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorQCDDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorQCDDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorQCDDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimePtEtaDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimePtEtaDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimePtEtaDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC1Down[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC1Down[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC1Down[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJEREC2Down[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJEREC2Down[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJEREC2Down[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeJERHFDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeJERHFDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeJERHFDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtBBDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtBBDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtBBDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC1Down[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC1Down[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC1Down[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtEC2Down[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtEC2Down[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtEC2Down[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativePtHFDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativePtHFDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativePtHFDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeBalDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeBalDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeBalDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeSampleDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeSampleDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeSampleDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeFSRDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeFSRDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeFSRDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatFSRDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatFSRDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatFSRDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatECDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatECDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatECDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesRelativeStatHFDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesRelativeStatHFDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesRelativeStatHFDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpDataMCDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpDataMCDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpDataMCDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtRefDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtRefDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtRefDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtBBDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtBBDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtBBDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC1Down[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC1Down[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC1Down[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtEC2Down[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtEC2Down[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtEC2Down[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpPtHFDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpPtHFDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpPtHFDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpMuZeroDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpMuZeroDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpMuZeroDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesPileUpEnvelopeDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesPileUpEnvelopeDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesPileUpEnvelopeDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPileUpDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPileUpDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPileUpDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalRelativeDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalRelativeDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalRelativeDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalPtDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalPtDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalPtDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalScaleDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalScaleDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalScaleDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalAbsoluteDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalAbsoluteDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalAbsoluteDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesSubTotalMCDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesSubTotalMCDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesSubTotalMCDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoTimeDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoTimeDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoTimeDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTotalNoFlavorNoTimeDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTotalNoFlavorNoTimeDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorZJetDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorZJetDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorZJetDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPhotonJetDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPhotonJetDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPhotonJetDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureGluonDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureGluonDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureGluonDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureQuarkDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureQuarkDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureQuarkDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureCharmDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureCharmDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureCharmDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesFlavorPureBottomDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesFlavorPureBottomDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesFlavorPureBottomDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunADown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunADown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunADown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunBDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunBDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunBDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunCDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunCDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunCDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesTimeRunDDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesTimeRunDDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesTimeRunDDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupMPFInSituDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupMPFInSituDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupIntercalibrationDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupIntercalibrationDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupbJESDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupbJESDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupbJESDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupFlavorDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupFlavorDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupFlavorDown[5];   //[nFatJet]
   Float_t         FatJet_pt_jesCorrelationGroupUncorrelatedDown[5];   //[nFatJet]
   Float_t         FatJet_mass_jesCorrelationGroupUncorrelatedDown[5];   //[nFatJet]
   Float_t         FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown[5];   //[nFatJet]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_btagWeight_CSVV2;   //!
   TBranch        *b_btagWeight_DeepCSVB;   //!
   TBranch        *b_CaloMET_phi;   //!
   TBranch        *b_CaloMET_pt;   //!
   TBranch        *b_CaloMET_sumEt;   //!
   TBranch        *b_ChsMET_phi;   //!
   TBranch        *b_ChsMET_pt;   //!
   TBranch        *b_ChsMET_sumEt;   //!
   TBranch        *b_nCorrT1METJet;   //!
   TBranch        *b_CorrT1METJet_area;   //!
   TBranch        *b_CorrT1METJet_eta;   //!
   TBranch        *b_CorrT1METJet_muonSubtrFactor;   //!
   TBranch        *b_CorrT1METJet_phi;   //!
   TBranch        *b_CorrT1METJet_rawPt;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eCorr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetPtRelv2;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_Electron_mvaFall17V1Iso;   //!
   TBranch        *b_Electron_mvaFall17V1noIso;   //!
   TBranch        *b_Electron_mvaFall17V2Iso;   //!
   TBranch        *b_Electron_mvaFall17V2noIso;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_mvaTTH;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_cutBased_Fall17_V1;   //!
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_tightCharge;   //!
   TBranch        *b_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_Electron_vidNestedWPBitmapHEEP;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_lostHits;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WPL;   //!
   TBranch        *b_Electron_seedGain;   //!
   TBranch        *b_Flag_ecalBadCalibFilterV2;   //!
   TBranch        *b_nFatJet;   //!
   TBranch        *b_FatJet_area;   //!
   TBranch        *b_FatJet_btagCMVA;   //!
   TBranch        *b_FatJet_btagCSVV2;   //!
   TBranch        *b_FatJet_btagDDBvL;   //!
   TBranch        *b_FatJet_btagDDBvL_noMD;   //!
   TBranch        *b_FatJet_btagDDCvB;   //!
   TBranch        *b_FatJet_btagDDCvB_noMD;   //!
   TBranch        *b_FatJet_btagDDCvL;   //!
   TBranch        *b_FatJet_btagDDCvL_noMD;   //!
   TBranch        *b_FatJet_btagDeepB;   //!
   TBranch        *b_FatJet_btagHbb;   //!
   TBranch        *b_FatJet_deepTagMD_H4qvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_HbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_TvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_WvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHccvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_bbvsLight;   //!
   TBranch        *b_FatJet_deepTagMD_ccvsLight;   //!
   TBranch        *b_FatJet_deepTag_H;   //!
   TBranch        *b_FatJet_deepTag_QCD;   //!
   TBranch        *b_FatJet_deepTag_QCDothers;   //!
   TBranch        *b_FatJet_deepTag_TvsQCD;   //!
   TBranch        *b_FatJet_deepTag_WvsQCD;   //!
   TBranch        *b_FatJet_deepTag_ZvsQCD;   //!
   TBranch        *b_FatJet_eta;   //!
   TBranch        *b_FatJet_mass;   //!
   TBranch        *b_FatJet_msoftdrop;   //!
   TBranch        *b_FatJet_n2b1;   //!
   TBranch        *b_FatJet_n3b1;   //!
   TBranch        *b_FatJet_phi;   //!
   TBranch        *b_FatJet_pt;   //!
   TBranch        *b_FatJet_rawFactor;   //!
   TBranch        *b_FatJet_tau1;   //!
   TBranch        *b_FatJet_tau2;   //!
   TBranch        *b_FatJet_tau3;   //!
   TBranch        *b_FatJet_tau4;   //!
   TBranch        *b_FatJet_jetId;   //!
   TBranch        *b_FatJet_subJetIdx1;   //!
   TBranch        *b_FatJet_subJetIdx2;   //!
   TBranch        *b_nFsrPhoton;   //!
   TBranch        *b_FsrPhoton_dROverEt2;   //!
   TBranch        *b_FsrPhoton_eta;   //!
   TBranch        *b_FsrPhoton_phi;   //!
   TBranch        *b_FsrPhoton_pt;   //!
   TBranch        *b_FsrPhoton_relIso03;   //!
   TBranch        *b_FsrPhoton_muonIdx;   //!
   TBranch        *b_nGenJetAK8;   //!
   TBranch        *b_GenJetAK8_eta;   //!
   TBranch        *b_GenJetAK8_mass;   //!
   TBranch        *b_GenJetAK8_phi;   //!
   TBranch        *b_GenJetAK8_pt;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_nSubGenJetAK8;   //!
   TBranch        *b_SubGenJetAK8_eta;   //!
   TBranch        *b_SubGenJetAK8_mass;   //!
   TBranch        *b_SubGenJetAK8_phi;   //!
   TBranch        *b_SubGenJetAK8_pt;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_nGenVisTau;   //!
   TBranch        *b_GenVisTau_eta;   //!
   TBranch        *b_GenVisTau_mass;   //!
   TBranch        *b_GenVisTau_phi;   //!
   TBranch        *b_GenVisTau_pt;   //!
   TBranch        *b_GenVisTau_charge;   //!
   TBranch        *b_GenVisTau_genPartIdxMother;   //!
   TBranch        *b_GenVisTau_status;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_LHEWeight_originalXWGTUP;   //!
   TBranch        *b_nLHEPdfWeight;   //!
   TBranch        *b_LHEPdfWeight;   //!
   TBranch        *b_nLHEReweightingWeight;   //!
   TBranch        *b_LHEReweightingWeight;   //!
   TBranch        *b_nLHEScaleWeight;   //!
   TBranch        *b_LHEScaleWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nIsoTrack;   //!
   TBranch        *b_IsoTrack_dxy;   //!
   TBranch        *b_IsoTrack_dz;   //!
   TBranch        *b_IsoTrack_eta;   //!
   TBranch        *b_IsoTrack_pfRelIso03_all;   //!
   TBranch        *b_IsoTrack_pfRelIso03_chg;   //!
   TBranch        *b_IsoTrack_phi;   //!
   TBranch        *b_IsoTrack_pt;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_all;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_chg;   //!
   TBranch        *b_IsoTrack_fromPV;   //!
   TBranch        *b_IsoTrack_pdgId;   //!
   TBranch        *b_IsoTrack_isHighPurityTrack;   //!
   TBranch        *b_IsoTrack_isPFcand;   //!
   TBranch        *b_IsoTrack_isFromLostTrack;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagCMVA;   //!
   TBranch        *b_Jet_btagCSVV2;   //!
   TBranch        *b_Jet_btagDeepB;   //!
   TBranch        *b_Jet_btagDeepC;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_btagDeepFlavC;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_jercCHF;   //!
   TBranch        *b_Jet_jercCHPUF;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_Jet_bRegCorr;   //!
   TBranch        *b_Jet_bRegRes;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_LHE_HT;   //!
   TBranch        *b_LHE_HTIncoming;   //!
   TBranch        *b_LHE_Vpt;   //!
   TBranch        *b_LHE_Njets;   //!
   TBranch        *b_LHE_Nb;   //!
   TBranch        *b_LHE_Nc;   //!
   TBranch        *b_LHE_Nuds;   //!
   TBranch        *b_LHE_Nglu;   //!
   TBranch        *b_LHE_NpNLO;   //!
   TBranch        *b_LHE_NpLO;   //!
   TBranch        *b_nLHEPart;   //!
   TBranch        *b_LHEPart_pt;   //!
   TBranch        *b_LHEPart_eta;   //!
   TBranch        *b_LHEPart_phi;   //!
   TBranch        *b_LHEPart_mass;   //!
   TBranch        *b_LHEPart_pdgId;   //!
   TBranch        *b_GenMET_phi;   //!
   TBranch        *b_GenMET_pt;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_MET_covXX;   //!
   TBranch        *b_MET_covXY;   //!
   TBranch        *b_MET_covYY;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetPtRelv2;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_softMva;   //!
   TBranch        *b_Muon_tkRelIso;   //!
   TBranch        *b_Muon_tunepRelPt;   //!
   TBranch        *b_Muon_mvaLowPt;   //!
   TBranch        *b_Muon_mvaTTH;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_fsrPhotonIdx;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_eCorr;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_mass;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_mvaIDV1;   //!
   TBranch        *b_Photon_pfRelIso03_all;   //!
   TBranch        *b_Photon_pfRelIso03_chg;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_charge;   //!
   TBranch        *b_Photon_cutBasedBitmap;   //!
   TBranch        *b_Photon_cutBasedV1Bitmap;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
   TBranch        *b_Photon_pdgId;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_seedGain;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nGenDressedLepton;   //!
   TBranch        *b_GenDressedLepton_eta;   //!
   TBranch        *b_GenDressedLepton_mass;   //!
   TBranch        *b_GenDressedLepton_phi;   //!
   TBranch        *b_GenDressedLepton_pt;   //!
   TBranch        *b_GenDressedLepton_pdgId;   //!
   TBranch        *b_GenDressedLepton_hasTauAnc;   //!
   TBranch        *b_nSoftActivityJet;   //!
   TBranch        *b_SoftActivityJet_eta;   //!
   TBranch        *b_SoftActivityJet_phi;   //!
   TBranch        *b_SoftActivityJet_pt;   //!
   TBranch        *b_SoftActivityJetHT;   //!
   TBranch        *b_SoftActivityJetHT10;   //!
   TBranch        *b_SoftActivityJetHT2;   //!
   TBranch        *b_SoftActivityJetHT5;   //!
   TBranch        *b_SoftActivityJetNjets10;   //!
   TBranch        *b_SoftActivityJetNjets2;   //!
   TBranch        *b_SoftActivityJetNjets5;   //!
   TBranch        *b_nSubJet;   //!
   TBranch        *b_SubJet_btagCMVA;   //!
   TBranch        *b_SubJet_btagCSVV2;   //!
   TBranch        *b_SubJet_btagDeepB;   //!
   TBranch        *b_SubJet_eta;   //!
   TBranch        *b_SubJet_mass;   //!
   TBranch        *b_SubJet_n2b1;   //!
   TBranch        *b_SubJet_n3b1;   //!
   TBranch        *b_SubJet_phi;   //!
   TBranch        *b_SubJet_pt;   //!
   TBranch        *b_SubJet_rawFactor;   //!
   TBranch        *b_SubJet_tau1;   //!
   TBranch        *b_SubJet_tau2;   //!
   TBranch        *b_SubJet_tau3;   //!
   TBranch        *b_SubJet_tau4;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_Tau_chargedIso;   //!
   TBranch        *b_Tau_dxy;   //!
   TBranch        *b_Tau_dz;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_leadTkDeltaEta;   //!
   TBranch        *b_Tau_leadTkDeltaPhi;   //!
   TBranch        *b_Tau_leadTkPtOverTauPt;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_neutralIso;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_photonsOutsideSignalCone;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_puCorr;   //!
   TBranch        *b_Tau_rawAntiEle;   //!
   TBranch        *b_Tau_rawAntiEle2018;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSmu;   //!
   TBranch        *b_Tau_rawIso;   //!
   TBranch        *b_Tau_rawIsodR03;   //!
   TBranch        *b_Tau_rawMVAnewDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDM;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v1;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDMdR032017v2;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_jetIdx;   //!
   TBranch        *b_Tau_rawAntiEleCat;   //!
   TBranch        *b_Tau_rawAntiEleCat2018;   //!
   TBranch        *b_Tau_idAntiEle;   //!
   TBranch        *b_Tau_idAntiEle2018;   //!
   TBranch        *b_Tau_idAntiMu;   //!
   TBranch        *b_Tau_idDecayMode;   //!
   TBranch        *b_Tau_idDecayModeNewDMs;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSmu;   //!
   TBranch        *b_Tau_idMVAnewDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDM;   //!
   TBranch        *b_Tau_idMVAoldDM2017v1;   //!
   TBranch        *b_Tau_idMVAoldDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDMdR032017v2;   //!
   TBranch        *b_TkMET_phi;   //!
   TBranch        *b_TkMET_pt;   //!
   TBranch        *b_TkMET_sumEt;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_genTtbarId;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   TBranch        *b_GenJetAK8_partonFlavour;   //!
   TBranch        *b_GenJetAK8_hadronFlavour;   //!
   TBranch        *b_GenJet_partonFlavour;   //!
   TBranch        *b_GenJet_hadronFlavour;   //!
   TBranch        *b_Jet_genJetIdx;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_Photon_genPartIdx;   //!
   TBranch        *b_Photon_genPartFlav;   //!
   TBranch        *b_MET_fiducialGenPhi;   //!
   TBranch        *b_MET_fiducialGenPt;   //!
   TBranch        *b_Electron_cleanmask;   //!
   TBranch        *b_Jet_cleanmask;   //!
   TBranch        *b_Muon_cleanmask;   //!
   TBranch        *b_Photon_cleanmask;   //!
   TBranch        *b_Tau_cleanmask;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_Tau_genPartIdx;   //!
   TBranch        *b_Tau_genPartFlav;   //!
   TBranch        *b_L1simulation_step;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLT_DoubleEle25_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle27_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle33_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle24_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_Ele27_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_TkMu27;   //!
   TBranch        *b_HLT_DoubleMu4_3_Bs;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi;   //!
   TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu;   //!
   TBranch        *b_HLT_DoubleMu3_TkMu_DsTau3Mu;   //!
   TBranch        *b_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track2_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track3p5_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track7_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track2_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track3p5_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track7_Upsilon;   //!
   TBranch        *b_HLT_Mu3_L1SingleMu5orSingleMu7;   //!
   TBranch        *b_HLT_Ele20_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele15_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele17_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele20_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele20_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;   //!
   TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele28_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele30_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf_L1EGMT;   //!
   TBranch        *b_HLT_Ele38_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele40_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_L1DoubleEG;   //!
   TBranch        *b_HLT_IsoMu20;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_IsoMu30;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu25_TkMu0_Onia;   //!
   TBranch        *b_HLT_Mu30_TkMu0_Psi;   //!
   TBranch        *b_HLT_Mu30_TkMu0_Upsilon;   //!
   TBranch        *b_HLT_Mu20_TkMu0_Phi;   //!
   TBranch        *b_HLT_Mu25_TkMu0_Phi;   //!
   TBranch        *b_HLT_Mu12;   //!
   TBranch        *b_HLT_Mu15;   //!
   TBranch        *b_HLT_Mu20;   //!
   TBranch        *b_HLT_Mu27;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_OldMu100;   //!
   TBranch        *b_HLT_TkMu100;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL;   //!
   TBranch        *b_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_DoubleIsoMu20_eta2p1;   //!
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu17;   //!
   TBranch        *b_HLT_Mu19;   //!
   TBranch        *b_HLT_Mu18_Mu9;   //!
   TBranch        *b_HLT_Mu18_Mu9_DZ;   //!
   TBranch        *b_HLT_Mu20_Mu10;   //!
   TBranch        *b_HLT_Mu20_Mu10_DZ;   //!
   TBranch        *b_HLT_Mu23_Mu12;   //!
   TBranch        *b_HLT_Mu23_Mu12_DZ;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_Jet_pt_raw;   //!
   TBranch        *b_Jet_pt_nom;   //!
   TBranch        *b_Jet_mass_raw;   //!
   TBranch        *b_Jet_mass_nom;   //!
   TBranch        *b_Jet_corr_JEC;   //!
   TBranch        *b_Jet_corr_JER;   //!
   TBranch        *b_MET_pt_nom;   //!
   TBranch        *b_MET_phi_nom;   //!
   TBranch        *b_MET_pt_jer;   //!
   TBranch        *b_MET_phi_jer;   //!
   TBranch        *b_Jet_pt_jerUp;   //!
   TBranch        *b_Jet_mass_jerUp;   //!
   TBranch        *b_MET_pt_jerUp;   //!
   TBranch        *b_MET_phi_jerUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteStatUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteScaleUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteSampleUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteSampleUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteSampleUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteSampleUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_phi_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_Jet_pt_jesFragmentationUp;   //!
   TBranch        *b_Jet_mass_jesFragmentationUp;   //!
   TBranch        *b_MET_pt_jesFragmentationUp;   //!
   TBranch        *b_MET_phi_jesFragmentationUp;   //!
   TBranch        *b_Jet_pt_jesSinglePionECALUp;   //!
   TBranch        *b_Jet_mass_jesSinglePionECALUp;   //!
   TBranch        *b_MET_pt_jesSinglePionECALUp;   //!
   TBranch        *b_MET_phi_jesSinglePionECALUp;   //!
   TBranch        *b_Jet_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_Jet_mass_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_phi_jesSinglePionHCALUp;   //!
   TBranch        *b_Jet_pt_jesFlavorQCDUp;   //!
   TBranch        *b_Jet_mass_jesFlavorQCDUp;   //!
   TBranch        *b_MET_pt_jesFlavorQCDUp;   //!
   TBranch        *b_MET_phi_jesFlavorQCDUp;   //!
   TBranch        *b_Jet_pt_jesTimePtEtaUp;   //!
   TBranch        *b_Jet_mass_jesTimePtEtaUp;   //!
   TBranch        *b_MET_pt_jesTimePtEtaUp;   //!
   TBranch        *b_MET_phi_jesTimePtEtaUp;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC1Up;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC2Up;   //!
   TBranch        *b_Jet_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_phi_jesRelativeJERHFUp;   //!
   TBranch        *b_Jet_pt_jesRelativePtBBUp;   //!
   TBranch        *b_Jet_mass_jesRelativePtBBUp;   //!
   TBranch        *b_MET_pt_jesRelativePtBBUp;   //!
   TBranch        *b_MET_phi_jesRelativePtBBUp;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_phi_jesRelativePtEC1Up;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_phi_jesRelativePtEC2Up;   //!
   TBranch        *b_Jet_pt_jesRelativePtHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativePtHFUp;   //!
   TBranch        *b_MET_pt_jesRelativePtHFUp;   //!
   TBranch        *b_MET_phi_jesRelativePtHFUp;   //!
   TBranch        *b_Jet_pt_jesRelativeBalUp;   //!
   TBranch        *b_Jet_mass_jesRelativeBalUp;   //!
   TBranch        *b_MET_pt_jesRelativeBalUp;   //!
   TBranch        *b_MET_phi_jesRelativeBalUp;   //!
   TBranch        *b_Jet_pt_jesRelativeSampleUp;   //!
   TBranch        *b_Jet_mass_jesRelativeSampleUp;   //!
   TBranch        *b_MET_pt_jesRelativeSampleUp;   //!
   TBranch        *b_MET_phi_jesRelativeSampleUp;   //!
   TBranch        *b_Jet_pt_jesRelativeFSRUp;   //!
   TBranch        *b_Jet_mass_jesRelativeFSRUp;   //!
   TBranch        *b_MET_pt_jesRelativeFSRUp;   //!
   TBranch        *b_MET_phi_jesRelativeFSRUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_phi_jesRelativeStatFSRUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatECUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatECUp;   //!
   TBranch        *b_MET_pt_jesRelativeStatECUp;   //!
   TBranch        *b_MET_phi_jesRelativeStatECUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_phi_jesRelativeStatHFUp;   //!
   TBranch        *b_Jet_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_Jet_mass_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_phi_jesPileUpDataMCUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_phi_jesPileUpPtRefUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_phi_jesPileUpPtBBUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC1Up;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC2Up;   //!
   TBranch        *b_Jet_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_phi_jesPileUpPtHFUp;   //!
   TBranch        *b_Jet_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_Jet_mass_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_phi_jesPileUpMuZeroUp;   //!
   TBranch        *b_Jet_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_Jet_mass_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_phi_jesPileUpEnvelopeUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_phi_jesSubTotalPileUpUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_phi_jesSubTotalRelativeUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalPtUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalPtUp;   //!
   TBranch        *b_MET_pt_jesSubTotalPtUp;   //!
   TBranch        *b_MET_phi_jesSubTotalPtUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_phi_jesSubTotalScaleUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_phi_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalMCUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalMCUp;   //!
   TBranch        *b_MET_pt_jesSubTotalMCUp;   //!
   TBranch        *b_MET_phi_jesSubTotalMCUp;   //!
   TBranch        *b_Jet_pt_jesTotalUp;   //!
   TBranch        *b_Jet_mass_jesTotalUp;   //!
   TBranch        *b_MET_pt_jesTotalUp;   //!
   TBranch        *b_MET_phi_jesTotalUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_phi_jesTotalNoTimeUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_Jet_pt_jesFlavorZJetUp;   //!
   TBranch        *b_Jet_mass_jesFlavorZJetUp;   //!
   TBranch        *b_MET_pt_jesFlavorZJetUp;   //!
   TBranch        *b_MET_phi_jesFlavorZJetUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_phi_jesFlavorPhotonJetUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureGluonUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureQuarkUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureCharmUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_phi_jesFlavorPureBottomUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunAUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunAUp;   //!
   TBranch        *b_MET_pt_jesTimeRunAUp;   //!
   TBranch        *b_MET_phi_jesTimeRunAUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunBUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunBUp;   //!
   TBranch        *b_MET_pt_jesTimeRunBUp;   //!
   TBranch        *b_MET_phi_jesTimeRunBUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunCUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunCUp;   //!
   TBranch        *b_MET_pt_jesTimeRunCUp;   //!
   TBranch        *b_MET_phi_jesTimeRunCUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunDUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunDUp;   //!
   TBranch        *b_MET_pt_jesTimeRunDUp;   //!
   TBranch        *b_MET_phi_jesTimeRunDUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_pt_unclustEnUp;   //!
   TBranch        *b_MET_phi_unclustEnUp;   //!
   TBranch        *b_Jet_pt_jerDown;   //!
   TBranch        *b_Jet_mass_jerDown;   //!
   TBranch        *b_MET_pt_jerDown;   //!
   TBranch        *b_MET_phi_jerDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteStatDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteScaleDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteSampleDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteSampleDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteSampleDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteSampleDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_phi_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_Jet_pt_jesFragmentationDown;   //!
   TBranch        *b_Jet_mass_jesFragmentationDown;   //!
   TBranch        *b_MET_pt_jesFragmentationDown;   //!
   TBranch        *b_MET_phi_jesFragmentationDown;   //!
   TBranch        *b_Jet_pt_jesSinglePionECALDown;   //!
   TBranch        *b_Jet_mass_jesSinglePionECALDown;   //!
   TBranch        *b_MET_pt_jesSinglePionECALDown;   //!
   TBranch        *b_MET_phi_jesSinglePionECALDown;   //!
   TBranch        *b_Jet_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_Jet_mass_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_phi_jesSinglePionHCALDown;   //!
   TBranch        *b_Jet_pt_jesFlavorQCDDown;   //!
   TBranch        *b_Jet_mass_jesFlavorQCDDown;   //!
   TBranch        *b_MET_pt_jesFlavorQCDDown;   //!
   TBranch        *b_MET_phi_jesFlavorQCDDown;   //!
   TBranch        *b_Jet_pt_jesTimePtEtaDown;   //!
   TBranch        *b_Jet_mass_jesTimePtEtaDown;   //!
   TBranch        *b_MET_pt_jesTimePtEtaDown;   //!
   TBranch        *b_MET_phi_jesTimePtEtaDown;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC1Down;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_phi_jesRelativeJEREC2Down;   //!
   TBranch        *b_Jet_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_phi_jesRelativeJERHFDown;   //!
   TBranch        *b_Jet_pt_jesRelativePtBBDown;   //!
   TBranch        *b_Jet_mass_jesRelativePtBBDown;   //!
   TBranch        *b_MET_pt_jesRelativePtBBDown;   //!
   TBranch        *b_MET_phi_jesRelativePtBBDown;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_phi_jesRelativePtEC1Down;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_phi_jesRelativePtEC2Down;   //!
   TBranch        *b_Jet_pt_jesRelativePtHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativePtHFDown;   //!
   TBranch        *b_MET_pt_jesRelativePtHFDown;   //!
   TBranch        *b_MET_phi_jesRelativePtHFDown;   //!
   TBranch        *b_Jet_pt_jesRelativeBalDown;   //!
   TBranch        *b_Jet_mass_jesRelativeBalDown;   //!
   TBranch        *b_MET_pt_jesRelativeBalDown;   //!
   TBranch        *b_MET_phi_jesRelativeBalDown;   //!
   TBranch        *b_Jet_pt_jesRelativeSampleDown;   //!
   TBranch        *b_Jet_mass_jesRelativeSampleDown;   //!
   TBranch        *b_MET_pt_jesRelativeSampleDown;   //!
   TBranch        *b_MET_phi_jesRelativeSampleDown;   //!
   TBranch        *b_Jet_pt_jesRelativeFSRDown;   //!
   TBranch        *b_Jet_mass_jesRelativeFSRDown;   //!
   TBranch        *b_MET_pt_jesRelativeFSRDown;   //!
   TBranch        *b_MET_phi_jesRelativeFSRDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_phi_jesRelativeStatFSRDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatECDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatECDown;   //!
   TBranch        *b_MET_pt_jesRelativeStatECDown;   //!
   TBranch        *b_MET_phi_jesRelativeStatECDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_phi_jesRelativeStatHFDown;   //!
   TBranch        *b_Jet_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_Jet_mass_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_phi_jesPileUpDataMCDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_phi_jesPileUpPtRefDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_phi_jesPileUpPtBBDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC1Down;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_phi_jesPileUpPtEC2Down;   //!
   TBranch        *b_Jet_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_phi_jesPileUpPtHFDown;   //!
   TBranch        *b_Jet_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_Jet_mass_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_phi_jesPileUpMuZeroDown;   //!
   TBranch        *b_Jet_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_Jet_mass_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_phi_jesPileUpEnvelopeDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_phi_jesSubTotalPileUpDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_phi_jesSubTotalRelativeDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalPtDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalPtDown;   //!
   TBranch        *b_MET_pt_jesSubTotalPtDown;   //!
   TBranch        *b_MET_phi_jesSubTotalPtDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_phi_jesSubTotalScaleDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_phi_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalMCDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalMCDown;   //!
   TBranch        *b_MET_pt_jesSubTotalMCDown;   //!
   TBranch        *b_MET_phi_jesSubTotalMCDown;   //!
   TBranch        *b_Jet_pt_jesTotalDown;   //!
   TBranch        *b_Jet_mass_jesTotalDown;   //!
   TBranch        *b_MET_pt_jesTotalDown;   //!
   TBranch        *b_MET_phi_jesTotalDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_phi_jesTotalNoTimeDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_phi_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_Jet_pt_jesFlavorZJetDown;   //!
   TBranch        *b_Jet_mass_jesFlavorZJetDown;   //!
   TBranch        *b_MET_pt_jesFlavorZJetDown;   //!
   TBranch        *b_MET_phi_jesFlavorZJetDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_phi_jesFlavorPhotonJetDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureGluonDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureQuarkDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureCharmDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_phi_jesFlavorPureBottomDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunADown;   //!
   TBranch        *b_Jet_mass_jesTimeRunADown;   //!
   TBranch        *b_MET_pt_jesTimeRunADown;   //!
   TBranch        *b_MET_phi_jesTimeRunADown;   //!
   TBranch        *b_Jet_pt_jesTimeRunBDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunBDown;   //!
   TBranch        *b_MET_pt_jesTimeRunBDown;   //!
   TBranch        *b_MET_phi_jesTimeRunBDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunCDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunCDown;   //!
   TBranch        *b_MET_pt_jesTimeRunCDown;   //!
   TBranch        *b_MET_phi_jesTimeRunCDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunDDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunDDown;   //!
   TBranch        *b_MET_pt_jesTimeRunDDown;   //!
   TBranch        *b_MET_phi_jesTimeRunDDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_phi_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_pt_unclustEnDown;   //!
   TBranch        *b_MET_phi_unclustEnDown;   //!
   TBranch        *b_FatJet_pt_raw;   //!
   TBranch        *b_FatJet_pt_nom;   //!
   TBranch        *b_FatJet_mass_raw;   //!
   TBranch        *b_FatJet_mass_nom;   //!
   TBranch        *b_FatJet_corr_JEC;   //!
   TBranch        *b_FatJet_corr_JER;   //!
   TBranch        *b_FatJet_corr_JMS;   //!
   TBranch        *b_FatJet_corr_JMR;   //!
   TBranch        *b_FatJet_msoftdrop_raw;   //!
   TBranch        *b_FatJet_msoftdrop_nom;   //!
   TBranch        *b_FatJet_msoftdrop_corr_JMR;   //!
   TBranch        *b_FatJet_msoftdrop_corr_JMS;   //!
   TBranch        *b_FatJet_msoftdrop_corr_PUPPI;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_nom;   //!
   TBranch        *b_FatJet_pt_jerUp;   //!
   TBranch        *b_FatJet_mass_jerUp;   //!
   TBranch        *b_FatJet_mass_jmrUp;   //!
   TBranch        *b_FatJet_mass_jmsUp;   //!
   TBranch        *b_FatJet_msoftdrop_jerUp;   //!
   TBranch        *b_FatJet_msoftdrop_jmrUp;   //!
   TBranch        *b_FatJet_msoftdrop_jmsUp;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jerUp;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmrUp;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmsUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteStatUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteStatUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteScaleUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteScaleUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteSampleUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteSampleUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteSampleUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_FatJet_pt_jesFragmentationUp;   //!
   TBranch        *b_FatJet_mass_jesFragmentationUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFragmentationUp;   //!
   TBranch        *b_FatJet_pt_jesSinglePionECALUp;   //!
   TBranch        *b_FatJet_mass_jesSinglePionECALUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionECALUp;   //!
   TBranch        *b_FatJet_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_FatJet_mass_jesSinglePionHCALUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionHCALUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorQCDUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorQCDUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorQCDUp;   //!
   TBranch        *b_FatJet_pt_jesTimePtEtaUp;   //!
   TBranch        *b_FatJet_mass_jesTimePtEtaUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimePtEtaUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC1Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC1Up;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC2Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC2Up;   //!
   TBranch        *b_FatJet_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeJERHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJERHFUp;   //!
   TBranch        *b_FatJet_pt_jesRelativePtBBUp;   //!
   TBranch        *b_FatJet_mass_jesRelativePtBBUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtBBUp;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC1Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC1Up;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC2Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC2Up;   //!
   TBranch        *b_FatJet_pt_jesRelativePtHFUp;   //!
   TBranch        *b_FatJet_mass_jesRelativePtHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtHFUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeBalUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeBalUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeBalUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeSampleUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeSampleUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeSampleUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeFSRUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeFSRUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeFSRUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatFSRUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatFSRUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatECUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatECUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatECUp;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatHFUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpDataMCUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpDataMCUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtRefUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtRefUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtBBUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtBBUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC1Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC1Up;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC2Up;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC2Up;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtHFUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtHFUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpMuZeroUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpMuZeroUp;   //!
   TBranch        *b_FatJet_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_FatJet_mass_jesPileUpEnvelopeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpEnvelopeUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPileUpUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPileUpUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalRelativeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalRelativeUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPtUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPtUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPtUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalScaleUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalScaleUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_FatJet_pt_jesSubTotalMCUp;   //!
   TBranch        *b_FatJet_mass_jesSubTotalMCUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalMCUp;   //!
   TBranch        *b_FatJet_pt_jesTotalUp;   //!
   TBranch        *b_FatJet_mass_jesTotalUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalUp;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorUp;   //!
   TBranch        *b_FatJet_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_FatJet_mass_jesTotalNoTimeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoTimeUp;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorZJetUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorZJetUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorZJetUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPhotonJetUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPhotonJetUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureGluonUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureGluonUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureQuarkUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureQuarkUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureCharmUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureCharmUp;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureBottomUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureBottomUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunAUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunAUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunAUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunBUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunBUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunBUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunCUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunCUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunCUp;   //!
   TBranch        *b_FatJet_pt_jesTimeRunDUp;   //!
   TBranch        *b_FatJet_mass_jesTimeRunDUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunDUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_FatJet_pt_jerDown;   //!
   TBranch        *b_FatJet_mass_jerDown;   //!
   TBranch        *b_FatJet_mass_jmrDown;   //!
   TBranch        *b_FatJet_mass_jmsDown;   //!
   TBranch        *b_FatJet_msoftdrop_jerDown;   //!
   TBranch        *b_FatJet_msoftdrop_jmrDown;   //!
   TBranch        *b_FatJet_msoftdrop_jmsDown;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jerDown;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmrDown;   //!
   TBranch        *b_FatJet_msoftdrop_tau21DDT_jmsDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteStatDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteStatDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteScaleDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteScaleDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteSampleDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteSampleDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteSampleDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_FatJet_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_FatJet_mass_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_FatJet_pt_jesFragmentationDown;   //!
   TBranch        *b_FatJet_mass_jesFragmentationDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFragmentationDown;   //!
   TBranch        *b_FatJet_pt_jesSinglePionECALDown;   //!
   TBranch        *b_FatJet_mass_jesSinglePionECALDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionECALDown;   //!
   TBranch        *b_FatJet_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_FatJet_mass_jesSinglePionHCALDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSinglePionHCALDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorQCDDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorQCDDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorQCDDown;   //!
   TBranch        *b_FatJet_pt_jesTimePtEtaDown;   //!
   TBranch        *b_FatJet_mass_jesTimePtEtaDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimePtEtaDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC1Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC1Down;   //!
   TBranch        *b_FatJet_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_FatJet_mass_jesRelativeJEREC2Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJEREC2Down;   //!
   TBranch        *b_FatJet_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeJERHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeJERHFDown;   //!
   TBranch        *b_FatJet_pt_jesRelativePtBBDown;   //!
   TBranch        *b_FatJet_mass_jesRelativePtBBDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtBBDown;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC1Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC1Down;   //!
   TBranch        *b_FatJet_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_FatJet_mass_jesRelativePtEC2Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtEC2Down;   //!
   TBranch        *b_FatJet_pt_jesRelativePtHFDown;   //!
   TBranch        *b_FatJet_mass_jesRelativePtHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativePtHFDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeBalDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeBalDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeBalDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeSampleDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeSampleDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeSampleDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeFSRDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeFSRDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeFSRDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatFSRDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatFSRDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatECDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatECDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatECDown;   //!
   TBranch        *b_FatJet_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_FatJet_mass_jesRelativeStatHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesRelativeStatHFDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpDataMCDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpDataMCDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtRefDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtRefDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtBBDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtBBDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC1Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC1Down;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtEC2Down;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtEC2Down;   //!
   TBranch        *b_FatJet_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpPtHFDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpPtHFDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpMuZeroDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpMuZeroDown;   //!
   TBranch        *b_FatJet_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_FatJet_mass_jesPileUpEnvelopeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesPileUpEnvelopeDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPileUpDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPileUpDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalRelativeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalRelativeDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalPtDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalPtDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalPtDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalScaleDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalScaleDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_FatJet_pt_jesSubTotalMCDown;   //!
   TBranch        *b_FatJet_mass_jesSubTotalMCDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesSubTotalMCDown;   //!
   TBranch        *b_FatJet_pt_jesTotalDown;   //!
   TBranch        *b_FatJet_mass_jesTotalDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalDown;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorDown;   //!
   TBranch        *b_FatJet_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_FatJet_mass_jesTotalNoTimeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoTimeDown;   //!
   TBranch        *b_FatJet_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_FatJet_mass_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorZJetDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorZJetDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorZJetDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPhotonJetDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPhotonJetDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureGluonDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureGluonDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureQuarkDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureQuarkDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureCharmDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureCharmDown;   //!
   TBranch        *b_FatJet_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_FatJet_mass_jesFlavorPureBottomDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesFlavorPureBottomDown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunADown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunADown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunADown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunBDown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunBDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunBDown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunCDown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunCDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunCDown;   //!
   TBranch        *b_FatJet_pt_jesTimeRunDDown;   //!
   TBranch        *b_FatJet_mass_jesTimeRunDDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesTimeRunDDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_FatJet_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_FatJet_mass_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown;   //!

  NanoAOD_MC(TTree *tree=0) { if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root");
      if (!f || !f->IsOpen()) {
	f = new TFile("root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root");
      }
      f->GetObject("Events",tree);

    }
    Init(tree);
  };
  virtual ~NanoAOD_MC() { 
    //if (!fChain) return;
    //delete fChain->GetCurrentFile();
  };
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
};

#endif

//#ifdef NanoAOD_MC_cxx
//NanoAOD_MC::NanoAOD_MC(TTree *tree) : fChain(0) 
//{
//// if parameter tree is not specified (or zero), connect the file
//// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root");
//      }
//      f->GetObject("Events",tree);
//
//   }
//   Init(tree);
//}
//
//NanoAOD_MC::~NanoAOD_MC()
//{
//   if (!fChain) return;
//   delete fChain->GetCurrentFile();
//}

Int_t NanoAOD_MC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NanoAOD_MC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void NanoAOD_MC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2, &b_btagWeight_CSVV2);
   fChain->SetBranchAddress("btagWeight_DeepCSVB", &btagWeight_DeepCSVB, &b_btagWeight_DeepCSVB);
   fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   fChain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
   fChain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
   fChain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   fChain->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet, &b_nCorrT1METJet);
   fChain->SetBranchAddress("CorrT1METJet_area", CorrT1METJet_area, &b_CorrT1METJet_area);
   fChain->SetBranchAddress("CorrT1METJet_eta", CorrT1METJet_eta, &b_CorrT1METJet_eta);
   fChain->SetBranchAddress("CorrT1METJet_muonSubtrFactor", CorrT1METJet_muonSubtrFactor, &b_CorrT1METJet_muonSubtrFactor);
   fChain->SetBranchAddress("CorrT1METJet_phi", CorrT1METJet_phi, &b_CorrT1METJet_phi);
   fChain->SetBranchAddress("CorrT1METJet_rawPt", CorrT1METJet_rawPt, &b_CorrT1METJet_rawPt);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2, &b_Electron_jetPtRelv2);
   fChain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso", Electron_mvaFall17V1Iso, &b_Electron_mvaFall17V1Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso", Electron_mvaFall17V1noIso, &b_Electron_mvaFall17V1noIso);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_Fall17_V1", Electron_cutBased_Fall17_V1, &b_Electron_cutBased_Fall17_V1);
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmapHEEP", Electron_vidNestedWPBitmapHEEP, &b_Electron_vidNestedWPBitmapHEEP);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP80", Electron_mvaFall17V1Iso_WP80, &b_Electron_mvaFall17V1Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP90", Electron_mvaFall17V1Iso_WP90, &b_Electron_mvaFall17V1Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WPL", Electron_mvaFall17V1Iso_WPL, &b_Electron_mvaFall17V1Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", Electron_mvaFall17V1noIso_WP80, &b_Electron_mvaFall17V1noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", Electron_mvaFall17V1noIso_WP90, &b_Electron_mvaFall17V1noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", Electron_mvaFall17V1noIso_WPL, &b_Electron_mvaFall17V1noIso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80, &b_Electron_mvaFall17V2Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90, &b_Electron_mvaFall17V2Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL, &b_Electron_mvaFall17V2Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80, &b_Electron_mvaFall17V2noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90, &b_Electron_mvaFall17V2noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", Electron_mvaFall17V2noIso_WPL, &b_Electron_mvaFall17V2noIso_WPL);
   fChain->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilterV2", &Flag_ecalBadCalibFilterV2, &b_Flag_ecalBadCalibFilterV2);
   fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   fChain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   fChain->SetBranchAddress("FatJet_btagCMVA", FatJet_btagCMVA, &b_FatJet_btagCMVA);
   fChain->SetBranchAddress("FatJet_btagCSVV2", FatJet_btagCSVV2, &b_FatJet_btagCSVV2);
   fChain->SetBranchAddress("FatJet_btagDDBvL", FatJet_btagDDBvL, &b_FatJet_btagDDBvL);
   fChain->SetBranchAddress("FatJet_btagDDBvL_noMD", FatJet_btagDDBvL_noMD, &b_FatJet_btagDDBvL_noMD);
   fChain->SetBranchAddress("FatJet_btagDDCvB", FatJet_btagDDCvB, &b_FatJet_btagDDCvB);
   fChain->SetBranchAddress("FatJet_btagDDCvB_noMD", FatJet_btagDDCvB_noMD, &b_FatJet_btagDDCvB_noMD);
   fChain->SetBranchAddress("FatJet_btagDDCvL", FatJet_btagDDCvL, &b_FatJet_btagDDCvL);
   fChain->SetBranchAddress("FatJet_btagDDCvL_noMD", FatJet_btagDDCvL_noMD, &b_FatJet_btagDDCvL_noMD);
   fChain->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
   fChain->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
   fChain->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD", FatJet_deepTagMD_H4qvsQCD, &b_FatJet_deepTagMD_H4qvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD", FatJet_deepTagMD_HbbvsQCD, &b_FatJet_deepTagMD_HbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_TvsQCD", FatJet_deepTagMD_TvsQCD, &b_FatJet_deepTagMD_TvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_WvsQCD", FatJet_deepTagMD_WvsQCD, &b_FatJet_deepTagMD_WvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD", FatJet_deepTagMD_ZHbbvsQCD, &b_FatJet_deepTagMD_ZHbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD", FatJet_deepTagMD_ZHccvsQCD, &b_FatJet_deepTagMD_ZHccvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZbbvsQCD", FatJet_deepTagMD_ZbbvsQCD, &b_FatJet_deepTagMD_ZbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZvsQCD", FatJet_deepTagMD_ZvsQCD, &b_FatJet_deepTagMD_ZvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_bbvsLight", FatJet_deepTagMD_bbvsLight, &b_FatJet_deepTagMD_bbvsLight);
   fChain->SetBranchAddress("FatJet_deepTagMD_ccvsLight", FatJet_deepTagMD_ccvsLight, &b_FatJet_deepTagMD_ccvsLight);
   fChain->SetBranchAddress("FatJet_deepTag_H", FatJet_deepTag_H, &b_FatJet_deepTag_H);
   fChain->SetBranchAddress("FatJet_deepTag_QCD", FatJet_deepTag_QCD, &b_FatJet_deepTag_QCD);
   fChain->SetBranchAddress("FatJet_deepTag_QCDothers", FatJet_deepTag_QCDothers, &b_FatJet_deepTag_QCDothers);
   fChain->SetBranchAddress("FatJet_deepTag_TvsQCD", FatJet_deepTag_TvsQCD, &b_FatJet_deepTag_TvsQCD);
   fChain->SetBranchAddress("FatJet_deepTag_WvsQCD", FatJet_deepTag_WvsQCD, &b_FatJet_deepTag_WvsQCD);
   fChain->SetBranchAddress("FatJet_deepTag_ZvsQCD", FatJet_deepTag_ZvsQCD, &b_FatJet_deepTag_ZvsQCD);
   fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   fChain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   fChain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   fChain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   fChain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   fChain->SetBranchAddress("nFsrPhoton", &nFsrPhoton, &b_nFsrPhoton);
   fChain->SetBranchAddress("FsrPhoton_dROverEt2", FsrPhoton_dROverEt2, &b_FsrPhoton_dROverEt2);
   fChain->SetBranchAddress("FsrPhoton_eta", FsrPhoton_eta, &b_FsrPhoton_eta);
   fChain->SetBranchAddress("FsrPhoton_phi", FsrPhoton_phi, &b_FsrPhoton_phi);
   fChain->SetBranchAddress("FsrPhoton_pt", FsrPhoton_pt, &b_FsrPhoton_pt);
   fChain->SetBranchAddress("FsrPhoton_relIso03", FsrPhoton_relIso03, &b_FsrPhoton_relIso03);
   fChain->SetBranchAddress("FsrPhoton_muonIdx", FsrPhoton_muonIdx, &b_FsrPhoton_muonIdx);
   fChain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   fChain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   fChain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   fChain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   fChain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8, &b_nSubGenJetAK8);
   fChain->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta, &b_SubGenJetAK8_eta);
   fChain->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass, &b_SubGenJetAK8_mass);
   fChain->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi, &b_SubGenJetAK8_phi);
   fChain->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt, &b_SubGenJetAK8_pt);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   fChain->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
   fChain->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
   fChain->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
   fChain->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
   fChain->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
   fChain->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   fChain->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP);
   fChain->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight, &b_nLHEPdfWeight);
   fChain->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight, &b_LHEPdfWeight);
   fChain->SetBranchAddress("nLHEReweightingWeight", &nLHEReweightingWeight, &b_nLHEReweightingWeight);
   fChain->SetBranchAddress("LHEReweightingWeight", &LHEReweightingWeight, &b_LHEReweightingWeight);
   fChain->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight, &b_nLHEScaleWeight);
   fChain->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight, &b_LHEScaleWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   fChain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   fChain->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   fChain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   fChain->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   fChain->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   fChain->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
   fChain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   fChain->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_btagDeepFlavC", Jet_btagDeepFlavC, &b_Jet_btagDeepFlavC);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_jercCHF", Jet_jercCHF, &b_Jet_jercCHF);
   fChain->SetBranchAddress("Jet_jercCHPUF", Jet_jercCHPUF, &b_Jet_jercCHPUF);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("Jet_bRegCorr", Jet_bRegCorr, &b_Jet_bRegCorr);
   fChain->SetBranchAddress("Jet_bRegRes", Jet_bRegRes, &b_Jet_bRegRes);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);
   fChain->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming, &b_LHE_HTIncoming);
   fChain->SetBranchAddress("LHE_Vpt", &LHE_Vpt, &b_LHE_Vpt);
   fChain->SetBranchAddress("LHE_Njets", &LHE_Njets, &b_LHE_Njets);
   fChain->SetBranchAddress("LHE_Nb", &LHE_Nb, &b_LHE_Nb);
   fChain->SetBranchAddress("LHE_Nc", &LHE_Nc, &b_LHE_Nc);
   fChain->SetBranchAddress("LHE_Nuds", &LHE_Nuds, &b_LHE_Nuds);
   fChain->SetBranchAddress("LHE_Nglu", &LHE_Nglu, &b_LHE_Nglu);
   fChain->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO, &b_LHE_NpNLO);
   fChain->SetBranchAddress("LHE_NpLO", &LHE_NpLO, &b_LHE_NpLO);
   fChain->SetBranchAddress("nLHEPart", &nLHEPart, &b_nLHEPart);
   fChain->SetBranchAddress("LHEPart_pt", LHEPart_pt, &b_LHEPart_pt);
   fChain->SetBranchAddress("LHEPart_eta", LHEPart_eta, &b_LHEPart_eta);
   fChain->SetBranchAddress("LHEPart_phi", LHEPart_phi, &b_LHEPart_phi);
   fChain->SetBranchAddress("LHEPart_mass", LHEPart_mass, &b_LHEPart_mass);
   fChain->SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId, &b_LHEPart_pdgId);
   fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   fChain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   fChain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   fChain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_softMva", Muon_softMva, &b_Muon_softMva);
   fChain->SetBranchAddress("Muon_tkRelIso", Muon_tkRelIso, &b_Muon_tkRelIso);
   fChain->SetBranchAddress("Muon_tunepRelPt", Muon_tunepRelPt, &b_Muon_tunepRelPt);
   fChain->SetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt, &b_Muon_mvaLowPt);
   fChain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_fsrPhotonIdx", Muon_fsrPhotonIdx, &b_Muon_fsrPhotonIdx);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaId", Muon_mvaId, &b_Muon_mvaId);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_mvaIDV1", Photon_mvaIDV1, &b_Photon_mvaIDV1);
   fChain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
   fChain->SetBranchAddress("Photon_cutBasedBitmap", Photon_cutBasedBitmap, &b_Photon_cutBasedBitmap);
   fChain->SetBranchAddress("Photon_cutBasedV1Bitmap", Photon_cutBasedV1Bitmap, &b_Photon_cutBasedV1Bitmap);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   fChain->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
   fChain->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
   fChain->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
   fChain->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
   fChain->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
   fChain->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
   fChain->SetBranchAddress("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc, &b_GenDressedLepton_hasTauAnc);
   fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   fChain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   fChain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   fChain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   fChain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   fChain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   fChain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   fChain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   fChain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   fChain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   fChain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   fChain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   fChain->SetBranchAddress("SubJet_btagCMVA", SubJet_btagCMVA, &b_SubJet_btagCMVA);
   fChain->SetBranchAddress("SubJet_btagCSVV2", SubJet_btagCSVV2, &b_SubJet_btagCSVV2);
   fChain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   fChain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   fChain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   fChain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   fChain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   fChain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   fChain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   fChain->SetBranchAddress("SubJet_rawFactor", SubJet_rawFactor, &b_SubJet_rawFactor);
   fChain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   fChain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   fChain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   fChain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("Tau_rawAntiEle", Tau_rawAntiEle, &b_Tau_rawAntiEle);
   fChain->SetBranchAddress("Tau_rawAntiEle2018", Tau_rawAntiEle2018, &b_Tau_rawAntiEle2018);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSe", Tau_rawDeepTau2017v2p1VSe, &b_Tau_rawDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSjet", Tau_rawDeepTau2017v2p1VSjet, &b_Tau_rawDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSmu", Tau_rawDeepTau2017v2p1VSmu, &b_Tau_rawDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   fChain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   fChain->SetBranchAddress("Tau_rawMVAnewDM2017v2", Tau_rawMVAnewDM2017v2, &b_Tau_rawMVAnewDM2017v2);
   fChain->SetBranchAddress("Tau_rawMVAoldDM", Tau_rawMVAoldDM, &b_Tau_rawMVAoldDM);
   fChain->SetBranchAddress("Tau_rawMVAoldDM2017v1", Tau_rawMVAoldDM2017v1, &b_Tau_rawMVAoldDM2017v1);
   fChain->SetBranchAddress("Tau_rawMVAoldDM2017v2", Tau_rawMVAoldDM2017v2, &b_Tau_rawMVAoldDM2017v2);
   fChain->SetBranchAddress("Tau_rawMVAoldDMdR032017v2", Tau_rawMVAoldDMdR032017v2, &b_Tau_rawMVAoldDMdR032017v2);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   fChain->SetBranchAddress("Tau_rawAntiEleCat", Tau_rawAntiEleCat, &b_Tau_rawAntiEleCat);
   fChain->SetBranchAddress("Tau_rawAntiEleCat2018", Tau_rawAntiEleCat2018, &b_Tau_rawAntiEleCat2018);
   fChain->SetBranchAddress("Tau_idAntiEle", Tau_idAntiEle, &b_Tau_idAntiEle);
   fChain->SetBranchAddress("Tau_idAntiEle2018", Tau_idAntiEle2018, &b_Tau_idAntiEle2018);
   fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   fChain->SetBranchAddress("Tau_idDecayMode", Tau_idDecayMode, &b_Tau_idDecayMode);
   fChain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSe", Tau_idDeepTau2017v2p1VSe, &b_Tau_idDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSjet", Tau_idDeepTau2017v2p1VSjet, &b_Tau_idDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSmu", Tau_idDeepTau2017v2p1VSmu, &b_Tau_idDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_idMVAnewDM2017v2", Tau_idMVAnewDM2017v2, &b_Tau_idMVAnewDM2017v2);
   fChain->SetBranchAddress("Tau_idMVAoldDM", Tau_idMVAoldDM, &b_Tau_idMVAoldDM);
   fChain->SetBranchAddress("Tau_idMVAoldDM2017v1", Tau_idMVAoldDM2017v1, &b_Tau_idMVAoldDM2017v1);
   fChain->SetBranchAddress("Tau_idMVAoldDM2017v2", Tau_idMVAoldDM2017v2, &b_Tau_idMVAoldDM2017v2);
   fChain->SetBranchAddress("Tau_idMVAoldDMdR032017v2", Tau_idMVAoldDMdR032017v2, &b_Tau_idMVAoldDMdR032017v2);
   fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   fChain->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
   fChain->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
   fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
   fChain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
   fChain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
   fChain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
   fChain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
   fChain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
   fChain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
   fChain->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask);
   fChain->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
   fChain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
   fChain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW, &b_HLT_DoubleEle25_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW, &b_HLT_DoubleEle27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf, &b_HLT_DoubleEle24_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_Ele37_CaloIdL_MW", &HLT_Ele27_Ele37_CaloIdL_MW, &b_HLT_Ele27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27, &b_HLT_Mu37_TkMu27);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi", &HLT_DoubleMu4_3_Jpsi, &b_HLT_DoubleMu4_3_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi, &b_HLT_Mu7p5_Track2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi, &b_HLT_Mu7p5_Track3p5_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi, &b_HLT_Mu7p5_Track7_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon, &b_HLT_Mu7p5_Track2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon, &b_HLT_Mu7p5_Track3p5_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon, &b_HLT_Mu7p5_Track7_Upsilon);
   fChain->SetBranchAddress("HLT_Mu3_L1SingleMu5orSingleMu7", &HLT_Mu3_L1SingleMu5orSingleMu7, &b_HLT_Mu3_L1SingleMu5orSingleMu7);
   fChain->SetBranchAddress("HLT_Ele20_WPTight_Gsf", &HLT_Ele20_WPTight_Gsf, &b_HLT_Ele20_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele15_WPLoose_Gsf", &HLT_Ele15_WPLoose_Gsf, &b_HLT_Ele15_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele17_WPLoose_Gsf", &HLT_Ele17_WPLoose_Gsf, &b_HLT_Ele17_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele20_WPLoose_Gsf", &HLT_Ele20_WPLoose_Gsf, &b_HLT_Ele20_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf", &HLT_Ele20_eta2p1_WPLoose_Gsf, &b_HLT_Ele20_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG, &b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele28_WPTight_Gsf", &HLT_Ele28_WPTight_Gsf, &b_HLT_Ele28_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf, &b_HLT_Ele30_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_L1EGMT", &HLT_Ele35_WPTight_Gsf_L1EGMT, &b_HLT_Ele35_WPTight_Gsf_L1EGMT);
   fChain->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, &b_HLT_Ele38_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, &b_HLT_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_IsoMu30", &HLT_IsoMu30, &b_HLT_IsoMu30);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Onia", &HLT_Mu25_TkMu0_Onia, &b_HLT_Mu25_TkMu0_Onia);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Psi", &HLT_Mu30_TkMu0_Psi, &b_HLT_Mu30_TkMu0_Psi);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Upsilon", &HLT_Mu30_TkMu0_Upsilon, &b_HLT_Mu30_TkMu0_Upsilon);
   fChain->SetBranchAddress("HLT_Mu20_TkMu0_Phi", &HLT_Mu20_TkMu0_Phi, &b_HLT_Mu20_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi, &b_HLT_Mu25_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu12", &HLT_Mu12, &b_HLT_Mu12);
   fChain->SetBranchAddress("HLT_Mu15", &HLT_Mu15, &b_HLT_Mu15);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_OldMu100", &HLT_OldMu100, &b_HLT_OldMu100);
   fChain->SetBranchAddress("HLT_TkMu100", &HLT_TkMu100, &b_HLT_TkMu100);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL", &HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu19", &HLT_Mu19, &b_HLT_Mu19);
   fChain->SetBranchAddress("HLT_Mu18_Mu9", &HLT_Mu18_Mu9, &b_HLT_Mu18_Mu9);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_DZ", &HLT_Mu18_Mu9_DZ, &b_HLT_Mu18_Mu9_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10, &b_HLT_Mu20_Mu10);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ, &b_HLT_Mu20_Mu10_DZ);
   fChain->SetBranchAddress("HLT_Mu23_Mu12", &HLT_Mu23_Mu12, &b_HLT_Mu23_Mu12);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_DZ", &HLT_Mu23_Mu12_DZ, &b_HLT_Mu23_Mu12_DZ);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("Jet_pt_raw", Jet_pt_raw, &b_Jet_pt_raw);
   fChain->SetBranchAddress("Jet_pt_nom", Jet_pt_nom, &b_Jet_pt_nom);
   fChain->SetBranchAddress("Jet_mass_raw", Jet_mass_raw, &b_Jet_mass_raw);
   fChain->SetBranchAddress("Jet_mass_nom", Jet_mass_nom, &b_Jet_mass_nom);
   fChain->SetBranchAddress("Jet_corr_JEC", Jet_corr_JEC, &b_Jet_corr_JEC);
   fChain->SetBranchAddress("Jet_corr_JER", Jet_corr_JER, &b_Jet_corr_JER);
   fChain->SetBranchAddress("MET_pt_nom", &MET_pt_nom, &b_MET_pt_nom);
   fChain->SetBranchAddress("MET_phi_nom", &MET_phi_nom, &b_MET_phi_nom);
   fChain->SetBranchAddress("MET_pt_jer", &MET_pt_jer, &b_MET_pt_jer);
   fChain->SetBranchAddress("MET_phi_jer", &MET_phi_jer, &b_MET_phi_jer);
   fChain->SetBranchAddress("Jet_pt_jerUp", Jet_pt_jerUp, &b_Jet_pt_jerUp);
   fChain->SetBranchAddress("Jet_mass_jerUp", Jet_mass_jerUp, &b_Jet_mass_jerUp);
   fChain->SetBranchAddress("MET_pt_jerUp", &MET_pt_jerUp, &b_MET_pt_jerUp);
   fChain->SetBranchAddress("MET_phi_jerUp", &MET_phi_jerUp, &b_MET_phi_jerUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteStatUp", Jet_pt_jesAbsoluteStatUp, &b_Jet_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteStatUp", Jet_mass_jesAbsoluteStatUp, &b_Jet_mass_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteStatUp", &MET_pt_jesAbsoluteStatUp, &b_MET_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteStatUp", &MET_phi_jesAbsoluteStatUp, &b_MET_phi_jesAbsoluteStatUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteScaleUp", Jet_pt_jesAbsoluteScaleUp, &b_Jet_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteScaleUp", Jet_mass_jesAbsoluteScaleUp, &b_Jet_mass_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteScaleUp", &MET_pt_jesAbsoluteScaleUp, &b_MET_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteScaleUp", &MET_phi_jesAbsoluteScaleUp, &b_MET_phi_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteSampleUp", Jet_pt_jesAbsoluteSampleUp, &b_Jet_pt_jesAbsoluteSampleUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteSampleUp", Jet_mass_jesAbsoluteSampleUp, &b_Jet_mass_jesAbsoluteSampleUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteSampleUp", &MET_pt_jesAbsoluteSampleUp, &b_MET_pt_jesAbsoluteSampleUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteSampleUp", &MET_phi_jesAbsoluteSampleUp, &b_MET_phi_jesAbsoluteSampleUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapUp", Jet_pt_jesAbsoluteFlavMapUp, &b_Jet_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapUp", Jet_mass_jesAbsoluteFlavMapUp, &b_Jet_mass_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteFlavMapUp", &MET_pt_jesAbsoluteFlavMapUp, &b_MET_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteFlavMapUp", &MET_phi_jesAbsoluteFlavMapUp, &b_MET_phi_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasUp", Jet_pt_jesAbsoluteMPFBiasUp, &b_Jet_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasUp", Jet_mass_jesAbsoluteMPFBiasUp, &b_Jet_mass_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteMPFBiasUp", &MET_pt_jesAbsoluteMPFBiasUp, &b_MET_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteMPFBiasUp", &MET_phi_jesAbsoluteMPFBiasUp, &b_MET_phi_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("Jet_pt_jesFragmentationUp", Jet_pt_jesFragmentationUp, &b_Jet_pt_jesFragmentationUp);
   fChain->SetBranchAddress("Jet_mass_jesFragmentationUp", Jet_mass_jesFragmentationUp, &b_Jet_mass_jesFragmentationUp);
   fChain->SetBranchAddress("MET_pt_jesFragmentationUp", &MET_pt_jesFragmentationUp, &b_MET_pt_jesFragmentationUp);
   fChain->SetBranchAddress("MET_phi_jesFragmentationUp", &MET_phi_jesFragmentationUp, &b_MET_phi_jesFragmentationUp);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionECALUp", Jet_pt_jesSinglePionECALUp, &b_Jet_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionECALUp", Jet_mass_jesSinglePionECALUp, &b_Jet_mass_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_pt_jesSinglePionECALUp", &MET_pt_jesSinglePionECALUp, &b_MET_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_phi_jesSinglePionECALUp", &MET_phi_jesSinglePionECALUp, &b_MET_phi_jesSinglePionECALUp);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionHCALUp", Jet_pt_jesSinglePionHCALUp, &b_Jet_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionHCALUp", Jet_mass_jesSinglePionHCALUp, &b_Jet_mass_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_pt_jesSinglePionHCALUp", &MET_pt_jesSinglePionHCALUp, &b_MET_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_phi_jesSinglePionHCALUp", &MET_phi_jesSinglePionHCALUp, &b_MET_phi_jesSinglePionHCALUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorQCDUp", Jet_pt_jesFlavorQCDUp, &b_Jet_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorQCDUp", Jet_mass_jesFlavorQCDUp, &b_Jet_mass_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorQCDUp", &MET_pt_jesFlavorQCDUp, &b_MET_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorQCDUp", &MET_phi_jesFlavorQCDUp, &b_MET_phi_jesFlavorQCDUp);
   fChain->SetBranchAddress("Jet_pt_jesTimePtEtaUp", Jet_pt_jesTimePtEtaUp, &b_Jet_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("Jet_mass_jesTimePtEtaUp", Jet_mass_jesTimePtEtaUp, &b_Jet_mass_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_pt_jesTimePtEtaUp", &MET_pt_jesTimePtEtaUp, &b_MET_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_phi_jesTimePtEtaUp", &MET_phi_jesTimePtEtaUp, &b_MET_phi_jesTimePtEtaUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC1Up", Jet_pt_jesRelativeJEREC1Up, &b_Jet_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC1Up", Jet_mass_jesRelativeJEREC1Up, &b_Jet_mass_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC1Up", &MET_pt_jesRelativeJEREC1Up, &b_MET_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC1Up", &MET_phi_jesRelativeJEREC1Up, &b_MET_phi_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC2Up", Jet_pt_jesRelativeJEREC2Up, &b_Jet_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC2Up", Jet_mass_jesRelativeJEREC2Up, &b_Jet_mass_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC2Up", &MET_pt_jesRelativeJEREC2Up, &b_MET_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC2Up", &MET_phi_jesRelativeJEREC2Up, &b_MET_phi_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJERHFUp", Jet_pt_jesRelativeJERHFUp, &b_Jet_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJERHFUp", Jet_mass_jesRelativeJERHFUp, &b_Jet_mass_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeJERHFUp", &MET_pt_jesRelativeJERHFUp, &b_MET_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeJERHFUp", &MET_phi_jesRelativeJERHFUp, &b_MET_phi_jesRelativeJERHFUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtBBUp", Jet_pt_jesRelativePtBBUp, &b_Jet_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtBBUp", Jet_mass_jesRelativePtBBUp, &b_Jet_mass_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_pt_jesRelativePtBBUp", &MET_pt_jesRelativePtBBUp, &b_MET_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_phi_jesRelativePtBBUp", &MET_phi_jesRelativePtBBUp, &b_MET_phi_jesRelativePtBBUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC1Up", Jet_pt_jesRelativePtEC1Up, &b_Jet_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC1Up", Jet_mass_jesRelativePtEC1Up, &b_Jet_mass_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC1Up", &MET_pt_jesRelativePtEC1Up, &b_MET_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC1Up", &MET_phi_jesRelativePtEC1Up, &b_MET_phi_jesRelativePtEC1Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC2Up", Jet_pt_jesRelativePtEC2Up, &b_Jet_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC2Up", Jet_mass_jesRelativePtEC2Up, &b_Jet_mass_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC2Up", &MET_pt_jesRelativePtEC2Up, &b_MET_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC2Up", &MET_phi_jesRelativePtEC2Up, &b_MET_phi_jesRelativePtEC2Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtHFUp", Jet_pt_jesRelativePtHFUp, &b_Jet_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtHFUp", Jet_mass_jesRelativePtHFUp, &b_Jet_mass_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_pt_jesRelativePtHFUp", &MET_pt_jesRelativePtHFUp, &b_MET_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_phi_jesRelativePtHFUp", &MET_phi_jesRelativePtHFUp, &b_MET_phi_jesRelativePtHFUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeBalUp", Jet_pt_jesRelativeBalUp, &b_Jet_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeBalUp", Jet_mass_jesRelativeBalUp, &b_Jet_mass_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeBalUp", &MET_pt_jesRelativeBalUp, &b_MET_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeBalUp", &MET_phi_jesRelativeBalUp, &b_MET_phi_jesRelativeBalUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeSampleUp", Jet_pt_jesRelativeSampleUp, &b_Jet_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeSampleUp", Jet_mass_jesRelativeSampleUp, &b_Jet_mass_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeSampleUp", &MET_pt_jesRelativeSampleUp, &b_MET_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeSampleUp", &MET_phi_jesRelativeSampleUp, &b_MET_phi_jesRelativeSampleUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeFSRUp", Jet_pt_jesRelativeFSRUp, &b_Jet_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeFSRUp", Jet_mass_jesRelativeFSRUp, &b_Jet_mass_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeFSRUp", &MET_pt_jesRelativeFSRUp, &b_MET_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeFSRUp", &MET_phi_jesRelativeFSRUp, &b_MET_phi_jesRelativeFSRUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatFSRUp", Jet_pt_jesRelativeStatFSRUp, &b_Jet_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatFSRUp", Jet_mass_jesRelativeStatFSRUp, &b_Jet_mass_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatFSRUp", &MET_pt_jesRelativeStatFSRUp, &b_MET_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatFSRUp", &MET_phi_jesRelativeStatFSRUp, &b_MET_phi_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatECUp", Jet_pt_jesRelativeStatECUp, &b_Jet_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatECUp", Jet_mass_jesRelativeStatECUp, &b_Jet_mass_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatECUp", &MET_pt_jesRelativeStatECUp, &b_MET_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatECUp", &MET_phi_jesRelativeStatECUp, &b_MET_phi_jesRelativeStatECUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatHFUp", Jet_pt_jesRelativeStatHFUp, &b_Jet_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatHFUp", Jet_mass_jesRelativeStatHFUp, &b_Jet_mass_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatHFUp", &MET_pt_jesRelativeStatHFUp, &b_MET_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatHFUp", &MET_phi_jesRelativeStatHFUp, &b_MET_phi_jesRelativeStatHFUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpDataMCUp", Jet_pt_jesPileUpDataMCUp, &b_Jet_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpDataMCUp", Jet_mass_jesPileUpDataMCUp, &b_Jet_mass_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpDataMCUp", &MET_pt_jesPileUpDataMCUp, &b_MET_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpDataMCUp", &MET_phi_jesPileUpDataMCUp, &b_MET_phi_jesPileUpDataMCUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtRefUp", Jet_pt_jesPileUpPtRefUp, &b_Jet_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtRefUp", Jet_mass_jesPileUpPtRefUp, &b_Jet_mass_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtRefUp", &MET_pt_jesPileUpPtRefUp, &b_MET_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtRefUp", &MET_phi_jesPileUpPtRefUp, &b_MET_phi_jesPileUpPtRefUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtBBUp", Jet_pt_jesPileUpPtBBUp, &b_Jet_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtBBUp", Jet_mass_jesPileUpPtBBUp, &b_Jet_mass_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtBBUp", &MET_pt_jesPileUpPtBBUp, &b_MET_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtBBUp", &MET_phi_jesPileUpPtBBUp, &b_MET_phi_jesPileUpPtBBUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC1Up", Jet_pt_jesPileUpPtEC1Up, &b_Jet_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC1Up", Jet_mass_jesPileUpPtEC1Up, &b_Jet_mass_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC1Up", &MET_pt_jesPileUpPtEC1Up, &b_MET_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC1Up", &MET_phi_jesPileUpPtEC1Up, &b_MET_phi_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC2Up", Jet_pt_jesPileUpPtEC2Up, &b_Jet_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC2Up", Jet_mass_jesPileUpPtEC2Up, &b_Jet_mass_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC2Up", &MET_pt_jesPileUpPtEC2Up, &b_MET_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC2Up", &MET_phi_jesPileUpPtEC2Up, &b_MET_phi_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtHFUp", Jet_pt_jesPileUpPtHFUp, &b_Jet_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtHFUp", Jet_mass_jesPileUpPtHFUp, &b_Jet_mass_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtHFUp", &MET_pt_jesPileUpPtHFUp, &b_MET_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtHFUp", &MET_phi_jesPileUpPtHFUp, &b_MET_phi_jesPileUpPtHFUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpMuZeroUp", Jet_pt_jesPileUpMuZeroUp, &b_Jet_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpMuZeroUp", Jet_mass_jesPileUpMuZeroUp, &b_Jet_mass_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpMuZeroUp", &MET_pt_jesPileUpMuZeroUp, &b_MET_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpMuZeroUp", &MET_phi_jesPileUpMuZeroUp, &b_MET_phi_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpEnvelopeUp", Jet_pt_jesPileUpEnvelopeUp, &b_Jet_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpEnvelopeUp", Jet_mass_jesPileUpEnvelopeUp, &b_Jet_mass_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_pt_jesPileUpEnvelopeUp", &MET_pt_jesPileUpEnvelopeUp, &b_MET_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_phi_jesPileUpEnvelopeUp", &MET_phi_jesPileUpEnvelopeUp, &b_MET_phi_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPileUpUp", Jet_pt_jesSubTotalPileUpUp, &b_Jet_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPileUpUp", Jet_mass_jesSubTotalPileUpUp, &b_Jet_mass_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPileUpUp", &MET_pt_jesSubTotalPileUpUp, &b_MET_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPileUpUp", &MET_phi_jesSubTotalPileUpUp, &b_MET_phi_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalRelativeUp", Jet_pt_jesSubTotalRelativeUp, &b_Jet_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalRelativeUp", Jet_mass_jesSubTotalRelativeUp, &b_Jet_mass_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalRelativeUp", &MET_pt_jesSubTotalRelativeUp, &b_MET_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalRelativeUp", &MET_phi_jesSubTotalRelativeUp, &b_MET_phi_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPtUp", Jet_pt_jesSubTotalPtUp, &b_Jet_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPtUp", Jet_mass_jesSubTotalPtUp, &b_Jet_mass_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPtUp", &MET_pt_jesSubTotalPtUp, &b_MET_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPtUp", &MET_phi_jesSubTotalPtUp, &b_MET_phi_jesSubTotalPtUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalScaleUp", Jet_pt_jesSubTotalScaleUp, &b_Jet_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalScaleUp", Jet_mass_jesSubTotalScaleUp, &b_Jet_mass_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalScaleUp", &MET_pt_jesSubTotalScaleUp, &b_MET_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalScaleUp", &MET_phi_jesSubTotalScaleUp, &b_MET_phi_jesSubTotalScaleUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteUp", Jet_pt_jesSubTotalAbsoluteUp, &b_Jet_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteUp", Jet_mass_jesSubTotalAbsoluteUp, &b_Jet_mass_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalAbsoluteUp", &MET_pt_jesSubTotalAbsoluteUp, &b_MET_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalAbsoluteUp", &MET_phi_jesSubTotalAbsoluteUp, &b_MET_phi_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalMCUp", Jet_pt_jesSubTotalMCUp, &b_Jet_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalMCUp", Jet_mass_jesSubTotalMCUp, &b_Jet_mass_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_pt_jesSubTotalMCUp", &MET_pt_jesSubTotalMCUp, &b_MET_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_phi_jesSubTotalMCUp", &MET_phi_jesSubTotalMCUp, &b_MET_phi_jesSubTotalMCUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalUp", Jet_pt_jesTotalUp, &b_Jet_pt_jesTotalUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalUp", Jet_mass_jesTotalUp, &b_Jet_mass_jesTotalUp);
   fChain->SetBranchAddress("MET_pt_jesTotalUp", &MET_pt_jesTotalUp, &b_MET_pt_jesTotalUp);
   fChain->SetBranchAddress("MET_phi_jesTotalUp", &MET_phi_jesTotalUp, &b_MET_phi_jesTotalUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorUp", Jet_pt_jesTotalNoFlavorUp, &b_Jet_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorUp", Jet_mass_jesTotalNoFlavorUp, &b_Jet_mass_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorUp", &MET_pt_jesTotalNoFlavorUp, &b_MET_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorUp", &MET_phi_jesTotalNoFlavorUp, &b_MET_phi_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoTimeUp", Jet_pt_jesTotalNoTimeUp, &b_Jet_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoTimeUp", Jet_mass_jesTotalNoTimeUp, &b_Jet_mass_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_pt_jesTotalNoTimeUp", &MET_pt_jesTotalNoTimeUp, &b_MET_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_phi_jesTotalNoTimeUp", &MET_phi_jesTotalNoTimeUp, &b_MET_phi_jesTotalNoTimeUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeUp", Jet_pt_jesTotalNoFlavorNoTimeUp, &b_Jet_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeUp", Jet_mass_jesTotalNoFlavorNoTimeUp, &b_Jet_mass_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorNoTimeUp", &MET_pt_jesTotalNoFlavorNoTimeUp, &b_MET_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorNoTimeUp", &MET_phi_jesTotalNoFlavorNoTimeUp, &b_MET_phi_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorZJetUp", Jet_pt_jesFlavorZJetUp, &b_Jet_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorZJetUp", Jet_mass_jesFlavorZJetUp, &b_Jet_mass_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorZJetUp", &MET_pt_jesFlavorZJetUp, &b_MET_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorZJetUp", &MET_phi_jesFlavorZJetUp, &b_MET_phi_jesFlavorZJetUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPhotonJetUp", Jet_pt_jesFlavorPhotonJetUp, &b_Jet_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPhotonJetUp", Jet_mass_jesFlavorPhotonJetUp, &b_Jet_mass_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPhotonJetUp", &MET_pt_jesFlavorPhotonJetUp, &b_MET_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPhotonJetUp", &MET_phi_jesFlavorPhotonJetUp, &b_MET_phi_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureGluonUp", Jet_pt_jesFlavorPureGluonUp, &b_Jet_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureGluonUp", Jet_mass_jesFlavorPureGluonUp, &b_Jet_mass_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureGluonUp", &MET_pt_jesFlavorPureGluonUp, &b_MET_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureGluonUp", &MET_phi_jesFlavorPureGluonUp, &b_MET_phi_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureQuarkUp", Jet_pt_jesFlavorPureQuarkUp, &b_Jet_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureQuarkUp", Jet_mass_jesFlavorPureQuarkUp, &b_Jet_mass_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureQuarkUp", &MET_pt_jesFlavorPureQuarkUp, &b_MET_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureQuarkUp", &MET_phi_jesFlavorPureQuarkUp, &b_MET_phi_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureCharmUp", Jet_pt_jesFlavorPureCharmUp, &b_Jet_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureCharmUp", Jet_mass_jesFlavorPureCharmUp, &b_Jet_mass_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureCharmUp", &MET_pt_jesFlavorPureCharmUp, &b_MET_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureCharmUp", &MET_phi_jesFlavorPureCharmUp, &b_MET_phi_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureBottomUp", Jet_pt_jesFlavorPureBottomUp, &b_Jet_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureBottomUp", Jet_mass_jesFlavorPureBottomUp, &b_Jet_mass_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureBottomUp", &MET_pt_jesFlavorPureBottomUp, &b_MET_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureBottomUp", &MET_phi_jesFlavorPureBottomUp, &b_MET_phi_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunAUp", Jet_pt_jesTimeRunAUp, &b_Jet_pt_jesTimeRunAUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunAUp", Jet_mass_jesTimeRunAUp, &b_Jet_mass_jesTimeRunAUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunAUp", &MET_pt_jesTimeRunAUp, &b_MET_pt_jesTimeRunAUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunAUp", &MET_phi_jesTimeRunAUp, &b_MET_phi_jesTimeRunAUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunBUp", Jet_pt_jesTimeRunBUp, &b_Jet_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunBUp", Jet_mass_jesTimeRunBUp, &b_Jet_mass_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunBUp", &MET_pt_jesTimeRunBUp, &b_MET_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunBUp", &MET_phi_jesTimeRunBUp, &b_MET_phi_jesTimeRunBUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunCUp", Jet_pt_jesTimeRunCUp, &b_Jet_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunCUp", Jet_mass_jesTimeRunCUp, &b_Jet_mass_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunCUp", &MET_pt_jesTimeRunCUp, &b_MET_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunCUp", &MET_phi_jesTimeRunCUp, &b_MET_phi_jesTimeRunCUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunDUp", Jet_pt_jesTimeRunDUp, &b_Jet_pt_jesTimeRunDUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunDUp", Jet_mass_jesTimeRunDUp, &b_Jet_mass_jesTimeRunDUp);
   fChain->SetBranchAddress("MET_pt_jesTimeRunDUp", &MET_pt_jesTimeRunDUp, &b_MET_pt_jesTimeRunDUp);
   fChain->SetBranchAddress("MET_phi_jesTimeRunDUp", &MET_phi_jesTimeRunDUp, &b_MET_phi_jesTimeRunDUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituUp", Jet_pt_jesCorrelationGroupMPFInSituUp, &b_Jet_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituUp", Jet_mass_jesCorrelationGroupMPFInSituUp, &b_Jet_mass_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupMPFInSituUp", &MET_pt_jesCorrelationGroupMPFInSituUp, &b_MET_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupMPFInSituUp", &MET_phi_jesCorrelationGroupMPFInSituUp, &b_MET_phi_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationUp", Jet_pt_jesCorrelationGroupIntercalibrationUp, &b_Jet_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationUp", Jet_mass_jesCorrelationGroupIntercalibrationUp, &b_Jet_mass_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupIntercalibrationUp", &MET_pt_jesCorrelationGroupIntercalibrationUp, &b_MET_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupIntercalibrationUp", &MET_phi_jesCorrelationGroupIntercalibrationUp, &b_MET_phi_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESUp", Jet_pt_jesCorrelationGroupbJESUp, &b_Jet_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESUp", Jet_mass_jesCorrelationGroupbJESUp, &b_Jet_mass_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupbJESUp", &MET_pt_jesCorrelationGroupbJESUp, &b_MET_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupbJESUp", &MET_phi_jesCorrelationGroupbJESUp, &b_MET_phi_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorUp", Jet_pt_jesCorrelationGroupFlavorUp, &b_Jet_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorUp", Jet_mass_jesCorrelationGroupFlavorUp, &b_Jet_mass_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupFlavorUp", &MET_pt_jesCorrelationGroupFlavorUp, &b_MET_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupFlavorUp", &MET_phi_jesCorrelationGroupFlavorUp, &b_MET_phi_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedUp", Jet_pt_jesCorrelationGroupUncorrelatedUp, &b_Jet_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedUp", Jet_mass_jesCorrelationGroupUncorrelatedUp, &b_Jet_mass_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupUncorrelatedUp", &MET_pt_jesCorrelationGroupUncorrelatedUp, &b_MET_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupUncorrelatedUp", &MET_phi_jesCorrelationGroupUncorrelatedUp, &b_MET_phi_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_pt_unclustEnUp", &MET_pt_unclustEnUp, &b_MET_pt_unclustEnUp);
   fChain->SetBranchAddress("MET_phi_unclustEnUp", &MET_phi_unclustEnUp, &b_MET_phi_unclustEnUp);
   fChain->SetBranchAddress("Jet_pt_jerDown", Jet_pt_jerDown, &b_Jet_pt_jerDown);
   fChain->SetBranchAddress("Jet_mass_jerDown", Jet_mass_jerDown, &b_Jet_mass_jerDown);
   fChain->SetBranchAddress("MET_pt_jerDown", &MET_pt_jerDown, &b_MET_pt_jerDown);
   fChain->SetBranchAddress("MET_phi_jerDown", &MET_phi_jerDown, &b_MET_phi_jerDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteStatDown", Jet_pt_jesAbsoluteStatDown, &b_Jet_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteStatDown", Jet_mass_jesAbsoluteStatDown, &b_Jet_mass_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteStatDown", &MET_pt_jesAbsoluteStatDown, &b_MET_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteStatDown", &MET_phi_jesAbsoluteStatDown, &b_MET_phi_jesAbsoluteStatDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteScaleDown", Jet_pt_jesAbsoluteScaleDown, &b_Jet_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteScaleDown", Jet_mass_jesAbsoluteScaleDown, &b_Jet_mass_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteScaleDown", &MET_pt_jesAbsoluteScaleDown, &b_MET_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteScaleDown", &MET_phi_jesAbsoluteScaleDown, &b_MET_phi_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteSampleDown", Jet_pt_jesAbsoluteSampleDown, &b_Jet_pt_jesAbsoluteSampleDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteSampleDown", Jet_mass_jesAbsoluteSampleDown, &b_Jet_mass_jesAbsoluteSampleDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteSampleDown", &MET_pt_jesAbsoluteSampleDown, &b_MET_pt_jesAbsoluteSampleDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteSampleDown", &MET_phi_jesAbsoluteSampleDown, &b_MET_phi_jesAbsoluteSampleDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapDown", Jet_pt_jesAbsoluteFlavMapDown, &b_Jet_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapDown", Jet_mass_jesAbsoluteFlavMapDown, &b_Jet_mass_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteFlavMapDown", &MET_pt_jesAbsoluteFlavMapDown, &b_MET_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteFlavMapDown", &MET_phi_jesAbsoluteFlavMapDown, &b_MET_phi_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasDown", Jet_pt_jesAbsoluteMPFBiasDown, &b_Jet_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasDown", Jet_mass_jesAbsoluteMPFBiasDown, &b_Jet_mass_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_pt_jesAbsoluteMPFBiasDown", &MET_pt_jesAbsoluteMPFBiasDown, &b_MET_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_phi_jesAbsoluteMPFBiasDown", &MET_phi_jesAbsoluteMPFBiasDown, &b_MET_phi_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("Jet_pt_jesFragmentationDown", Jet_pt_jesFragmentationDown, &b_Jet_pt_jesFragmentationDown);
   fChain->SetBranchAddress("Jet_mass_jesFragmentationDown", Jet_mass_jesFragmentationDown, &b_Jet_mass_jesFragmentationDown);
   fChain->SetBranchAddress("MET_pt_jesFragmentationDown", &MET_pt_jesFragmentationDown, &b_MET_pt_jesFragmentationDown);
   fChain->SetBranchAddress("MET_phi_jesFragmentationDown", &MET_phi_jesFragmentationDown, &b_MET_phi_jesFragmentationDown);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionECALDown", Jet_pt_jesSinglePionECALDown, &b_Jet_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionECALDown", Jet_mass_jesSinglePionECALDown, &b_Jet_mass_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_pt_jesSinglePionECALDown", &MET_pt_jesSinglePionECALDown, &b_MET_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_phi_jesSinglePionECALDown", &MET_phi_jesSinglePionECALDown, &b_MET_phi_jesSinglePionECALDown);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionHCALDown", Jet_pt_jesSinglePionHCALDown, &b_Jet_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionHCALDown", Jet_mass_jesSinglePionHCALDown, &b_Jet_mass_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_pt_jesSinglePionHCALDown", &MET_pt_jesSinglePionHCALDown, &b_MET_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_phi_jesSinglePionHCALDown", &MET_phi_jesSinglePionHCALDown, &b_MET_phi_jesSinglePionHCALDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorQCDDown", Jet_pt_jesFlavorQCDDown, &b_Jet_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorQCDDown", Jet_mass_jesFlavorQCDDown, &b_Jet_mass_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorQCDDown", &MET_pt_jesFlavorQCDDown, &b_MET_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorQCDDown", &MET_phi_jesFlavorQCDDown, &b_MET_phi_jesFlavorQCDDown);
   fChain->SetBranchAddress("Jet_pt_jesTimePtEtaDown", Jet_pt_jesTimePtEtaDown, &b_Jet_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("Jet_mass_jesTimePtEtaDown", Jet_mass_jesTimePtEtaDown, &b_Jet_mass_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_pt_jesTimePtEtaDown", &MET_pt_jesTimePtEtaDown, &b_MET_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_phi_jesTimePtEtaDown", &MET_phi_jesTimePtEtaDown, &b_MET_phi_jesTimePtEtaDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC1Down", Jet_pt_jesRelativeJEREC1Down, &b_Jet_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC1Down", Jet_mass_jesRelativeJEREC1Down, &b_Jet_mass_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC1Down", &MET_pt_jesRelativeJEREC1Down, &b_MET_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC1Down", &MET_phi_jesRelativeJEREC1Down, &b_MET_phi_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC2Down", Jet_pt_jesRelativeJEREC2Down, &b_Jet_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC2Down", Jet_mass_jesRelativeJEREC2Down, &b_Jet_mass_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_pt_jesRelativeJEREC2Down", &MET_pt_jesRelativeJEREC2Down, &b_MET_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_phi_jesRelativeJEREC2Down", &MET_phi_jesRelativeJEREC2Down, &b_MET_phi_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJERHFDown", Jet_pt_jesRelativeJERHFDown, &b_Jet_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJERHFDown", Jet_mass_jesRelativeJERHFDown, &b_Jet_mass_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeJERHFDown", &MET_pt_jesRelativeJERHFDown, &b_MET_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeJERHFDown", &MET_phi_jesRelativeJERHFDown, &b_MET_phi_jesRelativeJERHFDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtBBDown", Jet_pt_jesRelativePtBBDown, &b_Jet_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtBBDown", Jet_mass_jesRelativePtBBDown, &b_Jet_mass_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_pt_jesRelativePtBBDown", &MET_pt_jesRelativePtBBDown, &b_MET_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_phi_jesRelativePtBBDown", &MET_phi_jesRelativePtBBDown, &b_MET_phi_jesRelativePtBBDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC1Down", Jet_pt_jesRelativePtEC1Down, &b_Jet_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC1Down", Jet_mass_jesRelativePtEC1Down, &b_Jet_mass_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC1Down", &MET_pt_jesRelativePtEC1Down, &b_MET_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC1Down", &MET_phi_jesRelativePtEC1Down, &b_MET_phi_jesRelativePtEC1Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC2Down", Jet_pt_jesRelativePtEC2Down, &b_Jet_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC2Down", Jet_mass_jesRelativePtEC2Down, &b_Jet_mass_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_pt_jesRelativePtEC2Down", &MET_pt_jesRelativePtEC2Down, &b_MET_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_phi_jesRelativePtEC2Down", &MET_phi_jesRelativePtEC2Down, &b_MET_phi_jesRelativePtEC2Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtHFDown", Jet_pt_jesRelativePtHFDown, &b_Jet_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtHFDown", Jet_mass_jesRelativePtHFDown, &b_Jet_mass_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_pt_jesRelativePtHFDown", &MET_pt_jesRelativePtHFDown, &b_MET_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_phi_jesRelativePtHFDown", &MET_phi_jesRelativePtHFDown, &b_MET_phi_jesRelativePtHFDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeBalDown", Jet_pt_jesRelativeBalDown, &b_Jet_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeBalDown", Jet_mass_jesRelativeBalDown, &b_Jet_mass_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeBalDown", &MET_pt_jesRelativeBalDown, &b_MET_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeBalDown", &MET_phi_jesRelativeBalDown, &b_MET_phi_jesRelativeBalDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeSampleDown", Jet_pt_jesRelativeSampleDown, &b_Jet_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeSampleDown", Jet_mass_jesRelativeSampleDown, &b_Jet_mass_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeSampleDown", &MET_pt_jesRelativeSampleDown, &b_MET_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeSampleDown", &MET_phi_jesRelativeSampleDown, &b_MET_phi_jesRelativeSampleDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeFSRDown", Jet_pt_jesRelativeFSRDown, &b_Jet_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeFSRDown", Jet_mass_jesRelativeFSRDown, &b_Jet_mass_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeFSRDown", &MET_pt_jesRelativeFSRDown, &b_MET_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeFSRDown", &MET_phi_jesRelativeFSRDown, &b_MET_phi_jesRelativeFSRDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatFSRDown", Jet_pt_jesRelativeStatFSRDown, &b_Jet_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatFSRDown", Jet_mass_jesRelativeStatFSRDown, &b_Jet_mass_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatFSRDown", &MET_pt_jesRelativeStatFSRDown, &b_MET_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatFSRDown", &MET_phi_jesRelativeStatFSRDown, &b_MET_phi_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatECDown", Jet_pt_jesRelativeStatECDown, &b_Jet_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatECDown", Jet_mass_jesRelativeStatECDown, &b_Jet_mass_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatECDown", &MET_pt_jesRelativeStatECDown, &b_MET_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatECDown", &MET_phi_jesRelativeStatECDown, &b_MET_phi_jesRelativeStatECDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatHFDown", Jet_pt_jesRelativeStatHFDown, &b_Jet_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatHFDown", Jet_mass_jesRelativeStatHFDown, &b_Jet_mass_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_pt_jesRelativeStatHFDown", &MET_pt_jesRelativeStatHFDown, &b_MET_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_phi_jesRelativeStatHFDown", &MET_phi_jesRelativeStatHFDown, &b_MET_phi_jesRelativeStatHFDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpDataMCDown", Jet_pt_jesPileUpDataMCDown, &b_Jet_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpDataMCDown", Jet_mass_jesPileUpDataMCDown, &b_Jet_mass_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpDataMCDown", &MET_pt_jesPileUpDataMCDown, &b_MET_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpDataMCDown", &MET_phi_jesPileUpDataMCDown, &b_MET_phi_jesPileUpDataMCDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtRefDown", Jet_pt_jesPileUpPtRefDown, &b_Jet_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtRefDown", Jet_mass_jesPileUpPtRefDown, &b_Jet_mass_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtRefDown", &MET_pt_jesPileUpPtRefDown, &b_MET_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtRefDown", &MET_phi_jesPileUpPtRefDown, &b_MET_phi_jesPileUpPtRefDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtBBDown", Jet_pt_jesPileUpPtBBDown, &b_Jet_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtBBDown", Jet_mass_jesPileUpPtBBDown, &b_Jet_mass_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtBBDown", &MET_pt_jesPileUpPtBBDown, &b_MET_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtBBDown", &MET_phi_jesPileUpPtBBDown, &b_MET_phi_jesPileUpPtBBDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC1Down", Jet_pt_jesPileUpPtEC1Down, &b_Jet_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC1Down", Jet_mass_jesPileUpPtEC1Down, &b_Jet_mass_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC1Down", &MET_pt_jesPileUpPtEC1Down, &b_MET_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC1Down", &MET_phi_jesPileUpPtEC1Down, &b_MET_phi_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC2Down", Jet_pt_jesPileUpPtEC2Down, &b_Jet_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC2Down", Jet_mass_jesPileUpPtEC2Down, &b_Jet_mass_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtEC2Down", &MET_pt_jesPileUpPtEC2Down, &b_MET_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtEC2Down", &MET_phi_jesPileUpPtEC2Down, &b_MET_phi_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtHFDown", Jet_pt_jesPileUpPtHFDown, &b_Jet_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtHFDown", Jet_mass_jesPileUpPtHFDown, &b_Jet_mass_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpPtHFDown", &MET_pt_jesPileUpPtHFDown, &b_MET_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpPtHFDown", &MET_phi_jesPileUpPtHFDown, &b_MET_phi_jesPileUpPtHFDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpMuZeroDown", Jet_pt_jesPileUpMuZeroDown, &b_Jet_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpMuZeroDown", Jet_mass_jesPileUpMuZeroDown, &b_Jet_mass_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpMuZeroDown", &MET_pt_jesPileUpMuZeroDown, &b_MET_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpMuZeroDown", &MET_phi_jesPileUpMuZeroDown, &b_MET_phi_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpEnvelopeDown", Jet_pt_jesPileUpEnvelopeDown, &b_Jet_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpEnvelopeDown", Jet_mass_jesPileUpEnvelopeDown, &b_Jet_mass_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_pt_jesPileUpEnvelopeDown", &MET_pt_jesPileUpEnvelopeDown, &b_MET_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_phi_jesPileUpEnvelopeDown", &MET_phi_jesPileUpEnvelopeDown, &b_MET_phi_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPileUpDown", Jet_pt_jesSubTotalPileUpDown, &b_Jet_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPileUpDown", Jet_mass_jesSubTotalPileUpDown, &b_Jet_mass_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPileUpDown", &MET_pt_jesSubTotalPileUpDown, &b_MET_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPileUpDown", &MET_phi_jesSubTotalPileUpDown, &b_MET_phi_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalRelativeDown", Jet_pt_jesSubTotalRelativeDown, &b_Jet_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalRelativeDown", Jet_mass_jesSubTotalRelativeDown, &b_Jet_mass_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalRelativeDown", &MET_pt_jesSubTotalRelativeDown, &b_MET_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalRelativeDown", &MET_phi_jesSubTotalRelativeDown, &b_MET_phi_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPtDown", Jet_pt_jesSubTotalPtDown, &b_Jet_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPtDown", Jet_mass_jesSubTotalPtDown, &b_Jet_mass_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalPtDown", &MET_pt_jesSubTotalPtDown, &b_MET_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalPtDown", &MET_phi_jesSubTotalPtDown, &b_MET_phi_jesSubTotalPtDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalScaleDown", Jet_pt_jesSubTotalScaleDown, &b_Jet_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalScaleDown", Jet_mass_jesSubTotalScaleDown, &b_Jet_mass_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalScaleDown", &MET_pt_jesSubTotalScaleDown, &b_MET_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalScaleDown", &MET_phi_jesSubTotalScaleDown, &b_MET_phi_jesSubTotalScaleDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteDown", Jet_pt_jesSubTotalAbsoluteDown, &b_Jet_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteDown", Jet_mass_jesSubTotalAbsoluteDown, &b_Jet_mass_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalAbsoluteDown", &MET_pt_jesSubTotalAbsoluteDown, &b_MET_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalAbsoluteDown", &MET_phi_jesSubTotalAbsoluteDown, &b_MET_phi_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalMCDown", Jet_pt_jesSubTotalMCDown, &b_Jet_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalMCDown", Jet_mass_jesSubTotalMCDown, &b_Jet_mass_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_pt_jesSubTotalMCDown", &MET_pt_jesSubTotalMCDown, &b_MET_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_phi_jesSubTotalMCDown", &MET_phi_jesSubTotalMCDown, &b_MET_phi_jesSubTotalMCDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalDown", Jet_pt_jesTotalDown, &b_Jet_pt_jesTotalDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalDown", Jet_mass_jesTotalDown, &b_Jet_mass_jesTotalDown);
   fChain->SetBranchAddress("MET_pt_jesTotalDown", &MET_pt_jesTotalDown, &b_MET_pt_jesTotalDown);
   fChain->SetBranchAddress("MET_phi_jesTotalDown", &MET_phi_jesTotalDown, &b_MET_phi_jesTotalDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorDown", Jet_pt_jesTotalNoFlavorDown, &b_Jet_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorDown", Jet_mass_jesTotalNoFlavorDown, &b_Jet_mass_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorDown", &MET_pt_jesTotalNoFlavorDown, &b_MET_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorDown", &MET_phi_jesTotalNoFlavorDown, &b_MET_phi_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoTimeDown", Jet_pt_jesTotalNoTimeDown, &b_Jet_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoTimeDown", Jet_mass_jesTotalNoTimeDown, &b_Jet_mass_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_pt_jesTotalNoTimeDown", &MET_pt_jesTotalNoTimeDown, &b_MET_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_phi_jesTotalNoTimeDown", &MET_phi_jesTotalNoTimeDown, &b_MET_phi_jesTotalNoTimeDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeDown", Jet_pt_jesTotalNoFlavorNoTimeDown, &b_Jet_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeDown", Jet_mass_jesTotalNoFlavorNoTimeDown, &b_Jet_mass_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_pt_jesTotalNoFlavorNoTimeDown", &MET_pt_jesTotalNoFlavorNoTimeDown, &b_MET_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_phi_jesTotalNoFlavorNoTimeDown", &MET_phi_jesTotalNoFlavorNoTimeDown, &b_MET_phi_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorZJetDown", Jet_pt_jesFlavorZJetDown, &b_Jet_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorZJetDown", Jet_mass_jesFlavorZJetDown, &b_Jet_mass_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorZJetDown", &MET_pt_jesFlavorZJetDown, &b_MET_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorZJetDown", &MET_phi_jesFlavorZJetDown, &b_MET_phi_jesFlavorZJetDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPhotonJetDown", Jet_pt_jesFlavorPhotonJetDown, &b_Jet_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPhotonJetDown", Jet_mass_jesFlavorPhotonJetDown, &b_Jet_mass_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPhotonJetDown", &MET_pt_jesFlavorPhotonJetDown, &b_MET_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPhotonJetDown", &MET_phi_jesFlavorPhotonJetDown, &b_MET_phi_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureGluonDown", Jet_pt_jesFlavorPureGluonDown, &b_Jet_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureGluonDown", Jet_mass_jesFlavorPureGluonDown, &b_Jet_mass_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureGluonDown", &MET_pt_jesFlavorPureGluonDown, &b_MET_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureGluonDown", &MET_phi_jesFlavorPureGluonDown, &b_MET_phi_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureQuarkDown", Jet_pt_jesFlavorPureQuarkDown, &b_Jet_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureQuarkDown", Jet_mass_jesFlavorPureQuarkDown, &b_Jet_mass_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureQuarkDown", &MET_pt_jesFlavorPureQuarkDown, &b_MET_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureQuarkDown", &MET_phi_jesFlavorPureQuarkDown, &b_MET_phi_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureCharmDown", Jet_pt_jesFlavorPureCharmDown, &b_Jet_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureCharmDown", Jet_mass_jesFlavorPureCharmDown, &b_Jet_mass_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureCharmDown", &MET_pt_jesFlavorPureCharmDown, &b_MET_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureCharmDown", &MET_phi_jesFlavorPureCharmDown, &b_MET_phi_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureBottomDown", Jet_pt_jesFlavorPureBottomDown, &b_Jet_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureBottomDown", Jet_mass_jesFlavorPureBottomDown, &b_Jet_mass_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_pt_jesFlavorPureBottomDown", &MET_pt_jesFlavorPureBottomDown, &b_MET_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_phi_jesFlavorPureBottomDown", &MET_phi_jesFlavorPureBottomDown, &b_MET_phi_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunADown", Jet_pt_jesTimeRunADown, &b_Jet_pt_jesTimeRunADown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunADown", Jet_mass_jesTimeRunADown, &b_Jet_mass_jesTimeRunADown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunADown", &MET_pt_jesTimeRunADown, &b_MET_pt_jesTimeRunADown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunADown", &MET_phi_jesTimeRunADown, &b_MET_phi_jesTimeRunADown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunBDown", Jet_pt_jesTimeRunBDown, &b_Jet_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunBDown", Jet_mass_jesTimeRunBDown, &b_Jet_mass_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunBDown", &MET_pt_jesTimeRunBDown, &b_MET_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunBDown", &MET_phi_jesTimeRunBDown, &b_MET_phi_jesTimeRunBDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunCDown", Jet_pt_jesTimeRunCDown, &b_Jet_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunCDown", Jet_mass_jesTimeRunCDown, &b_Jet_mass_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunCDown", &MET_pt_jesTimeRunCDown, &b_MET_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunCDown", &MET_phi_jesTimeRunCDown, &b_MET_phi_jesTimeRunCDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunDDown", Jet_pt_jesTimeRunDDown, &b_Jet_pt_jesTimeRunDDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunDDown", Jet_mass_jesTimeRunDDown, &b_Jet_mass_jesTimeRunDDown);
   fChain->SetBranchAddress("MET_pt_jesTimeRunDDown", &MET_pt_jesTimeRunDDown, &b_MET_pt_jesTimeRunDDown);
   fChain->SetBranchAddress("MET_phi_jesTimeRunDDown", &MET_phi_jesTimeRunDDown, &b_MET_phi_jesTimeRunDDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituDown", Jet_pt_jesCorrelationGroupMPFInSituDown, &b_Jet_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituDown", Jet_mass_jesCorrelationGroupMPFInSituDown, &b_Jet_mass_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupMPFInSituDown", &MET_pt_jesCorrelationGroupMPFInSituDown, &b_MET_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupMPFInSituDown", &MET_phi_jesCorrelationGroupMPFInSituDown, &b_MET_phi_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationDown", Jet_pt_jesCorrelationGroupIntercalibrationDown, &b_Jet_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationDown", Jet_mass_jesCorrelationGroupIntercalibrationDown, &b_Jet_mass_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupIntercalibrationDown", &MET_pt_jesCorrelationGroupIntercalibrationDown, &b_MET_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupIntercalibrationDown", &MET_phi_jesCorrelationGroupIntercalibrationDown, &b_MET_phi_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESDown", Jet_pt_jesCorrelationGroupbJESDown, &b_Jet_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESDown", Jet_mass_jesCorrelationGroupbJESDown, &b_Jet_mass_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupbJESDown", &MET_pt_jesCorrelationGroupbJESDown, &b_MET_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupbJESDown", &MET_phi_jesCorrelationGroupbJESDown, &b_MET_phi_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorDown", Jet_pt_jesCorrelationGroupFlavorDown, &b_Jet_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorDown", Jet_mass_jesCorrelationGroupFlavorDown, &b_Jet_mass_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupFlavorDown", &MET_pt_jesCorrelationGroupFlavorDown, &b_MET_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupFlavorDown", &MET_phi_jesCorrelationGroupFlavorDown, &b_MET_phi_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedDown", Jet_pt_jesCorrelationGroupUncorrelatedDown, &b_Jet_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedDown", Jet_mass_jesCorrelationGroupUncorrelatedDown, &b_Jet_mass_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_pt_jesCorrelationGroupUncorrelatedDown", &MET_pt_jesCorrelationGroupUncorrelatedDown, &b_MET_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_phi_jesCorrelationGroupUncorrelatedDown", &MET_phi_jesCorrelationGroupUncorrelatedDown, &b_MET_phi_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_pt_unclustEnDown", &MET_pt_unclustEnDown, &b_MET_pt_unclustEnDown);
   fChain->SetBranchAddress("MET_phi_unclustEnDown", &MET_phi_unclustEnDown, &b_MET_phi_unclustEnDown);
   fChain->SetBranchAddress("FatJet_pt_raw", FatJet_pt_raw, &b_FatJet_pt_raw);
   fChain->SetBranchAddress("FatJet_pt_nom", FatJet_pt_nom, &b_FatJet_pt_nom);
   fChain->SetBranchAddress("FatJet_mass_raw", FatJet_mass_raw, &b_FatJet_mass_raw);
   fChain->SetBranchAddress("FatJet_mass_nom", FatJet_mass_nom, &b_FatJet_mass_nom);
   fChain->SetBranchAddress("FatJet_corr_JEC", FatJet_corr_JEC, &b_FatJet_corr_JEC);
   fChain->SetBranchAddress("FatJet_corr_JER", FatJet_corr_JER, &b_FatJet_corr_JER);
   fChain->SetBranchAddress("FatJet_corr_JMS", FatJet_corr_JMS, &b_FatJet_corr_JMS);
   fChain->SetBranchAddress("FatJet_corr_JMR", FatJet_corr_JMR, &b_FatJet_corr_JMR);
   fChain->SetBranchAddress("FatJet_msoftdrop_raw", FatJet_msoftdrop_raw, &b_FatJet_msoftdrop_raw);
   fChain->SetBranchAddress("FatJet_msoftdrop_nom", FatJet_msoftdrop_nom, &b_FatJet_msoftdrop_nom);
   fChain->SetBranchAddress("FatJet_msoftdrop_corr_JMR", FatJet_msoftdrop_corr_JMR, &b_FatJet_msoftdrop_corr_JMR);
   fChain->SetBranchAddress("FatJet_msoftdrop_corr_JMS", FatJet_msoftdrop_corr_JMS, &b_FatJet_msoftdrop_corr_JMS);
   fChain->SetBranchAddress("FatJet_msoftdrop_corr_PUPPI", FatJet_msoftdrop_corr_PUPPI, &b_FatJet_msoftdrop_corr_PUPPI);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_nom", FatJet_msoftdrop_tau21DDT_nom, &b_FatJet_msoftdrop_tau21DDT_nom);
   fChain->SetBranchAddress("FatJet_pt_jerUp", FatJet_pt_jerUp, &b_FatJet_pt_jerUp);
   fChain->SetBranchAddress("FatJet_mass_jerUp", FatJet_mass_jerUp, &b_FatJet_mass_jerUp);
   fChain->SetBranchAddress("FatJet_mass_jmrUp", FatJet_mass_jmrUp, &b_FatJet_mass_jmrUp);
   fChain->SetBranchAddress("FatJet_mass_jmsUp", FatJet_mass_jmsUp, &b_FatJet_mass_jmsUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jerUp", FatJet_msoftdrop_jerUp, &b_FatJet_msoftdrop_jerUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmrUp", FatJet_msoftdrop_jmrUp, &b_FatJet_msoftdrop_jmrUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmsUp", FatJet_msoftdrop_jmsUp, &b_FatJet_msoftdrop_jmsUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jerUp", FatJet_msoftdrop_tau21DDT_jerUp, &b_FatJet_msoftdrop_tau21DDT_jerUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmrUp", FatJet_msoftdrop_tau21DDT_jmrUp, &b_FatJet_msoftdrop_tau21DDT_jmrUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmsUp", FatJet_msoftdrop_tau21DDT_jmsUp, &b_FatJet_msoftdrop_tau21DDT_jmsUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteStatUp", FatJet_pt_jesAbsoluteStatUp, &b_FatJet_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteStatUp", FatJet_mass_jesAbsoluteStatUp, &b_FatJet_mass_jesAbsoluteStatUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteStatUp", FatJet_msoftdrop_jesAbsoluteStatUp, &b_FatJet_msoftdrop_jesAbsoluteStatUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteScaleUp", FatJet_pt_jesAbsoluteScaleUp, &b_FatJet_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteScaleUp", FatJet_mass_jesAbsoluteScaleUp, &b_FatJet_mass_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteScaleUp", FatJet_msoftdrop_jesAbsoluteScaleUp, &b_FatJet_msoftdrop_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteSampleUp", FatJet_pt_jesAbsoluteSampleUp, &b_FatJet_pt_jesAbsoluteSampleUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteSampleUp", FatJet_mass_jesAbsoluteSampleUp, &b_FatJet_mass_jesAbsoluteSampleUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteSampleUp", FatJet_msoftdrop_jesAbsoluteSampleUp, &b_FatJet_msoftdrop_jesAbsoluteSampleUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteFlavMapUp", FatJet_pt_jesAbsoluteFlavMapUp, &b_FatJet_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteFlavMapUp", FatJet_mass_jesAbsoluteFlavMapUp, &b_FatJet_mass_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteFlavMapUp", FatJet_msoftdrop_jesAbsoluteFlavMapUp, &b_FatJet_msoftdrop_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteMPFBiasUp", FatJet_pt_jesAbsoluteMPFBiasUp, &b_FatJet_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteMPFBiasUp", FatJet_mass_jesAbsoluteMPFBiasUp, &b_FatJet_mass_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteMPFBiasUp", FatJet_msoftdrop_jesAbsoluteMPFBiasUp, &b_FatJet_msoftdrop_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("FatJet_pt_jesFragmentationUp", FatJet_pt_jesFragmentationUp, &b_FatJet_pt_jesFragmentationUp);
   fChain->SetBranchAddress("FatJet_mass_jesFragmentationUp", FatJet_mass_jesFragmentationUp, &b_FatJet_mass_jesFragmentationUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFragmentationUp", FatJet_msoftdrop_jesFragmentationUp, &b_FatJet_msoftdrop_jesFragmentationUp);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionECALUp", FatJet_pt_jesSinglePionECALUp, &b_FatJet_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionECALUp", FatJet_mass_jesSinglePionECALUp, &b_FatJet_mass_jesSinglePionECALUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionECALUp", FatJet_msoftdrop_jesSinglePionECALUp, &b_FatJet_msoftdrop_jesSinglePionECALUp);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionHCALUp", FatJet_pt_jesSinglePionHCALUp, &b_FatJet_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionHCALUp", FatJet_mass_jesSinglePionHCALUp, &b_FatJet_mass_jesSinglePionHCALUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionHCALUp", FatJet_msoftdrop_jesSinglePionHCALUp, &b_FatJet_msoftdrop_jesSinglePionHCALUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorQCDUp", FatJet_pt_jesFlavorQCDUp, &b_FatJet_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorQCDUp", FatJet_mass_jesFlavorQCDUp, &b_FatJet_mass_jesFlavorQCDUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorQCDUp", FatJet_msoftdrop_jesFlavorQCDUp, &b_FatJet_msoftdrop_jesFlavorQCDUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimePtEtaUp", FatJet_pt_jesTimePtEtaUp, &b_FatJet_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimePtEtaUp", FatJet_mass_jesTimePtEtaUp, &b_FatJet_mass_jesTimePtEtaUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimePtEtaUp", FatJet_msoftdrop_jesTimePtEtaUp, &b_FatJet_msoftdrop_jesTimePtEtaUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC1Up", FatJet_pt_jesRelativeJEREC1Up, &b_FatJet_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC1Up", FatJet_mass_jesRelativeJEREC1Up, &b_FatJet_mass_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC1Up", FatJet_msoftdrop_jesRelativeJEREC1Up, &b_FatJet_msoftdrop_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC2Up", FatJet_pt_jesRelativeJEREC2Up, &b_FatJet_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC2Up", FatJet_mass_jesRelativeJEREC2Up, &b_FatJet_mass_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC2Up", FatJet_msoftdrop_jesRelativeJEREC2Up, &b_FatJet_msoftdrop_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJERHFUp", FatJet_pt_jesRelativeJERHFUp, &b_FatJet_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJERHFUp", FatJet_mass_jesRelativeJERHFUp, &b_FatJet_mass_jesRelativeJERHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJERHFUp", FatJet_msoftdrop_jesRelativeJERHFUp, &b_FatJet_msoftdrop_jesRelativeJERHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtBBUp", FatJet_pt_jesRelativePtBBUp, &b_FatJet_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtBBUp", FatJet_mass_jesRelativePtBBUp, &b_FatJet_mass_jesRelativePtBBUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtBBUp", FatJet_msoftdrop_jesRelativePtBBUp, &b_FatJet_msoftdrop_jesRelativePtBBUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC1Up", FatJet_pt_jesRelativePtEC1Up, &b_FatJet_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC1Up", FatJet_mass_jesRelativePtEC1Up, &b_FatJet_mass_jesRelativePtEC1Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC1Up", FatJet_msoftdrop_jesRelativePtEC1Up, &b_FatJet_msoftdrop_jesRelativePtEC1Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC2Up", FatJet_pt_jesRelativePtEC2Up, &b_FatJet_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC2Up", FatJet_mass_jesRelativePtEC2Up, &b_FatJet_mass_jesRelativePtEC2Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC2Up", FatJet_msoftdrop_jesRelativePtEC2Up, &b_FatJet_msoftdrop_jesRelativePtEC2Up);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtHFUp", FatJet_pt_jesRelativePtHFUp, &b_FatJet_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtHFUp", FatJet_mass_jesRelativePtHFUp, &b_FatJet_mass_jesRelativePtHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtHFUp", FatJet_msoftdrop_jesRelativePtHFUp, &b_FatJet_msoftdrop_jesRelativePtHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeBalUp", FatJet_pt_jesRelativeBalUp, &b_FatJet_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeBalUp", FatJet_mass_jesRelativeBalUp, &b_FatJet_mass_jesRelativeBalUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeBalUp", FatJet_msoftdrop_jesRelativeBalUp, &b_FatJet_msoftdrop_jesRelativeBalUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeSampleUp", FatJet_pt_jesRelativeSampleUp, &b_FatJet_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeSampleUp", FatJet_mass_jesRelativeSampleUp, &b_FatJet_mass_jesRelativeSampleUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeSampleUp", FatJet_msoftdrop_jesRelativeSampleUp, &b_FatJet_msoftdrop_jesRelativeSampleUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeFSRUp", FatJet_pt_jesRelativeFSRUp, &b_FatJet_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeFSRUp", FatJet_mass_jesRelativeFSRUp, &b_FatJet_mass_jesRelativeFSRUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeFSRUp", FatJet_msoftdrop_jesRelativeFSRUp, &b_FatJet_msoftdrop_jesRelativeFSRUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatFSRUp", FatJet_pt_jesRelativeStatFSRUp, &b_FatJet_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatFSRUp", FatJet_mass_jesRelativeStatFSRUp, &b_FatJet_mass_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatFSRUp", FatJet_msoftdrop_jesRelativeStatFSRUp, &b_FatJet_msoftdrop_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatECUp", FatJet_pt_jesRelativeStatECUp, &b_FatJet_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatECUp", FatJet_mass_jesRelativeStatECUp, &b_FatJet_mass_jesRelativeStatECUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatECUp", FatJet_msoftdrop_jesRelativeStatECUp, &b_FatJet_msoftdrop_jesRelativeStatECUp);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatHFUp", FatJet_pt_jesRelativeStatHFUp, &b_FatJet_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatHFUp", FatJet_mass_jesRelativeStatHFUp, &b_FatJet_mass_jesRelativeStatHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatHFUp", FatJet_msoftdrop_jesRelativeStatHFUp, &b_FatJet_msoftdrop_jesRelativeStatHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpDataMCUp", FatJet_pt_jesPileUpDataMCUp, &b_FatJet_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpDataMCUp", FatJet_mass_jesPileUpDataMCUp, &b_FatJet_mass_jesPileUpDataMCUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpDataMCUp", FatJet_msoftdrop_jesPileUpDataMCUp, &b_FatJet_msoftdrop_jesPileUpDataMCUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtRefUp", FatJet_pt_jesPileUpPtRefUp, &b_FatJet_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtRefUp", FatJet_mass_jesPileUpPtRefUp, &b_FatJet_mass_jesPileUpPtRefUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtRefUp", FatJet_msoftdrop_jesPileUpPtRefUp, &b_FatJet_msoftdrop_jesPileUpPtRefUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtBBUp", FatJet_pt_jesPileUpPtBBUp, &b_FatJet_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtBBUp", FatJet_mass_jesPileUpPtBBUp, &b_FatJet_mass_jesPileUpPtBBUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtBBUp", FatJet_msoftdrop_jesPileUpPtBBUp, &b_FatJet_msoftdrop_jesPileUpPtBBUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC1Up", FatJet_pt_jesPileUpPtEC1Up, &b_FatJet_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC1Up", FatJet_mass_jesPileUpPtEC1Up, &b_FatJet_mass_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC1Up", FatJet_msoftdrop_jesPileUpPtEC1Up, &b_FatJet_msoftdrop_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC2Up", FatJet_pt_jesPileUpPtEC2Up, &b_FatJet_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC2Up", FatJet_mass_jesPileUpPtEC2Up, &b_FatJet_mass_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC2Up", FatJet_msoftdrop_jesPileUpPtEC2Up, &b_FatJet_msoftdrop_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtHFUp", FatJet_pt_jesPileUpPtHFUp, &b_FatJet_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtHFUp", FatJet_mass_jesPileUpPtHFUp, &b_FatJet_mass_jesPileUpPtHFUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtHFUp", FatJet_msoftdrop_jesPileUpPtHFUp, &b_FatJet_msoftdrop_jesPileUpPtHFUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpMuZeroUp", FatJet_pt_jesPileUpMuZeroUp, &b_FatJet_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpMuZeroUp", FatJet_mass_jesPileUpMuZeroUp, &b_FatJet_mass_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpMuZeroUp", FatJet_msoftdrop_jesPileUpMuZeroUp, &b_FatJet_msoftdrop_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpEnvelopeUp", FatJet_pt_jesPileUpEnvelopeUp, &b_FatJet_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpEnvelopeUp", FatJet_mass_jesPileUpEnvelopeUp, &b_FatJet_mass_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpEnvelopeUp", FatJet_msoftdrop_jesPileUpEnvelopeUp, &b_FatJet_msoftdrop_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPileUpUp", FatJet_pt_jesSubTotalPileUpUp, &b_FatJet_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPileUpUp", FatJet_mass_jesSubTotalPileUpUp, &b_FatJet_mass_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPileUpUp", FatJet_msoftdrop_jesSubTotalPileUpUp, &b_FatJet_msoftdrop_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalRelativeUp", FatJet_pt_jesSubTotalRelativeUp, &b_FatJet_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalRelativeUp", FatJet_mass_jesSubTotalRelativeUp, &b_FatJet_mass_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalRelativeUp", FatJet_msoftdrop_jesSubTotalRelativeUp, &b_FatJet_msoftdrop_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPtUp", FatJet_pt_jesSubTotalPtUp, &b_FatJet_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPtUp", FatJet_mass_jesSubTotalPtUp, &b_FatJet_mass_jesSubTotalPtUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPtUp", FatJet_msoftdrop_jesSubTotalPtUp, &b_FatJet_msoftdrop_jesSubTotalPtUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalScaleUp", FatJet_pt_jesSubTotalScaleUp, &b_FatJet_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalScaleUp", FatJet_mass_jesSubTotalScaleUp, &b_FatJet_mass_jesSubTotalScaleUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalScaleUp", FatJet_msoftdrop_jesSubTotalScaleUp, &b_FatJet_msoftdrop_jesSubTotalScaleUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalAbsoluteUp", FatJet_pt_jesSubTotalAbsoluteUp, &b_FatJet_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalAbsoluteUp", FatJet_mass_jesSubTotalAbsoluteUp, &b_FatJet_mass_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalAbsoluteUp", FatJet_msoftdrop_jesSubTotalAbsoluteUp, &b_FatJet_msoftdrop_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalMCUp", FatJet_pt_jesSubTotalMCUp, &b_FatJet_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalMCUp", FatJet_mass_jesSubTotalMCUp, &b_FatJet_mass_jesSubTotalMCUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalMCUp", FatJet_msoftdrop_jesSubTotalMCUp, &b_FatJet_msoftdrop_jesSubTotalMCUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalUp", FatJet_pt_jesTotalUp, &b_FatJet_pt_jesTotalUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalUp", FatJet_mass_jesTotalUp, &b_FatJet_mass_jesTotalUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalUp", FatJet_msoftdrop_jesTotalUp, &b_FatJet_msoftdrop_jesTotalUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorUp", FatJet_pt_jesTotalNoFlavorUp, &b_FatJet_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorUp", FatJet_mass_jesTotalNoFlavorUp, &b_FatJet_mass_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorUp", FatJet_msoftdrop_jesTotalNoFlavorUp, &b_FatJet_msoftdrop_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoTimeUp", FatJet_pt_jesTotalNoTimeUp, &b_FatJet_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoTimeUp", FatJet_mass_jesTotalNoTimeUp, &b_FatJet_mass_jesTotalNoTimeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoTimeUp", FatJet_msoftdrop_jesTotalNoTimeUp, &b_FatJet_msoftdrop_jesTotalNoTimeUp);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorNoTimeUp", FatJet_pt_jesTotalNoFlavorNoTimeUp, &b_FatJet_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorNoTimeUp", FatJet_mass_jesTotalNoFlavorNoTimeUp, &b_FatJet_mass_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp", FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp, &b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorZJetUp", FatJet_pt_jesFlavorZJetUp, &b_FatJet_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorZJetUp", FatJet_mass_jesFlavorZJetUp, &b_FatJet_mass_jesFlavorZJetUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorZJetUp", FatJet_msoftdrop_jesFlavorZJetUp, &b_FatJet_msoftdrop_jesFlavorZJetUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPhotonJetUp", FatJet_pt_jesFlavorPhotonJetUp, &b_FatJet_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPhotonJetUp", FatJet_mass_jesFlavorPhotonJetUp, &b_FatJet_mass_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPhotonJetUp", FatJet_msoftdrop_jesFlavorPhotonJetUp, &b_FatJet_msoftdrop_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureGluonUp", FatJet_pt_jesFlavorPureGluonUp, &b_FatJet_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureGluonUp", FatJet_mass_jesFlavorPureGluonUp, &b_FatJet_mass_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureGluonUp", FatJet_msoftdrop_jesFlavorPureGluonUp, &b_FatJet_msoftdrop_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureQuarkUp", FatJet_pt_jesFlavorPureQuarkUp, &b_FatJet_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureQuarkUp", FatJet_mass_jesFlavorPureQuarkUp, &b_FatJet_mass_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureQuarkUp", FatJet_msoftdrop_jesFlavorPureQuarkUp, &b_FatJet_msoftdrop_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureCharmUp", FatJet_pt_jesFlavorPureCharmUp, &b_FatJet_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureCharmUp", FatJet_mass_jesFlavorPureCharmUp, &b_FatJet_mass_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureCharmUp", FatJet_msoftdrop_jesFlavorPureCharmUp, &b_FatJet_msoftdrop_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureBottomUp", FatJet_pt_jesFlavorPureBottomUp, &b_FatJet_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureBottomUp", FatJet_mass_jesFlavorPureBottomUp, &b_FatJet_mass_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureBottomUp", FatJet_msoftdrop_jesFlavorPureBottomUp, &b_FatJet_msoftdrop_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunAUp", FatJet_pt_jesTimeRunAUp, &b_FatJet_pt_jesTimeRunAUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunAUp", FatJet_mass_jesTimeRunAUp, &b_FatJet_mass_jesTimeRunAUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunAUp", FatJet_msoftdrop_jesTimeRunAUp, &b_FatJet_msoftdrop_jesTimeRunAUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunBUp", FatJet_pt_jesTimeRunBUp, &b_FatJet_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunBUp", FatJet_mass_jesTimeRunBUp, &b_FatJet_mass_jesTimeRunBUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunBUp", FatJet_msoftdrop_jesTimeRunBUp, &b_FatJet_msoftdrop_jesTimeRunBUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunCUp", FatJet_pt_jesTimeRunCUp, &b_FatJet_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunCUp", FatJet_mass_jesTimeRunCUp, &b_FatJet_mass_jesTimeRunCUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunCUp", FatJet_msoftdrop_jesTimeRunCUp, &b_FatJet_msoftdrop_jesTimeRunCUp);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunDUp", FatJet_pt_jesTimeRunDUp, &b_FatJet_pt_jesTimeRunDUp);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunDUp", FatJet_mass_jesTimeRunDUp, &b_FatJet_mass_jesTimeRunDUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunDUp", FatJet_msoftdrop_jesTimeRunDUp, &b_FatJet_msoftdrop_jesTimeRunDUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupMPFInSituUp", FatJet_pt_jesCorrelationGroupMPFInSituUp, &b_FatJet_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupMPFInSituUp", FatJet_mass_jesCorrelationGroupMPFInSituUp, &b_FatJet_mass_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp", FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp, &b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupIntercalibrationUp", FatJet_pt_jesCorrelationGroupIntercalibrationUp, &b_FatJet_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupIntercalibrationUp", FatJet_mass_jesCorrelationGroupIntercalibrationUp, &b_FatJet_mass_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp", FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp, &b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupbJESUp", FatJet_pt_jesCorrelationGroupbJESUp, &b_FatJet_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupbJESUp", FatJet_mass_jesCorrelationGroupbJESUp, &b_FatJet_mass_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupbJESUp", FatJet_msoftdrop_jesCorrelationGroupbJESUp, &b_FatJet_msoftdrop_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupFlavorUp", FatJet_pt_jesCorrelationGroupFlavorUp, &b_FatJet_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupFlavorUp", FatJet_mass_jesCorrelationGroupFlavorUp, &b_FatJet_mass_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupFlavorUp", FatJet_msoftdrop_jesCorrelationGroupFlavorUp, &b_FatJet_msoftdrop_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupUncorrelatedUp", FatJet_pt_jesCorrelationGroupUncorrelatedUp, &b_FatJet_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupUncorrelatedUp", FatJet_mass_jesCorrelationGroupUncorrelatedUp, &b_FatJet_mass_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp", FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp, &b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("FatJet_pt_jerDown", FatJet_pt_jerDown, &b_FatJet_pt_jerDown);
   fChain->SetBranchAddress("FatJet_mass_jerDown", FatJet_mass_jerDown, &b_FatJet_mass_jerDown);
   fChain->SetBranchAddress("FatJet_mass_jmrDown", FatJet_mass_jmrDown, &b_FatJet_mass_jmrDown);
   fChain->SetBranchAddress("FatJet_mass_jmsDown", FatJet_mass_jmsDown, &b_FatJet_mass_jmsDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jerDown", FatJet_msoftdrop_jerDown, &b_FatJet_msoftdrop_jerDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmrDown", FatJet_msoftdrop_jmrDown, &b_FatJet_msoftdrop_jmrDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jmsDown", FatJet_msoftdrop_jmsDown, &b_FatJet_msoftdrop_jmsDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jerDown", FatJet_msoftdrop_tau21DDT_jerDown, &b_FatJet_msoftdrop_tau21DDT_jerDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmrDown", FatJet_msoftdrop_tau21DDT_jmrDown, &b_FatJet_msoftdrop_tau21DDT_jmrDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_tau21DDT_jmsDown", FatJet_msoftdrop_tau21DDT_jmsDown, &b_FatJet_msoftdrop_tau21DDT_jmsDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteStatDown", FatJet_pt_jesAbsoluteStatDown, &b_FatJet_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteStatDown", FatJet_mass_jesAbsoluteStatDown, &b_FatJet_mass_jesAbsoluteStatDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteStatDown", FatJet_msoftdrop_jesAbsoluteStatDown, &b_FatJet_msoftdrop_jesAbsoluteStatDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteScaleDown", FatJet_pt_jesAbsoluteScaleDown, &b_FatJet_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteScaleDown", FatJet_mass_jesAbsoluteScaleDown, &b_FatJet_mass_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteScaleDown", FatJet_msoftdrop_jesAbsoluteScaleDown, &b_FatJet_msoftdrop_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteSampleDown", FatJet_pt_jesAbsoluteSampleDown, &b_FatJet_pt_jesAbsoluteSampleDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteSampleDown", FatJet_mass_jesAbsoluteSampleDown, &b_FatJet_mass_jesAbsoluteSampleDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteSampleDown", FatJet_msoftdrop_jesAbsoluteSampleDown, &b_FatJet_msoftdrop_jesAbsoluteSampleDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteFlavMapDown", FatJet_pt_jesAbsoluteFlavMapDown, &b_FatJet_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteFlavMapDown", FatJet_mass_jesAbsoluteFlavMapDown, &b_FatJet_mass_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteFlavMapDown", FatJet_msoftdrop_jesAbsoluteFlavMapDown, &b_FatJet_msoftdrop_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("FatJet_pt_jesAbsoluteMPFBiasDown", FatJet_pt_jesAbsoluteMPFBiasDown, &b_FatJet_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("FatJet_mass_jesAbsoluteMPFBiasDown", FatJet_mass_jesAbsoluteMPFBiasDown, &b_FatJet_mass_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesAbsoluteMPFBiasDown", FatJet_msoftdrop_jesAbsoluteMPFBiasDown, &b_FatJet_msoftdrop_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("FatJet_pt_jesFragmentationDown", FatJet_pt_jesFragmentationDown, &b_FatJet_pt_jesFragmentationDown);
   fChain->SetBranchAddress("FatJet_mass_jesFragmentationDown", FatJet_mass_jesFragmentationDown, &b_FatJet_mass_jesFragmentationDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFragmentationDown", FatJet_msoftdrop_jesFragmentationDown, &b_FatJet_msoftdrop_jesFragmentationDown);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionECALDown", FatJet_pt_jesSinglePionECALDown, &b_FatJet_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionECALDown", FatJet_mass_jesSinglePionECALDown, &b_FatJet_mass_jesSinglePionECALDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionECALDown", FatJet_msoftdrop_jesSinglePionECALDown, &b_FatJet_msoftdrop_jesSinglePionECALDown);
   fChain->SetBranchAddress("FatJet_pt_jesSinglePionHCALDown", FatJet_pt_jesSinglePionHCALDown, &b_FatJet_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("FatJet_mass_jesSinglePionHCALDown", FatJet_mass_jesSinglePionHCALDown, &b_FatJet_mass_jesSinglePionHCALDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSinglePionHCALDown", FatJet_msoftdrop_jesSinglePionHCALDown, &b_FatJet_msoftdrop_jesSinglePionHCALDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorQCDDown", FatJet_pt_jesFlavorQCDDown, &b_FatJet_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorQCDDown", FatJet_mass_jesFlavorQCDDown, &b_FatJet_mass_jesFlavorQCDDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorQCDDown", FatJet_msoftdrop_jesFlavorQCDDown, &b_FatJet_msoftdrop_jesFlavorQCDDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimePtEtaDown", FatJet_pt_jesTimePtEtaDown, &b_FatJet_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimePtEtaDown", FatJet_mass_jesTimePtEtaDown, &b_FatJet_mass_jesTimePtEtaDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimePtEtaDown", FatJet_msoftdrop_jesTimePtEtaDown, &b_FatJet_msoftdrop_jesTimePtEtaDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC1Down", FatJet_pt_jesRelativeJEREC1Down, &b_FatJet_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC1Down", FatJet_mass_jesRelativeJEREC1Down, &b_FatJet_mass_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC1Down", FatJet_msoftdrop_jesRelativeJEREC1Down, &b_FatJet_msoftdrop_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJEREC2Down", FatJet_pt_jesRelativeJEREC2Down, &b_FatJet_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJEREC2Down", FatJet_mass_jesRelativeJEREC2Down, &b_FatJet_mass_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJEREC2Down", FatJet_msoftdrop_jesRelativeJEREC2Down, &b_FatJet_msoftdrop_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeJERHFDown", FatJet_pt_jesRelativeJERHFDown, &b_FatJet_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeJERHFDown", FatJet_mass_jesRelativeJERHFDown, &b_FatJet_mass_jesRelativeJERHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeJERHFDown", FatJet_msoftdrop_jesRelativeJERHFDown, &b_FatJet_msoftdrop_jesRelativeJERHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtBBDown", FatJet_pt_jesRelativePtBBDown, &b_FatJet_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtBBDown", FatJet_mass_jesRelativePtBBDown, &b_FatJet_mass_jesRelativePtBBDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtBBDown", FatJet_msoftdrop_jesRelativePtBBDown, &b_FatJet_msoftdrop_jesRelativePtBBDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC1Down", FatJet_pt_jesRelativePtEC1Down, &b_FatJet_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC1Down", FatJet_mass_jesRelativePtEC1Down, &b_FatJet_mass_jesRelativePtEC1Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC1Down", FatJet_msoftdrop_jesRelativePtEC1Down, &b_FatJet_msoftdrop_jesRelativePtEC1Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtEC2Down", FatJet_pt_jesRelativePtEC2Down, &b_FatJet_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtEC2Down", FatJet_mass_jesRelativePtEC2Down, &b_FatJet_mass_jesRelativePtEC2Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtEC2Down", FatJet_msoftdrop_jesRelativePtEC2Down, &b_FatJet_msoftdrop_jesRelativePtEC2Down);
   fChain->SetBranchAddress("FatJet_pt_jesRelativePtHFDown", FatJet_pt_jesRelativePtHFDown, &b_FatJet_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativePtHFDown", FatJet_mass_jesRelativePtHFDown, &b_FatJet_mass_jesRelativePtHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativePtHFDown", FatJet_msoftdrop_jesRelativePtHFDown, &b_FatJet_msoftdrop_jesRelativePtHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeBalDown", FatJet_pt_jesRelativeBalDown, &b_FatJet_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeBalDown", FatJet_mass_jesRelativeBalDown, &b_FatJet_mass_jesRelativeBalDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeBalDown", FatJet_msoftdrop_jesRelativeBalDown, &b_FatJet_msoftdrop_jesRelativeBalDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeSampleDown", FatJet_pt_jesRelativeSampleDown, &b_FatJet_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeSampleDown", FatJet_mass_jesRelativeSampleDown, &b_FatJet_mass_jesRelativeSampleDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeSampleDown", FatJet_msoftdrop_jesRelativeSampleDown, &b_FatJet_msoftdrop_jesRelativeSampleDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeFSRDown", FatJet_pt_jesRelativeFSRDown, &b_FatJet_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeFSRDown", FatJet_mass_jesRelativeFSRDown, &b_FatJet_mass_jesRelativeFSRDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeFSRDown", FatJet_msoftdrop_jesRelativeFSRDown, &b_FatJet_msoftdrop_jesRelativeFSRDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatFSRDown", FatJet_pt_jesRelativeStatFSRDown, &b_FatJet_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatFSRDown", FatJet_mass_jesRelativeStatFSRDown, &b_FatJet_mass_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatFSRDown", FatJet_msoftdrop_jesRelativeStatFSRDown, &b_FatJet_msoftdrop_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatECDown", FatJet_pt_jesRelativeStatECDown, &b_FatJet_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatECDown", FatJet_mass_jesRelativeStatECDown, &b_FatJet_mass_jesRelativeStatECDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatECDown", FatJet_msoftdrop_jesRelativeStatECDown, &b_FatJet_msoftdrop_jesRelativeStatECDown);
   fChain->SetBranchAddress("FatJet_pt_jesRelativeStatHFDown", FatJet_pt_jesRelativeStatHFDown, &b_FatJet_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesRelativeStatHFDown", FatJet_mass_jesRelativeStatHFDown, &b_FatJet_mass_jesRelativeStatHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesRelativeStatHFDown", FatJet_msoftdrop_jesRelativeStatHFDown, &b_FatJet_msoftdrop_jesRelativeStatHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpDataMCDown", FatJet_pt_jesPileUpDataMCDown, &b_FatJet_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpDataMCDown", FatJet_mass_jesPileUpDataMCDown, &b_FatJet_mass_jesPileUpDataMCDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpDataMCDown", FatJet_msoftdrop_jesPileUpDataMCDown, &b_FatJet_msoftdrop_jesPileUpDataMCDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtRefDown", FatJet_pt_jesPileUpPtRefDown, &b_FatJet_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtRefDown", FatJet_mass_jesPileUpPtRefDown, &b_FatJet_mass_jesPileUpPtRefDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtRefDown", FatJet_msoftdrop_jesPileUpPtRefDown, &b_FatJet_msoftdrop_jesPileUpPtRefDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtBBDown", FatJet_pt_jesPileUpPtBBDown, &b_FatJet_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtBBDown", FatJet_mass_jesPileUpPtBBDown, &b_FatJet_mass_jesPileUpPtBBDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtBBDown", FatJet_msoftdrop_jesPileUpPtBBDown, &b_FatJet_msoftdrop_jesPileUpPtBBDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC1Down", FatJet_pt_jesPileUpPtEC1Down, &b_FatJet_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC1Down", FatJet_mass_jesPileUpPtEC1Down, &b_FatJet_mass_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC1Down", FatJet_msoftdrop_jesPileUpPtEC1Down, &b_FatJet_msoftdrop_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtEC2Down", FatJet_pt_jesPileUpPtEC2Down, &b_FatJet_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtEC2Down", FatJet_mass_jesPileUpPtEC2Down, &b_FatJet_mass_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtEC2Down", FatJet_msoftdrop_jesPileUpPtEC2Down, &b_FatJet_msoftdrop_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpPtHFDown", FatJet_pt_jesPileUpPtHFDown, &b_FatJet_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpPtHFDown", FatJet_mass_jesPileUpPtHFDown, &b_FatJet_mass_jesPileUpPtHFDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpPtHFDown", FatJet_msoftdrop_jesPileUpPtHFDown, &b_FatJet_msoftdrop_jesPileUpPtHFDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpMuZeroDown", FatJet_pt_jesPileUpMuZeroDown, &b_FatJet_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpMuZeroDown", FatJet_mass_jesPileUpMuZeroDown, &b_FatJet_mass_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpMuZeroDown", FatJet_msoftdrop_jesPileUpMuZeroDown, &b_FatJet_msoftdrop_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("FatJet_pt_jesPileUpEnvelopeDown", FatJet_pt_jesPileUpEnvelopeDown, &b_FatJet_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("FatJet_mass_jesPileUpEnvelopeDown", FatJet_mass_jesPileUpEnvelopeDown, &b_FatJet_mass_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesPileUpEnvelopeDown", FatJet_msoftdrop_jesPileUpEnvelopeDown, &b_FatJet_msoftdrop_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPileUpDown", FatJet_pt_jesSubTotalPileUpDown, &b_FatJet_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPileUpDown", FatJet_mass_jesSubTotalPileUpDown, &b_FatJet_mass_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPileUpDown", FatJet_msoftdrop_jesSubTotalPileUpDown, &b_FatJet_msoftdrop_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalRelativeDown", FatJet_pt_jesSubTotalRelativeDown, &b_FatJet_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalRelativeDown", FatJet_mass_jesSubTotalRelativeDown, &b_FatJet_mass_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalRelativeDown", FatJet_msoftdrop_jesSubTotalRelativeDown, &b_FatJet_msoftdrop_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalPtDown", FatJet_pt_jesSubTotalPtDown, &b_FatJet_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalPtDown", FatJet_mass_jesSubTotalPtDown, &b_FatJet_mass_jesSubTotalPtDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalPtDown", FatJet_msoftdrop_jesSubTotalPtDown, &b_FatJet_msoftdrop_jesSubTotalPtDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalScaleDown", FatJet_pt_jesSubTotalScaleDown, &b_FatJet_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalScaleDown", FatJet_mass_jesSubTotalScaleDown, &b_FatJet_mass_jesSubTotalScaleDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalScaleDown", FatJet_msoftdrop_jesSubTotalScaleDown, &b_FatJet_msoftdrop_jesSubTotalScaleDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalAbsoluteDown", FatJet_pt_jesSubTotalAbsoluteDown, &b_FatJet_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalAbsoluteDown", FatJet_mass_jesSubTotalAbsoluteDown, &b_FatJet_mass_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalAbsoluteDown", FatJet_msoftdrop_jesSubTotalAbsoluteDown, &b_FatJet_msoftdrop_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("FatJet_pt_jesSubTotalMCDown", FatJet_pt_jesSubTotalMCDown, &b_FatJet_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("FatJet_mass_jesSubTotalMCDown", FatJet_mass_jesSubTotalMCDown, &b_FatJet_mass_jesSubTotalMCDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesSubTotalMCDown", FatJet_msoftdrop_jesSubTotalMCDown, &b_FatJet_msoftdrop_jesSubTotalMCDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalDown", FatJet_pt_jesTotalDown, &b_FatJet_pt_jesTotalDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalDown", FatJet_mass_jesTotalDown, &b_FatJet_mass_jesTotalDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalDown", FatJet_msoftdrop_jesTotalDown, &b_FatJet_msoftdrop_jesTotalDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorDown", FatJet_pt_jesTotalNoFlavorDown, &b_FatJet_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorDown", FatJet_mass_jesTotalNoFlavorDown, &b_FatJet_mass_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorDown", FatJet_msoftdrop_jesTotalNoFlavorDown, &b_FatJet_msoftdrop_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoTimeDown", FatJet_pt_jesTotalNoTimeDown, &b_FatJet_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoTimeDown", FatJet_mass_jesTotalNoTimeDown, &b_FatJet_mass_jesTotalNoTimeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoTimeDown", FatJet_msoftdrop_jesTotalNoTimeDown, &b_FatJet_msoftdrop_jesTotalNoTimeDown);
   fChain->SetBranchAddress("FatJet_pt_jesTotalNoFlavorNoTimeDown", FatJet_pt_jesTotalNoFlavorNoTimeDown, &b_FatJet_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("FatJet_mass_jesTotalNoFlavorNoTimeDown", FatJet_mass_jesTotalNoFlavorNoTimeDown, &b_FatJet_mass_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown", FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown, &b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorZJetDown", FatJet_pt_jesFlavorZJetDown, &b_FatJet_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorZJetDown", FatJet_mass_jesFlavorZJetDown, &b_FatJet_mass_jesFlavorZJetDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorZJetDown", FatJet_msoftdrop_jesFlavorZJetDown, &b_FatJet_msoftdrop_jesFlavorZJetDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPhotonJetDown", FatJet_pt_jesFlavorPhotonJetDown, &b_FatJet_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPhotonJetDown", FatJet_mass_jesFlavorPhotonJetDown, &b_FatJet_mass_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPhotonJetDown", FatJet_msoftdrop_jesFlavorPhotonJetDown, &b_FatJet_msoftdrop_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureGluonDown", FatJet_pt_jesFlavorPureGluonDown, &b_FatJet_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureGluonDown", FatJet_mass_jesFlavorPureGluonDown, &b_FatJet_mass_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureGluonDown", FatJet_msoftdrop_jesFlavorPureGluonDown, &b_FatJet_msoftdrop_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureQuarkDown", FatJet_pt_jesFlavorPureQuarkDown, &b_FatJet_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureQuarkDown", FatJet_mass_jesFlavorPureQuarkDown, &b_FatJet_mass_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureQuarkDown", FatJet_msoftdrop_jesFlavorPureQuarkDown, &b_FatJet_msoftdrop_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureCharmDown", FatJet_pt_jesFlavorPureCharmDown, &b_FatJet_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureCharmDown", FatJet_mass_jesFlavorPureCharmDown, &b_FatJet_mass_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureCharmDown", FatJet_msoftdrop_jesFlavorPureCharmDown, &b_FatJet_msoftdrop_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("FatJet_pt_jesFlavorPureBottomDown", FatJet_pt_jesFlavorPureBottomDown, &b_FatJet_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("FatJet_mass_jesFlavorPureBottomDown", FatJet_mass_jesFlavorPureBottomDown, &b_FatJet_mass_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesFlavorPureBottomDown", FatJet_msoftdrop_jesFlavorPureBottomDown, &b_FatJet_msoftdrop_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunADown", FatJet_pt_jesTimeRunADown, &b_FatJet_pt_jesTimeRunADown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunADown", FatJet_mass_jesTimeRunADown, &b_FatJet_mass_jesTimeRunADown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunADown", FatJet_msoftdrop_jesTimeRunADown, &b_FatJet_msoftdrop_jesTimeRunADown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunBDown", FatJet_pt_jesTimeRunBDown, &b_FatJet_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunBDown", FatJet_mass_jesTimeRunBDown, &b_FatJet_mass_jesTimeRunBDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunBDown", FatJet_msoftdrop_jesTimeRunBDown, &b_FatJet_msoftdrop_jesTimeRunBDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunCDown", FatJet_pt_jesTimeRunCDown, &b_FatJet_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunCDown", FatJet_mass_jesTimeRunCDown, &b_FatJet_mass_jesTimeRunCDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunCDown", FatJet_msoftdrop_jesTimeRunCDown, &b_FatJet_msoftdrop_jesTimeRunCDown);
   fChain->SetBranchAddress("FatJet_pt_jesTimeRunDDown", FatJet_pt_jesTimeRunDDown, &b_FatJet_pt_jesTimeRunDDown);
   fChain->SetBranchAddress("FatJet_mass_jesTimeRunDDown", FatJet_mass_jesTimeRunDDown, &b_FatJet_mass_jesTimeRunDDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesTimeRunDDown", FatJet_msoftdrop_jesTimeRunDDown, &b_FatJet_msoftdrop_jesTimeRunDDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupMPFInSituDown", FatJet_pt_jesCorrelationGroupMPFInSituDown, &b_FatJet_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupMPFInSituDown", FatJet_mass_jesCorrelationGroupMPFInSituDown, &b_FatJet_mass_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown", FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown, &b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupIntercalibrationDown", FatJet_pt_jesCorrelationGroupIntercalibrationDown, &b_FatJet_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupIntercalibrationDown", FatJet_mass_jesCorrelationGroupIntercalibrationDown, &b_FatJet_mass_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown", FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown, &b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupbJESDown", FatJet_pt_jesCorrelationGroupbJESDown, &b_FatJet_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupbJESDown", FatJet_mass_jesCorrelationGroupbJESDown, &b_FatJet_mass_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupbJESDown", FatJet_msoftdrop_jesCorrelationGroupbJESDown, &b_FatJet_msoftdrop_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupFlavorDown", FatJet_pt_jesCorrelationGroupFlavorDown, &b_FatJet_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupFlavorDown", FatJet_mass_jesCorrelationGroupFlavorDown, &b_FatJet_mass_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupFlavorDown", FatJet_msoftdrop_jesCorrelationGroupFlavorDown, &b_FatJet_msoftdrop_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("FatJet_pt_jesCorrelationGroupUncorrelatedDown", FatJet_pt_jesCorrelationGroupUncorrelatedDown, &b_FatJet_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("FatJet_mass_jesCorrelationGroupUncorrelatedDown", FatJet_mass_jesCorrelationGroupUncorrelatedDown, &b_FatJet_mass_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown", FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown, &b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown);
}
//#endif // #ifdef Events_cxx
