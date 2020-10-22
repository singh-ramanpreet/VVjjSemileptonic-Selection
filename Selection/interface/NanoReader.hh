#ifndef NanoReader_h
#define NanoReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class NanoReader {
public :
  TTree *fChain;
  int fCurrent_;

  int era_;
  int nanoVersion_;
  bool isMC_;

  NanoReader(TTree *tree, int era, int nanoVersion, int isMC) {

    if (era!=2016 && era!=2017 && era!=2018) {
      //std::cout << "unknown era, assuming 2018" << std::endl;
      era_ = 2018;
    }
    else { era_ = era; }
    
    if (nanoVersion!=6 && nanoVersion!=7) {
      nanoVersion_ = 7;
    }
    else { nanoVersion_ = nanoVersion; }

    if (isMC==1) { isMC_ = true; }
    else { isMC_ = false; }

    Init(tree);

    if (isMC_) InitMC();
    else InitData();

  }

  virtual ~NanoReader() {};

  virtual int GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);

  virtual void Init(TTree *tree);
  virtual void InitMC();
  virtual void InitData();

  // Variables

  UInt_t          run;
  UInt_t          luminosityBlock;
  ULong64_t       event;
  Float_t         btagWeight_CSVV2;
  Float_t         btagWeight_CMVA;
  Float_t         btagWeight_DeepCSVB;
  Float_t         L1PreFiringWeight_Dn;
  Float_t         L1PreFiringWeight_Nom;
  Float_t         L1PreFiringWeight_Up;

  TBranch        *b_run;   //!                                                                                       
  TBranch        *b_luminosityBlock;   //!                                                                           
  TBranch        *b_event;   //!                                                                                     
  TBranch        *b_btagWeight_CSVV2;   //!                                                                          
  TBranch        *b_btagWeight_CMVA;   //!                                                                          
  TBranch        *b_btagWeight_DeepCSVB;   //!
  TBranch        *b_L1PreFiringWeight_Dn;
  TBranch        *b_L1PreFiringWeight_Nom;
  TBranch        *b_L1PreFiringWeight_Up;

  Bool_t          Flag_ecalBadCalibFilterV2;

  TBranch        *b_Flag_ecalBadCalibFilterV2;   //!

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

  // Electrons

  UInt_t          nElectron;
  Float_t         Electron_deltaEtaSC[10];   //[nElectron]                                                                                          
  Float_t         Electron_dr03EcalRecHitSumEt[10];   //[nElectron]                                                                                 
  Float_t         Electron_dr03HcalDepth1TowerSumEt[10];   //[nElectron]                                                                            
  Float_t         Electron_dr03TkSumPt[10];   //[nElectron]                                                                                         
  Float_t         Electron_dr03TkSumPtHEEP[10];   //[nElectron]                                                                                     
  Float_t         Electron_dxy[10];   //[nElectron]                                                                                                 
  Float_t         Electron_dxyErr[10];   //[nElectron]                                                                                              
  Float_t         Electron_dz[10];   //[nElectron]                                                                                                  
  Float_t         Electron_dzErr[10];   //[nElectron]                                                                                               
  Float_t         Electron_eCorr[10];   //[nElectron]                                                                                               
  Float_t         Electron_eInvMinusPInv[10];   //[nElectron]                                                                                       
  Float_t         Electron_energyErr[10];   //[nElectron]                                                                                           
  Float_t         Electron_eta[10];   //[nElectron]                                                                                                 
  Float_t         Electron_hoe[10];   //[nElectron]                                                                                                 
  Float_t         Electron_ip3d[10];   //[nElectron]                                                                                                
  Float_t         Electron_jetPtRelv2[10];   //[nElectron]                                                                                          
  Float_t         Electron_jetRelIso[10];   //[nElectron]                                                                                           
  Float_t         Electron_mass[10];   //[nElectron]                                                                                                
  Float_t         Electron_miniPFRelIso_all[10];   //[nElectron]                                                                                    
  Float_t         Electron_miniPFRelIso_chg[10];   //[nElectron]                                                                                    
  Float_t         Electron_mvaFall17V1Iso[10];   //[nElectron]                                                                                      
  Float_t         Electron_mvaFall17V1noIso[10];   //[nElectron]                                                                                    
  Float_t         Electron_mvaFall17V2Iso[10];   //[nElectron]                                                                                      
  Float_t         Electron_mvaFall17V2noIso[10];   //[nElectron]                                                                                    
  Float_t         Electron_pfRelIso03_all[10];   //[nElectron]                                                                                      
  Float_t         Electron_pfRelIso03_chg[10];   //[nElectron]                                                                                      
  Float_t         Electron_phi[10];   //[nElectron]                                                                                                 
  Float_t         Electron_pt[10];   //[nElectron]                                                                                                  
  Float_t         Electron_r9[10];   //[nElectron]                                                                                                  
  Float_t         Electron_sieie[10];   //[nElectron]                                                                                               
  Float_t         Electron_sip3d[10];   //[nElectron]                                                                                               
  Float_t         Electron_mvaTTH[10];   //[nElectron]
  Int_t           Electron_charge[10];   //[nElectron]                                                                                              
  Int_t           Electron_cutBased[10];   //[nElectron]                                                                                            
  Int_t           Electron_cutBased_Fall17_V1[10];   //[nElectron]                                                                                  
  Int_t           Electron_jetIdx[10];   //[nElectron]                                                                                              
  Int_t           Electron_pdgId[10];   //[nElectron]                                                                                               
  Int_t           Electron_photonIdx[10];   //[nElectron]                                                                                           
  Int_t           Electron_tightCharge[10];   //[nElectron]                                                                                         
  Int_t           Electron_vidNestedWPBitmap[10];   //[nElectron]                                                                                   
  Int_t           Electron_vidNestedWPBitmapHEEP[10];   //[nElectron]                                                                               
  Bool_t          Electron_convVeto[10];   //[nElectron]                                                                                            
  Bool_t          Electron_cutBased_HEEP[10];   //[nElectron]                                                                                       
  Bool_t          Electron_isPFcand[10];   //[nElectron]                                                                                            
  UChar_t         Electron_lostHits[10];   //[nElectron]                                                                                            
  Bool_t          Electron_mvaFall17V1Iso_WP80[10];   //[nElectron]                                                                                 
  Bool_t          Electron_mvaFall17V1Iso_WP90[10];   //[nElectron]                                                                                 
  Bool_t          Electron_mvaFall17V1Iso_WPL[10];   //[nElectron]                                                                                  
  Bool_t          Electron_mvaFall17V1noIso_WP80[10];   //[nElectron]                                                                               
  Bool_t          Electron_mvaFall17V1noIso_WP90[10];   //[nElectron]                                                                               
  Bool_t          Electron_mvaFall17V1noIso_WPL[10];   //[nElectron]                                                                                
  Bool_t          Electron_mvaFall17V2Iso_WP80[10];   //[nElectron]                                                                                 
  Bool_t          Electron_mvaFall17V2Iso_WP90[10];   //[nElectron]                                                                                 
  Bool_t          Electron_mvaFall17V2Iso_WPL[10];   //[nElectron]                                                                                  
  Bool_t          Electron_mvaFall17V2noIso_WP80[10];   //[nElectron]                                                                               
  Bool_t          Electron_mvaFall17V2noIso_WP90[10];   //[nElectron]                                                                               
  Bool_t          Electron_mvaFall17V2noIso_WPL[10];   //[nElectron]                                                                                
  UChar_t         Electron_seedGain[10];   //[nElectron]

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

  // FatJets

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

   // Jets

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
 
  // MET

  Float_t         MET_MetUnclustEnUpDeltaX;
  Float_t         MET_MetUnclustEnUpDeltaY;
  Float_t         MET_covXX;
  Float_t         MET_covXY;
  Float_t         MET_covYY;
  Float_t         MET_phi;
  Float_t         MET_pt;
  Float_t         MET_significance;
  Float_t         MET_sumEt;
  Float_t         METFixEE2017_phi;
  Float_t         METFixEE2017_pt;
  Float_t         METFixEE2017_significance;
  Float_t         METFixEE2017_sumEt;

  TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!                                                                  
  TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!                                                                  
  TBranch        *b_MET_covXX;   //!                                                                                 
  TBranch        *b_MET_covXY;   //!                                                                                 
  TBranch        *b_MET_covYY;   //!                                                                                 
  TBranch        *b_MET_phi;   //!                                                                                   
  TBranch        *b_MET_pt;   //!                                                                                    
  TBranch        *b_MET_significance;   //!                                                                          
  TBranch        *b_MET_sumEt;   //!
  TBranch        *b_METFixEE2017_phi;   //!                                                                              
  TBranch        *b_METFixEE2017_pt;   //!                                                                               
  TBranch        *b_METFixEE2017_significance;   //!                                                                     
  TBranch        *b_METFixEE2017_sumEt;   //!                                                                            

  Float_t         CaloMET_phi;
  Float_t         CaloMET_pt;
  Float_t         CaloMET_sumEt;
  Float_t         ChsMET_phi;
  Float_t         ChsMET_pt;
  Float_t         ChsMET_sumEt;

  TBranch        *b_CaloMET_phi;   //!                                                                               
  TBranch        *b_CaloMET_pt;   //!                                                                                
  TBranch        *b_CaloMET_sumEt;   //!                                                                             
  TBranch        *b_ChsMET_phi;   //!                                                                                
  TBranch        *b_ChsMET_pt;   //!                                                                                 
  TBranch        *b_ChsMET_sumEt;   //!

  Float_t         PuppiMET_phi;
  Float_t         PuppiMET_pt;
  Float_t         PuppiMET_sumEt;
  Float_t         RawMET_phi;
  Float_t         RawMET_pt;
  Float_t         RawMET_sumEt;

  TBranch        *b_PuppiMET_phi;   //!                                                                              
  TBranch        *b_PuppiMET_pt;   //!                                                                               
  TBranch        *b_PuppiMET_sumEt;   //!                                                                            
  TBranch        *b_RawMET_phi;   //!                                                                                
  TBranch        *b_RawMET_pt;   //!                                                                                 
  TBranch        *b_RawMET_sumEt;   //!

  Float_t         TkMET_phi;
  Float_t         TkMET_pt;
  Float_t         TkMET_sumEt;

  TBranch        *b_TkMET_phi;   //!                                                                                 
  TBranch        *b_TkMET_pt;   //!                                                                                  
  TBranch        *b_TkMET_sumEt;   //!

  // Muons

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

  //Photons

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

  // Pileup

  Float_t         Pileup_nTrueInt;
  Float_t         Pileup_pudensity;
  Float_t         Pileup_gpudensity;
  Int_t           Pileup_nPU;
  Int_t           Pileup_sumEOOT;
  Int_t           Pileup_sumLOOT;
  Float_t         fixedGridRhoFastjetAll;
  Float_t         fixedGridRhoFastjetCentral;
  Float_t         fixedGridRhoFastjetCentralCalo;
  Float_t         fixedGridRhoFastjetCentralChargedPileUp;
  Float_t         fixedGridRhoFastjetCentralNeutral;

  TBranch        *b_Pileup_nTrueInt;   //!                                                                           
  TBranch        *b_Pileup_pudensity;   //!                                                                          
  TBranch        *b_Pileup_gpudensity;   //!                                                                         
  TBranch        *b_Pileup_nPU;   //!                                                                                
  TBranch        *b_Pileup_sumEOOT;   //!                                                                            
  TBranch        *b_Pileup_sumLOOT;   //!
  TBranch        *b_fixedGridRhoFastjetAll;   //!                                                                    
  TBranch        *b_fixedGridRhoFastjetCentral;   //!                                                                
  TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!                                                            
  TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!                                                   
  TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!

  // Subjets

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

  // PVs

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

  // Triggers

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

  // SingleMuon triggers

  Bool_t HLT_IsoMu22;
  Bool_t HLT_IsoTkMu22;
  Bool_t HLT_IsoMu24;
  Bool_t HLT_IsoTkMu24;
  Bool_t HLT_IsoMu27;
  Bool_t HLT_IsoMu30;
  Bool_t HLT_Mu50;

  TBranch  *b_HLT_IsoMu22;
  TBranch  *b_HLT_IsoTkMu22;
  TBranch  *b_HLT_IsoMu24;
  TBranch  *b_HLT_IsoTkMu24;
  TBranch  *b_HLT_IsoMu27;
  TBranch  *b_HLT_IsoMu30;
  TBranch  *b_HLT_Mu50;

  // DoubleMuon triggers

  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  Bool_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Bool_t HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
  Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
  Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;

  TBranch  *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  TBranch  *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  TBranch  *b_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  TBranch  *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  TBranch  *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
  TBranch  *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
  TBranch  *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;

  // SingleElectron triggers

  Bool_t HLT_Ele25_eta2p1_WPTight_Gsf;
  Bool_t HLT_Ele27_eta2p1_WPLoose_Gsf;
  Bool_t HLT_Ele27_WPTight_Gsf;
  Bool_t HLT_Ele28_WPTight_Gsf;
  Bool_t HLT_Ele32_WPTight_Gsf;
  Bool_t HLT_Ele35_WPLooose_Gsf;
  Bool_t HLT_Ele35_WPTight_Gsf;
  Bool_t HLT_Ele38_WPTight_Gsf;
  Bool_t HLT_Ele40_WPTight_Gsf;

  TBranch  *b_HLT_Ele25_eta2p1_WPTight_Gsf;
  TBranch  *b_HLT_Ele27_eta2p1_WPLoose_Gsf;
  TBranch  *b_HLT_Ele27_WPTight_Gsf;
  TBranch  *b_HLT_Ele28_WPTight_Gsf;
  TBranch  *b_HLT_Ele32_WPTight_Gsf;
  TBranch  *b_HLT_Ele35_WPLooose_Gsf;
  TBranch  *b_HLT_Ele35_WPTight_Gsf;
  TBranch  *b_HLT_Ele38_WPTight_Gsf;
  TBranch  *b_HLT_Ele40_WPTight_Gsf;

  // DoubleElectron triggers

  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  Bool_t HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
  Bool_t HLT_DoubleEle25_CaloIdL_MW;
  Bool_t HLT_DoubleEle27_CaloIdL_MW;
  Bool_t HLT_DoubleEle33_CaloIdL_MW;

  TBranch  *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  TBranch  *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  TBranch  *b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
  TBranch  *b_HLT_DoubleEle25_CaloIdL_MW;
  TBranch  *b_HLT_DoubleEle27_CaloIdL_MW;
  TBranch  *b_HLT_DoubleEle33_CaloIdL_MW;

  // Corrected objects

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

  // MC only variables

  Float_t         Generator_binvar;
  Float_t         Generator_scalePDF;
  Float_t         Generator_weight;
  Float_t         Generator_x1;
  Float_t         Generator_x2;
  Float_t         Generator_xpdf1;
  Float_t         Generator_xpdf2;
  Int_t           Generator_id1;
  Int_t           Generator_id2;
  Float_t         genWeight;
  Float_t         LHEWeight_originalXWGTUP;
  UInt_t          nLHEPdfWeight;
  Float_t         LHEPdfWeight[200];   //[nLHEPdfWeight]                                                               
  UInt_t          nLHEReweightingWeight;
  Float_t         LHEReweightingWeight[1000];   //[nLHEReweightingWeight]                                               
  UInt_t          nLHEScaleWeight;
  Float_t         LHEScaleWeight[50];   //[nLHEScaleWeight]                                                           
  UInt_t          nPSWeight;
  Float_t         PSWeight[4];   //[nPSWeight]

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

  TBranch  *b_Jet_pt_jerUp;                                                                       
  TBranch  *b_Jet_mass_jerUp;                                                                     
  TBranch  *b_MET_pt_jerUp;
  TBranch  *b_MET_phi_jerUp;
  TBranch  *b_Jet_pt_jesAbsoluteStatUp;                                                           
  TBranch  *b_Jet_mass_jesAbsoluteStatUp;                                                         
  TBranch  *b_MET_pt_jesAbsoluteStatUp;
  TBranch  *b_MET_phi_jesAbsoluteStatUp;
  TBranch  *b_Jet_pt_jesAbsoluteScaleUp;                                                          
  TBranch  *b_Jet_mass_jesAbsoluteScaleUp;                                                        
  TBranch  *b_MET_pt_jesAbsoluteScaleUp;
  TBranch  *b_MET_phi_jesAbsoluteScaleUp;
  TBranch  *b_Jet_pt_jesAbsoluteSampleUp;                                                         
  TBranch  *b_Jet_mass_jesAbsoluteSampleUp;                                                       
  TBranch  *b_MET_pt_jesAbsoluteSampleUp;
  TBranch  *b_MET_phi_jesAbsoluteSampleUp;
  TBranch  *b_Jet_pt_jesAbsoluteFlavMapUp;                                                        
  TBranch  *b_Jet_mass_jesAbsoluteFlavMapUp;                                                      
  TBranch  *b_MET_pt_jesAbsoluteFlavMapUp;
  TBranch  *b_MET_phi_jesAbsoluteFlavMapUp;
  TBranch  *b_Jet_pt_jesAbsoluteMPFBiasUp;                                                        
  TBranch  *b_Jet_mass_jesAbsoluteMPFBiasUp;                                                      
  TBranch  *b_MET_pt_jesAbsoluteMPFBiasUp;
  TBranch  *b_MET_phi_jesAbsoluteMPFBiasUp;
  TBranch  *b_Jet_pt_jesFragmentationUp;                                                          
  TBranch  *b_Jet_mass_jesFragmentationUp;                                                        
  TBranch  *b_MET_pt_jesFragmentationUp;
  TBranch  *b_MET_phi_jesFragmentationUp;
  TBranch  *b_Jet_pt_jesSinglePionECALUp;                                                         
  TBranch  *b_Jet_mass_jesSinglePionECALUp;
  TBranch  *b_MET_pt_jesSinglePionECALUp;
  TBranch  *b_MET_phi_jesSinglePionECALUp;
  TBranch  *b_Jet_pt_jesSinglePionHCALUp;                                                         
  TBranch  *b_Jet_mass_jesSinglePionHCALUp;                                                       
  TBranch  *b_MET_pt_jesSinglePionHCALUp;
  TBranch  *b_MET_phi_jesSinglePionHCALUp;
  TBranch  *b_Jet_pt_jesFlavorQCDUp;                                                              
  TBranch  *b_Jet_mass_jesFlavorQCDUp;                                                            
  TBranch  *b_MET_pt_jesFlavorQCDUp;
  TBranch  *b_MET_phi_jesFlavorQCDUp;
  TBranch  *b_Jet_pt_jesTimePtEtaUp;                                                              
  TBranch  *b_Jet_mass_jesTimePtEtaUp;                                                            
  TBranch  *b_MET_pt_jesTimePtEtaUp;
  TBranch  *b_MET_phi_jesTimePtEtaUp;
  TBranch  *b_Jet_pt_jesRelativeJEREC1Up;                                                         
  TBranch  *b_Jet_mass_jesRelativeJEREC1Up;                                                       
  TBranch  *b_MET_pt_jesRelativeJEREC1Up;
  TBranch  *b_MET_phi_jesRelativeJEREC1Up;
  TBranch  *b_Jet_pt_jesRelativeJEREC2Up;                                                         
  TBranch  *b_Jet_mass_jesRelativeJEREC2Up;                                                       
  TBranch  *b_MET_pt_jesRelativeJEREC2Up;
  TBranch  *b_MET_phi_jesRelativeJEREC2Up;
  TBranch  *b_Jet_pt_jesRelativeJERHFUp;                                                          
  TBranch  *b_Jet_mass_jesRelativeJERHFUp;
  TBranch  *b_MET_pt_jesRelativeJERHFUp;
  TBranch  *b_MET_phi_jesRelativeJERHFUp;
  TBranch  *b_Jet_pt_jesRelativePtBBUp;                                                           
  TBranch  *b_Jet_mass_jesRelativePtBBUp;                                                         
  TBranch  *b_MET_pt_jesRelativePtBBUp;
  TBranch  *b_MET_phi_jesRelativePtBBUp;
  TBranch  *b_Jet_pt_jesRelativePtEC1Up;                                                          
  TBranch  *b_Jet_mass_jesRelativePtEC1Up;                                                        
  TBranch  *b_MET_pt_jesRelativePtEC1Up;
  TBranch  *b_MET_phi_jesRelativePtEC1Up;
  TBranch  *b_Jet_pt_jesRelativePtEC2Up;                                                          
  TBranch  *b_Jet_mass_jesRelativePtEC2Up;                                                        
  TBranch  *b_MET_pt_jesRelativePtEC2Up;
  TBranch  *b_MET_phi_jesRelativePtEC2Up;
  TBranch  *b_Jet_pt_jesRelativePtHFUp;                                                           
  TBranch  *b_Jet_mass_jesRelativePtHFUp;                                                         
  TBranch  *b_MET_pt_jesRelativePtHFUp;
  TBranch  *b_MET_phi_jesRelativePtHFUp;
  TBranch  *b_Jet_pt_jesRelativeBalUp;                                                            
  TBranch  *b_Jet_mass_jesRelativeBalUp;                                                          
  TBranch  *b_MET_pt_jesRelativeBalUp;
  TBranch  *b_MET_phi_jesRelativeBalUp;
  TBranch  *b_Jet_pt_jesRelativeSampleUp;                                                         
  TBranch  *b_Jet_mass_jesRelativeSampleUp;                                                       
  TBranch  *b_MET_pt_jesRelativeSampleUp;
  TBranch  *b_MET_phi_jesRelativeSampleUp;
  TBranch  *b_Jet_pt_jesRelativeFSRUp;                                                            
  TBranch  *b_Jet_mass_jesRelativeFSRUp;                                                          
  TBranch  *b_MET_pt_jesRelativeFSRUp;
  TBranch  *b_MET_phi_jesRelativeFSRUp;
  TBranch  *b_Jet_pt_jesRelativeStatFSRUp;                                                        
  TBranch  *b_Jet_mass_jesRelativeStatFSRUp;                                                      
  TBranch  *b_MET_pt_jesRelativeStatFSRUp;
  TBranch  *b_MET_phi_jesRelativeStatFSRUp;
  TBranch  *b_Jet_pt_jesRelativeStatECUp;                                                         
  TBranch  *b_Jet_mass_jesRelativeStatECUp;                                                       
  TBranch  *b_MET_pt_jesRelativeStatECUp;
  TBranch  *b_MET_phi_jesRelativeStatECUp;
  TBranch  *b_Jet_pt_jesRelativeStatHFUp;                                                         
  TBranch  *b_Jet_mass_jesRelativeStatHFUp;                                                       
  TBranch  *b_MET_pt_jesRelativeStatHFUp;
  TBranch  *b_MET_phi_jesRelativeStatHFUp;
  TBranch  *b_Jet_pt_jesPileUpDataMCUp;                                                           
  TBranch  *b_Jet_mass_jesPileUpDataMCUp;                                                         
  TBranch  *b_MET_pt_jesPileUpDataMCUp;
  TBranch  *b_MET_phi_jesPileUpDataMCUp;
  TBranch  *b_Jet_pt_jesPileUpPtRefUp;                                                            
  TBranch  *b_Jet_mass_jesPileUpPtRefUp;                                                          
  TBranch  *b_MET_pt_jesPileUpPtRefUp;
  TBranch  *b_MET_phi_jesPileUpPtRefUp;
  TBranch  *b_Jet_pt_jesPileUpPtBBUp;                                                             
  TBranch  *b_Jet_mass_jesPileUpPtBBUp;                                                           
  TBranch  *b_MET_pt_jesPileUpPtBBUp;
  TBranch  *b_MET_phi_jesPileUpPtBBUp;
  TBranch  *b_Jet_pt_jesPileUpPtEC1Up;                                                            
  TBranch  *b_Jet_mass_jesPileUpPtEC1Up;                                                          
  TBranch  *b_MET_pt_jesPileUpPtEC1Up;
  TBranch  *b_MET_phi_jesPileUpPtEC1Up;
  TBranch  *b_Jet_pt_jesPileUpPtEC2Up;                                                            
  TBranch  *b_Jet_mass_jesPileUpPtEC2Up;                                                          
  TBranch  *b_MET_pt_jesPileUpPtEC2Up;
  TBranch  *b_MET_phi_jesPileUpPtEC2Up;
  TBranch  *b_Jet_pt_jesPileUpPtHFUp;                                                             
  TBranch  *b_Jet_mass_jesPileUpPtHFUp;                                                           
  TBranch  *b_MET_pt_jesPileUpPtHFUp;
  TBranch  *b_MET_phi_jesPileUpPtHFUp;
  TBranch  *b_Jet_pt_jesPileUpMuZeroUp;                                                           
  TBranch  *b_Jet_mass_jesPileUpMuZeroUp; 
  TBranch  *b_MET_pt_jesPileUpMuZeroUp;
  TBranch  *b_MET_phi_jesPileUpMuZeroUp;
  TBranch  *b_Jet_pt_jesPileUpEnvelopeUp;                                                         
  TBranch  *b_Jet_mass_jesPileUpEnvelopeUp;                                                       
  TBranch  *b_MET_pt_jesPileUpEnvelopeUp;
  TBranch  *b_MET_phi_jesPileUpEnvelopeUp;
  TBranch  *b_Jet_pt_jesSubTotalPileUpUp;                                                         
  TBranch  *b_Jet_mass_jesSubTotalPileUpUp;                                                       
  TBranch  *b_MET_pt_jesSubTotalPileUpUp;
  TBranch  *b_MET_phi_jesSubTotalPileUpUp;
  TBranch  *b_Jet_pt_jesSubTotalRelativeUp;                                                       
  TBranch  *b_Jet_mass_jesSubTotalRelativeUp;                                                     
  TBranch  *b_MET_pt_jesSubTotalRelativeUp;
  TBranch  *b_MET_phi_jesSubTotalRelativeUp;
  TBranch  *b_Jet_pt_jesSubTotalPtUp;                                                             
  TBranch  *b_Jet_mass_jesSubTotalPtUp;                                                           
  TBranch  *b_MET_pt_jesSubTotalPtUp;
  TBranch  *b_MET_phi_jesSubTotalPtUp;
  TBranch  *b_Jet_pt_jesSubTotalScaleUp;                                                          
  TBranch  *b_Jet_mass_jesSubTotalScaleUp;                                                        
  TBranch  *b_MET_pt_jesSubTotalScaleUp;
  TBranch  *b_MET_phi_jesSubTotalScaleUp;
  TBranch  *b_Jet_pt_jesSubTotalAbsoluteUp;                                                       
  TBranch  *b_Jet_mass_jesSubTotalAbsoluteUp;                                                     
  TBranch  *b_MET_pt_jesSubTotalAbsoluteUp;
  TBranch  *b_MET_phi_jesSubTotalAbsoluteUp;
  TBranch  *b_Jet_pt_jesSubTotalMCUp;                                                             
  TBranch  *b_Jet_mass_jesSubTotalMCUp;                                                           
  TBranch  *b_MET_pt_jesSubTotalMCUp;
  TBranch  *b_MET_phi_jesSubTotalMCUp;
  TBranch  *b_Jet_pt_jesTotalUp;                                                                  
  TBranch  *b_Jet_mass_jesTotalUp;
  TBranch  *b_MET_pt_jesTotalUp;
  TBranch  *b_MET_phi_jesTotalUp;
  TBranch  *b_Jet_pt_jesTotalNoFlavorUp;                                                          
  TBranch  *b_Jet_mass_jesTotalNoFlavorUp;                                                        
  TBranch  *b_MET_pt_jesTotalNoFlavorUp;
  TBranch  *b_MET_phi_jesTotalNoFlavorUp;
  TBranch  *b_Jet_pt_jesTotalNoTimeUp;                                                            
  TBranch  *b_Jet_mass_jesTotalNoTimeUp;                                                          
  TBranch  *b_MET_pt_jesTotalNoTimeUp;
  TBranch  *b_MET_phi_jesTotalNoTimeUp;
  TBranch  *b_Jet_pt_jesTotalNoFlavorNoTimeUp;                                                    
  TBranch  *b_Jet_mass_jesTotalNoFlavorNoTimeUp;                                                  
  TBranch  *b_MET_pt_jesTotalNoFlavorNoTimeUp;
  TBranch  *b_MET_phi_jesTotalNoFlavorNoTimeUp;
  TBranch  *b_Jet_pt_jesFlavorZJetUp;                                                             
  TBranch  *b_Jet_mass_jesFlavorZJetUp;                                                           
  TBranch  *b_MET_pt_jesFlavorZJetUp;
  TBranch  *b_MET_phi_jesFlavorZJetUp;
  TBranch  *b_Jet_pt_jesFlavorPhotonJetUp;                                                        
  TBranch  *b_Jet_mass_jesFlavorPhotonJetUp;                                                      
  TBranch  *b_MET_pt_jesFlavorPhotonJetUp;
  TBranch  *b_MET_phi_jesFlavorPhotonJetUp;
  TBranch  *b_Jet_pt_jesFlavorPureGluonUp;                                                        
  TBranch  *b_Jet_mass_jesFlavorPureGluonUp;                                                      
  TBranch  *b_MET_pt_jesFlavorPureGluonUp;
  TBranch  *b_MET_phi_jesFlavorPureGluonUp;
  TBranch  *b_Jet_pt_jesFlavorPureQuarkUp;                                                        
  TBranch  *b_Jet_mass_jesFlavorPureQuarkUp;                                                      
  TBranch  *b_MET_pt_jesFlavorPureQuarkUp;
  TBranch  *b_MET_phi_jesFlavorPureQuarkUp;
  TBranch  *b_Jet_pt_jesFlavorPureCharmUp;                                                        
  TBranch  *b_Jet_mass_jesFlavorPureCharmUp;
  TBranch  *b_MET_pt_jesFlavorPureCharmUp;
  TBranch  *b_MET_phi_jesFlavorPureCharmUp;
  TBranch  *b_Jet_pt_jesFlavorPureBottomUp;                                                       
  TBranch  *b_Jet_mass_jesFlavorPureBottomUp;                                                     
  TBranch  *b_MET_pt_jesFlavorPureBottomUp;
  TBranch  *b_MET_phi_jesFlavorPureBottomUp;
  TBranch  *b_Jet_pt_jesTimeRunAUp;                                                               
  TBranch  *b_Jet_mass_jesTimeRunAUp;                                                             
  TBranch  *b_MET_pt_jesTimeRunAUp;
  TBranch  *b_MET_phi_jesTimeRunAUp;
  TBranch  *b_Jet_pt_jesTimeRunBUp;                                                               
  TBranch  *b_Jet_mass_jesTimeRunBUp;                                                             
  TBranch  *b_MET_pt_jesTimeRunBUp;
  TBranch  *b_MET_phi_jesTimeRunBUp;
  TBranch  *b_Jet_pt_jesTimeRunCUp;                                                               
  TBranch  *b_Jet_mass_jesTimeRunCUp;                                                             
  TBranch  *b_MET_pt_jesTimeRunCUp;
  TBranch  *b_MET_phi_jesTimeRunCUp;
  TBranch  *b_Jet_pt_jesTimeRunDUp;                                                               
  TBranch  *b_Jet_mass_jesTimeRunDUp;                                                             
  TBranch  *b_MET_pt_jesTimeRunDUp;
  TBranch  *b_MET_phi_jesTimeRunDUp;
  TBranch  *b_Jet_pt_jesCorrelationGroupMPFInSituUp;                                              
  TBranch  *b_Jet_mass_jesCorrelationGroupMPFInSituUp;                                            
  TBranch  *b_MET_pt_jesCorrelationGroupMPFInSituUp;
  TBranch  *b_MET_phi_jesCorrelationGroupMPFInSituUp;
  TBranch  *b_Jet_pt_jesCorrelationGroupIntercalibrationUp;                                       
  TBranch  *b_Jet_mass_jesCorrelationGroupIntercalibrationUp;                                     
  TBranch  *b_MET_pt_jesCorrelationGroupIntercalibrationUp;
  TBranch  *b_MET_phi_jesCorrelationGroupIntercalibrationUp;
  TBranch  *b_Jet_pt_jesCorrelationGroupbJESUp;                                                   
  TBranch  *b_Jet_mass_jesCorrelationGroupbJESUp;                                                 
  TBranch  *b_MET_pt_jesCorrelationGroupbJESUp;
  TBranch  *b_MET_phi_jesCorrelationGroupbJESUp;
  TBranch  *b_Jet_pt_jesCorrelationGroupFlavorUp;                                                 
  TBranch  *b_Jet_mass_jesCorrelationGroupFlavorUp;                                               
  TBranch  *b_MET_pt_jesCorrelationGroupFlavorUp;
  TBranch  *b_MET_phi_jesCorrelationGroupFlavorUp;
  TBranch  *b_Jet_pt_jesCorrelationGroupUncorrelatedUp;                                           
  TBranch  *b_Jet_mass_jesCorrelationGroupUncorrelatedUp;                                         
  TBranch  *b_MET_pt_jesCorrelationGroupUncorrelatedUp;
  TBranch  *b_MET_phi_jesCorrelationGroupUncorrelatedUp;
  TBranch  *b_MET_pt_unclustEnUp;
  TBranch  *b_MET_phi_unclustEnUp;
  TBranch  *b_Jet_pt_jerDown;                                                                     
  TBranch  *b_Jet_mass_jerDown;                                                                   
  TBranch  *b_MET_pt_jerDown;
  TBranch  *b_MET_phi_jerDown;
  TBranch  *b_Jet_pt_jesAbsoluteStatDown;                                                         
  TBranch  *b_Jet_mass_jesAbsoluteStatDown;                                                       
  TBranch  *b_MET_pt_jesAbsoluteStatDown;
  TBranch  *b_MET_phi_jesAbsoluteStatDown;
  TBranch  *b_Jet_pt_jesAbsoluteScaleDown;                                                        
  TBranch  *b_Jet_mass_jesAbsoluteScaleDown;                                                      
  TBranch  *b_MET_pt_jesAbsoluteScaleDown;
  TBranch  *b_MET_phi_jesAbsoluteScaleDown;
  TBranch  *b_Jet_pt_jesAbsoluteSampleDown;                                                       
  TBranch  *b_Jet_mass_jesAbsoluteSampleDown;                                                     
  TBranch  *b_MET_pt_jesAbsoluteSampleDown;
  TBranch  *b_MET_phi_jesAbsoluteSampleDown;
  TBranch  *b_Jet_pt_jesAbsoluteFlavMapDown;                                                      
  TBranch  *b_Jet_mass_jesAbsoluteFlavMapDown;                                                    
  TBranch  *b_MET_pt_jesAbsoluteFlavMapDown;
  TBranch  *b_MET_phi_jesAbsoluteFlavMapDown;
  TBranch  *b_Jet_pt_jesAbsoluteMPFBiasDown;                                                      
  TBranch  *b_Jet_mass_jesAbsoluteMPFBiasDown;                                                    
  TBranch  *b_MET_pt_jesAbsoluteMPFBiasDown;
  TBranch  *b_MET_phi_jesAbsoluteMPFBiasDown;
  TBranch  *b_Jet_pt_jesFragmentationDown;                                                        
  TBranch  *b_Jet_mass_jesFragmentationDown;                                                      
  TBranch  *b_MET_pt_jesFragmentationDown;
  TBranch  *b_MET_phi_jesFragmentationDown;
  TBranch  *b_Jet_pt_jesSinglePionECALDown;                                                       
  TBranch  *b_Jet_mass_jesSinglePionECALDown;                                                     
  TBranch  *b_MET_pt_jesSinglePionECALDown;
  TBranch  *b_MET_phi_jesSinglePionECALDown;
  TBranch  *b_Jet_pt_jesSinglePionHCALDown;                                                       
  TBranch  *b_Jet_mass_jesSinglePionHCALDown;                                                     
  TBranch  *b_MET_pt_jesSinglePionHCALDown;
  TBranch  *b_MET_phi_jesSinglePionHCALDown;
  TBranch  *b_Jet_pt_jesFlavorQCDDown;                                                            
  TBranch  *b_Jet_mass_jesFlavorQCDDown;                                                          
  TBranch  *b_MET_pt_jesFlavorQCDDown;
  TBranch  *b_MET_phi_jesFlavorQCDDown;
  TBranch  *b_Jet_pt_jesTimePtEtaDown;                                                            
  TBranch  *b_Jet_mass_jesTimePtEtaDown;                                                          
  TBranch  *b_MET_pt_jesTimePtEtaDown;
  TBranch  *b_MET_phi_jesTimePtEtaDown;
  TBranch  *b_Jet_pt_jesRelativeJEREC1Down;                                                       
  TBranch  *b_Jet_mass_jesRelativeJEREC1Down;                                                     
  TBranch  *b_MET_pt_jesRelativeJEREC1Down;
  TBranch  *b_MET_phi_jesRelativeJEREC1Down;
  TBranch  *b_Jet_pt_jesRelativeJEREC2Down;                                                       
  TBranch  *b_Jet_mass_jesRelativeJEREC2Down;                                                     
  TBranch  *b_MET_pt_jesRelativeJEREC2Down;
  TBranch  *b_MET_phi_jesRelativeJEREC2Down;
  TBranch  *b_Jet_pt_jesRelativeJERHFDown;                                                        
  TBranch  *b_Jet_mass_jesRelativeJERHFDown;                                                      
  TBranch  *b_MET_pt_jesRelativeJERHFDown;
  TBranch  *b_MET_phi_jesRelativeJERHFDown;
  TBranch  *b_Jet_pt_jesRelativePtBBDown;                                                         
  TBranch  *b_Jet_mass_jesRelativePtBBDown;                                                       
  TBranch  *b_MET_pt_jesRelativePtBBDown;
  TBranch  *b_MET_phi_jesRelativePtBBDown;
  TBranch  *b_Jet_pt_jesRelativePtEC1Down;                                                        
  TBranch  *b_Jet_mass_jesRelativePtEC1Down;                                                      
  TBranch  *b_MET_pt_jesRelativePtEC1Down;
  TBranch  *b_MET_phi_jesRelativePtEC1Down;
  TBranch  *b_Jet_pt_jesRelativePtEC2Down;                                                        
  TBranch  *b_Jet_mass_jesRelativePtEC2Down;                                                      
  TBranch  *b_MET_pt_jesRelativePtEC2Down;
  TBranch  *b_MET_phi_jesRelativePtEC2Down;
  TBranch  *b_Jet_pt_jesRelativePtHFDown;                                                         
  TBranch  *b_Jet_mass_jesRelativePtHFDown;                                                       
  TBranch  *b_MET_pt_jesRelativePtHFDown;
  TBranch  *b_MET_phi_jesRelativePtHFDown;
  TBranch  *b_Jet_pt_jesRelativeBalDown;                                                          
  TBranch  *b_Jet_mass_jesRelativeBalDown;                                                        
  TBranch  *b_MET_pt_jesRelativeBalDown;
  TBranch  *b_MET_phi_jesRelativeBalDown;
  TBranch  *b_Jet_pt_jesRelativeSampleDown;                                                       
  TBranch  *b_Jet_mass_jesRelativeSampleDown;                                                     
  TBranch  *b_MET_pt_jesRelativeSampleDown;
  TBranch  *b_MET_phi_jesRelativeSampleDown;
  TBranch  *b_Jet_pt_jesRelativeFSRDown;
  TBranch  *b_Jet_mass_jesRelativeFSRDown;                                                        
  TBranch  *b_MET_pt_jesRelativeFSRDown;
  TBranch  *b_MET_phi_jesRelativeFSRDown;
  TBranch  *b_Jet_pt_jesRelativeStatFSRDown;                                                      
  TBranch  *b_Jet_mass_jesRelativeStatFSRDown;                                                    
  TBranch  *b_MET_pt_jesRelativeStatFSRDown;
  TBranch  *b_MET_phi_jesRelativeStatFSRDown;
  TBranch  *b_Jet_pt_jesRelativeStatECDown;                                                       
  TBranch  *b_Jet_mass_jesRelativeStatECDown;                                                     
  TBranch  *b_MET_pt_jesRelativeStatECDown;
  TBranch  *b_MET_phi_jesRelativeStatECDown;
  TBranch  *b_Jet_pt_jesRelativeStatHFDown;                                                       
  TBranch  *b_Jet_mass_jesRelativeStatHFDown;                                                     
  TBranch  *b_MET_pt_jesRelativeStatHFDown;
  TBranch  *b_MET_phi_jesRelativeStatHFDown;
  TBranch  *b_Jet_pt_jesPileUpDataMCDown;                                                         
  TBranch  *b_Jet_mass_jesPileUpDataMCDown;                                                       
  TBranch  *b_MET_pt_jesPileUpDataMCDown;
  TBranch  *b_MET_phi_jesPileUpDataMCDown;
  TBranch  *b_Jet_pt_jesPileUpPtRefDown;                                                          
  TBranch  *b_Jet_mass_jesPileUpPtRefDown;                                                        
  TBranch  *b_MET_pt_jesPileUpPtRefDown;
  TBranch  *b_MET_phi_jesPileUpPtRefDown;
  TBranch  *b_Jet_pt_jesPileUpPtBBDown;                                                           
  TBranch  *b_Jet_mass_jesPileUpPtBBDown;                                                         
  TBranch  *b_MET_pt_jesPileUpPtBBDown;
  TBranch  *b_MET_phi_jesPileUpPtBBDown;
  TBranch  *b_Jet_pt_jesPileUpPtEC1Down;                                                          
  TBranch  *b_Jet_mass_jesPileUpPtEC1Down;                                                        
  TBranch  *b_MET_pt_jesPileUpPtEC1Down;
  TBranch  *b_MET_phi_jesPileUpPtEC1Down;
  TBranch  *b_Jet_pt_jesPileUpPtEC2Down;                                                          
  TBranch  *b_Jet_mass_jesPileUpPtEC2Down;                                                        
  TBranch  *b_MET_pt_jesPileUpPtEC2Down;
  TBranch  *b_MET_phi_jesPileUpPtEC2Down;
  TBranch  *b_Jet_pt_jesPileUpPtHFDown;                                                           
  TBranch  *b_Jet_mass_jesPileUpPtHFDown;                                                         
  TBranch  *b_MET_pt_jesPileUpPtHFDown;
  TBranch  *b_MET_phi_jesPileUpPtHFDown;
  TBranch  *b_Jet_pt_jesPileUpMuZeroDown;                                                         
  TBranch  *b_Jet_mass_jesPileUpMuZeroDown;                                                       
  TBranch  *b_MET_pt_jesPileUpMuZeroDown;
  TBranch  *b_MET_phi_jesPileUpMuZeroDown;
  TBranch  *b_Jet_pt_jesPileUpEnvelopeDown;                                                       
  TBranch  *b_Jet_mass_jesPileUpEnvelopeDown;                                                     
  TBranch  *b_MET_pt_jesPileUpEnvelopeDown;
  TBranch  *b_MET_phi_jesPileUpEnvelopeDown;
  TBranch  *b_Jet_pt_jesSubTotalPileUpDown;                                                       
  TBranch  *b_Jet_mass_jesSubTotalPileUpDown;                                                     
  TBranch  *b_MET_pt_jesSubTotalPileUpDown;
  TBranch  *b_MET_phi_jesSubTotalPileUpDown;
  TBranch  *b_Jet_pt_jesSubTotalRelativeDown;                                                     
  TBranch  *b_Jet_mass_jesSubTotalRelativeDown;                                                   
  TBranch  *b_MET_pt_jesSubTotalRelativeDown;
  TBranch  *b_MET_phi_jesSubTotalRelativeDown;
  TBranch  *b_Jet_pt_jesSubTotalPtDown;                                                           
  TBranch  *b_Jet_mass_jesSubTotalPtDown;                                                         
  TBranch  *b_MET_pt_jesSubTotalPtDown;
  TBranch  *b_MET_phi_jesSubTotalPtDown;
  TBranch  *b_Jet_pt_jesSubTotalScaleDown;                                                        
  TBranch  *b_Jet_mass_jesSubTotalScaleDown;                                                      
  TBranch  *b_MET_pt_jesSubTotalScaleDown;
  TBranch  *b_MET_phi_jesSubTotalScaleDown;
  TBranch  *b_Jet_pt_jesSubTotalAbsoluteDown;                                                     
  TBranch  *b_Jet_mass_jesSubTotalAbsoluteDown;                                                   
  TBranch  *b_MET_pt_jesSubTotalAbsoluteDown;
  TBranch  *b_MET_phi_jesSubTotalAbsoluteDown;
  TBranch  *b_Jet_pt_jesSubTotalMCDown;                                                           
  TBranch  *b_Jet_mass_jesSubTotalMCDown;                                                         
  TBranch  *b_MET_pt_jesSubTotalMCDown;
  TBranch  *b_MET_phi_jesSubTotalMCDown;
  TBranch  *b_Jet_pt_jesTotalDown;                                                                
  TBranch  *b_Jet_mass_jesTotalDown;                                                              
  TBranch  *b_MET_pt_jesTotalDown;
  TBranch  *b_MET_phi_jesTotalDown;
  TBranch  *b_Jet_pt_jesTotalNoFlavorDown;                                                        
  TBranch  *b_Jet_mass_jesTotalNoFlavorDown;                                                      
  TBranch  *b_MET_pt_jesTotalNoFlavorDown;
  TBranch  *b_MET_phi_jesTotalNoFlavorDown;
  TBranch  *b_Jet_pt_jesTotalNoTimeDown;                                                          
  TBranch  *b_Jet_mass_jesTotalNoTimeDown;                                                        
  TBranch  *b_MET_pt_jesTotalNoTimeDown;
  TBranch  *b_MET_phi_jesTotalNoTimeDown;
  TBranch  *b_Jet_pt_jesTotalNoFlavorNoTimeDown;                                                  
  TBranch  *b_Jet_mass_jesTotalNoFlavorNoTimeDown;                                                
  TBranch  *b_MET_pt_jesTotalNoFlavorNoTimeDown;
  TBranch  *b_MET_phi_jesTotalNoFlavorNoTimeDown;
  TBranch  *b_Jet_pt_jesFlavorZJetDown;                                                           
  TBranch  *b_Jet_mass_jesFlavorZJetDown;                                                         
  TBranch  *b_MET_pt_jesFlavorZJetDown;
  TBranch  *b_MET_phi_jesFlavorZJetDown;
  TBranch  *b_Jet_pt_jesFlavorPhotonJetDown;                                                      
  TBranch  *b_Jet_mass_jesFlavorPhotonJetDown;                                                    
  TBranch  *b_MET_pt_jesFlavorPhotonJetDown;
  TBranch  *b_MET_phi_jesFlavorPhotonJetDown;
  TBranch  *b_Jet_pt_jesFlavorPureGluonDown;                                                      
  TBranch  *b_Jet_mass_jesFlavorPureGluonDown;                                                    
  TBranch  *b_MET_pt_jesFlavorPureGluonDown;
  TBranch  *b_MET_phi_jesFlavorPureGluonDown;
  TBranch  *b_Jet_pt_jesFlavorPureQuarkDown;                                                      
  TBranch  *b_Jet_mass_jesFlavorPureQuarkDown;                                                    
  TBranch  *b_MET_pt_jesFlavorPureQuarkDown;
  TBranch  *b_MET_phi_jesFlavorPureQuarkDown;
  TBranch  *b_Jet_pt_jesFlavorPureCharmDown;                                                      
  TBranch  *b_Jet_mass_jesFlavorPureCharmDown;                                                    
  TBranch  *b_MET_pt_jesFlavorPureCharmDown;
  TBranch  *b_MET_phi_jesFlavorPureCharmDown;
  TBranch  *b_Jet_pt_jesFlavorPureBottomDown;                                                     
  TBranch  *b_Jet_mass_jesFlavorPureBottomDown;                                                   
  TBranch  *b_MET_pt_jesFlavorPureBottomDown;
  TBranch  *b_MET_phi_jesFlavorPureBottomDown;
  TBranch  *b_Jet_pt_jesTimeRunADown;                                                             
  TBranch  *b_Jet_mass_jesTimeRunADown;                                                           
  TBranch  *b_MET_pt_jesTimeRunADown;
  TBranch  *b_MET_phi_jesTimeRunADown;
  TBranch  *b_Jet_pt_jesTimeRunBDown;                                                             
  TBranch  *b_Jet_mass_jesTimeRunBDown;                                                           
  TBranch  *b_MET_pt_jesTimeRunBDown;
  TBranch  *b_MET_phi_jesTimeRunBDown;
  TBranch  *b_Jet_pt_jesTimeRunCDown;                                                             
  TBranch  *b_Jet_mass_jesTimeRunCDown;                                                           
  TBranch  *b_MET_pt_jesTimeRunCDown;
  TBranch  *b_MET_phi_jesTimeRunCDown;
  TBranch  *b_Jet_pt_jesTimeRunDDown;                                                             
  TBranch  *b_Jet_mass_jesTimeRunDDown;                                                           
  TBranch  *b_MET_pt_jesTimeRunDDown;
  TBranch  *b_MET_phi_jesTimeRunDDown;
  TBranch  *b_Jet_pt_jesCorrelationGroupMPFInSituDown;                                            
  TBranch  *b_Jet_mass_jesCorrelationGroupMPFInSituDown;                                          
  TBranch  *b_MET_pt_jesCorrelationGroupMPFInSituDown;
  TBranch  *b_MET_phi_jesCorrelationGroupMPFInSituDown;
  TBranch  *b_Jet_pt_jesCorrelationGroupIntercalibrationDown;                                     
  TBranch  *b_Jet_mass_jesCorrelationGroupIntercalibrationDown;                                   
  TBranch  *b_MET_pt_jesCorrelationGroupIntercalibrationDown;
  TBranch  *b_MET_phi_jesCorrelationGroupIntercalibrationDown;
  TBranch  *b_Jet_pt_jesCorrelationGroupbJESDown;                                                 
  TBranch  *b_Jet_mass_jesCorrelationGroupbJESDown;                                               
  TBranch  *b_MET_pt_jesCorrelationGroupbJESDown;
  TBranch  *b_MET_phi_jesCorrelationGroupbJESDown;
  TBranch  *b_Jet_pt_jesCorrelationGroupFlavorDown;                                               
  TBranch  *b_Jet_mass_jesCorrelationGroupFlavorDown;                                             
  TBranch  *b_MET_pt_jesCorrelationGroupFlavorDown;
  TBranch  *b_MET_phi_jesCorrelationGroupFlavorDown;
  TBranch  *b_Jet_pt_jesCorrelationGroupUncorrelatedDown;                                         
  TBranch  *b_Jet_mass_jesCorrelationGroupUncorrelatedDown;                                       
  TBranch  *b_MET_pt_jesCorrelationGroupUncorrelatedDown;
  TBranch  *b_MET_phi_jesCorrelationGroupUncorrelatedDown;
  TBranch  *b_MET_pt_unclustEnDown;
  TBranch  *b_MET_phi_unclustEnDown;
  TBranch  *b_FatJet_pt_jerUp;                                                                  
  TBranch  *b_FatJet_mass_jerUp;                                                                
  TBranch  *b_FatJet_mass_jmrUp;                                                                
  TBranch  *b_FatJet_mass_jmsUp;                                                                
  TBranch  *b_FatJet_msoftdrop_jerUp;                                                           
  TBranch  *b_FatJet_msoftdrop_jmrUp;                                                           
  TBranch  *b_FatJet_msoftdrop_jmsUp;                                                           
  TBranch  *b_FatJet_msoftdrop_tau21DDT_jerUp;                                                  
  TBranch  *b_FatJet_msoftdrop_tau21DDT_jmrUp;                                                  
  TBranch  *b_FatJet_msoftdrop_tau21DDT_jmsUp;                                                  
  TBranch  *b_FatJet_pt_jesAbsoluteStatUp;                                                      
  TBranch  *b_FatJet_mass_jesAbsoluteStatUp;                                                    
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteStatUp;                                               
  TBranch  *b_FatJet_pt_jesAbsoluteScaleUp;                                                     
  TBranch  *b_FatJet_mass_jesAbsoluteScaleUp;                                                   
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteScaleUp;                                              
  TBranch  *b_FatJet_pt_jesAbsoluteSampleUp;                                                    
  TBranch  *b_FatJet_mass_jesAbsoluteSampleUp;                                                  
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteSampleUp;                                             
  TBranch  *b_FatJet_pt_jesAbsoluteFlavMapUp;                                                   
  TBranch  *b_FatJet_mass_jesAbsoluteFlavMapUp;                                                 
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteFlavMapUp;                                            
  TBranch  *b_FatJet_pt_jesAbsoluteMPFBiasUp;                                                   
  TBranch  *b_FatJet_mass_jesAbsoluteMPFBiasUp;                                                 
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteMPFBiasUp;                                            
  TBranch  *b_FatJet_pt_jesFragmentationUp;                                                     
  TBranch  *b_FatJet_mass_jesFragmentationUp;                                                   
  TBranch  *b_FatJet_msoftdrop_jesFragmentationUp;                                              
  TBranch  *b_FatJet_pt_jesSinglePionECALUp;                                                    
  TBranch  *b_FatJet_mass_jesSinglePionECALUp;
  TBranch  *b_FatJet_msoftdrop_jesSinglePionECALUp;                                             
  TBranch  *b_FatJet_pt_jesSinglePionHCALUp;                                                    
  TBranch  *b_FatJet_mass_jesSinglePionHCALUp;                                                  
  TBranch  *b_FatJet_msoftdrop_jesSinglePionHCALUp;                                             
  TBranch  *b_FatJet_pt_jesFlavorQCDUp;                                                         
  TBranch  *b_FatJet_mass_jesFlavorQCDUp;                                                       
  TBranch  *b_FatJet_msoftdrop_jesFlavorQCDUp;                                                  
  TBranch  *b_FatJet_pt_jesTimePtEtaUp;                                                         
  TBranch  *b_FatJet_mass_jesTimePtEtaUp;                                                       
  TBranch  *b_FatJet_msoftdrop_jesTimePtEtaUp;                                                  
  TBranch  *b_FatJet_pt_jesRelativeJEREC1Up;                                                    
  TBranch  *b_FatJet_mass_jesRelativeJEREC1Up;                                                  
  TBranch  *b_FatJet_msoftdrop_jesRelativeJEREC1Up;                                             
  TBranch  *b_FatJet_pt_jesRelativeJEREC2Up;                                                    
  TBranch  *b_FatJet_mass_jesRelativeJEREC2Up;                                                  
  TBranch  *b_FatJet_msoftdrop_jesRelativeJEREC2Up;                                             
  TBranch  *b_FatJet_pt_jesRelativeJERHFUp;                                                     
  TBranch  *b_FatJet_mass_jesRelativeJERHFUp;                                                   
  TBranch  *b_FatJet_msoftdrop_jesRelativeJERHFUp;                                              
  TBranch  *b_FatJet_pt_jesRelativePtBBUp;                                                      
  TBranch  *b_FatJet_mass_jesRelativePtBBUp;                                                    
  TBranch  *b_FatJet_msoftdrop_jesRelativePtBBUp;                                               
  TBranch  *b_FatJet_pt_jesRelativePtEC1Up;                                                     
  TBranch  *b_FatJet_mass_jesRelativePtEC1Up;                                                   
  TBranch  *b_FatJet_msoftdrop_jesRelativePtEC1Up;                                              
  TBranch  *b_FatJet_pt_jesRelativePtEC2Up;                                                     
  TBranch  *b_FatJet_mass_jesRelativePtEC2Up;                                                   
  TBranch  *b_FatJet_msoftdrop_jesRelativePtEC2Up;                                              
  TBranch  *b_FatJet_pt_jesRelativePtHFUp;                                                      
  TBranch  *b_FatJet_mass_jesRelativePtHFUp;                                                    
  TBranch  *b_FatJet_msoftdrop_jesRelativePtHFUp;                                               
  TBranch  *b_FatJet_pt_jesRelativeBalUp;                                                       
  TBranch  *b_FatJet_mass_jesRelativeBalUp;                                                     
  TBranch  *b_FatJet_msoftdrop_jesRelativeBalUp; 
  TBranch  *b_FatJet_pt_jesRelativeSampleUp;                                                    
  TBranch  *b_FatJet_mass_jesRelativeSampleUp;                                                  
  TBranch  *b_FatJet_msoftdrop_jesRelativeSampleUp;                                             
  TBranch  *b_FatJet_pt_jesRelativeFSRUp;                                                       
  TBranch  *b_FatJet_mass_jesRelativeFSRUp;                                                     
  TBranch  *b_FatJet_msoftdrop_jesRelativeFSRUp;                                                
  TBranch  *b_FatJet_pt_jesRelativeStatFSRUp;                                                   
  TBranch  *b_FatJet_mass_jesRelativeStatFSRUp;                                                 
  TBranch  *b_FatJet_msoftdrop_jesRelativeStatFSRUp;                                            
  TBranch  *b_FatJet_pt_jesRelativeStatECUp;                                                    
  TBranch  *b_FatJet_mass_jesRelativeStatECUp;                                                  
  TBranch  *b_FatJet_msoftdrop_jesRelativeStatECUp;                                             
  TBranch  *b_FatJet_pt_jesRelativeStatHFUp;                                                    
  TBranch  *b_FatJet_mass_jesRelativeStatHFUp;                                                  
  TBranch  *b_FatJet_msoftdrop_jesRelativeStatHFUp;                                             
  TBranch  *b_FatJet_pt_jesPileUpDataMCUp;                                                      
  TBranch  *b_FatJet_mass_jesPileUpDataMCUp;                                                    
  TBranch  *b_FatJet_msoftdrop_jesPileUpDataMCUp;                                               
  TBranch  *b_FatJet_pt_jesPileUpPtRefUp;                                                       
  TBranch  *b_FatJet_mass_jesPileUpPtRefUp;                                                     
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtRefUp;                                                
  TBranch  *b_FatJet_pt_jesPileUpPtBBUp;                                                        
  TBranch  *b_FatJet_mass_jesPileUpPtBBUp;                                                      
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtBBUp;                                                 
  TBranch  *b_FatJet_pt_jesPileUpPtEC1Up;                                                       
  TBranch  *b_FatJet_mass_jesPileUpPtEC1Up;                                                     
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtEC1Up;                                                
  TBranch  *b_FatJet_pt_jesPileUpPtEC2Up;                                                       
  TBranch  *b_FatJet_mass_jesPileUpPtEC2Up;                                                     
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtEC2Up;                                                
  TBranch  *b_FatJet_pt_jesPileUpPtHFUp;                                                        
  TBranch  *b_FatJet_mass_jesPileUpPtHFUp;                                                      
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtHFUp;                                                 
  TBranch  *b_FatJet_pt_jesPileUpMuZeroUp;
  TBranch  *b_FatJet_mass_jesPileUpMuZeroUp;                                                    
  TBranch  *b_FatJet_msoftdrop_jesPileUpMuZeroUp;                                               
  TBranch  *b_FatJet_pt_jesPileUpEnvelopeUp;                                                    
  TBranch  *b_FatJet_mass_jesPileUpEnvelopeUp;                                                  
  TBranch  *b_FatJet_msoftdrop_jesPileUpEnvelopeUp;                                             
  TBranch  *b_FatJet_pt_jesSubTotalPileUpUp;                                                    
  TBranch  *b_FatJet_mass_jesSubTotalPileUpUp;                                                  
  TBranch  *b_FatJet_msoftdrop_jesSubTotalPileUpUp;                                             
  TBranch  *b_FatJet_pt_jesSubTotalRelativeUp;                                                  
  TBranch  *b_FatJet_mass_jesSubTotalRelativeUp;                                                
  TBranch  *b_FatJet_msoftdrop_jesSubTotalRelativeUp;                                           
  TBranch  *b_FatJet_pt_jesSubTotalPtUp;                                                        
  TBranch  *b_FatJet_mass_jesSubTotalPtUp;                                                      
  TBranch  *b_FatJet_msoftdrop_jesSubTotalPtUp;                                                 
  TBranch  *b_FatJet_pt_jesSubTotalScaleUp;                                                     
  TBranch  *b_FatJet_mass_jesSubTotalScaleUp;                                                   
  TBranch  *b_FatJet_msoftdrop_jesSubTotalScaleUp;                                              
  TBranch  *b_FatJet_pt_jesSubTotalAbsoluteUp;                                                  
  TBranch  *b_FatJet_mass_jesSubTotalAbsoluteUp;                                                
  TBranch  *b_FatJet_msoftdrop_jesSubTotalAbsoluteUp;                                           
  TBranch  *b_FatJet_pt_jesSubTotalMCUp;                                                        
  TBranch  *b_FatJet_mass_jesSubTotalMCUp;                                                      
  TBranch  *b_FatJet_msoftdrop_jesSubTotalMCUp;                                                 
  TBranch  *b_FatJet_pt_jesTotalUp;                                                             
  TBranch  *b_FatJet_mass_jesTotalUp;                                                           
  TBranch  *b_FatJet_msoftdrop_jesTotalUp;                                                      
  TBranch  *b_FatJet_pt_jesTotalNoFlavorUp;                                                     
  TBranch  *b_FatJet_mass_jesTotalNoFlavorUp;
  TBranch  *b_FatJet_msoftdrop_jesTotalNoFlavorUp;                                              
  TBranch  *b_FatJet_pt_jesTotalNoTimeUp;                                                       
  TBranch  *b_FatJet_mass_jesTotalNoTimeUp;                                                     
  TBranch  *b_FatJet_msoftdrop_jesTotalNoTimeUp;                                                
  TBranch  *b_FatJet_pt_jesTotalNoFlavorNoTimeUp;                                               
  TBranch  *b_FatJet_mass_jesTotalNoFlavorNoTimeUp;                                             
  TBranch  *b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeUp;                                        
  TBranch  *b_FatJet_pt_jesFlavorZJetUp;                                                        
  TBranch  *b_FatJet_mass_jesFlavorZJetUp;                                                      
  TBranch  *b_FatJet_msoftdrop_jesFlavorZJetUp;                                                 
  TBranch  *b_FatJet_pt_jesFlavorPhotonJetUp;                                                   
  TBranch  *b_FatJet_mass_jesFlavorPhotonJetUp;                                                 
  TBranch  *b_FatJet_msoftdrop_jesFlavorPhotonJetUp;                                            
  TBranch  *b_FatJet_pt_jesFlavorPureGluonUp;                                                   
  TBranch  *b_FatJet_mass_jesFlavorPureGluonUp;                                                 
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureGluonUp;                                            
  TBranch  *b_FatJet_pt_jesFlavorPureQuarkUp;                                                   
  TBranch  *b_FatJet_mass_jesFlavorPureQuarkUp;                                                 
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureQuarkUp;                                            
  TBranch  *b_FatJet_pt_jesFlavorPureCharmUp;                                                   
  TBranch  *b_FatJet_mass_jesFlavorPureCharmUp;                                                 
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureCharmUp;                                            
  TBranch  *b_FatJet_pt_jesFlavorPureBottomUp;                                                  
  TBranch  *b_FatJet_mass_jesFlavorPureBottomUp;                                                
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureBottomUp;                                           
  TBranch  *b_FatJet_pt_jesTimeRunAUp;                                                          
  TBranch  *b_FatJet_mass_jesTimeRunAUp;                                                        
  TBranch  *b_FatJet_msoftdrop_jesTimeRunAUp;                                                   
  TBranch  *b_FatJet_pt_jesTimeRunBUp;                                                          
  TBranch  *b_FatJet_mass_jesTimeRunBUp;                                                        
  TBranch  *b_FatJet_msoftdrop_jesTimeRunBUp;                                                   
  TBranch  *b_FatJet_pt_jesTimeRunCUp;                                                          
  TBranch  *b_FatJet_mass_jesTimeRunCUp;                                                        
  TBranch  *b_FatJet_msoftdrop_jesTimeRunCUp;                                                   
  TBranch  *b_FatJet_pt_jesTimeRunDUp;                                                          
  TBranch  *b_FatJet_mass_jesTimeRunDUp;                                                        
  TBranch  *b_FatJet_msoftdrop_jesTimeRunDUp;                                                   
  TBranch  *b_FatJet_pt_jesCorrelationGroupMPFInSituUp;                                         
  TBranch  *b_FatJet_mass_jesCorrelationGroupMPFInSituUp;
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituUp;                                  
  TBranch  *b_FatJet_pt_jesCorrelationGroupIntercalibrationUp;                                  
  TBranch  *b_FatJet_mass_jesCorrelationGroupIntercalibrationUp;                                
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationUp;                           
  TBranch  *b_FatJet_pt_jesCorrelationGroupbJESUp;                                              
  TBranch  *b_FatJet_mass_jesCorrelationGroupbJESUp;                                            
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupbJESUp;                                       
  TBranch  *b_FatJet_pt_jesCorrelationGroupFlavorUp;                                            
  TBranch  *b_FatJet_mass_jesCorrelationGroupFlavorUp;                                          
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupFlavorUp;                                     
  TBranch  *b_FatJet_pt_jesCorrelationGroupUncorrelatedUp;                                      
  TBranch  *b_FatJet_mass_jesCorrelationGroupUncorrelatedUp;                                    
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedUp;                               
  TBranch  *b_FatJet_pt_jerDown;                                                                
  TBranch  *b_FatJet_mass_jerDown;                                                              
  TBranch  *b_FatJet_mass_jmrDown;                                                              
  TBranch  *b_FatJet_mass_jmsDown;                                                              
  TBranch  *b_FatJet_msoftdrop_jerDown;                                                         
  TBranch  *b_FatJet_msoftdrop_jmrDown;                                                         
  TBranch  *b_FatJet_msoftdrop_jmsDown;                                                         
  TBranch  *b_FatJet_msoftdrop_tau21DDT_jerDown;                                                
  TBranch  *b_FatJet_msoftdrop_tau21DDT_jmrDown;                                                
  TBranch  *b_FatJet_msoftdrop_tau21DDT_jmsDown;                                                
  TBranch  *b_FatJet_pt_jesAbsoluteStatDown;                                                    
  TBranch  *b_FatJet_mass_jesAbsoluteStatDown;                                                  
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteStatDown;                                             
  TBranch  *b_FatJet_pt_jesAbsoluteScaleDown;                                                   
  TBranch  *b_FatJet_mass_jesAbsoluteScaleDown;                                                 
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteScaleDown;                                            
  TBranch  *b_FatJet_pt_jesAbsoluteSampleDown;                                                  
  TBranch  *b_FatJet_mass_jesAbsoluteSampleDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteSampleDown;                                           
  TBranch  *b_FatJet_pt_jesAbsoluteFlavMapDown;                                                 
  TBranch  *b_FatJet_mass_jesAbsoluteFlavMapDown;                                               
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteFlavMapDown;                                          
  TBranch  *b_FatJet_pt_jesAbsoluteMPFBiasDown;
  TBranch  *b_FatJet_mass_jesAbsoluteMPFBiasDown;                                               
  TBranch  *b_FatJet_msoftdrop_jesAbsoluteMPFBiasDown;                                          
  TBranch  *b_FatJet_pt_jesFragmentationDown;                                                   
  TBranch  *b_FatJet_mass_jesFragmentationDown;                                                 
  TBranch  *b_FatJet_msoftdrop_jesFragmentationDown;                                            
  TBranch  *b_FatJet_pt_jesSinglePionECALDown;                                                  
  TBranch  *b_FatJet_mass_jesSinglePionECALDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesSinglePionECALDown;                                           
  TBranch  *b_FatJet_pt_jesSinglePionHCALDown;                                                  
  TBranch  *b_FatJet_mass_jesSinglePionHCALDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesSinglePionHCALDown;                                           
  TBranch  *b_FatJet_pt_jesFlavorQCDDown;                                                       
  TBranch  *b_FatJet_mass_jesFlavorQCDDown;                                                     
  TBranch  *b_FatJet_msoftdrop_jesFlavorQCDDown;                                                
  TBranch  *b_FatJet_pt_jesTimePtEtaDown;                                                       
  TBranch  *b_FatJet_mass_jesTimePtEtaDown;                                                     
  TBranch  *b_FatJet_msoftdrop_jesTimePtEtaDown;                                                
  TBranch  *b_FatJet_pt_jesRelativeJEREC1Down;                                                  
  TBranch  *b_FatJet_mass_jesRelativeJEREC1Down;                                                
  TBranch  *b_FatJet_msoftdrop_jesRelativeJEREC1Down;                                           
  TBranch  *b_FatJet_pt_jesRelativeJEREC2Down;                                                  
  TBranch  *b_FatJet_mass_jesRelativeJEREC2Down;                                                
  TBranch  *b_FatJet_msoftdrop_jesRelativeJEREC2Down;                                           
  TBranch  *b_FatJet_pt_jesRelativeJERHFDown;                                                   
  TBranch  *b_FatJet_mass_jesRelativeJERHFDown;                                                 
  TBranch  *b_FatJet_msoftdrop_jesRelativeJERHFDown;                                            
  TBranch  *b_FatJet_pt_jesRelativePtBBDown;                                                    
  TBranch  *b_FatJet_mass_jesRelativePtBBDown;                                                  
  TBranch  *b_FatJet_msoftdrop_jesRelativePtBBDown;                                             
  TBranch  *b_FatJet_pt_jesRelativePtEC1Down;                                                   
  TBranch  *b_FatJet_mass_jesRelativePtEC1Down;                                                 
  TBranch  *b_FatJet_msoftdrop_jesRelativePtEC1Down;                                            
  TBranch  *b_FatJet_pt_jesRelativePtEC2Down;                                                   
  TBranch  *b_FatJet_mass_jesRelativePtEC2Down;
  TBranch  *b_FatJet_msoftdrop_jesRelativePtEC2Down;                                            
  TBranch  *b_FatJet_pt_jesRelativePtHFDown;                                                    
  TBranch  *b_FatJet_mass_jesRelativePtHFDown;                                                  
  TBranch  *b_FatJet_msoftdrop_jesRelativePtHFDown;                                             
  TBranch  *b_FatJet_pt_jesRelativeBalDown;                                                     
  TBranch  *b_FatJet_mass_jesRelativeBalDown;                                                   
  TBranch  *b_FatJet_msoftdrop_jesRelativeBalDown;                                              
  TBranch  *b_FatJet_pt_jesRelativeSampleDown;                                                  
  TBranch  *b_FatJet_mass_jesRelativeSampleDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesRelativeSampleDown;                                           
  TBranch  *b_FatJet_pt_jesRelativeFSRDown;                                                     
  TBranch  *b_FatJet_mass_jesRelativeFSRDown;                                                   
  TBranch  *b_FatJet_msoftdrop_jesRelativeFSRDown;                                              
  TBranch  *b_FatJet_pt_jesRelativeStatFSRDown;                                                 
  TBranch  *b_FatJet_mass_jesRelativeStatFSRDown;                                               
  TBranch  *b_FatJet_msoftdrop_jesRelativeStatFSRDown;                                          
  TBranch  *b_FatJet_pt_jesRelativeStatECDown;                                                  
  TBranch  *b_FatJet_mass_jesRelativeStatECDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesRelativeStatECDown;                                           
  TBranch  *b_FatJet_pt_jesRelativeStatHFDown;                                                  
  TBranch  *b_FatJet_mass_jesRelativeStatHFDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesRelativeStatHFDown;                                           
  TBranch  *b_FatJet_pt_jesPileUpDataMCDown;                                                    
  TBranch  *b_FatJet_mass_jesPileUpDataMCDown;                                                  
  TBranch  *b_FatJet_msoftdrop_jesPileUpDataMCDown;                                             
  TBranch  *b_FatJet_pt_jesPileUpPtRefDown;                                                     
  TBranch  *b_FatJet_mass_jesPileUpPtRefDown;                                                   
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtRefDown;                                              
  TBranch  *b_FatJet_pt_jesPileUpPtBBDown;                                                      
  TBranch  *b_FatJet_mass_jesPileUpPtBBDown;                                                    
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtBBDown;                                               
  TBranch  *b_FatJet_pt_jesPileUpPtEC1Down;
  TBranch  *b_FatJet_mass_jesPileUpPtEC1Down;                                                   
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtEC1Down;                                              
  TBranch  *b_FatJet_pt_jesPileUpPtEC2Down;                                                     
  TBranch  *b_FatJet_mass_jesPileUpPtEC2Down;                                                   
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtEC2Down;                                              
  TBranch  *b_FatJet_pt_jesPileUpPtHFDown;                                                      
  TBranch  *b_FatJet_mass_jesPileUpPtHFDown;                                                    
  TBranch  *b_FatJet_msoftdrop_jesPileUpPtHFDown;                                               
  TBranch  *b_FatJet_pt_jesPileUpMuZeroDown;                                                    
  TBranch  *b_FatJet_mass_jesPileUpMuZeroDown;                                                  
  TBranch  *b_FatJet_msoftdrop_jesPileUpMuZeroDown;                                             
  TBranch  *b_FatJet_pt_jesPileUpEnvelopeDown;                                                  
  TBranch  *b_FatJet_mass_jesPileUpEnvelopeDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesPileUpEnvelopeDown;                                           
  TBranch  *b_FatJet_pt_jesSubTotalPileUpDown;                                                  
  TBranch  *b_FatJet_mass_jesSubTotalPileUpDown;                                                
  TBranch  *b_FatJet_msoftdrop_jesSubTotalPileUpDown;                                           
  TBranch  *b_FatJet_pt_jesSubTotalRelativeDown;                                                
  TBranch  *b_FatJet_mass_jesSubTotalRelativeDown;                                              
  TBranch  *b_FatJet_msoftdrop_jesSubTotalRelativeDown;                                         
  TBranch  *b_FatJet_pt_jesSubTotalPtDown;                                                      
  TBranch  *b_FatJet_mass_jesSubTotalPtDown;                                                    
  TBranch  *b_FatJet_msoftdrop_jesSubTotalPtDown;                                               
  TBranch  *b_FatJet_pt_jesSubTotalScaleDown;                                                   
  TBranch  *b_FatJet_mass_jesSubTotalScaleDown;                                                 
  TBranch  *b_FatJet_msoftdrop_jesSubTotalScaleDown;                                            
  TBranch  *b_FatJet_pt_jesSubTotalAbsoluteDown;                                                
  TBranch  *b_FatJet_mass_jesSubTotalAbsoluteDown;                                              
  TBranch  *b_FatJet_msoftdrop_jesSubTotalAbsoluteDown;                                         
  TBranch  *b_FatJet_pt_jesSubTotalMCDown;                                                      
  TBranch  *b_FatJet_mass_jesSubTotalMCDown;                                                    
  TBranch  *b_FatJet_msoftdrop_jesSubTotalMCDown;
  TBranch  *b_FatJet_pt_jesTotalDown;                                                           
  TBranch  *b_FatJet_mass_jesTotalDown;                                                         
  TBranch  *b_FatJet_msoftdrop_jesTotalDown;                                                    
  TBranch  *b_FatJet_pt_jesTotalNoFlavorDown;                                                   
  TBranch  *b_FatJet_mass_jesTotalNoFlavorDown;                                                 
  TBranch  *b_FatJet_msoftdrop_jesTotalNoFlavorDown;                                            
  TBranch  *b_FatJet_pt_jesTotalNoTimeDown;                                                     
  TBranch  *b_FatJet_mass_jesTotalNoTimeDown;                                                   
  TBranch  *b_FatJet_msoftdrop_jesTotalNoTimeDown;                                              
  TBranch  *b_FatJet_pt_jesTotalNoFlavorNoTimeDown;                                             
  TBranch  *b_FatJet_mass_jesTotalNoFlavorNoTimeDown;                                           
  TBranch  *b_FatJet_msoftdrop_jesTotalNoFlavorNoTimeDown;                                      
  TBranch  *b_FatJet_pt_jesFlavorZJetDown;                                                      
  TBranch  *b_FatJet_mass_jesFlavorZJetDown;                                                    
  TBranch  *b_FatJet_msoftdrop_jesFlavorZJetDown;                                               
  TBranch  *b_FatJet_pt_jesFlavorPhotonJetDown;                                                 
  TBranch  *b_FatJet_mass_jesFlavorPhotonJetDown;                                               
  TBranch  *b_FatJet_msoftdrop_jesFlavorPhotonJetDown;                                          
  TBranch  *b_FatJet_pt_jesFlavorPureGluonDown;                                                 
  TBranch  *b_FatJet_mass_jesFlavorPureGluonDown;                                               
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureGluonDown;                                          
  TBranch  *b_FatJet_pt_jesFlavorPureQuarkDown;                                                 
  TBranch  *b_FatJet_mass_jesFlavorPureQuarkDown;                                               
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureQuarkDown;                                          
  TBranch  *b_FatJet_pt_jesFlavorPureCharmDown;                                                 
  TBranch  *b_FatJet_mass_jesFlavorPureCharmDown;                                               
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureCharmDown;                                          
  TBranch  *b_FatJet_pt_jesFlavorPureBottomDown;
  TBranch  *b_FatJet_mass_jesFlavorPureBottomDown;                                              
  TBranch  *b_FatJet_msoftdrop_jesFlavorPureBottomDown;                                         
  TBranch  *b_FatJet_pt_jesTimeRunADown;                                                        
  TBranch  *b_FatJet_mass_jesTimeRunADown;                                                      
  TBranch  *b_FatJet_msoftdrop_jesTimeRunADown;                                                 
  TBranch  *b_FatJet_pt_jesTimeRunBDown;                                                        
  TBranch  *b_FatJet_mass_jesTimeRunBDown;                                                      
  TBranch  *b_FatJet_msoftdrop_jesTimeRunBDown;                                                 
  TBranch  *b_FatJet_pt_jesTimeRunCDown;                                                        
  TBranch  *b_FatJet_mass_jesTimeRunCDown;                                                      
  TBranch  *b_FatJet_msoftdrop_jesTimeRunCDown;                                                 
  TBranch  *b_FatJet_pt_jesTimeRunDDown;                                                        
  TBranch  *b_FatJet_mass_jesTimeRunDDown;                                                      
  TBranch  *b_FatJet_msoftdrop_jesTimeRunDDown;                                                 
  TBranch  *b_FatJet_pt_jesCorrelationGroupMPFInSituDown;                                       
  TBranch  *b_FatJet_mass_jesCorrelationGroupMPFInSituDown;                                     
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupMPFInSituDown;                                
  TBranch  *b_FatJet_pt_jesCorrelationGroupIntercalibrationDown;                                
  TBranch  *b_FatJet_mass_jesCorrelationGroupIntercalibrationDown;                              
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupIntercalibrationDown;                         
  TBranch  *b_FatJet_pt_jesCorrelationGroupbJESDown;                                            
  TBranch  *b_FatJet_mass_jesCorrelationGroupbJESDown;                                          
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupbJESDown;                                     
  TBranch  *b_FatJet_pt_jesCorrelationGroupFlavorDown;                                          
  TBranch  *b_FatJet_mass_jesCorrelationGroupFlavorDown;                                        
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupFlavorDown;                                   
  TBranch  *b_FatJet_pt_jesCorrelationGroupUncorrelatedDown;                                    
  TBranch  *b_FatJet_mass_jesCorrelationGroupUncorrelatedDown;                                  
  TBranch  *b_FatJet_msoftdrop_jesCorrelationGroupUncorrelatedDown;
  
};

#endif
