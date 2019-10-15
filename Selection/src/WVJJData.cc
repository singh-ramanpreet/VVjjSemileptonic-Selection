#include "WVJJAna/Selection/interface/WVJJData.hh"

void WVJJData::init(TTree *t) {

  clearVars();

  fTree = t;
  
  //metadata and weights
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("ls",&ls,"ls/I");
  fTree->Branch("evt",&evt,"evt/I");
  fTree->Branch("genWeight",&genWeight,"genWeight/f");
  fTree->Branch("triggerEffWeight",&triggerEffWeight,"triggerEffWeight/f");
  fTree->Branch("L1PFWeight",&L1PFWeight,"L1PFWeight/f");
  //lepton 1
  fTree->Branch("lep1_pt",&lep1_pt,"lep1_pt/f");
  fTree->Branch("lep1_eta",&lep1_eta,"lep1_eta/f");
  fTree->Branch("lep1_phi",&lep1_phi,"lep1_phi/f");
  fTree->Branch("lep1_m",&lep1_m,"lep1_m/f");
  fTree->Branch("lep1_q",&lep1_q,"lep1_q/f");
  fTree->Branch("lep1_iso",&lep1_iso,"lep1_iso/f");
  fTree->Branch("lep1_idEffWeight",&lep1_idEffWeight,"lep1_idEffWeight/f");
  //lepton 1 scale variations
  fTree->Branch("lep1_pt_scaleUp",&lep1_pt_scaleUp,"lep1_pt_scaleUp/f");
  fTree->Branch("lep1_pt_scaleDn",&lep1_pt_scaleDn,"lep1_pt_scaleDn/f");
  //lepton 2
  fTree->Branch("lep2_pt",&lep2_pt,"lep2_pt/f");
  fTree->Branch("lep2_eta",&lep2_eta,"lep2_eta/f");
  fTree->Branch("lep2_phi",&lep2_phi,"lep2_phi/f");
  fTree->Branch("lep2_m",&lep2_m,"lep2_m/f");
  fTree->Branch("lep2_q",&lep2_q,"lep2_q/f");
  fTree->Branch("lep2_iso",&lep2_iso,"lep2_iso/f");
  fTree->Branch("lep2_idEffWeight",&lep2_idEffWeight,"lep2_idEffWeight/f");  
  //lepton 2 scale variations
  fTree->Branch("lep2_pt_scaleUp",&lep2_pt_scaleUp,"lep2_pt_scaleUp/f");
  fTree->Branch("lep2_pt_scaleDn",&lep2_pt_scaleDn,"lep2_pt_scaleDn/f");
  //dilepton final state
  fTree->Branch("dilep_m",&dilep_m,"dilep_m/f");
  fTree->Branch("dilep_pt",&dilep_pt,"dilep_pt/f");
  fTree->Branch("dilep_eta",&dilep_eta,"dilep_eta/f");
  fTree->Branch("dilep_phi",&dilep_phi,"dilep_phi/f");
  //dilepton scale variations
  fTree->Branch("dilep_m_scaleUp",&dilep_m_scaleUp,"dilep_m_scaleUp/f");
  fTree->Branch("dilep_m_scaleDn",&dilep_m_scaleDn,"dilep_m_scaleDn/f");
  fTree->Branch("dilep_pt_scaleUp",&dilep_pt_scaleUp,"dilep_pt_scaleUp/f");
  fTree->Branch("dilep_pt_scaleDn",&dilep_pt_scaleDn,"dilep_pt_scaleDn/f");
  //VBF jet 1
  fTree->Branch("vbf1_AK4_pt",&vbf1_AK4_pt,"vbf1_AK4_pt/f");
  fTree->Branch("vbf1_AK4_eta",&vbf1_AK4_eta,"vbf1_AK4_eta/f");
  fTree->Branch("vbf1_AK4_phi",&vbf1_AK4_phi,"vbf1_AK4_phi/f");
  fTree->Branch("vbf1_AK4_m",&vbf1_AK4_m,"vbf1_AK4_m/f");
  fTree->Branch("vbf1_AK4_e",&vbf1_AK4_e,"vbf1_AK4_e/f");
  //VBF jet 1 variations
  fTree->Branch("vbf1_AK4_pt_scaleUp",&vbf1_AK4_pt_scaleUp,"vbf1_AK4_pt_scaleUp/f");
  fTree->Branch("vbf1_AK4_pt_scaleDn",&vbf1_AK4_pt_scaleDn,"vbf1_AK4_pt_scaleDn/f");
  fTree->Branch("vbf1_AK4_m_scaleUp",&vbf1_AK4_m_scaleUp,"vbf1_AK4_m_scaleUp/f");
  fTree->Branch("vbf1_AK4_m_scaleDn",&vbf1_AK4_m_scaleDn,"vbf1_AK4_m_scaleDn/f");
  fTree->Branch("vbf1_AK4_e_scaleUp",&vbf1_AK4_e_scaleUp,"vbf1_AK4_e_scaleUp/f");
  fTree->Branch("vbf1_AK4_e_scaleDn",&vbf1_AK4_e_scaleDn,"vbf1_AK4_e_scaleDn/f");
  //VBF jet 2
  fTree->Branch("vbf2_AK4_pt",&vbf2_AK4_pt,"vbf2_AK4_pt/f");
  fTree->Branch("vbf2_AK4_eta",&vbf2_AK4_eta,"vbf2_AK4_eta/f");
  fTree->Branch("vbf2_AK4_phi",&vbf2_AK4_phi,"vbf2_AK4_phi/f");
  fTree->Branch("vbf2_AK4_m",&vbf2_AK4_m,"vbf2_AK4_m/f");
  fTree->Branch("vbf2_AK4_e",&vbf2_AK4_e,"vbf2_AK4_e/f");
  //VBF jet 2 variations
  fTree->Branch("vbf2_AK4_pt_scaleUp",&vbf2_AK4_pt_scaleUp,"vbf2_AK4_pt_scaleUp/f");
  fTree->Branch("vbf2_AK4_pt_scaleDn",&vbf2_AK4_pt_scaleDn,"vbf2_AK4_pt_scaleDn/f");
  fTree->Branch("vbf2_AK4_m_scaleUp",&vbf2_AK4_m_scaleUp,"vbf2_AK4_m_scaleUp/f");
  fTree->Branch("vbf2_AK4_m_scaleDn",&vbf2_AK4_m_scaleDn,"vbf2_AK4_m_scaleDn/f");
  fTree->Branch("vbf2_AK4_e_scaleUp",&vbf2_AK4_e_scaleUp,"vbf2_AK4_e_scaleUp/f");
  fTree->Branch("vbf2_AK4_e_scaleDn",&vbf2_AK4_e_scaleDn,"vbf2_AK4_e_scaleDn/f");
  //VBF dijet object
  fTree->Branch("vbf_AK4AK4_pt", &vbf_AK4AK4_pt, "vbf_AK4AK4_pt/f");
  fTree->Branch("vbf_AK4AK4_eta",&vbf_AK4AK4_eta,"vbf_AK4AK4_eta/f");
  fTree->Branch("vbf_AK4AK4_phi",&vbf_AK4AK4_phi,"vbf_AK4AK4_phi/f");
  fTree->Branch("vbf_AK4AK4_m",  &vbf_AK4AK4_m,  "vbf_AK4AK4_m/f");
  //Boson AK8 jet
  fTree->Branch("bos_PuppiAK8_m_sd0",&bos_PuppiAK8_m_sd0,"bos_PuppiAK8_m_sd0/f");
  fTree->Branch("bos_PuppiAK8_m_sd0_corr",&bos_PuppiAK8_m_sd0_corr,"bos_PuppiAK8_m_sd0_corr/f");
  fTree->Branch("bos_PuppiAK8_pt_ungroomed",&bos_PuppiAK8_pt_ungroomed,"bos_PuppiAK8_pt_ungroomed/f");
  fTree->Branch("bos_PuppiAK8_eta_ungroomed",&bos_PuppiAK8_eta_ungroomed,"bos_PuppiAK8_eta_ungroomed/f");
  fTree->Branch("bos_PuppiAK8_phi_ungroomed",&bos_PuppiAK8_phi_ungroomed,"bos_PuppiAK8_phi_ungroomed/f");
  fTree->Branch("bos_PuppiAK8_e_ungroomed",&bos_PuppiAK8_e_ungroomed,"bos_PuppiAK8_e_ungroomed/f");
  fTree->Branch("bos_PuppiAK8_tau2tau1",&bos_PuppiAK8_tau2tau1,"bos_PuppiAK8_tau2tau1/f");
  //Boson AK8 jet variations
  fTree->Branch("bos_PuppiAK8_m_sd0_corr_scaleUp",&bos_PuppiAK8_m_sd0_corr_scaleUp,"bos_PuppiAK8_m_sd0_corr_scaleUp/f");
  fTree->Branch("bos_PuppiAK8_m_sd0_corr_scaleDn",&bos_PuppiAK8_m_sd0_corr_scaleDn,"bos_PuppiAK8_m_sd0_corr_scaleDn/f");
  fTree->Branch("bos_PuppiAK8_pt_ungroomed_scaleUp",&bos_PuppiAK8_pt_ungroomed_scaleUp,"bos_PuppiAK8_pt_ungroomed_scaleUp/f");
  fTree->Branch("bos_PuppiAK8_pt_ungroomed_scaleDn",&bos_PuppiAK8_pt_ungroomed_scaleDn,"bos_PuppiAK8_pt_ungroomed_scaleDn/f");
  fTree->Branch("bos_PuppiAK8_e_ungroomed_scaleUp",&bos_PuppiAK8_e_ungroomed_scaleUp,"bos_PuppiAK8_e_ungroomed_scaleUp/f");
  fTree->Branch("bos_PuppiAK8_e_ungroomed_scaleDn",&bos_PuppiAK8_e_ungroomed_scaleDn,"bos_PuppiAK8_e_ungroomed_scaleDn/f");
  //Boson AK4 jet 1
  fTree->Branch("bos_j1_AK4_pt",&bos_j1_AK4_pt,"bos_j1_AK4_pt/f");
  fTree->Branch("bos_j1_AK4_eta",&bos_j1_AK4_eta,"bos_j1_AK4_eta/f");
  fTree->Branch("bos_j1_AK4_phi",&bos_j1_AK4_phi,"bos_j1_AK4_phi/f");
  fTree->Branch("bos_j1_AK4_m",&bos_j1_AK4_m,"bos_j1_AK4_m/f");
  fTree->Branch("bos_j1_AK4_e",&bos_j1_AK4_e,"bos_j1_AK4_e/f");
  //Boson AK4 jet 1 variations
  fTree->Branch("bos_j1_AK4_pt_scaleUp",&bos_j1_AK4_pt_scaleUp,"bos_j1_AK4_pt_scaleUp/f");
  fTree->Branch("bos_j1_AK4_pt_scaleDn",&bos_j1_AK4_pt_scaleDn,"bos_j1_AK4_pt_scaleDn/f");
  fTree->Branch("bos_j1_AK4_m_scaleUp",&bos_j1_AK4_m_scaleUp,"bos_j1_AK4_m_scaleUp/f");
  fTree->Branch("bos_j1_AK4_m_scaleDn",&bos_j1_AK4_m_scaleDn,"bos_j1_AK4_m_scaleDn/f");
  fTree->Branch("bos_j1_AK4_e_scaleUp",&bos_j1_AK4_e_scaleUp,"bos_j1_AK4_e_scaleUp/f");
  fTree->Branch("bos_j1_AK4_e_scaleDn",&bos_j1_AK4_e_scaleDn,"bos_j1_AK4_e_scaleDn/f");
  //Boson AK4 jet 2
  fTree->Branch("bos_j2_AK4_pt",&bos_j2_AK4_pt,"bos_j2_AK4_pt/f");
  fTree->Branch("bos_j2_AK4_eta",&bos_j2_AK4_eta,"bos_j2_AK4_eta/f");
  fTree->Branch("bos_j2_AK4_phi",&bos_j2_AK4_phi,"bos_j2_AK4_phi/f");
  fTree->Branch("bos_j2_AK4_m",&bos_j2_AK4_m,"bos_j2_AK4_m/f");
  fTree->Branch("bos_j2_AK4_e",&bos_j2_AK4_e,"bos_j2_AK4_e/f");
  //Boson AK4 jet 2 variations
  fTree->Branch("bos_j2_AK4_pt_scaleUp",&bos_j2_AK4_pt_scaleUp,"bos_j2_AK4_pt_scaleUp/f");
  fTree->Branch("bos_j2_AK4_pt_scaleDn",&bos_j2_AK4_pt_scaleDn,"bos_j2_AK4_pt_scaleDn/f");
  fTree->Branch("bos_j2_AK4_m_scaleUp",&bos_j2_AK4_m_scaleUp,"bos_j2_AK4_m_scaleUp/f");
  fTree->Branch("bos_j2_AK4_m_scaleDn",&bos_j2_AK4_m_scaleDn,"bos_j2_AK4_m_scaleDn/f");
  fTree->Branch("bos_j2_AK4_e_scaleUp",&bos_j2_AK4_e_scaleUp,"bos_j2_AK4_e_scaleUp/f");
  fTree->Branch("bos_j2_AK4_e_scaleDn",&bos_j2_AK4_e_scaleDn,"bos_j2_AK4_e_scaleDn/f");
  //Boson dijet object
  fTree->Branch("bos_AK4AK4_pt", &bos_AK4AK4_pt, "bos_AK4AK4_pt/f");
  fTree->Branch("bos_AK4AK4_eta",&bos_AK4AK4_eta,"bos_AK4AK4_eta/f");
  fTree->Branch("bos_AK4AK4_phi",&bos_AK4AK4_phi,"bos_AK4AK4_phi/f");
  fTree->Branch("bos_AK4AK4_m",  &bos_AK4AK4_m,  "bos_AK4AK4_m/f");
  
};

void WVJJData::clearVars() {

  //------------------------------------//
  //       METADATA AND EVENT WEIGHTS   //
  //------------------------------------//
  
  run = 0;
  ls = 0;
  evt = 0;
  
  genWeight = 1.0;
  triggerEffWeight = 1.0;
  L1PFWeight = 1.0;
  
  //------------------------------------//
  //       LEPTONS                      //
  //------------------------------------//

  //lepton 1
  lep1_pt = -999.0;
  lep1_eta = -999.0;
  lep1_phi = -999.0;
  lep1_m = -999.0;
  lep1_q = -999.0;
  lep1_iso = -999.0;
  lep1_idEffWeight = 1.0;

  //lepton 1 scale variations
  lep1_pt_scaleUp = -999.0;
  lep1_pt_scaleDn = -999.0;
  
  //lepton 2
  lep2_pt = -999.0;
  lep2_eta = -999.0;
  lep2_phi = -999.0;
  lep2_m = -999.0;
  lep2_q = -999.0;
  lep2_iso = -999.0;
  lep2_idEffWeight = 1.0;

  //lepton 2 scale variations
  lep2_pt_scaleUp = -999.0;
  lep2_pt_scaleDn = -999.0;

  //dilepton final state
  dilep_m = -999.0;
  dilep_pt = -999.0;
  dilep_eta = -999.0;
  dilep_phi = -999.0;

  //dilepton scale variations
  dilep_m_scaleUp = -999.0;
  dilep_m_scaleDn = -999.0;
  dilep_pt_scaleUp = -999.0;
  dilep_pt_scaleDn = -999.0;

  //------------------------------------//
  //       VBF/TAGGING JETS             //
  //------------------------------------//

  //VBF jet 1
  vbf1_AK4_pt = -999.0;
  vbf1_AK4_eta = -999.0;
  vbf1_AK4_phi = -999.0;
  vbf1_AK4_m = -999.0;
  vbf1_AK4_e = -999.0;

  //VBF jet 1 variations
  vbf1_AK4_pt_scaleUp = -999.0;
  vbf1_AK4_pt_scaleDn = -999.0;
  vbf1_AK4_m_scaleUp = -999.0;
  vbf1_AK4_m_scaleDn = -999.0;
  vbf1_AK4_e_scaleUp = -999.0;
  vbf1_AK4_e_scaleDn = -999.0;

  //VBF jet 2
  vbf2_AK4_pt = -999.0;
  vbf2_AK4_eta = -999.0;
  vbf2_AK4_phi = -999.0;
  vbf2_AK4_m = -999.0;
  vbf2_AK4_e = -999.0;

  //VBF jet 2 variations
  vbf2_AK4_pt_scaleUp = -999.0;
  vbf2_AK4_pt_scaleDn = -999.0;
  vbf2_AK4_m_scaleUp = -999.0;
  vbf2_AK4_m_scaleDn = -999.0;
  vbf2_AK4_e_scaleUp = -999.0;
  vbf2_AK4_e_scaleDn = -999.0;

  //------------------------------------//
  //       HADRONIC BOOSTED OBJECTS     //
  //------------------------------------//

  //Boson AK8 jet
  bos_PuppiAK8_m_sd0 = -999.0;
  bos_PuppiAK8_m_sd0_corr = -999.0;
  bos_PuppiAK8_pt_ungroomed = -999.0;
  bos_PuppiAK8_eta_ungroomed = -999.0;
  bos_PuppiAK8_phi_ungroomed = -999.0;
  bos_PuppiAK8_e_ungroomed = -999.0;
  bos_PuppiAK8_tau2tau1 = -999.0;

  //Boson AK8 jet variations
  bos_PuppiAK8_m_sd0_corr_scaleUp = -999.0;
  bos_PuppiAK8_m_sd0_corr_scaleDn = -999.0;
  bos_PuppiAK8_pt_ungroomed_scaleUp = -999.0;
  bos_PuppiAK8_pt_ungroomed_scaleDn = -999.0;
  bos_PuppiAK8_e_ungroomed_scaleUp = -999.0;
  bos_PuppiAK8_e_ungroomed_scaleDn = -999.0;

  //------------------------------------//
  //       HADRONIC RESOLVED OBJECTS    //
  //------------------------------------//

  //Boson AK4 jet 1
  bos_j1_AK4_pt = -999.0;
  bos_j1_AK4_eta = -999.0;
  bos_j1_AK4_phi = -999.0;
  bos_j1_AK4_m = -999.0;
  bos_j1_AK4_e = -999.0;

  //Boson AK4 jet 1 variations
  bos_j1_AK4_pt_scaleUp = -999.0;
  bos_j1_AK4_pt_scaleDn = -999.0;
  bos_j1_AK4_m_scaleUp = -999.0;
  bos_j1_AK4_m_scaleDn = -999.0;
  bos_j1_AK4_e_scaleUp = -999.0;
  bos_j1_AK4_e_scaleDn = -999.0;

  //Boson AK4 jet 2
  bos_j2_AK4_pt = -999.0;
  bos_j2_AK4_eta = -999.0;
  bos_j2_AK4_phi = -999.0;
  bos_j2_AK4_m = -999.0;
  bos_j2_AK4_e = -999.0;

  //Boson AK4 jet 2 variations
  bos_j2_AK4_pt_scaleUp = -999.0;
  bos_j2_AK4_pt_scaleDn = -999.0;
  bos_j2_AK4_m_scaleUp = -999.0;
  bos_j2_AK4_m_scaleDn = -999.0;
  bos_j2_AK4_e_scaleUp = -999.0;
  bos_j2_AK4_e_scaleDn = -999.0;

  //Boson dijet object
  bos_AK4AK4_pt = -999.0;
  bos_AK4AK4_eta = -999.0;
  bos_AK4AK4_phi = -999.0;
  bos_AK4AK4_m = -999.0;

};

