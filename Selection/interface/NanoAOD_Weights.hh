//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 12 03:34:06 2020 by ROOT version 6.12/07
// from TTree Runs/Runs
// found on file: root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root
//////////////////////////////////////////////////////////

#ifndef NanoAOD_Weights_h
#define NanoAOD_Weights_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class NanoAOD_Weights {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxgenEventCount = 1;
   static constexpr Int_t kMaxgenEventSumw = 1;
   static constexpr Int_t kMaxgenEventSumw2 = 1;
   static constexpr Int_t kMaxnLHEScaleSumw = 200;
   static constexpr Int_t kMaxLHEScaleSumw = 1;
   static constexpr Int_t kMaxnLHEPdfSumw = 200;
   static constexpr Int_t kMaxLHEPdfSumw = 1;

   // Declaration of leaf types
   UInt_t          run;
   Long64_t        genEventCount_;
   Double_t        genEventSumw_;
   Double_t        genEventSumw2_;
   UInt_t          nLHEScaleSumw_;
   Double_t        LHEScaleSumw_[kMaxnLHEScaleSumw];   //[nLHEScaleSumw_]
   UInt_t          nLHEPdfSumw_;
   Double_t        LHEPdfSumw_[kMaxnLHEPdfSumw];   //[nLHEPdfSumw_]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_genEventCount_;   //!
   TBranch        *b_genEventSumw_;   //!
   TBranch        *b_genEventSumw2_;   //!
   TBranch        *b_nLHEScaleSumw_;   //!
   TBranch        *b_LHEScaleSumw_;   //!
   TBranch        *b_nLHEPdfSumw_;   //!
   TBranch        *b_LHEPdfSumw_;   //!

  NanoAOD_Weights(TTree *tree=0) {
    if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root");
      if (!f || !f->IsOpen()) {
	f = new TFile("root://cmseos.fnal.gov//store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6/WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/200408_094108/FB0B3C91-4F33-9D41-BF4C-D6677CEFD73A_Skim.root");
      }
      f->GetObject("Runs",tree);
      
    }
    Init(tree);
  };
  virtual ~NanoAOD_Weights() {
    if (!fChain) return;
    delete fChain->GetCurrentFile();
  };
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
};

#endif
Int_t NanoAOD_Weights::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NanoAOD_Weights::LoadTree(Long64_t entry)
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

void NanoAOD_Weights::Init(TTree *tree)
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
   fChain->SetBranchAddress("genEventCount_", &genEventCount_, &b_genEventCount_);
   fChain->SetBranchAddress("genEventSumw_", &genEventSumw_, &b_genEventSumw_);
   fChain->SetBranchAddress("genEventSumw2_", &genEventSumw2_, &b_genEventSumw2_);
   fChain->SetBranchAddress("nLHEScaleSumw_", &nLHEScaleSumw_, &b_nLHEScaleSumw_);
   fChain->SetBranchAddress("LHEScaleSumw_", &LHEScaleSumw_, &b_LHEScaleSumw_);
   fChain->SetBranchAddress("nLHEPdfSumw_", &nLHEPdfSumw_, &b_nLHEPdfSumw_);
   fChain->SetBranchAddress("LHEPdfSumw_", &LHEPdfSumw_, &b_LHEPdfSumw_);
}
