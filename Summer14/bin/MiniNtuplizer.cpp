#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "../include/MuonLoader.hh"
#include "../include/VTaggingVariables.h"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/Selector.hh"

#include "CMSTopTagger.hh"
#include "fastjet/tools/JHTopTagger.hh"
#include "SubstrAna/Summer14/src/HEPTopTaggerWrapper.hh"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "TRandom3.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace fastjet;
using namespace contrib;

typedef vector<float> vfloat;
typedef vector<bool> vbool;

//Object Processors
GenLoader       *fGen      = 0; 
MuonLoader      *fMuon     = 0; 
PFLoader        *fPFCand   = 0; 

TClonesArray *fJet;
TBranch      *fJetBr;

TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}

// jet clustering R size
double jetR ;

// object for VTagging evaluation
VTaggingVariables vtagger;

// SoftDrop paramenters
double beta, symmetry_cut, R0;

//trimming paramenters
double R_trimming, PtFraction;

//pruning paramenters
double R_Cut, z_cut, R_jet_def_pruning;

//matching thresholds 
double dRMatching ;

//random seed
TRandom3 randNumber ;

fastjet::JetAlgorithm algorithm_Trimming, algorithm_Pruning;

///////// Structure used in order to fill output tree branches
class GenJetInfo {

 public :

  GenJetInfo(){};
  ~GenJetInfo(){};

  int npu ;

  //ungroomed jet variables
  vector<float> pt;     
  vector<float> ptcorr; // JEC applied
  vector<float> ptraw;  // raw pt 
  vector<float> ptunc;  // uncertainty on the pt
  vector<float> eta;
  vector<float> phi;
  vector<float> m;      // mass after JEC
  vector<float> mraw;
  
  vector<float> ptconst;
  vector<float> mconst;
  vector<float> ptclean;
  vector<float> mclean;
  vector<float> pttrim;
  vector<float> mtrim;
  vector<float> pttrimsafe;
  vector<float> mtrimsafe;
  vector<float> ptpruned;
  vector<float> mpruned;
  vector<float> ptprunedsafe;
  vector<float> mprunedsafe;
  vector<float> ptsoftdrop;
  vector<float> msoftdrop;
  vector<float> ptsoftdropsafe;
  vector<float> msoftdropsafe;

  vector<float> ptCA ;
  vector<float> ptCAcorr ;
  vector<float> ptCAraw ;
  vector<float> msoftdropCA ;
  vector<float> msoftdropCAsafe;

  vector<float> sdsymmetry ;
  vector<float> sddeltar ;
  vector<float> sdmu ;
  vector<float> sdenergyloss ;
  vector<float> sdarea ;
  vector<float> sdnconst ;
  vector<float> mfiltsoftdropCA ;

  vector<float> hepmass ;
  vector<float> hepmasscorr ;
  vector<float> hepwmass ;
  vector<float> hepm01 ;
  vector<float> hepm02 ;
  vector<float> hepm12 ;
  vector<float> hepm12m012 ;
  vector<float> hepatanm02m01;

  // V-tagging base variables
  vector<float> tau1;
  vector<float> tau2;
  vector<float> tau3;
  vector<float> tau4;
  vector<float> tau5;

  vector<float> tau1_pr;
  vector<float> tau2_pr;
  vector<float> tau3_pr;
  vector<float> tau4_pr;
  vector<float> tau5_pr;

  vector<float> tau1_softdrop;
  vector<float> tau2_softdrop;
  vector<float> tau3_softdrop;
  vector<float> tau4_softdrop;
  vector<float> tau5_softdrop;

  vector<float> Qjets;

  vector<float> charge_k05;
  vector<float> charge_k07;
  vector<float> charge_k10;

  vector<float> ecf_b05;
  vector<float> ecf_b10;
  vector<float> ecf_b15;
  vector<float> ecf_b20;

  // constituents info
  vector<int>   nparticles;
  vector<int>   nneutrals;
  vector<int>   ncharged;

  vector<float> cmsmass ;
  vector<float> cmsminmass ;
  vector<float> cmshelicity ;
  vector<float> cmsnsubjets ;
  vector<float> cmsarea ;
  vector<float> cmsmasscorr ;
  vector<float> cmsminmasscorr;


};

class JetInfo : public GenJetInfo {

 public:

  JetInfo(){};
  ~JetInfo(){};

  // gen level info
  vector<float> ptgen;
  vector<float> etagen;
  vector<float> phigen;
  vector<float> mgen;
  vector<float> mrawgen;
  vector<float> mtrimgen;
  vector<float> mtrimsafegen;
  vector<float> mcleangen; //needed?
  vector<float> mconstgen;//needed?
  vector<int>   imatch;

  vector<float> ptCAgen ;
  vector<float> msoftdropCAgen ;
  vector<float> msoftdropCAsafegen;
  vector<float> mfiltsoftdropCAgen;
  
  //matching to the Boson
  vector <bool> is_MatchedToBoson;

};

///////// divide jet particles after ClusterSequence in neutrals, charged from PV and PU charged
void getConstitsForCleansing(vector<PseudoJet> inputs, vector<PseudoJet> &oNeutrals, vector<PseudoJet> &oChargedLV, vector<PseudoJet> &oChargedPU){
  for (unsigned int i = 0; i < inputs.size(); i++){
    if (inputs[i].user_index() <= 1) oNeutrals.push_back(inputs[i]);
    if (inputs[i].user_index() == 2) oChargedLV.push_back(inputs[i]);
    if (inputs[i].user_index() == 3) oChargedPU.push_back(inputs[i]);
  }
}

//// Selector for charged particles acting on pseudojets
class SW_IsPupCharged : public SelectorWorker {
public:
  SW_IsPupCharged(){}
  virtual bool pass(const PseudoJet & jet) const {
    return (jet.user_index() > 1);
  }
};

Selector SelectorIsPupCharged(){
  return Selector(new SW_IsPupCharged());
}


//// Selector for charged particles from pile up vertexes acting on pseudojets
class SW_IsPupVertex : public SelectorWorker {
public:
  SW_IsPupVertex(){}
  virtual bool pass(const PseudoJet & jet) const {
    return (jet.user_index() == 2);
  }
};


Selector SelectorIsPupVertex(){
  return Selector(new SW_IsPupVertex());
}

//// function to get JEC from the factorize corrector
double correction( PseudoJet &iJet,FactorizedJetCorrector *iJetCorr,double iRho) { 
  iJetCorr->setJetPt (iJet.pt());
  iJetCorr->setJetEta(iJet.eta());
  iJetCorr->setJetPhi(iJet.phi());
  iJetCorr->setJetE  (iJet.e());
  iJetCorr->setJetA  (iJet.area());
  iJetCorr->setRho(iRho);
  iJetCorr->setJetEMF(-99.0);
  double jetcorr= iJetCorr->getCorrection();
  return jetcorr;
}


//// function to get JEC unncertainty
double unc( PseudoJet &iJet,JetCorrectionUncertainty *iJetUnc) { 
  if(fabs(iJet.eta()) > 5. || fabs(iJet.pt()) < 10.) return 1.;
  iJetUnc->setJetPt ( iJet.pt()  );
  iJetUnc->setJetEta( iJet.eta() );
  double jetunc = iJetUnc->getUncertainty(true);
  return jetunc;
}


//// function to match a jet in a collection of other jets --> dR = 0.3 set by default
int matchingIndex(const PseudoJet & jet, const vector<PseudoJet> & genjets) {
  float rmin = 9999.;  
  int imatch = -1;
  for(unsigned int i = 0; i < genjets.size(); i++) {
    float rtemp = jet.delta_R(genjets[i]);
    if ( rtemp > dRMatching ) continue;
    if ( rtemp < rmin ){
      rmin =  rtemp;
      imatch = i;
    }
  }
  return (imatch);  
}



//// function to match a jet in a collection of other jets --> dR = 0.3 set by default
int matchingIndexFromJetInfo(const PseudoJet& jet, const GenJetInfo& jetInfo) {
  float rmin = 9999.;  
  int imatch = -1;
  for(unsigned int i = 0; i < (jetInfo.eta).size(); i++) {
    double dEta = fabs(jet.eta() - (jetInfo.eta)[i]);
    double dPhi = fabs(jet.phi() - (jetInfo.phi)[i]);
    if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
    float rtemp = sqrt(dEta*dEta+dPhi*dPhi);
    if ( rtemp > dRMatching ) continue;
    if ( rtemp < rmin ){
      rmin =  rtemp;
      imatch = i;
    }
  }
  return (imatch);  
}

//// function to match a jet with a vector boson (generic 4V given phi and eta)
bool IsMatchedToGenBoson(const vfloat& eta, const vfloat& phi, const PseudoJet& Jet) {
 bool IsMatched=false;
   
  for (unsigned int iGen =0; iGen < eta.size(); ++iGen){
      double dEta = fabs(eta.at(iGen) - (Jet.eta()));
      double dPhi = fabs(phi.at(iGen) - (Jet.phi()));
      if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
      float rtemp = sqrt(dEta*dEta+dPhi*dPhi);
      if ( rtemp < dRMatching ){
	IsMatched = true;
      }
  }  
  return (IsMatched);  
}



//// Set the output tree structure
void setupGenTree(TTree *iTree, GenJetInfo &iJet, std::string iName) {
  iTree->Branch((iName+"npu"       ).c_str(),&iJet.npu       );

  iTree->Branch((iName+"pt"        ).c_str(),&iJet.pt        );   
  iTree->Branch((iName+"ptcorr"    ).c_str(),&iJet.ptcorr    );
  iTree->Branch((iName+"ptraw"     ).c_str(),&iJet.ptraw     );
  iTree->Branch((iName+"ptunc"     ).c_str(),&iJet.ptunc     );
  iTree->Branch((iName+"eta"       ).c_str(),&iJet.eta       );
  iTree->Branch((iName+"phi"       ).c_str(),&iJet.phi       );
  iTree->Branch((iName+"m"         ).c_str(),&iJet.m         );
  iTree->Branch((iName+"mraw"      ).c_str(),&iJet.mraw      );

  iTree->Branch((iName+"ptclean"   ).c_str(),&iJet.ptclean   );
  iTree->Branch((iName+"mclean"    ).c_str(),&iJet.mclean    );

  iTree->Branch((iName+"pttrim"    ).c_str(),&iJet.pttrim    );
  iTree->Branch((iName+"mtrim"     ).c_str(),&iJet.mtrim     );

  iTree->Branch((iName+"pttrimsafe").c_str(),&iJet.pttrimsafe);
  iTree->Branch((iName+"mtrimsafe" ).c_str(),&iJet.mtrimsafe );

  iTree->Branch((iName+"ptconst"   ).c_str(),&iJet.ptconst   );
  iTree->Branch((iName+"mconst"    ).c_str(),&iJet.mconst    );

  iTree->Branch((iName+"ptpruned" ).c_str(),&iJet.ptpruned );
  iTree->Branch((iName+"mpruned" ).c_str(),&iJet.mpruned );

  iTree->Branch((iName+"ptprunedsafe" ).c_str(),&iJet.ptprunedsafe );
  iTree->Branch((iName+"mprunedsafe" ).c_str(),&iJet.mprunedsafe );

  iTree->Branch((iName+"ptsoftdrop" ).c_str(),&iJet.msoftdrop);
  iTree->Branch((iName+"msoftdrop" ).c_str(),&iJet.msoftdrop);

  iTree->Branch((iName+"ptsoftdropsafe" ).c_str(),&iJet.msoftdropsafe);
  iTree->Branch((iName+"msoftdropsafe" ).c_str(),&iJet.msoftdropsafe);

  iTree->Branch((iName+"nparticles").c_str(),&iJet.nparticles);
  iTree->Branch((iName+"nneutrals" ).c_str(),&iJet.nneutrals);
  iTree->Branch((iName+"ncharged"  ).c_str(),&iJet.ncharged);

  iTree->Branch((iName+"ptCA" ).c_str(),&iJet.ptCA );
  iTree->Branch((iName+"ptCAcorr" ).c_str(),&iJet.ptCAcorr );
  iTree->Branch((iName+"ptCAraw" ).c_str(),&iJet.ptCAraw );
  iTree->Branch((iName+"msoftdropCA" ).c_str(),&iJet.msoftdropCA );
  iTree->Branch((iName+"msoftdropCAsafe" ).c_str(),&iJet.msoftdropCAsafe );

  iTree->Branch((iName+"sdsymmetry" ).c_str(),&iJet.sdsymmetry );
  iTree->Branch((iName+"sddeltar" ).c_str(),&iJet.sddeltar );
  iTree->Branch((iName+"sdmu" ).c_str(),&iJet.sdmu );
  iTree->Branch((iName+"sdenergyloss" ).c_str(),&iJet.sdenergyloss );
  iTree->Branch((iName+"sdarea" ).c_str(),&iJet.sdarea );
  iTree->Branch((iName+"sdnconst" ).c_str(),&iJet.sdnconst );
  iTree->Branch((iName+"mfiltsoftdropCA" ).c_str(),&iJet.mfiltsoftdropCA );


  iTree->Branch((iName+"tau1"  ).c_str(),&iJet.tau1);
  iTree->Branch((iName+"tau2"  ).c_str(),&iJet.tau2);
  iTree->Branch((iName+"tau3"  ).c_str(),&iJet.tau3);
  iTree->Branch((iName+"tau4"  ).c_str(),&iJet.tau4);
  iTree->Branch((iName+"tau5"  ).c_str(),&iJet.tau5);

  iTree->Branch((iName+"tau1_pr"  ).c_str(),&iJet.tau1_pr);
  iTree->Branch((iName+"tau2_pr"  ).c_str(),&iJet.tau2_pr);
  iTree->Branch((iName+"tau3_pr"  ).c_str(),&iJet.tau3_pr);
  iTree->Branch((iName+"tau4_pr"  ).c_str(),&iJet.tau4_pr);
  iTree->Branch((iName+"tau5_pr"  ).c_str(),&iJet.tau5_pr);

  iTree->Branch((iName+"tau1_softdrop"  ).c_str(),&iJet.tau1_softdrop);
  iTree->Branch((iName+"tau2_softdrop"  ).c_str(),&iJet.tau2_softdrop);
  iTree->Branch((iName+"tau3_softdrop"  ).c_str(),&iJet.tau3_softdrop);
  iTree->Branch((iName+"tau4_softdrop"  ).c_str(),&iJet.tau4_softdrop);
  iTree->Branch((iName+"tau5_softdrop"  ).c_str(),&iJet.tau5_softdrop);

  iTree->Branch((iName+"Qjets"  ).c_str(),&iJet.Qjets);

  iTree->Branch((iName+"charge_k05"  ).c_str(),&iJet.charge_k05);
  iTree->Branch((iName+"charge_k07"  ).c_str(),&iJet.charge_k07);
  iTree->Branch((iName+"charge_k10"  ).c_str(),&iJet.charge_k10);

  iTree->Branch((iName+"ecf_b05"  ).c_str(),&iJet.ecf_b05);
  iTree->Branch((iName+"ecf_b10"  ).c_str(),&iJet.ecf_b10);
  iTree->Branch((iName+"ecf_b15"  ).c_str(),&iJet.ecf_b15);
  iTree->Branch((iName+"ecf_b20"  ).c_str(),&iJet.ecf_b20);

  iTree->Branch((iName+"hepmass" ).c_str(),&iJet.hepmass );
  iTree->Branch((iName+"hepmasscorr" ).c_str(),&iJet.hepmasscorr );
  iTree->Branch((iName+"hepwmass" ).c_str(),&iJet.hepwmass );
  iTree->Branch((iName+"hepm01" ).c_str(),&iJet.hepm01 );
  iTree->Branch((iName+"hepm02" ).c_str(),&iJet.hepm02 );
  iTree->Branch((iName+"hepm12" ).c_str(),&iJet.hepm12 );
  iTree->Branch((iName+"hepm12m012" ).c_str(),&iJet.hepm12m012 );
  iTree->Branch((iName+"hepatanm02m01" ).c_str(),&iJet.hepatanm02m01 );

  iTree->Branch((iName+"cmsmass" ).c_str(),&iJet.cmsmass );
  iTree->Branch((iName+"cmsminmass" ).c_str(),&iJet.cmsminmass );
  iTree->Branch((iName+"cmshelicity" ).c_str(),&iJet.cmshelicity );
  iTree->Branch((iName+"cmsnsubjets" ).c_str(),&iJet.cmsnsubjets );
  iTree->Branch((iName+"cmsarea" ).c_str(),&iJet.cmsarea );
  iTree->Branch((iName+"cmsmasscorr" ).c_str(),&iJet.cmsmasscorr );
  iTree->Branch((iName+"cmsminmasscorr" ).c_str(),&iJet.cmsminmasscorr );



}

void setupTree(TTree *iTree, JetInfo &iJet, std::string iName) {

  setupGenTree(iTree,iJet,iName); 

  // gen info
  iTree->Branch((iName+"ptgen"       ).c_str(),&iJet.ptgen       );
  iTree->Branch((iName+"etagen"      ).c_str(),&iJet.etagen      );
  iTree->Branch((iName+"phigen"      ).c_str(),&iJet.phigen      );
  iTree->Branch((iName+"mgen"        ).c_str(),&iJet.mgen        );
  iTree->Branch((iName+"mrawgen"     ).c_str(),&iJet.mrawgen     );//needed?
  iTree->Branch((iName+"mtrimgen"    ).c_str(),&iJet.mtrimgen    );//needed?
  iTree->Branch((iName+"mtrimsafegen").c_str(),&iJet.mtrimsafegen);//needed?
  iTree->Branch((iName+"mcleangen"   ).c_str(),&iJet.mcleangen   );//needed?
  iTree->Branch((iName+"mconstgen"   ).c_str(),&iJet.mconstgen   );//needed?
  iTree->Branch((iName+"imatch"      ).c_str(),&iJet.imatch      );
  
  //matched to the boson
  iTree->Branch((iName+"is_MatchedToBoson"      ).c_str(),&iJet.is_MatchedToBoson      );
 
  iTree->Branch((iName+"ptCAgen" ).c_str(),&iJet.ptCAgen );
  iTree->Branch((iName+"msoftdropCAgen" ).c_str(),&iJet.msoftdropCAgen );
  iTree->Branch((iName+"msoftdropCAsafegen" ).c_str(),&iJet.msoftdropCAsafegen );
  iTree->Branch((iName+"mfiltsoftdropCAgen" ).c_str(),&iJet.mfiltsoftdropCAgen );

  
}


// clear tree structure content at the beginning of each event
void clear(GenJetInfo &iJet) {
  iJet.npu  = -1;

  iJet.pt         .clear();
  iJet.ptcorr     .clear();
  iJet.ptraw      .clear();
  iJet.eta        .clear();
  iJet.phi        .clear();
  iJet.m          .clear();
  iJet.mraw       .clear();

  iJet.mtrim      .clear();
  iJet.pttrim     .clear();
  iJet.mtrimsafe  .clear();
  iJet.pttrimsafe .clear();
  iJet.mpruned    .clear();
  iJet.ptpruned    .clear();
  iJet.mprunedsafe.clear();
  iJet.ptprunedsafe.clear();
  iJet.msoftdrop  .clear();
  iJet.ptsoftdrop  .clear();
  iJet.msoftdropsafe.clear();
  iJet.ptsoftdropsafe.clear();
  iJet.mclean     .clear();
  iJet.ptclean    .clear();

  iJet.mconst     .clear();
  iJet.nparticles .clear();
  iJet.nneutrals  .clear();
  iJet.ncharged   .clear();

  iJet.ptCA .clear();
  iJet.ptCAcorr .clear();
  iJet.ptCAraw .clear();
  iJet.msoftdropCA .clear();
  iJet.msoftdropCAsafe.clear();
  iJet.sdsymmetry .clear();
  iJet.sddeltar .clear();
  iJet.sdmu .clear();
  iJet.sdenergyloss .clear();
  iJet.sdarea .clear();
  iJet.sdnconst .clear();
  iJet.mfiltsoftdropCA.clear();


  iJet.tau1.clear();
  iJet.tau2.clear();
  iJet.tau3.clear();
  iJet.tau4.clear();
  iJet.tau5.clear();

  iJet.tau1_pr.clear();
  iJet.tau2_pr.clear();
  iJet.tau3_pr.clear();
  iJet.tau4_pr.clear();
  iJet.tau5_pr.clear();

  iJet.tau1_softdrop.clear();
  iJet.tau2_softdrop.clear();
  iJet.tau3_softdrop.clear();
  iJet.tau4_softdrop.clear();
  iJet.tau5_softdrop.clear();

  iJet.Qjets.clear();

  iJet.charge_k05.clear();
  iJet.charge_k07.clear();
  iJet.charge_k10.clear();

  iJet.ecf_b05.clear();
  iJet.ecf_b10.clear();
  iJet.ecf_b15.clear();
  iJet.ecf_b20.clear();

  iJet.hepmass .clear();
  iJet.hepmasscorr .clear();
  iJet.hepwmass .clear();
  iJet.hepm01 .clear();
  iJet.hepm02 .clear();
  iJet.hepm12 .clear();
  iJet.hepm12m012 .clear();
  iJet.hepatanm02m01 .clear();
  iJet.cmsmass .clear();
  iJet.cmsminmass .clear();
  iJet.cmshelicity .clear();
  iJet.cmsnsubjets .clear();
  iJet.cmsarea .clear();
  iJet.cmsmasscorr .clear();
  iJet.cmsminmasscorr .clear();


}

void clear(JetInfo &iJet) {

  clear((GenJetInfo &)iJet);

  iJet.ptgen       .clear();
  iJet.etagen      .clear();
  iJet.phigen      .clear();
  iJet.mgen        .clear();
  iJet.mrawgen     .clear();
  iJet.mtrimgen    .clear();
  iJet.mtrimsafegen.clear();
  iJet.mcleangen   .clear();
  iJet.mconstgen   .clear();

  iJet.imatch      .clear();
  iJet.is_MatchedToBoson.clear();

  iJet.ptCAgen .clear();
  iJet.msoftdropCAgen .clear();
  iJet.msoftdropCAsafegen.clear();
  iJet.mfiltsoftdropCAgen.clear();

}


void setJet(PseudoJet &iJet, JetInfo &iJetI, JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, 
	    bool isGEN, bool isCHS, FactorizedJetCorrector *iJetCorr, JetCorrectionUncertainty *iJetUnc, JetCleanser &gsn_cleanser, 
	    bool doGenMatching, vector<PseudoJet> genJets, vfloat eta_Boson, vfloat phi_Boson) {

  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  if(!isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  if( isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom,SelectorIsPupCharged(),SelectorIsPupVertex());
  PseudoJet lCorr =  (*area_subtractor)(iJet);
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet);

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU;
  getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
  PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
  
  // -- trimming
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(algorithm_Trimming, R_trimming), fastjet::SelectorPtFractionMin(PtFraction)));
  
  PseudoJet lTrim     = (trimmer)(iJet);
  trimmer.set_subtractor(area_subtractor);
  
  PseudoJet lTrimSafe = (trimmer)(iJet);
 
   //pruning
  JetDefinition jet_def_Pruning(algorithm_Pruning, R_jet_def_pruning);
  Pruner pruner(jet_def_Pruning, z_cut,R_Cut);
  PseudoJet lPruned = pruner(iJet);
  PseudoJet lPrunedSafe = pruner(lCorr);
 
  //softdrop
  contrib::SoftDrop softdrop(beta, symmetry_cut, R0);
  PseudoJet lSoftDropped = softdrop(iJet);
  softdrop.set_subtractor(area_subtractor);
  PseudoJet lSoftDroppedSafe = softdrop(iJet);
  
  
  // -- apply the JEC
  double lJEC = 1.;
  double lUnc = 0 ;
  if (!isGEN){
    lJEC = correction(iJet,iJetCorr,bge_rho.rho());  
    lUnc = unc       (iJet,iJetUnc);
  }

  // -- find the gen jet matched to this reco jet
  int imatch = -1;
  bool matched = IsMatchedToGenBoson(eta_Boson, phi_Boson, iJet);
  if (doGenMatching) imatch = matchingIndex(iJet,genJets);

  // -- fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());
  (iJetI.ptcorr    ).push_back(iJet      .pt()*lJEC);
  (iJetI.ptraw     ).push_back(iJet      .pt());
  (iJetI.ptclean   ).push_back(lClean    .pt());
  (iJetI.pttrim    ).push_back(lTrim     .pt());
  (iJetI.pttrimsafe).push_back(lTrimSafe .pt());
  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.ptunc     ).push_back(lUnc);
  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());
  (iJetI.mclean    ).push_back(lClean    .m());
  (iJetI.mtrim     ).push_back(lTrim     .m());
  (iJetI.mtrimsafe ).push_back(lTrimSafe .m());
  (iJetI.mpruned   ).push_back(lPruned   .m());
  (iJetI.mprunedsafe).push_back(lPrunedSafe.m());
  (iJetI.msoftdrop).push_back(lSoftDropped.m());
  (iJetI.msoftdropsafe).push_back(lSoftDroppedSafe.m());
  (iJetI.mconst    ).push_back(lConstit  .m());
  (iJetI.nparticles).push_back((iJet.constituents()).size());
  (iJetI.nneutrals ).push_back(neutrals.size());
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());
  (iJetI.is_MatchedToBoson ).push_back(matched);

  if (imatch > -1){
    (iJetI.imatch  ).push_back(imatch);
    (iJetI.ptgen    ).push_back(genJets[imatch].pt());
    (iJetI.etagen   ).push_back(genJets[imatch].eta());
    (iJetI.phigen   ).push_back(genJets[imatch].phi());
    (iJetI.mgen     ).push_back(genJets[imatch].m());

  }
  else {
    (iJetI.imatch   ).push_back(imatch);
    (iJetI.ptgen    ).push_back(-999.);
    (iJetI.etagen   ).push_back(-999.);
    (iJetI.phigen   ).push_back(-999.);
    (iJetI.mgen     ).push_back(-999.);

  }
  
}


void setRecoJet(PseudoJet &iJet, JetInfo &iJetI, GenJetInfo& iGenJetI, JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, bool isCHS, FactorizedJetCorrector *iJetCorr, JetCorrectionUncertainty *iJetUnc, JetCleanser &gsn_cleanser, vfloat eta_Boson, vfloat phi_Boson) {

  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  if(!isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  if( isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom,SelectorIsPupCharged(),SelectorIsPupVertex());
  PseudoJet lCorr =  (*area_subtractor)(iJet);
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet);

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU;
  getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
  PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
  
   
  // -- trimming
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(algorithm_Trimming, R_trimming), fastjet::SelectorPtFractionMin(PtFraction)));
  
  PseudoJet lTrim     = (trimmer)(iJet);
  trimmer.set_subtractor(area_subtractor);
  PseudoJet lTrimSafe = (trimmer)(iJet);
 
  //pruning
  JetDefinition jet_def_Pruning(algorithm_Pruning, R_jet_def_pruning);
  Pruner pruner(jet_def_Pruning, z_cut, R_Cut);
  PseudoJet lPruned = pruner(iJet);
  //PseudoJet lPrunedSafe = pruner(lCorr);
  PseudoJet lPrunedSafe = (*area_subtractor)(lPruned);

  //softdrop
  contrib::SoftDrop softdrop(beta, symmetry_cut, R0);
  PseudoJet lSoftDropped = softdrop(iJet);
  softdrop.set_subtractor(area_subtractor);
  PseudoJet lSoftDroppedSafe = softdrop(iJet);
  
  // -- apply the JEC
  double lJEC = correction(iJet,iJetCorr,bge_rho.rho());  
  double lUnc = unc       (iJet,iJetUnc);

  // -- recluster jet CA
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  JetDefinition jet_def_CA (fastjet::cambridge_algorithm, jetR*10); //large R to cluster all constituents of original jet
  fastjet::ClusterSequenceArea cs_Recluster (iJet.constituents(), jet_def_CA, area_def);
  vector<fastjet::PseudoJet> jets_Recluster = sorted_by_pt(cs_Recluster.inclusive_jets());
  fastjet::PseudoJet iJetCA = jets_Recluster[0];
  PseudoJet lCorrCA = (*area_subtractor)(iJetCA);

  double lJEC_CA = correction(iJetCA,iJetCorr,bge_rho.rho());

  // -- softdrop reclustered jet
  contrib::SoftDrop softdropCA(beta,symmetry_cut,R0);
  PseudoJet lSoftDroppedCA = softdropCA(iJetCA);
  softdropCA.set_subtractor(area_subtractor);
  PseudoJet lSoftDroppedCASafe = softdrop(iJetCA);
 

  double SoftDropedSymmetry = -1.0;
  double SoftDropedDR = -1.0;
  double SoftDropedMassDrop = -1.0;
  double SoftDropedEnergyLoss = -1.0;
  double SoftDropedArea = -1.0;
  double SoftDropedNconst = -1.0;
  PseudoJet filtered_softdropped_jet;

  if (lSoftDroppedCASafe!=0 and lSoftDroppedCASafe.m()>0.0){

    SoftDropedSymmetry = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().symmetry();
    SoftDropedDR = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().delta_R();
    SoftDropedMassDrop = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().mu();
    SoftDropedEnergyLoss = 1-lSoftDroppedCASafe.pt()/iJetCA.pt();
    SoftDropedArea = lSoftDroppedCASafe .area() ;
    SoftDropedNconst = lSoftDroppedCASafe .constituents().size() ;

    // filter jet dynamically based on deltaR between subjets (arXiv:0802.2470)
    double dyn_Rfilt = min(0.3, SoftDropedDR*0.5);
    int dyn_nfilt = 3;
    Filter filtersoft(dyn_Rfilt, SelectorNHardest(dyn_nfilt));
    filtered_softdropped_jet = filtersoft(lSoftDroppedCASafe);
  }

  // -- find the gen jet matched to this reco jet
  int imatch = matchingIndexFromJetInfo(iJet,iGenJetI);
  bool matched = IsMatchedToGenBoson( eta_Boson, phi_Boson, iJet);


  // Top taggers 
  double mass_drop_threshold=0.8;
  double max_subjet_mass=30;
  bool use_subjet_mass_cuts=false;
  HEPTopTagger hep_top_tagger(mass_drop_threshold, max_subjet_mass, use_subjet_mass_cuts);
  
  PseudoJet hep_top_candidate = hep_top_tagger( iJetCA );
  double hepttJetMass = -1;
  double hepttJetMassCorr= -1;
  double hepttWMass = -1;
  double hepttM01 = -1;
  double hepttM02 = -1;
  double hepttM12 = -1;
  double hepttM12M012 = -1;
  double hepttAtanM02M01 = -1;
  if (hep_top_candidate != 0){
      PseudoJet W = hep_top_candidate.structure_of<HEPTopTagger>().W();
      PseudoJet W1 = hep_top_candidate.structure_of<HEPTopTagger>().W1();
      PseudoJet W2 = hep_top_candidate.structure_of<HEPTopTagger>().W2();
      PseudoJet non_W = hep_top_candidate.structure_of<HEPTopTagger>().non_W();

      vector<PseudoJet> all_subjets;
      all_subjets.push_back(W1);
      all_subjets.push_back(W2);
      all_subjets.push_back(non_W);
      all_subjets = sorted_by_pt(all_subjets);

      PseudoJet sum012 = all_subjets[0]+all_subjets[1]+all_subjets[2];
      PseudoJet sum01 = all_subjets[0]+all_subjets[1];
      PseudoJet sum02 = all_subjets[0]+all_subjets[2];
      PseudoJet sum12 = all_subjets[1]+all_subjets[2];

      hepttJetMass = hep_top_candidate.m();
      hepttWMass = W.m();
      hepttM01 = sum01.m();
      hepttM02 = sum02.m();
      hepttM12 = sum12.m();
      if ( sum012.m()!=0 ) hepttM12M012 = sum12.m() / sum012.m() ;
      if ( sum01.m()!=0 ) hepttAtanM02M01 = atan( sum02.m() / sum01.m() ) ;

      //PseudoJet lCorrHEP = (*area_subtractor)(hep_top_candidate);
      hepttJetMassCorr = hep_top_candidate.m();

    }


  // -- CMS Top Tagger
  double cms_delta_p = 0.05;
  double cms_delta_r=0.4;
  double A=0.0004;

  CMSTopTagger cms_top_tagger(cms_delta_p, cms_delta_r, A);
  PseudoJet cms_top_candidate = cms_top_tagger( iJetCA );

  double cmsttJetMass = -1;
  double cmsttMinMass = -1;
  double cmsttHelicity = -1;
  double cmsttNsubjets = -1;
  double cmsttArea = -1;
  double cmsttJetMassCorr = -1;
  double cmsttMinMassCorr = -1;

  if (cms_top_candidate != 0){
      vector<PseudoJet> kept_subjets0 = cms_top_candidate.structure_of<CMSTopTagger>().W().pieces();
      vector<PseudoJet> kept_subjets1 = cms_top_candidate.structure_of<CMSTopTagger>().non_W().pieces();
      vector<PseudoJet> all_subjets = kept_subjets0;
      all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

      PseudoJet lCorrCMS = (*area_subtractor)(cms_top_candidate);

      cmsttJetMass = cms_top_candidate.m();
      cmsttMinMass = cms_top_candidate.structure_of<CMSTopTagger>().W().m();
      cmsttHelicity = cms_top_candidate.structure_of<CMSTopTagger>().cos_theta_W();
      cmsttNsubjets = all_subjets.size();
      cmsttArea = cms_top_candidate.area();
      cmsttJetMassCorr = lCorrCMS.m();
      cmsttMinMassCorr = cmsttMinMass * lJEC;
  }

 
  // -- fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());
  (iJetI.ptcorr    ).push_back(iJet      .pt()*lJEC);  
  (iJetI.ptraw     ).push_back(iJet      .pt());
  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());
  (iJetI.ptunc     ).push_back(lUnc);

  (iJetI.ptclean   ).push_back(lClean    .pt());
  (iJetI.mclean    ).push_back(lClean    .m());

  (iJetI.pttrim    ).push_back(lTrim     .pt());
  (iJetI.mtrim     ).push_back(lTrim     .m());

  (iJetI.pttrimsafe).push_back(lTrimSafe .pt());
  (iJetI.mtrimsafe ).push_back(lTrimSafe .m());

  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.mconst    ).push_back(lConstit  .m());

  (iJetI.mpruned   ).push_back(lPruned   .m());
  (iJetI.ptpruned  ).push_back(lPruned   .pt());

  (iJetI.mprunedsafe).push_back(lPrunedSafe.m());
  (iJetI.ptprunedsafe).push_back(lPrunedSafe.pt());

  (iJetI.msoftdrop).push_back(lSoftDropped.m());
  (iJetI.ptsoftdrop).push_back(lSoftDropped.pt());

  (iJetI.msoftdropsafe).push_back(lSoftDroppedSafe.m());
  (iJetI.ptsoftdropsafe).push_back(lSoftDroppedSafe.pt());

  (iJetI.nparticles).push_back((iJet.constituents()).size());
  (iJetI.nneutrals ).push_back(neutrals.size());
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());
  (iJetI.is_MatchedToBoson ).push_back(matched);

  (iJetI.ptCA ).push_back(lCorrCA .pt());
  (iJetI.ptCAcorr ).push_back(iJetCA .pt()*lJEC_CA);
  (iJetI.ptCAraw ).push_back(iJetCA .pt());
  (iJetI.msoftdropCA ).push_back(lSoftDroppedCA.m());
  (iJetI.msoftdropCAsafe).push_back(lSoftDroppedCASafe.m());

  (iJetI.sdsymmetry ).push_back( SoftDropedSymmetry );
  (iJetI.sddeltar ).push_back( SoftDropedDR );
  (iJetI.sdmu ).push_back( SoftDropedMassDrop );
  (iJetI.sdenergyloss ).push_back( SoftDropedEnergyLoss );
  (iJetI.sdarea ).push_back( SoftDropedArea );
  (iJetI.sdnconst ).push_back( SoftDropedNconst );
  (iJetI.mfiltsoftdropCA ).push_back( filtered_softdropped_jet.m() );

  (iJetI.hepmass ).push_back( hepttJetMass );
  (iJetI.hepmasscorr ).push_back( hepttJetMassCorr);
  (iJetI.hepwmass ).push_back( hepttWMass );
  (iJetI.hepm01 ).push_back( hepttM01 );
  (iJetI.hepm02 ).push_back( hepttM02 );
  (iJetI.hepm12 ).push_back( hepttM12 );
  (iJetI.hepm12m012 ).push_back( hepttM12M012 );
  (iJetI.hepatanm02m01).push_back( hepttAtanM02M01 );

  (iJetI.cmsmass ).push_back(cmsttJetMass );
  (iJetI.cmsminmass ).push_back(cmsttMinMass );
  (iJetI.cmshelicity ).push_back(cmsttHelicity );
  (iJetI.cmsnsubjets ).push_back(cmsttNsubjets );
  (iJetI.cmsarea ).push_back(cmsttArea );
  (iJetI.cmsmasscorr ).push_back(cmsttJetMassCorr );
  (iJetI.cmsminmasscorr ).push_back(cmsttMinMassCorr );
 
  // V-tagging variables 
  vtagger.setInputJet(iJet); 
  (iJetI.tau1 ).push_back(vtagger.computeNSubJettines(1,1.,jetR,jetR));
  (iJetI.tau2 ).push_back(vtagger.computeNSubJettines(2,1.,jetR,jetR));
  (iJetI.tau3 ).push_back(vtagger.computeNSubJettines(3,1.,jetR,jetR));
  (iJetI.tau4 ).push_back(vtagger.computeNSubJettines(4,1.,jetR,jetR));
  (iJetI.tau5 ).push_back(vtagger.computeNSubJettines(5,1.,jetR,jetR));

  (iJetI.Qjets).push_back(vtagger.computeQjets(35,25,randNumber.Uniform(0.,10000000)));
  
  (iJetI.ecf_b05).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,0.5));
  (iJetI.ecf_b10).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,1.0));
  (iJetI.ecf_b15).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,1.5));
  (iJetI.ecf_b20).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,2.0));
    
  (iJetI.charge_k05).push_back(vtagger.computeJetChargeReco(0.5));
  (iJetI.charge_k07).push_back(vtagger.computeJetChargeReco(0.7));
  (iJetI.charge_k10).push_back(vtagger.computeJetChargeReco(1.0));
 

  vtagger.setInputJet(lPruned); 
  (iJetI.tau1_pr ).push_back(vtagger.computeNSubJettines(1,1.,jetR,jetR));
  (iJetI.tau2_pr ).push_back(vtagger.computeNSubJettines(2,1.,jetR,jetR));
  (iJetI.tau3_pr ).push_back(vtagger.computeNSubJettines(3,1.,jetR,jetR));
  (iJetI.tau4_pr ).push_back(vtagger.computeNSubJettines(4,1.,jetR,jetR));
  (iJetI.tau5_pr ).push_back(vtagger.computeNSubJettines(5,1.,jetR,jetR));

  vtagger.setInputJet(lSoftDropped); 
  (iJetI.tau1_softdrop ).push_back(vtagger.computeNSubJettines(1,1.,jetR,jetR));
  (iJetI.tau2_softdrop ).push_back(vtagger.computeNSubJettines(2,1.,jetR,jetR));
  (iJetI.tau3_softdrop ).push_back(vtagger.computeNSubJettines(3,1.,jetR,jetR));
  (iJetI.tau4_softdrop ).push_back(vtagger.computeNSubJettines(4,1.,jetR,jetR));
  (iJetI.tau5_softdrop ).push_back(vtagger.computeNSubJettines(5,1.,jetR,jetR));
  

  if (imatch > -1){
    (iJetI.imatch).push_back(imatch);
    (iJetI.ptgen    ).push_back((iGenJetI.pt)[imatch]);
    (iJetI.etagen   ).push_back((iGenJetI.eta)[imatch]);
    (iJetI.phigen   ).push_back((iGenJetI.phi)[imatch]);
    (iJetI.mgen     ).push_back((iGenJetI.m)[imatch]);
    (iJetI.mrawgen     ).push_back((iGenJetI.mraw)[imatch]);
    (iJetI.mtrimgen    ).push_back((iGenJetI.mtrim)[imatch]);
    (iJetI.mtrimsafegen).push_back((iGenJetI.mtrimsafe)[imatch]);
    (iJetI.mcleangen   ).push_back((iGenJetI.mclean)[imatch]);
    (iJetI.mconstgen   ).push_back((iGenJetI.mconst)[imatch]);
  }
  else { 
    (iJetI.imatch).push_back(imatch);
    (iJetI.ptgen    ).push_back(-999.);
    (iJetI.etagen   ).push_back(-999.);
    (iJetI.phigen   ).push_back(-999.);
    (iJetI.mgen     ).push_back(-999.);
    (iJetI.mrawgen     ).push_back(-999.);
    (iJetI.mtrimgen    ).push_back(-999.);
    (iJetI.mtrimsafegen).push_back(-999.);
    (iJetI.mcleangen   ).push_back(-999.);
    (iJetI.mconstgen   ).push_back(-999.);
  }

}

//set the gen jet info in the output tree
void setGenJet(PseudoJet &iJet, GenJetInfo &iJetI,  JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, JetCleanser &gsn_cleanser) {

  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  PseudoJet lCorr =  (*area_subtractor)(iJet); // correct the jet for safe are subtraction
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet); // correct the jet for constituent subtraction

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU;
  getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
  PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU); // use cleansing
  
  // -- trimming
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(algorithm_Trimming, R_trimming), fastjet::SelectorPtFractionMin(PtFraction))); // apply trimming with the global par
  
  PseudoJet lTrim     = (trimmer)(iJet);
  
  trimmer.set_subtractor(area_subtractor); // apply safe subtraction on top of trimmed jets 
  
  PseudoJet lTrimSafe = (trimmer)(iJet); // safe subtracted trimmed jets 
  
  JetDefinition jet_def_Pruning(algorithm_Pruning, R_jet_def_pruning); // pruning and pruning safe subtracted jets 
  Pruner pruner(jet_def_Pruning, z_cut, R_Cut);
  PseudoJet lPruned = pruner(iJet);
  PseudoJet lPrunedSafe = pruner(lCorr);
   
  //softdrop
  contrib::SoftDrop softdrop(beta, symmetry_cut, R0); 
  PseudoJet lSoftDropped = softdrop(iJet);
 
  softdrop.set_subtractor(area_subtractor);
 
  PseudoJet lSoftDroppedSafe = softdrop(iJet);


  // -- recluster jet CA
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  JetDefinition jet_def_CA (fastjet::cambridge_algorithm, jetR);
  fastjet::ClusterSequenceArea cs_Recluster ( iJet.constituents(), jet_def_CA, area_def);
  vector<fastjet::PseudoJet> jets_Recluster = sorted_by_pt(cs_Recluster .inclusive_jets());
  fastjet::PseudoJet iJetCA = jets_Recluster [0];
  PseudoJet lCorrCA = (*area_subtractor)(iJetCA);

  

  // -- softdrop reclustered jet
  contrib::SoftDrop softdropCA(beta, symmetry_cut, R0);
  PseudoJet lSoftDroppedCA = softdropCA(iJetCA);
  softdropCA.set_subtractor(area_subtractor);
  PseudoJet lSoftDroppedCASafe = softdrop(iJetCA);
 
  double SoftDropedSymmetry = -1.0;
  double SoftDropedDR = -1.0;
  double SoftDropedMassDrop = -1.0;
  double SoftDropedEnergyLoss = -1.0;
  double SoftDropedArea = -1.0;
  double SoftDropedNconst = -1.0;
  PseudoJet filtered_softdropped_jet;

  if (lSoftDroppedCASafe!=0 and lSoftDroppedCASafe.m()>0.0){
      SoftDropedSymmetry = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().symmetry();
      SoftDropedDR = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().delta_R();
      SoftDropedMassDrop = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().mu();
      SoftDropedEnergyLoss = 1-lSoftDroppedCASafe.pt()/iJetCA.pt();
      SoftDropedArea = lSoftDroppedCASafe .area() ;
      SoftDropedNconst = lSoftDroppedCASafe .constituents().size() ;



      // filter jet dynamically based on deltaR between subjets (arXiv:0802.2470)
      double dyn_Rfilt = min(0.3, SoftDropedDR*0.5);
      int dyn_nfilt = 3;
      Filter filtersoft(dyn_Rfilt, SelectorNHardest(dyn_nfilt));
      filtered_softdropped_jet = filtersoft(lSoftDroppedCASafe);
  }


  // -- HEP Top Tagger
  double mass_drop_threshold=0.8;
  double max_subjet_mass=30;
  bool use_subjet_mass_cuts=false;
  HEPTopTagger hep_top_tagger(mass_drop_threshold, max_subjet_mass, use_subjet_mass_cuts);
  
  PseudoJet hep_top_candidate = hep_top_tagger( iJet );
  double hepttJetMass = -1;
  double hepttWMass = -1;
  double hepttM01 = -1;
  double hepttM02 = -1;
  double hepttM12 = -1;
  double hepttM12M012 = -1;
  double hepttAtanM02M01 = -1;
  if (hep_top_candidate != 0){
      PseudoJet W = hep_top_candidate.structure_of<HEPTopTagger>().W();
      PseudoJet W1 = hep_top_candidate.structure_of<HEPTopTagger>().W1();
      PseudoJet W2 = hep_top_candidate.structure_of<HEPTopTagger>().W2();
      PseudoJet non_W = hep_top_candidate.structure_of<HEPTopTagger>().non_W();

      vector<PseudoJet> all_subjets;
      all_subjets.push_back(W1);
      all_subjets.push_back(W2);
      all_subjets.push_back(non_W);
      all_subjets = sorted_by_pt(all_subjets);

      PseudoJet sum012 = all_subjets[0]+all_subjets[1]+all_subjets[2];
      PseudoJet sum01 = all_subjets[0]+all_subjets[1];
      PseudoJet sum02 = all_subjets[0]+all_subjets[2];
      PseudoJet sum12 = all_subjets[1]+all_subjets[2];

      hepttJetMass = hep_top_candidate.m();
      hepttWMass = W.m();
      hepttM01 = sum01.m();
      hepttM02 = sum02.m();
      hepttM12 = sum12.m();
      if ( sum012.m()!=0 ) hepttM12M012 = sum12.m() / sum012.m() ;
      if ( sum01.m()!=0 ) hepttAtanM02M01 = atan( sum02.m() / sum01.m() ) ;
  }


  // -- CMS Top Tagger
  double cms_delta_p = 0.05;
  double cms_delta_r=0.4;
  double A=0.0004;

  CMSTopTagger cms_top_tagger(cms_delta_p, cms_delta_r, A);
  PseudoJet cms_top_candidate = cms_top_tagger( iJetCA );

  double cmsttJetMass = -1;
  double cmsttMinMass = -1;
  double cmsttHelicity = -1;
  double cmsttNsubjets = -1;
  double cmsttArea = -1;
  double cmsttJetMassCorr = -1;
  double cmsttMinMassCorr = -1;

  if (cms_top_candidate != 0){
      vector<PseudoJet> kept_subjets0 = cms_top_candidate.structure_of<CMSTopTagger>().W().pieces();
      vector<PseudoJet> kept_subjets1 = cms_top_candidate.structure_of<CMSTopTagger>().non_W().pieces();
      vector<PseudoJet> all_subjets = kept_subjets0;
      all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

      PseudoJet lCorrCA = (*area_subtractor)(cms_top_candidate);
      //double lJEC = correction(cms_top_candidate,iJetCorr,bge_rho.rho());

      cmsttJetMass = cms_top_candidate.m();
      cmsttMinMass = cms_top_candidate.structure_of<CMSTopTagger>().W().m();
      cmsttHelicity = cms_top_candidate.structure_of<CMSTopTagger>().cos_theta_W();
      cmsttNsubjets = all_subjets.size();
      cmsttArea = cms_top_candidate.area();
      cmsttJetMassCorr = lCorrCA.m();
      cmsttMinMassCorr = cmsttMinMass;
  }

  // -- fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());  
  (iJetI.ptcorr    ).push_back(iJet      .pt());
  (iJetI.ptraw     ).push_back(iJet      .pt());
  (iJetI.ptunc     ).push_back(0.);
  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());

  (iJetI.pttrim    ).push_back(lTrim     .pt());
  (iJetI.mtrim     ).push_back(lTrim     .m());

  (iJetI.pttrimsafe).push_back(lTrimSafe .pt());
  (iJetI.mtrimsafe ).push_back(lTrimSafe .m());

  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.mconst    ).push_back(lConstit  .m());

  (iJetI.ptclean   ).push_back(lClean    .pt());
  (iJetI.mclean    ).push_back(lClean    .m());

  (iJetI.mpruned   ).push_back(lPruned   .m());
  (iJetI.ptpruned  ).push_back(lPruned   .pt());

  (iJetI.mprunedsafe).push_back(lPrunedSafe.m());
  (iJetI.ptprunedsafe).push_back(lPrunedSafe.pt());

  (iJetI.ptsoftdrop).push_back(lSoftDropped.pt());
  (iJetI.msoftdrop).push_back(lSoftDropped.m());

  (iJetI.ptsoftdropsafe).push_back(lSoftDroppedSafe.pt());
  (iJetI.msoftdropsafe).push_back(lSoftDroppedSafe.m());

  (iJetI.nparticles).push_back((iJet.constituents()).size());
  (iJetI.nneutrals ).push_back(neutrals.size());
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());


  (iJetI.ptCA ).push_back(lCorrCA .pt());
  (iJetI.ptCAcorr ).push_back(iJetCA .pt());
  (iJetI.ptCAraw ).push_back(iJetCA .pt());
  (iJetI.msoftdropCA ).push_back(lSoftDroppedCA.m());
  (iJetI.msoftdropCAsafe).push_back(lSoftDroppedCASafe.m());

  (iJetI.sdsymmetry ).push_back( SoftDropedSymmetry );
  (iJetI.sddeltar ).push_back( SoftDropedDR );
  (iJetI.sdmu ).push_back( SoftDropedMassDrop );
  (iJetI.sdenergyloss ).push_back( SoftDropedEnergyLoss );
  (iJetI.sdarea ).push_back( SoftDropedArea );
  (iJetI.sdnconst ).push_back( SoftDropedNconst );
  (iJetI.mfiltsoftdropCA ).push_back( filtered_softdropped_jet.m() );


  (iJetI.hepmass ).push_back( hepttJetMass );
  (iJetI.hepwmass ).push_back( hepttWMass );
  (iJetI.hepm01 ).push_back( hepttM01 );
  (iJetI.hepm02 ).push_back( hepttM02 );
  (iJetI.hepm12 ).push_back( hepttM12 );
  (iJetI.hepm12m012 ).push_back( hepttM12M012 );
  (iJetI.hepatanm02m01).push_back( hepttAtanM02M01 );

  (iJetI.cmsmass ).push_back(cmsttJetMass );
  (iJetI.cmsminmass ).push_back(cmsttMinMass );
  (iJetI.cmshelicity ).push_back(cmsttHelicity );
  (iJetI.cmsnsubjets ).push_back(cmsttNsubjets );
  (iJetI.cmsarea ).push_back(cmsttArea );
  (iJetI.cmsmasscorr ).push_back(cmsttJetMassCorr );
  (iJetI.cmsminmasscorr ).push_back(cmsttMinMassCorr );

  vtagger.setInputJet(iJet); 

  (iJetI.tau1 ).push_back(vtagger.computeNSubJettines(1,1.,jetR,jetR));
  (iJetI.tau2 ).push_back(vtagger.computeNSubJettines(2,1.,jetR,jetR));
  (iJetI.tau3 ).push_back(vtagger.computeNSubJettines(3,1.,jetR,jetR));
  (iJetI.tau4 ).push_back(vtagger.computeNSubJettines(4,1.,jetR,jetR));
  (iJetI.tau5 ).push_back(vtagger.computeNSubJettines(5,1.,jetR,jetR));

  (iJetI.Qjets).push_back(vtagger.computeQjets(35,25,randNumber.Uniform(0.,10000000)));

  (iJetI.ecf_b05).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,0.5));
  (iJetI.ecf_b10).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,1.0));
  (iJetI.ecf_b15).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,1.5));
  (iJetI.ecf_b20).push_back(vtagger.computeECF(fastjet::antikt_algorithm,2.0,2,2.0));
  
  (iJetI.charge_k05).push_back(vtagger.computeJetCharge(0.5));
  (iJetI.charge_k07).push_back(vtagger.computeJetCharge(0.7));
  (iJetI.charge_k10).push_back(vtagger.computeJetCharge(1.0));
  
  vtagger.setInputJet(lPruned); 
  (iJetI.tau1_pr ).push_back(vtagger.computeNSubJettines(1,1.,jetR,jetR));
  (iJetI.tau2_pr ).push_back(vtagger.computeNSubJettines(2,1.,jetR,jetR));
  (iJetI.tau3_pr ).push_back(vtagger.computeNSubJettines(3,1.,jetR,jetR));
  (iJetI.tau4_pr ).push_back(vtagger.computeNSubJettines(4,1.,jetR,jetR));
  (iJetI.tau5_pr ).push_back(vtagger.computeNSubJettines(5,1.,jetR,jetR));

  vtagger.setInputJet(lSoftDropped); 
  (iJetI.tau1_softdrop ).push_back(vtagger.computeNSubJettines(1,1.,jetR,jetR));
  (iJetI.tau2_softdrop ).push_back(vtagger.computeNSubJettines(2,1.,jetR,jetR));
  (iJetI.tau3_softdrop ).push_back(vtagger.computeNSubJettines(3,1.,jetR,jetR));
  (iJetI.tau4_softdrop ).push_back(vtagger.computeNSubJettines(4,1.,jetR,jetR));
  (iJetI.tau5_softdrop ).push_back(vtagger.computeNSubJettines(5,1.,jetR,jetR));
    
}


// ------------------------------------------------------------------------------------------
void fillGenJetsInfo(vector<PseudoJet> &iJets, vector<PseudoJet> &iParticles, GenJetInfo &iJetInfo, JetCleanser &gsn_cleanser, int nPU){

  // -- Compute rho, rho_m for SafeAreaSubtraction
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4); // jet clustering for background estimation
  Selector rho_range =  SelectorAbsRapMax(5.0); // apply just an eta acceptance cut
  ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def); // cluster the initial particles with the k_t + active ghost for bkg determination
  JetMedianBackgroundEstimator  bge_rho(rho_range, clust_seq_rho);
  JetMedianBackgroundEstimator  bge_rhom(rho_range, clust_seq_rho); // get mediam background estimator
  BackgroundJetPtMDensity m_density;
  bge_rhom.set_jet_density_class(&m_density);
  
  //// use GridMedianBackgroundEstimator (faster), however doesn't work with SafeAreaSubtraction
  //GridMedianBackgroundEstimator bge_rho(5.0,0.8);
  //bge_rho.set_particles(iParticles);
  //GridMedianBackgroundEstimator bge_rhom(5.0,0.8);
  //bge_rhom.set_particles(iParticles);
  
  // -- Background estimator for constituents subtractor
  JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
  bge_rhoC.set_jet_density_class(scalarPtDensity);
  bge_rhoC.set_particles(iParticles);
    
  // -- Clear jet info for each event                                                                                                                                     
  clear(iJetInfo);  
  iJetInfo.npu = nPU;
  // -- Loop over jets in the event and set jets variables                                                                                                           
  for (unsigned int j = 0; j < iJets.size(); j++){
    setGenJet( iJets[j], iJetInfo,  bge_rho, bge_rhom, bge_rhoC, gsn_cleanser); // give the original clustered jets, the background estimations and cleansing
    //cout << iTree.GetName() << "  " << (iJetInfo.pt)[j] << "  "<< (iJetInfo.ptcorr)[j] <<endl;                                                                               
  }

}

// ------------------------------------------------------------------------------------------
void fillRecoJetsInfo(vector<PseudoJet> &iJets,  vector<PseudoJet> &iParticles, JetInfo &iJetInfo, GenJetInfo iGenJetInfo, bool isCHS, FactorizedJetCorrector *jetCorr, JetCorrectionUncertainty *ijetUnc, JetCleanser &gsn_cleanser, int nPU, vfloat eta_Boson, vfloat phi_Boson){
  
  // -- Compute rho, rho_m for SafeAreaSubtraction -> same procedure is used for GenJets
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  Selector rho_range =  SelectorAbsRapMax(5.0);
  ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def);
  JetMedianBackgroundEstimator bge_rho(rho_range, clust_seq_rho);
  JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_rho);
  BackgroundJetPtMDensity m_density;
  bge_rhom.set_jet_density_class(&m_density);
  
  //// use GridMedianBackgroundEstimator (faster), however doesn't work with SafeAreaSubtraction
  //GridMedianBackgroundEstimator bge_rho(5.0,0.8);
  //bge_rho.set_particles(iParticles);
  //GridMedianBackgroundEstimator bge_rhom(5.0,0.8);
  //bge_rhom.set_particles(iParticles);
  
  // -- Background estimator for constituents subtractor
  JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
  bge_rhoC.set_jet_density_class(scalarPtDensity);
  bge_rhoC.set_particles(iParticles);
    
  // -- Clear jet info for each event                                                                                                                                           
  clear(iJetInfo);  
  iJetInfo.npu = nPU;

  // -- Loop over jets in the event and set jets variables                                                                                                                      
  for (unsigned int j = 0; j < iJets.size(); j++){
    setRecoJet( iJets[j], iJetInfo, iGenJetInfo,bge_rho, bge_rhom, bge_rhoC, isCHS, jetCorr, ijetUnc, gsn_cleanser, eta_Boson, phi_Boson);
    //cout << iTree.GetName() << "  " << (iJetInfo.pt)[j] << "  "<< (iJetInfo.ptcorr)[j] <<endl;                                                                                   
  }
  
}

// -------- Method to setbranch address for TJet in the input file, which have been produced (clustered) inside CMSSW
void setupCMSSWJetReadOut(TTree *iTree, float R ) {  
  cout << "Setting up to read jet collection : " << Form("Jet0%d",int(R*10)) << endl;
  fJet  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress(Form("Jet0%d",int(R*10)), &fJet);
  fJetBr  = iTree->GetBranch(Form("Jet0%d",int(R*10)));
}


// ------------------------------------------------------------------------------------------
void readCMSSWJet(int entry, TTree *iTree, TTree &oTree,  std::vector<fastjet::PseudoJet> genJets, JetInfo &iJetI) {

  // -- Clear jet info for each event
  clear(iJetI);

  // -- Read event and fill jet info
  iTree->GetEntry(entry);

  for (int i = 0; i < fJet->GetEntriesFast(); i++){
    TJet *pJet = (TJet*)((*fJet)[i]);
    
    // -- fill jet info                                                                                                                                                                                                       
    (iJetI.pt        ).push_back(pJet->pt);
    (iJetI.ptcorr    ).push_back(pJet->pt);
    (iJetI.ptraw     ).push_back(pJet->ptRaw);
    (iJetI.eta       ).push_back(pJet->eta);
    (iJetI.phi       ).push_back(pJet->phi);
    (iJetI.m         ).push_back(pJet->mass);
    (iJetI.nparticles).push_back(pJet->nParticles);
    (iJetI.nneutrals ).push_back(pJet->nNeutrals);
    (iJetI.ncharged  ).push_back(pJet->nCharged);
    
    // for now fill this branches with dummy value
    (iJetI.mraw      ).push_back(-999.);
    (iJetI.ptunc     ).push_back(-999.);
    (iJetI.ptclean   ).push_back(-999.);
    (iJetI.mclean   ).push_back(-999.);
    (iJetI.pttrim    ).push_back(-999.);
    (iJetI.mtrim     ).push_back(-999.);
    (iJetI.pttrimsafe).push_back(-999.);
    (iJetI.mtrimsafe ).push_back(-999.);
    (iJetI.ptconst   ).push_back(-999.);
    (iJetI.mconst    ).push_back(-999.);
    (iJetI.ptpruned   ).push_back(-999.);
    (iJetI.mpruned   ).push_back(-999.);
    (iJetI.ptprunedsafe).push_back(-999.);
    (iJetI.mprunedsafe).push_back(-999.);
    (iJetI.ptsoftdrop).push_back(-999.);
    (iJetI.msoftdrop).push_back(-999.);
    (iJetI.ptsoftdropsafe).push_back(-999.);
    (iJetI.msoftdropsafe).push_back(-999.);
  
    //-- gen matching
    int imatch = -1;
    float mindr = dRMatching;
    TLorentzVector *recojet = new TLorentzVector();
    recojet->SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    for (unsigned int ig = 0; ig < genJets.size(); ig++){
      TLorentzVector *genjet = new TLorentzVector();
      genjet->SetPtEtaPhiM(genJets[ig].pt(), genJets[ig].eta(), genJets[ig].phi(), genJets[ig].m());
      double dr = recojet->DeltaR(*genjet);
      if (dr < mindr){
	mindr = dr;
	imatch = ig;
      }
      delete genjet;
    }
    
    delete recojet;

    if (imatch > -1){
      (iJetI.imatch   ).push_back(imatch);
      (iJetI.ptgen    ).push_back(genJets[imatch].pt());
      (iJetI.etagen   ).push_back(genJets[imatch].eta());
      (iJetI.phigen   ).push_back(genJets[imatch].phi());
      (iJetI.mgen     ).push_back(genJets[imatch].m());
      (iJetI.mrawgen     ).push_back(-999.);// dummy val
      (iJetI.mtrimgen    ).push_back(-999.);// dummy val
      (iJetI.mtrimsafegen).push_back(-999.);// dummy val
      (iJetI.mcleangen   ).push_back(-999.);// dummy val
      (iJetI.mconstgen   ).push_back(-999.);// dummy val
    }
    else {
      (iJetI.imatch   ).push_back(imatch);
      (iJetI.ptgen    ).push_back(-999.);
      (iJetI.etagen   ).push_back(-999.);
      (iJetI.phigen   ).push_back(-999.);
      (iJetI.mgen     ).push_back(-999.);
      (iJetI.mrawgen     ).push_back(-999.);// dummy val
      (iJetI.mtrimgen    ).push_back(-999.);// dummy val
      (iJetI.mtrimsafegen).push_back(-999.);// dummy val
      (iJetI.mcleangen   ).push_back(-999.);// dummy val
      (iJetI.mconstgen   ).push_back(-999.);// dummy val
    }
  }

  // --- fill tree 
  oTree.Fill();
  
}
// ------------------------------------------------------------------------------------------

fastjet::JetAlgorithm get_algo(string algo)
{
  fastjet::JetAlgorithm jetalgo;
  if (algo=="kt") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="ca") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="ak") jetalgo = fastjet::antikt_algorithm ;
  else if (algo=="KT") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="CA") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="AK") jetalgo = fastjet::antikt_algorithm ;
  else if (algo=="kt_algorithm") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="cambridge_algorithm") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="antikt_algorithm") jetalgo = fastjet::antikt_algorithm ;
  else if (algo=="0") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="1") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="2") jetalgo = fastjet::antikt_algorithm ;
  else jetalgo = fastjet::antikt_algorithm ;
  return jetalgo;
}


// function to fill a TChain with the list of input files to be processed 
bool FillChain(TChain& chain, const std::string& inputFileList){

  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open()){
      std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
      return false;
  }
  
  while(1){
      inFile >> buffer;
      if(!inFile.good()) break;
      chain.Add(buffer.c_str());
  }

  return true;
}


//---------------------------------------------------------------------------------------------------------------
//--- MAIN PROGRAM
//---------------------------------------------------------------------------------------------------------------

int main (int argc, char ** argv) {

  // --- args
  if (argc<3){
    cout << "Missing arguments!!!" <<endl;
    cout << "Usage: MiniNtuplizer <config> <input files list> <output file>" <<endl;
  }

  // args 
  std::string inputFilesList = argv[2]; // input file name
  std::string fOut           = argv[3]; // output name

  // --- Read configurable parameters from config                                                                                                                        
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");
  int maxEvents              = Options.getParameter<int>("maxEvents");        // max num of events to analyze
  double jetPtCut            = Options.getParameter<double>("jetPtCut"); //pT cut applied when getting jets from cluster sequence
 
  jetR                       = Options.getParameter<double>("jetR");          // jet cone size  
  std::string jetAlgo        = Options.getParameter<std::string>("jetAlgo"); // jet clustering algorithm
  fastjet::JetAlgorithm fatjet_algo = get_algo("jetAlgo");

  bool doCMSSWJets           = Options.getParameter<bool>("doCMSSWJets");     // analyze also default CMSSW PF jets
  std::string puppiConfig    = Options.getParameter<std::string>("puppiConfig"); // Puppi congiguration file

  std::string L1FastJetJEC    = Options.getParameter<std::string>("L1FastJetJEC");  // L1 JEC 
  std::string L2RelativeJEC   = Options.getParameter<std::string>("L2RelativeJEC"); // L2
  std::string L3AbsoluteJEC   = Options.getParameter<std::string>("L3AbsoluteJEC"); // L3
  std::string L2L3ResidualJEC = Options.getParameter<std::string>("L2L3ResidualJEC"); // L2L3 residual (for data only)
  std::string JECUncertainty  = Options.getParameter<std::string>("JECUncertainty"); // Uncertainty

  bool DoMatchingToBoson      = Options.getParameter<bool>("DoMatchingToBoson"); // this is relevant for the WW, ttbar etc. samples
  int pdgIdBoson              = Options.getParameter<int>("pdgIdBoson"); // absolute value of pdgId of the boson. Can be used only if the DoMatchingToBoson is set to true.
  dRMatching                  = Options.getParameter<double>("dRMatiching");   // dR matching thresholds with the truth
  
  //softdrop parameters
  beta         = Options.getParameter<double>("beta");
  symmetry_cut = Options.getParameter<double>("symmetry_cut");
  R0           = Options.getParameter<double>("R0");

  //trimming paramenters
  R_trimming = Options.getParameter<double>("R_trimming");
  PtFraction = Options.getParameter<double>("PtFraction");
  std::string trimAlgo = Options.getParameter<std::string>("trimAlgo"); // jet cone size
  algorithm_Trimming = get_algo("trimAlgo");

  
  //pruning paramenters
  z_cut = Options.getParameter<double>("z_cut");
  R_Cut = Options.getParameter<double>("R_Cut");
  std::string pruneAlgo = Options.getParameter<std::string>("pruneAlgo"); // jet cone size
  algorithm_Pruning = get_algo("pruneAlgo");
  R_jet_def_pruning = Options.getParameter<double>("R_jet_def_pruning");

  // --- Read list of files to be analyzed and fill TChain
  TChain* lTree = new TChain("Events");
  FillChain(*lTree, inputFilesList);
  if (lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 

  cout << "This analysis will run on "<< maxEvents << " events" <<endl; 

  // --- Load branches from the input tree -->  only the one related to gen particles and PFcandidates
  fPFCand = new PFLoader (lTree,puppiConfig.c_str());
  fGen    = new GenLoader(lTree);
  if (doCMSSWJets) setupCMSSWJetReadOut(lTree, jetR);

  TEventInfo *eventInfo = new TEventInfo();
  lTree->SetBranchAddress("Info",&eventInfo);

  // --- Setup JEC on the fly  
  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters(L1FastJetJEC.c_str()));  
  corrParams.push_back(JetCorrectorParameters(L2RelativeJEC.c_str()));  
  corrParams.push_back(JetCorrectorParameters(L3AbsoluteJEC.c_str()));  
  if (L2L3ResidualJEC!="") corrParams.push_back(JetCorrectorParameters(L2L3ResidualJEC.c_str())); // 
  JetCorrectorParameters param(JECUncertainty.c_str());      

      
  FactorizedJetCorrector   *jetCorr = new FactorizedJetCorrector(corrParams);
  JetCorrectionUncertainty *jetUnc  = new JetCorrectionUncertainty(param);
 

  // --- Setup JetAlgos for basic clustering of the event
  JetDefinition jet_def(fatjet_algo,jetR);
  //JetDefinition jet_def_Pruning(antikt_algorithm,0.3);//this is a jet algorithm for pruning. Smaller radius to be used
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0))); // real ghosts in the PseudoJet list 
  
  // --- Setup cleansing
  JetDefinition subjet_def(kt_algorithm,0.2);
  JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);

  // --- Setup soft-killer
  SoftKiller soft_killer (2.5,0.4);

  // --- Setup output trees -> one tree for each jet collection type: GenJets, PFJets, PFCHS, Puppi, cmssw and softkiller
  TFile *fout = new TFile(fOut.c_str(),"RECREATE");

  TTree *genTree   = new TTree("gen"  , "gen"  );
  TTree *pfTree    = new TTree("pf"   , "pf"   );
  TTree *chsTree   = new TTree("chs"  , "chs"  );
  TTree *puppiTree = new TTree("puppi", "puppi");
  TTree *softkillerTree    = new TTree("softkiller", "softkiller");
  TTree *cmsswTree = new TTree("cmsswpf", "cmsswpf");
  
  GenJetInfo JGenInfo;
  JetInfo JPFInfo, JCHSInfo, JPuppiInfo, JSoftKillerInfo, JCMSSWPFInfo; // declare structures to fill the output tree information + make branches

  setupGenTree(genTree,   JGenInfo    , "" );
  setupTree(pfTree,    JPFInfo     , "" );
  setupTree(chsTree,   JCHSInfo    , "" );
  setupTree(puppiTree, JPuppiInfo  , "" );
  setupTree(softkillerTree, JSoftKillerInfo  , "" );
  if (doCMSSWJets) setupTree(cmsswTree, JCMSSWPFInfo, "" );
   

  // --- start loop over events
  for(int ientry = 0; ientry < maxEvents; ientry++) { 

    // -- For each event build collections of particles (gen, puppi, etc..) to cluster as a first step
    Long64_t localEntry = lTree->LoadTree(ientry);
    fPFCand->load(localEntry); // load pF information
    fGen   ->load(localEntry); // load gen information  
    
    vector<PseudoJet> gen_event       = fGen   ->genFetch();  //gen particles: only status 1 (ME) and user_index set 2
    vector<PseudoJet> pf_event        = fPFCand->pfFetch();   //return all the particles
    vector<PseudoJet> chs_event       = fPFCand->pfchsFetch(-1); //only chs particles -> user_index set to 1(neutrals) or 2 (chaged from PV)
    vector<PseudoJet> puppi_event     = fPFCand->puppiFetch();   // puppi particles from all pf with puppi weights 
    vector<PseudoJet> soft_event      = soft_killer(pf_event);   //retun the list from soft_killer contructor given all pf and the input parameters


    // -- Cluster jets -> make the clustering
    ClusterSequenceArea pGen    (gen_event    , jet_def, area_def);
    ClusterSequenceArea pPup    (puppi_event  , jet_def, area_def);
    ClusterSequenceArea pPF     (pf_event     , jet_def, area_def);
    ClusterSequenceArea pCHS    (chs_event    , jet_def, area_def);
    ClusterSequenceArea pSoft   (soft_event   , jet_def, area_def);

    // -- Order in decreasing pt the final jet collection with an inclusive cut on jets of 25GeV
    vector<PseudoJet> genJets     = sorted_by_pt(pGen    .inclusive_jets(jetPtCut));
    vector<PseudoJet> puppiJets   = sorted_by_pt(pPup    .inclusive_jets(jetPtCut));
    vector<PseudoJet> pfJets      = sorted_by_pt(pPF     .inclusive_jets(jetPtCut));
    vector<PseudoJet> chsJets     = sorted_by_pt(pCHS    .inclusive_jets(jetPtCut));
    vector<PseudoJet> softJets    = sorted_by_pt(pSoft   .inclusive_jets(jetPtCut));

    lTree->GetEntry(ientry);
    int nPU = eventInfo->nPU;
      
    vfloat eta_Boson, phi_Boson; // vector of eta and phi of all the vector bosons at gen level
    
    if (DoMatchingToBoson){
      fGen -> selectBoson(pdgIdBoson);
      eta_Boson = fGen -> eta_Boson;
      phi_Boson = fGen -> phi_Boson;
    }


    cout << "\r" ;
    cout << "===> Processed " << ientry << " - Done : " << (float(ientry)/float(maxEvents))*100 << "%"  ;
    
    // save jet info in a tree
    fillGenJetsInfo(genJets, gen_event, JGenInfo, gsn_cleanser, nPU); // take as input clustered jets, intial  list of particles, output structure, cleansing and pileup
    fillRecoJetsInfo(puppiJets, puppi_event, JPuppiInfo, JGenInfo, false, jetCorr, jetUnc, gsn_cleanser,nPU, eta_Boson, phi_Boson);
    fillRecoJetsInfo(pfJets , pf_event   , JPFInfo   , JGenInfo, false, jetCorr, jetUnc, gsn_cleanser,nPU, fGen -> eta_Boson,fGen -> phi_Boson);
    fillRecoJetsInfo(chsJets,  chs_event  , JCHSInfo  , JGenInfo, true , jetCorr, jetUnc, gsn_cleanser,nPU, fGen -> eta_Boson,fGen -> phi_Boson);
    fillRecoJetsInfo(softJets, soft_event  , JSoftKillerInfo  , JGenInfo, true , jetCorr, jetUnc, gsn_cleanser,nPU, fGen -> eta_Boson,fGen -> phi_Boson);
       

    // fill the output trees
    genTree->Fill();
    puppiTree->Fill();
    pfTree->Fill();
    chsTree->Fill();
    softkillerTree->Fill();

    if (doCMSSWJets)
      readCMSSWJet(ientry, lTree, *cmsswTree, genJets, JCMSSWPFInfo);        
    
    fGen->reset();         
    fPFCand->reset();
  }


  // --- Write trees 
  fout->cd();
  genTree  ->Write();
  pfTree   ->Write();
  chsTree  ->Write();
  puppiTree->Write();
  softkillerTree->Write();
  if (doCMSSWJets)  cmsswTree->Write();
}  

 
 
