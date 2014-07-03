#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "../include/GenLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/PFLoader.hh"
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
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/Selector.hh"

#include "CMSTopTagger.hh"
#include <fastjet/tools/JHTopTagger.hh>
#include "SubstrAna/Summer14/src/QjetsPlugin.h"
#include "SubstrAna/Summer14/src/Qjets.h"
#include "fastjet/contrib/EnergyCorrelator.hh"
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

// Jet parameters
double jetR;

// SoftDrop parameters
double beta, symmetry_cut, R0;

//trimming parameters
double R_trimming, PtFraction;

//pruning parameters
double R_Cut, z_cut, R_jet_def_pruning;

fastjet::JetAlgorithm algorithm_Trimming, algorithm_Pruning;

struct JetInfo {
  int npu ;
  /////////////////////////pt /////
  vector<float> pt;
  vector<float> ptcorr;
  vector<float> ptraw;
  vector<float> ptclean;
  vector<float> pttrim;
  vector<float> pttrimsafe;
  vector<float> ptconst;
  vector<float> ptunc;
  vector<float> eta;
  vector<float> phi;
  vector<float> m;
  vector<float> mraw;
  vector<float> mclean;
  vector<float> mpruned;
  vector<float> mprunedsafe;
  vector<float> mtrim;
  vector<float> mtrimsafe;
  vector<float> msoftdrop;
  vector<float> msoftdropsafe;
  vector<float> mconst;
  vector<int>   nparticles;
  vector<int>   nneutrals;
  vector<int>   ncharged;

  vector<float> ptCA           ;
  vector<float> ptCAcorr       ;
  vector<float> ptCAraw        ;
  vector<float> msoftdropCA    ;
  vector<float> msoftdropCAsafe;

  vector<float> sdsymmetry      ; 
  vector<float> sddeltar        ; 
  vector<float> sdmu            ; 
  vector<float> sdenergyloss    ; 
  vector<float> sdarea          ; 
  vector<float> sdnconst        ; 
  vector<float> mfiltsoftdropCA ; 

  vector<float> qjetsvol;

  vector<float> tau1;
  vector<float> tau2;
  vector<float> tau3;
  vector<float> tau4;
  vector<float> tau5;

  vector<float> tau1pruned;
  vector<float> tau2pruned;
  vector<float> tau3pruned;
  vector<float> tau4pruned;
  vector<float> tau5pruned;

  vector<float> tau1softdrop;
  vector<float> tau2softdrop;
  vector<float> tau3softdrop;
  vector<float> tau4softdrop;
  vector<float> tau5softdrop;

  vector<float> hepmass      ;
  vector<float> hepmasscorr  ;
  vector<float> hepwmass     ;
  vector<float> hepm01       ;
  vector<float> hepm02       ;
  vector<float> hepm12       ;
  vector<float> hepm12m012   ;
  vector<float> hepatanm02m01;

  vector<float> cmsmass       ;
  vector<float> cmsminmass    ;
  vector<float> cmshelicity   ;
  vector<float> cmsnsubjets   ;
  vector<float> cmsarea       ;
  vector<float> cmsmasscorr   ;
  vector<float> cmsminmasscorr;

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
  
  vector<float> ptCAgen           ;
  vector<float> msoftdropCAgen    ;
  vector<float> msoftdropCAsafegen;
  vector<float> mfiltsoftdropCAgen;
  //matching to the Boson
  vector <bool> is_MatchedToBoson;
};


void getConstitsForCleansing(vector<PseudoJet> inputs, vector<PseudoJet> &oNeutrals, vector<PseudoJet> &oChargedLV, vector<PseudoJet> &oChargedPU){
  for (unsigned int i = 0; i < inputs.size(); i++){
    if (inputs[i].user_index() <= 1) oNeutrals.push_back(inputs[i]);
    if (inputs[i].user_index() == 2) oChargedLV.push_back(inputs[i]);
    if (inputs[i].user_index() == 3) oChargedPU.push_back(inputs[i]);
  }
}


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


//Q jets stuff
float FindRMS( std::vector< float > qjetmasses ){
    
    float total = 0.;
    float ctr = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        total = total + qjetmasses[i];
        ctr++;
    }
    float mean = total/ctr;
    
    float totalsquared = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
    }
    float RMS = sqrt( totalsquared/ctr );
    return RMS;
}

float FindMean( std::vector< float > qjetmasses ){
    float total = 0.;
    float ctr = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        total = total + qjetmasses[i];
        ctr++;
    }
    return total/ctr;
}
double getQjetVolatility(std::vector < fastjet::PseudoJet > constits, int QJetsN, int seed){
    
    std::vector< float > qjetmasses;
    
    double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.01);
    
    for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
        QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
        qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
        fastjet::JetDefinition qjet_def(&qjet_plugin);
        fastjet::ClusterSequence qjet_seq(constits, qjet_def);
        vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));
        
        if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }
        
    }
    
    // find RMS of a vector
    float qjetsRMS = FindRMS( qjetmasses );
    // find mean of a vector
    float qjetsMean = FindMean( qjetmasses );
    float qjetsVolatility = qjetsRMS/qjetsMean;
    return qjetsVolatility;
}
double run_qjets_get_vol (PseudoJet &iJet, int rand)
{
  int QJetsPreclustering = 35;
  std::vector<fastjet::PseudoJet> constits;
  unsigned int nqjetconstits = iJet.constituents().size();
  if (nqjetconstits < (unsigned int) QJetsPreclustering) constits = iJet.constituents();
  else constits = iJet.associated_cluster_sequence()->exclusive_subjets_up_to(iJet,QJetsPreclustering);
  double qjet_vol = getQjetVolatility(constits, 25, rand*25) ;
  constits.clear();
  return qjet_vol;
}



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


double unc( PseudoJet &iJet,JetCorrectionUncertainty *iJetUnc) { 
  if(fabs(iJet.eta()) > 5. || fabs(iJet.pt()) < 10.) return 1.;
  iJetUnc->setJetPt ( iJet.pt()  );
  iJetUnc->setJetEta( iJet.eta() );
  double jetunc = iJetUnc->getUncertainty(true);
  return jetunc;
}


int matchingIndex(PseudoJet jet, vector<PseudoJet> genjets) {
  float rmin = 9999.;  
  int imatch = -1;
  for(unsigned int i = 0; i < genjets.size(); i++) {
    float rtemp = jet.delta_R(genjets[i]);
    if ( rtemp > 0.3 ) continue;
    if ( rtemp < rmin ){
      rmin =  rtemp;
      imatch = i;
    }
  }
  return (imatch);  
}



int matchingIndexFromJetInfo(PseudoJet jet, JetInfo jetInfo) {
  float rmin = 9999.;  
  int imatch = -1;
  for(unsigned int i = 0; i < (jetInfo.eta).size(); i++) {
    double dEta = fabs(jet.eta() - (jetInfo.eta)[i]);
    double dPhi = fabs(jet.phi() - (jetInfo.phi)[i]);
    if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
    float rtemp = sqrt(dEta*dEta+dPhi*dPhi);
    if ( rtemp > 0.3 ) continue;
    if ( rtemp < rmin ){
      rmin =  rtemp;
      imatch = i;
    }
  }
  return (imatch);  
}


bool IsMatchedToGenBoson(vfloat eta, vfloat phi, PseudoJet Jet) {
 bool IsMatched=false;
   
  for (unsigned int iGen =0; iGen < eta.size(); ++iGen)
  {
      double dEta = fabs(eta.at(iGen) - (Jet.eta()));
      double dPhi = fabs(phi.at(iGen) - (Jet.phi()));
      if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
      float rtemp = sqrt(dEta*dEta+dPhi*dPhi);
      if ( rtemp < 0.3 ){
	IsMatched = true;
      }
  }
  
  return (IsMatched);  
}



void setupTree(TTree *iTree, JetInfo &iJet, std::string iName) {
  iTree->Branch((iName+"npu"       ).c_str(),&iJet.npu       );
  iTree->Branch((iName+"pt"        ).c_str(),&iJet.pt        );
  
  iTree->Branch((iName+"ptcorr"    ).c_str(),&iJet.ptcorr    );
  iTree->Branch((iName+"ptraw"     ).c_str(),&iJet.ptraw     );
  iTree->Branch((iName+"ptclean"   ).c_str(),&iJet.ptclean   );
  iTree->Branch((iName+"pttrim"    ).c_str(),&iJet.pttrim    );
  iTree->Branch((iName+"pttrimsafe").c_str(),&iJet.pttrimsafe);
  iTree->Branch((iName+"ptconst"   ).c_str(),&iJet.ptconst   );
  iTree->Branch((iName+"ptunc"     ).c_str(),&iJet.ptunc     );
  iTree->Branch((iName+"eta"       ).c_str(),&iJet.eta       );
  iTree->Branch((iName+"phi"       ).c_str(),&iJet.phi       );
  iTree->Branch((iName+"m"         ).c_str(),&iJet.m         );
  iTree->Branch((iName+"mraw"      ).c_str(),&iJet.mraw      );
  iTree->Branch((iName+"mtrim"     ).c_str(),&iJet.mtrim     );
  iTree->Branch((iName+"mtrimsafe" ).c_str(),&iJet.mtrimsafe );
  iTree->Branch((iName+"mpruned" ).c_str(),&iJet.mpruned );
  iTree->Branch((iName+"mprunedsafe" ).c_str(),&iJet.mprunedsafe );
  iTree->Branch((iName+"msoftdrop" ).c_str(),&iJet.msoftdrop);
  iTree->Branch((iName+"msoftdropsafe" ).c_str(),&iJet.msoftdropsafe);
  iTree->Branch((iName+"mclean"    ).c_str(),&iJet.mclean    );
  iTree->Branch((iName+"mconst"    ).c_str(),&iJet.mconst    );
  iTree->Branch((iName+"nparticles").c_str(),&iJet.nparticles);
  iTree->Branch((iName+"nneutrals" ).c_str(),&iJet.nneutrals);
  iTree->Branch((iName+"ncharged"  ).c_str(),&iJet.ncharged);

  iTree->Branch((iName+"ptCA"             ).c_str(),&iJet.ptCA                    );
  iTree->Branch((iName+"ptCAcorr"         ).c_str(),&iJet.ptCAcorr                );
  iTree->Branch((iName+"ptCAraw"          ).c_str(),&iJet.ptCAraw                 );
  iTree->Branch((iName+"msoftdropCA"      ).c_str(),&iJet.msoftdropCA             );
  iTree->Branch((iName+"msoftdropCAsafe"  ).c_str(),&iJet.msoftdropCAsafe         );

  iTree->Branch((iName+"sdsymmetry"      ).c_str(),&iJet.sdsymmetry      );
  iTree->Branch((iName+"sddeltar"        ).c_str(),&iJet.sddeltar        );
  iTree->Branch((iName+"sdmu"            ).c_str(),&iJet.sdmu            );
  iTree->Branch((iName+"sdenergyloss"    ).c_str(),&iJet.sdenergyloss    );
  iTree->Branch((iName+"sdarea"          ).c_str(),&iJet.sdarea          );
  iTree->Branch((iName+"sdnconst"        ).c_str(),&iJet.sdnconst        );
  iTree->Branch((iName+"mfiltsoftdropCA" ).c_str(),&iJet.mfiltsoftdropCA );

  iTree->Branch((iName+"qjetsvol"         ).c_str(),&iJet.qjetsvol         );

  iTree->Branch((iName+"tau1"         ).c_str(),&iJet.tau1         );
  iTree->Branch((iName+"tau2"         ).c_str(),&iJet.tau2         );
  iTree->Branch((iName+"tau3"         ).c_str(),&iJet.tau3         );
  iTree->Branch((iName+"tau4"         ).c_str(),&iJet.tau4         );
  iTree->Branch((iName+"tau5"         ).c_str(),&iJet.tau5         );

  iTree->Branch((iName+"tau1pruned"         ).c_str(),&iJet.tau1pruned         );
  iTree->Branch((iName+"tau2pruned"         ).c_str(),&iJet.tau2pruned         );
  iTree->Branch((iName+"tau3pruned"         ).c_str(),&iJet.tau3pruned         );
  iTree->Branch((iName+"tau4pruned"         ).c_str(),&iJet.tau4pruned         );
  iTree->Branch((iName+"tau5pruned"         ).c_str(),&iJet.tau5pruned         );

  iTree->Branch((iName+"tau1softdrop"         ).c_str(),&iJet.tau1softdrop         );
  iTree->Branch((iName+"tau2softdrop"         ).c_str(),&iJet.tau2softdrop         );
  iTree->Branch((iName+"tau3softdrop"         ).c_str(),&iJet.tau3softdrop         );
  iTree->Branch((iName+"tau4softdrop"         ).c_str(),&iJet.tau4softdrop         );
  iTree->Branch((iName+"tau5softdrop"         ).c_str(),&iJet.tau5softdrop         );

  iTree->Branch((iName+"hepmass"               ).c_str(),&iJet.hepmass               );
  iTree->Branch((iName+"hepmasscorr"           ).c_str(),&iJet.hepmasscorr           );
  iTree->Branch((iName+"hepwmass"              ).c_str(),&iJet.hepwmass              );
  iTree->Branch((iName+"hepm01"                ).c_str(),&iJet.hepm01                );
  iTree->Branch((iName+"hepm02"                ).c_str(),&iJet.hepm02                );
  iTree->Branch((iName+"hepm12"                ).c_str(),&iJet.hepm12                );
  iTree->Branch((iName+"hepm12m012"            ).c_str(),&iJet.hepm12m012            );
  iTree->Branch((iName+"hepatanm02m01"         ).c_str(),&iJet.hepatanm02m01         );

  iTree->Branch((iName+"cmsmass"        ).c_str(),&iJet.cmsmass        );
  iTree->Branch((iName+"cmsminmass"     ).c_str(),&iJet.cmsminmass     );
  iTree->Branch((iName+"cmshelicity"    ).c_str(),&iJet.cmshelicity    );
  iTree->Branch((iName+"cmsnsubjets"    ).c_str(),&iJet.cmsnsubjets    );
  iTree->Branch((iName+"cmsarea"        ).c_str(),&iJet.cmsarea        );
  iTree->Branch((iName+"cmsmasscorr"    ).c_str(),&iJet.cmsmasscorr    );
  iTree->Branch((iName+"cmsminmasscorr" ).c_str(),&iJet.cmsminmasscorr );

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
  
  iTree->Branch((iName+"ptCAgen"                 ).c_str(),&iJet.ptCAgen                 );
  iTree->Branch((iName+"msoftdropCAgen"          ).c_str(),&iJet.msoftdropCAgen          );
  iTree->Branch((iName+"msoftdropCAsafegen"      ).c_str(),&iJet.msoftdropCAsafegen      );
  iTree->Branch((iName+"mfiltsoftdropCAgen"      ).c_str(),&iJet.mfiltsoftdropCAgen      );




  
  //matched to the boson
  iTree->Branch((iName+"is_MatchedToBoson"      ).c_str(),&iJet.is_MatchedToBoson      );
  
}


void clear(JetInfo &iJet) {
  iJet.npu  = -1;

  iJet.pt             .clear();
  iJet.ptcorr         .clear();
  iJet.ptraw          .clear();
  iJet.ptclean        .clear();
  iJet.pttrim         .clear();
  iJet.pttrimsafe     .clear();
  iJet.eta            .clear();
  iJet.phi            .clear();
  iJet.m              .clear();
  iJet.mraw           .clear();
  iJet.mtrim          .clear();
  iJet.mtrimsafe      .clear();
  iJet.mpruned        .clear();
  iJet.mprunedsafe    .clear();
  iJet.msoftdrop      .clear();
  iJet.msoftdropsafe  .clear();
  iJet.mclean         .clear();
  iJet.mconst         .clear();
  iJet.nparticles     .clear();
  iJet.nneutrals      .clear();
  iJet.ncharged       .clear();
  iJet.ptCA           .clear();
  iJet.ptCAcorr       .clear();
  iJet.ptCAraw        .clear();
  iJet.msoftdropCA    .clear();
  iJet.msoftdropCAsafe.clear();
  iJet.sdsymmetry     .clear(); 
  iJet.sddeltar       .clear(); 
  iJet.sdmu           .clear(); 
  iJet.sdenergyloss   .clear(); 
  iJet.sdarea         .clear(); 
  iJet.sdnconst       .clear(); 
  iJet.mfiltsoftdropCA.clear(); 
  iJet.qjetsvol       .clear();
  iJet.tau1           .clear();
  iJet.tau2           .clear();
  iJet.tau3           .clear();
  iJet.tau4           .clear();
  iJet.tau5           .clear();
  iJet.tau1pruned     .clear();
  iJet.tau2pruned     .clear();
  iJet.tau3pruned     .clear();
  iJet.tau4pruned     .clear();
  iJet.tau5pruned     .clear();
  iJet.tau1softdrop   .clear();
  iJet.tau2softdrop   .clear();
  iJet.tau3softdrop   .clear();
  iJet.tau4softdrop   .clear();
  iJet.tau5softdrop   .clear();
  iJet.hepmass        .clear();
  iJet.hepmasscorr    .clear();
  iJet.hepwmass       .clear();
  iJet.hepm01         .clear();
  iJet.hepm02         .clear();
  iJet.hepm12         .clear();
  iJet.hepm12m012     .clear();
  iJet.hepatanm02m01  .clear();
  iJet.cmsmass        .clear();
  iJet.cmsminmass     .clear();
  iJet.cmshelicity    .clear();
  iJet.cmsnsubjets    .clear();
  iJet.cmsarea        .clear();
  iJet.cmsmasscorr    .clear();
  iJet.cmsminmasscorr .clear();
  iJet.ptgen          .clear();
  iJet.etagen         .clear();
  iJet.phigen         .clear();
  iJet.mgen           .clear();
  iJet.mrawgen        .clear();
  iJet.mtrimgen       .clear();
  iJet.mtrimsafegen   .clear();
  iJet.mcleangen      .clear();
  iJet.mconstgen      .clear();
  iJet.imatch         .clear();
  iJet.ptCAgen           .clear();
  iJet.msoftdropCAgen    .clear();
  iJet.msoftdropCAsafegen.clear();
  iJet.mfiltsoftdropCAgen.clear();
  iJet.is_MatchedToBoson .clear();

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
  // ALternate: PseudoJet lCorr =  (*area_subtractor)(lTrim);

   //pruning
  JetDefinition jet_def_Pruning(algorithm_Pruning, R_jet_def_pruning);
  Pruner pruner(jet_def_Pruning, z_cut,R_Cut);
  PseudoJet lPruned = pruner(iJet);
  PseudoJet lPrunedSafe = pruner(lCorr);
  //PseudoJet lCorr =  (*area_subtractor)();

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
  (iJetI.pt           ).push_back(lCorr     .pt());
  (iJetI.ptcorr       ).push_back(iJet      .pt()*lJEC);
  (iJetI.ptraw        ).push_back(iJet      .pt());
  (iJetI.ptclean      ).push_back(lClean    .pt());
  (iJetI.pttrim       ).push_back(lTrim     .pt());
  (iJetI.pttrimsafe   ).push_back(lTrimSafe .pt());
  (iJetI.ptconst      ).push_back(lConstit  .pt());
  (iJetI.ptunc        ).push_back(lUnc);
  (iJetI.eta          ).push_back(iJet      .eta());
  (iJetI.phi          ).push_back(iJet      .phi());
  (iJetI.mraw         ).push_back(iJet      .m());
  (iJetI.m            ).push_back(lCorr     .m());
  (iJetI.mclean       ).push_back(lClean    .m());
  (iJetI.mtrim        ).push_back(lTrim     .m());
  (iJetI.mtrimsafe    ).push_back(lTrimSafe .m());
  (iJetI.mpruned      ).push_back(lPruned   .m());
  (iJetI.mprunedsafe  ).push_back(lPrunedSafe.m());
  (iJetI.msoftdrop    ).push_back(lSoftDropped.m());
  (iJetI.msoftdropsafe).push_back(lSoftDroppedSafe.m());
  (iJetI.mconst       ).push_back(lConstit  .m());
  (iJetI.nparticles   ).push_back((iJet.constituents()).size());
  (iJetI.nneutrals    ).push_back(neutrals.size());
  (iJetI.ncharged     ).push_back(chargedLV.size()+chargedPU.size());

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


void setRecoJet(PseudoJet &iJet, JetInfo &iJetI, JetInfo iGenJetI, JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, bool isCHS, FactorizedJetCorrector *iJetCorr, JetCorrectionUncertainty *iJetUnc, JetCleanser &gsn_cleanser, vfloat eta_Boson, vfloat phi_Boson) {

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
 
  // -- pruning
  JetDefinition jet_def_Pruning(algorithm_Pruning, R_jet_def_pruning);
  Pruner pruner(jet_def_Pruning, z_cut, R_Cut);
  PseudoJet lPruned = pruner(iJet);
  PseudoJet lPrunedSafe = pruner(lCorr);

  // -- softdrop
  contrib::SoftDrop softdrop(beta, symmetry_cut, R0);
  PseudoJet lSoftDropped = softdrop(iJet);
  softdrop.set_subtractor(area_subtractor);
  PseudoJet lSoftDroppedSafe = softdrop(iJet);

  // -- recluster jet CA
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  JetDefinition jet_def_CA (fastjet::cambridge_algorithm, jetR);  //MAKE infiintie
  fastjet::ClusterSequenceArea cs_Recluster ( iJet.constituents(), jet_def_CA, area_def); 
  vector<fastjet::PseudoJet> jets_Recluster  = sorted_by_pt(cs_Recluster .inclusive_jets()); 
  fastjet::PseudoJet  iJetCA =  jets_Recluster [0];
  PseudoJet lCorrCA =  (*area_subtractor)(iJetCA);

  // -- softdrop reclustered jet
  contrib::SoftDrop softdropCA(beta, symmetry_cut, R0);  
  PseudoJet lSoftDroppedCA = softdropCA(iJetCA);         
  softdropCA.set_subtractor(area_subtractor);            
  PseudoJet lSoftDroppedCASafe = softdrop(iJetCA);       
 
  double  SoftDropedSymmetry   = -1.0;
  double  SoftDropedDR         = -1.0;
  double  SoftDropedMassDrop   = -1.0;
  double  SoftDropedEnergyLoss = -1.0;
  double  SoftDropedArea       = -1.0;
  double  SoftDropedNconst     = -1.0;
  PseudoJet filtered_softdropped_jet; 

  if (lSoftDroppedCASafe!=0 and lSoftDroppedCASafe.m()>0.0)
  {
    SoftDropedSymmetry   = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().symmetry(); 
    SoftDropedDR         = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().delta_R();  
    SoftDropedMassDrop   = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().mu();       
    SoftDropedEnergyLoss = 1-lSoftDroppedCASafe.pt()/iJetCA.pt();                           
    SoftDropedArea       = lSoftDroppedCASafe .area() ;                                     
    SoftDropedNconst     = lSoftDroppedCASafe .constituents().size() ;                      

    // filter jet dynamically based on deltaR between subjets (arXiv:0802.2470)
    double dyn_Rfilt = min(0.3, SoftDropedDR*0.5);
    int    dyn_nfilt = 3;
    Filter filtersoft(dyn_Rfilt, SelectorNHardest(dyn_nfilt));
    filtered_softdropped_jet = filtersoft(lSoftDroppedCASafe);
  }


  // -- apply the JEC
  double lJEC = correction(iJet,iJetCorr,bge_rho.rho());  
  double lUnc = unc       (iJet,iJetUnc);
  double lJEC_CA = correction(iJetCA,iJetCorr,bge_rho.rho());  
  //double lUnc_CA = unc       (iJetCA,iJetUnc);

  // -- N-subjettiness
  double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
  //double R0 = 0.8; // Characteristic jet radius for normalization
  double Rcut = 10000.0; // maximum R particles can be from axis to be included in jet (large value for no cutoff)   
  NsubParameters paraNsub(beta, jetR, Rcut);
  Njettiness nSubKT      (Njettiness::kt_axes,paraNsub);
  Njettiness nSubMin     (Njettiness::min_axes,paraNsub);
  Njettiness nSubOnePass (Njettiness::onepass_kt_axes,paraNsub);

  vector<PseudoJet> jet_constituents  =  iJet.constituents();      
  double nSubOnePassTau1  = nSubOnePass.getTau(1,jet_constituents );  
  double nSubOnePassTau2  = nSubOnePass.getTau(2,jet_constituents );
  double nSubOnePassTau3  = nSubOnePass.getTau(3,jet_constituents );
  double nSubOnePassTau4  = nSubOnePass.getTau(4,jet_constituents );
  double nSubOnePassTau5  = nSubOnePass.getTau(5,jet_constituents );
  
  // -- pruned N-subjettiness
  double   nSubOnePassTau1Pruned = 9999;
  double   nSubOnePassTau2Pruned = 9999;
  double   nSubOnePassTau3Pruned = 9999;
  double   nSubOnePassTau4Pruned = 9999;
  double   nSubOnePassTau5Pruned = 9999;

  if (lPruned!=0)
  {
    vector<PseudoJet> jet_constituents_pruned  =  lPruned.constituents();      
    nSubOnePassTau1Pruned  = nSubOnePass.getTau(1,jet_constituents_pruned );  
    nSubOnePassTau2Pruned  = nSubOnePass.getTau(2,jet_constituents_pruned );
    nSubOnePassTau3Pruned  = nSubOnePass.getTau(3,jet_constituents_pruned );
    nSubOnePassTau4Pruned  = nSubOnePass.getTau(4,jet_constituents_pruned );
    nSubOnePassTau5Pruned  = nSubOnePass.getTau(5,jet_constituents_pruned );
  }
  
  // -- soft drop CA N-subjettiness
  double   nSubOnePassTau1SoftDrop = 9999;
  double   nSubOnePassTau2SoftDrop = 9999;
  double   nSubOnePassTau3SoftDrop = 9999;
  double   nSubOnePassTau4SoftDrop = 9999;
  double   nSubOnePassTau5SoftDrop = 9999;

  if (lSoftDroppedCA!=0)
  {
    vector<PseudoJet> jet_constituents_softdrop  =  lSoftDroppedCA.constituents();      
    nSubOnePassTau1SoftDrop  = nSubOnePass.getTau(1,jet_constituents_softdrop );  
    nSubOnePassTau2SoftDrop  = nSubOnePass.getTau(2,jet_constituents_softdrop );
    nSubOnePassTau3SoftDrop  = nSubOnePass.getTau(3,jet_constituents_softdrop );
    nSubOnePassTau4SoftDrop  = nSubOnePass.getTau(4,jet_constituents_softdrop );
    nSubOnePassTau5SoftDrop  = nSubOnePass.getTau(5,jet_constituents_softdrop );
  }  

  // -- HEP Top Tagger 
  double mass_drop_threshold=0.8;
  double max_subjet_mass=30;
  bool use_subjet_mass_cuts=false;
  HEPTopTagger hep_top_tagger(mass_drop_threshold, max_subjet_mass, use_subjet_mass_cuts);
  
  PseudoJet hep_top_candidate   = hep_top_tagger( iJetCA );
  double hepttJetMass    = -1; 
  double hepttJetMassCorr= -1; 
  double hepttWMass      = -1; 
  double hepttM01        = -1; 
  double hepttM02        = -1; 
  double hepttM12        = -1; 
  double hepttM12M012    = -1; 
  double hepttAtanM02M01 = -1; 
  if (hep_top_candidate != 0)
  {
    PseudoJet W =     hep_top_candidate.structure_of<HEPTopTagger>().W();
    PseudoJet W1 =    hep_top_candidate.structure_of<HEPTopTagger>().W1();
    PseudoJet W2 =    hep_top_candidate.structure_of<HEPTopTagger>().W2();
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

    hepttJetMass       = hep_top_candidate.m();
    hepttWMass         = W.m();
    hepttM01           = sum01.m();
    hepttM02           = sum02.m();
    hepttM12           = sum12.m();
    if ( sum012.m()!=0 ) hepttM12M012     = sum12.m() / sum012.m() ;
    if ( sum01.m()!=0 )  hepttAtanM02M01  = atan( sum02.m() / sum01.m() ) ;

    //PseudoJet lCorrHEP =  (*area_subtractor)(hep_top_candidate);
    hepttJetMassCorr       = hep_top_candidate.m();

  }


  // -- CMS Top Tagger 
  double cms_delta_p = 0.05;
  double cms_delta_r=0.4;
  double A=0.0004;

  CMSTopTagger cms_top_tagger(cms_delta_p, cms_delta_r, A);
  PseudoJet cms_top_candidate  = cms_top_tagger( iJetCA );

  double cmsttJetMass     = -1;
  double cmsttMinMass     = -1;
  double cmsttHelicity    = -1;
  double cmsttNsubjets    = -1;
  double cmsttArea        = -1;
  double cmsttJetMassCorr = -1;
  double cmsttMinMassCorr = -1;

  if (cms_top_candidate != 0)
  {
    vector<PseudoJet> kept_subjets0 = cms_top_candidate.structure_of<CMSTopTagger>().W().pieces();
    vector<PseudoJet> kept_subjets1 = cms_top_candidate.structure_of<CMSTopTagger>().non_W().pieces();
    vector<PseudoJet> all_subjets = kept_subjets0;
    all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

    PseudoJet lCorrCMS =  (*area_subtractor)(cms_top_candidate);

    cmsttJetMass      = cms_top_candidate.m();
    cmsttMinMass      = cms_top_candidate.structure_of<CMSTopTagger>().W().m();
    cmsttHelicity     = cms_top_candidate.structure_of<CMSTopTagger>().cos_theta_W();
    cmsttNsubjets     = all_subjets.size();
    cmsttArea         = cms_top_candidate.area();
    cmsttJetMassCorr  = lCorrCMS.m();
    cmsttMinMassCorr  = cmsttMinMass  * lJEC;
  } 

  // -- Qjets 
  double qjetVol  = run_qjets_get_vol( iJet, 2 ) ;  // FIXME seed
    

  // -- find the gen jet matched to this reco jet
  //int imatch = matchingIndex(iJet,genJets);
  int imatch = matchingIndexFromJetInfo(iJet, iGenJetI);
  bool matched = IsMatchedToGenBoson( eta_Boson, phi_Boson, iJet);
  //cout << matched << endl;

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

  (iJetI.ptCA           ).push_back(lCorrCA     .pt());
  (iJetI.ptCAcorr       ).push_back(iJetCA      .pt()*lJEC_CA);
  (iJetI.ptCAraw        ).push_back(iJetCA      .pt());
  (iJetI.msoftdropCA    ).push_back(lSoftDroppedCA.m());
  (iJetI.msoftdropCAsafe).push_back(lSoftDroppedCASafe.m());
             
  (iJetI.sdsymmetry       ).push_back( SoftDropedSymmetry           );
  (iJetI.sddeltar         ).push_back( SoftDropedDR                 );
  (iJetI.sdmu             ).push_back( SoftDropedMassDrop           );
  (iJetI.sdenergyloss     ).push_back( SoftDropedEnergyLoss         );
  (iJetI.sdarea           ).push_back( SoftDropedArea               );
  (iJetI.sdnconst         ).push_back( SoftDropedNconst             );
  (iJetI.mfiltsoftdropCA  ).push_back( filtered_softdropped_jet.m() );        

  (iJetI.qjetsvol).push_back( qjetVol );

  (iJetI.tau1).push_back( nSubOnePassTau1 );
  (iJetI.tau2).push_back( nSubOnePassTau2 );
  (iJetI.tau3).push_back( nSubOnePassTau3 );
  (iJetI.tau4).push_back( nSubOnePassTau4 );
  (iJetI.tau5).push_back( nSubOnePassTau5 );

  (iJetI.tau1pruned).push_back( nSubOnePassTau1Pruned );
  (iJetI.tau2pruned).push_back( nSubOnePassTau2Pruned );
  (iJetI.tau3pruned).push_back( nSubOnePassTau3Pruned );
  (iJetI.tau4pruned).push_back( nSubOnePassTau4Pruned );
  (iJetI.tau5pruned).push_back( nSubOnePassTau5Pruned );

  (iJetI.tau1softdrop).push_back( nSubOnePassTau1SoftDrop );
  (iJetI.tau2softdrop).push_back( nSubOnePassTau2SoftDrop );
  (iJetI.tau3softdrop).push_back( nSubOnePassTau3SoftDrop );
  (iJetI.tau4softdrop).push_back( nSubOnePassTau4SoftDrop );
  (iJetI.tau5softdrop).push_back( nSubOnePassTau5SoftDrop );

  (iJetI.hepmass      ).push_back( hepttJetMass    );
  (iJetI.hepmasscorr  ).push_back( hepttJetMassCorr);
  (iJetI.hepwmass     ).push_back( hepttWMass      );
  (iJetI.hepm01       ).push_back( hepttM01        );
  (iJetI.hepm02       ).push_back( hepttM02        );
  (iJetI.hepm12       ).push_back( hepttM12        );
  (iJetI.hepm12m012   ).push_back( hepttM12M012    );
  (iJetI.hepatanm02m01).push_back( hepttAtanM02M01 );

  (iJetI.cmsmass        ).push_back(cmsttJetMass     );
  (iJetI.cmsminmass     ).push_back(cmsttMinMass     );
  (iJetI.cmshelicity    ).push_back(cmsttHelicity    );
  (iJetI.cmsnsubjets    ).push_back(cmsttNsubjets    );
  (iJetI.cmsarea        ).push_back(cmsttArea        );
  (iJetI.cmsmasscorr    ).push_back(cmsttJetMassCorr );
  (iJetI.cmsminmasscorr ).push_back(cmsttMinMassCorr );

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

    (iJetI.ptCAgen               ).push_back((iGenJetI.ptCA)[imatch]);
    (iJetI.msoftdropCAgen        ).push_back((iGenJetI.msoftdropCA)[imatch]);
    (iJetI.msoftdropCAsafegen    ).push_back((iGenJetI.msoftdropCAsafe)[imatch]);
    (iJetI.mfiltsoftdropCAgen    ).push_back((iGenJetI.mfiltsoftdropCA)[imatch]);

  }
  else { //dummy values ?
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


    (iJetI.ptCAgen               ).push_back( -999.);
    (iJetI.msoftdropCAgen        ).push_back( -999.);
    (iJetI.msoftdropCAsafegen    ).push_back( -999.);
    (iJetI.mfiltsoftdropCAgen    ).push_back( -999.);

  }

}


void setGenJet(PseudoJet &iJet, JetInfo &iJetI,  JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, JetCleanser &gsn_cleanser) {

  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
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
  JetDefinition jet_def_Pruning(algorithm_Pruning, R_jet_def_pruning); 
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
  vector<fastjet::PseudoJet> jets_Recluster  = sorted_by_pt(cs_Recluster .inclusive_jets());
  fastjet::PseudoJet  iJetCA =  jets_Recluster [0];
  PseudoJet lCorrCA =  (*area_subtractor)(iJetCA);

  // -- softdrop reclustered jet
  contrib::SoftDrop softdropCA(beta, symmetry_cut, R0);
  PseudoJet lSoftDroppedCA = softdropCA(iJetCA);
  softdropCA.set_subtractor(area_subtractor);
  PseudoJet lSoftDroppedCASafe = softdrop(iJetCA);
 
  double  SoftDropedSymmetry   = -1.0;
  double  SoftDropedDR         = -1.0;
  double  SoftDropedMassDrop   = -1.0;
  double  SoftDropedEnergyLoss = -1.0;
  double  SoftDropedArea       = -1.0;
  double  SoftDropedNconst     = -1.0;
  PseudoJet filtered_softdropped_jet; 

  if (lSoftDroppedCASafe!=0 and lSoftDroppedCASafe.m()>0.0)
  {
    SoftDropedSymmetry   = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().symmetry(); 
    SoftDropedDR         = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().delta_R();  
    SoftDropedMassDrop   = lSoftDroppedCASafe.structure_of<contrib::SoftDrop>().mu();       
    SoftDropedEnergyLoss = 1-lSoftDroppedCASafe.pt()/iJetCA.pt();                           
    SoftDropedArea       = lSoftDroppedCASafe .area() ;                                     
    SoftDropedNconst     = lSoftDroppedCASafe .constituents().size() ;                      



    // filter jet dynamically based on deltaR between subjets (arXiv:0802.2470)
    double dyn_Rfilt = min(0.3, SoftDropedDR*0.5);
    int    dyn_nfilt = 3;
    Filter filtersoft(dyn_Rfilt, SelectorNHardest(dyn_nfilt));
    filtered_softdropped_jet = filtersoft(lSoftDroppedCASafe);
  }



  // -- N-subjettiness
    //N-subjettiness parameters
  double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
  //double R0 = 0.8; // Characteristic jet radius for normalization
  double Rcut = 10000.0; // maximum R particles can be from axis to be included in jet (large value for no cutoff)   
  NsubParameters paraNsub(beta, jetR, Rcut);
  Njettiness nSubKT      (Njettiness::kt_axes,paraNsub);
  Njettiness nSubMin     (Njettiness::min_axes,paraNsub);
  Njettiness nSubOnePass (Njettiness::onepass_kt_axes,paraNsub);

  vector<PseudoJet> jet_constituents  =  iJet.constituents();      
  double nSubOnePassTau1  = nSubOnePass.getTau(1,jet_constituents );  
  double nSubOnePassTau2  = nSubOnePass.getTau(2,jet_constituents );
  double nSubOnePassTau3  = nSubOnePass.getTau(3,jet_constituents );
  double nSubOnePassTau4  = nSubOnePass.getTau(4,jet_constituents );
  double nSubOnePassTau5  = nSubOnePass.getTau(5,jet_constituents );

  // ADD SHAPE SUTBRACTION

  
  // -- pruned N-subjettiness
  double   nSubOnePassTau1Pruned = 9999;
  double   nSubOnePassTau2Pruned = 9999;
  double   nSubOnePassTau3Pruned = 9999;
  double   nSubOnePassTau4Pruned = 9999;
  double   nSubOnePassTau5Pruned = 9999;

  if (lPruned!=0){
    vector<PseudoJet> jet_constituents_pruned  =  lPruned.constituents();      
    nSubOnePassTau1Pruned  = nSubOnePass.getTau(1,jet_constituents_pruned );  
    nSubOnePassTau2Pruned  = nSubOnePass.getTau(2,jet_constituents_pruned );
    nSubOnePassTau3Pruned  = nSubOnePass.getTau(3,jet_constituents_pruned );
    nSubOnePassTau4Pruned  = nSubOnePass.getTau(4,jet_constituents_pruned );
    nSubOnePassTau5Pruned  = nSubOnePass.getTau(5,jet_constituents_pruned );
  }

  // -- soft drop CA N-subjettiness
  double   nSubOnePassTau1SoftDrop = 9999;
  double   nSubOnePassTau2SoftDrop = 9999;
  double   nSubOnePassTau3SoftDrop = 9999;
  double   nSubOnePassTau4SoftDrop = 9999;
  double   nSubOnePassTau5SoftDrop = 9999;

  if (lSoftDroppedCA!=0){
    vector<PseudoJet> jet_constituents_softdrop  =  lSoftDroppedCA.constituents();      
    nSubOnePassTau1SoftDrop  = nSubOnePass.getTau(1,jet_constituents_softdrop );  
    nSubOnePassTau2SoftDrop  = nSubOnePass.getTau(2,jet_constituents_softdrop );
    nSubOnePassTau3SoftDrop  = nSubOnePass.getTau(3,jet_constituents_softdrop );
    nSubOnePassTau4SoftDrop  = nSubOnePass.getTau(4,jet_constituents_softdrop );
    nSubOnePassTau5SoftDrop  = nSubOnePass.getTau(5,jet_constituents_softdrop );
  }  


  // -- HEP Top Tagger 
  double mass_drop_threshold=0.8;
  double max_subjet_mass=30;
  bool use_subjet_mass_cuts=false;
  HEPTopTagger hep_top_tagger(mass_drop_threshold, max_subjet_mass, use_subjet_mass_cuts);
  
  PseudoJet hep_top_candidate   = hep_top_tagger( iJet );
  double hepttJetMass    = -1; 
  double hepttWMass      = -1; 
  double hepttM01        = -1; 
  double hepttM02        = -1; 
  double hepttM12        = -1; 
  double hepttM12M012    = -1; 
  double hepttAtanM02M01 = -1; 
  if (hep_top_candidate != 0)
  {
    PseudoJet W =     hep_top_candidate.structure_of<HEPTopTagger>().W();
    PseudoJet W1 =    hep_top_candidate.structure_of<HEPTopTagger>().W1();
    PseudoJet W2 =    hep_top_candidate.structure_of<HEPTopTagger>().W2();
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

    hepttJetMass       = hep_top_candidate.m();
    hepttWMass         = W.m();
    hepttM01           = sum01.m();
    hepttM02           = sum02.m();
    hepttM12           = sum12.m();
    if ( sum012.m()!=0 ) hepttM12M012     = sum12.m() / sum012.m() ;
    if ( sum01.m()!=0 )  hepttAtanM02M01  = atan( sum02.m() / sum01.m() ) ;
  }


  // -- CMS Top Tagger 
  double cms_delta_p = 0.05;
  double cms_delta_r=0.4;
  double A=0.0004;

  CMSTopTagger cms_top_tagger(cms_delta_p, cms_delta_r, A);
  PseudoJet cms_top_candidate  = cms_top_tagger( iJetCA );

  double cmsttJetMass     = -1;
  double cmsttMinMass     = -1;
  double cmsttHelicity    = -1;
  double cmsttNsubjets    = -1;
  double cmsttArea        = -1;
  double cmsttJetMassCorr = -1;
  double cmsttMinMassCorr = -1;

  if (cms_top_candidate != 0)
  {
    vector<PseudoJet> kept_subjets0 = cms_top_candidate.structure_of<CMSTopTagger>().W().pieces();
    vector<PseudoJet> kept_subjets1 = cms_top_candidate.structure_of<CMSTopTagger>().non_W().pieces();
    vector<PseudoJet> all_subjets = kept_subjets0;
    all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

    PseudoJet lCorrCA =  (*area_subtractor)(cms_top_candidate);
    //double lJEC = correction(cms_top_candidate,iJetCorr,bge_rho.rho());  

    cmsttJetMass      = cms_top_candidate.m();
    cmsttMinMass      = cms_top_candidate.structure_of<CMSTopTagger>().W().m();
    cmsttHelicity     = cms_top_candidate.structure_of<CMSTopTagger>().cos_theta_W();
    cmsttNsubjets     = all_subjets.size();
    cmsttArea         = cms_top_candidate.area();
    cmsttJetMassCorr  = lCorrCA.m();
    cmsttMinMassCorr  = cmsttMinMass;
  } 


  // -- Qjets 
  double qjetVol  =  run_qjets_get_vol( iJet, 2 ) ;  // FIXME seed
    

  // -- fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());
  

  (iJetI.ptcorr    ).push_back(iJet      .pt());
  (iJetI.ptraw     ).push_back(iJet      .pt());
  (iJetI.ptclean   ).push_back(lClean    .pt());
  (iJetI.pttrim    ).push_back(lTrim     .pt());
  (iJetI.pttrimsafe).push_back(lTrimSafe .pt());
  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.ptunc     ).push_back(0.);
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
    
  (iJetI.ptCA           ).push_back(lCorrCA     .pt());
  (iJetI.ptCAcorr       ).push_back(iJetCA      .pt());
  (iJetI.ptCAraw        ).push_back(iJetCA      .pt());
  (iJetI.msoftdropCA    ).push_back(lSoftDroppedCA.m());
  (iJetI.msoftdropCAsafe).push_back(lSoftDroppedCASafe.m());

  (iJetI.sdsymmetry       ).push_back( SoftDropedSymmetry           );
  (iJetI.sddeltar         ).push_back( SoftDropedDR                 );
  (iJetI.sdmu             ).push_back( SoftDropedMassDrop           );
  (iJetI.sdenergyloss     ).push_back( SoftDropedEnergyLoss         );
  (iJetI.sdarea           ).push_back( SoftDropedArea               );
  (iJetI.sdnconst         ).push_back( SoftDropedNconst             );
  (iJetI.mfiltsoftdropCA  ).push_back( filtered_softdropped_jet.m() );        

  (iJetI.qjetsvol).push_back( qjetVol );

  (iJetI.tau1).push_back( nSubOnePassTau1 );
  (iJetI.tau2).push_back( nSubOnePassTau2 );
  (iJetI.tau3).push_back( nSubOnePassTau3 );
  (iJetI.tau4).push_back( nSubOnePassTau4 );
  (iJetI.tau5).push_back( nSubOnePassTau5 );

  (iJetI.tau1pruned).push_back( nSubOnePassTau1Pruned );
  (iJetI.tau2pruned).push_back( nSubOnePassTau2Pruned );
  (iJetI.tau3pruned).push_back( nSubOnePassTau3Pruned );
  (iJetI.tau4pruned).push_back( nSubOnePassTau4Pruned );
  (iJetI.tau5pruned).push_back( nSubOnePassTau5Pruned );

  (iJetI.tau1softdrop).push_back( nSubOnePassTau1SoftDrop );
  (iJetI.tau2softdrop).push_back( nSubOnePassTau2SoftDrop );
  (iJetI.tau3softdrop).push_back( nSubOnePassTau3SoftDrop );
  (iJetI.tau4softdrop).push_back( nSubOnePassTau4SoftDrop );
  (iJetI.tau5softdrop).push_back( nSubOnePassTau5SoftDrop );

  (iJetI.hepmass      ).push_back( hepttJetMass    );
  (iJetI.hepwmass     ).push_back( hepttWMass      );
  (iJetI.hepm01       ).push_back( hepttM01        );
  (iJetI.hepm02       ).push_back( hepttM02        );
  (iJetI.hepm12       ).push_back( hepttM12        );
  (iJetI.hepm12m012   ).push_back( hepttM12M012    );
  (iJetI.hepatanm02m01).push_back( hepttAtanM02M01 );

  (iJetI.cmsmass        ).push_back(cmsttJetMass     );
  (iJetI.cmsminmass     ).push_back(cmsttMinMass     );
  (iJetI.cmshelicity    ).push_back(cmsttHelicity    );
  (iJetI.cmsnsubjets    ).push_back(cmsttNsubjets    );
  (iJetI.cmsarea        ).push_back(cmsttArea        );
  (iJetI.cmsmasscorr    ).push_back(cmsttJetMassCorr );
  (iJetI.cmsminmasscorr ).push_back(cmsttMinMassCorr );

  // some duplicate infos.... should be fixed.
  (iJetI.imatch      ).push_back((iJetI.pt).size()-1);
  (iJetI.ptgen       ).push_back(lCorr     .pt());
  (iJetI.etagen      ).push_back(iJet      .eta());
  (iJetI.phigen      ).push_back(iJet      .phi());
  (iJetI.mgen        ).push_back(lCorr     .m());
  (iJetI.mrawgen     ).push_back(iJet      .m());
  (iJetI.mtrimgen    ).push_back(lTrim     .m());
  (iJetI.mtrimsafegen).push_back(lTrimSafe .m());
  (iJetI.mcleangen   ).push_back(lClean    .m());
  (iJetI.mconstgen   ).push_back(lConstit  .m());

  (iJetI.ptCAgen             ).push_back(iJetCA      .pt() );
  (iJetI.msoftdropCAgen      ).push_back(lSoftDroppedCA.m());
  (iJetI.msoftdropCAsafegen  ).push_back(lSoftDroppedCASafe.m());
  (iJetI.mfiltsoftdropCAgen  ).push_back(filtered_softdropped_jet.m());





}


// ------------------------------------------------------------------------------------------
void fillGenJetsInfo(vector<PseudoJet> &iJets, vector<PseudoJet> &iParticles, JetInfo &iJetInfo, JetCleanser &gsn_cleanser, int nPU ){

  // -- Compute rho, rho_m for SafeAreaSubtraction
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
    setGenJet( iJets[j], iJetInfo,  bge_rho, bge_rhom, bge_rhoC, gsn_cleanser);
    //cout << iTree.GetName() << "  " << (iJetInfo.pt)[j] << "  "<< (iJetInfo.ptcorr)[j] <<endl;                                                                                                                                      
  }

}
// ------------------------------------------------------------------------------------------


void fillRecoJetsInfo(vector<PseudoJet> &iJets,  vector<PseudoJet> &iParticles, JetInfo &iJetInfo, JetInfo iGenJetInfo, bool isCHS, FactorizedJetCorrector *jetCorr, JetCorrectionUncertainty *ijetUnc,
		      JetCleanser &gsn_cleanser, int nPU, vfloat eta_Boson, vfloat phi_Boson ){
  
  // -- Compute rho, rho_m for SafeAreaSubtraction
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

// ------------------------------------------------------------------------------------------
void setupCMSSWJetReadOut(TTree *iTree, float R ) {
  
  cout << "Setting up to read jet collection : " << Form("Jet0%d",int(R*10)) << endl;
  fJet  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress(Form("Jet0%d",int(R*10)), &fJet);
  fJetBr  = iTree->GetBranch(Form("Jet0%d",int(R*10)));
}
// ------------------------------------------------------------------------------------------


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
    (iJetI.ptclean   ).push_back(-999.);
    (iJetI.pttrim    ).push_back(-999.);
    (iJetI.pttrimsafe).push_back(-999.);
    (iJetI.ptconst   ).push_back(-999.);
    (iJetI.ptunc     ).push_back(-999.);
    (iJetI.mraw      ).push_back(-999.);
    (iJetI.mclean    ).push_back(-999.);
    (iJetI.mtrim     ).push_back(-999.);
    (iJetI.mtrimsafe ).push_back(-999.);
    (iJetI.mpruned   ).push_back(-999.);
    (iJetI.mprunedsafe).push_back(-999.);
    (iJetI.msoftdrop).push_back(-999.);
    (iJetI.msoftdropsafe).push_back(-999.);
    (iJetI.mconst    ).push_back(-999.);
  
    //-- gen matching
    int imatch = -1;
    double mindr = 0.3;
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



// ------------------------------------------------------------------------------------------
bool FillChain(TChain& chain, const std::string& inputFileList)
{
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open())
    {
      std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
      return false;
    }
  
  while(1)
    {
      inFile >> buffer;
      if(!inFile.good()) break;
      chain.Add(buffer.c_str());
    }

  return true;
}
// ------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------------
//--- MAIN PROGRAM
//---------------------------------------------------------------------------------------------------------------
int main (int argc, char ** argv) {

  cout << "dfgasf"  << endl;
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
  int maxEvents              = Options.getParameter<int>("maxEvents"); // max num of events to analyze
  jetR                = Options.getParameter<double>("jetR"); // jet cone size  
  bool doCMSSWJets           = Options.getParameter<bool>("doCMSSWJets"); // analyze also default CMSSW PF jets
  std::string puppiConfig    = Options.getParameter<std::string>("puppiConfig"); // Puppi congiguration file

  std::string L1FastJetJEC   = Options.getParameter<std::string>("L1FastJetJEC");  // L1 JEC 
  std::string L2RelativeJEC  = Options.getParameter<std::string>("L2RelativeJEC"); // L2
  std::string L3AbsoluteJEC  = Options.getParameter<std::string>("L3AbsoluteJEC"); // L3
  std::string L2L3ResidualJEC= Options.getParameter<std::string>("L2L3ResidualJEC"); // L2L3 residual (for data only)
  std::string JECUncertainty = Options.getParameter<std::string>("JECUncertainty"); // Uncertainty
  bool DoMatchingToBoson =  Options.getParameter<bool>("DoMatchingToBoson"); // this is relevant for the WW, ttbar etc. samples
  int pdgIdBoson =  Options.getParameter<int>("pdgIdBoson"); // absolute value of pdgId of the boson. Can be used only if the DoMatchingToBoson is set to true.
  
  //softdrop parameters
  beta = Options.getParameter<double>("beta");
  symmetry_cut = Options.getParameter<double>("symmetry_cut");
  R0 = Options.getParameter<double>("R0");

  //trimming parameters
  R_trimming = Options.getParameter<double>("R_trimming");
  PtFraction = Options.getParameter<double>("PtFraction");
  algorithm_Trimming = fastjet::kt_algorithm;
  
  //pruning parameters
  z_cut = Options.getParameter<double>("z_cut");
  R_Cut = Options.getParameter<double>("R_Cut");
  algorithm_Pruning = fastjet::antikt_algorithm;
  R_jet_def_pruning = Options.getParameter<double>("R_jet_def_pruning");


  
  // --- Read list of files to be analyzed and fill TChain
  TChain* lTree = new TChain("Events");
  FillChain(*lTree, inputFilesList);
  if (lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 

  cout << "This analysis will run on "<< maxEvents << " events" <<endl; 

  // --- Load branches
  fPFCand = new PFLoader (lTree,puppiConfig.c_str());
  fGen    = new GenLoader(lTree);
  if (doCMSSWJets) setupCMSSWJetReadOut(lTree, jetR);

  TEventInfo *eventInfo = new TEventInfo();
  lTree->SetBranchAddress("Info",&eventInfo);

  // --- Setup JEC on the fly
  /*
  // std::string cmsenv = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_7_patch2/src/";
  std::string cmsenv = "/afs/cern.ch/user/b/bmahakud/public/JEC/";
  
  std::vector<JetCorrectorParameters> corrParams;
  // corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"POSTLS162_V5_L1FastJet_AK7PF.txt"));
  // corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"POSTLS162_V5_L2Relative_AK7PF.txt"));
  // corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"POSTLS162_V5_L3Absolute_AK7PF.txt"));
  
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  // JetCorrectorParameters     param(cmsenv+"BaconProd/Utils/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt");
  JetCorrectorParameters     param(cmsenv+"POSTLS162_V5_Uncertainty_AK7PF.txt");
  */
 
  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters(L1FastJetJEC.c_str()));  
  corrParams.push_back(JetCorrectorParameters(L2RelativeJEC.c_str()));  
  corrParams.push_back(JetCorrectorParameters(L3AbsoluteJEC.c_str()));  
  if (L2L3ResidualJEC!="") corrParams.push_back(JetCorrectorParameters(L2L3ResidualJEC.c_str())); // 
  JetCorrectorParameters param(JECUncertainty.c_str());      

  FactorizedJetCorrector   *jetCorr = new FactorizedJetCorrector(corrParams);
  JetCorrectionUncertainty *jetUnc  = new JetCorrectionUncertainty(param);

  // --- Setup JetAlgos
  JetDefinition jet_def(antikt_algorithm,jetR);         // the jet definition....
  //JetDefinition jet_def_Pruning(antikt_algorithm,0.3);//this is a jet algorithm for pruning. Smaller radius to be used
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  
  // --- Setup cleansing
  JetDefinition subjet_def(kt_algorithm,0.2);
  JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);

  // --- Setup soft-killer
  SoftKiller soft_killer   (2.5,0.4);

  // --- Setup output trees
  TFile *fout = new TFile(fOut.c_str(),"RECREATE");

  TTree *genTree   = new TTree("gen"  , "gen"  );
  TTree *pfTree    = new TTree("pf"   , "pf"   );
  TTree *chsTree   = new TTree("chs"  , "chs"  );
  TTree *puppiTree = new TTree("puppi", "puppi");
  TTree *softkillerTree    = new TTree("softkiller", "softkiller");
  TTree *cmsswTree = new TTree("cmsswpf", "cmsswpf");

  JetInfo JGenInfo, JPFInfo, JCHSInfo, JPuppiInfo, JSoftKillerInfo, JCMSSWPFInfo;
  setupTree(genTree,   JGenInfo    , "" );
  setupTree(pfTree,    JPFInfo     , "" );
  setupTree(chsTree,   JCHSInfo    , "" );
  setupTree(puppiTree, JPuppiInfo  , "" );
  setupTree(softkillerTree, JSoftKillerInfo  , "" );
  if (doCMSSWJets) setupTree(cmsswTree, JCMSSWPFInfo, "" );
   

  // --- start loop over events
  for(int ientry = 0; ientry < maxEvents; ientry++) { 
    cout<<"event "<<ientry<<endl;
    // -- For each event build collections of particles (gen, puppi, etc..) to cluster
    Long64_t localEntry = lTree->LoadTree(ientry);
    fPFCand->load(localEntry);
    fGen   ->load(localEntry);  
    
    vector<PseudoJet> gen_event       = fGen   ->genFetch();
    vector<PseudoJet> pf_event        = fPFCand->pfFetch();
    vector<PseudoJet> chs_event       = fPFCand->pfchsFetch(-1);
    vector<PseudoJet> puppi_event     = fPFCand->puppiFetch();
    vector<PseudoJet> soft_event      = soft_killer(pf_event);

    // -- Cluster jets
    ClusterSequenceArea pGen    (gen_event    , jet_def, area_def);
    ClusterSequenceArea pPup    (puppi_event  , jet_def, area_def);
    ClusterSequenceArea pPF     (pf_event     , jet_def, area_def);
    ClusterSequenceArea pCHS    (chs_event    , jet_def, area_def);
    ClusterSequenceArea pSoft   (soft_event   , jet_def, area_def);

    vector<PseudoJet> genJets     = sorted_by_pt(pGen    .inclusive_jets(25.));
    vector<PseudoJet> puppiJets   = sorted_by_pt(pPup    .inclusive_jets(25.));
    vector<PseudoJet> pfJets      = sorted_by_pt(pPF     .inclusive_jets(25.));
    vector<PseudoJet> chsJets     = sorted_by_pt(pCHS    .inclusive_jets(25.));
    vector<PseudoJet> softJets    = sorted_by_pt(pSoft   .inclusive_jets(25.));

    lTree->GetEntry(ientry);
    int nPU = eventInfo->nPU;
  
    
    vfloat eta_Boson, phi_Boson;
    
    if (DoMatchingToBoson)
    {
      fGen -> selectBoson(pdgIdBoson);
      eta_Boson = fGen -> eta_Boson;
      phi_Boson = fGen -> phi_Boson;
    }


    cout << "\r" ;
    cout << "===> Processed " << ientry << " - Done : " << (float(ientry)/float(maxEvents))*100 << "%"  ;

    // save jet info in a tree
    fillGenJetsInfo(genJets, gen_event, JGenInfo, gsn_cleanser, nPU);  
    fillRecoJetsInfo(puppiJets, puppi_event, JPuppiInfo       , JGenInfo, false, jetCorr, jetUnc, gsn_cleanser,nPU, eta_Boson, phi_Boson ); 
    fillRecoJetsInfo(pfJets   , pf_event   , JPFInfo          , JGenInfo, false, jetCorr, jetUnc, gsn_cleanser,nPU, fGen -> eta_Boson,fGen -> phi_Boson );
    fillRecoJetsInfo(chsJets  , chs_event  , JCHSInfo         , JGenInfo, true , jetCorr, jetUnc, gsn_cleanser,nPU, fGen -> eta_Boson,fGen -> phi_Boson );
    fillRecoJetsInfo(softJets , soft_event , JSoftKillerInfo  , JGenInfo, true , jetCorr, jetUnc, gsn_cleanser,nPU, fGen -> eta_Boson,fGen -> phi_Boson );

    genTree->Fill();
    puppiTree->Fill();
    pfTree->Fill();
    chsTree->Fill();
    softkillerTree->Fill();

    if (doCMSSWJets)
      readCMSSWJet(ientry, lTree, *cmsswTree, genJets, JCMSSWPFInfo);
    
    
    fGen -> reset();
    
 
  }

  cout<<"done event loop"<<endl;

  // --- Write trees 
  fout->cd();

  genTree  ->Write();
  pfTree   ->Write();
  chsTree  ->Write();
  puppiTree->Write();
  softkillerTree->Write();

  if (doCMSSWJets)  cmsswTree->Write();
  cout<<"done write trees"<<endl;

}  

 
 
