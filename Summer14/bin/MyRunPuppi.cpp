//#include "../include/puppiContainer.hh"
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
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/Selector.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"


using namespace std;
using namespace fastjet;
using namespace contrib;

//Object Processors
GenLoader       *fGen      = 0; 
MuonLoader      *fMuon     = 0; 
PFLoader        *fPFCand   = 0; 

TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}

struct JetInfo {
  float pt;
  float ptcorr;
  float ptraw;
  float ptclean;
  float pttrim;
  float pttrimsafe;
  float ptconst;
  float ptunc;
  float eta;
  float phi;
  float m;
  float mraw;
  float mclean;
  float mtrim;
  float mtrimsafe;
  float mconst;
  // gen level info
  float ptgen;
  float etagen;
  float phigen;
  float mgen;
  int ismatched;
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



void setJet(PseudoJet &iJet, JetInfo &iJetI,std::vector<PseudoJet> &iParticles, bool iCHS,
	    FactorizedJetCorrector *iJetCorr,JetCorrectionUncertainty *iJetUnc,JetCleanser &gsn_cleanser, 
	    bool doGenMatching, vector<PseudoJet> genJets) {

    vector<PseudoJet> neutrals,chargedLV,chargedPU;
    getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
    PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
    
    // define safeAreaSub (PF)
    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
    JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
    Selector rho_range =  SelectorAbsRapMax(5.0);
    ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def);
    // the two background estimators
    JetMedianBackgroundEstimator bge_rho (rho_range, clust_seq_rho);
    JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_rho);
    BackgroundJetPtMDensity m_density;
    bge_rhom.set_jet_density_class(&m_density);
    
    // declare an area-median subtractor from this
    contrib::SafeAreaSubtractor *area_subtractor = 0;
    if(!iCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
    if( iCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom,SelectorIsPupCharged(),SelectorIsPupVertex());
    //iGMBE->set_particles(iParticles);
    PseudoJet lCorr =  (*area_subtractor)(iJet);
    fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
    PseudoJet lTrim     = (trimmer)(iJet);
    trimmer.set_subtractor(area_subtractor);
    PseudoJet lTrimSafe = (trimmer)(iJet);
    JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
    BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
    bge_rhoC.set_jet_density_class(scalarPtDensity);
    bge_rhoC.set_particles(iParticles);
    contrib::ConstituentSubtractor subtractor(&bge_rhoC);
    subtractor.use_common_bge_for_rho_and_rhom(true);
    PseudoJet lConstit = subtractor(iJet);

    // Find the gen jet matched to this reco jet
    int imatch = -1;
    if (doGenMatching) imatch = matchingIndex(iJet,genJets);
    
    //Finally apply the JEC
    double lJEC = correction(iJet,iJetCorr,bge_rho.rho());  
    double lUnc = unc       (iJet,iJetUnc);

    iJetI.pt          = lCorr     .pt();
    iJetI.ptcorr      = iJet      .pt()*lJEC;
    iJetI.ptraw       = iJet      .pt();
    iJetI.ptclean     = lClean    .pt();
    iJetI.pttrim      = lTrim     .pt();
    iJetI.pttrimsafe  = lTrimSafe .pt();
    iJetI.ptconst     = lConstit  .pt();
    iJetI.ptunc       = lUnc;
    iJetI.eta         = iJet      .eta();
    iJetI.phi         = iJet      .phi();
    iJetI.mraw        = iJet      .m();
    iJetI.m           = lCorr     .m();
    iJetI.mclean      = lClean    .m();
    iJetI.mtrim       = lTrim     .m();
    iJetI.mtrimsafe   = lTrimSafe .m();
    iJetI.mconst      = lConstit  .m();
    
    if (imatch > -1){
      iJetI.ptgen        = genJets[imatch].pt();
      iJetI.etagen       = genJets[imatch].eta();
      iJetI.phigen       = genJets[imatch].phi();
      iJetI.mgen         = genJets[imatch].m();
      iJetI.ismatched    = 1;
    }
    else {
      iJetI.ismatched    = 0;
    }

}


void setupTree(TTree *iTree, JetInfo &iJet, std::string iName) {
    iTree->Branch((iName+"pt"        ).c_str(),&iJet.pt        ,(iName+"pt/F"        ).c_str());
    iTree->Branch((iName+"ptcorr"    ).c_str(),&iJet.ptcorr    ,(iName+"ptcorr/F"    ).c_str());
    iTree->Branch((iName+"ptraw"     ).c_str(),&iJet.ptraw     ,(iName+"ptraw/F"     ).c_str());
    iTree->Branch((iName+"ptclean"   ).c_str(),&iJet.ptclean,   (iName+"ptclean/F"   ).c_str());
    iTree->Branch((iName+"pttrim"    ).c_str(),&iJet.pttrim    ,(iName+"pttrim/F"    ).c_str());
    iTree->Branch((iName+"pttrimsafe").c_str(),&iJet.pttrimsafe,(iName+"pttrimsafe/F").c_str());
    iTree->Branch((iName+"ptconst"   ).c_str(),&iJet.ptconst   ,(iName+"ptconst/F"   ).c_str());
    iTree->Branch((iName+"ptunc"     ).c_str(),&iJet.ptunc     ,(iName+"ptunc/F"     ).c_str());
    iTree->Branch((iName+"eta"       ).c_str(),&iJet.eta       ,(iName+"eta/F"       ).c_str());
    iTree->Branch((iName+"phi"       ).c_str(),&iJet.phi       ,(iName+"phi/F"       ).c_str());
    iTree->Branch((iName+"m"         ).c_str(),&iJet.m         ,(iName+"m/F"         ).c_str());
    iTree->Branch((iName+"mraw"      ).c_str(),&iJet.mraw      ,(iName+"mraw/F"      ).c_str());
    iTree->Branch((iName+"mtrim"     ).c_str(),&iJet.mtrim     ,(iName+"mtrim/F"     ).c_str());
    iTree->Branch((iName+"mtrimsafe" ).c_str(),&iJet.mtrimsafe ,(iName+"mtrimsafe/F" ).c_str());
    iTree->Branch((iName+"mclean"    ).c_str(),&iJet.mclean    ,(iName+"mclean/F"    ).c_str());
    iTree->Branch((iName+"mconst"    ).c_str(),&iJet.mconst    ,(iName+"mconst/F"    ).c_str());
    // gen info
    iTree->Branch((iName+"ptgen"     ).c_str(),&iJet.ptgen     ,(iName+"ptgen/F"     ).c_str());
    iTree->Branch((iName+"etagen"    ).c_str(),&iJet.etagen    ,(iName+"etagen/F"    ).c_str());
    iTree->Branch((iName+"phigen"    ).c_str(),&iJet.phigen    ,(iName+"phigen/F"    ).c_str());
    iTree->Branch((iName+"mgen"      ).c_str(),&iJet.mgen      ,(iName+"mgen/F"      ).c_str());
    iTree->Branch((iName+"ismatched" ).c_str(),&iJet.ismatched ,(iName+"ismatched/I" ).c_str());
}

//vector<PseudoJet> threeHardest(vector<PseudoJet> &iParts, JetDefinition &iJetDef, Selector &iSelector,std::vector<ClusterSequence> &iCSs) {
// cluster full event (hard + pileup)
//  vector<PseudoJet> threehardest = iSelector(sorted_by_pt(cs.inclusive_jets()));
//  iCSs.push_back(cs);
//  return threehardest;
//}

//PseudoJet match(PseudoJet &iJet,vector<PseudoJet> &iJets) {
//    for(unsigned int i0 = 0; i0 < iJets.size(); i0++) {
//        double pEta = fabs(iJet.eta()-iJets[i0].eta());
//        double pPhi = fabs(iJet.phi() - iJets[i0].phi());
//        if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
//        if(sqrt(pEta*pEta+pPhi*pPhi) > 0.3) continue;
//        return iJets[i0];
//    }
//    return PseudoJet();
//}

void clear(JetInfo &iJet) {
    iJet.pt         = -1;
    iJet.ptraw      = -1;
    iJet.ptclean    = -1;
    iJet.pttrim     = -1;
    iJet.pttrimsafe = -1;
    iJet.eta        = -1;
    iJet.phi        = -1;
    iJet.m          = -1;
    iJet.mraw       = -1;
    iJet.mtrim      = -1;
    iJet.mtrimsafe  = -1;
    iJet.mclean     = -1;
    iJet.mconst     = -1;

    iJet.ptgen      = -1;
    iJet.etagen     = -1;
    iJet.phigen     = -1;
    iJet.mgen       = -1;
    iJet.ismatched  = -1;
}

void fillTree(vector<PseudoJet> &iJets, vector<PseudoJet> &iParticles, JetInfo &iJetInfo, bool isCHS, 
	      FactorizedJetCorrector *jetCorr, JetCorrectionUncertainty *ijetUnc, JetCleanser &gsn_cleanser, 
	      bool doGenMatching, vector<PseudoJet> genJets, TTree &iTree)
{

  for (unsigned int j = 0; j < iJets.size(); j++){
    clear(iJetInfo);
    setJet( iJets[j], iJetInfo, iParticles, isCHS, jetCorr, ijetUnc, gsn_cleanser, doGenMatching, genJets);
    cout << iTree.GetName() << "  " << iJetInfo.pt << "  "<< iJetInfo.eta <<endl;
    // fill tree for each event
  }

  iTree.Fill();    


}


//---------------------------------------------------------------------------------------------------------------
//--- MAIN PROGRAM
//---------------------------------------------------------------------------------------------------------------
int main (int argc, char ** argv) {

  // args 
  std::string fIn   = argv[1];        // input file name
  int maxEvents     = atoi(argv[2]);  // max events
  std::string fOut  = argv[3];        // output name
  float jetR        = atof(argv[4]);  // jet cone size      
  bool        lGen  = atoi(argv[5]);  // analyze gen

  // --- Read input file
  TTree *lTree = load(fIn); 
  if (lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 
  fPFCand       = new PFLoader (lTree,"Puppi_cff.py");
  if(lGen) fGen = new GenLoader(lTree);


  // --- Setup JEC on the fly
  std::string cmsenv = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_7_patch2/src/";
  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  JetCorrectorParameters     param(cmsenv+"BaconProd/Utils/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt");
  FactorizedJetCorrector   *jetCorr = new FactorizedJetCorrector(corrParams);
  JetCorrectionUncertainty *jetUnc  = new JetCorrectionUncertainty(param);

  // --- Setup JetAlgos
  JetDefinition jet_def(antikt_algorithm,jetR);         // the jet definition....
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  //Selector selector = SelectorNHardest(3);   // definition of a selector for the three hardest jets
  
  // --- Setup cleansing
  JetDefinition subjet_def(kt_algorithm,0.2);
  JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);


  // --- Setup output trees
  TTree *genTree   = new TTree("gen"  , "gen"  );
  TTree *pfTree    = new TTree("pf"   , "pf"   );
  TTree *chsTree   = new TTree("chs"  , "chs"  );
  TTree *puppiTree = new TTree("puppi", "puppi");

  JetInfo JGenInfo, JPFInfo, JCHSInfo, JPuppiInfo;    
  setupTree(genTree,   JGenInfo  , "" );
  setupTree(pfTree,    JPFInfo   , "" );
  setupTree(chsTree,   JCHSInfo  , "" );
  setupTree(puppiTree, JPuppiInfo, "" );

  
  // --- start loop over events
  for(int ientry = 0; ientry < maxEvents; ientry++) { 

    if(ientry % 2 == 0) 
      std::cout << "===> Processed " << ientry << " - Done : " << (float(ientry)/float(maxEvents)) << std::endl;
    
    // -- For each event build collections of particles (gen, puppi, etc..) to cluster
    fPFCand->load(ientry);
    fGen   ->load(ientry); 
    vector<PseudoJet> gen_event       = fGen   ->genFetch();
    vector<PseudoJet> pf_event        = fPFCand->pfFetch();
    vector<PseudoJet> chs_event       = fPFCand->pfchsFetch(-1);
    vector<PseudoJet> puppi_event     = fPFCand->puppiFetch();

    // -- Cluster jets
    ClusterSequenceArea pGen    (gen_event    , jet_def, area_def);
    ClusterSequenceArea pPup    (puppi_event  , jet_def, area_def);
    ClusterSequenceArea pPF     (pf_event     , jet_def, area_def);
    ClusterSequenceArea pCHS    (chs_event    , jet_def, area_def);

    //vector<PseudoJet> genJets     = selector(sorted_by_pt(pGen    .inclusive_jets()));
    //vector<PseudoJet> puppiJets   = selector(sorted_by_pt(pPup    .inclusive_jets()));
    //vector<PseudoJet> pfJets      = selector(sorted_by_pt(pPF     .inclusive_jets()));
    //vector<PseudoJet> chsJets     = selector(sorted_by_pt(pCHS    .inclusive_jets()));

    vector<PseudoJet> genJets     = sorted_by_pt(pGen    .inclusive_jets(20.));
    vector<PseudoJet> puppiJets   = sorted_by_pt(pPup    .inclusive_jets(20.));
    vector<PseudoJet> pfJets      = sorted_by_pt(pPF     .inclusive_jets(20.));
    vector<PseudoJet> chsJets     = sorted_by_pt(pCHS    .inclusive_jets(20.));

    // save jet info in a tree
    bool doGenMatching = true;
    fillTree(genJets  , gen_event  , JGenInfo  , false, jetCorr, jetUnc, gsn_cleanser, false       ,  genJets, *genTree);
    fillTree(puppiJets, puppi_event, JPuppiInfo, false, jetCorr, jetUnc, gsn_cleanser, doGenMatching, genJets, *puppiTree);
    fillTree(pfJets   , pf_event   , JPFInfo   , false, jetCorr, jetUnc, gsn_cleanser, doGenMatching, genJets, *pfTree);
    fillTree(chsJets  , chs_event  , JCHSInfo  , false, jetCorr, jetUnc, gsn_cleanser, doGenMatching, genJets, *chsTree);
    
  }

  TFile *fout = new TFile(fOut.c_str(),"RECREATE");
  fout->cd();
  genTree->Write();
  pfTree->Write();
  chsTree->Write();
  puppiTree->Write();
}  
