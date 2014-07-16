#include "../include/JetTreeAnalyzer.h"


// --- RecoToGenMatching -----------------------------------------------------------------
int RecoToGenMatching(float eta, float phi, vector<float> *etagen, vector<float> *phigen){
  float dRMatching = 0.3;
  float rmin = 9999.;
  int imatch = -1;
  for(unsigned int i = 0; i < etagen->size(); i++) {
    double dEta = fabs(eta - etagen->at(i));
    double dPhi = fabs(phi - phigen->at(i));
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


// --- constructor ------------------------------------------------------
JetTreeAnalyzer::JetTreeAnalyzer(TTree *tree, TTree *gentree, string treetype){
  treetype_ = treetype;
  Init(tree, gentree);
}

// ----------------------------------------------------------------------
JetTreeAnalyzer::~JetTreeAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

// --- Init tree --------------------------------------------------------
void JetTreeAnalyzer::Init(TTree *tree, TTree *gentree)
{
  // The Init() function is called when the selector needs to initialize                                                                                                   
  // a new tree or chain. Typically here the branch addresses and branch                                                    
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated     
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF 
  // (once per file to be processed).                                                                                                                                            
  
  pt = 0;
  ptcorr = 0;
  ptraw = 0;
  ptunc = 0;
  eta = 0;
  phi = 0;
  m = 0;
  mraw = 0;
  ptclean = 0;
  mclean = 0;
  pttrim_Rtrim_020_Ptfrac_005 = 0;
  mtrim_Rtrim_020_Ptfrac_005 = 0;
  pttrimsafe_Rtrim_020_Ptfrac_005 = 0;
  mtrimsafe_Rtrim_020_Ptfrac_005 = 0;
  pttrim_Rtrim_010_Ptfrac_003 = 0;
  mtrim_Rtrim_010_Ptfrac_003 = 0;
  pttrimsafe_Rtrim_010_Ptfrac_003 = 0;
  mtrimsafe_Rtrim_010_Ptfrac_003 = 0;
  pttrim_Rtrim_020_Ptfrac_003 = 0;
  mtrim_Rtrim_020_Ptfrac_003 = 0;
  pttrimsafe_Rtrim_020_Ptfrac_003 = 0;
  mtrimsafe_Rtrim_020_Ptfrac_003 = 0;
  pttrim_Rtrim_030_Ptfrac_003 = 0;
  mtrim_Rtrim_030_Ptfrac_003 = 0;
  pttrimsafe_Rtrim_030_Ptfrac_003 = 0;
  mtrimsafe_Rtrim_030_Ptfrac_003 = 0;
  ptconst = 0;
  mconst = 0;
  ptpruned_zcut_010_R_cut_050 = 0;
  mpruned_zcut_010_R_cut_050 = 0;
  ptprunedsafe_zcut_010_R_cut_050 = 0;
  mprunedsafe_zcut_010_R_cut_050 = 0;
  QGLikelihood_pr_zcut_010_R_cut_050 = 0;
  QGLikelihood_pr_sub1_zcut_010_R_cut_050 = 0;
  QGLikelihood_pr_sub2_zcut_010_R_cut_050 = 0;
  ptpruned_zcut_005_R_cut_050 = 0;
  mpruned_zcut_005_R_cut_050 = 0;
  ptprunedsafe_zcut_005_R_cut_050 = 0;
  mprunedsafe_zcut_005_R_cut_050 = 0;
  QGLikelihood_pr_zcut_005_R_cut_050 = 0;
  QGLikelihood_pr_sub1_zcut_005_R_cut_050 = 0;
  QGLikelihood_pr_sub2_zcut_005_R_cut_050 = 0;
  ptpruned_zcut_005_R_cut_075 = 0;
  mpruned_zcut_005_R_cut_075 = 0;
  ptprunedsafe_zcut_005_R_cut_075 = 0;
  mprunedsafe_zcut_005_R_cut_075 = 0;
  QGLikelihood_pr_zcut_005_R_cut_075 = 0;
  QGLikelihood_pr_sub1_zcut_005_R_cut_075 = 0;
  QGLikelihood_pr_sub2_zcut_005_R_cut_075 = 0;
  ptpruned_zcut_010_R_cut_075 = 0;
  mpruned_zcut_010_R_cut_075 = 0;
  ptprunedsafe_zcut_010_R_cut_075 = 0;
  mprunedsafe_zcut_010_R_cut_075 = 0;
  QGLikelihood_pr_zcut_010_R_cut_075 = 0;
  QGLikelihood_pr_sub1_zcut_010_R_cut_075 = 0;
  QGLikelihood_pr_sub2_zcut_010_R_cut_075 = 0;
  ptsoftdrop_beta20 = 0;
  msoftdrop_beta20 = 0;
  ptsoftdropsafe_beta20 = 0;
  msoftdropsafe_beta20 = 0;
  ptsoftdrop_beta00 = 0;
  msoftdrop_beta00 = 0;
  ptsoftdropsafe_beta00 = 0;
  msoftdropsafe_beta00 = 0;
  ptsoftdrop_beta10 = 0;
  msoftdrop_beta10 = 0;
  ptsoftdropsafe_beta10 = 0;
  msoftdropsafe_beta10 = 0;
  ptsoftdrop_betam10 = 0;
  msoftdrop_betam10 = 0;
  ptsoftdropsafe_betam10 = 0;
  msoftdropsafe_betam10 = 0;
  nparticles = 0;
  nneutrals = 0;
  ncharged = 0;
  sdsymmetry = 0;
  sddeltar = 0;
  sdmu = 0;
  sdenergyloss = 0;
  sdarea = 0;
  sdnconst = 0;
  mfiltsoftdrop = 0;
  tau1 = 0;
  tau2 = 0;
  tau3 = 0;
  tau4 = 0;
  tau5 = 0;
  tau1_pr = 0;
  tau2_pr = 0;
  tau3_pr = 0;
  tau4_pr = 0;
  tau5_pr = 0;
  tau1_softdrop = 0;
  tau2_softdrop = 0;
  tau3_softdrop = 0;
  tau4_softdrop = 0;
  tau5_softdrop = 0;
  Qjets = 0;
  charge_k05 = 0;
  charge_k07 = 0;
  charge_k10 = 0;
  ecf_beta_05 = 0;
  ecf_beta_10 = 0;
  ecf_beta_15 = 0;
  ecf_beta_20 = 0;
  hepmass = 0;
  hepwmass = 0;
  hepm01 = 0;
  hepm02 = 0;
  hepm12 = 0;
  hepm12m012 = 0;
  hepatanm02m01 = 0;
  cmsmass = 0;
  cmsminmass = 0;
  cmshelicity = 0;
  cmsnsubjets = 0;
  cmsnsubjets = 0;
  // gen vars
  ptgen = 0;
  etagen = 0;
  phigen = 0;
  mgen = 0;
  mrawgen = 0;
  mtrimgen = 0;
  mtrimsafegen = 0;
  mcleangen = 0;
  mconstgen = 0;
  pttrimgen = 0;
  pttrimsafegen = 0;
  ptcleangen = 0;
  ptconstgen = 0;
  imatch = 0;
  flavourgen = 0;
  msoftdropgen = 0;
  msoftdropsafegen = 0;
  mfiltsoftdropgen = 0;
  is_MatchedToBoson = 0;


  // Set branch addresses and branch pointers                                                                                                                                                                   
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  // disable all branches 
  fChain->SetBranchStatus("*",0);

  // enable only branches that are used
  fChain->SetBranchStatus("npu", 1);

  fChain->SetBranchStatus("eta", 1);
  fChain->SetBranchStatus("phi", 1);

  fChain->SetBranchStatus("pt", 1);
  fChain->SetBranchStatus("ptcorr", 1);
  fChain->SetBranchStatus("ptraw", 1);
  fChain->SetBranchStatus("ptclean", 1);
  fChain->SetBranchStatus("ptconst", 1);
  fChain->SetBranchStatus("pttrim_Rtrim_020_Ptfrac_005", 1);
  fChain->SetBranchStatus("pttrimsafe_Rtrim_020_Ptfrac_005", 1);
  
  fChain->SetBranchStatus("m", 1);
  fChain->SetBranchStatus("mraw", 1);
  fChain->SetBranchStatus("mclean", 1);
  fChain->SetBranchStatus("mconst", 1);
  fChain->SetBranchStatus("mtrim_Rtrim_020_Ptfrac_005", 1);
  fChain->SetBranchStatus("mtrimsafe_Rtrim_020_Ptfrac_005", 1);
  fChain->SetBranchStatus("msoftdrop_beta20", 1);
  fChain->SetBranchStatus("msoftdropsafe_beta20", 1);

  fChain->SetBranchStatus("nparticles", 1);
  fChain->SetBranchStatus("nneutrals", 1);
  fChain->SetBranchStatus("ncharged", 1);

  fChain->SetBranchStatus("tau1", 1);
  fChain->SetBranchStatus("tau2", 1);
  fChain->SetBranchStatus("tau1_softdrop", 1);
  fChain->SetBranchStatus("tau2_softdrop", 1);

  if ( treetype_ != "gen")
    fChain->SetBranchStatus("imatch", 1);
  
  // set branch addresses
  fChain->SetBranchAddress("npu", &npu, &b_npu);

  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("phi", &phi, &b_phi);

  fChain->SetBranchAddress("pt", &pt, &b_pt);
  fChain->SetBranchAddress("ptcorr", &ptcorr, &b_ptcorr);
  fChain->SetBranchAddress("ptraw", &ptraw, &b_ptraw);
  fChain->SetBranchAddress("ptclean", &ptclean, &b_ptclean);
  fChain->SetBranchAddress("ptconst", &ptconst, &b_ptconst);
  fChain->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_005", &pttrim_Rtrim_020_Ptfrac_005, &b_pttrim_Rtrim_020_Ptfrac_005);
  fChain->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_005", &pttrimsafe_Rtrim_020_Ptfrac_005, &b_pttrimsafe_Rtrim_020_Ptfrac_005);
  
  fChain->SetBranchAddress("m", &m, &b_m);
  fChain->SetBranchAddress("mraw", &mraw, &b_mraw);
  fChain->SetBranchAddress("mclean", &mclean, &b_mclean);
  fChain->SetBranchAddress("mconst", &mconst, &b_mconst);
  fChain->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_005", &mtrim_Rtrim_020_Ptfrac_005, &b_mtrim_Rtrim_020_Ptfrac_005);
  fChain->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_005", &mtrimsafe_Rtrim_020_Ptfrac_005, &b_mtrimsafe_Rtrim_020_Ptfrac_005);
  fChain->SetBranchAddress("msoftdrop_beta20", &msoftdrop_beta20, &b_msoftdrop_beta20);
  fChain->SetBranchAddress("msoftdropsafe_beta20", &msoftdropsafe_beta20, &b_msoftdropsafe_beta20);

  fChain->SetBranchAddress("nparticles", &nparticles, &b_nparticles);
  fChain->SetBranchAddress("nneutrals", &nneutrals, &b_nneutrals);
  fChain->SetBranchAddress("ncharged", &ncharged, &b_ncharged);

  fChain->SetBranchAddress("tau1", &tau1, &b_tau1);
  fChain->SetBranchAddress("tau2", &tau2, &b_tau2);
  fChain->SetBranchAddress("tau1_softdrop", &tau1_softdrop, &b_tau1_softdrop);
  fChain->SetBranchAddress("tau2_softdrop", &tau2_softdrop, &b_tau2_softdrop);

  if (treetype_ != "gen"){
    fChain->SetBranchAddress("imatch", &imatch, &b_imatch);
  }


  // gen tree

  if (treetype_ != "gen"){
    // Set branch addresses and branch pointers                                                                                                                                                                   
    if (!gentree) return;
    fChain2 = gentree;
    fCurrent2 = -1;
    fChain2->SetMakeClass(1);
    
    // disable all branches 
    fChain2->SetBranchStatus("*",0);

    // enable only branches that are used
    fChain2->SetBranchStatus("eta", 1);
    fChain2->SetBranchStatus("phi", 1);
    fChain2->SetBranchStatus("pt", 1);
    fChain2->SetBranchStatus("ptclean", 1);
    fChain2->SetBranchStatus("ptconst", 1);
    fChain2->SetBranchStatus("pttrim_Rtrim_020_Ptfrac_005", 1);
    fChain2->SetBranchStatus("pttrimsafe_Rtrim_020_Ptfrac_005", 1);
    
    fChain2->SetBranchStatus("m", 1);
    fChain2->SetBranchStatus("mraw", 1);
    fChain2->SetBranchStatus("mclean", 1);
    fChain2->SetBranchStatus("mconst", 1);
    fChain2->SetBranchStatus("mtrim_Rtrim_020_Ptfrac_005", 1);
    fChain2->SetBranchStatus("mtrimsafe_Rtrim_020_Ptfrac_005", 1);
    fChain2->SetBranchStatus("msoftdrop_beta20", 1);
    fChain2->SetBranchStatus("msoftdropsafe_beta20", 1);
    
    // set branch addresses
    fChain2->SetBranchAddress("eta", &etagen, &b_etagen);
    fChain2->SetBranchAddress("phi", &phigen, &b_phigen);
    
    fChain2->SetBranchAddress("pt", &ptgen, &b_ptgen);
    fChain2->SetBranchAddress("ptclean", &ptcleangen, &b_ptcleangen);
    fChain2->SetBranchAddress("ptconst", &ptconstgen, &b_ptconstgen);
    fChain2->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_005", &pttrimgen, &b_pttrimgen);
    fChain2->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_005", &pttrimsafegen, &b_pttrimsafegen);
    
    fChain2->SetBranchAddress("m", &mgen, &b_mgen);
    fChain2->SetBranchAddress("mraw", &mrawgen, &b_mrawgen);
    fChain2->SetBranchAddress("mclean", &mcleangen, &b_mcleangen);
    fChain2->SetBranchAddress("mconst", &mconstgen, &b_mconstgen);
    fChain2->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_005", &mtrimgen, &b_mtrimgen);
    fChain2->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_005", &mtrimsafegen, &b_mtrimsafegen);
    fChain2->SetBranchAddress("msoftdrop_beta20", &msoftdropgen, &b_msoftdropgen);
    fChain2->SetBranchAddress("msoftdropsafe_beta20", &msoftdropsafegen, &b_msoftdropsafegen);
  }
}


// --- get Tree entry ----------------------------------------------------------------
Int_t JetTreeAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.                                                                                                                                                                           
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}


// --- Book histograms ---------------------------------------------------------------
void JetTreeAnalyzer::bookHistograms(std::string suffix){

  std::cout << "Booking histograms for " << suffix.c_str() << " tree" << std::endl;

  hnjets = new TH1F(("hnjets"+suffix).c_str(), ("hnjets"+suffix).c_str(), 50, 0, 50 );

  // all jets
  hptgen = new TH1F(("hptgen"+suffix).c_str(), ("hptgen"+suffix).c_str(), 2000, 0, 2000 );
  hptgen_pu = new TH1F(("hptgen_pu"+suffix).c_str(), ("hptgen_pu"+suffix).c_str(), 2000, 0, 2000 );
  hptgen_good = new TH1F(("hptgen_good"+suffix).c_str(), ("hptgen_good"+suffix).c_str(), 2000, 0, 2000 );

  hptraw = new TH1F(("hptraw"+suffix).c_str(), ("hptraw"+suffix).c_str(), 2000, 0, 2000 );
  hptraw_pu = new TH1F(("hptraw_pu"+suffix).c_str(), ("hptraw_pu"+suffix).c_str(), 2000, 0, 2000 );
  hptraw_good = new TH1F(("hptraw_good"+suffix).c_str(), ("hptraw_good"+suffix).c_str(), 2000, 0, 2000 );
  hptraw_response = new TH1F(("hptraw_response"+suffix).c_str(), ("hptraw_response"+suffix).c_str(), 200, -100, 100 );

  hpt = new TH1F(("hpt"+suffix).c_str(), ("hpt"+suffix).c_str(), 2000, 0, 2000 );
  hpt_pu = new TH1F(("hpt_pu"+suffix).c_str(), ("hpt_pu"+suffix).c_str(), 2000, 0, 2000 );
  hpt_good = new TH1F(("hpt_good"+suffix).c_str(), ("hpt_good"+suffix).c_str(), 2000, 0, 2000 );
  hpt_response = new TH1F(("hpt_response"+suffix).c_str(), ("hpt_response"+suffix).c_str(), 200, -100, 100 );

  hptcorr = new TH1F(("hptcorr"+suffix).c_str(), ("hptcorr"+suffix).c_str(), 2000, 0, 2000 );
  hptcorr_pu = new TH1F(("hptcorr_pu"+suffix).c_str(), ("hptcorr_pu"+suffix).c_str(), 2000, 0, 2000 );
  hptcorr_good = new TH1F(("hptcorr_good"+suffix).c_str(), ("hptcorr_good"+suffix).c_str(), 2000, 0, 2000 );
  hptcorr_response = new TH1F(("hptcorr_response"+suffix).c_str(), ("hptcorr_response"+suffix).c_str(), 200, -100, 100 );

  hpttrim = new TH1F(("hpttrim"+suffix).c_str(), ("hpttrim"+suffix).c_str(), 2000, 0, 2000 );
  hpttrim_pu = new TH1F(("hpttrim_pu"+suffix).c_str(), ("hpttrim_pu"+suffix).c_str(), 2000, 0, 2000 );
  hpttrim_good = new TH1F(("hpttrim_good"+suffix).c_str(), ("hpttrim_good"+suffix).c_str(), 2000, 0, 2000 );
  hpttrim_response = new TH1F(("hpttrim_response"+suffix).c_str(), ("hpttrim_response"+suffix).c_str(), 200, -100, 100 );

  hpttrimsafe = new TH1F(("hpttrimsafe"+suffix).c_str(), ("hpttrimsafe"+suffix).c_str(), 2000, 0, 2000 );
  hpttrimsafe_pu = new TH1F(("hpttrimsafe_pu"+suffix).c_str(), ("hpttrimsafe_pu"+suffix).c_str(), 2000, 0, 2000 );
  hpttrimsafe_good = new TH1F(("hpttrimsafe_good"+suffix).c_str(), ("hpttrimsafe_good"+suffix).c_str(), 2000, 0, 2000 );
  hpttrimsafe_response = new TH1F(("hpttrimsafe_response"+suffix).c_str(), ("hpttrimsafe_response"+suffix).c_str(), 200, -100, 100 );

  hptsoftdrop = new TH1F(("hptsoftdrop"+suffix).c_str(), ("hptsoftdrop"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdrop_pu = new TH1F(("hptsoftdrop_pu"+suffix).c_str(), ("hptsoftdrop_pu"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdrop_good = new TH1F(("hptsoftdrop_good"+suffix).c_str(), ("hptsoftdrop_good"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdrop_response = new TH1F(("hptsoftdrop_response"+suffix).c_str(), ("hptsoftdrop_response"+suffix).c_str(), 200, -100, 100 );

  hptsoftdropsafe = new TH1F(("hptsoftdropsafe"+suffix).c_str(), ("hptsoftdropsafe"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdropsafe_pu = new TH1F(("hptsoftdropsafe_pu"+suffix).c_str(), ("hptsoftdropsafe_pu"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdropsafe_good = new TH1F(("hptsoftdropsafe_good"+suffix).c_str(), ("hptsoftdropsafe_good"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdropsafe_response = new TH1F(("hptsoftdropsafe_response"+suffix).c_str(), ("hptsoftdropsafe_response"+suffix).c_str(), 200, -100, 100 );

  hptconst = new TH1F(("hptconst"+suffix).c_str(), ("hptconst"+suffix).c_str(), 2000, 0, 2000 );
  hptconst_pu = new TH1F(("hptconst_pu"+suffix).c_str(), ("hptconst_pu"+suffix).c_str(), 2000, 0, 2000 );
  hptconst_good = new TH1F(("hptconst_good"+suffix).c_str(), ("hptconst_good"+suffix).c_str(), 2000, 0, 2000 );
  hptconst_response = new TH1F(("hptconst_response"+suffix).c_str(), ("hptconst_response"+suffix).c_str(), 200, -100, 100 );

  hptclean = new TH1F(("hptclean"+suffix).c_str(), ("hptclean"+suffix).c_str(), 2000, 0, 2000 );
  hptclean_pu = new TH1F(("hptclean_pu"+suffix).c_str(), ("hptclean_pu"+suffix).c_str(), 2000, 0, 2000 );
  hptclean_good = new TH1F(("hptclean_good"+suffix).c_str(), ("hptclean_good"+suffix).c_str(), 2000, 0, 2000 );
  hptclean_response = new TH1F(("hptclean_response"+suffix).c_str(), ("hptclean_response"+suffix).c_str(), 200, -100, 100 );

  heta = new TH1F(("heta"+suffix).c_str(), ("heta"+suffix).c_str(), 100, -5, 5 );
  heta_pu = new TH1F(("heta_pu"+suffix).c_str(), ("heta_pu"+suffix).c_str(), 100, -5, 5 );
  heta_good = new TH1F(("heta_good"+suffix).c_str(), ("heta_good"+suffix).c_str(), 100, -5, 5 );

  hnpu = new TH1F(("hnpu"+suffix).c_str(), ("hnpu"+suffix).c_str(), 100, 0, 100 );
  hnpu_pu = new TH1F(("hnpu_pu"+suffix).c_str(), ("hnpu_pu"+suffix).c_str(), 100, 0, 100 );
  hnpu_good = new TH1F(("hnpu_good"+suffix).c_str(), ("hnpu_good"+suffix).c_str(), 100, 0, 100 );

  hm = new TH1F(("hm"+suffix).c_str(), ("hm"+suffix).c_str(), 200, 0, 200 );
  hm_response = new TH1F(("hm_response"+suffix).c_str(), ("hm_response"+suffix).c_str(), 200, -100, 100 );

  hmraw = new TH1F(("hmraw"+suffix).c_str(), ("hmraw"+suffix).c_str(), 200, 0, 200 );
  hmraw_response = new TH1F(("hmraw_response"+suffix).c_str(), ("hmraw_response"+suffix).c_str(), 200, -100, 100 );
  
  hmtrim = new TH1F(("hmtrim"+suffix).c_str(), ("hmtrim"+suffix).c_str(), 200, 0, 200 );
  hmtrim_response = new TH1F(("hmtrim_response"+suffix).c_str(), ("hmtrim_response"+suffix).c_str(), 200, -100, 100 );

  hmtrimsafe = new TH1F(("hmtrimsafe"+suffix).c_str(), ("hmtrimsafe"+suffix).c_str(), 200, 0, 200 );
  hmtrimsafe_response = new TH1F(("hmtrimsafe_response"+suffix).c_str(), ("hmtrimsafe_response"+suffix).c_str(), 200, -100, 100 );

  hmclean = new TH1F(("hmclean"+suffix).c_str(), ("hmclean"+suffix).c_str(), 200, 0, 200 );
  hmclean_response = new TH1F(("hmclean_response"+suffix).c_str(), ("hmclean_response"+suffix).c_str(), 200, -100, 100 );

  hmconst = new TH1F(("hmconst"+suffix).c_str(), ("hmconst"+suffix).c_str(), 200, 0, 200 );
  hmconst_response = new TH1F(("hmconst_response"+suffix).c_str(), ("hmconst_response"+suffix).c_str(), 200, -100, 100 );

  hmsoftdrop = new TH1F(("hmsoftdrop"+suffix).c_str(), ("hmsoftdrop"+suffix).c_str(), 200, 0, 200 );
  hmsoftdrop_response = new TH1F(("hmsoftdrop_response"+suffix).c_str(), ("hmsoftdrop_response"+suffix).c_str(), 200, -100, 100 );

  hmsoftdropsafe = new TH1F(("hmsoftdropsafe"+suffix).c_str(), ("hmsoftdropsafe"+suffix).c_str(), 200, 0, 200 );
  hmsoftdropsafe_response = new TH1F(("hmsoftdropsafe_response"+suffix).c_str(), ("hmsoftdropsafe_response"+suffix).c_str(), 200, -100, 100 );

  hnparticles = new TH1F(("hnparticles"+suffix).c_str(), ("hnparticles"+suffix).c_str(), 1000, 0, 1000 );
  hnneutrals = new TH1F(("hnneutrals"+suffix).c_str(), ("hnneutrals"+suffix).c_str(), 1000, 0, 1000 );
  hncharged = new TH1F(("hncharged"+suffix).c_str(), ("hncharged"+suffix).c_str(), 1000, 0, 1000 );

  htau21          = new TH1F(("htau21"+suffix).c_str(), ("htau21"+suffix).c_str(), 1000, 0, 1 );
  htau21_softdrop = new TH1F(("htau21_softdrop"+suffix).c_str(), ("htau21_softdrop"+suffix).c_str(), 1000, 0, 1 );

  // leading jet 
  hptraw_leadjet = new TH1F(("hptraw_leadjet"+suffix).c_str(), ("hptraw_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptraw_pu_leadjet = new TH1F(("hptraw_pu_leadjet"+suffix).c_str(), ("hptraw_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptraw_good_leadjet = new TH1F(("hptraw_good_leadjet"+suffix).c_str(), ("hptraw_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptraw_response_leadjet = new TH1F(("hptraw_response_leadjet"+suffix).c_str(), ("hptraw_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hpt_leadjet = new TH1F(("hpt_leadjet"+suffix).c_str(), ("hpt_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpt_pu_leadjet = new TH1F(("hpt_pu_leadjet"+suffix).c_str(), ("hpt_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpt_good_leadjet = new TH1F(("hpt_good_leadjet"+suffix).c_str(), ("hpt_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpt_response_leadjet = new TH1F(("hpt_response_leadjet"+suffix).c_str(), ("hpt_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hptcorr_leadjet = new TH1F(("hptcorr_leadjet"+suffix).c_str(), ("hptcorr_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptcorr_pu_leadjet = new TH1F(("hptcorr_pu_leadjet"+suffix).c_str(), ("hptcorr_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptcorr_good_leadjet = new TH1F(("hptcorr_good_leadjet"+suffix).c_str(), ("hptcorr_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptcorr_response_leadjet = new TH1F(("hptcorr_response_leadjet"+suffix).c_str(), ("hptcorr_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hpttrim_leadjet = new TH1F(("hpttrim_leadjet"+suffix).c_str(), ("hpttrim_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpttrim_pu_leadjet = new TH1F(("hpttrim_pu_leadjet"+suffix).c_str(), ("hpttrim_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpttrim_good_leadjet = new TH1F(("hpttrim_good_leadjet"+suffix).c_str(), ("hpttrim_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpttrim_response_leadjet = new TH1F(("hpttrim_response_leadjet"+suffix).c_str(), ("hpttrim_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hpttrimsafe_leadjet = new TH1F(("hpttrimsafe_leadjet"+suffix).c_str(), ("hpttrimsafe_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpttrimsafe_pu_leadjet = new TH1F(("hpttrimsafe_pu_leadjet"+suffix).c_str(), ("hpttrimsafe_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpttrimsafe_good_leadjet = new TH1F(("hpttrimsafe_good_leadjet"+suffix).c_str(), ("hpttrimsafe_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpttrimsafe_response_leadjet = new TH1F(("hpttrimsafe_response_leadjet"+suffix).c_str(), ("hpttrimsafe_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hptsoftdrop_leadjet = new TH1F(("hptsoftdrop_leadjet"+suffix).c_str(), ("hptsoftdrop_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdrop_pu_leadjet = new TH1F(("hptsoftdrop_pu_leadjet"+suffix).c_str(), ("hptsoftdrop_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdrop_good_leadjet = new TH1F(("hptsoftdrop_good_leadjet"+suffix).c_str(), ("hptsoftdrop_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdrop_response_leadjet = new TH1F(("hptsoftdrop_response_leadjet"+suffix).c_str(), ("hptsoftdrop_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hptsoftdropsafe_leadjet = new TH1F(("hptsoftdropsafe_leadjet"+suffix).c_str(), ("hptsoftdropsafe_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdropsafe_pu_leadjet = new TH1F(("hptsoftdropsafe_pu_leadjet"+suffix).c_str(), ("hptsoftdropsafe_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdropsafe_good_leadjet = new TH1F(("hptsoftdropsafe_good_leadjet"+suffix).c_str(), ("hptsoftdropsafe_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptsoftdropsafe_response_leadjet = new TH1F(("hptsoftdropsafe_response_leadjet"+suffix).c_str(), ("hptsoftdropsafe_response_leadjet"+suffix).c_str(), 200, -100, 100 );
  
  hptclean_leadjet = new TH1F(("hptclean_leadjet"+suffix).c_str(), ("hptclean_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptclean_pu_leadjet = new TH1F(("hptclean_pu_leadjet"+suffix).c_str(), ("hptclean_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptclean_good_leadjet = new TH1F(("hptclean_good_leadjet"+suffix).c_str(), ("hptclean_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptclean_response_leadjet = new TH1F(("hptclean_response_leadjet"+suffix).c_str(), ("hptclean_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hptconst_leadjet = new TH1F(("hptconst_leadjet"+suffix).c_str(), ("hptconst_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptconst_pu_leadjet = new TH1F(("hptconst_pu_leadjet"+suffix).c_str(), ("hptconst_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptconst_good_leadjet = new TH1F(("hptconst_good_leadjet"+suffix).c_str(), ("hptconst_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hptconst_response_leadjet = new TH1F(("hptconst_response_leadjet"+suffix).c_str(), ("hptconst_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  heta_leadjet = new TH1F(("heta_leadjet"+suffix).c_str(), ("heta_leadjet"+suffix).c_str(), 100, -5, 5 );
  heta_pu_leadjet = new TH1F(("heta_pu_leadjet"+suffix).c_str(), ("heta_pu_leadjet"+suffix).c_str(), 100, -5, 5 );
  heta_good_leadjet = new TH1F(("heta_good_leadjet"+suffix).c_str(), ("heta_good_leadjet"+suffix).c_str(), 100, -5, 5 );

  hmraw_leadjet = new TH1F(("hmraw_leadjet"+suffix).c_str(), ("hmraw_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmraw_response_leadjet = new TH1F(("hmraw_response_leadjet"+suffix).c_str(), ("hmraw_response_leadjet"+suffix).c_str(), 200, -100, 100 );
  
  hm_leadjet = new TH1F(("hm_leadjet"+suffix).c_str(), ("hm_leadjet"+suffix).c_str(), 200, 0, 200 );
  hm_response_leadjet = new TH1F(("hm_response_leadjet"+suffix).c_str(), ("hm_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmtrim_leadjet = new TH1F(("hmtrim_leadjet"+suffix).c_str(), ("hmtrim_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmtrim_response_leadjet = new TH1F(("hmtrim_response_leadjet"+suffix).c_str(), ("hmtrim_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmtrimsafe_leadjet = new TH1F(("hmtrimsafe_leadjet"+suffix).c_str(), ("hmtrimsafe_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmtrimsafe_response_leadjet = new TH1F(("hmtrimsafe_response_leadjet"+suffix).c_str(), ("hmtrimsafe_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmclean_leadjet          = new TH1F(("hmclean_leadjet"+suffix).c_str(), ("hmclean_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmclean_response_leadjet = new TH1F(("hmclean_response_leadjet"+suffix).c_str(), ("hmclean_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  //for (int i = 0; i < 7; i++){
  //  hmclean_leadjet[i]          = new TH1F((Form("hmclean%d_leadjet",i)+suffix).c_str(),(Form("hmclean%d_leadjet",i)+suffix).c_str(), 200, 0, 200 );
  //  hmclean_response_leadjet[i] = new TH1F((Form("hmclean%d_response_leadjet",i)+suffix).c_str(), (Form("hmclean%d_response_leadjet",i)+suffix).c_str(), 200, -100, 100 );
  //}

  hmconst_leadjet = new TH1F(("hmconst_leadjet"+suffix).c_str(), ("hmconst_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmconst_response_leadjet = new TH1F(("hmconst_response_leadjet"+suffix).c_str(), ("hmconst_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmsoftdrop_leadjet = new TH1F(("hmsoftdrop_leadjet"+suffix).c_str(), ("hmsoftdrop_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmsoftdrop_response_leadjet = new TH1F(("hmsoftdrop_response_leadjet"+suffix).c_str(), ("hmsoftdrop_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmsoftdropsafe_leadjet = new TH1F(("hmsoftdropsafe_leadjet"+suffix).c_str(), ("hmsoftdropsafe_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmsoftdropsafe_response_leadjet = new TH1F(("hmsoftdropsafe_response_leadjet"+suffix).c_str(), ("hmsoftdropsafe_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hnparticles_leadjet = new TH1F(("hnparticles_leadjet "+suffix).c_str(), ("hnparticles_leadjet "+suffix).c_str(), 100, 0, 100 );
  hnneutrals_leadjet  = new TH1F(("hnneutrals_leadjet "+suffix).c_str(), ("hnneutrals_leadjet "+suffix).c_str(), 100, 0, 100 );
  hncharged_leadjet   = new TH1F(("hncharged_leadjet "+suffix).c_str(), ("hncharged_leadjet "+suffix).c_str(), 100, 0, 100 );

  htau21_leadjet          = new TH1F(("htau21_leadjet"+suffix).c_str(), ("htau21_leadjet"+suffix).c_str(), 1000, 0, 1 );
  htau21_softdrop_leadjet = new TH1F(("htau21_softdrop_leadjet"+suffix).c_str(), ("htau21_softdrop_leadjet"+suffix).c_str(), 1000, 0, 1 );

  // 2d histograms

  //hptraw_response_vs_pt     = new TH2F(("hptraw_response_vs_pt"+suffix).c_str(), ("hptraw_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  //hpt_response_vs_pt        = new TH2F(("hpt_response_vs_pt"+suffix).c_str(), ("hpt_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  //hptcorr_response_vs_pt        = new TH2F(("hptcorr_response_vs_pt"+suffix).c_str(), ("hptcorr_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hptraw_response_vs_pt     = new TH2F(("hptraw_response_vs_pt"+suffix).c_str(), ("hptraw_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -1, 1 );
  hpt_response_vs_pt        = new TH2F(("hpt_response_vs_pt"+suffix).c_str(), ("hpt_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -1, 1 );
  hptcorr_response_vs_pt    = new TH2F(("hptcorr_response_vs_pt"+suffix).c_str(), ("hptcorr_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -1, 1 );
  hmraw_response_vs_pt      = new TH2F(("hmraw_response_vs_pt"+suffix).c_str(), ("hmraw_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hm_response_vs_pt         = new TH2F(("hm_response_vs_pt"+suffix).c_str(), ("hm_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmtrim_response_vs_pt     = new TH2F(("hmtrim_response_vs_pt"+suffix).c_str(), ("hmtrim_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmtrimsafe_response_vs_pt = new TH2F(("hmtrimsafe_response_vs_pt"+suffix).c_str(), ("hmtrimsafe_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmclean_response_vs_pt    = new TH2F(("hmclean_response_vs_pt"+suffix).c_str(), ("hmclean_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmconst_response_vs_pt    = new TH2F(("hmconst_response_vs_pt"+suffix).c_str(), ("hmconst_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmsoftdrop_response_vs_pt    = new TH2F(("hmsoftdrop_response_vs_pt"+suffix).c_str(), ("hmsoftdrop_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmsoftdropsafe_response_vs_pt    = new TH2F(("hmsoftdropsafe_response_vs_pt"+suffix).c_str(), ("hmsoftdropsafe_response_vs_pt"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );

  //hptraw_response_vs_eta   = new TH2F(("hptraw_response_vs_eta"+suffix).c_str(), ("hptraw_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  //hpt_response_vs_eta      = new TH2F(("hpt_response_vs_eta"+suffix).c_str(), ("hpt_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  //hptcorr_response_vs_eta  = new TH2F(("hptcorr_response_vs_eta"+suffix).c_str(), ("hptcorr_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hptraw_response_vs_eta     = new TH2F(("hptraw_response_vs_eta"+suffix).c_str(), ("hptraw_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -1, 1 );
  hpt_response_vs_eta        = new TH2F(("hpt_response_vs_eta"+suffix).c_str(), ("hpt_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -1, 1 );
  hptcorr_response_vs_eta    = new TH2F(("hptcorr_response_vs_eta"+suffix).c_str(), ("hptcorr_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -1, 1 );
  hmraw_response_vs_eta      = new TH2F(("hmraw_response_vs_eta"+suffix).c_str(), ("hmraw_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hm_response_vs_eta         = new TH2F(("hm_response_vs_eta"+suffix).c_str(), ("hm_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmtrim_response_vs_eta     = new TH2F(("hmtrim_response_vs_eta"+suffix).c_str(), ("hmtrim_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmtrimsafe_response_vs_eta = new TH2F(("hmtrimsafe_response_vs_eta"+suffix).c_str(), ("hmtrimsafe_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmclean_response_vs_eta    = new TH2F(("hmclean_response_vs_eta"+suffix).c_str(), ("hmclean_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmconst_response_vs_eta    = new TH2F(("hmconst_response_vs_eta"+suffix).c_str(), ("hmconst_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmsoftdrop_response_vs_eta    = new TH2F(("hmsoftdrop_response_vs_eta"+suffix).c_str(), ("hmsoftdrop_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmsoftdropsafe_response_vs_eta    = new TH2F(("hmsoftdropsafe_response_vs_eta"+suffix).c_str(), ("hmsoftdropsafe_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );

  //hptraw_response_vs_npu   = new TH2F(("hptraw_response_vs_npu"+suffix).c_str(), ("hptraw_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  //hpt_response_vs_npu      = new TH2F(("hpt_response_vs_npu"+suffix).c_str(), ("hpt_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  //hptcorr_response_vs_npu  = new TH2F(("hptcorr_response_vs_npu"+suffix).c_str(), ("hptcorr_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hptraw_response_vs_npu     = new TH2F(("hptraw_response_vs_npu"+suffix).c_str(), ("hptraw_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -1, 1 );
  hpt_response_vs_npu        = new TH2F(("hpt_response_vs_npu"+suffix).c_str(), ("hpt_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -1, 1 );
  hptcorr_response_vs_npu    = new TH2F(("hptcorr_response_vs_npu"+suffix).c_str(), ("hptcorr_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -1, 1 );
  hmraw_response_vs_npu      = new TH2F(("hmraw_response_vs_npu"+suffix).c_str(), ("hmraw_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hm_response_vs_npu         = new TH2F(("hm_response_vs_npu"+suffix).c_str(), ("hm_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmtrim_response_vs_npu     = new TH2F(("hmtrim_response_vs_npu"+suffix).c_str(), ("hmtrim_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmtrimsafe_response_vs_npu = new TH2F(("hmtrimsafe_response_vs_npu"+suffix).c_str(), ("hmtrimsafe_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmclean_response_vs_npu    = new TH2F(("hmclean_response_vs_npu"+suffix).c_str(), ("hmclean_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmconst_response_vs_npu    = new TH2F(("hmconst_response_vs_npu"+suffix).c_str(), ("hmconst_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmsoftdrop_response_vs_npu    = new TH2F(("hmsoftdrop_response_vs_npu"+suffix).c_str(), ("hmsoftdrop_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmsoftdropsafe_response_vs_npu    = new TH2F(("hmsoftdropsafe_response_vs_npu"+suffix).c_str(), ("hmsoftdropsafe_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );

  hptraw_vs_npu     = new TH2F(("hptraw_vs_npu"+suffix).c_str(), ("hptraw_vs_npu"+suffix).c_str(), 100, 0, 100, 3000, 0, 3000 );
  hpt_vs_npu        = new TH2F(("hpt_vs_npu"+suffix).c_str(), ("hpt_vs_npu"+suffix).c_str(), 100, 0, 100, 3000, 0, 3000 );
  hptcorr_vs_npu    = new TH2F(("hptcorr_vs_npu"+suffix).c_str(), ("hptcorr_vs_npu"+suffix).c_str(), 100, 0, 100, 3000, 0, 3000 );
  hmraw_vs_npu      = new TH2F(("hmraw_vs_npu"+suffix).c_str(), ("hmraw_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hm_vs_npu         = new TH2F(("hm_vs_npu"+suffix).c_str(), ("hm_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmtrim_vs_npu     = new TH2F(("hmtrim_vs_npu"+suffix).c_str(), ("hmtrim_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmtrimsafe_vs_npu = new TH2F(("hmtrimsafe_vs_npu"+suffix).c_str(), ("hmtrimsafe_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmclean_vs_npu    = new TH2F(("hmclean_vs_npu"+suffix).c_str(), ("hmclean_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmconst_vs_npu    = new TH2F(("hmconst_vs_npu"+suffix).c_str(), ("hmconst_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmsoftdrop_vs_npu    = new TH2F(("hmsoftdrop_vs_npu"+suffix).c_str(), ("hmsoftdrop_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmsoftdropsafe_vs_npu    = new TH2F(("hmsoftdropsafe_vs_npu"+suffix).c_str(), ("hmsoftdropsafe_vs_npu"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );


  // leading jet only

  //hptraw_response_vs_pt_leadjet     = new TH2F(("hptraw_response_vs_pt_leadjet"+suffix).c_str(), ("hptraw_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  //hpt_response_vs_pt_leadjet        = new TH2F(("hpt_response_vs_pt_leadjet"+suffix).c_str(), ("hpt_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  //hptcorr_response_vs_pt_leadjet        = new TH2F(("hptcorr_response_vs_pt_leadjet"+suffix).c_str(), ("hptcorr_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hptraw_response_vs_pt_leadjet     = new TH2F(("hptraw_response_vs_pt_leadjet"+suffix).c_str(), ("hptraw_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -1, 1 );
  hpt_response_vs_pt_leadjet        = new TH2F(("hpt_response_vs_pt_leadjet"+suffix).c_str(), ("hpt_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -1, 1 );
  hptcorr_response_vs_pt_leadjet    = new TH2F(("hptcorr_response_vs_pt_leadjet"+suffix).c_str(), ("hptcorr_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -1, 1 );
  hmraw_response_vs_pt_leadjet      = new TH2F(("hmraw_response_vs_pt_leadjet"+suffix).c_str(), ("hmraw_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hm_response_vs_pt_leadjet         = new TH2F(("hm_response_vs_pt_leadjet"+suffix).c_str(), ("hm_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmtrim_response_vs_pt_leadjet     = new TH2F(("hmtrim_response_vs_pt_leadjet"+suffix).c_str(), ("hmtrim_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmtrimsafe_response_vs_pt_leadjet = new TH2F(("hmtrimsafe_response_vs_pt_leadjet"+suffix).c_str(), ("hmtrimsafe_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmclean_response_vs_pt_leadjet    = new TH2F(("hmclean_response_vs_pt_leadjet"+suffix).c_str(), ("hmclean_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmconst_response_vs_pt_leadjet    = new TH2F(("hmconst_response_vs_pt_leadjet"+suffix).c_str(), ("hmconst_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmsoftdrop_response_vs_pt_leadjet    = new TH2F(("hmsoftdrop_response_vs_pt_leadjet"+suffix).c_str(), ("hmsoftdrop_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );
  hmsoftdropsafe_response_vs_pt_leadjet    = new TH2F(("hmsoftdropsafe_response_vs_pt_leadjet"+suffix).c_str(), ("hmsoftdropsafe_response_vs_pt_leadjet"+suffix).c_str(), 2000, 0, 2000, 200, -100, 100 );

  //hptraw_response_vs_eta_leadjet   = new TH2F(("hptraw_response_vs_eta_leadjet"+suffix).c_str(), ("hptraw_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  //hpt_response_vs_eta_leadjet      = new TH2F(("hpt_response_vs_eta_leadjet"+suffix).c_str(), ("hpt_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  //hptcorr_response_vs_eta_leadjet  = new TH2F(("hptcorr_response_vs_eta_leadjet"+suffix).c_str(), ("hptcorr_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hptraw_response_vs_eta_leadjet     = new TH2F(("hptraw_response_vs_eta_leadjet"+suffix).c_str(), ("hptraw_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -1, 1 );
  hpt_response_vs_eta_leadjet        = new TH2F(("hpt_response_vs_eta_leadjet"+suffix).c_str(), ("hpt_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -1, 1 );
  hptcorr_response_vs_eta_leadjet    = new TH2F(("hptcorr_response_vs_eta_leadjet"+suffix).c_str(), ("hptcorr_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -1, 1 );
  hmraw_response_vs_eta_leadjet      = new TH2F(("hmraw_response_vs_eta_leadjet"+suffix).c_str(), ("hmraw_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hm_response_vs_eta_leadjet         = new TH2F(("hm_response_vs_eta_leadjet"+suffix).c_str(), ("hm_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmtrim_response_vs_eta_leadjet     = new TH2F(("hmtrim_response_vs_eta_leadjet"+suffix).c_str(), ("hmtrim_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmtrimsafe_response_vs_eta_leadjet = new TH2F(("hmtrimsafe_response_vs_eta_leadjet"+suffix).c_str(), ("hmtrimsafe_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmclean_response_vs_eta_leadjet    = new TH2F(("hmclean_response_vs_eta_leadjet"+suffix).c_str(), ("hmclean_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmconst_response_vs_eta_leadjet    = new TH2F(("hmconst_response_vs_eta_leadjet"+suffix).c_str(), ("hmconst_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmsoftdrop_response_vs_eta_leadjet    = new TH2F(("hmsoftdrop_response_vs_eta_leadjet"+suffix).c_str(), ("hmsoftdrop_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmsoftdropsafe_response_vs_eta_leadjet    = new TH2F(("hmsoftdropsafe_response_vs_eta_leadjet"+suffix).c_str(), ("hmsoftdropsafe_response_vs_eta_leadjet"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );

  //hptraw_response_vs_npu_leadjet   = new TH2F(("hptraw_response_vs_npu_leadjet"+suffix).c_str(), ("hptraw_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  //hpt_response_vs_npu_leadjet      = new TH2F(("hpt_response_vs_npu_leadjet"+suffix).c_str(), ("hpt_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  //hptcorr_response_vs_npu_leadjet  = new TH2F(("hptcorr_response_vs_npu_leadjet"+suffix).c_str(), ("hptcorr_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hptraw_response_vs_npu_leadjet     = new TH2F(("hptraw_response_vs_npu_leadjet"+suffix).c_str(), ("hptraw_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -1, 1 );
  hpt_response_vs_npu_leadjet        = new TH2F(("hpt_response_vs_npu_leadjet"+suffix).c_str(), ("hpt_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -1, 1 );
  hptcorr_response_vs_npu_leadjet    = new TH2F(("hptcorr_response_vs_npu_leadjet"+suffix).c_str(), ("hptcorr_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -1, 1 );
  hmraw_response_vs_npu_leadjet      = new TH2F(("hmraw_response_vs_npu_leadjet"+suffix).c_str(), ("hmraw_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hm_response_vs_npu_leadjet         = new TH2F(("hm_response_vs_npu_leadjet"+suffix).c_str(), ("hm_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmtrim_response_vs_npu_leadjet     = new TH2F(("hmtrim_response_vs_npu_leadjet"+suffix).c_str(), ("hmtrim_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmtrimsafe_response_vs_npu_leadjet = new TH2F(("hmtrimsafe_response_vs_npu_leadjet"+suffix).c_str(), ("hmtrimsafe_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmclean_response_vs_npu_leadjet    = new TH2F(("hmclean_response_vs_npu_leadjet"+suffix).c_str(), ("hmclean_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmconst_response_vs_npu_leadjet    = new TH2F(("hmconst_response_vs_npu_leadjet"+suffix).c_str(), ("hmconst_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmsoftdrop_response_vs_npu_leadjet    = new TH2F(("hmsoftdrop_response_vs_npu_leadjet"+suffix).c_str(), ("hmsoftdrop_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmsoftdropsafe_response_vs_npu_leadjet    = new TH2F(("hmsoftdropsafe_response_vs_npu_leadjet"+suffix).c_str(), ("hmsoftdropsafe_response_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );

  hptraw_vs_npu_leadjet     = new TH2F(("hptraw_vs_npu_leadjet"+suffix).c_str(), ("hptraw_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 3000, 0, 3000 );
  hpt_vs_npu_leadjet        = new TH2F(("hpt_vs_npu_leadjet"+suffix).c_str(), ("hpt_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 3000, 0, 3000 );
  hptcorr_vs_npu_leadjet    = new TH2F(("hptcorr_vs_npu_leadjet"+suffix).c_str(), ("hptcorr_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 3000, 0, 3000 );
  hmraw_vs_npu_leadjet      = new TH2F(("hmraw_vs_npu_leadjet"+suffix).c_str(), ("hmraw_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hm_vs_npu_leadjet         = new TH2F(("hm_vs_npu_leadjet"+suffix).c_str(), ("hm_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmtrim_vs_npu_leadjet     = new TH2F(("hmtrim_vs_npu_leadjet"+suffix).c_str(), ("hmtrim_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmtrimsafe_vs_npu_leadjet = new TH2F(("hmtrimsafe_vs_npu_leadjet"+suffix).c_str(), ("hmtrimsafe_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmclean_vs_npu_leadjet    = new TH2F(("hmclean_vs_npu_leadjet"+suffix).c_str(), ("hmclean_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmconst_vs_npu_leadjet    = new TH2F(("hmconst_vs_npu_leadjet"+suffix).c_str(), ("hmconst_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmsoftdrop_vs_npu_leadjet    = new TH2F(("hmsoftdrop_vs_npu_leadjet"+suffix).c_str(), ("hmsoftdrop_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );
  hmsoftdropsafe_vs_npu_leadjet    = new TH2F(("hmsoftdropsafe_vs_npu_leadjet"+suffix).c_str(), ("hmsoftdropsafe_vs_npu_leadjet"+suffix).c_str(), 100, 0, 100, 1000, 0, 1000 );

}


// --- Fill histograms ---------------------------------------------------------------
void JetTreeAnalyzer::fillHistograms(int maxEntries, float minPt, float maxPt, float minAbsEta, float maxAbsEta){

  std::cout << "Filling histograms..." << std::endl;
  
  if (fChain==0){
    std::cout<<"Error: cannot open " << fChain->GetName() << std::endl;
    exit(0);
  }

  if (maxEntries == -1)
    maxEntries = fChain->GetEntries();


  int iclean  = 3;

  for (int entry = 0; entry < maxEntries; entry++){

    fChain->GetEntry(entry);
    if (treetype_ != "gen")    
      fChain2->GetEntry(entry);
    
    if (entry%100==0) std::cout << "Analyzing entry : " << entry << "\r" << std::flush;
    //if (entry%1==0) std::cout << "Analyzing entry : " << entry << "\r" << std::endl;
    
    // --- Loop over jets in this event                                                                                                                                      
    int nj = 0;

    float thispt = 0;

    //  all jets
    for (unsigned int j = 0; j < ptraw->size(); j++){
    

      thispt = pt->at(j); // use pt 
      //float thispt = ptcorr->at(j); // use ptcorr 
      //float thispt = ptraw->at(j); // use pt raw
      
      if (thispt < minPt)  continue;
      if (thispt > maxPt)  continue;

      if (fabs(eta->at(j)) < minAbsEta) continue;
      if (fabs(eta->at(j)) > maxAbsEta) continue;

      nj++;

      // gen-reco matching
      int matchInd = -1;
      if (treetype_ != "gen")
	matchInd = RecoToGenMatching(eta->at(j), phi->at(j), etagen, phigen);
      //matchInd = imatch->at(j); // this is not working now...

      float genpt = -999;
      float genptconst = -999;
      float genptclean = -999;
      float genm = -999;
      float genmraw = -999; 
      float genmconst = -999;
      float genmclean = -999;

      // get gen level quantities
      if (matchInd > -1){
	genpt        = ptgen->at(matchInd);
	genptconst   = ptconst->at(matchInd);
	genptclean   = ptclean->at(matchInd);
	genm         = mgen->at(matchInd);
	genmraw      = mrawgen->at(matchInd);
	genmconst    = mconstgen->at(matchInd);
	genmclean    = mgen->at(matchInd); // cleansed variables are saved only for the leading jet!
      }

      // fill some basic distributions
      hptraw    -> Fill(ptraw->at(j));
      hpt       -> Fill(pt->at(j));
      hptcorr   -> Fill(ptcorr->at(j));
      hptconst  -> Fill(ptconst->at(j));
      if (j==0)       
	hptclean  -> Fill(ptclean->at(iclean)); // cleansing only for leading jet
      heta      -> Fill(eta->at(j));
      hnpu      -> Fill(npu);
      hmraw     -> Fill(mraw->at(j));
      hm        -> Fill(m->at(j));
      hmconst   -> Fill(mconst->at(j));
      if (j==0)       
	hmclean  -> Fill(mclean->at(iclean)); // cleansing only for leading jet

      hnparticles -> Fill(nparticles->at(j));
      hnneutrals  -> Fill(nneutrals->at(j));
      hncharged   -> Fill(ncharged->at(j));

      if (tau1->at(j)!=999 && tau2->at(j)!=999) htau21          ->Fill( tau2->at(j)/tau1->at(j) );
      if (tau1_softdrop->at(j)!=999 && tau2_softdrop->at(j)!=999) htau21_softdrop ->Fill( tau2_softdrop->at(j)/tau1_softdrop->at(j));

      // split matched ad unmatched jets
      if (matchInd == -1 ) {
	hptgen    -> Fill(genpt);
	hptgen_pu -> Fill(genpt);
	hptraw_pu -> Fill(ptraw->at(j));
	hpt_pu    -> Fill(pt->at(j));
	hptcorr_pu-> Fill(ptcorr->at(j));
	hptconst_pu->Fill(ptconst->at(j));
	if (j == 0) hptclean_pu->Fill(ptclean->at(iclean));
	heta_pu   -> Fill(eta->at(j));
	hnpu_pu   -> Fill(npu);
	if (j == 0){
	  hptraw_pu_leadjet -> Fill(ptraw->at(j));
	  hpt_pu_leadjet    -> Fill(pt->at(j));
	  hptcorr_pu_leadjet-> Fill(ptcorr->at(j));
	  hptconst_pu_leadjet->Fill(ptconst->at(j));
	  hptclean_pu_leadjet->Fill(ptclean->at(iclean));
	  heta_pu_leadjet   -> Fill(eta->at(j));
	}
      }
      else {
	hptgen      -> Fill(genpt);
	hptgen_good -> Fill(genpt);
	hptraw_good -> Fill(ptraw->at(j));
	hpt_good    -> Fill(pt->at(j));
	hptcorr_good-> Fill(ptcorr->at(j));
	hptconst_good->Fill(ptconst->at(j));
	if (j == 0) hptclean_good->Fill(ptclean->at(iclean));
	heta_good   -> Fill(eta->at(j));
	hnpu_good   -> Fill(npu);

	if (j == 0){
	  hptraw_good_leadjet -> Fill(ptraw->at(j));
	  hpt_good_leadjet    -> Fill(pt->at(j));
	  hptcorr_good_leadjet-> Fill(ptcorr->at(j));
	  hptconst_good_leadjet-> Fill(ptconst->at(j));
	  hptclean_good_leadjet->Fill(ptclean->at(iclean));
	  heta_good_leadjet   -> Fill(eta->at(j));
	}
      }

      // -- leading jet only
      if (j == 0){
	hptraw_leadjet    -> Fill(ptraw->at(j));
	hpt_leadjet       -> Fill(pt->at(j));
	hptcorr_leadjet   -> Fill(ptcorr->at(j));
	hptconst_leadjet  -> Fill(ptconst->at(j));
	hptclean_leadjet  -> Fill(ptclean->at(iclean)); // linear cleanser (mclean only for leading jet)
	heta_leadjet      -> Fill(eta->at(j));
	hmraw_leadjet     -> Fill(mraw->at(j));
	hm_leadjet        -> Fill(m->at(j));
	hmconst_leadjet   -> Fill(mconst->at(j));
	hmclean_leadjet   -> Fill(mclean->at(iclean)); // linear cleanser (mclean only for leading jet)

	hnparticles_leadjet -> Fill(nparticles->at(j));
	hnneutrals_leadjet  -> Fill(nneutrals->at(j));
	hncharged_leadjet   -> Fill(ncharged->at(j));
	
	htau21_leadjet          ->Fill( tau2->at(j)/tau1->at(j) );
	htau21_softdrop_leadjet ->Fill( tau2_softdrop->at(j)/tau1_softdrop->at(j));
      }


      // -- response plots
      if (matchInd > -1){

	hptraw_response     -> Fill(ptraw->at(j)-genpt);
	hpt_response        -> Fill(pt->at(j)-genpt);
	hptcorr_response    -> Fill(ptcorr->at(j)-genpt);
	hptconst_response   -> Fill(ptconst->at(j)-genptconst);
	if (j==0) hptclean_response -> Fill(ptclean->at(iclean)-genptclean);
	hmraw_response      -> Fill(mraw->at(j)-genmraw);
	hm_response         -> Fill(m->at(j)-genm);
	hmconst_response    -> Fill(mconst->at(j)-genmconst);
	if (j==0) hmclean_response -> Fill(mclean->at(iclean)-genmclean);

	// 2d plots
	//hptraw_response_vs_pt     -> Fill(genpt,ptraw->at(j)-genpt);
	//hpt_response_vs_pt        -> Fill(genpt,pt->at(j)-genpt);
	//hptcorr_response_vs_pt    -> Fill(genpt,ptcorr->at(j)-genpt);
	hptraw_response_vs_pt     -> Fill(genpt,ptraw->at(j)/genpt-1); // fill with (pt-ptgen)/ptgen
	hpt_response_vs_pt        -> Fill(genpt,pt->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	hptcorr_response_vs_pt    -> Fill(genpt,ptcorr->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	hmraw_response_vs_pt      -> Fill(genpt,mraw->at(j)-genmraw);
	hm_response_vs_pt         -> Fill(genpt,m->at(j)-genm);
	hmconst_response_vs_pt    -> Fill(genpt,mconst->at(j)-genmconst);

	//hptraw_response_vs_eta     -> Fill(eta->at(j),ptraw->at(j)-genpt);
	//hpt_response_vs_eta        -> Fill(eta->at(j),pt->at(j)-genpt);
	//hptcorr_response_vs_eta    -> Fill(eta->at(j),ptcorr->at(j)-genpt);
	hptraw_response_vs_eta     -> Fill(eta->at(j),ptraw->at(j)/genpt-1); // fill with (pt-ptgen)/ptgen
	hpt_response_vs_eta        -> Fill(eta->at(j),pt->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	hptcorr_response_vs_eta    -> Fill(eta->at(j),ptcorr->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	hmraw_response_vs_eta      -> Fill(eta->at(j),mraw->at(j)-genmraw);
	hm_response_vs_eta         -> Fill(eta->at(j),m->at(j)-genm);
	hmconst_response_vs_eta    -> Fill(eta->at(j),mconst->at(j)-genmconst);

	//hptraw_response_vs_npu     -> Fill(npu,ptraw->at(j)-genpt);
	//hpt_response_vs_npu        -> Fill(npu,pt->at(j)-genpt);
	//hptcorr_response_vs_npu    -> Fill(npu,ptcorr->at(j)-genpt);
	hptraw_response_vs_npu     -> Fill(npu,ptraw->at(j)/genpt-1); // fill with (pt-ptgen)/ptgen
	hpt_response_vs_npu        -> Fill(npu,pt->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	hptcorr_response_vs_npu    -> Fill(npu,ptcorr->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	hmraw_response_vs_npu      -> Fill(npu,mraw->at(j)-genmraw);
	hm_response_vs_npu         -> Fill(npu,m->at(j)-genm);
	hmconst_response_vs_npu    -> Fill(npu,mconst->at(j)-genmconst);

	hptraw_vs_npu     -> Fill(npu,ptraw->at(j)); // fill with (pt-ptgen)/ptgen
	hpt_vs_npu        -> Fill(npu,pt->at(j));    // fill with (pt-ptgen)/ptgen
	hptcorr_vs_npu    -> Fill(npu,ptcorr->at(j));    // fill with (pt-ptgen)/ptgen
	hmraw_vs_npu      -> Fill(npu,mraw->at(j));
	hm_vs_npu         -> Fill(npu,m->at(j));
	hmconst_vs_npu    -> Fill(npu,mconst->at(j));

	// leading jet
	if (j == 0){
	  hptraw_response_leadjet     -> Fill(ptraw->at(j)-genpt);
	  hpt_response_leadjet        -> Fill(pt->at(j)-genpt);
	  hptcorr_response_leadjet    -> Fill(ptcorr->at(j)-genpt);
	  hptclean_response_leadjet   -> Fill(ptclean->at(iclean)-genptclean);
	  hmraw_response_leadjet      -> Fill(mraw->at(j)-genmraw);
	  hm_response_leadjet         -> Fill(m->at(j)-genm);
	  hmconst_response_leadjet    -> Fill(mconst->at(j)-genmconst);
	  if (j==0)
	    hmclean_response_leadjet  -> Fill(mclean->at(iclean)-genmclean);

	  // 2d plots
	  //hptraw_response_vs_pt_leadjet     -> Fill(genpt,ptraw->at(j)-genpt);
	  //hpt_response_vs_pt_leadjet        -> Fill(genpt,pt->at(j)-genpt);
	  //hptcorr_response_vs_pt_leadjet    -> Fill(genpt,ptcorr->at(j)-genpt);
	  hptraw_response_vs_pt_leadjet     -> Fill(genpt,ptraw->at(j)/genpt-1); // fill with (pt-ptgen)/ptgen
	  hpt_response_vs_pt_leadjet        -> Fill(genpt,pt->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	  hptcorr_response_vs_pt_leadjet    -> Fill(genpt,ptcorr->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	  hmraw_response_vs_pt_leadjet      -> Fill(genpt,mraw->at(j)-genmraw);
	  hm_response_vs_pt_leadjet         -> Fill(genpt,m->at(j)-genm);
	  hmconst_response_vs_pt_leadjet    -> Fill(genpt,mconst->at(j)-genmconst);
	  if (j==0) 
	    hmclean_response_vs_pt_leadjet  -> Fill(genpt,mclean->at(iclean)-genmclean);

	  //hptraw_response_vs_eta_leadjet     -> Fill(eta->at(j),ptraw->at(j)-genpt);
	  //hpt_response_vs_eta_leadjet        -> Fill(eta->at(j),pt->at(j)-genpt);
	  //hptcorr_response_vs_eta_leadjet    -> Fill(eta->at(j),ptcorr->at(j)-genpt);
	  hptraw_response_vs_eta_leadjet     -> Fill(eta->at(j),ptraw->at(j)/genpt-1); // fill with (pt-ptgen)/ptgen
	  hpt_response_vs_eta_leadjet        -> Fill(eta->at(j),pt->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	  hptcorr_response_vs_eta_leadjet    -> Fill(eta->at(j),ptcorr->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	  hmraw_response_vs_eta_leadjet      -> Fill(eta->at(j),mraw->at(j)-genmraw);
	  hm_response_vs_eta_leadjet         -> Fill(eta->at(j),m->at(j)-genm);
	  hmconst_response_vs_eta_leadjet    -> Fill(eta->at(j),mconst->at(j)-genmconst);
	  if (j==0)  
	    hmclean_response_vs_eta_leadjet  -> Fill(eta->at(j),mclean->at(iclean)-genmclean);

	  //hptraw_response_vs_npu_leadjet     -> Fill(npu,ptraw->at(j)-genpt);
	  //hpt_response_vs_npu_leadjet        -> Fill(npu,pt->at(j)-genpt);
	  //hptcorr_response_vs_npu_leadjet    -> Fill(npu,ptcorr->at(j)-genpt);
	  hptraw_response_vs_npu_leadjet     -> Fill(npu,ptraw->at(j)/genpt-1); // fill with (pt-ptgen)/ptgen
	  hpt_response_vs_npu_leadjet        -> Fill(npu,pt->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	  hptcorr_response_vs_npu_leadjet    -> Fill(npu,ptcorr->at(j)/genpt-1);    // fill with (pt-ptgen)/ptgen
	  hmraw_response_vs_npu_leadjet      -> Fill(npu,mraw->at(j)-genmraw);
	  hm_response_vs_npu_leadjet         -> Fill(npu,m->at(j)-genm);
	  hmconst_response_vs_npu_leadjet    -> Fill(npu,mconst->at(j)-genmconst);
	  if (j==0)  
	    hmclean_response_vs_npu_leadjet  -> Fill(npu,mclean->at(iclean)-genmclean);


	hptraw_vs_npu_leadjet     -> Fill(npu,ptraw->at(j)); // fill with (pt-ptgen)/ptgen
	hpt_vs_npu_leadjet        -> Fill(npu,pt->at(j));    // fill with (pt-ptgen)/ptgen
	hptcorr_vs_npu_leadjet    -> Fill(npu,ptcorr->at(j));    // fill with (pt-ptgen)/ptgen
	hmraw_vs_npu_leadjet      -> Fill(npu,mraw->at(j));
	hm_vs_npu_leadjet         -> Fill(npu,m->at(j));
	hmconst_vs_npu_leadjet    -> Fill(npu,mconst->at(j));
	if (j==0)  
	  hmclean_vs_npu_leadjet  -> Fill(npu,mclean->at(iclean));

	} // end loop over leading jet

      } // end loop over matched jets

    }// end loop over jets

    hnjets->Fill(nj);

    // -----  grooming is done only for jets  with pT>100 GeV
    for (unsigned int j = 0; j < msoftdrop_beta20->size(); j++){
    
      thispt = pt->at(j); // use pt 
      //float thispt = ptcorr->at(j); // use ptcorr 
      //float thispt = ptraw->at(j); // use pt raw
      
      if (thispt < minPt)  continue;
      if (thispt > maxPt)  continue;

      if (fabs(eta->at(j)) < minAbsEta) continue;
      if (fabs(eta->at(j)) > maxAbsEta) continue;

      int matchInd = -1;
      if (treetype_ != "gen")
	matchInd = RecoToGenMatching(eta->at(j), phi->at(j), etagen, phigen);
      //matchInd = imatch->at(j); // this is not working now for groomed quantities...

      float mtrim         = mtrim_Rtrim_020_Ptfrac_005->at(j);
      float mtrimsafe     = mtrimsafe_Rtrim_020_Ptfrac_005->at(j);
      float msoftdrop     = msoftdrop_beta20->at(j);
      float msoftdropsafe = msoftdropsafe_beta20->at(j);

      float genpt = -999;
      float genmtrim = -999;
      float genmtrimsafe = -999;
      float genmsoftdrop= -999;
      float genmsoftdropsafe= -999;
      if (matchInd > -1){
	genpt        = ptgen->at(matchInd);
	genmtrim     = mtrimgen->at(matchInd);
	genmtrimsafe = mtrimsafegen->at(matchInd);
	genmsoftdrop = msoftdropgen->at(matchInd);
	genmsoftdropsafe = msoftdropsafegen->at(matchInd);
      }

      hmtrim    -> Fill(mtrim);
      hmtrimsafe-> Fill(mtrimsafe);
      hmsoftdrop-> Fill(msoftdrop);
      hmsoftdropsafe-> Fill(msoftdropsafe);
      
      if (j == 0){
	hmtrim_leadjet    -> Fill(mtrim);
	hmtrimsafe_leadjet-> Fill(mtrimsafe);
	hmsoftdrop_leadjet-> Fill(msoftdrop);
	hmsoftdropsafe_leadjet-> Fill(msoftdropsafe);
      }

      // -- response plots
      if (matchInd >- 1){
	hmtrim_response     -> Fill(mtrim-genmtrim);
	hmtrimsafe_response -> Fill(mtrimsafe- genmtrimsafe);
	hmsoftdrop_response -> Fill(msoftdrop-genmsoftdrop);
	hmsoftdropsafe_response -> Fill(msoftdropsafe-genmsoftdropsafe);

	// 2d plots
	hmtrim_response_vs_pt     -> Fill(genpt,mtrim-genmtrim);
	hmtrimsafe_response_vs_pt -> Fill(genpt,mtrimsafe- genmtrimsafe);
	hmsoftdrop_response_vs_pt     -> Fill(genpt,msoftdrop-genmsoftdrop);
	hmsoftdropsafe_response_vs_pt -> Fill(genpt,msoftdropsafe- genmsoftdropsafe);

	hmtrim_response_vs_eta     -> Fill(eta->at(j),mtrim-genmtrim);
	hmtrimsafe_response_vs_eta -> Fill(eta->at(j),mtrimsafe- genmtrimsafe);
	hmsoftdrop_response_vs_eta     -> Fill(eta->at(j),msoftdrop-genmsoftdrop);
	hmsoftdropsafe_response_vs_eta -> Fill(eta->at(j),msoftdropsafe- genmsoftdropsafe);

	hmtrim_response_vs_npu     -> Fill(npu,mtrim-genmtrim);
	hmtrimsafe_response_vs_npu -> Fill(npu,mtrimsafe- genmtrimsafe);
	hmsoftdrop_response_vs_npu     -> Fill(npu,msoftdrop-genmsoftdrop);
	hmsoftdropsafe_response_vs_npu -> Fill(npu,msoftdropsafe- genmsoftdropsafe);

	if (j == 0){
	  hmtrim_response_leadjet     -> Fill(mtrim-genmtrim);
	  hmtrimsafe_response_leadjet -> Fill(mtrimsafe- genmtrimsafe);
	  hmsoftdrop_response_leadjet     -> Fill(msoftdrop-genmsoftdrop);
	  hmsoftdropsafe_response_leadjet -> Fill(msoftdropsafe- genmsoftdropsafe);
	}
      }
 
    }// end loop over jets 


  }

}


// --- Save histograms ---------------------------------------------------------------
void JetTreeAnalyzer::saveHistograms(TFile *outfile, std::string dir){

  std::cout << "Saving histograms ... " << std::endl;
  
  outfile->cd();
  TDirectory *thisdir = outfile->mkdir(dir.c_str());
  thisdir->cd();    // make the "thisdir" directory
  
  hnjets->Write();

  hptgen->Write();
  hptgen_pu->Write();
  hptgen_good->Write();

  hptraw->Write();
  hptraw_pu->Write();
  hptraw_good->Write();
  hptraw_response->Write();

  hpt->Write();
  hpt_pu->Write();
  hpt_good->Write();
  hpt_response->Write();

  hptcorr->Write();
  hptcorr_pu->Write();
  hptcorr_good->Write();
  hptcorr_response->Write();

  heta->Write();
  heta_pu->Write();
  heta_good->Write();

  hnpu->Write();
  hnpu_pu->Write();
  hnpu_good->Write();

  hnparticles->Write();
  hnneutrals->Write();
  hncharged->Write();

  htau21->Write();
  htau21_softdrop->Write(); 

  hmraw->Write();
  hmraw_response->Write();
  hm->Write();
  hm_response->Write();
  hmtrim->Write();
  hmtrim_response->Write();
  hmtrimsafe->Write();
  hmtrimsafe_response->Write();
  hmclean->Write();
  hmclean_response->Write();
  hmconst->Write();
  hmconst_response->Write();
  hmsoftdrop->Write();
  hmsoftdrop_response->Write();
  hmsoftdropsafe->Write();
  hmsoftdropsafe_response->Write();

  // leading jet
  hptraw_leadjet->Write();
  hptraw_pu_leadjet->Write();
  hptraw_good_leadjet->Write();
  hptraw_response_leadjet->Write();

  hpt_leadjet->Write();
  hpt_pu_leadjet->Write();
  hpt_good_leadjet->Write();
  hpt_response_leadjet->Write();
  hptcorr_leadjet->Write();
  hptcorr_pu_leadjet->Write();
  hptcorr_good_leadjet->Write();
  hptcorr_response_leadjet->Write();

  heta_leadjet->Write();
  heta_pu_leadjet->Write();
  heta_good_leadjet->Write();

  hnparticles_leadjet->Write();
  hnneutrals_leadjet->Write();
  hncharged_leadjet->Write();

  htau21_leadjet->Write();
  htau21_softdrop_leadjet->Write(); 

  hmraw_leadjet->Write();
  hmraw_response_leadjet->Write();
  hm_leadjet->Write();
  hm_response_leadjet->Write();
  hmtrim_leadjet->Write();
  hmtrim_response_leadjet->Write();
  hmtrimsafe_leadjet->Write();
  hmtrimsafe_response_leadjet->Write();
  hmclean_leadjet->Write();
  hmclean_response_leadjet->Write();
  hmconst_leadjet->Write();
  hmconst_response_leadjet->Write();
  hmsoftdrop_leadjet->Write();
  hmsoftdrop_response_leadjet->Write();
  hmsoftdropsafe_leadjet->Write();
  hmsoftdropsafe_response_leadjet->Write();

  // 2d
  hptraw_response_vs_pt->Write();
  hpt_response_vs_pt->Write();
  hptcorr_response_vs_pt->Write();
  hmraw_response_vs_pt->Write();
  hm_response_vs_pt->Write();
  hmtrim_response_vs_pt->Write();
  hmtrimsafe_response_vs_pt->Write();
  hmclean_response_vs_pt->Write();
  hmconst_response_vs_pt->Write();
  hmsoftdrop_response_vs_pt->Write();
  hmsoftdropsafe_response_vs_pt->Write();

  hptraw_response_vs_eta->Write();
  hpt_response_vs_eta->Write();
  hptcorr_response_vs_eta->Write();
  hm_response_vs_eta->Write();
  hmraw_response_vs_eta->Write();
  hmtrim_response_vs_eta->Write();
  hmtrimsafe_response_vs_eta->Write();
  hmclean_response_vs_eta->Write();
  hmconst_response_vs_eta->Write();
  hmsoftdrop_response_vs_eta->Write();
  hmsoftdropsafe_response_vs_eta->Write();

  hptraw_response_vs_npu->Write();
  hpt_response_vs_npu->Write();
  hptcorr_response_vs_npu->Write();
  hmraw_response_vs_npu->Write();
  hm_response_vs_npu->Write();
  hmtrim_response_vs_npu->Write();
  hmtrimsafe_response_vs_npu->Write();
  hmclean_response_vs_npu->Write();
  hmconst_response_vs_npu->Write();
  hmsoftdrop_response_vs_npu->Write();
  hmsoftdropsafe_response_vs_npu->Write();

}
