#ifndef JetTreeAnalyzer_h
#define JetTreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TStyle.h"

#include "TMath.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

using namespace std;

class JetTreeAnalyzer{

 public :

  JetTreeAnalyzer(TTree *tree, TTree *gentree, string treetype);
  virtual ~JetTreeAnalyzer();

  virtual void Init(TTree *tree, TTree *gentree);
  virtual int  GetEntry(Long64_t entry);

  virtual void bookHistograms(std::string suffix="");
  virtual void fillHistograms(int maxEntries, float minPt, float maxPt, float minAbsEta, float maxAbsEta);
  virtual void saveHistograms(TFile *file, std::string dir);

  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain                                              
  TTree          *fChain2;  //!pointer to the analyzed TTree or TChain (gen)
  Int_t           fCurrent; //!current Tree number in a TChain
  Int_t           fCurrent2; //!current Tree number in a TChain

  // Declaration of leaf types
  int npu;
  int npv;
  vector<float>   *pt;
  vector<float>   *ptcorr;
  vector<float>   *ptraw;
  vector<float>   *ptunc;
  vector<float>   *eta;
  vector<float>   *phi;
  vector<float>   *m;
  vector<float>   *mraw;
  vector<float>   *ptclean;
  vector<float>   *mclean;
  vector<float>   *pttrim_Rtrim_020_Ptfrac_005;
  vector<float>   *mtrim_Rtrim_020_Ptfrac_005;
  vector<float>   *pttrimsafe_Rtrim_020_Ptfrac_005;
  vector<float>   *mtrimsafe_Rtrim_020_Ptfrac_005;
  vector<float>   *pttrim_Rtrim_010_Ptfrac_003;
  vector<float>   *mtrim_Rtrim_010_Ptfrac_003;
  vector<float>   *pttrimsafe_Rtrim_010_Ptfrac_003;
  vector<float>   *mtrimsafe_Rtrim_010_Ptfrac_003;
  vector<float>   *pttrim_Rtrim_020_Ptfrac_003;
  vector<float>   *mtrim_Rtrim_020_Ptfrac_003;
  vector<float>   *pttrimsafe_Rtrim_020_Ptfrac_003;
  vector<float>   *mtrimsafe_Rtrim_020_Ptfrac_003;
  vector<float>   *pttrim_Rtrim_030_Ptfrac_003;
  vector<float>   *mtrim_Rtrim_030_Ptfrac_003;
  vector<float>   *pttrimsafe_Rtrim_030_Ptfrac_003;
  vector<float>   *mtrimsafe_Rtrim_030_Ptfrac_003;
  vector<float>   *ptconst;
  vector<float>   *mconst;
  vector<float>   *ptpruned_zcut_010_R_cut_050;
  vector<float>   *mpruned_zcut_010_R_cut_050;
  vector<float>   *ptprunedsafe_zcut_010_R_cut_050;
  vector<float>   *mprunedsafe_zcut_010_R_cut_050;
  vector<float>   *QGLikelihood_pr_zcut_010_R_cut_050;
  vector<float>   *QGLikelihood_pr_sub1_zcut_010_R_cut_050;
  vector<float>   *QGLikelihood_pr_sub2_zcut_010_R_cut_050;
  vector<float>   *ptpruned_zcut_005_R_cut_050;
  vector<float>   *mpruned_zcut_005_R_cut_050;
  vector<float>   *ptprunedsafe_zcut_005_R_cut_050;
  vector<float>   *mprunedsafe_zcut_005_R_cut_050;
  vector<float>   *QGLikelihood_pr_zcut_005_R_cut_050;
  vector<float>   *QGLikelihood_pr_sub1_zcut_005_R_cut_050;
  vector<float>   *QGLikelihood_pr_sub2_zcut_005_R_cut_050;
  vector<float>   *ptpruned_zcut_005_R_cut_075;
  vector<float>   *mpruned_zcut_005_R_cut_075;
  vector<float>   *ptprunedsafe_zcut_005_R_cut_075;
  vector<float>   *mprunedsafe_zcut_005_R_cut_075;
  vector<float>   *QGLikelihood_pr_zcut_005_R_cut_075;
  vector<float>   *QGLikelihood_pr_sub1_zcut_005_R_cut_075;
  vector<float>   *QGLikelihood_pr_sub2_zcut_005_R_cut_075;
  vector<float>   *ptpruned_zcut_010_R_cut_075;
  vector<float>   *mpruned_zcut_010_R_cut_075;
  vector<float>   *ptprunedsafe_zcut_010_R_cut_075;
  vector<float>   *mprunedsafe_zcut_010_R_cut_075;
  vector<float>   *QGLikelihood_pr_zcut_010_R_cut_075;
  vector<float>   *QGLikelihood_pr_sub1_zcut_010_R_cut_075;
  vector<float>   *QGLikelihood_pr_sub2_zcut_010_R_cut_075;
  vector<float>   *ptsoftdrop_beta20;
  vector<float>   *msoftdrop_beta20;
  vector<float>   *ptsoftdropsafe_beta20;
  vector<float>   *msoftdropsafe_beta20;
  vector<float>   *ptsoftdrop_beta00;
  vector<float>   *msoftdrop_beta00;
  vector<float>   *ptsoftdropsafe_beta00;
  vector<float>   *msoftdropsafe_beta00;
  vector<float>   *ptsoftdrop_beta10;
  vector<float>   *msoftdrop_beta10;
  vector<float>   *ptsoftdropsafe_beta10;
  vector<float>   *msoftdropsafe_beta10;
  vector<float>   *ptsoftdrop_betam10;
  vector<float>   *msoftdrop_betam10;
  vector<float>   *ptsoftdropsafe_betam10;
  vector<float>   *msoftdropsafe_betam10;
  vector<int>     *nparticles;
  vector<int>     *nneutrals;
  vector<int>     *ncharged;
  vector<float>   *sdsymmetry;
  vector<float>   *sddeltar;
  vector<float>   *sdmu;
  vector<float>   *sdenergyloss;
  vector<float>   *sdarea;
  vector<float>   *sdnconst;
  vector<float>   *mfiltsoftdrop;
  vector<float>   *tau1;
  vector<float>   *tau2;
  vector<float>   *tau3;
  vector<float>   *tau4;
  vector<float>   *tau5;
  vector<float>   *tau1_pr;
  vector<float>   *tau2_pr;
  vector<float>   *tau3_pr;
  vector<float>   *tau4_pr;
  vector<float>   *tau5_pr;
  vector<float>   *tau1_softdrop;
  vector<float>   *tau2_softdrop;
  vector<float>   *tau3_softdrop;
  vector<float>   *tau4_softdrop;
  vector<float>   *tau5_softdrop;
  vector<float>   *Qjets;
  vector<float>   *charge_k05;
  vector<float>   *charge_k07;
  vector<float>   *charge_k10;
  vector<float>   *ecf_beta_05;
  vector<float>   *ecf_beta_10;
  vector<float>   *ecf_beta_15;
  vector<float>   *ecf_beta_20;
  vector<float>   *hepmass;
  vector<float>   *hepwmass;
  vector<float>   *hepm01;
  vector<float>   *hepm02;
  vector<float>   *hepm12;
  vector<float>   *hepm12m012;
  vector<float>   *hepatanm02m01;
  vector<float>   *cmsmass;
  vector<float>   *cmsminmass;
  vector<float>   *cmshelicity;
  vector<float>   *cmsnsubjets;
  vector<float>   *ptgen;
  vector<float>   *etagen;
  vector<float>   *phigen;
  vector<float>   *mgen;
  vector<float>   *mrawgen;
  vector<float>   *mtrimgen;
  vector<float>   *mtrimsafegen;
  vector<float>   *mcleangen;
  vector<float>   *mconstgen;
  vector<float>   *pttrimgen;
  vector<float>   *pttrimsafegen;
  vector<float>   *ptcleangen;
  vector<float>   *ptconstgen;
  vector<int>     *imatch;
  vector<int>     *flavourgen;
  vector<float>   *msoftdropgen;
  vector<float>   *msoftdropsafegen;
  vector<float>   *mfiltsoftdropgen;
  vector<bool>    *is_MatchedToBoson;

  // List of branches
  TBranch        *b_npu;
  TBranch        *b_npv;
  TBranch        *b_pt;
  TBranch        *b_ptcorr;
  TBranch        *b_ptraw;
  TBranch        *b_ptunc;
  TBranch        *b_eta;
  TBranch        *b_phi;
  TBranch        *b_m;
  TBranch        *b_mraw;
  TBranch        *b_ptclean;
  TBranch        *b_mclean;
  TBranch        *b_pttrim_Rtrim_020_Ptfrac_005;
  TBranch        *b_mtrim_Rtrim_020_Ptfrac_005;
  TBranch        *b_pttrimsafe_Rtrim_020_Ptfrac_005;
  TBranch        *b_mtrimsafe_Rtrim_020_Ptfrac_005;
  TBranch        *b_pttrim_Rtrim_010_Ptfrac_003;
  TBranch        *b_mtrim_Rtrim_010_Ptfrac_003;
  TBranch        *b_pttrimsafe_Rtrim_010_Ptfrac_003;
  TBranch        *b_mtrimsafe_Rtrim_010_Ptfrac_003;
  TBranch        *b_pttrim_Rtrim_020_Ptfrac_003;
  TBranch        *b_mtrim_Rtrim_020_Ptfrac_003;
  TBranch        *b_mtrimsafe_Rtrim_030_Ptfrac_003;
  TBranch        *b_ptconst;
  TBranch        *b_mconst;
  TBranch        *b_ptpruned_zcut_010_R_cut_050;
  TBranch        *b_mpruned_zcut_010_R_cut_050;
  TBranch        *b_ptprunedsafe_zcut_010_R_cut_050;
  TBranch        *b_mprunedsafe_zcut_010_R_cut_050;
  TBranch        *b_QGLikelihood_pr_zcut_010_R_cut_050;
  TBranch        *b_QGLikelihood_pr_sub1_zcut_010_R_cut_050;
  TBranch        *b_QGLikelihood_pr_sub2_zcut_010_R_cut_050;
  TBranch        *b_ptpruned_zcut_005_R_cut_050;
  TBranch        *b_mpruned_zcut_005_R_cut_050;
  TBranch        *b_ptprunedsafe_zcut_005_R_cut_050;
  TBranch        *b_mprunedsafe_zcut_005_R_cut_050;
  TBranch        *b_QGLikelihood_pr_zcut_005_R_cut_050;
  TBranch        *b_QGLikelihood_pr_sub1_zcut_005_R_cut_050;
  TBranch        *b_QGLikelihood_pr_sub2_zcut_005_R_cut_050;
  TBranch        *b_ptpruned_zcut_005_R_cut_075;
  TBranch        *b_mpruned_zcut_005_R_cut_075;
  TBranch        *b_ptprunedsafe_zcut_005_R_cut_075;
  TBranch        *b_mprunedsafe_zcut_005_R_cut_075;
  TBranch        *b_QGLikelihood_pr_zcut_005_R_cut_075;
  TBranch        *b_QGLikelihood_pr_sub1_zcut_005_R_cut_075;
  TBranch        *b_QGLikelihood_pr_sub2_zcut_010_R_cut_075;
  TBranch        *b_ptsoftdrop_beta20;
  TBranch        *b_msoftdrop_beta20;
  TBranch        *b_ptsoftdropsafe_beta20;
  TBranch        *b_msoftdropsafe_beta20;
  TBranch        *b_ptsoftdrop_beta00;
  TBranch        *b_msoftdrop_beta00;
  TBranch        *b_ptsoftdropsafe_beta00;
  TBranch        *b_msoftdropsafe_beta00;
  TBranch        *b_ptsoftdrop_beta10;
  TBranch        *b_msoftdrop_beta10;
  TBranch        *b_ptsoftdropsafe_beta10;
  TBranch        *b_msoftdropsafe_beta10;
  TBranch        *b_ptsoftdrop_betam10;
  TBranch        *b_msoftdrop_betam10;
  TBranch        *b_ptsoftdropsafe_betam10;
  TBranch        *b_msoftdropsafe_betam10;
  TBranch        *b_nparticles;
  TBranch        *b_nneutrals;
  TBranch        *b_ncharged;
  TBranch        *b_sdsymmetry;
  TBranch        *b_sddeltar;
  TBranch        *b_sdmu;
  TBranch        *b_sdenergyloss;
  TBranch        *b_sdarea;
  TBranch        *b_sdnconst;
  TBranch        *b_mfiltsoftdrop;
  TBranch        *b_tau1;
  TBranch        *b_tau2;
  TBranch        *b_tau3;
  TBranch        *b_tau4;
  TBranch        *b_tau5;
  TBranch        *b_tau1_pr;
  TBranch        *b_tau2_pr;
  TBranch        *b_tau3_pr;
  TBranch        *b_tau4_pr;
  TBranch        *b_tau5_pr;
  TBranch        *b_tau1_softdrop;
  TBranch        *b_tau2_softdrop;
  TBranch        *b_tau3_softdrop;
  TBranch        *b_tau4_softdrop;
  TBranch        *b_tau5_softdrop;
  TBranch        *b_Qjets;
  TBranch        *b_charge_k05;
  TBranch        *b_charge_k07;
  TBranch        *b_charge_k10;
  TBranch        *b_ecf_beta_05;
  TBranch        *b_ecf_beta_10;
  TBranch        *b_ecf_beta_15;
  TBranch        *b_ecf_beta_20;
  TBranch        *b_hepmass;
  TBranch        *b_hepwmass;
  TBranch        *b_hepm01;
  TBranch        *b_hepm02;
  TBranch        *b_hepm12;
  TBranch        *b_hepm12m012;
  TBranch        *b_hepatanm02m01;
  TBranch        *b_cmsmass;
  TBranch        *b_cmsminmass;
  TBranch        *b_cmshelicity;
  TBranch        *b_cmsnsubjets;
  TBranch        *b_ptgen;
  TBranch        *b_etagen;
  TBranch        *b_phigen;
  TBranch        *b_mgen;
  TBranch        *b_mrawgen;
  TBranch        *b_mtrimgen;
  TBranch        *b_mtrimsafegen;
  TBranch        *b_mcleangen;
  TBranch        *b_mconstgen;
  TBranch        *b_pttrimgen;
  TBranch        *b_pttrimsafegen;
  TBranch        *b_ptcleangen;
  TBranch        *b_ptconstgen;
  TBranch        *b_imatch;
  TBranch        *b_flavourgen;
  TBranch        *b_msoftdropgen;
  TBranch        *b_msoftdropsafegen;
  TBranch        *b_mfiltsoftdropgen;
  TBranch        *b_is_MatchedToBoson;


  // histograms declaration
  TH1F *hnjets;

  TH1F* hptgen;
  TH1F* hptgen_pu;
  TH1F* hptgen_good;

  TH1F* hptraw;
  TH1F* hptraw_pu;
  TH1F* hptraw_good;
  TH1F* hptraw_response;
  
  TH1F* hpt;
  TH1F* hpt_pu;
  TH1F* hpt_good;
  TH1F* hpt_response;

  TH1F* hptcorr;
  TH1F* hptcorr_pu;
  TH1F* hptcorr_good;
  TH1F* hptcorr_response;

  TH1F* hpttrim;
  TH1F* hpttrim_pu;
  TH1F* hpttrim_good;
  TH1F* hpttrim_response;

  TH1F* hpttrimsafe;
  TH1F* hpttrimsafe_pu;
  TH1F* hpttrimsafe_good;
  TH1F* hpttrimsafe_response;

  TH1F* hptsoftdrop;
  TH1F* hptsoftdrop_pu;
  TH1F* hptsoftdrop_good;
  TH1F* hptsoftdrop_response;

  TH1F* hptsoftdropsafe;
  TH1F* hptsoftdropsafe_pu;
  TH1F* hptsoftdropsafe_good;
  TH1F* hptsoftdropsafe_response;

  TH1F* hptconst;
  TH1F* hptconst_pu;
  TH1F* hptconst_good;
  TH1F* hptconst_response;

  TH1F* hptclean;
  TH1F* hptclean_pu;
  TH1F* hptclean_good;
  TH1F* hptclean_response;

  TH1F* heta;
  TH1F* heta_pu;
  TH1F* heta_good;

  TH1F* hnpu;
  TH1F* hnpu_pu;
  TH1F* hnpu_good;

  TH1F* hmraw;
  TH1F* hmraw_response;

  TH1F* hm;
  TH1F* hm_response;  

  TH1F* hmtrim;
  TH1F* hmtrim_response;

  TH1F* hmtrimsafe;
  TH1F* hmtrimsafe_response;

  TH1F* hmsoftdrop;
  TH1F* hmsoftdrop_response;

  TH1F* hmsoftdropsafe;
  TH1F* hmsoftdropsafe_response;

  TH1F* hmconst;
  TH1F* hmconst_response;

  TH1F* hmclean;
  TH1F* hmclean_response;

  TH1F* hnparticles;
  TH1F* hnneutrals;
  TH1F* hncharged;

  TH1F* htau21;
  TH1F* htau21_softdrop;

  TH2F* hpt_response_vs_pt;
  TH2F* hptraw_response_vs_pt;
  TH2F* hptcorr_response_vs_pt;
  TH2F* hmraw_response_vs_pt;
  TH2F* hm_response_vs_pt;
  TH2F* hmtrim_response_vs_pt;
  TH2F* hmtrimsafe_response_vs_pt;
  TH2F* hmclean_response_vs_pt;
  TH2F* hmconst_response_vs_pt;
  TH2F* hmsoftdrop_response_vs_pt;
  TH2F* hmsoftdropsafe_response_vs_pt;

  TH2F* hpt_response_vs_eta;
  TH2F* hptraw_response_vs_eta;
  TH2F* hptcorr_response_vs_eta;
  TH2F* hmraw_response_vs_eta;
  TH2F* hm_response_vs_eta;
  TH2F* hmtrim_response_vs_eta;
  TH2F* hmtrimsafe_response_vs_eta;
  TH2F* hmclean_response_vs_eta;
  TH2F* hmconst_response_vs_eta;
  TH2F* hmsoftdrop_response_vs_eta;
  TH2F* hmsoftdropsafe_response_vs_eta;

  TH2F* hpt_response_vs_npu;
  TH2F* hptraw_response_vs_npu;
  TH2F* hptcorr_response_vs_npu;
  TH2F* hmraw_response_vs_npu;
  TH2F* hm_response_vs_npu;
  TH2F* hmtrim_response_vs_npu;
  TH2F* hmtrimsafe_response_vs_npu;
  TH2F* hmclean_response_vs_npu;
  TH2F* hmconst_response_vs_npu;
  TH2F* hmsoftdrop_response_vs_npu;
  TH2F* hmsoftdropsafe_response_vs_npu;

  TH2F* hpt_vs_npu;
  TH2F* hptraw_vs_npu;
  TH2F* hptcorr_vs_npu;
  TH2F* hmraw_vs_npu;
  TH2F* hm_vs_npu;
  TH2F* hmtrim_vs_npu;
  TH2F* hmtrimsafe_vs_npu;
  TH2F* hmclean_vs_npu;
  TH2F* hmconst_vs_npu;
  TH2F* hmsoftdrop_vs_npu;
  TH2F* hmsoftdropsafe_vs_npu;
  TH2F* htau21_vs_npu;
  TH2F* htau21_softdrop_vs_npu;

  // leading jet
  TH1F* hptraw_leadjet;
  TH1F* hptraw_pu_leadjet;
  TH1F* hptraw_good_leadjet;
  TH1F* hptraw_response_leadjet;

  TH1F* hpt_leadjet;
  TH1F* hpt_pu_leadjet;
  TH1F* hpt_good_leadjet;
  TH1F* hpt_response_leadjet;

  TH1F* hptcorr_leadjet;
  TH1F* hptcorr_pu_leadjet;
  TH1F* hptcorr_good_leadjet;
  TH1F* hptcorr_response_leadjet;

  TH1F* hpttrim_leadjet;
  TH1F* hpttrim_pu_leadjet;
  TH1F* hpttrim_good_leadjet;
  TH1F* hpttrim_response_leadjet;

  TH1F* hpttrimsafe_leadjet;
  TH1F* hpttrimsafe_pu_leadjet;
  TH1F* hpttrimsafe_good_leadjet;
  TH1F* hpttrimsafe_response_leadjet;

  TH1F* hptsoftdrop_leadjet;
  TH1F* hptsoftdrop_pu_leadjet;
  TH1F* hptsoftdrop_good_leadjet;
  TH1F* hptsoftdrop_response_leadjet;

  TH1F* hptsoftdropsafe_leadjet;
  TH1F* hptsoftdropsafe_pu_leadjet;
  TH1F* hptsoftdropsafe_good_leadjet;
  TH1F* hptsoftdropsafe_response_leadjet;

  TH1F* hptclean_leadjet;
  TH1F* hptclean_pu_leadjet;
  TH1F* hptclean_good_leadjet;
  TH1F* hptclean_response_leadjet;

  TH1F* hptconst_leadjet;
  TH1F* hptconst_pu_leadjet;
  TH1F* hptconst_good_leadjet;
  TH1F* hptconst_response_leadjet;
 
  TH1F* heta_leadjet;
  TH1F* heta_pu_leadjet;
  TH1F* heta_good_leadjet;

  TH1F* hmraw_leadjet;
  TH1F* hmraw_response_leadjet;

  TH1F* hm_leadjet;
  TH1F* hm_response_leadjet;

  TH1F* hmtrim_leadjet;
  TH1F* hmtrim_response_leadjet;

  TH1F* hmtrimsafe_leadjet;
  TH1F* hmtrimsafe_response_leadjet;

  //TH1F* hmclean_leadjet[7];
  //TH1F* hmclean_response_leadjet[7];
  TH1F* hmclean_leadjet;
  TH1F* hmclean_response_leadjet;

  TH1F* hmconst_leadjet;
  TH1F* hmconst_response_leadjet;

  TH1F* hmsoftdrop_leadjet;
  TH1F* hmsoftdrop_response_leadjet;

  TH1F* hmsoftdropsafe_leadjet;
  TH1F* hmsoftdropsafe_response_leadjet;

  TH2F* hpt_response_vs_pt_leadjet;
  TH2F* hptraw_response_vs_pt_leadjet;
  TH2F* hptcorr_response_vs_pt_leadjet;
  TH2F* hmraw_response_vs_pt_leadjet;
  TH2F* hm_response_vs_pt_leadjet;
  TH2F* hmtrim_response_vs_pt_leadjet;
  TH2F* hmtrimsafe_response_vs_pt_leadjet;
  TH2F* hmclean_response_vs_pt_leadjet;
  TH2F* hmconst_response_vs_pt_leadjet;
  TH2F* hmsoftdrop_response_vs_pt_leadjet;
  TH2F* hmsoftdropsafe_response_vs_pt_leadjet;

  TH2F* hpt_response_vs_eta_leadjet;
  TH2F* hptraw_response_vs_eta_leadjet;
  TH2F* hptcorr_response_vs_eta_leadjet;
  TH2F* hmraw_response_vs_eta_leadjet;
  TH2F* hm_response_vs_eta_leadjet;
  TH2F* hmtrim_response_vs_eta_leadjet;
  TH2F* hmtrimsafe_response_vs_eta_leadjet;
  TH2F* hmclean_response_vs_eta_leadjet;
  TH2F* hmconst_response_vs_eta_leadjet;
  TH2F* hmsoftdrop_response_vs_eta_leadjet;
  TH2F* hmsoftdropsafe_response_vs_eta_leadjet;

  TH2F* hpt_response_vs_npu_leadjet;
  TH2F* hptraw_response_vs_npu_leadjet;
  TH2F* hptcorr_response_vs_npu_leadjet;
  TH2F* hmraw_response_vs_npu_leadjet;
  TH2F* hm_response_vs_npu_leadjet;
  TH2F* hmtrim_response_vs_npu_leadjet;
  TH2F* hmtrimsafe_response_vs_npu_leadjet;
  TH2F* hmclean_response_vs_npu_leadjet;
  TH2F* hmconst_response_vs_npu_leadjet;
  TH2F* hmsoftdrop_response_vs_npu_leadjet;
  TH2F* hmsoftdropsafe_response_vs_npu_leadjet;

  TH2F* hpt_vs_npu_leadjet;
  TH2F* hptraw_vs_npu_leadjet;
  TH2F* hptcorr_vs_npu_leadjet;
  TH2F* hmraw_vs_npu_leadjet;
  TH2F* hm_vs_npu_leadjet;
  TH2F* hmtrim_vs_npu_leadjet;
  TH2F* hmtrimsafe_vs_npu_leadjet;
  TH2F* hmclean_vs_npu_leadjet;
  TH2F* hmconst_vs_npu_leadjet;
  TH2F* hmsoftdrop_vs_npu_leadjet;
  TH2F* hmsoftdropsafe_vs_npu_leadjet;
  TH2F* htau21_vs_npu_leadjet;
  TH2F* htau21_softdrop_vs_npu_leadjet;

  TH1F* hnparticles_leadjet;
  TH1F* hnneutrals_leadjet;
  TH1F* hncharged_leadjet;

  TH1F* htau21_leadjet;
  TH1F* htau21_softdrop_leadjet;

 private:
  string treetype_;
  
};
#endif
