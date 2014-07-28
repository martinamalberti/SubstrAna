#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"
#include "TList.h"

#include "TKey.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TH2F.h"
#include "TMatrixDSym.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// fill a TList with all the methods                                                                                                                                                       
int GetListOfMethods( TList & methods, TDirectory *dir = 0);
int GetListOfTitles(  TDirectory*rfdir,TList & titles);
TKey *NextKey(TIter & keyIter, TString className);

struct MatrixEntry{

  std::string binXName ; 
  std::string binYName ;
  double value;

};


/// Main programme 
int main (int argc, char **argv){
  if(argc<2){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit "<<std::endl; return -1; }

  // Load TTree Lybrary                                                                                                                                                                   
  gSystem->Load("libTree.so");

  // Set Root style from global enviroment path                                                                                                                                           
  std::string ROOTStyle;
  if(getenv ("ROOTStyle")!=NULL){
    ROOTStyle = getenv ("ROOTStyle");
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetErrorX(0.5);
  gStyle->SetPaintTextFormat("3g");

  // parse config file parameter                                                                                                                                                          
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");


  std::vector<edm::ParameterSet> InputInformationParamLowPU ;
  if(Options.existsAs<std::vector<edm::ParameterSet>>("InputLowPUFiles"))
    InputInformationParamLowPU = Options.getParameter<std::vector<edm::ParameterSet>>("InputLowPUFiles");
  else{ std::cout<<" Exit from code, no input set found for low pile-up files"<<std::endl; return -1; }

  std::vector<edm::ParameterSet> InputInformationParamHighPU ;
  if(Options.existsAs<std::vector<edm::ParameterSet>>("InputHighPUFiles"))
    InputInformationParamHighPU = Options.getParameter<std::vector<edm::ParameterSet>>("InputHighPUFiles");
  else{ std::cout<<" Exit from code, no input set found for high pile-up files"<<std::endl; return -1; }


  std::string outputDirectory;
  if(Options.existsAs<std::string>("outputDirectory"))
    outputDirectory = Options.getParameter<std::string>("outputDirectory");
  else{ outputDirectory = "output/"; }

  double signalEfficiencyTarget;
  if(Options.existsAs<double>("signalEfficiencyTarget"))
    signalEfficiencyTarget = Options.getParameter<double>("signalEfficiencyTarget");
  else{ signalEfficiencyTarget = 0.5; }
 
  system(("mkdir -p "+outputDirectory).c_str());

  std::vector<MatrixEntry> performanceValue ;

  TFile* inputFile = NULL ;

  // take all the TH1 outputs for S and B for  each variable in input-> low Pile Up
  TList TrainingMethods;
  TList Titles;

  std::vector<edm::ParameterSet>::const_iterator itLowPileUp = InputInformationParamLowPU.begin();
  for( ; itLowPileUp != InputInformationParamLowPU.end() ; ++itLowPileUp){
    inputFile = TFile::Open((*itLowPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL){
      MatrixEntry entry ; 
      entry.binXName = (*itLowPileUp).getParameter<std::string>("variableNameX");
      entry.binYName = (*itLowPileUp).getParameter<std::string>("variableNameY");
      entry.value = 0;   
      performanceValue.push_back(entry);         
      continue;
    }
    // take the background efficiency related to the target
    inputFile->cd();
    TrainingMethods.Clear();
    int res = GetListOfMethods(TrainingMethods);
    if(res == 0) std::cout<<" No methods found "<<std::endl ;
    TIter next(&TrainingMethods);
    TKey *key = 0, *hkey = 0;
    while ((key = (TKey*)next())) {
      TDirectory *myDir = (TDirectory*)key->ReadObj();
      Titles.Clear();
      int nTitles = GetListOfTitles(myDir,Titles);
      if(nTitles == 0) std::cout<<" No titles found "<<std::endl ;
      TIter nextTitle(&Titles);
      TKey *titkey = 0;
      TDirectory *titDir = 0;
      while ((titkey = NextKey(nextTitle,"TDirectory"))) {
	titDir = (TDirectory*)titkey->ReadObj(); // read each object and take again the method title for each element of the list                                                         
	TString methodTitle;
        methodTitle = titDir->GetName();
	TIter nextKey( titDir->GetListOfKeys() ); // loop and the list of keys                                                                                                
	while ((hkey = NextKey(nextKey,"TH1"))) { // take only the TH1 object type                                                                                                   
	  TH1F* h = (TH1F*) hkey->ReadObj();
	  TString hname = h->GetName();    // only the one which are called rejBvsS            
	  if (hname.Contains("effBvsS") && hname.BeginsWith("MVA_") && not hname.Contains("effBvsSLocal")) {
            MatrixEntry entry ; 
            entry.binXName = (*itLowPileUp).getParameter<std::string>("variableNameX");
            entry.binYName = (*itLowPileUp).getParameter<std::string>("variableNameY");
	    entry.value = (1./h->GetBinContent(h->FindBin(signalEfficiencyTarget)));   
            performanceValue.push_back(entry);         
	  }
	}
      }
    }
  }

  //for getting the number of bins is enough to cycle on all the entry and count the diagonal terms
  int numberOfBins = 0;
  for(unsigned int ientry = 0 ; ientry < performanceValue.size(); ientry++){
    if(performanceValue.at(ientry).binXName == performanceValue.at(ientry).binYName) numberOfBins++;
  }

  TH2F* performanceBDT_lowPileUP  = new TH2F("backgroundMatrix_lowPileUP","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
  for(int iBinX = 0; iBinX < performanceBDT_lowPileUP->GetNbinsX(); iBinX ++){
   for(int iBinY = 0; iBinY < performanceBDT_lowPileUP->GetNbinsY(); iBinY ++){
     performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,0.);
   }
  }

  int vecPos = 0;
  for(int iBinX = 0; iBinX < performanceBDT_lowPileUP->GetNbinsX(); iBinX ++){  
   for(int iBinY = iBinX; iBinY < performanceBDT_lowPileUP->GetNbinsY(); iBinY ++){
     if(iBinY < (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_lowPileUP->GetNbinsX()-1) and vecPos < int(performanceValue.size()-1)){
       performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
       if(std::string(performanceBDT_lowPileUP->GetYaxis()->GetBinLabel(iBinY+1)) == "")
	 performanceBDT_lowPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       if(std::string(performanceBDT_lowPileUP->GetXaxis()->GetBinLabel(iBinY+1)) == "")
          performanceBDT_lowPileUP->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       vecPos++; 
     }    
     else if (iBinY == (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_lowPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY < (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_lowPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY == (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_lowPileUP->GetNbinsX()-1)){
       performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value));
       performanceBDT_lowPileUP->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_lowPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
     }
   }
  }

  TCanvas* cPerformance_lowPU = new TCanvas("","",180,52,500,550);

  cPerformance_lowPU->SetGrid();
  cPerformance_lowPU->SetTicks();
  cPerformance_lowPU->cd();

  performanceBDT_lowPileUP->SetMarkerSize(1.5);
  performanceBDT_lowPileUP->SetMarkerColor(0);
  performanceBDT_lowPileUP->GetXaxis()->SetLabelSize(0.035);
  performanceBDT_lowPileUP->GetYaxis()->SetLabelSize(0.035);
  performanceBDT_lowPileUP->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     
  performanceBDT_lowPileUP->SetMaximum(performanceBDT_lowPileUP->GetMaximum());
  performanceBDT_lowPileUP->SetMinimum(performanceBDT_lowPileUP->GetMinimum());

  performanceBDT_lowPileUP->Draw("colz");
  performanceBDT_lowPileUP->Draw("textsame");
  performanceBDT_lowPileUP->GetXaxis()->LabelsOption("v");
  performanceBDT_lowPileUP->GetYaxis()->LabelsOption("h");


  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(21); // align right                                                                                                                                                  
  latex.SetTextSize(0.033);
  latex.DrawLatex(0.547,0.92,Form("CMS Preliminary Simulation, #sqrt{s} = 13 TeV"));

  cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.pdf").c_str(),"pdf");
  cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.png").c_str(),"png");
  cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.root").c_str(),"root");

  cPerformance_lowPU->SetLogz();
 
  cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.pdf").c_str(),"pdf");
  cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.png").c_str(),"png");
  cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.root").c_str(),"root");

  // high pt part
  performanceValue.clear();
  std::vector<edm::ParameterSet>::const_iterator itHighPileUp = InputInformationParamHighPU.begin();
  for( ; itHighPileUp != InputInformationParamHighPU.end() ; ++itHighPileUp){
    inputFile = TFile::Open((*itHighPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL){
      MatrixEntry entry ; 
      entry.binXName = (*itHighPileUp).getParameter<std::string>("variableNameX");
      entry.binYName = (*itHighPileUp).getParameter<std::string>("variableNameY");
      entry.value = 0;   
      performanceValue.push_back(entry);         
      continue;
    }
    // take the background efficiency related to the target
    inputFile->cd();
    TrainingMethods.Clear();
    int res = GetListOfMethods(TrainingMethods);
    if(res == 0) std::cout<<" No methods found "<<std::endl ;
    TIter next(&TrainingMethods);
    TKey *key = 0, *hkey = 0;
    while ((key = (TKey*)next())) {
      TDirectory *myDir = (TDirectory*)key->ReadObj();
      Titles.Clear();
      int nTitles = GetListOfTitles(myDir,Titles);
      if(nTitles == 0) std::cout<<" No titles found "<<std::endl ;
      TIter nextTitle(&Titles);
      TKey *titkey = 0;
      TDirectory *titDir = 0;
      while ((titkey = NextKey(nextTitle,"TDirectory"))) {
	titDir = (TDirectory*)titkey->ReadObj(); // read each object and take again the method title for each element of the list                                                         
	TString methodTitle;
        methodTitle = titDir->GetName();
	TIter nextKey( titDir->GetListOfKeys() ); // loop and the list of keys                                                                                                
	while ((hkey = NextKey(nextKey,"TH1"))) { // take only the TH1 object type                                                                                                   
	  TH1F* h = (TH1F*) hkey->ReadObj();
	  TString hname = h->GetName();    // only the one which are called rejBvsS            
	  if (hname.Contains("effBvsS") && hname.BeginsWith("MVA_") && not hname.Contains("effBvsSLocal")) {
            MatrixEntry entry ; 
            entry.binXName = (*itHighPileUp).getParameter<std::string>("variableNameX");
            entry.binYName = (*itHighPileUp).getParameter<std::string>("variableNameY");
	    entry.value = (1./h->GetBinContent(h->FindBin(signalEfficiencyTarget)));   
            performanceValue.push_back(entry);         
	  }
	}
      }
    }
  }

  //for getting the number of bins is enough to cycle on all the entry and count the diagonal terms
  numberOfBins = 0;
  for(unsigned int ientry = 0 ; ientry < performanceValue.size(); ientry++){
    if(performanceValue.at(ientry).binXName == performanceValue.at(ientry).binYName) numberOfBins++;
  }

  TH2F* performanceBDT_highPileUP  = new TH2F("backgroundMatrix_highPileUP","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
  for(int iBinX = 0; iBinX < performanceBDT_highPileUP->GetNbinsX(); iBinX ++){
   for(int iBinY = 0; iBinY < performanceBDT_highPileUP->GetNbinsY(); iBinY ++){
     performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,0.);
   }
  }

  vecPos = 0;
  for(int iBinX = 0; iBinX < performanceBDT_highPileUP->GetNbinsX(); iBinX ++){
   for(int iBinY = iBinX; iBinY < performanceBDT_highPileUP->GetNbinsY(); iBinY ++){
     if(iBinY < (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_highPileUP->GetNbinsX()-1) and vecPos < int(performanceValue.size()-1)){
       performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
       if(std::string(performanceBDT_highPileUP->GetYaxis()->GetBinLabel(iBinY+1)) == "") 
          performanceBDT_highPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       if(std::string(performanceBDT_highPileUP->GetXaxis()->GetBinLabel(iBinY+1)) == "")
          performanceBDT_highPileUP->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       vecPos++;
     }    
     else if (iBinY == (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_highPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY < (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_highPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY == (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_highPileUP->GetNbinsX()-1)){
       performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value));
       performanceBDT_highPileUP->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_highPileUP->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
     }

   }
  }

  TCanvas* cPerformance_highPU = new TCanvas("","",180,52,500,550);

  cPerformance_highPU->SetGrid();
  cPerformance_highPU->SetTicks();
  cPerformance_highPU->cd();

  performanceBDT_highPileUP->SetMarkerSize(1.5);
  performanceBDT_highPileUP->SetMarkerColor(0);
  performanceBDT_highPileUP->GetXaxis()->SetLabelSize(0.035);
  performanceBDT_highPileUP->GetYaxis()->SetLabelSize(0.035);
  performanceBDT_highPileUP->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     
  performanceBDT_highPileUP->SetMaximum(performanceBDT_highPileUP->GetMaximum());
  performanceBDT_highPileUP->SetMinimum(performanceBDT_highPileUP->GetMinimum());

  performanceBDT_highPileUP->Draw("colz");
  performanceBDT_highPileUP->Draw("textsame");
  performanceBDT_highPileUP->GetXaxis()->LabelsOption("v");
  performanceBDT_highPileUP->GetYaxis()->LabelsOption("h");

  latex.DrawLatex(0.547,0.92,Form("CMS Preliminary Simulation, #sqrt{s} = 13 TeV"));

  cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.pdf").c_str(),"pdf");
  cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.png").c_str(),"png");
  cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.root").c_str(),"root");

  cPerformance_highPU->SetLogz();
 
  cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.pdf").c_str(),"pdf");
  cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.png").c_str(),"png");
  cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.root").c_str(),"root");
  
  return 0 ;

}

//////////////////////////////

int GetListOfMethods( TList & methods, TDirectory *dir){

  if (dir==0) dir = gDirectory;
  TIter mnext(dir->GetListOfKeys());
  TKey *mkey;
  methods.Clear();
  methods.SetOwner(kFALSE);
  UInt_t ni=0;
  while ((mkey = (TKey*)mnext())) { // make sure, that we only look at TDirectory with name Method_<xxx>                                                                      
    TString name = mkey->GetClassName();
    TClass *cl = gROOT->GetClass(name);
    if (cl->InheritsFrom("TDirectory")) {
      if (TString(mkey->GetName()).BeginsWith("Method_")) {
	methods.Add(mkey);
	ni++;
      }
    }
  }
  return ni;
}

// get a list of titles (i.e TDirectory) given a method dir                                                                                                           
int GetListOfTitles( TDirectory *rfdir, TList & titles ){
 UInt_t ni=0;
 if (rfdir==0) return 0;
 TList *keys = rfdir->GetListOfKeys();
 if(keys==0) {
    std::cout << "+++ Directory '" << rfdir->GetName() << "' contains no keys" << std::endl;
    return 0;
 }

 TIter rfnext(rfdir->GetListOfKeys());
 TKey *rfkey;
 titles.Clear();
 titles.SetOwner(kFALSE);
 while ((rfkey = (TKey*)rfnext())) { // make sure, that we only look at histograms                                                                                                       
   TClass *cl = gROOT->GetClass(rfkey->GetClassName());
   if (cl->InheritsFrom("TDirectory")) {
   titles.Add(rfkey);
   ni++;
   }
 }
 return ni;
}

// Next key iterator matching the className                                                                                                                    
TKey *NextKey(TIter & keyIter, TString className) {
 TKey *key  = (TKey *) keyIter.Next();
 TKey *rkey = 0;
 Bool_t loop = (key!=0);

 while (loop) {
  TClass *cl = gROOT->GetClass(key->GetClassName());
  if (cl->InheritsFrom(className.Data())) { loop = kFALSE;
       rkey = key;
  }
  else {
     key = (TKey *)keyIter.Next();
     if (key==0) loop = kFALSE;
  }
 }
 return rkey;

}
