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
  gStyle->SetErrorX(0.5);


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
 
  system(("mkdir -p "+outputDirectory).c_str());

  std::vector<TTree*> inputLowPileUpTrees ;
  std::vector<std::string> reducedNameLowPileUp ;
  TFile* inputFile = NULL ;
  // take all the TH1 outputs for S and B for  each variable in input-> low Pile Up
  std::vector<edm::ParameterSet>::const_iterator itLowPileUp = InputInformationParamLowPU.begin();
  for( ; itLowPileUp != InputInformationParamLowPU.end() ; ++itLowPileUp){
    inputFile = TFile::Open((*itLowPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL) continue;
    inputLowPileUpTrees.push_back((TTree*) inputFile->Get("TestTree"));
    reducedNameLowPileUp.push_back((*itLowPileUp).getParameter<std::string>("variableName"));
  }

  // Build Up the correlation matrix
  TMatrixDSym correlationMatrixS_lowPileUP(reducedNameLowPileUp.size()) ;
  TMatrixDSym correlationMatrixB_lowPileUP(reducedNameLowPileUp.size()) ;
  for( unsigned int irow = 0 ; irow < reducedNameLowPileUp.size() ; irow++){
   for( unsigned int icolum = 0 ; icolum < reducedNameLowPileUp.size() ; icolum++){
     correlationMatrixS_lowPileUP(irow,icolum) = 0;
     correlationMatrixB_lowPileUP(irow,icolum) = 0;
   }
  }

  std::vector<float> MeanValue_S ; MeanValue_S.resize(inputLowPileUpTrees.size());
  std::vector<float> MeanValue_B ; MeanValue_B.resize(inputLowPileUpTrees.size());

  std::vector<int> ientrySignal ; ientrySignal.resize(inputLowPileUpTrees.size());
  std::vector<int> ientryBackground ; ientryBackground.resize(inputLowPileUpTrees.size());

  int classID          = 0;
  float BDTG_NoPruning = 0;

  // compute mean values
  for(unsigned int iTree = 0 ; iTree<inputLowPileUpTrees.size(); iTree++){   
    inputLowPileUpTrees.at(iTree)->SetBranchAddress("classID",&classID);
    inputLowPileUpTrees.at(iTree)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning);
    for(int iEntries = 0; iEntries < inputLowPileUpTrees.at(iTree)->GetEntries(); iEntries++){
      inputLowPileUpTrees.at(iTree)->GetEntry(iEntries);
      if(classID == 0){ MeanValue_S.at(iTree) += BDTG_NoPruning; ientrySignal.at(iTree)++;}
      else{ MeanValue_B.at(iTree) += BDTG_NoPruning; ientryBackground.at(iTree)++;}
   }
  }

  for( unsigned int iVec = 0 ; iVec < MeanValue_S.size(); iVec++){
    MeanValue_S.at(iVec) = MeanValue_S.at(iVec)/ientrySignal.at(iVec);
    MeanValue_B.at(iVec) = MeanValue_B.at(iVec)/ientryBackground.at(iVec);
  }

  // compute distances
  int classID_X          = 0;
  float BDTG_NoPruning_X = 0;
  int classID_Y          = 0;
  float BDTG_NoPruning_Y = 0;

  for(unsigned int iTreeX = 0 ; iTreeX<inputLowPileUpTrees.size(); iTreeX++){   
    for(unsigned int iTreeY = 0 ; iTreeY<inputLowPileUpTrees.size(); iTreeY++){    
     int iEntriesX = 0; int iEntriesY = 0;
     for( ; iEntriesX < inputLowPileUpTrees.at(iTreeX)->GetEntries() and iEntriesY < inputLowPileUpTrees.at(iTreeY)->GetEntries(); iEntriesX++, iEntriesY++){
       inputLowPileUpTrees.at(iTreeX)->SetBranchAddress("classID",&classID_X);
       inputLowPileUpTrees.at(iTreeX)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning_X);
       inputLowPileUpTrees.at(iTreeX)->GetEntry(iEntriesX);
       inputLowPileUpTrees.at(iTreeY)->SetBranchAddress("classID",&classID_Y);
       inputLowPileUpTrees.at(iTreeY)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning_Y);
       inputLowPileUpTrees.at(iTreeY)->GetEntry(iEntriesY);
       if(classID_X == 0 and classID_Y == 0)
         correlationMatrixS_lowPileUP(iTreeX,iTreeY) += double((BDTG_NoPruning_X-MeanValue_S.at(iTreeX))*(BDTG_NoPruning_Y-MeanValue_S.at(iTreeY)));
       else if(classID_X == 1 and classID_Y == 1) correlationMatrixB_lowPileUP(iTreeX,iTreeY) += double((BDTG_NoPruning_X-MeanValue_B.at(iTreeX))*(BDTG_NoPruning_Y-MeanValue_B.at(iTreeY)));
     }
    }
  }
  
  for( unsigned int binX = 0; binX < reducedNameLowPileUp.size(); binX ++){
   for( unsigned int binY = 0; binY < reducedNameLowPileUp.size(); binY ++){
     correlationMatrixS_lowPileUP(binX,binY) /= double(ientrySignal.at(binX));
     correlationMatrixB_lowPileUP(binX,binY) /= double(ientryBackground.at(binX));
   }
  }
 
  std::vector<float> Sigma_S ; Sigma_S.resize(inputLowPileUpTrees.size());
  std::vector<float> Sigma_B ; Sigma_B.resize(inputLowPileUpTrees.size());
  for( unsigned int binX = 0; binX < reducedNameLowPileUp.size(); binX ++) {
    Sigma_S.at(binX) = sqrt(correlationMatrixS_lowPileUP(binX,binX));
    Sigma_B.at(binX) = sqrt(correlationMatrixB_lowPileUP(binX,binX));
  }
  
  for( unsigned int binX = 0; binX < reducedNameLowPileUp.size(); binX ++){
   for( unsigned int binY = 0; binY < reducedNameLowPileUp.size(); binY ++){
     correlationMatrixS_lowPileUP(binX,binY) /= double(Sigma_S.at(binX)*Sigma_S.at(binY));
     correlationMatrixB_lowPileUP(binX,binY) /= double(Sigma_B.at(binX)*Sigma_B.at(binY));
   }     
  }

  TH2D* Matrix_S_lowPileUp = new TH2D(correlationMatrixS_lowPileUP);
  Matrix_S_lowPileUp->SetName("Matrix_S_lowPileUp");
  TH2D* Matrix_B_lowPileUp = new TH2D(correlationMatrixB_lowPileUP);
  Matrix_B_lowPileUp->SetName("Matrix_B_lowPileUp");

  Matrix_S_lowPileUp->Scale(100);
  Matrix_B_lowPileUp->Scale(100);

  for( int iBinX = 0 ; iBinX < Matrix_S_lowPileUp->GetNbinsX() ; iBinX++){
   for( int iBinY = 0 ; iBinY < Matrix_S_lowPileUp->GetNbinsX() ; iBinY++){
     Matrix_S_lowPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_S_lowPileUp->GetBinContent(iBinX+1,iBinY+1)));
     Matrix_B_lowPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_B_lowPileUp->GetBinContent(iBinX+1,iBinY+1)));
     if(iBinX+1 == iBinY+1 ) {
      Matrix_S_lowPileUp->SetBinContent(iBinX+1,iBinY+1,1);
      Matrix_B_lowPileUp->SetBinContent(iBinX+1,iBinY+1,1);
     }
   }
  }


  TCanvas* cCorrelationSignal = new TCanvas("",Form("Correlation Matrix Signal"),180,52,550,550);

  float newMargin1 = 0.13;
  float newMargin2 = 0.15;
  float newMargin3 = 0.20;

  cCorrelationSignal->SetGrid();
  cCorrelationSignal->SetTicks();
  cCorrelationSignal->SetLeftMargin(newMargin3);
  cCorrelationSignal->SetBottomMargin(newMargin2);
  cCorrelationSignal->SetRightMargin(newMargin1);
  cCorrelationSignal->cd();
  gStyle->SetPaintTextFormat("3g");

  Matrix_S_lowPileUp->SetMarkerSize(1.5);
  Matrix_S_lowPileUp->SetMarkerColor(0);
  Matrix_S_lowPileUp->GetXaxis()->SetLabelSize(0.035);
  Matrix_S_lowPileUp->GetYaxis()->SetLabelSize(0.035);
  Matrix_S_lowPileUp->SetLabelOffset(0.011);// label offset on x axis                                                                                                                   

  Matrix_S_lowPileUp->SetMaximum(100);         
  Matrix_S_lowPileUp->SetMinimum(-100);         

  Matrix_S_lowPileUp->Draw("colz");
  for( int binX = 0 ; binX < Matrix_S_lowPileUp->GetNbinsX(); binX++){
    Matrix_S_lowPileUp->GetXaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
    Matrix_S_lowPileUp->GetYaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
  }
  Matrix_S_lowPileUp->Draw("textsame");
  Matrix_S_lowPileUp->GetXaxis()->LabelsOption("v");
  Matrix_S_lowPileUp->GetYaxis()->LabelsOption("h");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(21); // align right                                                                                                                                                  
  latex.SetTextSize(0.033);
  latex.DrawLatex(0.547,0.92,Form("CMS Preliminary Simulation, #sqrt{s} = 13 TeV"));

  cCorrelationSignal->Print((outputDirectory+"/CorrelationBDT_S.pdf").c_str(),"pdf");
  cCorrelationSignal->Print((outputDirectory+"/CorrelationBDT_S.png").c_str(),"png");
  cCorrelationSignal->Print((outputDirectory+"/CorrelationBDT_S.root").c_str(),"root");


  ////////////////////////
  TCanvas* cCorrelationBackground = new TCanvas("",Form("Correlation Matrix Background"),180,52,550,550);

  cCorrelationBackground->SetGrid();
  cCorrelationBackground->SetTicks();
  cCorrelationBackground->SetLeftMargin(newMargin3);
  cCorrelationBackground->SetBottomMargin(newMargin2);
  cCorrelationBackground->SetRightMargin(newMargin1);
  cCorrelationBackground->cd();

  Matrix_B_lowPileUp->SetMarkerSize(1.5);
  Matrix_B_lowPileUp->SetMarkerColor(0);
  Matrix_B_lowPileUp->GetXaxis()->SetLabelSize(0.035);
  Matrix_B_lowPileUp->GetYaxis()->SetLabelSize(0.035);
  Matrix_B_lowPileUp->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     

  Matrix_B_lowPileUp->SetMaximum(100);         
  Matrix_B_lowPileUp->SetMinimum(-100);         
           
  Matrix_B_lowPileUp->Draw("colz");
  for( int binX = 0 ; binX < Matrix_B_lowPileUp->GetNbinsX(); binX++){
    Matrix_B_lowPileUp->GetXaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
    Matrix_B_lowPileUp->GetYaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
  }
  Matrix_B_lowPileUp->Draw("textsame");
  Matrix_B_lowPileUp->GetXaxis()->LabelsOption("v");
  Matrix_B_lowPileUp->GetYaxis()->LabelsOption("h");
  
  latex.DrawLatex(0.547,0.92,Form("CMS Preliminary Simulation, #sqrt{s} = 13 TeV"));

  cCorrelationBackground->Print((outputDirectory+"/CorrelationBDT_B.pdf").c_str(),"pdf");
  cCorrelationBackground->Print((outputDirectory+"/CorrelationBDT_B.png").c_str(),"png");
  cCorrelationBackground->Print((outputDirectory+"/CorrelationBDT_B.root").c_str(),"root");

  //////////////////////////////////////////////////////////////////////
  inputLowPileUpTrees.clear() ;
  reducedNameLowPileUp.clear() ;
  // take all the TH1 outputs for S and B for  each variable in input-> low Pile Up
  itLowPileUp = InputInformationParamHighPU.begin();
  for( ; itLowPileUp != InputInformationParamHighPU.end() ; ++itLowPileUp){
    inputFile = TFile::Open((*itLowPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL) continue;
    inputLowPileUpTrees.push_back((TTree*) inputFile->Get("TestTree"));
    reducedNameLowPileUp.push_back((*itLowPileUp).getParameter<std::string>("variableName"));
  }

  // Build Up the correlation matrix
  TMatrixDSym correlationMatrixS_highPileUP(reducedNameLowPileUp.size()) ;
  TMatrixDSym correlationMatrixB_highPileUP(reducedNameLowPileUp.size()) ;
  for( unsigned int irow = 0 ; irow < reducedNameLowPileUp.size() ; irow++){
   for( unsigned int icolum = 0 ; icolum < reducedNameLowPileUp.size() ; icolum++){
     correlationMatrixS_highPileUP(irow,icolum) = 0;
     correlationMatrixB_highPileUP(irow,icolum) = 0;
   }
  }

  MeanValue_S.clear() ; MeanValue_S.resize(inputLowPileUpTrees.size());
  MeanValue_B.clear() ; MeanValue_B.resize(inputLowPileUpTrees.size());

  ientrySignal.clear() ; ientrySignal.resize(inputLowPileUpTrees.size());
  ientryBackground.clear() ; ientryBackground.resize(inputLowPileUpTrees.size());

  classID          = 0;
  BDTG_NoPruning = 0;

  // compute mean values
  for(unsigned int iTree = 0 ; iTree<inputLowPileUpTrees.size(); iTree++){   
    inputLowPileUpTrees.at(iTree)->SetBranchAddress("classID",&classID);
    inputLowPileUpTrees.at(iTree)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning);
    for(int iEntries = 0; iEntries < inputLowPileUpTrees.at(iTree)->GetEntries(); iEntries++){
      inputLowPileUpTrees.at(iTree)->GetEntry(iEntries);
      if(classID == 0){ MeanValue_S.at(iTree) += BDTG_NoPruning; ientrySignal.at(iTree)++;}
      else{ MeanValue_B.at(iTree) += BDTG_NoPruning; ientryBackground.at(iTree)++;}
   }
  }

  for( unsigned int iVec = 0 ; iVec < MeanValue_S.size(); iVec++){
    MeanValue_S.at(iVec) = MeanValue_S.at(iVec)/ientrySignal.at(iVec);
    MeanValue_B.at(iVec) = MeanValue_B.at(iVec)/ientryBackground.at(iVec);
  }

  // compute distances
  classID_X          = 0;
  BDTG_NoPruning_X = 0;
  classID_Y          = 0;
  BDTG_NoPruning_Y = 0;

  for(unsigned int iTreeX = 0 ; iTreeX<inputLowPileUpTrees.size(); iTreeX++){   
    for(unsigned int iTreeY = 0 ; iTreeY<inputLowPileUpTrees.size(); iTreeY++){    
     int iEntriesX = 0; int iEntriesY = 0;
     for( ; iEntriesX < inputLowPileUpTrees.at(iTreeX)->GetEntries() and iEntriesY < inputLowPileUpTrees.at(iTreeY)->GetEntries(); iEntriesX++, iEntriesY++){
       inputLowPileUpTrees.at(iTreeX)->SetBranchAddress("classID",&classID_X);
       inputLowPileUpTrees.at(iTreeX)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning_X);
       inputLowPileUpTrees.at(iTreeX)->GetEntry(iEntriesX);
       inputLowPileUpTrees.at(iTreeY)->SetBranchAddress("classID",&classID_Y);
       inputLowPileUpTrees.at(iTreeY)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning_Y);
       inputLowPileUpTrees.at(iTreeY)->GetEntry(iEntriesY);
       if(classID_X == 0 and classID_Y == 0) correlationMatrixS_highPileUP(iTreeX,iTreeY) += double((BDTG_NoPruning_X-MeanValue_S.at(iTreeX))*(BDTG_NoPruning_Y-MeanValue_S.at(iTreeY)));
       else if(classID_X == 1 and classID_Y == 1) correlationMatrixB_highPileUP(iTreeX,iTreeY) += double((BDTG_NoPruning_X-MeanValue_B.at(iTreeX))*(BDTG_NoPruning_Y-MeanValue_B.at(iTreeY)));
     }
    }
  }
  
  for( unsigned int binX = 0; binX < reducedNameLowPileUp.size(); binX ++){
   for( unsigned int binY = 0; binY < reducedNameLowPileUp.size(); binY ++){
     correlationMatrixS_highPileUP(binX,binY) /= double(ientrySignal.at(binX));
     correlationMatrixB_highPileUP(binX,binY) /= double(ientryBackground.at(binX));
   }
  }
 
  Sigma_S.clear() ; Sigma_S.resize(inputLowPileUpTrees.size());
  Sigma_B.clear() ; Sigma_B.resize(inputLowPileUpTrees.size());
  for( unsigned int binX = 0; binX < reducedNameLowPileUp.size(); binX ++) {
    Sigma_S.at(binX) = sqrt(correlationMatrixS_highPileUP(binX,binX));
    Sigma_B.at(binX) = sqrt(correlationMatrixB_highPileUP(binX,binX));
  }
  
  for( unsigned int binX = 0; binX < reducedNameLowPileUp.size(); binX ++){
   for( unsigned int binY = 0; binY < reducedNameLowPileUp.size(); binY ++){
     correlationMatrixS_highPileUP(binX,binY) /= double(Sigma_S.at(binX)*Sigma_S.at(binY));
     correlationMatrixB_highPileUP(binX,binY) /= double(Sigma_B.at(binX)*Sigma_B.at(binY));
   }     
  }

  TH2D* Matrix_S_highPileUp = new TH2D(correlationMatrixS_highPileUP);
  Matrix_S_highPileUp->SetName("Matrix_S_lowPileUp");
  TH2D* Matrix_B_highPileUp = new TH2D(correlationMatrixB_highPileUP);
  Matrix_B_highPileUp->SetName("Matrix_B_lowPileUp");

  Matrix_S_highPileUp->Scale(100);
  Matrix_B_highPileUp->Scale(100);

  for( int iBinX = 0 ; iBinX < Matrix_S_highPileUp->GetNbinsX() ; iBinX++){
   for( int iBinY = 0 ; iBinY < Matrix_S_highPileUp->GetNbinsX() ; iBinY++){
     Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_S_highPileUp->GetBinContent(iBinX+1,iBinY+1)));
     Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_B_highPileUp->GetBinContent(iBinX+1,iBinY+1)));
     if(iBinX+1 == iBinY+1 ) {
      Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,1);
      Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,1);
     }
   }
  }


  cCorrelationSignal = new TCanvas("",Form("Correlation Matrix Signal highPU"),180,52,500,550);
  cCorrelationSignal->SetGrid();
  cCorrelationSignal->SetTicks();
  cCorrelationSignal->SetLeftMargin(newMargin3);
  cCorrelationSignal->SetBottomMargin(newMargin2);
  cCorrelationSignal->SetRightMargin(newMargin1);
  cCorrelationSignal->cd();

  Matrix_S_highPileUp->SetMarkerSize(1.5);
  Matrix_S_highPileUp->SetMarkerColor(0);
  Matrix_S_highPileUp->GetXaxis()->SetLabelSize(0.035);
  Matrix_S_highPileUp->GetYaxis()->SetLabelSize(0.035);
  Matrix_S_highPileUp->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     
  Matrix_S_highPileUp->SetMaximum(100);         
  Matrix_S_highPileUp->SetMinimum(-100);         
           
  Matrix_S_highPileUp->Draw("colz");
  for( int binX = 0 ; binX < Matrix_S_highPileUp->GetNbinsX(); binX++){
    Matrix_S_highPileUp->GetXaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
    Matrix_S_highPileUp->GetYaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
  }
  Matrix_S_highPileUp->GetXaxis()->LabelsOption("v");
  Matrix_S_highPileUp->GetYaxis()->LabelsOption("h");
  Matrix_S_highPileUp->Draw("textsame");

  latex.DrawLatex(0.547,0.92,Form("CMS Preliminary Simulation, #sqrt{s} = 13 TeV"));
  cCorrelationSignal->Print((outputDirectory+"/CorrelationBDT_S_highPU.pdf").c_str(),"pdf");
  cCorrelationSignal->Print((outputDirectory+"/CorrelationBDT_S_highPU.png").c_str(),"png");
  cCorrelationSignal->Print((outputDirectory+"/CorrelationBDT_S_highPU.root").c_str(),"root");


  ////////////////////////
  cCorrelationBackground = new TCanvas("",Form("Correlation Matrix Background"),180,52,550,550);

  cCorrelationBackground->SetGrid();
  cCorrelationBackground->SetTicks();
  cCorrelationBackground->SetLeftMargin(newMargin3);
  cCorrelationBackground->SetBottomMargin(newMargin2);
  cCorrelationBackground->SetRightMargin(newMargin1);
  cCorrelationBackground->cd();

  Matrix_B_highPileUp->SetMarkerSize(1.5);
  Matrix_B_highPileUp->SetMarkerColor(0);
  Matrix_B_highPileUp->GetXaxis()->SetLabelSize(0.035);
  Matrix_B_highPileUp->GetYaxis()->SetLabelSize(0.035);
  Matrix_B_highPileUp->SetLabelOffset(0.011);// label offset on x axis                                                                                                                    
  Matrix_B_highPileUp->SetMaximum(100);         
  Matrix_B_highPileUp->SetMinimum(-100);         
           
            
  Matrix_B_highPileUp->Draw("colz");
  for( int binX = 0 ; binX < Matrix_B_highPileUp->GetNbinsX(); binX++){
    Matrix_B_highPileUp->GetXaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
    Matrix_B_highPileUp->GetYaxis()->SetBinLabel(binX+1,reducedNameLowPileUp.at(binX).c_str());
  }
  Matrix_B_highPileUp->GetXaxis()->LabelsOption("v");
  Matrix_B_highPileUp->GetYaxis()->LabelsOption("h");
  Matrix_B_highPileUp->Draw("textsame");

  latex.DrawLatex(0.547,0.92,Form("CMS Preliminary Simulation, #sqrt{s} = 13 TeV"));
  cCorrelationBackground->Print((outputDirectory+"/CorrelationBDT_B_highPU.pdf").c_str(),"pdf");
  cCorrelationBackground->Print((outputDirectory+"/CorrelationBDT_B_highPU.png").c_str(),"png");
  cCorrelationBackground->Print((outputDirectory+"/CorrelationBDT_B_highPU.root").c_str(),"root");

  return 0 ;

}
