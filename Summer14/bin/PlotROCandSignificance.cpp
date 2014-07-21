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

#include "../include/PlottingMVAResults.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"


// method to get the string cut for the related preselection used in the trainint --> implementation at bottom of the code
TString GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType, const double & pTJetMin_, const double & pTJetMax_, const double & npu_Min, const double & npu_Max, const std::string & TreeName);

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
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetErrorX(0.5);


  // parse config file parameter                                                                                                                                                          
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");


  std::vector<edm::ParameterSet> InputInformationParam ;
  if(Options.existsAs<std::vector<edm::ParameterSet>>("InputInformationParam"))
    InputInformationParam = Options.getParameter<std::vector<edm::ParameterSet>>("InputInformationParam");
  else{ std::cout<<" Exit from code, no input set found for MVA training output"<<std::endl; return -1; }

  std::string TreeName ;
  if(Options.existsAs<std::string>("TreeName"))
    TreeName = Options.getParameter<std::string>("TreeName");
  else{ std::cout<<" Exit from code, no TreeName found "<<std::endl; return -1; }

  std::string Label;
  if(Options.existsAs<std::string>("Label"))
    Label = Options.getParameter<std::string>("Label");
  else{ std::cout<<" Label set to Test "<<std::endl; Label = "Test"; }

  std::string LeptonType;
  if(Options.existsAs<std::string>("LeptonType"))
    LeptonType  = Options.getParameter<std::string>("LeptonType");
  else{ std::cout<<" Lepton type set by default to muon "<<std::endl; LeptonType = "Muon"; }

  std::cout<<std::endl;
  std::cout<<" Input TreeName       = "<<TreeName<<std::endl;
  std::cout<<" Input LeptonType     = "<<LeptonType<<std::endl;
  std::cout<<std::endl;


  std::string SignalggHName ;
  if(Options.existsAs<std::string>("SignalggHName"))
    SignalggHName  = Options.getParameter<std::string>("SignalggHName");
  else std::cout<<" SignalggHName --> not found --> empty string "<<std::endl;

  std::cout<<" Option Signal ggH Name = "<<SignalggHName<<std::endl;
  std::cout<<std::endl;

  std::string SignalqqHName ;
  if(Options.existsAs<std::string>("SignalqqHName"))
    SignalqqHName  = Options.getParameter<std::string>("SignalqqHName");
  else std::cout<<" SignalqqHName --> not found --> empty string "<<std::endl;

  std::cout<<" Option Signal qqH Name = "<<SignalqqHName<<std::endl;
  std::cout<<std::endl;

  std::string EventWeight ;
  if(Options.existsAs<std::string>("EventWeight"))
    EventWeight  = Options.getParameter<std::string>("EventWeight");
  else std::cout<<" EventWeight --> not found --> empty string "<<std::endl;

  std::cout<<" Option Event Weight = "<<EventWeight<<std::endl;
  std::cout<<std::endl;

  int useTypeOfSignal = 0;
  if(Options.existsAs<int>("useTypeOfSignal"))
    useTypeOfSignal  = Options.getParameter<int>("useTypeOfSignal");
  else std::cout<<" type of signal to be considered not found --> set to 0 "<<std::endl;

  std::cout<<" Option  useTypeOfSignal = "<<useTypeOfSignal<<std::endl;
  std::cout<<std::endl;


  double Lumi = 0;
  if(Options.existsAs<double>("Lumi"))
    Lumi  = Options.getParameter<double>("Lumi");
  else std::cout<<" Lumi --> not found --> set to 0 "<<std::endl;

  std::cout<<" Option Lumi = "<<Lumi<<std::endl;
  std::cout<<std::endl;


  std::string PreselectionCutType ;
  if(Options.existsAs<std::string>("PreselectionCutType"))
    PreselectionCutType  = Options.getParameter<std::string>("PreselectionCutType");
  else{ std::cout<<" PreselectionCutType --> not found --> empty string "<<std::endl; PreselectionCutType = "none"; }


  std::cout<<" Option Preselection Cut = "<<PreselectionCutType<<std::endl;
  std::cout<<std::endl;


  std::string outputPlotDirectory ;
  if(Options.existsAs<std::string>("outputPlotDirectory"))
    outputPlotDirectory = Options.getParameter<std::string>("outputPlotDirectory");
  else{std::cout<<" Default output directory "<<std::endl;
    outputPlotDirectory = "output/TMVATrainingResult";
  }

  std::cout<<"                      "<<std::endl;
  std::cout<<" Outout Directory     "<<outputPlotDirectory<<std::endl;
  std::cout<<"                      "<<std::endl;

  std::string command;
  command = "if [ ! -e "+outputPlotDirectory+" ] ; then mkdir "+outputPlotDirectory+" ; fi";
  std::cout<<" command = "<<command<<std::endl;
  std::cout<<"           "<<std::endl;

  system(command.c_str());

  command = "if [ ! -f "+outputPlotDirectory+" ] ; then rm "+outputPlotDirectory+"/* ; fi";
  std::cout<<" command = "<<command<<std::endl;
  std::cout<<"           "<<std::endl;

  system(command.c_str());

    
  // Declare the object for the manipolation of the TMVA ROOT file
  TMVAGlob* TMVATraining = new TMVAGlob();
  std::vector<std::string> fileName ;
  std::vector<std::string> inputVariableName ;
  std::vector<std::string> inputVariableReducedName ;
  std::vector<std::pair<double,double>> ptBin  ;
  std::vector<std::pair<double,double>> puBin  ;
  
  std::vector<edm::ParameterSet>::const_iterator itFile = InputInformationParam.begin();
  for( ; itFile != InputInformationParam.end() ; ++itFile){
   // O the set of inputFiles and get back to the main code  
   fileName.push_back((*itFile).getParameter<std::string>("fileName"));
   for( unsigned int iName = 0 ; iName < (*itFile).getParameter<std::vector<std::string> >("inputVariableOrMethodName").size(); iName++){
     inputVariableName.push_back((*itFile).getParameter<std::vector<std::string> >("inputVariableOrMethodName").at(iName));
     std::cout<<" method name "<<(*itFile).getParameter<std::vector<std::string> >("inputVariableOrMethodName").at(iName)<<std::endl;
   }
   inputVariableReducedName.push_back((*itFile).getParameter<std::string>("ReducedName"));
 
   std::pair<double,double> pairtemp ; 
   pairtemp.first = (*itFile).getParameter<std::vector<double> >("JetPtBinOfTraining").at(0);
   pairtemp.second =(*itFile).getParameter<std::vector<double> >("JetPtBinOfTraining").at(1);
   ptBin.push_back(pairtemp);
   pairtemp.first = (*itFile).getParameter<std::vector<double> >("PileUpBinOfTraining").at(0);
   pairtemp.second =(*itFile).getParameter<std::vector<double> >("PileUpBinOfTraining").at(1);
   puBin.push_back(pairtemp);
  }

  TMVATraining->openFileInput(fileName);
  TMVATraining->SetMethodName(inputVariableName);     
  std::vector<TFile*> inputFile = TMVATraining->GetInputFile();
  TMVATraining->plotEfficiency(inputFile,gDirectory,ptBin.at(0).first,ptBin.at(0).second,puBin.at(0).first,puBin.at(0).second); // call the plot efficiency function 
  TMVATraining->PrintImageROC(gDirectory,outputPlotDirectory);
  for(unsigned int iFile = 0 ; iFile < inputFile.size(); iFile++){
    TMVATraining->plotCorrelationMatrix(inputFile.at(iFile),inputVariableReducedName.at(iFile),outputPlotDirectory);
    TMVATraining->plotMVAs(inputFile.at(iFile),inputVariableReducedName.at(iFile),0,outputPlotDirectory);
    TMVATraining->plotMVAs(inputFile.at(iFile),inputVariableReducedName.at(iFile),3,outputPlotDirectory);
    TMVATraining->plotSignificance(inputFile.at(iFile),inputVariableReducedName.at(iFile),1,1,1,true,true,outputPlotDirectory);

  }
  return 0 ;
}


TString GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType, const double & pTJetMin_, const double & pTJetMax_, const double & npu_Min, const double & npu_Max, const std::string & TreeName){

  //--------------------------                                                                                                                                                            
  // Basic preselection CSA14                                                                                                                                                             
  //--------------------------                                                                                                                                                            
 
  if( preselectionCutType == "basicJetsCutCSA14" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "electron" || LeptonType == "El" || LeptonType== "el" || LeptonType == "Electron" || LeptonType == "Jets" || LeptonType == "jets") and  TreeName !="gen")
    return Form("ptraw[0]>200 && fabs(eta[0])<2.5 && imatch[0] >= 0 && (ptraw[0] > %f  && ptraw[0] < %f ) && (npu > %f && npu < %f)",pTJetMin_,pTJetMax_,npu_Min,npu_Max);

  else if ( preselectionCutType == "basicJetsCutCSA14" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "electron" || LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "Jets" || LeptonType == "jets") and  TreeName =="gen")
    return Form("ptraw[0] > 200 && abs(eta[0])<2.5 && (ptraw[0] > %f  && ptraw[0] < %f )  && (npu > %f && npu < %f)",pTJetMin_,pTJetMax_,npu_Min,npu_Max);

  else return Form("v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


}
