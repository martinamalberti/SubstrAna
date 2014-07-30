#include "../src/JetTreeAnalyzer.cc"
#include "TChain.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>


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
 
//-------------------------------------------------------                                                                                                              
// MAIN                                                                                                                                                                       
//-------------------------------------------------------                                                                                                                                  
  
int main( int argc, char **argv ) {

  gROOT->ProcessLine("#include <vector>");

  if (argc<8){
    cout << "Usage: TestJetTreeAnalysis <ntuple name> <output name> <min pt> <max pt> <min eta> <max eta> <analyze also CMSSW jets>"<<endl;
    exit(0);
  }


  //TFile *inputFile = TFile::Open(argv[1]);
  //if (inputFile==0){
  //  std::cout<<"Error: cannot open " << inputFile->GetName() << std::endl;
  //  exit(0);
  //}
  
  //std::string inputname = argv[1];
  std::string inputFilesList = argv[1];
  std::string outname = argv[2];

  int maxEntries = -1;
  float minpt = atof(argv[3]);
  float maxpt  = atof(argv[4]);
  float minAbsEta = atof(argv[5]);
  float maxAbsEta = atof(argv[6]);
  bool doCMSSWJets = atoi(argv[7]);

  TChain* lTree = new TChain("Events");
  FillChain(*lTree, inputFilesList);


  // -- gen
  //TTree *tree_gen   = (TTree *)inputFile->Get("gen");
  TChain* tree_gen = new TChain("gen");
  FillChain(*tree_gen, inputFilesList);
  JetTreeAnalyzer *genAnalyzer = new JetTreeAnalyzer(tree_gen, tree_gen, "gen");
  genAnalyzer->bookHistograms("_gen");
  genAnalyzer->fillHistograms(maxEntries,minpt,maxpt,minAbsEta,maxAbsEta);
  //delete tree_gen;

  // -- pf
  //TTree *tree_pf    = (TTree *)inputFile->Get("pf");
  TChain* tree_pf = new TChain("pf");
  FillChain(*tree_pf, inputFilesList);
  JetTreeAnalyzer *pfAnalyzer = new JetTreeAnalyzer(tree_pf, tree_gen, "");
  pfAnalyzer->bookHistograms("_pf");
  pfAnalyzer->fillHistograms(maxEntries,minpt,maxpt,minAbsEta,maxAbsEta);
  delete tree_pf;

  // -- pfchs
  //  TTree *tree_pfchs = (TTree *)inputFile->Get("chs");
  TChain* tree_pfchs = new TChain("chs");
  FillChain(*tree_pfchs, inputFilesList);
  JetTreeAnalyzer *pfchsAnalyzer = new JetTreeAnalyzer(tree_pfchs, tree_gen, "");
  pfchsAnalyzer->bookHistograms("_pfchs");
  pfchsAnalyzer->fillHistograms(maxEntries,minpt,maxpt,minAbsEta,maxAbsEta);

  // -- puppi
  //TTree *tree_puppi = (TTree *)inputFile->Get("puppi");
  TChain* tree_puppi = new TChain("puppi");
  FillChain(*tree_puppi, inputFilesList);
  JetTreeAnalyzer *puppiAnalyzer = new JetTreeAnalyzer(tree_puppi, tree_gen, "");
  puppiAnalyzer->bookHistograms("_puppi");
  puppiAnalyzer->fillHistograms(maxEntries,minpt,maxpt,minAbsEta,maxAbsEta);
  delete tree_puppi;

  // -- softkiller
  //TTree *tree_softkiller = (TTree *)inputFile->Get("softkiller");
  TChain* tree_softkiller = new TChain("softkiller");
  FillChain(*tree_softkiller, inputFilesList);
  JetTreeAnalyzer *softkillerAnalyzer = new JetTreeAnalyzer(tree_softkiller, tree_gen, "");
  softkillerAnalyzer->bookHistograms("_softkiller");
  softkillerAnalyzer->fillHistograms(maxEntries,minpt,maxpt,minAbsEta,maxAbsEta);
  delete tree_softkiller;

  // -- pf cmssw
  //TTree *tree_pfcmssw;
  TChain *tree_pfcmssw;
  JetTreeAnalyzer *pfcmsswAnalyzer = 0;
  if (doCMSSWJets){
    //tree_pfcmssw = (TTree *)inputFile->Get("cmsswpf");
    tree_pfcmssw = new TChain("pfcmssw");
    FillChain(*tree_pfcmssw, inputFilesList);
    pfcmsswAnalyzer = new JetTreeAnalyzer(tree_pfcmssw, tree_gen, "");
    pfcmsswAnalyzer->bookHistograms("_pfcmssw");
    pfcmsswAnalyzer->fillHistograms(maxEntries, minpt,maxpt,minAbsEta,maxAbsEta);
    delete tree_pfcmssw;
  }

  // save results in file
  TFile *outfile = new TFile(outname.c_str(),"RECREATE");
  genAnalyzer->saveHistograms(outfile,"gen");
  pfAnalyzer->saveHistograms(outfile,"pf");
  pfchsAnalyzer->saveHistograms(outfile,"pfchs");
  puppiAnalyzer->saveHistograms(outfile,"puppi");
  softkillerAnalyzer->saveHistograms(outfile,"softkiller");
  if (doCMSSWJets) pfcmsswAnalyzer->saveHistograms(outfile,"pfcmssw");
  

}
