
#include "../include/TrainingMVAClass.h"

// constructor giving files 
TrainingMVAClass::TrainingMVAClass(const std::vector<TFile*> & signalFileList, const std::vector<TFile*> & backgroundFileList, const std::string & TreeName,
				   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label, const std::string &  transformation){

 
  (*this).SetSignalTree (signalFileList,TreeName) ;

  (*this).SetBackgroundTree (backgroundFileList,TreeName) ;

  (*this).SetLabel (Label);

  (*this).SetOutputFile (outputFilePath , outputFileName) ;

  factory_ = new TMVA::Factory (TreeName_+"_"+Label_,outputFile_, Form("!V:!Silent:%sColor:DrawProgressBar:AnalysisType=Classification%s",gROOT->IsBatch()?"!":"",transformation.c_str()));

}

// constructor giving tree
TrainingMVAClass::TrainingMVAClass(const std::vector<TTree*> & signalTreeList, const std::vector<TTree*> & backgroundTreeList,  const std::string & TreeName,
                                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label,const std::string &  transformation){
   

  (*this).SetTreeName (TreeName) ;

  (*this).SetSignalTree (signalTreeList) ;

  (*this).SetBackgroundTree (backgroundTreeList) ;

  (*this).SetLabel (Label);

  (*this).SetOutputFile (outputFilePath , outputFileName) ;

  factory_ = new TMVA::Factory (TreeName_+"_"+Label_,outputFile_, Form("!V:!Silent:%sColor:DrawProgressBar:AnalysisType=Classification%s",gROOT->IsBatch()?"!":"",transformation.c_str()));

}

// Deconstructor
TrainingMVAClass::~TrainingMVAClass(){

  for(size_t iTree = 0; iTree < signalTreeList_.size() ; iTree++) { 
   if(signalTreeList_.at(iTree)!=0)  delete signalTreeList_.at(iTree) ; 
  }

  for(size_t iTree = 0; iTree < backgroundTreeList_.size() ; iTree++) {
   if(backgroundTreeList_.at(iTree)!=0)  delete backgroundTreeList_.at(iTree) ;
  }

  if(preselectionCut_!=0) preselectionCut_->Delete() ;

  if(outputFile_!=0) outputFile_->Close() ;
 
  if(factory_!=0) factory_->Delete() ;

}

// AddTrainingVariables in the MVA
void TrainingMVAClass::AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables){

  (*this).SetTrainingVariables(mapTrainingVariables);
  (*this).SetSpectatorVariables(mapSpectatorVariables);

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ){
    std::cout<<" train " <<mapTrainingVariables_.at(iVar)<<std::endl;
    factory_->AddVariable(mapTrainingVariables_.at(iVar)+" := "+mapTrainingVariables_.at(iVar),'F');
  }

  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ){
    std::cout<<" spectator " <<mapSpectatorVariables_.at(iVar)<<std::endl;
    factory_->AddSpectator(mapSpectatorVariables_.at(iVar),'F');
  }    
}

// AddTrainingVariables in the MVA
void TrainingMVAClass::AddTrainingVariables ( const std::string & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables){

  
  mapTrainingVariables_.push_back(mapTrainingVariables);
  (*this).SetSpectatorVariables(mapSpectatorVariables);

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ){
    std::cout<<" train " <<mapTrainingVariables_.at(iVar)<<std::endl;
    factory_->AddVariable(mapTrainingVariables_.at(iVar)+" := "+mapTrainingVariables_.at(iVar),'F');
  }

  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ){
    std::cout<<" spectator " <<mapSpectatorVariables_.at(iVar)<<std::endl;
    factory_->AddSpectator(mapSpectatorVariables_.at(iVar),'F');
  }    
}


// Book MVA Training Variables 
void TrainingMVAClass::BookMVATrees (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight){
  
  (*this).SetGlobalSampleWeight(signalGlobalWeight,backgroundGlobalWeight);
 
  if(signalGlobalWeight.size() == signalTreeList_.size()){

    for(size_t iTree = 0; iTree<signalTreeList_.size(); iTree ++) 
      factory_->AddSignalTree (signalTreeList_.at(iTree),signalGlobalWeight.at(iTree)) ;
  }
  else{
        
    for(size_t iTree = 0; iTree<signalTreeList_.size(); iTree ++) 
      factory_->AddSignalTree (signalTreeList_.at(iTree),1.0) ;
  }

  if(backgroundGlobalWeight.size() == backgroundTreeList_.size()){

    for(size_t iTree = 0; iTree<backgroundTreeList_.size(); iTree ++) 
      factory_->AddBackgroundTree (backgroundTreeList_.at(iTree),backgroundGlobalWeight.at(iTree)) ;
  }
  else{
        
    for(size_t iTree = 0; iTree<backgroundTreeList_.size(); iTree ++) 
      factory_->AddBackgroundTree (backgroundTreeList_.at(iTree),1.0) ;
  }
  
}


// Prepare the training addinf the preselection cut and the train/testing composition
void TrainingMVAClass::AddPrepareTraining ( const std::string & LeptonType, const std::string & preselectionCutType, const std::string & weightString,
                                            std::vector<double> * JetPtBinOfTraining, const int & pTBin, std::vector<double> * PileUpBinOfTraining, const int & puBin,
                                            const int & nTraining, const int & nTesting, const std::string & splitMode, const std::string & NormMode){



  if(JetPtBinOfTraining!=NULL || !JetPtBinOfTraining ){
     pTJetMin_ = JetPtBinOfTraining->at(pTBin) ;
     pTJetMax_ = JetPtBinOfTraining->at(pTBin+1) ;
  }
  else{
       pTJetMin_ = 0;
       pTJetMax_ = 2000 ;
  }

  if(PileUpBinOfTraining!=NULL || !PileUpBinOfTraining ){
     npuMin_ = PileUpBinOfTraining->at(puBin) ;
     npuMax_ = PileUpBinOfTraining->at(puBin+1) ;
  }
  else{
       npuMin_ = 0;
       npuMax_ = 100 ;
  }


  preselectionCut_ = new TCut (GetPreselectionCut(LeptonType,preselectionCutType).Data()) ;

  std::cout<<" Get Preselection Cut "<<GetPreselectionCut(LeptonType,preselectionCutType).Data()<<std::endl;

  TString Option = Form("nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:SplitMode=%s:NormMode=%s:!V",
                         nTraining,nTraining,nTesting,nTesting,splitMode.c_str(),NormMode.c_str());

  (*this).SetEventWeight (weightString);
  factory_->PrepareTrainingAndTestTree( *(preselectionCut_),*(preselectionCut_), Option.Data() );

}

// Train Rectangular Cuts
void TrainingMVAClass::BookandTrainRectangularCuts (const std::string & FitMethod, std::string variable){

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 
  std::replace(variable.begin(),variable.end(),'/', '_');
  std::replace(variable.begin(),variable.end(),'[', '_');
  std::replace(variable.begin(),variable.end(),']', '_');
  std::replace(variable.begin(),variable.end(),'(', '_');
  std::replace(variable.begin(),variable.end(),')', '_');
  std::replace(variable.begin(),variable.end(),':', '_');

  // Set Name of the Weight file for TMVA evaluating procedure
  outputFileWeightName_["Cuts"+FitMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_Cuts"+FitMethod+"_"+Label_+"_"+variable;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["Cuts"+FitMethod+"_"+Label_];

  if(FitMethod!=""){ 
     TString Option = Form("!H:!V:FitMethod=%s:EffSel:SampleSize=200000:VarProp=FSmart%s", FitMethod.c_str(),transformations_.c_str());
     TString Name = Form("Cuts%s",FitMethod.c_str());
     if(TString(Name).Contains("CutsGA"))  
       Option = Option+":CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95";
     else if (TString(Name).Contains("CutsSA"))
       Option = Option+":MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale";
     factory_->BookMethod( TMVA::Types::kCuts, Name.Data(),Option.Data());
  }

  else{
        TString Option = Form("!H:!V:FitMethod=MC:EffSel%s",transformations_.c_str());
        factory_->BookMethod( TMVA::Types::kCuts, "CutsMC"+FitMethod,Option.Data());
        Option = Form("!H:!V:FitMethod=GA::CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95%s",transformations_.c_str());
        factory_->BookMethod( TMVA::Types::kCuts, "CutsGA"+FitMethod,Option.Data());
        Option = Form("!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale%s",transformations_.c_str());
        factory_->BookMethod( TMVA::Types::kCuts, "CutsSA"+FitMethod,Option.Data());

  }

  factory_->TrainAllMethods();
  factory_->TestAllMethods();
  factory_->EvaluateAllMethods();
  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout<< "==> TMVAClassification is done!" << std::endl;

}

// Train Likelihood
void TrainingMVAClass::BookandTrainLikelihood ( const std::string & LikelihoodType ){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure
  outputFileWeightName_[LikelihoodType+"_"+Label_] = outputFilePath_+"/TMVAWeight_"+LikelihoodType+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_[LikelihoodType+"_"+Label_];

  TString Option ;

  if( LikelihoodType == "LikelihoodKDE") { Option = Form("LikelihoodKDE");
    factory_->BookMethod(TMVA::Types::kLikelihood, Option.Data(),"!H:!V:IgnoreNegWeightsInTraining:!TransformOutput:"
                                                                 "PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:CreateMVAPdfs:KDEFineFactor=0.3:KDEborder=None%s"+transformations_);
  }
  else if( LikelihoodType == "PDERS") { Option = Form("%s",LikelihoodType.c_str());
      factory_->BookMethod(TMVA::Types::kPDERS, Option.Data(),
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:CreateMVAPdfs:DeltaFrac=4:GaussSigma=0.3:NormTree=T"+transformations_);
  }
  else if( LikelihoodType == "PDEFoam") { Option = Form("%s",LikelihoodType.c_str());
       factory_->BookMethod(TMVA::Types::kPDEFoam, Option.Data(),"!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:SigBgSeparate=F:TailCut=0.001"
                                                                 ":VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T"+transformations_);
  }
  else if( LikelihoodType == "PDEFoamBoost") { Option = Form("%s",LikelihoodType.c_str());
      factory_->BookMethod(TMVA::Types::kPDEFoam, Option.Data(),
                           "!H:!V:IgnoreNegWeightsInTraining:Boost_Num=30:CreateMVAPdfs:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4"
                           ":UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=300:nBin=20:Nmin=300:Kernel=None:Compress=T"+transformations_);
  }
  else{ Option = Form("%s",LikelihoodType.c_str());
        factory_->BookMethod( TMVA::Types::kLikelihood, Option.Data(),"!H:!V:!TransformOutput:CreateMVAPdfs:IgnoreNegWeightsInTraining:PDFInterpol=Spline2"
			                                              ":NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50"+transformations_);
  }


  factory_->TrainAllMethods();
  factory_->TestAllMethods();
  factory_->EvaluateAllMethods();
  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

// Train Fisher Discriminant
void TrainingMVAClass::BookandTrainFisherDiscriminant(){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure                                                              

  outputFileWeightName_["Fisher"+Label_] = outputFilePath_+"/TMVAWeight_Fisher_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["Fisher"+Label_];


  factory_->BookMethod( TMVA::Types::kFisher, "Fisher",
                        "!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10:Fisher"+transformations_ );

  factory_->TrainAllMethods();
  factory_->TestAllMethods();
  factory_->EvaluateAllMethods();
  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

// Train Linear Discriminant
void TrainingMVAClass::BookandTrainLinearDiscriminant(){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 
  // Set Name of the Weight file for TMVA evaluating procedure

  outputFileWeightName_["LD"+Label_] = outputFilePath_+"/TMVAWeight_LD_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["LD"+Label_];

  // Training Testing and Evaluating   
  outputFile_->cd();
  factory_->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=I,N,D,P:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10"+transformations_);

  factory_->TrainAllMethods();
  factory_->TestAllMethods();
  factory_->EvaluateAllMethods();
  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

// Train MLP
void TrainingMVAClass::BookandTrainMLP(const int & nCycles, const std::string & HiddenLayers, const std::string & NeuronType,
				       const std::string & TrainingMethod, const int & TestRate, const int & ConvergenceTests){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                                
  outputFileWeightName_["MLP_"+NeuronType+"_"+TrainingMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_MLP_"+NeuronType+"_"+TrainingMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["MLP_"+NeuronType+"_"+TrainingMethod+"_"+Label_];

  TString Option = Form ("!H:!V:NCycles=%d:CalculateErrors:HiddenLayers=%s:NeuronType=%s:CreateMVAPdfs:TrainingMethod=%s:TestRate=%d"
			 ":ConvergenceTests=%d:UseRegulator%s",nCycles,HiddenLayers.c_str(),NeuronType.c_str(),TrainingMethod.c_str(),TestRate,ConvergenceTests,transformations_.c_str());

  factory_->BookMethod( TMVA::Types::kMLP, "MLP_"+NeuronType+"_"+TrainingMethod, Option.Data());
  
  factory_->TrainAllMethods();
  factory_->TestAllMethods();
  factory_->EvaluateAllMethods();
  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

// Train Clemont Ferrand ANN
void TrainingMVAClass::BookandTrainCFMlpANN ( const int & nCycles, const std::string & HiddenLayers){

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                           
  outputFileWeightName_["CFMlpANN_"+Label_] = outputFilePath_+"/TMVAWeight_CFMlpANN_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["CFMlpANN_"+Label_];

  TString Option = Form ("!H:!V:NCycles=%d:HiddenLayers=%s:CreateMVAPdfs%s",nCycles,HiddenLayers.c_str(),transformations_.c_str());

  factory_->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN",Option.Data());

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}


// Train TMVA ANN
void TrainingMVAClass::BookandTrainTMlpANN  ( const int & nCycles, const std::string & HiddenLayers,  const std::string & TrainingMethod, const float & ValidationFraction){

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                            
  outputFileWeightName_["TMlpANN_"+TrainingMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_TMlpANN_"+TrainingMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["TMlpANN_"+TrainingMethod+"_"+Label_];

  TString Option = Form ("!H:!V:NCycles=%d:HiddenLayers=%s:LearningMethod=%s:ValidationFraction=%f:CreateMVAPdfs%s",
			 nCycles,HiddenLayers.c_str(),TrainingMethod.c_str(),ValidationFraction,transformations_.c_str());

  factory_->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN_"+TrainingMethod,Option.Data());

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

// Train BDT
void TrainingMVAClass::BookandTrainBDT ( const int & NTrees, const bool & optimizeMethods, const std::string & BoostType, const float & AdaBoostBeta,
				         const std::string & PruneMethod, const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure                                                          
  outputFileWeightName_["BDT_"+BoostType+"_"+PruneMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_BDT_"+BoostType+"_"+PruneMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["BDT_"+BoostType+"_"+PruneMethod+"_"+Label_];
  
  TString Option = Form ("!H:!V:CreateMVAPdfs:NTrees=%d:BoostType=%s:AdaBoostBeta=%f:PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s:Shrinkage=0.1:NNodesMax=100000:UseYesNoLeaf=F:nEventsMin=200:nCuts=200%s",NTrees,BoostType.c_str(),AdaBoostBeta,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str(),transformations_.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDT_"+BoostType+"_"+PruneMethod, Option.Data());

  if(optimizeMethods) factory_->OptimizeAllMethods();                                                                                                                                                            
  factory_->TrainAllMethods();
  factory_->TestAllMethods();
  factory_->EvaluateAllMethods();
  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

// Train BDTG
void TrainingMVAClass::BookandTrainBDTG ( const int & NTrees, const bool & optimizeMethods, const float & GradBaggingFraction, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                           
  outputFileWeightName_["BDTG_"+PruneMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_BDTG_"+PruneMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["BDTG_"+PruneMethod+"_"+Label_];

  TString Option = Form ("CreateMVAPdfs:NTrees=%d:BoostType=Grad:!UseBaggedGrad:GradBaggingFraction=%f:PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s:Shrinkage=0.1:NNodesMax=100000:UseYesNoLeaf=F:nEventsMin=200:nCuts=200%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str(),transformations_.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDTG_"+PruneMethod, Option.Data());
  
  if(optimizeMethods) factory_->OptimizeAllMethods();                                                                                                                                                           
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

// Train BDTF
void TrainingMVAClass::BookandTrainBDTF ( const int & NTrees, const bool & optimizeMethods, const std::string & BoostType, const float & AdaBoostBeta, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  int result = system(command.c_str());
  if(result) std::cout<<"Directory created "<<outputFilePath_<<std::endl; 

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                           
  outputFileWeightName_["BDTF_"+PruneMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_BDTF_"+PruneMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["BDTF_"+PruneMethod+"_"+Label_];

  TString Option = Form ("!H:!V:CreateMVAPdfs:UseFisherCuts:NTrees=%d:BoostType=%s:AdaBoostBeta=%f:PruneMethod=%s:"
                         "PruneStrength=%d:MaxDepth=%d:SeparationType=%s:Shrinkage=0.10:nCuts=20%s",NTrees,BoostType.c_str(),
                          AdaBoostBeta,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str(),transformations_.c_str());

  factory_->BookMethod( TMVA::Types::kBDT,"BDTF_"+PruneMethod+"_", Option.Data());

  if(optimizeMethods) factory_->OptimizeAllMethods();
                                                                                                                                                             
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
}


//Set Signal Tree Name
void TrainingMVAClass::SetSignalTree (const std::vector<TFile*> & signalFileList, const std::string & TreeName){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;

   for(size_t iFile = 0 ; iFile < signalFileList.size() ; iFile ++){
     if(signalFileList.at(iFile)!=0)  signalTreeList_.push_back((TTree*) signalFileList.at(iFile)->Get(TreeName_.c_str()));
   }

   return ;

}


void TrainingMVAClass::SetSignalTree (const std::vector<TTree*> & signalTreeList){
  
  if(signalTreeList.size()!=0) signalTreeList_ = signalTreeList ; 

  return ;

}

// Set Background Tree Name
void TrainingMVAClass::SetBackgroundTree (const std::vector<TFile*> & backgroundFileList, const std::string & TreeName){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;
     
  for(size_t iFile = 0 ; iFile < backgroundFileList.size() ; iFile ++){
     if(backgroundFileList.at(iFile)!=0) backgroundTreeList_.push_back((TTree*) backgroundFileList.at(iFile)->Get(TreeName_.c_str()));
  }

  return ;


}

void TrainingMVAClass::SetBackgroundTree (const std::vector<TTree*> & backgroundTreeList){

   if(backgroundTreeList.size()!=0) backgroundTreeList_ = backgroundTreeList ;

   return ;

}

// Set Training Variables 
void TrainingMVAClass::SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables){

   if(mapTrainingVariables.size()!=0) mapTrainingVariables_=mapTrainingVariables;

   return ;

}

// Set Spectator Variables
void TrainingMVAClass::SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables){

   if(mapSpectatorVariables.size()!=0) mapSpectatorVariables_=mapSpectatorVariables;

   return ;
}

// Set label
void TrainingMVAClass::SetLabel (const std::string & Label ){

  Label_ = Label ;

  return ;

}

// Set Tree Name
void TrainingMVAClass::SetTreeName (const std::string & TreeName ){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;

  return ;
}

// set output file
void TrainingMVAClass::SetOutputFile ( const std::string & outputFilePath , const std::string & outputFileName ){
  
  if( !outputFilePath.empty() && !outputFileName.empty()) { 
 
   outputFilePath_=outputFilePath; 
   outputFileName_=outputFileName;

   std::replace(outputFileName_.begin(),outputFileName_.end(),'/','_');
   std::replace(outputFileName_.begin(),outputFileName_.end(),'[','_');
   std::replace(outputFileName_.begin(),outputFileName_.end(),']','_');
   std::replace(outputFileName_.begin(),outputFileName_.end(),'(','_');
   std::replace(outputFileName_.begin(),outputFileName_.end(),')','_');
   std::replace(outputFileName_.begin(),outputFileName_.end(),':','_');

   outputFileNameComplete_ = outputFilePath_+"/"+outputFileName_+"_"+Label_+".root" ;

   outputFile_ = new TFile((outputFilePath_+"/"+outputFileName_+"_"+Label_+".root").c_str(),"RECREATE");
   
   outputFile_->cd();
 }

  return ;

}

// set global event weight
void TrainingMVAClass::SetGlobalSampleWeight (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight){

  signalGlobalWeight_ =  signalGlobalWeight ;
  backgroundGlobalWeight_ = backgroundGlobalWeight ;

  return ;

}


void TrainingMVAClass::SetEventWeight (const std::string & weightString){

  factory_->SetSignalWeightExpression(weightString);
  factory_->SetBackgroundWeightExpression(weightString);

  return ;

}

void TrainingMVAClass::SetTransformations (const std::string & transformations){

  transformations_ = transformations;  

  return ;

}

// Take Preselection Selection 

TString TrainingMVAClass::GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType){

  //-----------------------------------------------
  // Basic preselection used in the exotica analysis
  //-----------------------------------------------

  if(preselectionCutType == "basicPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);
	     
  else if(preselectionCutType == "basicPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSBPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 ) ) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 ) ) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 ) ) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSRPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form(" issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=60 && jet_mass_pr <= 100 ) && (jet_GeneralizedECF >0 && jet_GeneralizedECF <1) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=60 && jet_mass_pr <= 100 ) && nbjets_csvm_veto == 0 && (jet_GeneralizedECF >0 && jet_GeneralizedECF < 1 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 50 && l_pt > 50 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=60 && jet_mass_pr <= 100 ) && nbjets_csvm_veto == 0 &&(jet_GeneralizedECF > 0 && jet_GeneralizedECF < 1 ) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  //-----------------------------------------------
  // Basic preselection used in the exotica analysis using RawMass
  //-----------------------------------------------

  else if(preselectionCutType == "basicSRPreselectionCutEXO_RawMass" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form(" issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( TMath::Sqrt(ungroomed_jet_e*ungroomed_jet_e-ungroomed_jet_pt*ungroomed_jet_pt-TMath::SinH(ungroomed_jet_eta)*TMath::SinH(ungroomed_jet_eta)*ungroomed_jet_pt*ungroomed_jet_pt) >=75 && TMath::Sqrt(ungroomed_jet_e*ungroomed_jet_e-ungroomed_jet_pt*ungroomed_jet_pt-TMath::SinH(ungroomed_jet_eta)*TMath::SinH(ungroomed_jet_eta)*ungroomed_jet_pt*ungroomed_jet_pt) <= 115 ) && (jet_GeneralizedECF >0 && jet_GeneralizedECF <1) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f ) ",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutEXO_RawMass" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( TMath::Sqrt(ungroomed_jet_e*ungroomed_jet_e-ungroomed_jet_pt*ungroomed_jet_pt-TMath::SinH(ungroomed_jet_eta)*TMath::SinH(ungroomed_jet_eta)*ungroomed_jet_pt*ungroomed_jet_pt) >=75 && TMath::Sqrt(ungroomed_jet_e*ungroomed_jet_e-ungroomed_jet_pt*ungroomed_jet_pt-TMath::SinH(ungroomed_jet_eta)*TMath::SinH(ungroomed_jet_eta)*ungroomed_jet_pt*ungroomed_jet_pt) <= 115 ) && nbjets_csvm_veto == 0 && (jet_GeneralizedECF >0 && jet_GeneralizedECF < 1 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutEXO_RawMass" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 50 && l_pt > 50 && ungroomed_jet_pt > 200 && ( TMath::Sqrt(ungroomed_jet_e*ungroomed_jet_e-ungroomed_jet_pt*ungroomed_jet_pt-TMath::SinH(ungroomed_jet_eta)*TMath::SinH(ungroomed_jet_eta)*ungroomed_jet_pt*ungroomed_jet_pt) >=75 && TMath::Sqrt(ungroomed_jet_e*ungroomed_jet_e-ungroomed_jet_pt*ungroomed_jet_pt-TMath::SinH(ungroomed_jet_eta)*TMath::SinH(ungroomed_jet_eta)*ungroomed_jet_pt*ungroomed_jet_pt) <= 115 ) && nbjets_csvm_veto == 0 &&(jet_GeneralizedECF > 0 && jet_GeneralizedECF < 1 ) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);
    

  //-----------------------------------------------
  // Basic preselection used in the Higgs analysis
  //-----------------------------------------------

  if(preselectionCutType == "basicPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBin < 2 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);
	      
 else if(preselectionCutType == "basicPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

 else if(preselectionCutType == "basicPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && numberJetBin < 2 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSBPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && nbjets_csvm_veto == 0 && ( (jet_mass_pr > 40 && jet_mass_pr < 65) || (jet_mass_pr > 105 && jet_mass_pr < 130) ) && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 60 ) || ( jet_mass_pr >=100 && jet_mass_pr <= 130 ) ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 ) ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSRPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && nbjets_csvm_veto == 0 && (jet_mass_pr > 65 && jet_mass_pr < 105) && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f)",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=65 && jet_mass_pr <= 105 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=65 && jet_mass_pr <= 105 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSRSBPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && nbjets_csvm_veto == 0 && (jet_mass_pr > 40 && jet_mass_pr < 130) && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);
    
  else if(preselectionCutType == "basicSRSBPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  //-----------------------------------------------
  // Basic preselection used in the Higgs analysis + VBF cuts
  //-----------------------------------------------

  else if( preselectionCutType == "basicVBFPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if( preselectionCutType == "basicVBFSBHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ((jet_mass_pr > 40 && jet_mass_pr <65) || (jet_mass_pr > 105 && jet_mass_pr < 130)) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ((jet_mass_pr > 40 && jet_mass_pr <65) || (jet_mass_pr > 105 && jet_mass_pr < 130)) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ((jet_mass_pr > 40 && jet_mass_pr <65) || (jet_mass_pr > 105 && jet_mass_pr < 130)) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if( preselectionCutType == "basicVBFSRHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSRHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSRHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if( preselectionCutType == "basicVBFSRHiggsTop" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f ) && mass_ungroomedjet_closerjet > 200 && mass_leptonic_closerjet > 200",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSRHiggsTop" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f ) && mass_ungroomedjet_closerjet > 200 && mass_leptonic_closerjet > 200",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSRHiggsTop" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f ) && mass_ungroomedjet_closerjet > 200 && mass_leptonic_closerjet > 200",pTJetMin_,pTJetMax_);

      
  else if( preselectionCutType == "basicVBFSBSRHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 40 && jet_mass_pr <130)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBSRHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 40 && jet_mass_pr <130)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBSRHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 40 && jet_mass_pr <130) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);
      
  //--------------------------
  // Basic preselection CSA14
  //--------------------------

  else if( preselectionCutType == "basicJetsCutCSA14" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "electron" || LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "Jets" || LeptonType == "jets") and  TreeName_ !="gen")
    return Form("pt[0]>200 && abs(eta[0])<2.5 && imatch[0] >= 0 && (pt[0] > %f  && pt[0] < %f ) && (npu > %f && npu < %f)",pTJetMin_,pTJetMax_,npuMin_,npuMax_);
  else if ( preselectionCutType == "basicJetsCutCSA14" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "electron" || LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "Jets" || LeptonType == "jets") and  TreeName_ =="gen")
    return Form("pt[0] > 200 && abs(eta[0])<2.5 && (pt[0] > %f  && pt[0] < %f )  && (npu > %f && npu < %f)",pTJetMin_,pTJetMax_,npuMin_,npuMax_);
 
  else return Form("v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);




}



// print Training results 

void  TrainingMVAClass::PrintTrainingResults (){

  std::string command = " if [ ! -e plots ] ; then mkdir plots ; fi";
  int result = system(command.c_str());
  if(!result) return ;
  std::cout << "******************************************************* "<<std::endl;  
  std::cout << "==> Print Output Plots For: " << outputFile_->GetName() << std::endl;
  std::cout << "******************************************************* "<<std::endl;  

  std::string ROOTStyle =  getenv ("ROOTStyle");

  gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());

  std::string PWD = getenv("PWD");

  // Run macros for plots

  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<" #################################################### "<<std::endl;
  std::cout<<" Run Variables Macro : /macros/TMVAMacro/variables.C  "<<std::endl;
  std::cout<<" #################################################### "<<std::endl;
  std::cout<<std::endl;
 
  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/variables.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  result = system(command.c_str());
  if(!result) return ;


  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<" ############################################################## "<<std::endl;
  std::cout<<" Run Variables Macro : /macros/TMVAMacro/correlationscatters.C  "<<std::endl;
  std::cout<<" ############################################################## "<<std::endl;
  std::cout<<std::endl;
 
  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/correlationscatters.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  result = system(command.c_str());
  if(!result) return ;

  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<" ####################################################### "<<std::endl;
  std::cout<<" Run Variables Macro : /macros/TMVAMacro/correlations.C  "<<std::endl;
  std::cout<<" ####################################################### "<<std::endl;
  std::cout<<std::endl;
 
  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/correlations.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  result = system(command.c_str());
  if(!result) return ;

  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<" ############################################### "<<std::endl;
  std::cout<<" Run Variables Macro : /macros/TMVAMacro/mvas.C  "<<std::endl;
  std::cout<<" ############################################### "<<std::endl;
  std::cout<<std::endl;
 
  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/mvas.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  result = system(command.c_str());
  if(!result) return ;

  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<" ####################################################### "<<std::endl;
  std::cout<<" Run Variables Macro : /macros/TMVAMacro/efficiencies.C  "<<std::endl;
  std::cout<<" ####################################################### "<<std::endl;
  std::cout<<std::endl;
 
  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/efficiencies.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  result = system(command.c_str());
  if(!result) return ;

  // Output directory for plots

  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<" ############################################ "<<std::endl;
  std::cout<<"  Move Files in the correct output directory  "<<std::endl;
  std::cout<<" ############################################ "<<std::endl;
  std::cout<<std::endl;
 
  command = " if [ ! -e "+outputFilePath_+"/trainingPlots ] ; then mkdir "+outputFilePath_+"/trainingPlots ; fi";
  result = system(command.c_str());
  if(!result) return ;
  std::cout<<" command : "<<command<<std::endl;
  std::cout<<std::endl;

  command = " if [ ! -e "+outputFilePath_+"/trainingPlots/"+outputFileName_+"_"+Label_+" ] ; then mkdir "+outputFilePath_+"/trainingPlots/"+
                         outputFileName_+"_"+Label_+" ; fi";
  result = system(command.c_str());
  if(!result) return ;
  std::cout<<" command : "<<command<<std::endl;
  std::cout<<std::endl;
  
  command = "mv ./plots/* "+outputFilePath_+"/trainingPlots/"+outputFileName_+"_"+Label_+"/" ;
  result = system(command.c_str());
  if(!result) return ;
  std::cout<<" command : "<<command<<std::endl;
  std::cout<<std::endl;

  // Output Directory for weights -> one directory for all the methods

  command = " if [ ! -e "+outputFilePath_+"/weights ] ; then mkdir "+outputFilePath_+"/weights ; fi";
  result = system(command.c_str());
  if(!result) return ;
  std::cout<<" command : "<<command<<std::endl;
  std::cout<<std::endl;

  std::map<std::string,std::string>::const_iterator itMap = outputFileWeightName_.begin();

  for( ; itMap!= outputFileWeightName_.end(); ++itMap){

  command = "scp "+itMap->second+"/* "+outputFilePath_+"/weights/";
  result = system(command.c_str());
  if(!result) return ;
  std::cout<<" command : "<<command<<std::endl;
  std::cout<<std::endl;

  command = "rm -rf "+itMap->second ;
  result = system(command.c_str());
  if(!result) return ;
  std::cout<<" command : "<<command<<std::endl;
  std::cout<<std::endl;

  }

  command = "rm -rf ./plots" ;
  result = system(command.c_str());
  if(!result) return ;
  std::cout<<" command : "<<command<<std::endl;
  std::cout<<std::endl;

  std::cout << "==> Wrote image files: " << outputFilePath_+"/trainingPlots/"+outputFileName_+"_"+Label_ << std::endl;
  std::cout << "==> TMVA Plots are done!" << std::endl;
  std::cout<<std::endl;

  return ;
   
}
