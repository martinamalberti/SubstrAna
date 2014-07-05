#include "../include/VTaggingVariables.h"

///////////////////////                                                                                                                                                                  
VTaggingVariables::VTaggingVariables(const PseudoJet & inputJet){

  inputJet_ = inputJet ;

}

///////////////////////                                                                                                                                                                  

double VTaggingVariables::computeJetCharge(const double & jetChargeKappa){

  jetChargeKappa_ = jetChargeKappa ;

  double  chargeVal = 0.;
  std::vector<fastjet::PseudoJet> particles ;
  std::vector<fastjet::PseudoJet> ghosts ;

  SelectorIsPureGhost().sift(inputJet_.constituents(), ghosts, particles);
  std::vector<fastjet::PseudoJet>::const_iterator itConstituents =  particles.begin();

  for ( ; itConstituents != particles.end(); ++itConstituents){
    if(fabs((*itConstituents).user_index()) > 2) continue ;
    chargeVal += (*itConstituents).user_index()*pow((*itConstituents).pt(),jetChargeKappa_);
  }

  chargeVal /= pow(inputJet_.pt(),jetChargeKappa_ );
  return chargeVal;

}

///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeJetChargeReco(const double & jetChargeKappa){

  jetChargeKappa_ = jetChargeKappa ;

  double  chargeVal = 0.;
  std::vector<fastjet::PseudoJet> particles ;
  std::vector<fastjet::PseudoJet> ghosts ;

  SelectorIsPureGhost().sift(inputJet_.constituents(), ghosts, particles);
  std::vector<fastjet::PseudoJet>::const_iterator itConstituents =  particles.begin();
  for ( ; itConstituents != particles.end(); ++itConstituents){
    double charge = 0.;        
    
    if(fabs((*itConstituents).user_index()) > 7) continue ;
    if((*itConstituents).user_index() == 0) charge = 0;
    else if(fabs((*itConstituents).user_index() )<=2) charge = (*itConstituents).user_index() ;
    else if((fabs((*itConstituents).user_index() ) > 2 )) charge = (*itConstituents).user_index() -5 ;
    else charge = 0 ;
    
    chargeVal += charge*pow((*itConstituents).pt(),jetChargeKappa_);
  }

  chargeVal /= pow(inputJet_.pt(),jetChargeKappa_ );
  return chargeVal;

}


///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeNSubJettines(const int & nJettines, const double & beta, const double & R0, const double & Rcut){

 if(inputJet_.structure_ptr() == 0 or inputJet_.structure_ptr()== NULL or inputJet_ == 0) return -99;
 contrib::Nsubjettiness nSub(nJettines,contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
 double nSubjettinessValues = nSub(inputJet_);

 return nSubjettinessValues ;
}
  

///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeECF(JetAlgorithm jetAlgoforECF, const double & Rparameter, const int & nPoint, const double & beta){

   double ECFValues = 0 ;
   fastjet::JetDefinition jet_def_forECF(jetAlgoforECF,Rparameter); // definition for ECF
   fastjet::ClusterSequence clust_seq_forECF(inputJet_.constituents(),jet_def_forECF);       // cluster sequence for ECF computation
   std::vector<fastjet::PseudoJet> incluisve_jets_forECF = clust_seq_forECF.inclusive_jets(0);       // take all the set of pseudojets from clustering
   contrib::EnergyCorrelatorRatio C2beta(nPoint,beta,contrib::EnergyCorrelator::pt_R); // calculate the ECF
   ECFValues = C2beta(incluisve_jets_forECF[0]);
  
   return ECFValues ;

}


///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeQjets(const int & QJetsPreclustering, const int & QJetsN, const int & seed){

  std::vector<fastjet::PseudoJet> constits;
  unsigned int nqjetconstits = inputJet_.constituents().size(); // take jet constituent size

  if (nqjetconstits < (unsigned int) QJetsPreclustering) constits = inputJet_.constituents(); // use of subset of particles
  else constits = inputJet_.associated_cluster_sequence()->exclusive_subjets_up_to(inputJet_,QJetsPreclustering); // take the exclusive subjets
  double qjet_vol = getQjetVolatility(constits, QJetsN, seed*QJetsN) ;
  constits.clear();
  return qjet_vol;
}


double VTaggingVariables::getQjetVolatility(std::vector < fastjet::PseudoJet > constits, const int & QJetsN, const int & seed){

    std::vector<double> qjetmasses;
    double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.01);

    for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
      QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
      qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed                                                                                                   
      JetDefinition qjet_def(&qjet_plugin);
      ClusterSequence qjet_seq(constits, qjet_def);
      std::vector<PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));

      if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }

    }

    // find RMS of a vector                                                                                                                                                               
    double qjetsRMS = FindRMS( qjetmasses );
    // find mean of a vector                                                                                                                                                              
    double qjetsMean = FindMean( qjetmasses );
    double qjetsVolatility = qjetsRMS/qjetsMean;
    return qjetsVolatility;
  
}
  
//Q jets stuff                                                                                                                                                                            
double VTaggingVariables::FindRMS( const std::vector< double > & qjetmasses ){

  double total = 0.;
  double ctr = 0.;

  for (unsigned int i = 0; i < qjetmasses.size(); i++){
      total = total + qjetmasses[i];
      ctr++;
  }

  double mean = total/ctr;
  double totalsquared = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
      totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
  }
  double RMS = sqrt( totalsquared/ctr );
  return RMS;
}

double VTaggingVariables::FindMean(const std::vector< double > & qjetmasses){

 double total = 0.;
 double ctr = 0.;
 for (unsigned int i = 0; i < qjetmasses.size(); i++){
  total = total + qjetmasses[i];
  ctr++;
 }
  
 return total/ctr;  
}

///////////////////////                                                                                                                                                                  
void VTaggingVariables::setInputJet(const PseudoJet & inputJet){

  inputJet_ = inputJet ;

}

 
