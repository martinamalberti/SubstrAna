#ifndef QGLikelihoodCalculator_h
#define QGLikelihoodCalculator_h

#include <string>

#include "TFile.h"
#include "TH1F.h"
#include <map>



class QGLikelihoodCalculator {

 public:
  QGLikelihoodCalculator( const std::string dataDir, Bool_t chs = false);
  ~QGLikelihoodCalculator();

  float computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );
  float computeQGLikelihood2012( float pt, float eta, float rho, int nPFCandidates_QC, float ptD_QC, float axis2_QC ); //new
  Float_t QGvalue(std::map<std::string, Float_t>);

  float likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand);

 private:

  TFile* histoFile_;
  std::map<std::string,TH1F*> plots_;
  unsigned int nPtBins_;
  unsigned int nRhoBins_;

};


#endif
