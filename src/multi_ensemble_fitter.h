// -*- C++ -*-
#ifndef __MULTI_ENSEM_FITTER_H__
#define __MULTI_ENSEM_FITTER_H__

#include "jackFitter/jackFitter.h"
  
// need to read in ensem values from an ensemble
// compute the covariance

// do this for multiple ensembles
// build a block diagonal covariance

class MultiEnsemFitter{
 public:
  MultiEnsemFitter(map<string, EnsemData> data_, Handle<FitFunction> ff_);
  
  bool fit(); // fit all the data
  //  bool fitEnsemble(string ensem); // fit a single ensemble
  //  bool fitEnsembles(vector<string> ensems); // fit some ensembles
  
  double getFitParValue(string name) {return fitParValue[name];};
  double getFitParError(string name) {return fitParError[name];};
  
  double getChisq() const {return chisq;};
  int getNDoF() const {return nDoF;};   
  
  bool getFitSuccess() const {return fitSuccess;};
  string getFitReport() const {return fitReport;};
  
  string makeFitPlotAxis(double xmin, double xmax, string label); // plot should discriminate between ensembles
  
 private:
  map<string, EnsemData> data; // string names the ensemble
  Handle<FitFunction> ff;
  
  vector<double> vFitParValue;
  vector<double> vFitParError;

  map<string, double> fitParValue;
  map<string, double> fitParError;
  
  double chisq;
  int nDoF;
  
  bool fitSuccess;
  string fitReport;
  
};

// do a chisq minimisation for some fit function

// the minimisation is not jackknife

// can re-use class ChiSquare
// and base class FitFunction

// just don't call them from JackknifeFitter



#endif
