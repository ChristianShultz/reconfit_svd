#ifndef __MULTI_PRINCORR_FIT_H__
#define __MULTI_PRINCORR_FIT_H__


#include <iostream>
#include <iomanip>

#include "jackknife_fitter.h"
#include "ensem_data.h"


//pars are
//pars[0] = common mass
//pars[1,3,5...] = mass' for corr_0,1,2
//pars[2,4,6...] = A for corr_0,1,2
 

class FitMultiplePrincipalCorrelator{
 public:
  FitMultiplePrincipalCorrelator(EnsemData data_, int nCorrs_, vector<int> t0_list_, int tLen,  vector<int> tmin_,  double noiseRatioCutoff);
  
  void saveFitPlot(string filename);
  string getFitPlotString() const {return axis_plot;};
  
  EnsemReal getMass0() const {return mass_0;};
  
  vector<EnsemReal> getMass1() const {return mass_1;};
  vector<EnsemReal> getA()       const {return A;};
  
  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};

  
  
 private:
  int nCorrs;
  vector<int> t0_list;
  vector<int> tmin;

  int tLen;
  

  EnsemData data;
  
  double chisq;
  double nDoF;
  string axis_plot;
  
  EnsemReal mass_0;
  vector<EnsemReal> mass_1;
  vector<EnsemReal> A;
  
};



//fit function:
class MultiPrinCorrTwoExp : public FitFunction {
 public:
 MultiPrinCorrTwoExp(int nCorrs_, vector<int> t0_list_, int tLen_) : nCorrs(nCorrs_), t0_list(t0_list_), tLen(tLen_), FitFunction( 2*nCorrs_ + 1 )
    {
      //      cout << "constructing MultPrinCorrTwoExp" << endl;
      setParName(0, "mass_0");
      for(int n=0; n < nCorrs; n++){
	stringstream dum; dum <<  "m_1_corr" << n;
	setParName(1 + 2*n , dum.str());
	
	stringstream dum2; dum2 << "A_corr" << n;
	setParName(2 + 2*n , dum2.str());
      }

      //      cout << "constructed MultiPrinCorrTwoExp" << endl;
      
    }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "MultiPrinCorrTwoExp";}

 private:
  int nCorrs;
  vector<int> t0_list;
  int tLen;
};


double MultiPrinCorrTwoExp::operator()(const vector<double>& pars, double x) const{
  
  double m = pars[0];
  double f;
  
  for(int n = 0; n < nCorrs; n++){
    if( (x < (n + 1)*tLen) && (x >= n*tLen) ){
      int t0 = t0_list[n];
      double m2 = pars[2*n + 1];
      double A = pars[2*n + 2];
      
      f = ( 1.0 - A ) * exp( - m * (x - double(t0) ) ) + A * exp( - m2 * (x - double(t0) ) );
      //   cout << "f(" << x <<") = " << f << endl;
    }
  } 
  return f; 
};

class MultiCorrWeight : public EnsemFunction{
 public:
 MultiCorrWeight(EnsemReal mass_, int nCorrs_, vector<int> t0_list_, int tLen_) : nCorrs(nCorrs_), mass(mass_), t0_list(t0_list_), tLen(tLen_), EnsemFunction(){}
  
  EnsemReal operator()(double(x)) const{
      EnsemReal f = Real(0.0)*mass;
      for(int n = 0; n < nCorrs; n++){
	if( (x < (n + 1)*tLen) && (x >= n*tLen) ){
	  int t0 = t0_list[n];
	  f = exp( mass * Real(x - double(t0) ) ); 
	}
      } 
      return f;
    } 
    
    //    EnsemReal operator()(double x) const{ return weight * exp( mass * Real(x) );}
    
 private:
  int nCorrs, tLen;
  vector<int> t0_list;
  EnsemReal mass;
};



#endif
