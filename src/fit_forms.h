#ifndef __FIT_FORMS_H__
#define __FIT_FORMS_H__

#include <iostream>
#include <iomanip>

#include "jackknife_fitter.h"
#include "ensem_data.h"


//**********************************
// CONSTANT
//**********************************

class Const : public FitFunction {
 public:
 Const() : FitFunction(1)
    {
      setParName(0, "a");
    }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "constant";}
};

double Const::operator()(const vector<double>& pars, double x) const{
  return pars[0];
};

//**********************************
// CONSTANT + EXP
//**********************************

class ConstPlusExp : public FitFunction {
 public:
 ConstPlusExp() : FitFunction(3)
    {
      setParName(0, "a");
      setParName(1, "b");
      setParName(2, "exp_mass");
    }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "constant_plus_exp";}
};

double ConstPlusExp::operator()(const vector<double>& pars, double x) const{
  return pars[0] + pars[1] * exp( - pars[2] * x );
};


//*****************************************************************************
// THREE POINT FIT FORMS

//constant + one exp at source (mass supplied) + one exp at sink (mass fitted)
/*class OneExpLinearSourceOneExpSink : public FitFunction {
 public:
 LinearSourceOneExpSink(int tsnk, EnsemReal W0, EnsemReal W1, EnsemReal E1) : FitFunction(4)
    {
      setParName(0, "f0");
      setParName(1, "f1_src");
      setParName(2, "f1_snk");
      setParName(3, "snk_mass");
    }
  inline double operator()(const vector<double
*/

//*****************************************************************************

//this is an example where Q*NDoF of the fit is returned (to be maximised)
class CompareFitsByQN : public FitComparator{
 public:
 CompareFitsByQN() : FitComparator(){};
  
  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const;
};
double CompareFitsByQN::operator()(const FitDescriptor& fitDesc, const JackFit& fit) const{
  return fit.getNDoF()* statQ( fit.getAvgChisq(), fit.getNDoF() );
};


class CompareZFits : public FitComparator{
public:
  CompareZFits() : FitComparator(){};

  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const;
};

double CompareZFits::operator()(const FitDescriptor& fitDesc, const JackFit& fit) const{
  double chisq = fit.getAvgChisq();
  int nDoF = fit.getNDoF();
  double Q = statQ(chisq, nDoF);

  double Z = fit.getAvgFitParValue("a");
  
  double out = Q * nDoF;

  if( (*(fitDesc.ff)).getFitType() == "constant_plus_exp" ){
    double mass = fit.getAvgFitParValue("exp_mass");
    double mass_err = fit.getAvgFitParError("exp_mass");

    double z2 = fit.getAvgFitParValue("b");
    double z2_err = fit.getAvgFitParError("b");

    // a cooked up function that varies between 1 for 1 sigma and 2 for infinite sigma
    double x0 = 2.0;
    double c = 3.14159/2.0 - 2*atan(1.0 - x0);
    double alpha = 2.0/(0.5*3.14159 + c);

    double e =  mass/mass_err;
    double f = alpha * ( atan( e - x0 ) + c );

    out *= f;

    if( mass < mass_err ){ out = 0.0; }
    if( fabs(z2) < z2_err ){ out = 0.0; }

    // introduce an ABSOLUTE scale for the exponential mass - just HARDWIRE here
    if( mass < 0.1 ){ out = 0.0; }
  }
  return out;
}


//********************************************************************************
class FitZ{
 public:
  FitZ(EnsemData data_, int t0, Handle<FitComparator> fitComp_, int minTSlices);
  
  void saveFitPlot(string filename);
  string getFitPlotString() const {return axis_plot;};
  string getFitSummary() const {return fit_summary;};

  string getFitType() const {return fit_type;};

  EnsemReal getZ() const {return Z;};
  EnsemReal getZ2() const 
  {
    if(fit_type != "constant"){return Z2;}
    else{cerr << "constant fit was best" << endl; exit(1);}; 
  };
  EnsemReal getExpMass()  const 
  {
    if(fit_type != "constant"){return exp_mass;}
    else{cerr << "constant fit was best" << endl; exit(1);}; 
  };

  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};
  string getFitName() const {return best_fit_name;};

  
 private:
  Handle<FitComparator> fitComp;
  JackFitLog fits;

  string fit_summary;
  string axis_plot;
  string fit_type;

  double chisq;
  double nDoF;
  string best_fit_name;

  EnsemReal Z;
  EnsemReal Z2;
  EnsemReal exp_mass;

};

//********************************************************************************
class FitVersusT0{
public:
  FitVersusT0(EnsemData data_, Handle<FitComparator> fitComp_ );
	
  void saveFitPlot(string filename);
  string getFitPlotString() const {return axis_plot;};
  string getFitSummary() const {return fit_summary;};
  
  string getFitType() const {return fit_type;};
  
  EnsemReal getConst() const {return a;};
  
  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};
  string getFitName() const {return best_fit_name;};
  
	
private:
  Handle<FitComparator> fitComp;
  JackFitLog fits;
  
  string fit_summary;
  string axis_plot;
  string fit_type;
  
  double chisq;
  double nDoF;
  string best_fit_name;
  
  EnsemReal a;	
};




#endif
