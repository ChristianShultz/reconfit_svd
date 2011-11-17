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
