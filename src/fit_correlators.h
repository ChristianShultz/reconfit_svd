#ifndef __FIT_CORRELATORS_H__
#define __FIT_CORRELATORS_H__

#include <iostream>
#include <iomanip>

#include "jackknife_fitter.h"
#include "ensem_data.h"

//////////////////////////////////////////////////
// define some concrete versions of FitFunction
// relavant for fitting timeslice data
//////////////////////////////////////////////////


//******************************
// ONE EXP - PRIN CORR
//******************************
class PrinCorrOneExp : public FitFunction {
 public:
 PrinCorrOneExp(int t0_) : t0(t0_), FitFunction(1)
    {
      setParName(0, "mass_0");
    }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "PrinCorrOneExp";}

 private:
  int t0;
};

double PrinCorrOneExp::operator()(const vector<double>& pars, double x) const{
  return exp( - pars[0] * ( x - double(t0) ) );
};

//******************************
// TWO EXP - PRIN CORR
//******************************
class PrinCorrTwoExp : public FitFunction {
 public:
 PrinCorrTwoExp(int t0_) : t0(t0_), FitFunction(3)
    {
      setParName(0, "mass_0");
      setParName(1, "mass_1");
      setParName(2, "A");
    }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "PrinCorrTwoExp";}

 private:
  int t0;
};

double PrinCorrTwoExp::operator()(const vector<double>& pars, double x) const{
  //params are stored as {mass_0, mass_1, amp}
  // (1-A) e^{- m_0 * (t-t0)} + Ae^{- m_1 * (t-t0)}
  return (1.0 - pars[2]) * exp( - pars[0] * (x - double(t0) ) ) 
    + pars[2] * exp( - pars[1] * (x - double(t0) ) );
};

//******************************
// ONE EXP - CORR
//******************************
class CorrOneExp : public FitFunction {
 public:
 CorrOneExp() : FitFunction(2)
  {
    setParName(0, "mass_0");
    setParName(1, "amp_0");
  }
 inline double operator()(const vector<double>& pars, double x) const;
 string getFitType() const {return "CorrOneExp";}
};

double CorrOneExp::operator()(const vector<double>& pars, double x) const{
  //params are stored as {mass_0, amp}
  // A e^{- m_0 * t} 
  return pars[1] * exp( - pars[0] * x );
};

//******************************
// TWO EXP - CORR
//******************************
class CorrTwoExp : public FitFunction {
 public:
 CorrTwoExp() : FitFunction(4)
  {
    setParName(0, "mass_0");
    setParName(1, "amp_0");
    setParName(2, "mass_1");
    setParName(3, "amp_1");
  }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "CorrTwoExp";}
};

double CorrTwoExp::operator()(const vector<double>& pars, double x) const{
  //params are stored as {mass_0, amp_0, mass_1, amp_1}
  // A_0 e^{- m_0 * t} + A_1 e^{- m_1 * t}
  return pars[1] * exp( - pars[0] * x ) + pars[3] * exp( - pars[2] * x );
};

//******************************
// ONE COSH - CORR
//******************************
class CorrOneCosh : public FitFunction {
 public:
 CorrOneCosh(int T_) : T(T_) , FitFunction(2)
  {
    setParName(0, "mass_0");
    setParName(1, "amp_0");
  }
 inline double operator()(const vector<double>& pars, double x) const;
 string getFitType() const {return "CorrOneCosh";}

 private:
  int T;
};

double CorrOneCosh::operator()(const vector<double>& pars, double x) const{
  //params are stored as {mass_0, amp}
  // A * exp( - m * T / 2 ) * 2 * cosh( m * (t - T/2) ) 
  return pars[1] * exp( - pars[0] * double(T / 2) ) * 2 * cosh( pars[0] * ( x - double(T/2) ) );
}; //NB as T-> infty this form tends to A * exp( - m * t)

//******************************
// TWO COSH - CORR
//******************************
class CorrTwoCosh : public FitFunction {
 public:
 CorrTwoCosh(int T_) : T(T_) , FitFunction(4)
  {
    setParName(0, "mass_0");
    setParName(1, "amp_0");
    setParName(2, "mass_1");
    setParName(3, "amp_1");
  }
 inline double operator()(const vector<double>& pars, double x) const;
 string getFitType() const {return "CorrTwoCosh";}

 private:
  int T;
};

double CorrTwoCosh::operator()(const vector<double>& pars, double x) const{
  //params are stored as {mass_0, amp_0, mass_1, amp_1}

  // A0 * exp( - m0 * T / 2 ) * 2 * cosh( m0 * (t - T/2) ) + A1 * ... 
  return pars[1] * exp( - pars[0] * double(T / 2) ) * 2 * cosh( pars[0] * ( x - double(T/2) ) ) 
    + pars[3] * exp( - pars[2] * double(T / 2) ) * 2 * cosh( pars[2] * ( x - double(T/2) ) );
};

//******************************
// ONE EXP & CONSTANT - CORR
//******************************
class CorrOneExpAndConst : public FitFunction {
 public:
 CorrOneExpAndConst() : FitFunction(3)
  {
    setParName(0, "const");
    setParName(1, "mass_0");
    setParName(2, "amp_0");
  }
 inline double operator()(const vector<double>& pars, double x) const;
 string getFitType() const {return "CorrOneExpAndConst";}
};

double CorrOneExpAndConst::operator()(const vector<double>& pars, double x) const{
  //params are stored as {const, mass_0, amp}
  // c + A e^{- m_0 * t} 
  return pars[0] + pars[2] * exp( - pars[1] * x );
};

//******************************
// TWO EXP & CONSTANT - CORR
//******************************
class CorrTwoExpAndConst : public FitFunction {
 public:
 CorrTwoExpAndConst() : FitFunction(5)
  {
    setParName(0, "const");
    setParName(1, "mass_0");
    setParName(2, "amp_0");
    setParName(3, "mass_1");
    setParName(4, "amp_1");
  }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "CorrTwoExpAndConst";}
};

double CorrTwoExpAndConst::operator()(const vector<double>& pars, double x) const{
  //params are stored as {const, mass_0, amp_0, mass_1, amp_1}
  // c + A_0 e^{- m_0 * t} + A_1 e^{- m_1 * t}
  return pars[0] + pars[2] * exp( - pars[1] * x ) + pars[4] * exp( - pars[3] * x );
};

//******************************
// ONE COSH & CONST- CORR
//******************************
class CorrOneCoshAndConst : public FitFunction {
 public:
 CorrOneCoshAndConst(int T_) : T(T_) , FitFunction(3)
  {
    setParName(0, "const");
    setParName(1, "mass_0");
    setParName(2, "amp_0");
  }
 inline double operator()(const vector<double>& pars, double x) const;
 string getFitType() const {return "CorrOneCoshAndConst";}

 private:
  int T;
};

double CorrOneCoshAndConst::operator()(const vector<double>& pars, double x) const{
  //params are stored as {const, mass_0, amp}
  // c + A * exp( - m * T / 2 ) * 2 * cosh( m * (t - T/2) ) 
  return pars[0] + pars[2] * exp( - pars[1] * double(T / 2) ) * 2 * cosh( pars[1] * ( x - double(T/2) ) );
}; //NB as T-> infty this form tends to A * exp( - m * t)

//******************************
// TWO COSH - CORR
//******************************
class CorrTwoCoshAndConst : public FitFunction {
 public:
 CorrTwoCoshAndConst(int T_) : T(T_) , FitFunction(5)
  {
    setParName(0, "const");  
    setParName(1, "mass_0");
    setParName(2, "amp_0");
    setParName(3, "mass_1");
    setParName(4, "amp_1");
  }
 inline double operator()(const vector<double>& pars, double x) const;
 string getFitType() const {return "CorrTwoCoshAndConst";}

 private:
  int T;
};

double CorrTwoCoshAndConst::operator()(const vector<double>& pars, double x) const{
  //params are stored as {const, mass_0, amp_0, mass_1, amp_1}

  // c + A0 * exp( - m0 * T / 2 ) * 2 * cosh( m0 * (t - T/2) ) + A1 * ... 
  return pars[0] + pars[2] * exp( - pars[1] * double(T / 2) ) * 2 * cosh( pars[1] * ( x - double(T/2) ) ) 
    + pars[4] * exp( - pars[3] * double(T / 2) ) * 2 * cosh( pars[3] * ( x - double(T/2) ) );
};






//********************************************************************
class FitPrincipalCorrelator{
 public:
  FitPrincipalCorrelator(EnsemData data_, int t0_,  Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices);
  
  void saveFitPlot(string filename);
  string getFitPlotString() const {return axis_plot;};
  string getFitSummary() const {return fit_summary;};

  int getNExp() const {return nExp;};

  EnsemReal getMass0() const {return mass_0;};
  EnsemReal getMass1() const {if(nExp==2){return mass_1;}else{cerr << "oneExp fit was best" << endl; exit(1);}; };
  EnsemReal getA()  const {if(nExp==2){return A;}else{cerr << "oneExp fit was best" << endl; exit(1);}; };

  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};
  string getFitName() const {return best_fit_name;};

  
 private:
  int t0;
  Handle<FitComparator> fitComp;
  JackFitLog fits;

  string fit_summary;
  string axis_plot;
  int nExp;
  double chisq;
  double nDoF;
  string best_fit_name;

  EnsemReal mass_0;
  EnsemReal mass_1;
  EnsemReal A;
};
//********************************************************************
class FitCorrelatorExp{
 public:
  FitCorrelatorExp(EnsemData data_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices);
  
  void saveFitPlot(string filename){ofstream out; out.open(filename.c_str());  out << axis_plot;  out.close(); }
  string getFitPlotString() const {return axis_plot;};
  string getFitSummary() const {return fit_summary;};

  int getNExp() const {return nExp;};

  EnsemReal getMass0() const {return mass_0;};
  EnsemReal getMass1() const {if(nExp==2){return mass_1;}else{cerr << "oneExp fit was best" << endl; exit(1);}; };
  EnsemReal getAmp0()  const {return amp_0;};
  EnsemReal getAmp1()  const {if(nExp==2){return amp_1;}else{cerr << "oneExp fit was best" << endl; exit(1);}; };

  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};
  string getFitName() const {return best_fit_name;};

  
 private:
  Handle<FitComparator> fitComp;
  JackFitLog fits;

  string fit_summary;
  string axis_plot;
  int nExp;
  double chisq;
  double nDoF;
  string best_fit_name;

  EnsemReal mass_0;
  EnsemReal mass_1;
  EnsemReal amp_0;
  EnsemReal amp_1;

};
//**************************************************
class FitCorrelatorCosh{
 public:
  FitCorrelatorCosh(EnsemData data_, int T_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices);
  
  void saveFitPlot(string filename);
  string getFitPlotString() const {return axis_plot;};
  string getFitSummary() const {return fit_summary;};

  int getNExp() const {return nExp;}; //obviously actually the number of coshes

  EnsemReal getMass0() const {return mass_0;};
  EnsemReal getMass1() const {if(nExp==2){return mass_1;}else{cerr << "oneCosh fit was best" << endl; exit(1);}; };
  EnsemReal getAmp0()  const {return amp_0;};
  EnsemReal getAmp1()  const {if(nExp==2){return amp_1;}else{cerr << "oneCosh fit was best" << endl; exit(1);}; };

  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};
  string getFitName() const {return best_fit_name;};

  
 private:
  int T;  // the time extent of the lattice

  Handle<FitComparator> fitComp;
  JackFitLog fits;

  string fit_summary;
  string axis_plot;
  int nExp;
  double chisq;
  double nDoF;
  string best_fit_name;

  EnsemReal mass_0;
  EnsemReal mass_1;
  EnsemReal amp_0;
  EnsemReal amp_1;

};
//**********************************************************
class FitCorrelatorExpAndConst{
 public:
  FitCorrelatorExpAndConst(EnsemData data_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices);
  
  void saveFitPlot(string filename){ofstream out; out.open(filename.c_str());  out << axis_plot;  out.close(); }
  string getFitPlotString() const {return axis_plot;};
  string getFitSummary() const {return fit_summary;};

  int getNExp() const {return nExp;};

  EnsemReal getConst() const {return constant;};
  EnsemReal getMass0() const {return mass_0;};
  EnsemReal getMass1() const {if(nExp==2){return mass_1;}else{cerr << "oneExp fit was best" << endl; exit(1);}; };
  EnsemReal getAmp0()  const {return amp_0;};
  EnsemReal getAmp1()  const {if(nExp==2){return amp_1;}else{cerr << "oneExp fit was best" << endl; exit(1);}; };

  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};
  string getFitName() const {return best_fit_name;};

  
 private:
  Handle<FitComparator> fitComp;
  JackFitLog fits;

  string fit_summary;
  string axis_plot;
  int nExp;
  double chisq;
  double nDoF;
  string best_fit_name;

  EnsemReal constant;
  EnsemReal mass_0;
  EnsemReal mass_1;
  EnsemReal amp_0;
  EnsemReal amp_1;

};

//**********************************************************
class FitCorrelatorCoshAndConst{
 public:
  FitCorrelatorCoshAndConst(EnsemData data_, int T_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices);
  
  void saveFitPlot(string filename){ofstream out; out.open(filename.c_str());  out << axis_plot;  out.close(); }
  string getFitPlotString() const {return axis_plot;};
  string getFitSummary() const {return fit_summary;};

  int getNExp() const {return nExp;};

  EnsemReal getConst() const {return constant;};
  EnsemReal getMass0() const {return mass_0;};
  EnsemReal getMass1() const {if(nExp==2){return mass_1;}else{cerr << "oneCosh fit was best" << endl; exit(1);}; };
  EnsemReal getAmp0()  const {return amp_0;};
  EnsemReal getAmp1()  const {if(nExp==2){return amp_1;}else{cerr << "oneCosh fit was best" << endl; exit(1);}; };

  double getChisq() const{ return chisq;};
  double getNDoF() const{ return nDoF;};
  string getFitName() const {return best_fit_name;};

  
 private:
  Handle<FitComparator> fitComp;
  JackFitLog fits;

  int T;

  string fit_summary;
  string axis_plot;
  int nExp;
  double chisq;
  double nDoF;
  string best_fit_name;

  EnsemReal constant;
  EnsemReal mass_0;
  EnsemReal mass_1;
  EnsemReal amp_0;
  EnsemReal amp_1;

};





//****************************************
//   FIT COMPARATORS
//****************************************

//this is an example where the "splitN" of the fit is returned (to be maximised) - should be able to use this for Prin & regular corrs
class CompareFitsBySplitN : public FitComparator{
 public:
 CompareFitsBySplitN() : FitComparator(){};
  
  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const;
};
double CompareFitsBySplitN::operator()(const FitDescriptor& fitDesc, const JackFit& fit) const{
  double dum =  fit.getNDoF()* statQ( fit.getAvgChisq(), fit.getNDoF() );
  
  vector<double> pars = fit.getAvgFitParValues();
  vector<double> errs = fit.getAvgFitParErrors();

  dum *= pars[0]/errs[0];
  if( (*(fitDesc.ff)).getFitType() == "PrinCorrTwoExp" ){
    dum *= (pars[1] - pars[0]) / errs[0];
  }
  // ADD MORE OPTIONS FOR REGULAR CORRS

  return dum;
};

//this is an example where a generic fitCrit of the fit is returned (to be maximised) - should be able to use this for Prin & regular corrs
class CompareFitsByGeneric : public FitComparator{
 public:
 CompareFitsByGeneric() : FitComparator(){};
  
  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const;
};
double CompareFitsByGeneric::operator()(const FitDescriptor& fitDesc, const JackFit& fit) const{
  CompareFitsBySplitN splitN;
  double dum = splitN(fitDesc, fit);

  vector<double> pars = fit.getAvgFitParValues();
  vector<double> errs = fit.getAvgFitParErrors();
  
  if( (*(fitDesc.ff)).getFitType() == "PrinCorrTwoExp" ){
    if(pars[1] < pars[0]){dum = 0.0;}; //second exp mass better be higher than the first exp mass
    if(errs[1] > 2.0*pars[1]){dum = 0.0;} //very noisy second exps are not allowed
    if((pars[2] > 0.5)||(pars[2] < 0.0)){dum = 0.0;}; //no large or negative second exp amplitudes
  }
  return dum;
};







#endif
