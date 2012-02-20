#ifndef __ENSEM_DATA_H__
#define __ENSEM_DATA_H__

#include <vector>
#include <algorithm>
#include "ensem/ensem.h"
#include <cassert>
#include <itpp/itbase.h>
#include <iostream>

#include "linear_algebra.h"

using namespace ENSEM;
using namespace std;


//this class takes a double argument and returns an EnsemReal
//the constructor takes some EnsemReal params
//this is the base class
class EnsemFunction{  
 public: 
  EnsemFunction(){};
  ~EnsemFunction(){};

  virtual EnsemReal operator()(double x) const = 0;
};
//here's a concrete example
class ConstTimesExp : public EnsemFunction{
 public:
 ConstTimesExp(EnsemReal weight_, EnsemReal mass_) : weight(weight_), mass(mass_), EnsemFunction(){}
  
  EnsemReal operator()(double x) const{ return weight * exp( mass * Real(x) );}
    
 private:
  EnsemReal weight, mass;
};
//here's a rather trivial example
class One : public EnsemFunction{
 public:
 One(EnsemReal anything_) : anything(anything_), EnsemFunction(){};
  
  EnsemReal operator()(double x) const { 
    EnsemReal one; one.resize(anything.size()); one = Real(1.0); return one; }

 private:
  EnsemReal anything;
};
//
class ConstDivideCosh : public EnsemFunction{
 public:
  ConstDivideCosh(EnsemReal weight_, EnsemReal mass_, int T_) : weight(weight_), mass(mass_), T(T_), EnsemFunction(){}
  
    EnsemReal operator()(double x) const{ 
      //      return weight / ( exp(- mass * Real( double(T/2) ) ) * cosh( mass * Real(x - T/2.0 ) )  ) ;
      EnsemReal cosh = Real(0.5)*(  exp(  mass * Real(x - T/2.0 ) ) + exp( -  mass * Real(x - T/2.0 ) ) );
      return weight / ( Real(2.0) * exp(- mass * Real( double(T/2) )) * cosh );
    }
    
 private:
  EnsemReal weight, mass;
  int T;
};


//****************************************
//some search functions
bool double_predicate(double a, double b);
int find_nearest(vector<double> v, double a);
int find_exact(vector<double> v, double a);
vector<int> find_below(vector<double> v, double a);
vector<int> find_above(vector<double> v, double a);
vector<int> find_in_range(vector<double> v, double low, double high);
//****************************************


//****************************************
itpp::mat makeEnsemCov(EnsemVectorReal y);
//****************************************

//****************************************
// Ensem y values, double x values
class EnsemData{
 public:
  EnsemData()
    {totalNData = 0; nBins = 0; initP = false; initCov = false; initInvCov = false; SVCutoff = 1.0e-6;};

  EnsemData(vector<double> x, EnsemVectorReal y);

  //  EnsemData(const EnsemData& e);
  // EnsemData& operator=(const EnsemData& e); copy constructors should be there automatically

  void addDatum(double x, EnsemReal y);

  //show/hide
  void hideDatumByIndex(int i);
  void hideDataByIndex(vector<int> list);
  
  bool hideDatumByX(double x);
  bool showDatumByX(double x);

  int hideDataByXRange(double xlow, double xhigh);
  int showDataByXRange(double xlow, double xhigh);
  int hideDataBelowX(double x);
  int showDataBelowX(double x);
  int hideDataAboveX(double x);
  int showDataAboveX(double x);

  int hideDataByYRange(double ylow, double yhigh);
  int showDataByYRange(double ylow, double yhigh);
  int hideDataBelowY(double y);
  int showDataBelowY(double y);
  int hideDataAboveY(double y);
  int showDataAboveY(double y);

  int hideDataAboveYErr(double err);
  int showDataAboveYErr(double err);
  int hideDataBelowYErr(double err); 
  int showDataBelowYErr(double err);
  
  int hideDataAboveYErrRat(double rat);
  int showDataAboveYErrRat(double rat);
  int hideDataBelowYErrRat(double rat);
  int showDataBelowYErrRat(double rat);

  void hideAll();
  void showAll();

  //get at innards
  int getNData() const;
  int getTotalNData() const {return totalNData;}
  int getNBins() const {return nBins;}

  //active data calls
  vector<double> getXData() const;  
  EnsemVectorReal getYData() const; 
  EnsemVectorReal getScaledYData() const;
  itpp::mat getCov() const;

  //total data calls
  vector<double> getAllXData() const {return x_data;}; 
  EnsemVectorReal getAllYData() const {return y_data;}; 
  vector<bool> getActiveDataList() const {return active_data;};
  itpp::mat getAllCov() const;
  
  //pull out single values
  EnsemReal getYUsingNearestX(double& x);

  //covariance calls
  void makeCov() const {
    if(!initCov){cov = makeEnsemCov(y_data); initCov = true;};
    //    cout << "making covariance" << endl;
  };
  bool isCovMade() const {return initCov;};
  bool setSVCutoff(double cutoff){SVCutoff = cutoff;};
  double getSVCutoff() const {return SVCutoff;};

  //inverse covariance
  void makeInvCov() const;
  void makeInvCov(double cutoff) const;
  itpp::mat getInvCov() const;
  int getNResetCovSingVals() const;

 private:
  EnsemVectorReal y_data; //all the ensemData
  vector<double> x_data; //all the x-values
  vector<bool> active_data; //which data points are switched on
  
  int nBins;
  int totalNData;  //total num data points stored 

  bool initP; //has any data been inserted yet
  
  mutable itpp::mat cov; // the big matrix - mutable should mean I can generate this by a 'getCov()' call
  mutable itpp::mat invCov; //inverse covarince for the active data
  mutable int nResetCovSingVals;

  void augmentCov(EnsemReal y); 
  mutable bool initCov; // has the big cov been made?
  mutable bool initInvCov; //has the covariance been inverted on this active data set
  double SVCutoff;
};  
//*****************************************************************

EnsemData concat(EnsemData a, EnsemData b);
ostream& operator<<(ostream& output, const EnsemData& e);


//******************************************************************
double ensemDataChisq(const EnsemData& data, const EnsemData& theory);


#endif


















