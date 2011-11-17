#ifndef __LUESCHER_H__
#define __LUESCHER_H__

#include "recipes/nr.h"
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "math.h"
#include "ensem/ensem.h"

using namespace std;
using namespace ENSEM;

/*    LUESCHER W FUNCTIONS   */
/*    NPB354 pp572 eqn (E.1) */

// these functions are divergent at each integer value of qsq,
// use cubic splines between neighbouring integers, more sampled points near the integers
// should be very accurate except within 0.001 of an integer
// might consider implementing more accurate method close to integers

class LuescherW00{
 public:
  LuescherW00(); //constructor sets up the splines

  double operator()(double qsq);

 private:
  vector<NR::Vec_DP> w_qsq; //vector over integer regions
  vector<NR::Vec_DP> w_w;
  vector<NR::Vec_DP> w_ddw; //spline
};

class LuescherW40{
 public:
  LuescherW40(); //constructor sets up the splines

  double operator()(double qsq);

 private:
  vector<NR::Vec_DP> w_qsq; //vector over integer regions
  vector<NR::Vec_DP> w_w;
  vector<NR::Vec_DP> w_ddw; //spline
};

class LuescherW60{
 public:
  LuescherW60(); //constructor sets up the splines

  double operator()(double qsq);

 private:
  vector<NR::Vec_DP> w_qsq; //vector over integer regions
  vector<NR::Vec_DP> w_w;
  vector<NR::Vec_DP> w_ddw; //spline
};

class LuescherW80{
 public:
  LuescherW80(); //constructor sets up the splines

  double operator()(double qsq);

 private:
  vector<NR::Vec_DP> w_qsq; //vector over integer regions
  vector<NR::Vec_DP> w_w;
  vector<NR::Vec_DP> w_ddw; //spline
};



class LuescherPhiSigma{
 public:
  LuescherPhiSigma();

  double phi(string irrep, double q);
  double sigma(string irrep, double q);

 private:
  LuescherW00 w00;
  LuescherW40 w40;
  LuescherW60 w60;
  LuescherW80 w80;

  map<string,vector<double> > mll_weights;
  map<string,vector<double> > ml4_weights;
  map<string,vector<double> > m44_weights; //not currently used
};

class LuescherPhaseShift{
  public:
   LuescherPhaseShift(string irrep, EnsemReal mA, EnsemReal mB, double L); //spatial volume in temporal lattice units

   EnsemReal k(EnsemReal E);
   EnsemReal q(EnsemReal E);
   EnsemReal phi(EnsemReal E);
   EnsemReal sigma(EnsemReal E);

   EnsemReal delta(EnsemReal E);
   
   EnsemReal delta_corrected(EnsemReal E, EnsemReal tan_delta4); //jackknife tan_d4
   EnsemReal delta_corrected(EnsemReal E, double tan_delta4); //single value of tan_d4
   vector<EnsemReal> delta_corrected(EnsemReal E, double tan_delta4, double tan_delta4_error); //tan_d4 and error uncorrelated with ensemble

   EnsemReal k2Lp1_cot_delta(EnsemReal E){return pow(k(E), 2*ang_mom + 1) / tan(delta(E));};

   EnsemReal k2Lp1_cot_delta_corrected(EnsemReal E, EnsemReal tan_delta4){return pow(k(E), 2*ang_mom + 1) / tan(delta_corrected(E, tan_delta4)); };
   EnsemReal k2Lp1_cot_delta_corrected(EnsemReal E, double tan_delta4){return pow(k(E), 2*ang_mom + 1) / tan(delta_corrected(E, tan_delta4)); };
   
   vector<EnsemReal> k2Lp1_cot_delta_corrected(EnsemReal E, double tan_delta4, double tan_delta4_error){
	vector<EnsemReal> d = delta_corrected(E, tan_delta4, tan_delta4_error);
	vector<EnsemReal> dum;	
	for(int i=0; i < 3; i++){ dum.push_back( pow(k(E), 2*ang_mom + 1) / tan( d[i] ) ); };	
	return dum;
   }

  private:
   string irrep;
   EnsemReal mA, mB;
   double L; 	
   int ang_mom;
   LuescherPhiSigma phi_sigma;
};



//*******************************

/*class LuescherPhi{
 public:
  LuescherPhi(); //constructor makes the cubic spline interpolation

  double phi(double q); //returns the value of phi
  string type(double q); //returns how phi is calculated ("interpolate", "out of range", "arctan")

  double sigma(double q); //return the sensitivity

 private:
  NR::Vec_DP phi_data_qsq;
  NR::Vec_DP phi_data_phi;
  NR::Vec_DP phi_data_ddphi; //cubic spline 

  NR::Vec_DP sigma_data_qsq;
  NR::Vec_DP sigma_data_sigma;
  NR::Vec_DP sigma_data_ddsigma; //cubic spline

};

class LuescherPhaseShift{
 public:
  LuescherPhaseShift(double L_, double mA_, double mB_); // L = lattice length in TEMPORAL lattice units (e.g. 16*3.5)

  double delta(double E, int n); //energy and level number

 private:
  LuescherPhi luPhi;
  double L, mA, mB;  
};


class LuescherPhaseShiftEnsem{
 public:
  LuescherPhaseShiftEnsem(double L_, EnsemReal mA_, EnsemReal mB_);

  EnsemReal k(EnsemReal E);
  EnsemReal delta(EnsemReal E, int n);
  EnsemReal delta(EnsemReal E);
  EnsemReal k2Lp1cotdelta(EnsemReal E, int ang_mom);

  EnsemReal phi(EnsemReal E);
  EnsemReal sigma(EnsemReal E);		

 private:
  LuescherPhi luPhi;
  double L;
  EnsemReal mA, mB;
};
*/

#endif
