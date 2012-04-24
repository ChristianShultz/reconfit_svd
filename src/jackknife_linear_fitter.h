#ifndef __JACKKNIFE_LINEAR_FITTER_H__
#define __JACKKNIFE_LINEAR_FITTER_H__

#include <vector>
#include <algorithm>
#include <sstream>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>

#include <itpp/itbase.h>
#include "ensem/ensem.h"
#include "adat/handle.h"

#include "ensem_data.h"


using namespace ENSEM;
using namespace ADAT;
using namespace std;

class JackLinearFit
{
 public:
  JackLinearFit(EnsemData data, vector< Handle<EnsemFunction> > basis_f, double sv_cutoff);
  //fit is done in the constructor

  //fit description
  double getChisq(){return chisq;}
  int getNDoF(){return nDoF;}
  double getChisqPerNDoF(){return chisq / nDoF;}
  double getQ(){return Q;}

  //fit results
  vector<EnsemReal> getFitParValues(){return par_vec;}
  itpp::mat getFitParCov(){return cov_par;}


 private:
  double chisq;
  int nDoF;
  double Q;
  vector<EnsemReal> par_vec;
  itpp::mat cov_par;

};



#endif
