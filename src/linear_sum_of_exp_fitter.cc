#include "jackknife_linear_fitter.h"
#include "ensem_data.h"
#include "plot.h"

#include <itpp/itbase.h>
#include "ensem/ensem.h"
#include "adat/handle.h"

#include <vector>
#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;
using namespace ENSEM;

int main(int argc, char *argv[]){
  if(argc != 6){
    cerr << "Usage:" << argv[0] << ": <data_filename> <nExp> <tmin> <tmax> <sv_cutoff>" << endl;
    exit(1);
  }

  string datafilename; {istringstream val(argv[1]); val >> datafilename;}
  //  int tSrc; {istringstream val(argv[2]); val >> tSrc;}
  int nExp; {istringstream val(argv[2]); val >> nExp;}
  int tmin; {istringstream val(argv[3]); val >> tmin;}
  int tmax; {istringstream val(argv[4]); val >> tmax;}
  double sv_cutoff; {istringstream val(argv[5]); val >> sv_cutoff;}
  //  string argsfilename; {istringstream val(argv[6]); val >> argsfilename;}

  //load the data
  EnsemVectorReal corr;
  {ostringstream filename;  filename << datafilename;  read(filename.str(), corr);}
  //  cout << "read file: " << filename << endl;

  //cout << "loaded correlator" << endl;

  vector<double> tvals;
  for(int t = 0; t < corr.numElem(); t++){
    tvals.push_back(double(t));}

  EnsemData data(tvals, corr);
  data.hideDataBelowX( tmin - 0.1 );
  data.hideDataAboveX( tmax + 0.1 );
  data.makeCov();
  data.makeInvCov();

  //  cout << "made covariance and inverted" << endl;


  //read the weights and "masses"
  //constant weight 
  EnsemReal W0; 
  {ostringstream filename;  filename << "W_0";  read(filename.str(), W0);}

  vector<EnsemReal> W; vector<EnsemReal> E;
  for(int i = 1; i <= nExp; i++){
    stringstream dum; dum << "W_" << i;
    EnsemReal w;  {ostringstream filename;  filename << dum.str();  read(filename.str(), w);}
    stringstream dum2; dum2 << "E_" << i;
    EnsemReal e;  {ostringstream filename;  filename << dum2.str();  read(filename.str(), e);}
    W.push_back(w); E.push_back(- e);
  }

  // cout << "loaded masses and weights" << endl;

  //make the functions
  vector< Handle<EnsemFunction> > basis;
  ConstTimesExp* dum = new ConstTimesExp(W0, Real(0.0)*W0); Handle<EnsemFunction> constant(dum); 
  basis.push_back(constant);
  //  cout << "made the constant function" << endl;

  for(int i = 1; i <= nExp; i++){
    ConstTimesExp* dum1 = new ConstTimesExp(W[i-1], E[i-1]); Handle<EnsemFunction> exp(dum1);
    basis.push_back(exp);
    //   cout << "added an exponential function" << endl;
  }

  JackLinearFit lin_fit(data, basis, sv_cutoff);
  
  cout << "** linear fit complete **" << endl;
  cout << "chisq/nDoF = " << lin_fit.getChisq() << "/" << lin_fit.getNDoF() << endl;
  cout << "Q = " << lin_fit.getQ() << endl << endl;

  vector<EnsemReal> fit_pars = lin_fit.getFitParValues();
  cout << "const = " << mean(fit_pars[0]) << " +/- " << sqrt(variance(fit_pars[0])) << endl;
  for(int i = 1; i <= nExp; i++){
    cout << "a_" << i << " = " << mean(fit_pars[i]) << " +/- " << sqrt(variance(fit_pars[i])) << endl;
  }

  cout << endl << "parameter correlation:" << endl;
  itpp::mat cov_par = lin_fit.getFitParCov();
  itpp::mat d(nExp +1, nExp +1);
  for(int i=0; i < nExp+1; i++){d(i,i) = 1.0 / sqrt(cov_par(i,i));}
  cout << d*cov_par*d << endl;

  // write out ensem files

  // make a plot
  AxisPlot plot;
  plot.addEnsemData(data);
  
  vector<int> tslices;
  vector<double> avg, avgplus, avgminus;
  for(int t=0; t <= tmin; t++){
    tslices.push_back(t);
    EnsemReal dum = W0 * fit_pars[0];
    for(int n = 1; n <= nExp; n++){ dum += fit_pars[n] * basis[n]->operator()(t); }
    double av = toDouble(mean(dum));
    double err = toDouble(sqrt(variance(dum)));
    avg.push_back(av); avgplus.push_back(av+err); avgminus.push_back(av-err);
  }
  plot.addLineData(tslices, avg, 3);
  plot.addLineData(tslices, avgplus, 3);
  plot.addLineData(tslices, avgminus, 3);

  tslices.clear(); avg.clear(); avgplus.clear(); avgminus.clear();
  for(int t=tmin; t <= tmax; t++){
    tslices.push_back(t);
    EnsemReal dum = W0 * fit_pars[0];
    for(int n = 1; n <= nExp; n++){ dum += fit_pars[n] * basis[n]->operator()(t); }
    double av = toDouble(mean(dum));
    double err = toDouble(sqrt(variance(dum)));
    avg.push_back(av); avgplus.push_back(av+err); avgminus.push_back(av-err);
  }
  plot.addLineData(tslices, avg, 2);
  plot.addLineData(tslices, avgplus, 2);
  plot.addLineData(tslices, avgminus, 2);
  
  tslices.clear(); avg.clear(); avgplus.clear(); avgminus.clear();
  for(int t=0; t <= int(tvals.back()); t++){
    tslices.push_back(t);
    EnsemReal dum = W0 * fit_pars[0];
    double av = toDouble(mean(dum));
    double err = toDouble(sqrt(variance(dum)));
    avg.push_back(av); avgplus.push_back(av+err); avgminus.push_back(av-err);
  }
  plot.addLineData(tslices, avg, 5);
  plot.addLineData(tslices, avgplus, 5);
  plot.addLineData(tslices, avgminus, 5);


  plot.sendToFile("plot.ax");

}
