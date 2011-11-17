#include "luescher.h"
#include "ensem_data.h"

#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

int main(int argc, char** argv){
  
  if(argc < 6){ cerr << " <L/a_t> <m_A filename> <m_B filename> <irrep> <E_1 filename> [<E_2 filename> ....]" << endl; exit(1);}

  double L; {istringstream a(argv[1]); a >> L;};
  string irrep; {istringstream a(argv[4]); a >> irrep;};

  //load files
  EnsemReal mA;
  {ostringstream fname; fname << argv[2]; read(fname.str(), mA);}
  EnsemReal mB;
  {ostringstream fname; fname << argv[3]; read(fname.str(), mB);}
  vector<EnsemReal> E;
  for(int n = 5; n < argc; n++){
    EnsemReal dum;
    {ostringstream fname; fname << argv[n]; read(fname.str(), dum);}
    E.push_back(dum);
  }

  //check ensems are compatible
  int nBins = mA.size();
  if(mB.size() != nBins){cerr << "ensembles don't match" << endl; exit(1);}
  for(int n=0; n < E.size(); n++){
    if(E[n].size() != nBins){cerr << "ensembles don't match" << endl; exit(1);}
  }

  //compute phase shifts etc...
  LuescherPhaseShift luDelta(irrep, mA, mB, L);

  vector<EnsemReal> k;
  vector<EnsemReal> ksq;
  vector<EnsemReal> q;
  vector<EnsemReal> qsq;
  vector<EnsemReal> phi;
  vector<EnsemReal> sigma;
  vector<EnsemReal> delta;
  vector<EnsemReal> k2Lp1cotdelta;

  EnsemData eff_range_fit_data; //nb can't deal with ensem values of x 

  for(int n=0; n < E.size(); n++){
    EnsemReal e = E[n];
    EnsemReal K = luDelta.k(e); k.push_back(K); ksq.push_back(K*K);
    EnsemReal Q = luDelta.q(e); q.push_back(Q); qsq.push_back(Q*Q);	

    EnsemReal Phi = luDelta.phi(e); phi.push_back(Phi);
    EnsemReal Sigma = luDelta.sigma(e); sigma.push_back(Sigma);

    EnsemReal Delta = luDelta.delta(e); delta.push_back(Delta);
    EnsemReal Kcot = luDelta.k2Lp1_cot_delta(e); k2Lp1cotdelta.push_back(Kcot);

    eff_range_fit_data.addDatum(toDouble(mean(K*K)), Kcot);
  }

  //write out results
  //need a header line
  {char buffer[1000];
    sprintf(buffer, "        E*a_t \t\t         (k*a_t)^2  \t     q^2 \t     phi \t     sigma \t         delta/deg \t   (k*a_t)^(2L+1) cot(delta) ");
    cout << buffer << endl;
  }

  for(int n=0; n < E.size(); n++){
    double e = toDouble(mean(E[n]));
    double e_err = toDouble(sqrt(variance(E[n])));

    double K2 = toDouble(mean(ksq[n]));
    double K2_err = toDouble(sqrt(variance(ksq[n])));

    double Q2 = toDouble(mean(qsq[n]));
    double Q2_err = toDouble(sqrt(variance(qsq[n])));

    double Phi = toDouble(mean(phi[n]));
    double Phi_err = toDouble(sqrt(variance(phi[n])));	

    double Sigma = toDouble(mean(sigma[n]));
    double Sigma_err = toDouble(sqrt(variance(sigma[n])));

    double pi = 2.0*atan2(1.0,0.0);
    double angle = (180/pi)*toDouble(mean(delta[n]));
    double angle_err = (180/pi)*toDouble(sqrt(variance(delta[n])));

    double eff_range = toDouble(mean(k2Lp1cotdelta[n]));
    double eff_range_err = toDouble(sqrt(variance(k2Lp1cotdelta[n])));


    //nicely formatted output
    char buffer[1000];
    sprintf(buffer, "%8.5g +/- %8.3g\t%8.5g +/- %8.3g\t%5.3g +/- %5.3g\t%5.3g +/- %5.3g  %5.3g +/- %5.3g\t%5.3g +/- %5.3g\t   %5.3g +/ %5.3g", e, e_err, 
K2, K2_err, Q2, 
Q2_err, Phi, Phi_err, Sigma, Sigma_err, angle, angle_err, eff_range, eff_range_err);
    cout << buffer << endl;
 
  }

  //output data covariance
  itpp::mat cov = eff_range_fit_data.getCov();
  cout << cov << endl;
};

