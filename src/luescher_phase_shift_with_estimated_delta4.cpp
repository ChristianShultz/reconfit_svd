#include "luescher.h"
#include "ensem_data.h"

#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

int main(int argc, char** argv){
  
  if(argc < 8){ cerr << " <L/a_t> <m_A filename> <m_B filename> <irrep> <est a_4> <est a_4 error> <E_1 filename> [<E_2 filename> ....]" << endl; exit(1);}

  double L; {istringstream a(argv[1]); a >> L;};
  string irrep; {istringstream a(argv[4]); a >> irrep;};
  double a4; {istringstream a(argv[5]); a >> a4;};
  double a4_err; {istringstream a(argv[6]); a >> a4_err;}; 


  //load files
  EnsemReal mA;
  {ostringstream fname; fname << argv[2]; read(fname.str(), mA);}
  EnsemReal mB;
  {ostringstream fname; fname << argv[3]; read(fname.str(), mB);}

  vector<EnsemReal> E;
  for(int n = 7; n < argc; n++){
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

  vector< vector<EnsemReal> > delta_corr;
  vector< vector<EnsemReal> > kcot_corr;

  EnsemData eff_range_fit_data; //nb can't deal with ensem values of x 

  for(int n=0; n < E.size(); n++){
    EnsemReal e = E[n];
    EnsemReal K = luDelta.k(e); k.push_back(K); ksq.push_back(K*K);
    EnsemReal Q = luDelta.q(e); q.push_back(Q); qsq.push_back(Q*Q);	

    EnsemReal Phi = luDelta.phi(e); phi.push_back(Phi);
    EnsemReal Sigma = luDelta.sigma(e); sigma.push_back(Sigma);

    EnsemReal Delta = luDelta.delta(e); delta.push_back(Delta);
    EnsemReal Kcot = luDelta.k2Lp1_cot_delta(e); k2Lp1cotdelta.push_back(Kcot);

    double t4 = a4 * pow( toDouble(mean(K)) , 9 );
    double t4_err = a4_err *  pow( toDouble(mean(K)) , 9 );	

    vector<EnsemReal> DeltaCorr = luDelta.delta_corrected(e, t4, t4_err); delta_corr.push_back(DeltaCorr);
    vector<EnsemReal> KcotCorr = luDelta.k2Lp1_cot_delta_corrected(e, t4, t4_err);kcot_corr.push_back(KcotCorr);	

    eff_range_fit_data.addDatum(toDouble(mean(K*K)), KcotCorr[1]);
  }

  //write out results
  //need a header line
/*  {char buffer[1000];
    sprintf(buffer, "        E*a_t \t\t         (k*a_t)^2  \t     q^2 \t     phi \t     sigma \t         delta/deg \t   (k*a_t)^(2L+1) cot(delta)\t delta corr \t kcot corr ");
    cout << buffer << endl;
  }*/

  for(int n=0; n < E.size(); n++){
    double e = toDouble(mean(E[n]));
    double e_err = toDouble(sqrt(variance(E[n])));
    cout << "at*E = " << e << " +/- " << e_err<< endl;	

    double K2 = toDouble(mean(ksq[n]));
    double K2_err = toDouble(sqrt(variance(ksq[n])));
    cout << "(at*k)^2 = " << K2 << " +/- " << K2_err << endl;

    double Q2 = toDouble(mean(qsq[n]));
    double Q2_err = toDouble(sqrt(variance(qsq[n])));
    cout << "q^2 = " << Q2 << " +/- " << Q2_err << endl;

    double Phi = toDouble(mean(phi[n]));
    double Phi_err = toDouble(sqrt(variance(phi[n])));	
    cout << "phi = " << Phi << " +/- " << Phi_err << endl;

    double Sigma = toDouble(mean(sigma[n]));
    double Sigma_err = toDouble(sqrt(variance(sigma[n])));
    cout << "sigma = " << Sigma << " +/- " << Sigma_err << endl;

    double pi = 2.0*atan2(1.0,0.0);
    double angle = (180/pi)*toDouble(mean(delta[n]));
    double angle_err = (180/pi)*toDouble(sqrt(variance(delta[n])));
    cout << "delta(uncorr) = " << angle << " +/- " << angle_err << endl;	

    double eff_range = toDouble(mean(k2Lp1cotdelta[n]));
    double eff_range_err = toDouble(sqrt(variance(k2Lp1cotdelta[n])));
    cout << "k^(2L+1) cot(delta(uncorr)) = " << eff_range << " +/- " << eff_range_err << endl; 

    vector<EnsemReal> delta = delta_corr[n];
    cout << "corrected phase shift : " << endl;
    cout << (180/pi)*toDouble(mean(delta[0])) << " +/- " << (180/pi)*toDouble(sqrt(variance(delta[0]))) << endl;
    cout << (180/pi)*toDouble(mean(delta[1])) << " +/- " << (180/pi)*toDouble(sqrt(variance(delta[1])))	<< endl;
    cout << (180/pi)*toDouble(mean(delta[2])) << " +/- " << (180/pi)*toDouble(sqrt(variance(delta[2])))	<< endl;	   

    vector<EnsemReal> kcot = kcot_corr[n];
    cout << "corrected kcot : " << endl;
    cout << toDouble(mean(kcot[0])) << " +/- " << toDouble(sqrt(variance(kcot[0])))	<< endl;
    cout << toDouble(mean(kcot[1])) << " +/- " << toDouble(sqrt(variance(kcot[1])))     << endl;
    cout << toDouble(mean(kcot[2])) << " +/- " << toDouble(sqrt(variance(kcot[2])))     << endl;

    cout << endl << endl;
    //nicely formatted output
/*    char buffer[1000];
    sprintf(buffer, "%8.5g +/- %8.3g\t%8.5g +/- %8.3g\t%5.3g +/- %5.3g\t%5.3g +/- %5.3g  %5.3g +/- %5.3g\t%5.3g +/- %5.3g\t   %5.3g +/ %5.3g", e, e_err, 
K2, K2_err, Q2, 
Q2_err, Phi, Phi_err, Sigma, Sigma_err, angle, angle_err, eff_range, eff_range_err);
    cout << buffer << endl; */
 
  }

  //output data covariance
  itpp::mat cov = eff_range_fit_data.getCov();
  cout << cov << endl;
};

