#include "fit_correlators.h"
#include "jackknife_fitter.h"
#include "fit_forms.h"
#include "semble_typedefs.h"

#include "ensem/ensem.h"
#include "adat/handle.h"

#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>



using namespace std;
using namespace ENSEM;
using namespace ADAT;

int main(int argc, char *argv[]){

  if(argc != 9){ cerr << "usage: <filename> <time extent of lattice> <tmin> <tmax> <noise_cutoff> <fitCrit> <min_tslices> <SVcut>" << 
endl; exit(1); }
  string filen; {istringstream val(argv[1]); val >> filen;}
  int T; {istringstream val(argv[2]); val >> T;}
  int tmin; {istringstream val(argv[3]); val >> tmin;}
  int tmax; {istringstream val(argv[4]); val >> tmax;}
  double cutoff; {istringstream val(argv[5]); val >> cutoff;}
  string fitCrit; {istringstream val(argv[6]); val >> fitCrit;}
  int minTSlices; {istringstream val(argv[7]); val >> minTSlices;}
  double svcut; {istringstream val(argv[8]); val >> svcut;}	

  if((fitCrit != "chisq_per_dof") && (fitCrit != "Q") && (fitCrit != "splitN") && (fitCrit != "generic") && (fitCrit != "QN"))
    { cerr << "fit criterion " << fitCrit << " is not known" << endl; 
      cerr << "known fitCrit values are : chisq_per_dof, Q, splitN, generic, QN" << endl;  
      exit(1);
    }
 
  //load the data
  EnsemVectorReal corr;  {ostringstream filename;   filename << filen;   read(filename.str(), corr);}
  cout << "read file: " << filen << endl;
  int bins = peekObs(corr, 0).size();  cout << "bins = " << bins << endl;

  //make the EnsemData object
  EnsemVectorReal lambda; lambda.resize(bins); lambda.resizeObs(tmax - tmin + 1);
  vector<double> tslices;
  for(int t = tmin; t <= tmax; t++){
    pokeObs(lambda, peekObs(corr, t), t - tmin); //correct use of peekObs ?
    tslices.push_back(t);
  }
  EnsemData corrData(tslices, lambda);

  corrData.setSVCutoff(svcut);

  //make the possible fit comparators
  map<string, Handle<FitComparator> > dum;
  
  CompareFitsByChisqPerNDoF* comp1 = new CompareFitsByChisqPerNDoF; Handle<FitComparator> comp1H(comp1);
  dum.insert( make_pair("chisq_per_dof", comp1H) );
  CompareFitsByQ* comp2 = new CompareFitsByQ; Handle<FitComparator> comp2H(comp2);
  dum.insert( make_pair("Q", comp2H) );
  CompareFitsBySplitN* comp3 = new CompareFitsBySplitN; Handle<FitComparator> comp3H(comp3);
  dum.insert( make_pair("splitN", comp3H) );
  CompareFitsByGeneric* comp4 = new CompareFitsByGeneric; Handle<FitComparator> comp4H(comp4);
  dum.insert( make_pair("generic", comp4H) );
  CompareFitsByQN* comp5 = new CompareFitsByQN; Handle<FitComparator> comp5H(comp5);
  dum.insert( make_pair("QN", comp5H) );
  

  //do the fit
  FitCorrelatorCosh fitCorr(corrData, T, dum[fitCrit], cutoff, minTSlices); 

  //output some business
  //fit summary 
  cout << endl << fitCorr.getFitSummary() << endl << endl;

  //the final result
  cout << "**********************************************************" << endl;
  cout << "   BEST FIT uses " << fitCorr.getNExp() << " coshes" << endl;
  cout << "   chisq/ndof = " << setprecision(2) << fitCorr.getChisq() << "/" << fitCorr.getNDoF() << " = " << double(fitCorr.getChisq() /  fitCorr.getNDoF() ) << endl;
  cout << "   " << fitCorr.getFitName() << endl;
  cout << "   m0 = " << setprecision(4) << toDouble(mean(fitCorr.getMass0())) << " +/- " << toDouble(sqrt(variance(fitCorr.getMass0()))) << endl;
  cout << "   a0 = " << toDouble(mean(fitCorr.getAmp0())) << " +/- " << toDouble(sqrt(variance(fitCorr.getAmp0()))) << endl;
  if( fitCorr.getNExp() == 2 ){
     cout << "   m1 = " << toDouble(mean(fitCorr.getMass1())) << " +/- " << toDouble(sqrt(variance(fitCorr.getMass1()))) << endl;
     cout << "   a1 = " << toDouble(mean(fitCorr.getAmp1())) << " +/- " << toDouble(sqrt(variance(fitCorr.getAmp1()))) << endl;
  }
  cout << "**********************************************************" << endl;

  //assess importance of cosh versus exp
  double m = toDouble( mean( fitCorr.getMass0() ) );
  int tt = min( tmax, T/2);
  double cosht = 2.0 * exp( - m * T / 2 ) * cosh( m * (tt - T/2) );
  double expt = exp( - m * tt );
  double diff = (cosht - expt) / cosht;
  cout << " at t= " << tt << ", cosh and exp differ at the " << 100*diff << "% level" << endl; 

  //write the log to a file
  {
    stringstream s; s << filen << "_cosh_fit.log"; 
    ofstream out; out.open(s.str().c_str());
    out << fitCorr.getFitSummary();
    out.close();
  }
  //write the plot to a file
  {  
    stringstream ss; ss << filen << "_cosh_fit.ax"; 
    ofstream out; out.open(ss.str().c_str());
    out << fitCorr.getFitPlotString();
    out.close();
  }

  //write the jackknife fit files
  EnsemReal out = fitCorr.getMass0();
  {  ostringstream filename;  filename << filen << "_cosh_fit_mass0.jack";  write(filename.str(), out );  } //cout << "wrote " << filename.str() << endl; }
  out = fitCorr.getAmp0();
  {  ostringstream filename;  filename << filen << "_cosh_fit_amp0.jack";  write(filename.str(), out );   }//cout << "wrote " << filename.str() << endl; }

  if(fitCorr.getNExp() == 2){
    out = fitCorr.getMass1();
    {  ostringstream filename;  filename << filen << "_cosh_fit_mass1.jack";  write(filename.str(), out );  }//cout << "wrote " << filename.str() << endl; }
    out = fitCorr.getAmp1();
    {  ostringstream filename;  filename << filen << "_cosh_fit_amp1.jack";  write(filename.str(), out );  }//cout << "wrote " << filename.str() << endl; }
  }



}
