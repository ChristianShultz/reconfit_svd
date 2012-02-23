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

  if(argc != 8){ cerr << "usage: <filename> <t0> <tmin> <tmax> <fitCrit> <min_tslices> <covar_SV_cutoff>" << endl; exit(1); }
  string filen; {istringstream val(argv[1]); val >> filen;}
  int t0; {istringstream val(argv[2]); val >> t0;}
  int tmin; {istringstream val(argv[3]); val >> tmin;}
  int tmax; {istringstream val(argv[4]); val >> tmax;}
  string fitCrit; {istringstream val(argv[5]); val >> fitCrit;}
  int minTSlices; {istringstream val(argv[6]); val >> minTSlices;}
  double SVcutoff;  {istringstream val(argv[7]); val >> SVcutoff;}

  if((fitCrit != "chisq_per_dof") && (fitCrit != "Q") && (fitCrit != "Zfit") && (fitCrit != "QN"))
    { cerr << "fit criterion " << fitCrit << " is not known" << endl; 
      cerr << "known fitCrit values are : chisq_per_dof, Q, fitZ, QN" << endl;  
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
    pokeObs(lambda, peekObs(corr, t), t - tmin);
    tslices.push_back(t);
  }
  EnsemData corrData(tslices, lambda);

  corrData.setSVCutoff(SVcutoff);
  cout << "covariance will be inverted cutting off SV below " << corrData.getSVCutoff() << endl; 

  //cout << corrData << endl; //(loading routine does work)

  //make the possible fit comparators
  map<string, Handle<FitComparator> > dum;
  
  CompareFitsByChisqPerNDoF* comp1 = new CompareFitsByChisqPerNDoF; Handle<FitComparator> comp1H(comp1);
  dum.insert( make_pair("chisq_per_dof", comp1H) );
  CompareFitsByQ* comp2 = new CompareFitsByQ; Handle<FitComparator> comp2H(comp2);
  dum.insert( make_pair("Q", comp2H) );
  CompareZFits* comp3 = new CompareZFits; Handle<FitComparator> comp3H(comp3);
  dum.insert( make_pair("Zfit", comp3H) );
  CompareFitsByQN* comp5 = new CompareFitsByQN; Handle<FitComparator> comp5H(comp5);
  dum.insert( make_pair("QN", comp5H) );
  
  //do the fit
  FitZ fitCorr(corrData, t0, dum[fitCrit], minTSlices); 

  //output some business
  //fit summary 
  cout << endl << fitCorr.getFitSummary() << endl << endl;
  
  //the final result
  cout << "**********************************************************" << endl;
  cout << "   chisq/ndof = " << setprecision(2) << fitCorr.getChisq() << "/" << fitCorr.getNDoF() << " = " << double(fitCorr.getChisq() /  fitCorr.getNDoF() ) << endl;
  cout << "   " << fitCorr.getFitName() << endl;
  cout << "   c = " << setprecision(4) << toDouble(mean(fitCorr.getZ())) << " +/- " << toDouble(sqrt(variance(fitCorr.getZ()))) << endl;
  if( fitCorr.getFitType() != "constant" ){
     cout << "   E = " << toDouble(mean(fitCorr.getExpMass())) << " +/- " << toDouble(sqrt(variance(fitCorr.getExpMass()))) << endl;
     cout << "   Z' = " << toDouble(mean(fitCorr.getZ2())) << " +/- " << toDouble(sqrt(variance(fitCorr.getZ2()))) << endl;
  }
  cout << "**********************************************************" << endl;

   //write the log to a file
  {
    stringstream s; s << filen << "_const_and_exp_fit.log"; 
    ofstream out; out.open(s.str().c_str());
    out << fitCorr.getFitSummary();
    out.close();
  }
  //write the plot to a file
  {  
    stringstream ss; ss << filen << "_const_and_exp_fit.ax"; 
    ofstream out; out.open(ss.str().c_str());
    out << fitCorr.getFitPlotString();
    out.close();
  }

  //write the jackknife fit files
  EnsemReal out = fitCorr.getZ();
  {  ostringstream filename;  filename << filen << "_const_and_exp_fit_Z.jack";  write(filename.str(), out );  } //cout << "wrote " << filename.str() << endl; }
  if(fitCorr.getFitType() != "constant"){
    out = fitCorr.getExpMass();
    {  ostringstream filename;  filename << filen << "_const_and_exp_fit_E.jack";  write(filename.str(), out );  }//cout << "wrote " << filename.str() << endl; }
    out = fitCorr.getZ2();
    {  ostringstream filename;  filename << filen << "_const_and_exp_fit_Z2.jack";  write(filename.str(), out );  }//cout << "wrote " << filename.str() << endl; }
  }






}
