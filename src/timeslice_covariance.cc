#include "ensem/ensem.h"
#include "adat/handle.h"
#include "itpp/itbase.h"
#include "ensem_data.h"

#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>



using namespace std;
using namespace ENSEM;
using namespace ADAT;

int main(int argc, char *argv[]){

  if(argc != 4){ cerr << "usage: <filename> <tmin> <tmax>" << endl; exit(1); }
  string filen; {istringstream val(argv[1]); val >> filen;}
  int tmin; {istringstream val(argv[2]); val >> tmin;}
  int tmax; {istringstream val(argv[3]); val >> tmax;}
 
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

  itpp::mat cov = corrData.getCov();
  itpp::vec v = 1.0 / itpp::sqrt( itpp::diag(cov) );
  itpp::mat norm = itpp::diag(v);
  itpp::mat norm_cov = norm*cov*norm;

  cout << "normalised data covariance = " << endl << norm_cov << endl;

  itpp::mat U = cov; itpp::mat V = cov; itpp::vec s; s.set_size(cov.rows());
  itpp::svd(norm_cov, U, s, V);
  cout << "singular values : " << endl << s << endl;
  cout << "condition number = " << s(0) / (s.right(1))(0) << endl;

  cout << "unnormalised covariance = " << endl << cov << endl;

}
