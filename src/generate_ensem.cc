#include <vector>
#include <algorithm>
#include "ensem/ensem.h"
#include <cassert>
#include <itpp/itbase.h>
#include <itpp/base/random.h>
#include <iostream>

using namespace ENSEM;
using namespace std;

int main(int argc, char** argv){
  
  if( argc < 5 ){ cerr << "generate_ensem <mean> <error> <Ncfgs> <filename> \n "; exit(1); }

  double mean; {istringstream a(argv[1]); a >> mean;};
  double error; {istringstream a(argv[2]); a >> error;};
  int N; {istringstream a(argv[3]); a >> N;};
  string filename; {istringstream a(argv[4]); a >> filename;};

  itpp::RNG_randomize();

  itpp::Normal_RNG rng(mean, (N-1)*pow(error,2) ); //check defn of error on the mean versus variance

  EnsemReal ensem; ensem.resize(N);

  for(int n = 0; n < N; n++){
    pokeEnsem(ensem, Real(rng()), n);
  }
  
  write(filename, ensem );  

  exit(0);
};
