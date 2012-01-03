#include "multi_princorr_fit.h"

#include "ensem_data.h"
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

int main(int argc, char** argv){

  if(argc < 2){cerr << " arguments are <tLen> <noiseCut> <filename1> <t0_1> <tmin1> [<filename2> <t0_2> <tmin2> ...]" << endl; exit(1);}

  vector<int> tmin;
  vector<int> t0;
  vector<EnsemVectorReal> data_list;

  int tLen; {istringstream a(argv[1]); a >> tLen;};
  double noiseCut; {istringstream a(argv[2]); a >> noiseCut;};

  //load the files
  for(int n = 3; n < argc; n +=3){
    EnsemVectorReal corr;
    {ostringstream filename;  filename << argv[n];  read(filename.str(), corr);}
    data_list.push_back(corr);

    int aa; {istringstream a(argv[n+2]); a >> aa;};
    int bb; {istringstream b(argv[n+1]); b >> bb;};

    tmin.push_back(aa);
    t0.push_back(bb);

    if(corr.numElem() >= tLen){ cerr << "length of a correlator is bigger than tLen" << endl; exit(1);}
  }


  //concatenate the data, padding out to tLen (specified in command line args), place into an EnsemData object
  int nCorrs = tmin.size();

  vector<EnsemData> x;

  int nBins = data_list[0].size();

  EnsemReal neg; neg.resize(nBins); neg = Real(-1.0);

  for(int n = 0; n < nCorrs; n++){
    vector<double> tvals;
    for(int t = 0; t < (data_list[n]).numElem(); t++){tvals.push_back(double(t + n*tLen));}
    EnsemData data_dum(tvals, data_list[n]);
    for(int t = (data_list[n]).numElem(); t < tLen; t++){ data_dum.addDatum(double(t + n*tLen), neg);}
    x.push_back(data_dum);
  }

  EnsemData data = x[0];
  for(int n=1; n < nCorrs; n++){
    EnsemData tmp = data;
    data = concat( tmp, x[n]);
  }
   
  //run the fit
  FitMultiplePrincipalCorrelator fit(data, nCorrs, t0, tLen, tmin, noiseCut);
  fit.saveFitPlot("plot.ax");
  
  //write out the results
  //add write-out routines here if you like
  EnsemReal mult = fit.getMass0();
  write("multi.jack", mult);


};
