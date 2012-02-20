#include "ensem_data.h"

//some search functions
bool double_predicate(double a, double b){
  double eps = 1.0e-15; //hard-coded comparison level
  return (abs(a - b) < eps);
}

int find_nearest(vector<double> v, double a){
  double lowest_diff = abs(v[0] - a);
  int pos = 0;
  for(int i=1; i < v.size(); i++){
    double temp =  abs(v[i] - a) ;
    if(temp < lowest_diff ){ lowest_diff = temp; pos = i;}
  }
  return pos;
}

int find_exact(vector<double> v, double a){
  vector<double>::iterator iter;
  double find[1]; find[0] = a;
  iter = search(v.begin(), v.end(), find, find + 1, double_predicate);    
  int pos =  iter - v.begin() ;
  if(iter == v.end()){pos = -1;}
  return pos;
}

vector<int> find_below(vector<double> v, double a){
  vector<int> below;
  for(int i=0; i < v.size(); i++){
    if( v[i] < a ){below.push_back(i);}
  }
  return below;
}

vector<int> find_above(vector<double> v, double a){
  vector<int> above;
  for(int i=0; i < v.size(); i++){
    if( v[i] > a ){above.push_back(i);}
  }
  return above;
}

vector<int> find_in_range(vector<double> v, double low, double high){
  vector<int> inrange;
  for(int i=0; i < v.size(); i++){
    if( (v[i] <= high) && (v[i] >= low) ){inrange.push_back(i);}
  }
  return inrange;
}

//****************************************

itpp::mat makeEnsemCov(EnsemVectorReal y){
  int size = y.numElem(); int bins = y.size();
  itpp::mat cov;
  cov.set_size(size, size);

  for(int i=0; i < size; i++){
    Real mi = mean( peekObs(y , i) );
    EnsemReal ei = peekObs(y , i);

    for(int j=i; j < size; j++){
      Real mj = mean( peekObs(y , j) );
      EnsemReal ej = peekObs(y , j);

      double covar = 0.;
      for(int bin = 0; bin < bins; bin++ ){
	covar += toDouble( (peekEnsem(ei , bin) - mi ) * (peekEnsem(ej , bin) - mj ) );
      }
      cov(i,j) = covar;
      cov(j,i) = covar;
    }
  }
  cov /= bins * (bins - 1. );
  return cov;
}



//****************************************
// ENSEMDATA CLASS
//****************************************
EnsemData::EnsemData(vector<double> x, EnsemVectorReal y){
  if(x.size() != y.numElem()){cerr << "sizes don't match" << endl; exit(1);}

  totalNData = x.size();  nBins = y.size(); SVCutoff = 1.0e-6;
  initP = true; initCov = false; initInvCov = false;

  x_data = x;
  y_data = y;
  for(int i=0; i < totalNData; i++){ active_data.push_back(true);} //all data active
};

/*EnsemData::EnsemData(const EnsemData& e){
  totalNData = e.getTotalNData(); nBins = e.getNBins(); SVCutoff = e.getSVCutoff();
  x_data = e.getAllXData();
  y_data = e.getAllYData();
  active_data = e.getActiveData();
  initP = true; initCov = e.isCovMade();
  if(initCov){cov = e.getAllCov();} 
};

EnsemData& EnsemData::operator=(const EnsemData& e){
  if(this != &e){
    totalNData = e.getTotalNData(); nBins = e.getNBins(); SVCutoff = e.getSVCutoff();
    x_data = e.getAllXData();
    y_data = e.getAllYData();
    active_data = e.getActiveData();
    initP = true; initCov = e.isCovMade();
    if(initCov){cov = e.getAllCov();} 
  }
  return *this;
  };
*/

//*****************
void EnsemData::addDatum(double x, EnsemReal y){
  if(initP){
    if(y.size() != nBins){cerr << "wrong binsize" << endl; exit(1);}  
    //check there isn't already a datum at this x-values
    int find = find_exact(x_data, x);
    
    if(find = -1){
      EnsemVectorReal temp = y_data; y_data.resizeObs(totalNData + 1);
      for(int i =0; i < totalNData; i++){pokeObs(y_data, peekObs(temp , i), i);}    
      pokeObs(y_data, y, totalNData);
      x_data.push_back(x);
      active_data.push_back(true);
      if(initCov){augmentCov(y);}
      totalNData++;
    }
    else{ 
      cerr << "WARNING: already a data point here, overwriting" << endl;
      x_data[find] = x; pokeObs(y_data, y, find);
      if(initCov){makeCov();} //remake the covariance with the new point - slow
    }

  }
  else{ //this is the first data-point    
    nBins = y.size();
    y_data.resize(nBins); y_data.resizeObs(1);
    pokeObs(y_data, y, 0);
    initP = true;
    x_data.clear(); active_data.clear();
    x_data.push_back(x);
    active_data.push_back(true);
    totalNData = 1; 
    if(initCov){makeCov();}
  }
  initInvCov = false; //inverse covariance no longer appropriate
};


void EnsemData::hideDatumByIndex(int i){
  if( (i >= totalNData)||(i < 0) ){cerr << "index outside range" << endl; exit(1);} 
  active_data[i] = false;
  initInvCov = false;
};

void EnsemData::hideDataByIndex(vector<int> list){
  for(int j =0; j < list.size(); j++){
    int i = list[j];
    hideDatumByIndex(i);
  }
  initInvCov = false;
};

bool EnsemData::hideDatumByX(double x){
  bool temp = false;
  int pos = find_exact(x_data, x);
  if(pos != -1){ active_data[pos] = false; temp = true; initInvCov = false;}
  return temp;
};

int EnsemData::hideDataByXRange(double xlow, double xhigh){
  vector<int> found = find_in_range(x_data, xlow, xhigh);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};

int EnsemData::hideDataBelowX(double x){
  vector<int> found = find_below(x_data, x);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};

int EnsemData::hideDataAboveX(double x){
  vector<int> found = find_above(x_data, x);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};

bool EnsemData::showDatumByX(double x){
  bool temp = false;
  int pos = find_exact(x_data, x);
  if(pos != -1){ active_data[pos] = true; temp = true;  initInvCov = false;}
  return temp;
};

int EnsemData::showDataByXRange(double xlow, double xhigh){
  vector<int> found = find_in_range(x_data, xlow, xhigh);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = true;  initInvCov = false;
  }
  return found.size();
};

int EnsemData::showDataBelowX(double x){
  vector<int> found = find_below(x_data, x);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = true;  initInvCov = false;
  }
  return found.size();
};

int EnsemData::showDataAboveX(double x){
  vector<int> found = find_above(x_data, x);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = true;  initInvCov = false;
  }
  return found.size();
};

int EnsemData::hideDataByYRange(double ylow, double yhigh){
  vector<double> y_mean;
  for(int i=0; i < totalNData; i++){
    y_mean.push_back(toDouble(mean( peekObs(y_data , i) )));
  } 
  vector<int> found = find_in_range(y_mean, ylow, yhigh);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};

int EnsemData::hideDataBelowY(double y){
  vector<double> y_mean;
  for(int i=0; i < totalNData; i++){
    y_mean.push_back(toDouble(mean( peekObs(y_data , i) )));
  } 
  vector<int> found = find_below(y_mean, y);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};


int EnsemData::hideDataAboveY(double y){
  vector<double> y_mean;
  for(int i=0; i < totalNData; i++){
    y_mean.push_back(toDouble(mean( peekObs(y_data , i) )));
  } 
  vector<int> found = find_above(y_mean, y);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};


int EnsemData::hideDataAboveYErr(double err){
  vector<double> y_err;
  for(int i=0; i < totalNData; i++){
    y_err.push_back(toDouble(sqrt(variance( peekObs(y_data,i) ))));
  } 
  vector<int> found = find_above(y_err, err);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};


int EnsemData::showDataByYRange(double ylow, double yhigh){
  vector<double> y_mean;
  for(int i=0; i < totalNData; i++){
    y_mean.push_back(toDouble(mean( peekObs(y_data , i) )));
  } 
  vector<int> found = find_in_range(y_mean, ylow, yhigh);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = true;  initInvCov = false;
  }
  return found.size();
};


int EnsemData::showDataBelowY(double y){
  vector<double> y_mean;
  for(int i=0; i < totalNData; i++){
    y_mean.push_back(toDouble(mean( peekObs(y_data , i) )));
  } 
  vector<int> found = find_below(y_mean, y);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = true; initInvCov = false;
  }
  return found.size();
};


int EnsemData::showDataAboveY(double y){
  vector<double> y_mean;
  for(int i=0; i < totalNData; i++){
    y_mean.push_back(toDouble(mean( peekObs(y_data , i) )));
  } 
  vector<int> found = find_above(y_mean, y);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = true;  initInvCov = false;
  }
  return found.size();
};


int EnsemData::showDataAboveYErr(double err){
  vector<double> y_err;
  for(int i=0; i < totalNData; i++){
    y_err.push_back(toDouble(sqrt(variance(peekObs(y_data,i)))));
  } 
  vector<int> found = find_above(y_err, err);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = true;  initInvCov = false;
  }
  return found.size();
};

//IMPLEMENT
/*
  int hideDataBelowYErr(double err); 
  int showDataBelowYErr(double err);
*/  
int EnsemData::hideDataAboveYErrRat(double rat){
  vector<double> y_rat;
  for(int i=0; i < totalNData; i++){
    y_rat.push_back( abs( toDouble(sqrt(variance(peekObs(y_data,i)))) / toDouble(mean( peekObs(y_data,i) ) ) ) );
  } 
  vector<int> found = find_above(y_rat, rat);
  for(int i = 0; i < found.size(); i++){
    active_data[found[i]] = false;  initInvCov = false;
  }
  return found.size();
};


/*
  int showDataAboveYErrRat(double rat);
  int hideDataBelowYErrRat(double rat);
  int showDataBelowYErrRat(double rat);
*/


void EnsemData::hideAll(){
  for(int i=0; i< totalNData; i++){active_data[i] = false;}
  initInvCov = false; 
}
void EnsemData::showAll(){
  for(int i=0; i< totalNData; i++){active_data[i] = true;}
   initInvCov = false;
} 


int EnsemData::getNData() const{
  int count = 0;
  for(int i=0; i < totalNData; i++){
    if(active_data[i]){count++;}
  }
  return count;
}


vector<double> EnsemData::getXData() const{
  vector<double> out;
  for(int i = 0; i < totalNData; i++){
    if(active_data[i]){out.push_back(x_data[i]);}
  }
  return out;
};

  
EnsemVectorReal EnsemData::getYData() const{
  EnsemVectorReal out; out.resize(nBins); out.resizeObs(getNData());
  int j = 0;
  for(int i = 0; i < totalNData; i++){
    if(active_data[i]){pokeObs(out, peekObs(y_data, i), j); j++;}
  }
  return out;
};

EnsemVectorReal EnsemData::getScaledYData() const{
  EnsemVectorReal out = getYData();
  return rescaleEnsemDown(out);
};

  
itpp::mat EnsemData::getCov() const{
  if(!initCov){makeCov();} //cov is mutable so this should be OK
  
  itpp::mat out(1, totalNData);
  for(int i = 0; i < totalNData; i++){
    if(active_data[i]){
      itpp::vec row = cov.get_row(i);
      out.append_row(row);}
  }
  out.del_row(0);

  for(int i = totalNData - 1; i > -1; i--){
    if(!active_data[i]){
      out.del_col(i);
    }
  }

  if( (out.rows() != getNData()) && (out.cols() != getNData() ) )
    {cerr << "cov matrix isn't the right size" << endl; exit(1);}

  return out;
}

itpp::mat EnsemData::getAllCov() const{
  if(!initCov){makeCov();}
  return cov;
};


EnsemReal EnsemData::getYUsingNearestX(double& x){
  int pos = find_nearest(x_data, x);
  x = x_data[pos];
  EnsemReal out =  peekObs(y_data, pos);
  return out;
}


// INVERSE COVARIANCE

void EnsemData::makeInvCov(double cutoff) const{
  //  cout << "inverting the covariance" << endl;

  itpp::Real_Timer timer; timer.tic();

  itpp::mat activeCov = getCov();
  // cout << "dimension " << activeCov.rows() << endl;


  invCov = invertSVDNorm(activeCov, cutoff, nResetCovSingVals);
  if(nResetCovSingVals > 0)
    {cout << "     " << __PRETTY_FUNCTION__ << " : inverted data covariance with " << nResetCovSingVals << " singular values reset from " << getNData() << endl;} 
  initInvCov = true;

  //  cout << "inverted the data covariance " << endl;

  double time = timer.toc();
  if(time > 10.0){
    cerr << "TIMING: inverse of covariance dimension " << getNData() << " took " << time << " seconds" << endl;} 
}

void EnsemData::makeInvCov() const {
  makeInvCov(SVCutoff);
}

itpp::mat EnsemData::getInvCov() const {
  if(!initInvCov){
    if(!initCov){makeCov();}
    makeInvCov();
  }
  return invCov;
}

int EnsemData::getNResetCovSingVals() const {
  //if you ask for this and it doesn't exist you'll get the inversion done
  if(!initInvCov){
    if(!initCov){makeCov();}
    makeInvCov();
  }
  return nResetCovSingVals;
}



void EnsemData::augmentCov(EnsemReal y){
  if(!initCov){cerr << " if you're reading this Jo screwed up the initCov flags" << endl;}  

  if(y.size() != nBins){cerr << "binsize doesn't match" << endl; exit(1);}
  itpp::vec row(totalNData);
  for(int i=0; i < totalNData; i++){
    Real mi = mean( peekObs(y_data , i) ); EnsemReal ei = peekObs(y_data , i);
    double covar = 0.;
    for(int bin = 0; bin < nBins; bin++){
      covar += toDouble( (peekEnsem(ei , bin) - mi ) * (peekEnsem(y , bin) - mean(y) ) );
    }
    row(i) = covar;
  }
  row /= nBins * (nBins - 1. );
  cov.append_row(row);
  
  row.set_size(totalNData + 1, true);
  double covar = 0.;
  for(int bin = 0; bin < nBins; bin++){
    covar += toDouble( (peekEnsem(y , bin) - mean(y) ) * (peekEnsem(y , bin) - mean(y) ) );
  }
  row(totalNData) = covar / (nBins * (nBins - 1.));
  cov.append_col(row);	
};

//**********************************************
// helper functions
//*********************************************

//concatenate data sets
 EnsemData concat(EnsemData a, EnsemData b){
  if(a.getNBins() != b.getNBins()){cerr << "binsize doesn't match" << endl; exit(1);}
  EnsemData c = a; //nb asymmetric - copies the covariance state of a, not b

  vector<double> x_b = b.getAllXData();
  EnsemVectorReal y_b = b.getAllYData();
  vector<bool> active_b = b.getActiveDataList();

  for(int i = 0; i < b.getTotalNData(); i++){
    double x = x_b[i]; 
    EnsemReal y = peekObs(y_b, i);  
    bool active = active_b[i];
    c.addDatum(x, y);
    if(!active){c.hideDatumByX(x);}
  }
  return c;
 };

ostream& operator<<(ostream& output, const EnsemData& e){
  output.precision(3); output << scientific;
  int nData = e.getTotalNData();
  vector<double> x_data = e.getAllXData();
  EnsemVectorReal y_data = e.getAllYData();
  vector<bool> active = e.getActiveDataList();

  output << "index\tx\t\tmean(y)\t\terror(y)\tratio" << endl;
  for(int i = 0; i < nData; i++){
    string act = "*"; if(active[i]){act = "";};
    output << act << i << "\t" << act << x_data[i] << "\t" << act <<toDouble(mean( peekObs(y_data , i) ) ) 
	   << "\t" << act << toDouble(sqrt(variance( peekObs(y_data , i) ) ) ) 
	   <<"\t" << act <<  toDouble(sqrt(variance( peekObs(y_data , i) ) )  / mean( peekObs(y_data , i) ) ) << endl;
  }
  
  return output;
}


double ensemDataChisq(const EnsemData& data, const EnsemData& theory){
  //data and theory should be evaluated at the same places
  //naughty - assume the data are ordered the same !
  if(data.getXData() != theory.getXData()){cerr << "data sets not evaluated at the same x-values" << endl; exit(1);}
  
  vector<double> x = data.getXData();
  EnsemVectorReal y_data = data.getYData();
  EnsemVectorReal y_theory = theory.getYData();

  //DONT use the covariance of the difference
  // EnsemVectorReal y_diff = y_data - y_theory;
  // EnsemData diff(x, y_diff);

  //  itpp::mat invcov = diff.getInvCov();
  itpp::mat invcov = data.getInvCov();

  //chisq computed using the mean data
  double chisq = 0.0;
  vector<double> dif; dif.resize(x.size());

  for(int i = 0; i < x.size(); i++){
    dif[i] = toDouble(mean( peekObs(y_data , i) - peekObs(y_theory , i) ) );
  }

  for(int i = 0; i < x.size(); i++){ 
    chisq += dif[i] * invcov(i,i) * dif[i];
    for(int j = i + 1; j < x.size(); j++){
      chisq += 2.0 * dif[i] * invcov(i,j) * dif[j];
    }
  }

  return chisq;
  
};

