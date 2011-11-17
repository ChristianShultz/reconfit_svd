#include "jackknife_fitter.h"

//*******************************************************************
//Fit Function
//*******************************************************************
FitFunction::FitFunction(int nPars_) : nPars(nPars_){
  

  //default par names
  parNames.resize(nPars);
  for(int i = 0; i < nPars; i++){
    stringstream ss; ss << i; 
    parNames[i] = ss.str();
  }
  
  //default default values & errors
  defaultParValues.resize(nPars);
  defaultParErrors.resize(nPars);
  for(int i = 0; i < nPars; i++){
    defaultParValues[i] = 1.0;
    defaultParErrors[i] = 1.0;
  }

  //by default nothing is fixed
  parFixed.resize(nPars);
  parUpperLimited.resize(nPars);
  parLowerLimited.resize(nPars);
  parLowLimits.resize(nPars);
  parHighLimits.resize(nPars);
  for(int i = 0; i < nPars; i++){
    parFixed[i] = false;
    parUpperLimited[i] = false;
    parLowerLimited[i] = false;    
  }
};

//default to a zero function
//double FitFunction::operator()(const vector<double>& pars, double x) const { return 0.0; }

int FitFunction::getNUnfixedPars() const {
  int dum = 0;
  for(int i=0; i < nPars; i++){
    if(!isParamFixed(i)){dum++;}
  }
  return dum;
};

void FitFunction::setParNames(vector<string> names){
  if(names.size() != nPars){cerr << "too many names" << endl; exit(1);}
  parNames = names;

  //check for names that differ only after the tenth char
  vector<string> shortname = parNames;
  for(int i = 0; i < parNames.size(); i++){
    shortname[i].resize(10);
  }
  sort( shortname.begin(), shortname.end() );
  vector<string> compare = shortname;
  vector<string>::iterator pos = unique(compare.begin(), compare.end());
  compare.erase( pos, compare.end() );
  if(compare != shortname){cerr << "some parameter names do not differ in the first ten chars - Minuit can't cope with this" << endl; exit(1);}

}
void FitFunction::setParName(int parNum, string name){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  parNames[parNum] = name;

  //check for names that differ only after the tenth char
  vector<string> shortname = parNames;
  for(int i = 0; i < parNames.size(); i++){
    shortname[i].resize(10);
  }
  sort( shortname.begin(), shortname.end() );
  vector<string> compare = shortname;
  vector<string>::iterator pos = unique(compare.begin(), compare.end());
  compare.erase( pos, compare.end() );
  if(compare != shortname){cerr << "some parameter names do not differ in the first ten chars - Minuit can't cope with this" << endl; exit(1);}
}


string FitFunction::getParName(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return parNames[parNum];
}
int FitFunction::getParNum(string name) const {
  bool found = false;
  int pos;
  for(int i =0 ; i < nPars; i++){
    if(parNames[i] == name){found = true; pos = i; break;}
  }
  if(!found){ cerr << "no such parameter" << endl; exit(1);}
  return pos;
}


void FitFunction::setDefaultParValues(vector<double> values){
  if(values.size() != nPars){cerr << "too many names" << endl; exit(1);}
  defaultParValues = values;
}
void FitFunction::setDefaultParValue(int parNum, double value){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  defaultParValues[parNum] = value;
}
void FitFunction::setDefaultParValue(string name, double value){
  int parNum = getParNum(name);
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  defaultParValues[parNum] = value;
}
double FitFunction::getDefaultParValue(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return defaultParValues[parNum];
}
double FitFunction::getDefaultParValue(string name) const {
  int parNum = getParNum(name);
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return defaultParValues[parNum];
}


void FitFunction::setDefaultParErrors(vector<double> errors){
  if(errors.size() != nPars){cerr << "too many names" << endl; exit(1);}
  defaultParErrors = errors;
}
void FitFunction::setDefaultParError(int parNum, double error){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  defaultParErrors[parNum] = error;
}
void FitFunction::setDefaultParError(string name, double error){
  int parNum = getParNum(name);
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  defaultParErrors[parNum] = error;
}
double FitFunction::getDefaultParError(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return defaultParErrors[parNum];
}
double FitFunction::getDefaultParError(string name) const {
  int parNum = getParNum(name);
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return defaultParErrors[parNum];
}


void FitFunction::fixParam(int parNum){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  parFixed[parNum] = true;
}
void FitFunction::fixParam(string name){
  int parNum = getParNum(name);
  fixParam(parNum);
}
void FitFunction::releaseParam(int parNum){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  parFixed[parNum] = false;
}
void FitFunction::releaseParam(string name){
  int parNum = getParNum(name);
  releaseParam(parNum);
}
bool FitFunction::isParamFixed(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return parFixed[parNum];
}
bool FitFunction::isParamFixed(string name) const {
  int parNum = getParNum(name);
  return isParamFixed(parNum);
}


void FitFunction::setParamUpperLimit(int parNum, double value){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  parHighLimits[parNum] = value;
  parUpperLimited[parNum] = true;
}
void FitFunction::setParamUpperLimit(string name, double value){
  int parNum = getParNum(name);
  parHighLimits[parNum] = value;
  parUpperLimited[parNum] = true;
}
void FitFunction::setParamLowerLimit(int parNum, double value){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  parLowLimits[parNum] = value;
  parLowerLimited[parNum] = true;
}
void FitFunction::setParamLowerLimit(string name, double value){
  int parNum = getParNum(name);
  parLowLimits[parNum] = value;
  parLowerLimited[parNum] = true;
}
void FitFunction::setParamLimits(int parNum, double low, double high){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  parLowLimits[parNum] = low;
  parHighLimits[parNum] = high;
  parUpperLimited[parNum] = true; parLowerLimited[parNum] = true;
}
void FitFunction::setParamLimits(string name, double low, double high){
  int parNum = getParNum(name);
  parLowLimits[parNum] = low;
  parHighLimits[parNum] = high;
  parUpperLimited[parNum] = true; parLowerLimited[parNum] = true;
}
void FitFunction::releaseParamLimits(int parNum){
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  parUpperLimited[parNum] = false; parLowerLimited[parNum] = false;
}
void FitFunction::releaseParamLimits(string name){
  int parNum = getParNum(name);
  parUpperLimited[parNum] = false; parLowerLimited[parNum] = false;
}

double FitFunction::getParamUpperLimit(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return parHighLimits[parNum];
}
double FitFunction::getParamUpperLimit(string name) const {
  int parNum = getParNum(name);
  return parHighLimits[parNum];
}
double FitFunction::getParamLowerLimit(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return parLowLimits[parNum];
}
double FitFunction::getParamLowerLimit(string name) const {
  int parNum = getParNum(name);
  return parLowLimits[parNum];
}


bool FitFunction::isParamUpperLimited(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return parUpperLimited[parNum];
}
bool FitFunction::isParamUpperLimited(string name) const {
  int parNum = getParNum(name);
  return parUpperLimited[parNum];
}
bool FitFunction::isParamLowerLimited(int parNum) const {
  if((parNum > nPars)||(parNum < 0)){cerr << "this param number out of range" << endl; exit(1);}
  return parLowerLimited[parNum];
}
bool FitFunction::isParamLowerLimited(string name) const {
  int parNum = getParNum(name);
  return parLowerLimited[parNum];
}

//******************************************************************************************

//*******************************************************************************************
double ChiSquare::operator()(const vector<double>& pars) const{
  if(pars.size() != ff->getNPars()){cerr << "chisquared can't be computed, wrong number of params" << endl; exit(1);}
  
  int n = x_data.size();
  double chisq = 0;
  vector<double> y_theory;
  for(int i = 0; i < n; i++){
    y_theory.push_back( (*ff)(pars, x_data[i])   );
  }

  for(int i = 0; i < n; i++){
    double diff_i = y_data[i] - y_theory[i];
    for(int j = 0; j < n; j++){
      double diff_j = y_data[j] - y_theory[j];
      chisq += diff_i * invcov(i,j) * diff_j;
    }
  }

  return chisq;
};



//********************************************************************************************
//           JACKKNIFE FITTER
//*******************************************************************************************
JackFit::JackFit(EnsemData data_, Handle<FitFunction> ff_): data(data_), ff(ff_){ 
  //data is a fixed copy - if you want to change it, construct another fit object
  //ff is virtual so we use a pointer

  int nBins = data.getNBins();
  jackChisqs.resize(nBins);
  jackFitSuccess.resize(nBins);
  jackFitReports.resize(nBins);

  initJack = false;
  avgFitSuccess = false;

  //fill the jackknife vectors with zeroes
  EnsemReal zero; zero.resize(nBins); zero = Real(0.0);
  vector<EnsemReal> veczero(ff->getNPars(), zero);
  scaledJackParValues = veczero;    
};



bool JackFit::runBinFit(int bin, vector<double> startValues, vector<double> startErrors){
  if( (bin > data.getNBins() ) || (bin < -1) ){cerr << "bin out of range" << endl; exit(1);}

  //prepare the data on this bin
  vector<double> y;

  if(bin == -1){ //average data
    EnsemVectorReal yy = data.getYData();
    for(int i=0; i < data.getNData(); i++){
      y.push_back( toDouble( mean( peekObs(yy, i) ) ) );
    }
  }
  else{ //pull out one bin
    if(!initJack){ cerr << "run a jack fit first before you play with individual bins" << endl; exit(1);}
    EnsemVectorReal yscaled = data.getScaledYData();
    for(int i=0; i < data.getNData(); i++){
      y.push_back( toDouble( peekEnsem( peekObs(yscaled, i) , bin ) ) );
    }
  }

  //construct the chisq object for Minuit
  ChiSquare chisquare(data.getXData(), y, data.getInvCov(), ff);

  //set up the pars
  MnUserParameters upar;
  for(int i = 0; i < ff->getNPars(); i++)
    { upar.Add(  (ff->getParName(i)).c_str(), startValues[i], startErrors[i]);}
  //fix if required
  for(int i = 0; i < ff->getNPars(); i++){ if(ff->isParamFixed(i)){upar.Fix(i);}; }
  //limit if required
  for(int i = 0; i < ff->getNPars(); i++){ 
    if(ff->isParamUpperLimited(i)){upar.SetUpperLimit(i, ff->getParamUpperLimit(i));};
    if(ff->isParamLowerLimited(i)){upar.SetLowerLimit(i, ff->getParamLowerLimit(i));};
  }

  //  cout << "MnUserParameters = " << endl << upar << endl;

  //cout << "in runBinFit, set all the parameter ranges etc...." << endl;

  //create the minimiser
  MnMigrad mini(chisquare, upar);
  //MINIMIZE
  FunctionMinimum min = mini();
  bool fitSuccess = min.IsValid();
 
  /*  cout << "in runBinFit, ran a fit, was ";
  if(fitSuccess){cout << "successful";}
  else{cout << "unsuccessful";}
  cout << endl;
  */
  //dump the result
  if(bin == -1){
    avgParValues.resize(ff->getNPars());
    avgParErrors.resize(ff->getNPars());
    for(int i=0; i < ff->getNPars(); i++){ 
      avgParValues[i] = min.UserState().Value(i);
      avgParErrors[i] = min.UserState().Error(i);
    }
    avgChisq = min.Fval();
    avgFitSuccess = fitSuccess;
    stringstream report; report << min;
    avgFitReport = report.str();

    //debug
    //    cout << "AVERAGE CHISQ = " << avgChisq << endl;
    // cout << "mass_0 = " << avgParValues[0] << endl;
    //
  }
  else{
    for(int i = 0; i < ff->getNPars(); i++)
      {pokeEnsem(scaledJackParValues[i] , Real(min.UserState().Value(i)) , bin);}
    jackChisqs[bin] = min.Fval();
    jackFitSuccess[bin] = fitSuccess;
    stringstream report; report << min;
    jackFitReports[bin] = report.str();
  }
    
  return fitSuccess;
}


bool JackFit::runAvgFit(){ 
  return runBinFit(-1, ff->getDefaultParValues(), ff->getDefaultParErrors()); 
}

vector<bool> JackFit::runJackFit(){
  if(!avgFitSuccess){runAvgFit();} //if no average fit has been done successfully - do one now
  initJack = true;

  //run over every bin using the avg fit results as the start guesses for the fit
  vector<bool> success;
  for(int bin = 0; bin < data.getNBins(); bin++){ success.push_back( runBinFit(bin, avgParValues, avgParErrors) );}
  return success;
}

double JackFit::getAvgFitParValue(int n) const {
  if((n > ff->getNPars() )||(n < 0)){ cerr << "parnum out of range" << endl; exit(1);}
  return avgParValues[n];
}

double JackFit::getAvgFitParValue(string name) const {
  int n = ff->getParNum(name);
  return avgParValues[n];
}

double JackFit::getAvgFitParError(int n) const {
  if((n >  ff->getNPars() )||(n < 0)){ cerr << "parnum out of range" << endl; exit(1);}
  return avgParErrors[n];
}

double JackFit::getAvgFitParError(string name) const {
  int n = ff->getParNum(name);
  return avgParErrors[n];
}

//investigate jack fit results
vector<EnsemReal> JackFit::getJackFitParValues() const {
  //rescale now
  vector<EnsemReal> jackParValues;
  for(int i = 0; i <  ff->getNPars(); i++){
    EnsemReal temp = rescaleEnsemUp(scaledJackParValues[i]);
    jackParValues.push_back(temp);
  }
  return jackParValues;
}

EnsemReal JackFit::getJackFitParValue(int n) const {
  return rescaleEnsemUp(scaledJackParValues[n]);
}

EnsemReal JackFit::getJackFitParValue(string name) const {
  int n = ff->getParNum(name);
  return rescaleEnsemUp(scaledJackParValues[n]);
}

double JackFit::getJackChisq() const{
  double mean = 0.0;
  for(int bin = 0; bin < data.getNBins(); bin++){mean += jackChisqs[bin];}
  mean /= data.getNBins();
  return mean;
}

int JackFit::getNFailedFits() const{
  int count = 0;
  for(int bin = 0; bin < data.getNBins(); bin++){ if(!jackFitSuccess[bin]){count++;};}
  return count;
}

vector<int> JackFit::getFailedFitBins() const{
  vector<int> failedList;
  for(int bin = 0; bin < data.getNBins(); bin++){ if(!jackFitSuccess[bin]){failedList.push_back(bin);};}
  return failedList;
}

string JackFit::getJackFitReport(int bin) const{
  return jackFitReports[bin];
}

bool JackFit::getJackFitSuccess() const{
  bool temp = false;
  if(getNFailedFits() == 0){temp = true;}
  return temp;
}

string JackFit::makeJackFitPlotAxis(double xmin, double xmax, string label){
  int nBins = data.getNBins();
  EnsemReal dum; dum.resize(nBins);
  One one(dum); //EnsemFunction weighting that does nothing
  return makeJackFitPlotAxis(one, xmin, xmax, label);
}

string JackFit::makeJackFitPlotAxis(EnsemFunction& weightFn, double xmin, double xmax, string label){
  vector<double> x; 
  vector<double> fit_mean; 
  vector<double> fit_mean_plus_err;
  vector<double> fit_mean_minus_err;
  
  double dx = (xmax - xmin) / 100.; //hardwired 100 point resolution
  
  for(int ix = 0; ix < 101; ix++){
    double xx = xmin + (dx * double(ix));
    x.push_back(xx);

    EnsemReal yscaled; yscaled.resize(data.getNBins());
    EnsemReal y = yscaled;

    for(int bin = 0; bin < data.getNBins(); bin++){
      vector<double> pars;
      for(int i = 0; i < ff->getNPars(); i++){ pars.push_back( toDouble( scaledJackParValues[i].elem(bin) ) );}
      double f = (*ff)(pars, xx);
      pokeEnsem(yscaled, Real(f), bin);
    }
    
    y = rescaleEnsemUp(yscaled);
    EnsemReal w = weightFn(xx);
    y *= w;

    double m = toDouble(mean(y)); double e = toDouble(sqrt(variance(y)));
    fit_mean.push_back( m );
    fit_mean_plus_err.push_back( m + e );
    fit_mean_minus_err.push_back( m - e );    
  }

  //  cout << "DEBUG:: made the data" << endl;	

  //make the plot
  AxisPlot plot;
  plot.setXRange(xmin, xmax);
    
  //find the lowest/highest x in the fit
  vector<double> x_data = data.getXData();
  double xfit_min = *min_element(x_data.begin(), x_data.end());
  double xfit_max = *max_element(x_data.begin(), x_data.end());
  
  //write the fit plot
  //below the fit region
  vector<double> xdum;
  vector<double> d, dp, dm; 
  int ix = 0;
  while (x[ix] <= xfit_min){
    xdum.push_back(x[ix]);
    d.push_back(fit_mean[ix]); dp.push_back(fit_mean_plus_err[ix]); dm.push_back(fit_mean_minus_err[ix]); 
    ix++;
  }
  plot.addLineData(xdum, d, 3); plot.addLineData(xdum, dp, 3); plot.addLineData(xdum, dm, 3);  
  ix--;

  //in the fit region
  xdum.clear(); d.clear(); dp.clear(); dm.clear();
  while (x[ix] <= xfit_max){
    xdum.push_back(x[ix]);
    d.push_back(fit_mean[ix]); dp.push_back(fit_mean_plus_err[ix]); dm.push_back(fit_mean_minus_err[ix]); 
    ix++;
  }
  plot.addLineData(xdum, d, 2); plot.addLineData(xdum, dp, 2); plot.addLineData(xdum, dm, 2);
  ix--;

  //above the fit region
  xdum.clear(); d.clear(); dp.clear(); dm.clear();
  while (ix < 101){
    xdum.push_back(x[ix]);
    d.push_back(fit_mean[ix]); dp.push_back(fit_mean_plus_err[ix]); dm.push_back(fit_mean_minus_err[ix]); 
    ix++;
  }
  plot.addLineData(xdum, d, 3); plot.addLineData(xdum, dp, 3); plot.addLineData(xdum, dm, 3);

  //  cout << "DEBUG: added the fit data" << endl;

  //add the raw data
  vector<EnsemReal> weighted_data;
  for(int i = 0; i < data.getNData(); i++){
    weighted_data.push_back( weightFn(x_data[i]) * peekObs(data.getYData() , i) );
  }
  plot.addEnsemData(x_data, weighted_data, "\\sq", 1);
  
  //add the hidden data
  if(data.getTotalNData() > data.getNData()){
    vector<EnsemReal> off_weighted_data;
    vector<double> off_x_data;
    for(int i = 0; i < data.getTotalNData(); i++){
      if( !((data.getActiveDataList())[i]) ){
	off_x_data.push_back( (data.getAllXData())[i] );
	off_weighted_data.push_back( weightFn( (data.getAllXData())[i] ) * peekObs(data.getAllYData() , i) );
      }
    }
    plot.addEnsemData(off_x_data, off_weighted_data, "\\di", 3);
  }

  //set the yrange using the active data
  vector<double> yp;
  vector<double> ym;
  for(int i=0; i < data.getNData(); i++){
    EnsemReal tmp = weighted_data[i];
    double m = toDouble(mean(tmp));
    double e = toDouble( sqrt( variance( tmp ) ) );
    yp.push_back( m + e );
    ym.push_back( m - e );
  }
  double ymax = *max_element(yp.begin(), yp.end());
  if(ymax > 0){ ymax *= 1.10; }else{ ymax *= 0.9; }
  
  double ymin = *min_element(ym.begin(), ym.end());
  if(ymin > 0){ ymin *= 0.9; }else{ ymin *= 1.10; }

  plot.setYRange(ymin, ymax);

  //add a label
  plot.addLabel(xmin + 10*dx , ymax - 0.05*(ymax - ymin) , label, 1 , 1.0);
 
  //  plot.sendToFile(filename);
  return plot.getAxisPlotString();

}



//*************************************************************
// JACKFITLOG
//*************************************************************


JackFitLog::JackFitLog(EnsemData data_) : data(data_){}

EnsemData& JackFitLog::getEnsemData(){
  return data;
}

void JackFitLog::addFit(string fitname, Handle<FitFunction> ff){
  JackFit fit(data, ff);
  fit.runAvgFit();

  FitDescriptor fitDesc(ff, data.getActiveDataList(), fitname);
  if(fit.getAvgFitSuccess()){ 
    keys.push_back(fitDesc);
    fits.push_back(fit); 
  } 
}

bool fitDescriptorPredicate(const FitDescriptor& a, const FitDescriptor& b){
  if( ( a.ff.operator->() == b.ff.operator->() ) && (a.activeData == b.activeData) && (a.fitname == b.fitname) ){ return true;}
  else{return false;}
}

void JackFitLog::removeFit(FitDescriptor fitDesc){

  vector<FitDescriptor>::iterator p;
  vector<FitDescriptor> dum; dum.push_back(fitDesc);
  
  p = search(keys.begin(), keys.end(), dum.begin(), dum.end(), fitDescriptorPredicate);
  int pos = p - keys.begin();
  keys.erase(p);
  
  vector<JackFit>::iterator pp;
  pp = fits.begin() + pos;
  fits.erase(pp);
}


map<double, FitDescriptor> JackFitLog::getFitList(FitComparator& fitComp) const{
  
  map<double, FitDescriptor> list;
  
  for(int i = 0; i < keys.size(); i++){
    double fitCrit = fitComp(keys[i] , fits[i]);
    list.insert( make_pair(fitCrit, keys[i]) );
  }
  
  return list;
}

JackFit& JackFitLog::getFit(FitDescriptor fitDesc){
  vector<FitDescriptor>::iterator p;
  vector<FitDescriptor> dum; dum.push_back(fitDesc);
  
  p = search(keys.begin(), keys.end(), dum.begin(), dum.end(), fitDescriptorPredicate);
  int pos = p - keys.begin();
  vector<JackFit>::iterator pp;
  pp = fits.begin() + pos;

  if(p != keys.end()){ return *pp; }
  else{ cerr << "this fit not present" << endl; exit(1);}
}

FitDescriptor JackFitLog::getBestFit(FitComparator& fitComp){
  //at least one fit must have suceeded

 map<double, FitDescriptor> list;
 list = getFitList(fitComp);
 
 //this use to be a failure mode.. 
 if(list.end() == list.begin())
   {
     cerr << "no fits present" << endl;
     exit(1);
   }
 map<double, FitDescriptor>::const_iterator p = list.end(); p--;

 return p->second;
}

FitDescriptor JackFitLog::getBestJackFit(FitComparator& fitComp, int& rank){

  map<double, FitDescriptor> list;
  list = getFitList(fitComp);

  map<double, FitDescriptor>::const_iterator p = list.end();
  bool success = false; rank = 0;

  while((!success) && (rank <= list.size()) ){
    p--;
    FitDescriptor thisFit = p->second;
    JackFit& bestFit = getFit(thisFit);
    bestFit.runJackFit();
    if(bestFit.getNFailedFits() > 0){} //cout << thisFit.fitname << "  " << bestFit.getNFailedFits()  <<" bins failed" << endl;} //write a log message ?
    else if(bestFit.getJackChisq() > 0){success = true;} //write a log message?
    //would like to also check that the jack fit is "compatible" with the average fit we started with ?
    
    rank++;
  }
  
  if(success){ return p->second;}
  else{ 
    cerr << "found no acceptable jackknife fits" << endl; 
    FitDescriptor fake = (list.begin())->second; //fake - assumes one entry at least
    fake.fitname = "FAILED";
    return fake;
  } 

}












