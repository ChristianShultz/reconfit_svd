#include "multi_ensemble_fitter.h"


MultiEnsemFitter::MultiEnsemFitter(map<string, EnsemData> data_, Handle<FitFunction> ff_) : data(data_), ff(ff_){};
  
bool MultiEnsemFitter::fit(){
  //convert the (active) data to a long list with a block diagonal covariance
  vector<double> x;
  vector<double> y;
  vector<itpp::mat> inv_covs;

  int nResetSingVal = 0;

  //loop over the ensembles
  for(map<string, EnsemData>::iterator it = data.begin(); it != data.end(); it++){
    EnsemVectorReal ytmp = ((*it).second).getYData();
    for(int i=0; i < ytmp.numElem(); i++){ y.push_back( toDouble( mean( peekObs(ytmp, i) ) ) ) ; }

    vector<double> xtmp = ((*it).second).getXData();
    for(int i=0; i < xtmp.size(); i++){ x.push_back( xtmp[i] ) ; }

    itpp::mat covtmp = ((*it).second).getInvCov();
    inv_covs.push_back(covtmp);

    nResetSingVal += ((*it).second).getNResetCovSingVals();
  }

  int data_size = x.size();
  itpp::sparse_mat sparse(data_size, data_size); 

  //build the block diagonal inverse covariance:
  int dum = 0;
  for(int i=0; i < inv_covs.size(); i++){  
    sparse.set_submatrix(dum, dum, inv_covs[i]);
    dum += inv_covs[i].rows();
  }
  itpp::mat big_inv_cov = sparse.full();

  //  cout << "big_inv_cov = " << endl << big_inv_cov << endl;

  //construct the chisq object for Minuit
  ChiSquare chisquare(x, y, big_inv_cov, ff);

  //set up the pars
  MnUserParameters upar;
  for(int i = 0; i < ff->getNPars(); i++)
    { upar.Add(  (ff->getParName(i)).c_str(), ff->getDefaultParValue(i), ff->getDefaultParError(i));} 
  //fix if required
  for(int i = 0; i < ff->getNPars(); i++){ if(ff->isParamFixed(i)){upar.Fix(i);}; }
  //limit if required
  for(int i = 0; i < ff->getNPars(); i++){ 
    if(ff->isParamUpperLimited(i)){upar.SetUpperLimit(i, ff->getParamUpperLimit(i));};
    if(ff->isParamLowerLimited(i)){upar.SetLowerLimit(i, ff->getParamLowerLimit(i));};
  }

  //  cout << "MnUserParameters = " << endl << upar << endl;

  //create the minimiser
  MnMigrad mini(chisquare, upar);
  //MINIMIZE
  FunctionMinimum min = mini();
  fitSuccess = min.IsValid();

  for(int i = 0; i < ff->getNPars(); i++){
    vFitParValue.push_back( min.UserState().Value(i) );
    vFitParError.push_back( min.UserState().Error(i) );

    fitParValue[ff->getParName(i)] = min.UserState().Value(i);
    fitParError[ff->getParName(i)] = min.UserState().Error(i);
  }
  chisq = min.Fval();
  stringstream report; report << min;
  fitReport = report.str();
  
  nDoF = x.size() - nResetSingVal - ff->getNUnfixedPars();

  return fitSuccess;

};

string MultiEnsemFitter::makeFitPlotAxis(double xmin, double xmax, string label){

  AxisPlot plot;

  int colour = 1;
  //add the raw data
  for(map<string, EnsemData>::iterator it = data.begin(); it != data.end(); it++){
    EnsemVectorReal y = ((*it).second).getYData();
    vector<double> x = ((*it).second).getXData();
    plot.addEnsemData(x, y, "\\sq", colour);
    colour++;
  }

  double dx = (xmax - xmin) / 101;
  
  //add the fit 
  vector<double> xx;
  vector<double> yy, ype, yme;

  //computing the error bands will need the covariance between the parameters - ADD LATER
  vector<double> pars;
  for(int i = 0; i < vFitParValue.size(); i++ ){
    pars.push_back( vFitParValue[i] );
  }

  for(int i=0; i < 101; i++){
    double x = xmin + i*dx;
    xx.push_back(x);
    yy.push_back( (*ff)(pars, x) );
  }

  plot.addLineData(xx, yy, 1);

  //setting the ranges
  plot.setXRange(xmin, xmax);

  double ymin = *min_element( yy.begin(), yy.end() );
  double ymax = *max_element( yy.begin(), yy.end() );
  plot.setYRange(ymin, ymax);

  plot.addLabel(xmin + 10*dx, ymax - 0.05*(ymax-ymin), label, 1, 1.0);

  return plot.getAxisPlotString();
};
