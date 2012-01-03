#include "multi_princorr_fit.h"


FitMultiplePrincipalCorrelator::FitMultiplePrincipalCorrelator(EnsemData data_, int nCorrs_, vector<int> t0_list_, int tLen_, vector<int> tmin_,  double noiseRatioCutoff) : nCorrs(nCorrs_), t0_list(t0_list_), tLen(tLen_), data(data_), tmin(tmin_){

  //turn off t0's
  for(int n=0; n < nCorrs; n++){
    int t0 = t0_list[n] + n*tLen;
    data.hideDataByXRange( t0 - 0.1 ,  t0 + 0.1 );
  }
  
  //exclude all points with negative mean value
  data.hideDataBelowY(0.0);

  //exclude noisy points
  data.hideDataAboveYErrRat(noiseRatioCutoff); 


  //shift the tmins & t0's
  vector<int> t0_list_shifted;
  vector<int> tmin_list;
  for(int n=0; n < nCorrs; n++){
    tmin_list.push_back(tmin[n] + n*tLen);
    t0_list_shifted.push_back(t0_list[n] + n*tLen);
  }

  //find the largest t-value being considered for each corr
  vector<double> tslices = data.getXData();

  cout << "tmin, tmax with shifted times:" << endl;
  vector<int> tmax_list;
  for(int n=0; n < nCorrs; n++){
    int tmax = n*tLen;
    for(int i = 0; i < tslices.size(); i++){
      if( (int(tslices[i]) < (n+1)*tLen) && (int(tslices[i]) >= n*tLen) && (int(tslices[i]) > tmax) ){ tmax = int(tslices[i]); } 
    }
    cout << "for corr " << n << "   tmin = " << tmin_list[n] << "  t0 = " << t0_list_shifted[n] <<  " , tmax = " << tmax << endl;
    tmax_list.push_back(tmax);
  }
  

  MultiPrinCorrTwoExp* tmp2 = new MultiPrinCorrTwoExp(nCorrs, t0_list_shifted, tLen);
  Handle<FitFunction> twoExp(tmp2);

  twoExp->setDefaultParValue("mass_0", 0.5 ); 
  twoExp->setDefaultParError("mass_0", 0.2);
  //limit the mass to positive values
  twoExp->setParamLowerLimit("mass_0", 0.0);
  
  for(int n=0; n < nCorrs; n++){
    stringstream dum; dum <<  "m_1_corr" << n;
    twoExp->setDefaultParValue(dum.str(), 1.0);
    twoExp->setDefaultParError(dum.str(), 0.3); 
    twoExp->setParamLowerLimit(dum.str(), 0.0);
    
    stringstream dum2; dum2 << "A_corr" << n;
    twoExp->setDefaultParValue(dum2.str(), 0.1);
    twoExp->setDefaultParError(dum2.str(), 0.1);
    twoExp->setParamLowerLimit(dum2.str(), 0.0);
  }
  
  //turn off the data below the tmin values
  for(int i = 0; i < tslices.size(); i++){
    int t = int(tslices[i]);
    for(int n=0; n < nCorrs; n++){
      if( (t < tmin_list[n]) && (t >= n*tLen) ){data.hideDatumByX(double(t));}
    }
  }

//  cout << "this is the data\n";
//  cout << data << endl;
  
  JackFit theFit(data, twoExp);
  theFit.runAvgFit();
  theFit.runJackFit();

  cout << endl << endl;
  cout << "jack fit ran with chisq = " << theFit.getJackChisq() << " for nDoF = " << theFit.getNDoF() <<  "    chisq/ndof = " <<      theFit.getJackChisq() / theFit.getNDoF()  <<  endl;

  chisq = theFit.getJackChisq();
  nDoF = theFit.getNDoF();


  mass_0 = theFit.getJackFitParValue(0);
  cout << "mass_0 = " << toDouble(mean(mass_0)) << " +/- " << toDouble(sqrt(variance(mass_0))) << endl << endl;
  
  for(int n=0; n < nCorrs; n++){
    stringstream dum; dum <<  "mass_1_corr" << n;
    EnsemReal a = theFit.getJackFitParValue(1 + 2*n);
    cout << dum.str() << " = " <<  toDouble(mean(a)) << " +/- " << toDouble(sqrt(variance(a))) << endl;
    mass_1.push_back(a);

    stringstream dum2; dum2 <<  "A_corr" << n;
    EnsemReal b = theFit.getJackFitParValue(2 + 2*n);
    cout << dum2.str() << " = " <<  toDouble(mean(b)) << " +/- " << toDouble(sqrt(variance(b))) << endl << endl;
    A.push_back(b);
  }
 
  //make the plot
  MultiCorrWeight weight(mass_0, nCorrs, t0_list_shifted, tLen);  

  stringstream lab; 
  lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << theFit.getJackChisq() << "/" << theFit.getNDoF(); 
  lab << "; m=" << fixed << setprecision(4) << toDouble(mean(mass_0)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(mass_0)));
  
  axis_plot = theFit.makeJackFitPlotAxis(weight, 0.0, double(nCorrs*tLen), lab.str() );  

};

  void FitMultiplePrincipalCorrelator::saveFitPlot(string filename){
  ofstream out; out.open(filename.c_str());
  out << axis_plot;
  out.close();
}

