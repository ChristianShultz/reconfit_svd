#include "fit_correlators.h"

FitPrincipalCorrelator::FitPrincipalCorrelator(EnsemData data_, int t0_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices) : t0(t0_), fits(data_), fitComp(fitComp_){

  EnsemVectorReal y = (fits.getEnsemData()).getYData();
  vector<double> t_temp = (fits.getEnsemData()).getXData();

  //turn off t0
  (fits.getEnsemData()).hideDataByXRange( t0 - 0.1 ,  t0 + 0.1 );
 
  //exclude all points with negative mean value
  (fits.getEnsemData()).hideDataBelowY(0.0);

  //exclude noisy points
  (fits.getEnsemData()).hideDataAboveYErrRat(noiseRatioCutoff); 
  // alternatively can exclude all points after the time where the data is first noisy - not implemented

  vector<double> all_tslices = (fits.getEnsemData()).getAllXData();
  int tlow = all_tslices[0];

  if(all_tslices.size() < minTSlices){ cerr << "** not enough timeslices to have minTSlices = " << minTSlices << endl; exit(1); }

  //ensure enough timeslices to do a fit
  if((fits.getEnsemData()).getNData() < minTSlices ){
    (fits.getEnsemData()).hideAll();
    for(int t = tlow; t <= tlow + minTSlices; t++){
      (fits.getEnsemData()).showDatumByX(t);
    }
  }

  //find the largest t-value being considered 
  vector<double> tslices = (fits.getEnsemData()).getXData();
  int tmax = int( *max_element(tslices.begin(), tslices.end() ) );
    
  //=======================
  //oneExp fits
  //=======================
  PrinCorrOneExp* tmp1 = new PrinCorrOneExp(t0);
  Handle<FitFunction> oneExp(tmp1);

  oneExp->setDefaultParValue("mass_0", 0.2); //totally arbitrary - could use meff to guess
  oneExp->setDefaultParError("mass_0", 0.2); //totally arbitrary

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( (tmax - int(tslices[i]) + 1.1 ) >= minTSlices){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "oneExp tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), oneExp);
    } 
  }// next t min

  //find the best oneExp fit
  FitDescriptor bestOneExp = fits.getBestFit(*fitComp);
  JackFit bestOneExpFit = fits.getFit( bestOneExp );
 

  //since timeslice data will be ordered
  //tmin is given by the lowest true in active data
  int tminOneExp;
  vector<bool> active = bestOneExp.activeData;
  for(int i = 0; i < active.size(); i++){ if(active[i]){tminOneExp = int(all_tslices[i]); break;}; }

  //reinstate the data
  for(int i = 0; i < tslices.size(); i++){(fits.getEnsemData()).showDatumByX( tslices[i] );}

  //=======================
  //twoExp fits
  //=======================
  PrinCorrTwoExp* tmp2 = new PrinCorrTwoExp(t0);
  Handle<FitFunction> twoExp(tmp2);

  twoExp->setDefaultParValue("mass_0", bestOneExpFit.getAvgFitParValue("mass_0") ); 
  twoExp->setDefaultParError("mass_0", 5.0*bestOneExpFit.getAvgFitParError("mass_0"));
  
  twoExp->setDefaultParValue("mass_1", bestOneExpFit.getAvgFitParValue("mass_0") + 0.5);
  twoExp->setDefaultParError("mass_1", 10.0*bestOneExpFit.getAvgFitParError("mass_0"));
  twoExp->setDefaultParValue("A", 0.1);
  twoExp->setDefaultParError("A", 0.1);

  //limit the masses to positive values
  twoExp->setParamLowerLimit("mass_0", 0.0);
  twoExp->setParamLowerLimit("mass_1", 0.0);


  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( ((tmax - int(tslices[i]) + 1.1 ) >= minTSlices ) && (tslices[i] + 0.1 < tminOneExp ) ){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "twoExp tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), twoExp);
    } 
  }// next t min

  //get the best jack fit
  int rank;
  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);

  if(best.fitname == "FAILED"){
    //no jackknife fits found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED";
    fit_summary = "FAILED -  set mass_0 to 1.0";
    //no plot !!!
  }
  else if(toDouble(mean(fits.getFit(best).getJackFitParValue("mass_0"))) <= 0.0)
    {
    //negative mass found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED (negative mass)";
    fit_summary = "FAILED: fit with negative mass_0 - set mass_0 to 1.0";
    //no plot !!!
    }
  else{
  
    if( (best.ff).operator->() == twoExp.operator->() ){nExp = 2;}else{nExp = 1;};
    JackFit& bestFit = fits.getFit(best);
    
    mass_0 = bestFit.getJackFitParValue("mass_0");
    if(nExp ==2){
      mass_1 = bestFit.getJackFitParValue("mass_1");
      A = bestFit.getJackFitParValue("A");
    }
    
    chisq = bestFit.getJackChisq();
    nDoF = bestFit.getNDoF();
    best_fit_name = best.fitname;
    
    
    //write out a summary of the fits
    int count = 1;
    map<double, FitDescriptor> list = fits.getFitList(*fitComp);
    stringstream ss; 
    ss << "                         | chisq/nDoF |     Q      |  fitCrit   | " << endl;
    for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
      JackFit& thisFit = fits.getFit(p->second);
      double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
      double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
      ss << setw(25) <<(p->second).fitname << "|";
      ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
      ss << setw(12) << fixed << setprecision(3) << Q <<"|";
      ss << setw(12) << scientific << setprecision(3) << p->first << "|";
      
      if(count == rank){ ss << "*";}else{ ss << " ";}
      
      ss << " m=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("mass_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("mass_0");
      
      if( ((p->second).ff).operator->() == twoExp.operator->() ){
	ss << ", m'=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("mass_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("mass_1");
	ss << ", A=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("A") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("A");
      }
      ss << endl;
      count++;
    } 
    fit_summary = ss.str();
    
    //make the plot
    ConstTimesExp expWeight(exp(- mass_0 * Real(t0)) , mass_0);
    
    stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; m=" << fixed << setprecision(4) << toDouble(mean(mass_0)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(mass_0)));
    
    axis_plot = bestFit.makeJackFitPlotAxis(expWeight, 0.0, double(tmax + 5.5), lab.str() );
  }

};

void FitPrincipalCorrelator::saveFitPlot(string filename){
  ofstream out; out.open(filename.c_str());
  out << axis_plot;
  out.close();
}



///========================================================================================================================
FitCorrelatorExp::FitCorrelatorExp(EnsemData data_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices) : fits(data_), fitComp(fitComp_){

  EnsemVectorReal y = (fits.getEnsemData()).getYData();
  vector<double> t_temp = (fits.getEnsemData()).getXData();

  //exclude noisy points
  (fits.getEnsemData()).hideDataAboveYErrRat(noiseRatioCutoff); 
 
  vector<double> all_tslices = (fits.getEnsemData()).getAllXData();
  int tlow = all_tslices[0];

  if(all_tslices.size() < minTSlices){ cerr << "** not enough timeslices to have minTSlices = " << minTSlices << endl; exit(1); }

  //ensure enough timeslices to do a fit
  if((fits.getEnsemData()).getNData() < minTSlices ){
    (fits.getEnsemData()).hideAll();
    for(int t = tlow; t <= tlow + minTSlices; t++){
      (fits.getEnsemData()).showDatumByX(t);
    }
  }

  //find the largest t-value being considered 
  vector<double> tslices = (fits.getEnsemData()).getXData();
  int tmax = int( *max_element(tslices.begin(), tslices.end() ) );

  //  cout << "about to start oneExp fits" << endl;
    
  //=======================
  //oneExp fits
  //=======================
  CorrOneExp* tmp1 = new CorrOneExp();
  Handle<FitFunction> oneExp(tmp1);

  oneExp->setDefaultParValue("mass_0", 0.2); //totally arbitrary - could use meff to guess
  oneExp->setDefaultParError("mass_0", 0.2); //totally arbitrary

  oneExp->setDefaultParValue("amp_0", 1.0); //totally arbitrary
  oneExp->setDefaultParError("amp_0", 1.0); //totally arbitrary

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( (tmax - int(tslices[i]) + 1.1 ) >= minTSlices){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "oneExp tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), oneExp);
    } 
  }// next t min

  //find the best oneExp fit
  FitDescriptor bestOneExp = fits.getBestFit(*fitComp);
  //  cout << "found the best oneExp fit" << endl;
  JackFit bestOneExpFit = fits.getFit( bestOneExp );



  //since timeslice data will be ordered
  //tmin is given by the lowest true in active data
  int tminOneExp;
  vector<bool> active = bestOneExp.activeData;
  for(int i = 0; i < active.size(); i++){ if(active[i]){tminOneExp = int(all_tslices[i]); break;}; }

  //reinstate the data
  for(int i = 0; i < tslices.size(); i++){(fits.getEnsemData()).showDatumByX( tslices[i] );}

  //=======================
  //twoExp fits
  //=======================
  CorrTwoExp* tmp2 = new CorrTwoExp();
  Handle<FitFunction> twoExp(tmp2);

  twoExp->setDefaultParValue("mass_0", bestOneExpFit.getAvgFitParValue("mass_0") ); 
  twoExp->setDefaultParError("mass_0", 5.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_0", bestOneExpFit.getAvgFitParValue("amp_0") ); 
  twoExp->setDefaultParError("amp_0", 5.0*bestOneExpFit.getAvgFitParError("amp_0"));
  
  twoExp->setDefaultParValue("mass_1", bestOneExpFit.getAvgFitParValue("mass_0") + 0.5);
  twoExp->setDefaultParError("mass_1", 10.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_1", bestOneExpFit.getAvgFitParValue("amp_0") );   //assume the excited amplitude the same 
  twoExp->setDefaultParError("amp_1", 5.0*bestOneExpFit.getAvgFitParError("amp_0")); //big error

  //limit the masses to positive values - TURNED OFF !
  //  twoExp->setParamLowerLimit("mass_0", 0.0);
  //  twoExp->setParamLowerLimit("mass_1", 0.0);

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( ((tmax - int(tslices[i]) + 1.1 ) >= minTSlices ) && (tslices[i] + 0.1 < tminOneExp ) ){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "twoExp tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), twoExp);
    } 
  }// next t min

  //get the best jack fit
  int rank;
  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);

  if(best.fitname == "FAILED"){
    //no jackknife fits found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED";
    fit_summary = "FAILED -  set mass_0 to 1.0";
    //no plot !!!
  }
  /*  else if(toDouble(mean(fits.getFit(best).getJackFitParValue("mass_0"))) <= 0.0)
    {
    //negative mass found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED (negative mass)";
    fit_summary = "FAILED: fit with negative mass_0 - set mass_0 to 1.0";
    //no plot !!!
    }
  */
  else{
    if( (best.ff).operator->() == twoExp.operator->() ){nExp = 2;}else{nExp = 1;};
    JackFit& bestFit = fits.getFit(best);
    
    mass_0 = bestFit.getJackFitParValue("mass_0");
    amp_0 = bestFit.getJackFitParValue("amp_0");
    if(nExp == 2){
      mass_1 = bestFit.getJackFitParValue("mass_1");
      amp_1 = bestFit.getJackFitParValue("amp_1");
    }
    
    chisq = bestFit.getJackChisq();
    nDoF = bestFit.getNDoF();
    best_fit_name = best.fitname;
    
    
    //write out a summary of the fits
    int count = 1;
    map<double, FitDescriptor> list = fits.getFitList(*fitComp);
    stringstream ss; 
    ss << "                         | chisq/nDoF |     Q      |  fitCrit   | " << endl;
    for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
      JackFit& thisFit = fits.getFit(p->second);
      double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
      double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
      ss << setw(25) <<(p->second).fitname << "|";
      ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
      ss << setw(12) << fixed << setprecision(3) << Q <<"|";
      ss << setw(12) << scientific << setprecision(3) << p->first << "|";
      
      if(count == rank){ ss << "*";}else{ ss << " ";}
      
      ss << " m0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("mass_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("mass_0");
      ss << " a0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("amp_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("amp_0");
      
      if( ((p->second).ff).operator->() == twoExp.operator->() ){
	ss << ", m1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("mass_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("mass_1");
	ss << ", a1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("amp_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("amp_1");
      }
      ss << endl;
      count++;
    } 
    fit_summary = ss.str();
    
    //make the plot
    ConstTimesExp expWeight(mass_0 / mass_0 , mass_0);
    
    stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; m=" << fixed << setprecision(4) << toDouble(mean(mass_0)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(mass_0)));
    
    axis_plot = bestFit.makeJackFitPlotAxis(expWeight, 0.0, double(tmax + 5.5), lab.str() );
  }


};


//========================================================================================================================
FitCorrelatorCosh::FitCorrelatorCosh(EnsemData data_, int T_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices)
 : T(T_), fits(data_), fitComp(fitComp_){

  EnsemVectorReal y = (fits.getEnsemData()).getYData();
  vector<double> t_temp = (fits.getEnsemData()).getXData();

  //exclude noisy points
  (fits.getEnsemData()).hideDataAboveYErrRat(noiseRatioCutoff);   // alternatively can exclude all points after the time where the data is first noisy - not implemented

  vector<double> all_tslices = (fits.getEnsemData()).getAllXData();
  int tlow = all_tslices[0];
  if(all_tslices.size() < minTSlices){ cerr << "** not enough timeslices to have minTSlices = " << minTSlices << endl; exit(1); }

  //ensure enough timeslices to do a fit
  if((fits.getEnsemData()).getNData() < minTSlices ){
    (fits.getEnsemData()).hideAll();
    for(int t = tlow; t <= tlow + minTSlices; t++){
      (fits.getEnsemData()).showDatumByX(t);
    }
  }

  //find the largest t-value being considered 
  vector<double> tslices = (fits.getEnsemData()).getXData();
  int tmax = int( *max_element(tslices.begin(), tslices.end() ) );
    
  //=======================
  //oneCosh fits
  //=======================
  CorrOneCosh* tmp1 = new CorrOneCosh(T);
  Handle<FitFunction> oneExp(tmp1);

  oneExp->setDefaultParValue("mass_0", 0.2); //totally arbitrary - could use meff to guess
  oneExp->setDefaultParError("mass_0", 0.2); //totally arbitrary

  oneExp->setDefaultParValue("amp_0", 1.0); //totally arbitrary
  oneExp->setDefaultParError("amp_0", 1.0); //totally arbitrary

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( (tmax - int(tslices[i]) + 1.1 ) >= minTSlices){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "oneCosh tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), oneExp);
    } 
  }// next t min

  //find the best oneExp fit
  FitDescriptor bestOneExp = fits.getBestFit(*fitComp);
  JackFit bestOneExpFit = fits.getFit( bestOneExp );
  
  //since timeslice data will be ordered
  //tmin is given by the lowest true in active data
  int tminOneExp;
  vector<bool> active = bestOneExp.activeData;
  for(int i = 0; i < active.size(); i++){ if(active[i]){tminOneExp = int(all_tslices[i]); break;}; }

  //reinstate the data
  for(int i = 0; i < tslices.size(); i++){(fits.getEnsemData()).showDatumByX( tslices[i] );}

  //=======================
  //twoExp fits
  //=======================
  CorrTwoCosh* tmp2 = new CorrTwoCosh(T);
  Handle<FitFunction> twoExp(tmp2);

  twoExp->setDefaultParValue("mass_0", bestOneExpFit.getAvgFitParValue("mass_0") ); 
  twoExp->setDefaultParError("mass_0", 5.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_0", bestOneExpFit.getAvgFitParValue("amp_0") ); 
  twoExp->setDefaultParError("amp_0", 5.0*bestOneExpFit.getAvgFitParError("amp_0"));
  
  twoExp->setDefaultParValue("mass_1", bestOneExpFit.getAvgFitParValue("mass_0") + 0.5   );
  twoExp->setDefaultParError("mass_1", 10.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_1", bestOneExpFit.getAvgFitParValue("amp_0") );   //assume the excited amplitude the same 
  twoExp->setDefaultParError("amp_1", 5.0*bestOneExpFit.getAvgFitParError("amp_0")); //big error

  //limit the masses to positive values
  twoExp->setParamLowerLimit("mass_0", 0.0);
  twoExp->setParamLowerLimit("mass_1", 0.0);

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( ((tmax - int(tslices[i]) + 1.1 ) >= minTSlices ) && (tslices[i] + 0.1 < tminOneExp ) ){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "twoCosh tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), twoExp);
    } 
  }// next t min

   //get the best jack fit
  int rank;
  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);

  if(best.fitname == "FAILED"){
    //no jackknife fits found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;   chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED";
    fit_summary = "FAILED -  set mass_0 to 1.0";
    //no plot !!!
  }
  /*  else if(toDouble(mean(fits.getFit(best).getJackFitParValue("mass_0"))) <= 0.0)
    {
    //negative mass found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;
    chisq = 1.0e10;
    nDoF = 1;
    nExp = 1;
    best_fit_name = "FAILED (negative mass)";
    fit_summary = "FAILED: fit with negative mass_0 - set mass_0 to 1.0";
    //no plot !!!
    }
  */
  else{
  
    if( (best.ff).operator->() == twoExp.operator->() ){nExp = 2;}else{nExp = 1;};
    JackFit& bestFit = fits.getFit(best);
    
    mass_0 = bestFit.getJackFitParValue("mass_0");
    amp_0 = bestFit.getJackFitParValue("amp_0");
    if(nExp == 2){
      mass_1 = bestFit.getJackFitParValue("mass_1");
      amp_1 = bestFit.getJackFitParValue("amp_1");
    }
    
    chisq = bestFit.getJackChisq();
    nDoF = bestFit.getNDoF();
    best_fit_name = best.fitname;
  
    //write out a summary of the fits
    int count = 1;
    map<double, FitDescriptor> list = fits.getFitList(*fitComp);
    stringstream ss; 
    ss << "                         | chisq/nDoF |     Q      |  fitCrit   | " << endl;
    for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
      JackFit& thisFit = fits.getFit(p->second);
      double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
      double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
      ss << setw(25) <<(p->second).fitname << "|";
      ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
      ss << setw(12) << fixed << setprecision(3) << Q <<"|";
      ss << setw(12) << scientific << setprecision(3) << p->first << "|";
      
      if(count == rank){ ss << "*";}else{ ss << " ";}
      
      ss << " m0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("mass_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("mass_0");
      ss << " a0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("amp_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("amp_0");
      
      if( ((p->second).ff).operator->() == twoExp.operator->() ){
        ss << ", m1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("mass_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("mass_1");
        ss << ", a1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("amp_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("amp_1");
      }
      ss << endl;
      count++;
    } 
    fit_summary = ss.str();
    
    

    //make the plot  -- this may need a rethink !!
    //better to divide the whole cosh
    ConstDivideCosh weight(mass_0 / mass_0 , mass_0, T);//Real(0.0)*mass_0);
    
    stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; m=" << fixed << setprecision(4) << toDouble(mean(mass_0)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(mass_0
)));
    
    axis_plot = bestFit.makeJackFitPlotAxis(weight, 0.0, double(tmax + 5.5), lab.str() );
  }


};
///========================================================================================================================
FitCorrelatorExpAndConst::FitCorrelatorExpAndConst(EnsemData data_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices) : fits(data_), fitComp(fitComp_){

  EnsemVectorReal y = (fits.getEnsemData()).getYData();
  vector<double> t_temp = (fits.getEnsemData()).getXData();

  //exclude noisy points
  (fits.getEnsemData()).hideDataAboveYErrRat(noiseRatioCutoff); 
 
  vector<double> all_tslices = (fits.getEnsemData()).getAllXData();
  int tlow = all_tslices[0];

  if(all_tslices.size() < minTSlices){ cerr << "** not enough timeslices to have minTSlices = " << minTSlices << endl; exit(1); }

  //ensure enough timeslices to do a fit
  if((fits.getEnsemData()).getNData() < minTSlices ){
    (fits.getEnsemData()).hideAll();
    for(int t = tlow; t <= tlow + minTSlices; t++){
      (fits.getEnsemData()).showDatumByX(t);
    }
  }

  //find the largest t-value being considered 
  vector<double> tslices = (fits.getEnsemData()).getXData();
  int tmax = int( *max_element(tslices.begin(), tslices.end() ) );

  //  cout << "about to start oneExp fits" << endl;
    
  //=======================
  //oneExp fits
  //=======================
  CorrOneExpAndConst* tmp1 = new CorrOneExpAndConst();
  Handle<FitFunction> oneExp(tmp1);

  //use the value of the correlator at the largest time as the start value for the correlator
  EnsemVectorReal ytmp = (fits.getEnsemData()).getYData();
  int size = ytmp.numElem(); 
  double c = toDouble(mean( peekObs(y, size - 1) ));
  double c_err = toDouble(sqrt(variance( peekObs(y, size - 1) )));

  oneExp->setDefaultParValue("const", c); 
  oneExp->setDefaultParError("const", 3.0*c_err); 

  oneExp->setDefaultParValue("mass_0", 0.2); //totally arbitrary - could use meff to guess
  oneExp->setDefaultParError("mass_0", 0.2); //totally arbitrary

  oneExp->setDefaultParValue("amp_0", 1.0); //totally arbitrary
  oneExp->setDefaultParError("amp_0", 1.0); //totally arbitrary

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( (tmax - int(tslices[i]) + 1.1 ) >= minTSlices){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "oneExpAndConst tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), oneExp);
    } 
  }// next t min

  //find the best oneExp fit
  FitDescriptor bestOneExp = fits.getBestFit(*fitComp);
  //  cout << "found the best oneExp fit" << endl;
  JackFit bestOneExpFit = fits.getFit( bestOneExp );



  //since timeslice data will be ordered
  //tmin is given by the lowest true in active data
  int tminOneExp;
  vector<bool> active = bestOneExp.activeData;
  for(int i = 0; i < active.size(); i++){ if(active[i]){tminOneExp = int(all_tslices[i]); break;}; }

  //reinstate the data
  for(int i = 0; i < tslices.size(); i++){(fits.getEnsemData()).showDatumByX( tslices[i] );}

  //=======================
  //twoExp fits
  //=======================
  CorrTwoExpAndConst* tmp2 = new CorrTwoExpAndConst();
  Handle<FitFunction> twoExp(tmp2);

  twoExp->setDefaultParValue("const", bestOneExpFit.getAvgFitParValue("const") ); 
  twoExp->setDefaultParError("const", 5.0*bestOneExpFit.getAvgFitParError("const"));

  twoExp->setDefaultParValue("mass_0", bestOneExpFit.getAvgFitParValue("mass_0") ); 
  twoExp->setDefaultParError("mass_0", 5.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_0", bestOneExpFit.getAvgFitParValue("amp_0") ); 
  twoExp->setDefaultParError("amp_0", 5.0*bestOneExpFit.getAvgFitParError("amp_0"));
  
  twoExp->setDefaultParValue("mass_1", bestOneExpFit.getAvgFitParValue("mass_0") + 0.5);
  twoExp->setDefaultParError("mass_1", 10.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_1", bestOneExpFit.getAvgFitParValue("amp_0") );   //assume the excited amplitude the same 
  twoExp->setDefaultParError("amp_1", 5.0*bestOneExpFit.getAvgFitParError("amp_0")); //big error

  //limit the masses to positive values - TURNED OFF !
  //  twoExp->setParamLowerLimit("mass_0", 0.0);
  //  twoExp->setParamLowerLimit("mass_1", 0.0);

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( ((tmax - int(tslices[i]) + 1.1 ) >= minTSlices ) && (tslices[i] + 0.1 < tminOneExp ) ){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "twoExpAndConst tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), twoExp);
    } 
  }// next t min

  //get the best jack fit
  int rank;
  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);

  if(best.fitname == "FAILED"){
    //no jackknife fits found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    constant = dum; mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED";
    fit_summary = "FAILED -  set mass_0 to 1.0";
    //no plot !!!
  }
  /*  else if(toDouble(mean(fits.getFit(best).getJackFitParValue("mass_0"))) <= 0.0)
    {
    //negative mass found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED (negative mass)";
    fit_summary = "FAILED: fit with negative mass_0 - set mass_0 to 1.0";
    //no plot !!!
    }
  */
  else{
    if( (best.ff).operator->() == twoExp.operator->() ){nExp = 2;}else{nExp = 1;};
    JackFit& bestFit = fits.getFit(best);
    
    constant =  bestFit.getJackFitParValue("const");
    mass_0 = bestFit.getJackFitParValue("mass_0");
    amp_0 = bestFit.getJackFitParValue("amp_0");
    if(nExp == 2){
      mass_1 = bestFit.getJackFitParValue("mass_1");
      amp_1 = bestFit.getJackFitParValue("amp_1");
    }
    
    chisq = bestFit.getJackChisq();
    nDoF = bestFit.getNDoF();
    best_fit_name = best.fitname;
    
    
    //write out a summary of the fits
    int count = 1;
    map<double, FitDescriptor> list = fits.getFitList(*fitComp);
    stringstream ss; 
    ss << "                               | chisq/nDoF |     Q      |  fitCrit   | " << endl;
    for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
      JackFit& thisFit = fits.getFit(p->second);
      double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
      double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
      ss << setw(25) <<(p->second).fitname << "|";
      ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
      ss << setw(12) << fixed << setprecision(3) << Q <<"|";
      ss << setw(12) << scientific << setprecision(3) << p->first << "|";
      
      if(count == rank){ ss << "*";}else{ ss << " ";}
      //      ss << " c=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("const") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("const");   
      ss << " c=" << setw(8)  << thisFit.getAvgFitParValue("const") << " +/-" <<  setw(8) << thisFit.getAvgFitParError("const");   
      ss << " m0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("mass_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("mass_0");
      ss << " a0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("amp_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("amp_0");
      
      if( ((p->second).ff).operator->() == twoExp.operator->() ){
	ss << ", m1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("mass_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("mass_1");
	ss << ", a1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("amp_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("amp_1");
      }
      ss << endl;
      count++;
    } 
    fit_summary = ss.str();
    
    //make the plot
    ConstTimesExp expWeight(mass_0 / mass_0 , mass_0);
    //One one(mass_0);

    stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; m=" << fixed << setprecision(4) << toDouble(mean(mass_0)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(mass_0)));
    
    axis_plot = bestFit.makeJackFitPlotAxis(expWeight, 0.0, double(tmax + 5.5), lab.str() );
    //axis_plot = bestFit.makeJackFitPlotAxis(one, 0.0, double(tmax + 5.5), lab.str() );
  }


};


///========================================================================================================================
FitCorrelatorCoshAndConst::FitCorrelatorCoshAndConst(EnsemData data_, int T_, Handle<FitComparator> fitComp_, double noiseRatioCutoff, int minTSlices) : fits(data_), fitComp(fitComp_), T(T_) {

  EnsemVectorReal y = (fits.getEnsemData()).getYData();
  vector<double> t_temp = (fits.getEnsemData()).getXData();

  //exclude noisy points
  (fits.getEnsemData()).hideDataAboveYErrRat(noiseRatioCutoff); 
 
  vector<double> all_tslices = (fits.getEnsemData()).getAllXData();
  int tlow = all_tslices[0];

  if(all_tslices.size() < minTSlices){ cerr << "** not enough timeslices to have minTSlices = " << minTSlices << endl; exit(1); }

  //ensure enough timeslices to do a fit
  if((fits.getEnsemData()).getNData() < minTSlices ){
    (fits.getEnsemData()).hideAll();
    for(int t = tlow; t <= tlow + minTSlices; t++){
      (fits.getEnsemData()).showDatumByX(t);
    }
  }

  //find the largest t-value being considered 
  vector<double> tslices = (fits.getEnsemData()).getXData();
  int tmax = int( *max_element(tslices.begin(), tslices.end() ) );

  //  cout << "about to start oneExp fits" << endl;
    
  //=======================
  //oneExp fits
  //=======================
  CorrOneCoshAndConst* tmp1 = new CorrOneCoshAndConst(T);
  Handle<FitFunction> oneExp(tmp1);

  //use the value of the correlator at the largest time as the start value for the correlator
  EnsemVectorReal ytmp = (fits.getEnsemData()).getYData();
  int size = ytmp.numElem(); 
  double c = toDouble(mean( peekObs(y, size - 1) ));
  double c_err = toDouble(sqrt(variance( peekObs(y, size - 1) )));

  oneExp->setDefaultParValue("const", c); 
  oneExp->setDefaultParError("const", 3.0*c_err); 

  oneExp->setDefaultParValue("mass_0", 0.2); //totally arbitrary - could use meff to guess
  oneExp->setDefaultParError("mass_0", 0.2); //totally arbitrary

  oneExp->setDefaultParValue("amp_0", 1.0); //totally arbitrary
  oneExp->setDefaultParError("amp_0", 1.0); //totally arbitrary

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( (tmax - int(tslices[i]) + 1.1 ) >= minTSlices){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "oneCoshAndConst tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), oneExp);
    } 
  }// next t min

  //find the best oneExp fit
  FitDescriptor bestOneExp = fits.getBestFit(*fitComp);
  //  cout << "found the best oneExp fit" << endl;
  JackFit bestOneExpFit = fits.getFit( bestOneExp );



  //since timeslice data will be ordered
  //tmin is given by the lowest true in active data
  int tminOneExp;
  vector<bool> active = bestOneExp.activeData;
  for(int i = 0; i < active.size(); i++){ if(active[i]){tminOneExp = int(all_tslices[i]); break;}; }

  //reinstate the data
  for(int i = 0; i < tslices.size(); i++){(fits.getEnsemData()).showDatumByX( tslices[i] );}

  //=======================
  //twoExp fits
  //=======================
  CorrTwoCoshAndConst* tmp2 = new CorrTwoCoshAndConst(T);
  Handle<FitFunction> twoExp(tmp2);

  twoExp->setDefaultParValue("const", bestOneExpFit.getAvgFitParValue("const") ); 
  twoExp->setDefaultParError("const", 5.0*bestOneExpFit.getAvgFitParError("const"));

  twoExp->setDefaultParValue("mass_0", bestOneExpFit.getAvgFitParValue("mass_0") ); 
  twoExp->setDefaultParError("mass_0", 5.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_0", bestOneExpFit.getAvgFitParValue("amp_0") ); 
  twoExp->setDefaultParError("amp_0", 5.0*bestOneExpFit.getAvgFitParError("amp_0"));
  
  twoExp->setDefaultParValue("mass_1", bestOneExpFit.getAvgFitParValue("mass_0") + 0.5);
  twoExp->setDefaultParError("mass_1", 10.0*bestOneExpFit.getAvgFitParError("mass_0"));

  twoExp->setDefaultParValue("amp_1", bestOneExpFit.getAvgFitParValue("amp_0") );   //assume the excited amplitude the same 
  twoExp->setDefaultParError("amp_1", 5.0*bestOneExpFit.getAvgFitParError("amp_0")); //big error

  //limit the masses to positive values - TURNED OFF !
  //  twoExp->setParamLowerLimit("mass_0", 0.0);
  //  twoExp->setParamLowerLimit("mass_1", 0.0);

  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( ((tmax - int(tslices[i]) + 1.1 ) >= minTSlices ) && (tslices[i] + 0.1 < tminOneExp ) ){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
      stringstream s; s << "twoCoshAndConst tmin= " << int(tslices[i]) << " tmax= " << tmax;
      fits.addFit(s.str(), twoExp);
    } 
  }// next t min

  //get the best jack fit
  int rank;
  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);

  if(best.fitname == "FAILED"){
    //no jackknife fits found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    constant = dum; mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED";
    fit_summary = "FAILED -  set mass_0 to 1.0";
    //no plot !!!
  }
  /*  else if(toDouble(mean(fits.getFit(best).getJackFitParValue("mass_0"))) <= 0.0)
    {
    //negative mass found - default to the cutoff - HORRIBLE !!!
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(1.0);
    mass_0 = dum;    chisq = 1.0e10;    nDoF = 1;    nExp = 1;
    best_fit_name = "FAILED (negative mass)";
    fit_summary = "FAILED: fit with negative mass_0 - set mass_0 to 1.0";
    //no plot !!!
    }
  */
  else{
    if( (best.ff).operator->() == twoExp.operator->() ){nExp = 2;}else{nExp = 1;};
    JackFit& bestFit = fits.getFit(best);
    
    constant =  bestFit.getJackFitParValue("const");
    mass_0 = bestFit.getJackFitParValue("mass_0");
    amp_0 = bestFit.getJackFitParValue("amp_0");
    if(nExp == 2){
      mass_1 = bestFit.getJackFitParValue("mass_1");
      amp_1 = bestFit.getJackFitParValue("amp_1");
    }
    
    chisq = bestFit.getJackChisq();
    nDoF = bestFit.getNDoF();
    best_fit_name = best.fitname;
    
    
    //write out a summary of the fits
    int count = 1;
    map<double, FitDescriptor> list = fits.getFitList(*fitComp);
    stringstream ss; 
    ss << "                               | chisq/nDoF |     Q      |  fitCrit   | " << endl;
    for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
      JackFit& thisFit = fits.getFit(p->second);
      double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
      double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
      ss << setw(25) <<(p->second).fitname << "|";
      ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
      ss << setw(12) << fixed << setprecision(3) << Q <<"|";
      ss << setw(12) << scientific << setprecision(3) << p->first << "|";
      
      if(count == rank){ ss << "*";}else{ ss << " ";}
      //ss << " c=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("const") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("const");   
      ss << " c=" << setw(8) << thisFit.getAvgFitParValue("const") << " +/-" <<  setw(8) << thisFit.getAvgFitParError("const");   
      ss << " m0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("mass_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("mass_0");
      ss << " a0=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("amp_0") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("amp_0");
      
      if( ((p->second).ff).operator->() == twoExp.operator->() ){
	ss << ", m1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("mass_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("mass_1");
	ss << ", a1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("amp_1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("amp_1");
      }
      ss << endl;
      count++;
    } 
    fit_summary = ss.str();
    
    //make the plot
    ConstTimesExp expWeight(mass_0 / mass_0 , mass_0);
    //One one(mass_0);

    stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; m=" << fixed << setprecision(4) << toDouble(mean(mass_0)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(mass_0)));
    
    axis_plot = bestFit.makeJackFitPlotAxis(expWeight, 0.0, double(tmax + 5.5), lab.str() );
    //axis_plot = bestFit.makeJackFitPlotAxis(one, 0.0, double(tmax + 5.5), lab.str() );
  }


};
