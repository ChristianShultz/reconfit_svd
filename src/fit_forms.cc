#include "fit_forms.h"

// JJD - attempt ot improve the fitting of Z(t)
// this method seems to reduce the number of completely spurious fits
// it still can't do anything if the variation over time is incompatible with being a constant
// suggests that matching timeslice-timeslice of eigenvectors is still imperfect






FitZ::FitZ(EnsemData data_, int t0, Handle<FitComparator> fitComp_, int minTSlices) : fits(data_), fitComp(fitComp_){

  double t0plus1 = double(t0 + 1);
  EnsemReal start = data_.getYUsingNearestX( t0plus1 );
  if(fabs( t0plus1 - double( t0 + 1 ) ) > 0.1){ cerr << "somethings gone wrong here - failing" << endl; }

  vector<double> tslices = data_.getAllXData(); 
  int tend = int( *max_element(tslices.begin(), tslices.end() ) );

  //  cout << "tend = " << tend << endl;

  // try constant fits, expanding "outwards" from t0+1
  Const* tmp1 = new Const; Handle<FitFunction> constant(tmp1);

  constant->setDefaultParValue("a", toDouble( mean(start) ) );
  constant->setDefaultParValue("a", toDouble( sqrt(variance(start)) ) );
  
  // need more intelligent handling if t0+1 is a baddy - NOT CURRENTLY PRESENT !!!!!!!!!!!

  for(int tmax = t0+1; tmax <= tend; tmax++){
    for(int tmin = tmax-3; tmin >= 1; tmin--){
      (fits.getEnsemData()).showAll();
      (fits.getEnsemData()).hideDatumByX(t0);
      (fits.getEnsemData()).hideDataAboveX(tmax +0.1);
      (fits.getEnsemData()).hideDataBelowX(tmin -0.1);
      stringstream s; s << "constant tmin= " << tmin << " tmax= " << tmax;

      if( (fits.getEnsemData()).getNData() >= minTSlices ){
	fits.addFit(s.str(), constant);
      }

    }
  }
  
  FitDescriptor best_const = fits.getBestFit(*fitComp);
  JackFit best_const_fit = fits.getFit( best_const );
  
  //need to find the tmin, tmax for the best constant fit
  int tmin_best_const, tmax_best_const;
  vector<bool> active = best_const.activeData;
  for(int i = 0; i < active.size(); i++){ if(active[i]){tmin_best_const = int(tslices[i]); break;}; } 
  for(int i = 0; i < active.size(); i++){ if(active[i]){tmax_best_const = int(tslices[i]);}; }

  //try adding an exponential
  ConstPlusExp* tmp2 = new ConstPlusExp(); Handle<FitFunction> constant_plus_exp(tmp2);

  constant_plus_exp->setDefaultParValue("a", best_const_fit.getAvgFitParValue("a") ); 
  constant_plus_exp->setDefaultParError("a", 5.0*best_const_fit.getAvgFitParError("a"));
  constant_plus_exp->setDefaultParValue("b", 0.0);
  constant_plus_exp->setDefaultParError("b", 10.0*best_const_fit.getAvgFitParError("a") );
  constant_plus_exp->setDefaultParValue("exp_mass", 0.1);
  constant_plus_exp->setDefaultParError("exp_mass", 0.1);
  //limit the exponential mass to positive values
  constant_plus_exp->setParamLowerLimit("exp_mass", 0.0);

  //keep the best tmax from the constant fit - don't change that
  //walk backwards from the constant fit tmin
  for(int tmin = tmin_best_const; tmin >= 1; tmin--){
    (fits.getEnsemData()).showAll();
    (fits.getEnsemData()).hideDatumByX(t0);
    (fits.getEnsemData()).hideDataAboveX(tmax_best_const + 0.1 );
    (fits.getEnsemData()).hideDataBelowX(tmin - 0.1 );
    stringstream s; s << "constant_plus_exp tmin= " << tmin << " tmax= " << tmax_best_const;

    if( (fits.getEnsemData()).getNData() >= minTSlices ){
      fits.addFit(s.str(), constant_plus_exp);
    }
  }

  //get the best jack fit
  int rank;  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);
  if( best.fitname == "FAILED" ){
    //no jackknife fits found - zero out the Z
    EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(0.0);
    Z = dum;
    chisq = 1.0e10;
    nDoF = 1;
    best_fit_name = "FAILED";
    fit_summary = "FAILED -  set Z to 0.0";
    //no plot !!!
  }
  else{
    if( (best.ff).operator->() == constant_plus_exp.operator->() ){fit_type = "constant_plus_exp";}
    else{fit_type = "constant";};
    JackFit& bestFit = fits.getFit(best);
    
    Z = bestFit.getJackFitParValue("a");
    if(fit_type != "constant"){
      Z2 = bestFit.getJackFitParValue("b");
      exp_mass = bestFit.getJackFitParValue("exp_mass");
    }
    
    chisq = bestFit.getJackChisq();
    nDoF = bestFit.getNDoF();
    best_fit_name = best.fitname;
    
    //write out a summary of the fits

    // MAKE THIS PRETTIER WITH printf

    int count = 1;
    map<double, FitDescriptor> list = fits.getFitList(*fitComp);
    stringstream ss; 
    ss << "                                   | chisq/nDoF |     Q      |  fitCrit   | " << endl;
    for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
      JackFit& thisFit = fits.getFit(p->second);
      double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
      double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
      ss << setw(35) <<(p->second).fitname << "|";
      ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
      ss << setw(12) << fixed << setprecision(3) << Q <<"|";
      ss << setw(12) << scientific << setprecision(3) << p->first << "|";
      
      if(count == rank){ ss << "*";}else{ ss << " ";}
      ss << " Z=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("a") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("a");
      if( ((p->second).ff).operator->() == constant_plus_exp.operator->() ){
	ss << ", Z'=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("b") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("b");
	ss << ", m=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("exp_mass") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("exp_mass");
      }
      ss << endl;count++;
    } 
    //    cout << ss.str();
    fit_summary = ss.str();
    
    //make the plot
    stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; Z=" << fixed << setprecision(4) << toDouble(mean(Z)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(Z)));
    axis_plot = bestFit.makeJackFitPlotAxis(0.0, double(tend + 5.5), lab.str() );    
  }
  
};






// 'original' Z(t) fitting 

// FitZ::FitZ(EnsemData data_, int t0, Handle<FitComparator> fitComp_, int minTSlices) : fits(data_), fitComp(fitComp_){

//   (fits.getEnsemData()).hideDatumByX(t0);

//   // DON'T EXCLUDE ANY DATA BY NOISE ?


//   //exclude noisy points - but not if they are consistent with zero
//   EnsemVectorReal y = (fits.getEnsemData()).getYData();
//   vector<double> x = (fits.getEnsemData()).getXData();
//   vector<double> keep;

//   // set an ABSOLUTE noise scale
//   // exclude points with noise far above this ?


//   /*  int t0_index = find_nearest(x, t0);
//   vector<int> above = find_above(x, t0);
  
//   for(int i = 0; i < above.size(); i++){
//     if( ( x[i] > t0 ) && (x[i] << 
//   EnsemReal yy = peekObs(y, t0_index + 1);
//   double noise_scale = toDouble
//   */

//   /*
//   if(y.numElem() > 2){
//     //set a scale using the first few points
//     EnsemReal scale = ( peekObs(y, 1) + peekObs(y, 2) ) / Real(2.0) ;
//     //if this is consistent with zero, the fit may be too
//     if(abs( toDouble(mean(scale)) ) < 5.0 * toDouble( sqrt(variance(scale)))){
//       //      Zfitlog << "  may be consistent with zero" << endl;
//       for(int i = 0; i < y.numElem(); i++){
// 	EnsemReal yy = peekObs(y, i);
// 	if( //abs(toDouble( mean(yy) )) < 1.5 * toDouble( sqrt(variance(yy)) ) &&
// 	    ( toDouble(sqrt(variance(yy))) < 7.0 * toDouble(sqrt(variance(scale))) )    ){ keep.push_back(x[i]); } 
//       }
//     }
//   }

//   (fits.getEnsemData()).hideDataAboveYErrRat(noiseRatioCutoff); 


//   for(int i = 0; i < keep.size(); i++){
//     (fits.getEnsemData()).showDatumByX(keep[i]);
//   }
//   */

//   //look for jumps as we go through t0
//   // if there's a large jump then hide below t0
//   y = (fits.getEnsemData()).getAllYData();
//   x.clear();  x = data_.getAllXData();
//   int t0pos = find_exact(x, t0);
//   EnsemReal y_diff = peekObs(y, t0pos - 1) - peekObs(y, t0pos + 1);;
//   if( toDouble(mean(y_diff)) > 5.0 * toDouble(sqrt(variance(y_diff)))){
//     //Zfitlog << "  large jump at t0 - will fits only above t0" << endl;
//     (fits.getEnsemData()).hideDataBelowX(t0);
//   }



//   /*
//   //keep increasing the error cutoff until there are some data points
//   int nData = (fits.getEnsemData()).getNData();
//   double rat = noiseRatioCutoff;
//   while( nData < minTSlices ){
//     //Zfitlog << "  increasing the ratio by 25% to capture some points" << endl;
//     rat *= 1.25;
//     y = data_.getYData();
//     x.clear(); x = data_.getXData();
//     for(int i = 0; i < x.size(); i++){
//       if(( abs( toDouble( sqrt(variance(peekObs(y,i))) / mean(peekObs(y,i)) )) < rat) && !( (x[i] > t0 - 0.1)&&(x[i] < t0 + 0.1) ) ){
// 	(fits.getEnsemData()).showDatumByX(x[i]);
//       }
//     }
//     nData = (fits.getEnsemData()).getNData();
//   }
//   */

//   //find the largest t-value being considered 
//   vector<double> tslices = (fits.getEnsemData()).getXData();
//   int tmax = int( *max_element(tslices.begin(), tslices.end() ) );
  
//   vector<double> all_tslices = (fits.getEnsemData()).getAllXData();
//   //cout << "tmax is " << tmax << endl;

//   //*******************
//   // constant fits
//   //*******************
//   Const* tmp1 = new Const;
//   Handle<FitFunction> constant(tmp1);
 
//   EnsemReal start = peekObs( (fits.getEnsemData()).getYData() , (fits.getEnsemData()).getNData() - 1 );

//   constant->setDefaultParValue("a", toDouble( mean(start) ) );
//   constant->setDefaultParValue("a", toDouble( sqrt(variance(start)) ) );
  
//   //loop over tmins
//   for(int i = 0; i < tslices.size(); i++){
//     if( (tmax - int(tslices[i]) + 1.1 ) >= minTSlices){ 
//       (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
//       stringstream s; s << "constant tmin= " << int(tslices[i]) << " tmax= " << tmax;
//       fits.addFit(s.str(), constant);
//     } 
//   }

//   //find the best oneExp fit
//   FitDescriptor bestConst = fits.getBestFit(*fitComp);

//   JackFit bestConstFit = fits.getFit( bestConst );
  
//   //since timeslice data will be ordered, tmin is given by the lowest true in active data
//   int tminConst;
//   vector<bool> active = bestConst.activeData;
//   for(int i = 0; i < active.size(); i++){ if(active[i]){tminConst = int(all_tslices[i]); break;}; }

//   //reinstate the data
//   for(int i = 0; i < tslices.size(); i++){(fits.getEnsemData()).showDatumByX( tslices[i] );}

//   //=======================
//   // const_plus_exp fits
//   //=======================
//   ConstPlusExp* tmp2 = new ConstPlusExp();
//   Handle<FitFunction> constPlusExp(tmp2);

//   constPlusExp->setDefaultParValue("a", bestConstFit.getAvgFitParValue("a") ); 
//   constPlusExp->setDefaultParError("a", 5.0*bestConstFit.getAvgFitParError("a"));
  
//   constPlusExp->setDefaultParValue("b", 0.0);
//   constPlusExp->setDefaultParError("b", 10.0*bestConstFit.getAvgFitParError("a") );
//   constPlusExp->setDefaultParValue("exp_mass", 0.1);
//   constPlusExp->setDefaultParError("exp_mass", 0.1);

//   //limit the exponential mass to positive values
//   constPlusExp->setParamLowerLimit("exp_mass", 0.0);

//   //loop over tmin values
//   // const plus exp only allowed for tmin < t0
//   for(int i = 0; i < tslices.size(); i++){
//     if( ((tmax - int(tslices[i]) + 1.1 ) >= minTSlices ) && (tslices[i] + 0.1 < tminConst ) && (tslices[i] + 0.1 < t0) ){ 
//       (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   //timeslices are ordered, so this is safe
      
//       stringstream s; s << "constant_plus_exp tmin= " << int(tslices[i]) << " tmax= " << tmax;
//       fits.addFit(s.str(), constPlusExp);
//     } 
//   }// next t min



//   //get the best jack fit
//   int rank;  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);
//   if( best.fitname == "FAILED" ){
//     //no jackknife fits found - zero out the Z
//     EnsemReal dum; dum.resize( (fits.getEnsemData()).getNBins() ); dum = Real(0.0);
//     Z = dum;
//     chisq = 1.0e10;
//     nDoF = 1;
//     best_fit_name = "FAILED";
//     fit_summary = "FAILED -  set Z to 0.0";
//     //no plot !!!
//   }
//   else{

//     if( (best.ff).operator->() == constPlusExp.operator->() ){fit_type = "constant_plus_exp";}
//     else{fit_type = "constant";};
//     JackFit& bestFit = fits.getFit(best);
    
//     Z = bestFit.getJackFitParValue("a");
//     if(fit_type != "constant"){
//       Z2 = bestFit.getJackFitParValue("b");
//       exp_mass = bestFit.getJackFitParValue("exp_mass");
//     }
    
//     chisq = bestFit.getJackChisq();
//     nDoF = bestFit.getNDoF();
//     best_fit_name = best.fitname;
    
//     //write out a summary of the fits

//     // MAKE THIS PRETTIER WITH printf

//     int count = 1;
//     map<double, FitDescriptor> list = fits.getFitList(*fitComp);
//     stringstream ss; 
//     ss << "                                   | chisq/nDoF |     Q      |  fitCrit   | " << endl;
//     for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
//       JackFit& thisFit = fits.getFit(p->second);
//       double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
//       double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
//       ss << setw(35) <<(p->second).fitname << "|";
//       ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
//       ss << setw(12) << fixed << setprecision(3) << Q <<"|";
//       ss << setw(12) << scientific << setprecision(3) << p->first << "|";
      
//       if(count == rank){ ss << "*";}else{ ss << " ";}
//       ss << " Z=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("a") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("a");
//       if( ((p->second).ff).operator->() == constPlusExp.operator->() ){
// 	ss << ", Z'=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("b") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("b");
// 	ss << ", m=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("exp_mass") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("exp_mass");
//       }
//       ss << endl;count++;
//     } 
//     fit_summary = ss.str();
    
//     //make the plot
//     stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
//     lab << "; Z=" << fixed << setprecision(4) << toDouble(mean(Z)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(Z)));
//     axis_plot = bestFit.makeJackFitPlotAxis(0.0, double(tmax + 5.5), lab.str() );
    
//   }
  
// };



void FitZ::saveFitPlot(string filename){
  ofstream out; out.open(filename.c_str());
  out << axis_plot;
  out.close();
}

//***********************************************************************

FitVersusT0::FitVersusT0(EnsemData data_, Handle<FitComparator> fitComp_) : fits(data_), fitComp(fitComp_){
	
  vector<double> tslices = (fits.getEnsemData()).getXData();
  int tmax = *tslices.rbegin();
  //  cout << "tmax = " << tmax << endl;
  vector<double> all_tslices = (fits.getEnsemData()).getAllXData();
  
  //*******************
  // constant fits
  //*******************
  Const* tmp1 = new Const;
  Handle<FitFunction> constant(tmp1);
  
  EnsemReal start = peekObs( (fits.getEnsemData()).getYData() , (fits.getEnsemData()).getNData() - 1 );

  constant->setDefaultParValue("a", toDouble( mean(start) ) );
  constant->setDefaultParValue("a", toDouble( sqrt(variance(start)) ) );
  
  //loop over tmins
  for(int i = 0; i < tslices.size(); i++){
    if( (tmax - int(tslices[i]) + 1.1 ) >= 2){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 );   
      
      stringstream s; s << "constant tmin= " << int(tslices[i]) ;
      fits.addFit(s.str(), constant);
    } 
  }
  
  //find the best constant fit
  FitDescriptor bestConst = fits.getBestFit(*fitComp);
  JackFit bestConstFit = fits.getFit( bestConst );
  
  int tminConst;
  vector<bool> active = bestConst.activeData;
  for(int i = 0; i < active.size(); i++)
    { if(active[i]){tminConst = int(all_tslices[i]); break;}; }
  
  //  cout << "did the constant fits" << endl;
  // cout << "tminConst = " << tminConst << endl;

  //reinstate the data
  for(int i = 0; i < tslices.size(); i++)
    {(fits.getEnsemData()).showDatumByX( tslices[i] );}
  
  //cout << "turned all the data back on" << endl;

  //=======================
  // const_plus_exp fits
  //=======================
  ConstPlusExp* tmp2 = new ConstPlusExp();
  Handle<FitFunction> constPlusExp(tmp2);
  
  constPlusExp->setDefaultParValue("a", bestConstFit.getAvgFitParValue("a") ); 
  constPlusExp->setDefaultParError("a", 5.0*bestConstFit.getAvgFitParError("a"));
  
  constPlusExp->setDefaultParValue("b", 0.0);
  constPlusExp->setDefaultParError("b", 10.0*bestConstFit.getAvgFitParError("a") );
  constPlusExp->setDefaultParValue("exp_mass", 0.1);
  constPlusExp->setDefaultParError("exp_mass", 0.1);
  
  //limit the exponential mass to positive values
  constPlusExp->setParamLowerLimit("exp_mass", 0.0);
  
  //loop over tmin values
  for(int i = 0; i < tslices.size(); i++){
    if( ((tmax - int(tslices[i]) + 1.1 ) >= 4 ) && (tslices[i] + 0.1 < tminConst ) ){ 
      (fits.getEnsemData()).hideDataBelowX( tslices[i] - 0.1 ); 
      
      stringstream s; s << "constant_plus_exp tmin= " << int(tslices[i]);
      fits.addFit(s.str(), constPlusExp);
    } 
  }// next t min
  
  //get the best jack fit
  int rank;  FitDescriptor best = fits.getBestJackFit(*fitComp, rank);
  
  //cout << "found the best fit" << endl;

  if( (best.ff).operator->() == constPlusExp.operator->() ){fit_type = "constant_plus_exp";}
  else{fit_type = "constant";};
  JackFit& bestFit = fits.getFit(best);
  
  //cout << "got a ref to the best fit" << endl;

  a = bestFit.getJackFitParValue("a");
  
  //  cout << "extracted the fit value a= " << toDouble(mean(a)) << endl;

  chisq = bestFit.getJackChisq();
  nDoF = bestFit.getNDoF();
  best_fit_name = best.fitname;
  


  //write out a summary of the fits
  int count = 1;
  map<double, FitDescriptor> list = fits.getFitList(*fitComp);
  stringstream ss; 
  ss << "                                   | chisq/nDoF |     Q      |  fitCrit   | " << endl;
  for( map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++){
    JackFit& thisFit = fits.getFit(p->second);
    double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
    double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
    ss << setw(35) <<(p->second).fitname << "|";
    ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
    ss << setw(12) << fixed << setprecision(3) << Q <<"|";
    ss << setw(12) << scientific << setprecision(3) << p->first << "|";
    
    if(count == rank){ ss << "*";}else{ ss << " ";}
    ss << " a=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("a") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("a");
    if( ((p->second).ff).operator->() == constPlusExp.operator->() ){
      ss << ", b=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("b") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("b");
      ss << ", exp_mass=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("exp_mass") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("exp_mass");
    }
    ss << endl;count++;
  } 
  fit_summary = ss.str();


  //make the plot
  stringstream lab; lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
  lab << "; a=" << fixed << setprecision(4) << toDouble(mean(a)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(a)));
  
  axis_plot = bestFit.makeJackFitPlotAxis(0.0, double(tmax + 1.5), lab.str() );


};


void FitVersusT0::saveFitPlot(string filename){
  ofstream out; out.open(filename.c_str());
  out << axis_plot;
  out.close();
}

