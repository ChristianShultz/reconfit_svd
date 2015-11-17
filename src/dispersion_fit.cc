#include "multi_ensemble_fitter.h"

#include "ensem/ensem.h"
#include "adat/handle.h"
#include "jackFitter/jackknife_fitter.h"

#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;
using namespace ENSEM;
using namespace ADAT;

class DispnP2 : public FitFunction {
public:
  DispnP2() : FitFunction(2){
    setParName(0, "m");
    setParName(1, "xi");
  }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "dispnP2";}
};

double DispnP2::operator()(const vector<double>& pars, double x) const{
  return sqrt( pars[0]*pars[0] + x / (pars[1]*pars[1]) );
};

class DispnP2P4 : public FitFunction {
public:
  DispnP2P4() : FitFunction(3){
    setParName(0, "m");
    setParName(1, "xi");
    setParName(2, "rho");
  }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "dispnP2P4";}
};

double DispnP2P4::operator()(const vector<double>& pars, double x) const{
  return sqrt( pars[0]*pars[0] + x / ( pars[1]*pars[1] ) + pars[2]*x*x );
};

class DispnLatBoson : public FitFunction {
public:
  DispnLatBoson() : FitFunction(2){
    setParName(0, "m");
    setParName(1, "xi");
  }
  inline double operator()(const vector<double>& pars, double x) const;
  string getFitType() const {return "dispnLatBoson";}
};

double DispnLatBoson::operator()(const vector<double>& pars, double x) const{
  double y = 0.5 * sqrt( pars[0]*pars[0] + x / (pars[1]*pars[1]) );
  y =  2.0 * log( y + sqrt(y*y + 1.0) );
  //cout << "m=" << pars[0] << ", xi=" << pars[1] << "; x=" << x << ", E = " << y << endl;
  return y;
};

//##################################################################################################


int main(int argc, char *argv[])
{
  if(argc < 3){ cerr << "usage: <V1> <listfile1> <V2> <listfile2> ... " << endl; exit(1); }
  
  double pi = 2.0*atan2(1.0,0.0);

  map<string, EnsemData> data;

  for(int n = 1; n < argc; n += 2){
    double V; {istringstream val(argv[n]); val >> V;}
    string filen; {istringstream val(argv[n+1]); val >> filen;}
    
    //    cout << "V= " << V << ", listfile = " << filen << endl;

    // listfiles are formatted as
    // int(psq) Epi_ensem_file
    EnsemData ensemble;
    string line, buf; 
    ifstream listfile(filen.c_str());

    if(listfile.is_open()){
      cout << "reading from " << filen << endl;

      while( !(listfile.eof()) ){
	if( !getline(listfile, line) ) break;
	// cout << "  " << line <<" : ";
	vector<string> tokens;
	stringstream ss(line); while(ss >> buf){tokens.push_back(buf);}
	
	int psq; istringstream sss(tokens[0].c_str()); sss >> psq;
	string jackfile = tokens[1]; EnsemReal ensem;
	read(jackfile, ensem);

	// cout << psq << "  " << toDouble(mean(ensem)) << endl;

	ensemble.addDatum(double(psq) * pow( 2.0 * pi / V , 2) , ensem);
      }

      // cout << "finished reading from " << filen << endl;
      listfile.close();
    }
    else
    {
      std::cerr << __func__ << ": error reading " << filen << std::endl;
      exit(1);
    }

    //    cout << "V = " << V << endl;
    //   cout << ensemble << endl << endl;

    stringstream s; s << (n-1)/2 << "_" << V;
    data.insert( pair<string, EnsemData>(s.str(), ensemble) );

  }//next listfile

  //  cout << endl << endl;
  for(map<string, EnsemData>::iterator it = data.begin(); it != data.end(); it++){
    cout << "ensemble = " <<  (*it).first << endl << (*it).second << endl << endl;
  }
  // DATA LOADED

  DispnP2* tmp = new DispnP2();
  Handle<FitFunction> dispnP2(tmp);

  dispnP2->setDefaultParValue("m", 0.07);
  dispnP2->setDefaultParError("m", 0.02);

  dispnP2->setDefaultParValue("xi", 3.5);
  dispnP2->setDefaultParError("xi", 0.3);

  MultiEnsemFitter fitter(data, dispnP2);

  if( ! fitter.fit() ){
    cout << " FIT FAILED " << endl; exit(1);
  }
  
  cout << "DISPERSION P2" << endl;

  cout << "chisq/nDoF = " << fixed << setprecision(2) << fitter.getChisq() << "/" << fitter.getNDoF() << endl;
  cout << "m = " << fixed << setprecision(5) << fitter.getFitParValue("m") << " +/- " << fitter.getFitParError("m") << endl;
  cout << "xi = " << fixed << fitter.getFitParValue("xi") << " +/- " << fitter.getFitParError("xi") << endl;

  stringstream label; label << "chisq/nDoF=" << fixed << setprecision(0) << fitter.getChisq() << "/" << fitter.getNDoF(); 
  label << " xi=" << setprecision(4) << fitter.getFitParValue("xi") << "+/-" << fitter.getFitParError("xi");
  string plot = fitter.makeFitPlotAxis(0, .5, label.str());

  // cout << fitter.getFitReport() << endl;

  //write the plot to a file
  {  
    stringstream ss; ss << "dispnP2_fit.ax"; 
    ofstream out; out.open(ss.str().c_str());
    out << plot;
    out.close();
  }

  //*********************************************************

  DispnP2P4* tmp2 = new DispnP2P4();
  Handle<FitFunction> dispnP2P4(tmp2);

  dispnP2P4->setDefaultParValue("m", fitter.getFitParValue("m") );
  dispnP2P4->setDefaultParError("m", fitter.getFitParValue("m") );

  dispnP2P4->setDefaultParValue("xi", fitter.getFitParValue("xi") );
  dispnP2P4->setDefaultParError("xi", fitter.getFitParError("xi") );

  dispnP2P4->setDefaultParValue("rho", 0.1);
  dispnP2P4->setDefaultParError("rho", 0.1);

  MultiEnsemFitter fitter2(data, dispnP2P4);

  if( ! fitter2.fit() ){
    cout << " FIT FAILED " << endl; exit(1);
  }
  
  cout << endl << "********************************" << endl <<"DISPERSION P2P4" << endl;

  cout << "chisq/nDoF = " << fixed << setprecision(2) << fitter2.getChisq() << "/" << fitter2.getNDoF() << endl;
  cout << "m = " << fixed << setprecision(5) << fitter2.getFitParValue("m") << " +/- " << fitter2.getFitParError("m") << endl;
  cout << "xi = " << fixed << fitter2.getFitParValue("xi") << " +/- " << fitter2.getFitParError("xi") << endl;
  cout << "rho = " << fixed << fitter2.getFitParValue("rho") << " +/- " << fitter2.getFitParError("rho") << endl;

  {
  stringstream label2; label2 << "chisq/nDoF=" << fixed << setprecision(0) << fitter2.getChisq() << "/" << fitter2.getNDoF(); 
  label2 << " xi=" << setprecision(4) << fitter2.getFitParValue("xi") << "+/-" << fitter2.getFitParError("xi");
  string plot2 = fitter2.makeFitPlotAxis(0, .5, label2.str());
  //write the plot to a file
  {  
    stringstream ss; ss << "dispnP2P4_fit.ax"; 
    ofstream out; out.open(ss.str().c_str());
    out << plot2;
    out.close();
  }
  }
  //*******************************************************************
  DispnLatBoson* tmp3 = new DispnLatBoson();
  Handle<FitFunction> dispnLatBoson(tmp3);

  dispnLatBoson->setDefaultParValue("m", 0.07);
  dispnLatBoson->setDefaultParError("m", 0.02);

  dispnLatBoson->setDefaultParValue("xi", 3.5);
  dispnLatBoson->setDefaultParError("xi", 0.3);

  MultiEnsemFitter fitter3(data, dispnLatBoson);

  if( ! fitter3.fit() ){
    cout << " FIT FAILED " << endl; exit(1);
  }

  cout << endl << "********************************" << endl;
  cout << "DISPERSION LATTICE BOSON" << endl;

  cout << "chisq/nDoF = " << fixed << setprecision(2) << fitter3.getChisq() << "/" << fitter3.getNDoF() << endl;
  cout << "m = " << fixed << setprecision(5) << fitter3.getFitParValue("m") << " +/- " << fitter3.getFitParError("m") << endl;
  cout << "xi = " << fixed << fitter3.getFitParValue("xi") << " +/- " << fitter3.getFitParError("xi") << endl;
  {
    stringstream label; label << "chisq/nDoF=" << fixed << setprecision(0) << fitter3.getChisq() << "/" << fitter3.getNDoF(); 
    label << " xi=" << setprecision(4) << fitter3.getFitParValue("xi") << "+/-" << fitter3.getFitParError("xi");
    string plot = fitter3.makeFitPlotAxis(0, .5, label.str());
    //write the plot to a file
    {  
      stringstream ss; ss << "dispnLatBoson_fit.ax"; 
      ofstream out; out.open(ss.str().c_str());
      out << plot;
      out.close();
    }
  }

}
