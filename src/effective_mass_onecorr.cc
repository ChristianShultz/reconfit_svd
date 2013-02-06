// effective_mass_onecorr.cc -
//
// Thursday, June 28 2012
//

#include "semble/semble_semble.h"
#include "semble_load_correlators.h"
#include "semble_fit_ini_xml.h"
#include "semble_multi_t0fit.h"
#include <itpp/itbase.h>
#include <vector>
#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;
using namespace ENSEM;
using namespace SEMBLE;

int main(int argc, char *argv[])
{
  if(argc != 5)
    {
      cerr << "usage: " << argv[0] << ": <xmlinifile>  <row> <col> <delta_t>" << endl;
      exit(1);
    }


  string xmlinifile;
  { istringstream val(argv[1]); val >> xmlinifile; }
  cout << "Loading xmlinifile: " << xmlinifile << endl;

  int row; {istringstream val(argv[2]); val >> row;}
  int col; {istringstream val(argv[3]); val >> col;}
  int deltat; {istringstream val(argv[4]); val >> deltat;}

  // Read parameters from xml ini file
  FitIniProps_t inikeys;

  try
    {
      XMLReader xml_in(xmlinifile);
      read(xml_in, "/FitIniParams", inikeys);
    }
  catch(const std::string &e)
    {
      cerr << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << endl;
      exit(1);
    }

  //run a sainity check on the ini, this may overwrite some parameters such as the global tmax
  check_ini(inikeys);

  // Load correlation matrix from appropriate source
  SembleRCorrs twoPoints;
  loadCorr(twoPoints, inikeys);

  cout << "LOADED " << endl;

  EnsemVectorReal oneCorr = twoPoints.getCij(row,col);
  const int loopt = twoPoints.getLt() - deltat;

  ENSEM::write(std::string("oneCorr.jack"),oneCorr); 

  EnsemVectorReal effMass;
  effMass.resize(oneCorr.size());
  effMass.resizeObs(loopt);
  
  for(int i = 0; i < loopt; i++)
    pokeObs(effMass, log(peekObs(oneCorr,i)/peekObs(oneCorr,i+deltat))/Real(double(deltat)), i);
      
  AxisPlot effMassPlot;
  effMassPlot.addEnsemData(effMass,"\\sq",1);

  effMassPlot.sendToFile(std::string("effMass.ax"));

  ENSEM::write(std::string("effMass.jack"),effMass); 
  
  EnsemVectorReal dt_effMass;
  dt_effMass.resize(effMass.size());
  dt_effMass.resizeObs(loopt-1);

  // forward derivative of effective mass -- should be zero
  for(int i = 0; i < loopt-1; i++)
    pokeObs(dt_effMass,peekObs(effMass,i+1) - peekObs(effMass,i),i);


  AxisPlot dt_effMassPlot;
  dt_effMassPlot.addEnsemData(dt_effMass,"\\sq",1);
  
  dt_effMassPlot.sendToFile(std::string("dt_effMass.ax"));

  return 0;
}
