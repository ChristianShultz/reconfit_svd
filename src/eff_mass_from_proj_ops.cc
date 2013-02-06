/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : eff_mass_from_proj_ops.cc

* Purpose :

* Creation Date : 06-02-2013

* Last Modified : Wed Feb  6 10:34:38 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "semble/semble_semble.h"
#include "semble_load_correlators.h"
#include "semble_fit_ini_xml.h"
#include "semble_multi_t0fit.h"
#include <itpp/itbase.h>
#include <vector>
#include <ostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>


struct ProjFileLine
{
  ProjFileLine(const std::string _proj_op, std::string _op, std::string _weight)
  : proj_op(_proj_op) , op(_op) , weight(_weight) {}

  double get_weight(void) const
  {
    return atof(weight.c_str()); 
  }

  std::string proj_op;
  std::string op;
  std::string weight; 
};

using namespace std;
using namespace ENSEM;
using namespace SEMBLE;

int main(int argc, char *argv[])
{
  if(argc != 3)
    {
      cerr << "usage: " << argv[0] << ": <xmlinifile> <proj_op_file.list> " << endl;
      exit(1);
    }

  string xmlinifile;
  istringstream val(argv[1]);
  val >> xmlinifile;
  cout << "Loading xmlinifile: " << xmlinifile << endl;

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


  std::string proj_op_listf;
  std::istringstream val2(argv[2]); 
  val2 >> proj_op_listf; 
  std::vector<ProjFileLine> proj_ops;
   
  std::ifstream file;
  file.open(proj_op_listf.c_str()); 
  std::string aline;
  while(std::getline(file,aline))
  {
    std::istringstream foo(aline); 
    std::string proj, op, weight; 

    if(!!!(foo >> proj >> op >> weight))
    {
      std::cerr << "Error reading " << proj_op_listf << std::endl;
      exit(1);
    }

    proj_ops.push_back(ProjFileLine(proj,op,weight)); 
  }
  
  std::cout << "read the following..." << std::endl;

  std::vector<double> v(proj_ops.size()); 
  for(int i =0; i < proj_ops.size(); ++i)
  {
    std::cout << proj_ops[i].proj_op << " " << proj_ops[i].op << " " << proj_ops[i].get_weight() << std::endl;
    v[i] = proj_ops[i].get_weight();
  }


  ENSEM::EnsemVectorReal corr = twoPoints.getCij(0,0);
  corr = SEMBLE::toScalar(double(0.)); 

  for(int i = 0; i < v.size(); ++i)
    for(int j = 0; j < v.size(); ++j)
      corr += SEMBLE::toScalar(v[i]*v[j])*twoPoints.getCij(i,j); 

  ENSEM::write(proj_ops[0].proj_op + std::string("_corr.jack"), corr); 
  
  int Lt = corr.numElem(); 
  int deltat = 1; 
  EnsemVectorReal effMass;
  effMass.resize(corr.size()); 
  effMass.resizeObs(Lt - deltat);

  for(int t = 0; t < Lt - deltat; ++t)
    ENSEM::pokeObs(effMass,ENSEM::log(ENSEM::peekObs(corr,t)/ENSEM::peekObs(corr,t+deltat))/SEMBLE::toScalar(double(deltat)),t);


  AxisPlot effMassPlot;
  effMassPlot.addEnsemData(effMass,"\\sq",1);

  effMassPlot.sendToFile(proj_ops[0].proj_op + std::string("__effMass.ax"));

  ENSEM::write(proj_ops[0].proj_op + std::string("_effMass.jack"),effMass); 

  
  EnsemVectorReal dt_effMass;
  dt_effMass.resize(effMass.size());
  dt_effMass.resizeObs(Lt-deltat-1);

  // forward derivative of effective mass -- should be zero
  for(int i = 0; i <Lt-deltat-1; i++)
    pokeObs(dt_effMass,peekObs(effMass,i+1) - peekObs(effMass,i),i);


  AxisPlot dt_effMassPlot;
  dt_effMassPlot.addEnsemData(dt_effMass,"\\sq",1);
  
  dt_effMassPlot.sendToFile(proj_ops[0].proj_op + std::string("dt_effMass.ax"));

  return 0;
  

}

