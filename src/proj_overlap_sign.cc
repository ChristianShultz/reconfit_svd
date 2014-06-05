/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : eff_mass_from_proj_ops.cc

 * Purpose :

 * Creation Date : 06-02-2013

 * Last Modified : Thu 05 Jun 2014 11:23:07 AM EDT

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


double SVcutoff = 1e-6; 
double cutoff = 0.2; 
int minTSlices = 6; 


struct ProjFileLine
{
  ProjFileLine(const std::string _proj_op, 
      std::string _op, 
      std::string _weight)
    : proj_op(_proj_op) , op(_op) , weight(_weight) 
  {}

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

  int
main(int argc, char *argv[])
{
  if(argc != 4)
  {
    cerr << "usage: " << argv[0] << ": <xmlinifile> <proj_op_file.list> <comp_op>" << std::endl;
    exit(1);
  }

  string xmlinifile;
  istringstream val(argv[1]);
  val >> xmlinifile;
  std::cout << "Loading xmlinifile: " << xmlinifile << std::endl;

  // Read parameters from xml ini file
  FitIniProps_t inikeys;

  try
  {
    XMLReader xml_in(xmlinifile);
    read(xml_in, "/FitIniParams", inikeys);
  }
  catch(const std::string &e)
  {
    cerr << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }

  //run a sainity check on the ini, this may overwrite some parameters such as the global tmax
  check_ini(inikeys);

  // Load correlation matrix from appropriate source
  SembleRCorrs twoPoints;
  loadCorr(twoPoints, inikeys);

  std::cout << "LOADED " << std::endl;


  std::string proj_op_listf;
  std::istringstream val2(argv[2]); 
  val2 >> proj_op_listf; 
  std::vector<ProjFileLine> proj_ops;

  std::istringstream val3(argv[3]);
  int comp_op_index; 
  val3 >> comp_op_index; 

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
    corr += SEMBLE::toScalar(v[i])*twoPoints.getCij(i,comp_op_index); 

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


  vector<double> tslices;
  for(int t = 0; t < Lt; t++)
    tslices.push_back(t);

  EnsemData corrData(tslices, corr);

  corrData.setSVCutoff(SVcutoff);
  std::cout << "covariance will be inverted cutting off SV below " << corrData.getSVCutoff() << std::endl; 


  Handle<FitComparator> comp1H(new CompareFitsByChisqPerNDoF); 

  //do the fit
  FitCorrelatorExp fitCorr(corrData, comp1H , cutoff, minTSlices); 

  stringstream ss;
  ss << "exp_fit.ax"; 
  ofstream out; 
  out.open(ss.str().c_str());
  out << fitCorr.getFitPlotString();
  out.close();

  //write the jackknife fit files
  EnsemReal outcorr = fitCorr.getMass0();
  stringstream filename; 
  filename << "exp_fit_mass0.jack";  
  write(filename.str(), outcorr );

  outcorr = fitCorr.getAmp0();
  stringstream filename2; 
  filename2 << "exp_fit_amp0.jack"; 
  write(filename2.str(), outcorr );   

  return 0;
}

