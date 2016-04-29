/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : proj_op_util.cc

 * Purpose :

 * Creation Date : 18-06-2014

 * Last Modified : Thu 19 Jun 2014 01:29:11 PM EDT

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
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace ENSEM;
using namespace SEMBLE;



typedef std::pair<SembleRCorrs,std::vector<std::string> > usage_t;


SembleRCorrs grab_correlation_matrix(const std::string &ini)
{
  // Read parameters from xml ini file
  FitIniProps_t inikeys;

  try
  {
    XMLReader xml_in(ini);
    read(xml_in, "/FitIniParams", inikeys);
  }
  catch(const std::string &e)
  {
    cerr << __func__ << ": ERROR: can't read xmlinifile (" << ini << "): " << e << endl;
    exit(1);
  }

  //run a sainity check on the ini, this may overwrite some parameters such as the global tmax
  check_ini(inikeys);

  // Load correlation matrix from appropriate source
  SembleRCorrs twoPoints;
  loadCorr(twoPoints, inikeys);

  return twoPoints;
}

usage_t usage(int argc, char *argv[])
{
  if(argc < 2)
  {
    cerr << "usage: " << argv[0] << ": <xmlinifile> <operation> <args>" << endl;
    exit(1);
  }

  string xmlinifile;
  std::vector<std::string> commands;

  {
    istringstream val(argv[1]);
    val >> xmlinifile;
  }


  for(int i = 2; i < argc; ++i)
  {

    istringstream val(argv[i]);
    std::string tmp; 
    val >> tmp;
    commands.push_back(tmp);
  }

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

  return usage_t(twoPoints,commands);
}

//**********************************************
//**********************************************
//**********************************************
//**********************************************


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

void check_exit_pars( const int np, const std::vector<std::string> &p, const std::string &msg, const char *f)
{
  if( np != p.size()) 
  {
    std::cout << "error: usage " << f << " <sfit.ini.xml> <operation> " << msg << std::endl;  
    exit(1); 
  }
}

using namespace std;
using namespace ENSEM;
using namespace SEMBLE;

std::vector<ProjFileLine> read_proj_ops(const std::string &proj_op_listf)
{
  std::vector<ProjFileLine> proj_ops;

  std::ifstream file;
  file.open(proj_op_listf.c_str()); 
  std::string aline;
  while(std::getline(file,aline))
  {
    std::istringstream foo(aline); 
    std::string proj, op, weight; 

    if (aline.find("version") != std::string::npos)
      continue;

    if(!!!(foo >> proj >> op >> weight))
    {
      std::cerr << "Error reading " << proj_op_listf << std::endl;
      exit(1);
    }

    proj_ops.push_back(ProjFileLine(proj,op,weight)); 
  }

  return proj_ops; 
}

//**********************************************

void  effective_mass(const usage_t &u)
{
  std::vector<std::string> options = u.second; 

  check_exit_pars(2,options," <proj.list> ",__func__); 

  SembleRCorrs twoPoints = u.first;

  std::vector<ProjFileLine> proj_ops = read_proj_ops(options[1]); 

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

  ENSEM::write(proj_ops[0].proj_op + std::string("_corr.dat"), corr); 

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

  ENSEM::write(proj_ops[0].proj_op + std::string("_effMass.dat"),effMass); 

}

//**********************************************

void  check_ortho(const usage_t &u)
{
  std::vector<std::string> options = u.second; 

  check_exit_pars(3,options," <proj_left.list> <proj_right.list> ",__func__); 

  SembleRCorrs twoPoints = u.first;

  std::vector<ProjFileLine> proj_ops_l = read_proj_ops(options[1]); 
  std::vector<ProjFileLine> proj_ops_r = read_proj_ops(options[2]); 


  std::vector<double> vl(proj_ops_l.size()); 
  for(int i =0; i < proj_ops_l.size(); ++i)
  {
    vl[i] = proj_ops_l[i].get_weight();
    std::cout << proj_ops_l[i].proj_op << " " << proj_ops_l[i].op << " " << proj_ops_l[i].get_weight() << std::endl;
  }

  std::vector<double> vr(proj_ops_r.size()); 
  for(int i =0; i < proj_ops_r.size(); ++i)
  {
    vr[i] = proj_ops_r[i].get_weight();
    std::cout << proj_ops_r[i].proj_op << " " << proj_ops_r[i].op << " " << proj_ops_r[i].get_weight() << std::endl;
  }


  std::cout << "vl.size() " << vl.size() << " vr.size() " << vr.size() << std::endl;

  ENSEM::EnsemVectorReal corr = twoPoints.getCij(0,0);
  corr = SEMBLE::toScalar(double(0.)); 

  for(int i = 0; i < vl.size(); ++i)
    for(int j = 0; j < vr.size(); ++j)
      corr += SEMBLE::toScalar(vl[i]*vr[j])*twoPoints.getCij(i,j); 

  ENSEM::write(std::string("check_ortho_") + proj_ops_l[0].proj_op + std::string("_") + proj_ops_r[0].proj_op + std::string(".dat"), corr); 

  int Lt = corr.numElem(); 
  int deltat = 1; 

}

//**********************************************
void project_against(const usage_t &u)
{ 
  std::vector<std::string> options = u.second; 
  check_exit_pars(3,options," <proj.list> <comp_op> ",__func__); 

  // Load correlation matrix from appropriate source
  SembleRCorrs twoPoints = u.first; 

  std::vector<ProjFileLine> proj_ops = read_proj_ops(options[1]); 

  std::istringstream val3(options[2]);
  int comp_op_index; 
  val3 >> comp_op_index; 

  std::vector<double> v(proj_ops.size()); 
  for(int i =0; i < proj_ops.size(); ++i)
    v[i] = proj_ops[i].get_weight();


  ENSEM::EnsemVectorReal corr = twoPoints.getCij(0,0);
  corr = SEMBLE::toScalar(double(0.)); 

  // make the correlator 
  for(int i = 0; i < v.size(); ++i)
    corr += SEMBLE::toScalar(v[i])*twoPoints.getCij(i,comp_op_index); 

  ENSEM::write("project_against_" + proj_ops[0].proj_op + std::string("_op") + options[2] + std::string(".dat"), corr); 

  double SVcutoff = 1e-6; 
  double cutoff = 0.2; 
  int minTSlices = 6; 
  int maxT = 25; 

  vector<double> tslices;
  ENSEM::EnsemVectorReal corr_fit;
  corr_fit.resize( corr.size() ); 
  corr_fit.resizeObs( maxT ); 
  for(int t = 0; t < maxT; t++)
  {
    tslices.push_back(t);
    ENSEM::pokeObs(corr_fit,ENSEM::peekObs(corr,t),t); 
  }

  EnsemData corrData(tslices, corr_fit);

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
  filename << "exp_fit_mass0.dat";  
  write(filename.str(), outcorr );

  outcorr = fitCorr.getAmp0();
  stringstream filename2; 
  filename2 << "exp_fit_amp0.dat"; 
  write(filename2.str(), outcorr );   
}


//**********************************************
//**********************************************

typedef void (*fptr)(const usage_t &);  
std::map<std::string,fptr> ptr_map; 

void init_map(void)
{
  ptr_map.insert(std::pair<std::string,fptr>("effective_mass",
        &effective_mass)); 
  ptr_map.insert(std::pair<std::string,fptr>("check_ortho",
        &check_ortho)); 
  ptr_map.insert(std::pair<std::string,fptr>("project_against",
        &project_against)); 
}


// dispatch 
void work_handler(const usage_t &u)
{
  init_map(); 

  if(ptr_map.find(u.second.front()) == ptr_map.end())
  {
    std::cerr << "unrecognized command" << u.second.front() << std::endl; 
    std::map<std::string,fptr>::const_iterator it;
    for(it = ptr_map.begin(); it != ptr_map.end(); ++it)
      std::cerr << it->first << std::endl;  
    exit (1); 
  }

  fptr foo = ptr_map[u.second.front()]; 

  foo(u); // yup, I just made a foo of u 
}



  int
main(int argc, char *argv[])
{
  usage_t fred = usage(argc,argv); 

  work_handler(fred); 

  return 0;
}





