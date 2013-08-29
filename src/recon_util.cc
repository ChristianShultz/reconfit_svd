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

void matrix_plot(const usage_t &u)
{
  SembleRCorrs tp = u.first; 
  std::vector<std::string> opts = u.second; 

  if ( opts.size() != 2 ) 
  {
    std::cerr << __func__ << ": expected <option=matrix_plot> <t0_extract> " << std::endl; 
    exit (1);
  } 

  int t = atoi( opts[1].c_str() ); 


  SEMBLE::SembleMatrix<double> s = tp.getCt(t); 

  itpp::Mat<double> m = s.mean();

  std::stringstream f; 
  f << "matrix_plot_" << opts[1] << ".dat"; 
  std::ofstream out; 
  out.open(f.str().c_str()); 

  for (int i = 0; i < m.rows(); ++i)
  {
    for(int j = 0; j < m.cols(); ++j) 
      out << m(i,j) << ", "; 

    out << "\n";
  }

  out.close(); 
}


void norm_matrix_plot(const usage_t &u)
{
  SembleRCorrs tp = u.first; 
  std::vector<std::string> opts = u.second; 

  if ( opts.size() != 2 ) 
  {
    std::cerr << __func__ << ": expected <option=matrix_plot> <t0_extract> " << std::endl; 
    exit (1);
  } 

  int t = atoi( opts[1].c_str() ); 


  SEMBLE::SembleMatrix<double> s = tp.getCt(t); 

  itpp::Mat<double> m = s.mean();

  std::stringstream f; 
  f << "matrix_plot_" << opts[1] << ".dat"; 
  std::ofstream out; 
  out.open(f.str().c_str()); 

  for (int i = 0; i < m.rows(); ++i)
  {
    for(int j = 0; j < m.cols(); ++j) 
      out << m(i,j)/sqrt(fabs(m(i,i)*m(j,j))) << ", "; 

    out << "\n";
  }

  out.close(); 
}

namespace 
{
  int do_round(const double &d) {return int(d+0.5);}
}

void sign_matrix_plot(const usage_t &u)
{
  SembleRCorrs tp = u.first; 
  std::vector<std::string> opts = u.second; 

  if ( opts.size() != 2 ) 
  {
    std::cerr << __func__ << ": expected <option=matrix_plot> <t0_extract> " << std::endl; 
    exit (1);
  } 

  int t = atoi( opts[1].c_str() ); 


  SEMBLE::SembleMatrix<double> s = tp.getCt(t); 

  itpp::Mat<double> m = s.mean();

  std::stringstream f; 
  f << "matrix_plot_" << opts[1] << ".dat"; 
  std::ofstream out; 
  out.open(f.str().c_str()); 

  for (int i = 0; i < m.rows(); ++i)
  {
    for(int j = 0; j < m.cols(); ++j) 
      out << do_round(m(i,j)/sqrt(fabs(m(i,i)*m(j,j)))) << ", "; 

    out << "\n";
  }

  out.close(); 
}

typedef void (*fptr)(const usage_t &);  
std::map<std::string,fptr> ptr_map; 

void init_map(void)
{
  ptr_map.insert(std::pair<std::string,fptr>("matrix_plot",&matrix_plot)); 
  ptr_map.insert(std::pair<std::string,fptr>("norm_matrix_plot",&norm_matrix_plot)); 
  ptr_map.insert(std::pair<std::string,fptr>("sign_matrix_plot",&sign_matrix_plot)); 
}

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





