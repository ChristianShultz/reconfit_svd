#include "semble_semble.h"
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
  if(argc != 2)
    {
      cerr << "usage: " << argv[0] << ": <xmlinifile> " << endl;
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

  // do fits at all desired t0 values
  SMT0Fit<double> foo(twoPoints, inikeys);

  return 0;

}




