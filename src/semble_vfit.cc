#include "semble/semble_semble.h"
#include "correlator_reader.h"
#include "correlator_reader_factory.h"
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
  if(argc != 2)
  {
    std::cerr << "usage: " << argv[0] << ": <xmlinifile> " << std::endl;
    exit(1);
  }

  bool doReg = CorrReaderEnv::registerAll();

  std::string xmlinifile;
  std::istringstream val(argv[1]);
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
    std::cerr << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }

  //run a sainity check on the ini, this may overwrite some parameters such as the global tmax
  check_ini(inikeys);

  // Load correlation matrix from appropriate source
  std::cout << "Correlator reader successfully initialized" << std::endl;
  std::vector< SEMBLE::SembleMatrix<double> > twoPoints = CorrReaderEnv::getCorrs(inikeys.inputProps);
  std::cout << "LOADED " << std::endl;

  // do fits at all desired t0 values
  SMT0Fit<double> foo(twoPoints, inikeys);

  return 0;

}




