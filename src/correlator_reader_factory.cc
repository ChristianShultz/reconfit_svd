/*! \file
 *  \brief Factories of correlator readers
 */

#include <adat/handle.h>
#include "correlator_reader.h"
#include "correlator_reader_factory.h"
#include "correlator_util.h"
//#include "redstar_reader.h"
#include "redstarSUN_reader.h"

#include <string>

namespace CorrReaderEnv
{
  //----------------------------------------------------------------------------------
  //! Return the complex correlators
  std::vector< SEMBLE::SembleMatrix<double> > ComplexCorrReader::getCt()
  {
    // Actually do the rephasing
    return rephaseCorrs(this->getComplexCorrs(), getRephaseMode(), getOpsList(), 20); //stupid hardwire of tmax
  }


  //----------------------------------------------------------------------------------
  //! For the lazy, this is a quick way to return correlators
  std::vector< SEMBLE::SembleMatrix<double> > getCorrs(const GroupXML_t& inputProps)
  {
    try
    {
      std::cout << "Initalize correlator reader" << std::endl;
  
      std::istringstream  xml_l(inputProps.xml);
      XMLReader  linktop(xml_l);
      std::cout << "Correlator type = " << inputProps.id << std::endl;
	
      Handle<RealCorrReader> corrs(TheCorrReaderFactory::Instance().createObject(inputProps.id,
										 linktop, 
										 inputProps.path));
    
      // Read correlators
      return corrs->getCt();
    }
    catch(const std::string& e) 
    {
      std::cerr << __func__ << ": Caught Exception correlator reader: " << e << std::endl;
      exit(1);
    }
  }


  //---------------------------------------------------------------------------
  namespace
  {
    //! Local registration flag
    bool registered = false;
  }

  //---------------------------------------------------------------------------
  // Register all the possible readers
  bool registerAll() 
  {
    if (registered) {return true;}

    bool success = true;
    success &= RedstarSUNReaderEnv::registerAll();

    return success &= (registered = true);
  }

}  // end namespace CorrReaderEnv
