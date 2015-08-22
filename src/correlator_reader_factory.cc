/*! \file
 *  \brief Factories of correlator readers
 */

#include "correlator_reader.h"
#include "correlator_util.h"
//#include "redstar_reader.h"
#include "redstarSUN_reader.h"

#include <string>

namespace CorrReaderEnv
{
  //----------------------------------------------------------------------------------
  //! Return the complex correlators
  std::vector< SEMBLE::SembleMatrix<double> > ComplexCorrReader::getCt() const
  {
    // Actually do the rephasing
    return rephaseCorrs(this->getComplexCorrs(), getRephaseMode(), 20); //stupid hardwire of tmax
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
