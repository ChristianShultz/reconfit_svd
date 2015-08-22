// -*- C++ -*-
/*! \file
 * \brief Base class for reading correlators from disk
 */

#ifndef __correlator_reader_h__
#define __correlator_reader_h__

#include "io/adat_io.h"
#include "io/adat_xmlio.h"

#include "ensem/ensem.h"

#include "semble/semble_matrix.h"
//#include "semble/semble_algebra.h"
//#include "semble_fit_ini_xml.h"
#include <vector>

namespace CorrReaderEnv
{
  using namespace ENSEM;
  using namespace ADATIO;
  using namespace ADATXML;
  using namespace ADAT;

  //----------------------------------------------------------------------------------
  //! Read correlators matrix
  class RealCorrReader
  {
  public:
    //! Virtual destructor
    virtual ~RealCorrReader() {}

    //! Only one thing - read the correlators
    virtual std::vector< SEMBLE::SembleMatrix<double> > getCt() = 0;
  };




  //----------------------------------------------------------------------------------
  //! Read complex correlators matrix
  class ComplexCorrReader : public RealCorrReader
  {
  public:
    //! Virtual destructor
    virtual ~ComplexCorrReader() {}

    //! Can provide a default version
    virtual std::vector< SEMBLE::SembleMatrix<double> > getCt();

  protected:
    //! Read the complex corrs
    virtual std::vector< SEMBLE::SembleMatrix< std::complex<double> > > getComplexCorrs() = 0;

    //! Rephasing mode
    virtual std::string getRephaseMode() const = 0;
  
    //! This is "required" as a convenience - it is nice to see bad operators, whether or not they rephase
    virtual std::vector<std::string> getOpsList() const = 0;
};

} // namespace SEMBLE

#endif
