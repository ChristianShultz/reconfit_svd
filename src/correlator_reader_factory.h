// -*- C++ -*-
/*! \file
 *  \brief Factories of correlator readers
 */

#ifndef __correlator_reader_factory_h__
#define __correlator_reader_factory_h__

#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "correlator_reader.h"
#include <string>

namespace CorrReaderEnv
{
  using namespace Util;

  //---------------------------------------------------------------------------
  //! Hadron flavor factory (foundry)
  typedef SingletonHolder< 
    ObjectFactory<RealCorrReader, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  RealCorrReader* (*)(XMLReader&, const std::string&),
		  StringFactoryError> >
  TheCorrReaderFactory;

  bool registerAll();

}

#endif
