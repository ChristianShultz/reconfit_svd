// -*- C++ -*-
/*! \file
 *  \brief Routines to support loading correlators
 */

#ifndef __correlator_util_h__
#define __correlator_util_h__

#include "AllConfStoreDB.h"
#include "io/key_val_db.h"
#include "semble/semble_matrix.h"
#include "semble/semble_algebra.h"
#include "semble_fit_ini_xml.h"

#include <map>
#include <vector>


namespace CorrReaderEnv
{
  //----------------------------------------------------------------------------------
  // filedb (redstar) database specific functions
  //! Get the metadata
  std::string getMetaData(const std::string& dbase);

  //! Read the dbtype from the metadata
  std::string getDBType(const std::string& dbfile);

  //----------------------------------------------------------------------------------
  //! Read ops list
  std::vector<std::string> readOpsList(const std::string& opslistfile);
  

  //----------------------------------------------------------------------------------
  // Rephase complex correlators
  std::vector< SEMBLE::SembleMatrix<double> > rephaseCorrs(const std::vector< SEMBLE::SembleMatrix<std::complex<double> > >& ComplexCorrs,
							   const std::string& rephaseMode, int tmax);


  //----------------------------------------------------------------------------------
  // Various db utils (a lot of these are modified version of dbutils.cc)
  //! Get a key/value
  // New database format
  template<typename K, typename V>
  V printKeyValue(const K& ky,
                  FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<K>,  ADATIO::SerialDBData<typename EnsemScalar<V>::Type_t> >& database)
  {
    typedef typename EnsemScalar<V>::Type_t SV;

    ADATIO::SerialDBKey<K> key;
    key.key() = ky;

    std::vector< ADATIO::SerialDBData<SV> > vals;
    int ret;

    if((ret = database.get(key, vals)) != 0)
    {
      std::cerr << __func__ << ": key not found\n" << ky;
      exit(1);
    }

    V eval;
    eval.resize(vals.size());
    eval.resizeObs(vals[0].data().numElem());

    for(int i = 0; i < vals.size(); ++i)
    {
      SV sval = vals[i].data();
      pokeEnsem(eval, sval, i);
    }

    return eval;
  }


  //----------------------------------------------------------------------------------
  // Read operator operator map from XML
  template<typename T>
  std::map<std::string, T> readOpsMap(const std::vector<std::string>& opsxmlfiles)
  {
    if (opsxmlfiles.empty())
    {
      std::cerr << __func__ << ": ERROR: there must be at least one opsxmlfile" << std::endl;
      exit(1);
    }

    // Ops in a map
    std::map<std::string, T> opsmap;

    for(auto xml = opsxmlfiles.begin(); xml != opsxmlfiles.end(); ++xml)
    {
      // Continually read into the same map
      try
      {
	XMLReader xml_in(*xml);
	read(xml_in, "/OpsList", opsmap);
      }
      catch(const std::string &e)
      {
	std::cerr << __func__ << ": ERROR: can't read opsxmlfile (" << *xml << "): " << e << std::endl;
	exit(1);
      }

      std::cout << __func__ << ": currently at " << opsmap.size() << " elems in after reading " << *xml << std::endl;
    }

    return opsmap;
  }


  //----------------------------------------------------------------------------------
  // Find operator xml-s
  template<typename T>
  std::vector<T> findOpsXml(const std::map<std::string, T> opsmap, const std::vector<std::string>& opsList)
  {
    // Find the ops in the list
    std::vector<T> keys;
	
    for(auto op = opsList.begin(); op != opsList.end(); ++op)
    {
      auto ptr = opsmap.find(*op);
      if (ptr == opsmap.end())
      {
	std::cerr << __func__ << ": ERROR: can't find op = " << *op << "  in opsxmlfiles\n";
	exit(1);
      }

      keys.push_back(ptr->second);
    }

    return keys;
  }


  //----------------------------------------------------------------------------------
  // Get correlators from database, averaging over lorentz/spin if required
  template<typename T>
  std::vector< SEMBLE::SembleMatrix<std::complex<double> > > loadCorrs(FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<T>, ADATIO::SerialDBData<EnsemScalar<EnsemVectorComplex>::Type_t> >& database,
								       const Array2d<T>& keys)
  {
    if ((keys.nrows() != keys.ncols()) || keys.nrows() == 0)
    {
      std::cerr << __func__ << ": matrix of keys not square\n";
      exit(1);
    }
    
    int dim = keys.nrows();
    EnsemVectorComplex Test = printKeyValue<T, EnsemVectorComplex>(keys(0,0), database);
    int Lt = Test.numElem();
    int nbins = peekObs(Test, 0).size();

    SEMBLE::SembleMatrix< std::complex<double> > dum(nbins, dim, dim);
    std::vector< SEMBLE::SembleMatrix<std::complex<double> > > ComplexCorrs;
    ComplexCorrs.resize(Lt, dum);

    for(int j_src = 0; j_src < dim; j_src++)
    {
      for(int j_snk = 0; j_snk < dim; j_snk++)
      {
	EnsemVectorComplex TempCorr = printKeyValue<T, EnsemVectorComplex>(keys(j_src,j_snk), database);

	for(int t = 0; t < Lt; t++)
	{
	  ComplexCorrs[t].loadEnsemElement(j_src, j_snk, peekObs(TempCorr, t));
	}   // loop over t
      }   // loop over j_snk
    }   // loop over j_src

    return ComplexCorrs;
  }


}

#endif
