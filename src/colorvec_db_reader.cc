/*! \file
 *  \brief ColorVec Hadron2PtCorr reader
 */

#include "correlator_reader_factory.h"
#include "correlator_util.h"
#include "colorvec_db_reader.h"
#include <formfac/hadron_2pt_corr.h>
#include <adat/handle.h>
#include <xml_array2d.h>
#include <string>
#include <map>

//#error "NOT SUPPORTING THIS FORMAT ANYMORE..."

namespace CorrReaderEnv
{
  namespace ColorVecDBReaderEnv
  { 
    //---------------------------------------------------------------------------
    namespace
    {
      //! Parameter structure
      struct Params
      {
	Params() {}
	Params(XMLReader& xml_in, const std::string& path);

	std::vector<std::string>   dbFnames;
	std::string                opsListFname;
	std::string                rephaseMode;
	FF::KeyHadron2PtCorr_t     KeyParams;
      };

      //----------------------------------------------------------------------------
      // Param stuff
      Params::Params(XMLReader& xml_in, const std::string& path) 
      {
	XMLReader ptop(xml_in, path);

	// Read program parameters
	read(ptop, "dbFnames", dbFnames);
	read(ptop, "opsListFname", opsListFname);
	read(ptop, "rephaseMode", rephaseMode);
	read(ptop, "KeyParams", KeyParams);
      }


      //----------------------------------------------------------------------------------
      //! Create db keys from operator list
      Array2d<FF::KeyHadron2PtCorr_t> createKeys(const std::vector<std::string>& opsList,
						 const FF::KeyHadron2PtCorr_t& DefaultKey)
      {
	int dim = opsList.size();
	
	Array2d<FF::KeyHadron2PtCorr_t> keys(dim, dim);

	for(int j_src =  0; j_src < dim; j_src++)
	{
	  for(int j_snk = 0; j_snk < dim; j_snk++)
	  {
   	    FF::KeyHadron2PtCorr_t key = DefaultKey;
	    
	    // This is "spin" in the old code
	    int i = 1;
	    
	    key.src_name = opsList[j_src];
	    key.snk_name = opsList[j_snk];
	    key.src_spin = i;
	    key.snk_spin = i;

	    keys(j_src,j_snk) = key;
	    
	  } // loop over j_snk
	} // loop over j_src
	
	return keys;
      }

 

      //----------------------------------------------------------------------------------
      //----------------------------------------------------------------------------------
      //----------------------------------------------------------------------------------
      //! Reader
      class CorrReaderImpl : public ComplexCorrReader
      {
      public:
	//! Constructor
	CorrReaderImpl(const Params& p_);

	//! Virtual destructor
	virtual ~CorrReaderImpl() {}

      protected:
	//! Read the complex version of the correlators
	virtual std::vector< SEMBLE::SembleMatrix<std::complex<double> > > getComplexCorrs();

	//! Rephasing mode
	virtual std::string getRephaseMode() const {return params.rephaseMode;}
  
	//! This is "required" as a convenience - it is nice to see bad operators, whether or not they rephase
	virtual std::vector<std::string> getOpsList() const {return opsList;}

      private:
	Params                                            params;
	std::vector<std::string>                          opsList;
	FILEDB::AllConfStoreMultipleDB< ADATIO::SerialDBKey<FF::KeyHadron2PtCorr_t>,  ADATIO::SerialDBData<EnsemScalar<EnsemVectorComplex>::Type_t> > database;
      };


      //----------------------------------------------------------------------------------
      //! Constructor
      CorrReaderImpl::CorrReaderImpl(const Params& p_) : params(p_)
      {
	if (params.dbFnames.empty())
	{
	  std::cerr << __func__ << ": error - expect at least one ensemble file\n";
	  exit(1);
	}

	try
	{
	  // Read the dbtype from the metadata of the first file
	  std::string dbtype = getDBType(params.dbFnames[0]);

	  if (dbtype != "hadron2Pt")
	  {
	    std::cerr << "Error - corr edb dbFname = " << params.dbFnames[0] << "  not appropriate type, found type = " << dbtype << std::endl;
	    exit(1);
	  }

	  // Open DB
	  if (database.open(params.dbFnames) != 0)
	  {
	    std::cerr << __func__ << ": error opening databases " << params.dbFnames[0] << std::endl;
	    exit(1);
	  }
	  
	  // Read the operator list file
	  opsList = readOpsList(params.opsListFname);
	}
	catch(const std::string &e)
	{
	  std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
	  exit(1);
	}
	catch(std::exception &e)
	{
	  std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
	  exit(1);
	}

	std::cout << __func__ << ": finished loading the complex correlators" << std::endl;
      }


      //----------------------------------------------------------------------------------
      //! Return the complex correlators
      std::vector< SEMBLE::SembleMatrix<std::complex<double> > > CorrReaderImpl::getComplexCorrs()
      {
	try
	{
	  // And construct the keys
	  Array2d<FF::KeyHadron2PtCorr_t> keys = createKeys(opsList, params.KeyParams);
	  std::cout << __func__ << ": Will extract nkeys = " << keys.nrows() << std::endl;

	  // Check reading the Db
	  EnsemVectorComplex Test = printKeyValue<FF::KeyHadron2PtCorr_t, EnsemVectorComplex>(keys(0,0), database);
	  int Lt = Test.numElem();
	  int nbins = peekObs(Test, 0).size();
	  std::cout << __func__ << ": filedb database (" << params.dbFnames[0] << ") has Lt = " << Lt << ", nbins = " << nbins << std::endl;
      
	  // Read the correlators
	  return loadCorrs(database, keys);
	}
	catch(const std::string &e)
	{
	  std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
	  exit(1);
	}
	catch(std::exception &e)
	{
	  std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
	  exit(1);
	}

	std::cout << __func__ << ": finished loading the complex correlators" << std::endl;
      }


    }


    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    namespace
    {
      RealCorrReader* createReader(XMLReader& xml_in, const std::string& path) 
      {
	return new CorrReaderImpl(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "dbnew";
    }

    //---------------------------------------------------------------------------
    //! Register all the factories
    bool registerAll() 
    {
      if (registered) {return true;}
      bool success = true; 

      success &= TheCorrReaderFactory::Instance().registerObject(name, createReader);

      return success &= (registered = true);
    }

  }  // namespace ColorVecDBReaderEnv

} // namespace CorrReaderEnv
