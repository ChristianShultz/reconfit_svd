/*! \file
 *  \brief Factory for reading redstar SU(2)
 */

#include "correlator_reader_factory.h"
#include "correlator_util.h"
#include "redstar_reader.h"
#include "hadron/su2_corr/hadron_npart_npt_corr.h"
#include "hadron/cgc_irrep_mom.h"
#include <adat/handle.h>
#include <xml_array2d.h>
#include <string>
#include <map>

namespace CorrReaderEnv
{
  using namespace Hadron;
  
  namespace RedstarReaderEnv
  { 
    //---------------------------------------------------------------------------
    namespace
    {
      struct InputKeys_t
      {
	int                        source_tslice;
	int                        twoI_z;         /*!< Target isospin component */
	KeyCGCIrrepMom_t           irrep_mom;      /*!< Target irrep row and D-1 momentum of this N-particle op */
      };


      //! Parameter structure
      struct Params
      {
	Params() {}
	Params(XMLReader& xml_in, const std::string& path);

	std::vector<std::string>   dbFnames;
	std::string                opsListFname;
	std::vector<std::string>   opsXMLFiles;
	std::string                rephaseMode;
	InputKeys_t                KeyParams;
      };

      void read(XMLReader &xml, const std::string &path, InputKeys_t& prop)
      {
	XMLReader ptop(xml, path);

	read(ptop, "source_tslice", prop.source_tslice);
	read(ptop, "twoI_z", prop.twoI_z);
	read(ptop, "irmom", prop.irrep_mom);
      }

      //----------------------------------------------------------------------------
      // Param stuff
      Params::Params(XMLReader& xml_in, const std::string& path) 
      {
	XMLReader ptop(xml_in, path);

	// Read program parameters
	read(ptop, "dbFnames", dbFnames);
	read(ptop, "opsListFname", opsListFname);
	read(ptop, "opsXMLFiles", opsXMLFiles);
	read(ptop, "rephaseMode", rephaseMode);
	read(ptop, "KeyParams", KeyParams);
      }


      //----------------------------------------------------------------------------------
      //! Create db keys from operator list
      Array2d<KeyHadronNPartNPtCorr_t> createKeys(const std::vector<KeyHadronNPartIrrepOp_t>& opsxml,
						  int twoI_z, const KeyCGCIrrepMom_t& irmom,
						  int t_source)
      {
	int dim = opsxml.size();

	Array2d<KeyHadronNPartNPtCorr_t> keys(dim, dim);

	for(int j_src =  0; j_src < dim; j_src++)
	{
	  for(int j_snk = 0; j_snk < dim; j_snk++)
	  {
	    KeyHadronNPartNPtCorr_t key;
	    key.npoint.resize(2);
	    
	    // The sink op
	    key.npoint[1].t_slice           = -2;
	    key.npoint[1].irrep.twoI_z      = twoI_z;
	    key.npoint[1].irrep.row         = irmom.row;
	    key.npoint[1].irrep.mom         = irmom.mom;
	    key.npoint[1].irrep.creation_op = false;
	    key.npoint[1].irrep.smearedP    = true;
	    key.npoint[1].irrep.op          = opsxml[j_snk];

	    // The source op
	    key.npoint[2].t_slice           = t_source;
	    key.npoint[2].irrep.twoI_z      = twoI_z;
	    key.npoint[2].irrep.row         = irmom.row;
	    key.npoint[2].irrep.mom         = irmom.mom;
	    key.npoint[2].irrep.creation_op = true;
	    key.npoint[2].irrep.smearedP    = true;
	    key.npoint[2].irrep.op          = opsxml[j_src];

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
	std::map<std::string, KeyHadronNPartIrrepOp_t>    opsMap;
 	std::vector<KeyHadronNPartIrrepOp_t>              opsxml;
	FILEDB::AllConfStoreMultipleDB< ADATIO::SerialDBKey<KeyHadronNPartNPtCorr_t>,  ADATIO::SerialDBData<EnsemScalar<EnsemVectorComplex>::Type_t> > database;
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

	  if (dbtype != "hadronNPartNPtCorr")
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

	  // Read the operator maps
	  opsMap = readOpsMap<KeyHadronNPartIrrepOp_t>(params.opsXMLFiles);

	  // Find the xml for the desired operator list
	  opsxml = findOpsXml<KeyHadronNPartIrrepOp_t>(opsMap, opsList);
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
	  Array2d<KeyHadronNPartNPtCorr_t> keys = createKeys(opsxml, params.KeyParams.twoI_z, params.KeyParams.irrep_mom, params.KeyParams.source_tslice);

	  std::cout << __func__ << ": Will extract nkeys = " << keys.nrows() << std::endl;

	  // Check reading the Db
	  EnsemVectorComplex Test = printKeyValue<KeyHadronNPartNPtCorr_t, EnsemVectorComplex>(keys(0,0), database);
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

      const std::string name = "redstar";
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

  }  // namespace RedstarReaderEnv

} // namespace CorrReaderEnv
