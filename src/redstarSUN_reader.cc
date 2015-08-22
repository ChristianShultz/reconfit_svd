/*! \file
 *  \brief Factory for reading redstarSUN
 */

#include "correlator_reader_factory.h"
#include "correlator_util.h"
#include "redstarSUN_reader.h"
//#include "hadron/su2_corr/hadron_npart_npt_corr.h"
#include "hadron/hadron_sun_npart_npt_corr.h"
#include <adat/handle.h>
#include <xml_array2d.h>
#include <string>
#include <map>

namespace CorrReaderEnv
{
  using namespace Hadron;
  
  namespace RedstarSUNReaderEnv
  { 
    //---------------------------------------------------------------------------
    namespace
    {
      struct InputPropsRedstarKeys_t
      {
	int                        source_tslice;
	KeyCGCSU3_t                flavor;         /*!< Target flavor component */
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
	std::string                badList;
	InputPropsRedstarKeys_t    KeyParams;
      };

      void read(XMLReader &xml, const std::string &path, InputPropsRedstarKeys_t& prop)
      {
	XMLReader ptop(xml, path);

	read(ptop, "flavor", prop.flavor);
	read(ptop, "irmom", prop.irrep_mom);
      }

      void read(XMLReader& xml, const std::string& path, Params& prop)
      {
	XMLReader ptop(xml, path);

	read(ptop, "dbFnames", prop.dbFnames);
	read(ptop, "opsListFname", prop.opsListFname);
	read(ptop, "opsXMLFiles", prop.opsXMLFiles);
	read(ptop, "rephaseMode", prop.rephaseMode);
	read(ptop, "badList", prop.badList);
	read(ptop, "KeyParams", prop.KeyParams);
      }




      //----------------------------------------------------------------------------
      // Param stuff
      Params::Params(XMLReader& xml_in, const std::string& path) 
      {
	try 
	{
	  XMLReader paramtop(xml_in, path);

	  // Read program parameters
	  read(paramtop, "Param", *this);
	}
	catch(const std::string& e) 
	{
	  std::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	  exit(1);
	}
      }


      //----------------------------------------------------------------------------------
      //! Create db keys from operator list
      Array2d<KeyHadronSUNNPartNPtCorr_t> createKeys(const std::vector<KeyHadronSUNNPartIrrep_t>& opsxml,
						     const KeyCGCSU3_t& flavor, const KeyCGCIrrepMom_t& irmom,
						     int t_source)
      {
	int dim = opsxml.size();

	Array2d<KeyHadronSUNNPartNPtCorr_t> keys(dim, dim);

	for(int j_src =  0; j_src < dim; j_src++)
	{
	  for(int j_snk = 0; j_snk < dim; j_snk++)
	  {
	    KeyHadronSUNNPartNPtCorr_t key;
	    key.npoint.resize(2);
	    
	    // The sink op
	    key.npoint[1].t_slice           = -2;
	    key.npoint[1].irrep.flavor      = flavor;
	    key.npoint[1].irrep.irrep_mom   = irmom;
	    key.npoint[1].irrep.creation_op = false;
	    key.npoint[1].irrep.smearedP    = true;
	    key.npoint[1].irrep.op          = opsxml[j_snk].op;

	    // The source op
	    key.npoint[2].t_slice           = t_source;
	    key.npoint[2].irrep.flavor      = flavor;
	    key.npoint[2].irrep.irrep_mom   = irmom;
	    key.npoint[2].irrep.creation_op = true;
	    key.npoint[2].irrep.smearedP    = true;
	    key.npoint[2].irrep.op          = opsxml[j_src].op;

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
	Params                                          params;
	std::vector<std::string>                        opsList;
	std::map<std::string, KeyHadronSUNNPartIrrep_t> opsMap;
 	std::vector<KeyHadronSUNNPartIrrep_t>           opsxml;
	FILEDB::AllConfStoreMultipleDB< ADATIO::SerialDBKey<KeyHadronSUNNPartNPtCorr_t>,  ADATIO::SerialDBData<EnsemScalar<EnsemVectorComplex>::Type_t> > database;
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

	  if (dbtype != "hadronSUNNPartNPtCorr")
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
	  opsMap = readOpsMap<KeyHadronSUNNPartIrrep_t>(params.opsXMLFiles);

	  // Find the xml for the desired operator list
	  opsxml = findOpsXml<KeyHadronSUNNPartIrrep_t>(opsMap, opsList);
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
	  Array2d<KeyHadronSUNNPartNPtCorr_t> keys = createKeys(opsxml, params.KeyParams.flavor, params.KeyParams.irrep_mom, params.KeyParams.source_tslice);

	  std::cout << __func__ << ": Will extract nkeys = " << keys.nrows() << std::endl;

	  // Check reading the Db
	  EnsemVectorComplex Test = printKeyValue<KeyHadronSUNNPartNPtCorr_t, EnsemVectorComplex>(keys(0,0), database);
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

      const std::string name = "REDSTAR_SUN";
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

  }  // namespace RedstarSUNReaderEnv

} // namespace CorrReaderEnv
