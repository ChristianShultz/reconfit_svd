/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : recon_dbutil.cc

 * Purpose :

 * Creation Date : 05-10-2012

 * Last Modified : Fri Oct  5 17:12:21 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>
#include <vector>
#include <utility>
#include <map>
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"
#include "semble_fit_ini_xml.h"
#include "hadron/hadron_npart_irrep.h"
#include "ensem/ensem.h"
#include "semble/semble_key_val_db.h"
#include "semble/semble_meta.h"
#include "AllConfStoreDB.h"


// useful object
struct NameListXML
{
  ADATXML::Array<std::string> stateNameList;
  ADATXML::Array<std::string> operatorNameList;

  void read(ADATXML::XMLReader &xml, const std::string &path)
  {
    ADATXML::XMLReader ptop(xml,path);
    if(ptop.count("stateNameList") > 0)
      ADATXML::read(ptop,"stateNameList",stateNameList);
    else
    {
      std::cerr << __func__ << ": need to provide a valid name list file" << std::endl;
      exit(1);
    }

    if(ptop.count("operatorNameList") > 0)
      ADATXML::read(ptop,"operatorNameList",operatorNameList);
    else
    {
      std::cerr << __func__ << ": need to provide a valid name list file" << std::endl;
      exit(1);
    }
  }
};

// read the ops phases to get the number associated with each operator
std::map<std::string, int> readOpsPhases(void)
{
  std::map<std::string, int> ret;
  std::ifstream infile;
  infile.open("ops_phases");

  while(!infile.eof())
  {
    std::string line;
    getline(infile,line);
    std::stringstream ss(line);
    std::string buf;
    std::vector<std::string> tokens;

    while(ss >> buf)
      tokens.push_back(buf);

    if(!tokens.empty()) // guard against empty lines
      ret[tokens[1]] = ::atoi(tokens[0].c_str());
  }

  infile.close();

  return ret;
}


// opsxml input structure
struct opsxml_t
{
  std::string opnamekey;                     // operator name
  Hadron::KeyHadronNPartIrrep_t irrep;  // operator xml structure
};

void read(XMLReader &xml, const std::string &path, opsxml_t &param)
{
  XMLReader paramtop(xml, path);

  if(paramtop.count("OpNameKey") > 0)
    read(paramtop, "OpNameKey", param.opnamekey);
  else
    param.opnamekey = "";

  if(paramtop.count("Irrep") > 0)
    read(paramtop, "Irrep", param.irrep);
}

// get the keys we care about
std::map<std::string, Hadron::KeyHadronNPartIrrep_t> readOpsXML(const ADATXML::Array<std::string> &opsxmlfiles, const NameListXML &list)
{
  if(opsxmlfiles.size() < 1)
  {
    std::cerr << __func__ << ": ERROR: there must be at least one opsxmlfile" << std::endl;
    exit(1);
  }


  ADATXML::Array<opsxml_t> opsxml;

  for(int i = 0; i < opsxmlfiles.size(); i++)
  {
    ADATXML::Array<opsxml_t> opsxml_temp;

    try
    {
      ADATXML::XMLReader xml_in(opsxmlfiles[i]);
      read(xml_in, "/OpsList", opsxml_temp);
    }
    catch(const std::string &e)
    {
      std::cerr << __func__ << ": ERROR: can't read opsxmlfile (" << opsxmlfiles[i] << "): " << e << std::endl;
      exit(1);
    }

    std::cout << __func__ << ": read " << opsxml_temp.size() << " elems from " << opsxmlfiles[i] << std::endl;

    if(i == 0)
    {
      opsxml = opsxml_temp;
    }
    else
    {
      ADATXML::Array<opsxml_t> opsxml_orig = opsxml;
      opsxml = concat(opsxml_orig, opsxml_temp);
    }
  }

  std::map<std::string,Hadron::KeyHadronNPartIrrep_t> keymap;

  for(int i = 0; i < list.operatorNameList.size(); ++i)
  {
    std::string op = list.operatorNameList[i];

    bool found = false;
    for(int j = 0; j < opsxml.size(); ++j)
    {
      if(opsxml[j].opnamekey == op)
      {
        found = true;
        keymap[op] = opsxml[j].irrep;
        break;
      }
    }

    if(!!!found)
    { 
      std::cerr << __func__ << ": ERROR: couldnt find op " << op << " exiting" << std::endl;
      exit(1);
    }

  }
  return keymap;
}

// another useful intermediary object
struct stupid_datum
{
  std::string particle_id;
  Hadron::KeyHadronNPartIrrep_t key;
  int op;
  int state;
  int t0;
};

// stream it for checking
std::ostream& operator<<(std::ostream& o, const stupid_datum &d)
{
  std::string n("\n");
  o << "particle_id " << d.particle_id << n;
  o << "key " << d.key ;
  o << "op_num " << d.op << n;
  o << "state_num " << d.state << n;
  o << "t0 " << d.t0 << n;
  return o;
}

// eat all the stuff and make a datum
std::vector<stupid_datum> chomp(const NameListXML &namelist, 
    const std::map<std::string,Hadron::KeyHadronNPartIrrep_t> &keymap, 
    const std::map<std::string,int> &opmap,
    const int t)
{

  std::vector<stupid_datum> ret;

  for(int particle = 0; particle < namelist.stateNameList.size(); ++particle)
  {
    std::string pid = namelist.stateNameList[particle];
    if(pid == "skip")
      continue;

    for(int operator_ = 0; operator_ < namelist.operatorNameList.size(); ++operator_)
    {
      stupid_datum dat;
      dat.particle_id = pid;
      if(keymap.find(namelist.operatorNameList[operator_]) != keymap.end())
        dat.key = keymap.find(namelist.operatorNameList[operator_])->second;
      else
      {
        std::cerr << __func__ << ": ERROR: something bad happened" << std::endl;
        exit(1);
      }

      if(opmap.find(namelist.operatorNameList[operator_]) != opmap.end())
        dat.op = opmap.find(namelist.operatorNameList[operator_])->second;
      else
      {
        std::cerr << __func__ << ": ERROR: something bad happened" << std::endl;
        exit(1);
      }

      dat.state = particle;
      dat.t0 = t;
      ret.push_back(dat);
    }
  }

  return ret;
}


struct dbInterface
{
  typedef SEMBLE::SembleExtendedKeyHadronNPartIrrep_t K;
  typedef SEMBLE::SembleMassOverlapData_t D;
  typedef ADATIO::SerialDBKey<K> SK;
  typedef ADATIO::SerialDBData<D> SD;

  dbInterface(void); // hide ctor
  dbInterface(const std::string &dbname) : m_dbname(dbname) , m_db(NULL) {}
  ~dbInterface(void) {if (m_db); delete m_db; m_db = NULL;}

  bool alloc(const stupid_datum &d)
  {
  //  std::stringstream ss;
  //  ss << "t0" << d.t0 << "/MassJackFiles/mass_t0_" << d.t0 << "_reorder_state" << d.state << ".jack";
  //  ENSEM::EnsemReal E;
  //  ENSEM::read(ss.str(),E);
  //  const int ncfg = E.size();
  //  std::vector<int> cfg(ncfg);
  // 
  // for(int i = 0; i < 1; ++i)
  //    cfg[i] = i;

    try
    {
      m_db = new FILEDB::AllConfStoreDB<SK,SD>(std::vector<int>(1,1));
    }
    catch(...)
    {
      std::cerr << __func__ << ": ERROR: couldn't alloc a database" << std::endl;
      delete m_db;
      m_db = NULL;
      exit(1);
    }

    if(m_db->open(m_dbname,  O_RDWR | O_TRUNC | O_CREAT, 0664) != 0)
    {
      std::cerr << __func__ << ": error opening dbase= " << m_dbname << std::endl;
      exit(1);
    }
    return true;
  }


  void insert(const stupid_datum &d)
  {
    if(!!!m_db)
      if(!!!alloc(d))
      {
        std::cerr << __func__ << ": ERROR: something bad happened" << std::endl;
        exit(1);
      }

    SK key;
    SD data;

    key.key() = K(d.particle_id,d.key);
    std::stringstream E,Z;
    E << "t0" << d.t0 << "/MassJackFiles/mass_t0_" << d.t0 << "_reorder_state" << d.state << ".jack";
    Z << "t0" << d.t0 << "/ZJackFiles/Z_t0_" << d.t0 << "_reorder_state" << d.state << "_op" << d.op << ".jack";
    ENSEM::read(E.str(),data.data().E());
    ENSEM::read(Z.str(),data.data().Z());

    if(m_db->insert(key,std::vector<SD>(1,data)) != 0)
    {
      std::cerr << __func__ << ": could not insert key \n" << key.key() << std::endl;
      exit(1);
    }
  }

  std::string m_dbname;
  FILEDB::AllConfStoreDB< SK , SD > *m_db;
};




/////////////////////////////////////////////////////////////////////////



int main(int argc, char *argv[])
{

  if(argc != 5)
  {
    std::cerr << "usage: " << argv[0] << ": <sfit.ini.xml> <t0 extract> <name_list.xml> <dbname>" << std::endl;
    exit(1);
  }

  std::string xmlinifile;
  std::istringstream val(argv[1]);
  val >> xmlinifile;
  std::cout << "Loading xmlinifile: " << xmlinifile << std::endl;

  // Read parameters from xml ini file
  SEMBLE::FitIniProps_t inikeys;

  try
  {
    ADATXML::XMLReader xml_in(xmlinifile);
    SEMBLE::read(xml_in, "/FitIniParams", inikeys);
  }
  catch(const std::string &e)
  {
    std::cerr << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }

  // sanity
  if(inikeys.dbInputType != "redstar")
  {
    std::cerr << __func__ << ": ERROR: only redstar is supported, you chose "
      << inikeys.dbInputType << ", better luck next time.." << std::endl;
    exit(1);
  }

  // Read t0 extract
  int t0 = ::atoi(argv[2]);
  std::cout << "Extracting on t0 = " << t0 << std::endl;


  // Read the name list file
  std::string namelistfile;
  std::istringstream namelistfilestream(argv[3]);
  namelistfilestream >> namelistfile; 
  std::cout << "Loading namelistfile: " << namelistfile << std::endl;

  NameListXML namelist;
  try
  {
    ADATXML::XMLReader xml_in(namelistfile);
    namelist.read(xml_in,"/NameListXML");
  }
  catch(const std::string &e)
  { 
    std::cerr << __func__ << ": ERROR: can't read namelistfile (" << namelistfile << "): " << e << std::endl;
    exit(1);
  }

  // Read the dbname
  std::string dbname;
  std::istringstream dbnamestream(argv[4]);
  dbnamestream >> dbname;


  // get the data into a useful form
  std::map<std::string, int> opmap = readOpsPhases();
  std::map<std::string, Hadron::KeyHadronNPartIrrep_t> keymap = readOpsXML(inikeys.inputPropsRedstar.opsXMLFiles, namelist);
  std::vector<stupid_datum> data = chomp(namelist, keymap, opmap , t0);
  std::vector<stupid_datum>::const_iterator it;

  std::cout << __func__  << ": extracting the following: " << std::endl;
  for(it = data.begin(); it != data.end(); ++it)
    std::cout << *it << std::endl;


  // now lets make this thing
  dbInterface dbase(dbname);
  for(it = data.begin(); it != data.end(); ++it)
    dbase.insert(*it);




  return 0;
}
