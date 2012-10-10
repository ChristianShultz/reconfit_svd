/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : recondb_util.cc

 * Purpose :

 * Creation Date : 05-10-2012

 * Last Modified : Wed Oct 10 11:47:08 2012

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


typedef SEMBLE::SembleExtendedKeyHadronNPartIrrep_t K;
typedef SEMBLE::SembleMassOverlapData_t D;
typedef ADATIO::SerialDBKey<K> SK;
typedef ADATIO::SerialDBData<D> SD;


template<typename KEY, typename DATA> 
  typename std::vector< ADATIO::SerialDBKey<KEY> >
keys(typename FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > &db)
{
  typename std::vector< ADATIO::SerialDBKey<KEY> > skeys;
  db.keys(skeys);
  return skeys;
}

template<typename KEY, typename DATA> 
  void
print_keys(typename FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > &db)
{
  typename std::vector< ADATIO::SerialDBKey<KEY> > m_keys = keys(db);
  typename std::vector< ADATIO::SerialDBKey<KEY> >::const_iterator it;

  for(it = m_keys.begin(); it != m_keys.end(); ++it)
    std::cout << it->key() << std::endl;
}

  template<typename KEY, typename DATA>
int work_handler(const std::string &dbname, const std::string &op)
{

  FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > db;
  if(db.open(dbname, O_RDONLY , 0400) != 0)
  {
    std::cerr << __func__ << ": error opening database " << dbname << std::endl;
    exit(1);
  }

  if(op == "keys")
    print_keys(db);

  return 0;
}


int main(int argc , char *argv[])
{

  if(argc != 3)
  {
    std::cerr << "usage; " << argv[0] << ": <database> <operation>" << std::endl;
    exit(1);
  }

  // Read the dbname
  std::string dbname;
  std::istringstream dbnamestream(argv[1]);
  dbnamestream >> dbname;

  std::string op;
  std::istringstream opstream(argv[2]);
  opstream >> op;
  return work_handler<K,D>(dbname,op);
} 





