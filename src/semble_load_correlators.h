#ifndef __LOAD_CORRELATORS_H__
#define __LOAD_CORRELATORS_H__

#include "formfac/hadron_2pt_corr.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "AllConfStoreDB.h"


#include "ensem/ensem.h"
#include "semble_matrix.h"
#include "semble_algebra.h"
#include "semble_fit_ini_xml.h"
#include <vector>

using namespace std;
using namespace ENSEM;
using namespace ADATXML;
using namespace ADATIO;

namespace SEMBLE
{

//*******************************************************************
// Complex Correlators
  class SembleCCorrs
  {
  public:
    SembleCCorrs();

    void loadFromDB(const string &dbfile, vector<string> opsList_, vector<int> opsListCConj_, int irrepdim_, FF::KeyHadron2PtCorr_t DefaultKeys, const string &avgMode, double avgTol, const string &badlistfile, bool avgMom, vector< Array<int> > momList_);

    void loadFromDB(const string &dbfile, vector<string> opsList_, vector<int> opsListCConj_, const Array<string>& opsxmlfiles, Array<int> avgRows, double avgTol, const string &badlistfile, bool avgMom, vector< Array<int> > momList_, InputPropsRedstarKeys_t keyParams);

    vector<SembleMatrix<double> > rephaseCorrs(const string &rephaseMode, const string &foldTimeReversal, int tmax); //tmax is new

    SembleMatrix<std::complex<double> > getCt(int t) const;
    EnsemVectorComplex getCij(int i, int j) const;
    int getLt() const;
    int getDim() const;
    int getBins() const;

  private:
    typedef EnsemScalar<EnsemVectorComplex>::Type_t SV;

    vector<SembleMatrix<std::complex<double> > > Ct;   //correlator matrix on each timeslice

    int dim;   // Dimension of matrix to invert, i.e. number of operators
    int Lt;   // Number of time slices
    int irrepdim;   // Dimension of irrep
    int nbins;   // Number of configurations
    vector<string> opsList;   // List of operators to use
    vector<int> opsListCConj;   // List of operators' charge conjugation
    vector< Array<int> > momList;   // List of momenta to use

    // Various utils and internal functions
    Array<FF::KeyHadron2PtCorr_t> createKeys(FF::KeyHadron2PtCorr_t DefaultKeys, const string &avgMode, bool avgMom);
    Array<Hadron::KeyHadronNPartNPtCorr_t> createRedstarKeys(const Array<string>& opsxmlfiles, InputPropsRedstarKeys_t keyParams, const Array<int> avgRows, bool avgMom);

    vector<SembleMatrix<std::complex<double> > > loadCorrs(FILEDB::AllConfStoreDB< SerialDBKey<FF::KeyHadron2PtCorr_t>,  SerialDBData<SV> >& database,  Array<FF::KeyHadron2PtCorr_t> keys, const string &avgMode, double avgTol, const string &badlistfile, bool avgMom);
    vector<SembleMatrix<std::complex<double> > > loadCorrs(FILEDB::AllConfStoreDB< SerialDBKey<Hadron::KeyHadronNPartNPtCorr_t>,  SerialDBData<SV> >& database,  Array<Hadron::KeyHadronNPartNPtCorr_t> keys, const Array<int> avgRows, double avgTol, const string &badlistfile, bool avgMom);

    Array<Hadron::KeyHadronNPartIrrep_t> readOpsxml(const Array<string>& opsxmlfiles);
  };


//*******************************************************************
// Real Correlators
  class SembleRCorrs
  {
  public:
    SembleRCorrs();

    void loadRephaseComplexCorrs(const string &dbfile, const string &opslistfile, int irrepdim_, FF::KeyHadron2PtCorr_t DefaultKeys,  const string &avgMode, double avgTol, string const &badlistfile, const string &foldTimeReversal_, const string &rephaseMode_, bool avgMom, const string &momListFile);

    void loadRephaseComplexCorrs(const string &dbfile, const string &opslistfile, const Array<string>& opsxmlfiles, Array<int> avgRows, double avgTol, const string &badlistfile, const string &foldTimeReversal_, const string &rephaseMode_, const bool avgMom, const string &momListFile, InputPropsRedstarKeys_t keyParams);

    void loadFromEnsemFiles(int dim_, const string &filepath);
    void loadFromEnsemFile(const string &ensemfilename);
    void changeNumCfgs(int newbins);
    void useEnsemMean();

    void skipTimeslices(int nt);
    void shiftCorrs(int dt);
    void expWeightCorrs(double E);

    SembleMatrix<double> getCt(int t) const;
    int getLt() const;
    int getDim() const;
    int getBins() const;

    EnsemVectorReal getCij(int i, int j) const;

  private:
    typedef EnsemScalar<EnsemVectorComplex>::Type_t SV;

    vector<SembleMatrix<double> > Ct;   //correlator matrix on each timeslice

    int dim;   // Dimension of matrix to invert, i.e. number of operators
    int Lt;   // Number of time slices
    int irrepdim;   // Dimension of irrep
    int nbins;   // Number of configurations
    vector<string> opsList;   // List of operators to use
    vector<int> opsListCConj;   // List of operators' charge conjugation
    string foldTimeReversal;   // Fold time reversal option
    string rephaseMode;   // Method to use for rephasing correlators
    vector< Array<int> > momList;   // List of momenta to use

    // Various utils
    void readOpsList(const string &opslistfile);
    void readMomList(const string &momListFile);
    Array<int> stringToIntArray3(const string &input);
    void doRephasing(SembleCCorrs &ComplexCorrs);

  };

  //************************************************************************
  void loadCorr(SembleRCorrs &tp, const FitIniProps_t &inikeys);


// **************************************************************************
// Various utils

  bool phaseCorr(SembleMatrix<std::complex<double> > C, int row, int &phase_re, int  &phase_im, int refrow, int ref_phase_re, int ref_phase_im);
  bool test_done(vector<bool> in);

// opsxml input structure
  struct opsxml_t
  {
    string opnamekey;                     // operator name
    Hadron::KeyHadronNPartIrrep_t irrep;  // operator xml structure
  };
  void read(XMLReader &xml, const std::string &path, opsxml_t &param);


//*******************************************************************
// Various db utils (a lot of these are modified version of dbutils.cc)

//! Get a key/value
  template<typename K, typename V>
  V printKeyValue(const K &ky,
                  FILEDB::AllConfStoreDB< SerialDBKey<K>,  SerialDBData<typename EnsemScalar<V>::Type_t> >& database);

//cook up a promotion scheme to use with template types
///////////////////////////////////////////////////////////////////////////////
  template<class T>
  struct PromoteCorr
  {
    typedef T Type;
  };

  template<>
  struct PromoteCorr<double>
  {
    typedef SembleRCorrs Type;
  };

  template<>
  struct PromoteCorr<std::complex<double> >
  {
    typedef SembleCCorrs Type;
  };


  /////////////////////////////////////////////////////////////////////////////////////////

}

#endif
