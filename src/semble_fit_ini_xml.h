#ifndef __FIT_INI_XML_H__
#define __FIT_INI_XML_H__

#include "formfac/hadron_2pt_corr.h"
#include "io/adat_xmlio.h"

#include <iostream>
#include <string>
#include <sstream>

using namespace ADATXML;
using namespace ADATIO;


namespace SEMBLE
{
  //forward
  struct SortingProps_t;
  struct GenEigProps_t;
  struct PrinCorrProps_t;
  struct ZProps_t;
  struct T0FitProps_t;
  struct T0Props_t;
  struct ReconProps_t;
  struct OutputProps_t;
  struct InputPropsEnsem_t;
  struct InputPropsDB_t;
  struct InputPropsRedstar_t;
  struct InputPropsRedstarKeys_t;
  struct ShiftProps_t;
  struct WeightProps_t;
  struct FixedCoeffProps_t;
  struct GlobalProps_t;
  struct FitIniProps_t;

  //a parser for sanity, ie printing z w/o fitting being true, should eliminate stupid crashes from ill-conceived ini files
  void check_ini(FitIniProps_t &ini);

  //the associated readers
  void read(XMLReader &xml, const std::string &path, SortingProps_t &prop);
  void read(XMLReader &xml, const std::string &path, GenEigProps_t &prop);
  void read(XMLReader &xml, const std::string &path, PrinCorrProps_t &prop);
  void read(XMLReader &xml, const std::string &path, ZProps_t &prop);
  void read(XMLReader &xml, const std::string &path, T0FitProps_t &prop);
  void read(XMLReader &xml, const std::string &path, T0Props_t &prop);
  void read(XMLReader &xml, const std::string &path, ReconProps_t &prop);
  void read(XMLReader &xml, const std::string &path, OutputProps_t &prop);
  void read(XMLReader &xml, const std::string &path, InputPropsEnsem_t &prop);
  void read(XMLReader &xml, const std::string &path, InputPropsDB_t &prop);
  void read(XMLReader &xml, const std::string &path, InputPropsRedstar_t &prop);
  void read(XMLReader &xml, const std::string &path, InputPropsRedstarKeys_t &prop);
  void read(XMLReader &xml, const std::string &path, ShiftProps_t &prop);
  void read(XMLReader &xml, const std::string &path, FixedCoeffProps_t &prop);
  void read(XMLReader &xml, const std::string &path, WeightProps_t &prop);
  void read(XMLReader &xml, const std::string &path, GlobalProps_t &prop);
  void read(XMLReader &xml, const std::string &path, FitIniProps_t &prop);

  //the associated writers
  std::string write_params(const SortingProps_t &prop);
  std::string write_params(const GenEigProps_t &prop);
  std::string write_params(const PrinCorrProps_t &prop);
  std::string write_params(const ZProps_t &prop);
  std::string write_params(const T0FitProps_t &prop);
  std::string write_params(const T0Props_t &prop);
  std::string write_params(const ReconProps_t &prop);
  std::string write_params(const OutputProps_t &prop);
  std::string write_params(const InputPropsEnsem_t &prop);
  std::string write_params(const InputPropsDB_t &prop);
  std::string write_params(const InputPropsRedstar_t &prop);
  std::string write_params(const InputPropsRedstarKeys_t &prop);
  std::string write_params(const ShiftProps_t &prop);
  std::string write_params(const FixedCoeffProps_t &prop);
  std::string write_params(const WeightProps_t &prop);
  std::string write_params(const GlobalProps_t &prop);
  std::string write_params(const FitIniProps_t &prop);

  //overload <<
  std::ostream &operator<<(std::ostream &o, const SortingProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const GenEigProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const PrinCorrProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const ZProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const T0FitProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const T0Props_t &prop);
  std::ostream &operator<<(std::ostream &o, const ReconProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const OutputProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const InputPropsEnsem_t &prop);
  std::ostream &operator<<(std::ostream &o, const InputPropsDB_t &prop);
  std::ostream &operator<<(std::ostream &o, const InputPropsRedstar_t &prop);
  std::ostream &operator<<(std::ostream &o, const InputPropsRedstarKeys_t &prop);
  std::ostream &operator<<(std::ostream &o, const ShiftProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const FixedCoeffProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const WeightProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const GlobalProps_t &prop);
  std::ostream &operator<<(std::ostream &o, const FitIniProps_t &prop);


  //the actual structs
  struct SortingProps_t
  {
    std::string sortEvecsCfg;
    std::string sortEvecsTimeslice;
    int deltaRef;
  };

  struct GenEigProps_t
  {
    std::string type;
    double thresh;
    double sigma;
    bool svdHisto;
    int nHistoBins;
  };

  struct PrinCorrProps_t
  {
    int tmax;
    int minTSlices;
    double noiseCutoff;
    std::string fitCrit;
    double accChisq;
  };

  struct ZProps_t
  {
    bool fit;
    int tmax;
    int minTSlices;
    std::string fitCrit;
    double accChisq;
  };

  struct T0FitProps_t
  {
    bool ZT0;
    bool MT0;
    std::string fitCrit;
    double accChisq;
  };

  struct T0Props_t
  {
    int t0low;
    int t0high;
    int t0ref;
  };

  struct ReconProps_t
  {
    bool recon;
    std::string type;
    int tmax;
    double accChisq;
    std::string selectT0;
  };

  struct OutputProps_t
  {
    bool mass;
    bool Z_t;
    bool V_t;
    bool pcorrFiles;
    bool zFitFiles;
    bool pcorrPlots;
    bool zFitPlots;
    bool mT0Files;
    bool zT0Files;
    bool mT0Plots;
    bool zT0Plots;
    bool reconPlots;
    std::string reconType;
    bool logs;
  };

  struct InputPropsEnsem_t
  {
    std::string dbFname;
    int dim;
  };

  struct InputPropsDB_t
  {
    std::string dbFname;
    std::string opsListFname;
    int irrepDim;
    std::string rephaseMode;
    std::string foldTimeReversal;
    std::string avgMode;
    double avgTol;
    std::string badList;
    bool avgMom;
    std::string momListFname;
    FF::KeyHadron2PtCorr_t keys;
  };

  struct InputPropsRedstarKeys_t
  {
    std::string ensemble;
    Array<int> mom;
    int twoIz;
    int sourceTSlice;
  };

  struct InputPropsRedstar_t
  {
    std::string dbFname;
    std::string opsListFname;
    Array<std::string> opsXMLFiles;
    std::string rephaseMode;
    std::string foldTimeReversal;
    Array<int> avgRows;
    double avgTol;
    std::string badList;
    bool avgMom;
    std::string momListFname;
    InputPropsRedstarKeys_t redKeys;
  };

  struct ShiftProps_t
  {
    bool shift;
    int dt;
  };

  struct FixedCoeffProps_t
  {
    bool fixed;
    int tstar;
  };

  struct WeightProps_t
  {
    bool weight;
    double E;
  };

  struct GlobalProps_t
  {
    int tmin;
    int tmax;
    double SVCut;
    bool verbose;
  };

  struct FitIniProps_t
  {
    int version;
    std::string dbInputType;
    SortingProps_t sortingProps;
    GenEigProps_t genEigProps;
    PrinCorrProps_t prinCorrProps;
    ZProps_t zProps;
    T0FitProps_t t0FitProps;
    T0Props_t t0Props;
    ReconProps_t reconProps;
    OutputProps_t outputProps;
    InputPropsEnsem_t inputPropsEnsem;
    InputPropsDB_t inputPropsDB;
    InputPropsRedstar_t inputPropsRedstar;
    ShiftProps_t shiftProps;
    FixedCoeffProps_t fixedCoeffProps;
    WeightProps_t weightProps;
    GlobalProps_t globalProps;
  };


}//end namespace

#endif
