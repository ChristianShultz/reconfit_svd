#include "semble_fit_ini_xml.h"

using namespace SEMBLE;

//a parser for sanity, ie printing z w/o fitting being true, should eliminate stupid crashes from ill-conceived ini files
void SEMBLE::check_ini(FitIniProps_t &ini)
{
  OutputProps_t out = ini.outputProps;

  //check t0ref
  ini.t0Props.t0ref = (ini.t0Props.t0ref < ini.t0Props.t0low) ? ini.t0Props.t0low + 1 : ini.t0Props.t0ref;
  ini.t0Props.t0ref = (ini.t0Props.t0ref > ini.t0Props.t0high) ? ini.t0Props.t0high : ini.t0Props.t0ref;

  //check output props
  ini.zProps.fit = (out.Z_t || out.zFitFiles || out.zFitPlots || out.zT0Files || out.zT0Plots) ? true : ini.zProps.fit;
  ini.t0FitProps.ZT0 = (out.zT0Plots || out.zT0Files) ? true : ini.t0FitProps.ZT0;
  ini.t0FitProps.MT0 = (out.mT0Files || out.mT0Plots) ? true : ini.t0FitProps.MT0;
  ini.reconProps.recon = (out.reconPlots) ? true : ini.reconProps.recon;

  //check global props
  ini.globalProps.tmax = (ini.globalProps.tmax < ini.prinCorrProps.tmax) ? ini.prinCorrProps.tmax : ini.globalProps.tmax;

  if(ini.zProps.fit)
    ini.globalProps.tmax = (ini.globalProps.tmax < ini.zProps.tmax) ? ini.zProps.tmax : ini.globalProps.tmax;

  if(ini.reconProps.recon)
    ini.globalProps.tmax = (ini.globalProps.tmax < ini.reconProps.tmax) ? ini.reconProps.tmax : ini.globalProps.tmax;

  if(out.reconPlots)
    ini.zProps.tmax = (ini.zProps.tmax < ini.reconProps.tmax) ? ini.reconProps.tmax : ini.zProps.tmax;

  //check shift
  if(ini.shiftProps.shift && ini.shiftProps.dt == 0)
    ini.shiftProps.dt = 1;
}

/* The XML readers*/

void SEMBLE::read(XMLReader &xml, const std::string &path, SortingProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("sortEvecsCfg") > 0)
    read(ptop, "sortEvecsCfg", prop.sortEvecsCfg);
  else
    prop.sortEvecsCfg = "Refvecs";

  if(ptop.count("sortEvecsTimeslice") > 0)
    read(ptop, "sortEvecsTimeslice", prop.sortEvecsTimeslice);
  else
    prop.sortEvecsTimeslice = "Refvecs";

  if(ptop.count("deltaRef") > 0)
    read(ptop, "deltaRef", prop.deltaRef);
  else
    prop.deltaRef = 1;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, GenEigProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("type") > 0)
    read(ptop, "type", prop.type);
  else
    prop.type = "SvdCond";

  if(ptop.count("type") > 0)
    read(ptop, "thresh", prop.thresh);
  else
    prop.thresh = 1e6;

  if(ptop.count("sigma") > 0)
    read(ptop, "sigma", prop.sigma);
  else
    prop.sigma = 3;

  if(ptop.count("svdHisto") > 0)
    read(ptop, "svdHisto", prop.svdHisto);
  else
    prop.svdHisto = false;

  if(ptop.count("nHistoBins") > 0)
    read(ptop, "nHistoBins", prop.nHistoBins);
  else
    prop.nHistoBins = 20;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, PrinCorrProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("tmax") > 0)
    read(ptop, "tmax", prop.tmax);
  else
    prop.tmax = 25;

  if(ptop.count("minTSlices") > 0)
    read(ptop, "minTSlices", prop.minTSlices);
  else
    prop.minTSlices = 4;

  if(ptop.count("noiseCutoff") > 0)
    read(ptop, "noiseCutoff", prop.noiseCutoff);
  else
    prop.noiseCutoff = 0.1;

  if(ptop.count("fitCrit") > 0)
    read(ptop, "fitCrit", prop.fitCrit);
  else
    prop.fitCrit = "";

  if(ptop.count("accChisq") > 0)
    read(ptop, "accChisq", prop.accChisq);
  else
    prop.accChisq = 3.0;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, ZProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("fit") > 0)
    read(ptop, "fit", prop.fit);
  else
    prop.fit = true;

  if(prop.fit)
    {
      if(ptop.count("tmax") > 0)
        read(ptop, "tmax", prop.tmax);
      else
        prop.tmax = 25;

      if(ptop.count("minTSlices") > 0)
        read(ptop, "minTSlices", prop.minTSlices);
      else
        prop.minTSlices = 4;

      if(ptop.count("fitCrit") > 0)
        read(ptop, "fitCrit", prop.fitCrit);
      else
        prop.fitCrit = "";

      if(ptop.count("accChisq") > 0)
        read(ptop, "accChisq", prop.accChisq);
      else
        prop.accChisq = 3.0;
    }
}

void SEMBLE::read(XMLReader &xml, const std::string &path, T0FitProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("ZT0") > 0)
    read(ptop, "ZT0", prop.ZT0);
  else
    prop.ZT0 = false;

  if(ptop.count("MT0") > 0)
    read(ptop, "MT0", prop.MT0);
  else
    prop.MT0 = false;

  if(ptop.count("fitCrit") > 0)
    read(ptop, "fitCrit", prop.fitCrit);
  else
    prop.fitCrit = "";

  if(ptop.count("accChisq") > 0)
    read(ptop, "accChisq", prop.accChisq);
  else
    prop.accChisq = 3.0;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, T0Props_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("t0low") > 0)
    read(ptop, "t0low", prop.t0low);
  else
    prop.t0low = 5;

  if(ptop.count("t0high") > 0)
    read(ptop, "t0high", prop.t0high);
  else
    prop.t0high = 11;

  if(ptop.count("t0ref") > 0)
    read(ptop, "t0ref", prop.t0ref);
  else
    prop.t0ref = prop.t0high;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, ReconProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("recon") > 0)
    read(ptop, "recon", prop.recon);
  else
    prop.recon = false;

  if(ptop.count("type") > 0)
    read(ptop, "type", prop.type);
  else
    prop.type = "";

  if(ptop.count("tmax") > 0)
    read(ptop, "tmax", prop.tmax);
  else
    prop.tmax = 25;

  if(ptop.count("accChisq") > 0)
    read(ptop, "accChisq", prop.accChisq);
  else
    prop.accChisq = 3.0;

  if(ptop.count("selectT0") > 0)
    read(ptop, "selectT0", prop.selectT0);
  else
    prop.selectT0 = "";
}

void SEMBLE::read(XMLReader &xml, const std::string &path, OutputProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("mass") > 0)
    read(ptop, "mass", prop.mass);
  else
    prop.mass = false;

  if(ptop.count("Z_t") > 0)
    read(ptop, "Z_t", prop.Z_t);
  else
    prop.Z_t = false;

  if(ptop.count("V_t") > 0)
    read(ptop, "V_t", prop.V_t);
  else
    prop.V_t = false;

  if(ptop.count("pcorrFiles") > 0)
    read(ptop, "pcorrFiles", prop.pcorrFiles);
  else
    prop.pcorrFiles = false;

  if(ptop.count("zFitFiles") > 0)
    read(ptop, "zFitFiles", prop.zFitFiles);
  else
    prop.zFitFiles = false;

  if(ptop.count("pcorrPlots") > 0)
    read(ptop, "pcorrPlots", prop.pcorrPlots);
  else
    prop.pcorrPlots = false;

  if(ptop.count("zFitPlots") > 0)
    read(ptop, "zFitPlots", prop.zFitPlots);
  else
    prop.zFitPlots = false;

  if(ptop.count("mT0Files") > 0)
    read(ptop, "mT0Files", prop.mT0Files);
  else
    prop.mT0Files = false;

  if(ptop.count("zT0Files") > 0)
    read(ptop, "zT0Files", prop.zT0Files);
  else
    prop.zT0Files = false;

  if(ptop.count("mT0Plots") > 0)
    read(ptop, "mT0Plots", prop.mT0Plots);
  else
    prop.mT0Plots = false;

  if(ptop.count("zT0Plots") > 0)
    read(ptop, "zT0Plots", prop.zT0Plots);
  else
    prop.zT0Plots = false;

  if(ptop.count("reconPlots") > 0)
    read(ptop, "reconPlots", prop.reconPlots);
  else
    prop.reconPlots = false;

  if(ptop.count("reconType") > 0)
    read(ptop, "reconType", prop.reconType);
  else
    prop.reconType = "diag";

  if(ptop.count("logs") > 0)
    read(ptop, "logs", prop.logs);
  else
    prop.logs = false;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, InputPropsEnsem_t &param)
{
  XMLReader paramtop(xml, path);

  if(paramtop.count("dbFname") > 0)
    read(paramtop, "dbFname", param.dbFname);
  else
    param.dbFname = "";

  if(paramtop.count("dim") > 0)
    read(paramtop, "dim", param.dim);
  else
    param.dim = -1;
}


void SEMBLE::read(XMLReader &xml, const std::string &path, InputPropsDB_t &param)
{
  XMLReader paramtop(xml, path);

  if(paramtop.count("dbFname") > 0)
    read(paramtop, "dbFname", param.dbFname);
  else
    param.dbFname = "";

  if(paramtop.count("opsListFname") > 0)
    read(paramtop, "opsListFname", param.opsListFname);
  else
    param.opsListFname = "";

  if(paramtop.count("irrepDim") > 0)
    read(paramtop, "irrepDim", param.irrepDim);
  else
    param.irrepDim = -1;

  if(paramtop.count("rephaseMode") > 0)
    read(paramtop, "rephaseMode", param.rephaseMode);
  else
    param.rephaseMode = "";

  if(paramtop.count("foldTimeReversal") > 0)
    read(paramtop, "foldTimeReversal", param.foldTimeReversal);
  else
    param.foldTimeReversal = "";

  if(paramtop.count("avgMode") > 0)
    read(paramtop, "avgMode", param.avgMode);
  else
    param.avgMode = "";

  if(paramtop.count("avgTol") > 0)
    read(paramtop, "avgTol", param.avgTol);
  else
    param.avgTol = 0;

  if(paramtop.count("badList") > 0)
    read(paramtop, "badList", param.badList);
  else
    param.badList = "badList";

  if(paramtop.count("avgMom") > 0)
    read(paramtop, "avgMom", param.avgMom);
  else
    param.avgMom = false;

  if(paramtop.count("momListFname") > 0)
    read(paramtop, "momListFname", param.momListFname);
  else
    param.momListFname = "";

  if(paramtop.count("keys") > 0)
    read(paramtop, "keys", param.keys);
}

void SEMBLE::read(XMLReader &xml, const std::string &path, InputPropsRedstar_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("dbFname") > 0)
    read(ptop, "dbFname", prop.dbFname);
  else
    prop.dbFname = "";

  if(ptop.count("opsListFname") > 0)
    read(ptop, "opsListFname", prop.opsListFname);
  else
    prop.opsListFname = "";

  if(ptop.count("opsXMLFiles") > 0)
    read(ptop, "opsXMLFiles", prop.opsXMLFiles);
  else
    prop.opsXMLFiles.resize(0);

  if(ptop.count("rephaseMode") > 0)
    read(ptop, "rephaseMode", prop.rephaseMode);
  else
    prop.rephaseMode = "";

  if(ptop.count("foldTimeReversal") > 0)
    read(ptop, "foldTimeReversal", prop.foldTimeReversal);
  else
    prop.foldTimeReversal = "";

  if(ptop.count("avdRows") > 0)
    read(ptop, "avgRows", prop.avgRows);
  else
    {
      prop.avgRows.resize(1);
      prop.avgRows[0] = 0;
    }

  if(ptop.count("avgTol") > 0)
    read(ptop, "avgTol", prop.avgTol);
  else
    prop.avgTol = 0.;

  if(ptop.count("badList") > 0)
    read(ptop, "badList", prop.badList);
  else
    prop.badList = "";

  if(ptop.count("avgMom") > 0)
    read(ptop, "avgMom", prop.avgMom);
  else
    prop.avgMom = false;

  if(ptop.count("momListFname") > 0)
    read(ptop, "momListFname", prop.momListFname);
  else
    prop.momListFname = "";

  if(ptop.count("redKeys") > 0)
    read(ptop, "redKeys", prop.redKeys);
  else
    {
      std::cout << __PRETTY_FUNCTION__ << __LINE__ << __FILE__ << "Error : redKeys required in xml in file, exiting." << std::endl;
      exit(1);
    }

}

void SEMBLE::read(XMLReader &xml, const std::string &path, InputPropsRedstarKeys_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("ensemble") > 0)
    read(ptop, "ensemble", prop.ensemble);
  else
    prop.ensemble = "";

  if(ptop.count("mom") > 0)
    read(ptop, "mom", prop.mom);
  else
    {
      prop.mom.resize(3);
      prop.mom[0] = 0;
      prop.mom[1] = 0;
      prop.mom[2] = 0;
    }

  if(ptop.count("twoIz") > 0)
    read(ptop, "twoIz", prop.twoIz);
  else
    prop.twoIz = 0;

  if(ptop.count("sourceTSlice") > 0)
    read(ptop, "sourceTSlice", prop.sourceTSlice);
  else
    prop.sourceTSlice = 0;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, ShiftProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("shift") > 0)
    read(ptop, "shift", prop.shift);
  else
    prop.shift = false;

  if(prop.shift)
    {
      if(ptop.count("dt") > 0)
        read(ptop, "dt", prop.dt);
      else
        prop.dt = 1;
    }
}

void SEMBLE::read(XMLReader &xml, const std::string &path, FixedCoeffProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("fixed") > 0)
    read(ptop, "fixed", prop.fixed);
  else
    prop.fixed = false;

  if(prop.fixed)
    if(ptop.count("tstar") > 0)
      read(ptop, "tstar", prop.tstar);
    else
      prop.tstar = 0;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, WeightProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("weight") > 0)
    read(ptop, "weight", prop.weight);
  else
    prop.weight = false;

  if(prop.weight)
    {
      if(ptop.count("E") > 0)
        read(ptop, "E", prop.E);
      else
        prop.E = 0.0;
    }
}

void SEMBLE::read(XMLReader &xml, const std::string &path, GlobalProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("tmin") > 0)
    read(ptop, "tmin", prop.tmin);
  else
    prop.tmin = 1;

  if(ptop.count("tmax") > 0)
    read(ptop, "tmax", prop.tmax);
  else
    prop.tmax = 25;

  if(ptop.count("SVCut") > 0)
    read(ptop, "SVCut", prop.SVCut);
  else
    prop.SVCut = 1e-6;

  if(ptop.count("verbose") > 0)
    read(ptop, "verbose", prop.verbose);
  else
    prop.verbose = false;
}

void SEMBLE::read(XMLReader &xml, const std::string &path, FitIniProps_t &prop)
{
  XMLReader ptop(xml, path);
  const int version = 0;

  //check version
  if(ptop.count("version") > 0)
    read(ptop, "version", prop.version);
  else
    {
      std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n Missing Version input, trying to continue \n" << std::endl;
      prop.version = version;
    }

  if(version != prop.version)
    {
      std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n" << "Error: version" << prop.version
                << " is not supported, trying to continue, supported version is " << version << std::endl;
    }

  //get db input
  if(ptop.count("dbInputType") > 0)
    read(ptop, "dbInputType", prop.dbInputType);
  else
    {
      std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n" << "Error: dbInputType required, exiting." << std::endl;
      exit(1);
    }

  //try to read db input
  if(prop.dbInputType == "ensem" || prop.dbInputType == "ensem_debug" || prop.dbInputType == "ensem_onecorr")
    {
      if(ptop.count("inputPropsEnsem") > 0)
        read(ptop, "inputPropsEnsem", prop.inputPropsEnsem);
      else
        {
          std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Error: no ensem input, exiting." << std::endl;
          exit(1);
        }
    }
  else if(prop.dbInputType == "dbnew" || prop.dbInputType == "dbnew_debug")
    {
      if(ptop.count("inputPropsDB") > 0)
        read(ptop, "inputPropsDB", prop.inputPropsDB);
      else
        {
          std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Error: no dbnew input, exiting." << std::endl;
          exit(1);
        }
    }
  else if(prop.dbInputType == "redstar" || prop.dbInputType == "redstar_debug")
    {
      if(ptop.count("inputPropsRedstar") > 0)
        read(ptop, "inputPropsRedstar", prop.inputPropsRedstar);
      else
        {
          std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Error: no redstar input, exiting." << std::endl;
          exit(1);
        }
    }
  else
    {
      std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n Error: dbInputType," << prop.dbInputType << ", not recognized,  exiting." << std::endl;
      exit(1);
    }

  //assume they are there and try to read, if not there it should default to the else values
  read(ptop, "sortingProps", prop.sortingProps);
  read(ptop, "genEigProps", prop.genEigProps);
  read(ptop, "prinCorrProps", prop.prinCorrProps);
  read(ptop, "zProps", prop.zProps);
  read(ptop, "t0FitProps", prop.t0FitProps);
  read(ptop, "t0Props", prop.t0Props);
  read(ptop, "reconProps", prop.reconProps);
  read(ptop, "outputProps", prop.outputProps);
  read(ptop, "shiftProps", prop.shiftProps);
  read(ptop, "fixedCoeffProps", prop.fixedCoeffProps);
  read(ptop, "weightProps", prop.weightProps);
  read(ptop, "globalProps", prop.globalProps);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


/*the write param functions*/
std::string SEMBLE::write_params(const SortingProps_t &prop)
{
  std::stringstream ss;

  ss << "sortingProps \n";
  ss << "sortEvecsCfg " << prop.sortEvecsCfg << "\n";
  ss << "sortEvecsTimeslice " << prop.sortEvecsTimeslice << "\n";
  ss << "deltaRef " << prop.deltaRef << "\n";
  ss << "\n\n";

  return ss.str();
}

std::string SEMBLE::write_params(const GenEigProps_t &prop)
{
  std::stringstream ss;

  ss << "genEigProps \n";
  ss << "type " << prop.type << "\n";
  ss << "thresh " << prop.thresh << "\n";
  ss << "sigma " << prop.sigma << "\n";
  ss << "histo " << prop.svdHisto << "\n";
  ss << "nHistoBins " << prop.nHistoBins << "\n";
  ss << "\n\n";

  return ss.str();
}

std::string SEMBLE::write_params(const PrinCorrProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "prinCorrProps \n";
  ss << "tmax " << prop.tmax << n;
  ss << "minTSlices " << prop.minTSlices << n;
  ss << "noiseCutoff " << prop.noiseCutoff << n;
  ss << "fitCrit " << prop.fitCrit << n;
  ss << "accChisq " << prop.accChisq << n;
  ss << "\n\n";

  return ss.str();
}

std::string SEMBLE::write_params(const ZProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "zProps \n";
  ss << "fit " << prop.fit << n;
  ss << "tmax " << prop.tmax << n;
  ss << "minTSlices " << prop.minTSlices << n;
  ss << "fitCrit " << prop.fitCrit << n;
  ss << "accChisq " << prop.accChisq << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const T0FitProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "t0FitProps" << n;
  ss << "ZT0 " << prop.ZT0 << n;
  ss << "MT0 " << prop.MT0 << n;
  ss << "fitCrit " << prop.fitCrit << n;
  ss << "accChisq " << prop.accChisq << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const T0Props_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "t0Props" << n;
  ss << "t0low " << prop.t0low << n;
  ss << "t0high " << prop.t0high << n;
  ss << "t0ref " << prop.t0ref << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const ReconProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "reconProps" << n;
  ss << "recon " << prop.recon << n;
  ss << "type " << prop.type << n;
  ss << "tmax " << prop.tmax << n;
  ss << "accChisq " << prop.accChisq << n;
  ss << "selectT0 " << prop.selectT0 << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const OutputProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "outputProps" << n;
  ss << "mass " << prop.mass << n;
  ss << "Z_t " << prop.Z_t << n;
  ss << "V_t " << prop.V_t << n;
  ss << "pcorrFiles " << prop.pcorrFiles << n;
  ss << "zFitFiles " << prop.zFitFiles << n;
  ss << "pcorrPlots " << prop.pcorrPlots << n;
  ss << "zFitPlots " << prop.zFitPlots << n;
  ss << "mT0Files " << prop.mT0Files << n;
  ss << "zT0Files " << prop.zT0Files << n;
  ss << "mT0Plots " << prop.mT0Plots << n;
  ss << "zT0Plots " << prop.zT0Plots << n;
  ss << "reconPlots " << prop.reconPlots << n;
  ss << "reconType " << prop.reconType << n;
  ss << "logs " << prop.logs << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const InputPropsEnsem_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "inputPropsEnsem" << n;
  ss << "dbFname " << prop.dbFname << n;
  ss << "dim " << prop.dim << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const InputPropsDB_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "inputPropsDB" << n;
  ss << "dbFname " << prop.dbFname << n;
  ss << "opsListFname " << prop.opsListFname << n;
  ss << "irrepDim " << prop.irrepDim << n;
  ss << "rephaseMode " << prop.rephaseMode << n;
  ss << "foldTimeReversal " << prop.foldTimeReversal << n;
  ss << "avgMode " << prop.avgMode << n;
  ss << "avgTol " << prop.avgTol << n;
  ss << "badList " << prop.badList << n;
  ss << "avgMom " << prop.avgMom << n;
  ss << "momListFname " << prop.momListFname << n;
  ss << "keys -- no printing" << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const InputPropsRedstar_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "inputPropsRedstar" << n;
  ss << "dbFname " << prop.dbFname << n;
  ss << "opsListFname " << prop.opsListFname << n;
  ss << "opsXMLFiles -- no printing support" << n;
  ss << "rephaseMode " << prop.rephaseMode << n;
  ss << "foldTimeReversal " << prop.foldTimeReversal << n;
  ss << "avgRows -- no printing support" << n;
  ss << "avgTol " << prop.avgTol << n;
  ss << "badList " << prop.badList << n;
  ss << "avgMom " << prop.avgMom << n;
  ss << "momListFname " << prop.momListFname << n;
  ss << "redKeys " << prop.redKeys << n;
  ss << n << n;

  return ss.str();
}
std::string SEMBLE::write_params(const InputPropsRedstarKeys_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "inputPropsRedstarKeys" << n;
  ss << "ensemble " << prop.ensemble << n;
  ss << "mom -- no printing support " << n;
  ss << "twoIZ " << prop.twoIz << n;
  ss << "sourceTSlice " << prop.sourceTSlice << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const ShiftProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "shiftProps" << n;
  ss << "shift " << prop.shift << n;
  ss << "dt " << prop.dt;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const FixedCoeffProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "fixedCoeffProps" << n;
  ss << "fixed " << prop.fixed << n;
  ss << "tstar " << prop.tstar << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const WeightProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "weightProps" << n;
  ss << "weight " << prop.weight << n;
  ss << "E " << prop.E << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const GlobalProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "globalProps" << n;
  ss << "tmin " << prop.tmin << n;
  ss << "tmax " << prop.tmax << n;
  ss << "SVCut " << prop.SVCut << n;
  ss << "verbose " << prop.verbose << n;
  ss << n << n;

  return ss.str();
}

std::string SEMBLE::write_params(const FitIniProps_t &prop)
{
  std::stringstream ss;
  std::string n("\n");

  ss << "fitIniProps" << n;
  ss << "version " << prop.version << n;
  ss << "dbInputType " << prop.dbInputType << n << n;
  ss << "sortingProps " << prop.sortingProps << n;
  ss << "genEigProps " << prop.genEigProps << n;
  ss << "prinCorrProps " << prop.prinCorrProps << n;
  ss << "zProps " << prop.zProps << n;
  ss << "t0FitProps " << prop.t0FitProps << n;
  ss << "t0Props " << prop.t0Props << n;
  ss << "reconProps " << prop.reconProps << n;
  ss << "outputProps " << prop.outputProps << n;
  ss << "inputPropsEnsem " << prop.inputPropsEnsem << n;
  ss << "inputPropsDB " << prop.inputPropsDB << n;
  ss << "inputPropsRedstar " << prop.inputPropsRedstar << n;
  ss << "shiftProps " << prop.shiftProps << n;
  ss << "fixedCoeffProps " << prop.fixedCoeffProps << n;
  ss << "weightProps " << prop.weightProps << n;
  ss << "globalProps " << prop.globalProps;
  ss << n << n;

  return ss.str();
}



//overload <<
std::ostream &SEMBLE::operator<<(std::ostream &o, const SortingProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const GenEigProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const PrinCorrProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const ZProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const T0FitProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const T0Props_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const ReconProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const OutputProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const InputPropsEnsem_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const InputPropsDB_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const InputPropsRedstar_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const InputPropsRedstarKeys_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const ShiftProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const FixedCoeffProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const WeightProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const GlobalProps_t &prop)
{
  o << write_params(prop);
  return o;
}

std::ostream &SEMBLE::operator<<(std::ostream &o, const FitIniProps_t &prop)
{
  o << write_params(prop);
  return o;
}
