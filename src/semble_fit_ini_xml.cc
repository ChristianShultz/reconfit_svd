#include "semble_fit_ini_xml.h"
#include "io/adat_xml_group_reader.h"

using namespace SEMBLE;

//a parser for sanity, ie printing z w/o fitting being true, should eliminate stupid crashes from ill-conceived ini files
void SEMBLE::check_ini(FitIniProps_t &ini)
{

  OutputProps_t out = ini.outputProps;

  //check t0ref
  ini.t0Props.t0ref = (ini.t0Props.t0ref < ini.t0Props.t0low) ? ini.t0Props.t0low + 1 : ini.t0Props.t0ref;
  ini.t0Props.t0ref = (ini.t0Props.t0ref > ini.t0Props.t0high) ? ini.t0Props.t0high : ini.t0Props.t0ref;

  //check output props

  // NB: give full control of z fitting to <zfit>
  if(!!!ini.zProps.fit)
    {

      if(out.Z_t || out.zFitFiles || out.zFitPlots || out.zT0Files || out.zT0Plots)
	std::cout << "Notification:  zFit was set to false so no Z output files will be produced" << std::endl;
	
      ini.outputProps.Z_t = false;
      ini.outputProps.zFitFiles = false;
      ini.outputProps.zFitPlots = false;
      ini.outputProps.zT0Files = false;
      ini.outputProps.zT0Plots = false;
      ini.t0FitProps.ZT0 = false;
    } 

  //  ini.zProps.fit = (out.Z_t || out.zFitFiles || out.zFitPlots || out.zT0Files || out.zT0Plots) ? true : ini.zProps.fit;
  ini.t0FitProps.ZT0 = (out.zT0Plots || out.zT0Files) ? true : ini.t0FitProps.ZT0;
  ini.t0FitProps.MT0 = (out.mT0Files || out.mT0Plots) ? true : ini.t0FitProps.MT0;
  ini.reconProps.recon = (out.reconPlots) ? true : ini.reconProps.recon;

  //check global props
  ini.globalProps.tmax = (ini.globalProps.tmax < ini.prinCorrProps.tmax) ? ini.prinCorrProps.tmax : ini.globalProps.tmax;



  if(ini.zProps.fit)
     ini.globalProps.tmax = (ini.globalProps.tmax < ini.zProps.tmax) ? ini.zProps.tmax : ini.globalProps.tmax;


  if(ini.reconProps.recon)
    ini.globalProps.tmax = (ini.globalProps.tmax < ini.reconProps.tmax) ? ini.reconProps.tmax : ini.globalProps.tmax;

  //if(out.reconPlots)
  //  ini.zProps.tmax = (ini.zProps.tmax < ini.reconProps.tmax) ? ini.reconProps.tmax : ini.zProps.tmax;

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


      if(ptop.count("tmax") > 0)
        read(ptop, "tmax", prop.tmax);
      else
        prop.tmax = 25;

      //////////////////////////////////////
      //std::cout << "read Zprops.tmax = " << prop.tmax << std::endl;
      ////////////////////////////////

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

  if(ptop.count("nThreads") > 0)
    read(ptop, "nThreads", prop.nThreads);
  else
    prop.nThreads = -1; // will not try to set the number of threads

  //  if(ptop.count("svdPhop") > 0)
  //    read(ptop,"svdPhop", prop.svdPhop);
  //  else
  //    prop.svdPhop = 0; //no hop
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
  {
    if(ptop.count("tstar") > 0)
      read(ptop, "tstar", prop.tstar);
    else
      prop.tstar = 0;
  }
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





void SEMBLE::read(XMLReader &xml, const std::string &path, WeightShiftCorrectProps_t &prop)
{
  XMLReader ptop(xml, path);

  if(ptop.count("E_dt") > 0){
    
    Array<GroupXML_t> params = readXMLArrayGroup(ptop, "E_dt", "weight_energy");
    for(int p = 0; p < params.size(); p++){
      std::pair<double, int> tmp;
      
      const GroupXML_t& param_xml = params[p];
      std::istringstream xml_s(param_xml.xml);
      XMLReader paramtop(xml_s);
      
      read( paramtop, "weight_energy", tmp.first );
      read( paramtop, "shift_tslices", tmp.second );
      
      prop.E_dt.push_back(tmp); 
      
    }
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

  if(ptop.count("skip_nt") > 0)
    read(ptop, "skip_nt", prop.skip_nt);
  else
    prop.skip_nt = 0;

}







void SEMBLE::read(XMLReader &xml, const std::string &path, FitIniProps_t &prop)
{
  XMLReader ptop(xml, path);
  const int version = 2;

  //check version
  read(ptop, "version", prop.version);

  // Read
  switch (version) 
  {
  case 2:
    break;
    
  default:
    std::cerr << "Input parameter version " << version << " unsupported." << std::endl;
    exit(1);
  }

  //get db input
  prop.inputProps  = readXMLGroup(ptop, "inputProps", "dbInputType");

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
  read(ptop, "weightShiftCorrectProps", prop.weightShiftCorrectProps);
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
  ss << "nThreads" << prop.nThreads << n;
  //ss << "svdPhop " << prop.svdPhop << n;
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
  ss << "skip_nt " << prop.skip_nt << n;
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
  ss << "inputProps " << prop.inputProps.xml << n;
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
