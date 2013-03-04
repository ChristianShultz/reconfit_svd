#ifndef SEMBLE_T0FIT_H_H_GUARD
#define SEMBLE_T0FIT_H_H_GUARD

#include"ensem/ensem.h"
#include"adat/handle.h"
#include"semble_fit_ini_xml.h"
#include"jackFitter/fit_correlators.h"
#include"jackFitter/fit_forms.h"
#include"jackFitter/jackknife_fitter.h"
#include"jackFitter/plot.h"
#include"semble/semble_file_management.h"
#include"semble/semble_matrix.h"
#include"semble/semble_vector.h"
#include"semble_t0fit_prim.h"
#include"semble_load_correlators.h"
#include"semble_functions.h"
#include"semble/semble_meta.h"
#include<map>
#include<vector>
#include<iostream>
#include<string>
#include<utility>

/*
  A class that solves the GEVP for a single value of t0 then uses the 
  generalized eigenvectors and generalized eigenvalues to solve for the 
  spectrum and the overlaps (Z)
  
  The reorder maps used in multi-t0 fitting obey the following conventions
  1) map keys are the number assigned to the ordered states, data is the 
     state on this t0
  2) When using SVD we could potentially have different numbers of 
     states on each t0, the sorting of any states with index larger than 
     dim(t0_ref) is arbitrary
*/


namespace SEMBLE
{

  template<class T>
  struct ST0Fit
  {
  public: // constructor, destructor, copy, assignment
    ST0Fit(void);
    ST0Fit(int t0_, const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_);
    ST0Fit(int t0_, const typename std::vector<SembleMatrix<T> > &fake, const FitIniProps_t &inikeys_);
    ST0Fit(const ST0Fit<T> &o);
    ~ST0Fit(void) {} //nothing fancy
    ST0Fit<T>& operator=(const ST0Fit<T> &o);


  public: // load methods
    void load(int t0_, const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_);
    void load(int t0_, const typename std::vector<SembleMatrix<T> > &fake, const FitIniProps_t &inikeys_);

    // convert from a typename PromoteCorr<T>::Type to a std::vector<SembleMatrix<T> >
    typename std::vector<SembleMatrix<T> > maketp(const typename PromoteCorr<T>::Type &tp_);

  public:
    //principal correlator methods
    void fitPrinCorr(void);
    void printSpectrum(std::ostream &out) const;                           // print the fit data (m,m',A) to an ostream
    void printPrinCorrFits(void) const;                                    // print the unordered plots to a file
    void printPrinCorrReorder(const std::map<int, int> &mapp) const;       // print the reordered plots to a file (multit0 fitting)
    void printMassFile(void) const;                                        // print the unordered jack files
    void printMassFileReorder(const std::map<int, int> &mapp) const;       // print the reordered jack files
    void printPrinCorrFiles(void) const;                                   // print the unordered pcorr files
    void printPrinCorrFilesReorder(const std::map<int, int> &mapp) const;  // print the reordered pcorr files
    void printPrinCorrFitLog(void) const;                                  // print the unordered pcorr fit log
    void printPrinCorrFitLogReorder(const std::map<int,int> &mapp) const;  // print the reordered pcorr fit log

  public:
    //overlap methods
    void makeZ(void);                                                      // construct overlaps
    void fitZ(void);                                                       // fit overlaps
    void printZFits(void) const;                                           // unordered plot
    void printZFitsReorder(const std::map<int, int> &mapp) const;          // reordered plot
    void printZFile(void) const;                                           // unordered jack file
    void printZFileReorder(const std::map<int, int> &mapp) const;          // reordered jack file
    void printZt(void) const;                                              // print the z at each time
    void printZtReorder(const std::map<int, int> &mapp) const;             // print the reordered z at each time
    void printZFitLog(void) const;                                         // print the unordered fit log
    void printZFitLogReorder(const std::map<int,int> &mapp) const;         // print the reordered fit log

  public:
    //eigenvector methods
    void printVt(void) const;                                               // print the eigenvectors at each time
    void printVtReorder(const std::map<int, int> &mapp) const;              // print the reordered eigenvectors at each time

  public:
    //reconstruction methods
    void reconCorr(void);                                                   // do the reconstruction
    double reconChisq(const int tZ, const ReconProps_t &rec);               // this is very slow
    double reconChisqFast(const int tz, const ReconProps_t &rec);           // chisq w/o inverted covariance
    void clearRecon(void)                                                   // clean up the big recon obj
    {
      recon_total.clear();
      recon = false;
    }
    pair<int, double> findBestTZ(const ReconProps_t &rec);                  // find the best set of overlaps by reconstruction
    void makeReconPlots(std::string &mode);
    //void makeReconPlotsReorder(cost std::map<int,int> &mapp);
    void printReconPlots(void);
    void clearReconPlots(void)
    {
      recon_plots.clear();
      recon_flat = false;
    }

  public:// return the fit log
    std::string getPCorrFitLog(const std::map<int,int> &mapp);              // return a string of the fit log in the old reconfit
                                                                            // format, the map is the reorder map for indexing 
  public:
    //fixed coefficient method
    void fixedCoeffMethod(const int tstar) const;

  public:
    //rephase the states for multit0fitting -- key is state, data is  new phase
    void rephase(const typename std::map<int, typename PromoteScalar<T>::Type> &mapp);  

    //peek methods
  public:
    const SembleVector<double>& peekVals(int i) const
    {
      initChk();
      return eVals[i];
    }
    const SembleMatrix<T>&  peekVecs(int i) const
    {
      initChk();
      return eVecs[i];
    }
    const SembleMatrix<T>& peekZ(int i) const
    {
      zinitChk();  //look at on at a single time
      return Z[i];
    }
    const EnsemReal& peekMass0(int i) const
    {
      pcorrChk();
      return mass_0[i];
    }
    const EnsemReal& peekMass1(int i) const
    {
      pcorrChk();
      return mass_1[i];
    }
    const EnsemReal& peekA(int i) const
    {
      pcorrChk();
      return A[i];
    }
    typename PromoteEnsem<T>::Type const peekZ(int s, int o) const
    {
      zfitChk();  //pull the ensemble from the fit
      return zFit.getEnsemElement(s, o);          
    }
    const SembleMatrix<T>& getZFit(void) const
    {
      zfitChk();
      return zFit;
    }
    const SembleMatrix<T>& getCt0(void) const
    {
      initChk();
      return Ct0;
    }
    int getN(void) const
    {
      initChk();
      return N;
    }
    int getM(void) const
    {
      initChk();
      return M;
    }
    int getB(void) const
    {
      initChk();
      return B;
    }

    std::pair<double,double> getCondNum(void) const
    {
      initChk();
      return condNum;
    }

    double getPCorrChiSq(int i) const
    {
      pcorrChk();
      return double(pCorrChiSqPDoF[i]);
    }
    double getZFitChiSq(int state, int op) const
    {
      zfitChk();
      return double(zFitChiSqPDof[state][op]);
    }

  private: //bool checks;
    void initChk(void) const;
    void pcorrChk(void) const;
    void zfitChk(void) const;
    void zinitChk(void) const;

  private: //internal functions
    Handle<FitComparator> fitComp(std::string &s) const;
    void cleanUp(void);            //clean up to load a new problem

  private: //data store
    int N, M, B, t0, tz_best;               //nrows,ncols,nbins,t0
    FitIniProps_t inikeys;
    SembleMatrix<T> Ct0;
    ST0FitPrim<T> primFit;
    std::vector<SembleVector<double> > eVals;
    typename std::vector<SembleMatrix<T> > eVecs, Z;
    std::pair<double,double> condNum;       //mean,var

    //init bools
    bool init, pcorr, zfit, zinit, recon, recon_flat, recon_tz;

    //pcorr data
    std::vector<EnsemReal> mass_0, mass_1, A;
    std::stringstream pCorrLog;
    std::vector<std::string> pCorrFitName, pCorrFitPlot, pCorrFitLog;
    std::vector<double> pCorrChiSqPDoF;

    //overlap data
    std::vector<std::vector<std::string> > zFitPlot, zFitResults;
    std::vector<std::vector<double> > zFitChiSqPDof;
    SembleMatrix<T> zFit;

    //reconstruction data
    std::map<int, typename std::vector<SembleMatrix<T> > > recon_total;         // int gives the set of Z used to construct,
    EnsemData flatData;                                                         //the inner vector is the timeslice, massive object
    std::map<std::pair<int, int>, std::pair<std::string, std::string> > recon_plots;           //key is _ij, data is fname and plot string

    //two point correlator matrix at each time slice
    typename std::vector<SembleMatrix<T> > tp;

  };

  template<class T>
  ST0Fit<T>::ST0Fit(void)
  {
    cleanUp();
  }

  template<class T>
  ST0Fit<T>::ST0Fit(int t0_, const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_)
  {
    load(t0_, tp_, inikeys_);
  }

  template<class T>
  ST0Fit<T>::ST0Fit(int t0_, const typename std::vector<SembleMatrix<T> > &fake, const FitIniProps_t &inikeys_)
  {
    load(t0_, fake, inikeys_);
  }

  template<class T>
  ST0Fit<T>::ST0Fit(const ST0Fit<T> &o)
    : N(o.N) , M(o.M) , B(o.B) , t0(o.t0) 
  {
    init = o.init;
    pcorr = o.pcorr;
    zinit = o.zinit;
    zfit = o.zfit;
    recon = o.recon;
    recon_flat = o.recon_flat;
    recon_tz = o.recon_tz;

    if(o.init)
      {
        init = true;
        primFit = o.primFit;
        inikeys = o.inikeys;
        Ct0 = o.Ct0;
        eVals = o.eVals;
        eVecs = o.eVecs;
        tp = o.tp;
      }

    if(o.pcorr)
      {
        mass_0 = o.mass_0;
        mass_1 = o.mass_1;
        A = o.A;
        //pCorrLog = o.pCorrLog;    //can't copy a buffer
        pCorrFitName = o.pCorrFitName;
        pCorrFitPlot = o.pCorrFitPlot;
        pCorrChiSqPDoF = o.pCorrChiSqPDoF;
        pcorr = true;
      }

    if(o.zinit)
      {
        Z = o.Z;
      }

    if(o.zfit)
      {
        zFitPlot = o.zFitPlot;
        zFitResults = o.zFitResults;
        zFitChiSqPDof = o.zFitChiSqPDof;
        zfit = true;
      }

    if(o.recon)
      {
        recon_total = o.recon_total;
        recon = true;
      }

    if(o.recon_flat)
      {
	recon_flat = true;
	flatData = o.flatData;
      }

    if(o.recon_tz)
      tz_best = o.tz_best;
  }

  template<class T>
  ST0Fit<T>& ST0Fit<T>::operator=(const ST0Fit<T> &o)
  {
    if(this != &o)
      {
        cleanUp();

        N = o.N;
        M = o.M;
        B = o.B;
        t0 = o.t0;


        init = o.init;
        pcorr = o.pcorr;
        zinit = o.zinit;
        zfit = o.zfit;
        recon = o.recon;
        recon_flat = o.recon_flat;
        recon_tz = o.recon_tz;

        if(o.init)
          {
            init = true;
            primFit = o.primFit;
            inikeys = o.inikeys;
            Ct0 = o.Ct0;
            eVals = o.eVals;
            eVecs = o.eVecs;
            tp = o.tp;
          }

        if(o.pcorr)
          {
            mass_0 = o.mass_0;
            mass_1 = o.mass_1;
            A = o.A;
            //pCorrLog = o.pCorrLog;  //can't copy the buffer
            pCorrFitName = o.pCorrFitName;
            pCorrFitPlot = o.pCorrFitPlot;
            pCorrChiSqPDoF = o.pCorrChiSqPDoF;
            pcorr = true;
          }

        if(o.zinit)
          {
            Z = o.Z;
          }

        if(o.zfit)
          {
            zFitPlot = o.zFitPlot;
            zFitResults = o.zFitResults;
            zFitChiSqPDof = o.zFitChiSqPDof;
            zfit = true;
          }

        if(o.recon)
          {
            recon_total = o.recon_total;
            recon = true;
          }

        if(o.recon_flat)
          {
	    recon_flat = true;
	    flatData = o.flatData;
	  }

        if(o.recon_tz)
          tz_best = o.tz_best;

      }

    return *this;
  }

  template<class T>
  void ST0Fit<T>::load(int t0_, const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_)
  {
    load(t0_,maketp(tp_),inikeys_);
  }

  template<class T>
  void ST0Fit<T>::load(int t0_, const typename std::vector<SembleMatrix<T> > &fake, const FitIniProps_t &inikeys_)
  {
    cleanUp();
    t0 = t0_;
    inikeys = inikeys_;
    tp = fake;
    primFit.clear();
    primFit.load(t0_, fake, inikeys_);
    primFit.solve();
    primFit.sort_solved();
    Ct0 = primFit.getCt0();
    N = primFit.getN();
    M = primFit.getM();
    B = primFit.getB();
    eVals = primFit.getEvals();
    eVecs = primFit.getEvecs();
    primFit.clear();
    init = true;

    SembleMatrix<T> U, V;
    SembleVector<double> s;
    
    svd(Ct0,U,s,V);
    condNum = std::pair<double,double>(toScalar(mean(s(0)/s(N-1))),std::sqrt(toScalar(variance(s(0)/s(N-1)))));

    std::cout << "ST0Fit constructed, Nbins = " << B << ", colRank = " << M << ", rowRank = " << N << ", t0 = " << t0 << std::endl;
  }

  template<class T>
  typename std::vector<SembleMatrix<T> > ST0Fit<T>::maketp(const typename PromoteCorr<T>::Type &tp_)
  {
    typename std::vector<SembleMatrix<T> > dum;

    for(int t = 0; t <= inikeys.globalProps.tmax; ++t)
      {
        dum.push_back(tp_.getCt(t));
      }

    return dum;
  }


//fit the principal correlators
  template<class T>
  void ST0Fit<T>::fitPrinCorr(void)
  {
    initChk();

    SembleVector<double> dum(B, 1);
    EnsemReal one,zero;

    dum.ones();
    one = dum.getEnsemElement(0);
    dum.zeros();
    zero = dum.getEnsemElement(0);

    //clean out anything thats old
    mass_0.clear();
    mass_1.clear();
    A.clear();
    pCorrFitName.clear();
    pCorrFitPlot.clear();
    pCorrChiSqPDoF.clear();

    //reserve the space we need
    mass_0.resize(M, zero);
    mass_1.resize(M, zero);
    A.resize(M, zero);
    pCorrFitName.resize(M);
    pCorrFitPlot.resize(M);
    pCorrChiSqPDoF.resize(M);

    std::cout << "Constructing Principal Correlator Fits.." << std::endl;
    pCorrLog << "Constructing Principal Correlator Fits for t0 = " << t0 << " rowRank(Ct0) = " << N << " colRank(Ct0) = " << M << std::endl;
    std::stringstream pfitlog;

    //loop on states
    for(int state = 0; state < M; ++state)
      {

	pfitlog.str(std::string());
	pfitlog << "fitting unordered state #" << state << std::endl;
        pCorrLog << "fitting unordered state #" << state << std::endl;

        EnsemVectorReal lambda;
        std::vector<double> tslice;
        lambda.resize(B);
        lambda.resizeObs(inikeys.prinCorrProps.tmax + 1 - inikeys.globalProps.tmin);

        for(int t = inikeys.globalProps.tmin; t <= inikeys.prinCorrProps.tmax; ++t)
          {
            pokeObs(lambda, eVals[t].getEnsemElement(state), t - inikeys.globalProps.tmin);
            tslice.push_back(t);
          }

        EnsemData pCorrData(tslice, lambda);
        Handle<FitComparator> fitC = fitComp(inikeys.prinCorrProps.fitCrit);
        FitPrincipalCorrelator pCorrFit(pCorrData, t0, fitC, inikeys.prinCorrProps.noiseCutoff, inikeys.prinCorrProps.minTSlices);

        mass_0[state] = pCorrFit.getMass0();

        if(toScalar(mean(mass_0[state])) <= 0.0)
          mass_0[state] = one;

        if(pCorrFit.getNExp() == 2)
          {
            mass_1[state] = pCorrFit.getMass1();
            A[state] = pCorrFit.getA();

            if(toScalar(mean(mass_1[state] - mass_0[state])) <= 0.0)
              {
                pCorrLog << "Zeroing m1 and A b/c m1 < m0" << std::endl;
                mass_0[state] = mass_0[state] / mass_0[state];
                mass_1[state] = toScalar(0.0) * mass_1[state];
                A[state] = toScalar(0.0) * A[state];
              }
          }
        else
          {
            mass_1[state] = zero;
            A[state] = zero;
          }

        pCorrFitName[state] = pCorrFit.getFitName();
        pCorrFitPlot[state] = pCorrFit.getFitPlotString();
        pCorrChiSqPDoF[state] = pCorrFit.getChisq() / pCorrFit.getNDoF();
	pfitlog << pCorrFit.getFitSummary() << std::endl;
	pfitlog << "** ensem fit value ** m= " << toScalar(mean(mass_0[state])) 
		<< " +/- " << toScalar(variance(mass_0[state])) << std::endl;
	pCorrFitLog.push_back(pfitlog.str());

      }//next pcorr




    pcorr = true;
  }
  
//*** write out the spectrum in a nice format -  whatever order the genEig solver chose  ***
  template<class T>
  void ST0Fit<T>::printSpectrum(ostream &output) const
  {
    if(!!!pcorr)
      {
        output << "Need to first fit principal correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    output << endl << endl;
    output << "**********************" << endl;
    output << "* UNORDERED SPECTRUM *" << endl;
    output << "**********************" << endl;

    output << "state|        mass         |           fit           |chisq/nDoF |  t0 = " << t0 << "|" << endl;

    //** MAKE THIS PRETTIER WITH prinf **

    for(int i = 0; i < M; ++i)
      {
        output << setw(5) << i << "|" << setw(8) << fixed << setprecision(5) << toDouble(mean(mass_0[i]))
               << " +/- " <<  setw(8) << fixed << setprecision(5) << toDouble(sqrt(variance(mass_0[i]))) << "|";
        output << setw(25) << pCorrFitName[i] << "|" << setw(10) << fixed << setprecision(2) << pCorrChiSqPDoF[i] << " |";

        if(toDouble(mean(mass_1[i])) > 0.0)
          {
            output << "  m'=" << setw(6) << fixed << setprecision(3) << toDouble(mean(mass_1[i]))
                   << " +/- " <<  setw(6) << fixed << setprecision(3) << toDouble(sqrt(variance(mass_1[i])));
            output << ", A="  << setw(6) << fixed << setprecision(3) << toDouble(mean(A[i]))
                   << " +/- " <<  setw(6) << fixed << setprecision(3) << toDouble(sqrt(variance(A[i])));
          }

        output << endl;
      }
  }

//prinCorrPlots
  template<class T>
  void ST0Fit<T>::printPrinCorrFits(void) const
  {

    if(!!!pcorr)
      {
        std::cout << "Need to first fit principal correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/PrinCorrPlots");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int state = 0; state < M; ++state)
      {
        std::stringstream s;
        s <<  path << "prin_corr_fit_t0" << t0 << "_state" << state << ".ax" ;
        std::ofstream out;
        out.open(s.str().c_str());
        out << pCorrFitPlot[state];
        out.close();
      }

    std::stringstream s;
    s << path << "pcorr_log";
    std::ofstream out;
    out.open(s.str().c_str());
    out << pCorrLog.str();
    out.close();
  }

  template<class T>
  void ST0Fit<T>::printPrinCorrReorder(const std::map<int, int> &mapp) const
  {

    if(!!!pcorr)
      {
        std::cout << "Need to first fit principal correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/PrinCorrPlots");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int, int>::const_iterator it;

    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	std::stringstream s;
	s <<  path << "prin_corr_fit_t0" << t0 << "_reorder_state" << it->first << ".ax" ;
	std::ofstream out;
	out.open(s.str().c_str());
	out << pCorrFitPlot[it->second];
	out.close();
      }
  }

  template<class T>
  void ST0Fit<T>::printMassFile(void) const
  {
    if(!!!pcorr)
      {
        std::cout << "Need to first fit principal correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/MassJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int state = 0; state < M; ++state)
      {
        std::ostringstream file;
        file << path <<  "mass_t0_" << t0 << "_state" << state << ".jack";
        write(file.str(), mass_0[state]);
      }
  }

  template<class T>
  void ST0Fit<T>::printMassFileReorder(const std::map<int, int> &mapp) const
  {

    if(!!!pcorr)
      {
        std::cout << "Need to first fit principal correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/MassJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int, int>::const_iterator it;

    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	std::ostringstream file;
	file << path <<  "mass_t0_" << t0 << "_reorder_state" << it->first << ".jack";
	write(file.str(), mass_0[it->second]);
      }
  }

  template<class T>
  void ST0Fit<T>::printPrinCorrFiles(void) const
  {
    if(!!!init)
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " not initialized, exiting.." << std::endl;
        exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/PrinCorrFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");


    for(int state = 0; state < M; ++state)
      {
        EnsemVectorReal lambdat;
        lambdat.resize(B);
        lambdat.resizeObs(inikeys.prinCorrProps.tmax + 1);

        for(int t = 0; t < inikeys.globalProps.tmin; ++t)
          pokeObs(lambdat, eVals[0].getEnsemElement(state) / eVals[0].getEnsemElement(state), t);

        for(int t = inikeys.globalProps.tmin; t < inikeys.prinCorrProps.tmax; ++t)
          pokeObs(lambdat, eVals[t].getEnsemElement(state), t);

        std::ostringstream file;
        file << path << "prin_corr_t0" << t0 << "_state" << state << ".jack";
        write(file.str(), lambdat);
      }

  }

  template<class T>
  void ST0Fit<T>::printPrinCorrFilesReorder(const std::map<int, int> &mapp) const
  {
    if(!!!init)
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " not initialized, exiting.." << std::endl;
        exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/PrinCorrFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");


    std::map<int, int>::const_iterator it;

    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	EnsemVectorReal lambdat;
	lambdat.resize(B);
	lambdat.resizeObs(inikeys.prinCorrProps.tmax + 1);
	
	EnsemReal one;
	one.resize(B);
	one = toScalar(1.);
	
	for(int t = 0; t < inikeys.globalProps.tmin; ++t)
	  pokeObs(lambdat, one, t);
	
	for(int t = inikeys.globalProps.tmin; t <= inikeys.prinCorrProps.tmax; ++t)
	  pokeObs(lambdat, eVals[t].getEnsemElement(it->second), t);
	
	std::ostringstream file;
	file << path << "prin_corr_t0" << t0 << "_state" << it->first << "_reordered.jack";
	write(file.str(), lambdat);
	
      }
  }

// print pcorr fit logs unordered
  template<class T>
  void ST0Fit<T>::printPrinCorrFitLog(void) const
  {
    if(!!!pcorr)
      {
	std::cout << "Need to first fit principal correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
	exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/PrinCorrFitLogs");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int state = 0; state < M; ++state)
      {
	std::stringstream s;
	s << path << "prin_corr_fit_log_t0" << t0 << "_state" << state << ".log";
	std::ofstream out;
	out.open(s.str().c_str());
	out << pCorrFitLog[state];
	out.close();
      }
  }

// print pcorr fit logs reordered 
  template<class T>
  void ST0Fit<T>::printPrinCorrFitLogReorder(const std::map<int, int> &mapp) const
  {
    
    if(!!!pcorr)
      {
	std::cout << "Need to first fit principal correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
	exit(1);
      }

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/PrinCorrFitLogs");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int,int>::const_iterator it;
    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	if(it->second >= M)
	  {
	    std::cout << "remapping error, out of bounds" << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	std::stringstream s;
	s << path << "prin_corr_fit_log_t0" << t0 << "_reorder_state" << it->first << ".log";
	std::ofstream out;
	out.open(s.str().c_str());
	out << pCorrFitLog[it->second];
	out.close();
      }
  }

//make the overlaps
  template<class T>
  void ST0Fit<T>::makeZ(void)
  {
    Z.clear();
    initChk();
    pcorrChk();

    if(inikeys.globalProps.verbose)
      std::cout << "Making Z.." << std::endl;


    //make it for everything
    //    for(int t = 0; t <= inikeys.zProps.tmax; ++t)
    for(int t = 0; t <= std::max(inikeys.zProps.tmax, inikeys.reconProps.tmax); ++t)
      //jjd change to make zProps belong to fitting
      {
        SembleVector<T> dum_vals(B, M);
        dum_vals.zeros();

        for(int index = 0; index < M; ++ index)
          dum_vals.loadEnsemElement(index, (sqrt(Real(2.0) * mass_0[index]) * exp(mass_0[index] * Real(t0) / Real(2.0))));

        Z.push_back(diag(dum_vals)*adj(eVecs[t])*Ct0);
      }

    zinit = true;
  }

//fit the overlaps
  template<class T>
  void ST0Fit<T>::fitZ(void)
  {
    std::cout << "Fitting Z.." << std::endl;
    initChk();
    pcorrChk();

    if(!!!zinit)
      makeZ();

    zFit.reDim(B, M, N);
    zFit.zeros();

    zFitPlot.clear();
    zFitResults.clear();
    zFitChiSqPDof.clear();

    zFitPlot.resize(M, std::vector<std::string>(N, ""));
    zFitResults = zFitPlot;
    zFitChiSqPDof.resize(M, std::vector<double>(N, 0.));

    for(int state = 0; state < M; ++state)
      {

        if(inikeys.globalProps.verbose){
	  std::cout << "**********************************" << std::endl;
          std::cout << "t0 = " << t0 << " - fitting Z(state=" << state << ")" << std::endl;
	  std::cout << "**********************************" << std::endl;
	}

        for(int op = 0; op < N; ++op)
          {
	    if(inikeys.globalProps.verbose)
	      std::cout << "**** t0 = " << t0 << "   -> fitting state = " << state << ", op = " << op << std::endl;

            EnsemVectorReal z;
            std::vector<double> tslice;
            z.resize(B);
            z.resizeObs(inikeys.zProps.tmax + 1 - inikeys.globalProps.tmin);

	    /////////////////////////////////////////////////////////////////////////////
	    //std::cout << "DEBUG : zProps.tmax = " <<  inikeys.zProps.tmax << std::endl;
	    ////////////////////////////////////////////////////////////////////////////

            for(int t = inikeys.globalProps.tmin; t <= inikeys.zProps.tmax; ++ t)
              {
                pokeObs(z, Z[t].getEnsemElement(state, op), t - inikeys.globalProps.tmin);
                tslice.push_back(t);
              }

            EnsemData zData(tslice, z);
	    Handle<FitComparator> fitC = fitComp(inikeys.zProps.fitCrit);
	    //Handle<FitComparator> fitC = new CompareZFits;
            FitZ fitZ(zData, t0, fitC, inikeys.zProps.minTSlices);

	    /*	    if(inikeys.globalProps.verbose){
	      string tmp = fitZ.getFitSummary();
	      tmp.resize(std::min(int(tmp.size()), 1000));
	      std::cout << tmp << "    ..." << std::endl << std::endl ;
	      }*/

            zFit.loadEnsemElement(state, op, fitZ.getZ());
            zFitPlot[state][op] = fitZ.getFitPlotString();
            zFitResults[state][op] = fitZ.getFitSummary();
            zFitChiSqPDof[state][op] = fitZ.getChisq() / fitZ.getNDoF();
          }
      }

    zfit = true;
  }

//print the fits
  template<class T>
  void ST0Fit<T>::printZFits(void) const
  {

    zfitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZFitPlots");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int index = 0; index < M; ++index)
      {
        for(int op = 0; op < N; ++op)
          {
            std::stringstream s;
            s << path << "Z_fit_t0" << t0 << "_state" << index << "_op" << op << ".ax";
            std::ofstream out;
            out.open(s.str().c_str());
            out << zFitPlot[index][op];
            out.close();
          }
      }
  }

  template<class T>
  void ST0Fit<T>::printZFitsReorder(const std::map<int, int> &mapp) const
  {

    zfitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZFitPlots");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int, int>::const_iterator it;

    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	for(int op = 0; op < N; ++op)
	  {
	    std::stringstream s;
	    s << path << "Z_fit_t0" << t0 << "_reorder_state" << it->first << "_op" << op << ".ax";
	    std::ofstream out;
	    out.open(s.str().c_str());
	    out << zFitPlot[it->second][op];
	    out.close();
	  }
      }
  }

  template<class T>
  void ST0Fit<T>::printZFile(void) const
  {
    zfitChk();
    zfitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int state = 0; state < M; ++state)
      for(int op = 0; op < N; ++op)
        {
          std::ostringstream file;
          file << path << "Z_t0_" << t0 << "_state" << state << "_op" << op << ".jack";
          write(file.str(), zFit.getEnsemElement(state, op));
        }
  }

  template<class T>
  void ST0Fit<T>::printZFileReorder(const std::map<int, int> &mapp) const
  {
    zfitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int, int>::const_iterator it;

    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	for(int op = 0; op < N; ++op)
	  {
	    std::ostringstream file;
	    file << path << "Z_t0_" << t0 << "_reorder_state" << it->first << "_op" << op << ".jack";
	    write(file.str(), zFit.getEnsemElement(it->second, op));
	  }
      }
  }

  template<class T>
  void ST0Fit<T>::printZt(void) const
  {
    zinitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/Z_tJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    typename PromoteEnsem<T>::Type one;
    one.resize(B);
    one = toScalar(T(1.0));


    for(int state = 0; state < M; ++state)
      for(int op = 0; op < N; ++op)
        {
          typename PromoteEnsemVec<T>::Type zt;
          zt.resize(B);
          zt.resizeObs(inikeys.zProps.tmax);

          for(int t = 0; t < inikeys.globalProps.tmin; ++t)
            pokeObs(zt, one, t);

          for(int t = inikeys.globalProps.tmin; t < inikeys.zProps.tmax; ++t)
            pokeObs(zt, Z[t].getEnsemElement(state, op), t);

          std::ostringstream file;
          file << path << "Zt_t0_" << t0 << "_state" << state << "_op" << op << ".jack";
          write(file.str(), zt);
        }
  }

  template<class T>
  void ST0Fit<T>::printZtReorder(const std::map<int, int> &mapp) const
  {
    zinitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/Z_tJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int, int>::const_iterator it;

    typename PromoteEnsem<T>::Type one;
    one.resize(B);
    one = toScalar(T(1.0));

    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	for(int op = 0; op < N; ++op)
	  {
	    
	    typename PromoteEnsemVec<T>::Type zt;
	    zt.resize(B);
	    zt.resizeObs(inikeys.zProps.tmax);

	    for(int t = 0; t < inikeys.globalProps.tmin; ++t)
	      pokeObs(zt, one, t);
	    
	    for(int t = inikeys.globalProps.tmin; t < inikeys.zProps.tmax; ++t)
	      pokeObs(zt, Z[t].getEnsemElement(it->second, op), t);
	    
	    
	    std::ostringstream file;
	    file << path << "Zt_t0_" << t0 << "_reorder_state" << it->first << "_op" << op << ".jack";
                write(file.str(), zt);
	  }
      }    
  }


//print the fits
  template<class T>
  void ST0Fit<T>::printZFitLog(void) const
  {

    zfitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZFitLogs");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int index = 0; index < M; ++index)
      {
        for(int op = 0; op < N; ++op)
          {
            std::stringstream s;
            s << path << "Z_fit_log_t0" << t0 << "_state" << index << "_op" << op << ".log";
            std::ofstream out;
            out.open(s.str().c_str());
            out << zFitResults[index][op];
            out.close();
          }
      }
  }

  template<class T>
  void ST0Fit<T>::printZFitLogReorder(const std::map<int, int> &mapp) const
  {

    zfitChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZFitLogs");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int, int>::const_iterator it;

    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	for(int op = 0; op < N; ++op)
	  {
	    std::stringstream s;
	    s << path << "Z_fit_log_t0" << t0 << "_reorder_state" << it->first << "_op" << op << ".log";
	    std::ofstream out;
	    out.open(s.str().c_str());
	    out << zFitResults[it->second][op];
	    out.close();
	  }
      }
  }


  template<class T>
  void ST0Fit<T>::printVt(void) const
  {
    initChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/V_tJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int state = 0; state < M; ++state)
      for(int op = 0; op < N; ++op)
        {
          typename PromoteEnsemVec<T>::Type evec;
          evec.resize(B);
          evec.resizeObs(inikeys.globalProps.tmax + 1);
	      
	      
          for(int t = 0; t <= inikeys.globalProps.tmax; ++t)
            pokeObs(evec, eVecs[t].getEnsemElement(op, state), t);

          std::ostringstream file;
          file << path << "V_t0_" << t0 << "_state" << state << "_op" << op << ".jack";
          write(file.str(), evec);
        }

  }

  template<class T>
  void ST0Fit<T>::printVtReorder(const std::map<int, int> &mapp) const
  {
    initChk();

    std::stringstream ss;
    ss << "t0" << t0;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/V_tJackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    std::map<int, int>::const_iterator it;


    for(it = mapp.begin(); it != mapp.end(); ++it)
      {
	if(it->second >= M) //in bounds check
	  {
	    std::cout << "remapping error, out of bounds " << __PRETTY_FUNCTION__ << std::endl;
	    exit(1);
	  }
	
	for(int op = 0; op < N; ++op)
	  {
	    typename PromoteEnsemVec<T>::Type evec;
	    evec.resize(B);
	    evec.resizeObs(inikeys.globalProps.tmax + 1);
	    
	    for(int t = 0; t <= inikeys.globalProps.tmax; ++t)
	      pokeObs(evec, eVecs[t].getEnsemElement(op, it->second), t);
	    
	    std::ostringstream file;
	    file << path << "V_t0_" << t0 <<  "_reordered_state" << it->first << "_op" << op << ".jack";
	    write(file.str(), evec);
	    
	  }
      }    
  }



//reconstruction
  template<class T>
  void ST0Fit<T>::reconCorr(void)
  {
    initChk();
    pcorrChk();

    if(!!!zinit)
      makeZ();

    std::cout << "Reconstructing Correlators.. " << std::endl;
    itpp::Real_Timer rt;
    rt.tic();

    recon_total.clear();

    SembleMatrix<T> Zt = Z[t0 + 1];
    typename std::vector<SembleMatrix<T> > recon_vec, M;
    SembleMatrix<T> massM(B, mass_0.size(), mass_0.size());
    massM.zeros();

    for(int t = t0 + 1; t <= inikeys.reconProps.tmax; ++t)
      {
        for(int state = 0; state < mass_0.size(); ++state)
          {
            massM.loadEnsemElement(state, state, exp(-mass_0[state]*Real(t)) / (Real(2.0)*mass_0[state]));
          }

        M.push_back(massM);
        recon_vec.push_back(adj(Zt)*massM * Zt);
      }

    recon_total[t0 + 1] = recon_vec;

    for(int tz = t0 + 2; tz <= inikeys.reconProps.tmax; ++tz)
      {
        Zt = Z[tz];
        recon_vec.clear();

        typename std::vector<SembleMatrix<T> >::const_iterator it;  //M runs from recon t0 + 1 to tmax

        for(it = M.begin(); it != M.end(); ++it)
          recon_vec.push_back(adj(Zt) * (*it)*Zt);

        recon_total[tz] = recon_vec;                                //runs from t0 +1 to recon_props.tmax
      }

    recon = true;

    std::cout << "Correlator Reconstruction took " << rt.toc() << " seconds" << std::endl;
  }

  template<class T>
  double ST0Fit<T>::reconChisq(const int tz, const ReconProps_t &rec)
  {
    initChk();
    pcorrChk();
    zinitChk();

    if(!!!recon)
      reconCorr();


    if(!!!recon_flat)
      {
	recon_flat = true;
	typename std::vector<SembleMatrix<T> > d;
	typename std::vector<SembleMatrix<T> >::const_iterator it;
	for(it = tp.begin() +t0 +1; it != tp.begin() + rec.tmax +1; it++)
	  {
	    if(it == tp.end())
	      {
		std::cout << "BARFING" << std::endl;
		exit(1);
	      }
	    d.push_back(*it);
	  }

	flatData = flattenSym(d);   //tp goes [0,inikeys.globalProps.tmax]
      }

    flatData.setSVCutoff(inikeys.globalProps.SVCut);

    typename std::vector<SembleMatrix<T> > dd;
    typename std::vector<SembleMatrix<T> >::const_iterator dit;
    for(dit = recon_total[tz].begin(); dit != recon_total[tz].end(); dit++)
      {
	if(dit == tp.end())
	  {
	    std::cout << "BARFING" << std::endl;
	    exit(1);
	  }
	dd.push_back(*dit);
      }

    EnsemData flatRecon = flattenSym(dd); //recon_total goes [t0+1,rec.tmax]

    return ensemDataChisq(flatData, flatRecon) / (flatData.getNData() - flatData.getNResetCovSingVals() - (N + N * (N + 1) / 2));
  }

  template<class T>
  double ST0Fit<T>::reconChisqFast(const int tz, const ReconProps_t &rec)
  {
    initChk();
    pcorrChk();
    zinitChk();

    if(!!!recon)
      reconCorr();

    int nDoF = -(N + N * (N + 1) / 2); // a 'guess' at the number of variable pars you'd have if you did a fit
    double chisq = 0.0;

    for(int i = 0; i < N; ++i)
      {
        for(int j = i; j < N; ++j)
          {
            std::vector<double> tslice;
            EnsemVectorReal data, recon;
            data.resize(B);
            data.resizeObs(rec.tmax - t0);
            recon = data;

            for(int t = t0 + 1; t <= rec.tmax; ++t)
              {
                tslice.push_back(t);
                pokeObs(data, (tp[t].getEnsemElement(i, j) + tp[t].getEnsemElement(j, i))*Real(0.5), t - t0 - 1); //symmetrize tp                  //there should be a cc here, want C_ij + C_ji* for complex types 
                pokeObs(recon, (recon_total[tz])[t - t0 - 1].getEnsemElement(i, j), t - t0 - 1); //this was 'made' so its symmetric to round off
              }

            EnsemData edata(tslice, data), erecon(tslice, recon);
            edata.setSVCutoff(inikeys.globalProps.SVCut);

            chisq += ensemDataChisq(edata, erecon) * ((i == j) ? 1.0 : 2.0);
            nDoF += (rec.tmax - t0 - edata.getNResetCovSingVals()) * ((i == j) ? 1 : 2);
          }
      }

    return chisq / double(nDoF);
  }

//find the best set of overlaps by reconstruction
  template<class T>
  pair<int, double> ST0Fit<T>::findBestTZ(const ReconProps_t &rec)
  {
    initChk();
    pcorrChk();

    if(!!!zinit)
      makeZ();

    if(!!!recon)
      reconCorr();

    double(ST0Fit<T>::*ptr)(const int, const ReconProps_t &) = NULL;   //member-function pointer
    int best = t0 + 1;
    double prev, next;
    std::stringstream ss; 
  

    ptr = (rec.type == "fast") ? &ST0Fit<T>::reconChisqFast : &ST0Fit<T>::reconChisq; //maybe this should be a private member function?
    prev = (*this.*ptr)(best, rec);

    if(inikeys.globalProps.verbose)
      std::cout << "tz = " << best << " with chisq/ndof = " << prev << std::endl;

    ss << "tz = " << best << " chisq/dof = " << prev << "\n";

    for(int tz = t0 + 2; tz <= rec.tmax; ++tz)                         //find best by minimizing the chisq
    {
      next = (*this.*ptr)(tz, rec);

      ss << "tz = " << tz << " chisq/dof = " << next << "\n";

      if(inikeys.globalProps.verbose)
        std::cout << "tz = " << tz << " with chisq/ndof = " << next << std::endl;

      if(next < prev)
      {
        best = tz;
        prev = next;
      }
    }

    ptr = NULL;

    if(inikeys.globalProps.verbose)
      std::cout << "tz best = " << best << " with chisq/ndof = " << prev << std::endl;

    tz_best = best;
    recon_tz = true;


    // we just did a whole bunch of work.. might as well keep track of it
    std::stringstream name;
    name << "t0" << t0 ;
    std::string out = SEMBLE::SEMBLEIO::getPath() + name.str();
    SEMBLE::SEMBLEIO::makeDirectoryPath(out);
    out += std::string("/tz_chisq.log");
    std::ofstream fout(out.c_str());
    fout << ss.str() << std::endl;
    fout.close(); 

    return pair<int, double>(best, prev);
  }

  template<class T>
    void ST0Fit<T>::makeReconPlots(std::string &mode)
    {
      //do a recon if we haven't
      if(!!!recon_tz)
        std::pair<int, double> dum = findBestTZ(inikeys.reconProps);

      recon_plots.clear();

      SembleMatrix<T> Zt = Z[tz_best];
      SembleMatrix<T> Zero = 0.0 * Ct0;
      typename std::vector<std::vector<itpp::Mat<T> > > summed_by_state;
      typename std::vector<SembleMatrix<T> > total;
      total.resize(inikeys.reconProps.tmax + 1, Zero);

      for(int state = 0; state < M; ++state)
      {
        SembleVector<T> z = getRow(Zt, state);
        SembleMatrix<T> zz = outterProduct(z, z);
        std::vector<itpp::Mat<T> > corr_t_mean;

        for(int t = 0; t <= inikeys.reconProps.tmax; ++t)
        {
          EnsemReal e = exp(- mass_0[state] * Real(t)) / (Real(2.0) * mass_0[state]);
          total[t] += zz * e;
          corr_t_mean.push_back(mean(total[t]));
        }

        summed_by_state.push_back(corr_t_mean);
      }

      //now make the plots for tz_best
      //suppose the states have been reordered, we want to multiply by the lowest mass..
      int lightest = 0;
      double light = toScalar(mean(mass_0[lightest]));

      for(int i = 1; i < M; ++i)
      {
        double m = toScalar(mean(mass_0[i]));

        if(m < light)
        {
          lightest = i;
          light = m;
        }
      }

      ConstTimesExp weight(exp(- mass_0[lightest] * Real(t0)), mass_0[lightest]);
      std::vector<Real> dweight;
      for(int t = 0; t <= inikeys.reconProps.tmax; ++t)
        dweight.push_back(toScalar(exp(light*double(t-t0))));


      for(int i = 0; i < N; ++i)
      {
        int bound = (mode == "diag") ? i + 1 : N;

        for(int j = i; j < bound; ++j)
        {
          AxisPlot plot;
          std::vector<std::vector<double> > sum_ij; //outter is state, inner time
          std::vector<double> time(inikeys.reconProps.tmax + 1, 0.0);
          sum_ij.resize(M, time);


          for(int state = 0; state < M; ++state)
            for(int t = inikeys.globalProps.tmin; t <= inikeys.reconProps.tmax; ++t)
              sum_ij[state][t] = summed_by_state[state][t](i, j);

          std::vector<double> dat, tslices;

          for(int t = inikeys.globalProps.tmin; t <= inikeys.reconProps.tmax; ++t)
          {
            dat.push_back(toScalar(mean(weight(t)))*sum_ij[0][t]);
            tslices.push_back(t);
          }

          plot.addLineData(tslices, dat, 1);
          dat.clear();

          //plot data line if it is significant
          for(int state = 1; state < M - 1; ++state)
          {
            if(fabs(sum_ij[state][t0] - sum_ij[state - 1][t0]) < sum_ij[M - 1][t0] / M)
              continue;

            for(int t = inikeys.globalProps.tmin; t <= inikeys.reconProps.tmax; ++t)
              dat.push_back(toScalar(mean(weight(t)))*sum_ij[state][t]);

            plot.addLineData(tslices, dat, (state) % 5 + 1);
            dat.clear();
          }

          //add last state and total

          //add the final state & plot the total
          std::vector<double> m, mpe, mme;
          tslices.clear();

          for(int t = t0; t <= inikeys.reconProps.tmax; t++)
          {
            EnsemReal d = weight(t) * total[t](i, j) ;
            //	cout << "read an element of total at t=" << t << endl;
            double me = toDouble(mean(d));
            double err = toDouble(sqrt(variance(d)));
            m.push_back(me);
            mpe.push_back(me + err);
            mme.push_back(me - err);
            tslices.push_back(t);
          }

          plot.addLineData(tslices, m, 1);
          plot.addLineData(tslices, mpe, 1);
          plot.addLineData(tslices, mme, 1);

          m.clear();
          mpe.clear();
          mme.clear();
          tslices.clear();

          for(int t = inikeys.globalProps.tmin; t <= t0; t++)
          {
            EnsemReal d = weight(t) * total[t](i, j) ;
            double me = toDouble(mean(d));
            double err = toDouble(sqrt(variance(d)));
            m.push_back(me);
            mpe.push_back(me + err);
            mme.push_back(me - err);
            tslices.push_back(t);
          }

          plot.addLineData(tslices, m, 3);
          plot.addLineData(tslices, mpe, 3);
          plot.addLineData(tslices, mme, 3);

          //add the ensem data
          EnsemVectorReal raw;
          raw.resize(B);
          raw.resizeObs(inikeys.reconProps.tmax - inikeys.globalProps.tmin + 1);

          mme.clear();
          mpe.clear();
          tslices.clear();

          for(int t = inikeys.globalProps.tmin; t <= inikeys.reconProps.tmax; t++)
          {
            pokeObs(raw, dweight[t]*tp[t](i, j), t - inikeys.globalProps.tmin);
            double me = toDouble(mean(tp[t](i, j)));
            double err = toDouble(sqrt(variance(tp[t](i, j))));

            if(t >= t0)
            {
              mpe.push_back(me + err);
              mme.push_back(me - err);
            }

            tslices.push_back(t);
          }

          plot.addEnsemData(tslices, raw, "\\sq", 1);

          //set ranges & add labels
          double low = *min_element(mme.begin(), mme.end());
          double high = *max_element(mpe.begin(), mpe.end()) * 1.5 ;
          plot.setYRange(low, high);
          plot.setXRange(-0.5, inikeys.reconProps.tmax + 0.5);

          double height = high - low;
          double unit = height / (M + 1);

          for(int state = 0; state < M - 1; state++)
          {
            int colour = (state) % 5 + 1;
            char c[10];
            int n = sprintf(c, "%.3f", toDouble(mean(mass_0[state])));
            string str(c);
            plot.addLabel(inikeys.reconProps.tmax - 1, low + (state + 1)*unit, str, colour, 0.8);
          }

          char c[10];
          int n = sprintf(c, "%.3f", toDouble(mean(mass_0[M - 1 ])));
          std::string str(c);
          plot.addLabel(inikeys.reconProps.tmax - 1, low + M * unit, str, 1, 0.8);

          std::pair<int, int> key(std::pair<int, int>(i, j));
          std::pair<std::string, std::string> data;
          std::string plot_string = plot.getAxisPlotString();
          std::stringstream ss;
          ss <<  "recon_t0" << t0 << "_" << i << "_" << j << ".ax";
          data = std::pair<std::string, std::string>(ss.str(), plot_string);

          recon_plots[key] = data;

        }//end j
      }//end i
    }

  template<class T>
    void ST0Fit<T>::printReconPlots(void)
    {
      if(recon_plots.begin() == recon_plots.end())
        makeReconPlots(inikeys.outputProps.reconType);

      std::stringstream ss;
      ss << "t0" << t0;
      std::string path = SEMBLEIO::getPath() += ss.str();
      SEMBLEIO::makeDirectoryPath(path);
      path += std::string("/ReconPlots");
      SEMBLEIO::makeDirectoryPath(path);
      path += std::string("/");

      std::map<std::pair<int, int>, std::pair<std::string, std::string> >::const_iterator it;

      for(it = recon_plots.begin(); it != recon_plots.end(); ++it)
      {
        std::stringstream ss;
        std::ofstream out;
        ss << path << it->second.first.c_str();
        out.open(ss.str().c_str());
        out << it->second.second;
        out.close();
      }
    }


  template<class T>
    std::string ST0Fit<T>::getPCorrFitLog(const std::map<int,int> &mapp)
    {
      std::stringstream ss;
      std::string n("\n");
      std::map<int,int>::const_iterator it;

      ss << n << n;
      ss << "**********************" << n;
      ss << "* REORDERED SPECTRUM *" << n;
      ss << "**********************" << n;

      ss << "state|unord|        mass       |           fit           |chisq/nDoF|   fit parameters " << endl;


      for(it = mapp.begin(); it != mapp.end(); ++it)
      {
        ss << setw(5) << it->first << "|" << setw(5) << it->second << "|"                               // states
          << setw(8) << fixed << setprecision(5) << toScalar(mean(mass_0[it->second])) << "+/-"        // m_0
          << setw(8) << fixed << setprecision(5) << std::sqrt(toScalar(variance(mass_0[it->second])))  // var(m_0)
          << "|" << setw(25) << pCorrFitName[it->second] << "|"                                        // fit name
          << setw(10) << fixed << setprecision(2) << pCorrChiSqPDoF[it->second] << "|";                // fit chisq

        if(toScalar(mean(mass_0[it->second])) > 0.)
          ss << " m'=" << setw(6) << fixed << setprecision(3) << toScalar(mean(mass_1[it->second]))   // second exp mass
            << " +/- " << setw(6) << fixed << setprecision(3) << std::sqrt(toScalar(variance(mass_1[it->second])))
            << ", A=" << setw(6) << fixed << setprecision(3) << toScalar(mean(A[it->second]))        // coeff
            << " +/- " << setw(6) << fixed << setprecision(3) << std::sqrt(toScalar(variance(A[it->second]))) << n;
      }

      return ss.str();
    }



  template<class T>
    void ST0Fit<T>::fixedCoeffMethod(const int tstar) const
    {
      if(!!!init)
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " cant solve without init, exiting" << std::endl;
        exit(1);
      }

      if(inikeys.globalProps.verbose)
        std::cout << __func__ << ": Trying the fixed coeff method at t0 = " << t0 << ": solving at t* = " << tstar << std::endl;

      SembleMatrix<T> refEvecs = eVecs[tstar];
      typename std::vector<SembleMatrix<T> > rCorr;
      typename std::vector<SembleVector<T> > rPCorr;

      for(int t = inikeys.globalProps.tmin; t <= inikeys.prinCorrProps.tmax; ++t)
      {
        rCorr.push_back(adj(refEvecs)*tp[t]*refEvecs);
        SembleVector<T> tmp(B, M);

        for(int elem = 0; elem < M; ++elem)
          tmp.loadEnsemElement(elem, (rCorr.back())->getEnsemElement(elem, elem));

        rPCorr.push_back(tmp);
      }

      //now write them out to a file
      if(inikeys.globalProps.verbose)
        std::cout << __func__ << ": Writing fixedcoeff_prin_corr files" << std::endl;

      //make a directory
      std::stringstream ss;
      ss << "t0" << t0;
      std::string path = SEMBLEIO::getPath() += ss.str();
      SEMBLEIO::makeDirectoryPath(path);
      path += std::string("/FixedCoeffPrinCorr");
      SEMBLEIO::makeDirectoryPath(path);
      path += std::string("/");

      //write out each state
      for(int state = 0; state < M; ++state)
      {
        typename PromoteEnsemVec<T>::Type pt;
        pt.resize(B);
        pt.resizeObs(inikeys.prinCorrProps.tmax + 1);
        typename PromoteEnsem<T>::Type one;
        one.resize(B);
        one = toScalar(T(1.0));

        for(int t = 0; t < inikeys.globalProps.tmin; ++t)
          pokeObs(pt, one, t);

        for(int t = inikeys.globalProps.tmin; t <= inikeys.prinCorrProps.tmax; ++t)
          pokeObs(pt, rPCorr[t].getEnsemElement(state), t);

        std::ostringstream file;
        file << path << "fixedCoeffPrinCorr_t0" << t0 << "_ts" << tstar << "_state" << state << ".jack";
        write(file.str(), pt);
      }

    }


  template<class T> //this is expensive!!! it should be called before everything else, it assumes the states are already ordered within this t0 fit
    void ST0Fit<T>::rephase(const typename std::map<int, typename PromoteScalar<T>::Type> &mapped)
    {
      initChk();
      //sanity
      std::vector<bool> used(M, false);
      typename std::map<int, typename PromoteScalar<T>::Type>::const_iterator it;
      bool check = true;

      for(it = mapped.begin(); it != mapped.end(); ++it)
      {
        if(it->first >= M)
        {
          check = false;
          std::cout << "t0 = " << t0 << " it->first = " << it->first << " " << M << " states possible in" << std::endl;
          std::cout << __PRETTY_FUNCTION__ << std::endl;
          break;
        }

        used[it->first] = true;
      }

      if(check)
      {
        std::vector<bool>::const_iterator bit;

        for(bit = used.begin(); bit != used.end(); ++bit)
          if(!!!*bit)
          {
            std::cout << "unphased state in" << std::endl;
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            check = false;
            break;
          }
      }

      if(!!!check)
      {
        std::cout << "Rephasing error in " << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

      for(int t = inikeys.globalProps.tmin; t < eVecs.size(); ++t)
      {
        SembleMatrix<T> cp = eVecs[t];

        for(it = mapped.begin(); it != mapped.end(); ++it)
        {
          for(int row = 0; row < N; ++row)
            eVecs[t].loadEnsemElement(row, it->first, it->second*cp.getEnsemElement(row, it->first));
        }
      }

      //make sure any Z's we already made have the correct phase
      if(zinit)
        makeZ();

      if(zfit)
        fitZ();
    }

  //Internal Functions
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class T>
    Handle<FitComparator> ST0Fit<T>::fitComp(std::string &in) const
    {
      if(in == "chisq_per_dof")
        return Handle<FitComparator>(new CompareFitsByChisqPerNDoF);

      if(in == "Q")
        return Handle<FitComparator>(new CompareFitsByQ);

      if(in == "splitN")
        return Handle<FitComparator>(new CompareFitsBySplitN);

      if(in == "generic")
        return Handle<FitComparator>(new CompareFitsByGeneric);

      if(in == "QN")
        return Handle<FitComparator>(new CompareFitsByQN);

      if(in == "Zfit")
        return Handle<FitComparator>(new CompareZFits);

      std::cout << __PRETTY_FUNCTION__ << "Fit Criterion " << in << " unknown, defaulting to chisq_per_dof" << std::endl;
      return Handle<FitComparator>(new CompareFitsByChisqPerNDoF);
    }

  template<class T>
    void ST0Fit<T>::cleanUp(void)
    {
      N = -1;
      M = -1;
      B = -1;
      t0 = -1;
      inikeys = FitIniProps_t();
      Ct0 = SembleMatrix<T>();
      primFit.clear();
      primFit = ST0FitPrim<T>();
      eVals.clear();
      eVecs.clear();
      Z.clear();
      mass_0.clear();
      mass_1.clear();
      A.clear();
      init = pcorr = zfit = zinit = recon = recon_flat = recon_tz = false;
      pCorrFitName.clear();
      pCorrFitPlot.clear();
      pCorrChiSqPDoF.clear();
      zFitPlot.clear();
      zFitResults.clear();
      zFitChiSqPDof.clear();
      zFit = SembleMatrix<T>();
      recon_total.clear();
      tp.clear();

    }

  //bool checks
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<class T>
    void ST0Fit<T>::initChk(void) const
    {
      if(!!!init)
      {
        std::cout << "Need to Initialize in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(100);
      }
    }

  template<class T>
    void ST0Fit<T>::pcorrChk(void) const
    {
      if(!!!pcorr)
      {
        std::cout << "Need to construct Principal Correlators in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(101);
      }
    }

  template<class T>
    void ST0Fit<T>::zfitChk(void) const
    {
      if(!!!zfit)
      {
        std::cout << "Need to fit Z in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(102);
      }
    }

  template<class T>
    void ST0Fit<T>::zinitChk(void) const
    {
      if(!!!zinit)
      {
        std::cout << "Need to construct Z in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(102);
      }
    }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
