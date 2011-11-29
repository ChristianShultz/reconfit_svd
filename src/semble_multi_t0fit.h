#ifndef SEMBLE_MULTI_T0FIT_H_H_GUARD
#define SEMBLE_MULTI_T0FIT_H_H_GUARD


#include"semble_load_correlators.h"
#include"semble_matrix.h"
#include"semble_vector.h"
#include"semble_file_management.h"
#include"semble_t0fit.h"
#include"ensem/ensem.h"
#include"adat/handle.h"
#include"semble_fit_ini_xml.h"
#include"fit_forms.h"
#include"fit_correlators.h"
#include"jackknife_fitter.h"
#include<map>
#include<vector>
#include<iostream>
#include<string>

namespace SEMBLE
{

  template<class T>
  struct SMT0Fit
  {

  public: //constructor, destructor, copy, assignment
    SMT0Fit(void);
    SMT0Fit(const SMT0Fit<T> &o);
    SMT0Fit(const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_);
    SMT0Fit(const typename std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_);
    ~SMT0Fit(void) {}
    SMT0Fit<T>& operator=(const SMT0Fit<T> &o);

  public:
    typename std::vector<SembleMatrix<T> > maketp(const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_);
    void load(const typename std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_);

  public:
    void rephaseAndReorder(void);

    void fitMassT0(void);                                //fit mass t0
    void printMass(void);                                //princorr plot
    void printMassT0(void);                              //mass t0 plot
    void printMassFiles(void);                           //jack files
    void printMassT0Files(void);                         //jack file of mt0 fit
    void printEvalFiles(void);                           //the principal correlator in jack format

    void fitZT0(void);                                   //fit Zt0
    void printZ(void);                                   //z plots
    void printZT0(void);                                 //zt0 plot
    void printZT0Files(void);                            //jack files of Z_ij(t0) fit
    void printZFiles(void);                              //jack files
    void printZtFiles(void);                             //z_ij(t) in jack form

    void printVtFiles(void);                             //eigenvectors

    void printReorderLog(void);                          //how the remaping worked
    void printNResetLog(void);                           //the number of reset singular values at each t0

    void printReconPlots(void);                          //print the recon plots

    void reorderStates(void);                            //order the states according to t0ref

  private:
    void clear(void);
    Handle<FitComparator> fitComp(const std::string &in) const;
    void printer(void);

  private:  //data store
    typename std::map<int, Handle<ST0Fit<T> > > t0_fits;   //key is t0 values, data is the ST0Fit(t0)
    std::map<int, std::pair<int, double> > tz_chisq;       //key is t0, data is tz_best and tz_chisq
    std::map<int, std::map<int, int> > reorder;            //key is t0, data is the reorder map, key of inner map is ref state, data of inner map is t0state
    FitIniProps_t inikeys;
    typename std::vector<SembleMatrix<T> > tp;             //the two point correlation matrix at each tslice
    int max_states;                                        //the maximum number of states considered

    //bool check
    bool init, mass_fit, mass_fit_t0, z_fit, z_fit_t0, reordered;   //bool checks

    //prin corr
    std::vector<EnsemReal> mass;                            //extracted mass from multi t0 fits
    std::vector<double> mass_chisq;                         //chisq on multi t0 mass fits
    std::vector<std::string> mass_summary, mass_plots;      //plot strings and summary for same
    std::stringstream fitlog;                               //a fit log

    //overlaps
    std::vector<std::vector<typename PromoteEnsem<T>::Type > > Z;   //extracted Z from the multi t0 fit
    std::vector<std::vector<double> > Z_chisq;                      //chisq on multi t0 Z fits
    std::vector<std::vector<std::string> > Z_summary, Z_plots;      //plots strings and summary
  };


  //constructors
  template<class T>
  SMT0Fit<T>::SMT0Fit(void)
    : init(false) , mass_fit(false) , mass_fit_t0(false) , z_fit(false) , z_fit_t0(false) , reordered(false) , max_states(-1)
  {}

  template<class T>
  SMT0Fit<T>::SMT0Fit(const SMT0Fit<T> &o)
    :  init(false) , mass_fit(false) , mass_fit_t0(false) , z_fit(false) , z_fit_t0(false) , reordered(false) , max_states(-1)
  {

    if(o.init)
      {
        t0_fits = o.t0_fits;
        inikeys = o.inikeys;
        tp = o.tp;
        tz_chisq = o.tz_chisq;
        max_states = o.max_states;
        init = true;
      }

    if(o.reordered)
      {
        reordered = true;
        reorder = o.reorder;
      }

    if(o.mass_fit)
      mass_fit = true;

    if(o.mass_fit_t0)
      {
        mass = o.mass;
        mass_chisq = o.mass_chisq;
        mass_summary = o.mass_summary;
        mass_plots = o.mass_plots;
        mass_fit_t0 = true;
      }

    if(o.z_fit)
      z_fit = true;

    if(o.z_fit_t0)
      {
        Z = o.Z;
        Z_chisq = o.Z_chisq;
        Z_summary = o.Z_summary;
        Z_plots = o.Z_plots;
        z_fit_t0 = true;
      }

  }

  template<class T>
  SMT0Fit<T>::SMT0Fit(const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_)
  {
    load(maketp(tp_, inikeys_), inikeys_);
  }

  template<class T>
  SMT0Fit<T>::SMT0Fit(const typename std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_)
  {
    load(tp_, inikeys_);
  }

  //assignment
  template<class T>
  SMT0Fit<T>& SMT0Fit<T>::operator=(const SMT0Fit<T> &o)
  {
    if(this != &o)
      {
        if(o.init)
          {
            t0_fits = o.t0_fits;
            inikeys = o.inikeys;
            tp = o.tp;
            tz_chisq = o.tz_chisq;
            max_states = o.max_states;
            init = true;
          }

        if(o.reordered)
          {
            reordered = true;
            reorder = o.reorder;
          }

        if(o.mass_fit)
          mass_fit = true;

        if(o.mass_fit_t0)
          {
            mass = o.mass;
            mass_chisq = o.mass_chisq;
            mass_summary = o.mass_summary;
            mass_plots = o.mass_plots;
            mass_fit_t0 = true;
          }

        if(o.z_fit)
          z_fit = true;

        if(o.z_fit_t0)
          {
            Z = o.Z;
            Z_chisq = o.Z_chisq;
            Z_summary = o.Z_summary;
            Z_plots = o.Z_plots;
            z_fit_t0 = true;
          }
      }

    return *this;
  }

  template<class T>
  typename std::vector<SembleMatrix<T> > SMT0Fit<T>::maketp(const typename PromoteCorr<T>::Type &tp_, const FitIniProps_t &inikeys_)
  {
    typename std::vector<SembleMatrix<T> > dum;

    for(int t = 0; t <= inikeys_.globalProps.tmax; ++t)
      dum.push_back(tp_.getCt(t));

    return dum;
  }

  //this is where most of the heavy lifting is done
  template<class T>
  void SMT0Fit<T>::load(const typename std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_)
  {
    clear();
    tp = tp_;
    inikeys = inikeys_;
    init = true;
    int t0;

    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        ST0Fit<T> *ptr = new ST0Fit<T>(t0, tp, inikeys);
        Handle<ST0Fit<T> > fit(ptr);
        t0_fits.insert(make_pair(t0, fit));

        t0_fits[t0]->fitPrinCorr();

        if(inikeys.globalProps.verbose)
          t0_fits[t0]->printSpectrum(std::cout);

        if(inikeys.outputProps.logs)
          {
            #pragma omp critical
            {
              fitlog << "\n\n t0 = " << t0 << "\n";
              t0_fits[t0]->printSpectrum(fitlog);
            }
          }

        if(inikeys.reconProps.recon)
          {
            tz_chisq.insert(make_pair(t0, t0_fits[t0]->findBestTZ(inikeys.reconProps)));
            t0_fits[t0]->clearRecon();
          }
        else //fake a map
          {
            tz_chisq.insert(make_pair(t0, std::pair<int, double>(inikeys.t0Props.t0ref + 1, 1.e7)));
          }
      }

    mass_fit = true;

    //pick a ref, most states or best chisq
    int max = t0_fits[inikeys.t0Props.t0low]->getM(), t0_ref = inikeys.t0Props.t0ref;

    for(int t = inikeys.t0Props.t0low + 1; t <= inikeys.t0Props.t0high; ++t)
      if(t0_fits[t]->getM() > max)
        max = t0_fits[t]->getM();

    max_states = max;

    //rephase and reorder the states
    SembleMatrix<T> Ct0_half, U, V;
    SembleVector<T> s;

    svd(t0_fits[t0_ref]->getCt0(), U, s, V);
    Ct0_half = U * sqrt(s);
    V = t0_fits[t0_ref]->peekVecs(tz_chisq[t0_ref].first);

    #pragma omp parallel for private(t0) shared(std::cout,Ct0_half,U,V)

    for(int t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        SembleMatrix<T> metric, Up, Vp;
        SembleVector<T> sp;

        svd(t0_fits[t0]->getCt0(), Up, sp, Vp);
        metric = Ct0_half * sqrt(sp) * adj(Up);
        Vp = metric * t0_fits[t0]->peekVecs(tz_chisq[t0].first);

        t0_fits[t0]->rephase(rephaseEigenVectorsEnsembleMap(V, Vp));
      }

    reorderStates();

    if(inikeys.t0FitProps.MT0)
      fitMassT0();

    if(inikeys.zProps.fit)
      {
        #pragma omp parallel for private(t0)

        for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
          {
            t0_fits[t0]->makeZ();
            t0_fits[t0]->fitZ();
          }

        z_fit = true;


        if(inikeys.t0FitProps.ZT0)
          fitZT0();

      }

    if(inikeys.reconProps.recon && inikeys.outputProps.logs)   //write out a recon_tz log
      {
        std::stringstream ss;
        ss << "multi_t0_fits";
        std::string path = SEMBLEIO::getPath() += ss.str();
        SEMBLEIO::makeDirectoryPath(path);
        path += std::string("/tz_chisq.log");
        std::ofstream out;
        out.open(path.c_str());
        out << "t_0 tz_best tz_chisq\n";

        for(int t = inikeys.t0Props.t0low; t <= inikeys.t0Props.t0high; ++t)
          out << t << " " << tz_chisq[t].second << "\n";

        out.close();
      }

    printer();
  }

//fit masses across t0s
  template<class T>
  void SMT0Fit<T>::fitMassT0(void)
  {

    if(mass_fit_t0)
      return;

    if(!!!mass_fit)
      {
	load(tp,inikeys);
      }

    std::cout << "fitting m(t0).." << std::endl;

    EnsemReal zero;
    zero.resize(t0_fits[inikeys.t0Props.t0low]->getB());
    zero = Real(0.);

    mass.resize(max_states, zero);
    mass_chisq.resize(max_states, 0.);
    mass_summary.resize(max_states, std::string(""));
    mass_plots = mass_summary;


    int state;
    #pragma omp parallel for private(state)

    for(state = 0; state < max_states; ++state)
      {
        std::vector<double> t;
        std::vector<EnsemReal> _data;
        int count = 0;

        for(int t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
          {
            int mapped;

            if(reorder[t0][state] == -1) //doesnt exist
              continue;
            else
              mapped = reorder[t0][state];

            if(t0_fits[t0]->getPCorrChiSq(mapped) < inikeys.prinCorrProps.accChisq)
              {
                t.push_back(t0);
                _data.push_back(t0_fits[t0]->peekMass0(mapped));
                ++count;
              }
          }

        if(count <= 1)
          {
            if(count == 0)
              mass_summary[state] = "no data";
            else
              {
                mass[state] = _data[0];
                mass_summary[state] = "one point";
              }
          }
        else
          {
            EnsemVectorReal data;
            data.resize(t0_fits[inikeys.t0Props.t0low]->getB());
            data.resizeObs(count);

            for(int index = 0; index < count; ++index)
              pokeObs(data, _data[index], index);

            FitVersusT0 mfit(EnsemData(t, data), fitComp(inikeys.t0FitProps.fitCrit));

            mass[state] = mfit.getConst();
            mass_summary[state] = mfit.getFitSummary();
            mass_chisq[state] = mfit.getChisq() / mfit.getNDoF();
            mass_plots[state] = mfit.getFitPlotString();
          }
      }//end parallel

    mass_fit_t0 = true;
  }

  template<class T>
  void SMT0Fit<T>::printMass(void)
  {
    if(!!!mass_fit)
      {
        std::cout << "no masses have been fit " << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

    int t0;

    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      if(reordered)
        t0_fits[t0]->printPrinCorrReorder(reorder[t0]);
      else
        t0_fits[t0]->printPrinCorrFits();
  }

  template<class T>
  void SMT0Fit<T>::printMassFiles(void)
  {
    if(!!!mass_fit)
      {
        std::cout << "no masses have been fit " << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

    int t0;

    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      if(reordered)
        t0_fits[t0]->printMassFileReorder(reorder[t0]);
      else
        t0_fits[t0]->printMassFile();

  }

  template<class T>
  void SMT0Fit<T>::printMassT0Files(void)
  {
    if(!!!mass_fit_t0)
      fitMassT0();

    std::stringstream ss;
    ss << "multi_t0_fits";
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/Mt0JackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int state = 0; state < mass.size(); ++state)
      {
        std::ostringstream file;
        file << path << "mass_fit_state" << state << ".jack";
        write(file.str(), mass[state]);
      }
  }


  template<class T>
  void SMT0Fit<T>::printEvalFiles(void)
  {
    if(!!!init)
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " not initialized, exiting.." << std::endl;
        exit(1);
      }

    int t0;
    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        if(reordered)
          t0_fits[t0]->printPrinCorrFilesReorder(reorder[t0]);
        else
          t0_fits[t0]->printPrinCorrFiles();
      }

  }


  template<class T>
  void SMT0Fit<T>::printMassT0(void)
  {
    if(!!!mass_fit_t0)
      fitMassT0();

    std::stringstream ss;
    ss << "multi_t0_fits";
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/PrinCorrPlots");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    for(int state = 0; state < max_states; ++state)
      {
        std::stringstream s;
        s <<  path << "multi_t0_prin_corr_fit_state" << state << ".ax" ;
        std::ofstream out;
        out.open(s.str().c_str());
        out << mass_plots[state];
        out.close();
      }
  }

  //fit Z across t0s
  template<class T>
  void SMT0Fit<T>::fitZT0(void)
  {

    if(z_fit_t0)
      return;

    if(!!!z_fit)
      {
        int t0;
        #pragma omp parallel for private(t0)

        for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
          {
            t0_fits[t0]->makeZ();
            t0_fits[t0]->fitZ();
          }

        z_fit = true;
        printZ();
      }

    std::cout << "fitting z(t0).." << std::endl;

    typename PromoteEnsem<T>::Type zero;
    zero.resize(t0_fits[inikeys.t0Props.t0low]->getB());
    zero = toScalar(0.);

    int nops = t0_fits[inikeys.t0Props.t0low]->getN();
    Z.resize(max_states, std::vector<typename PromoteEnsem<T>::Type>(nops, zero));
    Z_chisq.resize(max_states, std::vector<double>(nops, 0.));
    Z_summary.resize(max_states, std::vector<std::string>(nops, std::string("")));
    Z_plots = Z_summary;

    int state = 0;

    #pragma omp parallel for private(state) shared(std::cout)

    for(state = 0; state < max_states; ++state)
      {
        if(inikeys.globalProps.verbose)
          std::cout << "fitting state " << state << " of " << max_states << std::endl;

        for(int op = 0; op < nops; ++op)
          {
            std::vector<double> t;
            std::vector<typename PromoteEnsem<T>::Type> _data;
            int count = 0;

            for(int t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
              {
                int mapped;

                if(reorder[t0][state] == -1)
                  continue;
                else
                  mapped = reorder[t0][state];

                if(t0_fits[t0]->getZFitChiSq(mapped, op) < inikeys.zProps.accChisq)
                  {
                    t.push_back(t0);
                    _data.push_back(t0_fits[t0]->peekZ(mapped, op));
                    ++count;
                  }
              }

            if(count <= 1)
              {
                if(count == 0)
                  Z_summary[state][op] = "no data";
                else
                  {
                    Z[state][op] = _data[0];
                    Z_summary[state][op] = "one point";
                  }
              }
            else
              {
                typename PromoteEnsemVec<T>::Type data;
                data.resize(t0_fits[inikeys.t0Props.t0low]->getB());
                data.resizeObs(count);

                for(int index = 0; index < count; ++index)
                  pokeObs(data, _data[index], index);

                FitVersusT0 zfit(EnsemData(t, data), fitComp(inikeys.t0FitProps.fitCrit));

                Z[state][op] = zfit.getConst();
                Z_summary[state][op] = zfit.getFitSummary();
                Z_chisq[state][op]  = zfit.getChisq() / zfit.getNDoF();
                Z_plots[state][op] = zfit.getFitPlotString();
              }
          }
      }//end parallel

    z_fit_t0 = true;
  }

  template<class T>
  void SMT0Fit<T>::printZ(void)
  {
    if(!!!z_fit)
      {
        std::cout << "no overlaps have been fit" << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

    int t0;
    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        if(reordered)
          t0_fits[t0]->printZFitsReorder(reorder[t0]);
        else
          t0_fits[t0]->printZFits();
      }
  }

  template<class T>
  void SMT0Fit<T>::printZFiles(void)
  {
    if(!!!z_fit)
      {
        std::cout << "no overlaps have been fit" << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

    int t0;
    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        if(reordered)
          t0_fits[t0]->printZFileReorder(reorder[t0]);
        else
          t0_fits[t0]->printZFile();
      }
  }


  template<class T>
  void SMT0Fit<T>::printZtFiles(void)
  {
    if(!!!z_fit)
      {
        std::cout << "no overlaps have been fit" << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

    int t0;

    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        if(reordered)
          t0_fits[t0]->printZtReorder(reorder[t0]);
        else
          t0_fits[t0]->printZt();
      }
  }


  template<class T>
  void SMT0Fit<T>::printZT0(void)
  {
    if(!!!z_fit_t0)
      fitZT0();

    std::stringstream ss;
    ss << "multi_t0_fits";
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZFitPlots");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    int nops = t0_fits[inikeys.t0Props.t0low]->getN();

    for(int index = 0; index < max_states; ++index)
      {
        for(int op = 0; op < nops; ++op)
          {
            std::stringstream s;
            s << path << "multi_t0_Z_fit_state" << index << "_op" << op << ".ax";
            std::ofstream out;
            out.open(s.str().c_str());
            out << Z_plots[index][op];
            out.close();
          }
      }
  }


  template<class T>
  void SMT0Fit<T>::printZT0Files(void)
  {
    if(!!!z_fit_t0)
      fitZT0();

    std::stringstream ss;
    ss << "multi_t0_fits";
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/ZT0JackFiles");
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    int nops = t0_fits[inikeys.t0Props.t0low]->getN();

    for(int index = 0; index < max_states; ++index)
      for(int op = 0; op < nops; ++op)
        {
          std::ostringstream file;
          file << path << "z_fit_state" << index << "_op" << op << ".jack";
          write(file.str(), Z[index][op]);
        }
  }


  template<class T>
  void SMT0Fit<T>::printVtFiles(void)
  {
    if(!!!init)
      {
        std::cout << "not initialized" << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

    int t0;
    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        if(reordered)
          t0_fits[t0]->printVtReorder(reorder[t0]);
        else
          t0_fits[t0]->printVt();
      }
  }

  template<class T>
  void SMT0Fit<T>::printReorderLog(void)
  {
    std::stringstream sp;
    sp << "multi_t0_fits";
    std::string path = SEMBLEIO::getPath() += sp.str();
    SEMBLEIO::makeDirectoryPath(path);

    std::stringstream ss;

    std::map<int, std::map<int, int> >::const_iterator mapit;
    std::map<int, int>::const_iterator it;

    SembleMatrix<T> Ct0_half, U, V;
    SembleVector<T> s;

    svd(t0_fits[inikeys.t0Props.t0ref]->getCt0(), U, s, V);
    Ct0_half = U * sqrt(s);
    V = t0_fits[inikeys.t0Props.t0ref]->peekVecs(tz_chisq[inikeys.t0Props.t0ref].first);

    for(mapit = reorder.begin(); mapit != reorder.end(); ++mapit)
      {
        ss << "t0ref = " << inikeys.t0Props.t0ref << " t0 = " << mapit->first << "\n";

        for(it = mapit->second.begin(); it != mapit->second.end(); ++it)
          {
            ss << "ref state " << it->first << " is t0 state " << it->second << "\n";
          }

        ss << "\n\n mean ovelaps adj(ref)*metric*t0state \n\n";

        SembleMatrix<T> metric, Up, Vp;
        SembleVector<T> sp;

        svd(t0_fits[mapit->first]->getCt0(), Up, sp, Vp);
        metric = Ct0_half * sqrt(sp) * adj(Up);
        Vp = metric * t0_fits[mapit->first]->peekVecs(tz_chisq[mapit->first].first);

        ss << mean(adj(V)*Vp);

        ss << "\n \n \n \n";
      }

    path += std::string("/reorder_log");
    std::ofstream out;
    out.open(path.c_str());
    out << ss.str();
    out.close();
  }

  template<class T>
  void SMT0Fit<T>::printNResetLog(void)
  {
    std::stringstream ss;
    ss << "multi_t0_fits";
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);

    std::stringstream s;
    typename std::map<int, Handle<ST0Fit<T> > >::const_iterator it;

    for(it = t0_fits.begin(); it != t0_fits.end(); ++it)
      s << (*it).first << " " << (*it).second->getN() - (*it).second->getM() << "\n";

    path += std::string("/nreset_log");
    std::ofstream out;
    out.open(path.c_str());
    out << s.str();
    out.close();
  }

  template<class T>
  void SMT0Fit<T>::printReconPlots(void)
  {
    if(!!!z_fit)
      {
std:
        cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " can't recon without overlap fitting" << std::endl;
        exit(1);
      }

    if(inikeys.globalProps.verbose)
      std::cout << "Making Recon Plots.." << std::endl;

    int t0;
    #pragma omp parallel for private(t0)

    for(t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        t0_fits[t0]->printReconPlots();
      }
  }

  template<class T>
  void SMT0Fit<T>::reorderStates(void)
  {
    if(!!!init)
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " can't reorder w/o init" << std::endl;
        exit(1);
      }

    reorder.clear();

    const int nops = tp.end()->getN();
    const int B = tp.end()->getB();
    const int nvecs = max_states;
    const int t0_ref = inikeys.t0Props.t0ref;

    typename std::map<int, SembleMatrix<T> > metric;

    //determine the metric for each t0 once
    for(int t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
        SembleMatrix<T> U, V;
        SembleVector<T> s;
        svd(t0_fits[t0]->getCt0(), U, s, V);
        metric.insert(std::make_pair(t0, sqrt(s)*adj(U)));
      }

    std::vector<int> zero(nvecs,0);
    std::vector<std::vector<int> > dum(nvecs,zero);
    std::vector<std::vector<std::vector<int> > > counts(inikeys.t0Props.t0high - inikeys.t0Props.t0low + 1, dum);

    //use counts to keep a tally for a maximum likelyhood overlap map
    //outter vector is t0, inner vector is the ref state, innermost vector is the state on that t0

    //now lets loop over each t0 value
    for(int t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
      {
	if(t == t0_ref) //no eigenvectors, don't want to contaminate map with junk
	  continue;

	const SembleMatrix<T> W_t_ref = metric[t0_ref]*t0_fits[t0_ref]->peekVecs(t);

	for(int t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
	  {
	    if(t == t0) //no eigenvectors, don't want to contaminate map with junk
	      continue;
	    
	    std::map<int,int> mapp = makeRemap(W_t_ref,metric[t0]*t0_fits[t0]->peekVecs(t),nvecs);
	    std::map<int,int>::const_iterator it;

	    for(it = mapp.begin(); it != mapp.end(); ++it)
	      if(it->second == -1)    
		continue;                             //no mapping
	      else
		++counts[t0 - inikeys.t0Props.t0low][it->first][it->second];  //increment
	  }
      }

    const int refdim = t0_fits[t0_ref]->getM();

    for(int t0 = inikeys.t0Props.t0low; t0 <= inikeys.t0Props.t0high; ++t0)
      {
	const int vecdim = t0_fits[t0]->getM();
	const int bound = (refdim < vecdim) ? refdim : vecdim;
	std::vector<bool> ur(nvecs,false),uv(nvecs,false);
	int mr,mv,m;
	std::map<int,int> rmap_t0;

	dum = counts[t0 - inikeys.t0Props.t0low];

	for(int state = 0; state < bound; ++state)
	  {
	    mr = 0;
	    mv = 0;
	    m = 0;

	    for(int r = 0; r < nvecs; ++r)
	      if(ur[r]) 
		continue;
	      else
		for(int v = 0; v < nvecs; ++v)
		  {
		    if(uv[v])
		      continue;
		    
		    if(dum[r][v] > m)
		      {
			m = dum[r][v];
			mr = r;
			mv = v;
		      }
		  }

	    ur[mr] = true;
	    uv[mv] = true;
	    rmap_t0[mr] = mv;
	  }


	/*
	    CASES
	           1) refdim = vecdim - do nothing
		   2) refdim < vecdim - pile anything that didn't map on top with a lowest to lowest convention
		   3) refdim > vecdim - don't need to do anything, the map takes care of itself
	 */

	if(refdim < vecdim)
	  {
	    ur.clear();
	    uv.clear();

	    ur.resize(vecdim,false);
	    uv.resize(vecdim,false);
	    
	    std::map<int,int>::const_iterator it;
	    
	    for(it = rmap_t0.begin(); it != rmap_t0.end(); ++it)
	      {
		ur[it->first] = true;
		uv[it->second] = true;
	      }

	    std::vector<int> rr,vv;

	    for(int i = 0; i < vecdim; ++i)
	      {
		if(ur[i])
		  rr.push_back(i);
		if(uv[i])
		  vv.push_back(i);
	      }

	    if(rr.size() != vv.size())
	      {
		std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
		std::cout << " If you see this something went horribly wrong" << std::endl;
		exit(1);
	      }

	    for(int i = 0; i < rr.size(); ++i)
	      rmap_t0[rr[i]] = vv[i];
	  }

	reorder[t0] = rmap_t0;

      }//end t0 loop

    reordered = true;
  }


  template<class T>
  void SMT0Fit<T>::clear(void)
  {
    tz_chisq.clear();
    t0_fits.clear();
    reorder.clear();
    tp.clear();
    mass.clear();
    mass_chisq.clear();
    mass_summary.clear();
    mass_plots.clear();
    Z.clear();
    Z_chisq.clear();
    Z_summary.clear();
    Z_plots.clear();
    max_states = -1;
    init = false;
    mass_fit = false;
    mass_fit_t0 = false;
    z_fit = false;
    z_fit_t0 = false;
    reordered = false;
  }



  template<class T>
  Handle<FitComparator> SMT0Fit<T>::fitComp(const std::string &in) const
  {
    if(in == "chisq_per_dof")
      return Handle<FitComparator>(new CompareFitsByChisqPerNDoF);

    if(in == "Q")
      return Handle<FitComparator>(new CompareFitsByQ);

    if(in == "QN")
      return Handle<FitComparator>(new CompareFitsByQN);

    std::cout << __PRETTY_FUNCTION__ << "Fit Criterion unknown, defaulting to chisq_per_dof" << std::endl;
    return Handle<FitComparator>(new CompareFitsByChisqPerNDoF);
  }

  template<class T>
  void SMT0Fit<T>::printer(void)
  {
    OutputProps_t out = inikeys.outputProps;

    if(out.mass)
      printMassFiles();

    if(out.Z_t)
      printZtFiles();

    if(out.V_t)
      printVtFiles();

    if(out.pcorrFiles)
      printEvalFiles();

    if(out.zFitFiles)
      printZFiles();

    if(out.pcorrPlots)
      printMass();

    if(out.zFitPlots)
      printZ();

    if(out.mT0Files)
      printMassT0Files();

    if(out.zT0Files)
      printZT0Files();

    if(out.mT0Plots)
      printMassT0();

    if(out.zT0Plots)
      printZT0();

    if(out.reconPlots)
      printReconPlots();

    if(out.logs)
      {
        printReorderLog();
        printNResetLog();
        std::ofstream out;
        out.open("ini_file");
        out << inikeys;
        out.close();
        out.open("fitlog");
        out << fitlog.str();
        out.close();
      }

  }


}
#endif
