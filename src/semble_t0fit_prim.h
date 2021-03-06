#ifndef SEMBLE_T0FIT_PRIM_H_H_GUARD
#define SEMBLE_T0FIT_PRIM_H_H_GUARD

#include"ensem/ensem.h"
#include"semble_fit_ini_xml.h"
#include"semble/semble_file_management.h"
#include"semble/semble_meta.h"
#include"semble/semble_matrix.h"
#include"semble/semble_vector.h"
#include"semble/semble_linear_algebra.h"
#include"semble_load_correlators.h"
#include"semble_t0fit_base.h"
#include"semble/semble_histogram.h"
#include"semble/semble_histogram_aux.h"
#include<string>
#include<vector>
#include<iostream>
#include<sstream>


//this header contains the ST0FitPrim class that hides the choice of generalized eigen problem solution and is
//a member of the ST0Fit class, it only knows how to take a vector of correlation matricies and give back
//sorted generalized eigenvalues/vectors

namespace SEMBLE
{

//enumerate the gen eig solution choices
/////////////////////////////////////////////////////////////////////////////////////////

  enum eload_gen_eig_string //put this string into inikeys.gen_eigen_props.inversionType
  {
    eCho,           //Cholesky
    eSvdCond,       //Svd on condition number
    eSvdValue,      //Svd on value
    eSvdSigma,      //Svd on num sigma away from zero
    eSvdSigmaValue,
    eSvdSigmaCond,
    eGenErr
  };


  //enumerate the sorting props
  enum eload_sorting_props_stringCfg
  {
    eNoneCfg,      //sort by eigenvalues
    eRefvecsCfg,   //sort by orthogonality of eigenvectors
  };


  //enumerate the sort evecs tslice choices
  enum eload_sort_evecs_tslice
  {
    eNone,                    //no sorting across timeslice
    eRefvecs_Moving,          //sort t+1 by the eigenvectors of t
    eRefvecs_Fixed,           //sort t_ref+n by the eigenvectors of t_ref
    eRefvecs_Fixed_Auto       //sort by the auto ref with vixed vectors
  };


//ST0FitPrim<T> class definition
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class T>
  struct ST0FitPrim
  {
  public: //constructors, destructor, copy assignment
    ST0FitPrim(void);
    ST0FitPrim(const int t0_, const typename std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_);
    ST0FitPrim(const ST0FitPrim<T> &o);
    ~ST0FitPrim(void); // needs to delete base pointers
    ST0FitPrim<T>& operator=(const ST0FitPrim<T> &o);

  public: //const access
    int getN(void) const                 //return the number of rows(ops)
    {
      initChk();
      return N;
    }
    int getM(void) const                 //return the number of cols(states)
    {
      initChk();
      return M;
    }
    int getB(void) const
    {
      initChk();
      return B;
    }
    inline int rows(void) const
    {
      return getN();
    }
    inline int cols(void) const
    {
      return getM();
    }
    inline int bins(void) const
    {
      return getB();
    }
    SembleMatrix<T> getCt0(void) const
    {
      initChk();
      return Ct0;
    }

  public: //non const access
    void load(const int t0_, const std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_);
    void solve(void);                                                    //solve the problem
    void sort_solved(void);                                              //sort the evecs/vals
    void clear(void);                                                    //clean up to load a new problem

  private:
    void clear_data(void);                                               //delete pointers then clear vector _data
    void initChk(void) const;
    void solveChk(void) const;
    bool sortEvecsCfg(void) const;
    int findRefT(void) const;
    void makeSVDHistos(const SembleVector<double> &svals) const;
    eload_gen_eig_string  lges_hash(const std::string &in) const;
    eload_sorting_props_stringCfg lspscfg_hash(const std::string &in) const;
    eload_sort_evecs_tslice lset_hash(const std::string &in) const; 

  public:                                                                 //will solve if initialized but not solved
    std::vector<SembleVector<double> > getEvals(void);                    //indexed [0,inikeys.globalProps.tmax]
    typename std::vector<SembleMatrix<T> > getEvecs(void);

  private: //data store
    int N, M, B, t0;
    typename std::vector<ST0Base<T>* > _data;
    bool init, solved;
    SembleMatrix<T> Ct0;
    FitIniProps_t inikeys;
    std::stringstream svdRematchingLog;
  };

  //ST0FitPrim<T> Implementation
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<class T>
  ST0FitPrim<T>::ST0FitPrim(void)
    : N(0), M(0), B(0), t0(-1), init(false) , solved(false)
  {}

  template<class T>
  ST0FitPrim<T>::ST0FitPrim(const int t0_, const typename std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_)
    : N(0), M(0), B(0), t0(-1), init(false) , solved(false)
  {
    load(t0_, tp_, inikeys_);
    solve();
    sort_solved();
  }

  template<class T>
  ST0FitPrim<T>::ST0FitPrim(const ST0FitPrim<T> &o)
    : N(o.N) , M(o.M), B(o.B), t0(o.t0), init(o.init), solved(o.solved)
  {
    if(init)
      {
        Ct0 = o.Ct0;
        inikeys = o.inikeys;
        clear_data();
        typename std::vector<ST0Base<T>* >::const_iterator it;

        for(it = o._data.begin(); it != o._data.end(); ++it)
          _data.push_back((*it)->clone());
      }
  }

  template<class T>   //we need to clear the vector of pointers to avoid mem leak
  ST0FitPrim<T>::~ST0FitPrim(void)
  {
    clear_data();
  }

  template<class T>
  ST0FitPrim<T>& ST0FitPrim<T>::operator=(const ST0FitPrim<T> &o)
  {
    if(this != &o)
      {
        N = o.N;
        M = o.M;
        B = o.B;
        t0 = o.t0;
        init = o.init;
        solved = o.solved;

        if(init)
          {
            Ct0 = o.Ct0;
            inikeys = o.inikeys;
            clear_data();
            typename std::vector<ST0Base<T>* >::const_iterator it;

            for(it = o._data.begin(); it != o._data.end(); ++it)
              _data.push_back((*it)->clone());
          }
      }

    return *this;
  }


  //public non const
  template<class T>
  void ST0FitPrim<T>::load(const int t0_, const std::vector<SembleMatrix<T> > &tp_, const FitIniProps_t &inikeys_)
  {
    if(init)
      clear_data();

    init = true;
    solved = false;

    t0 = t0_;
    inikeys = inikeys_;
    bool sort = sortEvecsCfg();
    SembleMatrix<T> Cts, F, U, Ud, sCt,Ur,Udr,Ct0r;
    SembleMatrix<double> RSP;
    SembleVector<double> rsp,sr;
    double thresh, sigma;

    Cts = tp_[t0];
    Ct0 = symmetrize(Cts);

    N = Ct0.getN();
    M = Ct0.getM();
    B = Ct0.getB();


    std::ofstream out;
    std::stringstream ss;
    std::string path;

    //fill with the base class to solve the generalized eigen problem

    //NB for svd if we have degenerate/near degenerate states then there could be an 
    //ordering ambiguity across multiple values of t0, there also will always be a 
    //multi-t0 phasing ambiguity in the matrix U, to enforce a consistent phase
    //and ordering we can just match the U and Sigma to t0_ref's ordering
    //and phasing conventions, we do this for each svd case below before we 
    //do our resetting and move into the reduced subspace

    _data.resize(inikeys.globalProps.tmax - inikeys.globalProps.tmin + 1); 
    int t;
    int toff = inikeys.globalProps.tmin; // the jlab gcc is stupid and doesnt like const vars in omp 

    switch(lges_hash(inikeys.genEigProps.type))
    {
      case eCho:
        if(!!!cholesky(Ct0, F))
        {
          std::cout << "WARNING CHOLESKY FACTORIZATION FAILED ON t0 = " << t0 << "\n EXITING..." << std::endl;
          exit(1);
        }
#pragma omp parallel for shared(t,toff) schedule(dynamic)
        for(t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
        {
          Cts = tp_[t];
          sCt = symmetrize(Cts);
          _data[t - toff] = new ChoF<double>(sCt, F, t, sort);
        }

        break;

      case eSvdCond:

        thresh = inikeys.genEigProps.thresh;
        svdRematchingLog << "--\n--\n--\ndoing svd on Ct0\n" << svd(Ct0, U, rsp, Ud);

        path = SEMBLEIO::getPath() + std::string("/SVDLogs/");
        SEMBLEIO::makeDirectoryPath(path);

        ss << path << "svdRematchingLog_t0" << t0 << ".log";
        out.open(ss.str().c_str());
        out << svdRematchingLog.str();
        out.close();

        if(inikeys.genEigProps.svdHisto)
          makeSVDHistos(rsp);

        M = svdResetAverageCond(rsp, thresh);
        svdResetPseudoInvertRoot(rsp, M);
        RSP = diag(rsp);
        RSP.rows(M);

#pragma omp parallel for shared(t,toff) schedule(dynamic)
        for(t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
        {
          Cts = tp_[t];
          sCt = symmetrize(Cts);
          _data[t-toff] = new SvdF<double>(sCt, RSP, U, t, sort);
        }

        break;

      case eSvdValue:
        thresh = inikeys.genEigProps.thresh;
        svdRematchingLog << "--\n--\n--\ndoing svd on Ct0\n" << svd(Ct0, U, rsp, Ud);

        path = SEMBLEIO::getPath() + std::string("/SVDLogs/");
        SEMBLEIO::makeDirectoryPath(path);

        ss << path << "svdRematchingLog_t0" << t0 << ".log";
        out.open(ss.str().c_str());
        out << svdRematchingLog.str();
        out.close();

        if(inikeys.genEigProps.svdHisto)
          makeSVDHistos(rsp);

        M = svdResetAverageValue(rsp, thresh);
        svdResetPseudoInvertRoot(rsp, M);
        RSP = diag(rsp);
        RSP.rows(M);

#pragma omp parallel for shared(t,toff) schedule(dynamic)
        for(t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
        {
          Cts = tp_[t];
          sCt = symmetrize(Cts);
          _data[t-toff] = new SvdF<double>(sCt, RSP, U, t, sort);
        }

        break;

      case eSvdSigma:
        sigma = inikeys.genEigProps.sigma;
        svdRematchingLog << "--\n--\n--\ndoing svd on Ct0\n" << svd(Ct0, U, rsp, Ud);

        path = SEMBLEIO::getPath() + std::string("/SVDLogs/");
        SEMBLEIO::makeDirectoryPath(path);

        ss << path << "svdRematchingLog_t0" << t0 << ".log";
        out.open(ss.str().c_str());
        out << svdRematchingLog.str();
        out.close();

        if(inikeys.genEigProps.svdHisto)
          makeSVDHistos(rsp);

        M = svdResetSigma(rsp, U, sigma);
        svdResetPseudoInvertRoot(rsp, M);
        RSP = diag(rsp);
        RSP.rows(M);

#pragma omp parallel for shared(t,toff) schedule(dynamic)
        for(t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
        {
          Cts = tp_[t];
          sCt = symmetrize(Cts);
          _data[t-toff] = new SvdF<double>(sCt, RSP, U, t, sort);
        }

        break;

      case eSvdSigmaValue:
        sigma = inikeys.genEigProps.sigma;
        thresh = inikeys.genEigProps.thresh;
        svdRematchingLog << "--\n--\n--\ndoing svd on Ct0\n" << svd(Ct0, U, rsp, Ud);

        path = SEMBLEIO::getPath() + std::string("/SVDLogs/");
        SEMBLEIO::makeDirectoryPath(path);

        ss << path << "svdRematchingLog_t0" << t0 << ".log";
        out.open(ss.str().c_str());
        out << svdRematchingLog.str();
        out.close();

        if(inikeys.genEigProps.svdHisto)
          makeSVDHistos(rsp);

        M = svdResetAvgValueAndSigma(rsp, U, thresh, sigma);
        svdResetPseudoInvertRoot(rsp, M);
        RSP = diag(rsp);
        RSP.rows(M);

#pragma omp parallel for shared(t,toff) schedule(dynamic)
        for(t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
        {
          Cts = tp_[t];
          sCt = symmetrize(Cts);
          _data[t-toff] = new SvdF<double>(sCt, RSP, U, t, sort);
        }

        break;


      case eSvdSigmaCond:
        sigma = inikeys.genEigProps.sigma;
        thresh = inikeys.genEigProps.thresh;
        svdRematchingLog << "--\n--\n--\ndoing svd on Ct0\n" << svd(Ct0, U, rsp, Ud);

        path = SEMBLEIO::getPath() + std::string("/SVDLogs/");
        SEMBLEIO::makeDirectoryPath(path);

        ss << path << "svdRematchingLog_t0" << t0 << ".log";
        out.open(ss.str().c_str());
        out << svdRematchingLog.str();
        out.close();

        if(inikeys.genEigProps.svdHisto)
          makeSVDHistos(rsp);

        M = svdResetAvgCondAndSigma(rsp, U, thresh, sigma);
        svdResetPseudoInvertRoot(rsp, M);
        RSP = diag(rsp);
        RSP.rows(M);

#pragma omp parallel for shared(t,toff) schedule(dynamic)
        for(t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
        {
          Cts = tp_[t];
          sCt = symmetrize(Cts);
          _data[t-toff] = new SvdF<double>(sCt, RSP, U, t, sort);
        }

        break;


      case eGenErr:
        std::cout << "Unsupported generalized eigen problem solution type, defaulting to Cholesky" << std::endl;
        cholesky(Ct0, F);

#pragma omp parallel for shared(t,toff) schedule(dynamic)
        for(t = inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax; ++t)
        {
          Cts = tp_[t];
          sCt = symmetrize(Cts);
          _data[t-toff] = new ChoF<double>(sCt, F, t, sort);
        }

    }//end switch



  }//end load

  template<class T>
    void ST0FitPrim<T>::solve(void)
    {
      initChk();

      if(inikeys.globalProps.verbose)
        std::cout << "Solving the generalized eigen problem using " << inikeys.genEigProps.type << ".." << std::endl;

      if(!!!solved)
      {

        int t;
#pragma omp parallel for shared(t)
        for(t = 0; t < _data.size(); ++t)
          _data[t]->eval();

        solved = true;
      }


      typename std::vector<ST0Base<T>* >::const_iterator it; 

      // pump out log files -- single thread so we dont piss off the file system too much
      for(it = _data.begin(); it != _data.end(); ++it)
      {
        std::stringstream path;  
        path << SEMBLEIO::getPath() + std::string("/SVDLogs/");
        SEMBLEIO::makeDirectoryPath(path.str());
        path << "t_logs_t0" << t0 << "/";
        SEMBLEIO::makeDirectoryPath(path.str());
        path << "svdRematchingLog_t0" << t0 << "__t" << (*it)->get_t() << ".log";
        std::ofstream out;
        out.open(path.str().c_str());
        out << (*it)->get_log();
        out.close();
      }
    }


  // this looks like it should be parallelized -- someday CJS should do it
  template<class T>
    void ST0FitPrim<T>::sort_solved(void)   //NB _data[index], index runs [0,inikeys.globalProps.tmax] so need to include possible offsets
    {
      initChk();

      if(!!!solved)
        solve();

      if(inikeys.globalProps.verbose)
        std::cout << "sorting on timeslices using " << inikeys.sortingProps.sortEvecsTimeslice << std::endl;
    
      std::stringstream matchlog; 

      typename std::vector<ST0Base<T>* >::iterator it;
      ST0Base<T> *ptr, *ptrRef;

      switch(lset_hash(inikeys.sortingProps.sortEvecsTimeslice))
      {
        case eNone: 
          for(it = _data.begin(); it != _data.end(); ++it)
          {
            if((*it)->getT() >= t0)
              reorderEigenValues((*it)->evals(), (*it)->w(), (*it)->evecs(), false);
            else
              reorderEigenValues((*it)->evals(), (*it)->w(), (*it)->evecs(), true);
          }
          break;

        case eRefvecs_Moving: //phase starting with t0 -1 and carry it forward (backwards)

          //enforce the convention that the lowest state is indexed by 0
          ptrRef = _data[t0 - 1 - inikeys.globalProps.tmin];
          reorderEigenValues(ptrRef->evals(),ptrRef->w(),ptrRef->evecs(),true);


          //rephase/reoder from t0-1 to inikeys.globalProps.tmin
          for(int t = t0 - 2 - inikeys.globalProps.tmin; t >= 0; --t)
          {
            ptrRef = _data[t + 1];
            ptr = _data[t];
            matchlog << "--\n--\n--   tr = " << t + 1 + inikeys.globalProps.tmin
              << " t = " << t + inikeys.globalProps.tmin 
              << matchEigenVectorsEnsembleMetric(Ct0, ptrRef->evecs(), ptr->evecs(), ptr->evals());
          }

          //the gen eig solution is junk on t = t0 so skip it, put phase/order from t0-1 to t0+1
          matchlog << "--\n--\n--   tr = " << t0 - 1 << " t = " << t0 + 1  
            << matchEigenVectorsEnsembleMetric(Ct0,
                (_data[t0 - 1 - inikeys.globalProps.tmin])->evecs(),
                (_data[t0 + 1 - inikeys.globalProps.tmin])->evecs(),
                (_data[t0 + 1 - inikeys.globalProps.tmin])->evals());

          //carry the phase from t0+1 forward to inikeys.globalProps.tmax
          for(int t = t0 + 2 - inikeys.globalProps.tmin; t <= inikeys.globalProps.tmax - inikeys.globalProps.tmin; ++t)
          {
            ptrRef = _data[t - 1];
            ptr = _data[t];
            matchlog << "--\n--\n--   tr = " << t - 1 + inikeys.globalProps.tmin
              << " t = " << t + inikeys.globalProps.tmin 
              <<matchEigenVectorsEnsembleMetric(Ct0, ptrRef->evecs(), ptr->evecs(), ptr->evals());
          }

          break;

        case eRefvecs_Fixed:

          //enforce the convention that the lowest state is indexed by 0
          ptrRef = _data[t0 + inikeys.sortingProps.deltaRef - inikeys.globalProps.tmin];
          if(ptrRef->getT() >= t0)
            reorderEigenValues(ptrRef->evals(),ptrRef->w(),ptrRef->evecs(),false);
          else
            reorderEigenValues(ptrRef->evals(),ptrRef->w(),ptrRef->evecs(),true);

          for(int t = 0; t <= inikeys.globalProps.tmax  - inikeys.globalProps.tmin; ++t)
          {
            if(t == t0 - inikeys.globalProps.tmin) //skip b/c its junk
              continue;

            if(t == t0 + inikeys.sortingProps.deltaRef - inikeys.globalProps.tmin)  //t = tref, no need
              continue;

            ptr = _data[t];
            matchlog << "--\n--\n--   tr = " << t0 + inikeys.sortingProps.deltaRef
              << " t = " << t + inikeys.globalProps.tmin 
              << matchEigenVectorsEnsembleMetric(Ct0, ptrRef->evecs(), ptr->evecs(), ptr->evals());
          }

          break;

        case eRefvecs_Fixed_Auto:
          {
            //enforce the convention that the lowest state is indexed by 0
            int refT = findRefT(); 
            ptrRef = _data[refT - inikeys.globalProps.tmin];
            if(ptrRef->getT() >= t0)
              reorderEigenValues(ptrRef->evals(),ptrRef->w(),ptrRef->evecs(),false);
            else
              reorderEigenValues(ptrRef->evals(),ptrRef->w(),ptrRef->evecs(),true);

            for(int t = 0; t <= inikeys.globalProps.tmax - inikeys.globalProps.tmin; ++t)
            {
              if(t == t0  - inikeys.globalProps.tmin)
                continue;

              ptr = _data[t];
              matchlog << "--\n--\n--   tr = " << refT << " t = " << t + inikeys.globalProps.tmin 
                << matchEigenVectorsEnsembleMetric(Ct0, ptrRef->evecs(), ptr->evecs(), ptr->evals());
            }

            break;
          }
        default:
          std::cerr << "you open a door to a roomfull of orcs and are slain" 
            << "\n seriously though what did you do to get this error" << std::endl;
          exit(1); 
      }//end switch


        std::stringstream path;  
        path << SEMBLEIO::getPath() + std::string("/t0");
        path << t0 << "/";
        SEMBLEIO::makeDirectoryPath(path.str());
        path << "t_2_t_eigenvector_reordering" << ".log";
        std::ofstream out;
        out.open(path.str().c_str());
        out << matchlog.str();
        out.close();


    }

  template<class T>
    void ST0FitPrim<T>::clear(void)
    {
      clear_data();
      Ct0 = SembleMatrix<T>(1,1,1);
      N = M = B = t0 = 0;
      init = solved = false;
      svdRematchingLog.str(std::string());
    }

  //private

  template<class T>
    void ST0FitPrim<T>::clear_data(void)   //need to delete the pointers before clearing the vector..
    {
      typename std::vector<ST0Base<T>* >::iterator it;

      for(it = _data.begin(); it != _data.end(); ++it)
        delete *it;

      _data.clear();
    }

  template<class T>
    void ST0FitPrim<T>::initChk(void) const
    {
      if(!!!init)
      {
        std::cout << __PRETTY_FUNCTION__ << " need to initialize" << std::endl;
        exit(1);
      }
    }

  template<class T>
    void ST0FitPrim<T>::solveChk(void) const
    {
      if(!!!solved)
      {
        std::cout << __PRETTY_FUNCTION__ << " need to solve" << std::endl;
        exit(1);
      }
    }

  template<class T>
    bool ST0FitPrim<T>::sortEvecsCfg(void) const
    {

      if(!!!init)
      {
        std::cout << "Problem not initialized in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

      bool r = true;

      switch(lspscfg_hash(inikeys.sortingProps.sortEvecsCfg))
      {
        case eNoneCfg:
          r = false;
          break;
        case eRefvecsCfg: //this is the err condition in the 'hash'
          r = true;
          break;
      }

      return r;
    }

  template<class T>
    int ST0FitPrim<T>::findRefT(void) const
    {
      if(inikeys.globalProps.verbose)
        std::cout << "Finding a suitable reference t close to t0" << std::endl;

      typename std::map<int, SembleMatrix<T> > tvecs;
      typename std::vector<ST0Base<T>* >::const_iterator it;

      //now we have a map of where the key is time and the data is the generalized eigenvector
      for(it = _data.begin(); it != _data.end(); ++it)
        tvecs[(*it)->getT()] = (*it)->getEvecs();

      //sort all of the generalized eigenvectors that we have access too
      typename std::map<int, SembleMatrix<T> >::iterator mapit;
      typename std::map<int, SembleMatrix<T> >::const_iterator refit;
      SembleVector<T> dum = (*_data[0]).getEvals();
      refit = tvecs.begin();

      for(mapit = tvecs.begin(); mapit != tvecs.end(); ++mapit)
        matchEigenVectorsEnsembleMetric(Ct0, refit->second, mapit->second, dum);

      //now we have the overlaps and variance for each vector onto itself
      typename std::map<int, std::pair<itpp::Vec<T>, itpp::Vec<T> > > mandv;

      for(int t = inikeys.globalProps.tmin; t < inikeys.globalProps.tmax; ++t)
      {
        SembleVector<T> dum = diag(adj(tvecs[t]) * Ct0 * tvecs[t + 1]);
        mandv[t] = std::pair<itpp::Vec<T>, itpp::Vec<T> >(mean(dum), variance(dum));
      }

      //now lets find tbest
      int tbest = t0 + 1;
      double bcomp = itpp::sum(mandv[tbest - 2].second) + itpp::sum(mandv[tbest].second);

      // search +/- 7 tslices of t0
      for(int t = t0 - 7; t < t0 + 7; ++t)
      {
        if((t < inikeys.globalProps.tmin) || (t == t0) || (t > inikeys.globalProps.tmax))
          continue;

        double comp;
        int tm = t - 1;

        if(tm == t0)
          --tm;

        comp = itpp::sum(mandv[tm].second) + itpp::sum(mandv[t].second);

        if(comp < bcomp)
        {
          bcomp = comp;
          tbest = t;
        }
      }

      return tbest;
    }

  template<class T>
    void ST0FitPrim<T>::makeSVDHistos(const SembleVector<double> &svals) const
    {


      std::pair<std::vector<double>, std::vector<double> > range = findRange(svals);
      std::vector<int> bins(svals.getN(), inikeys.genEigProps.nHistoBins);

      SembleMultiHisto histo(range.first, range.second, bins);
      histo.Add(svals);

      std::vector<std::string> histograms = histo.genHisto();

      std::string path = SEMBLEIO::getPath();
      std::stringstream ss;
      ss << "t0" << t0 << "/";
      path += ss.str();
      SEMBLEIO::makeDirectoryPath(path);
      path += std::string("SVDHistos/");
      SEMBLEIO::makeDirectoryPath(path);

      for(int elem = 0; elem < histograms.size(); ++elem)
      {
        std::stringstream fname;
        std::ofstream out;
        fname << path << "svdHisto_t0" << t0 << "_singularvalue" << elem << ".txt";
        out.open(fname.str().c_str());
        out << histograms[elem];
        out.close();

        std::ostringstream file;
        file << path << "svdSingluarValue_" << elem << "_t0" << t0 << ".jack";
        write(file.str(), svals.getEnsemElement(elem));
      }
    }

  template<class T>   //hash the result for something useful in a switch
    eload_gen_eig_string ST0FitPrim<T>::lges_hash(const std::string &in) const
    {
      if(in == "Cho") return eCho;

      if(in == "Cholesky") return eCho;

      if(in == "SvdCond") return eSvdCond;

      if(in == "SvdValue") return eSvdValue;

      if(in == "SvdSigma") return eSvdSigma;

      if(in == "SvdSigmaValue") return eSvdSigmaValue;

      if(in == "SvdSigmaCond") return eSvdSigmaCond;

      std::cout << __PRETTY_FUNCTION__ << "key: " << in << " is not a supported type" << std::endl;

      return eGenErr;
    }

  template<class T> //hash the result
    eload_sorting_props_stringCfg ST0FitPrim<T>::lspscfg_hash(const std::string &in) const
    {
      if(in == "None") return eNoneCfg;

      if(in == "Refvecs") return eRefvecsCfg;

      std::cout << in << " is not a supported type, defaulting to Refvecs" << std::endl;

      return eRefvecsCfg; //default
    }

  template<class T> //hash the result
    eload_sort_evecs_tslice ST0FitPrim<T>::lset_hash(const std::string &in) const
    {
      if(in == "None") return eNone;

      if(in == "Refvecs_Moving") return eRefvecs_Moving;

      if(in == "Refvecs_Fixed") return eRefvecs_Fixed;

      if(in == "Refvecs_Fixed_Auto") return eRefvecs_Fixed_Auto;

      std::cout << in << " is not a supported type, defaulting to Refvecs_Fixed" << std::endl;

      return eRefvecs_Fixed; //default
    }

  template<class T>
    std::vector<SembleVector<double> > ST0FitPrim<T>::getEvals(void)
    {
      if(!!!solved)
      {
        initChk();
        solve();
        sort_solved();
      }

      std::vector<SembleVector<double> > dum;

      SembleVector<double> one(B, M);
      one.ones();

      for(int i = 0; i < inikeys.globalProps.tmin; ++i)
        dum.push_back(one);


      for(int t = 0; t <= inikeys.globalProps.tmax - inikeys.globalProps.tmin; ++t)
        dum.push_back(_data[t]->getEvals());

      return dum;
    }

  template<class T>
    typename std::vector<SembleMatrix<T> > ST0FitPrim<T>::getEvecs(void)
    {
      if(!!!solved)
      {
        initChk();
        solve();
        sort_solved();
      }

      typename std::vector<SembleMatrix<T> > dum;

      SembleMatrix<T> one(B,N, M);
      one.ones();

      for(int i = 0; i < inikeys.globalProps.tmin; ++i)
        dum.push_back(one);

      for(int t = 0; t <= inikeys.globalProps.tmax - inikeys.globalProps.tmin; ++t)
        dum.push_back(_data[t]->getEvecs());

      return dum;
    }


}
#endif
