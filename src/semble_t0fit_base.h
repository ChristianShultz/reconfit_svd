#ifndef SEMBLE_T0FIT_BASE_H_H_GUARD
#define SEMBLE_T0FIT_BASE_H_H_GUARD

#include"semble_matrix.h"
#include"semble_vector.h"
#include"semble_linear_algebra.h"
#include<string>
#include<iostream>

//this header contains template classes that encapsulate the solution to the generalized eigen-problem,
//they derive from a base so that the primitive ST0Fit class doesnt need to care how the problem was solved,
//all derived class data members should be public so that they can be accessed via pointers in ST0Fit_prim

namespace SEMBLE
{
  template<class T>
  struct ST0Base   //this is the abstract base class, its never actually used but I'm writing a full definition so
  {                //I can cut and paste, this and the derived class are just pieces of ST0Fit
  public:
    ST0Base(void)
      : init(false) , solved(false) , sort(true) , t(-1) {}
    ST0Base(const ST0Base &o)
    {
      init = o.init;
      solved = o.solved;
      sort = o.sort;
      t = o.t;

      if(o.init)
        {
          Ct = o.Ct;
          Ct0 = o.Ct0;
          F = o.F;
          RSP = o.RSP;
          U = o.U;

          if(o.solved)
            {
              V = o.V;
              W = o.W;
              eVals = o.eVals;
            }
        }
    }
    virtual ~ST0Base(void) {}
    virtual ST0Base<T>& operator=(const ST0Base<T> &o)
    {
      if(this != &o)
        {
          init = o.init;
          solved = o.solved;
          sort = o.sort;
          t = o.t;

          if(o.init)
            {
              Ct = o.Ct;
              Ct0 = o.Ct0;
              F = o.F;
              RSP = o.RSP;
              U = o.U;

              if(o.solved)
                {
                  V = o.V;
                  W = o.W;
                  eVals = o.eVals;
                }
            }
        }

      return *this;
    }
    virtual SembleVector<double> getEvals(void) const = 0;
    virtual SembleVector<double>& evals(void) = 0;
    virtual const SembleVector<double>& evals(void) const = 0;
    virtual SembleMatrix<T> getEvecs(void) const = 0;
    virtual SembleMatrix<T>& evecs(void) = 0;
    virtual const SembleMatrix<T>& evecs(void) const = 0;
    virtual SembleMatrix<T> getW(void) const = 0;
    virtual SembleMatrix<T>& w(void) = 0;
    virtual int getT(void) const = 0;
    virtual std::string echo(void) const = 0;
    virtual void eval(void) = 0;
    virtual ST0Base<T>* clone(void) const = 0;

  public:                                       //data members need to be public b/c this is effectively a fancy container class
    SembleMatrix<T> Ct, Ct0, W, V, F, U;        //possible ensembles
    SembleVector<double> eVals;                 //doing it this way will make t0fitting look cleaner and easier to implement
    SembleMatrix<double> RSP;
    int t;
    bool sort, init, solved;
  };


//a factorized cholesky base class
  template<class T>
  struct ChoF : public ST0Base<T>
  {
  public:
    ChoF(void)
      : init(false), solved(false) , sort(true) , t(-1) {}
    ChoF(const ChoF<T> &o)
    {
      init = o.init;
      solved = o.solved;
      sort = o.sort;
      t = o.t;

      if(o.init)
        {
          Ct = o.Ct;
          F = o.F;

          if(o.solved)
            {
              W = o.W;
              V = o.V;
              eVals = o.eVals;
            }
        }
    }
    ChoF(const SembleMatrix<T> &Ct_, const SembleMatrix<T> &F_, int t_, bool s = true)
      : Ct(Ct_) , F(F_) , t(t_) , init(true) , sort(s) {}
    ~ChoF(void) {}
    ChoF<T>& operator=(const ChoF<T> &o)
    {
      if(this != &o)
        {
          init = o.init;
          solved = o.solved;
          sort = o.sort;
          t = o.t;

          if(o.init)
            {
              Ct = o.Ct;
              F = o.F;

              if(o.solved)
                {
                  W = o.W;
                  V = o.V;
                  eVals = o.eVals;
                }
            }
        }

      return *this;
    }
    SembleVector<double> getEvals(void) const
    {
      if(solved)
        return eVals;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleVector<double>& evals(void)
    {
      if(solved)
        return eVals;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    const SembleVector<double>& evals(void) const
    {
      if(solved)
        return eVals;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }

    SembleMatrix<T> getEvecs(void) const
    {
      if(solved)
        return V;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleMatrix<T>& evecs(void)
    {
      if(solved)
        return V;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    const SembleMatrix<T>& evecs(void) const
    {
      if(solved)
        return V;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleMatrix<T> getW(void) const
    {
      if(solved)
        return W;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleMatrix<T>& w(void)
    {
      if(solved)
        return W;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    int getT(void) const
    {
      if(init)
        return t;

      std::cout << "Need to initialize the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    std::string echo(void) const
    {
      return std::string("ChoF");
    }
    void eval(void)
    {
      if(!!!init)
        {
          std::cout << "Need to initialize the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
          exit(1);
        }

      solved = true;
      genEigChoF(Ct, F, V, W, eVals, sort);
    }
    ChoF<T>* clone(void) const
    {
      return new ChoF<T>(*this);
    }

  public:
    SembleMatrix<T> Ct, W, V, F;
    SembleVector<double> eVals;
    int t;
    bool sort, init, solved;
  };


//a factorized svd base class
  template<class T>
  struct SvdF : public ST0Base<T>
  {
  public:
    SvdF(void)
      : init(false) , solved(false) , sort(true) , t(-1) {}
    SvdF(const SvdF<T> &o)
    {
      init = o.init;
      solved = o.solved;
      sort = o.sort;
      t = o.t;

      if(o.init)
        {
          Ct = o.Ct;
          RSP = o.RSP;
          U = o.U;

          if(o.solved)
            {
              W = o.W;
              V = o.V;
              eVals = o.eVals;
            }
        }
    }
    SvdF(const SembleMatrix<T> &Ct_, const SembleMatrix<double> &RSP_, const SembleMatrix<T> &U_, int t_, bool s = true)
      : Ct(Ct_) , RSP(RSP_) , U(U_) , init(true) , t(t_) , sort(s) {}
    ~SvdF(void) {}
    SvdF<T>& operator=(const SvdF<T> &o)
    {
      if(this != &o)
        {
          init = o.init;
          solved = o.solved;
          sort = o.sort;
          t = o.t;

          if(o.init)
            {
              Ct = o.Ct;
              RSP = o.RSP;
              U = o.U;

              if(o.solved)
                {
                  W = o.W;
                  V = o.V;
                  eVals = o.eVals;
                }
            }
        }

      return *this;
    }
    SembleVector<double> getEvals(void) const
    {
      if(solved)
        return eVals;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleVector<double>& evals(void)
    {
      if(solved)
        return eVals;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    const SembleVector<double>& evals(void) const
    {
      if(solved)
        return eVals;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleMatrix<T> getEvecs(void) const
    {
      if(solved)
        return V;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleMatrix<T>& evecs(void)
    {
      if(solved)
        return V;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    const SembleMatrix<T>& evecs(void) const
    {
      if(solved)
        return V;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleMatrix<T> getW(void) const
    {
      if(solved)
        return W;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    SembleMatrix<T>& w(void)
    {
      if(solved)
        return W;

      std::cout << "Need to solve the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    int getT(void) const
    {
      if(init)
        return t;

      std::cout << "Need to initialize the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
      exit(1);
    }
    std::string echo(void) const
    {
      return std::string("SvdF");
    }
    void eval(void)
    {
      if(!!!init)
        {
          std::cout << "Need to initialize the problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
          exit(1);
        }

      solved = true;
      genEigSvdF(Ct, RSP, U, V, W, eVals, sort);
    }
    SvdF<T>* clone(void) const
    {
      return new SvdF<T>(*this);
    }

  public:
    SembleMatrix<T> Ct, W, V, U;
    SembleMatrix<double> RSP;
    SembleVector<double> eVals;
    int t;
    bool sort, init, solved;
  };

}
#endif
