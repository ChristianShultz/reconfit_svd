#ifndef SEMBLE_MATRIX_H_H_GUARD
#define SEMBLE_MATRIX_H_H_GUARD

#include "semble_meta.h"
#include "semble_vector.h"
#include "itpp/itbase.h"
#include <iostream>

using namespace ENSEM;

namespace SEMBLE
{
//forward declaration
  template<class T>
  struct SembleMatrix;

//forward dec of unary minus friend
  template<class T>
  SembleMatrix<T> operator-(const SembleMatrix<T> &plus);

//no asymmetric operators are supported b/c they arent supported int itpp

  template <class T>
  struct SembleMatrix
  {

  public:
    //constructors,destructors,assignment
    SembleMatrix(void);
    SembleMatrix(const int bins_, const int rows_, const int cols_);
    SembleMatrix(const SembleMatrix<T> &other);
    SembleMatrix(const typename std::vector<itpp::Mat<T> > &in);
    ~SembleMatrix() {}
    SembleMatrix<T>& operator=(const SembleMatrix<T> &other);
    SembleMatrix<T>& operator=(const itpp::Mat<T> &other);

  public:
    //peek data
    itpp::Mat<T> operator[](const int bin_) const;
    typename PromoteEnsem<T>::Type operator()(const int row_, const int col_) const
    {
      return getEnsemElement(row_, col_);
    }
    typename PromoteEnsem<T>::Type getEnsemElement(const int row_, const int col_) const;
    typename PromoteScalar<T>::Type operator()(const int bin_, const int row_, const int col_) const;
    SembleVector<T> getRow(const int row) const;
    SembleVector<T> getCol(const int col) const;

  public:
    //poke data--requires modification
    itpp::Mat<T>& operator[](const int bin);
    void loadEnsemElement(const int row_, const int col_, const typename PromoteEnsem<T>::Type &ensem);
    void setElement(const int bin_, const int row_, const int col_, const T elem_);
    void setElement(const int bin_, const int row_, const int col_, const typename PromoteScalar<T>::Type &scalar);
    void zeros(void);
    void ones(void);

  public:
    //dimension functions--require modification
    int getN(void) const
    {
      return N;
    }
    int getM(void) const
    {
      return M;
    }
    int getB(void) const
    {
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
    void setN(const int n);
    void setM(const int m);
    void setB(const int b);
    void del_rows(const int r1, const int r2);
    void del_cols(const int c1, const int c2);
    inline void rows(const int r)
    {
      setN(r);
    }
    inline void cols(const int c)
    {
      setM(c);
    }
    inline void bins(const int b)
    {
      setB(b);
    }
    inline void reDim(const int b, const int c)
    {
      setN(c);
      setM(c);
      setB(b);
    }
    inline void reDim(const int b, const int r, const int c)
    {
      setB(b);
      setN(r);
      setM(c);
    }

  public:
    //statistics
    itpp::Mat<T> mean(void) const;
    itpp::Mat<T> variance(void) const;

  public:
    //rescaling--requires modification
    void rescaleSembleDown(void);
    void rescaleSembleUp(void);
    inline void rescaleEnsemDown(void)
    {
      rescaleSembleDown();
    }
    inline void rescaleEnsemUp(void)
    {
      rescaleSembleUp();
    }

  public:
    //algebra--requires modification
    SembleMatrix<T>& operator+=(const SembleMatrix<T> &rhs);
    SembleMatrix<T>& operator+=(const itpp::Mat<T> &rhs);
    SembleMatrix<T>& operator-=(const SembleMatrix<T> &rhs);
    SembleMatrix<T>& operator-=(const itpp::Mat<T> &rhs);
    SembleMatrix<T>& operator*=(const T &rhs);
    SembleMatrix<T>& operator*=(const typename PromoteScalar<T>::Type &rhs);
    SembleMatrix<T>& operator*=(const typename PromoteEnsem<T>::Type &rhs);   //requires rescaling
    SembleMatrix<T>& operator*=(const SembleMatrix<T> &rhs);
    SembleMatrix<T>& operator*=(const itpp::Mat<T> &rhs);
    SembleMatrix<T>& operator/=(const T &rhs);
    SembleMatrix<T>& operator/=(const typename PromoteScalar<T>::Type &rhs);
    SembleMatrix<T>& operator/=(const typename PromoteEnsem<T>::Type &rhs);   //requires rescaling

    //unary minus--requires modification
    friend SembleMatrix<T> operator-<>(const SembleMatrix<T> &plus);

    //hermitian conjugate--requires modification
    SembleMatrix<T>& hermitianConjugate(void);
    inline SembleMatrix<T>& H(void)
    {
      return hermitianConjugate();
    }
    inline SembleMatrix<T>& adj(void)
    {
      return hermitianConjugate();
    }

    //transpose--requires modification
    SembleMatrix<T>& transpose(void);
    inline SembleMatrix<T>& tran(void)
    {
      return transpose();
    }

  private:         //data
    int N, M, B;
    typename std::vector<itpp::Mat<T> > semble;

  };

//Implementation
////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<class T>
  SembleMatrix<T>::SembleMatrix(void) : N(0), M(0), B(0)
  {
  }


//reserve constructor
  template<class T>
  SembleMatrix<T>::SembleMatrix(const int b, const int r, const int c) : B(b), N(r), M(c)
  {
    itpp::Mat<T> zero(N, M);
    zero.zeros();
    semble = std::vector<itpp::Mat<T> >(B, zero);
  }

  template<class T>
  SembleMatrix<T>::SembleMatrix(const SembleMatrix<T> &other)
  {
    N = other.N;
    M = other.M;
    B = other.B;
    semble = other.semble;
  }

//create from a vector of itpp mats
  template<class T>
  SembleMatrix<T>::SembleMatrix(const typename std::vector<itpp::Mat<T> > &in)
  {

    semble = typename std::vector<itpp::Mat<T> >(in);

    B = semble.size();

    if(B != 0)
      {
        N = semble[0].rows();
        M = semble[0].cols();
      }
    else
      {
        N = 0;
        M = 0;
      }

  }

  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator=(const SembleMatrix<T> &other)
  {
    if(this != &other)
      {
        N = other.N;
        M = other.M;
        B = other.B;
        semble = other.semble;
      }

    return *this;
  }

  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator=(const itpp::Mat<T> &other)
  {

    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it = other;

    N = other.rows();
    M = other.cols();

    return *this;
  }


//const methods -
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//peek

//get a mat on a bin
  template<class T>
  itpp::Mat<T> SembleMatrix<T>::operator[](const int bin_) const
  {
    if(bin_ >= B)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    return semble.at(bin_);
  }

//get the ensemble on the element
  template<class T>
  typename PromoteEnsem<T>::Type SembleMatrix<T>::getEnsemElement(const int row_, const int col_) const
  {
    if((row_ >= N) || (col_ >= M))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    typename PromoteEnsem<T>::Type dum;
    dum.resize(B);

    for(int bin = 0; bin < B; ++bin)
      pokeEnsem(dum, toScalar((semble.at(bin)).get(row_, col_)), bin);

    return dum;
  }

//get a number
  template<class T>
  typename PromoteScalar<T>::Type SembleMatrix<T>::operator()(const int bin_, const int row_, const int col_) const
  {
    if((bin_ >= B) || (row_ >= N) || (col_ >= M))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    return toScalar((semble.at(bin_)).get(row_, col_));
  }

  //get an ensemble vector from a row
  template<class T>
  SembleVector<T> SembleMatrix<T>::getRow(const int row) const
  {
    SembleVector<T> ret(B, M);

    for(int col = 0; col < M; ++col)
      ret.loadEnsemElement(col, this->getEnsemElement(row, col));

    return ret;
  }

  //get an ensemble vector from a col
  template<class T>
  SembleVector<T> SembleMatrix<T>::getCol(const int col) const
  {
    SembleVector<T> ret(B, N);

    for(int row = 0; row < N; ++row)
      ret.loadEnsemElement(row, this->getEnsemElement(row, col));

    return ret;
  }


//methods requiring
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//poke data

//ref to Mat on bin
  template<class T>
  itpp::Mat<T>& SembleMatrix<T>::operator[](const int bin_)
  {
    if(bin_ >= B)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    return semble.at(bin_); //technically at checks the bound for us but its nice to check
  }

//load an ensemble into the element
  template<class T>
  void SembleMatrix<T>::loadEnsemElement(const int row_, const int col_, const typename PromoteEnsem<T>::Type &ensem)
  {
    if((row_ >= N) || (col_ >= M) || (ensem.size() != B))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    for(int bin = 0; bin < B; ++bin)
      (semble.at(bin))(row_, col_) = toScalar(ensem.elem(bin));
  }

  template<class T>
  void SembleMatrix<T>::setElement(const int bin_, const int row_, const int col_, const T elem_)
  {
    if((bin_ >= B) || (row_ >= N) || (col_ >= N))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    semble.at(bin_).set(row_, col_, elem_);
  }

  template<class T>
  void SembleMatrix<T>::setElement(const int bin_, const int row_, const int col_, const typename PromoteScalar<T>::Type &scalar)
  {
    setElement(bin_, row_, col_, toScalar(scalar));
  }

  template<class T>
  void SembleMatrix<T>::zeros(void)
  {

    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      it->zeros();
  }

  template<class T>
  void SembleMatrix<T>::ones(void)
  {

    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      it->ones();
  }

//
//dimension functions

//change n rows
  template<class T>
  void SembleMatrix<T>::setN(const int n)
  {
    if(N != n)
      {

        typename std::vector<itpp::Mat<T> >::iterator it;

        for(it = semble.begin(); it != semble.end(); ++it)
          it->set_size(n, M, true);

        N = n;
      }
  }

//change n cols
  template<class T>
  void SembleMatrix<T>::setM(const int m)
  {
    if(M != m)
      {

        typename std::vector<itpp::Mat<T> >::iterator it;

        for(it = semble.begin(); it != semble.end(); ++it)
          it->set_size(N, m, true);

        M = m;
      }
  }

//change n bins
  template<class T>
  void SembleMatrix<T>::setB(const int b)
  {
    if(B != b)
      {

        semble.resize(b);
        B = b;
      }
  }

//delete rows r1 to r2
  template<class T>
  void SembleMatrix<T>::del_rows(const int r1, const int r2)
  {
    if(r1 > r2)
      {
        int tmp = r1;
        r1 = r2;
        r2 = tmp;
      }

    if(r2 >= N)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      it->del_rows(r1, r2);
  }

//delete cols c1 to c2
  template<class T>
  void SembleMatrix<T>::del_cols(const int c1, const int c2)
  {
    if(c1 > c2)
      {
        int tmp = c1;
        c1 = c2;
        c2 = tmp;
      }

    if(c2 >= N)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      it->del_cols(c1, c2);
  }

//const methods -- no modification
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//statistics

//mean
  template<class T>
  itpp::Mat<T> SembleMatrix<T>::mean(void) const
  {
    itpp::Mat<T> sum;
    sum.set_size(N, M);
    sum.zeros();

    typename std::vector<itpp::Mat<T> >::const_iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      sum += *it;

    return sum / double(B);
  }

//variance
  template<class T>
  itpp::Mat<T> SembleMatrix<T>::variance(void) const
  {
    itpp::Mat<T> var, ave = mean();
    var.set_size(N, M);
    var.zeros();

    typename std::vector<itpp::Mat<T> >::const_iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      var += itpp::elem_mult((*it - ave), (*it - ave));

    return var / double(B * (B - 1));
  }


//methods requiring modification
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//rescale ensembles
  template<class T>
  void SembleMatrix<T>::rescaleSembleDown(void)
  {

    itpp::Mat<T> average = mean();
    double scale = 1.0 / double(B - 1);
    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it = average - scale * (*it - average);
  }

  template<class T>
  void SembleMatrix<T>::rescaleSembleUp(void)
  {

    itpp::Mat<T> average = mean();
    double scale = double(B - 1);
    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it = average - scale * (*it - average);
  }

//
//algebra methods

//add semble to semble
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator+=(const SembleMatrix<T> &rhs)
  {
    if((rhs.getN() != N) || (rhs.getM() != M) || (rhs.getB() != B))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] += rhs[bin];

    return *this;
  }

//add itpp::Mat<T> to semble
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator+=(const itpp::Mat<T> &rhs)
  {
    if((rhs.rows() != N) || (rhs.cols() != M))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] += rhs;

    return *this;
  }

//subtract a semble from a semble
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator-=(const SembleMatrix<T> &rhs)
  {
    if((rhs.getN() != N) || (rhs.getM() != M) || (rhs.getB() != B))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] -= rhs[bin];

    return *this;
  }

//subtract itpp::Mat<T> to semble
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator-=(const itpp::Mat<T> &rhs)
  {
    if((rhs.rows() != N) || (rhs.cols() != M))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] -= rhs;

    return *this;
  }

//multiply by a scalar T
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator*=(const T &rhs)
  {


    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it *= rhs;

    return *this;
  }

//multiply by an ensem type scalar T
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator*=(const typename PromoteScalar<T>::Type &rhs)
  {


    T dum = toScalar(rhs);

    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it *= dum;

    return *this;
  }

//multiply by an ensemble of scalars -- requires rescaling
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator*=(const typename PromoteEnsem<T>::Type &ensem)
  {
    if(ensem.size() != B)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    rescaleSembleDown();
    typename PromoteEnsem<T>::Type rhsprime = ::rescaleEnsemDown(ensem);

    for(int bin = 0; bin < B; ++bin)
      semble[bin] *= toScalar(rhsprime.elem(bin));

    rescaleSembleUp();
    return *this;
  }

//multiply by a semble
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator*=(const SembleMatrix<T> &rhs)
  {

    if(this == &rhs)
      {
        if(N != M)
          {
            std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
            exit(1);
          }


        rescaleSembleDown();
        typename std::vector<itpp::Mat<T> >::iterator it;

        for(it = semble.begin(); it != semble.end(); ++it)
          *it *= *it;

        rescaleSembleUp();
        return *this;
      }

    if(rhs.getN() != M)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    if(rhs.getB() != B)
      {
        std::cout << "Bin mismatch in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    rescaleSembleDown();

    SembleMatrix<T> dum(rhs);
    dum.rescaleSembleDown();

    for(int bin = 0; bin < B; ++bin)
      semble[bin] *= dum[bin];

    M = rhs.getM();

    rescaleSembleUp();

    return *this;
  }

  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator*=(const itpp::Mat<T> &rhs)
  {
    rescaleSembleDown();

    typename std::vector<itpp::Mat<T> >::iterator it;
    for(it = semble.begin(); it != semble.end(); ++it)
      *it *= rhs;

    rescaleSembleUp();
    
    return *this;
  }


//divide by a scalar T
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator/=(const T &rhs)
  {


    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it /= rhs;

    return *this;
  }

//divide by an ensem type scalar T
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator/=(const typename PromoteScalar<T>::Type &rhs)
  {


    T dum = toScalar(rhs);

    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it /= dum;

    return *this;
  }

//divide by an ensemble of scalars -- requires rescaling
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::operator/=(const typename PromoteEnsem<T>::Type &ensem)
  {
    if(ensem.size() != B)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    rescaleSembleDown();
    typename PromoteEnsem<T>::Type rhsprime = ::rescaleEnsemDown(ensem);

    for(int bin = 0; bin < B; ++bin)
      semble[bin] /= toScalar(rhsprime.elem(bin));

    rescaleSembleUp();
    return *this;
  }

//conjugate transpose
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::hermitianConjugate(void)
  {


    int tmp = M;
    M = N;
    N = tmp;

    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it = itpp::hermitian_transpose(*it);

    return *this;
  }

//transpose
  template<class T>
  SembleMatrix<T>& SembleMatrix<T>::transpose(void)
  {

    int tmp = M;
    M = N;
    N = tmp;

    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it = itpp::transpose(*it);

    return *this;
  }

//friends, only one -- apparently SembleMatrix isnt very popular
  template<class T>
  SembleMatrix<T> operator-(const SembleMatrix<T> &plus)
  {
    SembleMatrix<T> minus(plus);
    typename std::vector<itpp::Mat<T> >::iterator it;

    for(it = minus.semble.begin(); it != minus.semble.end(); ++it)
      *it = -*it;

    return minus;
  }

}
#endif
