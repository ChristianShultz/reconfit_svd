#ifndef SEMBLE_VECTOR_H_H_GUARD
#define SEMBLE_VECTOR_H_H_GUARD

#include "semble_meta.h"
#include "itpp/itbase.h"
#include <iostream>

using namespace ENSEM;

namespace SEMBLE
{

//forward declaration of vector
  template<class T>
  struct SembleVector;


//unary minus friend
  template<class T>
  SembleVector<T> operator-(const SembleVector<T> &plus);

//no asymmetric opperations are supported -- ie cant add a semble<double> to semble<complex<double> >

  template <class T>
  struct SembleVector
  {

  public:
    //constructors,destructor,assignment
    SembleVector(void);
    SembleVector(int bins_, int num_elem);
    SembleVector(const SembleVector<T> &other);
    SembleVector(const typename std::vector<itpp::Vec<T> > &in);
    virtual ~SembleVector();
    SembleVector<T>& operator=(const SembleVector<T> &other);
    SembleVector<T>& operator=(const typename PromoteEnsemVec<T>::Type &other);
    SembleVector<T>& operator=(const itpp::Vec<T> &other);

  private:
    //lazy copy -> deep copy
    void modify(void);

  public:
    //peek data
    itpp::Vec<T> operator[](const int bin) const;                                                 //return the vector on a bin
    typename PromoteEnsem<T>::Type operator()(const int elem) const
    {
      return getEnsemElement(elem); //return the ensemble of an elem
    }
    typename PromoteEnsem<T>::Type getEnsemElement(const int elem) const;                         //return the ensemble on an elem
    typename PromoteScalar<T>::Type operator()(const int bin, const int elem) const;              //return a number

  public:
    //poke data--requires modification
    itpp::Vec<T>& operator[](const int bin);
    void loadEnsemElement(int elem, const typename PromoteEnsem<T>::Type &ensem);
    void setElement(const int bin_, const int row_, const T elem_);
    void setElement(const int bin_, const int row_, const typename PromoteScalar<T>::Type &scalar);
    void zeros(void);
    void ones(void);

  public:
    //dimension functions--require modification
    int getN(void) const
    {
      return N;
    }
    int getB(void) const
    {
      return B;
    }
    inline int bins(void) const
    {
      return getB();
    }
    inline int rows(void) const
    {
      return getN();
    }
    inline int elem(void) const
    {
      return getN();
    }
    void setN(const int n);                                             //drop the higher elements by index
    void setB(const int b);                                             //drop the higher bins by index
    void del(const int e1, const int e2);                               //drop elements e1 to e2
    inline void rows(const int r)
    {
      setN(r);
    }
    inline void bins(const int b)
    {
      setB(b);
    }
    inline void reDim(const int b, const int r)
    {
      setB(b);
      setN(r);
    }

  public:
    //statistics
    itpp::Vec<T> mean(void) const;
    itpp::Vec<T> variance(void) const;

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
    SembleVector<T>& operator+=(const SembleVector<T> &rhs);
    SembleVector<T>& operator+=(const itpp::Vec<T> &rhs);
    SembleVector<T>& operator-=(const SembleVector<T> &rhs);
    SembleVector<T>& operator-=(const itpp::Vec<T> &rhs);
    SembleVector<T>& operator*=(const T &rhs);
    SembleVector<T>& operator*=(const typename PromoteScalar<T>::Type &scalar);
    SembleVector<T>& operator*=(const typename PromoteEnsem<T>::Type &ensem);    //requires rescaling
    SembleVector<T>& operator/=(const T &rhs);
    SembleVector<T>& operator/=(const typename PromoteScalar<T>::Type &scalar);
    SembleVector<T>& operator/=(const typename PromoteEnsem<T>::Type &ensem);    //requires rescaling

    //unary minus
    friend SembleVector<T> operator-<>(const SembleVector<T> &plus);

  public:
    //dot product
    typename PromoteEnsem<T>::Type dot(const SembleVector<T> &other);


  private:

  private:      //data
    int N, B;
    typename std::vector<itpp::Vec<T> > semble; //shared data for lazy copies
  };


//Implementation
/////////////////////////////////////////////////////////////

//default constructor--put a single zero in it, make it one bin
  template<class T>
  SembleVector<T>::SembleVector(void) : N(0), B(0)
  {
  }

//reserve constructor
  template<class T>
  SembleVector<T>::SembleVector(int b, int num_elem) : B(b) , N(num_elem)
  {

    itpp::Vec<T> zero(N);
    zero.zeros();
    semble = typename std::vector<itpp::Vec<T> >(B, zero);

  }

  template<class T>
  SembleVector<T>::SembleVector(const SembleVector<T> &other) : N(other.N) , B(other.B), semble(other.semble)
  {
  }

//create from vector of itpp vecs
  template<class T>
  SembleVector<T>::SembleVector(const typename std::vector<itpp::Vec<T> > &in)
  {

    semble = typename std::vector<itpp::Vec<T> >(in);

    B = semble.size();

    if(B != 0)
      N = semble[0].length();
    else
      N = 0;
  }

//cleanup but keep track of refs
  template<class T>
  SembleVector<T>::~SembleVector(void)
  {
  }

//assignment operator
  template<class T>
  SembleVector<T>& SembleVector<T>::operator=(const SembleVector<T> &other)
  {
    if(this != &other)                                              //self assignment
      {
        N = other.N;
        B = other.B;
        semble = other.semble;
      }

    return *this;
  }

  //convert from EnsemVector
  template<class T>
  SembleVector<T>& SembleVector<T>::operator=(const typename PromoteEnsemVec<T>::Type &other)
  {
    B = other.size();
    N = other.numElem();
    semble.resize(B);
    itpp::Vec<T> v(N);

    for(int bin_ = 0; bin_ < B; ++bin_)
      {
	for(int elem_ = 0; elem_ < N; ++elem_)
	  v[elem_] = toScalar(peekObs(other,elem_).elem(bin_));
	semble[bin_] = v;
      }
    return *this;
  }

//copy the vector to every bin
  template<class T>
  SembleVector<T>& SembleVector<T>::operator=(const itpp::Vec<T> &other)
  {
    semble.resize(B, other);
    N = other.size();
    return *this;
  }

//const methods
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//peek

//get the vector on the bin
  template<class T>
  itpp::Vec<T> SembleVector<T>::operator[](const int bin) const
  {
    if(bin >= B)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    return semble.at(bin); //technically stl at(int) checks bounds but the err message is nice and we dont care about speed
  }

//get the ensemble element
  template<class T>
  typename PromoteEnsem<T>::Type SembleVector<T>::getEnsemElement(const int elem) const
  {
    if(elem >= N)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    typename PromoteEnsem<T>::Type dum;
    dum.resize(B);

    for(int bin = 0; bin < B; ++bin)
      pokeEnsem(dum, toScalar((semble.at(bin)).get(elem)), bin); //technically stl at(int) checks bounds but the err message is nice and we dont care about speed

    return dum;
  }

//get a number
  template<class T>
  typename PromoteScalar<T>::Type SembleVector<T>::operator()(const int bin, const int elem) const
  {
    if((elem >= N) || (bin >= B))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    return toScalar((semble.at(bin)).get(elem));  //technically stl at(int) checks bounds but the err message is nice and we dont care about speed
  }

//non const methods
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//poke data

//ref to vec on bin
  template<class T>
  itpp::Vec<T>& SembleVector<T>::operator[](const int bin)
  {
    if(bin >= B)
      {
        std::cout << "Out of bounds in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    return semble.at(bin); //technically at checks the bound for us but its nice to check
  }

//load an ensemble into the element
  template<class T>
  void SembleVector<T>::loadEnsemElement(const int elem, const typename PromoteEnsem<T>::Type &ensem)
  {
    if((ensem.size() != B) || (elem >= N))
      {
        std::cout << "Out of Bounds in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      (semble.at(bin)).set(elem, toDouble(ensem.elem(bin)));
  }

  template<class T> //put something in the specifed slot
  void SembleVector<T>::setElement(const int bin_, const int row_, const T elem_)
  {
    if((bin_ >= B) || (row_ >= N))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    semble.at(bin_).set(row_, elem_);
  }

  template<class T>
  void SembleVector<T>::setElement(const int bin_, const int row_, const typename PromoteScalar<T>::Type &scalar)
  {
    setElement(bin_, row_, toScalar(scalar));
  }

  template<class T>
  void SembleVector<T>::zeros(void)
  {

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      it->zeros();
  }

  template<class T>
  void SembleVector<T>::ones(void)
  {

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      it->ones();
  }

//
//dimension functions

//change num elems
  template<class T>
  void SembleVector<T>::setN(const int n)
  {
    if(N != n)
      {

        typename std::vector<itpp::Vec<T> >::iterator it;

        for(it = semble.begin(); it != semble.end(); ++it)
          it->set_size(n, true);

        N = n;
      }
  }

//change num bins, delete the higher index guys if b < B, copy the last element if b > B
  template<class T>
  void SembleVector<T>::setB(const int b)
  {
    if(B != b)
      {

        semble.resize(b);
        B = b;
      }
  }

//delete elements e1 to e2
  template<class T>
  void SembleVector<T>::del(const int e1, const int e2)
  {
    if(e2 >= N)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      it->del(e1, e2);

    N -= e2 - e1 + 1;
  }

//const methods -- no modification
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//statistics

//mean
  template<class T>
  itpp::Vec<T> SembleVector<T>::mean(void) const
  {
    itpp::Vec<T> sum;
    sum.set_size(N);
    sum.zeros();

    typename std::vector<itpp::Vec<T> >::const_iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      sum += *it;

    return sum / double(B);
  }

//variance--elementwise
  template<class T>
  itpp::Vec<T> SembleVector<T>::variance(void) const
  {
    itpp::Vec<T> var, ave = mean();
    var.set_size(N);
    var.zeros();

    typename std::vector<itpp::Vec<T> >::const_iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      var += itpp::elem_mult((*it - ave), (*it - ave));

    return var / (double(B * (B - 1)));
  }

//methods requiring modification
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//rescale ensembles
  template<class T>
  void SembleVector<T>::rescaleSembleDown(void)
  {

    itpp::Vec<T> average = mean();
    double scale = 1.0 / double(B - 1);

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it = average - scale * (*it - average);
  }

  template<class T>
  void SembleVector<T>::rescaleSembleUp(void)
  {

    itpp::Vec<T> average = mean();
    double scale = double(B - 1);

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it = average - scale * (*it - average);
  }

//
//algebra functions

//add two sembles
  template<class T>
  SembleVector<T>& SembleVector<T>::operator+=(const SembleVector<T> &rhs)
  {
    if((rhs.getN() != N) || (rhs.getB() != B))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] += rhs[bin];

    return *this;
  }

//add a semble and an itpp
  template<class T>
  SembleVector<T>& SembleVector<T>::operator+=(const itpp::Vec<T> &rhs)
  {
    if(rhs.size() != N)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] += rhs;

    return *this;
  }

//subtract semble
  template<class T>
  SembleVector<T>& SembleVector<T>::operator-=(const SembleVector<T> &rhs)
  {
    if((rhs.getN() != N) || (rhs.getB() != B))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] -= rhs[bin];

    return *this;
  }

//subtract a semble and an itpp
  template<class T>
  SembleVector<T>& SembleVector<T>::operator-=(const itpp::Vec<T> &rhs)
  {
    if(rhs.size() != N)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }



    for(int bin = 0; bin < B; ++bin)
      semble[bin] -= rhs;

    return *this;
  }

//multiply by a T scalar
  template<class T>
  SembleVector<T>& SembleVector<T>::operator*=(const T &rhs)
  {

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it *= rhs;

    return *this;
  }

//multiply by an Ensem type scalar
  template<class T>
  SembleVector<T>& SembleVector<T>::operator*=(const typename PromoteScalar<T>::Type &rhs)
  {

    T dum = toScalar(rhs);

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it *= dum;

    return *this;
  }

//multiply by an ensemble of scalars --rescaling operation
  template<class T>
  SembleVector<T>& SembleVector<T>::operator*=(const typename PromoteEnsem<T>::Type &rhs)
  {
    if(rhs.size() != B)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    rescaleSembleDown();

    typename PromoteEnsem<T>::Type rhsprime = ::rescaleEnsemDown(rhs);

    for(int bin = 0; bin < B; ++bin)
      semble[bin] *= toScalar(rhsprime.elem(bin));

    rescaleSembleUp();

    return *this;
  }

//divide by a scalar
  template<class T>
  SembleVector<T>& SembleVector<T>::operator/=(const T &rhs)
  {


    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it /= rhs;

    return *this;
  }

//divide by an ensem type scalar
  template<class T>
  SembleVector<T>& SembleVector<T>::operator/=(const typename PromoteScalar<T>::Type &rhs)
  {

    T dum = toScalar(rhs);

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = semble.begin(); it != semble.end(); ++it)
      *it /= dum;

    return *this;
  }

//divide by an ensemble of scalars
  template<class T>
  SembleVector<T>& SembleVector<T>::operator/=(const typename PromoteEnsem<T>::Type &rhs)
  {
    if(rhs.size() != B)
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    rescaleSembleDown();

    typename PromoteEnsem<T>::Type rhsprime = ::rescaleEnsemDown(rhs);

    for(int bin = 0; bin < B; ++bin)
      semble[bin] /= toScalar(rhsprime.elem(bin));

    rescaleSembleUp();

    return *this;
  }

  template<class T>
  typename PromoteEnsem<T>::Type SembleVector<T>::dot(const SembleVector<T> &V)
  {

    if(this == &V)
      {

        typename PromoteEnsem<T>::Type dum2;
        dum2.resize(B);
        rescaleSembleDown();

        for(int bin = 0; bin < B; ++bin)
          pokeEnsem(dum2, toScalar(itpp::dot(semble[bin], semble[bin])), bin);

        rescaleSembleUp();
        return dum2;
      }

    if((V.bins() != B) || (V.rows() != N))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    rescaleSembleDown();

    SembleVector<T> dum(V);
    dum.rescaleSembleDown();

    typename PromoteEnsem<T>::Type dum2;
    dum2.resize(B);

    for(int bin = 0; bin < B; ++bin)
      pokeEnsem(dum2, toScalar(itpp::dot(semble[bin], dum[bin])), bin);

    rescaleSembleUp();

    return ::rescaleEnsemUp(dum2);
  }

//friends
  template<class T>
  SembleVector<T> operator-(const SembleVector<T> &plus)
  {
    SembleVector<T> minus(plus);

    typename std::vector<itpp::Vec<T> >::iterator it;

    for(it = minus.semble.begin(); it != minus.semble.end(); ++it)
      *it = -*it;

    return minus;
  }


}

#endif
