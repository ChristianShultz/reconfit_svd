#ifndef SEMBLE_META_H_H_GUARD
#define SEMBLE_META_H_H_GUARD

#include <complex>
#include"ensem/ensem.h"

using namespace ENSEM;

namespace SEMBLE
{
  inline
  Real toScalar(const double val)
  {
    return Real(val);
  }

  inline
  Complex toScalar(const std::complex<double> val)
  {
    return cmplx(Real(val.real()), Real(val.imag()));
  }

  inline
  double toScalar(const Real val)
  {
    return toDouble(val);
  }

  //so I needed to include these to get it to compile, something about the innerworkings
  //of the EnsemReal, EnsemVec classes force me to need to also provide definitions for these
  //which I find to be a bit strange..
  inline
  double toScalar(const OScalar<PScalar<PScalar<RScalar<REAL> > > > val)
  {
    return toDouble(val);
  }

  inline
  std::complex<double> toScalar(const OScalar<PScalar<PScalar<RComplex<REAL> > > > val)
  {
    return (std::complex<double>(toDouble(real(val)), toDouble(imag(val))));
  }

  inline
  std::complex<double> toScalar(Complex val)
  {
    return (std::complex<double>(toDouble(real(val)), toDouble(imag(val))));
  }

//////////////////////////////////////////////////////////////////////////
  template<class T1>
  struct PromoteScalar
  {
    typedef T1 Type;
  };

  template<>
  struct PromoteScalar<double>
  {
    typedef Real Type;
  };

  template<>
  struct PromoteScalar<std::complex<double> >
  {
    typedef Complex Type;
  };


//////////////////////////////////////////////////////////////////////////
  template<class T1>
  struct PromoteEnsem
  {
    typedef T1 Type;
  };

  template<>
  struct PromoteEnsem<double>
  {
    typedef EnsemReal Type;
  };

  template<>
  struct PromoteEnsem<std::complex<double> >
  {
    typedef EnsemComplex Type;
  };


/////////////////////////////////////////////////////////////////////////
  template<class T1>
  struct PromoteEnsemVec
  {
    typedef T1 Type;
  };

  template<>
  struct PromoteEnsemVec<double>
  {
    typedef EnsemVectorReal Type;
  };

  template<>
  struct PromoteEnsemVec<std::complex<double> >
  {
    typedef EnsemVectorComplex Type;
  };


/////////////////////////////////////////////////////////////////////////
  template<class T1, class T2>
  struct PromoteStl
  {
    typedef T1 Type;
  };

  template<>
  struct PromoteStl<double, double>
  {
    typedef double Type;
  };

  template<>
  struct PromoteStl<double, std::complex<double> >
  {
    typedef std::complex<double> Type;
  };

  template<>
  struct PromoteStl<std::complex<double>, double>
  {
    typedef std::complex<double>  Type;
  };

  template<>
  struct PromoteStl<std::complex<double>, std::complex<double> >
  {
    typedef std::complex<double> Type;
  };

}
#endif
