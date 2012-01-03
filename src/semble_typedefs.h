#ifndef SEMBLE_TYPEDEFS_H_H_GUARD
#define SEMBLE_TYPEDEFS_H_H_GUARD

#include "semble_matrix.h"
#include "semble_vector.h"
#include <complex>



  //typedef the Semble templates so that they can slot into legacy code that uses the EnsemMatrix and EnsemVec classes
  //the interface is the same so we can just use the typedef to save a bit of work

typedef SEMBLE::SembleMatrix<double> EnsemMatrixReal;
typedef SEMBLE::SembleMatrix<std::complex<double> > EnsemMatrixComplex;
typedef SEMBLE::SembleVector<double> EnsemVecReal;
typedef SEMBLE::SembleVector<std::complex<double> > EnsemVecComplex;

namespace SEMBLE
{
  //for people who don't like template arguments
  typedef SembleMatrix<double> SembleMatrixReal;
  typedef SembleMatrix<std::complex<double> > SembleMatrixComplex;
  typedef SembleVector<double> SembleVectorReal;
  typedef SembleVector<std::complex<double> > SembleVectorComplex;
}

#endif
