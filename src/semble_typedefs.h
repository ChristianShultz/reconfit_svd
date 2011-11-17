#ifndef SEMBLE_TYPEDEFS_H_H_GUARD
#define SEMBLE_TYPEDEFS_H_H_GUARD

#include "semble_matrix.h"
#include "semble_vector.h"
#include <complex>

namespace SEMBLE
{

  //typedef the Semble templates so that they can slot into legacy code that uses the EnsemMatrix and EnsemVec classes
  //the interface is the same so we can just use the typedef to save a bit of work

  typedef SembleMatrix<double> EnsemMatrixReal;
  typedef SembleMatrix<std::complex<double> > EnsemMatrixComplex;
  typedef SembleVector<double> EnsemVecReal;
  typedef SembleVector<std::complex<double> > EnsemVecComplex;


  //for people who don't like template arguments
  typedef SembleMatrix<double> SembleMatrixReal;
  typedef SembleMatrix<std::complex<double> > SembleMatrixComplex;
  typedef SembleVector<double> SembleVectorReal;
  typedef SembleVector<std::complex<double> > SembleVectorComplex;

}

#endif
