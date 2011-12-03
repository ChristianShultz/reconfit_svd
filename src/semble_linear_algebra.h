#ifndef SEMBLE_LINEAR_ALGEBRA_H_H_GUARD
#define SEMBLE_LINEAR_ALGEBRA_H_H_GUARD

#include"semble_meta.h"
#include"semble_vector.h"
#include"semble_matrix.h"
#include"semble_algebra.h"
#include<math.h>
#include<stdlib.h>
#include<map>
#include<utility>


namespace SEMBLE
{
//forward declarations
  template<class T>  //return T  if real, T* if complex, calls overloaded eConjPhase(T)
  typename PromoteScalar<T>::Type evecConjPhase(const T &in);
 
  //return toScalar(-T)
  PromoteScalar<double>::Type eConjPhase(const double &in);
  
  //return toScalar(T*)
  PromoteScalar<std::complex<double> >::Type eConjPhase(const std::complex<double> &in);

  template<class T>//solves R*P = V, P = R^dagger * V, will include the evecConjPhase in the return matrix
  itpp::Mat<T> constructPermutationMatrix(const itpp::Mat<T> &ref, const itpp::Mat<T> &vec);

//row/col operations
  template<class T>
  SembleVector<T> getRow(const SembleMatrix<T> &in, const int row);

  template<class T>
  SembleVector<T> getCol(const SembleMatrix<T> &in, const int col);

//semble linear algebra
  template<class T>//hermitian conjugate, transpose for reals
  SembleMatrix<T> hermitianConjugate(const SembleMatrix<T> &M);

  template<class T>//same
  SembleMatrix<T> adj(const SembleMatrix<T> &M);

  template<class T>//same
  SembleMatrix<T> adjoint(const SembleMatrix<T> &M);

  template<class T>//transpose
  SembleMatrix<T> transpose(const SembleMatrix<T> &M);

  template<class T>//same
  SembleMatrix<T> tran(const SembleMatrix<T> &M);

  template<class T> //loads into first row
  SembleMatrix<T> toMatrix(const SembleVector<T> &V);

  template<class T> //loads the corresponding row or col
  SembleVector<T> toVector(const SembleMatrix<T> &M);

  template<class T> //assumes a column vector
  SembleMatrix<T> adjoint(const SembleVector<T> &V);

  template<class T> //assumes a column vector
  SembleMatrix<T> adj(const SembleVector<T> &V);

  template<class T> //assumes a column vector
  SembleMatrix<T> transpose(const SembleVector<T> &V);

  template<class T> //assumes a column vector
  SembleMatrix<T> tran(const SembleVector<T> &V);

  template<class T> //symmetrize/hermitize? a matrix
  void symmetrize(const SembleMatrix<T> &in, SembleMatrix<T> &out);

  template<class T> //symmetrize/hermitize(sp?) a matrix
  SembleMatrix<T> symmetrize(const SembleMatrix<T> &in);

  template<class T> //create a diagonal matrix from a list of values
  void diag(const SembleVector<T> &in, SembleMatrix<T> &D);

  template<class T>//in place diagonal matrix from a list of values
  SembleMatrix<T> diag(const SembleVector<T> &in);

  template<class T>//take the diagonal elements and stick them in a vector
  void diag(const SembleMatrix<T> &in, SembleVector<T> &D);

  template<class T> //make a vector from the diagonal elements
  SembleVector<T> diag(const SembleMatrix<T> &in);

  template<class T> //return the diagonal matrix with the square root of the vector elements on the diagonal
  SembleMatrix<T> sqrt(const SembleVector<T> &in);

  template<class T> //return a diagonal matrix with the sqrt of in, in MUST be diagonal
  SembleMatrix<T> sqrt(const SembleMatrix<T> &in);

//mixed vector matrix
  template<class T>
  SembleVector<T> operator*(const SembleMatrix<T> &lhs, const SembleVector<T> &rhs);

  template<class T>
  SembleVector<T> operator*(const SembleVector<T>&lhs, const SembleMatrix<T> &rhs);

  template<class T>
  SembleMatrix<T> outterProduct(const SembleVector<T> &lhs, const SembleVector<T> &rhs);

//itpp extension
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//match unsorted eigenvectors by max overlap and resign so the sign of the overlap is positive
//for use within one timeslice to make sure all steps in the chain have eigenvectors ordered the same
  template<class T>
  void matchEigenVectors(const itpp::Mat<T> &ref, itpp::Mat<T> &in, itpp::Vec<double> &in2);

  template<class T> // will match the rescaled ensemble to the first bin and -- requires rescaling, defaults to rescaling
  void matchEigenVectorsEnsemble(SembleMatrix<T> &inM, SembleVector<double> &inV, bool rescale = true);

  template<class T>// order the eigenvalues by size
  void matchEigenValueSize(itpp::Mat<T> &inM, itpp::Vec<double> &inV, bool reverse = false);

  template<class T> //order the eigenvalues by size -- since we are comparing rescaled bins we need to move rescaled vectors too, defaults to rescaling
  void matchEigenValueSizeEnsemble(SembleMatrix<T> &inM, SembleVector<double> &inV, bool rescale = true);

  template<class T> //match eigen systems of the same size according to some set of ref vecs, its assumed that both ensembles are already sorted w/in themselves -- tslice to tslice
  void matchEigenVectorsEnsemble(const SembleMatrix<T> &ref, SembleMatrix<T> &vecs, SembleVector<double> &vals);

  template<class T> //match the eigen system on some metric to have positive overlap -- tslice to tslice
  void matchEigenVectorsEnsembleMetric(const SembleMatrix<T> &metric, const SembleMatrix<T> &ref, SembleMatrix<T> &vec, SembleVector<double> &val);

  template<class T> //order eigenvalues by magnitude, asc or dsc
  void reorderEigenValues(SembleVector<double> &eVals, SembleMatrix<T> &w, SembleMatrix<T> &eVec, bool orderAsc);

  template<class T> //attach a phase, largest element by modulus is real & positive
  void rephaseEigenVectors(itpp::Mat<T> &evecs);

  template<class T> //enforce a phase convenction on the mean -- calls overloaded rephaseEVectors
  void rephaseEigenVectors(SembleMatrix<T> &vecs);

//rephase doubles, need to overload for complex types
  void rephaseEVectors(SembleMatrix<double> &vecs);

  template<class T> //return a rephase map
  typename std::map<int, typename PromoteScalar<T>::Type> rephaseEigenVectorsEnsembleMap(const SembleMatrix<T> &ref, const SembleMatrix<T> &vecs);

  template<class T>
  typename std::map<int, typename PromoteScalar<T>::Type> rephaseEigenVectorsEnsembleMetricMap(const SembleMatrix<T> &metric, const SembleMatrix<T> &ref, const SembleMatrix<T> &vec);

  template<class T> //return a map, key is ref, data is vecnum, -1 means state didnt map, works for different sizes
  std::map<int, int> makeRemap(const SembleMatrix<T> &ref, const SembleMatrix<T> &vec);

  template<class T>
  std::map<int, int> makeRemapMetric(const SembleMatrix<T> &ref, const SembleMatrix<T> &vec, const SembleMatrix<T> &metric);

  template<class T> //return a map keys in [0,dim-1], data is state, -1 means state didnt map, works for different sizes
  std::map<int, int> makeRemap(const SembleMatrix<T> &ref, const SembleMatrix<T> &vec, const int dim);

  template<class T>
  std::map<int, int> makeRemapMetric(const SembleMatrix<T> &metric, const SembleMatrix<T> &ref, const SembleMatrix<T> &vec, const int dim);

  template<class T>//returns the upper triangular factorized matrix, rescales by default and returns unscaled, A = (F^t)*F 
  void chol(const SembleMatrix<T> &in, SembleMatrix<T> &factorized, bool rescale = true);

  template<class T> //another name for chol
  void cholesky(const SembleMatrix<T> &in, SembleMatrix<T> &factorized, bool rescale = true);

  template<class T>//return the eigenvectors and eigenvalues of a hermitian ensemble -- rescales by default
  void eig_sym(const SembleMatrix<T> &M, SembleVector<double> &vals, SembleMatrix<T> &V, bool rescale = true);

  template<class T> // A = US(V^t) -- rescales by default
  std::string svd(const SembleMatrix<T> &A, SembleMatrix<T> &U, SembleVector<double> &s, SembleMatrix<T> &V, bool rescale = true);

  template<class T> //another name for svd
  std::string singularValueDecomposition(const SembleMatrix<T> &A, SembleMatrix<T> &U, SembleVector<double> &s, SembleMatrix<T> &V, bool rescale = true);

  template<class T> //itpp inversion, I think it uses LU -> unstable for ill conditioned matricies, rescales by default
  void inv(const SembleMatrix<T> &M, SembleMatrix<T> &Minv, bool rescale = true);

  void pseudoInvert(SembleVector<double> &inout, const int dim, bool rescale = true);

  template<class T>//pseudo invert and zero low elements until cond < condMax
  void pseudoInvertCond(const SembleMatrix<T> &M, SembleMatrix<T> &Minv, double condMax = 1e7, bool rescale = true);

  template<class T>//pseudo invert and zero anything below thresh, thresh = 0 -> no resets
  void pseudoInvertValue(const SembleMatrix<T> &M, SembleMatrix<T> &Minv, double thresh = 1e-6, bool rescale = true);

  /*
    The routines invSVD, invSVDNorm and solveLinearSVD are the original routines from reconfit2, they seem to be 'cheating' in 
    that there could exist cases in which we reset different numbers of singular values on a bin by bin basis, the though behind doing this 
    is that they are jackkniffed down (very little variance for reasonable sized ensembles) and so we would only be resetting 
    extreme outliers, this may be true and the routines are left in place as they were in the original code.  
    
    New routines of using the same name with the suffix _H (for honest) are also below, instead of looking bin by bin they will 
    reset only on the mean but are oterwise identical.  It is left to the user to decide which they want to use. 
  */
  
 
  //do a moore penrose inverse (pseudo inverse)
  template<class T>
  SembleMatrix<T> invSVD(const SembleMatrix<T> &A, const double tol, double &avgReset);
  
  //invert using SVD scaling terms to O(1)
  template<class T>
  SembleMatrix<T> invSVDNorm(const SembleMatrix<T> &A, const double tol, double &avgReset);
  
  //solve A.x = b and compute the residual
  template<class T>
  void solveLinearSVD(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b, const double tol, double &res);

  //also compute the covariance for the solution -- BETTER
  template<class T>
  void solveLinearSVD(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b,
		      SembleMatrix<T> &cov_x, const double tol, double &res, int &nReset);

  template<class T>
  SembleMatrix<T> invSVD_H(const SembleMatrix<T> &A, const double tol, double &avgReset);
  
  template<class T>
  SembleMatrix<T> invSVDNorm_H(const SembleMatrix<T> &A, const double tol, double &avgReset);

  template<class T>
  void solveLinearSVD_H(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b, const double tol, double &res);

  template<class T>
  void solveLinearSVD_H(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b,
		      SembleMatrix<T> &cov_x, const double tol, double &res, int &nReset);



//generalized eigenproblem routines and helpers
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//A*V = B*V*Lambda

//***NB I have assumed A and B are unscaled, they are also explicitly symmetrized

  template<class T>
  void genEigCho(const SembleMatrix<T> &A,                    //cholesky solution using LU inversion
                 const SembleMatrix<T> &B,                    //defaults to sorting the system by the order of bin 0
                 SembleMatrix<T> &eVecs,                      //defaults to resigning to the sign convention of bin 0
                 SembleMatrix<T> &wVecs,                      //if sortByVecs == false the system will be ordered by
                 SembleVector<double> &eVals,                      //the size of the eigenValues
                 const bool sortByVecs = true);

  template<class T>                                           //the same but give it the upper triangular factorized version of B
  void genEigChoF(const SembleMatrix<T> &A,                   //this way one doesnt have to factorize the same SembleMatrix
                  const SembleMatrix<T> &F,                   //tmax - tmin times, instead just loop on A
                  SembleMatrix<T> &eVecs,
                  SembleMatrix<T> &wVecs,
                  SembleVector<double> &eVals,
                  const bool sortByVecs = true);

  template<class T>                                           //pass in the upper triangular factorized version of B
  void genEigChoFRef(const SembleMatrix<T> &A,                //first sort within the ensemble according to sortByVecs
                     const SembleMatrix<T> &F,                //then sort the ensemble to match refW
                     const SembleMatrix<T> &refW,             //nb assumed that everything is of the same rank
                     SembleMatrix<T> &eVecs,
                     SembleMatrix<T> &wVecs,
                     SembleVector<double> &eVals,
                     const bool sortByVecs = true);

  template<class T>                                           //svd solution
  void genEigSvdCond(const SembleMatrix<T> &A,                //defaults to sorting the system by order of bin 0
                     const SembleMatrix<T> &B,                //defaults to resigning by bin 0
                     SembleMatrix<T> &eVecs,                  //defaults to resetting until condition number is less than condMax
                     SembleMatrix<T> &wVecs,                  //condMax defaults to 1e7
                     SembleVector<double> &eVals,             //if sortByVecs=false the system is orderd by eigenvalue size
                     const bool sortByVecs = true,
                     const double condMax = 1e7);

  template<class T>
  void genEigSvdValue(const SembleMatrix<T> &A,               //svd solution
                      const SembleMatrix<T> &B,               //same but resets anything under thresh which defaults to 1e-6
                      SembleMatrix<T> &eVecs,
                      SembleMatrix<T> &wVecs,
                      SembleVector<double> &eVals,
                      const bool sortByVecs = true,
                      const double thresh = 1e-6);

  template<class T>                                           //**NB this assumes that rootSinvPlus is rectangular with ncols >= nrows!!, ie reDim rootSinvPlus before sticking it in
  void genEigSvdF(const SembleMatrix<T> &A,                   //assume B is symmetric -> B = U*S*(U^t)
                  const SembleMatrix<double> &rootSinvPlus,   //rootSinvPlus is sqrt(S^-1) with anything that was bad reset
                  const SembleMatrix<T> &U,                   //U is the column vector that diagonalizes B
                  SembleMatrix<T> &eVecs,
                  SembleMatrix<T> &wVecs,
                  SembleVector<double> &eVals,
                  const bool sortByVecs = true);

  template<class T>                                           //**NB this assumes that rootSinvPlus is rectangular with ncols >= nrows!!, ie reDim rootSinvPlus before sticking it in
  void genEigSvdFRef(const SembleMatrix<T> &A,                //assume B is symmetric -> B = U*S*(U^t)
                     const SembleMatrix<double> &rootSinvPlus,//rootSinvPlus is sqrt(S^-1) with anything that was bad reset
                     const SembleMatrix<T> &U,                //U is the column vector that diagonalizes B
                     const SembleMatrix<T> &refW,             //sort the ensemble to match refW
                     SembleMatrix<T> &eVecs,
                     SembleMatrix<T> &wVecs,
                     SembleVector<double> &eVals,
                     const bool sortByVecs = true);
//
//helpers
//change the elements to 1/sqrt(s(i)) and zero any element i > dim_
  void svdResetPseudoInvertRoot(SembleVector<double> &inout, const int dim_, const bool rescale = true);

//reset on the average value -- returns the dimension of the reduced sub space, zero element w/ index greater than ret val
  int svdResetAverageValue(const SembleVector<double> &inout, const double thresh = 1e-6);

//reset on the condition number -- returns the dimension of the reduced sub space, zero element w/ index greater than ret val
  int svdResetAverageCond(const SembleVector<double> &inout, const double condMax = 1e7);

  template<class T>//reset anything within factor*sigma of zero and then rearrange the the matrix so sing[i] > sing[i+1]
  int svdResetSigma(SembleVector<double> &sing, SembleMatrix<T> &U, const double factor = 3, const bool scaleUp = false); // assume jackknifed up data

  template<class T>//reset anything below thresh and then anything consistent with zero
  int svdResetAvgValueAndSigma(SembleVector<double> &inout, SembleMatrix<T> &U, const double thresh = 1e-6, const double sigma = 3, const bool scaleUp = false);

  template<class T>//reset on the condition number and then anything consistent with zero
  int svdResetAvgCondAndSigma(SembleVector<double> &inout, SembleMatrix<T> &U, const double thresh = 1e7, const double sigma = 3, const bool scaleUp = false);

  /*
    NB: The svd algorithm orders the singular values by size, largest first.  Say sigma_i and sigma_j's jackknife distributions overlap,
    because of low statistics or a noisy ensemble, then there could be cases in which the SVD solver (LAPACK:dgesvd) will flip their order
    relative to the rest of the ensemble.  We should undo the flip by performing a reordering according to the U or V unitary matricies.

    template<class T>
    std::string svdMatchSingValsByUandRephaseEnsemble(SembleMatrix<T> &inoutU, SembleVector<double> &inoutSigma, const bool rescale=true);

    does the reordering and keeps track of what it did, sends the result back as a string so that a higher level routine can use it in a log.
   */

  template<class T>//test to see if the svd algorithm moved the singular values around
  std::string svdMatchSingValsByUandRephaseEnsemble(SembleMatrix<T> &inoutU, SembleVector<double> &inoutSigma, const bool rescale = true);

//wrappers
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class T>
  void cholesky(const SembleMatrix<T> &in, SembleMatrix<T> &F, bool rescale = true)
  {
    chol(in, F, rescale);
  }

  template<class T>
  void singularValueDecomposition(const SembleMatrix<T> &A, SembleMatrix<T> &U, SembleVector<double> &S, SembleMatrix<T> &V, bool rescale = true)
  {
    return svd(A, U, S, V, rescale);
  }

  template<class T>
  SembleMatrix<T> inv(const SembleMatrix<T> &in, bool rescale = true)
  {
    SembleMatrix<T> out;
    inv(in, out, rescale);
    return out;
  }

//implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<class T>  //return -T  if real, T* if complex, calls overloaded eConjPhase(T)
  typename PromoteScalar<T>::Type evecConjPhase(const T &in)
  {
    return eConjPhase(in);
  }


  template<class T>//solves R*P = V, P = R^dagger * V, will include the evecConjPhase in the return matrix
  itpp::Mat<T> constructPermutationMatrix(const itpp::Mat<T> &ref, const itpp::Mat<T> &vec)
  {
    itpp::Mat<T> lap = (itpp::hermitian_transpose(ref)*vec);
    const int rdim = ref.cols();
    const int vdim = vec.cols();
    const int bound = (rdim < vdim) ? rdim : vdim;
    std::map<int,int> mapp;
    int mr,mv;
    std::vector<bool> ur(rdim,false), uv(vdim,false);
    T m;
    
    for(int p = 0; p < bound; ++p)
      {
	mr = 0; 
	mv = 0;
	m = T(0.);
	
	for(int v = 0; v < vdim; ++v)
	  if(uv[v])
	    continue;
	  else
	    for(int r = 0; r < rdim; ++r)
	      if(ur[r])
		continue;
	      else
		if(fabs(lap(r,v)) >= m)
		  {
		    m = fabs(lap(r,v));
		    mr = r;
		    mv = v;
		  }
	ur[mr] = true;
	uv[mv] = true;
	mapp[mr] = mv;      
      }
    
    itpp::Mat<T> r(lap);
    
    r.zeros();
  
    //fill in a generalized permutaion matrix, phases included
    std::map<int,int>::const_iterator it;
    for(it = mapp.begin(); it != mapp.end(); ++it)
      r(it->first,it->second) = toScalar(evecConjPhase(lap(it->first,it->second)/fabs(lap(it->first,it->second))));
    
    return r;  //this is rdim X vdim
  }



//semble linear algebra
  template<class T>
  SembleMatrix<T> hermitianConjugate(const SembleMatrix<T> &M)
  {
    SembleMatrix<T> dum(M);
    dum.hermitianConjugate();
    return dum;
  }

  template<class T>
  SembleMatrix<T> adj(const SembleMatrix<T> &M)
  {
    return hermitianConjugate(M);
  }

  template<class T>
  SembleMatrix<T> adjoint(const SembleMatrix<T> &M)
  {
    return hermitianConjugate(M);
  }

  template<class T>
  SembleMatrix<T> transpose(const SembleMatrix<T> &M)
  {
    SembleMatrix<T> dum(M);
    dum.transpose();
    return dum;
  }

  template<class T>
  SembleMatrix<T> tran(const SembleMatrix<T> &M)
  {
    return transpose(M);
  }

  template<class T> //loads into first row
  SembleMatrix<T> toMatrix(const SembleVector<T> &V)
  {
    SembleMatrix<T> dum(V.bins(), 1, V.elem());
    int size = V.elem();

    for(int elem = 0; elem < size; ++elem)
      dum.loadEnsemElement(0, elem, V(elem));

    return dum;
  }

  template<class T> //loads the corresponding row or col
  SembleVector<T> toVector(const SembleMatrix<T> &M)
  {
    if((M.getM() != 1) && (M.getN() != 1))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    int dim_ = 0, bin_ = M.getB(), index = 0, *row, *col;
    bool dum_dum = false;

    if(M.getM() == 1)
      {
        col = new int(0);
        row = &index;
        dim_ = M.getN();
      }
    else
      {
        row = new int(0);
        col = &index;
        dum_dum = true;
        dim_ = M.getM();
      }

    SembleVector<T> dum(bin_, dim_);

    for(index = 0; index < dim_; ++index)
      dum.loadEnsemElement(index, M(*row, *col));

    if(dum_dum)
      delete row;
    else
      delete col;

    row = NULL;
    col = NULL;

    return dum;
  }

  template<class T> //assumes a column vector
  SembleMatrix<T> adjoint(const SembleVector<T> &V)
  {
    SembleMatrix<T> dum(V.bins(), V.elem(), 1);
    int size = V.elem();

    for(int elem = 0; elem < size; ++elem)
      dum.loadEnsemElement(elem, 0, V(elem));

    return hermitianConjugate(dum);
  }

  template<class T> //assumes a column vector
  SembleMatrix<T> adj(const SembleVector<T> &V)
  {
    return adjoint(V);
  }

  template<class T> //assumes a column vector
  SembleMatrix<T> transpose(const SembleVector<T> &V)
  {
    SembleMatrix<T> dum(V.bins(), V.elem(), 1);
    int size = V.elem();

    for(int elem = 0; elem < size; ++elem)
      dum.loadEnsemElement(elem, 0, V(elem));

    return transpose(dum);
  }

  template<class T> //assumes a column vector
  SembleMatrix<T> tran(const SembleVector<T> &V)
  {
    return transpose(V);
  }

  template<class T> //explicitly symmetrize a matrix
  void symmetrize(const SembleMatrix<T> &in, SembleMatrix<T> &out)
  {
    out = 0.5 * (in + hermitianConjugate(in));
  }

  template<class T>
  SembleMatrix<T> symmetrize(const SembleMatrix<T> &in)
  {
    SembleMatrix<T> dum;
    symmetrize(in, dum);
    return dum;
  }

  template<class T>
  SembleVector<T> getRow(const SembleMatrix<T> &in, const int row)
  {
    int dim = in.getM();
    SembleVector<T> dum(in.getB(), dim);

    for(int col = 0; col < dim; ++col)
      dum.loadEnsemElement(col, in.getEnsemElement(row, col));

    return dum;
  }

  template<class T>
  SembleVector<T> getCol(const SembleMatrix<T> &in, const int col)
  {
    int dim = in.getN();
    SembleVector<T> dum(in.getB(), dim);

    for(int row = 0; row < dim; ++row)
      dum.loadEnsemElement(row, in.getEnsemElement(row, col));

    return dum;
  }

  template<class T>
  SembleMatrix<T> outterProduct(const SembleVector<T> &lhs, const SembleVector<T> &rhs)
  {
    if(lhs.getB() != rhs.getB())
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Don't have same num bins, exiting.." << std::endl;
        exit(1);
      }

    int N = lhs.getN();
    int M = rhs.getN();

    SembleMatrix<T> dum(lhs.getB(), N, M);

    for(int row = 0; row < N; ++row)
      for(int col = 0; col < M; ++col)
        dum.loadEnsemElement(row, col, lhs.getEnsemElement(row)*rhs.getEnsemElement(col));

    return dum;
  }

//mixed vector matrix
  template<class T>
  SembleVector<T> operator*(const SembleMatrix<T> &lhs, const SembleVector<T> &rhs)
  {
    if((rhs.getN() != lhs.getN()) || (rhs.getB() != lhs.getB()))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleVector<T> V = rescaleEnsemDown(rhs);
    SembleMatrix<T> M = rescaleEnsemDown(lhs);

    int bin_ = rhs.getB();

    for(int bin = 0; bin < bin_; ++bin)
      V[bin] = M[bin] * V[bin];

    V.rescaleEnsemUp();
    return V;
  }

  template<class T>
  SembleVector<T> operator*(const SembleVector<T>&lhs, const SembleMatrix<T> &rhs)
  {
    return toVector(toMatrix(lhs) * rhs);
  }

//itpp extension
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//match unsorted eigenvectors by max overlap and resign so the sign of the overlap is positive
//for use within one timeslice to make sure all steps in the chain have eigenvectors ordered the same
  template<class T>
  void matchEigenVectors(const itpp::Mat<T> &ref, itpp::Mat<T> &vecs, itpp::Vec<double> &vals)
  {
    if((ref.cols() != vecs.cols()) || (ref.rows() != vecs.rows()) || (vals.size() != vecs.cols()))
      {
        std::cout << "This matching only matches full space to full space" << std::endl;
      }


    itpp::Mat<T> laps, vec_cp(vecs), val_cp(vals);
    int maxr, maxv, dim = ref.cols();
    T max_lap;
    std::map<int, int> lap_map;
    std::vector<bool> ur(dim, false), uv(dim, false);

    laps = itpp::hermitian_transpose(ref) * vecs;

    for(int state = 0; state < dim; ++ state)
      {
        max_lap = 0.;
        maxr = 0;
        maxv = 0;

        for(int ref = 0; ref < dim; ++ref)
          {
            if(ur[ref])
              continue;

            for(int vec = 0; vec < dim; ++vec)
              {
                if(uv[vec])
                  continue;

                if(fabs(laps(ref, vec)) >= max_lap)
                  {
                    max_lap = fabs(laps(ref, vec));
                    maxr = ref;
                    maxv = vec;
                  }
              }
          }

        lap_map[maxr] = maxv;
        ur[maxr] = true;
        uv[maxv] = true;
      }

    std::map<int,int>::const_iterator it;
    for(it = lap_map.begin(); it != lap_map.end(); ++it)
      {
	if(it->first == it->second)
	  {
	    T factor = (laps(it->first,it->second)/fabs(laps(it->first, it->second)));
	    if(factor != T(1.))
	      vecs.set_col(it->first,factor*vec_cp.get_col(it->second));
	  }
	else
	  {
	    vals.set(it->first, val_cp.get(it->second));
	    T factor = (laps(it->first,it->second)/fabs(laps(it->first, it->second)));
	    vecs.set_col(it->first,factor*vec_cp.get_col(it->second));
	  }
      }
  }

  template<class T> // will match the rescaled ensemble to the first bin and -- requires rescaling, defaults to rescaling
  void matchEigenVectorsEnsemble(SembleMatrix<T> &inM, SembleVector<double> &inV, bool rescale = true)
  {
    if(inV.getB() != inM.getB())
      {
        std::cout << "Bins don't match in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    if(rescale)
      {
        inM.rescaleEnsemDown();
        inV.rescaleEnsemDown();
      }

    int bin_ = inV.getB();

    //set a phase on bin 0
    rephaseEigenVectors(inM[0]);

    //set an ordering on bin 0
    matchEigenValueSize(inM[0], inV[0], true);

    //enforce the phase/ordering
    for(int bin = 1; bin < bin_; ++bin)
      {
        matchEigenVectors(inM[0], inM[bin], inV[bin]);
      }

    if(rescale)
      {
        inM.rescaleEnsemUp();
        inV.rescaleEnsemUp();
      }
  }

  template<class T>// order the eigenvalues by size
  void matchEigenValueSize(itpp::Mat<T> &inM, itpp::Vec<double> &inV, bool reverse = false)
  {
    if(inV.size() != inM.cols())
      {
        std::cout << "This matching only matches full space to full space" << std::endl;
        exit(1);
      }

    itpp::Vec<int> sorted = itpp::sort_index(inV);

    if(reverse)
      sorted = itpp::reverse(sorted);

    int dim = sorted.size();
    itpp::Vec<T> inVc = inV;
    itpp::Mat<T> inMc = inM;

    for(int col = 0; col < dim; ++col)
      if(col != sorted(col))
        {
          inM.set_col(col, inMc.get_col(sorted(col)));
          inV.set(col, inVc.get(sorted(col)));
        }
  }

//NB, this sorts bin by bin so one would expect to see things like the average and variance change if you run this filter
  template<class T> //order the eigenvalues by size -- since we are comparing rescaled bins we need to move rescaled vectors too, defaults to rescaling
  void matchEigenValueSizeEnsemble(SembleMatrix<T> &inM, SembleVector<double> &inV, bool rescale = true)
  {
    int bin_ = inV.getB();

    if(bin_ != inM.getB())
      {
        std::cout << "Bins don't match in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    if(rescale)
      {
        inM.rescaleEnsemDown();
        inV.rescaleEnsemDown();
      }

    for(int bin = 0; bin < bin_; ++bin)
      matchEigenValueSize(inM[bin], inV[bin]);

    if(rescale)
      {
        inM.rescaleEnsemUp();
        inV.rescaleEnsemUp();
      }
  }

  template<class T> //match eigen systems of the same size  
  void matchEigenVectorsEnsemble(const SembleMatrix<T> &ref, SembleMatrix<T> &vecs, SembleVector<double> &vals)
  {

    int dim_ = ref.getM();
    itpp::Mat<T> laps = mean(adj(ref) * vecs); // this checks the dim
    int maxr, maxv, dim = laps.cols();
    std::vector<bool> ur(dim, false), uv(dim, false);
    std::map<int, int> lap_map;
    T max_lap;

    if(laps.rows() != laps.cols())
      {
        std::cout << "Dimensionally incorrect system in " << __PRETTY_FUNCTION__ << __LINE__ << __FILE__ << std::endl;
        exit(1);
      }


    for(int state = 0; state < dim; ++ state)
      {
        max_lap = 0.;
        maxr = 0;
        maxv = 0;

        for(int ref = 0; ref < dim; ++ref)
          {
            if(ur[ref])
              continue;

            for(int vec = 0; vec < dim; ++vec)
              {
                if(uv[vec])
                  continue;

                if(fabs(laps(ref, vec)) >= max_lap)
                  {
                    max_lap = fabs(laps(ref, vec));
                    maxr = ref;
                    maxv = vec;
                  }
              }
          }

        lap_map[maxr] = maxv;
        ur[maxr] = true;
        uv[maxv] = true;
      }

    int row_ = vecs.rows();
    SembleMatrix<T> vec_cp(vecs);
    SembleVector<double> val_cp(vals);
    std::map<int,int>::const_iterator it;

    //technically the factor is an ensemble object, the old code didn't use an ensemble so we won't either, we 
    //assume that since this is all done jackknife down that the ensemble factor is basically an ensemble of the mean..
    for(it = lap_map.begin(); it != lap_map.end(); ++it)
      {
	if(it->first == it->second)
	  {
	    typename PromoteScalar<T>::Type factor = toScalar(laps(it->first,it->second)/fabs(laps(it->first,it->second)));
	    if(toScalar(factor) != T(1.))
	      for(int row(0); row < row_; ++row)
		vecs.loadEnsemElement(row,it->first,factor*vec_cp.getEnsemElement(row,it->second));
	  }
	else
	  {
	    typename PromoteScalar<T>::Type factor = toScalar(laps(it->first,it->second)/fabs(laps(it->first,it->second)));
	    vals.loadEnsemElement(it->first,val_cp.getEnsemElement(it->second));
	    for(int row(0); row < row_; ++row)
		vecs.loadEnsemElement(row,it->first,factor*vec_cp.getEnsemElement(row,it->second));
	  }

      }
  }

  template<class T> //match the eigen system on some metric to have positive overlap V_A^t*M*V_B = delta_i_j
  void matchEigenVectorsEnsembleMetric(const SembleMatrix<T> &metric, const SembleMatrix<T> &ref, SembleMatrix<T> &evec, SembleVector<double> &eval)
  {
    matchEigenVectorsEnsemble(adj(metric)*ref, evec, eval);
  }

  template<class T>
  void reorderEigenValues(SembleVector<double> &eVals, SembleMatrix<T> &w, SembleMatrix<T> &eVec, bool orderAsc)
  {
    int bin_ = eVals.getB(), row_ = eVec.getN(), col_ = w.getM();

    if((bin_ != w.getB()) || (bin_ != eVec.getB()) || (w.getN() != col_) || (eVals.getN() != col_)
        || (col_ != eVec.getM()))
      {
        std::cout << "Dimensionally incorrect problem in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    itpp::Vec<T> vals = mean(eVals);
    itpp::Vec<int> index = itpp::sort_index(mean(eVals));
    SembleMatrix<T> eVecp = eVec, wcp = w;
    SembleVector<double> eValcp = eVals;

    if(!!!orderAsc)
      index = itpp::reverse(index);

    for(int col = 0; col < col_; ++col)
      {
        eVals.loadEnsemElement(col, eValcp.getEnsemElement(index(col)));

        for(int row = 0; row < row_; ++row)
          {
            eVec.loadEnsemElement(row, col, eVecp.getEnsemElement(row, index(col)));

            if(row < col_)
              w.loadEnsemElement(row, col, wcp.getEnsemElement(row, index(col)));
          }

      }
  }

  template<class T> //attach a phase, largest element by modulus is real & positive
  void rephaseEigenVectors(itpp::Mat<T> &evecs)
  {
    int nvecs = evecs.cols(); 
    int nelem = evecs.rows();

    for(int vec = 0; vec < nvecs; ++vec)
      {
	itpp::Vec<T> v = evecs.get_col(vec);
	T max = T(0.);

	for(int elem = 0; elem < nelem; ++elem)
	  if(abs(v(elem)) > abs(max))
	    max = v(elem);

	v *= max/abs(max);
	evecs.set_col(vec,v);
      }
  }

  template<class T>
  void rephaseEigenVectors(SembleMatrix<T> &vecs)
  {
    rephaseEVectors(vecs);
  }

//these can be of differing sizes, key is where it goes, data is phase
  template<class T>
  typename std::map<int, typename  PromoteScalar<T>::Type> rephaseEigenVectorsEnsembleMap(const SembleMatrix<T> &ref, const SembleMatrix<T> &vec)
  {

    typename std::map<int, typename PromoteScalar<T>::Type> ret;
    std::vector<bool> ur(ref.getM(), false), uv(vec.getM(), false);
    itpp::Mat<T> lap = mean(adj(ref) * vec);
    int maxr, maxv;
    T maxlap;

    int rdim = ref.getM(), vdim = vec.getM();
    int bound = (rdim < vdim) ? rdim : vdim;

    for(int state = 0; state < bound; ++state)
      {
        maxlap = T(0.);
        maxr = 0;
        maxv = 0;

        for(int r = 0; r < rdim; ++r)
          {
            if(ur[r])
              continue;

            for(int v = 0; v < vdim; ++v)
              {
                if(uv[v])
                  continue;

                if(fabs(lap(r, v)) >= maxlap)
                  {
                    maxr = r;
                    maxv = v;
                    maxlap = fabs(lap(r, v));
                  }
              }

            ur[maxr] = true;
            uv[maxv] = true;
	    
	    maxlap = toScalar(evecConjPhase(lap(maxr,maxv)/maxlap));
	    
            ret.insert(make_pair(maxv, toScalar(maxlap)));

            if(std::isinf(maxlap))
              ret[maxv] = toScalar(T((maxlap > 0) - (maxlap < 0))); //will fail if inf complex, pulls sign if inf
          }
      }

    //fill in a phase of one in the case that the dims arent the same
    for(int i = 0; i < vdim; ++i)
      if(!!!uv[i])
        ret[i] = toScalar(1.0);

    return ret;
  }

  template<class T>
  typename std::map<int, typename PromoteScalar<T>::Type> rephaseEigenVectorsEnsembleMetricMap(const SembleMatrix<T> &metric, const SembleMatrix<T> &ref, const SembleMatrix<T> &vecs)
  {
    return rephaseEigenVectorsEnsembleMap(adj(metric) * ref, vecs);
  }


  template<class T>
  std::map<int, int> makeRemap(const SembleMatrix<T> &ref, const SembleMatrix<T> &vec)
  {
    std::map<int, int> ret;
    itpp::Mat<T> lap = mean(adj(ref) * vec);
    std::vector<bool> ur(ref.getM(), false), uv(vec.getM(), false);
    int maxr, maxv;
    T maxlap;

    int rdim = ref.getM(), vdim = vec.getM();
    int bound = (rdim < vdim) ? rdim : vdim;

    for(int state = 0; state < bound; ++state)
      {
        maxr = 0;
        maxv = 0;
        maxlap = T(0.);

        for(int r = 0; r < rdim; ++r)
          {
            if(ur[r])
              continue;

            for(int v = 0; v < vdim; ++v)
              {
                if(uv[v])
                  continue;

                if(fabs(lap(r, v)) >= maxlap)
                  {
                    maxr = r;
                    maxv = v;
                    maxlap = fabs(lap(r, v));
                  }
              }
          }

        ur[maxr] = true;
        uv[maxv] = true;
        ret[maxr] = maxv;
      }

    /* Three cases
       rdim = vdim - do nothing
       rdim < vdim - match lowest to lowest, arbitrary
       rdim > vdim - match then fill rest of unused states with -1 as doesnt exist
    */

    if(rdim != vdim)
      {

        if(rdim < vdim)
          {
            std::vector<int> availablev, availabler;

            for(int index = 0; index < vdim; ++index)
              {
                if(!!!uv[index])
                  availablev.push_back(index);

                if(index < rdim)
                  {
                    if(!!!ur[index])
                      availabler.push_back(index);
                  }
                else
                  availabler.push_back(index);
              }

            for(int index = 0; index < availablev.size(); ++index)
              ret[availabler[index]] = availablev[index];

          }
        else
          {
            for(int index = 0; index < rdim; ++index)
              if(!!!ur[index])
                ret[index] = -1;
          }
      }

    return ret;
  }

  template<class T>
  std::map<int, int> makeRemapMetric(const SembleMatrix<T> &ref, const SembleMatrix<T> &vec, const SembleMatrix<T> &metric)
  {
    return makeRemap(adj(metric) * ref, vec);
  }

  template<class T> //return a map keys in [0,dim-1], data is state, -1 means state didnt map, works for different sizes
  std::map<int, int> makeRemap(const SembleMatrix<T> &ref, const SembleMatrix<T> &vec, const int dim)
  {
    std::map<int, int> ret = makeRemap(ref, vec);
    std::map<int, int>::const_iterator it;
    std::vector<bool> used(dim, false);
    bool sane = true;

    for(it = ret.begin(); it != ret.end(); ++it)
      {
        if((it->first >= dim) || (it->second >= dim))
          {
            sane = false;
            break;
          }

        used[it->first] = true;
      }

    if(!!!sane)
      {
        std::cout << "remapping error " << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
      }

    for(int index = 0; index < dim; ++index)
      if(!!!used[index])
        ret.insert(std::pair<int, int>(index, -1));

    return ret;
  }

  template<class T>
  std::map<int, int> makeRemapMetric(const SembleMatrix<T> &metric, const SembleMatrix<T> &ref, const SembleMatrix<T> &vec, const int dim)
  {
    return makeRemap(adj(metric)*ref,vec, dim);
  }



  template<class T>//returns the upper triangular factorized matrix, rescales by default and returns unscaled
  void chol(const SembleMatrix<T> &in, SembleMatrix<T> &factorized, bool rescale = true)
  {
    SembleMatrix<T> dum(in);
    bool success = true;

    if(rescale)
      dum.rescaleEnsemDown();

    factorized.reDim(dum.getB(), dum.getN(), dum.getM());
    factorized.zeros();

    int bin_ = dum.getB();

    for(int bin = 0; bin < bin_; ++bin)
      if(!!!itpp::chol(dum[bin], factorized[bin]))
        success = false;

    if(!!!success)
      {
        std::cout << "Cholesky Factorization Failed in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    if(rescale)
      factorized.rescaleEnsemUp();
  }

  template<class T>//return the eigenvectors and eigenvalues of a hermitian ensemble -- rescales by default, these should be matched bin to bin
  void eig_sym(const SembleMatrix<T> &M, SembleVector<double> &vals, SembleMatrix<T> &V, bool rescale = true)
  {
    if(M.getM() != M.getN())
      {
        std::cout << "eig_sym only works for square matricies " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> dum(M);
    bool success = true;

    if(rescale)
      dum.rescaleEnsemDown();

    int bin_ = dum.getB();

    vals.reDim(bin_, dum.getM());
    vals.zeros();
    V.reDim(bin_, dum.getM(), dum.getM());
    V.zeros();


    for(int bin = 0; bin < bin_; ++bin)
      if(!!!itpp::eig_sym(dum[bin], vals[bin], V[bin]))
        success = false;

    //match everything to bin 0 and put bin 0 phase convention
    matchEigenVectorsEnsemble(V,vals,false);


    //now do the largest element by modulus overall phase convention on the mean
    //this enforces a new phase convention, its a bit stupid since another new one gets enforced later for the GEVP
    //but eig_sym is a separate routine so it should have some phase convention attached to it for the sake of linear algebra
    rephaseEigenVectors(V); 


    if(!!!success)
      {
        std::cout << "eig_sym Failed in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    if(rescale)
      {
        vals.rescaleEnsemUp();
        V.rescaleEnsemUp();
      }
  }

  template<class T> // A = US(V^t) -- rescales by default
  std::string svd(const SembleMatrix<T> &A, SembleMatrix<T> &U, SembleVector<double> &S, SembleMatrix<T> &V, bool rescale = true)
  {
    if(A.getM() != A.getN())
      {
        std::cout << "svd only implemented for square matricies " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> dum(A);
    bool success = true;

    if(rescale)
      dum.rescaleEnsemDown();

    int bin_ = dum.getB();
    int dim_ = dum.getM();

    U.reDim(bin_, dim_, dim_);
    U.zeros();
    S.reDim(bin_, dim_);
    S.zeros();
    V.reDim(bin_, dim_, dim_);
    V.zeros();

    for(int bin = 0; bin < bin_; ++bin)
      if(!!!itpp::svd(dum[bin], U[bin], S[bin], V[bin]))
        {
          std:: cout << bin << std::endl;
          success = false;
          break;
        }

    if(!!!success)
      {
        std::cout << "svd failed in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    //fake a ref using the mean
    itpp::Mat<T> MM(mean(A)), MU, MV;
    itpp::Vec<double> MS;
    itpp::svd(MM, MU, MS, MV);
    std::string str;
    str = svdMatchSingValsByUandRephaseEnsemble(MU, U, S, false); //don't rescale, it already is/isn't

    if(rescale)
      {
        U.rescaleEnsemUp();
        S.rescaleEnsemUp();
        V.rescaleEnsemUp();
      }

    return str;
  }

  template<class T> //itpp inversion, I think it uses LU -> unstable for ill conditioned matricies, rescales by default
  void inv(const SembleMatrix<T> &M, SembleMatrix<T> &Minv, bool rescale = true)
  {
    int bin_ = M.getB();
    int dim_ = M.getN();

    if(M.getM() != dim_)
      {
        std::cout << "inv only defined for a square matrix in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> dum(M);
    bool success = true;

    if(rescale)
      dum.rescaleEnsemDown();

    Minv.reDim(bin_, dim_, dim_);

    for(int bin = 0; bin < bin_; ++bin)
      if(!!!itpp::inv(dum[bin], Minv[bin]))
        success = false;

    if(!!!success)
      {
        std::cout << "itpp::inv failed in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    if(rescale)
      Minv.rescaleEnsemUp();
  }

  template<class T>
  void pseudoInvertCond(const SembleMatrix<T> &M, SembleMatrix<T> &Minv, double condMax = 1e7, bool rescale = true)
  {
    int bin_ = M.getB();
    int dim_ = M.getN();

    if(M.getM() != dim_)
      {
        std::cout << "inv only defined for a square matrix in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> U, V;
    SembleVector<double> s;

    svd(M, U, s, V, rescale);

    itpp::Vec<T> ss = mean(s);
    int rdim = dim_;

    while(true)
      {
        if(ss(0) / ss(rdim - 1) < condMax)
          break;

        --rdim;

        if(rdim == 0)
          {
            std::cout << "All Null Space in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
            exit(1);
          }
      }

    pseudoInvert(s, rdim, rescale);

    Minv = adj(U * diag(s) * adj(V));
  }

  template<class T>
  void pseudoInvertValue(const SembleMatrix<T> &M, SembleMatrix<T> &Minv, double thresh = 1e-6, bool rescale = true)
  {
    int bin_ = M.getB();
    int dim_ = M.getN();

    if(M.getM() != dim_)
      {
        std::cout << "inv only defined for a square matrix in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> U, V;
    SembleVector<double> s;

    svd(M, U, s, V, rescale);

    itpp::Vec<T> ss = mean(s);
    int rdim = dim_;

    while(true)
      {
        if(ss(rdim - 1) > thresh)
          break;

        --rdim;

        if(rdim == 0)
          {
            std::cout << "All Null Space in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
            exit(1);
          }
      }

    pseudoInvert(s, rdim, rescale);

    Minv = adj(U * diag(s) * adj(V));
  }

  template<class T> //create a diagonal matrix from a list of values
  void diag(const SembleVector<T> &in, SembleMatrix<T> &D)
  {
    int bin_ = in.getB(), dim_ = in.getN();
    D.reDim(bin_, dim_, dim_);
    D.zeros();

    for(int bin = 0; bin < bin_; ++bin)
      D[bin] = itpp::diag(in[bin]);
  }

  template<class T>
  SembleMatrix<T> diag(const SembleVector<T> &in)
  {
    SembleMatrix<T> dum;
    diag(in, dum);
    return dum;
  }

  template<class T>
  void diag(const SembleMatrix<T> &in, SembleVector<T> &D)
  {
    int bin_ = in.getB(), dim_ = in.getM();

    if(dim_ != in.getN())
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " only supports square matrix" << std::endl;
        exit(1);
      }

    D.reDim(bin_, dim_);

    for(int elem = 0; elem < dim_; ++elem)
      D.loadEnsemElement(elem, in.getEnsemElement(elem, elem));
  }

  template<class T>
  SembleVector<T> diag(const SembleMatrix<T> &in)
  {
    SembleVector<T> dum;
    diag(in, dum);
    return dum;
  }

  template<class T>
  SembleMatrix<T> sqrt(const SembleVector<T> &in)
  {
    return sqrt(diag(in));
  }

  template<class T>
  SembleMatrix<T> sqrt(const SembleMatrix<T> &in)
  {
    SembleMatrix<T> dum(in);
    dum.rescaleSembleDown();

    int dim_ = (in.rows() < in.cols()) ? in.rows() : in.cols();
    int bin_ = in.bins();

    for(int bin = 0; bin < bin_; ++bin)
      for(int dim = 0; dim < dim_; ++dim)
        dum[bin](dim, dim) = std::sqrt(dum[bin](dim, dim));

    dum.rescaleSembleUp();

    return dum;
  }


  /*
    The routines invSVD, invSVDNorm and solveLinearSVD are the original routines from reconfit2, they seem to be 'cheating' in 
    that there could exist cases in which we reset different numbers of singular values on a bin by bin basis, the thought behind doing this 
    is that they are jackkniffed down (very little variance for reasonable sized ensembles) and so we would only be resetting 
    extreme outliers (or that it doesn't actually happen). This may be true and the routines are left in place as they were in the
    original code.  
    
    New routines of using the same name with the suffix _H (for honest) are also below, instead of looking bin by bin they will 
    reset only on the mean but are oterwise identical.  It is left to the user to decide which they want to use. 
  */

  //do a moore penrose inverse (pseudo inverse)
  template<class T>
  SembleMatrix<T> invSVD(const SembleMatrix<T> &A, const double tol, double &avgReset)
  {
    //assumes a square matrix

    SembleMatrix<T> U,V;
    SembleVector<double> s, sinv;
    int nreset(0);
    int bins = A.getB();
    int elems = A.getN();

    svd(A,U,s,V);
    s.rescaleEnsemDown();
    sinv = 0.*s;
    
    for(int bin = 0; bin < bins; ++bin)
      for(int elem = 0; elem < elems; ++elem)	
	if(s[bin][elem] > tol*s[bin][0])                  //this could introduce a bias
	  sinv[bin][elem] = 1.0/s[bin][elem];
	else
	  ++nreset;
    
    avgReset = double(nreset)/double(bins);
    sinv.rescaleEnsemUp();
    
    return V*diag(sinv)*adj(U);
  }
  
  //invert using SVD scaling terms to O(1)
  template<class T>
  SembleMatrix<T> invSVDNorm(const SembleMatrix<T> &A, const double tol, double &avgReset)
  {
    SembleMatrix<T> U,V,Normed;
    SembleVector<double> s;

    Normed = rescaleEnsemDown(A);
    int bins(A.getB()), elems(A.getN());
    if(elems != A.getM())
      {
	std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " non-square matrix, exiting" << std::endl;
	exit(1);
      }
    
    itpp::Mat<T> UU,VV;
    itpp::Vec<double> ss;

    for(int bin = 0; bin < bins; ++bin)
      {
	itpp::svd(Normed[bin],UU,ss,VV);
	for(int elem = 0; elem < elems; ++elem)
	  ss[elem] = 1./std::sqrt(fabs(ss[elem]));
      
	Normed[bin] = itpp::diag(ss)*Normed[bin]*itpp::diag(ss);
      }

    Normed.rescaleEnsemUp();

    return invSVD(Normed,tol,avgReset);
  }
  
  //solve A.x = b and compute the residual
  template<class T>
  void solveLinearSVD(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b, const double tol, double &res)
  {
    double avgResets;
    SembleMatrix<T> Apinv = invSVD(A,tol,avgResets);
    x.setB(A.getB());
    x.setN(A.getM());
    x = Apinv*b;
    SembleVector<T> rres;
    res = A*x - b;
    typename PromoteEnsem<T>::Type dum = rres.dot(rres);
    res = toScalar(mean(dum));
  }

  //also compute the covariance for the solution -- BETTER
  template<class T>
  void solveLinearSVD(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b,
		      SembleMatrix<T> &cov_x, const double tol, double &res, int &nReset)
  {
    double avgResets;
    SembleMatrix<T> Apinv = invSVD(A,tol,avgResets);
    x.setB(A.getB());
    x.setN(A.getM());
    x = Apinv*b;
    SembleVector<T> rres;
    res = A*x - b;
    typename PromoteEnsem<T>::Type dum = rres.dot(rres);
    res = std::sqrt(toScalar(mean(dum)));
    nReset = int(avgResets);
  }

  template<class T>
  SembleMatrix<T> invSVD_H(const SembleMatrix<T> &A, const double tol, double &avgReset);
  
  template<class T>
  SembleMatrix<T> invSVDNorm_H(const SembleMatrix<T> &A, const double tol, double &avgReset);

  template<class T>
  void solveLinearSVD_H(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b, const double tol, double &res);

  template<class T>
  void solveLinearSVD_H(const SembleMatrix<T> &A, SembleVector<T> &x, const SembleVector<T> &b,
		      SembleMatrix<T> &cov_x, const double tol, double &res, int &nReset);



//generalized eigenproblem routines and helpers
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//cho
  template<class T>
  void genEigCho(const SembleMatrix<T> &A, const SembleMatrix<T> &B, SembleMatrix<T> &eVecs,
                 SembleMatrix<T> &wVecs, SembleVector<double> &eVals, const bool sortByVecs = true)
  {
    int bin_ = A.bins(), row_ = A.rows(), col_ = A.cols();

    if((bin_ != B.bins()) || (row_ != B.cols()) || (row_ != B.rows()) || (row_ != col_))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    eVecs.reDim(bin_, row_, col_);
    eVecs.zeros();
    wVecs = eVecs;
    eVals.reDim(bin_, row_);
    eVals.zeros();

    SembleMatrix<T> F(eVecs), Finv(eVecs), AA, BB;

    BB = symmetrize(B);
    AA = symmetrize(A);

    chol(BB, F);
    inv(F, Finv);

    AA = symmetrize(inv(adj(F)) * AA * Finv);
    eig_sym(AA, eVals, wVecs);

    if(sortByVecs)
      matchEigenVectorsEnsemble(wVecs, eVals);
    else
      matchEigenValueSizeEnsemble(wVecs, eVals);

    eVecs = Finv * wVecs;

    //enforce a phase convention
    //rephaseEigenVectors(eVecs);
  }

//cho
  template<class T>
  void genEigChoF(const SembleMatrix<T> &A, const SembleMatrix<T> &F, SembleMatrix<T> &eVecs,
                  SembleMatrix<T> &wVecs, SembleVector<double> &eVals, const bool sortByVecs = true)
  {
    int bin_ = A.bins(), row_ = A.rows(), col_ = A.cols();

    if((bin_ != F.bins()) || (row_ != F.cols()) || (row_ != F.rows()) || (row_ != col_))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    eVecs.reDim(bin_, row_, col_);
    eVecs.zeros();
    wVecs = eVecs;
    eVals.reDim(bin_, row_);
    eVals.zeros();

    SembleMatrix<T> AA, Finv;

    AA = symmetrize(A);
    inv(F, Finv);
    AA = symmetrize(inv(adj(F)) * AA * Finv);
    eig_sym(AA, eVals, wVecs);

    if(sortByVecs)
      matchEigenVectorsEnsemble(wVecs, eVals);
    else
      matchEigenValueSizeEnsemble(wVecs, eVals);

    eVecs = Finv * wVecs;

    //enforce a phase convention
    //rephaseEigenVectors(eVecs);
  }

//cho
  template<class T>
  void genEigChoFRef(const SembleMatrix<T> &A, const SembleMatrix<T> &F, const SembleMatrix<T> &refW,
                     SembleMatrix<T> &eVecs, SembleMatrix<T> &wVecs, SembleVector<double> &eVals, const bool sortByVecs = true)
  {
    int bin_ = A.bins(), row_ = A.rows(), col_ = A.cols();

    if((bin_ != F.bins()) || (row_ != F.cols()) || (row_ != F.rows()) || (row_ != col_))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    eVecs.reDim(bin_, row_, col_);
    eVecs.zeros();
    wVecs = eVecs;
    eVals.reDim(bin_, row_);
    eVals.zeros();

    SembleMatrix<T> AA, Finv;

    AA = symmetrize(A);
    inv(F, Finv);

    AA = symmetrize(inv(adj(F)) * AA * Finv);
    eig_sym(AA, eVals, wVecs);

    //this orders within the system
    if(sortByVecs)
      matchEigenVectorsEnsemble(wVecs, eVals);
    else
      matchEigenValueSizeEnsemble(wVecs, eVals);

    //this orders the ensemble
    matchEigenVectorsEnsemble(refW, wVecs, eVals);

    eVecs = Finv * wVecs;

    //enforce a phase convention
    //rephaseEigenVectors(eVecs);
  }


  template<class T>
  void genEigSvdCond(const SembleMatrix<T> &A, const SembleMatrix<T> &B, SembleMatrix<T> &eVecs, SembleMatrix<T> &wVecs,
                     SembleVector<double> &eVals, const bool sortByVecs = true, const double condMax = 1e7)
  {
    int bin_ = A.bins(), row_ = A.rows(), col_ = A.cols();

    if((bin_ != B.bins()) || (row_ != B.cols()) || (row_ != B.rows()) || (row_ != col_))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> U, Udag, M, RSIP, AA, BB;
    SembleVector<double> RSinvP;

    BB = symmetrize(B);
    AA = symmetrize(A);

    svd(BB, U, RSinvP, Udag); // bb is symmetric, Udag is a placeholder here

    Udag = adj(U);   //true Udag

    int rdim = svdResetAverageCond(RSinvP, condMax);

    eVecs.reDim(bin_, row_, rdim);
    eVecs.zeros();
    wVecs.reDim(bin_, rdim, rdim);
    wVecs.zeros();
    eVals.reDim(bin_, rdim);
    eVals.zeros();

    svdResetPseudoInvertRoot(RSinvP, rdim);

    RSIP = diag(RSinvP);

    RSIP.rows(rdim);

    M = symmetrize(RSIP * Udag * AA * U * adj(RSIP));

    eig_sym(M, eVals, wVecs);

    if(sortByVecs)
      matchEigenVectorsEnsemble(wVecs, eVals);
    else
      matchEigenValueSizeEnsemble(wVecs, eVals);

    eVecs = U * adj(RSIP) * wVecs;

    //enforce a phase convention
    //rephaseEigenVectors(eVecs);
  }

  template<class T>
  void genEigSvdValue(const SembleMatrix<T> &A, const SembleMatrix<T> &B, SembleMatrix<T> &eVecs, SembleMatrix<T> &wVecs,
                      SembleVector<double> &eVals, const bool sortByVecs = true, const double condMax = 1e7)
  {
    int bin_ = A.bins(), row_ = A.rows(), col_ = A.cols();

    if((bin_ != B.bins()) || (row_ != B.cols()) || (row_ != B.rows()) || (row_ != col_))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }


    SembleMatrix<T> U, Udag, M, RSIP, AA, BB;
    SembleVector<double> RSinvP;

    BB = symmetrize(B);
    AA = symmetrize(A);

    svd(BB, U, RSinvP, Udag); // bb is symmetric, Udag is a placeholder here

    Udag = adj(U);   //true Udag

    int rdim = svdResetAverageValue(RSinvP, condMax);       //***this is the only different line from above

    eVecs.reDim(bin_, row_, rdim);
    eVecs.zeros();
    wVecs.reDim(bin_, rdim, rdim);
    wVecs.zeros();
    eVals.reDim(bin_, rdim);
    eVals.zeros();

    svdResetPseudoInvertRoot(RSinvP, rdim);

    RSIP = diag(RSinvP);

    RSIP.rows(rdim);

    M = symmetrize(RSIP * Udag * AA * U * adj(RSIP));

    eig_sym(M, eVals, wVecs);

    if(sortByVecs)
      matchEigenVectorsEnsemble(wVecs, eVals);
    else
      matchEigenValueSizeEnsemble(wVecs, eVals);

    eVecs = U * adj(RSIP) * wVecs;

    //enforce a phase convention
    //rephaseEigenVectors(eVecs);
  }

  template<class T> //nb only for hermitian matricies
  void genEigSvdF(const SembleMatrix<T> &A, const SembleMatrix<double> &rootSinvPlus, const SembleMatrix<T> &U, SembleMatrix<T> &eVecs,
                  SembleMatrix<T> &wVecs, SembleVector<double> &eVals, const bool sortByVecs = true)
  {
    int bin_ = A.bins(), row_ = A.rows(), col_ = A.cols();

    if(((bin_ != U.bins()) || (row_ != U.cols()) || (row_ != U.rows()) || (row_ != col_)) ||
        (rootSinvPlus.rows() > rootSinvPlus.cols()) || (rootSinvPlus.cols() != row_))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> M, AA;

    int ldim = rootSinvPlus.rows();

    eVecs.reDim(bin_, row_, ldim);
    eVecs.zeros();
    wVecs.reDim(bin_, ldim, ldim);
    wVecs.zeros();
    eVals.reDim(bin_, ldim);
    eVals.zeros();

    AA = symmetrize(A);
    M = symmetrize(rootSinvPlus * adj(U) * AA * U * adj(rootSinvPlus));
    eig_sym(M, eVals, wVecs);

    if(sortByVecs)
      matchEigenVectorsEnsemble(wVecs, eVals);
    else
      matchEigenValueSizeEnsemble(wVecs, eVals);

    eVecs = U * adj(rootSinvPlus) * wVecs;

    //enforce a phase convention
    //rephaseEigenVectors(eVecs);

  }

  template<class T> //nb only for hermitian matricies
  void genEigSvdFRef(const SembleMatrix<T> &A, const SembleMatrix<T> &rootSinvPlus, const SembleMatrix<T> &U, const SembleMatrix<T> &refW,
                     SembleMatrix<T> &eVecs, SembleMatrix<T> &wVecs, SembleVector<double> &eVals, const bool sortByVecs = true)
  {
    int bin_ = A.bins(), row_ = A.rows(), col_ = A.cols();

    if(((bin_ != U.bins()) || (row_ != U.cols()) || (row_ != U.rows()) || (row_ != col_)) ||
        (rootSinvPlus.rows() > rootSinvPlus.cols()) || (rootSinvPlus.cols() != row_))
      {
        std::cout << "Dimensionally Incorrect in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
        exit(1);
      }

    SembleMatrix<T> M, AA;

    int ldim = rootSinvPlus.rows();

    eVecs.reDim(bin_, row_, ldim);
    eVecs.zeros();
    wVecs.reDim(bin_, ldim, ldim);
    wVecs.zeros();
    eVals.reDim(bin_, ldim);
    eVals.zeros();

    AA = symmetrize(A);
    M = symmetrize(rootSinvPlus * adj(U) * AA * U * adj(rootSinvPlus));

    eig_sym(M, eVals, wVecs);

    if(sortByVecs)
      matchEigenVectorsEnsemble(wVecs, eVals);
    else
      matchEigenValueSizeEnsemble(wVecs, eVals);

    //this orders the ensemble
    matchEigenVectorsEnsemble(refW, wVecs, eVals);
 
    eVecs = U * adj(rootSinvPlus) * wVecs;

    //enforce a phase convention
    //rephaseEigenVectors(eVecs);

  }

//helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  template<class T>
  int svdResetSigma(SembleVector<double> &sing, SembleMatrix<T> &U, const double f = 3, const bool r = false)
  {
    int dim = sing.getN();

    if((dim != U.getN()) || (dim != U.getM()))
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " only implemented for square matricies, exiting.." << std::endl;
        exit(1);
      }

    if(r)
      sing.rescaleEnsemUp();

    itpp::Vec<double> m(mean(sing)), v(variance(sing));
    std::vector<std::pair<int, int> > vmap;
    std::vector<int> resets;

    int val(0);

    for(int index = 0; index < dim; ++index)
      {
        if(m[index] - std::sqrt(v[index])*f > 0)
          {
            vmap.push_back(std::pair<int, int>(val, index));
            ++val;
          }
        else
          resets.push_back(index);
      }

    //now lets reorder it so that we only have the things that are not compatible with zero on first block and everythign else on second
    SembleVector<double> cps(sing);
    SembleMatrix<T> cpU(U);

    std::vector<std::pair<int, int> >::const_iterator it;

    for(it = vmap.begin(); it != vmap.end(); ++it)
      {
        if(it->first == it->second)
          continue;

        sing.loadEnsemElement(it->first, cps.getEnsemElement(it->second));

        for(int row = 0; row < dim; ++row)
          U.loadEnsemElement(row, it->first, cpU.getEnsemElement(row, it->second));
      }

    //now put the rest of them in
    int off = vmap.size();

    for(int i = 0; i < dim - off; ++i)
      {
        sing.loadEnsemElement(off + i, cps.getEnsemElement(resets[i]));

        for(int row = 0; row < dim; ++row)
          U.loadEnsemElement(row, i + off, cpU.getEnsemElement(row, resets[i]));
      }

    if(r)
      sing.rescaleEnsemDown();

    return off; // this is the dimension of the good space
  }

  template<class T>
  int svdResetAvgValueAndSigma(SembleVector<double> &sing, SembleMatrix<T> &U, const double thresh = 1e-6, const double sigma = 3, const bool scaleUp = false)
  {
    //get rid of anything with variance overlapping onto zero
    int rdim = svdResetSigma(sing, U, sigma, scaleUp);

    itpp::Vec<double> s = mean(sing);

    //get rid of anything below the value cut
    while(true)
      {
        if(s(rdim - 1) > thresh)
          break;

        --rdim;

        if(rdim == 0)
          {
            std::cout << "All Null Space in" << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
            exit(1);
          }
      }

    return rdim;
  }

  template<class T>
  int svdResetAvgCondAndSigma(SembleVector<double> &sing, SembleMatrix<T> &U, const double thresh = 1e7, const double sigma = 3, const bool scaleUp = false)
  {
    //get rid of anything compatible with zero
    int rdim = svdResetSigma(sing, U, sigma, scaleUp);

    itpp::Vec<double> s = mean(sing);

    while(true)
      {
        if(s(0) / s(rdim - 1) < thresh)
          break;

        --rdim;

        if(rdim == 0)
          {
            std::cout << "All Null Space in" << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
            exit(1);
          }
      }

    return rdim;
  }

  template<class T>//test to see if the svd algorithm moved the singular values around, move them back
  std::string svdMatchSingValsByUandRephaseEnsemble(const itpp::Mat<T> &ref, SembleMatrix<T> &U, SembleVector<double> &S, const bool rescale = true)
  {
    if((U.getN() != S.getN()) || (U.getB() != S.getB()))
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Dimensionally incorrect" << std::endl;
        exit(1);
      }

    if(rescale)
      {
        U.rescaleEnsemDown();
        S.rescaleEnsemDown();
      }

    int sz = U.getB();
    int dim = U.getN();
    std::stringstream ss;
    bool swapped = false;
    bool phase = false;

    for(int i = 0; i < sz; ++i)
      {
        itpp::Vec<double> cpS = S[i];
        itpp::Mat<T> cpU = U[i];
        matchEigenVectors(ref, U[i], S[i]);

        if(cpS != S[i])
          {
            swapped = true;

            ss << "Bin[" << i << "]\n";

            std::vector<int> match;

            for(int elem = 0; elem < dim; ++elem)
              if(cpS[elem] != S[i][elem])
                match.push_back(elem);

            std::vector<int> reff(match);
            std::vector<int>::const_iterator itr;
            std::vector<int>::iterator itm;

            for(itr = reff.begin(); itr != reff.end(); ++itr)
              {
                for(itm = match.begin(); itm != match.end(); ++itm)
                  if(cpS[*itr] == S[i][*itm])
                    break;

                ss << *itr << " maps to " << *itm << "\n";
                match.erase(itm);
              }

            if(!!!match.empty())
              {
                ss << "Something ill-defined happened, matched has " << match.size() << "elements..\n";

                for(itm = match.begin(); itm != match.end(); ++itm)
                  ss << *itm << "\n";
              }

            ss << "\n\n";
          }
        else if(cpU != U[i])
          phase = true;
      }

    if(rescale)
      {
        U.rescaleEnsemUp();
        S.rescaleEnsemUp();
      }

    if(phase)
      ss << "\n there was a phasing ambuguity \n";

    if(swapped || phase)
      return ss.str();

    return std::string("No ordering/phasing ambiguity\n");
  }


}//namespace


#endif
