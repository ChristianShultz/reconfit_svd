#include "linear_algebra.h"

itpp::mat pseudoinvertSVD(const itpp::mat& A, double tol, int& nReset)
{
  int M = A.rows();
  int N = A.cols();
  
  itpp::mat U(M,M);
  itpp::mat V(N,N);
  itpp::vec s(M);
  
  itpp::svd(A, U, s, V);
  
  //  cout << "done the svd" << endl;

  double smax = itpp::max(s); 
  int count = 0;
  
  itpp::mat Sinv(N,M);
  Sinv.zeros();
  for(int i = 0; i < N; i++)
    { //test for small singular values - kills those smaller than tol of maximum
      if(s(i) > tol * smax){ Sinv(i,i) = 1./s(i);}
      else{Sinv(i,i) = 0; count++;}
    }
  
  nReset = count;
  // cout << "reset " << count << " sing values" << endl;
  itpp::mat inverse = V*Sinv*transpose(U);  
  return inverse;
}

itpp::cmat pseudoinvertSVD(const itpp::cmat& A, double tol, int& nReset)
{
  int M = A.rows();
  int N = A.cols();
  
  itpp::cmat U(M,M);
  itpp::cmat V(N,N);
  itpp::vec s(M);
  
  itpp::svd(A, U, s, V);
	
  double smax = itpp::max(s); 
  int count = 0;
  
  itpp::mat Sinv(N,M);
  Sinv.zeros();
  for(int i = 0; i < N; i++)
    { //test for small singular values - kills those smaller than tol of maximum
      if(s(i) > tol * smax){ Sinv(i,i) = 1./s(i);}
      else{Sinv(i,i) = 0; count++;}
    }
  
  nReset = count;
  // cout << "reset " << count << " sing values" << endl;
  
  itpp::cmat inverse = V*Sinv*itpp::hermitian_transpose(U);
  
  return inverse;
}


itpp::mat invertSVDNorm(const itpp::mat& A, double tol, int& nReset)
{
  //normalise the matrix to O(1) using the sqrt of the diag
  //only for square matrices now
  if(A.rows() != A.cols()){cerr << "square matrices only in invertSVDNorm" << endl; exit(1);}
  
  itpp::mat D(A.rows(), A.cols());
  D.zeros();
  for(int i=0; i < A.rows(); i++){D(i,i) = 1.0 / sqrt(fabs(A(i,i)));}
  
  itpp::mat DAD = D*A*D;
  
  // cout << "normalised the covariance" << endl;

  itpp::mat AA = pseudoinvertSVD(DAD, tol, nReset);
  
  return D*AA*D;	
}

itpp::cmat invertSVDNorm(const itpp::cmat& A, double tol, int& nReset)
{
  //normalise the matrix to O(1) using the sqrt of the diag
  //only for square matrices now
  if(A.rows() != A.cols()){cerr << "square matrices only in invertSVDNorm" << endl; exit(1);}
  
  itpp::cmat D(A.rows(), A.cols());
  D.zeros();
  for(int i=0; i < A.rows(); i++){D(i,i) = 1.0 / sqrt(abs(A(i,i)));}
  
  itpp::cmat DAD = D*A*D;
  
  itpp::cmat AA = pseudoinvertSVD(DAD, tol, nReset);
  
  return D*AA*D;	
}


void sortEigen(itpp::vec& evals, itpp::mat& evecs, bool asc)
{
  itpp::ivec indexes = itpp::sort_index(evals);
  
  int j=0;
  itpp::vec cpevals(evals.size());
  for(int i=0; i < evals.size(); i++){cpevals(i) = evals(i);}
  
  itpp::mat cpevecs = evecs;
  if(asc){
    //ascending
    for(int i=0;i<evals.size();i++){
      if(indexes[i]!=i){
	evals[i]=cpevals[indexes[i]];
	evecs.set_col(i,cpevecs.get_col(indexes[i]));	
      }
    }
  }
  else{
    //descending
    for(int i=evals.size()-1; i>-1;i--){
      if(indexes[i]!=j){
	evals[j]=cpevals[indexes[i]];
	evecs.set_col(j,cpevecs.get_col(indexes[i]));	
      }
      j++;
    }
  }
  
}

void sortEigen(itpp::vec& evals, itpp::cmat& evecs, bool asc)
{
  itpp::ivec indexes = itpp::sort_index(evals);
  
  int j=0;
  itpp::vec cpevals(evals.size());
  for(int i=0; i < evals.size(); i++){cpevals(i) = evals(i);}
  
  itpp::cmat cpevecs = evecs;
  if(asc){
    //ascending
    for(int i=0;i<evals.size();i++){
      if(indexes[i]!=i){
	evals[i]=cpevals[indexes[i]];
	evecs.set_col(i,cpevecs.get_col(indexes[i]));	
      }
    }
  }
  else{
    //descending
    for(int i=evals.size()-1; i>-1;i--){
      if(indexes[i]!=j){
	evals[j]=cpevals[indexes[i]];
	evecs.set_col(j,cpevecs.get_col(indexes[i]));	
      }
      j++;
    }
  }
  
}

void sortEigen(itpp::cvec& evals, itpp::cmat& evecs, bool asc)
{
  itpp::vec mods = abs(evals);
  
  itpp::ivec indexes = itpp::sort_index(mods);
  
  int j=0;
  itpp::cvec cpevals(evals.size());
  for(int i=0; i < evals.size(); i++){cpevals(i) = evals(i);}
  
  itpp::cmat cpevecs = evecs;
  if(asc){
    //ascending
    for(int i=0;i<evals.size();i++){
      if(indexes[i]!=i){
	evals[i]=cpevals[indexes[i]];
	evecs.set_col(i,cpevecs.get_col(indexes[i]));	
      }
    }
  }
  else{
    //descending
    for(int i=evals.size()-1; i>-1;i--){
      if(indexes[i]!=j){
	evals[j]=cpevals[indexes[i]];
	evecs.set_col(j,cpevecs.get_col(indexes[i]));	
      }
      j++;
    }
  }
  
}


bool sortEigen(itpp::vec& evals, itpp::mat& evecs, itpp::mat RefVecs)
{
  // Sort eigenvalues and eigenvectors by comparing with RefVecs
  
  int dim = evecs.cols();
  itpp::vec evals_old = evals;
  itpp::mat evecs_old = evecs;
  
  //assert ( (evecs.cols() == evecs.rows()) && (evals.rows() == evecs.rows()) );
  //assert ( (RefVecs.cols() == RefVecs.row()) && (RefVecs.rows() == evecs.rows()) );

  if(!(evecs.cols() == evecs.rows()) && (evals.size() == evecs.rows()) ){cerr << "evecs & evals don't match" << endl; exit(1);}
  if(!(RefVecs.cols() == RefVecs.rows()) && (RefVecs.rows() == evecs.rows()) ){cerr << "Refvecs & evecs don't match" << endl; exit(1);}

  
  // Find the ref vector that maximizes overlap (dot product) with eigenvector
  itpp::mat Overlaps = transpose(evecs_old) * RefVecs;
  vector<int> NewIndex(dim);
  vector<int> UsedEvecs(dim); for(int i=0; i<dim; i++) { UsedEvecs[i] = 0; }
  vector<int> UsedRefVecs(dim); for(int i=0; i<dim; i++) { UsedRefVecs[i] = 0; }
  
  for(int i = 0; i < dim; i++)
    {
      double MaxOverlap = 0.0;
      int MaxOverlapEvec = 0;
      int MaxOverlapRefVec = 0;
      // Find maximum overlap in matrix
      for(int j=0; j < dim; j++)
	{
	  for(int k=0; k < dim; k++)
	    {
	      double Overlap = Overlaps(j,k);
	      if ( (abs(Overlap) > abs(MaxOverlap)) && (UsedEvecs[j] == 0) && (UsedRefVecs[k] == 0) )
		{
		  MaxOverlap = Overlap;
		  MaxOverlapEvec = j;
		  MaxOverlapRefVec = k;
		}
	    }
	}
      NewIndex[MaxOverlapEvec] = MaxOverlapRefVec;
      UsedEvecs[MaxOverlapEvec] += 1;
      UsedRefVecs[MaxOverlapRefVec] += 1;
      if (abs(MaxOverlap) < 0.7)
	{
	  //	  ofstream reorderlog("reorder_log", ios::app);
	  // reorderlog << __func__ << ": WARNING: reasonably low MaxOverlap = " << MaxOverlap << endl;
	  // reorderlog.close();
	}
    }
  
  bool testflag = true;
  for(int i = 0; i < dim; i++)
    {
      if (NewIndex[i] != i)
	{
	  testflag = false;
	  break;
	}
    }
  
  // Reorder the eigenvectors and eigenvalues
  for(int state = 0; state < dim; state++)
    {
      evecs.set_col(NewIndex[state],evecs_old.get_col(state));
      evals[NewIndex[state]] = evals_old[state];
    }
  
  return testflag;
}

//resigns each eigenvector by forcing a positive overlap with RefVecs
void resignEigenByOverlap(itpp::mat& evecs, itpp::mat& RefVecs)
{
  //  assert ( (evecs.rows()==evecs.cols()) && (evecs.rows()==RefVecs.rows()) && (RefVecs.rows()==RefVecs.cols()) );
  if(!(evecs.cols() == evecs.rows()) && (evecs.rows() == RefVecs.rows()) ){cerr << "evecs & Refvecs don't match" << endl; exit(1);}

  int dim = evecs.rows();
  itpp::mat Overlaps = transpose(evecs) * RefVecs;
  
  for(int evecNum = 0; evecNum < dim; evecNum++)
    {
      double Overlap = Overlaps(evecNum, evecNum);
      double OverlapSign = Overlap / abs(Overlap);
      for(int i=0; i < dim; i++){evecs(i, evecNum) = evecs(i, evecNum) * OverlapSign;}
    }//next eigenvector
}

//resigns each eigenvector so that its largest element (by modulus) is positive
void resignEigen(itpp::mat& evecs)
{
  int dim = evecs.rows();
  for(int evecNum = 0; evecNum < dim; evecNum++)
    {
      itpp::vec evectorMod(dim);
      itpp::vec evectorSign(dim);
      for(int i=0; i < dim; i++){   
	double m = evecs(i, evecNum);
	evectorMod(i) = abs(m);
	evectorSign(i) = m / abs(m);
      }
      
      double signFactor = evectorSign(itpp::max_index(evectorMod)) ;
      for(int i=0; i < dim; i++){evecs(i, evecNum) = evecs(i, evecNum) * signFactor;}
    }//next eigenvector
}

//rephases each eigenvector so that its largest element (by modulus) is real and positive
void rephaseEigen(itpp::cmat& evecs)
{
  int dim = evecs.rows();
  for(int evecNum = 0; evecNum < dim; evecNum++)
    {
      itpp::vec evectorMod(dim);
      itpp::cvec evectorPhase(dim);
      for(int i=0; i < dim; i++){   
	complex<double> m = evecs(i, evecNum);
	evectorMod(i) = abs(m);
	evectorPhase(i) = m / abs(m);
      }
      
      complex<double> phaseFactor = conj( evectorPhase(itpp::max_index(evectorMod))  );
      for(int i=0; i < dim; i++){evecs(i, evecNum) = evecs(i, evecNum) * phaseFactor;}
    }//next eigenvector
}


double statQ(double chisq, int nDoF)
{
  NR::DP c = chisq / 2.0;
  NR::DP n = double(nDoF) / 2.0;
  
 if(n > 0.0){ return double(NR::gammq(n,c)); }
 else{ return 0.0;}	
}
