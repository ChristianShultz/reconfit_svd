#ifndef __LINEAR_ALGEBRA_H__
#define __LINEAR_ALGEBRA_H__

#include <vector>
#include <math.h>
#include <iostream>
#include <itpp/itbase.h>
#include "recipes/nr.h"

using namespace std;


itpp::mat pseudoinvertSVD(const itpp::mat& A, double tol, int& nReset);
itpp::cmat pseudoinvertSVD(const itpp::cmat& A, double tol, int& nReset);
itpp::mat invertSVDNorm(const itpp::mat& A, double tol, int& nReset);
itpp::cmat invertSVDNorm(const itpp::cmat& A, double tol, int& nReset);

void sortEigen(itpp::vec& evals, itpp::mat& evecs, bool asc);
void sortEigen(itpp::vec& evals, itpp::cmat& evecs, bool asc);
void sortEigen(itpp::cvec& evals, itpp::cmat& evecs, bool asc);

//void sortEigen(itpp::vec& evals, itpp::cmat& evecs, itpp::mat RefVecs);

bool sortEigen(itpp::vec& evals, itpp::mat& evecs, itpp::mat RefVecs);

void resignEigen(itpp::mat& evecs);
void resignEigenByOverlap(itpp::mat& evecs, itpp::mat& RefVecs);
void rephaseEigen(itpp::cmat& evecs);

double statQ(double chisq, int nDoF); // the statistical Q function


#endif
