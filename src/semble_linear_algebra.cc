// semble_linear_algebra.cc -
//
// Thursday, November 17 2011
//

#include"semble_linear_algebra.h"

namespace SEMBLE
{

  PromoteScalar<double>::Type eConjPhase(const double &in)
  {
    double dum = double((in > 0) - (in < 0))*in;
    return toScalar(dum);
  }

  PromoteScalar<std::complex<double> >::Type eConjPhase(const std::complex<double> &in)
  {
    return toScalar(std::conj(in));
  }

//overloaded, not templated
  void rephaseEVectors(SembleMatrix<double> &vecs)
  {
    int n = vecs.rows(), m = vecs.cols();
    double max;
    std::map<int,bool> rephase;
    itpp::Mat<double> mvec = mean(vecs);

    for(int vec = 0; vec < m; ++vec)
      {
	max = 0;
	itpp::Vec<double> dum = mvec.get_col(vec);

	for(int elem = 0; elem < n; ++elem)
	  if(fabs(dum(elem)) > fabs(max))
	    max = dum(elem);

	if(max < 0)
	  rephase.insert(std::pair<int,bool>(vec,true));
	else
	  rephase.insert(std::pair<int,bool>(vec,false));
      }
    
    std::map<int,bool>::const_iterator it;
    
    for(it = rephase.begin(); it != rephase.end(); ++it)
      {
	if(it->second)
	  for(int row = 0; row < n; ++row)
	    vecs.loadEnsemElement(row,it->first,toScalar(-1.0)*vecs.getEnsemElement(row,it->first));
      }
  }

  void rephaseEVectors(SembleMatrix<std::complex<double> > &vecs)
  {
    std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " not implemented, exiting" << std::endl;
    exit(1);
  }


  void pseudoInvert(SembleVector<double> &inout, const int dim_, bool rescale/*=true*/)
  {
    int bin_ = inout.getB();
    int bdim_ = inout.getN();

    if(rescale)
      inout.rescaleEnsemDown();

    for(int bin = 0; bin < bin_; ++bin)
      {
        for(int elem = 0; elem < bdim_; ++elem)
          {
            if(elem < dim_)
              inout[bin][elem] = 1. / inout[bin][elem];
            else
              inout[bin][elem] = 0.;
          }
      }

    if(rescale)
      inout.rescaleEnsemUp();
  }

//reset index greater than/eq dim and
  void svdResetPseudoInvertRoot(SembleVector<double> &inout, const int dim_, const bool rescale/*=true*/)
  {
    int bin_ = inout.getB();
    int bdim_ = inout.getN();

    if(rescale)
      inout.rescaleEnsemDown();

    for(int bin = 0; bin < bin_; ++bin)
      {
        for(int elem = 0; elem < bdim_; ++elem)
          {
            if(elem < dim_)
              inout[bin][elem] = 1. / std::sqrt(inout[bin][elem]);
            else
              inout[bin][elem] = 0.;
          }
      }

    if(rescale)
      inout.rescaleEnsemUp();
  }

  int svdResetAverageValue(const SembleVector<double> &in, const double thresh/*= 1e-6*/)
  {
    itpp::Vec<double> s = mean(in);
    int rdim = in.getN();

    while(true)
      {
        if(s(rdim - 1) > thresh)
          break;

        --rdim;

        if(rdim == 0)
          {
            std::cout << "All Null Space in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
            exit(1);
          }
      }

    return rdim;
  }

  int svdResetAverageCond(const SembleVector<double> &inout, const double condMax/*= 1e7*/)
  {
    itpp::Vec<double> s = mean(inout);
    int rdim = inout.getN();

    while(true)
      {
        if(s(0) / s(rdim - 1) < condMax) //nb this will obviously blow up on zero or negative (svd has positive sing vals..) condition numbers..
          break;

        --rdim;

        if(rdim == 0)
          {
            std::cout << "All Null Space in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
            exit(1);
          }
      }

    return rdim;
  }

}//namespace
