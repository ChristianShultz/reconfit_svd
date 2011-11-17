// semble_linear_algebra.cc -
//
// Thursday, November 17 2011
//

#include"semble_linear_algebra.h"

namespace SEMBLE
{


//overloaded, not templated
  void rephaseEVectors(SembleMatrix<double> &vecs)
  {
    int b = vecs.bins(), n = vecs.rows(), m = vecs.cols();
    bool rephase;
    double max;

    for(int bin = 0; bin < b; ++bin)
      {

        for(int vec = 0; vec < m; ++vec)
          {
            rephase = false;
            max = 0.;

            for(int elem = 0; elem < n; ++elem)
              {
                if(fabs(vecs[bin](elem, vec)) > max)
                  {
                    max = vecs[bin](elem, vec);

                    if(max < 0.)
                      {
                        max *= -1.;
                        rephase = true;
                      }
                    else
                      rephase = false;
                  }

                if(rephase)
                  (vecs[bin]).set_col(vec, -(vecs[bin]).get_col(vec));
              }
          }
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
