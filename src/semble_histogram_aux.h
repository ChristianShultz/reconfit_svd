#ifndef SEMBLE_HISTOGRAM_AUX_H_H_GUARD
#define SEMBLE_HISTOGRAM_AUX_H_H_GUARD

#include <utility>
#include <vector>
#include "semble_vector.h"

namespace SEMBLE
{

  template<class T>
  std::pair<std::vector<T>, std::vector<T> > findRange(const SembleVector<T> &v, int freq = 5)
  {
    std::vector<T> st, end;
    int sz = v.getB();
    int nelem = v.getN();

    if(sz < freq)
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " unable to histogram, exiting" << std::endl;
        exit(1);
      }

    st.resize(nelem);
    end.resize(nelem);

    for(int elem = 0; elem < nelem; ++elem)
      {
        T a = v[0][elem];
        T b = v[freq - 1][elem];

        if(a > b)
          {
            T swap = a;
            a = b;
            b = swap;
          }

        st[elem] = a;
        end[elem] = b;
      }

    for(int bin = freq; bin < sz; bin += freq)
      {
        const itpp::Vec<T> dum = v[bin];

        for(int elem = 0; elem < nelem; ++elem)
          {
            if(dum[elem] < st[elem])
              st[elem] = dum[elem];
            else if(end[elem] < dum[elem])
              end[elem] = dum[elem];
          }
      }

    return std::pair<std::vector<T>, std::vector<T> >(st, end);
  }


}
#endif
