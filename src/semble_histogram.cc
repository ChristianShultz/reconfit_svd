// semble_histogram.cc -
//
// Wednesday, October 26 2011
//




#include"semble_histogram.h"
#include <sstream>

namespace SEMBLE
{

//SembleHisto Implementation
///////////////////////////////////////////////////////////////////////////

  SembleHisto::SembleHisto(const double st, const double end, const int nbins)
    : nBins(nbins), start(st)
  {
    freq.resize(nBins + 2);
    bins_by_interval = double(nBins) / (end - start);
  }

  SembleHisto::SembleHisto(const SembleHisto &o)
    : bins_by_interval(o.bins_by_interval), start(o.start), freq(o.freq), nBins(o.nBins)
  {}

  SembleHisto &SembleHisto::operator=(const SembleHisto &o)
  {
    if(this != &o)
      {
        bins_by_interval = o.bins_by_interval;
        start = o.start;
        freq = o.freq;
        nBins = o.nBins;
      }

    return *this;
  }

  void SembleHisto::Add(const double x)
  {
    if(x - start < 0)
      ++freq[0];
    else
      {
        int i = int((x - start) * bins_by_interval);

        if(i > nBins)
          ++freq[nBins + 1];
        else
          ++freq[i];
      }
  }

  void SembleHisto::Add(const std::vector<double> &xx)
  {
    std::vector<double>::const_iterator it;

    for(it = xx.begin(); it != xx.end(); ++it)
      Add(*it);
  }

  int SembleHisto::numCounts(void)
  {
    int ct = 0;
    std::vector<int>::const_iterator it;

    for(it = freq.begin(); it != freq.end(); ++it)
      ct += *it;

    return ct;
  }

  void SembleHisto::clearCounts(void)
  {
    freq.clear();
  }

  std::string SembleHisto::genHisto(void) const
  {
    std::stringstream ss;
    double bin_len = 1. / bins_by_interval;
    double off = bin_len * 0.5;

    for(int i = 0; i < nBins + 2; ++i)
      ss << i *bin_len - off << " " << freq[i] << "\n";

    return ss.str();
  }

//End SembleHisto Implementation
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//SemblMultiHisto Implementation
////////////////////////////////////////////////////////////////////////
  SembleMultiHisto::SembleMultiHisto(const std::vector<double> &s, const std::vector<double> &e, const std::vector<int> &n)
  {
    if((s.size() != e.size()) || (e.size() != n.size()))
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Dimension mismatch, exiting.." << std::endl;
        exit(1);
      }

    sz = s.size();

    for(int i = 0; i < sz; ++i)
      hists.push_back(SembleHisto(s[i], e[i], n[i]));
  }

  SembleMultiHisto::SembleMultiHisto(const SembleMultiHisto &o)
    : hists(o.hists), sz(o.sz)
  {
  }

  SembleMultiHisto &SembleMultiHisto::operator=(const SembleMultiHisto &o)
  {
    if(this != &o)
      {
        hists = o.hists;
        sz = o.sz;
      }

    return *this;
  }

  void SembleMultiHisto::Add(const itpp::Vec<double> &xx)
  {
    if(xx.size() != sz)
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Dimension mismatch, exiting.." << std::endl;
        exit(1);
      }

    for(int i = 0; i < sz; ++i)
      hists[i].Add(xx[i]);
  }

  void SembleMultiHisto::fastAdd(const itpp::Vec<double> &xx)
  {
    for(int i = 0; i < sz; ++i)
      hists[i].Add(xx[i]);
  }

  void SembleMultiHisto::Add(const SembleVector<double> &xxx)
  {
    int sz = xxx.getB();
    int nelems = xxx.getN();

    if(nelems != hists.size())
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << " Dimension mismatch, exiting.." << std::endl;
        exit(1);
      }

    for(int i = 0; i < sz; ++ i)
      {
        itpp::Vec<double> dum = xxx[i];

        for(int elem = 0; elem < nelems; ++elem)
          hists[elem].Add(dum[elem]);
      }
  }

  std::vector<std::string> SembleMultiHisto::genHisto(void) const
  {
    std::vector<std::string> histograms;
    std::vector<SembleHisto>::const_iterator it;

    for(it = hists.begin(); it != hists.end(); ++it)
      histograms.push_back(it->genHisto());

    return histograms;
  }

}
