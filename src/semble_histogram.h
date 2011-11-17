#ifndef SEMBLE_HISTOGRAM_H_H_GUARD
#define SEMBLE_HISTOGRAM_H_H_GUARD

#include <vector>
#include <string>
#include "itpp/itbase.h"
#include "semble_vector.h"


namespace SEMBLE
{

  /*
     A simple histogram class to keep a tally of values withing a range
     specified during construction.  The range is arranged into some number of
     bins specified during construction. bin zero is underflow, bin nBins+1 is overflow.
   */
  struct SembleHisto
  {
    //constructors, destructors, copy assignment
    SembleHisto(const double st, const double end, const int nbins);
    SembleHisto(const SembleHisto &o);
    SembleHisto &operator=(const SembleHisto &o);
    ~SembleHisto(void) {}

    //add a point
    inline
    void Add(const double x);

    //add a bunch of points
    void Add(const std::vector<double> &xx);

    //get the total counts
    int numCounts(void);

    //clear the total counts
    void clearCounts(void);

    //generate a 'histogram'.. really a string with the middle of the bin and the num of counts on each line
    std::string genHisto(void) const;

  private: //data store
    double bins_by_interval, start;
    std::vector<int> freq;               //underflow is bin zero, overflow is bin nBins+1
    int nBins;
  };

  /*
    An intermediary class to histogram a SembleVector by element.
   */
  struct SembleMultiHisto
  {
    //constructors, destructor, copy assignment
    SembleMultiHisto(const std::vector<double> &starts, const std::vector<double> &ends, const std::vector<int> &nbins);
    SembleMultiHisto(const SembleMultiHisto &o);
    SembleMultiHisto &operator=(const SembleMultiHisto &o);
    ~SembleMultiHisto(void) {}

    //add, will check that size of xx is size of hists, histograms the ith element of xx into the ith histogram of hists
    inline
    void Add(const itpp::Vec<double> &xx);

    //same but don't check sizes.. no extra if statement
    inline
    void fastAdd(const itpp::Vec<double> &xx);

    //add the entire histogram, no checks
    void Add(const SembleVector<double> &xxx);

    std::vector<std::string> genHisto(void) const;

  private: //data store
    int sz;
    std::vector<SembleHisto> hists;

  };


}//end namespace

#endif
