#ifndef SEMBLE_FUNCTIONS_H_H_GUARD
#define SEMBLE_FUNCTIONS_H_H_GUARD

#include"semble_matrix.h"
#include"ensem_data.h"
#include<vector>

namespace SEMBLE
{
//flatten a vector of ensemble matrices --- ONLY SQUARE
  EnsemData flatten(const std::vector<SembleMatrix<double> > &M)
  {
    EnsemData flat;
    int Lt = M.size();
    int bins = M[0].getB();
    int dim = M[0].getN();

    for(int i = 0; i < dim; i++)
      for(int j = 0; j < dim; j++)
        {
          int k = i * dim + j; //element location - read right to left from top left

          for(int t = 0; t < Lt; t++)
            {
              double x = k * Lt + t; //timeslice location
              EnsemReal y = M[t].getEnsemElement(i, j);
              flat.addDatum(x, y);
            }
        }

    return flat;
  }

  EnsemData flattenSym(const std::vector<SembleMatrix<double> > &M)
  {
    EnsemData flat;
    int Lt = M.size();
    int bins = M[0].getB();
    int dim = M[0].getN();

    for(int i = 0; i < dim; i++)
      for(int j = i; j < dim; j++)
        {
          int k = i * dim + j; //element location - read right to left from top left

          for(int t = 0; t < Lt; t++)
            {
              double x = k * Lt + t; //timeslice location
              EnsemReal y = (M[t].getEnsemElement(i, j) + M[t].getEnsemElement(j, i)) * Real(0.5);
              flat.addDatum(x, y);
            }
        }

    return flat;
  }
}

#endif
