#include "jackknife_linear_fitter.h"

using namespace SEMBLE;

JackLinearFit::JackLinearFit(EnsemData data, vector< Handle<EnsemFunction> > basis_f, double sv_cutoff){

  //get the data x,y-values
  int nData = data.getNData();
  vector<double> x = data.getXData();
  EnsemVectorReal y = data.getYData();
  EnsemVecReal dat; dat = y;

  //make the basis
  int nPars = basis_f.size();
  EnsemMatrixReal basis(data.getNBins(), nData, nPars);

  for(int par=0; par < nPars; par++){
    for(int ix=0; ix < nData; ix++){
      basis.loadEnsemElement( ix , par,  basis_f[par]->operator()(x[ix])   );
    }
  }

  //get the data covariance
  itpp::mat inv_data_cov = data.getInvCov();
  int nResetSVDataCov = data.getNResetCovSingVals();
  if(nResetSVDataCov > 0){ cout << "inverse data covariance requires " << nResetSVDataCov << " reset singular values" << endl;}

  //form the design matrix
  EnsemMatrixReal A = transpose(basis)* (inv_data_cov * basis);
  //for the RHS
  EnsemVecReal b = (transpose(basis)*inv_data_cov) * dat;

  //solve the linear system
  double res; EnsemVecReal xx; EnsemMatrixReal cov_xx; int nResetSVLinSys;
  solveLinearSVD(A, xx, b, cov_xx, sv_cutoff, res, nResetSVLinSys);

  if(nResetSVLinSys > 0){ cout << "linear system inverted with " << nResetSVLinSys << " reset singular values" << endl;}

  //compute the chisq
  EnsemReal dum1 = dot(dat , (inv_data_cov * dat) );
  EnsemReal dum2 = dot(xx, b);
  EnsemReal ensemChisq = dum1 - dum2;
  
  chisq = toDouble(mean(ensemChisq));
  nDoF = nData - nPars - nResetSVDataCov - nResetSVLinSys;
  Q = statQ(chisq, nDoF);
  
  for(int par = 0; par < nPars; par++){
    par_vec.push_back(xx(par));
  }
  cov_par = mean(cov_xx);

}
