/*! \file
 *  \brief Routines to support loading correlators
 */

#include "correlator_util.h"
#include <ConfDataStoreDB.h>
#include <DBString.h>
#include <ensem/ensem.h>
#include <xml_array2d.h>
#include <string>
#include <map>

namespace CorrReaderEnv
{
  //----------------------------------------------------------------------------------
  namespace
  {
    bool test_done(const std::vector<bool>& in)
    {
      bool test = true;

      for(auto n = in.begin(); n != in.end(); ++n)
	test &= *n;

      return test;
    }
  }


  // **************************************************************************
  // filedb (redstar) database specific functions
  // Get the metadata
  std::string getMetaData(const std::string& dbase)
  {
    // Open DB solely to get the user data
    FILEDB::ConfDataStoreDB<FILEDB::StringKey,FILEDB::UserData> database;
    database.open(dbase, O_RDONLY, 0400);

    std::string meta_str;
    database.getUserdata(meta_str);

    return meta_str;
  }

  //! Read the dbtype from the metadata
  std::string getDBType(const std::string& dbfile)
  {
    std::string meta_str = getMetaData(dbfile);

    std::istringstream is(meta_str);
    XMLReader xml_in(is);
      
    std::string dbtype;
    read(xml_in, "/DBMetaData/id", dbtype); 

    return dbtype;
  }


  //----------------------------------------------------------------------------------
  //! Read ops list
  std::vector<std::string> readOpsList(const std::string& opslistfile)
  {
    // Read in list of operators to use at source and sink
    // Format is "irrep opname"
    std::ifstream opsListData(opslistfile.c_str());
    std::vector<std::string> opsList;

    if (! opsListData)
    {
      std::cerr << __func__ << ": opslistfile " << opslistfile << " read error" << std::endl;
      exit(1);
    }

    while(! opsListData.eof())
    {
      std::string irrep;
      std::string op;

      opsListData >> irrep >> op;

      opsList.push_back(op);
    } // end while

    int dim = opsList.size();
    std::cout << __func__ << ": using " << dim << " operators from " << opslistfile << ": ";

    for(int i = 0; i < dim; i++)
    {
      std::cout << opsList[i];

      if(i != dim - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    return opsList;
  }


  //----------------------------------------------------------------------------------
  // Rephase complex correlators
  std::vector< SEMBLE::SembleMatrix<double> > rephaseCorrs(const std::vector< SEMBLE::SembleMatrix<std::complex<double> > >& ComplexCorrs,
							   const std::string& rephaseMode,
							   const std::vector<std::string>& opsList,
							   int tmax) //NEEDS A TMAX
  {
    // Method to use to rephase correlators
    if(rephaseMode == "auto")
    {
      std::cout << __func__ << ": automatically determining correlator phases and rephasing" << std::endl;
    }
    else if(rephaseMode == "auto_positive")
    {
      std::cout << __func__ << ": automatically determining correlator phases (only +1 and +i allowed) and rephasing" << std::endl;
    }
    else if(rephaseMode == "real")
    {
      std::cout << __func__ << ": using real parts of complex correlators (not rephasing)" << std::endl;
    }
    else
    {
      std::cerr << __func__ << ": ERROR: Unknown rephaseMode " << rephaseMode
	   << " - rephaseMode must be one of \"auto\", \"auto_positive\", or \"real\" " << std::endl;
      exit(1);
    }

    int nbins = ComplexCorrs[0].getB();
    int dim   = ComplexCorrs[0].getN();
    int Lt    = ComplexCorrs.size();

    SEMBLE::SembleMatrix<double> phases(nbins, dim, dim);
    std::vector< std::vector<bool> > phase_set;
    std::vector<bool> tmp;
    tmp.resize(dim, false);

    for(int i = 0; i < dim; i++)
    {
      phase_set.push_back(tmp);
    }

    double numSigma = 5.0;

    Array<int> op_phase_re(dim), op_phase_im(dim);

    if((rephaseMode == "auto") || (rephaseMode == "auto_positive"))
    {
      // Determine the correlator phases
      // loop over elements
      for(int row = 0; row < dim; row++)
      {
	for(int col = 0; col < dim; col++)
	{

	  //try timeslices until either imag(C) signifcant or real(C) significant
	  bool signif = false;
	  int t = 1;
	  EnsemComplex c = ComplexCorrs[0].getEnsemElement(0, 0);

	  do
	  {
	    //find the raw phase of each correlation matrix element
	    c = ComplexCorrs[t].getEnsemElement(row, col);
	    //if both elements are constitent with zero, keep increasing the timeslice used until one isn't

	    if((abs(toDouble(mean(real(c)))) < numSigma * toDouble(sqrt(variance(real(c)))))
	       && (abs(toDouble(mean(imag(c)))) < numSigma * toDouble(sqrt(variance(imag(c))))))
	    {
	      signif = false;
	    }
	    else
	    {
	      //this might be OK

	      if((abs(toDouble(mean(real(c)))) >= numSigma * toDouble(sqrt(variance(real(c)))))
		 && (abs(toDouble(mean(imag(c)))) >= numSigma * toDouble(sqrt(variance(imag(c))))))
		std::cerr << "Warning: both real and imaginary parts of element C[t=" << t << "](" << row << "," << col << ") are significant" << std::endl;

	      EnsemReal phase =  atan2(imag(c) , real(c));

	      double mean_phase = toDouble(mean(phase));
	      double err_phase = toDouble(sqrt(variance(phase)));

	      //reject noisy phases
	      if(err_phase < 3.14159 / (2.0 * numSigma))
	      {
		phases.loadEnsemElement(row, col, phase);
		phase_set[row][col] = true;
		std::cout << "found the ACCEPTABLY noisy phase of element C[t=" << t << "](" << row << "," << col << ") = " << mean_phase / 3.14159 << " +/- " << err_phase / 3.14159 << std::endl;
		signif = true;
	      }
	      else
	      {
		std::cout << "found the UNACCEPTABLY noisy phase of element C[t=" << t << "](" << row << "," << col << ") = " << mean_phase / 3.14159 << " +/- " << err_phase / 3.14159 << std::endl;
		signif = false;
	      }
	    }

	    if(t == tmax)
	    {
	      std::cerr << "couldn't find any signal on any timeslice for C(" << row << "," << col << ")" << std::endl;
	      break;
	    };

	    t++;
	  }
	  while(!signif);

	  if(!signif)
	  {
	    phase_set[row][col] = false;
	  }

	}//next col
      }//next row

      for(int row = 0; row < dim; row++)
	for(int col = 0; col < dim; col++)
	{
	  if(!phase_set[row][col])
	  {
	    std::cout << "phase of element C(" << row << "," << col << ") is not determined" << std::endl;
	  }
	}

      //now have all the phases (not rounded off)

      std::vector<bool> op_phase_set(dim, false);
      op_phase_re[0] = 1;
      op_phase_im[0] = 0;
      op_phase_set[0] = true;
      // set phase of operator 0 to be real
      int col = 0;
      int count = 0;

      do
      {
	if(op_phase_set[col])  //can use this column since we know the ref phase
	{
	  std::cout << "trying column " << col << std::endl;

	  int ref_phase_re = op_phase_re[col];
	  int ref_phase_im = op_phase_im[col];

	  for(int row = 1; row < dim; row++)
	  {
	    if(!phase_set[row][col])
	    {
	      std::cerr << "  can't use row " << row << " since phase wasn't determined" << std::endl;  // can't use this, no reliable phase known
	      continue;
	    }

	    if(op_phase_set[row])
	    {
	      std::cout << "  already set phase for op " << row << std::endl;  // already set this op phase
	      continue;
	    }

	    // need to test if phase compatible with a pure phase WITHIN NOISE - IMPLEMENT
	    double phase = toDouble(mean(phases.getEnsemElement(row, col)));

	    int phase_re = int(round(cos(phase)) * ref_phase_re - round(sin(phase)) * ref_phase_im);
	    int phase_im = int(round(cos(phase)) * ref_phase_im + round(sin(phase)) * ref_phase_re);

	    if(((phase_re != 0) && (phase_im != 0)) || ((phase_re == 0) && (phase_im == 0)))
	    {
	      continue;  //this phase is not unique
	    }

	    //insert the found phase
	    std::cout << "setting the phase of operator #" << row << " using phase of C(" << row << "," << col << ")" << std::endl;
	    op_phase_re[row] = phase_re;
	    op_phase_im[row] = phase_im;
	    op_phase_set[row] = true;

	  }
	}

	col = (col + 1) % dim;
	count++;

	if(count == dim * dim)
	{
	  std::cerr << "looped over the matrix dim^2 times and still didn't set the phases, exiting" << std::endl;
	  exit(1);
	}

      }
      while (! test_done(op_phase_set));  //keep going until all the phases are set

      if(rephaseMode == "auto_positive")
      {
	// Set phase to be +1 or +i (not -1 or -i)
	for(int row = 0; row < dim; row++)
	{
	  op_phase_re[row] = abs(op_phase_re[row]);
	  op_phase_im[row] = abs(op_phase_im[row]);
	}
      }

    }
    else if(rephaseMode == "real")
    {
      // Set all the correlator phases to be real
      for(int row = 0; row < dim; row++)
      {
	op_phase_re[row] = 1;
	op_phase_im[row] = 0;
      }
    }
    else
    {
      std::cerr << __func__ << ": ERROR: Unknown rephaseMode " << rephaseMode
	   << " - rephaseMode must be one of \"auto\", \"auto_positive\", or \"real\" " << std::endl;
      exit(1);
    }


    
    // Write out phases to file
    std::ofstream opphases("ops_phases");
    for (int i = 0; i < dim; i++)
    {
      // std::cout << "i=  " << i << ": " << opsList[i] << " (" << phase_re[i] << " , " << phase_im[i] << ")" << std::endl;
      opphases << i << " " << opsList[i] << " " << op_phase_re[i] << " " << op_phase_im[i] << std::endl;
    }
    opphases.close();
    
    std::vector< SEMBLE::SembleMatrix<double> > Cttemp;
    SEMBLE::SembleMatrix<double> dum1(nbins, dim, dim);
    Cttemp.push_back(dum1);
    Cttemp.resize(Lt, dum1);

    for(int i = 0; i < dim; i++)
    {
      for(int j = 0; j < dim; j++)
      {
	// Only one of pp_re and pp_im should be non-zero
	int pp_re = int(op_phase_re[i] * op_phase_re[j]) + int(op_phase_im[i] * op_phase_im[j]);
	int pp_im = int(op_phase_im[i] * op_phase_re[j]) - int(op_phase_re[i] * op_phase_im[j]);
	int expected_time_rev_sign = 0;

	if(pp_re != 0 && pp_im == 0)
	{
	  // REAL ELEMENT
	  for(int t = 0; t < Lt; t++)
	  {
	    // add a check that this really is real - but only if significantly bigger than zero
	    double rr = abs(toDouble(mean(real(ComplexCorrs[t].getEnsemElement(i, j)))));
	    double ii = abs(toDouble(mean(imag(ComplexCorrs[t].getEnsemElement(i, j)))));
	    double ee = toDouble(sqrt(variance(imag(ComplexCorrs[t].getEnsemElement(i, j)))));

	    if((rr < ii) && (ii > 4.0 * ee) && (t > 0) && (t < tmax))
	    {
	      std::cerr << "WARNING : C[t=" << t << "](" << i << "," << j << ") doesn't have the expected phase - should be real, |imag| part is " << ii << " +/- " << ee  << ", " << ii / ee << " sigma discrepancy" << std::endl;
	    }

	    EnsemReal R = Real(pp_re) * real(ComplexCorrs[t].getEnsemElement(i, j));
	    Cttemp[t].loadEnsemElement(i, j, R);
	    expected_time_rev_sign = 1;
	  }
	}
	else if(pp_im != 0 && pp_re == 0)
	{
	  // IMAG ELEMENT
	  for(int t = 0; t < Lt; t++)
	  {
	    // add a check that this is really imag
	    double rr = abs(toDouble(mean(real(ComplexCorrs[t].getEnsemElement(i, j)))));
	    double ii = abs(toDouble(mean(imag(ComplexCorrs[t].getEnsemElement(i, j)))));
	    double ee = toDouble(sqrt(variance(real(ComplexCorrs[t].getEnsemElement(i, j)))));

	    if((rr > ii) && (rr > 4.0 * ee) && (t > 0) && (t < tmax))
	    {
	      std::cerr << "WARNING : C[t=" << t << "](" << i << "," << j << ") doesn't have the expected phase - should be imag, |real| part is " << rr << " +/- " << ee  << ", " << rr / ee << " sigma discrepancy" << std::endl;
	    }

	    EnsemReal R = -Real(pp_im) * imag(ComplexCorrs[t].getEnsemElement(i, j));
	    Cttemp[t].loadEnsemElement(i, j, R);
	    expected_time_rev_sign = -1;
	  }
	}
	else
	{
	  std::cerr << __func__ << ": Error: one and only one of pp_re and pp_im should be non-zero" << std::endl;
	  exit(1);
	}

      }
    }

    //DEBUG QUIT
    // exit(1);

    return Cttemp;
  }

} // namespace CorrReaderEnv
