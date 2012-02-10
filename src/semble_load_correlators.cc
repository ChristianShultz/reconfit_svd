#include "semble_load_correlators.h"

namespace SEMBLE
{
  SembleRCorrs::SembleRCorrs()
  {
    dim = 0;
    Lt = 0;
    irrepdim = -1;
    nbins = 0;
    foldTimeReversal = "";
    opsList.clear();
    opsListCConj.clear();
    Ct.clear();
  } //end constructor


// Load correlators from filedb (meson_2pt) database
  void SembleRCorrs::loadRephaseComplexCorrs(const string &dbfile, const string &opslistfile, int irrepdim_, FF::KeyHadron2PtCorr_t DefaultKeys, const string &avgMode, double avgTol, const string &badlistfile, const string &foldTimeReversal_, const string &rephaseMode_, const bool avgMom, const string &momListFile)
  {
    cout << __func__ << ": loading the real correlators from complex correlators in database" << endl;

    irrepdim = irrepdim_;
    foldTimeReversal = foldTimeReversal_;
    rephaseMode = rephaseMode_;
    readOpsList(opslistfile);

    if(avgMom)
      {
        readMomList(momListFile);
      }
    else
      {
        momList.clear();
        momList.push_back(DefaultKeys.mom);
      }

    SembleCCorrs ComplexCorrs;
    ComplexCorrs.loadFromDB(dbfile, opsList, opsListCConj, irrepdim, DefaultKeys, avgMode, avgTol, badlistfile, avgMom, momList);

    doRephasing(ComplexCorrs);

    nbins = ComplexCorrs.getBins();
    Lt = ComplexCorrs.getLt();

    cout << __func__ << ": finished loading the real correlators" << endl;
  }


// Load correlators from filedb (redstar) database
  void SembleRCorrs::loadRephaseComplexCorrs(const string &dbfile, const string &opslistfile, const Array<string>& opsxmlfiles, Array<int> avgRows, double avgTol, const string &badlistfile, const string &foldTimeReversal_, const string &rephaseMode_, const bool avgMom, const string &momListFile, InputPropsRedstarKeys_t keyParams)
  {
    cout << __func__ << ": loading the real correlators from complex correlators in database" << endl;

    if(avgRows.size() < 1)
      {
        cerr << __func__ << ": avgRows must contain at least 1 row" << endl;
        exit(1);
      }

    irrepdim = avgRows.size();
    foldTimeReversal = foldTimeReversal_;
    rephaseMode = rephaseMode_;
    readOpsList(opslistfile);

    if(avgMom)
      {
        readMomList(momListFile);
      }
    else
      {
        momList.clear();
        momList.push_back(keyParams.mom);
      }

    SembleCCorrs ComplexCorrs;
    ComplexCorrs.loadFromDB(dbfile, opsList, opsListCConj, opsxmlfiles, avgRows, avgTol, badlistfile, avgMom, momList, keyParams);

    doRephasing(ComplexCorrs);

    nbins = ComplexCorrs.getBins();
    Lt = ComplexCorrs.getLt();

    cout << __func__ << ": finished loading the real correlators" << endl;
  }


  void SembleRCorrs::doRephasing(SembleCCorrs &ComplexCorrs)
  {
    // Method to use to rephase correlators
    if(rephaseMode == "auto")
      {
        cout << __func__ << ": automatically determining correlator phases and rephasing" << endl;
      }
    else if(rephaseMode == "auto_positive")
      {
        cout << __func__ << ": automatically determining correlator phases (only +1 and +i allowed) and rephasing" << endl;
      }
    else if(rephaseMode == "real")
      {
        cout << __func__ << ": using real parts of complex correlators (not rephasing)" << endl;
      }
    else
      {
        cerr << __func__ << ": ERROR: Unknown rephaseMode " << rephaseMode
             << " - rephaseMode must be one of \"auto\", \"auto_positive\", or \"real\" " << endl;
        exit(1);
      }

    // Rephase complex correlators and put output in Ct;
    if(foldTimeReversal == "auto_majority")
      cout << __func__ << ": will fold in time reversal using majority decision to determine symmetry" << endl;
    else if(foldTimeReversal == "expected")
      cout << __func__ << ": will fold in time reversal using expected symmetry" << endl;
    else
      cout << __func__ << ": will NOT fold in time reversal" << endl;


    // Actually do the rephasing
    Ct = ComplexCorrs.rephaseCorrs(rephaseMode, foldTimeReversal, 20); //stupid hardwire of tmax
  }


  void SembleRCorrs::loadFromEnsemFiles(int dim_, const string &filepath)
  {
    cout << __func__ << ": loading the real correlators from ensemble files" << endl;

    dim = dim_;

    //size up Ct
    ostringstream real_file_name;
    real_file_name << filepath;
    real_file_name << "R_" << 0 << "_" << 0 ;
    EnsemVectorReal real;
    read(real_file_name.str(), real);
    nbins = peekObs(real, 0).size();
    Lt = real.numElem();
    cout << __func__ << ": dim =  " << dim << ", ensemble files have Lt = " << Lt << ", nbins = " << nbins << endl;

    SembleMatrix<double> dum(nbins, dim, dim);
    Ct.push_back(dum);
    Ct.resize(Lt, dum);

    //load the correlator files
    for(int i = 0; i < dim; ++i)
      {
        for(int j = 0; j < dim; ++j)
          {
            ostringstream real_file_name;
            real_file_name << filepath;
            real_file_name << "R_" << i << "_" << j ;
            EnsemVectorReal real;
            read(real_file_name.str(), real);
            cout << "read " << real_file_name.str() << endl;

            for(int t = 0; t < Lt; t++)
              {
                EnsemReal R = peekObs(real, t);
                Ct[t].loadEnsemElement(i, j, R);
              }

          }
      }

    cout << __func__ << ": finished loading the real correlators" << endl;
  }

  void SembleRCorrs::loadFromEnsemFile(const string &ensemfilename)
  {
    cout << __func__ << ": loading the real correlators from ensemble file: " << ensemfilename << endl;

    dim = 1;

    //size up Ct
    EnsemVectorReal real;
    read(ensemfilename, real);
    nbins = peekObs(real, 0).size();
    Lt = real.numElem();
    cout << __func__ << ": dim =  " << dim << ", ensemble file has Lt = " << Lt << ", nbins = " << nbins << endl;

    SembleMatrix<double> dum(nbins, dim, dim);
    Ct.push_back(dum);
    Ct.resize(Lt, dum);

    //load the correlator files
    for(int t = 0; t < Lt; t++)
      {
        EnsemReal R = peekObs(real, t);
        Ct[t].loadEnsemElement(0, 0, R);
      }

    cout << __func__ << ": finished loading the real correlators" << endl;
  }


  void SembleRCorrs::changeNumCfgs(int newbins)
  {
    // Decrease the number of configurations (nbins) in the correlators
    if(newbins >= nbins)
      {
        cerr << __func__ << ": can't decrease size if newbins (" << newbins << ") >= bins (" << nbins << ") " << endl;
        exit(1);
      }

    vector<SembleMatrix<double> > Ct_new;
    SembleMatrix<double> dum_new(newbins, dim, dim);
    Ct_new.push_back(dum_new);
    Ct_new.resize(Lt, dum_new);

    for(int t = 0; t < Lt; t++)
      {
        SembleMatrix<double> TempEnsemMat = Ct[t];

        for(int i = 0; i < dim; i++)
          {
            for(int j = 0; j < dim; j++)
              {
                EnsemReal TempEnsem = TempEnsemMat.getEnsemElement(i, j);
                EnsemReal NewTempEnsem;
                NewTempEnsem.resize(newbins);

                for(int k = 0; k < newbins; k++)
                  {
                    pokeEnsem(NewTempEnsem, peekEnsem(TempEnsem, k), k);
                  }

                Ct_new[t].loadEnsemElement(i, j, NewTempEnsem);
              }
          }
      }

    Ct = Ct_new;

    if(Ct[0].bins() != newbins)
      {
        cerr << __func__ << ": Error resizing ensembles (nbins = " << nbins << ", newbins = " << newbins << ")" << endl;
        exit(1);
      }

    cout << __func__ << ": Decreased nbins from " << nbins << " to " << newbins << endl;
    nbins = newbins;
  }


  void SembleRCorrs::useEnsemMean()
  {
    // Use the mean instead of an ensemble
    // Hack this by having an ensemble with one element, the mean of the original ensemble

    const int newbins = 1;
    vector<SembleMatrix<double> > Ct_new;
    SembleMatrix<double> dum_new(newbins, dim, dim);
    Ct_new.push_back(dum_new);
    Ct_new.resize(Lt, dum_new);

    for(int t = 0; t < Lt; t++)
      {
        SembleMatrix<double> TempEnsemMat = Ct[t];

        for(int i = 0; i < dim; i++)
          {
            for(int j = 0; j < dim; j++)
              {
                EnsemReal TempEnsem = TempEnsemMat.getEnsemElement(i, j);
                EnsemReal NewTempEnsem;
                NewTempEnsem.resize(newbins);
                pokeEnsem(NewTempEnsem, mean(TempEnsem), 0);
                Ct_new[t].loadEnsemElement(i, j, NewTempEnsem);
              }
          }
      }

    Ct = Ct_new;

    if(Ct[0].bins() != newbins)
      {
        cerr << __func__ << ": Error replacing the ensemble of correlators by the mean correlator" << endl;
        exit(1);
      }

    cout << __func__ << ": Replaced the ensemble (" << nbins << " bins) with the mean (an ensemble with " << newbins << " bins)" << endl;
    nbins = newbins;
  }


  void SembleRCorrs::skipTimeslices(int nt)
  {
    // Consider every n'th timeslice
    
    int Lt_new = floor( float(Lt) / float(nt) );
    vector<SembleMatrix<double> > Ct_new;
    SembleMatrix<double> dum_new(nbins, dim, dim);
    Ct_new.push_back(dum_new);
    Ct_new.resize(Lt_new, dum_new);
    
    for(int t=0; t<Lt_new; t++)
      {
	Ct_new[t] = Ct[t*nt];
      }
    
    Lt = Lt_new;
    Ct = Ct_new;
    
    //  cout << __func__ << ": Only using every n'th timeslice where n= " << nt << std::endl;
  }


  void SembleRCorrs::shiftCorrs(int dt)
  {
    // Consider C(t) - C(t+dt) instead of C(t)

    int Lt_new = Lt - dt;
    vector<SembleMatrix<double> > Ct_new;
    SembleMatrix<double> dum_new(nbins, dim, dim);
    Ct_new.push_back(dum_new);
    Ct_new.resize(Lt_new, dum_new);

    for(int t = 0; t < Lt_new; t++)
      {
        Ct_new[t] = Ct[t] - Ct[t + dt];
      }

    Lt = Lt_new;
    Ct = Ct_new;

    //  cout << __func__ << ": Using C(t) - C(t+dt) where dt= " << dt << " , new Lt = " << Lt << endl;
  }

  void SembleRCorrs::expWeightCorrs(double E)
  {
    // replace C(t) with exp(E*t)*C(t)
    vector<SembleMatrix<double> > Ct_new;
    SembleMatrix<double> dum_new(nbins, dim, dim);
    Ct_new.push_back(dum_new);
    Ct_new.resize(Lt, dum_new);

    for(int t = 0; t < Lt; t++)
      {
        Ct_new[t] = exp(Real(E * t)) * Ct[t];
      }

    Ct = Ct_new;
  }


  void SembleRCorrs::readOpsList(const string &opslistfile)
  {
    // Read in list of operators to use at source and sink
    opsList.clear();
    string oneline;
    ifstream opsListData(opslistfile.c_str());

    if(!opsListData)
      {
        cerr << __func__ << ": opslistfile " << opslistfile << " read error" << endl;
        exit(1);
      }

    while(getline(opsListData, oneline))
      {
        int sepIndex = oneline.find(" ");

        if(sepIndex < 0)
          {
            cerr << __func__ << ": Error: can't extract operator (" << opsList.size() << ") name and irrep from " << opslistfile << endl;
            exit(1);
          }

        string opIrrep = oneline.substr(0, sepIndex);
        string opName = oneline.substr(sepIndex + 1);
        string opParity = oneline.substr(opIrrep.size() - 2, 1);
        string opChargeConj = oneline.substr(opIrrep.size() - 1, 1);
        int opCConj = 0;

        if(opParity != "m" && opParity != "-" && opParity != "p" && opParity != "+")
          {
            if(foldTimeReversal != "none")
              {
                cerr << __func__ << ": Error: can't extract parity and/or charge conj for operator " << opsList.size() << " (" << opIrrep << " " << opName << ") from " << opslistfile << " and foldTimeReversal != none" << endl;
                exit(1);
              }
            else
              opCConj = 0;
          }
        else if((opChargeConj == "m") || (opChargeConj == "-"))
          opCConj = -1;
        else if((opChargeConj == "p") || (opChargeConj == "+"))
          opCConj = +1;
        else if(foldTimeReversal != "none")
          {
            cerr << __func__ << ": Error: can't extract parity and/or charge conj for operator " << opsList.size() << " (" << opIrrep << " " << opName << ") from " << opslistfile << " and foldTimeReversal != none" << endl;
            exit(1);
          }
        else
          opCConj = 0;

        opsList.push_back(opName);
        opsListCConj.push_back(opCConj);

        if(opCConj != 0)
          cout << __func__ << ": read opslist line: op = " << opName << ", irrep = " << opIrrep << ", Parity = " << opParity << ", CConj = " << opCConj << endl;
        else
          cout << __func__ << ": read opslist line: op = " << opName << ", irrep = " << opIrrep << ", couldn't extract parity and charge conj." << endl;

      }   // while(getlines ...)

    dim = opsList.size();
    cout << __func__ << ": using " << dim << " operators from " << opslistfile << ": ";

    for(int i = 0; i < dim; i++)
      {
        cout << opsList[i];

        if(i != dim - 1) cout << ", ";
      }

    cout << endl;
  }


  void SembleRCorrs::readMomList(const string &momListFile)
  {
    // Read in list of momenta to use
    momList.clear();
    string oneline;
    ifstream momListData(momListFile.c_str());

    if(!momListData)
      {
        cerr << __func__ << ": momListFile " << momListFile << " read error" << endl;
        exit(1);
      }

    while(getline(momListData, oneline))
      {
        momList.push_back(stringToIntArray3(oneline));
      }

    cout << __func__ << ": using the following (" << momList.size() << ") momenta from " << momListFile << ": ";

    for(int i = 0; i < momList.size(); i++)
      {
        for(int j = 0; j < momList[i].size(); j++)
          {
            cout << " " << (momList[i])[j];
          }

        cout << ", ";
      }

    cout << endl;
  }


// Split a string in to a 3D array of integers
  Array<int> SembleRCorrs::stringToIntArray3(const string &input)
  {

    int sepIndex1 = input.find(" ");

    if(sepIndex1 < 0)
      {
        cerr << __func__ << ": Error: can't extract array from string " << input << endl;
        exit(1);
      }

    int sepIndex2 = input.find(" ", sepIndex1 + 1);

    if(sepIndex2 < 0)
      {
        cerr << __func__ << ": Error: can't extract array from string " << input << endl;
        exit(1);
      }

    string e1 = input.substr(0, sepIndex1);
    string e2 = input.substr(sepIndex1 + 1, sepIndex2 - sepIndex1 - 1);
    string e3 = input.substr(sepIndex2 + 1, std::string::npos);

    Array<int> output(3);
    output[0] = atoi(e1.c_str());
    output[1] = atoi(e2.c_str());
    output[2] = atoi(e3.c_str());
    return output;
  }


  SembleMatrix<double> SembleRCorrs::getCt(int t) const
  {
    return Ct[t];
  }

  int SembleRCorrs::getLt() const
  {
    return Lt;
  }
  int SembleRCorrs::getDim() const
  {
    return dim;
  }
  int SembleRCorrs::getBins() const
  {
    return nbins;
  }


  EnsemVectorReal SembleRCorrs::getCij(int i, int j) const
  {
    EnsemVectorReal dum;
    dum.resize(Ct[0].bins());
    dum.resizeObs(Lt);

    for(int t = 0; t < Lt; t++)
      {
        pokeObs(dum, Ct[t].getEnsemElement(i, j), t);
      }

    return dum;
  }



// **************************************************************************
  SembleCCorrs::SembleCCorrs()
  {
    Lt = 0;
    nbins = 0;
    Ct.clear();
    irrepdim = -1;
    opsList.clear();
    opsListCConj.clear();
    dim = 0;
  }


// **************************************************************************
// filedb (meson_2pt) database specific functions

// Load from database
  void SembleCCorrs::loadFromDB(const string &dbfile, vector<string> opsList_, vector<int> opsListCConj_, int irrepdim_, FF::KeyHadron2PtCorr_t DefaultKeys, const string &avgMode, double avgTol, const string &badlistfile, bool avgMom, vector< Array<int> > momList_)
  {
    cout << __func__ << ": started loading the complex correlators" << endl;

    irrepdim = irrepdim_;
    opsList = opsList_;
    opsListCConj = opsListCConj_;
    dim = opsList.size();
    momList = momList_;

    // Open DB
    FILEDB::AllConfStoreDB< SerialDBKey<FF::KeyHadron2PtCorr_t>,  SerialDBData<SV> > database;

    if(database.open(dbfile, O_RDONLY, 0400) != 0)
      {
        cerr << __func__ << ": error opening database " << dbfile << endl;
        exit(1);
      }

    try
      {
        // Construct the keys to extract by iterating over operators
        Array<FF::KeyHadron2PtCorr_t> keys = createKeys(DefaultKeys, avgMode, avgMom);
        cout << __func__ << ": Will extract nkeys = " << keys.size() << endl;

        // Debugging
        // cout << __func__ << ": Keys are:" << endl;
        // for (int i = 0; i < keys.size(); i++)
        //   {
        //     cout << keys[i];
        //   }

        // Find number of time slices and check size of ensemble
        EnsemVectorComplex Test = printKeyValue<FF::KeyHadron2PtCorr_t, EnsemVectorComplex>(keys[0], database);
        Lt = Test.numElem();
        nbins = peekObs(Test, 0).size();
        cout << __func__ << ": filedb database (" << dbfile << ") has Lt = " << Lt << ", nbins = " << nbins << endl;

        // Get correlators from database (keys to match are specified in keys)
        // Average (lorentz/spin and/or momenta) if required
        Ct = loadCorrs(database, keys, avgMode, avgTol, badlistfile, avgMom);

        // Debugging
        // for (int j_src = 0; j_src < dim; j_src++)
        // {
        //   for (int j_snk = 0; j_snk < dim; j_snk++)
        //     {
        //       stringstream FileName;
        //       FileName << "corr_" << j_src << "," << j_snk << ".jack";
        //       write(FileName.str(),getCij(j_src,j_snk));
        //     }
        // }

      }
    catch(const std::string &e)
      {
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
      }
    catch(std::exception &e)
      {
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
      }

    cout << __func__ << ": finished loading the complex correlators" << endl;

  }


// Create db keys from operator list
  Array<FF::KeyHadron2PtCorr_t> SembleCCorrs::createKeys(FF::KeyHadron2PtCorr_t DefaultKeys, const string &avgMode, bool avgMom)
  {
    int nkeys = 0;

    if(avgMode == "none")
      nkeys = dim * dim;
    else if(avgMode == "lorentz")
      nkeys = irrepdim * dim * dim;
    else if(avgMode == "spin")
      nkeys =  irrepdim * dim * dim;
    else
      {
        cerr << __func__ << ": Error: unknown avgMode = " << avgMode << endl;
        exit(1);
      }

    if(avgMom) nkeys *= momList.size();

    Array<FF::KeyHadron2PtCorr_t> keys(nkeys);

    int count = 0;

    for(int j_src = 0; j_src < dim; j_src++)
      {
        for(int j_snk = 0; j_snk < dim; j_snk++)
          {

            if(avgMode == "none")
              {
                for(int i_mom = 0; i_mom < momList.size(); i_mom++)
                  {
                    keys[count].num_vecs = DefaultKeys.num_vecs;
                    keys[count].src_smear = DefaultKeys.src_smear;
                    keys[count].src_spin = DefaultKeys.src_spin;
                    keys[count].src_lorentz = DefaultKeys.src_lorentz;
                    keys[count].snk_smear = DefaultKeys.snk_smear;
                    keys[count].snk_spin = DefaultKeys.snk_spin;
                    keys[count].snk_lorentz = DefaultKeys.snk_lorentz;
                    keys[count].mom = momList[i_mom];
                    keys[count].mass = DefaultKeys.mass;
                    keys[count].ensemble = DefaultKeys.ensemble;

                    keys[count].src_name = opsList[j_src];
                    keys[count].snk_name = opsList[j_snk];
                    count++;
                  } // loop over i_mom
              } // avgMode == none

            else if(avgMode == "lorentz")
              {
                for(int i = 0; i < irrepdim; i++)
                  {
                    for(int i_mom = 0; i_mom < momList.size(); i_mom++)
                      {
                        keys[count].num_vecs = DefaultKeys.num_vecs;
                        keys[count].src_smear = DefaultKeys.src_smear;
                        keys[count].src_spin = DefaultKeys.src_spin;
                        keys[count].snk_smear = DefaultKeys.snk_smear;
                        keys[count].snk_spin = DefaultKeys.snk_spin;
                        keys[count].mom = momList[i_mom];
                        keys[count].mass = DefaultKeys.mass;
                        keys[count].ensemble = DefaultKeys.ensemble;

                        keys[count].src_name = opsList[j_src];
                        keys[count].snk_name = opsList[j_snk];

                        if(irrepdim != 1)
                          {
                            keys[count].src_lorentz.resize(1);
                            keys[count].snk_lorentz.resize(1);
                            keys[count].src_lorentz[0] = i;
                            keys[count].snk_lorentz[0] = i;
                          }

                        count++;
                      } // loop over i_mom
                  } // loop over i (irrepdim)
              } // avgMode == lorentz

            else if(avgMode == "spin")
              {
                for(int i = 1; i <= irrepdim; i++)
                  {
                    for(int i_mom = 0; i_mom < momList.size(); i_mom++)
                      {
                        keys[count].num_vecs = DefaultKeys.num_vecs;
                        keys[count].src_smear = DefaultKeys.src_smear;
                        keys[count].src_lorentz = DefaultKeys.src_lorentz;
                        keys[count].snk_smear = DefaultKeys.snk_smear;
                        keys[count].snk_lorentz = DefaultKeys.snk_lorentz;
                        keys[count].mom = momList[i_mom];
                        keys[count].mass = DefaultKeys.mass;
                        keys[count].ensemble = DefaultKeys.ensemble;

                        keys[count].src_name = opsList[j_src];
                        keys[count].snk_name = opsList[j_snk];
                        keys[count].src_spin = i;
                        keys[count].snk_spin = i;
                        count++;
                      } // loop over i_mom
                  } // loop over i (irrepdim)
              } // avgMode == spin

            else
              {
                cerr << __func__ << ": unknown avgMode = " << avgMode << endl;
                exit(1);
              }

          } // loop over j_snk
      } // loop over j_src

    if(count != nkeys)
      {
        cerr << __func__ << ": error: count not equal to " << nkeys << endl;
        exit(1);
      }

    return keys;
  }


// Get correlators from database, averaging over lorentz/spin if required
  vector<SembleMatrix<std::complex<double> > > SembleCCorrs::loadCorrs(FILEDB::AllConfStoreDB< SerialDBKey<FF::KeyHadron2PtCorr_t>,  SerialDBData<SV> >& database,  Array<FF::KeyHadron2PtCorr_t> keys, const string &avgMode, double avgTol, const string &badlistfile, bool avgMom)
  {
    SembleMatrix<std::complex<double> > dum(nbins, dim, dim);
    vector<SembleMatrix<std::complex<double> > > ComplexCorrs;
    ComplexCorrs.push_back(dum);
    ComplexCorrs.resize(Lt, dum);


    if((avgMode == "none") && (!avgMom))
      {
        cout << __func__ << ": extracting complex correlators from database and not averaging" << endl;
        int count = 0;

        for(int j_src = 0; j_src < dim; j_src++)
          {
            for(int j_snk = 0; j_snk < dim; j_snk++)
              {
                EnsemVectorComplex TempCorr;
                TempCorr = printKeyValue<FF::KeyHadron2PtCorr_t, EnsemVectorComplex>(keys[count], database);
                //  cout << __func__ << ": getting from db key (" << j_src << "," << j_snk << "): " << keys[count];
                count++;

                for(int t = 0; t < Lt; t++)
                  {
                    EnsemComplex TempC = peekObs(TempCorr, t);
                    ComplexCorrs[t].loadEnsemElement(j_src, j_snk, TempC);
                  }
              }
          }
      } // avgMode == none


    else if((avgMode == "lorentz") || (avgMode == "spin") || ((avgMode == "none") && avgMom))
      {
        int avgdim = 0;

        if(avgMom) avgdim = irrepdim * momList.size();
        else avgdim = irrepdim;

        if(avgMode == "lorentz")
          {
            if(avgMom)
              cout << __func__ << ": extracting complex correlators from database and averaging over momenta and diagonal lorentz indicies" << endl;
            else
              cout << __func__ << ": extracting complex correlators from database and averaging over diagonal lorentz indicies" << endl;
          }
        else if(avgMode == "spin")
          {
            if(avgMom)
              cout << __func__ << ": extracting complex correlators from database and averaging over momenta and spin indicies" << endl;
            else
              cout << __func__ << ": extracting complex correlators from database and averaging over spin indicies" << endl;
          }
        else if(avgMode == "none" && avgMom)
          {
            cout << __func__ << ": extracting complex correlators from database and averaging over momenta only" << endl;
            avgdim = momList.size();
          }

        ofstream badlist(badlistfile.c_str());

        if(!badlist)
          {
            cerr << __func__ << ": error: bad ops list file (" << badlistfile << ") write error" << endl;
            exit(1);
          }

        badlist << "These operator combinations fail the row average test:" << endl;

        int count = 0;

        for(int j_src = 0; j_src < dim; j_src++)
          {
            for(int j_snk = 0; j_snk < dim; j_snk++)
              {
                Array<EnsemVectorComplex> TempCorrs(avgdim);
                TempCorrs[0] = printKeyValue<FF::KeyHadron2PtCorr_t, EnsemVectorComplex>(keys[count], database);
                //  cout << __func__ << ": ----- start of avg for " << j_src << "," << j_snk << " -----" << endl;
                //  cout << __func__ << ": getting from db key: " << keys[count];
                EnsemVectorComplex TempCorr = TempCorrs[0];
                count++;

                // Average over spin/lorentz and/or momenta
                for(int i = 1; i < avgdim; i++)
                  {
                    TempCorrs[i] = printKeyValue<FF::KeyHadron2PtCorr_t, EnsemVectorComplex>(keys[count], database);
                    // cout << __func__ << ": getting from db key: " << keys[count];
                    TempCorr += TempCorrs[i];
                    count++;
                  }

                TempCorr = TempCorr / Real(avgdim);

                int flag = 0;

                for(int t = 0; t < Lt; t++)
                  {
                    EnsemComplex TempC = peekObs(TempCorr, t);
                    ComplexCorrs[t].loadEnsemElement(j_src, j_snk, TempC);

                    // Check that each term in the average is within tolerance of the average
                    for(int k = 0; k < avgdim; k++)
                      {
                        EnsemComplex TempElement = peekObs(TempCorrs[k], t);

                        if((toDouble(mean(real(TempElement))) < 3.0 * toDouble(sqrt(variance(TempElement))))
                            && (toDouble(mean(real(TempC))) < 3.0 * toDouble(sqrt(variance(TempC)))))
                          continue;   // Don't perform the test if both this particular term and the average are consistent (3 sigma) with zero

                        EnsemComplex Ratio = Real(avgdim) * TempC / TempElement;
                        Real TestRatio = real(mean(Ratio)) - Real(avgdim);
                        Real errorRatio2 = variance(Ratio);

                        if((toDouble(TestRatio * TestRatio) > avgTol * avgTol * toDouble(errorRatio2)) && (avgdim > 1) && (t > 0))
                          {
                            cout << __func__ << ": Correlator " << opsList[j_src] << " , " << opsList[j_snk] << " has avg outside tolerance at t = " << t << ", component = " << k << endl;
                            cout << "Ratio is " << mean(Ratio) << ", error is " << sqrt(errorRatio2) << endl;

                            if(flag == 0)
                              {
                                badlist << opsList[j_src] << " , " << opsList[j_snk] << endl;
                                flag = 1;
                              }
                          }
                      }   // loop over k
                  }   // loop over t

                flag = 0;

              }   // loop over j_snk
          }   // loop over j_src
      }   // avgMode == lorentz or spin


    else
      {
        cerr << __func__ << ": Error: unknown avgMode = " << avgMode << endl;
        exit(1);
      }

    return ComplexCorrs;
  }


// **************************************************************************
// filedb (redstar) database specific functions

// Load from database
  void SembleCCorrs::loadFromDB(const string &dbfile, vector<string> opsList_, vector<int> opsListCConj_, const Array<string>& opsxmlfiles, Array<int> avgRows, double avgTol, const string &badlistfile, bool avgMom, vector< Array<int> > momList_, InputPropsRedstarKeys_t keyParams)
  {
    cout << __func__ << ": started loading the complex correlators" << endl;

    if(avgRows.size() < 1)
      {
        cerr << __func__ << ": avgRows must contain at least 1 row" << endl;
        exit(1);
      }

    irrepdim = avgRows.size();
    opsList = opsList_;
    opsListCConj = opsListCConj_;
    dim = opsList.size();
    momList = momList_;

    // Open DB
    FILEDB::AllConfStoreDB< SerialDBKey<Hadron::KeyHadronNPartNPtCorr_t>,  SerialDBData<SV> > database;

    if(database.open(dbfile, O_RDONLY, 0400) != 0)
      {
        cerr << __func__ << ": error opening database " << dbfile << endl;
        exit(1);
      }

    try
      {
        // Construct the keys to extract by iterating over operators
        Array<Hadron::KeyHadronNPartNPtCorr_t> keys = createRedstarKeys(opsxmlfiles, keyParams, avgRows, avgMom);
        cout << __func__ << ": Will extract nkeys = " << keys.size() << endl;

        // Debugging
        cout << __func__ << ": Keys are:" << endl;

        for(int i = 0; i < keys.size(); i++)
          {
            cout << keys[i];
          }
       

        // Find number of time slices and check size of ensemble
        EnsemVectorComplex Test = printKeyValue<Hadron::KeyHadronNPartNPtCorr_t, EnsemVectorComplex>(keys[0], database);
        Lt = Test.numElem();
        nbins = peekObs(Test, 0).size();
        cout << __func__ << ": filedb database (" << dbfile << ") has Lt = " << Lt << ", nbins = " << nbins << endl;

        // Get correlators from database (keys to match are specified in keys)
        // Average over rows and/or momenta if required
        Ct = loadCorrs(database, keys, avgRows, avgTol, badlistfile, avgMom);

        // Debugging
        // for (int j_src = 0; j_src < dim; j_src++)
        // {
        //   for (int j_snk = 0; j_snk < dim; j_snk++)
        //     {
        //       stringstream FileName;
        //       FileName << "corr_" << j_src << "," << j_snk << ".jack";
        //       write(FileName.str(),getCij(j_src,j_snk));
        //     }
        // }

      }
    catch(const std::string &e)
      {
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
      }
    catch(std::exception &e)
      {
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
      }

    cout << __func__ << ": finished loading the complex correlators" << endl;

  }


// Create db keys from operator list
  Array<Hadron::KeyHadronNPartNPtCorr_t> SembleCCorrs::createRedstarKeys(const Array<string>& opsxmlfiles, InputPropsRedstarKeys_t keyParams, const Array<int> avgRows, bool avgMom)
  {
    int nkeys = dim * dim * irrepdim;

    if(avgMom) nkeys *= momList.size();

    Array<Hadron::KeyHadronNPartIrrep_t> opsxml = readOpsxml(opsxmlfiles);

    if(opsxml.size() != dim)
      {
        cerr << __func__ << ": ERROR: only " << opsxml.size() << " elems read from opsxml files but there should be " << dim << endl;
        exit(1);
      }

    Array<Hadron::KeyHadronNPartNPtCorr_t> keys(nkeys);
    int count = 0;

    for(int j_src = 0; j_src < dim; j_src++)
      {
        for(int j_snk = 0; j_snk < dim; j_snk++)
          {
            for(int i = 0; i < irrepdim; i++)
              {
                for(int i_mom = 0; i_mom < momList.size(); i_mom++)
                  {
                    // cout << count << ": " << j_src << ", " << j_snk << ", " << i << ", " << i_mom << endl;
                    keys[count].ensemble = keyParams.ensemble;
                    keys[count].npoint.resize(2);  // This is a 2-point correlator. N.B. THIS IS Array1d0 - A ONE BASED ARRAY

                    // The sink op
                    keys[count].npoint[1].t_slice           = -2;
                    keys[count].npoint[1].irrep.row         = avgRows[i];
                    keys[count].npoint[1].irrep.twoI_z      = keyParams.twoI_z;
                    keys[count].npoint[1].irrep.mom         = momList[i_mom];
                    keys[count].npoint[1].irrep.creation_op = false;
                    keys[count].npoint[1].irrep.smearedP    = opsxml[j_snk].smearedP;
                    keys[count].npoint[1].irrep.CGs         = opsxml[j_snk].CGs;
                    keys[count].npoint[1].irrep.ops         = opsxml[j_snk].ops;

                    // The source op
                    keys[count].npoint[2].t_slice           = keyParams.source_tslice;
                    keys[count].npoint[2].irrep.row         = avgRows[i];
                    keys[count].npoint[2].irrep.twoI_z      = keyParams.twoI_z;
                    keys[count].npoint[2].irrep.mom         = momList[i_mom];
                    keys[count].npoint[2].irrep.creation_op = true;
                    keys[count].npoint[2].irrep.smearedP    = opsxml[j_src].smearedP;
                    keys[count].npoint[2].irrep.CGs         = opsxml[j_src].CGs;
                    keys[count].npoint[2].irrep.ops         = opsxml[j_src].ops;

                    count++;
                  } // loop over i_mom
              } // loop over i (irrepdim)
          } // loop over j_snk
      } // loop over j_src

    if(count != nkeys)
      {
        cerr << __func__ << ": error: count not equal to " << nkeys << endl;
        exit(1);
      }

    return keys;
  }


  Array<Hadron::KeyHadronNPartIrrep_t> SembleCCorrs::readOpsxml(const Array<string>& opsxmlfiles)
  {
    if(opsxmlfiles.size() < 1)
      {
        cerr << __func__ << ": ERROR: there must be at least one opsxmlfile" << endl;
        exit(1);
      }

    Array<opsxml_t> opsxml;

    for(int i = 0; i < opsxmlfiles.size(); i++)
      {
        Array<opsxml_t> opsxml_temp;

        try
          {
            XMLReader xml_in(opsxmlfiles[i]);
            read(xml_in, "/OpsList", opsxml_temp);
          }
        catch(const std::string &e)
          {
            cerr << __func__ << ": ERROR: can't read opsxmlfile (" << opsxmlfiles[i] << "): " << e << endl;
            exit(1);
          }

        cout << __func__ << ": read " << opsxml_temp.size() << " elems from " << opsxmlfiles[i] << endl;

        if(i == 0)
          {
            opsxml = opsxml_temp;
          }
        else
          {
            Array<opsxml_t> opsxml_orig = opsxml;
            opsxml = concat(opsxml_orig, opsxml_temp);
          }
      }

    Array<Hadron::KeyHadronNPartIrrep_t> keys(dim);

    for(int i = 0; i < dim; i++)
      {
        bool found = false;

        for(int k = 0; k < opsxml.size(); k++)
          {
            if(opsxml[k].opnamekey == opsList[i])
              {
                keys[i] = opsxml[k].irrep;
                found = true;
                // cout << "Op " << i << " (" << opsList[i] << ") matches with xml elem " << k << endl;
                break;
              }
          }

        if(!found)
          {
            cerr << __func__ << ": ERROR: can't find op " << i << " (" << opsList[i] << ") ";
            cerr << "in opsxmlfiles\n";
            exit(1);
          }
      }

    return keys;
  }


// Get correlators from database, averaging over lorentz/spin if required
  vector<SembleMatrix<std::complex<double> > > SembleCCorrs::loadCorrs(FILEDB::AllConfStoreDB< SerialDBKey<Hadron::KeyHadronNPartNPtCorr_t>,  SerialDBData<SV> >& database,  Array<Hadron::KeyHadronNPartNPtCorr_t> keys, const Array<int> avgRows, double avgTol, const string &badlistfile, bool avgMom)
  {
    SembleMatrix<std::complex<double> > dum(nbins, dim, dim);
    vector<SembleMatrix<std::complex<double> > > ComplexCorrs;
    ComplexCorrs.push_back(dum);
    ComplexCorrs.resize(Lt, dum);

    int avgdim = 0;

    if(avgMom) avgdim = irrepdim * momList.size();
    else avgdim = irrepdim;

    if(avgMom)
      cout << __func__ << ": extracting complex correlators from database and averaging over momenta and diagonal rows  ";
    else
      cout << __func__ << ": extracting complex correlators from database and averaging over diagonal rows ";

    cout << avgRows[0];

    for(int i = 1; i < irrepdim; i++)
      cout << ", " << avgRows[i];

    cout << endl;

    ofstream badlist(badlistfile.c_str());

    if(!badlist)
      {
        cerr << __func__ << ": error: bad ops list file (" << badlistfile << ") write error" << endl;
        exit(1);
      }

    badlist << "These operator combinations fail the row average test:" << endl;

    int count = 0;

    for(int j_src = 0; j_src < dim; j_src++)
      {
        for(int j_snk = 0; j_snk < dim; j_snk++)
          {
            Array<EnsemVectorComplex> TempCorrs(avgdim);
            TempCorrs[0] = printKeyValue<Hadron::KeyHadronNPartNPtCorr_t, EnsemVectorComplex>(keys[count], database);
            //  cout << __func__ << ": ----- start of avg for " << j_src << "," << j_snk << " -----" << endl;
            //  cout << __func__ << ": getting from db key: " << keys[count];
            EnsemVectorComplex TempCorr = TempCorrs[0];
            count++;

            // Average over rows and/or momenta
            for(int i = 1; i < avgdim; i++)
              {
                TempCorrs[i] = printKeyValue<Hadron::KeyHadronNPartNPtCorr_t, EnsemVectorComplex>(keys[count], database);
                // cout << __func__ << ": getting from db key: " << keys[count];
                TempCorr += TempCorrs[i];
                count++;
              }

            TempCorr = TempCorr / Real(avgdim);

            int flag = 0;

            for(int t = 0; t < Lt; t++)
              {
                EnsemComplex TempC = peekObs(TempCorr, t);
                ComplexCorrs[t].loadEnsemElement(j_src, j_snk, TempC);

                // Check that each term in the average is within tolerance of the average
                for(int k = 0; k < avgdim; k++)
                  {
                    EnsemComplex TempElement = peekObs(TempCorrs[k], t);

                    if((toDouble(mean(real(TempElement))) < 3.0 * toDouble(sqrt(variance(TempElement))))
                        && (toDouble(mean(real(TempC))) < 3.0 * toDouble(sqrt(variance(TempC)))))
                      continue;   // Don't perform the test if both this particular term and the average are consistent (3 sigma) with zero

                    EnsemComplex Ratio = Real(avgdim) * TempC / TempElement;
                    Real TestRatio = real(mean(Ratio)) - Real(avgdim);
                    Real errorRatio2 = variance(Ratio);

                    if((toDouble(TestRatio * TestRatio) > avgTol * avgTol * toDouble(errorRatio2)) && (avgdim > 1) && (t > 0))
                      {
                        cout << __func__ << ": Correlator " << opsList[j_src] << " , " << opsList[j_snk] << " has avg outside tolerance at t = " << t << ", component = " << k << endl;
                        cout << "Ratio is " << mean(Ratio) << ", error is " << sqrt(errorRatio2) << endl;

                        if(flag == 0)
                          {
                            badlist << opsList[j_src] << " , " << opsList[j_snk] << endl;
                            flag = 1;
                          }
                      }
                  }   // loop over k
              }   // loop over t

            flag = 0;

          }   // loop over j_snk
      }   // loop over j_src

    return ComplexCorrs;
  }


// **************************************************************************
// Rephasing stuff

// Rephase complex correlators and possibly fold in time reversal;
  vector<SembleMatrix<double> > SembleCCorrs::rephaseCorrs(const string &rephaseMode, const string &foldTimeReversal, int tmax) //NEEDS A TMAX
  {
    vector<SembleMatrix<std::complex<double> > > ComplexCorrs = Ct;

    SembleMatrix<double> phases(nbins, dim, dim);
    vector< vector<bool> > phase_set;
    vector<bool> tmp;
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
                EnsemComplex c = Ct[0].getEnsemElement(0, 0);

                do
                  {
                    //find the raw phase of each correlation matrix element
                    c = Ct[t].getEnsemElement(row, col);
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
                          cerr << "Warning: both real and imaginary parts of element C[t=" << t << "](" << row << "," << col << ") are significant" << endl;

                        EnsemReal phase =  atan2(imag(c) , real(c));

                        double mean_phase = toDouble(mean(phase));
                        double err_phase = toDouble(sqrt(variance(phase)));

                        //reject noisy phases
                        if(err_phase < 3.14159 / (2.0 * numSigma))
                          {
                            phases.loadEnsemElement(row, col, phase);
                            phase_set[row][col] = true;
                            cout << "found the ACCEPTABLY noisy phase of element C[t=" << t << "](" << row << "," << col << ") = " << mean_phase / 3.14159 << " +/- " << err_phase / 3.14159 << endl;
                            signif = true;
                          }
                        else
                          {
                            cout << "found the UNACCEPTABLY noisy phase of element C[t=" << t << "](" << row << "," << col << ") = " << mean_phase / 3.14159 << " +/- " << err_phase / 3.14159 << endl;
                            signif = false;
                          }
                      }

                    if(t == tmax)
                      {
                        cerr << "couldn't find any signal on any timeslice for C(" << row << "," << col << ")" << endl;
                        break;
                      };

                    t++;
                  }
                while(!signif);

                if(!signif)
                  {
                    phase_set[row][col] = false;
                  }

                /*      if(signif){
                //if a significant signal was found, compute the phase
                EnsemReal phase =  atan2( imag( c ) , real( c ) );

                double mean_phase = toDouble(mean(phase));
                double err_phase = toDouble(sqrt(variance(phase)));

                //reject noisy phases
                if(err_phase < 3.14159 /( 2.0 * numSigma) ){
                phases.loadEnsemElement(row, col, phase);
                phase_set[row][col] = true;
                cout << "found the ACCEPTABLY noisy phase of element C[t=" << t - 1 <<"](" << row << "," << col <<") = " << mean_phase << " +/- " << err_phase << endl;
                }
                else{
                cout << "found the UNACCEPTABLLY noisy phase of element C[t=" << t - 1 <<"](" << row << "," << col <<") = " << mean_phase << " +/- " << err_phase << endl;
                }

                }
                else{
                phase_set[row][col] = false;
                }
                */

              }//next col
          }//next row

        for(int row = 0; row < dim; row++)
          for(int col = 0; col < dim; col++)
            {
              if(!phase_set[row][col])
                {
                  cout << "phase of element C(" << row << "," << col << ") is not determined" << endl;
                }
            }

        //now have all the phases (not rounded off)

        vector<bool> op_phase_set(dim, false);
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
                cout << "trying column " << col << endl;

                int ref_phase_re = op_phase_re[col];
                int ref_phase_im = op_phase_im[col];

                for(int row = 1; row < dim; row++)
                  {
                    if(!phase_set[row][col])
                      {
                        cerr << "  can't use row " << row << " since phase wasn't determined" << endl;  // can't use this, no reliable phase known
                        continue;
                      }

                    if(op_phase_set[row])
                      {
                        cout << "  already set phase for op " << row << endl;  // already set this op phase
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
                    cout << "setting the phase of operator #" << row << " using phase of C(" << row << "," << col << ")" << endl;
                    op_phase_re[row] = phase_re;
                    op_phase_im[row] = phase_im;
                    op_phase_set[row] = true;

                  }
              }

            col = (col + 1) % dim;
            count++;

            if(count == dim * dim)
              {
                cerr << "looped over the matrix dim^2 times and still didn't set the phases, exiting" << endl;
                exit(1);
              }

          }
        while(!test_done(op_phase_set));  //keep going until all the phases are set

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
        cerr << __func__ << ": ERROR: Unknown rephaseMode " << rephaseMode
             << " - rephaseMode must be one of \"auto\", \"auto_positive\", or \"real\" " << endl;
        exit(1);
      }



    // Write out phases to file
    ofstream opphases("ops_phases");

    for(int i = 0; i < dim; i++)
      {
        // cout << "i=  " << i << ": " << opsList[i] << " (" << phase_re[i] << " , " << phase_im[i] << ")" << endl;
        opphases << i << " " << opsList[i] << " " << op_phase_re[i] << " " << op_phase_im[i] << endl;
      }

    opphases.close();

    vector<SembleMatrix<double> > Cttemp;
    SembleMatrix<double> dum1(nbins, dim, dim);
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
                        cerr << "WARNING : C[t=" << t << "](" << i << "," << j << ") doesn't have the expected phase - should be real, |imag| part is " << ii << " +/- " << ee  << ", " << ii / ee << " sigma discrepancy" << endl;
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
                        cerr << "WARNING : C[t=" << t << "](" << i << "," << j << ") doesn't have the expected phase - should be imag, |real| part is " << rr << " +/- " << ee  << ", " << rr / ee << " sigma discrepancy" << endl;
                      }

                    EnsemReal R = -Real(pp_im) * imag(ComplexCorrs[t].getEnsemElement(i, j));
                    Cttemp[t].loadEnsemElement(i, j, R);
                    expected_time_rev_sign = -1;
                  }
              }
            else
              {
                cerr << __func__ << ": Error: one and only one of pp_re and pp_im should be non-zero" << endl;
                exit(1);
              }


            //#############################################################################################
            //Folding in time reversal
            if((foldTimeReversal == "auto_majority") || (foldTimeReversal == "expected"))
              {

                // If the correlator changes charge conj (or would if charge conj was a good sym) then there's an extra minus sign:
                if((opsListCConj[i] * opsListCConj[j]) == -1)
                  {
                    expected_time_rev_sign = - expected_time_rev_sign;
                    // cout << __func__ << ": Correlator (" << i << "," << j << ") changes charge conj and so has extra -ve sign in time reversal sym" << endl;
                  }
                else if(opsListCConj[i] * opsListCConj[j] != 1)
                  {
                    cerr << __func__ << ": Error: problem with charge conj of correlator (" << i << ", " << j << ") = (" << opsListCConj[i] << ", " << opsListCConj[j] << ")" << endl;
                    exit(1);
                  }

                // Take a majority decision on the time-reversal symmetry of the correlator based on 11 timeslices
                int time_rev_sign = 0;

                for(int tTest = 1; tTest <= 11; tTest++)
                  {
                    EnsemReal TestR_pos = Cttemp[tTest].getEnsemElement(i, j) - Cttemp[Lt - tTest].getEnsemElement(i, j);
                    EnsemReal TestR_neg = Cttemp[tTest].getEnsemElement(i, j) + Cttemp[Lt - tTest].getEnsemElement(i, j);

                    if(abs(toDouble(mean(TestR_pos))) < abs(toDouble(mean(TestR_neg))))
                      time_rev_sign++;
                    else
                      time_rev_sign--;
                  }

                if(time_rev_sign > 0) time_rev_sign = +1;
                else if(time_rev_sign < 0) time_rev_sign = -1;
                else if(foldTimeReversal == "auto_majority")
                  {
                    cerr << __func__ << ": Error: Majority decision to decide on time-reversal symmetry of correlator failed ";
                    cerr << "for correlator " << i << "(" << opsList[i] << "), " << j << "(" << opsList[j] << ") " << endl;
                    exit(1);
                  }

                // Check if majority decision is the same as the 'expected' sign
                if(time_rev_sign != expected_time_rev_sign)
                  {
                    cout << "Warning: correlator " << i << "(" << opsList[i] << "), " << j << "(" << opsList[j] << ") ";
                    cout << "symmetry from majority decision (" << time_rev_sign << ") is not the same as the expected time reversal symmetry (" << expected_time_rev_sign << ") ";

                    if(foldTimeReversal == "expected")
                      cout << " -- using expected sign (" << expected_time_rev_sign << ")" << endl;
                    else
                      cout << " -- using majority decision sign (" << time_rev_sign << ")" << endl;
                  }

                // Is the majority decision or the expected sign used
                if(foldTimeReversal == "expected")
                  {
                    time_rev_sign = expected_time_rev_sign;
                  }
                else if(foldTimeReversal == "auto_majority")
                  {
                    //  time_rev_sign = time_rev_sign;
                  }
                else
                  {
                    cerr << __func__ << ": Error: unknown foldTimeReversal = " << foldTimeReversal << endl;
                    exit(1);
                  }

                // Perform time reveral and fold in
                vector<EnsemReal> Ctcopy;

                for(int t = 0; t < Lt; t++)
                  {
                    Ctcopy.push_back((Real(0.5)*Cttemp[t].getEnsemElement(i, j)));
                  }

                EnsemReal R = Ctcopy[0];
                R += R;
                Cttemp[0].loadEnsemElement(i, j, R);

                for(int t = 1; t < Lt; t++)
                  {
                    EnsemReal Rtemp = Ctcopy[t];
                    Rtemp += Real(time_rev_sign) * Ctcopy[Lt - t];
                    Cttemp[t].loadEnsemElement(i, j, Rtemp);
                  }
              }  // if (foldTimeReversal == ...)

          }
      }

    //DEBUG QUIT
    //  exit(1);

    return Cttemp;
  }


  bool test_done(vector<bool> in)
  {
    bool test = true;

    for(int i = 0; i < in.size(); i++)
      {
        if(!in[i])
          {
            test = false;
            break;
          }
      }

    return test;
  }

  //*************************************************************************
  void loadCorr(SembleRCorrs &twoPoints, const FitIniProps_t &inikeys)
  {

    if((inikeys.dbInputType == "ensem") || (inikeys.dbInputType == "ensem_debug"))
      {
        cout << "Using ensem files from directory: " << inikeys.inputPropsEnsem.dbFname << endl;
        twoPoints.loadFromEnsemFiles(inikeys.inputPropsEnsem.dim, inikeys.inputPropsEnsem.dbFname);
      }
    else if(inikeys.dbInputType == "ensem_onecorr")
      {
        cout << "Using ensem file: " << inikeys.inputPropsEnsem.dbFname << endl;
        twoPoints.loadFromEnsemFile(inikeys.inputPropsEnsem.dbFname);
      }

    else if((inikeys.dbInputType == "dbnew") || (inikeys.dbInputType == "dbnew_debug"))
      {
        cout << "Using filedb (meson_2pt) database: " << inikeys.inputPropsDB.dbFname << endl;

        // Check foldTimeReversal value
        if((inikeys.inputPropsDB.foldTimeReversal != "none") && (inikeys.inputPropsDB.foldTimeReversal != "auto_majority") && (inikeys.inputPropsDB.foldTimeReversal != "expected"))
          {
            cerr << __func__ << ": ERROR: Unknown foldTimeReversal " << inikeys.inputPropsDB.foldTimeReversal
                 << " - foldTimeReversal must be \"none\", \"auto_majority\", or \"expected\" " << endl;
            exit(1);
          }

        // Check avgMode value
        if((inikeys.inputPropsDB.avgMode != "none") && (inikeys.inputPropsDB.avgMode != "lorentz") && (inikeys.inputPropsDB.avgMode != "spin"))
          {
            cerr << __func__ << ": ERROR: Unknown avgMode " << inikeys.inputPropsDB.avgMode
                 << " - avgMode must be \"none\", \"lorentz\", or \"spin\" " << endl;
            exit(1);
          }

        // Check rephaseMode value
        if((inikeys.inputPropsDB.rephaseMode != "auto") && (inikeys.inputPropsDB.rephaseMode != "auto_positive") && (inikeys.inputPropsDB.rephaseMode != "real"))
          {
            cerr << __func__ << ": ERROR: Unknown rephaseMode " << inikeys.inputPropsDB.rephaseMode
                 << " - rephaseMode must be one of \"auto\", \"auto_positive\", or \"real\" " << endl;
            exit(1);
          }

        twoPoints.loadRephaseComplexCorrs(inikeys.inputPropsDB.dbFname,
                                          inikeys.inputPropsDB.opsListFname,
                                          inikeys.inputPropsDB.irrepDim,
                                          inikeys.inputPropsDB.keys,
                                          inikeys.inputPropsDB.avgMode,
                                          inikeys.inputPropsDB.avgTol,
                                          inikeys.inputPropsDB.badList,
                                          inikeys.inputPropsDB.foldTimeReversal,
                                          inikeys.inputPropsDB.rephaseMode,
                                          inikeys.inputPropsDB.avgMom,
                                          inikeys.inputPropsDB.momListFname);
      }

    else if((inikeys.dbInputType == "redstar") || (inikeys.dbInputType == "redstar_debug"))
      {
        cout << "Using filedb (redstar) database: " << inikeys.inputPropsRedstar.dbFname << endl;

        // Check foldTimeReversal value
        if((inikeys.inputPropsRedstar.foldTimeReversal != "none") && (inikeys.inputPropsRedstar.foldTimeReversal != "auto_majority") && (inikeys.inputPropsRedstar.foldTimeReversal != "expected"))
          {
            cerr << __func__ << ": ERROR: Unknown foldTimeReversal " << inikeys.inputPropsRedstar.foldTimeReversal
                 << " - foldTimeReversal must be \"none\", \"auto_majority\", or \"expected\" " << endl;
            exit(1);
          }

        // Check rephaseMode value
        if((inikeys.inputPropsRedstar.rephaseMode != "auto") && (inikeys.inputPropsRedstar.rephaseMode != "auto_positive") && (inikeys.inputPropsRedstar.rephaseMode != "real"))
          {
            cerr << __func__ << ": ERROR: Unknown rephaseMode " << inikeys.inputPropsRedstar.rephaseMode
                 << " - rephaseMode must be one of \"auto\", \"auto_positive\", or \"real\" " << endl;
            exit(1);
          }

        twoPoints.loadRephaseComplexCorrs(inikeys.inputPropsRedstar.dbFname,
                                          inikeys.inputPropsRedstar.opsListFname,
                                          inikeys.inputPropsRedstar.opsXMLFiles,
                                          inikeys.inputPropsRedstar.avgRows,
                                          inikeys.inputPropsRedstar.avgTol,
                                          inikeys.inputPropsRedstar.badList,
                                          inikeys.inputPropsRedstar.foldTimeReversal,
                                          inikeys.inputPropsRedstar.rephaseMode,
                                          inikeys.inputPropsRedstar.avgMom,
                                          inikeys.inputPropsRedstar.momListFname,
                                          inikeys.inputPropsRedstar.KeyParams);
      }

    else
      {
        cerr << __func__ << ": ERROR: Unknown inputtype " << inikeys.dbInputType
             << " - inputtype must be \"dbnew\", \"redstar\", \"ensem\", \"ensem_onecorr\" or, \"..._debug\" " << endl;
        exit(1);
      }


    // If skip_nt != 0, only use every n'th timeslice [new Lt = floor(Lt/nt)]
    if (inikeys.globalProps.skip_nt != 0)
      {
	cout << "Only using every n'th timeslice where n= " << inikeys.globalProps.skip_nt << endl;
	twoPoints.skipTimeslices(inikeys.globalProps.skip_nt);
      }
    
    
    // Debugging
    if(inikeys.dbInputType == "ensem_debug" || inikeys.dbInputType == "dbnew_debug" || inikeys.dbInputType == "redstar_debug")
      {
        ofstream debug("debug.log");

        if(inikeys.dbInputType == "redstar_debug")
          {
            debug << "***** WARNING *****\n";
            debug << "These correlators have been divided by 12 before printing out to make them consistent with old corrs\n";
          }

        for(int ttest = 1; ttest < twoPoints.getLt(); ttest++)
          {
            SembleMatrix<double> testMat = twoPoints.getCt(ttest);

            if(inikeys.dbInputType == "redstar_debug") testMat = testMat / Real(12.0);  // HACK TO MAKE CONSISTENT WITH OLD FORMAT

            for(int i = 0; i < twoPoints.getDim(); i++)
              {
                for(int j = 0; j < twoPoints.getDim(); j++)
                  {
                    debug << "t = " << ttest << " (" << i << "," << j << "): " << mean(testMat.getEnsemElement(i, j)) << " , " << sqrt(variance(testMat.getEnsemElement(i, j))) << endl;
                  }
              }
          }

        debug.close();
        cout << "Debugging mode (" << inikeys.dbInputType << "): written out real correlators (mean and sd) to debug.log" << endl;
      }

    // If expWeight == true, consider exp(Et)*C(t)
    if(inikeys.weightProps.weight)
      {
        cout << "Using exp(Et)*C(t) where E= " << inikeys.weightProps.E << endl;
        twoPoints.expWeightCorrs(inikeys.weightProps.E);
      }
    else
      {
        cout << "Not exp Weighting correlators" << endl;
      }


    // If shift_dt != 0, consider C(t) - C(t+dt)
    if(inikeys.shiftProps.shift)
      {
        cout << "Using C(t) - C(t+dt) where dt= " << inikeys.shiftProps.dt << endl;
        twoPoints.shiftCorrs(inikeys.shiftProps.dt);
      }
    else
      {
        cout << "Not shifting correlators" << endl;
      }




  }

// **************************************************************************
// Various SembleCCorrs utils

  SembleMatrix<std::complex<double> > SembleCCorrs::getCt(int t) const
  {
    return Ct[t];
  }

  int SembleCCorrs::getLt() const
  {
    return Lt;
  }
  int SembleCCorrs::getDim() const
  {
    return dim;
  }
  int SembleCCorrs::getBins() const
  {
    return nbins;
  }

  EnsemVectorComplex SembleCCorrs::getCij(int i, int j) const
  {
    EnsemVectorComplex dum;
    dum.resize(Ct[0].bins());
    dum.resizeObs(Lt);

    for(int t = 0; t < Lt; t++)
      {
        pokeObs(dum, Ct[t].getEnsemElement(i, j), t);
      }

    return dum;
  }


// **************************************************************************
// Various utils

  bool phaseCorr(SembleMatrix<std::complex<double> > C, int row, int &phase_re, int  &phase_im, int refrow, int ref_phase_re, int ref_phase_im)
  {
    EnsemComplex c = C.getEnsemElement(refrow, row);
    EnsemComplex cT = C.getEnsemElement(row, refrow);
    double phase = toDouble(atan2(mean(imag(c - cT)) , mean(real(c + cT))));
    phase_re = int(round(cos(phase)) * ref_phase_re - round(sin(phase)) * ref_phase_im);
    phase_im = int(round(cos(phase)) * ref_phase_im + round(sin(phase)) * ref_phase_re);

    bool worked = true;

    if(((phase_re != 0) && (phase_im != 0)) || ((phase_re == 0) && (phase_im == 0)))
      {
        worked = false;
      }

    return worked;
  }

  void read(XMLReader &xml, const std::string &path, opsxml_t &param)
  {
    XMLReader paramtop(xml, path);

    if(paramtop.count("OpNameKey") > 0)
      read(paramtop, "OpNameKey", param.opnamekey);
    else
      param.opnamekey = "";

    if(paramtop.count("Irrep") > 0)
      read(paramtop, "Irrep", param.irrep);
  }


// **************************************************************************
// Various db utils (a lot of these are modified version of dbutils.cc)

//! Get a key/value
// New database format
  template<typename K, typename V>
  V printKeyValue(const K &ky,
                  FILEDB::AllConfStoreDB< SerialDBKey<K>,  SerialDBData<typename EnsemScalar<V>::Type_t> >& database)
  {
    typedef typename EnsemScalar<V>::Type_t SV;

    SerialDBKey<K> key;
    key.key() = ky;

    std::vector< SerialDBData<SV> > vals;
    int ret;

    if((ret = database.get(key, vals)) != 0)
      {
        std::cerr << __func__ << ": key not found\n" << ky;
        exit(1);
      }

    V eval;
    eval.resize(vals.size());
    eval.resizeObs(vals[0].data().numElem());

    for(int i = 0; i < vals.size(); ++i)
      {
        SV sval = vals[i].data();
        pokeEnsem(eval, sval, i);
      }

    return eval;
  }

}
