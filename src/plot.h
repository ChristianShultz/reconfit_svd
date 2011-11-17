#ifndef __PLOT_H__
#define __PLOT_H__

#include <ostream>
#include <iostream>
#include <vector>
#include "ensem/ensem.h"
#include "ensem_data.h"

using namespace std;
using namespace ENSEM;

class AxisPlot{	
 public:
  AxisPlot();
	
  void addEnsemData(vector<int> t, EnsemVectorReal data, string style, int colour);
  void addEnsemData(vector<double> x, EnsemVectorReal data, string style, int colour);
  void addEnsemData(EnsemVectorReal data, string style, int colour);
  void addEnsemData(vector<double> x, vector<EnsemReal> y, string style, int colour);
    
  void addLineData(vector<int> t, vector<double> data, int colour);
  void addLineData(vector<double> x, vector<double> y, int colour);
  
  void addErrorData(vector<int> t, vector<double> avg, vector<double> err, string style, int colour);
  void addErrorData(vector<double> x, vector<double> avg, vector<double> err, string style, int colour);
  void addErrorData(vector<double> avg, vector<double> err, string style, int colour);
	
  void setXRange(int tmin, int tmax);
  void setXRange(double xmin, double xmax);
  void setYRange(double ymin, double ymax);

  void addLabel(double xpos, double ypos, string label, int colour, double size);
    
  
  void addEnsemData(EnsemData& data);  //adds ensemData hidden and shown
 	
  string getAxisPlotString();
  void sendToFile(string filename);

private:
  stringstream plot;
  //  ofstream outfs;
  
};

#endif
