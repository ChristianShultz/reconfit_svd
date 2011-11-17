#include "plot.h"

AxisPlot::AxisPlot(){}

void AxisPlot::addEnsemData(vector<int> t, EnsemVectorReal data, string style, int colour){
  //write the ensem data
  plot << "#e" << endl << "#c " << style << endl << "#cs 0.5" << endl << "#m " << colour <<endl;
   
  for(int it=0; it < t.size(); it++){
    plot << t[it] << "  " << toDouble( mean( peekObs(data , it) ) ) 
	  << "  " <<  toDouble( sqrt( variance( peekObs(data , it) ) ) )  << endl;
  }
  plot << endl << endl << "#" << endl;
}

void AxisPlot::addEnsemData(vector<double> x, EnsemVectorReal data, string style, int colour){
  //write the ensem data
  plot << "#e" << endl << "#c " << style << endl << "#cs 0.5" << endl << "#m " << colour <<endl;
   
  for(int it=0; it < x.size(); it++){
    plot << x[it] << "  " << toDouble( mean( peekObs(data , it) ) ) 
	  << "  " <<  toDouble( sqrt( variance( peekObs(data , it) ) ) )  << endl;
  }
  plot << endl << endl << "#" << endl;
}

void AxisPlot::addEnsemData(EnsemVectorReal data, string style, int colour){
  vector<int> t;
  for(int it = 0; it < data.numElem(); it++){
    t.push_back(it);
  }
  addEnsemData(t, data, style, colour);
}

void AxisPlot::addEnsemData(vector<double> x, vector<EnsemReal> y, string style, int colour){
  EnsemVectorReal yy; yy.resize(y[0].size()); yy.resizeObs(y.size());
  for(int it = 0; it < y.size(); it++){ pokeObs(yy, y[it], it); }
  addEnsemData(x, yy, style, colour);
}

void AxisPlot::addLineData(vector<int> t, vector<double> data, int colour){
  plot << "#e0" << endl << "#c0" << endl << "#m " << colour << endl;
  for(int it=0; it < t.size(); it++){
    plot<< t[it] << "  " << data[it]  << endl;
  }
  plot << endl<< endl << "#" << endl;	
}

void AxisPlot::addLineData(vector<double> x, vector<double> y, int colour){
  plot << "#e0" << endl << "#c0" << endl << "#m " << colour << endl;
  for(int it=0; it < x.size(); it++){
    plot<< x[it] << "  " << y[it]  << endl;
  }
  plot << endl<< endl << "#" << endl;	
}

void AxisPlot::addErrorData(vector<int> t, vector<double> avg, vector<double> err, string style, int colour){
  plot << "#e" << endl << "#c " << style << endl << "#cs 0.5" << endl << "#m " << colour <<endl;
  
  for(int it=0; it < t.size(); it++){
    plot << t[it] << "  " << avg[it] << "  " << err[it] << endl;
  }
  plot << endl << endl <<"#" << endl;
}

void AxisPlot::addErrorData(vector<double> x, vector<double> avg, vector<double> err, string style, int colour){
  plot << "#e" << endl << "#c " << style << endl << "#cs 0.5" << endl << "#m " << colour <<endl;
  
  for(int it=0; it < x.size(); it++){
    plot << x[it] << "  " << avg[it] << "  " << err[it] << endl;
  }
  plot << endl << endl <<"#" << endl;
}

void AxisPlot::addErrorData(vector<double> avg, vector<double> err, string style, int colour){
  vector<int> t;
  for(int it = 0; it < avg.size(); it++){
    t.push_back(it);
  }
  addErrorData(t, avg, err, style, colour);
}


void AxisPlot::setXRange(int tmin, int tmax){plot << "#x " << tmin << " " << tmax << endl;}
void AxisPlot::setXRange(double xmin, double xmax){plot << "#x " << xmin << " " << xmax << endl;}
void AxisPlot::setYRange(double ymin, double ymax){plot << "#y " << ymin << " " << ymax << endl;}

void AxisPlot::addLabel(double xpos, double ypos, string label, int colour, double size)
{
  plot << "#e0" << endl << "#c0" << endl << "#m " << colour << endl << "#cs " << size << endl ;
  plot << xpos <<" " << ypos <<"\"" << label << "\"" << endl;
}

void AxisPlot::addEnsemData(EnsemData& data){
  //get shown data
  vector<double> x = data.getXData();
  EnsemVectorReal y = data.getYData();
  addEnsemData(x, y, "\\sq", 1);
  
  //get hidden data
  vector<double> all_x = data.getAllXData(); vector<double> hidden_x;
  EnsemVectorReal all_y = data.getAllYData(); vector<EnsemReal> hidden_y;
  vector<bool> active_data = data.getActiveDataList();
  for(int i=0; i < all_x.size(); i++){
    if(!active_data[i]){
      hidden_x.push_back(all_x[i]);
      hidden_y.push_back( peekObs(all_y , i) );
    }
  } 
  addEnsemData(hidden_x, hidden_y, "\\di", 3);
}

string AxisPlot::getAxisPlotString(){
  return plot.str();
} 

void AxisPlot::sendToFile(string filename){
  ofstream outfs; outfs.open(filename.c_str());
  outfs << plot.str();
  outfs.close();
}

