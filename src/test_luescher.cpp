#include "luescher.h"

using namespace std;

int main(int argc, char** argv){

  /*  LuescherW00 w00;

  for(int i=0; i < 100; i++){
    double qsq = 0.0005 + i/100.0;
    cout << "qsq = " << qsq << ", w00 = " << w00(qsq) << endl;
    }*/

  LuescherPhiSigma phi_sigma;

  string irrep = "A1p";
  for(int i=0; i < 200; i++){
    double ksq = 0.0001 + (0.045 / 200)*i;
    double qsq = ksq * pow(20*3.5/(2*3.14159)   , 2);
    double sigma = phi_sigma.sigma(irrep, sqrt(qsq));
    double t4 = -4.6e5 * pow( sqrt(ksq) , 9);
    cout << ksq << " " << (180/3.14159)*sigma*t4 << endl;
  }


/*
  string irrep = "A1p";

  for(int i=0; i < 10; i++){
    double qsq = 0.0001 + (0.1 / 10) * i;
    cout << "qsq = " << qsq << ", phi = " << phi_sigma.phi(irrep, sqrt(qsq)) << endl;		
  }

  cout << endl << endl; 	

  for(int i=0; i < 10; i++){
    double qsq = 0.1 + 0.1*i;
    double phi = phi_sigma.phi(irrep, sqrt(qsq));
    if(phi < 0){phi += 3.14159;}
    cout << "qsq = " << qsq << ", phi/(pi*qsq) = " << phi / (3.14159 * qsq) << endl;					
  } 	

  cout << endl; 	

  for(int i=0; i < 10; i++){
    double qsq = 1.1 + 0.1*i;
    double phi = phi_sigma.phi(irrep, sqrt(qsq));
    if(phi < 0){phi += 3.14159;}
    if(phi < 3.14159){phi += 3.14159;}	
    cout << "qsq = " << qsq << ", phi/(pi*qsq) = " << phi / (3.14159 * qsq) << endl;
  }

 cout << endl;

  for(int i=0; i < 10; i++){
    double qsq = 2.1 + 0.1*i;
    double phi = phi_sigma.phi(irrep, sqrt(qsq));
    if(phi < 0){phi += 3.14159;}
    if(phi < 3.14159){phi += 3.14159;}
    if(phi < 2*3.14159){phi += 3.14159;}	
    cout << "qsq = " << qsq << ", phi/(pi*qsq) = " << phi / (3.14159 * qsq) << endl;
  }

 cout << endl;

  for(int i=0; i < 10; i++){
    double qsq = 3.1 + 0.1*i;
    double phi = phi_sigma.phi(irrep, sqrt(qsq));                     
    if(phi < 0){phi += 3.14159;}
    if(phi < 3.14159){phi += 3.14159;}
    if(phi < 2*3.14159){phi += 3.14159;}
    if(phi < 3*3.14159){phi += 3.14159;}	
    cout << "qsq = " << qsq << ", phi/(pi*qsq) = " << phi / (3.14159 * qsq) << endl;
  }

 cout << endl;
    
  for(int i=0; i < 10; i++){
    double qsq = 4.1 + 0.1*i;
    double phi = phi_sigma.phi(irrep, sqrt(qsq));
    if(phi < 0){phi += 3.14159;}
    if(phi < 3.14159){phi += 3.14159;}
    if(phi < 2*3.14159){phi += 3.14159;}
    if(phi < 3*3.14159){phi += 3.14159;} 
    if(phi < 4*3.14159){phi += 3.14159;}
    cout << "qsq = " << qsq << ", phi/(pi*qsq) = " << phi / (3.14159 * qsq) << endl;
  }

 cout << endl;

  for(int i=0; i < 10; i++){
    double qsq = 5.1 + 0.1*i;
    double phi = phi_sigma.phi(irrep, sqrt(qsq));
    if(phi < 0){phi += 3.14159;}
    if(phi < 3.14159){phi += 3.14159;}
    if(phi < 2*3.14159){phi += 3.14159;}
    if(phi < 3*3.14159){phi += 3.14159;}
    if(phi < 4*3.14159){phi += 3.14159;}
    if(phi < 5*3.14159){phi += 3.14159;}
    cout << "qsq = " << qsq << ", phi/(pi*qsq) = " << phi / (3.14159 * qsq) << endl;
  }
*/

/*
  cout << endl << endl;
  for(int i=0; i < 500; i++){
    double qsq = 0.0001 + (6.0/500)*i;
    double sigma = phi_sigma.sigma(irrep, sqrt(qsq));
    cout << qsq << " " << sigma << endl;
  }
 */

/*for(int i=0; i < 600; i++){
  double qsq = 0.001 + i / 100.0;
  string A1p = "A1p";
  string T2p = "T2p";
  string Ep = "Ep";
  string T1p = "T1p";
//  cout << qsq << " " << phi_sigma.phi(A1p, sqrt(qsq)) << " " << phi_sigma.phi(T2p, sqrt(qsq)) << " " << phi_sigma.phi(Ep, sqrt(qsq)) ;
//  cout << " " << phi_sigma.phi(T1p, sqrt(qsq)) << endl; 
  cout << qsq << " " << phi_sigma.sigma(A1p, sqrt(qsq)) << " " << phi_sigma.sigma(T2p, sqrt(qsq)) << " " << phi_sigma.sigma(Ep, sqrt(qsq)) << endl ;
}
*/
};


  /*  if(argc != 2){ cerr << "test_luescher <q> returns phi(q) and sigma(q)" << endl; exit(1);}
  double q; {istringstream a(argv[1]); a>> q;}

  LuescherPhi luPhi;

  cout << "q = " << q << ",  phi(q) = " << luPhi.phi(q) << ",  sigma(q) = " << luPhi.sigma(q) << endl;
  */
  /*  {

    double pi = 2.0*atan2(1.0,0.0);

    cout << "high values of q, check versus tabulated values : " << endl;

    LuescherPhi luPhi;
    for(int i = 1; i < 90; i++){
      double q2 = 0.1*i;
      cout << "q^2 = " << q2 << " , phi/(pi*q^2) = " << luPhi.phi(sqrt(q2)) / (pi*q2) << " by " << luPhi.type(sqrt(q2)) <<  endl;
    }
    
    cout << "low values of q, check versus the formula : " << endl;
    for(int i = 1; i <= 10; i++){
      double q2 = 0.01*i;
      cout << "q^2 = " << q2 << " , phi  = " << luPhi.phi(sqrt(q2)) << " by " << luPhi.type(sqrt(q2)) <<  endl;
    }

    cout << endl << endl;
    cout << "check values of sigma" << endl;
    for(int i = 1; i <=10; i++){
	double q2 = 7.0 + 0.1 * i;
	cout << "q^2 = " << q2 << " , sigma = " << luPhi.sigma(sqrt(q2)) << endl;
    }			
    cout << endl;
    cout << "q = 1.864, phi = " << luPhi.phi(1.864) << " , sigma = " << luPhi.sigma(1.864) << endl;

    
  }
  cout << endl << endl;*/



  /*  {
    LuescherPhaseShift luDelta(3.5*16, 0.148333, 0.148333);
    cout << "delta(.300455) = " << (360/(2.0 * 3.14159)) * luDelta.delta(0.300455, 0) << endl;
    cout << "delta(.382271) = " << (360/(2.0 * 3.14159)) * luDelta.delta(0.382271, 1) << endl;
    cout << "delta(.448683) = " << (360/(2.0 * 3.14159)) * luDelta.delta(0.448683, 2) << endl;
    cout << "delta(.498682) = " << (360/(2.0 * 3.14159)) * luDelta.delta(0.498682, 3) << endl;
    cout << "delta(.543619) = " << (360/(2.0 * 3.14159)) * luDelta.delta(0.543619, 4) << endl;
    
  }
  
  cout << endl << endl;
  */

  //  double L = 20.0;

  /*    {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi0";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*L, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 0);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

    }*/   
  /*   {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi1";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*L, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 1);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

    }   
 {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi2";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*L, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 2);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

  }    
  /*  {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi3";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*L, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 3);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

    }*/   
  /*      {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi4";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*L, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 4);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

    }   
 


  /* {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi2";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*16, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 2);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

    } */   
 /* {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi3";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*16, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 3);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

    }*/    
  /* {
    //load some masses
    EnsemReal m_pi;
    {ostringstream filename;  filename << "m_pi";  read(filename.str(), m_pi);}
    EnsemReal E_pipi;
    {ostringstream filename;  filename << "E_pipi4";  read(filename.str(), E_pipi);}
    
    LuescherPhaseShiftEnsem luDelta(3.5*16, m_pi, m_pi);
    EnsemReal delta = Real(360/(2.0 * 3.14159)) * luDelta.delta(E_pipi, 2);

    cout << "delta(" << toDouble(mean(E_pipi)) << " +/- " << toDouble(sqrt(variance(E_pipi))) << ") = " << toDouble(mean(delta)) << " +/- " << toDouble(sqrt(variance(delta))) << endl;

    } */   
//};

