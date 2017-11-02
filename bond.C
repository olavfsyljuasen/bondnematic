#include<string>
#include<iostream>
#include<cassert>
using namespace std;

#ifndef M_PI
#define  M_PI 3.141592653589793238462643
#endif
// what model to consider:


// universal definitions
const double PI=M_PI;
const double TWOPI=2.*PI;

// Options for this run
#if defined NDEBUG 
const bool TRACE= false; 
#else
const bool TRACE= true; // debugging option
#endif


const int  MAXRUNS = 200; // the maximum number of runs in the read.in file  
const int  MAXPAR = 100;  // the maximum allowed input parameters


const string RESULTSNAME= "res.dat";
const string RESULTS2NAME= "DeltavsT.dat";
const string RESULTS3NAME= "DeltaoverVvsT.dat";

const string PARAMETERFILENAME = "Deltas.in";

const string READIN   = "read.in";

#include <fstream>
ofstream logfile("log.txt",ios::app);

//#include "rnddef.h" // the random number generator RAN


#define LATTICETYPE BravaisLattice


#include "globalheader.h"

#include "J123rules.h"

#include "RunParameter.h"

#include "bravaislattices.h"

#include "bond.h"



int main()
{
  if(TRACE) cout << "Starting Main()" << endl;
  Timer mytimer;  
  logfile << "\nStarting at: " << mytimer << endl;
  mytimer.Start();
  
  RunParameters rp;
  logfile << "parameters: " << rp;

  int nc = rp.GetNR();

  for(int ic=1; ic <= nc; ic++)
    {
      double* par = rp.GetPars(ic);
      Simulation<LATTICETYPE> sim(par,ic);

      logfile << "Starting calculation " << ic << endl;
      sim.Run();
    }
            
  mytimer.Stop();
  double time_spent = mytimer.GetTimeElapsed();
  logfile << "Time elapsed:  " << time_spent << " sec. :" 
	  << time_spent/3600. << " hours." << endl;
  logfile << "Ending at: " << mytimer;
  logfile << "---Done--- " << endl;
  logfile.close();
  
  ofstream terminationsign("end.exec");
  terminationsign << "execution ended at " << mytimer;
  terminationsign.close();
  
  cout << "program ended \n";
  exit(0);
}










