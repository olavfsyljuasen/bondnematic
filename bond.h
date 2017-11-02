#ifndef BOND_H
#define BOND_H

#include<vector>
#include<complex>
#include<fstream>
#include<iostream>
#include<fstream>

// Routine for pretty printing of matrices, assuming a square matrix
ostream& operator<<(ostream& os,vector<double>& M)
{
  const int L=int(sqrt(M.size()));
  const int pcols=(L>3 ? 3:L);
  const int prows=(L>2 ? 2:L);

  os << "(" << L << " X " << L << "):" << endl;
  for(int j=0; j<prows; j++)
    {
      for(int i=0; i<pcols; i++) os << M[j*L+i] << " ";
      if(L>pcols)
	os << "..." << M[j*L+L-1] << endl;
      else
	os << endl;
    }
  if(L > prows) os << "...." << endl;
  for(int i=0; i<pcols; i++) os << M[(L-1)*L+i] << " ";
  if(L>pcols)
    os << "..." << M[L*L-1];
  return os;
}

// Routine for printing fourier matrices
void FourierPrint(int N,vector<complex<double> >& M)
{
  const int L=int(sqrt(N));
  const int LX=L/2+1;
  const int pcols=(LX>20? 20:LX);
  const int prows=(L>3 ? 2:L);


  cout << "(" << L << " X " << LX << "):" << endl;
  for(int j=0; j<prows; j++)
    {
      for(int i=0; i<pcols; i++) cout << setprecision(14) << scientific << M[j*LX+i] << " ";
      if(LX>pcols)
	cout << "..." << M[j*LX+LX-1] << endl;
      else
	cout << endl;
    }
  if(L > prows) cout << "...." << endl;
  for(int i=0; i<pcols; i++) cout << M[(L-1)*LX+i] << " ";
  if(LX>pcols)
    cout << "..." << M[L*LX-1];

  cout << endl;
}



void SubtractMinimum(vector<double>& K)
{
  if(TRACE) cout << "Start SubtractMinimum" << endl;
  double min=1.E99; // a really big number
  
  for(int i=0; i<K.size(); i++){if(K[i]<min){min=K[i];}} // find minimum
  for(int i=0; i<K.size(); i++){K[i]-=min;}   // subtract it
  if(TRACE) cout << "Done SubtractMinimum" << endl;
}





#include "solver.h"


template<class LATTICE>
class Driver
{
 public:
  Driver(double*,LATTICE&,Rule&);
  double CalculateT();
  vector<double> CalculateOrderPars(double);
  void SolveSelfConsistentEquation(vector<double>& dK0,double Delta);

  void Solve(const double delta)
  {
    if(TRACE) cout << "Starting Solve with Delta= " << delta << endl;
    Delta=delta;
    vector<double>& dK0=rule.GetInitialState();
    SolveSelfConsistentEquation(dK0,Delta);
    if(TRACE) cout << "Done Solve" << endl;
  }



 private:
  double* param;
  LATTICE& lattice; 
  Rule& rule;

  int Vf; // number of Fourier sites
  int Vq; // number of q-space sites.

  double Delta; 

  vector<double> K0; // K0(q) ; energy function 
  vector<double> K1; // K1(q)
  vector<double> SigmaE; // SigmaE(q), the self-energy

};

template<class LATTICE>
Driver<LATTICE>::Driver(double* pars,LATTICE& la,Rule& r):param(pars),lattice(la),rule(r),Vq(la.SiterVol()),Vf(la.SiteqVol()),Delta(0), K0(Vq),K1(Vq),SigmaE(Vq)
{
  if(TRACE) cout << "Initializing solver " << endl;
  K0=rule.GetInteraction();
  SubtractMinimum(K0);  

  if(TRACE) cout << "Done initializing solver " << endl;
}



template<class LATTICE>
double Driver<LATTICE>::CalculateT()
{
  if(TRACE) cout << "Starting CalculateT " << endl;
  double sum=0.;
  for(int i=0; i<Vq; i++){sum+= 1./K1[i];}
  if(TRACE) cout << "Done CalculateT " << endl;
  return (2.*Vq)/(NS*sum);
}

template<class LATTICE>
vector<double> Driver<LATTICE>::CalculateOrderPars(double T)
{
  if(TRACE) cout << "Starting CalculateOrderPars" << endl;

  vector<double> opars(NOBSERVABLES);
  for(int j=0; j<NOBSERVABLES; j++)
    {
      double sum=0.;
      vector<double>& f=rule.GetIrrep(j);

      for(int i=0; i<Vq; i++)
	{
	  sum+= f[i]/K1[i];
	}
      opars[j]=0.5*NS*T*sum/Vq;
    }

  if(TRACE) cout << "Done CalculateOrderPars" << endl;
  return opars;
}




template<class LATTICE>
void Driver<LATTICE>::SolveSelfConsistentEquation(vector<double>& dK0, const double Delta)
{
  if(TRACE) cout << "Starting SolveSelfConsistentEquation with Delta=" << Delta << endl;

  if(TRACE) cout << "K0: " << K0 << endl;
  if(TRACE) cout << "dK0: " << dK0 << endl;


  for(int i=0; i<Vq; i++){K1[i]= K0[i]+dK0[i];}
  SubtractMinimum(K1);
  for(int i=0; i<Vq; i++){K1[i]+=Delta;}


  double oldT=-1.; // something unphysical
  double newT=-1.;
  int iter=0;
  bool converged=false;
  while(iter< param[MAXITER] && !converged)
    {
      int dim=lattice.D();
      vector<int> dims=lattice.SiterDims();
      ComputeSelfEnergy(K1,NS,dim,dims,SigmaE);

      if(TRACE) cout << "SigmaE " << SigmaE << endl;

      for(int i=0; i<Vq; i++){K1[i]= K0[i]+SigmaE[i];}      
      SubtractMinimum(K1);
      for(int i=0; i<Vq; i++){K1[i] +=Delta;}      

      if(TRACE) cout << "Final K1: " << K1 << endl;


      iter++;
      // convergence checks
      newT=CalculateT();
      if( fabs((newT-oldT)/oldT) < param[TOLERANCE]) converged=true;
      oldT=newT;

      vector<double> nobs=CalculateOrderPars(newT);


      logfile << iter << " T=" << newT << " ";
      for(int i=0; i<NOBSERVABLES; i++)
	{
	  logfile << nobs[i] << " ";
	} 
      logfile << " converged: " << (converged ? "true": "false") << endl;
    }
  if(!converged)
    {
      logfile << "reached MAXITER=" << param[MAXITER] << " iterations without converging, increase MAXITER!" << endl;
    }
  else
    {
      logfile << "Convergence reached after " << iter << " steps." << endl;
      vector<double> result=CalculateOrderPars(newT);
      ofstream outfile(RESULTSNAME.c_str(),ios::app);
      outfile << Delta << " " << Delta/lattice.Nsites() << " " << newT << " "; 
      for(int j=0; j<NOBSERVABLES; j++) outfile << result[j] << " ";
      outfile << endl;
      outfile.close();

      ofstream outfile2(RESULTS2NAME.c_str(),ios::app);
      outfile2 << Delta << " " << " " << newT << " "; 
      outfile2 << endl;
      outfile2.close();

      ofstream outfile3(RESULTS3NAME.c_str(),ios::app);
      outfile3 << Delta/lattice.Nsites() << " " << " " << newT << " "; 
      outfile3 << endl;
      outfile3.close();
    }

  if(TRACE) cout << "Done SolveSelfConsistentEquation " << endl;
}


template<class LATTICE>
class Simulation{
  friend ostream& operator<<(ostream& os,Simulation& s){
    os << endl; return os;}
 public:
  Simulation(double* pars,int i);
  void Run();
 private:
  double* param;
  int ic;
  LATTICE lattice;
  Rule rule;

  Driver<LATTICE> mysolver;
  vector<double> Deltalist;
};

template<class LATTICE>
Simulation<LATTICE>::Simulation(double* pars,int i): param(pars),ic(i),lattice(pars),rule(pars,lattice),mysolver(pars,lattice,rule),Deltalist(0)
{
  if(TRACE) cout << "Initializing Simulation" << endl;

  ifstream parameterfile(PARAMETERFILENAME.c_str());
  if(!parameterfile)
    {
      if(TRACE) 
	cout << "No file " << PARAMETERFILENAME << " found." 
	     << " Using Delta=" << param[DELTA] << endl;
      logfile << "No file " << PARAMETERFILENAME << " found." 
	      << " Using Delta=" << param[DELTA] << endl;
      Deltalist.push_back(param[DELTA]); 
    }
  else
    {
      if(TRACE) cout << "Reading " << PARAMETERFILENAME << " from disk" << endl;
      logfile << "Reading " << PARAMETERFILENAME << " from disk" << endl;
      while(parameterfile)
	{
	  double dvalue; 
	  parameterfile >> dvalue;
	  if(parameterfile)
	    Deltalist.push_back(dvalue);
	}
    }
  if(TRACE) cout << "Deltalist has " << Deltalist.size() << " entries" << endl; 
  if(TRACE) cout << "Done Initializing Simulation" << endl;
}

template<class LATTICE>
void Simulation<LATTICE>::Run()
{
  if(TRACE) cout << "Starting Run" << endl;
  for(int i=0; i< Deltalist.size(); i++)
    {
      mysolver.Solve(Deltalist[i]);
    }
  if(TRACE) cout << "Done Run" << endl;
}
  


#endif //BOND_H
