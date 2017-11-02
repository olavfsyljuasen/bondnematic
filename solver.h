#ifndef SOLVER_H
#define SOLVER_H

#include<vector>
#include<complex>
#include<iostream>

extern "C"
{
  //#g++ -I/data/sylju/include -L/data/sylju/lib FFTWtest.C -lfftw3 -lm;
#include <fftw3.h>
}

/*
void ComputeSelfEnergy(vector<double>& K1,int NS,int dim, vector<int>& dims,vector<double>& SigmaE)
{
  if(TRACE) cout << "Starting ComputeSelfEnergy Delta=" << endl;

  if(TRACE) cout << "K1: " << K1 << endl;

  int Vq=1;
  int Vf=1;
  for(int i=0; i<(dim-1); i++){Vq*=dims[i]; Vf*=dims[i];}
  Vq*= dims[dim-1];
  Vf*= dims[dim-1]/2+1;


  const double invVq=1./Vq;
  
  vector<double> K1inv(Vq);
  vector<complex<double> > K1invtilde(Vf);

  fftw_plan p1 = fftw_plan_dft_r2c(dim,&dims[0],&K1inv[0],reinterpret_cast<fftw_complex*>(&K1invtilde[0]),FFTW_ESTIMATE);

  for(int i=0; i<Vq; i++){K1inv[i]=invVq/K1[i];}
  fftw_execute(p1); // K1inv -> K1invtilde

  if(TRACE)
    { 
      cout << "K1invtilde" << endl;
      FourierPrint(Vq,K1invtilde);
    }

  vector<complex<double> > K1K1kernel(Vf);
  vector<double> D1(Vq);
  fftw_plan p2 = fftw_plan_dft_c2r(dim,&dims[0],reinterpret_cast<fftw_complex*>(&K1K1kernel[0]),&D1[0],FFTW_ESTIMATE);

  vector<complex<double> > D1tilde(Vf);
  fftw_plan p3 = fftw_plan_dft_r2c(dim,&dims[0],&D1[0],reinterpret_cast<fftw_complex*>(&D1tilde[0]),FFTW_ESTIMATE);

  for(int i=0; i<Vf; i++) K1K1kernel[i] = K1invtilde[i]*conj(K1invtilde[i]);

  if(TRACE)
    {
      cout << "K1K1kernel" << endl;
      FourierPrint(Vq,K1K1kernel);
    }

  fftw_execute(p2);// K1K1kernel->D1 

  const double factor=(2./NS)*invVq; // FFTW (unnormalized backtransform) 
  for(int i=0; i<Vq; i++) D1[i] = factor/D1[i];

  
  if(TRACE) cout << "D1: " << D1 << endl;
  

  fftw_execute(p3);// D1->D1tilde

  if(TRACE)
    {
      cout << "D1tilde" << endl;
      FourierPrint(Vq,D1tilde);
    }

  vector<complex<double> > D1K1kernel(Vf);
  fftw_plan p4 = fftw_plan_dft_c2r(dim,&dims[0],reinterpret_cast<fftw_complex*>(&D1K1kernel[0]),&SigmaE[0],FFTW_ESTIMATE);
  
  for(int i=0; i<Vf; i++) D1K1kernel[i]=(D1tilde[i]*conj(K1invtilde[i])).real();

  if(TRACE)
    {
      cout << "D1K1kernel" << endl;
      FourierPrint(Vq,D1K1kernel);
    }

  fftw_execute(p4);// D1K1kernel -> SigmaE

  if(TRACE) cout << "SigmaE: " << SigmaE << endl;

  if(TRACE) cout << "Done with ComputeSelfEnergy " << endl;
}
*/

/* saving some space */
void ComputeSelfEnergy(vector<double>& K1,int NS,int dim, vector<int>& dims,vector<double>& SigmaE)
{
  if(TRACE) cout << "Starting ComputeSelfEnergy2 Delta=" << endl;

  if(TRACE) cout << "K1: " << K1 << endl;

  int Vq=1;
  int Vf=1;
  for(int i=0; i<(dim-1); i++){Vq*=dims[i]; Vf*=dims[i];}
  Vq*= dims[dim-1];
  Vf*= dims[dim-1]/2+1;


  const double invVq=1./Vq;
  
  //construct the temporary arrays
  vector<double> T1(Vq);
  vector<complex<double> > T1f(Vf);
  vector<complex<double> > T2f(Vf);

  //construct the Fourier plans:
  fftw_plan T1T1f = fftw_plan_dft_r2c(dim,&dims[0],&T1[0],reinterpret_cast<fftw_complex*>(&T1f[0]),FFTW_ESTIMATE);
  fftw_plan T1T2f = fftw_plan_dft_r2c(dim,&dims[0],&T1[0],reinterpret_cast<fftw_complex*>(&T2f[0]),FFTW_ESTIMATE);
  fftw_plan T1fT1 = fftw_plan_dft_c2r(dim,&dims[0],reinterpret_cast<fftw_complex*>(&T1f[0]),&T1[0],FFTW_ESTIMATE);
  fftw_plan T1fSigmaE = fftw_plan_dft_c2r(dim,&dims[0],reinterpret_cast<fftw_complex*>(&T1f[0]),&SigmaE[0],FFTW_ESTIMATE);

  // fill in T1 with K1inv
  for(int i=0; i<Vq; i++){T1[i]=invVq/K1[i];}

  fftw_execute(T1T2f); // K1inv -> K1invtilde

  if(TRACE)
    { 
      cout << "K1invtilde" << endl;
      FourierPrint(Vq,T2f);
    }

  // fill in T1f with K1K1kernel
  for(int i=0; i<Vf; i++) T1f[i] = T2f[i]*conj(T2f[i]);

  if(TRACE)
    {
      cout << "K1K1kernel" << endl;
      FourierPrint(Vq,T1f);
    }


  fftw_execute(T1fT1); // K1K1kernel -> D1

  // invert D1
  const double factor=(2./NS)*invVq; 
  for(int i=0; i<Vq; i++) T1[i] = factor/T1[i];

  
  if(TRACE) cout << "D1: " << T1 << endl;

  fftw_execute(T1T1f);// D1->D1tilde

  if(TRACE)
    {
      cout << "D1tilde" << endl;
      FourierPrint(Vq,T1f);
    }

  // construct D1K1kernel in T1f
  for(int i=0; i<Vf; i++) T1f[i]=(T1f[i]*conj(T2f[i])).real();

  if(TRACE)
    {
      cout << "D1K1kernel" << endl;
      FourierPrint(Vq,T1f);
    }

  fftw_execute(T1fSigmaE);// D1K1kernel -> SigmaE

  if(TRACE) cout << "SigmaE: " << SigmaE << endl;

  if(TRACE) cout << "Done with ComputeSelfEnergy " << endl;
}


#endif //SOLVER_H
