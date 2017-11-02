#ifndef BRAVAISLATTICES_H
#define BRAVAISLATTICES_H


#include<complex>
#include<iomanip>
#include<iostream>


class BravaisLattice
{
 public:
  BravaisLattice(double* par);
  ~BravaisLattice(){if(TRACE) cout << "deleting BravaisLattice\n";}

  const int N1;
  const int N2;
  const int N3;

 private:
  //  const int A;
  const int v;
  const int d;

 public:

  int Nsites() const {return v;} 
  int D() const {return d;}


  Coord Pos(const int i){return site_r[i];}

  int SiterVol() const {return site_rvol;} 
  vector<int> SiterDims() const {return site_rdims;}
  int SiteqVol() const {return site_qvol;}
  vector<int> SiteqDims() const {return site_qdims;}

 private:
  const Coord a1;
  const Coord a2;
  const Coord a3;
  const Coord origo;

  Coord b1;
  Coord b2;
  Coord b3;

  //  vector<int> phase;

  const int site_rvol;
  vector<int> site_rdims;
  vector<Coord> site_r;

  const int site_qvol;
  vector<int> site_qdims;
  vector<Coord> site_q;

};



BravaisLattice::BravaisLattice(double* par):
  N1(static_cast<int>(par[NX])),
     N2(static_cast<int>(par[NY])),
     N3(static_cast<int>(par[NZ])),


     v(N1*N2*N3),
     d( (N1>1)+(N2>1)+(N3>1) ),

   // basis vectors:
 
#ifdef SIMPLECUBIC
    a1(1.,0.,0.),
    a2(0.,1.,0.),
    a3(0.,0.,1.),
#elif defined BODYCENTEREDCUBIC
    a1(1.,0.,0.),
    a2(0.,1.,0.),
    a3(0.5,0.5,0.5),
#elif defined FACECENTEREDCUBIC
    a1(0.5,0.5,0.0),
    a2(0.5,0.0,0.5),
    a3(0. ,0.5,0.5),
#elif defined HEXAGONALLAYERS
    a1(1.,0.,0.),
    a2(0.5,0.8660254038,0.),
    a3(0.,0.,1.),
#endif

    origo(0.,0.,0.),   
    
    site_rvol(v),
    site_rdims(d),
    site_r(site_rvol),
    
    site_qvol( N3*N2*(N1/2+1)),
    site_qdims(d),
    site_q(site_qvol)
{
  logfile << "initializing BravaisLattice" << endl;
  if(TRACE) cout << "d: " << d << endl;
  if(TRACE) cout << "v: " << v << endl;

  if(d==1){site_rdims[0]=N1;}
  if(d==2){site_rdims[0]=N2; site_rdims[1]=N1;}
  if(d==3){site_rdims[0]=N3; site_rdims[1]=N2; site_rdims[2]=N1;}

  // reciprocal basisvectors, but without the 2pi factor
  double invunitvolume=1./scalarproduct(a1,crossproduct(a2,a3));
  b1= invunitvolume*crossproduct(a2,a3);
  b2= invunitvolume*crossproduct(a3,a1);
  b3= invunitvolume*crossproduct(a1,a2);


  if(TRACE) cout << "Starting to record positions" << endl;
  if(TRACE) cout << "N1=" << N1 << " N2=" << N2 << " N3=" << N3 << endl;
  // set the coordinate positions and phase

  int s=0;
  for(int i3=0; i3<N3; i3++)
    for(int i2=0; i2<N2; i2++)
      for(int i1=0; i1<N1; i1++)
	{

#ifdef RECIPROCAL_LATTICE
	  Coord position=origo+TWOPI*((i1*1./N1)*a1+(i2*1./N2)*a2+(i3*1./N3)*a3);
#else
	  Coord position=origo+i1*a1+i2*a2+i3*a3;
#endif

	  site_r[s++]=position;

	  if(TRACE) cout << "site: " << s 
		      //			 << " phase: " << phase[s] 
			 << " int coord: (" 
			 << i1 << " " << i2 << " " << i3 << ") " << position << " \n"; 
	}

  //  ofstream file_siter( SITERPTS.c_str());
  //  for(int i=0; i<site_rvol; i++) file_siter << site_r[i] << endl; 


  // initializing the site q-coordinates
  if(d==1){site_qdims[0]=N1/2+1;}
  if(d==2){site_qdims[0]=N2; site_qdims[1]=N1/2+1;}
  if(d==3){site_qdims[0]=N3; site_qdims[1]=N2; site_qdims[2]=N1/2+1;}

  // specify the q-coordinates of the entries in the
  // fourier-transform matrix;

  int i=0;
  if(d==1)
    { 
      for(int q1=0; q1< site_qdims[0]; q1++) 
#ifdef RECIPROCAL_LATTICE
	site_q[i++]=q1*b1;
#else
	site_q[i++]=TWOPI*(q1*1./N1)*b1;
#endif
    }
  if(d==2)
    {
      for(int q2=0; q2< site_qdims[0]; q2++)
	for(int q1=0; q1< site_qdims[1]; q1++) 
#ifdef RECIPROCAL_LATTICE
      site_q[i++]= q1*b1+q2*b2;
#else
      site_q[i++]=TWOPI*((q1*1./N1)*b1+(q2*1./N2)*b2);
#endif
    }
  if(d==3)
    {
      for(int q3=0; q3< site_qdims[0]; q3++)
	for(int q2=0; q2< site_qdims[1]; q2++)
	  for(int q1=0; q1< site_qdims[2]; q1++) 
#ifdef RECIPROCAL_LATTICE
      site_q[i++]= q1*b1+q2*b2+q3*b3;
#else
      site_q[i++]=TWOPI*( (q1*1./N1)*b1+(q2*1./N2)*b2+(q3*1./N3)*b3 );
#endif
    }
  
  
  //  ofstream file_siteq( SITEQPTS.c_str());
  //  for(int i=0; i<site_qvol; i++) file_siteq << setprecision(16) << site_q[i] << endl; 

}

#endif
