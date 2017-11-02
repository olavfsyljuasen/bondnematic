#ifndef GLOBALHEADER_H
#define GLOBALHEADER_H

#include <fstream>
#include <iomanip>
#include<stdio.h>
#include<stdlib.h>
#include <string>
#include<vector>
#include<time.h>
#include "math.h"
using namespace std;





string bitstring(int s,int cycle)
{
  //  string str=(s==0 ? "0":"");
  string str="";
  //  for(int i=0; (s!=0 || i<cycle); s>>=1, i++){if(s & 1){str='1'+str;}else{str='0'+str;}}
  for(int i=0; i<cycle; s>>=1, i++){if(s & 1){str='1'+str;}else{str='0'+str;}}
  return str;
}


class Timer{
  friend ostream&
    operator<<(ostream& os, Timer timer)
    { time_t t;
      if( time(&t) == time_t(-1)){exit(1);}
      os << ctime(&t); return os;
    }  
 public:
  Timer(){};
  ~Timer(){if(TRACE) cout << "deleting timer\n";}
  void Start(){t1=time(0);}
  void Stop(){t2=time(0);}
  double GetTimeElapsed(){return difftime(t2,t1);}
 private:
  time_t t1;
  time_t t2;
};


class Coord
{  friend ostream& operator<<(ostream& os,Coord& c)
    { os << c.x << " " << c.y << " " << c.z; return os; }
 public:
  Coord(double xx=0.,double yy=0.,double zz=0.):x(xx),y(yy),z(zz){}
  double x;
  double y;
  double z;
  void clear(){x=0.;y=0.;z=0.;}
  double Norm(){return sqrt(x*x+y*y+z*z);}
};

inline bool operator==(const Coord& a,const Coord& b)
{ return bool( (a.x == b.x) && (a.y==b.y) && (a.z==b.z));}

Coord operator-(const Coord& l,const Coord& r){
  return Coord(l.x-r.x,l.y-r.y,l.z-r.z);
}

Coord operator+(const Coord& l,const Coord& r){
  return Coord(l.x+r.x,l.y+r.y,l.z+r.z);
}

double operator*(const Coord& l,const Coord& r){
  return double(l.x*r.x+l.y*r.y+l.z*r.z);
}

Coord operator*(const int& k,const Coord& a){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord operator*(const Coord& a,const int& k){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord operator*(const double& k,const Coord& a){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord operator*(const Coord& a,const double& k){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord& operator+=(Coord& a,const Coord& b){
  a.x += b.x;   a.y += b.y;   a.z += b.z; return a;
}

Coord& operator-=(Coord& a,const Coord& b){
  a.x -= b.x;   a.y -= b.y;   a.z -= b.z; return a;
}

double scalarproduct(const Coord& a,const Coord& b){
return a.x*b.x+a.y*b.y+a.z*b.z;
}

Coord crossproduct(const Coord& a,const Coord& b){
  Coord result(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); return result;
}




// a routine that converts integers to strings
const int MAXDIGITS = 20;
const int ASCII0 = 48;
string int2string(int a)
{
  int digit[MAXDIGITS];
  int d=0;
  string s("");
  if (a == 0){ s+=ASCII0; return s;}
  while( a > 0){ digit[d++] = a-10*(a/10); a /=10 ;}
  for(int i=d-1; i>=0; i--) s+= digit[i]+ASCII0;
  return s;
}


// this routine also handles leading zeros
string int2string(int a,int nsiffer)
{
  int digit[MAXDIGITS];
  string s("");
  for(int d=0; d<MAXDIGITS; d++){ digit[d] = a-10*(a/10); a /=10 ;}
  for(int i=nsiffer-1; i>=0; i--) s+= digit[i]+ASCII0;
  return s;
}




#endif
