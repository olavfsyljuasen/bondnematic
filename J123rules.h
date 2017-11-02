#ifndef J123RULES_Hq
#define J123RULES_H


using namespace std;

// The input parameter sequence
const int NPARAMS = 14;
enum params{LINEID,J1,J2,J3,JC,DELTA,NX,NY,NZ,DA,MAXITER,TOLERANCE,NBINS,EQFLAG};


const int NS=3; // the number of spin components


//OBSERVABLES:

const int NOBSERVABLES=2;
enum observables{SIGMAA,SIGMAD};


class Rule
{
  friend ostream& operator<<(ostream& os,Rule& c){ return c.Write(os);}
  friend istream& operator>>(istream& is,Rule& c){ c.Read(is); return is;}
 public:
  template<class LATTICE> 
    Rule(double*,LATTICE& la);

  
  vector<double>& GetInitialState(){return Kinit;}
  vector<double>& GetIrrep(const int i){return irrep[i];}
  vector<double>& GetInteraction(){return  Jq;}

  ostream& Write(ostream& os);
  istream& Read(istream& is);
  double* GetPars(){return param;}
 private:
  double* param;

  int Vq;
  
  vector<double> Kinit;
  vector<double> Jq;
  vector<vector<double> > irrep;

  void SetInteraction();  
  void SetIrreps();
  void SetInitialState(double);

};

template<class LATTICE>
Rule::Rule(double* par,LATTICE& lattice):param(par),Vq(lattice.SiterVol()),Kinit(Vq),Jq(Vq),irrep(NOBSERVABLES,vector<double>(Vq))
{
  if(TRACE) cout << "Initializing Rule " << endl; 
  for(int i=0; i<Vq; i++)
    {
      double qx=lattice.Pos(i).x; // should be momentum-space coordinates
      double qy=lattice.Pos(i).y;
      double qz=lattice.Pos(i).z;
      
      double da=param[DA];
      Kinit[i]=da*(cos(qx)-cos(qy)+cos(qx+qy)-cos(qx-qy)
		   +cos(4*qx)*sin(qy)*sin(qy)-cos(4*qy)*sin(qx)*sin(qx)
		   +cos(2*qx-4*qy)-sin(qx+3*qy)*sin(qx+3*qy)
		   +cos(qz) + sin(qz)*sin(qz));

      Jq[i]=
	//	-param[J1]*(cos(qx)+cos(qy))
	+param[J1]*(cos(qx)+cos(qy))
	+2.*param[J2]*cos(qx)*cos(qy)
	+param[J3]*(cos(2.*qx)+cos(2.*qy))
        +param[JC]*cos(qz);

      irrep[SIGMAA][i]=cos(qx)-cos(qy);       // the sigma_a order parameter
      irrep[SIGMAD][i]=cos(qx+qy)-cos(qx-qy); // the sigma_d order parameter      
    }

  if(TRACE) cout << "Done Initializing Rule " << endl; 
}



ostream& Rule::Write(ostream& os)
{
  //  os.write( (char*) &Nbonds,sizeof(Nbonds));
  //  os.write( (char*) &lambda_total,sizeof(lambda_total));
  return os;
}


istream& Rule::Read(istream& is)
{
  //  is.read( (char*) &Nbonds,sizeof(Nbonds));
  //  is.read( (char*) &lambda_total,sizeof(lambda_total));
  return is;
}

#endif
