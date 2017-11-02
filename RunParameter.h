class RunParameters{
  friend ostream& operator<<(ostream& os,RunParameters& r)
  {
    for(int i=1; i<=r.GetNR(); i++)
      {
	for(int j=0; j<r.GetN(); j++) os << r.GetPar(i,j) << " ";
	os << endl;
      }
    return os;
  }
public:
  RunParameters(string name=READIN);
  ~RunParameters(){if(TRACE) cout << "removing runparameters\n";}
  int GetN(){return N;}
  int GetNR(){return NR;}
  double GetPar(const int r,const int i){return param[r-1][i];}
  void   Decreasenr(const int r){param[r-1][NBINS]--;}
  void   SetNoEqFlag(const int r){param[r-1][EQFLAG]=1;}
  void   UnSetNoEqFlag(const int r){param[r-1][EQFLAG]=0;}
  void   Write(string=READIN){
    ofstream outfile(file.c_str()); outfile << *this;}
  double* GetPars(const int r){return &param[r-1][0];}
private:
  int NR;
  int N;
  string file;
  double param[MAXRUNS][MAXPAR];
};  

  
RunParameters::RunParameters(string name): 
  NR(0),N(0),file(name){
  
  // find number of lines in file
  ifstream infile(file.c_str());
  string d;
  NR = 0;
  while( getline(infile,d,'\n')){ if( d != "") NR++;}
  infile.close();


  // now read the file
  ifstream ifile(file.c_str());
  for(int r = 0; r < NR; r++)
    for(int j=0; j< NPARAMS; j++)
      { ifile >> param[r][j]; 
      //       cout << " " << param[r][j] << " ";
      }
  N = NPARAMS;   
  ifile.close();
}
