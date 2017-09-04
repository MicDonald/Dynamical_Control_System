#include "KernelMatrix.h"

using namespace Eigen;
using namespace std;

KernelMatrix::KernelMatrix(const VRAtomInfo& model):model(model){
    gdof = model.atomGID.size()*3;
    vdof = model.atomVirtual.size()*3;
    rdof = model.atomReal.size()*3;
}

void KernelMatrix::setK_mass(double k_m){
  k_mass=k_m;
}

double KernelMatrix::getK_mass(){
        return k_mass;
}
double KernelMatrix::length3D(vector<double> p1, vector<double> p2){
    if(p1.size() != 3 || p2.size() != 3){
        cout<<"Error: position must have 3 dimension"<<endl;
        return 0.;
    }
    else
        return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
}
Eigen::MatrixXd KernelMatrix::localKe(vector<double> p1, vector<double> p2){
    Matrix2d k;
    k<<1,-1,-1,1;
    double bondLength = length3D(p1, p2);
    double cx=(p2[0]-p1[0])/bondLength;
    double cy=(p2[1]-p1[1])/bondLength;
    double cz=(p2[2]-p1[2])/bondLength;
    MatrixXd J1(2,6);
    MatrixXd J2(2,6);
    MatrixXd J3(2,6);
    J1<<cx,cy,cz,0.,0.,0.,0.,0.,0.,cx,cy,cz;
    // J2<<cy,cz,cx,0.,0.,0.,0.,0.,0.,cy,cz,cx;
    // J3<<cz,cx,cy,0.,0.,0.,0.,0.,0.,cz,cx,cy;
    // MatrixXd ke = J1.transpose() * k * J1 +
    //               J2.transpose() * k * J2 +
    //               J3.transpose() * k * J3;
    MatrixXd ke=J1.transpose() * k * J1;
    return ke*k_mass;
}

void KernelMatrix::bondNeighborIdendifier(){
    double L = 10000.;
    for (int i = 0 ;i<model.atomGID.size();++i){
        for (int j = 0 ;j<model.atomGID.size();++j){
            if (i!=j){
              double temp = length3D(model.atomCoord[i],model.atomCoord[j]);
              if (temp<L) L = temp;
            }
        }
    }
    cout<<"bond Length: "<<L<<endl;
    for (int i = 0 ;i<model.atomGID.size();++i){
      for (int j = 0 ;j<model.atomGID.size();++j){        
        if(i!=j){
          double temp = length3D(model.atomCoord[i],model.atomCoord[j]);
          if (abs(temp-L)<1e-1){
            vector<int> tempV={j,i};
            if (count(bonds.begin(),bonds.end(),tempV) == 0)
              bonds.insert( { i, j } );
            if (count(model.atomVirtual.begin(),model.atomVirtual.end(),i) == 1
            && count(model.atomReal.begin(),model.atomReal.end(),j) == 1)
              VRbonds.insert( { i, j } );
          }
        }
      }
    }
    set<int> tempV,tempR;
    cout<<"Identify "<<bonds.size()<<" bonds:"<<endl;
    cout<<VRbonds.size()<<" VRbonds:"<<endl;
    for (const auto& e : VRbonds){
        tempV.insert(e[0]);
        tempR.insert(e[1]);
        cout<<"V :"<<e[0]<<" &R: "<<e[1]<<" ";
    }
    cout<<endl;
    model.atomV2r.assign(tempV.begin(),tempV.end());
    model.atomR2v.assign(tempR.begin(),tempR.end());
    for (const auto& e : model.atomV2r){
        cout<<"["<<e<<"]"<<model.atomGID[e]<<" ";
    }
    cout<<endl;
    for (const auto& e : model.atomR2v){
        cout<<"["<<e<<"]"<<model.atomGID[e]<<" ";
    }
    cout<<endl;
}

void KernelMatrix::angleNeighborIdendifier(){


}

vector<int> KernelMatrix::bondOrAtom2MatrixDof(vector<int> atomList){
    vector<int> dof;
    for (int i=0;i<atomList.size();++i)
        for (int ii=0;ii<3;++ii)
            dof.push_back(atomList[i]*3+ii);
    return dof;
}
vector<int> KernelMatrix::bondOrAtom2MatrixDof(int atomID){
    vector<int> dof;
    for (int ii=0;ii<3;++ii)
    	dof.push_back(atomID*3+ii);
    return dof;
}


void KernelMatrix::calculateEigen(){
  gdof = model.atomGID.size()*3;
  vdof = model.atomVirtual.size()*3;
  rdof = model.atomReal.size()*3;
  v2rdof = model.atomV2r.size()*3;
  r2vdof = model.atomR2v.size()*3;
  SparseMatrix<double> D;
  d.resize(vdof, vdof);
  X.resize(vdof, vdof);
  D.resize(gdof,gdof);
  DVRreduced.resize(vdof, r2vdof); 
  MatrixXd DV(vdof,vdof);

#pragma omp single
{
  for(const auto& e : bonds){
    auto dof = bondOrAtom2MatrixDof(e);
    int i=0,j=0;
    auto ke=localKe(model.atomCoord[e[0]],model.atomCoord[e[1]]);
    for (const auto& eRow : dof){
      for (const auto& eCol :dof){
        D.coeffRef(eRow,eCol)+=ke(i,j);
        ++j;
      }
      j=0;
      ++i;
    }
  }


  
  cout<<"D: "<<gdof<<"x"<<D.outerSize()<<endl;
  Vdof=bondOrAtom2MatrixDof(model.atomVirtual);
  Rdof=bondOrAtom2MatrixDof(model.atomReal);
  V2Rdof=bondOrAtom2MatrixDof(model.atomV2r);
  R2Vdof=bondOrAtom2MatrixDof(model.atomR2v);
  int i=0,j=0,k=0;
  for (const auto& eRow : Vdof){
    for (const auto& eCol : Vdof){
      DV(i,j)=D.coeff(eRow,eCol);
      ++j;
    }
    for (const auto& eCol : R2Vdof){
      DVRreduced.insert(i,k)=D.coeff(eRow,eCol);
      ++k;
    }
	++i;
	j=k=0;
  }
  cout<<"DV: "<<vdof<<"x"<<DV.size()/vdof<<endl;
  cout<<"DVRreduced:"<<vdof<<"x"<<r2vdof<<endl;
}
  bool succ1 = loadMatrix("evalue",this->d,vdof);
  bool succ2 = loadMatrix("evector",this->X,vdof);
  if (!(succ1 && succ2)){
  #pragma omp parallel
  {
  SelfAdjointEigenSolver<MatrixXd> solver(DV);
  #pragma omp single
  {
    d = solver.eigenvalues().asDiagonal();
    X = solver.eigenvectors(); 
    saveMatrix("evalue",d,true);
    saveMatrix("evector",X,true);
    cout<<"eigen calculation done"<<endl;
  }
  }
  }
  right=X.transpose()*DVRreduced;
}



//SparseMatrix<double> KernelMatrix::calculateKernelMatrix(double t){ 
MatrixXd KernelMatrix::calculateKernelMatrix(double t){ 
//#pragma omp critical
//{
  SparseMatrix<double> sinHF;
  sinHF.resize(vdof,vdof);

  //SparseMatrix<double> reducedKM;
  MatrixXd reducedKM;
  reducedKM.resize(v2rdof, r2vdof);
  clock_t Ttemp=clock();
  for (int i = 0;i<vdof; ++i){
    if (d(i,i)==d(i,i) && abs(d(i,i))>1e-9) {
      complex<double> wi = sqrt( complex<double>(uv0(i,0)*uv0(i,0)/4 - d(i,i) +vv0(i,0)));
      complex<double> sinH = sinh(t * wi);
      sinHF.insert(i,i)=exp(uv0(i,0)*t/2)*real(sinH/wi);
            //no uv0
      // complex<double> wi = sqrt( complex<double>(- d(i,i) ));
      // complex<double> sinH = sinh(t * wi);
      // sinHF.insert(i,i)=real(sinH/wi);
    }
  }
  t_w+=clock()-Ttemp;

  Ttemp=clock();

  MatrixXd left=-X*sinHF;
  
  // Matrix tempKM=left*right;
  t_calK+=clock()-Ttemp;
  // Ttemp=clock();
  // int j=0;

  #pragma omp parallel for
  for (int i = 0;i<v2rdof; ++i){
    // if(count(V2Rdof.begin(),V2Rdof.end(),Vdof[i])==1){
    auto iter = find( Vdof.begin(), Vdof.end(), V2Rdof[i]);
    int ii = distance( Vdof.begin(), iter);
    for (int j = 0;j<r2vdof; ++j){
      // cout<<"Vdof: "<<V2Rdof[i]<<" Rdof: "<<R2Vdof[j]<<"i: "<<i<<"j: "<<j<<"ii: "<<ii<<endl;
      // reducedKM.insert(i,j)=tempKM(ii,j);
      // Test(i,j);
      //reducedKM.insert(i,j)=left.row(ii)*right.col(j);
      reducedKM(i,j)=left.row(ii)*right.col(j);
    }
  }
  return reducedKM;
}



void KernelMatrix::set0state(MatrixXd uv0, MatrixXd vv0){
  if(this->uv0.size()==uv0.size()){
    cout<<"zero state error"<<endl;
    exit(1);
  }
  this->uv0=uv0;
  this->vv0=vv0;

}

 
 //template<typename int RowsAtCompileTime, int ColsAtCompileTime>
inline bool loadMatrix(string filename, MatrixXd& m,int dof)
{
   // General structure
   // 1. Read file contents into vector<double> and count number of lines
   // 2. Initialize matrix
   // 3. Put data in vector<double> into matrix
   
  ifstream input(filename.c_str());
  if (input.fail())
  {
    cout << "No '" << filename << "'." << endl;
    return false;
  }
  string line;
  double d;
   
  vector<double> v;
  while (getline(input, line))
  {
    stringstream input_line(line);
    while (!input_line.eof())
    {
      input_line >> d;
      v.push_back(d);
    }
  }
  input.close();
  if (v.size()!=dof*dof) {
    cerr << "ERROR. Please delete '"<< filename << "'." << endl;
    exit(1);
  }
  m.setZero(dof,dof);

  for (int i=0; i<dof; i++)
    for (int j=0; j<dof; j++)
      m(i,j) = v[i*dof + j];

  cout<<"MatrixXd of "<<filename<<" read done!"<<endl; 
  // cout<<"Matrix:"<<m; 
  return true;
}
 
 //template<typename int RowsAtCompileTime, int ColsAtCompileTime>
inline bool saveMatrix(string filename, MatrixXd matrix, bool overwrite)
{
  //if (boost::filesystem::exists(filename)){
    if (!overwrite)
    {
      // File exists, but overwriting is not allowed. Abort.
      cerr << "File '" << filename << "' already exists. Not saving data." << endl;
      return false;
    }
   //} 
  ofstream file;
  file.open(filename.c_str());
  if (!file.is_open())
  {
    cerr << "Couldn't open file '" << filename << "' for writing." << endl;
    return false;
  }
  // MatrixXd tempMatrix = matrix;
  file << fixed;
  file << matrix;
  file.close();
   
  return true;
   
}

inline void Test( int n, int m )
{
  printf( "<T:%d> - %d, %d\n", omp_get_thread_num(), n, m );
}
 




