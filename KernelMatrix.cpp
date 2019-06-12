#include "KernelMatrix.h"

using namespace Eigen;
using namespace std;

KernelMatrix::KernelMatrix(const VRAtomInfo& model): model(model) {
  gdof = model.atomGID.size() * 3;
  vdof = model.atomVirtual.size() * 3;
  rdof = model.atomReal.size() * 3;
}

void KernelMatrix::setK_mass(double k_m) {
  k_mass = k_m;
}
void KernelMatrix::setAngle(double ang) {
  angle = ang;
}

double KernelMatrix::getK_mass() {
  return k_mass;
}


vector<double> KernelMatrix::vectorBy2Points(vector<double> p1, vector<double> p2) {
  vector<double> v12 = {p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]};
  return v12;
}


double KernelMatrix::lengthOfVector(vector<double> v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

double KernelMatrix::dotBy2Vectors(vector<double> v1, vector<double> v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double KernelMatrix::cosBy2Vectors(vector<double> v1, vector<double> v2) {
  return dotBy2Vectors(v1, v2) / lengthOfVector(v1) / lengthOfVector(v2);
}

vector<int> KernelMatrix::bondOrAtom2MatrixDof(const vector<int> atomList) {
  vector<int> dof;
  for (int i = 0; i < atomList.size(); ++i)
    for (int ii = 0; ii < 3; ++ii)
      dof.push_back(atomList[i] * 3 + ii);
  return dof;
}

vector<int> KernelMatrix::bondOrAtom2MatrixDof(const pair<int, int> atomList) {
  vector<int> dof;
  for (int ii = 0; ii < 3; ++ii)
    dof.push_back(atomList.first * 3 + ii);
  for (int ii = 0; ii < 3; ++ii)
    dof.push_back(atomList.second * 3 + ii);
  return dof;
}

vector<int> KernelMatrix::bondOrAtom2MatrixDof(const int atomID) {
  vector<int> dof;
  for (int ii = 0; ii < 3; ++ii)
    dof.push_back(atomID * 3 + ii);
  return dof;
}

void KernelMatrix::bondIdendifier() {
  double L = 10000.;
  for (int i = 0 ; i < model.atomGID.size(); ++i) {
    for (int j = 0 ; j < model.atomGID.size(); ++j) {
      if (i != j) {
        double temp = lengthOfVector(vectorBy2Points(model.atomCoord[i], model.atomCoord[j]));
        if (temp < L) L = temp;
      }
    }
  }
  cout << "bond Length: " << L << endl;
  for (int i = 0 ; i < model.atomGID.size(); ++i) {
    for (int j = 0 ; j < model.atomGID.size(); ++j) {
      if (i != j) {
        double temp = lengthOfVector(vectorBy2Points(model.atomCoord[i], model.atomCoord[j]));
        if (abs(temp - L) < 1e-1) {
          bonds.insert( { i, j } );
          std::vector<int> tempV = {j, i};
          if (count(pairs.begin(), pairs.end(), tempV) == 0)
            pairs.insert( { i, j } );
          if (count(model.atomVirtual.begin(), model.atomVirtual.end(), i) == 1
              && count(model.atomReal.begin(), model.atomReal.end(), j) == 1)
            VRbonds.insert( { i, j } );
        }
      }
    }
  }

  set<int> tempV, tempR;
  cout << "Identify " << bonds.size() << " bonds." << endl;
  cout << "Identify " << pairs.size() << " pairs." << endl;
  cout << VRbonds.size() << " VRbonds:" << endl;
  for (const auto& e : VRbonds) {
    // tempV.insert(e[0]);
    // tempR.insert(e[1]);
    tempV.insert(e.first);
    tempR.insert(e.second);
    cout << "V :" << e.first << " &R: " << e.second << " ";
  }
  cout << endl;
  model.atomV2r.assign(tempV.begin(), tempV.end());
  model.atomR2v.assign(tempR.begin(), tempR.end());
  for (const auto& e : model.atomV2r) {
    cout << "[" << e << "]" << model.atomGID[e] << " ";
  }
  cout << endl;
  for (const auto& e : model.atomR2v) {
    cout << "[" << e << "]" << model.atomGID[e] << " ";
  }
  cout << endl;
  gdof = model.atomGID.size() * 3;
  vdof = model.atomVirtual.size() * 3;
  rdof = model.atomReal.size() * 3;
  v2rdof = model.atomV2r.size() * 3;
  r2vdof = model.atomR2v.size() * 3;
}

void KernelMatrix::angleIdendifier() {
  if (angle != 0) {
    for (int o = 0 ; o < model.atomGID.size(); ++o) {
      pair<multimap<int, int>::iterator, multimap<int, int>::iterator> allBonded = bonds.equal_range(o);
      for (auto a = allBonded.first; a != allBonded.second; ++a) {
        for (auto b = a; b != allBonded.second; ++b) {
          if (a == b) continue;
          int A = a->second, B = b->second;
          double OA = lengthOfVector(vectorBy2Points(model.atomCoord[o], model.atomCoord[A]));
          double OB = lengthOfVector(vectorBy2Points(model.atomCoord[o], model.atomCoord[B]));
          double AB = lengthOfVector(vectorBy2Points(model.atomCoord[A], model.atomCoord[B]));
          //cout<<" { "<< AB<<", "<< OA<<", "<< OB<<" } "<< acos( (OA*OA + OB*OB - AB*AB)/2/OA/OB )*180/3.14159265;
          if ( abs( acos( (OA * OA + OB * OB - AB * AB) / 2 / OA / OB ) * 180 / 3.14159265 - 60) < 1e-4 )
            angles.push_back(make_tuple(A, o, B));
        }
      }
    }
  }
}

Eigen::MatrixXd KernelMatrix::localKe(vector<double> p1, vector<double> p2) {
  Matrix2d k;
  k << 1, -1, -1, 1;
  double bondLength = lengthOfVector(vectorBy2Points(p1, p2));
  double cx = (p2[0] - p1[0]) / bondLength;
  double cy = (p2[1] - p1[1]) / bondLength;
  double cz = (p2[2] - p1[2]) / bondLength;
  MatrixXd J1(2, 6);
  MatrixXd J2(2, 6);
  MatrixXd J3(2, 6);
  J1 << cx, cy, cz, 0., 0., 0., 0., 0., 0., cx, cy, cz;
  // J2<<cy,cz,cx,0.,0.,0.,0.,0.,0.,cy,cz,cx;
  // J3<<cz,cx,cy,0.,0.,0.,0.,0.,0.,cz,cx,cy;
  // MatrixXd ke = J1.transpose() * k * J1 +
  // J2.transpose() * k * J2 +
  // J3.transpose() * k * J3;
  MatrixXd ke = J1.transpose() * k * J1;
  return ke * k_mass;
}

void KernelMatrix::calculateEigen() {
  // gdof = model.atomGID.size() * 3;
  // vdof = model.atomVirtual.size() * 3;
  // rdof = model.atomReal.size() * 3;
  // v2rdof = model.atomV2r.size() * 3;
  // r2vdof = model.atomR2v.size() * 3;
  SparseMatrix<double> D;
  d.resize(vdof, vdof);
  X.resize(vdof, vdof);
  D.resize(gdof, gdof);
  DVRreduced.resize(vdof, r2vdof);
  MatrixXd DV(vdof, vdof);

  #pragma omp single
  {
    for (const auto& e : pairs) {
      auto dof = bondOrAtom2MatrixDof(e);
      int i = 0, j = 0;
      auto ke = localKe(model.atomCoord[e[0]], model.atomCoord[e[1]]);
      //auto ke=localKe(model.atomCoord[e.first],model.atomCoord[e.second]);
      for (const auto& eRow : dof) {
        for (const auto& eCol : dof) {
          D.coeffRef(eRow, eCol) += ke(i, j);
          ++j;
        }
        j = 0;
        ++i;
      }
    }


    cout << "D: " << gdof << "x" << D.outerSize() << endl;
    Vdof = bondOrAtom2MatrixDof(model.atomVirtual);
    Rdof = bondOrAtom2MatrixDof(model.atomReal);
    V2Rdof = bondOrAtom2MatrixDof(model.atomV2r);
    R2Vdof = bondOrAtom2MatrixDof(model.atomR2v);
    int i = 0, j = 0, k = 0;
    for (const auto& eRow : Vdof) {
      for (const auto& eCol : Vdof) {
        DV(i, j) = D.coeff(eRow, eCol);
        ++j;
      }
      for (const auto& eCol : R2Vdof) {
        DVRreduced.insert(i, k) = D.coeff(eRow, eCol);
        ++k;
      }
      ++i;
      j = k = 0;
    }
    cout << "DV: " << vdof << "x" << DV.size() / vdof << endl;
    cout << "DVRreduced:" << vdof << "x" << r2vdof << endl;
  }
  if (vdof > 1000) {
    bool succ1 = loadMatrix("evalue", d, vdof);
    bool succ2 = loadMatrix("evector", X, vdof);
    if (!(succ1 && succ2)) {
      #pragma omp parallel
      {
        SelfAdjointEigenSolver<MatrixXd> solver(DV);
        #pragma omp single
        {
          d = solver.eigenvalues().asDiagonal();
          X = solver.eigenvectors();
          saveMatrix("evalue", d, true);
          saveMatrix("evector", X, true);
          cout << "eigen calculation done" << endl;
        }
      }
    }
  }
  else {
    #pragma omp parallel
    {
      SelfAdjointEigenSolver<MatrixXd> solver(DV);
      #pragma omp single
      {
        cout << "Vdof is small. Ignore save & load." << endl;
        d = solver.eigenvalues().asDiagonal();
        X = solver.eigenvectors();
        cout << "eigen calculation done" << endl;
      }
    }

  }
  right = X.transpose() * DVRreduced;
}



//SparseMatrix<double> KernelMatrix::calculateKernelMatrix(double t){
MatrixXd KernelMatrix::calculateKernelMatrix(double t, bool deri = false) {
// #pragma omp critical
// {
  SparseMatrix<double> sinHF; //, cosDHF;
  sinHF.resize(vdof, vdof);
  // cosDHF.resize(vdof, vdof);
  //SparseMatrix<double> reducedKM;
  MatrixXd reducedKM;
  reducedKM.resize(v2rdof, r2vdof);
  clock_t Ttemp = clock();
  for (int i = 0; i < vdof; ++i) {
    if (d(i, i) == d(i, i) && abs(d(i, i)) > 1e-9) {
      complex<double> wi = sqrt( complex<double>(uv0(i, 0) * uv0(i, 0) / 4 - d(i, i) + vv0(i, 0)));
      complex<double> sinH = sinh(t * wi);
      // complex<double> cosH = cosh(t * wi);
      sinHF.insert(i, i) = exp(uv0(i, 0) * t / 2) * real(sinH / wi);
      // cosDHF.insert(i, i) = exp(uv0(i, 0) * t / 2) * uv0(i, 0) / 2 * real(sinH / wi) 
      //                     + exp(uv0(i, 0) * t / 2) * real(cosH);
    }
  }
  t_w += clock() - Ttemp;

  Ttemp = clock();
  
  MatrixXd left;
  // if (!deri)
  left = -X * sinHF;
  // else
    // left = -X * cosDHF;

  // Matrix tempKM=left*right;
  t_calK += clock() - Ttemp;
  // Ttemp=clock();
  // int j=0;

  #pragma omp parallel for
  for (int i = 0; i < v2rdof; ++i) {
    // if(count(V2Rdof.begin(),V2Rdof.end(),Vdof[i])==1){
    auto iter = find( Vdof.begin(), Vdof.end(), V2Rdof[i]);
    int ii = distance( Vdof.begin(), iter);
    for (int j = 0; j < r2vdof; ++j) {
      // cout<<"Vdof: "<<V2Rdof[i]<<" Rdof: "<<R2Vdof[j]<<"i: "<<i<<"j: "<<j<<"ii: "<<ii<<endl;
      // reducedKM.insert(i,j)=tempKM(ii,j);
      // Test(i,j);
      //reducedKM.insert(i,j)=left.row(ii)*right.col(j);
      reducedKM(i, j) = left.row(ii) * right.col(j);
    }
  }
  return reducedKM;
}



void KernelMatrix::set0state(MatrixXd uv0, MatrixXd vv0) {
  if (this->uv0.size() == uv0.size()) {
    cout << "zero state error" << endl;
    exit(1);
  }
  this->uv0 = uv0;
  this->vv0 = vv0;

}
void KernelMatrix::printAngles() {
  cout << "Angles: (A, o, B)" << endl;
  for (auto i = angles.begin(); i != angles.end(); ++i)
  {
    cout << "{ " << get<0>(*i) << ", " << get<1>(*i) << ", " << get<2>(*i) << " }  ";
  }
  cout << endl;
}


//template<typename int RowsAtCompileTime, int ColsAtCompileTime>
inline bool loadMatrix(string filename, MatrixXd& m, int dof)
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
  if (v.size() != dof * dof) {
    cerr << "ERROR. Please delete '" << filename << "'." << endl;
    exit(1);
  }
  m.setZero(dof, dof);

  for (int i = 0; i < dof; i++)
    for (int j = 0; j < dof; j++)
      m(i, j) = v[i * dof + j];

  cout << "MatrixXd of " << filename << " read done!" << endl;
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





