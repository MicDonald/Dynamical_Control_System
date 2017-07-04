#include "KernelMatrix.h"
using namespace Eigen;
using namespace std;

KernelMatrix::KernelMatrix(const VRAtomInfo& model):model(model){
    gdof = model.atomGID.size()*3;
    vdof = model.atomVirtual.size()*3;
    rdof = model.atomReal.size()*3;
    DVR.resize(vdof, rdof);
    d.resize(vdof, vdof);
    X.resize(vdof, vdof);
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
    J2<<cy,cz,cx,0.,0.,0.,0.,0.,0.,cy,cz,cx;
    J3<<cz,cx,cy,0.,0.,0.,0.,0.,0.,cz,cx,cy;
    MatrixXd ke = J1.transpose() * k * J1 +
                  J2.transpose() * k * J2 +
                  J3.transpose() * k * J3;
    //MatrixXd ke=J1.transpose() * k * J1;
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
            if (abs(temp-L)<10e-3){
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

    MatrixXd D;
    d.setZero(vdof, vdof);
    X.setZero(vdof, vdof);
    D.setZero(gdof,gdof);
    DVR.setZero(vdof, rdof); 
    MatrixXd DV(vdof,vdof);
#pragma omp single
{
    for(const auto& e : bonds){
        auto dof = bondOrAtom2MatrixDof(e);
	int i=0,j=0;
	auto ke=localKe(model.atomCoord[e[0]],model.atomCoord[e[1]]);
        for (const auto& eRow : dof){
            for (const auto& eCol :dof){
                D(eRow,eCol)+=ke(i,j);
                ++j;
            }
            j=0;
            ++i;
        }
    }
    cout<<"D: "<<gdof<<"x"<<D.size()/gdof<<endl;
    Vdof=bondOrAtom2MatrixDof(model.atomVirtual);
    Rdof=bondOrAtom2MatrixDof(model.atomReal);
    V2Rdof=bondOrAtom2MatrixDof(model.atomV2r);
    R2Vdof=bondOrAtom2MatrixDof(model.atomR2v);
    int i=0,j=0,k=0;
    for (const auto& eRow : Vdof){
        for (const auto& eCol : Vdof){
            DV(i,j)=D(eRow,eCol);
            ++j;
        }
        for (const auto& eCol : Rdof){
            DVR(i,k)=D(eRow,eCol);
            ++k;
        }
	++i;
	j=k=0;
    }

    cout<<"DV: "<<vdof<<"x"<<DV.size()/vdof<<endl;
    cout<<"DVR:"<<vdof<<"x"<<rdof<<endl;
}
    loadMatrix("evalue",this->d,vdof);
    loadMatrix("evector",this->X,vdof);
if (d(0,0)<1e-10){
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
// cout<<"eigen value after loading"<<endl;
// cout<<d<<endl;
}

MatrixXd KernelMatrix::calculateKernelMatrix(double t){ 
//#pragma omp critical
//{
    MatrixXd sinF,reducedKM;
    sinF.setZero(vdof,vdof);
    reducedKM.setZero(v2rdof, r2vdof);
    for (int i = 0;i<vdof; ++i)
	for (int j =0;j<vdof; ++j){
        	if (i==j) sinF(i,j)=sin(sqrt(d(i,j))*t)/sqrt(d(i,j));
	}
    MatrixXd tempKM=(X*sinF*(X.transpose()))*(-DVR);
    //cout<<"tempKM\n"<<tempKM<<endl;
    //cout<<"Reduced KM size:"<<v2rdof<<"x"<<r2vdof<<endl;
    int k=0,l=0;
    for (int i = 0;i<vdof; ++i){
	if(count(V2Rdof.begin(),V2Rdof.end(),Vdof[i])==1){
            for (int j =0;j<rdof; ++j){
                if(count(R2Vdof.begin(),R2Vdof.end(),Rdof[j])==1){
	            //cout<<"Vdof: "<<Vdof[i]<<" Rdof: "<<Rdof[j]<<" "<<k<<" "<<l<<" "<<i<<" "<<j<<" | ";
		    reducedKM(k,l)=tempKM(i,j);
	    	    ++l;
	        }
	    }
	    ++k;
	    l=0;
        }
    }
    //cout<<"reducedKM:\n"<<reducedKM<<endl;
    return reducedKM;
//}
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
   
   file << fixed;
   file << matrix;
   file.close();
   
   return true;
   
 }
 




