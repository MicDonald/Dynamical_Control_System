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
    kernelMatrix.resize(vdof, rdof);
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
    gdof = model.atomGID.size()*3;
    vdof = model.atomVirtual.size()*3;
    rdof = model.atomReal.size()*3;
    DVR.resize(vdof, rdof);
    d.resize(vdof, vdof);
    X.resize(vdof, vdof);
    kernelMatrix.resize(vdof, rdof);    

    double L = 10000.;
    for (int i = 0 ;i!=model.atomGID.size();++i){
        for (int j = 0 ;j!=model.atomGID.size();++j){
            if (i!=j){
	    double temp = length3D(model.atomCoord[i],model.atomCoord[j]);
            if (temp<L) L = temp;
            }
	}
    }
    cout<<"bond Length: "<<L<<endl;
    for (int i = 0 ;i!=model.atomGID.size();++i){
        for (int j = 0 ;j!=model.atomGID.size();++j){        
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
    for (int i=0;i!=atomList.size();++i)
        for (int ii=0;ii!=3;++ii)
            dof.push_back(atomList[i]*3+ii);
    return dof;
}
vector<int> KernelMatrix::bondOrAtom2MatrixDof(int atomID){
    vector<int> dof;
    for (int ii=0;ii!=3;++ii)
    	dof.push_back(atomID*3+ii);
    return dof;
}


void KernelMatrix::calculateEigen(){
    D.setZero(gdof,gdof);
    for(const auto& e : bonds){
        auto dof = bondOrAtom2MatrixDof(e);
	int i=0,j=0;
	auto ke=localKe(model.atomCoord[e[0]],model.atomCoord[e[1]]);
        for (const auto& eRow : dof){
            for (const auto& eCol :dof){
                D(eRow,eCol)+=ke(i,j);
		//cout<<eRow<<" "<<eCol<<" "<<i<<" "<<j<<endl;
                ++j;
            }
            j=0;
            ++i;
        }
    }
    cout<<"D: "<<gdof<<"x"<<D.size()/gdof<<endl;
    MatrixXd DV(vdof,vdof);
    auto Vdof=bondOrAtom2MatrixDof(model.atomVirtual);
    //cout<<"Vdof: ";
    //for (auto e: Vdof) cout<<e<<" ";
    //cout<<"\nRdof: ";
    auto Rdof=bondOrAtom2MatrixDof(model.atomReal);
    //for (auto e: Rdof) cout<<e<<" ";
    //cout<<endl;
    int i=0,j=0,k=0;
    for (const auto& eRow : Vdof){
        for (const auto& eCol : Vdof){
            DV(i,j)=D(eRow,eCol);
            //cout<<i<<" "<<j<<" "<<eRow<<" "<<eCol<<" | ";
            ++j;
        }
        for (const auto& eCol : Rdof){
            DVR(i,k)=D(eRow,eCol);
	    //cout<<i<<" "<<k<<" "<<eRow<<" "<<eCol<<" | ";
            ++k;
        }
        ++i;
        j=k=0;
    }
    cout<<"DV: "<<vdof<<"x"<<DV.size()/vdof<<"\n"<<endl;
    cout<<"DVR:"<<vdof<<"x"<<rdof<<"\n"<<endl;
    Eigen::initParallel();
    Eigen::setNbThreads(8);
    //EigenSolver<MatrixXd> solver(DV);
    SelfAdjointEigenSolver<MatrixXd> solver(DV);
    d = solver.eigenvalues().real().asDiagonal();
    X = solver.eigenvectors().real();
    cout<<"eigenvalues:\n"<<d<<endl;
    cout<<"eigenvectors:\n"<<X<<endl;
}

MatrixXd KernelMatrix::calculateKernelMatrix(double t){
    MatrixXd sinF(vdof,vdof);
    for (int i = 0;i!=vdof; ++i)
	for (int j =0;j!=vdof; ++j){
        if (i==j)
		sinF(i,j)=sin(sqrt(d(i,j))*t)/sqrt(d(i,j));
	else
		sinF(i,j)=0;
	}
    return (X*sinF*(X.transpose()))*(-DVR);
}
void KernelMatrix::printKernelMatrix(){
    cout<<"Kernel Matrix\n"<<kernelMatrix<<endl;
}



