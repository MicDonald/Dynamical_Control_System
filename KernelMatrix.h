
#ifndef KERNELMATRIX_H
#define KERNELMATRIX_H
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "mkl.h"
#include "VRAtomInfo.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"
#include "omp.h"
#include <set>
class KernelMatrix {
public:
    KernelMatrix(){};
    KernelMatrix(const VRAtomInfo&);
    double length3D(std::vector<double> , std::vector<double>);
    Eigen::MatrixXd localKe(std::vector<double>, std::vector<double>);
    void bondNeighborIdendifier();
    std::vector<int> bondOrAtom2MatrixDof(std::vector<int>);
    std::vector<int> bondOrAtom2MatrixDof(int);
    void calculateEigen();
    Eigen::SparseMatrix<double> calculateKernelMatrix(double);
    void setK_mass(double);
    double getK_mass();
    int getvDof(){return vdof;};
    int getrDof(){return rdof;};
    int getgDof(){return gdof;};
    int getv2rDof(){return v2rdof;};
    int getr2vDof(){return r2vdof;};
    void set0state(Eigen::MatrixXd,Eigen::MatrixXd);
    VRAtomInfo model;
protected:
    int gdof,vdof,rdof,v2rdof,r2vdof;
    double k_mass;
    std::set<std::vector<int> > bonds;
    std::set<std::vector<int> > VRbonds;
    std::vector<int> Vdof,Rdof,V2Rdof,R2Vdof;
    Eigen::MatrixXd vv0;
    Eigen::MatrixXd uv0;
    Eigen::MatrixXd X;
    Eigen::MatrixXd d;
    // Eigen::MatrixXd DVR;
    Eigen::MatrixXd DVRreduced;
};
bool loadMatrix(std::string,Eigen::MatrixXd&,int);
bool saveMatrix(std::string,Eigen::MatrixXd,bool);
#endif
