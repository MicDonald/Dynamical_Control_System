
#ifndef KERNELMATRIX_H
#define KERNELMATRIX_H
//#define EIGEN_USE_MKL_ALL
#include <vector>
#include <iostream>
#include <algorithm>
#include "mkl.h"
#include "VRAtomInfo.h"
#include "Eigen/Dense"
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
    //void calculateKernelMatrix(double);
    Eigen::MatrixXd calculateKernelMatrix(double);
    void printKernelMatrix();
    void setK_mass(double);
    void setDimension(int d){dimension = d;}
    double getK_mass();
    int getvDof(){return vdof;};
    int getrDof(){return rdof;};
    int getgDof(){return gdof;};
    Eigen::MatrixXd kernelMatrix;
    Eigen::MatrixXd D;
    VRAtomInfo model;
protected:
    int dimension=3;
    int gdof;
    int vdof;
    int rdof;
    double k_mass;
    std::set<std::vector<int> > bonds;
    std::set<std::vector<int> > VRbonds;
    Eigen::MatrixXd X;
    Eigen::MatrixXd d;
    Eigen::MatrixXd DVR;
};
#endif
