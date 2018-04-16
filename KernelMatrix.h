
#ifndef KERNELMATRIX_H
#define KERNELMATRIX_H
//#define EIGEN_USE_MKL_ALL
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
//#include "mkl.h"
#include "AtomInfo.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"
#include <ctime>
#include "omp.h"
#include <set>
#include <map>
#include <tuple>
#include <iterator>
class KernelMatrix {
public:
    KernelMatrix(){};
    KernelMatrix(const VRAtomInfo&);
    std::vector<double> vectorBy2Points(std::vector<double>, std::vector<double>);
    double lengthOfVector(std::vector<double>);
    double dotBy2Vectors(std::vector<double> , std::vector<double>);
    double cosBy2Vectors(std::vector<double> , std::vector<double>);
    std::vector<int> bondOrAtom2MatrixDof(const std::vector<int>);
    std::vector<int> bondOrAtom2MatrixDof(const std::pair<int,int>);
    std::vector<int> bondOrAtom2MatrixDof(const int);
    void bondIdendifier();
    void angleIdendifier();
    Eigen::MatrixXd localKe(std::vector<double>, std::vector<double>);
    void calculateEigen();
    //Eigen::SparseMatrix<double> calculateKernelMatrix(double);
    Eigen::MatrixXd calculateKernelMatrix(double, bool);
    void setK_mass(double);
    void setAngle(double);
    double getK_mass();
    int getvDof(){return vdof;};
    int getrDof(){return rdof;};
    int getgDof(){return gdof;};
    int getv2rDof(){return v2rdof;};
    int getr2vDof(){return r2vdof;};
    void set0state(Eigen::MatrixXd,Eigen::MatrixXd);
    VRAtomInfo model;
    void printAngles();
    clock_t t_w=0,t_calK=0;
protected:
    int gdof,vdof,rdof,v2rdof,r2vdof;
    double k_mass, angle;  
    std::set<std::vector<int> > pairs;
    //std::set<std::vector<int> > VRbonds;
    std::multimap<int,int> bonds;
    std::multimap<int,int> VRbonds;
    std::vector<std::tuple<int, int, int>> angles;
    std::vector<int> Vdof,Rdof,V2Rdof,R2Vdof;
    Eigen::MatrixXd vv0;
    Eigen::MatrixXd uv0;
    Eigen::MatrixXd X;
    Eigen::MatrixXd d;
    Eigen::MatrixXd right;
    // Eigen::MatrixXd DVR;
    Eigen::SparseMatrix<double> DVRreduced;

};
void Test(int ,int);
bool loadMatrix(std::string,Eigen::MatrixXd&,int);
bool saveMatrix(std::string,Eigen::MatrixXd,bool);
#endif
