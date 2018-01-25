
#ifndef DYCTRLSYS_H
#define DYCTRLSYS_H
#include <map>
#include <vector>
#include <utility>
#include "Eigen/Dense"

struct DynamicalCtrlSystem {
  DynamicalCtrlSystem() {};
  DynamicalCtrlSystem(std::vector<int> atomGID,
                      std::vector< std::vector<double> > atomCoord,
                      std::vector<int> atomVirtual,
                      std::vector<int> atomReal):
    atomGID(atomGID),
    atomCoord(atomCoord),
    atomVirtual(atomVirtual),
    atomReal(atomReal) {};
  void operator()(std::vector<int> atomGID,
                  std::vector< std::vector<double> > atomCoord,
                  std::vector<int> atomVirtual,
                  std::vector<int> atomReal)
  {
    this->atomGID = atomGID;
    this->atomCoord = atomCoord;
    this->atomVirtual = atomVirtual;
    this->atomReal = atomReal;
  };

  std::vector<int> atomGID;
  std::vector< std::vector<double> > atomCoord;
  std::vector<int> atomVirtual;
  std::vector<int> atomReal;
  std::vector<int> atomV2r;
  std::vector<int> atomR2v;
  Eigen::MatrixXd atomPos;
};

#endif
