#include "Field.h"
#include <fstream>
#include <sstream>
#include <iostream>

// Dimensions of incoming field.  Defined by the matlab output field file.
int nR = 610;
int nZ = 1220;
double deltaX = 0.25e-3;
double deltaY = deltaX;
double deltaZ = 0.5e-3;

void Field::Init(std::string FieldFile){
    Potential.resize(nR,std::vector<double>(nZ,0));
    CoordinatesX.resize(nR, 0);
    CoordinatesY.resize(nR, 0);
    CoordinatesZ.resize(nZ, 0);

    std::ifstream ffile(FieldFile.c_str());
    std::string line;

    double X0 = -76.2e-3;
    double Z0 = -152.4e-3;

    int ir = 0;
    int iz = 0;
    while (std::getline(ffile,line)){
        std::stringstream ss(line);
        CoordinatesZ[iz] = Z0 + iz * deltaZ;

        double incomingval=0;
        ir=0;
        while (ss >> incomingval){
            CoordinatesX[ir] = X0 + ir * deltaX;
            CoordinatesY[ir] = X0 + ir * deltaX;
            Potential[ir][iz] = incomingval;
            ir++;
        }
        iz++;
    }
}


std::vector<double> Field::GetCoordinates(std::string _dir){
    if (_dir=="X" || _dir=="x") return CoordinatesX;
    else if (_dir=="Y" || _dir=="y") return CoordinatesY;
    else if (_dir=="Z" || _dir=="z") return CoordinatesZ;
}



Eigen::Vector3d Field::GetFieldAtPosition(double _x, double _y,  double _z){
    // Compare X
    int indexX = -1;
    int indexY = -1;
    for (int iX=0; iX<CoordinatesX.size(); iX++){
        if (indexX>0 && indexY>0) break;
        if (indexX==-1 && CoordinatesX[iX]>_x){
            indexX = iX;
        }
        if (indexY==-1 && CoordinatesY[iX]>_y){
            indexY = iX;
        }
    }

    // Do same for Z
    int indexZ = -1;
    for (int iZ=0; iZ<CoordinatesZ.size(); iZ++){
        if (CoordinatesZ[iZ]>_z){
            indexZ = iZ;
            break;
        }
    }

    Eigen::Vector3d Efield(0,0,0);
    Efield[0] = (Potential[indexX][indexZ] - Potential[indexX-1][indexZ])/deltaX;
    Efield[1] = (Potential[indexY][indexZ] - Potential[indexY-1][indexZ])/deltaY;
    Efield[2] = (Potential[indexX][indexZ] - Potential[indexX][indexZ-1])/deltaZ;
    return Efield;
}
