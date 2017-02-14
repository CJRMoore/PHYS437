#include "Field.h"
#include <fstream>
#include <sstream>
#include <iostream>

// Dimensions of incoming field.  Defined by the matlab output field file.
int nR = 0;//610;
int nZ = 0;//1220;
double deltaX = 0.25e-3;
double deltaY = deltaX;
double deltaZ = 0.5e-3;
double X0, Z0;

void Field::Init(std::string FieldFile){
    std::ifstream ffile(FieldFile.c_str());
    std::string line;
    //std::getline(ffile,line);
    //std::getline(ffile,line);

    // Get metadata from first two lines of the file.
    std::getline(ffile,line);
    std::stringstream sX(line);
    double xmax, zmax;
    sX >> X0 >> xmax >> nR;
    X0   *= 1e-3; // mm to metres
    xmax *= 1e-3;
    deltaX = (xmax-X0)/(nR-1);
    deltaY = deltaX;

    std::getline(ffile,line);
    std::stringstream sZ(line);
    sZ >> Z0 >> zmax >> nZ;
    Z0   *= 1e-3;
    zmax *= 1e-3;
    deltaZ = (zmax-Z0)/(nZ-1);

    Potential.resize(nR,std::vector<double>(nZ,0));
    CoordinatesX.resize(nR, 0);
    CoordinatesY.resize(nR, 0);
    CoordinatesZ.resize(nZ, 0);

//    double X0 = -76.2e-3;
//    double Z0 = -152.4e-3;

    int ir = 0;
    int iz = 0;
    while (std::getline(ffile,line)){
        std::stringstream ss(line);
        CoordinatesZ[iz] = Z0 + iz * deltaZ;

        double incomingval=0;
        ir=0;
        while (ss >> incomingval || !ss.eof()){
            if (ss.fail()){
                ss.clear();
                std::string dummy;
                ss >> dummy;
                incomingval = 0;
            }
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
    //int indexM = -1;
    //double _m  = pow(pow(_x,2) + pow(_y,2),0.5);
    for (int iX=0; iX<CoordinatesX.size(); iX++){
        if (indexX>0 && indexY>0)break;// && indexM>0) break;
        if (indexX==-1 && CoordinatesX[iX]>_x) indexX = iX;
        if (indexY==-1 && CoordinatesY[iX]>_y) indexY = iX;
        //if (indexM==-1 && CoordinatesX[iX]>_m) indexM = iX;
    }

    // Do same for Z
    int indexZ = -1;
    for (int iZ=0; iZ<CoordinatesZ.size(); iZ++){
        if (CoordinatesZ[iZ]>_z){
            indexZ = iZ;
            break;
        }
    }

    int indexM = (indexX + indexY) / 2;

    Eigen::Vector3d Efield(0,0,0);
    if (fabs(Potential[indexM][indexZ])>0 && fabs(Potential[indexM][indexZ-1])>0){
        Efield[0] = -(Potential[indexX][indexZ] - Potential[indexX-1][indexZ])/deltaX;
        Efield[1] = -(Potential[indexY][indexZ] - Potential[indexY-1][indexZ])/deltaY;
        Efield[2] = -(Potential[indexM][indexZ] - Potential[indexM][indexZ-1])/deltaZ;
    }
    return Efield;
}
