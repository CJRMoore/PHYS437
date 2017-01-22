/* Particle.h
* Author: Colin Moore (2016)
* Describes a particle to be used in the Coulomb Explosion Monte Carlo.
*/

#ifndef __PARTICLE__
#define __PARTICLE__
#include <vector>
#include <string>
#include <iostream>
#include "Atom.h"
#include <Eigen/Core>

//class Molecule:
//  Container for the geometry of a molecule and the description of the atoms within the molecule.
//  Responsible for keeping track of the atoms.
class Molecule{
  public:
    Molecule(){ Init("OCS",0); }
    ~Molecule();
    Molecule(std::string aMolecule, unsigned int seed=0){ Init(aMolecule,seed); }

    unsigned short GetNatoms(){ return nAtoms; };
    Atom* GetAtom(int iAtom){
        if (iAtom>=nAtoms){
            std::cout << "Atom requested outside range\n";
            return NULL;
        }
        else return Atoms[iAtom];
    }

    void AddAtom(std::string _atom, unsigned int seed=0);
    void Ionize(int I1=1, int I2=1, int I3=1);
    Eigen::Matrix3d GenerateRotation(unsigned seed=0);
    void Rotate(unsigned int seed=0);
    void GenerateVelocity(unsigned int seed=0);
    bool EventFinished();
    double GetKE();
    Eigen::Vector3d GetInitialVelocity(){ return InitialVelocity; };

    double GetAngle(){ return bondangle; };
    std::vector<double> GetBondLengths(){
        std::vector<double> bl(2,0);
        bl[0] = pow(Atoms[0]->GetPosition().dot(Atoms[0]->GetPosition()),.5);
        bl[1] = pow(Atoms[2]->GetPosition().dot(Atoms[2]->GetPosition()),.5);
        return bl;
    }

    void Randomize();

//    std::string GetName(){ return MoleculeName; };
    
  protected:   
    void Init(std::string aMolecule, unsigned int seed=0);

  private:
    std::string MoleculeName;
    std::vector<Atom*> Atoms;
    unsigned short nAtoms;
    double bondangle;
    Eigen::Vector3d InitialVelocity;

};
#endif
