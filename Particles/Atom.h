#ifndef __ATOM__
#define __ATOM__
#include <vector>
#include <string>
#include "Eigen/Core"

//class Atom:
//  Describes an atom within the initial molecule.
//  Keeps track of the physical properties and position/momentum of an atom
class Atom{
  public:
    Atom();
    Atom(std::string aName, double aAtomicMass, int aAtomicCharge,  Eigen::Vector3d pos, int aIndex)
        { Init(aName, aAtomicMass, aAtomicCharge, pos, aIndex); };

    std::string GetName(){ return AtomName; };

    Eigen::Vector3d GetPosition(){ return position; };
    void SetPosition(Eigen::Vector3d _position){ position = _position; };

    Eigen::Vector3d GetVelocity(){ return velocity; };
    void SetVelocity(Eigen::Vector3d _velocity){ velocity = _velocity; };

    double GetMass(){ return mass; };
    void SetMass(double _mass) { mass = _mass; };

    double GetCharge(){ return charge; };
    double GetTotalCharge(){ return charge-nElectrons*1.602e-19; };
    void SetCharge(double _charge){ charge = _charge; };

    unsigned int GetNelectrons(){ return nElectrons; };
    void SetNelectrons(unsigned int _electrons){ nElectrons = _electrons; };

    double GetTimeOfFlight(){ return TimeOfFlight; };
    void SetTimeOfFlight(double _tof){ TimeOfFlight=_tof; };
    
    unsigned int GetIndex(){ return index; };

    double GetChargeMassRatio(){ return qm_ratio; };

  protected:
    void Init(std::string aName, double aAtomicMass, int aAtomicCharge, Eigen::Vector3d pos, int aIndex);


  private:
    std::string AtomName;
    double TimeOfFlight;
    double mass;
    double charge;
    double qm_ratio;
    unsigned short nElectrons;
    unsigned short index;

    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
};
#endif
