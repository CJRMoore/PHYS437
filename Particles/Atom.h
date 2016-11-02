#ifndef __ATOM__
#define __ATOM__
#include <vector>
#include <string>

//class Atom:
//  Describes an atom within the initial molecule.
//  Keeps track of the physical properties and position/momentum of an atom
class Atom{
  public:
    Atom();
    Atom(std::string aName, double aAtomicMass, int aAtomicCharge, std::vector<double> pos, int aIndex){ Init(aName, aAtomicMass, aAtomicCharge, pos, aIndex); };

    std::string GetName(){ return AtomName; };

    std::vector<double> GetPosition(){ return position; };
    void SetPosition(std::vector<double> _position){ position = _position; };

    std::vector<double> GetVelocity(){ return velocity; };
    void SetVelocity(std::vector<double> _velocity){ velocity = _velocity; };

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
    void Init(std::string aName, double aAtomicMass, int aAtomicCharge, std::vector<double> pos, int aIndex);


  private:
    std::string AtomName;
    double TimeOfFlight;
    double mass;
    double charge;
    double qm_ratio;
    unsigned short nElectrons;
    unsigned short index;

    std::vector<double> position;
    std::vector<double> velocity;
};
#endif
