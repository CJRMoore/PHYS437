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
    Atom(std::string aName, double aAtomicMass, int aAtomicCharge, double posX, double posY, double posZ){ Init(aName, aAtomicMass, aAtomicCharge, posX, posY, posZ); };

    std::string GetName(){ return AtomName; };

    std::vector<double> GetPosition(){ return position; };
    void SetPosition(std::vector<double> _position){ position = _position; };

    std::vector<double> GetMomentum(){ return momentum; };
    void SetMomentum(std::vector<double> _momentum){ momentum = _momentum; };

    double GetMass(){ return mass; };
    void SetMass(double _mass) { mass = _mass; };

    double GetCharge(){ return charge; };
    void SetCharge(double _charge){ charge = _charge; };

    unsigned int GetNelectrons(){ return nelectrons; };
    void SetNelectrons(unsigned int _electrons){ nelectrons = _electrons; };

  protected:
    void Init(std::string aName, double aAtomicMass, int aAtomicCharge, double posX, double posY, double posZ);


  private:
    std::string AtomName;
    double mass;
    double charge;
    unsigned short nelectrons;

    std::vector<double> position;
    std::vector<double> momentum;
};
#endif
