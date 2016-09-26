#ifndef __EVENT__
#define __EVENT__
#include <vector>
#include <string>

class Field;
class Molecule;
class Atom;

class EventHandler{
  public:
    EventHandler(){ Init(NULL,NULL); };
    EventHandler(Field* aField, Molecule* aMolecule){ Init(aField, aMolecule); };

    unsigned int GetNiter(){ return nIter; };

    void SetField(Field* _Field){ mField = _Field; };
    void SetMolecule(Molecule* _molecule){ mMolecule = _molecule; };

    std::vector<double> EfieldFromCharge(int aAtom, double dr=0);
    void Run();
    void RungeKutta();
    void Reset();
    double UpdateDistance(double v0, double dt, std::vector<double> x, int direction, int whichatom, double acharge, double k = 0);


  protected:
    void Init(Field* aField, Molecule* aMolecule);

  private:
    Field* mField;
    Molecule* mMolecule;

    unsigned int nIter;
    double time;
    double timedelta;

};

#endif
