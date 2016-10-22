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

    std::vector<double> EfieldFromCharge(std::vector<double> k, double dt=0);
    double Run();
    void RungeKutta();
    void Reset();
    std::vector<double> UpdateDistance(std::vector<std::vector<double> > k, double dt=0);


  protected:
    void Init(Field* aField, Molecule* aMolecule);

  private:
    Field* mField;
    Molecule* mMolecule;
    Atom* mAtom;

    unsigned int mask;
    int nIter;
    double time;
    double timedelta;
    std::vector<double> errors;

    // For the adaptive step sizing.
    std::vector<unsigned int> fail;
    std::vector<double> a_ij;
    std::vector<std::vector<double> > b_ij;
    std::vector<double> c_ij;
    std::vector<double> cs_ij;

};

#endif
