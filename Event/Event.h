#ifndef __EVENT__
#define __EVENT__
#include <vector>
#include <string>

class Field;
class Molecule;
class Atom;

class EventHandler{
  public:
    EventHandler(){ Init(NULL); };
    EventHandler(Field* aField){ Init(aField); };

    void SetField(Field* _Field){ mField = _Field; };

    void SetMolecule(Molecule* _molecule){ mMolecule = _molecule; };

    void Run();

  protected:
    void Init(Field* aField);

  private:
    Field* mField;
    Molecule* mMolecule;
    Atom* mAtom;

    unsigned int nIter;
    double time;
    double timedelta;

};

#endif
