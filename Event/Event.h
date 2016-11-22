#ifndef __EVENT__
#define __EVENT__
#include <vector>
#include <string>
#include "Eigen/Core"

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

    Eigen::Vector3d EfieldFromCharge(Eigen::Vector3d atompos, double dt=0);
    double Run(int RunType);
    bool RungeKutta(int RunType);
    void Reset();
    Eigen::ArrayXXd UpdateDistance(std::vector<Eigen::ArrayXXd> k, double dt, int Runtype, int index);

    bool FinalCondition(int RunType, double Contidion=0.);


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
    std::vector<std::vector<double> > errors_vel;
    std::vector<std::vector<double> > errors_pos;

    double Energy_PreExplosion;

    // For the adaptive step sizing.
    std::vector<unsigned int> fail;
    Eigen::VectorXd a_ij;
    Eigen::MatrixXd b_ij;
    Eigen::VectorXd c_ij;
    Eigen::VectorXd cs_ij;

};

#endif
