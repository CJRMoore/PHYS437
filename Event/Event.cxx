#include <math.h>
#include "Event.h"
#include "Field.h"
#include "/home/colin/GoogleDrive/437A/Code/Particles/Molecule.h"
#include "/home/colin/GoogleDrive/437A/Code/Particles/Atom.h"
#define PI 3.14159265

double epsilon = 8.85419e-12;
double k = 1./(4. * PI * epsilon);

void EventHandler::Init(Field* aField, Molecule* aMolecule){
    mField = aField;
    mMolecule = aMolecule;

    nIter = 0;
    time = 0;
    timedelta = 1e-3;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// EnentHandler::Reset
// Prepare the Monte Carlo to be run again.
// TODO: reset and randomize the geometry of the molecule
//////////////////////////////////////////////////////////////////////////////////////////////////
void EventHandler::Reset(){
    nIter = 0;
    time = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::Run
// Loop the Runge-Kutta function until the molecule finished flag is returned (all atoms
// have reached the bottom of the detector)
//////////////////////////////////////////////////////////////////////////////////////////////////
void EventHandler::Run(){
    // Get initial conditions
    while (!mMolecule->EventFinished()){
        RungeKutta();
        nIter++;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::EfieldFromCharge
// Calculate the acceleration on the particles due to the other particles
//////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> EventHandler::EfieldFromCharge(int aAtom, double dr){
    std::vector<double> a_new(3,0);

    Atom *atom = mMolecule->GetAtom(aAtom);
    std::vector<double> atompos = atom->GetPosition();

    Atom *other = 0;
    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
        if (aAtom==iA) continue;
        other = mMolecule->GetAtom(iA);
        std::vector<double> otherpos = other->GetPosition();
        for (int iD=0; iD<2; iD++){
            double r = otherpos[iD] - atompos[iD] + dr; // distance between atoms along an axis
            if (r<1e-15) continue;
            double a = k * other->GetCharge() / pow(r,2); // updated accel.
            if (atompos[iD] < otherpos[iD]) a *= -1.;
            a_new[iD] += a;
        }
    }
    return a_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::RungeKutta
// Advance the atoms in time
// TODO fix z-position updater...it looks  like X & Y work correctly.
//////////////////////////////////////////////////////////////////////////////////////////////////
void EventHandler::RungeKutta(){
    Atom *atom = 0;
    
    std::vector<std::vector<double> > newpos(3,std::vector<double>(3,0));
    std::vector<std::vector<double> > oldpos(3,std::vector<double>(3,0));
    std::vector<std::vector<double> > newmom(3,std::vector<double>(3,0));

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        atom = mMolecule->GetAtom(iA);
        for (int iD=0; iD<3; iD++){
            newpos[iA][iD] = atom->GetPosition()[iD];
            std::cout << newpos[iA][iD] << std::endl;
            oldpos[iA][iD] = atom->GetPosition()[iD];
            newmom[iA][iD] = atom->GetMomentum()[iD];
            double v0 = atom->GetMomentum()[iD]/atom->GetMass();

            double q = atom->GetCharge();
            double k1 = UpdateDistance(v0, 0, oldpos[iA], iD, iA, q, 0);
            double k2 = UpdateDistance(v0, timedelta/2, oldpos[iA], iD, iA, q, k1);
            double k3 = UpdateDistance(v0, timedelta/2, oldpos[iA], iD, iA, q, k2);
            double k4 = UpdateDistance(v0, timedelta, oldpos[iA], iD, iA, q, k3);

            newpos[iA][iD] += timedelta/6 * (k1 + 2*k2 + 2*k3 + k4);
            newmom[iA][iD] += atom->GetMass() * 1./6 * (k1 + 2*k2 + 2*k3 + k4);
            //std::cout << iA << " " << iD << " " << newpos[iA][iD] << std::endl;
        }
    }

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
    std::cout << iA << " " << newpos[iA][0] << " " << newpos[iA][1] << " " << newpos[iA][1] << std::endl;
        atom = mMolecule->GetAtom(iA);
        atom->SetPosition(newpos[iA]);
        atom->SetMomentum(newmom[iA]);
    }
}


double EventHandler::UpdateDistance(double v0, double dt,  std::vector<double> x, int direction, int whichatom, double acharge, double k){
    x[direction] += k*dt;
    double E = EfieldFromCharge(whichatom, k*dt)[direction];
    double a = mField->GetFieldAtPosition(x)[direction] + E;

    a *= acharge;
    v0 += k;
    return v0 + a * dt;
}
