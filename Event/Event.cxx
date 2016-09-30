#include <math.h>
#include <limits>
#include "Event.h"
#include "Field.h"
#include "/home/colin/GoogleDrive/437A/Code/Particles/Molecule.h"
#include "/home/colin/GoogleDrive/437A/Code/Particles/Atom.h"
#define PI 3.14159265

double epsilon = 8.85419e-12;
double K_const = 1./(4. * PI * epsilon);

int X=0;
int Y=1;
int Z=2;

void EventHandler::Init(Field* aField, Molecule* aMolecule){
    mField = aField;
    mMolecule = aMolecule;
    mAtom = 0;

    nIter = 0;
    time = 0;
    timedelta = 1e-10;
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
double EventHandler::Run(){
    // Get initial conditions
    while (!mMolecule->EventFinished()){
        RungeKutta();
        for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
            mAtom = mMolecule->GetAtom(iA);
            if (mAtom->GetPosition()[Z] <= 0 && mAtom->GetTimeOfFlight() < 0) 
                mAtom->SetTimeOfFlight(nIter*timedelta);
        }
        nIter++;
    }
    return nIter*timedelta;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::EfieldFromCharge
// Calculate the acceleration on the particles due to the other particles
//////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> EventHandler::EfieldFromCharge(std::vector<double> k, double dt){
    std::vector<double> atompos = mAtom->GetPosition();
    for (int ik=0; ik<atompos.size(); ik++) atompos[ik] += k[ik]*dt;

    double machine_eps = std::numeric_limits<double>::epsilon(); // added to denominator

    Atom *other = 0;
    std::vector<double> Evector(0);
    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
        other = mMolecule->GetAtom(iA);
        if (mAtom->GetIndex()==other->GetIndex()) continue;
        std::vector<double> otherpos = other->GetPosition();

        // Find abs. distance
        double rho = pow(atompos[X]-otherpos[X],2)+pow(atompos[Y]-otherpos[Y],2)+
            pow(atompos[Z]-otherpos[Z],2);
        double E = K_const * other->GetTotalCharge() / rho;
        rho = pow(rho,0.5);

        double theta = acos(rho/(atompos[Z]-otherpos[Z] + machine_eps));
        if (theta!=theta) theta=PI/2;
        double phi   = atan((atompos[Y]-otherpos[Y])/(atompos[X]-otherpos[X]+machine_eps));
        std::cout << rho << " " << (atompos[Z]-otherpos[Z]) << " " << theta << " " << phi << std::endl;

        Evector.push_back(E*sin(theta)*cos(phi)); //E_x
        if (atompos[X] < otherpos[X]) Evector.back() *= -1;
        Evector.push_back(E*sin(theta)*sin(phi)); //E_y
        if (atompos[Y] < otherpos[Y]) Evector.back() *= -1;
        Evector.push_back(E*cos(theta)); //E_z
        if (atompos[Z] < otherpos[Z]) Evector.back() *= -1;
        std::cout << "\t" << Evector[0] << " " << Evector[1] << " " << Evector[2] << std::endl;
    }
    return Evector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::RungeKutta
// Advance the atoms in time
// TODO fix z-position updater...it looks  like X & Y work correctly.
//////////////////////////////////////////////////////////////////////////////////////////////////
void EventHandler::RungeKutta(){
    std::vector<std::vector<double> > newpos(mMolecule->GetNatoms(),std::vector<double>(3,0));
    std::vector<std::vector<double> > newmom(mMolecule->GetNatoms(),std::vector<double>(3,0));

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        mAtom = mMolecule->GetAtom(iA);
        newpos[iA] = mAtom->GetPosition();
        if (mAtom->GetPosition()[Z] <= 0) continue;

        newmom[iA] = mAtom->GetMomentum();
        std::vector<double> vel = newmom[iA];
        for (int i=0; i<vel.size(); i++) vel[i] /= mAtom->GetMass();

        std::vector<double> k1 = UpdateDistance(vel, 0);
        std::vector<double> k2 = UpdateDistance(k1, timedelta/2);
        std::vector<double> k3 = UpdateDistance(k2, timedelta/2);
        std::vector<double> k4 = UpdateDistance(k3, timedelta);

        for (int iD=0; iD<3; iD++){
            newpos[iA][iD] += timedelta/6 * (k1[iD]+2*k2[iD]+2*k3[iD]+k4[iD]);
            newmom[iA][iD] = mAtom->GetMass()/6 * (k1[iD]+2*k2[iD]+2*k3[iD]+k4[iD]);
        }
    }

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        mAtom = mMolecule->GetAtom(iA);
        mAtom->SetPosition(newpos[iA]);
        mAtom->SetMomentum(newmom[iA]);
    }
}


std::vector<double> EventHandler::UpdateDistance(std::vector<double> k, double dt){
    std::vector<double> position = mAtom->GetPosition();

    std::vector<double> E_q = EfieldFromCharge(k, dt);
    std::vector<double> E_s = mField->GetFieldAtPosition(position);
    std::vector<double> accel(0);
    std::vector<double> v_return(3,0);
    double mass = mAtom->GetMass();
    double qm_ratio = mAtom->GetTotalCharge()/mass;
    for (int i=0; i<3; i++) {
        accel.push_back(qm_ratio * (E_q[i] + E_s[i]));
        v_return[i] = k[i] + accel[i] * dt;
    }
    return v_return;
}
