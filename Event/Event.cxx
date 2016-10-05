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

    double referenceDistance = 1e-11;
    std::vector<std::vector<long double> > V_at_R(3,std::vector<long double>(2,0));
    std::vector<std::vector<double> > absR(2,std::vector<double>(0,0));
    std::vector<std::vector<double> > r0(3,std::vector<double>(2,0));
    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
        other = mMolecule->GetAtom(iA);
        if (mAtom->GetIndex()==other->GetIndex()) continue;
        std::vector<double> otherpos = other->GetPosition();
        
        for (int iD=0; iD<3; iD++){
        // Find abs. distance
            absR[0] = otherpos;
            absR[1] = otherpos;
            absR[0][iD] += referenceDistance;
            absR[1][iD] -= referenceDistance;
            for (int iD2=0; iD2<3; iD2++) {
                r0[iD][0] += pow(atompos[iD2] - absR[0][iD2],2);
                r0[iD][1] += pow(atompos[iD2] - absR[1][iD2],2);
            }
            r0[iD][0] = pow(r0[iD][0],.5);
            r0[iD][1] = pow(r0[iD][1],.5);

            V_at_R[iD][0] += K_const * other->GetCharge() / r0[iD][0];
            V_at_R[iD][1] += K_const * other->GetCharge() / r0[iD][1];
        }
    }
    std::vector<double> Evector(3,0);
    Evector[0] = (V_at_R[0][1] - V_at_R[0][0]) / (2*referenceDistance);
//        if (atompos[X] < otherpos[X]) Evector.back() *= -1;
    Evector[1] = (V_at_R[1][1] - V_at_R[1][0]) / (2*referenceDistance);
//        if (atompos[Y] < otherpos[Y]) Evector.back() *= -1;
    Evector[2] = (V_at_R[2][1] - V_at_R[2][0]) / (2*referenceDistance);
//        if (atompos[Z] < otherpos[Z]) Evector.back() *= -1;
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
     //   std::cout << std::endl;

        for (int iD=0; iD<3; iD++){
            newpos[iA][iD] += timedelta/6 * (k1[iD]+2*k2[iD]+2*k3[iD]+k4[iD]);
            newmom[iA][iD] = mAtom->GetMass()/6 * (k1[iD]+2*k2[iD]+2*k3[iD]+k4[iD]);
        }
    }

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        mAtom = mMolecule->GetAtom(iA);
        mAtom->SetPosition(newpos[iA]);
        mAtom->SetMomentum(newmom[iA]);
   //     std::cout << newpos[iA][0] << " " << newpos[iA][1] << " " << newpos[iA][2] << " | " << newmom[iA][0] << " " << newmom[iA][1] << " " << newmom[iA][2] << std::endl;
    }
//    std::cout << nIter << std::endl << std::endl;;
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
        //std::cout << " " << mAtom->GetIndex() << "\n " << E_q[i] << " " << E_s[i] << " " << qm_ratio << std::endl;
        accel.push_back(qm_ratio * E_s[i] + qm_ratio * E_q[i]);// (E_q[i] + E_s[i]));
        v_return[i] = k[i] + accel[i] * dt;
    }
    //std::cout << "\t" << mAtom->GetIndex() << " " << qm_ratio << "|" << E_q[0] << " " << accel[0] << "|" << E_q[1] << " " << accel[1] << "|" << E_q[2] << " " << accel[2] << std::endl;
    return v_return;
}
