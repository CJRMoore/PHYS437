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
    timedelta = std::numeric_limits<double>::epsilon();//1e-10;

    fail.resize(mMolecule->GetNatoms(),0);
    mask = 0;

    // resize Runge-Kutta constant vectors and set their values
    errors.resize(2,0);
    errors[0] = 1e-16;
    errors[1] = 1e-14;

    a_ij.resize(6,0);
    b_ij.resize(6,std::vector<double>(5,0));
    c_ij.resize(6,0);
    cs_ij.resize(6,0);

    a_ij[1] = 1./5;
    a_ij[2] = 3./10;
    a_ij[3] = 3./5;
    a_ij[4] = 1.;
    a_ij[5] = 7./8;

    b_ij[1][0] = 1./5;
    b_ij[2][0] = 3./40;     b_ij[2][1] = 9./40;
    b_ij[3][0] = 3./10;     b_ij[3][1] = -9./10;    b_ij[3][2] = 6./5;
    b_ij[4][0] = -11./54;   b_ij[4][1] = 5./2;      b_ij[4][2] = -70./27;   b_ij[4][3] = 35./27;
    b_ij[5][0] = 1631./55296;   b_ij[5][1] = 175./512;  b_ij[5][2] = 575./13824;
    b_ij[5][3] = 44275./110592; b_ij[5][4] = 253./4096;

    c_ij[0] = 37./378;
    c_ij[2] = 250./621;
    c_ij[3] = 125./594;
    c_ij[5] = 512./1771;

    cs_ij[0] = 2825./27648;
    cs_ij[2] = 18575./48384;
    cs_ij[3] = 13525./55296;
    cs_ij[4] = 277./14336;
    cs_ij[5] = 1./4;
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
    std::vector<double> momentum(3,0);
    
    // Get initial conditions
    while (!mMolecule->EventFinished()){
        for (int iA=0; iA<mMolecule->GetNatoms(); iA++) fail[iA] = 0;
        RungeKutta();
        nIter++;
        time += timedelta;
    }
    return time;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::EfieldFromCharge
// Calculate the acceleration on the particles due to the other particles
//////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> EventHandler::EfieldFromCharge(std::vector<double> atompos, double dt){
//    std::vector<double> atompos = mAtom->GetPosition();
//    for (int ik=0; ik<atompos.size(); ik++) atompos[ik] += k[ik]*dt;

    Atom *other = 0;

    double referenceDistance = 1e-11;
    std::vector<std::vector<long double> > V_at_R(3,std::vector<long double>(2,0));
    std::vector<std::vector<double> > absR(2,std::vector<double>(0,0)); // distance between atoms
    std::vector<std::vector<double> > r0(3,std::vector<double>(2,0));   // distance to reference point
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
    std::vector<std::vector<double> > newvel(mMolecule->GetNatoms(),std::vector<double>(3,0));

    mask = 0;
    double delta_max = 0;
    int which_delta = -1;
    std::vector<double> momentum(3,0);

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        mAtom = mMolecule->GetAtom(iA);
        newpos[iA] = mAtom->GetPosition();
        if (mAtom->GetPosition()[Z] <= 0) {
            fail[iA] = 2;
            continue;
        }
        
        newvel[iA] = mAtom->GetVelocity();
        std::vector<double> vel = newvel[iA];
        std::vector<std::vector<double> > k_vector(0);
        k_vector.push_back(vel);

        for (int i=1; i<6; i++) 
            k_vector.push_back(UpdateDistance(k_vector, timedelta*a_ij[i]));

        std::vector<std::vector<double> > y(3,std::vector<double>(2,0));
        std::vector<std::vector<double> > v(3,std::vector<double>(2,0));
        for (int iD=0; iD<3; iD++){
            y[iD][0] = newpos[iA][iD] + c_ij[0] * k_vector[0][iD] * timedelta
                                      + c_ij[1] * k_vector[1][iD] * timedelta 
                                      + c_ij[2] * k_vector[2][iD] * timedelta
                                      + c_ij[3] * k_vector[3][iD] * timedelta
                                      + c_ij[4] * k_vector[4][iD] * timedelta
                                      + c_ij[5] * k_vector[5][iD] * timedelta;
            y[iD][1] = newpos[iA][iD] + cs_ij[0] * k_vector[0][iD] * timedelta
                                      + cs_ij[1] * k_vector[1][iD] * timedelta
                                      + cs_ij[2] * k_vector[2][iD] * timedelta
                                      + cs_ij[3] * k_vector[3][iD] * timedelta
                                      + cs_ij[4] * k_vector[4][iD] * timedelta
                                      + cs_ij[5] * k_vector[5][iD] * timedelta;
            double delta_y = fabs(y[iD][0] - y[iD][1]);

            v[iD][0] = newvel[iA][iD] + c_ij[0] * k_vector[0][iD]
                     + c_ij[1] * k_vector[1][iD]
                     + c_ij[2] * k_vector[2][iD]
                     + c_ij[3] * k_vector[3][iD]
                     + c_ij[4] * k_vector[4][iD]
                     + c_ij[5] * k_vector[5][iD];
            v[iD][1] = newvel[iA][iD] +cs_ij[0] * k_vector[0][iD]
                     + cs_ij[1] * k_vector[1][iD]
                     + cs_ij[2] * k_vector[2][iD]
                     + cs_ij[3] * k_vector[3][iD]
                     + cs_ij[4] * k_vector[4][iD]
                     + cs_ij[5] * k_vector[5][iD];
            double delta_v = fabs(v[iD][0] - v[iD][1]);
            
            newpos[iA][iD] = y[iD][0];
            newvel[iA][iD] = mAtom->GetMass() * v[iD][0];
            if (delta_y>errors[0]) fail[0] = 1;
            if (delta_v>errors[1]) fail[1] = 1;

            if (delta_y>delta_max) {
                delta_max = delta_y;
                which_delta = 0;
            }
            if (delta_v>delta_max) {
                delta_max = delta_v;
                which_delta = 1;
            }
            std::cout << iD << " " << delta_y << " " << delta_v << std::endl;
        }
        std::cout << std::endl;
    }

    //std::cout << timedelta << " " << fabs(delta_max) << std::endl;
    //std::cout << delta_max << " " << which_delta << " " << errors[which_delta] << std::endl;
//    if (which_delta>=0 && delta_max>0) timedelta *= pow(errors[which_delta]/delta_max,0.2);
//    else if (delta_max>0) timedelta *= pow(errors[1]/fabs(delta_max),0.2);

    
    //for (int iA=0; iA<mMolecule->GetNatoms(); iA++) mask |= 1 << fail[iA];
/*    switch (mask){
        case 1:
            break;
        case 2:

        case 3:
            return;
    }*/

    // check for validity of update (both that (Z+dz)>Z and that errors are withiin bounds)
    // MUST be done before next loop to maintain atomicity
    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
        mAtom = mMolecule->GetAtom(iA);
        if ((mAtom->GetPosition()[2] - newpos[iA][2])==0) mask|=1<<1;
    }
    if ((mask&2)==2) nIter--;

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        if ((mask&2)==2) continue;
        mAtom = mMolecule->GetAtom(iA);
        mAtom->SetPosition(newpos[iA]);
        mAtom->SetVelocity(newvel[iA]);
        mAtom->SetTimeOfFlight(mAtom->GetTimeOfFlight() + timedelta);
        //std::cout << newpos[iA][0] << " " << newpos[iA][1] << " " << newpos[iA][2] << " | " << newvel[iA][0] << " " << newvel[iA][1] << " " << newvel[iA][2] << std::endl;
    }
//std::cout << std::endl;
    
    //std::cout << nIter << " " << timedelta << " " << momentum[0] << " " << momentum[1] << " " << momentum[2] << std::endl;
//    timedelta *= pow(pow(10,-24)/oneofthem,0.2);

//std::cout << std::endl << std::endl;
//    std::cout << nIter << std::endl << std::endl;;
}


std::vector<double> EventHandler::UpdateDistance(std::vector<std::vector<double> > k, double dt){
    std::vector<double> position = mAtom->GetPosition();
    std::vector<double> v_return = mAtom->GetVelocity();
    int k_index = k.size()-1;
    double normalization=0;
    for (int i=0; i<k.size(); i++){ 
        for (int j=0; j<position.size(); j++) v_return[j] += k[i][j] * b_ij[k_index][i];
        for (int j=0; j<position.size(); j++) position[j] += v_return[j] * timedelta;
    }

    std::vector<double> E_q = EfieldFromCharge(position, dt);
    std::vector<double> E_s = mField->GetFieldAtPosition(position);
    std::vector<double> accel(0);
    double mass = mAtom->GetMass();
    double qm_ratio = mAtom->GetTotalCharge()/mass;
    for (int i=0; i<3; i++) {
        accel.push_back(qm_ratio * E_q[i] + qm_ratio * E_s[i]);// (E_q[i] + E_s[i]));
        v_return[i] += accel[i] * dt;
    }
//    std::cout << "\t" << mAtom->GetIndex() << " " << qm_ratio << "|" << E_s[0] << " " << v_return[0] << " " << accel[0] << "|" << E_s[1] << " " << v_return[1] << " " << accel[1] << "|" << E_s[2] << " " << v_return[2] << " " << accel[2] << std::endl;
    return v_return;
}
