#include <math.h>
#include <limits>
#include <algorithm>
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
    errors_vel.resize(2,std::vector<double> (3,2.5e-3));
    errors_pos.resize(2,std::vector<double> (3,2.5e-14));

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
    
    // Run coulomb explosion 
    while (RungeKutta(0)){
        for (int iA=0; iA<mMolecule->GetNatoms(); iA++) fail[iA] = 0;
        //nIter++;
    }
    std::cout << "Finished explosion\n";

    while (!RungeKutta(1)){
        nIter++;
    }
    time -= timedelta;
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
            absR[0] = atompos;//otherpos;
            absR[1] = atompos;//otherpos;
            absR[0][iD] += referenceDistance;
            absR[1][iD] -= referenceDistance;
            for (int iD2=0; iD2<3; iD2++) {
                r0[iD][0] += pow(otherpos[iD2] - absR[0][iD2],2);
                r0[iD][1] += pow(otherpos[iD2] - absR[1][iD2],2);
            }
            r0[iD][0] = pow(r0[iD][0],.5);
            r0[iD][1] = pow(r0[iD][1],.5);

            V_at_R[iD][0] += K_const * other->GetTotalCharge() / r0[iD][0];
            V_at_R[iD][1] += K_const * other->GetTotalCharge() / r0[iD][1];
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
// EventHandler::FinalCondition
// Takes a string, "EQ" for coulomb explosion or "ES" for the static field, and tells the Runge-
// Kutta algorithm whether it is complete.
//////////////////////////////////////////////////////////////////////////////////////////////////
bool EventHandler::FinalCondition(int RunType, double Condition){
    if (RunType==0){
        return (Condition<=1e-3);
    }
    else if (RunType==1) return (!mMolecule->EventFinished());
    else std::cerr << "Missing condition type\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::RungeKutta
// Advance the atoms in time
// TODO fix z-position updater...it looks  like X & Y work correctly.
//////////////////////////////////////////////////////////////////////////////////////////////////
bool EventHandler::RungeKutta(int RunType){
    std::vector<std::vector<double> > newpos(mMolecule->GetNatoms(),std::vector<double>(3,0));
    std::vector<std::vector<double> > newvel(mMolecule->GetNatoms(),std::vector<double>(3,0));

    mask = 0;
    std::vector<double> momentum(3,0);

    // For recording the error values for each position, velocity (per atom, per direction)
    std::vector<std::vector<double> > delta_y(mMolecule->GetNatoms(),std::vector<double> (3,0));
    std::vector<std::vector<double> > delta_v(mMolecule->GetNatoms(),std::vector<double> (3,0));

    std::vector<bool> isFinished(mMolecule->GetNatoms(),0);
    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        mAtom = mMolecule->GetAtom(iA);
        newpos[iA] = mAtom->GetPosition();
        if (mAtom->GetPosition()[Z] <= 0) {
            isFinished[iA] = 1;
            continue;
        }
        
        newvel[iA] = mAtom->GetVelocity();
        std::vector<double> vel = newvel[iA];
        std::vector<std::vector<double> > k_vector(0);
        k_vector.push_back(vel);

        for (int i=1; i<6; i++) 
            k_vector.push_back(UpdateDistance(k_vector, timedelta*a_ij[i], RunType));

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
            delta_y[iA][iD] = fabs(y[iD][0] - y[iD][1]);

            v[iD][0] = c_ij[0] * k_vector[0][iD]
                     + c_ij[1] * k_vector[1][iD]
                     + c_ij[2] * k_vector[2][iD]
                     + c_ij[3] * k_vector[3][iD]
                     + c_ij[4] * k_vector[4][iD]
                     + c_ij[5] * k_vector[5][iD];
//            std::cout << k_vector[0][iD] << " " << k_vector[1][iD] << " " << k_vector[2][iD] << " " << k_vector[3][iD] << " " << k_vector[4][iD] << " " << k_vector[5][iD] << "\n";
            v[iD][1] = cs_ij[0] * k_vector[0][iD]
                     + cs_ij[1] * k_vector[1][iD]
                     + cs_ij[2] * k_vector[2][iD]
                     + cs_ij[3] * k_vector[3][iD]
                     + cs_ij[4] * k_vector[4][iD]
                     + cs_ij[5] * k_vector[5][iD];
            delta_v[iA][iD] = fabs(v[iD][0] - v[iD][1]);

            //std::cout << "\t" << v[iD][0]-newvel[iA][iD] << std::endl;
            
            newpos[iA][iD] = y[iD][0];
            newvel[iA][iD] = v[iD][0];

//            std::cout << "\t" << v[iD][0] << " " << v[iD][1] << std::endl;
//            std::cout << iD << " " << delta_y << " " << delta_v << std::endl;
        }
//        std::cout << std::endl;
//        std::cout << delta_y[iA][0] << " " << delta_y[iA][1] << " " << delta_y[iA][2] << " " << delta_v[iA][0] << " " << delta_v[iA][1] << " " << delta_v[iA][2] << "\n";
    }

    // check for validity of update (both that (Z+dz)>Z and that errors are withiin bounds)
    // MUST be done before next loop to maintain atomicity
    double maxnewvel = 0;
    double timefactor = 0;
    double maxPosErr = 0;
    double maxVelErr = 0;
    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
        double atomPosErr = *std::max_element(delta_y[iA].begin(),delta_y[iA].end());
        double atomVelErr = *std::max_element(delta_v[iA].begin(),delta_v[iA].end());
        if (atomPosErr/errors_pos[RunType][0] > maxPosErr/errors_pos[RunType][0]) 
            maxPosErr = atomPosErr/errors_pos[RunType][0];
        if (atomVelErr/errors_vel[RunType][0] > maxVelErr/errors_vel[RunType][0]) 
            maxVelErr = atomVelErr/errors_vel[RunType][0];

        mAtom = mMolecule->GetAtom(iA);
        for (int i=0; i<3; i++){
            if (maxnewvel < (fabs(mAtom->GetVelocity()[i]-newvel[iA][i]) / fabs(newvel[iA][i])))
                maxnewvel = (fabs(mAtom->GetVelocity()[i]-newvel[iA][i]) / fabs(newvel[iA][i]));
        }
    }
    if ((maxPosErr > 1) || (maxVelErr > 1)) mask |= 1 << 1;
    if ((mask&2)==2 && RunType==1) nIter--;
    bool finality = FinalCondition(RunType, maxnewvel);

    for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
        if ((mask&2)==2 || isFinished[iA]) continue;
        mAtom = mMolecule->GetAtom(iA);
        momentum[0] += mAtom->GetVelocity()[0] * mAtom->GetMass();
        momentum[1] += mAtom->GetVelocity()[1] * mAtom->GetMass();
        momentum[2] += mAtom->GetVelocity()[2] * mAtom->GetMass();

        mAtom->SetPosition(newpos[iA]);
        mAtom->SetVelocity(newvel[iA]);
        mAtom->SetTimeOfFlight(mAtom->GetTimeOfFlight() + timedelta);
//        std::cout << newpos[iA][0] << " " << newpos[iA][1] << " " << newpos[iA][2] << " | " << newvel[iA][0]*mAtom->GetMass() << " " << newvel[iA][1]*mAtom->GetMass() << " " << newvel[iA][2]*mAtom->GetMass() << std::endl;
    }
//std::cout << std::endl;
    
    //std::cout << nIter << " " << timedelta << " " << momentum[0] << " " << momentum[1] << " " << momentum[2] << std::endl;
//    std::cout << newvel[0][0] << " " << newvel[0][1] << " " << newvel[0][2] << std::endl;
//    std::cout << maxPosErr << " " << maxVelErr << std::endl;
    if ((mask&2)==0) time += timedelta;
    if (maxPosErr>0 || maxVelErr>0) timedelta *= pow(1./std::max(maxPosErr,maxVelErr),0.2);
//    std::cout << timedelta << std::endl;

//std::cout << std::endl << std::endl;
//    std::cout << nIter << std::endl << std::endl;;
    return !finality;
}


std::vector<double> EventHandler::UpdateDistance(std::vector<std::vector<double> > k, double dt, int RunType){
    std::vector<double> position = mAtom->GetPosition();
    std::vector<double> v_return = mAtom->GetVelocity();
    int k_index = k.size();
    for (int j=0; j<position.size(); j++){
        for (int i=0; i<k.size(); i++) {
            position[j] += k[i][j] * timedelta * b_ij[k_index][i];
        }
    }

    std::vector<double> E_Field;
    if (RunType==0) E_Field = EfieldFromCharge(position, dt);
    else if (RunType==1) E_Field = mField->GetFieldAtPosition(position);
    double mass = mAtom->GetMass();
    double qm_ratio = mAtom->GetTotalCharge()/mass;
    for (int i=0; i<3; i++) {
        v_return[i] += qm_ratio * E_Field[i] * dt;
    }
//    std::cout << E_Field[0] << " " << accel[0] << std::endl;
//    std::cout << accel[0] << std::endl;
//    std::cout << "\t" << mAtom->GetIndex() << " " << qm_ratio << "|" << E_s[0] << " " << v_return[0] << " " << accel[0] << "|" << E_s[1] << " " << v_return[1] << " " << accel[1] << "|" << E_s[2] << " " << v_return[2] << " " << accel[2] << std::endl;
    return v_return;
}
