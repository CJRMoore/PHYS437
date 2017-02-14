#include <math.h>
#include <limits>
#include <algorithm>

#include "Event.h"
#include "Field.h"
#include "/home/colin/GoogleDrive/437A/Code/Particles/Molecule.h"
#include "/home/colin/GoogleDrive/437A/Code/Particles/Atom.h"

#include "Eigen/Core"
#include "Eigen/Geometry"

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
    mParticleEndpoint = mMolecule->GetEndpoint();

    nIter = 0;
    time = 0;
    timedelta = std::numeric_limits<double>::epsilon();//1e-10;

    fail.resize(mMolecule->GetNatoms(),0);
    mask = 0;

    // resize Runge-Kutta constant vectors and set their values
    errors_pos.resize(2);
    errors_pos[0] = std::vector<double> (3,1e-26);
    errors_pos[1] = std::vector<double> (3,1e-11);

    errors_vel.resize(2);
    errors_vel[0] = std::vector<double> (3,3e-11);
    errors_vel[1] = std::vector<double> (3,3e-3);

    a_ij.resize(6);
    b_ij.resize(6,5);
    c_ij.resize(6);
    cs_ij.resize(6);

    // Initialize matrices
    a_ij << 0,0,0,0,0,0;
    c_ij << 0,0,0,0,0,0;
    cs_ij << 0,0,0,0,0,0;
    for (int i=0; i<6; i++) for (int j=0; j<5; j++) b_ij(i,j)=0;


    a_ij[1] = 1./5;
    a_ij[2] = 3./10;
    a_ij[3] = 3./5;
    a_ij[4] = 1.;
    a_ij[5] = 7./8;

    b_ij(1,0) = 1./5;
    b_ij(2,0) = 3./40;     b_ij(2,1) = 9./40;
    b_ij(3,0) = 3./10;     b_ij(3,1) = -9./10;         b_ij(3,2) = 6./5;
    b_ij(4,0) = -11./54;   b_ij(4,1) = 5./2;           b_ij(4,2) = -70./27;   b_ij(4,3) = 35./27;
    b_ij(5,0) = 1631./55296;   b_ij(5,1) = 175./512 ;  b_ij(5,2) = 575./13824;
    b_ij(5,3) = 44275./110592; b_ij(5,4) = 253./4096;

    c_ij[0] = 37./378;
    c_ij[2] = 250./621;
    c_ij[3] = 125./594;
    c_ij[5] = 512./1771;

    cs_ij[0] = 2825./27648;
    cs_ij[2] = 18575./48384;
    cs_ij[3] = 13525./55296;
    cs_ij[4] = 277./14336;
    cs_ij[5] = 1./4;

    beta_prev = Eigen::Matrix3d::Zero();
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
void EventHandler::Run(int RunType){
    std::vector<double> momentum(3,0);

    if (!RunType) {
        Energy_PreExplosion = 0;
        for (int j=0; j<mMolecule->GetNatoms(); j++){
            mAtom = mMolecule->GetAtom(j);
            for (int k=j; k<mMolecule->GetNatoms(); k++){
                Atom* oAtom = mMolecule->GetAtom(k);
                if (oAtom->GetIndex() <= mAtom->GetIndex()) continue;
                Eigen::Vector3d mpos = mAtom->GetPosition();
                Eigen::Vector3d opos = oAtom->GetPosition();
                //double r = pow(mpos[0]-opos[0],2) + pow(mpos[1]-opos[1],2) + pow(mpos[2]-opos[2],2);
                //r = pow(r,0.5);
                mpos -= opos;
                double r = pow(mpos.dot(mpos),0.5);
                Energy_PreExplosion+= K_const * mAtom->GetTotalCharge() * oAtom->GetTotalCharge() / r;
            }
        }
    }
    //std::cout << "Potential energy of the system before explosion:\t" << Energy_PreExplosion << std::endl;

    nIter=0;
    timedelta = std::numeric_limits<double>::epsilon();
    while (RungeKutta(RunType)){
        //if (RunType) 
        nIter++;
    }

/*    if (!RunType){
        double E = 0;
        for (int i=0; i<mMolecule->GetNatoms(); i++){
            mAtom = mMolecule->GetAtom(i);
            std::vector<double> mvel = mAtom->GetVelocity();
            double vsq = pow(mvel[0],2) + pow(mvel[1],2) + pow(mvel[2],2);
//            v = pow(v,.5);
            E += 0.5 * mAtom->GetMass() * vsq;
        }
        std::cout << "Total energy of the system after explosion:\t\t" << E << std::endl;
    }*/
    //return time;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::EfieldFromCharge
// Calculate the acceleration on the particles due to the other particles
//////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3d EventHandler::EfieldFromCharge(Eigen::Vector3d atompos, double dt){
//    std::vector<double> atompos = mAtom->GetPosition();
//    for (int ik=0; ik<atompos.size(); ik++) atompos[ik] += k[ik]*dt;

    Atom *other = 0;

//    for (int i=0; i<3; i++) apos(i) = atompos[i];

    Eigen::Vector3d opos;
    Eigen::Vector3d Evector = Eigen::Vector3d::Zero();
    for (int iAtom=0; iAtom < mMolecule->GetNatoms(); iAtom++){
        other = mMolecule->GetAtom(iAtom);
        if (mAtom->GetIndex() == other->GetIndex()) continue;
        opos = other->GetPosition();
        double q = other->GetTotalCharge();

        // Lienard-Wiechert potentials
        /*Eigen::Vector3d norm = (atompos - opos);
        norm /= norm.norm();
        Eigen::Vector3d beta = other->GetVelocity() / 3e8;
        Eigen::Vector3d betadot = beta - beta_prev.col(mAtom->GetIndex()-1);
        double phi = K_const*q * 1./( (1.-norm.dot(beta)) * (atompos-opos).norm() );
        double lorentz = 1./pow(1.-beta.squaredNorm(),.5);

        Eigen::Vector3d E_intermediate = q*K_const*
          ( (norm-beta) / 
                (pow(lorentz,2) * pow(1.-norm.dot(beta),3) * (atompos - opos).squaredNorm())
           + norm.cross((norm-beta).cross(betadot)) / 
                (3e8*pow(1.-norm.dot(beta),3)*(atompos - opos).norm())
          );
          Evector += E_intermediate;*/
        Eigen::Vector3d Rvector = atompos - opos;
        double Rscalar = pow( Rvector.dot(Rvector), 0.5 );

        double AbsE = K_const * other->GetTotalCharge() / pow(Rscalar,2);
        Rvector = Rvector / Rscalar; // unit vector*/
        Evector += AbsE * Rvector;
        
    }
    return Evector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::FinalCondition
// Takes an integer (0 for Coulomb explosion, 1 for the extraction), and tells the Runge-
// Kutta algorithm whether it is complete.
//////////////////////////////////////////////////////////////////////////////////////////////////
bool EventHandler::FinalCondition(int RunType, double Condition){
    if (RunType==0){
        //return (Condition<=1e-8);
        double E = 0;
        for (int i=0; i<mMolecule->GetNatoms(); i++){
            mAtom = mMolecule->GetAtom(i);
            Eigen::Vector3d mvel = mAtom->GetVelocity() - mMolecule->GetInitialVelocity();
            double vsq = mvel.dot(mvel);//pow(mvel[0],2) + pow(mvel[1],2) + pow(mvel[2],2);
//            v = pow(v,.5);
            E += 0.5 * mAtom->GetMass() * vsq;
        }

        double E_charge = 0;
        for (int j=0; j<mMolecule->GetNatoms(); j++){
            mAtom = mMolecule->GetAtom(j);
            for (int k=j; k<mMolecule->GetNatoms(); k++){
                Atom* oAtom = mMolecule->GetAtom(k);
                if (oAtom->GetIndex() <= mAtom->GetIndex()) continue;
                Eigen::Vector3d mpos = mAtom->GetPosition();
                Eigen::Vector3d opos = oAtom->GetPosition();
                //double r = pow(mpos[0]-opos[0],2) + pow(mpos[1]-opos[1],2) + pow(mpos[2]-opos[2],2);
                mpos -= opos;
                double r = pow(mpos.dot(mpos),0.5);
                E_charge += K_const * mAtom->GetTotalCharge() * oAtom->GetTotalCharge() / r;
            }
        }

//        if (((int)Condition&2)==0) printf("%12.10f %12.10f %12.10e\n",E/Energy_PreExplosion,E_charge/Energy_PreExplosion,time);
        return (fabs(E - Energy_PreExplosion)/Energy_PreExplosion) < 0.005; 
    }
    else if (RunType){
        //std::cout << mMolecule->EventFinished() << std::endl;
        return (mMolecule->EventFinished());
    }
    else std::cerr << "Missing condition type\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// EventHandler::RungeKutta
// Advance the atoms in time
// TODO fix z-position updater...it looks  like X & Y work correctly.
//////////////////////////////////////////////////////////////////////////////////////////////////
bool EventHandler::RungeKutta(int RunType){
    std::vector<Eigen::Vector3d > newpos(mMolecule->GetNatoms(),Eigen::Vector3d (0,0,0));
    std::vector<Eigen::Vector3d > newvel(mMolecule->GetNatoms(),Eigen::Vector3d (0,0,0));

    mask = 0;
    std::vector<double> momentum(3,0);

    // For recording the error values for each position, velocity (per atom, per direction)
    //std::vector<Eigen::Vector3d > delta_y(mMolecule->GetNatoms(),Eigen::Vector3d(0,0,0));
    //std::vector<Eigen::Vector3d > delta_v(mMolecule->GetNatoms(),Eigen::Vector3d(0,0,0));
    Eigen::VectorXd delta_y = Eigen::VectorXd(mMolecule->GetNatoms());
    Eigen::VectorXd delta_v = Eigen::VectorXd(mMolecule->GetNatoms());

    std::vector<Eigen::MatrixXd> k_vector(mMolecule->GetNatoms(), Eigen::MatrixXd(6,3));
    for (int i=0; i<mMolecule->GetNatoms(); i++) k_vector[i] = Eigen::MatrixXd::Zero(6,3);

    std::vector<bool> isFinished(mMolecule->GetNatoms(),0);

    for (int iK=0; iK < 6; iK++){
        for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
            mAtom = mMolecule->GetAtom(iA);
            newpos[iA] = mAtom->GetPosition();
            if (fabs(newpos[iA](0)) > 0.07 || fabs(newpos[iA](1)) > 0.07 || (newpos[iA](2))>1){
                Eigen::Vector3d br(-1, -1, -1);
                mAtom->SetPosition(br);
            }
            if (mAtom->GetPosition()[2] <= mParticleEndpoint) {
                isFinished[iA] = 1;
                continue;
            }
        }
        UpdateDistance(k_vector, timedelta * a_ij[iK], RunType, iK);
    }

    for (int iA = 0; iA < mMolecule->GetNatoms(); iA++){
        std::vector<Eigen::Vector3d > y(2,Eigen::Vector3d (0,0,0));
        std::vector<Eigen::Vector3d > v(2,Eigen::Vector3d (0,0,0));

        mAtom = mMolecule->GetAtom(iA);
        newpos[iA] = mAtom->GetPosition();
        newvel[iA] = mAtom->GetVelocity();
        Eigen::Vector3d vel = newvel[iA];
        Eigen::Vector3d Accel_fifth=Eigen::Vector3d::Zero(), Accel_fourth=Eigen::Vector3d::Zero();
        for (int iD=0; iD<3; iD++){
            Accel_fifth(iD) = c_ij.dot(k_vector[iA].col(iD));/*c_ij[0] * k_vector[0](iA,iD)
                               + c_ij[1] * k_vector[1](iA,iD)
                               + c_ij[2] * k_vector[2](iA,iD)
                               + c_ij[3] * k_vector[3](iA,iD)
                               + c_ij[4] * k_vector[4](iA,iD)
                               + c_ij[5] * k_vector[5](iA,iD);*/

            Accel_fourth(iD) = cs_ij.dot(k_vector[iA].col(iD));/*cs_ij[0] * k_vector[0](iA,iD)
                                + cs_ij[1] * k_vector[1](iA,iD)
                                + cs_ij[2] * k_vector[2](iA,iD)
                                + cs_ij[3] * k_vector[3](iA,iD)
                                + cs_ij[4] * k_vector[4](iA,iD)
                                + cs_ij[5] * k_vector[5](iA,iD);*/
        }
            y[0] = newpos[iA] + vel * timedelta + 
                0.5 * Accel_fifth * pow(timedelta,2);
            y[1] = newpos[iA] + vel * timedelta + 
                0.5 * Accel_fourth * pow(timedelta,2);
            Eigen::Vector3d tempdy = y[1]-y[0];
            delta_y[iA] = tempdy.maxCoeff();//pow(tempdy.dot(tempdy),0.5);
            //delta_y /= pow(3,.5);
            //delta_y[iA] = y[1] - y[0];
            //delta_y[iA][iD] = fabs(y[iD][0] - y[iD][1]);

            v[0] = vel + Accel_fifth * timedelta;
            v[1] = vel + Accel_fourth * timedelta;
            Eigen::Vector3d tempdv = v[1]-v[0];
            delta_v[iA] = tempdv.maxCoeff();//pow(tempdv.dot(tempdv),0.5);
            //delta_v /= pow(3,.5);
            //delta_v[iA] = v[1] - v[0];

            newpos[iA] = y[0];
            newvel[iA] = v[0];

            //printf("%2i\t%6.4e %6.4e\n",mAtom->GetIndex(),delta_y[iA],delta_v[iA]);      
//printf("%2i\t%6.4e %6.4e %6.4e\t%6.4e %6.4e %6.4e\n",mAtom->GetIndex(),delta_y[iA][0],delta_y[iA][1],delta_y[iA][2],delta_v[iA][0],delta_v[iA][1],delta_v[iA][2]);
    }

    // check for validity of update (both that (Z+dz)>Z and that errors are withiin bounds)
    // MUST be done before next loop to maintain atomicity
    double maxnewvel = 0;
    double timefactor = 0;
    double maxPosErr = 0;
    double maxVelErr = 0;
//    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
/*        Eigen::VectorXd tempPosErrVec(mMolecule->GetNatoms());
        Eigen::VectorXd tempVelErrVec(mMolecule->GetNatoms());
        for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
            tempPosErrVec(iA) = delta_y[iA].maxCoeff();
            tempVelErrVec(iA) = delta_v[iA].maxCoeff();
        }*/
    double atomPosErr = delta_y.maxCoeff();//tempPosErrVec.maxCoeff();
    double atomVelErr = delta_v.maxCoeff();//tempVelErrVec.maxCoeff();
    if (atomPosErr/errors_pos[RunType?1:0][0] > maxPosErr)///errors_pos[RunType][0]) 
        maxPosErr = atomPosErr/errors_pos[RunType?1:0][0];
    if (atomVelErr/errors_vel[RunType?1:0][0] > maxVelErr)///errors_vel[RunType][0]) 
        maxVelErr = atomVelErr/errors_vel[RunType?1:0][0];
    if ((maxPosErr > .9) || (maxVelErr > .9)) mask |= 1 << 1;
    if ((mask&2)==2 && RunType) nIter--;
    bool finality = FinalCondition(RunType, mask);

//    if ((mask&2)==0){
        for (int iA=0; iA < mMolecule->GetNatoms(); iA++){
            //mAtom = mMolecule->GetAtom(iA);
//        printf("%2i %6.4e %2i %d %6.4e %6.4e %6.4e\n",iA,mAtom->GetTimeOfFlight(),mask,(int)isFinished[iA],mAtom->GetMass()*newvel[iA][0],mAtom->GetMass()*newvel[iA][1],mAtom->GetMass()*newvel[iA][2]);
            mAtom = mMolecule->GetAtom(iA);
            if (RunType && isFinished[iA]) {
                mAtom->SetVelocity(Eigen::Vector3d::Zero());
                beta_prev.col(iA) = Eigen::Vector3d::Zero();
                continue;
            }

            beta_prev.col(iA) = mAtom->GetVelocity();
            //printf("%2i %6.4e %6.4e %6.4e %6.4e\n",iA,timedelta,newpos[iA][0],newpos[iA][1],newpos[iA][2]);
            mAtom->SetPosition(newpos[iA]);
            mAtom->SetVelocity(newvel[iA]);
            mAtom->SetTimeOfFlight(mAtom->GetTimeOfFlight() + timedelta);
        }
//std::cout << std::endl;
    
//     std::cout << nIter << " " << timedelta << " " << momentum[0] << " " << newvel[0][0]*mMolecule->GetAtom(0)->GetMass() << " " << newvel[1][0]*mMolecule->GetAtom(1)->GetMass() << std::endl;
//momentum[1] << " " << momentum[2] << std::endl;
//    std::cout << newvel[0][0] << " " << newvel[0][1] << " " << newvel[0][2] << std::endl;
    //std::cout << maxPosErr << " " << maxVelErr << " " << timedelta << std::endl;
        time += timedelta;
  //  }

    double Factor_Power = 0.2;
    double prefactor = 1;
    if ((mask&2)==2) prefactor = 0.9;
    if ((maxPosErr>0 || maxVelErr>0)) { 
        timedelta *= prefactor*pow(1./std::max(maxPosErr,maxVelErr),Factor_Power);
    }
    else timedelta *= 2*pow(2,.5);


//    std::cout << newpos[0](0) << " " << newpos[0](1) << " " << newpos[0](2) << "\t" << newvel[0](0) << " " << newvel[0](1) << " " << newvel[0](2) << std::endl;


//    double m1 = pow(pow(momentum[0],2) + pow(momentum[1],2) + pow(momentum[2],2),0.5);
//    mAtom = mMolecule->GetAtom(0);
//    double m2 = pow( pow(mAtom->GetVelocity()[0],2) + pow(mAtom->GetVelocity()[1],2) + pow(mAtom->GetVelocity()[2],2),0.5) * mAtom->GetMass();
//    if ((mask&2)==0) std::cout << m1 << " " << m2 << std::endl;
    //printf("%6.4e\t%10.9e\t%10.9e\n",momentum[0],momentum[1],momentum[2]);



//std::cout << std::endl << std::endl;
//    std::cout << nIter << std::endl << std::endl;;
    return !finality;
}


void EventHandler::UpdateDistance(std::vector<Eigen::MatrixXd> &k, const double &dt, const int &RunType, const int &index){
//(std::vector<Eigen::MatrixXd> &k, double dt, int RunType, int index){
    //Eigen::MatrixXd accel(mMolecule->GetNatoms(),3);
    std::vector<Eigen::Vector3d > position (mMolecule->GetNatoms());
    std::vector<Eigen::Vector3d > velocity (mMolecule->GetNatoms());
    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
        mAtom = mMolecule->GetAtom(iA);
        position[iA] = mAtom->GetPosition();
        if (fabs(position[iA](0)) > 0.07 || fabs(position[iA](1)) > 0.07) continue;
        velocity[iA] = mAtom->GetVelocity();
        position[iA] += velocity[iA] * dt;
        for (int iDir=0; iDir<position.size(); iDir++){
            position[iA][iDir] += k[iA].col(iDir).head(5).dot(b_ij.row(index)) * timedelta * dt;
            velocity[iA][iDir] += k[iA].col(iDir).head(5).dot(b_ij.row(index)) * timedelta;
        }
    }

    for (int iA=0; iA<mMolecule->GetNatoms(); iA++){
        if (fabs(position[iA][0]) > 0.07 || fabs(position[iA][1]) > 0.07) continue;
        mAtom = mMolecule->GetAtom(iA);
        Eigen::Vector3d E_Field = Eigen::Vector3d::Zero();
        if (RunType==0) {
//            double Z0 = 0.5 * (0.5 * (89.61+92.91) + 0.5 * (82.51+85.81)) * 1e-3;
            E_Field = EfieldFromCharge(position[iA], dt);// + mField->GetFieldAtPosition(ptemp);
        }
        else if (RunType) E_Field = mField->GetFieldAtPosition(position[iA]);// * (1001.-RunType)/1000;
        double mass = mAtom->GetMass();
        double qm_ratio = mAtom->GetTotalCharge()/mass;
        k[iA].row(index) = qm_ratio * E_Field.transpose();    
    }
    //return accel;
}
