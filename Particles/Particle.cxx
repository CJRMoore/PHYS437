#include "Molecule.h"
#include "Atom.h"
#include <cstdlib>
#include <random>
#include <ctime>
#include <math.h>
#include <chrono>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/SVD"

#include "prob.cpp"


// Initial position of molecule before coulomb explosion.  Currently approximated by
// half of the distance between the centers of rings 4 and 5 (paper is not specific about
// initial position of molecule/z-position of laser focal point).
//double Zinitial = 0.5 * (0.5 * (89.61+92.91) + 0.5 * (82.51+85.81)) * 1e-3;

// Physical constants
double mp = 1.6726219e-27;
double pi = 3.14159265358979323;
double Q  = 1.60217662e-19 ;
double Zinitial = 0.5 * (0.5 * (89.61+92.91) + 0.5 * (82.51+85.81)) * 1e-3;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule Functions
//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::Init
// Create the atoms in the molecule
//////////////////////////////////////////////////////////////////////////////////////////////////
void Molecule::Init(std::string aMolecule, unsigned int seed){
    MoleculeName = aMolecule;
    nAtoms = 0;
    Atoms.resize(0,0);

    // Add ability later for different molecules, for now just OCS
    AddAtom("O",seed);
    AddAtom("C",seed);
    AddAtom("S",seed);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::~Molecule
// Destroy the atoms in the Atoms vector
//////////////////////////////////////////////////////////////////////////////////////////////////
Molecule::~Molecule(){
    for (int iAtom=0; iAtom<Atoms.size(); iAtom++) delete Atoms[iAtom];
}

// AddAtom: perform lookup of atom description (mass/charge/etc.) from atomic symbol, append
//          to end of list.
void Molecule::AddAtom(std::string _atom, unsigned int seed){

    //time_t seed;
    //time(&seed);
    if (seed==0) seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::normal_distribution<double> angle_dist(1.73951e+02,3.31818e+00);
    std::normal_distribution<double> OClength_dist(1.13969e+00,1.24398e-01);
    std::normal_distribution<double> CSlength_dist(1.62649e+00,1.24398e-01);
    std::uniform_real_distribution<double> dist(0,1);

    Eigen::Vector3d position;
    position << 0, 0, 0;
    //position[2] = Zinitial;

    double m=0;
    int q=0;
    if (_atom=="O"){
        m = 15.9994;
        q = 8;
        position[2] = OClength_dist(generator)*1e-10;//115.78e-12;
    }
    else if (_atom=="C"){
        m = 12.0107;
        q = 6;
    }
    else if (_atom=="S"){
        m = 32.065;
        q = 16;
        double bondlength = CSlength_dist(generator)*1e-10;//156.01e-12;
        double theta = angle_dist(generator);
        double phi = 2.*pi*dist(generator);
//        if (theta>180) theta = 360. - theta;
        theta *= pi/180.;//175. * pi / 180;
        position[2]  = bondlength * cos(theta);
        position[0]  = bondlength * sin(theta);// * cos(phi);
//        position[1]  = bondlength * sin(theta) * sin(phi);
    }
    else{
        std::cerr << "Unitientified atom: " << _atom << ".\n";
    }
    Atoms.push_back(new Atom(_atom, m, q, position, nAtoms+1));
    nAtoms += 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::GetKE()
// Calculate and return a vector of kinetic energies
//////////////////////////////////////////////////////////////////////////////////////////////////
double Molecule::GetKE(){
    double KE=0;
    for (int i=0; i<nAtoms; i++){
        double v = Atoms[i]->GetVelocity().norm();
        KE += .5 * Atoms[i]->GetMass() * pow(v,2);
    }

    return KE;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// GenerateRotation
// Create a rotation matrix
//////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix3d Molecule::GenerateRotation(unsigned seed){
    if (seed==0) seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution (0,1);
    std::uniform_int_distribution<int> int_dist(0,1);
    std::normal_distribution<double> norm_dist(0.0,100.0);

    // Gram-Schmidt Process
    /*===========================
    Eigen::Matrix3d V = Eigen::Matrix3d::Zero();
    for (int i=0; i<3; i++) for(int j=0; j<3; j++) V(i,j) = norm_dist(generator);
    Eigen::Matrix3d U = Eigen::Matrix3d::Zero();
    U.col(0) = V.col(0)/pow(V.col(0).transpose()*V.col(0),.5);
    U.col(1) = V.col(1) - (double)(V.col(1).dot(U.col(0)))/(U.col(0).dot(U.col(0)))*U.col(0);
    U.col(1) /= pow(U.col(1).dot(U.col(1)),.5);
    U.col(2) = V.col(2) - (double)(V.col(2).dot(U.col(0)))/(U.col(0).dot(U.col(0)))*U.col(0)
             - (double)(V.col(2).dot(U.col(1)))/(U.col(1).dot(U.col(1)))*U.col(1);
    U.col(2) /= pow(U.col(2).dot(U.col(2)),.5);

    Eigen::Matrix3d R;
    R = Eigen::AngleAxisd(distribution(generator)*2*pi,U.col(0)) *
        Eigen::AngleAxisd(distribution(generator) * pi,U.col(1)) *
        Eigen::AngleAxisd(distribution(generator)*2*pi,U.col(2));
    return R;
    ===============================*/


    /*===============================================================
    // Uniform rotation generated from 
    // from https://www-preview.ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf
    double theta = 2.*pi*distribution(generator) - pi;
    double phi = acos(1. - 2.*distribution(generator)) + pi/2;
    if (int_dist(generator)==0){
        if (phi<pi) phi += pi;
        else phi -= pi;
    }
    double eta = 2. * distribution(generator) - pi;

    Eigen::Matrix3d Z = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d Y = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d X = Eigen::Matrix3d::Zero();

    Z(0,0) =  cos(theta); Z(0,1) = sin(theta);
    Z(1,0) = -sin(theta); Z(1,1) = cos(theta); Z(2,2) = 1;

    Y(0,0) = cos(phi); Y(0,2) = -sin(phi);
    Y(2,0) = sin(phi); Y(2,2) =  cos(phi);     Y(1,1) = 1;

    X(1,1) =  cos(eta); X(1,2) = sin(eta);
    X(2,1) = -sin(eta); X(2,2) = cos(eta);     X(0,0) = 1;


    Eigen::Matrix3d R = Z*Y*X;
    return R;
    ===================================================*/

/*
    //Arvo's Implementation
    //https://pdfs.semanticscholar.org/04f3/beeee1ce89b9adf17a6fabde1221a328dbad.pdf
    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
    double x1 = distribution(generator)/2;
    R(0,0) = cos(2.*pi*x1); R(0,1) = sin(2.*pi*x1);
    R(1,0) =-sin(2.*pi*x1); R(1,1) = cos(2.*pi*x1);
    R(2,2) = 1;

    Eigen::Vector3d v = Eigen::Vector3d::Zero();
    double x2 = distribution(generator), x3 = distribution(generator);
    v(0) = cos(2.*pi*x2) * pow(x3,.5);
    v(1) = sin(2.*pi*x2) * pow(x3,.5);
    v(2) = pow(1.-x3,.5);

    Eigen::Matrix3d H = Eigen::Matrix3d::Identity() - 2. * v * v.transpose();
    Eigen::Matrix3d M = -1. * H * R;


    return M;
    */

    /*=================================================
    
    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();

    Eigen::Vector3d X = Eigen::Vector3d::Zero();
    X << distribution(generator),distribution(generator),distribution(generator);
    double theta1 = 2.*pi*X(1), theta2 = 2.*pi*X(2), r1 = pow(1.-X(0),.5), r2 = pow(X(0),.5);
    X /= X.norm();

    Eigen::Quaterniond q(cos(theta2)*r2,sin(theta1)*r1,cos(theta1)*r1,sin(theta2)*r2);
    q.normalize();
    R = q;
    return R;
    ==================================================*/

    /*==================================================
    // from http://www.mech.utah.edu/~brannon/public/rotation.pdf
    Eigen::Vector3d v = Eigen::Vector3d::Zero();
    v << 1,1,1;
    while (true){
        v << distribution(generator)-.5,distribution(generator)-.5,distribution(generator)-.5;
        if (v.squaredNorm()<=1) break;
    }
    v /= v.norm();
    Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(Atoms[0]->GetPosition(),v);

    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
    R = Eigen::AngleAxisd(distribution(generator)*2*pi,Eigen::Vector3d::UnitZ());
    R = q*R;
    return R;
    ===================================================*/

    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
    Eigen::Quaterniond q;
    Eigen::Vector3d v = Eigen::Vector3d::Zero();

    double u1 = distribution(generator), u2 = distribution(generator), u3 = distribution(generator);
    v(0) = pow(1.-u1,.5)*cos(2.*pi*u2);
    v(1) = pow(u1,.5) *  sin(2.*pi*u3);
    v(2) = pow(u1,.5) *  cos(2.*pi*u3);
    q.vec() = v;
    q.w() = pow(1.-u1,.5)*sin(2.*pi*u2);

    Eigen::AngleAxisd Rz(distribution(generator)*2.*pi,Eigen::Vector3d::UnitZ());
    R = Rz;
    R = Eigen::Matrix3d(q)*R;

    return R;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::Rotate
// Rotate the molecule using Euler angles
//////////////////////////////////////////////////////////////////////////////////////////////////
void Molecule::Rotate(unsigned int seed){
    // If the angles are not set in the call then use three random numbers.
    //time_t seed;
    //time(&seed);
    if (seed==0) seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution (0,1);
    std::normal_distribution<double> norm_dist(0.0,1.0);

/*
    Eigen::Vector3d Euler_vector;

        Euler_vector << 0,0,0;
        double x = distribution(generator);
        double y = distribution(generator);
        double z = distribution(generator);
        Euler_vector << x, y, z;
        Euler_vector /= pow(Euler_vector.dot(Euler_vector),0.5); // Make unit vector
    
        Eigen::Vector3d ref_vector = Atoms[0]->GetPosition();
        ref_vector /= pow(ref_vector.dot(ref_vector),0.5);
        Eigen::Vector3d cross_vector = (Eigen::Vector3d) ref_vector.cross(Euler_vector);
        Eigen::Vector3d cross_vector2 = (Eigen::Vector3d) cross_vector.cross(Euler_vector);
    
        double phi   = distribution(generator) * pi;
        double theta = distribution(generator) * pi/2;
        double psi   = distribution(generator) * pi;
        for (int iAtom=0; iAtom<nAtoms; iAtom++){
            Atom* cAtom = Atoms[iAtom];
            Eigen::Vector3d mPos = cAtom->GetPosition();
            Eigen::Vector3d rot1 = cos(phi)*mPos + sin(phi)*(Euler_vector.cross(mPos)) + (1.-cos(phi)) * (Euler_vector.dot(mPos)) * Euler_vector;
            Eigen::Vector3d rot2 = cos(psi)*rot1 + sin(psi)*(cross_vector.cross(rot1)) + (1.-cos(psi)) * (cross_vector.dot(rot1)) * cross_vector;
            Eigen::Vector3d rot3 = cos(theta)*rot2 + sin(theta)*(cross_vector2.cross(rot2)) + (1.-cos(theta)) * (cross_vector2.dot(rot2)) * cross_vector2;

            cAtom->SetPosition(rot3);
        }

    for (int iAtom=0; iAtom<nAtoms; iAtom++){
        Atom* cAtom = Atoms[iAtom];
        Eigen::Vector3d mPos = cAtom->GetPosition();
        mPos(2) += Zinitial;
        cAtom->SetPosition(mPos);
    }
*/
/*
    alpha = distribution(generator) * 2. * pi;
    beta  = distribution(generator) * pi;
    gamma = distribution(generator) * pi;

    // Set up Euler matrices (rotate Z, then Y, then Z)
    Eigen::MatrixXd R_z_alpha(3,3);
    R_z_alpha << 0,0,0,0,0,0,0,0,0;
    Eigen::MatrixXd R_x_gamma(3,3);
    R_x_gamma << 0,0,0,0,0,0,0,0,0;
    Eigen::MatrixXd R_y_beta(3,3);
    R_y_beta << 0,0,0,0,0,0,0,0,0;
    Eigen::MatrixXd R(3,3);
    R << 0,0,0,0,0,0,0,0,0;
    //std::vector<Eigen::MatrixXd> Rabg(3);

    R_z_alpha(0,0) = cos(alpha); R_z_alpha(0,1) = -sin(alpha);
    R_z_alpha(1,0) = sin(alpha); R_z_alpha(1,1) = cos(alpha);
    R_z_alpha(2,2) = 1;

    R_y_beta(0,0) = cos(beta);  R_y_beta(0,2) = sin(beta);
    R_y_beta(1,1) = 1;
    R_y_beta(2,0) = -sin(beta); R_y_beta(2,2) = cos(beta);

    R_x_gamma(0,0) = 1;
    R_x_gamma(1,1) = cos(gamma); R_x_gamma(1,2) = -sin(gamma);
    R_x_gamma(2,1) = sin(gamma); R_x_gamma(2,2) = cos(gamma);

    // Have the matrices in order.
    R = R_z_alpha * R_y_beta * R_x_gamma;
*/
/*
    Eigen::Quaterniond q;
    q.w() = norm_dist(generator);
    Eigen::Vector3d v(norm_dist(generator),norm_dist(generator),norm_dist(generator));
    q.vec() = v;
    q.normalize();
*/


/*
    double angle = 2*pi*distribution(generator);
    Eigen::Vector3d v(norm_dist(generator),norm_dist(generator),norm_dist(generator));
    v /= pow(v.dot(v),.5);
    Eigen::AngleAxisd R (angle,v);
*/


    //Eigen::Matrix3d R = GenerateRotation();
/*
    double u1 = norm_dist(generator),u2 = norm_dist(generator),u3 = norm_dist(generator);
    Eigen::Quaterniond q;
    q.w() = pow(1.-u1,0.5)*sin(2.*pi*u2);
    Eigen::Vector3d v(pow(1.-u1,0.5)*cos(2.*pi*u2),
                      pow(u1,0.5)*cos(2.*pi*u3),
                      pow(u1,0.5)*cos(2.*pi*u3));
    q.vec() = v;
    q.normalize();
*/

/*
    Eigen::Vector3d a = Atoms[0]->GetPosition();
    a /= pow(a.dot(a),0.5);
    
    Eigen::Vector3d b(norm_dist(generator),norm_dist(generator),norm_dist(generator));
    b /= pow(b.dot(b),.5);
    Eigen::Matrix3d G = Eigen::Matrix3d::Zero();
    G(0,0) =   a.dot(b);
    G(0,1) = -(a.cross(b).norm());
    G(1,0) =  (a.cross(b).norm());
    G(1,1) =   a.dot(b);
    G(2,2) =   1;

    Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
    F.col(0) = a;
    F.col(1) = (b - (a.dot(b))*a) / ((b - (a.dot(b))*a).norm());
    F.col(2) = b.cross(a);

    Eigen::Matrix3d R = F * G * F.inverse();
    //std::cout << "\n\n" << R << "\t" << R.lpNorm<2>() << "\n\n";
*/

/*
    Eigen::Vector3d b(norm_dist(generator),norm_dist(generator),norm_dist(generator));
    b /= pow(b.dot(b),.5);
    Eigen::AngleAxisd aa (2.*pi*distribution(generator),b);

    Eigen::Matrix3d R = (Eigen::Matrix3d) aa;

*/

    Eigen::Matrix3d R = GenerateRotation(seed);
    for (int iAtom=0; iAtom<nAtoms; iAtom++){
        Atom* cAtom = Atoms[iAtom];

        Eigen::Vector3d mPos = cAtom->GetPosition();
        mPos = R*mPos;
/*
        Eigen::Quaterniond p;
        p.w() = 0;
        p.vec() = cAtom->GetPosition();
        Eigen::Quaterniond p_rotated = q * p * q.inverse();
        Eigen::Vector3d mPos = p_rotated.vec();
*/
        cAtom->SetPosition(mPos);
    }
    /*
    a = Atoms[2]->GetPosition();
    std::cout << "\n" << a.norm() << " " << a.transpose() << "\t\t\t";
    a /= pow(a.dot(a),0.5);
    g = Atoms[0]->GetPosition();
    std::cout << g.norm() << " " << g.transpose();
    g /= pow(g.dot(g),0.5);
    std::cout << "\n" << acos(a.dot(g))*180/pi << std::endl; 
    */
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::Ionize
// Currently removes all electrons from each atom. 
// TODO: find out from Benji how the electrons will be ejected and what is expected to happen.
//////////////////////////////////////////////////////////////////////////////////////////////////
void Molecule::Ionize(int I1, int I2, int I3){
    Atoms[0]->SetNelectrons(Atoms[0]->GetNelectrons()-I1);
    Atoms[1]->SetNelectrons(Atoms[1]->GetNelectrons()-I2);
    Atoms[2]->SetNelectrons(Atoms[2]->GetNelectrons()-I3);
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::EventFinished
// Checks if all of the atoms in the molecule have reached the bottom of the detector.
//////////////////////////////////////////////////////////////////////////////////////////////////
bool Molecule::EventFinished(){
    bool fin = true;
    for (int iA=0; iA<nAtoms; iA++) {
   //     printf("%6.4e\t",Atoms[iA]->GetPosition()(2));
        if (Atoms[iA]->GetPosition()(2) > 0.) fin=false;
    }
 //   std::cout << std::endl;
    return fin;
}


/************************************************************************************************
** Atom Functions
************************************************************************************************/
void Atom::Init(std::string aName, double aAtomicMass, int aAtomicCharge, Eigen::Vector3d pos, int aIndex){
    AtomName = aName;
    mass = aAtomicMass * mp;
    charge = aAtomicCharge * Q;
    nElectrons = aAtomicCharge;
    qm_ratio = charge/mass;
    TimeOfFlight = 0;

    velocity << 0., 0., 0.;
    position = pos;

    index = aIndex;
}
