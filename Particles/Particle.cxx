#include "Molecule.h"
#include "Atom.h"
#include <cstdlib>
#include <random>
#include <ctime>
#include <math.h>

#include "Eigen/Core"

// Initial position of molecule before coulomb explosion.  Currently approximated by
// half of the distance between the centers of rings 4 and 5 (paper is not specific about
// initial position of molecule/z-position of laser focal point).
double Zinitial = 0.5 * (0.5 * (89.61+92.91) + 0.5 * (82.51+85.81)) * 1e-3;

// Physical constants
double Q = 1.60217662e-19;
double mp = 1.6726219e-27;
double pi = 3.14159265;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule Functions
//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::Init
// Create the atoms in the molecule
//////////////////////////////////////////////////////////////////////////////////////////////////
void Molecule::Init(std::string aMolecule){
    MoleculeName = aMolecule;
    nAtoms = 0;
    Atoms.resize(0,0);

    // Add ability later for different molecules, for now just OCS
    AddAtom("O");
    AddAtom("C");
    AddAtom("S");
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
void Molecule::AddAtom(std::string _atom){

    time_t seed;
    time(&seed);

    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(6,2);

    std::vector<double> position(3,0);
    position[2] = Zinitial;

    double m=0;
    int q=0;
    if (_atom=="O"){
        m = 15.9994;
        q = 8;
        position[0] = 115.78e-12;
    }
    else if (_atom=="C"){
        m = 12.0107;
        q = 6;
    }
    else if (_atom=="S"){
        m = 32.065;
        q = 16;
        double bondlength = 156.01e-12;
        double theta = 175 * pi / 180;
        position[0] = bondlength * cos(theta);
        position[2] += bondlength * sin(theta);
    }
    else{
        std::cerr << "Unitientified atom: " << _atom << ".\n";
    }
    Atoms.push_back(new Atom(_atom, m, q, position, nAtoms+1));
    nAtoms += 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::Rotate
// Rotate the molecule using Euler angles
//////////////////////////////////////////////////////////////////////////////////////////////////
void Molecule::Rotate(double alpha, double beta, double gamma){
    // If the angles are not set in the call then use three random numbers.
    if (alpha==-1 && beta==-1 && gamma==-1){
        time_t seed;
        time(&seed);
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distribution (0.0, 2*pi);
        alpha = distribution(generator);
        beta  = distribution(generator);
        gamma = distribution(generator);
    }

    std::cout << "a = " << alpha << ", b = " << beta << ", g = " << gamma << std::endl << std::endl;

    // Set up Euler matrices (rotate Z, then Y, then Z)
    std::vector<std::vector<double> > R_z_alpha(3,std::vector<double>(3,0));
    std::vector<std::vector<double> > R_z_gamma(3,std::vector<double>(3,0));
    std::vector<std::vector<double> > R_y_beta(3,std::vector<double>(3,0));
    std::vector<std::vector<std::vector<double> > > Rabg(3);

    R_z_alpha[0][0] = cos(alpha); R_z_alpha[0][1] = -sin(alpha);
    R_z_alpha[1][0] = sin(alpha); R_z_alpha[1][1] = cos(alpha);
    R_z_alpha[2][2] = 1;

    R_y_beta[0][0] = cos(beta);  R_y_beta[0][2] = sin(beta);
    R_y_beta[1][1] = 1;
    R_y_beta[2][0] = -sin(beta); R_y_beta[2][2] = cos(beta);

    R_z_gamma[0][0] = cos(gamma); R_z_gamma[0][1] = -sin(gamma);
    R_z_gamma[1][0] = sin(gamma); R_z_gamma[1][1] = cos(gamma);
    R_z_gamma[2][2] = 1;

    // Have the matrices in order.
    Rabg[0] = R_z_alpha;
    Rabg[1] = R_y_beta;
    Rabg[2] = R_z_gamma;

    std::cout << "Rotating\n";
    for (int iAtom=0; iAtom<nAtoms; iAtom++){
        Atom* cAtom = Atoms[iAtom];
        std::vector<double> pos = cAtom->GetPosition();
        pos[2] -= Zinitial; // Z-position needs to be set back to relative to an origin
        std::cout << pos[0] << " " << pos[1] << " " << pos[2] << " " << pow(pow(pos[0],2)+pow(pos[1],2)+pow(pos[2],2),0.5) << std::endl;
        for (int i=0; i<3; i++){ // Rabg loop
            for (int j=0; j<3; j++){ // First index
                pos[j] = Rabg[i][j][0] * pos[0] + Rabg[i][j][1] * pos[1] + Rabg[i][j][2] * pos[2];
//                if (i==1) std::cout << Rabg[i][j][0] * pos[0] << " " << Rabg[i][j][1] * pos[1] << " " << Rabg[i][j][2] * pos[2] << std::endl;
            }
            std::cout << pos[0] << " " << pos[1] << " " << pos[2] << " " << pow(pow(pos[0],2)+pow(pos[1],2)+pow(pos[2],2),0.5) << std::endl;
        }
        std::cout << pos[0] << " " << pos[1] << " " << pos[2] << " " << pow(pow(pos[0],2)+pow(pos[1],2)+pow(pos[2],2),0.5) << std::endl << std::endl;
        pos[2] += Zinitial;
        cAtom->SetPosition(pos);
    }
    std::cout << "Done\n";
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::Ionize
// Currently removes all electrons from each atom. 
// TODO: find out from Benji how the electrons will be ejected and what is expected to happen.
//////////////////////////////////////////////////////////////////////////////////////////////////
void Molecule::Ionize(){
/*    Atoms[0]->SetNelectrons(Atoms[0]->GetNelectrons()-1);
    Atoms[1]->SetNelectrons(Atoms[1]->GetNelectrons()-1);
    Atoms[2]->SetNelectrons(Atoms[2]->GetNelectrons()-1);
    return;*/

    int total_e = 0;

    time_t seed;
    time(&seed);

    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(6,2);
    int nEjectedElectrons = (int)distribution(generator);
    if (nEjectedElectrons<nAtoms) nEjectedElectrons=nAtoms;
    
    std::cout << "Ejecting " << nEjectedElectrons << " electrons from the molecule.\n";
    std::vector<int> nE(nAtoms,0);

    for (int iAtom=0; iAtom<nAtoms-1; iAtom++) nE[iAtom] = int(nEjectedElectrons/nAtoms);
    nE.back() = int(nEjectedElectrons - nEjectedElectrons * (nAtoms-1)/nAtoms);
    for (int i=0; i<nAtoms; i++) Atoms[i]->SetNelectrons(Atoms[i]->GetNelectrons()-nE[i]);

/*    for (unsigned int iAtom=0; iAtom<nAtoms; iAtom++) total_e += Atoms[iAtom]->GetNelectrons();
    for (int iE=0; iE<nEjectedElectrons; iE++){
        int atom = generator()%(total_e);

        int ecount = 0;
        for (int i=0; i<nAtoms; i++){
            if (ecount + Atoms[i]->GetNelectrons() > atom && Atoms[i]->GetNelectrons()>0) {
                Atoms[i]->SetNelectrons(Atoms[i]->GetNelectrons()-1);
                nE[i]++;
                break;
            }
            ecount += Atoms[i]->GetNelectrons();
        }

        total_e--;
    }*/
    for (int i=0; i<nAtoms; i++){
        std::cout << "\tElectrons ejected from " << Atoms[i]->GetName() << ": " << nE[i] << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::EventFinished
// Checks if all of the atoms in the molecule have reached the bottom of the detector.
//////////////////////////////////////////////////////////////////////////////////////////////////
bool Molecule::EventFinished(){
    bool fin = true;
    for (int iA=0; iA<nAtoms; iA++) if (Atoms[iA]->GetPosition().back() > 0.) fin=false;
    return fin;
}


/************************************************************************************************
** Atom Functions
************************************************************************************************/
void Atom::Init(std::string aName, double aAtomicMass, int aAtomicCharge, std::vector<double> pos, int aIndex){
    AtomName = aName;
    mass = aAtomicMass * mp;
    charge = Q * aAtomicCharge;
    nElectrons = aAtomicCharge;
    qm_ratio = charge/mass;
    TimeOfFlight = 0;

    velocity.resize(3,0);
    position.resize(3,0);
    position[0] = pos[0];
    position[1] = pos[1];
    position[2] = pos[2];

    index = aIndex;
}
