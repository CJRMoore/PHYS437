#include "Molecule.h"
#include "Atom.h"
#include <cstdlib>
#include <random>
#include <ctime>

// Initial position of molecule before coulomb explosion.  Currently approximated by
// half of the distance between the centers of rings 4 and 5 (paper is not specific about
// initial position of molecule/z-position of laser focal point).
double Zinitial = 0.5 * (0.5 * (89.61+92.91) + 0.5 * (82.51+85.81)) * 1e-3;

// Physical constants
double Q = 1.602e-19;
double mp = 1.67e-27;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule Functions
//////////////////////////////////////////////////////////////////////////////////////////////////

void Molecule::Init(std::string aMolecule){
    MoleculeName = aMolecule;
    nAtoms = 0;
    Atoms.resize(0,0);

    // Add ability later for different molecules, for now just OCS
    AddAtom("O");
    AddAtom("C");
    //AddAtom("S");
}

// AddAtom: perform lookup of atom description (mass/charge/etc.) from atomic symbol, append
//          to end of list.
void Molecule::AddAtom(std::string _atom){
    double m=0;
    int q=0;
    double x0=0;
    if (_atom=="O"){
        m = 15.9994;
        q = 8;
        x0 = -115.78e-12;
    }
    else if (_atom=="C"){
        m = 12.0107;
        q = 6;
        x0 = 0;
    }
    else if (_atom=="S"){
        m = 32.065;
        q = 16;
        x0 = 156.01e-12;
    }
    else{
        std::cerr << "Unitientified atom: " << _atom << ".\n";
    }
    Atoms.push_back(new Atom(_atom, m, q, x0, 0, Zinitial, nAtoms+1));
    nAtoms += 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Molecule::Ionize
// Currently removes all electrons from each atom. 
// TODO: find out from Benji how the electrons will be ejected and what is expected to happen.
//////////////////////////////////////////////////////////////////////////////////////////////////
void Molecule::Ionize(){
    Atoms[0]->SetNelectrons(Atoms[0]->GetNelectrons()-1);
    Atoms[1]->SetNelectrons(Atoms[1]->GetNelectrons()-1);
    //Atoms[2]->SetNelectrons(Atoms[2]->GetNelectrons()-1);
    return;

    int total_e = 0;

    time_t seed;
    time(&seed);

    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(10,1);
    int nEjectedElectrons = (int)distribution(generator);

    std::cout << "Ejecting " << nEjectedElectrons << " electrons from the molecule.\n";
    std::vector<int> nE(nAtoms,0);

    for (unsigned int iAtom=0; iAtom<nAtoms; iAtom++) total_e += Atoms[iAtom]->GetNelectrons();
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
    }
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
void Atom::Init(std::string aName, double aAtomicMass, int aAtomicCharge, double posX, double posY, double posZ, int aIndex){
    AtomName = aName;
    mass = aAtomicMass * mp;
    charge = Q * aAtomicCharge;
    nElectrons = aAtomicCharge;
    qm_ratio = charge/mass;
    TimeOfFlight = 0;

    momentum.resize(3,0);
    position.resize(3,0);
    position[0] = posX;
    position[1] = posY;
    position[2] = posZ;

    index = aIndex;
}
