#include <iostream>
#include "Particles/Molecule.h"
#include "Event/Field.h"
#include "Event/Event.h"
#include "Event/Event.cxx"

#include "TFile.h"
#include "TTree.h"


int main(int argc, char** argv){
    int nIterations = 1;
    std::string outFile = "output.root";
    if (argc==3) {
        nIterations = atoi(argv[1]);
        outFile = argv[2];
    }
    if (argc==2) nIterations = atoi(argv[1]);

    /*TFile *file = new TFile(outFile.c_str(),"RECREATE");
    TTree *tree = new TTree("data","OCS explosion data");

    std::vector<Atom*> InitAtomVector;
    std::vector<Atom*> FinalAtomVector;
    tree->Branch("initial",&InitAtomVector);
    tree->Branch("final",&FinalAtomVector);*/
    
    for (int i=0; i<nIterations; i++){
        std::cout << "\rProgress: " << i+1 << std::flush;
//        std::cout << "================Initializing!================\n";
        Molecule *m = new Molecule();
        m->Rotate();
        m->Ionize();

        Field *f = new Field();
        EventHandler *e = new EventHandler(f, m);

        std::cout << "Potential energy of the system before explosion:\t";
        Atom*  mAtom = m->GetAtom(0);
        double E = 0;
        for (int j=0; j<m->GetNatoms(); j++){
            Atom* oAtom = m->GetAtom(j);
            if (mAtom->GetIndex() == oAtom->GetIndex()) continue;
            std::vector<double> mpos = mAtom->GetPosition();
            std::vector<double> opos = oAtom->GetPosition();
            double r = pow(mpos[0]-opos[0],2) + pow(mpos[1]-opos[1],2) + pow(mpos[2]-opos[2],2);
            r = pow(r,0.5);
            E += K_const * mAtom->GetTotalCharge() * oAtom->GetTotalCharge() / r;
        }   
        std::cout << E << std::endl; 

        double ExplosionTime = e->Run(0);

        std::cout << "Total energy of the system after explosion:\t\t";   
        E = 0;
        for (int i=0; i<m->GetNatoms(); i++){
            mAtom = m->GetAtom(i);
            std::vector<double> mvel = mAtom->GetVelocity();
            double v = pow(mvel[0],2) + pow(mvel[1],2) + pow(mvel[2],2);
            v = pow(v,.5);
            E += 0.5 * mAtom->GetMass() *pow(v,2);//* pow(v,2);
        }
        mAtom = m->GetAtom(0);
/*        for (int j=0; j<m->GetNatoms(); j++){
            Atom* oAtom = m->GetAtom(j);
            if (mAtom->GetIndex() == oAtom->GetIndex()) continue;
            std::vector<double> mpos = mAtom->GetPosition();
            std::vector<double> opos = oAtom->GetPosition();
            double r = pow(mpos[0]-opos[0],2) + pow(mpos[1]-opos[1],2) + pow(mpos[2]-opos[2],2);
            r = pow(r,0.5);
            E += K_const * mAtom->GetTotalCharge() * oAtom->GetTotalCharge() / r;
        }*/

        std::cout << E << std::endl;
continue;

        double ToF = e->Run(1);

/*        std::cout << "================Finished!================\nTime to finish (from explosion): " << ToF << " with " << e->GetNiter() << " iterations\nResults:\n";
        for (int i=0; i<m->GetNatoms(); i++){
            std::vector<double> pos = m->GetAtom(i)->GetPosition();
            char buf[500];
            sprintf(buf,"Particle: %s\tToF: %5e\tX: %5e Y:%5e Z:%5e",m->GetAtom(i)->GetName().c_str(),
                m->GetAtom(i)->GetTimeOfFlight(), pos[0], pos[1], pos[2]);
            std::cout << buf << std::endl;
        }*/

//        tree->Fill();
        delete m;
        delete f;
        delete e;
    }
//    file->Write();
}
