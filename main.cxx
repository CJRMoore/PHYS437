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

    TFile *file = new TFile(outFile.c_str(),"RECREATE");
    TTree *tree = new TTree("data","OCS explosion data");

    std::vector<Atom*> InitAtomVector;
    std::vector<Atom*> FinalAtomVector;
    tree->Branch("initial",&InitAtomVector);
    tree->Branch("final",&FinalAtomVector);
    
    for (int i=0; i<nIterations; i++){
        std::cout << "\rProgress: " << i+1 << std::flush;
//        std::cout << "================Initializing!================\n";
        Molecule *m = new Molecule();
        m->Rotate();
        m->Ionize();

        InitialAtomVector.resize(m->GetNatoms(),0);
        FinalAtomVector.resize(m->GetNatoms(),0);
        for (int i=0; i<m->GetNatoms(); i++) InitialAtomVector[i] = mAtom->GetAtom(i);

        Field *f = new Field();
        EventHandler *e = new EventHandler(f, m);
        ExplosionTime = e->Run(0);

        ToF = e->Run(1);

/*        std::cout << "================Finished!================\nTime to finish (from explosion): " << ToF << " with " << e->GetNiter() << " iterations\nResults:\n";
        for (int i=0; i<m->GetNatoms(); i++){
            std::vector<double> pos = m->GetAtom(i)->GetPosition();
            char buf[500];
            sprintf(buf,"Particle: %s\tToF: %5e\tX: %5e Y:%5e Z:%5e",m->GetAtom(i)->GetName().c_str(),
                m->GetAtom(i)->GetTimeOfFlight(), pos[0], pos[1], pos[2]);
            std::cout << buf << std::endl;
        }*/

        tree->Fill();
        delete m;
        delete f;
        delete e;
    }
    file->Write();
}
