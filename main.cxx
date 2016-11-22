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

    std::vector<double> mass(3,0);
    std::vector<double> charge(3,0);
    std::vector<double> x(3,0);
    std::vector<double> y(3,0);
    std::vector<double> tof(3,0);
    std::vector<double> px(3,0);
    std::vector<double> py(3,0);
    std::vector<double> pz(3,0);

    tree->Branch("mass",&mass);
    tree->Branch("charge",&charge);
    tree->Branch("x",&x);
    tree->Branch("y",&y);
    tree->Branch("tof",&tof);
    tree->Branch("px",&px);
    tree->Branch("py",&py);
    tree->Branch("pz",&pz);

    for (int i=0; i<nIterations; i++){
        std::cout << "\rProgress: " << i+1 << "\tExplosion" << std::flush;
//        std::cout << "================Initializing!================\n";
        Molecule *m = new Molecule();
        m->Rotate();
        m->Ionize();

        Field *f = new Field();
        EventHandler *e = new EventHandler(f, m);
        double ExplosionTime = e->Run(0);

        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            mass[j] = atom->GetMass();
            charge[j] = atom->GetTotalCharge();
            Eigen::Vector3d vel = atom->GetVelocity();
            px[j] = vel[0] * mass[j];
            py[j] = vel[1] * mass[j];
            pz[j] = vel[2] * mass[j];
        }

//        std::cout << "Finished explosion in " << ExplosionTime << " seconds.\n";
        std::cout << "\rProgress: " << i+1 << "\tExtraction" << std::flush;
        double ToF = e->Run(1);

/*        std::cout << "================Finished!================\nTime to finish (from explosion): " << ToF << " with " << e->GetNiter() << " iterations\nResults:\n";
        for (int i=0; i<m->GetNatoms(); i++){
            std::vector<double> pos = m->GetAtom(i)->GetPosition();
            char buf[500];
            sprintf(buf,"Particle: %s\tToF: %5e\tX: %5e Y:%5e Z:%5e",m->GetAtom(i)->GetName().c_str(),
                m->GetAtom(i)->GetTimeOfFlight(), pos[0], pos[1], pos[2]);
            std::cout << buf << std::endl;
        }*/

        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            x[j] = atom->GetPosition()[0];
            y[j] = atom->GetPosition()[1];
            tof[j] = atom->GetTimeOfFlight();
        }

        tree->Fill();
        delete m;
        delete f;
        delete e;
    }
    file->Write();
    std::cout <<std::endl;
}
