#include <iostream>
#include "Particles/Molecule.h"
#include "Event/Field.h"
#include "Event/Event.h"
#include "Event/Event.cxx"

int main(int argc, char** argv){
    int nIterations = 1;
    if (argc==2) nIterations = atoi(argv[1]);
    
    for (int i=0; i<nIterations; i++){
        std::cout << "================Initializing!================\n";
        Molecule *m = new Molecule();
        m->Ionize();

        Field *f = new Field();
        EventHandler *e = new EventHandler(f, m);
        double ToF = e->Run();

        std::cout << "================Finished!================\nTime to finish (from explosion): " << ToF << " with " << e->GetNiter() << " iterations\nResults:\n";
        for (int i=0; i<m->GetNatoms(); i++){
            std::vector<double> pos = m->GetAtom(i)->GetPosition();
            char buf[500];
            sprintf(buf,"Particle: %s\tToF: %5e\tX: %5e Y:%5e Z:%5e",m->GetAtom(i)->GetName().c_str(),
                m->GetAtom(i)->GetTimeOfFlight(), pos[0], pos[1], pos[2]);
            std::cout << buf << std::endl;
        }
    }

}
