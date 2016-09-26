#include <iostream>
#include "Particles/Molecule.h"
#include "Event/Field.h"
#include "Event/Event.h"
#include "Event/Event.cxx"


int main(int argc, char** argv){
    Molecule *m = new Molecule();
    m->Ionize();

    std::cout << m->GetAtom(0)->GetPosition()[2] << std::endl;
    Field *f = new Field();
    EventHandler *e = new EventHandler(f, m);
    e->Run();

    std::cout << e->GetNiter() << std::endl;
    std::cout << m->GetAtom(0)->GetPosition()[2] << std::endl;
}
