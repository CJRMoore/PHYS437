#include <iostream>
#include "Particles/Molecule.h"
#include "Event/Field.h"


int main(int argc, char** argv){
    Molecule *m = new Molecule();
    m->Ionize();

    Field *f = new Field();
    std::vector< std::vector<double> > field;
    field = f->GetField();
    std::cout << field[field.size()/2][field[0].size()/2] << std::endl;
    std::vector<double> efield = f->GetFieldAtPoint(0,0,50e-3);
    std::cout << efield[0] << " " << efield[1] << " " << efield[2] << std::endl;
    
}
