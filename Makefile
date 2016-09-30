R_LDFLAGS =     `root-config --ldflags --glibs`
R_LIBS    =     `root-config --libs`
R_CFLAGS  =     `root-config --cflags`
R_ALL     =     $(R_LDFLAGS) $(R_LIBS) $(R_CFLAGS)

FLAGS = -O0 -std=c++11
RFLAGS= ${FLAGS} ${R_ALL}

Run: Particle.o Field.o Event.o main.cxx
	g++ ${RFLAGS} Particle.o Field.o main.cxx -o Run

Particle.o: Particles/Atom.h Particles/Molecule.h Particles/Particle.cxx
	g++ -c ${FLAGS} Particles/Atom.h Particles/Molecule.h Particles/Particle.cxx
	
Field.o: Event/Field.h Event/Field.cxx
	g++ -c ${FLAGS} Event/Field.h Event/Field.cxx

Event.o: Event/Event.h Event/Event.cxx
	g++ -c ${FLAGS} Event/Event.h Event/Event.cxx

clean:
	rm *.o
