R_LDFLAGS =     `root-config --ldflags --glibs`
R_LIBS    =     `root-config --libs`
R_CFLAGS  =     `root-config --cflags`
R_ALL     =    $(R_LDFLAGS) $(R_LIBS) $(R_CFLAGS)

FLAGS = -O0 -std=c++11 -I/home/colin/GoogleDrive/437A/Code/  -g
RFLAGS= ${FLAGS} ${R_ALL}

Run: Particle.o Field.o Event.o main.cxx
	`root-config --cxx --cflags` -o Run Particle.o Field.o main.cxx ${FLAGS} `root-config --glibs`

Particle.o: Particles/Atom.h Particles/Molecule.h Particles/Particle.cxx
	`root-config --cxx --cflags` -c Particles/Atom.h Particles/Molecule.h Particles/Particle.cxx ${FLAGS} `root-config --glibs`
#g++ -c ${FLAGS} Particles/Atom.h Particles/Molecule.h Particles/Particle.cxx
	
Field.o: Event/Field.h Event/Field.cxx
	g++ -c ${FLAGS} Event/Field.h Event/Field.cxx

Event.o: Event/Event.h Event/Event.cxx
	g++ -c ${FLAGS} Event/Event.h Event/Event.cxx

clean:
	rm *.o
