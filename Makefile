Run: Particle.o Field.o Event.o main.cxx
	g++ -O3 Particle.o Field.o main.cxx -o Run

Particle.o: Particles/Atom.h Particles/Molecule.h Particles/Particle.cxx
	g++ -c -O3 Particles/Atom.h Particles/Molecule.h Particles/Particle.cxx
	
Field.o: Event/Field.h Event/Field.cxx
	g++ -c -O3 Event/Field.h Event/Field.cxx

Event.o: Event/Event.h Event/Event.cxx
	g++ -c -O3 Event/Event.h Event/Event.cxx

clean:
	rm *.o
