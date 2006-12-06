#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "Simbody.h"

class Molecule {
public:
	Molecule();
	~Molecule();

	Molecule & setPosition(SimTK::Real x, SimTK::Real y, SimTK::Real z);

private:
	class MoleculeRep * rep;
};

#endif /*MOLECULE_H_*/
