#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "Simbody.h"

class Molecule {
public:
	Molecule & setPosition(SimTK::Real x, SimTK::Real y, SimTK::Real z);
};

#endif /*MOLECULE_H_*/
