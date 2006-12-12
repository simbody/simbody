#ifndef ATOM_H_
#define ATOM_H_

#include "ChemicalElement.h"
#include "Simmatrix.h"

class Atom {
public:
	Atom(const ChemicalElement & element, const SimTK::String & atomName = "unknown atom name");

	const SimTK::Vec3 & getPosition() const;
	Atom & setPosition(SimTK::Vec3 pos);

	const SimTK::String & getName() const;
	Atom & setName(const SimTK::String & name);

	SimTK::Real getMass() const;
private:
	class AtomRep * rep;
};

#endif /*ATOM_H_*/
