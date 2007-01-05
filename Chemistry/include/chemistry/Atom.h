#ifndef ATOM_H_
#define ATOM_H_

#include "ChemicalElement.h"
#include "Simmatrix.h"

class Atom {
public:
	Atom(const ChemicalElement & element, const SimTK::String & atomName = "unknown atom name");
	Atom(const Atom & src);
	Atom & operator=(const Atom & src);
	~Atom();

	const SimTK::Vec3 & getDefaultPosition() const;
	Atom & setDefaultPosition(SimTK::Vec3 pos);

	const SimTK::String & getName() const;
	Atom & setName(const SimTK::String & name);

	SimTK::Real getMass() const;

	const ChemicalElement & getElement() const;
private:
	class AtomRep * rep;
};

#endif /*ATOM_H_*/
