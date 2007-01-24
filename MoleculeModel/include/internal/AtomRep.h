#ifndef _ATOMREP_H_
#define _ATOMREP_H_

#include "chemistry/Atom.h"

class AtomRep {
friend class Atom;
friend class AtomModel;
protected:

	const Atom * myHandle;

	virtual AtomRep* clone() const {
		// I am concerned about equals operator here
		AtomRep* dup = new AtomRep(*this);
		return dup;
	}

	AtomRep(const ChemicalElement & e, const SimTK::String & atomName, const Atom * handle) 
		: element(e), name(atomName), myHandle(handle) {}

	const ChemicalElement & element;
	SimTK::Vec3 defaultPosition;
	SimTK::String name;
};

#endif /* _ATOMREP_H_ */
