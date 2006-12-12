#include "Atom.h"

using namespace SimTK;

class AtomRep {
friend class Atom;
private:
	AtomRep(const ChemicalElement & e, const SimTK::String & atomName, const Atom & handle) 
		: element(e), name(atomName), myHandle(&handle) {}
	
	const Atom * myHandle;
	const ChemicalElement & element;
	Vec3 position;
	SimTK::String name;
};

Atom::Atom(const ChemicalElement & element, const SimTK::String & atomName) {
	rep = new AtomRep(element, atomName, *this);
}

Real Atom::getMass() const {
	return rep->element.getMass();
}

const SimTK::Vec3 & Atom::getPosition() const {
	return rep->position;
}
Atom & Atom::setPosition(SimTK::Vec3 pos) {
	rep->position = pos;
	return *this;
}
