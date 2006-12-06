#include "Atom.h"

using namespace SimTK;

class AtomRep {
friend class Atom;
private:
	AtomRep(const ChemicalElement & e, const Atom & handle) 
		: element(e), myHandle(&handle) {}
	
	const Atom * myHandle;
	const ChemicalElement & element;
	Vec3 position;
};

Atom::Atom(const ChemicalElement & element) {
	rep = new AtomRep(element, *this);
}

Real Atom::getMass() const {
	return rep->element.getMass();
}

const SimTK::Vec3 Atom::getPosition() const {
	return rep->position;
}
Atom & Atom::setPosition(SimTK::Vec3 pos) {
	rep->position = pos;
	return *this;
}
