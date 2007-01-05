#include "chemistry/Atom.h"

using namespace SimTK;

class AtomRep {
friend class Atom;
private:
	AtomRep(const ChemicalElement & e, const SimTK::String & atomName, const Atom * handle) 
		: element(e), name(atomName), myHandle(handle) {}
	
	AtomRep* clone(Atom * newHandle) const {
		AtomRep* dup = new AtomRep(*this);
		dup->myHandle = newHandle;
		return dup;
	}

	const Atom * myHandle;
	const ChemicalElement & element;
	Vec3 defaultPosition;
	SimTK::String name;
};

Atom::Atom(const ChemicalElement & element, const SimTK::String & atomName) {
	rep = new AtomRep(element, atomName, this);
}
Atom::Atom(const Atom & src) {
	if (this == &src) return;
	rep = NULL;
	*this = src;
}
Atom & Atom::operator=(const Atom & src) {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = src.rep->clone(this);
	return *this;
}
Atom::~Atom() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

Real Atom::getMass() const {
	return rep->element.getMass();
}

const SimTK::Vec3 & Atom::getDefaultPosition() const {
	return rep->defaultPosition;
}
Atom & Atom::setDefaultPosition(SimTK::Vec3 pos) {
	rep->defaultPosition = pos;
	return *this;
}

const ChemicalElement & Atom::getElement() const {return rep->element;}
