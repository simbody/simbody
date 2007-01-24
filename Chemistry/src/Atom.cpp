#include "internal/AtomRep.h"

using namespace SimTK;


Atom::Atom(const ChemicalElement & element, const SimTK::String & atomName) {
	rep = new AtomRep(element, atomName, this);
}
Atom::Atom(const Atom & src) {
	if (this == &src) return; // short circuit for self assignment
	rep = src.rep->clone(); // create new secret self
	rep->myHandle = this;
}

Atom & Atom::operator=(const Atom & src) {
	if (rep && (this == rep->myHandle)) delete rep; // destroy old secret self
	rep = src.rep->clone(); // create new secret self
	rep->myHandle = this;
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
	
const SimTK::String & Atom::getName() const {return rep->name;}
Atom & Atom::setName(const SimTK::String & name) {
	rep->name = name;
	return *this;
}
