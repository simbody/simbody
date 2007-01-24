#include "chemistry/Atom.h"
#include <iostream>

using namespace std;

Atom createAtom() {
	Atom atom(ChemicalElement::Silicon, "test atom");
	atom.setDefaultPosition(SimTK::Vec3(1,2,3));
	return atom;
}

int main() {
	Atom atom = createAtom();
	cout << "Position = " << atom.getDefaultPosition() << endl;
	cout << "Name = " << atom.getName() << endl;
	cout << "Mass = " << atom.getMass() << endl;
	return 0;
}
