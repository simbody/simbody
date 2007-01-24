#include "AtomModel.h"
#include <iostream>

using namespace std;

AtomModel createAtomModel() {
	MoleculeModeler modeler;
	Atom atom(ChemicalElement::Silicon, "test atom");
	atom.setDefaultPosition(SimTK::Vec3(1,2,3));
	cout << "Position1 = " << atom.getDefaultPosition() << endl;

	AtomModel atomModel(modeler, atom);
	return atomModel;
}

int main() {
	AtomModel atomModel = createAtomModel();
	cout << "Position = " << atomModel.getDefaultPosition() << endl;
	cout << "Name = " << atomModel.getName() << endl;
	cout << "Mass = " << atomModel.getMass() << endl;
	return 0;
}
