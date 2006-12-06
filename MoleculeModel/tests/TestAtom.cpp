#include <iostream>
#include "Atom.h"

using namespace std;

int main() {
	Atom atom(ChemicalElement::Oxygen);
	cout << "atom mass = " << atom.getMass() << endl;

	return 0;	
}
