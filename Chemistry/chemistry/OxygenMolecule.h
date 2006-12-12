#ifndef OXYGENMOLECULE_H_
#define OXYGENMOLECULE_H_

#include "Molecule.h"

class OxygenMolecule : public Molecule {
public:
	OxygenMolecule() {
		int o1 = addAtom(Atom(ChemicalElement::Oxygen, "O1"));
		int o2 = addAtom(Atom(ChemicalElement::Oxygen, "O2"));
		addBond(o1, o2, BondType::DoubleBond);
	}
};

#endif /*OXYGENMOLECULE_H_*/
