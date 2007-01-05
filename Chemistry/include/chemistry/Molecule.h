#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "Simbody.h"
#include "Atom.h"
#include "BondType.h"

class Molecule 
{
public:
	Molecule();
	~Molecule();

	Molecule & setDefaultPosition(SimTK::Real x, SimTK::Real y, SimTK::Real z);
	char getPdbChainId() const;

	int addAtom(const Atom & atom);
	Molecule & addBond(int atom1, int atom2, const BondType & bondType);

private:
	class MoleculeRep * rep;
};

#endif /*MOLECULE_H_*/
