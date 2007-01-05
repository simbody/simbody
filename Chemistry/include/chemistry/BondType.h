#ifndef BONDTYPE_H_
#define BONDTYPE_H_

#include "SimTKSimbody.h"

// Generic bond types found in proteins
class BondType {
public:
	// Potentially freely rotating sigma bond
	static const BondType SingleBond;
	
	// Restricted rotation higher order bonds
	static const BondType SingleBondWithDoubleBondCharacter; // single bond with restricted rotation
	static const BondType SingleDoubleBond; // like in aromatic rings
	static const BondType DoubleBond;
	static const BondType TripleBond;

private:
	BondType(SimTK::Real valence);
	class BondTypeRep * rep;
};

#endif /*BOND_H_*/
