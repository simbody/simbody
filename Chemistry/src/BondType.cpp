#include "chemistry/BondType.h"

const BondType BondType::SingleBond(1.0);
const BondType BondType::SingleBondWithDoubleBondCharacter(1.0);
const BondType BondType::SingleDoubleBond(1.5);
const BondType BondType::DoubleBond(2.0);
const BondType BondType::TripleBond(3.0);

class BondTypeRep {
friend class BondType;
private:
	BondTypeRep(SimTK::Real v, const BondType & handle) 
		: valence(v), myHandle(&handle)
	{}

	const BondType * myHandle;
	SimTK::Real valence;
};

BondType::BondType(SimTK::Real valence) {
	rep = new BondTypeRep(valence, *this);
}