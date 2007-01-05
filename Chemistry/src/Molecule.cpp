#include "chemistry/Molecule.h"
#include "chemistry/Atom.h"

using namespace SimTK;
using namespace std;

class BondRep {
public:
	BondRep(int a1, int a2, const BondType & t) 
		: atomIndex1(a1), atomIndex2(a2), type(&t) 
	{}
	BondRep(const BondRep & src) 
		: atomIndex1(src.atomIndex1), atomIndex2(src.atomIndex2), type(src.type)
	{}
	BondRep & operator=(const BondRep & src) {
		atomIndex1 = src.atomIndex1;
		atomIndex2 = src.atomIndex2;
		type = src.type;
		return *this;
	}

	int atomIndex1;
	int atomIndex2;
	const BondType * type;
};

class MoleculeRep {
friend class Molecule;
private:
	MoleculeRep(const Molecule & handle) :
	   myHandle(&handle), cachedCenterOfMass(0,0,0), centerOfMassIsDirty(true), chainId('?')
	{}

	MoleculeRep & setDefaultPosition(Real x, Real y, Real z) {
		updateCenterOfMass();
		Vec3 newCenter(x,y,z);

		if (newCenter == cachedCenterOfMass) return *this; // no change, do nothing

		Vec3 translation = newCenter - cachedCenterOfMass;

		// Update all coordinates
		for (vector<Atom>::iterator i = atomVec.begin();
			i != atomVec.end();
			++i) 
		{
			i->setDefaultPosition(i->getDefaultPosition() + translation);
		}

		cachedCenterOfMass = newCenter;
		centerOfMassIsDirty = false;

		return *this;
	}

	void updateCenterOfMass() const {

		if (centerOfMassIsDirty) {

			// Compute center of mass
			double totalMass = 0.0;
			Vec3 weightedCenter(0,0,0);
			
			for (vector<Atom>::const_iterator i = atomVec.begin();
				i != atomVec.end();
				++i) 
			{
				totalMass += i->getMass();
				weightedCenter += i->getDefaultPosition() * i->getMass();
			}
			if (totalMass > 0.0)
				cachedCenterOfMass = weightedCenter * 1.0/totalMass;

			centerOfMassIsDirty = false;
		}
		
	}

	const Molecule * myHandle;

	mutable Vec3 cachedCenterOfMass;
	mutable bool centerOfMassIsDirty;

	vector<Atom> atomVec;
	vector<BondRep> bondVec;
	char chainId; // as defined in PDB file
};

Molecule::Molecule() {
	rep = new MoleculeRep(*this);
}

Molecule::~Molecule() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

Molecule & Molecule::setDefaultPosition(Real x, Real y, Real z) {
	rep->setDefaultPosition(x,y,z);
	return *this;
}

char Molecule::getPdbChainId() const {
	return rep->chainId;
}

Molecule & Molecule::addBond(int atom1, int atom2, const BondType & bondType) {
	rep->bondVec.push_back(BondRep(atom1, atom2, bondType));
	return *this;
}

int Molecule::addAtom(const Atom & atom) {
	int index = rep->atomVec.size();
	rep->atomVec.push_back(atom);
	return index;
}
