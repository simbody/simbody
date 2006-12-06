#include "Molecule.h"
#include "Atom.h"

using namespace SimTK;
using namespace std;

class MoleculeRep {
friend class Molecule;
private:
	MoleculeRep(const Molecule & handle) :
	   myHandle(&handle), cachedCenterOfMass(0,0,0), centerOfMassIsDirty(true)
	{}

	MoleculeRep & setPosition(Real x, Real y, Real z) {
		updateCenterOfMass();
		Vec3 newCenter(x,y,z);

		if (newCenter == cachedCenterOfMass) return *this; // no change, do nothing

		Vec3 translation = newCenter - cachedCenterOfMass;

		// Update all coordinates
		for (vector<Atom>::iterator i = atomVec.begin();
			i != atomVec.end();
			++i) 
		{
			i->setPosition(i->getPosition() + translation);
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
				weightedCenter += i->getPosition() * i->getMass();
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
};

Molecule::Molecule() {
	rep = new MoleculeRep(*this);
}

Molecule::~Molecule() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

Molecule & Molecule::setPosition(Real x, Real y, Real z) {
	rep->setPosition(x,y,z);
	return *this;
}
