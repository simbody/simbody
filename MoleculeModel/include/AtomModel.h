#ifndef ATOMMODEL_H_
#define ATOMMODEL_H_

#include "SimTKSimbody.h"
#include "chemistry/Atom.h"
#include "MoleculeModeler.h"

/**
 * Model of a single argon atom.
 */
class AtomModel : public Atom {
public:
	AtomModel(const MoleculeModeler & modeler, const Atom & atom);
	AtomModel(const AtomModel & src);
	AtomModel & operator=(const AtomModel & src);
	~AtomModel();

	SimTK::Vec3 getPosition(const SimTK::State & state) const;

	AtomModel & setMmAtomId(int id);
	AtomModel & setMmClusterId(int id);
	AtomModel & setMmBodyId(int id);

	int getMmAtomId() const;
	int getMmClusterId() const;
	int getMmBodyId() const;

private:
	// class AtomModelRep * rep; // Don't do this
};

#endif /*ATOMMODEL_H_*/
