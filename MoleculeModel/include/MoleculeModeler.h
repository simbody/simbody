#ifndef MOLECULEMODELER_H_
#define MOLECULEMODELER_H_

#include "SimTKSimbody.h"
#include "chemistry/Atom.h"
#include "chemistry/Molecule.h"
#include "chemistry/AminoAcid.h"

class AtomModel;

class MoleculeModeler {
public:
	MoleculeModeler();
	~MoleculeModeler();

	AtomModel addAtomLike( const Atom & atom );
	AtomModel addAtomLike( const Atom & atom, const SimTK::Vec3 & position );

	SimTK::MolecularMechanicsSystem & realizeSystem(SimTK::State & state);
	const SimTK::SimbodyMatterSubsystem & MoleculeModeler::getMatter() const;
private:
	class MoleculeModelerRep * rep;
};

#endif /*MOLECULEMODELER_H_*/
