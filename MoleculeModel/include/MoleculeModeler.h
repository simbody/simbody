#ifndef MOLECULEMODELER_H_
#define MOLECULEMODELER_H_

#include "Simbody.h"
#include "Molecule.h"
#include "AminoAcid.h"

class MoleculeModeler {
public:
	MoleculeModeler();
	~MoleculeModeler();

	MoleculeModeler & addMolecule(const Molecule & molecule);
	SimTK::MultibodySystem & system();
	SimTK::Real getPsiTorsion(const SimTK::State & state, const AminoAcid & aminoAcid);
private:
	class MoleculeModelerRep * rep;
};

#endif /*MOLECULEMODELER_H_*/
