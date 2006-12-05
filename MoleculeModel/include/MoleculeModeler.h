#ifndef MOLECULEMODELER_H_
#define MOLECULEMODELER_H_

#include "Simbody.h"
#include "Molecule.h"
#include "AminoAcid.h"
#include "simbody/internal/MultibodySystem.h"

class MoleculeModeler {
public:
	MoleculeModeler & addMolecule(const Molecule & molecule);
	SimTK::MultibodySystem & system();
	SimTK::Real getPsiTorsion(const SimTK::State & state, const AminoAcid & aminoAcid);
};

#endif /*MOLECULEMODELER_H_*/
