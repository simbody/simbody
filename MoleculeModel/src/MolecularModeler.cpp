#include "SimTKcommon.h"
#include "Simmatrix.h"
#include "Simbody.h"

#include "MoleculeModeler.h"

using namespace SimTK;

class MoleculeModelerRep {
friend class MoleculeModeler;
private:
	MoleculeModelerRep(const MoleculeModeler & handle) 
		: myHandle(&handle) 
	{
		mbs.setMatterSubsystem(matter);
		mbs.setMolecularMechanicsForceSubsystem(mm);
		mbs.addForceSubsystem(forces);
	}

	MoleculeModelerRep & addMolecule(const Molecule & molecule);
	MultibodySystem & system() {return mbs;}
	Real getPsiTorsion(const State & state, const AminoAcid & aminoAcid);

	SimbodyMatterSubsystem   matter;
    DuMMForceFieldSubsystem  mm;
    GeneralForceElements     forces;

	MolecularMechanicsSystem mbs;

	const MoleculeModeler * myHandle;
};

MoleculeModeler::MoleculeModeler() {
	rep = new MoleculeModelerRep(*this);
}

MoleculeModeler::~MoleculeModeler() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

MoleculeModeler & MoleculeModeler::addMolecule(const Molecule & molecule) {
	rep->addMolecule(molecule);
	return *this;
}

MultibodySystem & MoleculeModeler::system() {
	return rep->system();
}

Real MoleculeModeler::getPsiTorsion(const State & state, const AminoAcid & aminoAcid) {
	return rep->getPsiTorsion(state, aminoAcid);
}
