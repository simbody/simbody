#include "SimTKcommon.h"
#include "Simmatrix.h"
#include "Simbody.h"

#include "chemistry/OxygenMolecule.h"
#include "MoleculeModeler.h"

using namespace SimTK;
using namespace std;

class MoleculeModelerRep {
friend class MoleculeModeler;
private:
	MoleculeModelerRep(const MoleculeModeler & handle) 
		: myHandle(&handle) 
	{
		mbs.setMatterSubsystem(matter);
		mbs.setMolecularMechanicsForceSubsystem(dumm);
		mbs.addForceSubsystem(forces);

		dumm.defineAtomClass_KA(25, "Amber99 O2", 8, 1, 1.6612, 0.2100); 
		dumm.defineChargedAtomType_KA(9999, "Sherm's O2", 25, 0); // must be neutral by symmetry
		dumm.defineBondStretch_KA(25,25, 570.0, 1.21); // bond length is right, stiffness is from C=O.

	}

	void loadArgonParameters() {
		dumm.defineAtomClass_KA(
			25, // class ID
			"Argon", // class name
			ChemicalElement::Argon.number(), 
			0, // "valence", i.e. bond count
			3.82198, // radius 
			0.23725 // well depth
			); 
		dumm.defineChargedAtomType_KA(
			25, // type ID
			"Argon", // type name
			25, // atom class ID
			0.0 // charge
			); // must be neutral by symmetry
	}
	

	AtomModel & addAtomLike(const Atom & atom, const SimTK::Vec3 & position);

	//MoleculeModelerRep & addMolecule(const Molecule & molecule) {
	//	const OxygenMolecule * oxy = dynamic_cast<const OxygenMolecule *>(&molecule);
	//	if (oxy != null) { // Is an oxygen molecule

	//		int atomTypeId = 9999;

	//		int dummO1Id = dumm.addAtom(atomTypeId);
	//		int dummO2Id = dumm.addAtom(atomTypeId);
	//		dumm.addBond(dummO1Id, dummO2Id);
	//		
	//		int cluster = dumm.createCluster("two oxygens");
	//		SimTK::Real oOBondLength = 1.21 * DuMMForceFieldSubsystem::Ang2Nm;
	//		dumm.placeAtomInCluster(dummO1Id, cluster, Vec3(0, 0, -oOBondLength * 0.5));
	//		dumm.placeAtomInCluster(dummO2Id, cluster, Vec3(0, 0,  oOBondLength * 0.5));

	//		bodies.push_back(
	//			matter.addRigidBody(
	//				MassProperties(0,Vec3(0),Inertia(0)),
	//				Transform(),            // inboard mobilizer frame
	//				parent, Transform(),    // parent mobilizer frame
	//				Mobilizer::Cartesian));
	//		// y
	//		bodies.push_back(
	//			matter.addRigidBody(
	//				MassProperties(0,Vec3(0),Inertia(0)),
	//				Transform(Rotation::aboutX(-90*Deg2Rad)),            // inboard mobilizer frame
	//				bodies.back(), Transform(Rotation::aboutX(-90*Deg2Rad)),    // parent mobilizer frame
	//				Mobilizer::Pin));
	//		// x
	//		MassProperties mprops = mm.calcClusterMassProperties(twoOxygens, Transform());
	//		cout << "Inertia:" << mprops.getInertia();
	//		cout << "inertia kludge:" << mprops.getInertia()+Inertia(0,0,.4);
	//		MassProperties mpropsKludge(mprops.getMass(), mprops.getCOM(), mprops.getInertia() + Inertia(0,0,.4));
	//		bodies.push_back(
	//			matter.addRigidBody(
	//				mpropsKludge,
	//				Transform(Rotation::aboutY(90*Deg2Rad)),            // inboard mobilizer frame
	//				bodies.back(), Transform(Rotation::aboutY(90*Deg2Rad)),    // parent mobilizer frame
	//				Mobilizer::Pin));
	//        
	//		/*
	//		bodies.push_back(
	//			matter.addRigidBody(
	//				mm.calcClusterMassProperties(twoOxygens, Transform()), 
	//				Transform(),            // inboard mobilizer frame
	//				parent, Transform(),    // parent mobilizer frame
	//				Mobilizer::FreeLine));  // TODO: DOESN'T WORK!
	//				*/
	//		mm.attachClusterToBody(twoOxygens, bodies.back(), Transform()); 

	//		// TODO
	//	}
	//	else {
	//		throw std::exception("unknown molecule type");
	//	}
	//}

	MolecularMechanicsSystem & system() {return mbs;}
	// Real getPsiTorsion(const State & state, const AminoAcid & aminoAcid);

	SimbodyMatterSubsystem   matter;
    DuMMForceFieldSubsystem  dumm;
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

//MoleculeModeler & MoleculeModeler::addMolecule(const Molecule & molecule) {
//	rep->addMolecule(molecule);
//	return *this;
//}

MultibodySystem & MoleculeModeler::getSystem() {
	return rep->system();
}
	
AtomModel & MoleculeModeler::addAtomLike(const Atom & atom) {
	return rep->addAtomLike(atom, atom.getDefaultPosition());
}

AtomModel & MoleculeModeler::addAtomLike(const Atom & atom, const SimTK::Vec3 & position) {
	return rep->addAtomLike(atom, position);
}

// Real MoleculeModeler::getPsiTorsion(const State & state, const AminoAcid & aminoAcid) {
// 	return rep->getPsiTorsion(state, aminoAcid);
// }
