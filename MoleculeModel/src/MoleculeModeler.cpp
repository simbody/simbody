#include "SimTKcommon.h"
#include "SimTKSimbody.h"

#include "MoleculeModeler.h"
#include "AtomModel.h"

#include <map>

using namespace SimTK;
using namespace std;

class MoleculeModelerRep {
friend class MoleculeModeler;
private:
	int nextUnusedAtomClassId;
	int nextUnusedAtomTypeId;
	bool systemIsRealized;
	// int argonTypeId;

	MoleculeModelerRep(const MoleculeModeler & handle) 
		: myHandle(&handle), nextUnusedAtomClassId(1), nextUnusedAtomTypeId(1), systemIsRealized(false)
	{
		mbs.setMatterSubsystem(matter);
		mbs.setMolecularMechanicsForceSubsystem(dumm);
		mbs.addForceSubsystem(forces);

		loadArgonParameters();
	}

	void loadArgonParameters() {
		std::string atomName = "Argon";

		// If atom is already defined, skip data entry
		if (atomStringIdMap.find(atomName) != atomStringIdMap.end())
			return;

		int classId = nextUnusedAtomClassId;
		int typeId = nextUnusedAtomTypeId;
		nextUnusedAtomClassId ++;
		nextUnusedAtomTypeId ++;

		dumm.defineAtomClass_KA(
			classId, // class ID
			atomName.c_str(), // class name
			ChemicalElement::Argon.number(), 
			0, // "valence", i.e. bond count
			3.82198, // radius 
			0.23725 // well depth
			); 
		dumm.defineChargedAtomType_KA(
			typeId, // type ID
			atomName.c_str(), // type name
			classId, // atom class ID
			0.0 // charge
			); // must be neutral by symmetry

		atomStringIdMap[atomName] = typeId;
	}
	
	int getAtomTypeId(const Atom & atom) {
		if (atom.getElement() == ChemicalElement::Argon) {
			return getAtomTypeId("Argon");
		} else {
			throw new exception("Unknown atom type");
		}
	}

	int getAtomTypeId(std::string atomName) {
		if (atomName == "Argon")
			loadArgonParameters();

		return atomStringIdMap[atomName];
	}

	AtomModel addAtomLike(const Atom & atom, const SimTK::Vec3 & position) {
		if (systemIsRealized) throw new exception("Cannot add atom after system is realized");

		AtomModel a(*myHandle, atom);
		a.setDefaultPosition(position);

		a.setMmAtomId( dumm.addAtom(getAtomTypeId(a)) ); // TODO - don't hard code on Argon
		a.setMmClusterId(dumm.createCluster(""));
		dumm.placeAtomInCluster(a.getMmAtomId(), a.getMmClusterId(), Vec3(0));

		int groundMmId = 0;
		a.setMmBodyId(
			matter.addRigidBody(
				MassProperties(atom.getMass(),Vec3(0),Inertia(0)),
				Transform(),            // inboard mobilizer frame
				groundMmId, Transform(),    // parent mobilizer frame
				Mobilizer::Cartesian
			)
		);
		dumm.attachClusterToBody(a.getMmClusterId(), a.getMmBodyId(), Transform());

		atomModels.push_back(a);

		return a;
	}

	MolecularMechanicsSystem & realizeSystem(State & state) {
		mbs.realize(state, Stage::Model);

		// Set positions
		vector<AtomModel>::iterator a = atomModels.begin();
		while (a != atomModels.end()) {
			matter.setMobilizerPosition(state, a->getMmBodyId(), a->getDefaultPosition());
			++a;
		}

		return mbs;
	}

	SimbodyMatterSubsystem   matter;
    DuMMForceFieldSubsystem  dumm;
    GeneralForceElements     forces;
	MolecularMechanicsSystem mbs;

	std::map<std::string, int> atomStringIdMap;
	std::vector<AtomModel> atomModels;

	const MoleculeModeler * myHandle;
};

MoleculeModeler::MoleculeModeler() {
	rep = new MoleculeModelerRep(*this);
}

MoleculeModeler::~MoleculeModeler() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

MolecularMechanicsSystem & MoleculeModeler::realizeSystem(State & state) {
	return rep->realizeSystem(state);
}
	
AtomModel MoleculeModeler::addAtomLike(const Atom & atom) {
	return rep->addAtomLike(atom, atom.getDefaultPosition());
}

AtomModel MoleculeModeler::addAtomLike(const Atom & atom, const SimTK::Vec3 & position) {
	return rep->addAtomLike(atom, position);
}

const SimbodyMatterSubsystem & MoleculeModeler::getMatter() const {return rep->matter;}