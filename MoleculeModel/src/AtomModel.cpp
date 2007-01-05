#include "AtomModel.h"

using namespace SimTK;

class AtomModelRep {
friend class AtomModel;
private:
	AtomModelRep(const AtomModel * handle) 
		: myHandle(handle) {}

	AtomModelRep* clone(AtomModel * newHandle) const {
		AtomModelRep* dup = new AtomModelRep(*this);
		dup->myHandle = newHandle;
		return dup;
	}

	Vec3 getPosition(const State & state) const {
		return modeler->getMatter().getMobilizerPosition(state, mmBodyId).T();
	}

	const AtomModel * myHandle;
	int mmAtomId;
	int mmClusterId;
	int mmBodyId;
	const MoleculeModeler * modeler;
};

// "rep" mechanism requires cooperation of constructor, copy constructor, destructor, and assignment operator
AtomModel::AtomModel(const MoleculeModeler & m, const Atom & atom) 
	: Atom(atom)
{
	rep = new AtomModelRep(this);
	rep->modeler = & m;
}
AtomModel::AtomModel(const AtomModel & src) 
	: Atom(src)
{
	if (this == &src) return;
	rep = NULL;
	*this = src;
}
AtomModel & AtomModel::operator=(const AtomModel & src) {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = src.rep->clone(this);
	return *this;
}
AtomModel::~AtomModel() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}


AtomModel & AtomModel::setMmAtomId(int id) {
	rep->mmAtomId = id;
	return * this;
}
AtomModel & AtomModel::setMmClusterId(int id) {
	rep->mmClusterId = id;
	return * this;
}
AtomModel & AtomModel::setMmBodyId(int id) {
	rep->mmBodyId = id;
	return * this;
}

int AtomModel::getMmAtomId() const {return rep->mmAtomId;}
int AtomModel::getMmClusterId() const {return rep->mmClusterId;}
int AtomModel::getMmBodyId() const {return rep->mmBodyId;}

Vec3 AtomModel::getPosition(const State & state) const {
	return rep->getPosition(state);
}
