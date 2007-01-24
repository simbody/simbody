#include "AtomModel.h"
#include "internal/AtomRep.h"

using namespace SimTK;

class AtomModelRep : public AtomRep {
friend class AtomModel;
private:
	AtomModelRep(const AtomModel * handle, const Atom & atom) 
		: AtomRep(atom.getElement(), atom.getName(), handle)
	{
		defaultPosition = atom.getDefaultPosition();
	}

	AtomModelRep* clone() const {
		AtomModelRep* dup = new AtomModelRep(*this);
		return dup;
	}

	Vec3 getPosition(const State & state) const {
		return modeler->getMatter().getMobilizerPosition(state, mmBodyId).T();
	}

	int mmAtomId;
	int mmClusterId;
	int mmBodyId;
	const MoleculeModeler * modeler;
};

// "rep" mechanism requires cooperation of constructor, copy constructor, destructor, and assignment operator
AtomModel::AtomModel(const MoleculeModeler & m, const Atom & atom) 
	: Atom(atom)
{
	rep = new AtomModelRep(this, atom);
	((AtomModelRep*)rep)->modeler = & m;
}
AtomModel::AtomModel(const AtomModel & src) 
	: Atom(src)
{
	if (this == &src) return;
	rep = src.rep->clone();
	rep->myHandle = this;
}
AtomModel & AtomModel::operator=(const AtomModel & src) {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = src.rep->clone();
	rep->myHandle = this;
	return *this;
}
AtomModel::~AtomModel() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}


AtomModel & AtomModel::setMmAtomId(int id) {
	((AtomModelRep*)rep)->mmAtomId = id;
	return * this;
}
AtomModel & AtomModel::setMmClusterId(int id) {
	((AtomModelRep*)rep)->mmClusterId = id;
	return * this;
}
AtomModel & AtomModel::setMmBodyId(int id) {
	((AtomModelRep*)rep)->mmBodyId = id;
	return * this;
}

int AtomModel::getMmAtomId() const {return ((AtomModelRep*)rep)->mmAtomId;}
int AtomModel::getMmClusterId() const {return ((AtomModelRep*)rep)->mmClusterId;}
int AtomModel::getMmBodyId() const {return ((AtomModelRep*)rep)->mmBodyId;}

Vec3 AtomModel::getPosition(const State & state) const {
	return ((AtomModelRep*)rep)->getPosition(state);
}
