#include "atom.h"
#include "simulation.h"

//
// this constant is the largest number representable in the pdb format
//
const float_type Atom::INVALID_COORD = 9999.999;

//
// constructors
//

Atom::Atom() //FIX: should not return valid atom!
  : sim(Simulation::currentSimulation()), index_(-200)
{}

Atom::Atom(Simulation* sim, int index)
  : sim(sim), index_(index)
{}

Atom::Atom(const Simulation* sim, int index)
  : sim( (Simulation*)sim ), index_(index)
{}


//**
//**//
//**// note that coordinate/velocity/deriv/mass/fric/radius accessors
//**// are defined inline in the atom.hh file
//**//

const char* Atom::segmentName() const { return sim->segmentName(index_); }
const char* Atom::residueName() const { return sim->residueName(index_); }
const char* Atom::atomName() const    { return sim->atomName(index_); }
const char* Atom::chemType() const    { return sim->chemType(index_); }
int   Atom::residueNum() const        { return sim->residueNum(index_); }
Vec3  Atom::pos() const               { return sim->atomPos(index_); }
Vec3  Atom::vel() const               { return sim->atomVel(index_); }
const float_type& Atom::mass() const  { return sim->atomMass(index_); }
const float_type& Atom::fric() const  { return sim->atomFric(index_); }
//const float_type& Atom::radius() const  { return sim->atomRadius(index_); }
const float_type& Atom::charge() const  { return sim->atomCharge(index_); }

void 
Atom::setSegmentName(const char* x)
    { sim->setSegmentName(index(),x); }

void 
Atom::setAtomName(const char* x)
    { sim->setAtomName(index(),x); }

void 
Atom::setChemType(const char* x)
    { sim->setChemType(index(),x); }

void 
Atom::setResidueName(const char* x)
    { sim->setResidueName(index(),x); }

void 
Atom::setResidueNum(int x) 
    { sim->setResidueNum(index(),x); }

void
Atom::setPos(const Vec3& v)
    { sim->setAtomPos(index(),v); }

void
Atom::setVel(const Vec3& v)
    { sim->setAtomVel(index(),v); }

void
Atom::setMass(const float_type& x)
    { sim->setAtomMass(index(),x); }

void
Atom::setFric(const float_type& x)
    { sim->setAtomFric(index(),x); }

//void
//Atom::setRadius(const float_type& x)
//{ sim->setAtomRadius(index(),x); }

void
Atom::setCharge(const float_type& x)
    { sim->setAtomCharge(index(),x); }

String
Atom::string() const 
    { return sim->atomString(index()); }

ostream&
operator<<(ostream& s, const Atom& x) {
    s << x.simulation()->atomString(x.index());
    return s;
}
