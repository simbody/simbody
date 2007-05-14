/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Implementation of MultibodySystem, a concrete System.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

#include "MultibodySystemRep.h"

namespace SimTK {

    /////////////////////
    // MultibodySystem //
    /////////////////////


/*static*/ bool 
MultibodySystem::isInstanceOf(const System& s) {
    return MultibodySystemRep::isA(s.getRep());
}
/*static*/ const MultibodySystem&
MultibodySystem::downcast(const System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MultibodySystem&>(s);
}
/*static*/ MultibodySystem&
MultibodySystem::updDowncast(System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MultibodySystem&>(s);
}

const MultibodySystemRep& 
MultibodySystem::getRep() const {
    return dynamic_cast<const MultibodySystemRep&>(*rep);
}
MultibodySystemRep&       
MultibodySystem::updRep() {
    return dynamic_cast<MultibodySystemRep&>(*rep);
}

// Create generic multibody system by default.
MultibodySystem::MultibodySystem() {
    rep = new MultibodySystemRep();
    rep->setMyHandle(*this);
    updRep().setGlobalSubsystem();
}

MultibodySystem::MultibodySystem(MatterSubsystem& m)
{
    rep = new MultibodySystemRep();
    rep->setMyHandle(*this);
    updRep().setGlobalSubsystem();
    setMatterSubsystem(m);
}

// This is a protected constructor for use by derived classes which
// allocate a more specialized MultibodySystemRep.
MultibodySystem::MultibodySystem(MultibodySystemRep* rp) {
    rep = rp;
    rep->setMyHandle(*this);
    updRep().setGlobalSubsystem();
}

bool MultibodySystem::project(State& s, Vector& y_err, 
             const Real& tol,
             const Real& dontProjectFac,
             const Real& targetTol
             ) const
{
    return MultibodySystemRep::downcast(*rep).project(
                s,y_err,tol,dontProjectFac,targetTol);
}


int MultibodySystem::setMatterSubsystem(MatterSubsystem& m) {
    return MultibodySystemRep::downcast(*rep).setMatterSubsystem(m);
}
int MultibodySystem::addForceSubsystem(ForceSubsystem& f) {
    return MultibodySystemRep::downcast(*rep).addForceSubsystem(f);
}
int MultibodySystem::setDecorationSubsystem(DecorationSubsystem& m) {
    return MultibodySystemRep::downcast(*rep).setDecorationSubsystem(m);
}

const MatterSubsystem&       
MultibodySystem::getMatterSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getMatterSubsystem();
}
MatterSubsystem&       
MultibodySystem::updMatterSubsystem() {
    return MultibodySystemRep::downcast(*rep).updMatterSubsystem();
}

const DecorationSubsystem&       
MultibodySystem::getDecorationSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getDecorationSubsystem();
}
DecorationSubsystem&       
MultibodySystem::updDecorationSubsystem() {
    return MultibodySystemRep::downcast(*rep).updDecorationSubsystem();
}

const Real&                
MultibodySystem::getPotentialEnergy(const State& s) const {
    return getRep().getPotentialEnergy(s);
}
const Real&                
MultibodySystem::getKineticEnergy(const State& s) const {
    return getRep().getKineticEnergy(s);
}

const Vector_<SpatialVec>& 
MultibodySystem::getRigidBodyForces(const State& s) const {
    return getRep().getRigidBodyForces(s);
}
const Vector_<Vec3>&       
MultibodySystem::getParticleForces(const State& s) const {
    return getRep().getParticleForces(s);
}
const Vector&              
MultibodySystem::getMobilityForces(const State& s) const {
    return getRep().getMobilityForces(s);
}

Real&                
MultibodySystem::updPotentialEnergy(const State& s) const {
    return getRep().updPotentialEnergy(s);
}
Real&                
MultibodySystem::updKineticEnergy(const State& s) const {
    return getRep().updKineticEnergy(s);
}

Vector_<SpatialVec>& 
MultibodySystem::updRigidBodyForces(const State& s) const {
    return getRep().updRigidBodyForces(s);
}
Vector_<Vec3>&       
MultibodySystem::updParticleForces(const State& s) const {
    return getRep().updParticleForces(s);
}
Vector&              
MultibodySystem::updMobilityForces(const State& s) const {
    return getRep().updMobilityForces(s);
}


    ////////////////////////////////////
    // MultibodySystemGlobalSubsystem //
    ////////////////////////////////////


/*static*/ bool 
MultibodySystemGlobalSubsystem::isInstanceOf(const Subsystem& s) {
    return MultibodySystemGlobalSubsystemRep::isA(s.getRep());
}
/*static*/ const MultibodySystemGlobalSubsystem&
MultibodySystemGlobalSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MultibodySystemGlobalSubsystem&>(s);
}
/*static*/ MultibodySystemGlobalSubsystem&
MultibodySystemGlobalSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MultibodySystemGlobalSubsystem&>(s);
}

const MultibodySystemGlobalSubsystemRep& 
MultibodySystemGlobalSubsystem::getRep() const {
    return dynamic_cast<const MultibodySystemGlobalSubsystemRep&>(*rep);
}
MultibodySystemGlobalSubsystemRep&       
MultibodySystemGlobalSubsystem::updRep() {
    return dynamic_cast<MultibodySystemGlobalSubsystemRep&>(*rep);
}

    //////////////////////////////
    // MolecularMechanicsSystem //
    //////////////////////////////

class DuMMForceFieldSubsystem;

/*static*/ bool 
MolecularMechanicsSystem::isInstanceOf(const System& s) {
    return MolecularMechanicsSystemRep::isA(s.getRep());
}
/*static*/ const MolecularMechanicsSystem&
MolecularMechanicsSystem::downcast(const System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MolecularMechanicsSystem&>(s);
}
/*static*/ MolecularMechanicsSystem&
MolecularMechanicsSystem::updDowncast(System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MolecularMechanicsSystem&>(s);
}

const MolecularMechanicsSystemRep& 
MolecularMechanicsSystem::getRep() const {
    return dynamic_cast<const MolecularMechanicsSystemRep&>(*rep);
}
MolecularMechanicsSystemRep&       
MolecularMechanicsSystem::updRep() {
    return dynamic_cast<MolecularMechanicsSystemRep&>(*rep);
}

MolecularMechanicsSystem::MolecularMechanicsSystem() 
  : MultibodySystem(new MolecularMechanicsSystemRep())
{
}

MolecularMechanicsSystem::MolecularMechanicsSystem
   (MatterSubsystem& matter, DuMMForceFieldSubsystem& mm)
  : MultibodySystem(new MolecularMechanicsSystemRep())
{
    setMatterSubsystem(matter);
    setMolecularMechanicsForceSubsystem(mm);
}

int MolecularMechanicsSystem::setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem& mm) {
    return MolecularMechanicsSystemRep::downcast(*rep).setMolecularMechanicsForceSubsystem(mm);
}

const DuMMForceFieldSubsystem&       
MolecularMechanicsSystem::getMolecularMechanicsForceSubsystem() const {
    return MolecularMechanicsSystemRep::downcast(*rep).getMolecularMechanicsForceSubsystem();
}

DuMMForceFieldSubsystem&       
MolecularMechanicsSystem::updMolecularMechanicsForceSubsystem() {
    return MolecularMechanicsSystemRep::downcast(*rep).updMolecularMechanicsForceSubsystem();
}

} // namespace SimTK

