/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-8 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**@file
 *
 * Implementation of MultibodySystem, a concrete System.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

#include "MultibodySystemRep.h"
#include "DecorationSubsystemRep.h"

namespace SimTK {


    //////////////////////
    // MULTIBODY SYSTEM //
    //////////////////////

/*static*/ bool 
MultibodySystem::isInstanceOf(const System& s) {
    return MultibodySystemRep::isA(s.getSystemGuts());
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
    return dynamic_cast<const MultibodySystemRep&>(getSystemGuts());
}
MultibodySystemRep&       
MultibodySystem::updRep() {
    return dynamic_cast<MultibodySystemRep&>(updSystemGuts());
}

// Create generic multibody system by default.
MultibodySystem::MultibodySystem() {
    adoptSystemGuts(new MultibodySystemRep());
    DefaultSystemSubsystem defsub(*this); // This invokes adoptSubsystem().
    updRep().setGlobalSubsystem();
}

MultibodySystem::MultibodySystem(SimbodyMatterSubsystem& m)
{
    adoptSystemGuts(new MultibodySystemRep());
    DefaultSystemSubsystem defsub(*this); // This invokes adoptSubsystem().
    updRep().setGlobalSubsystem();
    setMatterSubsystem(m);
}

// This is a protected constructor for use by derived classes which
// allocate a more specialized MultibodySystemRep.
MultibodySystem::MultibodySystem(MultibodySystemRep* rp) {
    adoptSystemGuts(rp);
    DefaultSystemSubsystem defsub(*this); // This invokes adoptSubsystem().
    updRep().setGlobalSubsystem();
}

int MultibodySystem::setMatterSubsystem(SimbodyMatterSubsystem& m) {
    return updRep().setMatterSubsystem(m);
}
int MultibodySystem::addForceSubsystem(ForceSubsystem& f) {
    return updRep().addForceSubsystem(f);
}
int MultibodySystem::setDecorationSubsystem(DecorationSubsystem& m) {
    return updRep().setDecorationSubsystem(m);
}
int MultibodySystem::setContactSubsystem(GeneralContactSubsystem& m) {
    return updRep().setContactSubsystem(m);
}

const SimbodyMatterSubsystem&       
MultibodySystem::getMatterSubsystem() const {
    return getRep().getMatterSubsystem();
}
SimbodyMatterSubsystem&       
MultibodySystem::updMatterSubsystem() {
    return updRep().updMatterSubsystem();
}
bool MultibodySystem::hasMatterSubsystem() const {
    return getRep().hasMatterSubsystem();
}

const DecorationSubsystem&       
MultibodySystem::getDecorationSubsystem() const {
    return getRep().getDecorationSubsystem();
}
DecorationSubsystem&       
MultibodySystem::updDecorationSubsystem() {
    return updRep().updDecorationSubsystem();
}
bool MultibodySystem::hasDecorationSubsystem() const {
    return getRep().hasDecorationSubsystem();
}

const GeneralContactSubsystem&       
MultibodySystem::getContactSubsystem() const {
    return getRep().getContactSubsystem();
}
GeneralContactSubsystem&       
MultibodySystem::updContactSubsystem() {
    return updRep().updContactSubsystem();
}
bool MultibodySystem::hasContactSubsystem() const {
    return getRep().hasContactSubsystem();
}

const Real
MultibodySystem::calcPotentialEnergy(const State& s) const {
    return getRep().calcPotentialEnergy(s);
}
const Real
MultibodySystem::calcKineticEnergy(const State& s) const {
    return getMatterSubsystem().getRep().calcKineticEnergy(s);
}

const Vector_<SpatialVec>& 
MultibodySystem::getRigidBodyForces(const State& s, Stage g) const {
    return getRep().getRigidBodyForces(s,g);
}
const Vector_<Vec3>&       
MultibodySystem::getParticleForces(const State& s, Stage g) const {
    return getRep().getParticleForces(s,g);
}
const Vector&              
MultibodySystem::getMobilityForces(const State& s, Stage g) const {
    return getRep().getMobilityForces(s,g);
}

Vector_<SpatialVec>& 
MultibodySystem::updRigidBodyForces(const State& s, Stage g) const {
    return getRep().updRigidBodyForces(s,g);
}
Vector_<Vec3>&       
MultibodySystem::updParticleForces(const State& s, Stage g) const {
    return getRep().updParticleForces(s,g);
}
Vector&              
MultibodySystem::updMobilityForces(const State& s, Stage g) const {
    return getRep().updMobilityForces(s,g);
}


    //////////////////////////
    // MULTIBODY SYSTEM REP //
    //////////////////////////

int MultibodySystemRep::realizeTopologyImpl(State& s) const {
    assert(globalSub.isValid());
    assert(matterSub.isValid());

    // We do Matter subsystem first here in case any of the GlobalSubsystem
    // topology depends on Matter topology. That's unlikely though since
    // we don't know sizes until Model stage.
    getMatterSubsystem().getRep().realizeSubsystemTopology(s);
    getGlobalSubsystem().getRep().realizeSubsystemTopology(s);
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemTopology(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemTopology(s);

    return 0;
}
int MultibodySystemRep::realizeModelImpl(State& s) const {

    // Here it is essential to do the Matter subsystem first because the
    // force accumulation arrays in the Global subsystem depend on the
    // Stage::Model dimensions of the Matter subsystem.
    getMatterSubsystem().getRep().realizeSubsystemModel(s);
    getGlobalSubsystem().getRep().realizeSubsystemModel(s);
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemModel(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemModel(s);

    return 0;
}
int MultibodySystemRep::realizeInstanceImpl(const State& s) const {
    getGlobalSubsystem().getRep().realizeSubsystemInstance(s);
    getMatterSubsystem().getRep().realizeSubsystemInstance(s);
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemInstance(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemInstance(s);

    return 0;
}
int MultibodySystemRep::realizeTimeImpl(const State& s) const {
    getGlobalSubsystem().getRep().realizeSubsystemTime(s);
    getMatterSubsystem().getRep().realizeSubsystemTime(s);
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemTime(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemTime(s);

    return 0;
}
int MultibodySystemRep::realizePositionImpl(const State& s) const {
    getGlobalSubsystem().getRep().realizeSubsystemPosition(s);
    getMatterSubsystem().getRep().realizeSubsystemPosition(s);
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemPosition(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemPosition(s);

    return 0;
}
int MultibodySystemRep::realizeVelocityImpl(const State& s) const {
    getGlobalSubsystem().getRep().realizeSubsystemVelocity(s);
    getMatterSubsystem().getRep().realizeSubsystemVelocity(s);
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemVelocity(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemVelocity(s);

    return 0;
}
int MultibodySystemRep::realizeDynamicsImpl(const State& s) const {
    getGlobalSubsystem().getRep().realizeSubsystemDynamics(s);
    if (hasContactSubsystem())
        getContactSubsystem().getSubsystemGuts().realizeSubsystemDynamics(s);
    // note order: forces first (TODO: does that matter?)
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemDynamics(s);
    getMatterSubsystem().getRep().realizeSubsystemDynamics(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemDynamics(s);

    return 0;
}
int MultibodySystemRep::realizeAccelerationImpl(const State& s) const {
    getGlobalSubsystem().getRep().realizeSubsystemAcceleration(s);
    // note order: forces first (TODO: does that matter?)
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemAcceleration(s);
    getMatterSubsystem().getRep().realizeSubsystemAcceleration(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemAcceleration(s);

    return 0;
}
int MultibodySystemRep::realizeReportImpl(const State& s) const {
    getGlobalSubsystem().getRep().realizeSubsystemReport(s);
    // note order: forces first (TODO: does that matter?)
    for (int i=0; i < (int)forceSubs.size(); ++i)
        getForceSubsystem(forceSubs[i]).getRep().realizeSubsystemReport(s);
    getMatterSubsystem().getRep().realizeSubsystemReport(s);

    if (hasDecorationSubsystem())
        getDecorationSubsystem().getGuts().realizeSubsystemReport(s);

    return 0;
}


    ///////////////////////////////////////
    // MULTIBODY SYSTEM GLOBAL SUBSYSTEM //
    ///////////////////////////////////////


/*static*/ bool 
MultibodySystemGlobalSubsystem::isInstanceOf(const Subsystem& s) {
    return MultibodySystemGlobalSubsystemRep::isA(s.getSubsystemGuts());
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
    return dynamic_cast<const MultibodySystemGlobalSubsystemRep&>(getSubsystemGuts());
}
MultibodySystemGlobalSubsystemRep&       
MultibodySystemGlobalSubsystem::updRep() {
    return dynamic_cast<MultibodySystemGlobalSubsystemRep&>(updSubsystemGuts());
}


} // namespace SimTK

