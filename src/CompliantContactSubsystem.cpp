/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-10 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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
Private implementation of CompliantContactSubsystem and the built in contact
patch analysis methods.
**/

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/ForceSubsystemGuts.h"
#include "simbody/internal/CompliantContactSubsystem.h"
#include "simbody/internal/ContactTrackerSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/MultibodySystem.h"

namespace SimTK {

//==============================================================================
//                    COMPLIANT CONTACT SUBSYSTEM IMPL
//==============================================================================

class CompliantContactSubsystemImpl : public ForceSubsystemRep {
typedef std::map<ContactTypeId, const ContactForceGenerator*> GeneratorMap;
public:

CompliantContactSubsystemImpl(const ContactTrackerSubsystem& tracker)
:   ForceSubsystemRep("CompliantContactSubsystem", "0.0.1"),
    m_tracker(tracker), m_transitionVelocity(0.01), 
    m_ooTransitionVelocity(1/m_transitionVelocity), 
    m_defaultGenerator(0) 
{   
}

Real getTransitionVelocity() const  {return m_transitionVelocity;}
Real getOOTransitionVelocity() const  {return m_ooTransitionVelocity;}
void setTransitionVelocity(Real vt) 
{   m_transitionVelocity=vt; m_ooTransitionVelocity=1/vt;}


int getNumContactForces(const State& s) const {
    ensureForceCacheValid(s);
    const Array_<ContactForce>& forces = getForceCache(s);
    return (int)forces.size();
}

const ContactForce& getContactForce(const State& s, int n) const {
    const int numContactForces = getNumContactForces(s); // realize if needed
    SimTK_ERRCHK2_ALWAYS(n < numContactForces,
        "CompliantContactSubsystem::getContactForce()",
        "There are currently only %d contact forces but force %d"
        " was requested (they are numbered from 0). Use getNumContactForces()"
        " first to determine how many are available.", numContactForces, n);
    const Array_<ContactForce>& forces = getForceCache(s);
    return forces[n];
}


const ContactForce& getContactForceById(const State& s, ContactId id) const {
    static const ContactForce invalidForce;
    const int numContactForces = getNumContactForces(s); // realize if needed
    const Array_<ContactForce>& forces = getForceCache(s);
    // TODO: use a map
    for (int i=0; i < numContactForces; ++i)
        if (forces[i].getContactId()==id)
            return forces[i];

    return invalidForce;
}

bool calcContactPatchDetailsById(const State&   state,
                                 ContactId      id,
                                 ContactPatch&  patch_G) const
{
    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Velocity,
        "CompliantContactSubystemImpl::ensureForceCacheValid()");

    const ContactSnapshot& active = m_tracker.getActiveContacts(state);
    const Contact& contact = active.getContactById(id);
    if (contact.isEmpty()) {
        patch_G.clear();
        return false;
    }

    const Transform& X_S1S2 = contact.getTransform();
    const ContactSurfaceIndex surf1(contact.getSurface1());
    const ContactSurfaceIndex surf2(contact.getSurface2());
    const MobilizedBody& mobod1 = m_tracker.getMobilizedBody(surf1);
    const MobilizedBody& mobod2 = m_tracker.getMobilizedBody(surf2);

    const Transform X_GS1 = mobod1.findFrameTransformInGround
        (state, m_tracker.getContactSurfaceTransform(surf1));
    const Transform X_GS2 = mobod2.findFrameTransformInGround
        (state, m_tracker.getContactSurfaceTransform(surf2));

    const SpatialVec V_GS1 = mobod1.findFrameVelocityInGround
        (state, m_tracker.getContactSurfaceTransform(surf1));
    const SpatialVec V_GS2 = mobod2.findFrameVelocityInGround
        (state, m_tracker.getContactSurfaceTransform(surf2));

    // Calculate the relative velocity of S2 in S1, expressed in S1.
    const SpatialVec V_S1S2 =
        findRelativeVelocity(X_GS1, V_GS1, X_GS2, V_GS2);

    // Get the right force generator to use for this kind of contact.
    const ContactForceGenerator& generator = 
        getForceGenerator(contact.getTypeId());

    // Calculate the contact patch measured and expressed in S1.
    generator.calcContactPatch(state, contact, V_S1S2, 
        patch_G); // it is actually in S1 at this point

    // Re-express the contact patch in Ground for later use.
    if (patch_G.isValid()) {
        patch_G.changeFrameInPlace(X_GS1); // switch from S1 to Ground
        return true;
    }

    return false;
}

const ContactTrackerSubsystem& getContactTrackerSubsystem() const
{   return m_tracker; }

~CompliantContactSubsystemImpl() {
    delete m_defaultGenerator;
    for (GeneratorMap::iterator p  = m_generators.begin(); 
                                p != m_generators.end(); ++p)
        delete p->second; // free the generator
    // The map itself gets freed by its destructor here.
}

// We're going to take over ownership of this object. If the generator map
// already contains a generator for this type of contact, the new one replaces
// it.
void adoptForceGenerator(CompliantContactSubsystem* subsys, 
                         ContactForceGenerator*     generator) {
    assert(generator);
    generator->setCompliantContactSubsystem(subsys);
    invalidateSubsystemTopologyCache();
    GeneratorMap::iterator p=m_generators.find(generator->getContactTypeId());
    if (p != m_generators.end()) {
        // We're replacing an existing generator.
        delete p->second;
        p->second = generator;  // just copying the pointer
    } else {
        // Insert the new generator.
        m_generators.insert(std::make_pair(generator->getContactTypeId(), 
                                           generator));
    }
}

// Set the default force generator, deleting the previous one. Null is OK.
void adoptDefaultForceGenerator(CompliantContactSubsystem* subsys, 
                                ContactForceGenerator* generator) {
    invalidateSubsystemTopologyCache();
    delete m_defaultGenerator;
    m_defaultGenerator = generator; // just copying the pointer
    if (m_defaultGenerator) 
        m_defaultGenerator->setCompliantContactSubsystem(subsys);
}

// Do we have a generator registered for this type? If not, we'll invoke our
// default generator.
bool hasForceGenerator(ContactTypeId type) const 
{   return m_generators.find(type) != m_generators.end(); }

bool hasDefaultForceGenerator() const 
{   return m_defaultGenerator != 0; }

// Get the generator registered for this type of contact or the default
// generator if nothing is registered. 
const ContactForceGenerator& getForceGenerator(ContactTypeId type) const 
{   GeneratorMap::const_iterator p=m_generators.find(type);
    return p != m_generators.end() ? *p->second : getDefaultForceGenerator(); }

// Get the default generator.
const ContactForceGenerator& getDefaultForceGenerator() const
{   assert(m_defaultGenerator); return *m_defaultGenerator; }


// These override default implementations of virtual methods in the 
// Subsystem::Guts class.

CompliantContactSubsystemImpl* cloneImpl() const 
{   return new CompliantContactSubsystemImpl(*this); }

// Allocate lazy evaluation cache entries to hold potential energy (calculable
// after position stage) and forces and power (calculable after velocity stage).
int realizeSubsystemTopologyImpl(State& s) const {
    // Get writability briefly to fill in the Topology cache.
    CompliantContactSubsystemImpl* wThis = 
        const_cast<CompliantContactSubsystemImpl*>(this);

    // Calculating forces includes calculating PE for each force.
    wThis->m_forceCacheIx = allocateLazyCacheEntry(s, 
        Stage::Velocity, new Value<Array_<ContactForce> >());
    // This usually just requires summing up the forceCache PE values, but
    // may also have to be calculated at Position stage in which case we can't
    // calculate the forces.
    wThis->m_potEnergyCacheIx = allocateLazyCacheEntry(s, 
        Stage::Position, new Value<Real>(NaN));

    // This state variable is used to integrate power to get dissipated
    // energy.
    Vector einit(1, Real(0));
    wThis->m_dissipatedEnergyIx = allocateZ(s,einit);

    return 0;
}

int realizeSubsystemModelImpl(State& s) const {
    return 0;
}

int realizeSubsystemInstanceImpl(const State& s) const {
    return 0;
}

int realizeSubsystemTimeImpl(const State& s) const {
    return 0;
}

int realizeSubsystemPositionImpl(const State& s) const {
    return 0;
}

int realizeSubsystemVelocityImpl(const State& s) const {
    return 0;
}

int realizeSubsystemDynamicsImpl(const State& s) const {
    ensureForceCacheValid(s);

    const MultibodySystem&        mbs    = getMultibodySystem(); // my owner
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // Get access to System-global force cache array.
    Vector_<SpatialVec>& rigidBodyForces =
        mbs.updRigidBodyForces(s, Stage::Dynamics);

    // Accumulate the values from the cache into the global arrays.
    const ContactSnapshot& contacts = m_tracker.getActiveContacts(s);
    const Array_<ContactForce>& forces = getForceCache(s);
    for (unsigned i=0; i < forces.size(); ++i) {
        const ContactForce& force = forces[i];
        const Contact& contact = contacts.getContactById(force.getContactId());
        const MobilizedBody& mobod1 = m_tracker.getMobilizedBody
                                                (contact.getSurface1());
        const MobilizedBody& mobod2 = m_tracker.getMobilizedBody
                                                (contact.getSurface2());
        const Vec3 r1 = force.getContactPoint() - mobod1.getBodyOriginLocation(s);
        const Vec3 r2 = force.getContactPoint() - mobod2.getBodyOriginLocation(s);
        const SpatialVec& F2cpt = force.getForceOnSurface2(); // at contact pt
        // Shift applied force to body origins.
        const SpatialVec F2( F2cpt[0] + r2 %  F2cpt[1],  F2cpt[1]);
        const SpatialVec F1(-F2cpt[0] + r1 % -F2cpt[1], -F2cpt[1]);
        mobod1.applyBodyForce(s, F1, rigidBodyForces);
        mobod2.applyBodyForce(s, F2, rigidBodyForces);
    }

    return 0;
}

// Potential energy is normally a side effect of force calculation done after
// Velocity stage. But if only positions are available, we
// have to calculate forces at zero velocity and then throw away everything
// but the PE.
Real calcPotentialEnergy(const State& state) const {
    ensurePotentialEnergyCacheValid(state);
    return getPotentialEnergyCache(state);
}

int realizeSubsystemAccelerationImpl(const State& state) const {
    ensureForceCacheValid(state);
    Real powerLoss = 0;
    const Array_<ContactForce>& forces = getForceCache(state);
    for (unsigned i=0; i < forces.size(); ++i)
        powerLoss += forces[i].getPowerDissipation();
    updDissipatedEnergyDeriv(state) = powerLoss;
    return 0;
}

int realizeSubsystemReportImpl(const State&) const {
    return 0;
}

const Real& getDissipatedEnergyVar(const State& s) const
{   return getZ(s)[m_dissipatedEnergyIx]; }
Real& updDissipatedEnergyVar(State& s) const
{   return updZ(s)[m_dissipatedEnergyIx]; }

//--------------------------------------------------------------------------
                                  private:

Real& updDissipatedEnergyDeriv(const State& s) const
{   return updZDot(s)[m_dissipatedEnergyIx]; }

const Real& getPotentialEnergyCache(const State& s) const
{   return Value<Real>::downcast(getCacheEntry(s,m_potEnergyCacheIx)); }
const Array_<ContactForce>& getForceCache(const State& s) const
{   return Value<Array_<ContactForce> >::downcast
                                    (getCacheEntry(s,m_forceCacheIx)); }

Real& updPotentialEnergyCache(const State& s) const
{   return Value<Real>::updDowncast(updCacheEntry(s,m_potEnergyCacheIx)); }
Array_<ContactForce>& updForceCache(const State& s) const
{   return Value<Array_<ContactForce> >::updDowncast
                                    (updCacheEntry(s,m_forceCacheIx)); }

bool isPotentialEnergyCacheValid(const State& s) const
{   return isCacheValueRealized(s,m_potEnergyCacheIx); }
bool isForceCacheValid(const State& s) const
{   return isCacheValueRealized(s,m_forceCacheIx); }


void markPotentialEnergyCacheValid(const State& s) const
{   markCacheValueRealized(s,m_potEnergyCacheIx); }
void markForceCacheValid(const State& s) const
{   markCacheValueRealized(s,m_forceCacheIx); }

void ensurePotentialEnergyCacheValid(const State&) const;
void ensureForceCacheValid(const State&) const;



    // TOPOLOGY "STATE"

// This is the ContactTracker that we will consult to obtain Contact objects
// to use as raw material for contact force generation.
const ContactTrackerSubsystem&      m_tracker;

// This is the value that should be used by generators for the friction
// transition velocity if they are using a Hollars-like friction model.
// Also precalculate the inverse to avoid expensive runtime division.
Real                                m_transitionVelocity;
Real                                m_ooTransitionVelocity; // 1/vt

// This map owns the generator objects; be sure to clean up on destruction or
// when a generator is replaced.
GeneratorMap                        m_generators;

// This is the generator to use for an unrecognized ContactTypeId. Typically
// this will either do nothing silently or throw an error.
ContactForceGenerator*              m_defaultGenerator;

    // TOPOLOGY "CACHE"

// These must be set during realizeTopology and treated as const thereafter.
ZIndex                              m_dissipatedEnergyIx;
CacheEntryIndex                     m_potEnergyCacheIx;
CacheEntryIndex                     m_forceCacheIx;
};

void CompliantContactSubsystemImpl::
ensurePotentialEnergyCacheValid(const State& state) const {
    if (isPotentialEnergyCacheValid(state)) return;

    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Position,
        "CompliantContactSubystemImpl::ensurePotentialEnergyCacheValid()");

    Real& pe = updPotentialEnergyCache(state);
    pe = 0;

    // If the State has been realized to Velocity stage, we can just sum up
    // the PEs from the ContactForce objects.
    if (getStage(state) >= Stage::Velocity) {
        ensureForceCacheValid(state);
        const Array_<ContactForce>& forces = getForceCache(state);
        for (unsigned i=0; i < forces.size(); ++i)
            pe += forces[i].getPotentialEnergy();
        markPotentialEnergyCacheValid(state);
        return;
    }

    // The State has been realized to Position stage, so we're going to have
    // to calculate forces at zero velocity and then throw away all the 
    // results except for the PE.
    const ContactSnapshot& active = m_tracker.getActiveContacts(state);
    const int nContacts = active.getNumContacts();
    for (int i=0; i<nContacts; ++i) {
        const Contact& contact = active.getContact(i);
        const ContactForceGenerator& generator = 
            getForceGenerator(contact.getTypeId());
        ContactForce force;
        generator.calcContactForce(state,contact,SpatialVec(Vec3(0)), force);
        pe += force.getPotentialEnergy();
    }

    markPotentialEnergyCacheValid(state);
}


void CompliantContactSubsystemImpl::
ensureForceCacheValid(const State& state) const {
    if (isForceCacheValid(state)) return;

    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Velocity,
        "CompliantContactSubystemImpl::ensureForceCacheValid()");

    Array_<ContactForce>& forces = updForceCache(state);
    forces.clear();


    const ContactSnapshot& active = m_tracker.getActiveContacts(state);
    const int nContacts = active.getNumContacts();
    for (int i=0; i<nContacts; ++i) {
        const Contact& contact = active.getContact(i);
        const Transform& X_S1S2 = contact.getTransform();
        const ContactSurfaceIndex surf1(contact.getSurface1());
        const ContactSurfaceIndex surf2(contact.getSurface2());
        const MobilizedBody& mobod1 = m_tracker.getMobilizedBody(surf1);
        const MobilizedBody& mobod2 = m_tracker.getMobilizedBody(surf2);

        // TODO: These two are expensive (63 flops each) and shouldn't have 
        // to be recalculated here since we must have used them in creating
        // the Contact and X_S1S2.
        const Transform X_GS1 = mobod1.findFrameTransformInGround
            (state, m_tracker.getContactSurfaceTransform(surf1));
        const Transform X_GS2 = mobod2.findFrameTransformInGround
            (state, m_tracker.getContactSurfaceTransform(surf2));

        const SpatialVec V_GS1 = mobod1.findFrameVelocityInGround
            (state, m_tracker.getContactSurfaceTransform(surf1));
        const SpatialVec V_GS2 = mobod2.findFrameVelocityInGround
            (state, m_tracker.getContactSurfaceTransform(surf2));

        // Calculate the relative velocity of S2 in S1, expressed in S1.
        const SpatialVec V_S1S2 =
            findRelativeVelocity(X_GS1, V_GS1, X_GS2, V_GS2);   // 51 flops

        const ContactForceGenerator& generator = 
            getForceGenerator(contact.getTypeId());
        forces.push_back(); // allocate a new garbage ContactForce
        // Calculate the contact force measured and expressed in S1.
        generator.calcContactForce(state, contact, V_S1S2, forces.back());
        // Re-express the contact force in Ground for later use.
        if (forces.back().isValid())
            forces.back().changeFrameInPlace(X_GS1); // switch to Ground
        else
            forces.pop_back(); // never mind ...
    }

    markForceCacheValid(state);
}


//==============================================================================
//                      COMPLIANT CONTACT SUBSYSTEM
//==============================================================================

/*static*/ bool 
CompliantContactSubsystem::isInstanceOf(const ForceSubsystem& s) 
{   return CompliantContactSubsystemImpl::isA(s.getRep()); }
/*static*/ const CompliantContactSubsystem&
CompliantContactSubsystem::downcast(const ForceSubsystem& s) 
{   assert(isInstanceOf(s));
    return reinterpret_cast<const CompliantContactSubsystem&>(s); }
/*static*/ CompliantContactSubsystem&
CompliantContactSubsystem::updDowncast(ForceSubsystem& s) 
{   assert(isInstanceOf(s));
    return reinterpret_cast<CompliantContactSubsystem&>(s); }

const CompliantContactSubsystemImpl& 
CompliantContactSubsystem::getImpl() const 
{   return dynamic_cast<const CompliantContactSubsystemImpl&>
                                            (ForceSubsystem::getRep()); }
CompliantContactSubsystemImpl&       
CompliantContactSubsystem::updImpl() 
{   return dynamic_cast<CompliantContactSubsystemImpl&>
                                            (ForceSubsystem::updRep()); }

CompliantContactSubsystem::CompliantContactSubsystem
   (MultibodySystem& mbs, const ContactTrackerSubsystem& tracker)
:   ForceSubsystem() 
{   adoptSubsystemGuts(new CompliantContactSubsystemImpl(tracker));
    mbs.addForceSubsystem(*this);  // mbs steals ownership

    // Register default contact force generators.
    adoptForceGenerator(new ContactForceGenerator::HertzCircular());
    adoptForceGenerator(new ContactForceGenerator::ElasticFoundation());
    adoptDefaultForceGenerator(new ContactForceGenerator::DoNothing());
}

Real CompliantContactSubsystem::getTransitionVelocity() const
{   return getImpl().getTransitionVelocity(); }
Real CompliantContactSubsystem::getOOTransitionVelocity() const
{   return getImpl().getOOTransitionVelocity(); }
void CompliantContactSubsystem::setTransitionVelocity(Real vt) {
    SimTK_ERRCHK1_ALWAYS(vt > 0, 
        "CompliantContactSubsystem::setTransitionVelocity()",
        "Friction transition velocity %g is illegal.", vt);
    updImpl().setTransitionVelocity(vt);
}

int CompliantContactSubsystem::getNumContactForces(const State& s) const
{   return getImpl().getNumContactForces(s); }

const ContactForce& CompliantContactSubsystem::
getContactForce(const State& s, int n) const
{   return getImpl().getContactForce(s,n); }


const ContactForce& CompliantContactSubsystem::
getContactForceById(const State& s, ContactId id) const
{   return getImpl().getContactForceById(s,id); }

bool CompliantContactSubsystem::
calcContactPatchDetailsById(const State&   state,
                            ContactId      id,
                            ContactPatch&  patch_G) const
{   return getImpl().calcContactPatchDetailsById(state,id,patch_G); }

Real CompliantContactSubsystem::
getDissipatedEnergy(const State& s) const
{   return getImpl().getDissipatedEnergyVar(s); }

void CompliantContactSubsystem::
setDissipatedEnergy(State& s, Real energy) const {
    SimTK_ERRCHK1_ALWAYS(energy >= 0,
        "CompliantContactSubsystem::setDissipatedEnergy()",
        "The initial value for the dissipated energy must be nonnegative"
        " but an attempt was made to set it to %g.", energy);
    getImpl().updDissipatedEnergyVar(s) = energy; 
}


void CompliantContactSubsystem::adoptForceGenerator
   (ContactForceGenerator* generator)
{   return updImpl().adoptForceGenerator(this, generator); }

void CompliantContactSubsystem::adoptDefaultForceGenerator
   (ContactForceGenerator* generator)
{   return updImpl().adoptDefaultForceGenerator(this, generator); }

bool CompliantContactSubsystem::hasForceGenerator(ContactTypeId type) const
{   return getImpl().hasForceGenerator(type); }
bool CompliantContactSubsystem::hasDefaultForceGenerator() const
{   return getImpl().hasDefaultForceGenerator(); }

const ContactForceGenerator& CompliantContactSubsystem::
getContactForceGenerator(ContactTypeId contactType) const
{   return getImpl().getForceGenerator(contactType); }

const ContactForceGenerator& CompliantContactSubsystem::
getDefaultForceGenerator() const
{   return getImpl().getDefaultForceGenerator(); }

const ContactTrackerSubsystem& CompliantContactSubsystem::
getContactTrackerSubsystem() const
{   return getImpl().getContactTrackerSubsystem(); }

const MultibodySystem& CompliantContactSubsystem::getMultibodySystem() const
{   return MultibodySystem::downcast(getSystem()); }


// Input x goes from 0 to 1; output goes 0 to 1 but smoothed with an S-shaped 
// quintic with two zero derivatives at either end. Cost is 7 flops.
inline static Real step5(Real x) {
    assert(0 <= x && x <= 1);
    const Real x3=x*x*x;
    return x3*(10+x*(6*x-15)); //10x^3-15x^4+6x^5
}

// Input x goes from 0 to 1; output goes from 0 to y. First derivative
// is 0 at the beginning and yd at the end. Second derivatives are zero
// at both ends:
//
//    y -               *
//                   *  | slope=yd
//                *------
//             *
//    0 -  * * 
//       x=0            1
//    
// Cost is 16 flops.
inline static Real step5d(Real y, Real yd, Real x) {
    assert(0 <= x && x <= 1);
    const Real a=6*y-3*yd, b=-15*y+7*yd, c=10*y-4*yd, x3=x*x*x;
    return x3*(c + x*(b + x*a));
}


// This is the sum of two curves:
// (1) a wet friction term mu_wet which is a linear function of velocity: 
//     mu_wet = uv*v
// (2) a dry friction term mu_dry which is a quintic spline with 4 segments:
//     mu_dry = 
//      (a) v=0..1: smooth interpolation from 0 to us
//      (b) v=1..3: smooth interp from us down to ud (Stribeck)
//      (c) v=3..Inf ud
// CAUTION: uv and v must be dimensionless in multiples of transition velocity.
// The function mu = mu_wet + mu_dry is zero at v=0 with 1st deriv (slope) uv
// and 2nd deriv (curvature) 0. At large velocities v>>0 the value is 
// ud+uv*v, again with slope uv and zero curvature. We want mu(v) to be c2
// continuous, so mu_wet(v) must have zero slope and curvature at v==0 and
// at v==3 where it takes on a constant value ud.
//
// Cost: stiction 12 flops
//       stribeck 14 flops
//       sliding 3 flops
// Curve looks like this:
//
//  us+uv     ***
//           *    *                     *
//           *     *               *____| slope = uv at Inf
//           *      *         *
// ud+3uv    *        *  *      
//          *          
//          *        
//         *
//  0  *____| slope = uv at 0
//
//     |    |           |
//   v=0    1           3 
//
// This calculates a composite coefficient of friction that you should use
// to scale the normal force to produce the friction force.
inline static Real stribeck(Real us, Real ud, Real uv, Real v) {
    const Real mu_wet = uv*v;
    Real mu_dry;
    if      (v >= 3) mu_dry = ud; // sliding
    else if (v >= 1) mu_dry = us - (us-ud)*step5((v-1)/2); // Stribeck
    else             mu_dry = us*step5(v); // 0 <= v < 1 (stiction)
    return mu_dry + mu_wet;
}

// CAUTION: uv and v must be dimensionless in multiples of transition velocity.
// Const 9 flops + 1 divide = approx 25 flops.
// This calculates a composite coefficient of friction that you should use
// to scale the normal force to produce the friction force.
// Curve is similar to Stribeck above but more violent with a discontinuous
// derivative at v==1.
inline static Real hollars(Real us, Real ud, Real uv, Real v) {
    const Real mu =   std::min(v, 1.0)*(ud + 2*(us-ud)/(1+v*v))
                    + uv*v;
    return mu;
}



//==============================================================================
//                         HERTZ CIRCULAR GENERATOR
//==============================================================================
void ContactForceGenerator::HertzCircular::calcContactForce
   (const State&            state,
    const Contact&          overlap,    // contains X_S1S2
    const SpatialVec&       V_S1S2,     // relative surface velocity, S2 in S1
    ContactForce&           contactForce_S1) const
{
    SimTK_ASSERT(CircularPointContact::isInstance(overlap),
        "ContactForceGenerator::HertzCircular::calcContactForce(): expected"
        " CircularPointContact.");

    const CircularPointContact& contact = CircularPointContact::getAs(overlap);
    const Real depth = contact.getDepth();

    if (depth <= 0) {
        contactForce_S1.clear(); // no contact; invalidate return result
        return;
    }

    const CompliantContactSubsystem& subsys = getCompliantContactSubsystem();
    const ContactTrackerSubsystem&   tracker = subsys.getContactTrackerSubsystem();
    
    const ContactSurfaceIndex surf1x = contact.getSurface1();
    const ContactSurfaceIndex surf2x = contact.getSurface2();
    const Transform&          X_S1S2 = contact.getTransform();

    // Abbreviations.
    const Rotation& R12 = X_S1S2.R(); // orientation of S2 in S1
    const Vec3&     p12 = X_S1S2.p(); // position of O2 in S1
    const Vec3&     w12 = V_S1S2[0];  // ang. vel. of S2 in S1
    const Vec3&     v12 = V_S1S2[1];  // vel. of O2 in S1

    const ContactSurface&  surf1  = tracker.getContactSurface(surf1x);
    const ContactSurface&  surf2  = tracker.getContactSurface(surf2x);

    const ContactMaterial& mat1   = surf1.getMaterial();
    const ContactMaterial& mat2   = surf2.getMaterial();

    // Use 2/3 power of stiffness for combining here.
    const Real k1=mat1.getStiffness23(), k2=mat2.getStiffness23();
    const Real c1=mat1.getDissipation(), c2=mat2.getDissipation();

    // Adjust the contact location based on the relative stiffness of the 
    // two materials. s1 is the fraction of the "squishing" (deformation) that
    // is done by surface1 -- if surface2 is very stiff then surface1 does most.
    
    const Real s1 = k2/(k1+k2); // 0..1
    const Real s2 = 1-s1;       // 1..0
    const Real x  = contact.getDepth();
    // normal points away from surf1, expressed in surf1's frame 
    const UnitVec3& normal_S1 = contact.getNormal();
    // origin is half way between the two surfaces as though they had
    // equal stiffness
    const Vec3& origin_S1 = contact.getOrigin();
    // Actual contact point moves closer to stiffer surface.
    const Vec3 contactPt_S1 = origin_S1 + (x*(0.5-s1))*normal_S1;
    
    // Calculate the Hertz force fH, which is conservative.
    const Real k = k1*s1; // (==k2*s2) == E^(2/3)
    const Real c = c1*s1 + c2*s2;
    const Real R = contact.getEffectiveRadius();
    const Real fH = (4./3.)*k*x*std::sqrt(R*k*x); // always >= 0
    
    // Calculate the relative velocity of the two bodies at the contact point.
    // We're considering S1 fixed, so we just need the velocity in S1 of
    // the station of S2 that is coincident with the contact point.
    const Vec3 contactPt2_S1 = contactPt_S1 - p12;  // S2 station, exp. in S1

    // All vectors are in S1; dropping the "_S1" notation now.

    // Velocity of surf2 at contact point is opposite direction of normal
    // when penetration is increasing.
    const Vec3 vel = v12 + w12 % contactPt2_S1;
    // Want xdot > 0 when x is increasing; normal points the other way
    const Real xdot = -dot(vel, normal_S1); // penetration rate (signed)
    const Vec3 velNormal  = -xdot*normal_S1;
    const Vec3 velTangent = vel-velNormal;
    
    // Calculate the Hunt-Crossley force, which is dissipative, and the 
    // total normal force including elasticity and dissipation.
    const Real fHC      = fH*1.5*c*xdot; // same sign as xdot
    const Real fNormal  = fH + fHC;      // < 0 means "sticking"; see below

    // Start filling out the contact force.
    contactForce_S1.setContactId(contact.getContactId());
    contactForce_S1.setContactPoint(contactPt_S1);

    // Total force can be negative under unusual circumstances ("yanking");
    // that means no force is generated and no stored PE will be recovered.
    // This will most often occur in to-be-rejected trial steps but can
    // occasionally be real.
    if (fNormal <= 0) {
        //std::cout << "YANKING!!!\n";
        contactForce_S1.setForceOnSurface2(SpatialVec(Vec3(0)));
        contactForce_S1.setPotentialEnergy(0);
        contactForce_S1.setPowerDissipation(0);
        return; // there is contact, but no force
    }

    const Vec3 forceH          = fH *normal_S1; // as applied to surf2
    const Vec3 forceHC         = fHC*normal_S1;
    const Real potentialEnergy = (2./5.)*fH*x;
    const Real powerHC         = fHC*xdot; // rate of energy loss, >= 0

    // Calculate the friction force.
    Vec3 forceFriction(0);
    Real powerFriction = 0;
    const Real vslipSq = velTangent.normSqr();
    if (vslipSq > square(SignificantReal)) {
        const Real vslip = std::sqrt(vslipSq); // expensive
        const Real us1=mat1.getStaticFriction(), us2=mat2.getStaticFriction();
        const Real ud1=mat1.getDynamicFriction(), ud2=mat2.getDynamicFriction();
        const Real uv1=mat1.getViscousFriction(), uv2=mat2.getViscousFriction();
        const Real vtrans = subsys.getTransitionVelocity();
        const Real ooVtrans = subsys.getOOTransitionVelocity(); // 1/vtrans

        // Calculate effective coefficients of friction, being careful not
        // to divide 0/0 if both are frictionless.
        Real us = 2*us1*us2; if (us!=0) us /= (us1+us2);
        Real ud = 2*ud1*ud2; if (ud!=0) ud /= (ud1+ud2);
        Real uv = 2*uv1*uv2; if (uv!=0) uv /= (uv1+uv2);
        assert(us >= ud);

        // Express slip velocity as unitless multiple of transition velocity.
        const Real v = vslip * ooVtrans;
        // Must scale viscous coefficient to match unitless velocity.
        const Real mu=stribeck(us,ud,uv*vtrans,v);
        //const Real mu=hollars(us,ud,uv*vtrans,v);
        const Real fFriction = fNormal * mu;
        // Force direction on S2 opposes S2's velocity.
        forceFriction = (-fFriction/vslip)*velTangent; // in S1
        powerFriction = fFriction * vslip; // >= 0
    }

    const Vec3 forceLoss  = forceHC + forceFriction;
    const Vec3 forceTotal = forceH + forceLoss;
    
    // Report the force (as applied at contactPt).
    contactForce_S1.setForceOnSurface2(SpatialVec(Vec3(0),forceTotal));
    contactForce_S1.setPotentialEnergy(potentialEnergy);
    // Don't include dot(forceH,velNormal) power due to conservative force
    // here. This way we don't double-count the energy on the way in as
    // integrated power and potential energy. Although the books would balance
    // again when the contact is broken, it makes continuous contact look as
    // though some energy has been lost. In the "yanking" case above, without
    // including the conservative power term we will actually lose energy
    // because the deformed material isn't allowed to push back on us so the
    // energy is lost to surface vibrations or some other unmodeled effect.
    contactForce_S1.setPowerDissipation(powerHC + powerFriction);
}

void ContactForceGenerator::HertzCircular::calcContactPatch
   (const State&      state,
    const Contact&    overlap,
    const SpatialVec& V_S1S2,
    ContactPatch&     patch_S1) const
{  
    patch_S1.m_elements.clear(); // TODO no details yet
    calcContactForce(state,overlap,V_S1S2,patch_S1.m_resultant);
}



//==============================================================================
//                         ELASTIC FOUNDATION GENERATOR
//==============================================================================
void ContactForceGenerator::ElasticFoundation::calcContactForce
   (const State&            state,
    const Contact&          overlap,    // contains X_S1S2
    const SpatialVec&       V_S1S2,     // relative surface velocity
    ContactForce&           contactForce_S1) const
{
    calcContactForceAndDetails(state,overlap,V_S1S2,contactForce_S1,0);
}



void ContactForceGenerator::ElasticFoundation::calcContactPatch
   (const State&      state,
    const Contact&    overlap,
    const SpatialVec& V_S1S2,
    ContactPatch&     patch_S1) const
{
    calcContactForceAndDetails(state,overlap,V_S1S2,
        patch_S1.m_resultant, &patch_S1.m_elements);
}



void ContactForceGenerator::ElasticFoundation::calcContactForceAndDetails
   (const State&            state,
    const Contact&          overlap,    // contains X_S1S2
    const SpatialVec&       V_S1S2,     // relative surface velocity
    ContactForce&           contactForce_S1,
    Array_<ContactDetail>*  contactDetails_S1) const
{
    SimTK_ASSERT(TriangleMeshContact::isInstance(overlap),
        "ContactForceGenerator::ElasticFoundation::calcContactForce(): expected"
        " TriangleMeshContact.");

    const bool wantDetails = (contactDetails_S1 != 0);
    if (wantDetails)
        contactDetails_S1->clear();

    const TriangleMeshContact& contact = TriangleMeshContact::getAs(overlap);

    const CompliantContactSubsystem& subsys = getCompliantContactSubsystem();
    const ContactTrackerSubsystem&   tracker = subsys.getContactTrackerSubsystem();
    
    const ContactSurfaceIndex surf1x = contact.getSurface1();
    const ContactSurfaceIndex surf2x = contact.getSurface2();
    const Transform&          X_S1S2 = contact.getTransform();

    const ContactSurface&  surf1  = tracker.getContactSurface(surf1x);
    const ContactSurface&  surf2  = tracker.getContactSurface(surf2x);

    // Compute combined material properties just once -- they are the same
    // for all triangles.
    const ContactMaterial& mat1   = surf1.getMaterial();
    const ContactMaterial& mat2   = surf2.getMaterial();

    // Stiffnesses combine linearly here, not as 2/3 power as for Hertz.
    const Real k1=mat1.getStiffness(), k2=mat2.getStiffness();
    const Real c1=mat1.getDissipation(), c2=mat2.getDissipation();

    // Adjust the contact location based on the relative stiffness of the 
    // two materials. s1 is the fraction of the "squishing" (deformation) that
    // is done by surface1 -- if surface2 is very stiff then surface1 does most.
    const Real s1 = k2/(k1+k2); // 0..1
    const Real s2 = 1-s1;       // 1..0
    
    // k is the effective stiffness, c the effective dissipation
    const Real k = k1*s1; // (==k2*s2) == E
    const Real c = c1*s1 + c2*s2;

    // Calculate effective coefficients of friction, being careful not
    // to divide 0/0 if both are frictionless.
    const Real us1=mat1.getStaticFriction(), us2=mat2.getStaticFriction();
    const Real ud1=mat1.getDynamicFriction(), ud2=mat2.getDynamicFriction();
    const Real uv1=mat1.getViscousFriction(), uv2=mat2.getViscousFriction();
    Real us = 2*us1*us2; if (us!=0) us /= (us1+us2);
    Real ud = 2*ud1*ud2; if (ud!=0) ud /= (ud1+ud2);
    Real uv = 2*uv1*uv2; if (uv!=0) uv /= (uv1+uv2);

    // Now generate forces using the meshed surfaces only (one or two).
    const ContactGeometry& shape1 = surf1.getShape();
    const ContactGeometry& shape2 = surf2.getShape();

    // We want both patches to accumulate forces at the same point in
    // space. For numerical reasons this should be near the center of the
    // patch.
    Vec3 weightedPatchCentroid1_S1(0), weightedPatchCentroid2_S1(0);
    Real patchArea1 = 0, patchArea2 = 0;
    if (shape1.getTypeId() == ContactGeometry::TriangleMesh::classTypeId()) {
        const ContactGeometry::TriangleMesh& mesh1 = 
            ContactGeometry::TriangleMesh::getAs(shape1);

        calcWeightedPatchCentroid(mesh1, contact.getSurface1Faces(),
                                  weightedPatchCentroid1_S1, patchArea1);
    }
    if (shape2.getTypeId() == ContactGeometry::TriangleMesh::classTypeId()) {
        const ContactGeometry::TriangleMesh& mesh2 = 
            ContactGeometry::TriangleMesh::getAs(shape2);
        Vec3 weightedPatchCentroid2_S2;

        calcWeightedPatchCentroid(mesh2, contact.getSurface2Faces(),
                                  weightedPatchCentroid2_S2, patchArea2);
        // Remeasure patch2's weighted centroid from surface1's frame;
        // be sure to weight the new offset also.
        weightedPatchCentroid2_S1 = X_S1S2.R()*weightedPatchCentroid2_S2
                                    + patchArea2*X_S1S2.p();
    }

    // At this point one or two patch centroids are known; if one is unused
    // it is (0,0,0) with 0 weight. Combine them into a single composite
    // centroid for the patch.
    const Real patchArea = patchArea1+patchArea2;
    const Vec3 patchCentroid_S1 = 
        patchArea > 0 ? (  weightedPatchCentroid1_S1
                         + weightedPatchCentroid2_S1) / patchArea
                      : Vec3(0);

    // The patch centroid as calculated from one or two meshes is now
    // in patchCentroid_S1, measured from and expressed in S1. Now we'll
    // calculate all the patch forces and accumulate them at the patch
    // centroid.


    // Forces must be as applied to surface 2 at the patch centroid, but
    // expressed in surface 1's frame.
    SpatialVec force1_S1(Vec3(0)), force2_S1(Vec3(0));
    Real potEnergy1=0, potEnergy2=0;
    Real powerLoss1=0, powerLoss2=0;

    // Center of pressure r_c is calculated like this:
    //           sum_i (r_i * |r_i X Fn_i|) 
    //     r_c = --------------------------
    //              sum_i |r_i X Fn_i|
    // where Fn is the normal force applied by element i at location r_i,
    // with all locations measured from the patch centroid calculated above.
    //
    // We're going to calculate the weighted-points numerator for each
    // mesh separately in weightedCOP1 and weightedCOP2, and the corresponding
    // denominators (sum of all pressure-moment magnitudes) in weightCOP1 and
    // weightCOP2. At the end we'll combine and divide to calculate the actual 
    // center of pressure where we'll ultimately apply the contact force.
    Vec3 weightedCOP1_PC_S1(0), weightedCOP2_PC_S1(0); // from patch centroid
    Real weightCOP1=0, weightCOP2=0;

    if (shape1.getTypeId() == ContactGeometry::TriangleMesh::classTypeId()) {
        const ContactGeometry::TriangleMesh& mesh = 
            ContactGeometry::TriangleMesh::getAs(shape1);

        processOneMesh(state, 
            mesh, contact.getSurface1Faces(),
            X_S1S2, V_S1S2, shape2,
            s1, k, c, us, ud, uv,
            patchCentroid_S1,
            force1_S1, potEnergy1, powerLoss1,
            weightedCOP1_PC_S1, weightCOP1, 
            contactDetails_S1); // details returned if this is non-null
    }

    if (shape2.getTypeId() == ContactGeometry::TriangleMesh::classTypeId()) {
        const ContactGeometry::TriangleMesh& mesh = 
            ContactGeometry::TriangleMesh::getAs(shape2);

        // Costs 120 flops to flip everything into S2 for processing and
        // then put the results back in S1.

        const Transform  X_S2S1 = ~X_S1S2;                      //  3 flops
        const SpatialVec V_S2S1 = 
            reverseRelativeVelocity(X_S1S2,V_S1S2);             // 51 flops

        const Vec3 patchCentroid_S2 = X_S2S1*patchCentroid_S1;  // 18 flops

        SpatialVec       force2_S2; 
        Vec3             weightedCOP2_PC_S2;

        const unsigned nDetailsBefore = 
            wantDetails ? contactDetails_S1->size() : 0;

        processOneMesh(state, 
            mesh, contact.getSurface2Faces(),
            X_S2S1, V_S2S1, shape1,
            s2, k, c, us, ud, uv,
            patchCentroid_S2,
            force2_S2, potEnergy2, powerLoss2,
            weightedCOP2_PC_S2, weightCOP2,
            contactDetails_S1); // details returned *in S2* if this is non-null

        // Switch contact details from S1 in S2 to S2 in S1.
        if (wantDetails) {
            for (unsigned i=nDetailsBefore; i < contactDetails_S1->size(); ++i)
                (*contactDetails_S1)[i].changeFrameAndSwitchSurfacesInPlace(X_S1S2);
        }

        // Returned force is as applied to surface 1 (at the patch centroid); 
        // negate to make it the force applied to surface 2. Also re-express 
        // the moment and force vectors in S1.
        force2_S1 = -(X_S1S2.R()*force2_S2);                    // 36 flops
        // Returned weighted COP is measured from patch centroid but 
        // expressed in S2; reexpress in S1.
        weightedCOP2_PC_S1 = X_S1S2.R()*weightedCOP2_PC_S2;     // 15 flops
    }

    const Real weightCOP = weightCOP1+weightCOP2;
    Vec3 centerOfPressure_PC_S1 = // offset from patch centroid
        weightCOP > 0 ? (weightedCOP1_PC_S1+weightedCOP2_PC_S1) / weightCOP
                      : Vec3(0);

    // Calculate total force applied to surface 2 at the patch centroid,
    // expressed in S1.
    const SpatialVec force_PC_S1 = force1_S1 + force2_S1;

    // Shift to center of pressure.
    const SpatialVec forceCOP_S1 =
        shiftForceBy(force_PC_S1, centerOfPressure_PC_S1);

    // Contact point is the center of pressure measured from the S1 origin.
    const Vec3 contactPt_S1 = patchCentroid_S1 + centerOfPressure_PC_S1;

    // Fill in contact force object.
    contactForce_S1.setTo(contact.getContactId(),
        contactPt_S1,  // contact point in S1
        forceCOP_S1,   // force on surf2 at contact pt exp. in S1
        potEnergy1 + potEnergy2,
        powerLoss1 + powerLoss2);
}



// Compute the approximate geometric center of the patch represented by the 
// set of "inside faces" of this mesh, weighted by the total patch area. 
// We'll use the inside faces' centroids, weighted by face area, to calculate 
// a point that will be somewhere near the center of the patch, times the total
// area. Dividing by the area yields the actual patch centroid point which can 
// be used as a numerically good nearby point for accumulating all the forces, 
// since the moments will all be small  and thus
// not suffer huge roundoff errors when combined. Afterwards, the numerically
// good resultant can be shifted to whatever point is convenient, such as the
// body frame. We return the total area represented by this patch so that
// it can be combined with another mesh's patch if needed.
// The weighted centroid is returned in the mesh surface's own frame.
void ContactForceGenerator::ElasticFoundation::
calcWeightedPatchCentroid
   (const ContactGeometry::TriangleMesh&    mesh,
    const std::set<int>&                    insideFaces,
    Vec3&                                   weightedPatchCentroid,
    Real&                                   patchArea) const
{
    weightedPatchCentroid = Vec3(0); patchArea = 0;
    for (std::set<int>::const_iterator iter = insideFaces.begin(); 
                                       iter != insideFaces.end(); ++iter)
    {   const int  face = *iter;
        const Real area = mesh.getFaceArea(face);
        weightedPatchCentroid   += area*mesh.findCentroid(face); 
        patchArea               += area; 
    }
}



// Private method that calculates the net contact force produced by a single triangle
// mesh in contact with some other object (which might be another
// mesh; we don't care). We are given the relative spatial pose and velocity of
// the two surface frames, and for the mesh we are given a specific list of
// the faces that are suspected of being at least partially inside the 
// other surface (in the undeformed geometry overlap).
// Normally this just computes the resultant force but it can optionally append
// contact patch details (one entry per element) as well, if the contactDetails
// argument is non-null.
void ContactForceGenerator::ElasticFoundation::
processOneMesh
   (const State&                            state,
    const ContactGeometry::TriangleMesh&    mesh,
    const std::set<int>&                    insideFaces,
    const Transform&                        X_MO, 
    const SpatialVec&                       V_MO,
    const ContactGeometry&                  other,
    Real                                    meshDeformationFraction, // 0..1
    Real k, Real c, Real us, Real ud, Real uv, // composite material props
    const Vec3&                 resultantPt_M, // where to apply forces
    SpatialVec&                 resultantForceOnOther_M, // at resultant pt
    Real&                       potentialEnergy,
    Real&                       powerLoss,
    Vec3&                       weightedCenterOfPressure_M,
    Real&                       sumOfAllPressureMoments,   // COP weight
    Array_<ContactDetail>*      contactDetails_M) const    // in/out if present 
{
    assert(!insideFaces.empty());
    const bool wantDetails = (contactDetails_M != 0);

    // Initialize all results to zero.
    resultantForceOnOther_M = SpatialVec(Vec3(0),Vec3(0));
    potentialEnergy = powerLoss = 0;
    weightedCenterOfPressure_M = Vec3(0);
    sumOfAllPressureMoments = 0;

    // Don't initialize contact details; we're going to append them.

    // Abbreviations.
    const Rotation& RMO = X_MO.R(); // orientation of O in M
    const Vec3&     pMO = X_MO.p(); // position of OO (origin of O) in M
    const Vec3&     wMO = V_MO[0];  // ang. vel. of O in M
    const Vec3&     vMO = V_MO[1];  // vel. of OO in M

    // Get friction model transition velocity vt and 1/vt.
    const CompliantContactSubsystem& subsys = getCompliantContactSubsystem();
    const Real vtrans   = subsys.getTransitionVelocity();
    const Real ooVtrans = subsys.getOOTransitionVelocity(); // 1/vtrans

    // Now loop over all the faces again, evaluate the force from each 
    // spring, and apply it at the patch centroid.
    // This costs roughly 300 flops per contacting face.
    for (std::set<int>::const_iterator iter = insideFaces.begin(); 
                                       iter != insideFaces.end(); ++iter) 
    {   const int   face        = *iter;
        const Vec3  springPos_M = mesh.findCentroid(face);
        const Real  faceArea    = mesh.getFaceArea(face);

        bool        inside;
        UnitVec3    normal_O; // not used
        const Vec3  nearestPoint_O = // 18 flops + cost of findNearestPoint
            other.findNearestPoint(~X_MO*springPos_M, inside, normal_O);
        if (!inside)
            continue;
        
        // Although the "spring" is associated with just one surface (the mesh M)
        // it is considered here to include the compression of both surfaces
        // together, using composite material properties for stiffness and 
        // dissipation properties of the spring. The total displacement vector
        // for both surfaces points from the nearest point on the undeformed other 
        // surface to the undeformed spring position (face centroid) on the 
        // mesh. Since these overlap (we checked above) the nearest point is 
        // *inside* the mesh thus the vector points towards the mesh exterior; 
        // i.e., in the  direction that the force will be applied to the 
        // "other" body. This is the same convention we use for the patch 
        // normal for Hertz contact.
        const Vec3 nearestPoint_M = X_MO*nearestPoint_O; // 18 flops
        const Vec3 overlap_M      = springPos_M - nearestPoint_M; // 3 flops
        const Real overlap        = overlap_M.norm(); // ~40 flops

        // If there is no overlap we can't generate forces.
        if (overlap == 0)
            continue;

        // The surfaces are compressed by total amount "overlap".
        const UnitVec3 normal_M(overlap_M/overlap, true); // ~15 flops

        // Calculate the contact point location based on the relative 
        // squishiness of the two surfaces. The mesh deformation fraction (0-1)
        // gives the fraction of the material squishing that is done by the
        // mesh; the rest is done by the other surface. At 0 (rigid mesh) the 
        // contact point will be at the undeformed mesh face centroid; at 1 
        // (other body rigid) it will be at the (undeformed) nearest point on 
        // the other body.
        const Real meshSquish = meshDeformationFraction*overlap; // mesh displacement
        // Remember that the normal points towards the exterior of this mesh.
        const Vec3 contactPt_M = springPos_M - meshSquish*normal_M; // 6 flops
        
        // Calculate the relative velocity of the two bodies at the contact 
        // point. We're considering the mesh M fixed, so we just need the 
        // velocity in M of the station of O that is coincident with the 
        // contact point.

        // O station, exp. in M
        const Vec3 contactPtO_M = contactPt_M - pMO;    // 3 flops 

        // All vectors are in M; dropping the "_M" notation now.

        // Velocity of other at contact point is opposite direction of normal
        // when penetration is increasing.
        const Vec3 vel = vMO + wMO % contactPtO_M;      // 12 flops

        // Want odot > 0 when overlap is increasing; normal points the 
        // other way. odot is signed penetration (overlap) rate.
        const Real odot = -dot(vel, normal_M);          // 6 flops
        const Vec3 velNormal  = -odot*normal_M;         // 4 flops
        const Vec3 velTangent = vel-velNormal;          // 3 flops
        
        // Calculate scalar normal force                  (5 flops)
        const Real fK = k*faceArea*overlap; // normal elastic force (conservative)
        const Real fC = fK*c*odot;          // normal dissipation force (loss)
        const Real fNormal = fK + fC;       // normal force

        // Total force can be negative under unusual circumstances ("yanking");
        // that means no force is generated and no stored PE will be recovered.
        // This will most often occur in to-be-rejected trial steps but can
        // occasionally be real.
        if (fNormal <= 0) {
            SimTK_DEBUG1("YANKING!!! (face %d)\n", face);
            //printf("YANKING!!! (face %d)\n", face);
            continue;
        }

        // 12 flops in this series.
        const Vec3 forceK          = fK*normal_M;   // as applied to other surf
        const Vec3 forceC          = fC*normal_M;
        const Real PE              = fK*overlap/2;  // 1/2 kAx^2
        const Real powerC          = fC*odot;       // rate of energy loss, >= 0
        const Vec3 forceNormal     = forceK + forceC;

        // This is the moment r X f about the resultant point produced by 
        // applying this pure force at the contact point. Cost ~60 flops.
        const Vec3 r = contactPt_M - resultantPt_M;
        const Real pressureMoment = (r % forceNormal).norm();
        weightedCenterOfPressure_M += pressureMoment*r;
        sumOfAllPressureMoments    += pressureMoment;
        
        // Calculate the friction force. Cost is about 60 flops.
        Vec3 forceFriction(0);
        Real powerFriction = 0;
        const Real vslipSq = velTangent.normSqr();  // 5 flops
        if (vslipSq > square(SignificantReal)) {
            const Real vslip = std::sqrt(vslipSq); // expensive: ~25 flops
            // Express slip velocity as unitless multiple of transition velocity.
            const Real v = vslip * ooVtrans;
            // Must scale viscous coefficient to match unitless velocity.
            const Real mu=stribeck(us,ud,uv*vtrans,v); // ~10 flops
            //const Real mu=hollars(us,ud,uv*vtrans,v);
            const Real fFriction = fNormal * mu;
            // Force direction on O opposes O's velocity.
            forceFriction = (-fFriction/vslip)*velTangent; // ~20 flops
            powerFriction = fFriction * vslip; // always >= 0
        }

        const Vec3 forceLoss  = forceC + forceFriction;     // 3 flops
        const Vec3 forceTotal = forceK + forceLoss;         // 3 flops

        // Accumulate the moment and force on the *other* surface as though 
        // applied at the point of O that is coincident with the resultant
        // point; we'll move it later.                      (15 flops)
        resultantForceOnOther_M += SpatialVec(r % forceTotal, forceTotal);

        // Accumulate potential energy stored in elastic displacement.
        potentialEnergy += PE;                  // 1 flop

        // Don't include dot(forceK,velNormal) power due to conservative force
        // here. This way we don't double-count the energy on the way in as
        // integrated power and potential energy. Although the books would 
        // balance again when the contact is broken, it makes continuous 
        // contact look as though some energy has been lost. In the "yanking" 
        // case above, without including the conservative power term we will 
        // actually lose energy because the deformed material isn't allowed to 
        // push back on us so the energy is lost to surface vibrations or some
        // other unmodeled effect.
        const Real powerLossThisElement = powerC + powerFriction; // 1 flop
        powerLoss += powerLossThisElement;                        // 1 flop

        if (wantDetails) {
            contactDetails_M->push_back();
            ContactDetail& detail = contactDetails_M->back();
            detail.m_contactPt          = contactPt_M;
            detail.m_patchNormal        = normal_M;
            detail.m_slipVelocity       = velTangent;
            detail.m_forceOnSurface2    = forceTotal;
            detail.m_deformation        = overlap;
            detail.m_deformationRate    = odot;
            detail.m_patchArea          = faceArea;
            detail.m_peakPressure       = (faceArea != 0 ? fNormal/faceArea : Real(0));
            detail.m_potentialEnergy    = PE;
            detail.m_powerLoss          = powerLossThisElement;
        }
    }
}

} // namespace SimTK

