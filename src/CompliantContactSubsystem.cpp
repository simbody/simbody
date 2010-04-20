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
    m_tracker(tracker), m_transitionVelocity(0.01), m_defaultGenerator(0) 
{   
}

Real getTransitionVelocity() const  {return m_transitionVelocity;}
void setTransitionVelocity(Real vt) {m_transitionVelocity=vt;}

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
        const Contact& contact = contacts.getContactById(force.m_contactId);
        const MobilizedBody& mobod1 = m_tracker.getMobilizedBody
                                                (contact.getSurface1());
        const MobilizedBody& mobod2 = m_tracker.getMobilizedBody
                                                (contact.getSurface2());
        const Vec3 r1 = force.m_centerOfPressureInG 
                        - mobod1.getBodyOriginLocation(s);
        const Vec3 r2 = force.m_centerOfPressureInG 
                        - mobod2.getBodyOriginLocation(s);
        const SpatialVec& F2cop = force.m_forceOnSurface2InG; // at C.O.P.
        // Shift applied force to body origins.
        const SpatialVec F2( F2cop[0] + r2 %  F2cop[1],  F2cop[1]);
        const SpatialVec F1(-F2cop[0] + r1 % -F2cop[1], -F2cop[1]);
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
    Real power = 0;
    const Array_<ContactForce>& forces = getForceCache(state);
    for (unsigned i=0; i < forces.size(); ++i)
        power += forces[i].m_power;
    updDissipatedEnergyDeriv(state) = power;
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
Real                                m_transitionVelocity;

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
            pe += forces[i].m_potentialEnergy;
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
        generator.calcContactForce(state,contact,SpatialVec(Vec3(0)),
            SpatialVec(Vec3(0)), force);
        pe += force.m_potentialEnergy;
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
        const ContactSurfaceIndex surf1(contact.getSurface1());
        const ContactSurfaceIndex surf2(contact.getSurface2());
        const MobilizedBody& mobod1 = m_tracker.getMobilizedBody(surf1);
        const MobilizedBody& mobod2 = m_tracker.getMobilizedBody(surf2);
        const SpatialVec V_GS1 = mobod1.findFrameVelocityInGround
            (state, m_tracker.getContactSurfaceTransform(surf1));
        const SpatialVec V_GS2 = mobod2.findFrameVelocityInGround
            (state, m_tracker.getContactSurfaceTransform(surf2));

        const ContactForceGenerator& generator = 
            getForceGenerator(contact.getTypeId());
        forces.push_back();
        generator.calcContactForce(state,contact,V_GS1,V_GS2, forces.back());
        if (!forces.back().isValid())
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
void CompliantContactSubsystem::setTransitionVelocity(Real vt) {
    SimTK_ERRCHK1_ALWAYS(vt > 0, 
        "CompliantContactSubsystem::setTransitionVelocity()",
        "Friction transition velocity %g is illegal.", vt);
    updImpl().setTransitionVelocity(vt);
}


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






//==============================================================================
//                         HERTZ CIRCULAR GENERATOR
//==============================================================================

void ContactForceGenerator::HertzCircular::calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_GS1,  // surface velocities
    const SpatialVec&       V_GS2,
    ContactForce&           contactForce) const
{
    SimTK_ASSERT(CircularPointContact::isInstance(overlapping),
        "ContactForceGenerator::HertzCircular::calcContactForce(): expected"
        " CircularPointContact.");

    const CircularPointContact& contact =  
        CircularPointContact::getAs(overlapping);
    const Real depth = contact.getDepth();

    if (depth <= 0) {
        contactForce.clear(); // invalidate
        return;
    }

    const CompliantContactSubsystem& subsys = getCompliantContactSubsystem();
    const ContactTrackerSubsystem&   tracker = subsys.getContactTrackerSubsystem();

    const MobilizedBody& mobod1 = tracker.getMobilizedBody(contact.getSurface1());
    const MobilizedBody& mobod2 = tracker.getMobilizedBody(contact.getSurface2());  
    const ContactMaterial& mat1 =
        tracker.getContactSurface(contact.getSurface1()).getMaterial();
    const ContactMaterial& mat2 =
        tracker.getContactSurface(contact.getSurface1()).getMaterial();

    // Use 2/3 power of stiffness for combining here.
    const Real k1=mat1.getStiffness23(), k2=mat2.getStiffness23();
    const Real c1=mat1.getDissipation(), c2=mat2.getDissipation();

    // Adjust the contact location based on the relative stiffness of the 
    // two materials. s1 is the fraction of the "squishing" (deformation) that
    // is done by surface1 -- if surface2 is very stiff then surface1 does most.
    
    const Real s1 = k2/(k1+k2); // 0..1
    const Real s2 = 1-s1;       // 1..0
    const Real x  = contact.getDepth();
    // normal points from surf1 to surf2 
    const UnitVec3& normal = contact.getNormal();
    // origin is half way between the two surfaces as though they had
    // equal stiffness
    const Vec3& origin = contact.getOrigin();
    // Actual contact point moves closer to stiffer surface.
    const Vec3 contactPt = origin+(x*(0.5-s1))*normal;
    
    // Calculate the Hertz force fH, which is conservative.

    const Real k = k1*s1; // (==k2*s2) == E^(2/3)
    const Real c = c1*s1 + c2*s2;
    const Real R = contact.getEffectiveRadius();
    const Real fH = (4./3.)*k*x*std::sqrt(R*k*x);
    
    // Calculate the relative velocity of the two bodies at the contact point.
    //TODO: wrong velocity
    const Vec3 station1 = mobod1.findStationAtGroundPoint(state, contactPt);
    const Vec3 station2 = mobod2.findStationAtGroundPoint(state, contactPt);
    const Vec3 vel1 = mobod1.findStationVelocityInGround(state, station1);
    const Vec3 vel2 = mobod2.findStationVelocityInGround(state, station2);
    const Vec3 vel = vel1-vel2;
    const Real vnormal = dot(vel, normal);
    const Vec3 velNormal = vnormal*normal;
    const Vec3 velTangent = vel-velNormal;
    
    // Calculate the Hunt-Crossley force, which is dissipative.
    const Real fHC             = fH*1.5*c*vnormal;
    // Total force can be negative under unusual circumstances ("yanking");
    // that means no force is generated and no stored PE will be recovered.
    const Real f               = fH + fHC;
    Vec3 forceH(0), forceHC(0), forceFriction(0);
    Real powerHC=0, fFriction=0, potentialEnergy=0, powerFriction=0;
    if (f > 0) {
        forceH          = fH*normal;
        forceHC         = fHC*normal;
        potentialEnergy = (2./5.)*fH*x;
        powerHC         = dot(forceHC,velNormal);
    
        // Calculate the friction force.
    
        const Real vslipSq = velTangent.normSqr();
        if (vslipSq > square(SignificantReal)) {
            const Real vslip = std::sqrt(vslipSq); // expensive
            const Real us1=mat1.getStaticFriction(), us2=mat2.getStaticFriction();
            const Real ud1=mat1.getDynamicFriction(), ud2=mat2.getDynamicFriction();
            const Real uv1=mat1.getViscousFriction(), uv2=mat2.getViscousFriction();
            const Real vtrans = subsys.getTransitionVelocity();

            const bool hasStatic = (us1 != 0 || us2 != 0);
            const bool hasDynamic= (ud1 != 0 || ud2 != 0);
            const bool hasViscous = (uv1 != 0 || uv2 != 0);
            const Real us = hasStatic ? 2*us1*us2/(us1+us2) : 0;
            const Real ud = hasDynamic ? 2*ud1*ud2/(ud1+ud2) : 0;
            const Real uv = hasViscous ? 2*uv1*uv2/(uv1+uv2) : 0;
            const Real vrel = vslip/vtrans;
            const Real fFriction = 
                f*(std::min(vrel, 1.0)*(ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
            forceFriction = fFriction*velTangent/vslip;
            powerFriction = dot(forceFriction,velTangent);
        }
    } else {
        //std::cout << "YANKING!!!\n";
    }

    const Vec3 forceLoss  = forceHC + forceFriction;
    const Vec3 forceTotal = forceH + forceLoss;
    
    // Report the force.
    contactForce.m_contactId = contact.getContactId();
    contactForce.m_centerOfPressureInG = contactPt;
    contactForce.m_forceOnSurface2InG = SpatialVec(Vec3(0),forceTotal);
    contactForce.m_potentialEnergy = potentialEnergy;
    // Don't include dot(forceH,velNormal) power due to conservative force
    // here. This way we don't double-count the energy on the way in as
    // integrated power and potential energy. Although the books will balance
    // again when the contact is broken, it makes continuous contact look as
    // though some energy has been lost. In the "yanking" case above, without
    // including the conservative power term we will actually lose energy
    // because the deformed material isn't allowed to push back on us so the
    // energy is lost to surface vibrations or some other unmodeled effect.
    contactForce.m_power = powerHC + powerFriction;
}

void ContactForceGenerator::HertzCircular::calcContactPatch
   (const State&      state,
    const Contact&    overlapping,
    const SpatialVec& V_GS1,  // surface velocities
    const SpatialVec& V_GS2,
    ContactPatch&     patch) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactForceGenerator::HertzCircular::calcContactPatch() not implemented yet."); }




} // namespace SimTK

