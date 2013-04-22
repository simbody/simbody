/* -------------------------------------------------------------------------- *
 *                          Simbody(tm) - Tim's Box                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/* Attempt to duplicate Tim's OpenSim Box to diagnose nonrepeatability
and excessive sticking problems. */

//#define NDEBUG 1

#include "Simbody.h"

#include <string>
#include <iostream>
#include <exception>

using std::cout;
using std::endl;

using namespace SimTK;

//#define USE_TIMS_PARAMS
//#define USE_COMPLIANT_CONTACT
//#define USE_ELASTIC_FOUNDATION // if compliant
#define ANIMATE // off to get more accurate CPU time (you can still playback)

// Define this to run the simulation NTries times, saving the states and
// comparing them bitwise to see if the simulations are perfectly repeatable
// as they should be. You should see nothing but exact zeroes print out for
// second and subsequent runs.
//#define TEST_REPEATABILITY
static const int NTries=2;

const Real ReportInterval=1./30;
   

//==============================================================================
//                           UNI POINT IN PLANE IMPL
//==============================================================================
// This is Simbody's PointInPlane Constraint with sliding friction added.
// TODO: not using sliding forces here yet -- still using separate force
// element MySlidingFrictionForce (below).
class UniPointInPlaneImpl : public Constraint::Custom::Implementation {
public:
UniPointInPlaneImpl
   (MobilizedBody& plane,    const UnitVec3& defPlaneNormal, Real defPlaneHeight,
    MobilizedBody& follower, const Vec3&     defFollowerPoint)
:   Implementation(plane.updMatterSubsystem(),1,0,0) 
{
    planeBody    = addConstrainedBody(plane);
    followerBody = addConstrainedBody(follower);
    defaultPlaneNormal   = defPlaneNormal;
    defaultPlaneHeight   = defPlaneHeight;
    defaultFollowerPoint = defFollowerPoint;
}


void setPlaneDisplayHalfWidth(Real h) {
    // h <= 0 means don't display plane
    invalidateTopologyCache();
    planeHalfWidth = h > 0 ? h : 0;
}
Real getPlaneDisplayHalfWidth() const {return planeHalfWidth;}

void setPointDisplayRadius(Real r) {
    // r <= 0 means don't display point
    invalidateTopologyCache();
    pointRadius= r > 0 ? r : 0;
}
Real getPointDisplayRadius() const {return pointRadius;}

// Allocate state variables:
//  - s, the expected sign of the normal multiplier (-1,0,1)
//  - e, the "previous slip direction" (a 2d unit vector or NaN if none)
//
// The auto update for s is s=sign(|lambda|>=SignificantReal ? lambda : 0) so 
// that there is a dead zone in which we consider lambda to be zero.
//
// The auto update for e is v/|v| if slip speed |v| is large enough to be 
// considered "reliable", and dot(e,v)>0. Otherwise stick with the current e.
// A stiction->sliding transition causes e to be set to the impending slip
// direction (opposite the last known stiction force direction).
//
// Let d=s for bilateral constraint, d=max(s,0) for unilateral. Then
// sliding force multiplier sigma=d*mu_d*lambda, which will be >= 0 if d was
// correct. The sliding friction force will be -sigma*e.
//


void realizeTopology(State& state) const OVERRIDE_11 {
    isSlidingIx = getMatterSubsystem().allocateDiscreteVariable
       (state, Stage::Acceleration, new Value<bool>(false));
    prevSignIx = getMatterSubsystem().allocateAutoUpdateDiscreteVariable
       (state, Stage::Acceleration, new Value<Real>(1), 
        Stage::Acceleration); // cache update depends on lambda
    prevSlipDirIx = getMatterSubsystem().allocateAutoUpdateDiscreteVariable
       (state, Stage::Acceleration, new Value<Vec2>(Vec2(NaN)), 
        Stage::Velocity);     // cache update depends only on slip velocity
}

bool isSliding(const State& state) const {
    const bool& sliding = Value<bool>::downcast
        (getMatterSubsystem().getDiscreteVariable(state, isSlidingIx));
    return sliding;
}
Real getPrevSign(const State& state) const {
    const Real& prevSign = Value<Real>::downcast
        (getMatterSubsystem().getDiscreteVariable(state, prevSignIx));
    return prevSign;
}
Vec2 getPrevSlipDir(const State& state) const {
    const Vec2& prevSlipDir = Value<Vec2>::downcast
        (getMatterSubsystem().getDiscreteVariable(state, prevSlipDirIx));
    return prevSlipDir;
}
/*
// If we're sliding, set the update value for the previous slip direction
// if the current slip velocity is usable.
void realizeVelocity(const State& state) const OVERRIDE_11 {
    if (!hasSlidingForce(state))
        return; // nothing to do 
    const Vec2 Vslip = m_contact.getSlipVelocity(state);
    const Vec2 prevVslipDir = getPrevSlipDir(state);

    if (shouldUpdate(Vslip, prevVslipDir, m_vtol)) {
        Vec2& prevSlipUpdate = Value<Vec2>::updDowncast
            (m_forces.updDiscreteVarUpdateValue(state, m_prevSlipDirIx));
        const Real v = Vslip.norm();
        const Vec2 slipDir = Vslip / v;
        prevSlipUpdate = slipDir;
        m_forces.markDiscreteVarUpdateValueRealized(state, m_prevSlipDirIx);

        #ifndef NDEBUG
        //printf("UPDATE %d: prevSlipDir=%g %g; now=%g %g; |v|=%g dot=%g vdot=%g\n",
        //    m_friction.getFrictionIndex(),
        //    prevVslipDir[0],prevVslipDir[1],slipDir[0],slipDir[1],
        //    v, ~slipDir*prevVslipDir, ~Vslip*prevVslipDir);
        #endif
    } else {
        #ifndef NDEBUG
        printf("NO UPDATE %d: prevSlipDir=%g %g; Vnow=%g %g; |v|=%g vdot=%g\n",
            m_friction.getFrictionIndex(),
            prevVslipDir[0],prevVslipDir[1],Vslip[0],Vslip[1],
            Vslip.norm(), ~Vslip*prevVslipDir);
        #endif
    }
}

// Regardless of whether we're sticking or sliding, as long as the master
// contact is active use its normal force scalar as the update for our
// saved normal force.
void realizeAcceleration(const State& state) const OVERRIDE_11 {
    if (!m_friction.isMasterActive(state))
        return; // nothing to save
    const Real N = m_contact.getForce(state); // normal force
    const Real prevN = getPrevN(state);
    if (N==prevN) return; // no need for an update

    Real& prevNupdate = Value<Real>::updDowncast
        (m_forces.updDiscreteVarUpdateValue(state, m_prevNix));

    #ifndef NDEBUG
    printf("UPDATE %d: N changing from %g -> %g (%.3g)\n",
        m_friction.getFrictionIndex(), 
        prevN, N, std::abs(N-prevN)/std::max(N,prevN));
    #endif
    prevNupdate = N;
    m_forces.markDiscreteVarUpdateValueRealized(state, m_prevNix); 
}
*/

// Implementation of virtuals required for holonomic constraints.
UniPointInPlaneImpl* clone() const OVERRIDE_11
{   return new UniPointInPlaneImpl(*this); }

void calcDecorativeGeometryAndAppend
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    OVERRIDE_11 {}

void calcPositionErrors      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const OVERRIDE_11;

void calcPositionDotErrors      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const OVERRIDE_11;

void calcPositionDotDotErrors      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const OVERRIDE_11;

// apply f=lambda*n to the follower point S of body F,
//       -f         to point C (coincident point) of body B
void addInPositionConstraintForces
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const OVERRIDE_11
{
    assert(multipliers.size()==1 && bodyForcesInA.size()==2 
           && qForces.size()==0);
    const Real lambda = multipliers[0];

    //TODO: should be able to get p info from State
    const Vec3&      p_FS    = defaultFollowerPoint; // measured & expressed in F
    const Vec3       p_AS    = findStationLocationFromState(s, followerBody, 
                                                            defaultFollowerPoint);
    const Transform& X_AB    = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC    = ~X_AB * p_AS;         // measured & expressed in B
    const Vec3       force_A = X_AB.R()*(lambda*defaultPlaneNormal);

    addInStationForce(s, followerBody, p_FS,  force_A, bodyForcesInA);
    addInStationForce(s, planeBody,    p_BC, -force_A, bodyForcesInA);
}

//------------------------------------------------------------------------------
                                    private:

ConstrainedBodyIndex    planeBody;    // B1
ConstrainedBodyIndex    followerBody; // B2

UnitVec3                defaultPlaneNormal;   // on body 1, exp. in B1 frame
Real                    defaultPlaneHeight;
Vec3                    defaultFollowerPoint; // on body 2, exp. in B2 frame

// These are just for visualization
Real                    planeHalfWidth;
Real                    pointRadius;

mutable DiscreteVariableIndex   isSlidingIx;   // should generate sliding force?
mutable DiscreteVariableIndex   prevSignIx;    // last sign of lambda 
mutable DiscreteVariableIndex   prevSlipDirIx; // previous slip direction e
};


//==============================================================================
//                           MY CONTACT ELEMENT
//==============================================================================
// This abstract class hides the details about which kind of contact constraint
// we're dealing with, while giving us enough to work with for deciding what's
// on and off and generating impulses.
//
// There is always a scalar associated with the constraint for making 
// decisions. There may be a friction element associated with this contact.
class MyFrictionElement;
class MyContactElement {
public:
    enum ImpulseType {Compression,Expansion,Capture};

    MyContactElement(Constraint uni, Real multSign, Real coefRest) 
    :   m_uni(uni), m_multSign(multSign), m_coefRest(coefRest), 
        m_index(-1), m_friction(0),
        m_velocityDependentCOR(NaN), m_restitutionDone(false) 
    {   m_uni.setDisabledByDefault(true); }

    virtual ~MyContactElement() {}
    
    // (Re)initialize base & concrete class. If overridden, be sure to
    // invoke base class first.
    virtual void initialize() {
        setRestitutionDone(false); 
        m_velocityDependentCOR = NaN;
        m_Ic = m_Ie = m_I = 0;
    }

    // Provide a human-readable string identifying the type of contact
    // constraint.
    virtual String getContactType() const = 0;

    // These must be constructed so that a negative value means the 
    // unilateral constraint condition is violated.
    virtual Real getPerr(const State& state) const = 0;
    virtual Real getVerr(const State& state) const = 0;
    virtual Real getAerr(const State& state) const = 0;

    // This returns a point in the ground frame at which you might want to
    // say the constraint is "located", for purposes of display. This should
    // return something useful even if the constraint is currently off.
    virtual Vec3 whereToDisplay(const State& state) const = 0;

    // This is used by some constraints to collect position information that
    // may be used later to set instance variables when enabling the underlying
    // Simbody constraint. All constraints zero impulses here.
    virtual void initializeForImpact(const State& state, Real captureVelocity) { 
        if (-captureVelocity <= getVerr(state) && getVerr(state) < 0) {
            m_velocityDependentCOR = 0;
            SimTK_DEBUG3("CAPTURING %d because %g <= v=%g < 0\n",
                m_index, -captureVelocity, getVerr(state));
        } else {
            m_velocityDependentCOR = m_coefRest;
        }
        
        setRestitutionDone(false);        
        m_Ic = m_Ie = m_I = 0; }

    // Returns zero if the constraint is not currently enabled. Otherwise 
    // return the signed constraint force, with a negative value indicating
    // that the unilateral force condition is violated.
    Real getForce(const State& s) const {
        if (isDisabled(s)) return 0;
        const Vector mult = m_uni.getMultipliersAsVector(s);
        assert(mult.size() == 1);
        if (isNaN(mult[0]))
            printf("*** getForce(): mult is NaN\n");
        return m_multSign*mult[0];
    }

    bool isProximal(const State& state, Real posTol) const
    {   return !isDisabled(state) || getPerr(state) <= posTol; }
    bool isCandidate(const State& state, Real posTol, Real velTol) const
    {   return isProximal(state, posTol) && getVerr(state) <= velTol; }


    void enable(State& state) const {m_uni.enable(state);}
    void disable(State& state) const {m_uni.disable(state);}
    bool isDisabled(const State& state) const {return m_uni.isDisabled(state);}

    void setMyDesiredDeltaV(const State&    s,
                            Vector&         desiredDeltaV) const
    {   Vector myDesiredDV(1); myDesiredDV[0] = m_multSign*getVerr(s);
        m_uni.setMyPartInConstraintSpaceVector(s, myDesiredDV, 
                                                   desiredDeltaV); }

    void recordImpulse(ImpulseType type, const State& state,
                               const Vector& lambda) {
        Vector myImpulse(1);
        m_uni.getMyPartFromConstraintSpaceVector(state, lambda, myImpulse);
        const Real I = myImpulse[0];
        if (type==Compression) m_Ic = I;
        else if (type==Expansion) m_Ie = I;
        m_I += I;
    }

    // Impulse is accumulated internally.
    Real getImpulse()            const {return -m_multSign*m_I;}
    Real getCompressionImpulse() const {return -m_multSign*m_Ic;}
    Real getExpansionImpulse()   const {return -m_multSign*m_Ie;}

    Real getMyValueFromConstraintSpaceVector(const State& state,
                                             const Vector& lambda) const
    {   Vector myValue(1);
        m_uni.getMyPartFromConstraintSpaceVector(state, lambda, myValue);
        return -m_multSign*myValue[0]; }

    void setMyExpansionImpulse(const State& state,
                               Real         coefRest,
                               Vector&      lambda) const
    {   const Real I = coefRest * m_Ic;
        Vector myImp(1); myImp[0] = I;
        m_uni.setMyPartInConstraintSpaceVector(state, myImp, lambda); }


    Real getMaxCoefRest() const {return m_coefRest;}
    Real getEffectiveCoefRest() const {return m_velocityDependentCOR;}
    void setRestitutionDone(bool isDone) {m_restitutionDone=isDone;}
    bool isRestitutionDone() const {return m_restitutionDone;}

    // Record position within the set of unilateral contact constraints.
    void setContactIndex(int index) {m_index=index;}
    int getContactIndex() const {return m_index;}
    // If there is a friction element for which this is the master contact,
    // record it here.
    void setFrictionElement(MyFrictionElement& friction)
    {   m_friction = &friction; }
    // Return true if there is a friction element associated with this contact
    // element.
    bool hasFrictionElement() const {return m_friction != 0;}
    // Get the associated friction element.
    const MyFrictionElement& getFrictionElement() const
    {   assert(hasFrictionElement()); return *m_friction; }
    MyFrictionElement& updFrictionElement() const
    {   assert(hasFrictionElement()); return *m_friction; }

protected:
    Constraint          m_uni;
    const Real          m_multSign; // 1 or -1
    const Real          m_coefRest;

    int                 m_index; // contact index in unilateral constraint set
    MyFrictionElement*  m_friction; // if any (just a reference, not owned)

    // Runtime -- initialized at start of impact handler.
    Real m_velocityDependentCOR; // Calculated at start of impact 
    bool m_restitutionDone;
    Real m_Ic, m_Ie, m_I; // impulses
};



//==============================================================================
//                           MY FRICTION ELEMENT
//==============================================================================
// A Coulomb friction element consists of both a sliding force and a stiction 
// constraint, at most one of which is active. There is a boolean state variable 
// associated with each element that says whether it is in sliding or stiction,
// and that state can only be changed during event handling.
//
// Generated forces depend on a scalar normal force N that comes from a 
// separate "normal force master", which might be one of the following:
//  - a unilateral constraint
//  - a bilateral constraint 
//  - a mobilizer
//  - a compliant force element 
// If the master is an inactive unilateral constraint, or if N=0, then no 
// friction forces are generated. In this example, we're only going to use
// a unilateral contact constraint as the "normal force master".
//
// For all but the compliant normal force master, the normal force N is 
// acceleration-dependent and thus may be coupled to the force produced by a
// sliding friction element. This may require iteration to ensure consistency
// between the sliding friction force and its master contact's normal force.
//
// A Coulomb friction element depends on a scalar slip speed defined by the
// normal force master (this might be the magnitude of a generalized speed or
// slip velocity vector). When the slip velocity goes to zero, the stiction 
// constraint is enabled if its constraint force magnitude can be kept to
// mu_s*|N| or less. Otherwise, or if the slip velocity is nonzero, the sliding
// force is enabled instead and generates a force of constant magnitude mu_d*|N| 
// that opposes the slip direction, or impending slip direction, as defined by 
// the master.
//
// There are two witness functions generated: (1) in slip mode, observes slip 
// velocity reversal and triggers stiction, and (2) in stiction mode, observes
// stiction force increase past mu_s*|N| and triggers switch to sliding.
class MyFrictionElement {
public:
    MyFrictionElement(Real mu_d, Real mu_s, Real mu_v)
    :   mu_d(mu_d), mu_s(mu_s), mu_v(mu_v), m_index(-1) {}

    virtual ~MyFrictionElement() {}

    // (Re)initialize base & concrete class. If overridden, be sure to
    // invoke base class first.
    virtual void initialize() {
    }

    Real getDynamicFrictionCoef() const {return mu_d;}
    Real getStaticFrictionCoef()  const {return mu_s;}
    Real getViscousFrictionCoef() const {return mu_v;}

    // Return true if the stiction constraint is enabled.
    virtual bool isSticking(const State&) const = 0;

    virtual void enableStiction(State&) const = 0;
    virtual void disableStiction(State&) const = 0;

    // When sticking, record -f/|f| as the previous slip direction, and 
    // max(N,0) as the previous normal force. Stiction
    // must be currently active and constraint multipliers available.
    virtual void recordImpendingSlipInfo(const State&) = 0;
    // When sliding, record current slip velocity as the previous slip 
    // direction.
    virtual void recordSlipDir(const State&) = 0;

    // In an event handler or at initialization only, set the last recorded slip
    // direction as the previous direction. This invalidates Velocity stage.
    virtual void updatePreviousSlipDirFromRecorded(State& state) const = 0;

    // This is the dot product of the current sliding velocity and the
    // saved previous slip direction. This changes sign when a sliding friction
    // force of mu_d*|N| would cause a reversal, meaning a switch to stiction is
    // in order. State must be realized to Velocity stage.
    virtual Real calcSlipSpeedWitness(const State&) const = 0;

    // When in stiction, this calculates mu_s*|N| - |f|, which is negative if
    // the stiction force exceeds its limit. (Not suitable for impacts where
    // the dynamic coefficient should be used.) State must be realized to
    // Acceleration stage.
    virtual Real calcStictionForceWitness(const State&) const = 0;

    // This is the magnitude of the current slip velocity. State must be 
    // realized to Velocity stage.
    virtual Real getActualSlipSpeed(const State&) const = 0;

    // This is the magnitude of the current friction force, whether sliding
    // or sticking. State must be realized to Acceleration stage.
    virtual Real getActualFrictionForce(const State&) const = 0;

    // Return the scalar normal force N being generated by the contact master
    // of this friction element. This may be negative if the master is a
    // unilateral constraint whose "no-stick" condition is violated. 
    virtual Real getMasterNormalForce(const State&) const = 0;

    // If the sliding force element uses an estimated normal forces to generate
    // an approximate sliding force, return the estimate it is using here.
    virtual Real getEstimatedNormalForce(const State&) const = 0;

    // Modify the estimated normal force in the state.
    virtual void setEstimatedNormalForce(State&, Real N) const = 0;

    // Return true if the normal force master *could* be involved in an 
    // impact event (because it is touching).
    virtual bool isMasterProximal(const State&, Real posTol) const = 0;
    // Return true if the normal force master *could* be involved in contact
    // force generation (because it is touching and not separating).
    virtual bool isMasterCandidate(const State&, Real posTol, Real velTol)
        const = 0;
    // Return true if the normal force master is currently generating a
    // normal force (or impulse) so that this friction element might be 
    // generating a force also.
    virtual bool isMasterActive(const State&) const = 0;


    // This is used by some stiction constraints to collect position information
    // that may be used later to set instance variables when enabling the 
    // underlying Simbody constraint. Recorded impulses should be zeroed.
    virtual void initializeForStiction(const State& state) = 0; 

    // If this friction element's stiction constraint is enabled, set its
    // constraint-space velocity entry(s) in desiredDeltaV to the current
    // slip velocity (which might be a scalar or 2-vector).
    virtual void setMyDesiredDeltaV(const State& s,
                                    Vector&      desiredDeltaV) const = 0;

    // We just applied constraint-space impulse lambda to all active 
    // constraints. If this friction element's stiction constraint is enabled,
    // save its part of the impulse internally for reporting.
    virtual void recordImpulse(MyContactElement::ImpulseType type, 
                               const State& state,
                               const Vector& lambda) = 0;

    // Output the status, friction force, slip velocity, prev slip direction
    // (scalar or vector) to the given ostream, indented as indicated and 
    // followed by a newline. May generate multiple lines.
    virtual std::ostream& writeFrictionInfo(const State& state,
                                            const String& indent,
                                            std::ostream& o) const = 0;

    // Optional: give some kind of visual representation for the friction force.
    virtual void showFrictionForce(const State& state, 
        Array_<DecorativeGeometry>& geometry) const {}


    void setFrictionIndex(int index) {m_index=index;}
    int getFrictionIndex() const {return m_index;}

private:
    Real mu_d, mu_s, mu_v;
    int  m_index; // friction index within unilateral constraint set
};



//==============================================================================
//                       MY UNILATERAL CONSTRAINT SET
//==============================================================================

// These are indices into the unilateral constraint set arrays.
struct MyElementSubset {
    void clear() {m_contact.clear();m_friction.clear();m_sliding.clear();}
    Array_<int> m_contact;
    Array_<int> m_friction; // friction elements that might stick
    Array_<int> m_sliding;  // friction elements that can only slide
};

class MyUnilateralConstraintSet {
public:
    // Capture velocity is used two ways: (1) if the normal approach velocity
    // is smaller, the coefficient of restitution is set to zero for the 
    // upcoming impact, and (2) if a slip velocity is smaller than this the
    // contact is a candidate for stiction.
    MyUnilateralConstraintSet(const MultibodySystem& mbs, Real captureVelocity)
    :   m_mbs(mbs), m_captureVelocity(captureVelocity) {}

    // This class takes over ownership of the heap-allocated contact element.
    int addContactElement(MyContactElement* contact) {
        const int index = (int)m_contact.size();
        m_contact.push_back(contact);
        contact->setContactIndex(index);
        return index;
    }
    // This class takes over ownership of the heap-allocated friction element.
    int addFrictionElement(MyFrictionElement* friction) {
        const int index = (int)m_friction.size();
        m_friction.push_back(friction);
        friction->setFrictionIndex(index);
        return index;
    }

    Real getCaptureVelocity() const {return m_captureVelocity;}
    void setCaptureVelocity(Real v) {m_captureVelocity=v;}

    int getNumContactElements() const {return (int)m_contact.size();}
    int getNumFrictionElements() const {return (int)m_friction.size();}
    const MyContactElement& getContactElement(int ix) const 
    {   return *m_contact[ix]; }
    const MyFrictionElement& getFrictionElement(int ix) const 
    {   return *m_friction[ix]; }

    // Allow writable access to elements from const set so we can record
    // runtime results (e.g. impulses).
    MyContactElement&  updContactElement(int ix) const {return *m_contact[ix];}
    MyFrictionElement& updFrictionElement(int ix) const {return *m_friction[ix];}

    // Initialize all runtime fields in the contact & friction elements.
    void initialize()
    {
        for (unsigned i=0; i < m_contact.size(); ++i)
            m_contact[i]->initialize();
        for (unsigned i=0; i < m_friction.size(); ++i)
            m_friction[i]->initialize();
    }

    // Return the contact and friction elements that might be involved in an
    // impact occurring in this configuration. They are the contact elements 
    // for which perr <= posTol, and friction elements whose normal force 
    // masters can be involved in the impact. State must be realized through 
    // Position stage.
    void findProximalElements(const State&      state,
                              Real              posTol,
                              MyElementSubset&  proximals) const
    {
        proximals.clear();
        for (unsigned i=0; i < m_contact.size(); ++i)
            if (m_contact[i]->isProximal(state,posTol)) 
                proximals.m_contact.push_back(i);
        for (unsigned i=0; i < m_friction.size(); ++i)
            if (m_friction[i]->isMasterProximal(state,posTol))
                proximals.m_friction.push_back(i);
        // Any friction elements might stick if they are proximal since
        // we'll be changing velocities, so no m_sliding entries in proximals.
    }

    // Return the contact and friction elements that might be involved in 
    // generating contact forces at the current state. Candidate contact
    // elements are those that are (a) already enabled, or (b) for which 
    // perr <= posTol and verr <= velTol. Candidate friction elements are those
    // whose normal force master is unconditional or a candidate and (a) which 
    // are already sticking, or (b) for which vslip <= velTol, or (c) for which
    // vslip opposes the previous slip direction, meaning it has reversed and 
    // must have passed through zero during the last step. These are the elements 
    // that can be activated without making any changes to the configuration or 
    // velocity state variables, except slightly for constraint projection. 
    //
    // We also record the friction elements that, if their masters are active, 
    // can only slide because they have a significant slip velocity. State must 
    // be realized through Velocity stage.
    void findCandidateElements(const State&     s,
                               Real             posTol,
                               Real             velTol,
                               MyElementSubset& candidates) const
    {
        candidates.clear();
        for (unsigned i=0; i < m_contact.size(); ++i)
            if (m_contact[i]->isCandidate(s,posTol,velTol)) 
                candidates.m_contact.push_back(i);
        for (unsigned i=0; i < m_friction.size(); ++i) {
            MyFrictionElement& fric = updFrictionElement(i);
            if (!fric.isMasterCandidate(s,posTol,velTol))
                continue;
            if (fric.isSticking(s) 
                || fric.getActualSlipSpeed(s) <= velTol
                || fric.calcSlipSpeedWitness(s) <= 0) 
            {
                fric.initializeForStiction(s);
                candidates.m_friction.push_back(i); // could stick or slide
            } else {
                fric.recordSlipDir(s);
                candidates.m_sliding.push_back(i);  // could only slide
            }
        }
    }

    // Look through the given constraint subset and enable any constraints
    // that are currently disabled. Returns true if any change was made.
    // If includeStiction==false, we'll only enable contact constraints.
    bool enableConstraintSubset(const MyElementSubset& subset,
                                bool                   includeStiction,
                                State&                 state) const
    {
        bool changedSomething = false;

        // Enable contact constraints.
        for (unsigned i=0; i < subset.m_contact.size(); ++i) {
            const int which = subset.m_contact[i];
            const MyContactElement& cont = getContactElement(which);
            if (cont.isDisabled(state)) {
                cont.enable(state);
                changedSomething = true;
            }
        }

        if (includeStiction) {
            // Enable all stiction constraints.
            for (unsigned i=0; i < subset.m_friction.size(); ++i) {
                const int which = subset.m_friction[i];
                const MyFrictionElement& fric = getFrictionElement(which);
                if (!fric.isSticking(state)) {
                    assert(fric.isMasterActive(state));
                    fric.enableStiction(state);
                    changedSomething = true;
                }
            }
        }

        m_mbs.realize(state, Stage::Instance);
        return changedSomething;
    }

    // All event handlers call this method before returning. Given a state for
    // which no (further) impulse is required, here we decide which contact and
    // stiction constraints are active, and ensure that they satisfy the 
    // required constraint tolerances to the given accuracy. For sliding 
    // contacts, we will have recorded the slip or impending slip direction and 
    // converged the normal forces.
    // TODO: in future this may return indicating that an impulse is required
    // after all, as in Painleve's paradox.
    void selectActiveConstraints(State& state, Real accuracy) const;

    // This is the inner loop of selectActiveConstraints(). Given a set of
    // candidates to consider, it finds an active subset and enables those
    // constraints.
    void findActiveCandidates(State&                 state, 
                              const MyElementSubset& candidates) const;

    // If there are sliding contacts, adjust normal force estimates and then
    // recalculate actual normal forces until they agree to a tolerance.
    void convergeNormalForcesForSlidingContacts(State& state, Real rtol) const;

    // In Debug mode, produce a useful summary of the current state of the
    // contact and friction elements.
    void showConstraintStatus(const State& state, const String& place) const;

    ~MyUnilateralConstraintSet() {
        for (unsigned i=0; i < m_contact.size(); ++i)
            delete m_contact[i];
        for (unsigned i=0; i < m_friction.size(); ++i)
            delete m_friction[i];
    }

    const MultibodySystem& getMultibodySystem() const {return m_mbs;}
private:
    const MultibodySystem&      m_mbs;
    Real                        m_captureVelocity;
    Array_<MyContactElement*>   m_contact;
    Array_<MyFrictionElement*>  m_friction;
};



//==============================================================================
//                               STATE SAVER
//==============================================================================
// This reporter is called now and again to save the current state so we can
// play back a movie at the end.
class StateSaver : public PeriodicEventReporter {
public:
    StateSaver(const MultibodySystem&                   mbs,
               const MyUnilateralConstraintSet&         unis,
               const Integrator&                        integ,
               Real                                     reportInterval)
    :   PeriodicEventReporter(reportInterval), 
        m_mbs(mbs), m_unis(unis), m_integ(integ) 
    {   m_states.reserve(2000); }

    ~StateSaver() {}

    void clear() {m_states.clear();}
    int getNumSavedStates() const {return (int)m_states.size();}
    const State& getState(int n) const {return m_states[n];}

    void handleEvent(const State& s) const {
        const SimbodyMatterSubsystem& matter=m_mbs.getMatterSubsystem();
        const SpatialVec PG = matter.calcSystemMomentumAboutGroundOrigin(s);
        m_mbs.realize(s, Stage::Acceleration);

#ifndef NDEBUG
        printf("%3d: %5g mom=%g,%g E=%g", m_integ.getNumStepsTaken(),
            s.getTime(),
            PG[0].norm(), PG[1].norm(), m_mbs.calcEnergy(s));
        cout << " Triggers=" << s.getEventTriggers() << endl;
        m_unis.showConstraintStatus(s, "STATE SAVER");
#endif

        m_states.push_back(s);
    }
private:
    const MultibodySystem&                  m_mbs;
    const MyUnilateralConstraintSet&        m_unis;
    const Integrator&                       m_integ;
    mutable Array_<State>                   m_states;
};

// A periodic event reporter that does nothing; useful for exploring the
// effects of interrupting the simulation.
class Nada : public PeriodicEventReporter {
public:
    explicit Nada(Real reportInterval)
    :   PeriodicEventReporter(reportInterval) {} 

    void handleEvent(const State& s) const {
#ifndef NDEBUG
        printf("%7g NADA\n", s.getTime());
#endif
    }
};


//==============================================================================
//                          CONTACT ON HANDLER
//==============================================================================
// Allocate three of these for each unilateral contact constraint, using
// a position, velocity, or acceleration witness function. When the associated
// contact constraint is inactive, the event triggers are:
// 1. separation distance goes from positive to negative
// 2. separation rate goes from positive to negative while distance is zero
// 3. separation acceleration goes from positive to negative while both 
//    distance and rate are zero
// The first two cases may require an impulse, since the velocities may have to
// change discontinuously to satisfy the constraints. Case 3 requires only
// recalculation of the active contacts. In any case the particular contact
// element that triggered the handler is irrelevant; all "proximal" contacts
// are solved simultaneously.
class ContactOn: public TriggeredEventHandler {
public:
    ContactOn(const MultibodySystem&            system,
              const MyUnilateralConstraintSet&  unis,
              unsigned                          which,
              Stage                             stage) 
    :   TriggeredEventHandler(stage), 
        m_mbs(system), m_unis(unis), m_which(which),
        m_stage(stage)
    { 
        // Trigger only as height goes from positive to negative.
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    // This is the witness function.
    Real getValue(const State& state) const {
        const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
        const MyContactElement& uni = m_unis.getContactElement(m_which);
        if (!uni.isDisabled(state)) 
            return 0; // already locked

        const Real height = uni.getPerr(state);
        //printf("getValue %d(%.17g) perr=%g\n", m_which, state.getTime(), height);

        if (m_stage == Stage::Position)
            return height;

        // Velocity and acceleration triggers are not needed if we're
        // above ground.
        if (height > 0) return 0;

        const Real dheight = uni.getVerr(state);
        //printf("... verr=%g\n", dheight);

        if (m_stage == Stage::Velocity)
            return dheight;

        // Acceleration trigger is not needed if velocity is positive.
        if (dheight > 0) return 0;

        const Real ddheight = uni.getAerr(state);
        //printf("... aerr=%g\n", ddheight);

        return ddheight;
    }

    // We're using Poisson's definition of the coefficient of 
    // restitution, relating impulses, rather than Newton's, 
    // relating velocities, since Newton's can produce non-physical 
    // results for a multibody system. For Poisson, calculate the impulse
    // that would bring the velocity to zero, multiply by the coefficient
    // of restitution to calculate the rest of the impulse, then apply
    // both impulses to produce changes in velocity. In most cases this
    // will produce the same rebound velocity as Newton, but not always.
    void handleEvent(State& s, Real accuracy, bool& shouldTerminate) const;

    // Given the set of proximal constraints, prevent penetration by applying
    // a nonnegative least squares impulse generating a step change in 
    // velocity. On return, the applied impulse and new velocities are recorded
    // in the proximal elements, and state is updated to the new velocities and 
    // realized through Velocity stage. Constraints that ended up in contact
    // are enabled, those that rebounded are disabled.
    void processCompressionPhase(MyElementSubset&   proximal,
                                 State&             state) const;

    // Given a solution to the compression phase, including the compression
    // impulse, the set of impacters (enabled) and rebounders (disabled and
    // with positive rebound velocity), apply an expansion impulse based on
    // the effective coefficients of restitution of the impacters. Wherever
    // restitution is applied, the effective coefficient is reset to zero so
    // that further restitution will not be done for that contact. Returns
    // true if any expansion was done; otherwise nothing has changed.
    // Expansion may result in some negative velocities, in which case it has
    // induced further compression so another compression phase is required.
    bool processExpansionPhase(MyElementSubset& proximal,
                               State&           state) const;

    // Given only the subset of proximal constraints that are active, calculate
    // the impulse that would eliminate all their velocity errors. No change is
    // made to the set of active constraints. Some of the resulting impulses
    // may be negative.
    void calcStoppingImpulse(const MyElementSubset& proximal,
                             const State&           state,
                             Vector&                lambda0) const;

    // Given the initial generalized speeds u0, and a constraint-space impulse
    // lambda, calculate the resulting step velocity change du, modify the
    // generalized speeds in state to u0+du, and realize Velocity stage.
    void updateVelocities(const Vector& u0, 
                          const Vector& lambda, 
                          State&        state) const;


private:
    const MultibodySystem&              m_mbs; 
    const MyUnilateralConstraintSet&    m_unis;
    const unsigned                      m_which;
    const Stage                         m_stage;
};



//==============================================================================
//                          CONTACT OFF HANDLER
//==============================================================================
// Allocate one of these for each unilateral contact constraint. This handler 
// is invoked when an active contact constraint's contact force crosses zero
// from positive to negative, meaning it has gone from pushing to sticking.
// This simply invokes recalculation of the active contacts; the particular
// source of the event trigger doesn't matter.
class ContactOff: public TriggeredEventHandler {
public:
    ContactOff(const MultibodySystem&       system,
        const MyUnilateralConstraintSet&    unis,
        unsigned                            which) 
    :   TriggeredEventHandler(Stage::Acceleration), 
        m_mbs(system), m_unis(unis), m_which(which)
    { 
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    // This is the witness function.
    Real getValue(const State& state) const {
        const MyContactElement& uni = m_unis.getContactElement(m_which);
        if (uni.isDisabled(state)) return 0;
        const Real f = uni.getForce(state);
        return f;
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const 
    {
        SimTK_DEBUG2("\nhandle %d liftoff@%.17g\n", m_which, s.getTime());
        SimTK_DEBUG("\n----------------------------------------------------\n");
        SimTK_DEBUG2("LIFTOFF triggered by constraint %d @t=%.15g\n", 
            m_which, s.getTime());
        m_mbs.realize(s, Stage::Acceleration);

        #ifndef NDEBUG
        cout << " triggers=" << s.getEventTriggers() << "\n";
        #endif

        m_unis.selectActiveConstraints(s, accuracy);

        SimTK_DEBUG("LIFTOFF DONE.\n");
        SimTK_DEBUG("----------------------------------------------------\n");
    }

private:
    const MultibodySystem&              m_mbs; 
    const MyUnilateralConstraintSet&    m_unis;
    const unsigned                      m_which; // one of the contact elements
};



//==============================================================================
//                             MY POINT CONTACT
//==============================================================================
// Define a unilateral constraint to represent contact of a point on a moving
// body with the ground plane. The ground normal is assumed to be +y.
class MyPointContact : public MyContactElement {
    typedef MyContactElement Super;
public:
    MyPointContact(MobilizedBody& body, const Vec3& point, 
                   Real coefRest)
    :   MyContactElement( 
            //Constraint::PointInPlane(updGround(body), UnitVec3(YAxis), Zero,
            //                         body, point),
            Constraint::Custom(new UniPointInPlaneImpl(
                                     updGround(body), UnitVec3(YAxis), Zero,
                                     body, point)),
             Real(-1), // multiplier sign
             coefRest),
        m_body(body), m_point(point)
    {
    }

    Real getPerr(const State& s) const OVERRIDE_11 {
        const Vec3 p = m_body.findStationLocationInGround(s, m_point);
        return p[YAxis];
    }
    Real getVerr(const State& s) const OVERRIDE_11 {
        const Vec3 v = m_body.findStationVelocityInGround(s, m_point);
        return v[YAxis];
    }
    Real getAerr(const State& s) const OVERRIDE_11 {
        const Vec3 a = m_body.findStationAccelerationInGround(s, m_point);
        return a[YAxis];
    }

    String getContactType() const OVERRIDE_11 {return "Point";}
    Vec3 whereToDisplay(const State& state) const OVERRIDE_11 {
        return m_body.findStationLocationInGround(state,m_point);
    }

    // Will be zero if the stiction constraints are on.
    Vec2 getSlipVelocity(const State& s) const {
        const Vec3 v = m_body.findStationVelocityInGround(s, m_point);
        return Vec2(v[XAxis],v[ZAxis]);
    }
    // Will be zero if the stiction constraints are on.
    Vec2 getSlipAcceleration(const State& s) const {
        const Vec3 a = m_body.findStationAccelerationInGround(s, m_point);
        return Vec2(a[XAxis],a[ZAxis]);
    }

    Vec3 getContactPointInPlaneBody(const State& s) const
    {   return m_body.findStationLocationInGround(s, m_point); }

    const MobilizedBody& getBody() const {return m_body;}
    MobilizedBody& updBody() {return m_body;}
    const Vec3& getBodyStation() const {return m_point;}

    const MobilizedBody& getPlaneBody() const  {
        const SimbodyMatterSubsystem& matter = m_body.getMatterSubsystem();
        return matter.getGround();
    }
    MobilizedBody& updPlaneBody() const {return updGround(m_body);}

private:
    // For use during construction before m_body is set.
    MobilizedBody& updGround(MobilizedBody& body) const {
        SimbodyMatterSubsystem& matter = body.updMatterSubsystem();
        return matter.updGround();
    }

    MobilizedBody&    m_body;
    const Vec3        m_point;
};




//==============================================================================
//               MY SLIDING FRICTION FORCE -- Declaration
//==============================================================================

// A nice handle for the sliding friction force. The real code is in the Impl
// class defined at the bottom of this file.
class MySlidingFrictionForce : public Force::Custom {
public:
    // Add a sliding friction force element to the given force subsystem,
    // and associate it with a particular contact point.
    MySlidingFrictionForce(GeneralForceSubsystem&               forces,
                           const class MyPointContactFriction&  ptFriction,
                           Real                                 vtol);

    void setPrevN(State& state, Real N) const;
    // This should be a unit vector.
    void setPrevSlipDir(State& state, const Vec2& slipDir) const;

    Real getPrevN(const State& state) const;
    Vec2 getPrevSlipDir(const State& state) const;

    bool hasPrevSlipDir(const State& state) const;

    Real calcSlidingForceMagnitude(const State& state) const; 
    Vec2 calcSlidingForce(const State& state) const;

private:
    const class MySlidingFrictionForceImpl& getImpl() const;
};


//==============================================================================
//                        MY POINT CONTACT FRICTION
//==============================================================================
// This friction element expects its master to be a unilateral point contact 
// constraint. It provides slipping forces or stiction constraint forces acting
// in the plane, based on the normal force being applied by the point contact 
// constraint.
class MyPointContactFriction : public MyFrictionElement {
    typedef MyFrictionElement Super;
public:
    // The constructor allocates two NoSlip1D constraints and a sliding
    // friction force element.
    MyPointContactFriction(MyPointContact& contact,
        Real mu_d, Real mu_s, Real mu_v, Real vtol, //TODO: shouldn't go here
        GeneralForceSubsystem& forces)
    :   MyFrictionElement(mu_d,mu_s,mu_v), m_contact(contact),
        m_noslipX(contact.updPlaneBody(), Vec3(0), UnitVec3(XAxis), 
                  contact.updPlaneBody(), contact.updBody()),
        m_noslipZ(contact.updPlaneBody(), Vec3(0), UnitVec3(ZAxis), 
                  contact.updPlaneBody(), contact.updBody())
    {
        assert((0 <= mu_d && mu_d <= mu_s) && (0 <= mu_v));
        contact.setFrictionElement(*this);
        m_noslipX.setDisabledByDefault(true);
        m_noslipZ.setDisabledByDefault(true);
        m_sliding = new MySlidingFrictionForce(forces, *this, vtol);
        initializeRuntimeFields();
    }

    ~MyPointContactFriction() {delete m_sliding;}

    void initialize() OVERRIDE_11 {
        Super::initialize();
        initializeRuntimeFields();
    }

    // The way we constructed the NoSlip1D constraints makes the multipliers be
    // the force on Ground; we negate here so we'll get the force on the sliding
    // body instead.
    Vec2 getStictionForce(const State& s) const {
        assert(isSticking(s));
        return Vec2(-m_noslipX.getMultiplier(s), -m_noslipZ.getMultiplier(s));
    }

    // Implement pure virtuals from MyFrictionElement base class.

    bool isSticking(const State& s) const OVERRIDE_11
    {   return !m_noslipX.isDisabled(s); } // X,Z always on or off together

    // Note that initializeForStiction() must have been called first.
    void enableStiction(State& s) const OVERRIDE_11
    {   m_noslipX.setContactPoint(s, m_contactPointInPlane);
        m_noslipZ.setContactPoint(s, m_contactPointInPlane);
        m_noslipX.enable(s); m_noslipZ.enable(s); }

    void disableStiction(State& s) const OVERRIDE_11
    {   m_sliding->setPrevN(s, std::max(m_prevN, Real(0)));
        m_sliding->setPrevSlipDir(s, m_prevSlipDir);
        m_noslipX.disable(s); m_noslipZ.disable(s); }

    // When sticking with stiction force f, record -f/|f| as the previous slip 
    // direction. If the force is zero we'll leave the direction unchanged.
    // Also record the master's normal force as the previous normal force
    // unless it is negative; in that case record zero.
    // State must be realized through Acceleration stage.
    void recordImpendingSlipInfo(const State& s) OVERRIDE_11 {
        const Vec2 f = getStictionForce(s);
        SimTK_DEBUG3("%d: RECORD IMPENDING, f=%g %g\n", 
            getFrictionIndex(), f[0], f[1]);
        const Real fmag = f.norm();
        if (fmag > 0) // TODO: could this ever be zero?
            m_prevSlipDir = -f/fmag;
        const Real N = getMasterNormalForce(s); // might be negative
        m_prevN = N;
    }
    // When sliding, record current slip velocity (normalized) as the previous 
    // slip direction. If there is no slip velocity we leave the slip direction
    // unchanged. State must be realized through Velocity stage.
    void recordSlipDir(const State& s) OVERRIDE_11 {
        assert(!isSticking(s));
        Vec2 v = m_contact.getSlipVelocity(s);
        Real vmag = v.norm();
        if (vmag > 0)
            m_prevSlipDir = v/vmag;
    }

    // Transfer the locally-stored previous slip direction to the state variable.
    void updatePreviousSlipDirFromRecorded(State& s) const OVERRIDE_11 {
        m_sliding->setPrevSlipDir(s, m_prevSlipDir);
    }

    Real calcSlipSpeedWitness(const State& s) const OVERRIDE_11 {
        if (isSticking(s)) return 0;
        const Vec2 vNow = m_contact.getSlipVelocity(s);
        if (!m_sliding->hasPrevSlipDir(s)) return vNow.norm();
        const Vec2 vPrev = m_sliding->getPrevSlipDir(s);
        return dot(vNow, vPrev);
    }

    Real calcStictionForceWitness(const State& s) const OVERRIDE_11 {
        if (!isSticking(s)) return 0;
        const Real mu_s = getStaticFrictionCoef();
        const Real N = getMasterNormalForce(s); // might be negative
        const Vec2 f = getStictionForce(s);
        const Real fmag = f.norm();
        return mu_s*N - fmag;
    }

    Real getActualSlipSpeed(const State& s) const OVERRIDE_11 {
        const Vec2 vNow = m_contact.getSlipVelocity(s); 
        return vNow.norm();
    }

    Real getActualFrictionForce(const State& s) const OVERRIDE_11 {
        const Real f = isSticking(s) ? getStictionForce(s).norm() 
                                     : m_sliding->calcSlidingForceMagnitude(s);
        return f;
    }

    Real getMasterNormalForce(const State& s) const OVERRIDE_11 {
        const Real N = m_contact.getForce(s); // might be negative
        return N;
    }

    Real getEstimatedNormalForce(const State& s) const OVERRIDE_11 {
        return m_sliding->getPrevN(s);
    }

    void setEstimatedNormalForce(State& s, Real N) const OVERRIDE_11 {
        m_sliding->setPrevN(s, N);
    }


    bool isMasterProximal(const State& s, Real posTol) const OVERRIDE_11
    {   return m_contact.isProximal(s, posTol); }
    bool isMasterCandidate(const State& s, Real posTol, Real velTol) const 
        OVERRIDE_11
    {   return m_contact.isCandidate(s, posTol, velTol); }
    bool isMasterActive(const State& s) const OVERRIDE_11
    {   return !m_contact.isDisabled(s); }


    // Set the friction application point to be the projection of the contact 
    // point onto the contact plane. This will be used the next time stiction
    // is enabled. Requires state realized to Position stage.
    void initializeForStiction(const State& s) OVERRIDE_11 {
        const Vec3 p = m_contact.getContactPointInPlaneBody(s);
        m_contactPointInPlane = p;
        m_tIc = m_tIe = m_tI = Vec2(0);
    }

    void recordImpulse(MyContactElement::ImpulseType type, const State& state,
                      const Vector& lambda) OVERRIDE_11
    {
        if (!isSticking(state)) return;

        // Record translational impulse.
        Vector myImpulseX(1), myImpulseZ(1);
        m_noslipX.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseX);
        m_noslipZ.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseZ);
        const Vec2 tI(myImpulseX[0], myImpulseZ[0]);
        if (type==MyContactElement::Compression) m_tIc = tI;
        else if (type==MyContactElement::Expansion) m_tIe = tI;
        m_tI += tI;
    }

    // Fill in deltaV to eliminate slip velocity using the stiction 
    // constraints.
    void setMyDesiredDeltaV(const State& s,
                            Vector& desiredDeltaV) const OVERRIDE_11
    {
        if (!isSticking(s)) return;

        const Vec2 dv = -m_contact.getSlipVelocity(s); // X,Z
        Vector myDesiredDV(1); // Nuke translational velocity also.
        myDesiredDV[0] = dv[0];
        m_noslipX.setMyPartInConstraintSpaceVector(s, myDesiredDV, desiredDeltaV);
        myDesiredDV[0] = dv[1];
        m_noslipZ.setMyPartInConstraintSpaceVector(s, myDesiredDV, desiredDeltaV);
    }

    Real getMyImpulseMagnitudeFromConstraintSpaceVector(const State& state,
                                                        const Vector& lambda) const
    {   Vector myImpulseX(1), myImpulseZ(1);
        m_noslipX.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseX);
        m_noslipZ.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseZ);
        const Vec2 tI(myImpulseX[0], myImpulseZ[0]);
        return tI.norm();
    }


    std::ostream& writeFrictionInfo(const State& s, const String& indent, 
                                    std::ostream& o) const OVERRIDE_11 
    {
        o << indent;
        if (!isMasterActive(s)) o << "OFF";
        else if (isSticking(s)) o << "STICK";
        else o << "SLIP";

        const Vec2 v = m_contact.getSlipVelocity(s);
        const Vec2 pd = m_sliding->getPrevSlipDir(s);
        const Vec2 f = isSticking(s) ? getStictionForce(s)
                                     : m_sliding->calcSlidingForce(s);
        o << " prevDir=" << ~pd << " V=" << ~v << " Vdot=" << ~v*pd 
          << " F=" << ~f << endl;
        return o;
    }


    void showFrictionForce(const State& s, Array_<DecorativeGeometry>& geometry) 
            const OVERRIDE_11
    {
        const Real Scale = 0.1;
        const Vec2 f = isSticking(s) ? getStictionForce(s)
                                     : m_sliding->calcSlidingForce(s);
        if (f.normSqr() < square(SignificantReal))
            return;
        const MobilizedBody& bodyB = m_contact.getBody();
        const Vec3& stationB = m_contact.getBodyStation();
        const Vec3 stationG = bodyB.getBodyTransform(s)*stationB;
        const Vec3 endG = stationG - Scale*Vec3(f[0], 0, f[1]);
        geometry.push_back(DecorativeLine(endG     + Vec3(0,.05,0),
                                          stationG + Vec3(0,.05,0))
                            .setColor(isSticking(s)?Green:Orange));
    }

    const MyPointContact& getMyPointContact() const {return m_contact;}
    const MySlidingFrictionForce& getMySlidingFrictionForce() const
    {   return *m_sliding; }
private:
    void initializeRuntimeFields() {
        m_contactPointInPlane = Vec3(0); 
        m_tIc = m_tIe = m_tI = Vec2(NaN);
        m_prevN = 0;
        m_prevSlipDir = Vec2(NaN);
    }
    const MyPointContact&   m_contact;

    Constraint::NoSlip1D    m_noslipX, m_noslipZ; // stiction
    MySlidingFrictionForce* m_sliding;  // sliding friction force element

    // Runtime
    Vec3 m_contactPointInPlane; // point on plane body where friction will act
    Vec2 m_tIc, m_tIe, m_tI; // impulses

    Real m_prevN;       // master's recorded normal force (could be negative)
    Vec2 m_prevSlipDir; // master's last recording slip or impending direction
};


//==============================================================================
//                            SHOW CONTACT
//==============================================================================
// For each visualization frame, generate ephemeral geometry to show just 
// during this frame, based on the current State.
class ShowContact : public DecorationGenerator {
public:
    ShowContact(const MyUnilateralConstraintSet& unis) 
    :   m_unis(unis) {}

    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        for (int i=0; i < m_unis.getNumContactElements(); ++i) {
            const MyContactElement& contact = m_unis.getContactElement(i);
            const Vec3 loc = contact.whereToDisplay(state);
            if (!contact.isDisabled(state)) {
                geometry.push_back(DecorativeSphere(.1)
                    .setTransform(loc)
                    .setColor(Red).setOpacity(.25));
                String text = "LOCKED";
                if (contact.hasFrictionElement()) {
                    const MyFrictionElement& friction = contact.getFrictionElement();
                    text = friction.isSticking(state) ? "STICKING"
                                                      : "CONTACT";
                    m_unis.getMultibodySystem().realize(state, Stage::Acceleration);
                    friction.showFrictionForce(state, geometry);
                }
                geometry.push_back(DecorativeText(String(i)+"-"+text)
                    .setColor(White).setScale(.1)
                    .setTransform(loc+Vec3(0,.04,0)));
            } else {
                geometry.push_back(DecorativeText(String(i))
                    .setColor(White).setScale(.1)
                    .setTransform(loc+Vec3(0,.02,0)));
            }
        }
    }
private:
    const MyUnilateralConstraintSet& m_unis;
};


//==============================================================================
//                            BODY WATCHER
//==============================================================================
// Prior to rendering each frame, point the camera at the given body's
// origin.
class BodyWatcher : public Visualizer::FrameController {
public:
    explicit BodyWatcher(const MobilizedBody& body) : m_body(body) {}

    void generateControls(const Visualizer&             viz, 
                          const State&                  state, 
                          Array_< DecorativeGeometry >& geometry) OVERRIDE_11
    {
        const Vec3 Bo = m_body.getBodyOriginLocation(state);
        const Vec3 p_GC = Bo + Vec3(0, 1, 5); // above and back
        const Rotation R_GC(UnitVec3(0,1,0), YAxis,
                            p_GC-Bo, ZAxis);
        viz.setCameraTransform(Transform(R_GC,p_GC));
        //viz.pointCameraAt(Bo, Vec3(0,1,0));
    }
private:
    const MobilizedBody m_body;
};

//==============================================================================
//                          STICTION ON HANDLER
//==============================================================================
// Allocate one of these for each contact constraint that has friction. This 
// handler takes care of turning on the stiction constraints when the sliding 
// velocity drops to zero.
class StictionOn: public TriggeredEventHandler {
public:
    StictionOn(const MultibodySystem&       system,
        const MyUnilateralConstraintSet&    unis,
        unsigned                            which) 
    :   TriggeredEventHandler(Stage::Velocity), 
        m_mbs(system), m_unis(unis), m_which(which)
    { 
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    // This is the witness function.
    Real getValue(const State& state) const {
        const MyFrictionElement& friction = m_unis.getFrictionElement(m_which);
        if (!friction.isMasterActive(state)) return 0;
        const Real signedSlipSpeed = friction.calcSlipSpeedWitness(state);
        return signedSlipSpeed;
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const 
    {
        SimTK_DEBUG2("\nhandle %d slide->stick@%.17g\n", m_which, s.getTime());
        SimTK_DEBUG("\n----------------------------------------------------\n");
        SimTK_DEBUG2("STICTION ON triggered by friction element %d @t=%.15g\n", 
            m_which, s.getTime());
        m_mbs.realize(s, Stage::Acceleration);

        #ifndef NDEBUG
        m_unis.showConstraintStatus(s, "ENTER STICTION ON");
        cout << " triggers=" << s.getEventTriggers() << "\n";
        #endif

        m_unis.selectActiveConstraints(s, accuracy);

        SimTK_DEBUG("STICTION ON done.\n");
        SimTK_DEBUG("----------------------------------------------------\n");
    }

private:
    const MultibodySystem&              m_mbs; 
    const MyUnilateralConstraintSet&    m_unis;
    const int                           m_which; // one of the friction elements
};



//==============================================================================
//                          STICTION OFF HANDLER
//==============================================================================
// Allocate one of these for each contact constraint. This handler takes
// care of turning off stiction constraints when the stiction force magnitude
// exceeds mu*N.
class StictionOff: public TriggeredEventHandler {
public:
    StictionOff(const MultibodySystem&      system,
        const MyUnilateralConstraintSet&    unis,
        unsigned                            which) 
    :   TriggeredEventHandler(Stage::Acceleration), 
        m_mbs(system), m_unis(unis), m_which(which)
    {
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    // This is the witness function. It is positive as long as mu_s*N is greater
    // than the friction force magnitude.
    Real getValue(const State& state) const {
        const MyFrictionElement& friction = m_unis.getFrictionElement(m_which);
        if (!friction.isMasterActive(state)) return 0;
        const Real forceMargin = friction.calcStictionForceWitness(state);
        return forceMargin; // how much stiction capacity is left
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const 
    {
        SimTK_DEBUG2("\nhandle %d stick->slide@%.17g\n", m_which, s.getTime());
        SimTK_DEBUG("\n----------------------------------------------------\n");
        SimTK_DEBUG2("STICTION OFF triggered by friction element %d @t=%.15g\n", 
            m_which, s.getTime());
        m_mbs.realize(s, Stage::Acceleration);

        #ifndef NDEBUG
        cout << " triggers=" << s.getEventTriggers() << "\n";
        #endif

        m_unis.selectActiveConstraints(s, accuracy);

        SimTK_DEBUG("STICTION OFF done.\n");
        SimTK_DEBUG("----------------------------------------------------\n");
    }

private:
    const MultibodySystem&              m_mbs; 
    const MyUnilateralConstraintSet&    m_unis;
    const int                           m_which; // one of the friction elements
};


//==============================================================================
//                            MY PUSH FORCE
//==============================================================================
// This is a force element that generates a constant force on a body for a
// specified time interval. It is useful to perturb the system to force it
// to transition from sticking to sliding, for example.
class MyPushForceImpl : public Force::Custom::Implementation {
public:
    MyPushForceImpl(const MobilizedBody& bodyB, 
                    const Vec3&          stationB,
                    const Vec3&          forceG, // force direction in Ground!
                    Real                 onTime,
                    Real                 offTime)
    :   m_bodyB(bodyB), m_stationB(stationB), m_forceG(forceG),
        m_on(onTime), m_off(offTime)
    {    }

    //--------------------------------------------------------------------------
    //                       Custom force virtuals
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   OVERRIDE_11
    {
        if (!(m_on <= state.getTime() && state.getTime() <= m_off))
            return;

        m_bodyB.applyForceToBodyPoint(state, m_stationB, m_forceG, bodyForces);

        //SimTK_DEBUG4("PUSHING @t=%g (%g,%g,%g)\n", state.getTime(),
        //    m_forceG[0], m_forceG[1], m_forceG[2]);
    }

    // No potential energy.
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {return 0;}

    void calcDecorativeGeometryAndAppend
       (const State& state, Stage stage, 
        Array_<DecorativeGeometry>& geometry) const OVERRIDE_11
    {
        const Real ScaleFactor = 0.1;
        if (stage != Stage::Time) return;
        if (!(m_on <= state.getTime() && state.getTime() <= m_off))
            return;
        const Vec3 stationG = m_bodyB.findStationLocationInGround(state, m_stationB);
        geometry.push_back(DecorativeLine(stationG-ScaleFactor*m_forceG, stationG)
                            .setColor(Yellow)
                            .setLineThickness(3));
    }
private:
    const MobilizedBody& m_bodyB; 
    const Vec3           m_stationB;
    const Vec3           m_forceG;
    Real                 m_on;
    Real                 m_off;
};



//==============================================================================
//                                   MAIN
//==============================================================================
int main(int argc, char** argv) {
  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem             mbs;
    SimbodyMatterSubsystem      matter(mbs);
    GeneralForceSubsystem       forces(mbs);
    Force::Gravity              gravity(forces, matter, -YAxis, 9.81);

    MobilizedBody& Ground = matter.updGround();


    // Predefine some handy rotations.
    const Rotation Z90(Pi/2, ZAxis); // rotate +90 deg about z

    const Vec3 BrickHalfDims(.1, .25, .5);
    const Real BrickMass = 10;
    #ifdef USE_TIMS_PARAMS
        const Real RunTime=16;  // Tim's time
        const Real Stiffness = 2e7;
        const Real Dissipation = 0.1;
        const Real CoefRest = 0; 
        // Painleve problem with these friction coefficients.
        //const Real mu_d = 1; /* compliant: .7*/
        //const Real mu_s = 1; /* compliant: .7*/
        const Real mu_d = .7;
        const Real mu_s = .7;
        const Real mu_v = /*0.05*/0;
        const Real CaptureVelocity = 0.01;
        const Real TransitionVelocity = 0.05;
        const Inertia brickInertia(.1,.1,.1);
    #else
        const Real RunTime=20;
        const Real Stiffness = 1e6;
        const Real CoefRest = 0.5; 
        const Real TargetVelocity = 3; // speed at which to match coef rest
//        const Real Dissipation = (1-CoefRest)/TargetVelocity;
        const Real Dissipation = 0.1;
        const Real mu_d = .5;
        const Real mu_s = .7;
        const Real mu_v = 0.05;
        const Real CaptureVelocity = 0.01;
        const Real TransitionVelocity = 0.05;
        const Inertia brickInertia(BrickMass*UnitInertia::brick(BrickHalfDims));
    #endif

    printf("\n******************** Tim's Box ********************\n");
    #ifdef USE_COMPLIANT_CONTACT
    printf("USING COMPLIANT CONTACT (with %s)\n",
        #ifdef USE_ELASTIC_FOUNDATION
        "elastic foundation mesh");
        #else
        "Hertz spheres at corners");
        #endif
    #else
    printf("USING RIGID CONTACT\n");
    #endif
    #ifdef USE_TIMS_PARAMS
    printf("Using Tim's parameters:\n");
    #else
    printf("Using Sherm's parameters:\n");
    #endif
    #ifdef USE_COMPLIANT_CONTACT
        printf("  stiffness=%g dissipation=%g\n", Stiffness, Dissipation);
    #else
        printf("  coef restitution=%g\n", CoefRest);
    #endif
    printf("  mu_d=%g mu_s=%g mu_v=%g\n", mu_d, mu_s, mu_v);
    printf("  transition velocity=%g\n", TransitionVelocity);
    printf("  brick inertia=%g %g %g\n",
        brickInertia.getMoments()[0], brickInertia.getMoments()[1], 
        brickInertia.getMoments()[2]); 
    printf("******************** Tim's Box ********************\n\n");

    #ifdef USE_COMPLIANT_CONTACT
    ContactTrackerSubsystem tracker(mbs);
    CompliantContactSubsystem contact(mbs, tracker);
    // Set stiction max slip velocity to make it less stiff.
    contact.setTransitionVelocity(0.05);

    const Real Radius = BrickHalfDims[0]/3;

    // Define a material to use for contact. This is not very stiff.
    ContactMaterial material(Stiffness,
                             Dissipation,
                             mu_s,  // mu_static
                             mu_d,  // mu_dynamic
                             mu_v); // mu_viscous

    // Make up an assumed thickness for the elastic layer in case we use
    // an ElasticFoundation-based compliant contact.
    const Real Thickness = min(BrickHalfDims);

    // Add a contact surface to represent the ground.
    // Half space normal is -x; must rotate about z to make it +y.
    Ground.updBody().addContactSurface(Rotation(-Pi/2,ZAxis),
       ContactSurface(ContactGeometry::HalfSpace(), material));
    #endif

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS


    MyUnilateralConstraintSet unis(mbs, CaptureVelocity);

    Body::Rigid brickBody = 
        Body::Rigid(MassProperties(BrickMass, Vec3(0), brickInertia));

    // First body: cube
    MobilizedBody::Cartesian loc(Ground, MassProperties(0,Vec3(0),Inertia(0)));
    MobilizedBody::Ball brick(loc, Vec3(0),
                             brickBody, Vec3(0));
    brick.addBodyDecoration(Transform(), DecorativeBrick(BrickHalfDims)
                                                .setColor(Red).setOpacity(.3));
/*
1) t= 0.5, dt = 2 sec, pt = (0.05, 0.2, 0.4), fdir = (1,0,0), mag = 50N
2) t= 4.0, dt = 0.5 sec, pt = (0.03, 0.06, 0.09), fdir = (0.2,0.8,0), mag = 300N
3) t= 0.9, dt = 2 sec, pt = (0,0,0), fdir = (0,1,0), mag = 49.333N (half the weight of the block)
4) t= 13.0, dt = 1 sec, pt = (0 0 0), fdir = (-1,0,0), mag = 200N
*/
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0.05,0.2,0.4),
                                                    50 * Vec3(1,0,0),
                                                    0.5, 0.5+2));
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0.03, 0.06, 0.09),
                                                    300 * UnitVec3(0.2,0.8,0),
                                                    //300 * Vec3(0.2,0.8,0),
                                                    4, 4+0.5));
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0),
                                                    49.333 * Vec3(0,1,0),
                                                    0.9, 0.9+2));
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0),
                                                    200 * Vec3(-1,0,0),
                                                    13, 13+1));

    #ifndef USE_TIMS_PARAMS
    // Extra late force.
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(.1, 0, .45),
                                                    20 * Vec3(-1,-1,.5),
                                                    15, Infinity));
    #endif

    #ifdef USE_COMPLIANT_CONTACT
    #ifndef USE_ELASTIC_FOUNDATION
    // Using spheres at the corners:
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(BrickHalfDims);
        const Vec3 shrink = 
            (abs(pt)-Radius).elementwiseDivide(abs(pt)).elementwiseMultiply(pt);
        brick.updBody().addContactSurface(shrink,
            ContactSurface(ContactGeometry::Sphere(Radius), material));
        brick.addBodyDecoration(shrink, DecorativeSphere(Radius));
    }
    #else
    // Using a brick mesh.
    //const PolygonalMesh mesh = 
    //    PolygonalMesh::createCylinderMesh(UnitVec3(1,1,1),BrickHalfDims[0],BrickHalfDims[2],1);
    //DecorativeMesh cylDec(cylMesh);
    //matter.updGround().addBodyDecoration(Vec3(-2,3,0), 
    //    cylDec
    //    );
    //matter.updGround().addBodyDecoration(Vec3(-2,3,0), 
    //    DecorativeMesh(cylMesh).setScale(1.02)
    //    .setRepresentation(DecorativeGeometry::DrawWireframe).setColor(Red)
    //    );

    const int rez = 2*(int)(max(BrickHalfDims)/min(BrickHalfDims) + 0.5);
    const PolygonalMesh mesh = 
        PolygonalMesh::createBrickMesh(BrickHalfDims,rez);

    const ContactGeometry::TriangleMesh triMesh(mesh);
    DecorativeMesh decMesh(triMesh.createPolygonalMesh());

    decMesh.setLineThickness(3)
           .setRepresentation(DecorativeGeometry::DrawWireframe)
           ;
    brick.updBody().addContactSurface(Vec3(0), 
        ContactSurface(triMesh,material,Thickness));
    brick.addBodyDecoration(Vec3(0), decMesh);
    #endif  
    #else
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(BrickHalfDims);
        MyPointContact* contact = new MyPointContact(brick, pt, CoefRest);
        unis.addContactElement(contact);
        unis.addFrictionElement(
            new MyPointContactFriction(*contact, mu_d, mu_s, mu_v, 
                                       CaptureVelocity, // TODO: vtol?
                                       forces));

    }
    #endif

    matter.setShowDefaultGeometry(false);
    Visualizer viz(mbs);
    //viz.setDesiredFrameRate(1000);
    viz.setShowSimTime(true);
    viz.setShowFrameNumber(true);
    viz.setShowFrameRate(true);
    viz.addDecorationGenerator(new ShowContact(unis));
#ifdef ANIMATE
    mbs.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));
#else
    // This does nothing but interrupt the simulation.
    mbs.addEventReporter(new Nada(ReportInterval));
#endif

    viz.addFrameController(new BodyWatcher(brick));

    Vec3 cameraPos(0, 1, 2);
    UnitVec3 cameraZ(0,0,1);
    viz.setCameraTransform(Transform(Rotation(cameraZ, ZAxis, UnitVec3(YAxis), YAxis), cameraPos));
    viz.pointCameraAt(Vec3(0,0,0), Vec3(0,1,0));

#ifdef USE_TIMS_PARAMS
    Real accuracy = 1e-4;
    RungeKuttaMersonIntegrator integ(mbs);
#else
    Real accuracy = 1e-2;
    RungeKutta3Integrator integ(mbs);
#endif
    integ.setAccuracy(accuracy);
    //integ.setMaximumStepSize(0.25);
    integ.setMaximumStepSize(0.1);
    //integ.setMaximumStepSize(0.001);

    StateSaver* stateSaver = new StateSaver(mbs,unis,integ,ReportInterval);
    mbs.addEventReporter(stateSaver);

    #ifndef USE_COMPLIANT_CONTACT    
    for (int i=0; i < unis.getNumContactElements(); ++i) {
        mbs.addEventHandler(new ContactOn(mbs, unis,i, Stage::Position));
        mbs.addEventHandler(new ContactOn(mbs, unis,i, Stage::Velocity));
        mbs.addEventHandler(new ContactOn(mbs, unis,i, Stage::Acceleration));
        mbs.addEventHandler(new ContactOff(mbs, unis,i));
    }

    for (int i=0; i < unis.getNumFrictionElements(); ++i) {
        mbs.addEventHandler(new StictionOn(mbs, unis, i));
        mbs.addEventHandler(new StictionOff(mbs, unis, i));
    }
    #endif

    State s = mbs.realizeTopology(); // returns a reference to the the default state
    
    //matter.setUseEulerAngles(s, true);
    
    mbs.realizeModel(s); // define appropriate states for this System
    mbs.realize(s, Stage::Instance); // instantiate constraints if any


/*
rX_q = 0.7 rad
rX_u = 1.0 rad/sec

rY_q = 0.6 rad
rY_u = 0.0 rad/sec

rZ_q = 0.5 rad
rZ_u = 0.2 rad/sec

tX_q = 0.0 m
tX_u = 10 m/s

tY_q = 1.4 m
tY_u = 0.0 m/s

tZ_q = 0.0 m
tZ_u = 0.0 m/s
*/

#ifdef USE_TIMS_PARAMS
    loc.setQToFitTranslation(s, Vec3(0,10,0));
    loc.setUToFitLinearVelocity(s, Vec3(0,0,0));
#else
    loc.setQToFitTranslation(s, Vec3(0,1.4,0));
    loc.setUToFitLinearVelocity(s, Vec3(10,0,0));
    const Rotation R_BC(SimTK::BodyRotationSequence,
                                0.7, XAxis, 0.6, YAxis, 0.5, ZAxis);
    brick.setQToFitRotation(s, R_BC);
    brick.setUToFitAngularVelocity(s, Vec3(1,0,.2));
#endif

    mbs.realize(s, Stage::Velocity);
    viz.report(s);

    Array_<int> enableTheseContacts;
    for (int i=0; i < unis.getNumContactElements(); ++i) {
        const Real perr = unis.getContactElement(i).getPerr(s);
        printf("contact constraint %d has perr=%g%s\n", i, perr,
            perr<=0?" (ENABLING CONTACT)":"");
        if (perr <= 0)
            enableTheseContacts.push_back(i);
    }

    cout << "Initial configuration shown. Next? ";
    getchar();

    for (unsigned i=0; i < enableTheseContacts.size(); ++i)
        unis.getContactElement(enableTheseContacts[i]).enable(s);

    Assembler(mbs).setErrorTolerance(1e-6).assemble(s);
    viz.report(s);
    cout << "Assembled configuration shown. Ready? ";
    getchar();

    // Now look for enabled contacts that aren't sliding; turn on stiction
    // for those.
    mbs.realize(s, Stage::Velocity);
    Array_<int> enableTheseStictions;
    for (int i=0; i < unis.getNumFrictionElements(); ++i) {
        break; //TODO
        MyFrictionElement& fric = unis.updFrictionElement(i);
        const Real vSlip = fric.getActualSlipSpeed(s);
        fric.initializeForStiction(s); // just in case
        printf("friction element %d has v_slip=%g%s\n", i, vSlip,
            vSlip==0?" (ENABLING STICTION)":"");
        if (vSlip == 0)
            enableTheseStictions.push_back(i);
    }
    for (unsigned i=0; i < enableTheseStictions.size(); ++i)
        unis.getFrictionElement(enableTheseStictions[i]).enableStiction(s);

    // Make sure Lapack gets initialized.
    Matrix M(1,1); M(0,0)=1.23;
    FactorLU Mlu(M);

    
    // Simulate it.

    integ.setReturnEveryInternalStep(true);
    //integ.setAllowInterpolation(false);
    TimeStepper ts(mbs, integ);
    ts.setReportAllSignificantStates(true);

    #ifdef TEST_REPEATABILITY
        const int tries = NTries;
    #else
        const int tries = 1;
    #endif
        
    Array_< Array_<State> > states(tries);
    Array_< Array_<Real> > timeDiff(tries-1);
    Array_< Array_<Vector> > yDiff(tries-1);
    cout.precision(18);
    for (int ntry=0; ntry < tries; ++ntry) {
        mbs.resetAllCountersToZero();
        unis.initialize(); // reinitialize
        ts.updIntegrator().resetAllStatistics();
        ts.initialize(s);
        int nStepsWithEvent = 0;

        const double startReal = realTime();
        const double startCPU = cpuTime();

        Integrator::SuccessfulStepStatus status;
        do {
            status=ts.stepTo(RunTime);
            #ifdef TEST_REPEATABILITY
                states[ntry].push_back(ts.getState());
            #endif 
            const int j = states[ntry].size()-1;
            if (ntry>0) {
                int prev = ntry-1;
                //prev=0;
                Real dt = states[ntry][j].getTime() 
                          - states[prev][j].getTime();
                volatile double vd;
                const int ny = states[ntry][j].getNY();
                Vector d(ny);
                for (int i=0; i<ny; ++i) {
                    vd = states[ntry][j].getY()[i] 
                           - states[prev][j].getY()[i];
                    d[i] = vd;
                }
                timeDiff[prev].push_back(dt);
                yDiff[prev].push_back(d);
                cout << "\n" << states[prev][j].getTime() 
                     << "+" << dt << ":" << d << endl;
            }

            const bool handledEvent = status==Integrator::ReachedEventTrigger;
            if (handledEvent) {
                ++nStepsWithEvent;
                SimTK_DEBUG3("EVENT %3d @%.17g status=%s\n\n", 
                    nStepsWithEvent, ts.getTime(),
                    Integrator::getSuccessfulStepStatusString(status).c_str());
            } else {
                SimTK_DEBUG3("step  %3d @%.17g status=%s\n", integ.getNumStepsTaken(),
                    ts.getTime(),
                    Integrator::getSuccessfulStepStatusString(status).c_str());
            }
            #ifndef NDEBUG
                viz.report(ts.getState());
            #endif
        } while (ts.getTime() < RunTime);


        const double timeInSec = realTime()-startReal;
        const double cpuInSec = cpuTime()-startCPU;
        const int evals = integ.getNumRealizations();
        cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
            timeInSec << "s for " << ts.getTime() << "s sim (avg step=" 
            << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) " 
            << (1000*ts.getTime())/evals << "ms/eval\n";
        cout << "CPUtime " << cpuInSec << endl;

        printf("Used Integrator %s at accuracy %g:\n", 
            integ.getMethodName(), integ.getAccuracyInUse());
        printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
        printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
        printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());
        printf("# EVENT STEPS/HANDLER CALLS = %d/%d\n", 
            nStepsWithEvent, mbs.getNumHandleEventCalls());
    }

    for (int i=0; i<tries; ++i)
        cout << "nstates " << i << " " << states[i].size() << endl;

    // Instant replay.
    while(true) {
        printf("Hit ENTER for replay (%d states) ...", 
                stateSaver->getNumSavedStates());
        getchar();
        for (int i=0; i < stateSaver->getNumSavedStates(); ++i) {
            mbs.realize(stateSaver->getState(i), Stage::Velocity);
            viz.report(stateSaver->getState(i));
        }
    }

  } 
  catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
  catch (...) {
    printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

}


//==============================================================================
//                        IMPACT HANDLING (CONTACT ON)
//==============================================================================

//------------------------------ HANDLE EVENT ----------------------------------
// There are three different witness functions that cause this handler to get
// invoked. The position- and velocity-level ones require an impulse. The
// acceleration-level one just requires recalculating the active set, in the
// same manner as liftoff or friction transition events.

void ContactOn::
handleEvent(State& s, Real accuracy, bool& shouldTerminate) const 
{
    SimTK_DEBUG3("\nhandle %d impact@%.17g (%s)\n", m_which, s.getTime(),
         m_stage.getName().c_str());

    if (m_stage == Stage::Acceleration) {
        SimTK_DEBUG("\n---------------CONTACT ON (ACCEL)--------------\n");
        SimTK_DEBUG2("CONTACT triggered by element %d @t=%.15g\n", 
            m_which, s.getTime());
        m_mbs.realize(s, Stage::Acceleration);

        #ifndef NDEBUG
        cout << " triggers=" << s.getEventTriggers() << "\n";
        #endif

        m_unis.selectActiveConstraints(s, accuracy);
        SimTK_DEBUG("---------------CONTACT ON (ACCEL) done.--------------\n");
        return;
    }

    MyElementSubset proximal;
    m_unis.findProximalElements(s, accuracy, proximal);

    // Zero out accumulated impulses and perform any other necessary 
    // initialization of contact and friction elements.
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        m_unis.updContactElement(which)
            .initializeForImpact(s, m_unis.getCaptureVelocity());
    }
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        m_unis.updFrictionElement(which).initializeForStiction(s);
    }

    SimTK_DEBUG("\n---------------------CONTACT ON---------------------\n");
    SimTK_DEBUG3("\nIMPACT (%s) for contact %d at t=%.16g\n", 
        m_stage.getName().c_str(), m_which, s.getTime());
    SimTK_DEBUG2("  %d/%d proximal contact/friction elements\n", 
        proximal.m_contact.size(), proximal.m_friction.size());

    m_unis.showConstraintStatus(s, "ENTER IMPACT (CONTACT ON)");

    bool needMoreCompression = true;
    while (needMoreCompression) {
        processCompressionPhase(proximal, s);
        needMoreCompression = false;
        if (processExpansionPhase(proximal, s)) {
            for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
                const int which = proximal.m_contact[i];
                const MyContactElement& contact = 
                    m_unis.getContactElement(which);
                if (contact.getVerr(s) < 0) {
                    needMoreCompression = true;
                    break;
                }
            }
        }
    }

    // Record new previous slip velocities for all the sliding friction
    // since velocities have changed. First loop collects the velocities.
    m_mbs.realize(s, Stage::Velocity);
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        MyFrictionElement& fric = m_unis.updFrictionElement(which);
        if (!fric.isMasterActive(s) || fric.isSticking(s)) continue;
        fric.recordSlipDir(s);
    }

    // Now update all the previous slip direction state variables from the
    // recorded values.
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        const MyFrictionElement& fric = m_unis.getFrictionElement(which);
        if (!fric.isMasterActive(s) || fric.isSticking(s)) continue;
        fric.updatePreviousSlipDirFromRecorded(s);
    }

    m_unis.selectActiveConstraints(s, accuracy);

    SimTK_DEBUG("\n-----------------END CONTACT ON---------------------\n");
}




//------------------------ PROCESS COMPRESSION PHASE ---------------------------
// Given a list of proximal unilateral constraints (contact and stiction),
// determine which ones are active in the least squares solution for the
// constraint impulses. Candidates are those constraints that meet the 
// kinematic proximity condition -- for contacts, position less than
// tolerance; for stiction, master contact is proximal. Also, any
// constraint that is currently active is a candidate, regardless of its
// kinematics.
//
// TODO: NOT ENABLING STICTION AT ALL; assuming stiction on initially leads
// to a lot of wasted time and solution difficulty. Must start with sliding.
// TODO: sliding friction impulses
//
// Algorithm
// ---------
// We're given a set of candidates including contacts and stiction. If any are
// inactive, activate them.
// -- at this point all verr==0, some impulses f might be < 0
//
// loop
// - Calculate impulses with the current active set
//     (iterate sliding impulses until f=mu_d*N to tol, where N is normal 
//      impulse)
// - Calculate f for active constraints, verr for inactive
// - If all f>=0, verr>=0 -> break loop
// - Check for verr < 0 [tol?]. Shouldn't happen but if it does must turn on the
//     associated constraint for the worst violation, then -> continue loop
// - Find worst (most negative) offender:
//    contact offense  = fc < 0 ? fc : 0
//    stiction offense = mu_d*max(0, fc) - |fs|
// - Choose constraint to deactivate:
//     worst is a stiction constraint: choose it
//     worst is a contact constraint: if it has stiction, chose that
//                                    otherwise choose the contact constraint
// - Inactivate chosen constraint
//     (if stiction, record impending slip direction & N for stick->slide)
// end loop 
//
void ContactOn::
processCompressionPhase(MyElementSubset&    proximal,
                        State&              s) const
{
    const int ReviseNormalNIters = 6;

    SimTK_DEBUG2("Entering processCompressionPhase(): "
        "%d/%d impact/stick candidates ...\n", proximal.m_contact.size(),
        proximal.m_friction.size());

    if (proximal.m_contact.empty()) {
        // Can't be any friction either, if there are no contacts.
        SimTK_DEBUG("EXIT processCompressionPhase: no proximal candidates.\n");
        return;
    }

    Vector lambda;
    const Vector u0 = s.getU(); // save presenting velocity

    // Assume at first that all proximal contacts will participate. This is 
    // necessary to ensure that we get a least squares solution for the impulse 
    // involving as many constraints as possible sharing the impulse. 
    // TODO: note stiction is unconditionally disabled
    m_unis.enableConstraintSubset(proximal, false, s); 

    int pass=0, nContactsDisabled=0, nStictionDisabled=0, nContactsReenabled=0;
    while (true) {
        ++pass; 
        SimTK_DEBUG1("processCompressionPhase(): pass %d\n", pass);

        // Given an active set, evaluate impulse multipliers & forces, and
        // evaluate resulting constraint velocity errors.
        calcStoppingImpulse(proximal, s, lambda);
        // TODO: ignoring sliding impacts; should converge here.
        updateVelocities(u0, lambda, s);

        // Scan all proximal contacts to find the active one that has the
        // most negative normal force, and the inactive one that has the 
        // most negative velocity error (hopefully none will).

        int worstActiveContact=-1; Real mostNegativeContactImpulse=0;
        int worstInactiveContact=-1; Real mostNegativeVerr=0;
        
        SimTK_DEBUG("analyzing impact constraints ...\n");
        for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
            const int which = proximal.m_contact[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyContactElement& cont = m_unis.getContactElement(which);
            if (cont.isDisabled(s)) {
                const Real verr = cont.getVerr(s);
                SimTK_DEBUG1("off verr=%g\n", verr);
                if (verr < mostNegativeVerr) {
                    worstInactiveContact = which;
                    mostNegativeVerr = verr;
                }
            } else {
                const Real f = 
                    cont.getMyValueFromConstraintSpaceVector(s, lambda);
                SimTK_DEBUG1("on impulse=%g\n", f);
                if (f < mostNegativeContactImpulse) {
                    worstActiveContact = which;
                    mostNegativeContactImpulse = f;
                }
            }
        }

        // This is bad and might cause cycling.
        if (mostNegativeVerr < 0) {
            SimTK_DEBUG2("  !!! Inactive contact %d violated, verr=%g\n", 
                worstInactiveContact, mostNegativeVerr);
            const MyContactElement& cont = 
                m_unis.getContactElement(worstInactiveContact);
            //TODO -- must use a tolerance or prevent looping
            //++nContactsReenabled;
            //cont.enable(s);
            //continue;
        }

        SimTK_DEBUG("  All inactive contacts are satisfied.\n");

        #ifndef NDEBUG
        if (mostNegativeContactImpulse == 0)
            printf("  All active contacts are satisfied.\n");
        else 
            printf("  Active contact %d was worst violator with f=%g\n",
                worstActiveContact, mostNegativeContactImpulse);
        #endif

        int worstActiveStiction=-1; Real mostNegativeStictionCapacity=0;     
        SimTK_DEBUG("analyzing stiction constraints ...\n");
        for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
            const int which = proximal.m_friction[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyFrictionElement& fric = m_unis.getFrictionElement(which);
            if (!fric.isSticking(s)) {
                SimTK_DEBUG("off\n");
                continue;
            }
            // TODO: Kludge -- must be point contact.
            const MyPointContactFriction& ptfric = 
                dynamic_cast<const MyPointContactFriction&>(fric);
            const MyPointContact& cont = ptfric.getMyPointContact();
            const Real N = cont.getMyValueFromConstraintSpaceVector(s, lambda);
            const Real fsmag = 
                ptfric.getMyImpulseMagnitudeFromConstraintSpaceVector(s, lambda);
            const Real mu_d = fric.getDynamicFrictionCoef();
            const Real capacity = mu_d*std::max(N,Real(0)) - fsmag;
            SimTK_DEBUG2("on capacity=%g (N=%g)\n", capacity, N);

            if (capacity < mostNegativeStictionCapacity) {
                worstActiveStiction = which;
                mostNegativeStictionCapacity = capacity;
            }
        }

        #ifndef NDEBUG
        if (mostNegativeStictionCapacity == 0)
            printf("  All active stiction constraints are satisfied.\n");
        else 
            printf("  Active stiction %d was worst violator with capacity=%g\n",
                worstActiveStiction, mostNegativeStictionCapacity);
        #endif

        if (mostNegativeContactImpulse==0 && mostNegativeStictionCapacity==0) {
            SimTK_DEBUG("DONE. Current active set is a winner.\n");
            break;
        }

        // Restore original velocity.
        s.updU() = u0;

        if (mostNegativeStictionCapacity <= mostNegativeContactImpulse) {
            SimTK_DEBUG1("  Disable stiction %d\n", worstActiveStiction);
            MyFrictionElement& fric = 
                m_unis.updFrictionElement(worstActiveStiction);
            // TODO: need the impulse version of this
            //fric.recordImpendingSlipInfo(s);
            ++nStictionDisabled;
            fric.disableStiction(s);
            continue;
        }

        // A contact constraint was the worst violator. If that contact
        // constraint has an active stiction constraint, we have to disable
        // the stiction constraint first.
        SimTK_DEBUG1("  Contact %d was the worst violator.\n", 
            worstActiveContact);
        const MyContactElement& cont = 
            m_unis.getContactElement(worstActiveContact);
        assert(!cont.isDisabled(s));

        if (cont.hasFrictionElement()) {
            MyFrictionElement& fric = cont.updFrictionElement();
            if (fric.isSticking(s)) {
                SimTK_DEBUG1("  ... but must disable stiction %d first.\n",
                    fric.getFrictionIndex());
                // TODO: need the impulse version of this
                //fric.recordImpendingSlipInfo(s);
                ++nStictionDisabled;
                fric.disableStiction(s);
                continue;
            }
        }

        SimTK_DEBUG1("  Disable contact %d\n", worstActiveContact); 
        ++nContactsDisabled;
        cont.disable(s);
        // Go back for another pass.
    }


    // Now update the entries for each proximal constraint to reflect the
    // compression impulse and post-compression velocity.
    SimTK_DEBUG("  Compression results:\n");
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        MyContactElement& cont = m_unis.updContactElement(which);
        if (!cont.isDisabled(s))
            cont.recordImpulse(MyContactElement::Compression, s, lambda);
        SimTK_DEBUG4("  %d %3s: Ic=%g, V=%g\n",
            which, cont.isDisabled(s) ? "off" : "ON", 
            cont.getCompressionImpulse(), cont.getVerr(s));
    }

    SimTK_DEBUG("... compression phase done.\n");
}

//------------------------- PROCESS EXPANSION PHASE ----------------------------
bool ContactOn::
processExpansionPhase(MyElementSubset&  proximal,
                      State&            s) const
{
    SimTK_DEBUG("Entering processExpansionPhase() ...\n");

    // Generate an expansion impulse if there were any active contacts that
    // still have some restitution remaining.
    Vector expansionImpulse;

    bool anyChange = false;
    for (unsigned i=0; i<proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        MyContactElement& uni = m_unis.updContactElement(which);
        if (uni.isDisabled(s)||uni.isRestitutionDone()
            ||uni.getEffectiveCoefRest()==0
            ||uni.getCompressionImpulse()<=0)
            continue;
        uni.setMyExpansionImpulse(s, uni.getEffectiveCoefRest(), 
                                  expansionImpulse);
        uni.recordImpulse(MyContactElement::Expansion,s,expansionImpulse);
        uni.setRestitutionDone(true);
        anyChange = true;
    }

    if (!anyChange) {
        SimTK_DEBUG("... no expansion impulse -- done.\n");
        return false;
    }

    // We generated an expansion impulse. Apply it and update velocities.
    updateVelocities(Vector(), expansionImpulse, s);

    // Release any constraint that now has a positive velocity.
    Array_<int> toDisable;
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        if (!uni.isDisabled(s) && uni.getVerr(s) > 0)
            toDisable.push_back(which);
    }

    // Now do the actual disabling (can't mix this with checking velocities)
    // because disabling invalidates Instance stage.
    for (unsigned i=0; i < toDisable.size(); ++i) {
        const int which = toDisable[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        uni.disable(s);
    }

    SimTK_DEBUG("  Expansion results:\n");
    m_mbs.realize(s, Stage::Velocity);
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        SimTK_DEBUG4("  %d %3s: Ie=%g, V=%g\n",
            which, uni.isDisabled(s) ? "off" : "ON", 
            uni.getExpansionImpulse(), uni.getVerr(s));
    }

    SimTK_DEBUG("... expansion phase done.\n");

    return true;
}




//-------------------------- CALC STOPPING IMPULSE -----------------------------
// Calculate the impulse that eliminates all residual velocity for the
// current set of enabled constraints.
// Note: if you have applied impulses also (like sliding friction), 
// convert to generalized impulse f, then to desired delta V in constraint
// space like this: deltaV = G*M\f; add that to the verrs to get the total
// velocity change that must be produced by the impulse.
void ContactOn::
calcStoppingImpulse(const MyElementSubset&  proximal,
                    const State&            s,
                    Vector&                 lambda0) const
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(s, Stage::Dynamics); // TODO: should only need Position
    Vector desiredDeltaV;  // in constraint space
    SimTK_DEBUG("  Entering calcStoppingImpulse() ...\n");
    bool gotOne = false;
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        if (uni.isDisabled(s))
            continue;
        SimTK_DEBUG2("    uni constraint %d enabled, v=%g\n",
            which, uni.getVerr(s));
        uni.setMyDesiredDeltaV(s, desiredDeltaV);
        gotOne = true;
    }
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        const MyFrictionElement& fric = m_unis.getFrictionElement(which);
        if (!fric.isSticking(s))
            continue;
        SimTK_DEBUG2("    friction constraint %d enabled, |v|=%g\n",
            which, fric.getActualSlipSpeed(s));
        fric.setMyDesiredDeltaV(s, desiredDeltaV);
        gotOne = true;
    }

    if (gotOne) matter.solveForConstraintImpulses(s, desiredDeltaV, lambda0);
    else lambda0.clear();
#ifndef NDEBUG
    cout << "  ... done. Stopping impulse=" << lambda0 << endl;
#endif
}



//---------------------------- UPDATE VELOCITIES -------------------------------
void ContactOn::
updateVelocities(const Vector& u0, const Vector& lambda, State& state) const {
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    Vector f, deltaU;
    assert(u0.size()==0 || u0.size() == state.getNU());

    m_mbs.realize(state, Stage::Dynamics); // TODO: Position
    matter.multiplyByGTranspose(state,lambda,f);
    matter.multiplyByMInv(state,f,deltaU);
    if (u0.size()) state.updU() = u0 + deltaU;
    else state.updU() += deltaU;
    m_mbs.realize(state, Stage::Velocity);
}



//==============================================================================
//                       MY UNILATERAL CONSTRAINT SET
//==============================================================================


//------------------------ SELECT ACTIVE CONSTRAINTS ---------------------------
void MyUnilateralConstraintSet::
selectActiveConstraints(State& state, Real accuracy) const {

    // Find all the contacts and stiction elements that might be active based
    // on kinematics.
    MyElementSubset candidates;

    bool needRestart;
    do {
        //TODO: this (mis)use of accuracy needs to be revisited.
        findCandidateElements(state, accuracy, accuracy, candidates);

        // Evaluate accelerations and reaction forces and check if 
        // any of the active contacts are generating negative ("pulling") 
        // forces; if so, inactivate them.
        findActiveCandidates(state, candidates);

        Array_<Real> slipVels(candidates.m_friction.size());
        for (unsigned i=0; i < candidates.m_friction.size(); ++i) {
            const int which = candidates.m_friction[i];
            const MyFrictionElement& fric = getFrictionElement(which);
            slipVels[i] = fric.getActualSlipSpeed(state);
        }

        // Finally, project active constraints to the constraint manifolds.
        const Real tol = accuracy/1000;
        SimTK_DEBUG1("projecting active constraints to tol=%g\n", tol);
        m_mbs.project(state, tol);

        // It is possible that the projection above took some marginally-sliding
        // friction and slowed it down enough to make it a stiction candidate.
        needRestart = false;
        for (unsigned i=0; i < candidates.m_sliding.size(); ++i) {
            const int which = candidates.m_sliding[i];
            const MyFrictionElement& fric = getFrictionElement(which);
            if (   fric.getActualSlipSpeed(state) <= accuracy
                || fric.calcSlipSpeedWitness(state) <= 0) 
            {
                SimTK_DEBUG3("***RESTART1** selectActiveConstraints() because "
                    "sliding velocity of friction %d is |v|=%g or witness=%g\n",
                    which, fric.getActualSlipSpeed(state),
                    fric.calcSlipSpeedWitness(state));
                needRestart = true;
                break;
            }
        }
        if (needRestart) continue;
        for (unsigned i=0; i < candidates.m_friction.size(); ++i) {
            const int which = candidates.m_friction[i];
            const MyFrictionElement& fric = getFrictionElement(which);
            if (fric.isSticking(state)) continue;
            if (slipVels[i] > accuracy
                && fric.getActualSlipSpeed(state) <= accuracy)
            {
                SimTK_DEBUG3("***RESTART2** selectActiveConstraints() because "
                    "sliding velocity of friction %d went from |v|=%g to |v|=%g\n",
                    which, slipVels[i], fric.getActualSlipSpeed(state));
                needRestart = true;
                break;
            }
        }

    } while (needRestart);
}




//-------------------------- FIND ACTIVE CANDIDATES ---------------------------
// Given a list of candidate unilateral constraints (contact and stiction),
// determine which ones are active in the least squares solution for the
// constraint multipliers. Candidates are those constraints that meet all 
// kinematic conditions -- for contacts, position and velocity errors less than
// tolerance; for stiction, sliding velocity less than tolerance. Also, any
// constraint that is currently active is a candidate, regardless of its
// kinematics.
//
// This method should be called only from within an event handler. For sliding
// friction it will have reset the "previous slip direction" to the current
// slip or impending slip direction, and converged the remembered normal force.
//
// Algorithm
// ---------
// We're given a set of candidates including contacts and stiction. If any are
// inactive, activate them.
// -- at this point all aerr==0, some ferr might be < 0
//
// loop
// - Realize(Acceleration) with the current active set
//     (iterate sliding forces until f=mu_d*N to tol)
// - Calculate ferr for active constraints, aerr for inactive
// - If all ferr>=0, aerr>=0 -> break loop
// - Check for aerr < 0 [tol?]. Shouldn't happen but if it does must turn on the
//     associated constraint for the worst violation, then -> continue loop
// - Find worst (most negative) offender:
//    contact offense  = fc < 0 ? fc : 0
//    stiction offense = mu_s*max(0, fc) - |fs|
// - Choose constraint to deactivate:
//     worst is a stiction constraint: choose it
//     worst is a contact constraint: if it has stiction, chose that
//                                    otherwise choose the contact constraint
// - Inactivate chosen constraint
//     (if stiction, record impending slip direction & N for stick->slide)
// end loop 
//
void MyUnilateralConstraintSet::
findActiveCandidates(State& s, const MyElementSubset& candidates) const
{
    const int ReviseNormalNIters = 6;
    showConstraintStatus(s, "ENTER findActiveCandidates()");
    if (candidates.m_contact.empty()) {
        // Can't be any friction either, if there are no contacts.
        SimTK_DEBUG("EXIT findActiveCandidates: no candidates.\n");
        m_mbs.realize(s, Stage::Acceleration);
        return;
    }

    SimTK_DEBUG3(
        "findActiveCandidates() for %d/%d/%d contact/stick/slip candidates ...\n",
        candidates.m_contact.size(), candidates.m_friction.size(),
        candidates.m_sliding.size());

    // Enable all candidate contact and stiction constraints if any are
    // currently disabled.
    enableConstraintSubset(candidates, true, s);

    int pass=0, nContactsDisabled=0, nStictionDisabled=0, nContactsReenabled=0;
    while (true) {
        ++pass; 
        SimTK_DEBUG1("\nfindActiveCandidates(): pass %d\n", pass);

        // Given an active set, evaluate multipliers and accelerations, and
        // converge sliding forces.
        m_mbs.realize(s, Stage::Acceleration);
        SimTK_DEBUG("findActiveCandidates(): CONVERGE NORMALS\n");
        convergeNormalForcesForSlidingContacts(s, SqrtEps);
 /*       for (int i=0; i < ReviseNormalNIters; ++i) {
            s.autoUpdateDiscreteVariables();
            s.invalidateAllCacheAtOrAbove(Stage::Dynamics);
            convergeNormalForcesForSlidingContacts(s, SqrtEps);
            m_mbs.realize(s, Stage::Acceleration);
        }*/
        m_mbs.realize(s, Stage::Acceleration);

        // Scan all candidate contacts to find the active one that has the
        // most negative normal force, and the inactive one that has the 
        // most negative acceleration error (hopefully none will).

        int worstActiveContact=-1; Real mostNegativeContactForce=0;
        int worstInactiveContact=-1; Real mostNegativeAerr=0;
        
        SimTK_DEBUG("analyzing contact constraints ...\n");
        for (unsigned i=0; i < candidates.m_contact.size(); ++i) {
            const int which = candidates.m_contact[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyContactElement& cont = getContactElement(which);
            if (cont.isDisabled(s)) {
                const Real aerr = cont.getAerr(s);
                SimTK_DEBUG1("off aerr=%g\n", aerr);
                if (aerr < mostNegativeAerr) {
                    worstInactiveContact = which;
                    mostNegativeAerr = aerr;
                }
            } else {
                const Real f = cont.getForce(s);
                SimTK_DEBUG1("on f=%g\n", f);
                if (f < mostNegativeContactForce) {
                    worstActiveContact = which;
                    mostNegativeContactForce = f;
                }
            }
        }

        // This is bad and might cause cycling.
        if (mostNegativeAerr < 0) {
            SimTK_DEBUG2("  !!! Inactive contact %d violated, aerr=%g\n", 
                worstInactiveContact, mostNegativeAerr);
            const MyContactElement& cont = getContactElement(worstInactiveContact);
            //TODO -- must use a tolerance or prevent looping
            //++nContactsReenabled;
            //cont.enable(s);
            //continue;
        }

        SimTK_DEBUG("  All inactive contacts are satisfied.\n");

        #ifndef NDEBUG
        if (mostNegativeContactForce == 0)
            printf("  All active contacts are satisfied.\n");
        else 
            printf("  Active contact %d was worst violator with f=%g\n",
                worstActiveContact, mostNegativeContactForce);
        #endif

        int worstActiveStiction=-1; Real mostNegativeStictionCapacity=0;     
        SimTK_DEBUG("analyzing stiction constraints ...\n");
        for (unsigned i=0; i < candidates.m_friction.size(); ++i) {
            const int which = candidates.m_friction[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyFrictionElement& fric = getFrictionElement(which);
            if (!fric.isSticking(s)) {
                SimTK_DEBUG("off\n");
                continue;
            }
            const Real mu_s = fric.getStaticFrictionCoef();
            const Real N = fric.getMasterNormalForce(s); // might be negative
            const Real fsmag = fric.getActualFrictionForce(s);
            const Real capacity = mu_s*std::max(N,Real(0)) - fsmag;
            SimTK_DEBUG2("on capacity=%g (N=%g)\n", capacity, N);

            if (capacity < mostNegativeStictionCapacity) {
                worstActiveStiction = which;
                mostNegativeStictionCapacity = capacity;
            }
        }

        #ifndef NDEBUG
        if (mostNegativeStictionCapacity == 0)
            printf("  All active stiction constraints are satisfied.\n");
        else 
            printf("  Active stiction %d was worst violator with capacity=%g\n",
                worstActiveStiction, mostNegativeStictionCapacity);
        #endif

        if (mostNegativeContactForce==0 && mostNegativeStictionCapacity==0) {
            SimTK_DEBUG("DONE. Current active set is a winner.\n");
            break;
        }

        if (mostNegativeStictionCapacity <= mostNegativeContactForce) {
            SimTK_DEBUG1("  Disable stiction %d\n", worstActiveStiction);
            MyFrictionElement& fric = updFrictionElement(worstActiveStiction);
            fric.recordImpendingSlipInfo(s);
            ++nStictionDisabled;
            fric.disableStiction(s);
            continue;
        }

        // A contact constraint was the worst violator. If that contact
        // constraint has an active stiction constraint, we have to disable
        // the stiction constraint first.
        SimTK_DEBUG1("  Contact %d was the worst violator.\n", worstActiveContact);
        const MyContactElement& cont = getContactElement(worstActiveContact);
        assert(!cont.isDisabled(s));

        if (cont.hasFrictionElement()) {
            MyFrictionElement& fric = cont.updFrictionElement();
            if (fric.isSticking(s)) {
                SimTK_DEBUG1("  ... but must disable stiction %d first.\n",
                    fric.getFrictionIndex());
                fric.recordImpendingSlipInfo(s);
                ++nStictionDisabled;
                fric.disableStiction(s);
                continue;
            }
        }

        SimTK_DEBUG1("  Disable contact %d\n", worstActiveContact); 
        ++nContactsDisabled;
        cont.disable(s);
        // Go back for another pass.
    }

    // Reset all the slip directions so that all slip->stick event witnesses 
    // will be positive when integration resumes.
    for (unsigned i=0; i < candidates.m_sliding.size(); ++i) {
        const int which = candidates.m_sliding[i];
        const MyFrictionElement& fric = getFrictionElement(which);
        if (!fric.isMasterActive(s)) continue;
        fric.updatePreviousSlipDirFromRecorded(s);
    }

    // Always leave at acceleration stage.
    m_mbs.realize(s, Stage::Acceleration);

    SimTK_DEBUG3("... Done; %d contacts, %d stictions broken; %d re-enabled.\n", 
        nContactsDisabled, nStictionDisabled, nContactsReenabled);

    showConstraintStatus(s, "EXIT findActiveCandidates()");
}

//---------------- CONVERGE NORMAL FORCES FOR SLIDING CONTACTS -----------------
namespace {
class NormalErrors : public Differentiator::JacobianFunction {
public:
    NormalErrors(const MyUnilateralConstraintSet& unis,
                 const Array_<int>&               sliders,
                 State&                           state) 
        : Differentiator::JacobianFunction(sliders.size(), sliders.size()), 
          m_unis(unis), m_state(state), m_sliders(sliders) { }

    void getEstimates(const State& s, Vector& estN) const {
        estN.resize(m_sliders.size());
        for (int i=0; i < (int)m_sliders.size(); ++i) {
            const MyFrictionElement& friction = 
                m_unis.getFrictionElement(m_sliders[i]);
            estN[i] = friction.getEstimatedNormalForce(s);
        }
    }

    // Alternate signature for convenience (extra copy).
    Vector getEstimates(const State& s) {
        Vector estN(m_sliders.size());
        getEstimates(s,estN);
        return estN;
    }

    void setEstimates(State& s, const Vector& estN) const {
        for (int i=0; i < (int)m_sliders.size(); ++i) {
            const MyFrictionElement& friction = 
                m_unis.getFrictionElement(m_sliders[i]);
            friction.setEstimatedNormalForce(m_state, estN[i]);
        }
    }

    void calcErrors(const State& s, Vector& errs) const {
        errs.resize(m_sliders.size());
        m_unis.getMultibodySystem().realize(s, Stage::Acceleration);
        for (int i=0; i < (int)m_sliders.size(); ++i) {
            const MyFrictionElement& friction = 
                m_unis.getFrictionElement(m_sliders[i]);
            const Real estN = friction.getEstimatedNormalForce(s);
            const Real N = friction.getMasterNormalForce(s);
            errs[i] = estN - N;
        }
    }

    // Alternate signature for convenience (extra copy).
    Vector calcErrors(const State& s) const {
        Vector errs(m_sliders.size()); calcErrors(s, errs);
        return errs;
    }

    // Must provide this pure virtual function.
    int f(const Vector& y, Vector& fy) const OVERRIDE_11 {
        assert(y.size() == m_sliders.size() && fy.size() == m_sliders.size());
        Vector prevEst(m_sliders.size());
        getEstimates(m_state, prevEst);
        setEstimates(m_state, y);
        calcErrors(m_state, fy);
        setEstimates(m_state, prevEst); // restore
        return 0;
    }
private:
    const MyUnilateralConstraintSet& m_unis;
    const Array_<int>&               m_sliders;
    State&                           m_state;
};
}


void MyUnilateralConstraintSet::
convergeNormalForcesForSlidingContacts(State& s, Real rtol) const {
    Array_<int> sliders;
    for (int i=0; i < getNumFrictionElements(); ++i) {
        const MyFrictionElement& friction = getFrictionElement(i);
        if (friction.isMasterActive(s) && !friction.isSticking(s))
            sliders.push_back(i);
    }
    SimTK_DEBUG1("convergeNormalForces: %d sliding contacts\n", sliders.size());
    if (sliders.empty()) return;

    NormalErrors normErrs(*this,sliders,s);
    Differentiator dnormErrs(normErrs);

    Vector estN(sliders.size()), err(sliders.size()), destN(sliders.size()),
           newEstN(sliders.size());
    normErrs.getEstimates(s, estN);
    normErrs.calcErrors(s, err);

    Real prevNorm = err.normInf();

    #ifndef NDEBUG
    cout << "====> Initial estN=" << estN << "-> err=" << err
         << " norm=" << prevNorm << endl;
    #endif

    if (prevNorm <= rtol)
        return;

    #ifndef NDEBUG
    cout << "\nNeed to converge sliding force and normal:\n";
    #endif

    Matrix J(sliders.size(), sliders.size());
    bool improved;
    const int MaxIters = 7;
    int iter = 0;
    do {
        ++iter;
        improved = false;
        dnormErrs.calcJacobian(estN, err, J, Differentiator::CentralDifference);
        FactorQTZ Jinv(J);
        #ifndef NDEBUG
        cout << "====> iter " << iter << " rank " << Jinv.getRank() 
             << " Jacobian\n" << J << endl;
        #endif  
        Jinv.solve(err, destN);
        Real fac = 2, norm;
        Real worseningAllowance = iter == 1 ? 1.2 : 1.01;
        do {
            fac /= 2;
            newEstN = estN - fac*destN;
            normErrs.setEstimates(s, newEstN);
            normErrs.calcErrors(s, err);
            norm = err.normInf();
            #ifndef NDEBUG
            printf("fac=%g new norm=%g", fac, norm);
            cout << " est=" << newEstN << endl;
            #endif  
        } while (norm > worseningAllowance*prevNorm);
        if (iter==1 || norm < prevNorm) {
            prevNorm = norm;
            estN = newEstN;
            improved = true;
        }
        #ifndef NDEBUG
        cout << "====> revised estN=" << estN << "-> err=" << err 
             << " norm=" << prevNorm << endl;
        #endif  
    } while (iter < MaxIters && improved && err.normInf() > rtol);
}


//-------------------------- SHOW CONSTRAINT STATUS ----------------------------
void MyUnilateralConstraintSet::
showConstraintStatus(const State& s, const String& place) const
{
#ifndef NDEBUG
    printf("\n%s: uni status @t=%.15g\n", place.c_str(), s.getTime());
    m_mbs.realize(s, Stage::Acceleration);
    for (int i=0; i < getNumContactElements(); ++i) {
        const MyContactElement& contact = getContactElement(i);
        const bool isActive = !contact.isDisabled(s);
        printf("  %6s %2d %s, h=%g dh=%g f=%g\n", 
                isActive?"ACTIVE":"off", i, contact.getContactType().c_str(), 
                contact.getPerr(s),contact.getVerr(s),
                isActive?contact.getForce(s):Zero);
    }
    for (int i=0; i < getNumFrictionElements(); ++i) {
        const MyFrictionElement& friction = getFrictionElement(i);
        if (!friction.isMasterActive(s))
            continue;
        const bool isSticking = friction.isSticking(s);
        printf("  %8s friction %2d, |v|=%g witness=%g\n", 
                isSticking?"STICKING":"sliding", i,
                friction.getActualSlipSpeed(s),
                isSticking?friction.calcStictionForceWitness(s)
                          :friction.calcSlipSpeedWitness(s));
        friction.writeFrictionInfo(s, "    ", std::cout);
    }
    printf("\n");
#endif
}

//==============================================================================
//               MY SLIDING FRICTION FORCE -- Implementation
//==============================================================================
// This is used for friction when slipping or when slip is impending, so it
// will generate a force even if there is no slip velocity. Whenever the slip
// speed |v|<=vtol, or if it has reversed direction, we consider it unreliable 
// and leave the applied force direction unchanged until the next transition 
// event. At that point activating the stiction constraint will be attempted. 
// If the stiction condition is violated, a new impending slip direction is 
// recorded opposing the direction of the resulting constraint force.
//
// The friction force is a 2-vector F calculated at Dynamics stage, applied 
// equal and opposite to the two contacting bodies at their mutual contact
// point:
//      F = -mu*N_est*d_eff
// d_eff is the effective slip direction that is to be opposed by the force.
//
// This is composed of several functions:
//      shouldUpdate(v) = ~v*d_prev > 0 && |v|>vtol
//      d_eff(v)   = shouldUpdate(v) ? v/|v| : d_prev
//      v_eff(v)   = ~v*d_eff
//      mu(v)      = mu_d + mu_v * max(v_eff,0)
//      Fmag(v; N) = mu(v)*N
//      F(v; N)    = -Fmag(v,N)*d_eff(v)
//      N_est      = N_prev
//
// mu_d  ... the dynamic coefficient of friction (a scalar constant >= 0)
// mu_v  ... viscous coefficient of friction (>= 0)
// N_prev... previous normal force magnitude (a discrete state variable)
// d_prev... previous or impending slip direction (a discrete state variable)
// d_eff ... the effective slip direction, a unit length 2-vector
// v_eff ... slip speed in d_eff direction, for viscous friction
//
// There is a sliding-to-stiction event witness function
//              e1(v)=dot(v, d_prev)    velocity reversal
//
// TODO: other possible witness functions:
//              e2(v)=|v|-vtol                slowdown (would often be missed)
//              e(v) = dot(v, d_prev) - vtol  [signed]
//
// N_prev is an auto-update variable whose update value is set at Acceleration
// stage from the actual normal force magnitude N of this friction element's 
// master contact element.  N_prev is also set manually whenever sliding is 
// enabled, to the current normal force. In general we have Ni=Ni(F) (where
// F={Fi} is a vector of all nf active sliding friction forces), and 
// Fi=Fi(Ni_est), so the error in the magnitude of the i'th applied friction 
// force is Fi_err=mu_i(v_i)*(Ni_est-Ni). If this is too large we have to 
// iterate until Ni_est is close enough to Ni for all i (this has to be done 
// simultaneously for the system as a whole).
//
// d_prev is an auto-update variable whose update value is set at Velocity
// stage, if shouldUpdate(v), otherwise it remains unchanged. It is also set 
// manually when transitioning from sticking to sliding, to -F/|F| where F was 
// the last stiction force vector.
// 
class MySlidingFrictionForceImpl : public Force::Custom::Implementation {
public:
    MySlidingFrictionForceImpl(const GeneralForceSubsystem&     forces,
                               const MyPointContactFriction&    ptFriction,
                               Real                             vtol)
    :   m_forces(forces), m_friction(ptFriction), 
        m_contact(ptFriction.getMyPointContact()), m_vtol(vtol)
    {}

    bool hasSlidingForce(const State& s) const 
    {   return m_friction.isMasterActive(s) && !m_friction.isSticking(s); }


    // Calculate d_eff, the direction to be opposed by the sliding force.
    Vec2 getEffectiveSlipDir(const State& s) const {
        const Vec2 Vslip = m_contact.getSlipVelocity(s);
        const Vec2 prevVslipDir = getPrevSlipDir(s);
        if (shouldUpdate(Vslip, prevVslipDir, m_vtol)) { // TODO: tol?
            const Real v = Vslip.norm();
            return Vslip/v;
        }
        return prevVslipDir;
    }

    // This is useful for reporting.
    Real calcSlidingForceMagnitude(const State& state) const {
        if (!hasSlidingForce(state)) return 0;
        const Real slipV = m_contact.getSlipVelocity(state).norm();
        return calcSlidingForceMagnitude(state, slipV);
    }

    // This is the force that will be applied next.
    Vec2 calcSlidingForce(const State& state) const {
        if (!hasSlidingForce(state))
            return Vec2(0);

        Vec2 dEff = getEffectiveSlipDir(state);
        if (dEff.isNaN()) {
            SimTK_DEBUG("NO SLIDING DIRECTION AVAILABLE\n");
            return Vec2(0);
        }

        const Vec2 Vslip = m_contact.getSlipVelocity(state);
        const Real vEff = ~Vslip*dEff;

        const Real FMag = calcSlidingForceMagnitude(state, std::max(vEff,Zero));
        assert(!isNaN(FMag));

        const Vec2 F = -FMag*dEff;
        return F;
    }

    // Return the related contact constraint's normal force value and slip
    // velocity as recorded at the end of the last time step. Will be zero if 
    // the contact was not active then.
    Real getPrevN(const State& state) const {
        const Real& prevN = Value<Real>::downcast
           (m_forces.getDiscreteVariable(state, m_prevNix));
        return prevN;
    }
    void setPrevN(State& state, Real N) const {
        Real& prevN = Value<Real>::updDowncast
           (m_forces.updDiscreteVariable(state, m_prevNix));
        if (isNaN(N))
            printf("*** setPrevN(): N is NaN\n");
        //SimTK_DEBUG3("STATE CHG %d: N from %g->%g\n",
        //    m_friction.getFrictionIndex(), prevN, N);
        prevN = N;

    }
    Vec2 getPrevSlipDir(const State& state) const {
        const Vec2& prevSlipDir = Value<Vec2>::downcast
           (m_forces.getDiscreteVariable(state, m_prevSlipDirIx));
        return prevSlipDir;
    }
    void setPrevSlipDir(State& state, const Vec2& slipDir) const {
        Vec2& prevSlipDir = Value<Vec2>::updDowncast
           (m_forces.updDiscreteVariable(state, m_prevSlipDirIx));
        prevSlipDir = slipDir;
        SimTK_DEBUG3("STATE CHG %d: prevDir to %g %g\n",
            m_friction.getFrictionIndex(), slipDir[0], slipDir[1]);
    }
    //--------------------------------------------------------------------------
    //                       Custom force virtuals
    // Apply the sliding friction force if this is enabled.
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   OVERRIDE_11
    {
        if (!hasSlidingForce(state))
            return; // nothing to do 

        const MobilizedBody& bodyB = m_contact.getBody();
        const MobilizedBody& bodyP = m_contact.getPlaneBody();
        const Vec3& stationB = m_contact.getBodyStation();
        const Vec3 stationP = bodyB.findStationLocationInAnotherBody
                                                    (state, stationB, bodyP);
        const Vec2 fSlip = calcSlidingForce(state);
        const Vec3 forceG(fSlip[0], 0, fSlip[1]); // only X,Z components
        bodyB.applyForceToBodyPoint(state, stationB,  forceG, bodyForces);
        bodyP.applyForceToBodyPoint(state, stationP, -forceG, bodyForces);
    }

    // Sliding friction does not store energy.
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {return 0;}

    // Allocate state variables for storing the previous normal force and
    // sliding direction.
    void realizeTopology(State& state) const OVERRIDE_11 {
        // The previous normal force N is used as a first estimate for the 
        // mu*N sliding friction force calculated at Dynamics stage. However,
        // the update value N cannot be determined until Acceleration stage.
        m_prevNix = m_forces.allocateAutoUpdateDiscreteVariable
           (state, Stage::Dynamics, new Value<Real>(0), Stage::Acceleration);
        // The previous sliding direction is used in an event witness that 
        // is evaluated at Velocity stage.
        m_prevSlipDirIx = m_forces.allocateAutoUpdateDiscreteVariable
           (state, Stage::Velocity, new Value<Vec2>(Vec2(NaN)), 
            Stage::Velocity);
    }

    // If we're sliding, set the update value for the previous slip direction
    // if the current slip velocity is usable.
    void realizeVelocity(const State& state) const OVERRIDE_11 {
        if (!hasSlidingForce(state))
            return; // nothing to do 
        const Vec2 Vslip = m_contact.getSlipVelocity(state);
        const Vec2 prevVslipDir = getPrevSlipDir(state);

        if (shouldUpdate(Vslip, prevVslipDir, m_vtol)) {
            Vec2& prevSlipUpdate = Value<Vec2>::updDowncast
               (m_forces.updDiscreteVarUpdateValue(state, m_prevSlipDirIx));
            const Real v = Vslip.norm();
            const Vec2 slipDir = Vslip / v;
            prevSlipUpdate = slipDir;
            m_forces.markDiscreteVarUpdateValueRealized(state, m_prevSlipDirIx);

            #ifndef NDEBUG
            //printf("UPDATE %d: prevSlipDir=%g %g; now=%g %g; |v|=%g dot=%g vdot=%g\n",
            //    m_friction.getFrictionIndex(),
            //    prevVslipDir[0],prevVslipDir[1],slipDir[0],slipDir[1],
            //    v, ~slipDir*prevVslipDir, ~Vslip*prevVslipDir);
            #endif
        } else {
            #ifndef NDEBUG
            printf("NO UPDATE %d: prevSlipDir=%g %g; Vnow=%g %g; |v|=%g vdot=%g\n",
                m_friction.getFrictionIndex(),
                prevVslipDir[0],prevVslipDir[1],Vslip[0],Vslip[1],
                Vslip.norm(), ~Vslip*prevVslipDir);
            #endif
        }
    }

    // Regardless of whether we're sticking or sliding, as long as the master
    // contact is active use its normal force scalar as the update for our
    // saved normal force.
    void realizeAcceleration(const State& state) const OVERRIDE_11 {
        if (!m_friction.isMasterActive(state))
            return; // nothing to save
        const Real N = m_contact.getForce(state); // normal force
        const Real prevN = getPrevN(state);
        if (N==prevN) return; // no need for an update

        Real& prevNupdate = Value<Real>::updDowncast
           (m_forces.updDiscreteVarUpdateValue(state, m_prevNix));

        #ifndef NDEBUG
        //printf("realize(Acc) CACHE UPDATE %d: N changing from %g -> %g (%.3g)\n",
        //    m_friction.getFrictionIndex(), 
        //    prevN, N, std::abs(N-prevN)/std::max(N,prevN));
        #endif
        prevNupdate = N;
        m_forces.markDiscreteVarUpdateValueRealized(state, m_prevNix);
    }

    //--------------------------------------------------------------------------

private:
    // Given the norm of the slip velocity already calculated, determine the
    // magnitude of the slip force. If there is no viscous friction you can
    // pass a zero vEff since it won't otherwise affect the force.
    // Don't call this unless you know there may be a sliding force.
    Real calcSlidingForceMagnitude(const State& state, Real vEff) const {
        assert(vEff >= 0);
        const Real prevN = getPrevN(state);
        if (prevN <= 0) return 0; // no normal force -> no friction force

        const Real mu_d = m_friction.getDynamicFrictionCoef();
        const Real mu_v = m_friction.getViscousFrictionCoef();
        const Real fMag = (mu_d + mu_v*vEff)*prevN;
        return fMag;
    }

    // Determine whether the current slip velocity is reliable enough that
    // we should use it to replace the previous slip velocity.
    static bool shouldUpdate(const Vec2& Vslip, const Vec2& prevVslipDir,
                             Real velTol) {
        const Real v2 = Vslip.normSqr();
        if (prevVslipDir.isNaN())
            return v2 > 0; // we'll take anything

        // Check for reversal.
        bool reversed = (~Vslip*prevVslipDir < 0);
        return !reversed && (v2 > square(velTol));
    }

    const GeneralForceSubsystem&    m_forces;
    const MyPointContactFriction&   m_friction;
    const MyPointContact&           m_contact;
    const Real                      m_vtol;

    mutable DiscreteVariableIndex   m_prevNix;       // previous normal force
    mutable DiscreteVariableIndex   m_prevSlipDirIx; // previous slip direction
};

// This is the force handle's constructor; it just creates the force
// implementation object.
MySlidingFrictionForce::MySlidingFrictionForce
   (GeneralForceSubsystem&          forces,
    const MyPointContactFriction&   friction,
    Real                            vtol) 
:   Force::Custom(forces, new MySlidingFrictionForceImpl(forces,friction,vtol)) 
{}

Real MySlidingFrictionForce::getPrevN(const State& state) const 
{   return getImpl().getPrevN(state); }
void MySlidingFrictionForce::setPrevN(State& state, Real N) const 
{   getImpl().setPrevN(state,N); }
Vec2 MySlidingFrictionForce::getPrevSlipDir(const State& state) const 
{   return getImpl().getPrevSlipDir(state); }
bool MySlidingFrictionForce::hasPrevSlipDir(const State& state) const 
{   return !getImpl().getPrevSlipDir(state).isNaN(); }
void MySlidingFrictionForce::
setPrevSlipDir(State& state, const Vec2& slipDir) const
{   getImpl().setPrevSlipDir(state, slipDir); }
Real MySlidingFrictionForce::calcSlidingForceMagnitude(const State& state) const 
{   return getImpl().calcSlidingForceMagnitude(state); }
Vec2 MySlidingFrictionForce::calcSlidingForce(const State& state) const 
{   return getImpl().calcSlidingForce(state); }


const MySlidingFrictionForceImpl& MySlidingFrictionForce::
getImpl() const {
    return dynamic_cast<const MySlidingFrictionForceImpl&>(getImplementation()); 
}


//    --------------------------------
//    perr = ~p_BS*n - h
//    --------------------------------
void UniPointInPlaneImpl::calcPositionErrors      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const 
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Vec3       p_AS = findStationLocation(allX_AB, followerBody, 
                                                defaultFollowerPoint);
    const Transform& X_AB = getBodyTransform(allX_AB, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                     // C is material pt of B coincident with S

    // We'll calculate this scalar using B-frame vectors, but any frame would 
    // have done.
    perr[0] = dot(p_BC, defaultPlaneNormal) - defaultPlaneHeight;
}

//    --------------------------------
//    verr = ~v_CS_A*n
//    --------------------------------
void UniPointInPlaneImpl::calcPositionDotErrors      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const
{
    assert(V_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);
    //TODO: should be able to get p info from State

    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Transform& X_AB = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                   // C is material point of B coincident with S
    const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

    const Vec3       v_AS = findStationVelocity(s, V_AB, followerBody, 
                                                defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocity(s, V_AB, planeBody, p_BC);

    // Calculate this scalar using A-frame vectors.
    pverr[0] = dot( v_AS-v_AC, n_A );
}

//    -------------------------------------
//    aerr = ~(a_CS_A - 2 w_AB X v_CS_A) * n
//    -------------------------------------
void UniPointInPlaneImpl::calcPositionDotDotErrors      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const
{
    assert(A_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size() == 1);
    //TODO: should be able to get p and v info from State
    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Transform& X_AB = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                   // C is material point of B coincident with S
    const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

    const Vec3&      w_AB = getBodyAngularVelocityFromState(s, planeBody);
    const Vec3       v_AS = findStationVelocityFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocityFromState(s, planeBody, p_BC);

    const Vec3       a_AS = findStationAcceleration(s, A_AB, followerBody, 
                                                    defaultFollowerPoint);;
    const Vec3       a_AC = findStationAcceleration(s, A_AB, planeBody, p_BC);

    paerr[0] = dot( (a_AS-a_AC) - 2.*w_AB % (v_AS-v_AC), n_A );
}


