/* -------------------------------------------------------------------------- *
 *         Simbody(tm) - UnilateralPointContactWithFriction Example           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

/*
This example shows a manual approach to dealing with unilateral constraints in
Simbody, which does not currently have built-in support but has sufficiently
general facilities. In this example we'll implement non-slipping point
contact, joint limit constraints, and a rope-like one-sided distance 
constraint. We'll use Simbody bilateral constraints turned on and off with 
manual switching conditions that are set by discrete event handlers.

For each designated contact point that is not in contact, we'll track the 
vertical height over the ground plane and its first and second time derivatives
and use those to construct switching ("witness") functions to trigger an event 
that may enable the constraint. For each enabled contact constraint, we'll 
track the sign of the normal reaction force and use it as a witness to disable 
the constraint.

Note that there are two separate conditions involving these constraints:
impact (collision) and contact. Impact occurs during an infinitesimal 
interval and involves impulses and velocities, while contact persists over time
and involves forces and accelerations. Contact between rigid objects is a 
simple, physically justifiable process in which contact constraints generate 
forces if necessary to prevent interpenetration. Impact of rigid objects, on 
the other hand, requires assumptions to be made about the non-modeled details
of collision behavior that is assumed to occur in an infinitesimal interval.
Just producing logically-consistent behavior during impact is very difficult;
justifying it physically even more so.

How we handle contact
---------------------
In this example each contact consists of a constraint that prevents penetration
of a point on a moving body normal to the ground plane, and constraints
that prevent slipping tangent to the plane. We implement non-penetration with 
Simbody's "PointInPlane" constraint. We enable this 
constraint when a contact begins, defined so 
that its multiplier is the y component of the reaction force, with +y
being the ground plane normal. We monitor the reaction force y component, and
declare the contact broken if that component is negative. The no-slip condition
is enforced with two of Simbody's "NoSlip1D" constraints, one in the x 
direction and one in the z direction.

A rope is implemented similarly using Simbody's "Rod" (distance) constraint,
and joint stops are implemented (somewhat inadequately) using the existing
ConstantSpeed constraint, with the speed set to zero.

How we handle impacts
---------------------
In this example, an impact is signaled by a contact point that reaches the 
ground plane with a negative vertical speed vy, with similar conditions for
the other constraints. This requires a step change to
the system velocities to avoid penetration or constraint violation. We achieve
this step change by applying a constraint-space
impulse to the system, representing constraint-space contact 
forces integrated over the assumed-infinitesimal impact interval. The system
equations of motion are used to ensure that the velocity changes produced by
the impulses satisfy Newton's laws. This can produce velocity changes anywhere 
in the system and may result in other impacts or breaking of 
existing contacts.

When an impact is signaled, we determine the subset of potential contacts that
may be involved in this event; those are called "proximal" contacts and are 
just those whose contact points are at zero height, within a small tolerance.
The rest are ignored during handling of the impact.

We use Poisson's interpretation of coefficient of restitution as a ratio of
impulses, rather than Newton's more commonly known but inconsistent 
interpretation as a ratio of velocities. To apply Poisson's interpretation, we 
divide the impact into two distinct phases: compression and expansion. 
During compression we determine what impulse is required to prevent any 
penetration at the proximal contacts, by eliminating any negative speeds vy. 
The task for expansion is to determine an expansion impulse, based
on the compression impulse, the coefficient of restitution e at each contact, 
and a "capture velocity" vc that says when a rebound velocity is so small we 
should consider a new persistent contact to be initiated. 

0) Initialize: initialize the effective coefficients of restitution e(i) for
each of the proximal constraints. These can be constants associated with the
contact parameters, or can be calculated from the initial velocities. Set
the total applied impulse I=0, determine current velocity V. Activate all
proximal constraints.

1) Compression phase: Determine the nonnegative least squares constraint-space
impulse Ic that brings any impacting proximal contacts ("impacters") to a stop, 
and leaves non-impacting ones ("rebounders") with a positive vertical speed 
(however small). Note that the set of impacters might not end up being the same 
ones as came in with negative vy; some new ones might be added and some of the
originals might turn out to be rebounders due to the effects of other impacts.
At the end we have for the i'th proximal constraint a compression impulse 
Ic(i)>=0 and a post-compression velocity Vc(i)>=0, with contact constraint i 
active (impacter) if Ic(i)>0 or Ic(i)==0 && Vc(i)==0 and inactive 
(rebounder) otherwise, with Ic(i)==0 and Vc(i)>0. Although it may take a few
iterations to figure out what's going on, we consider everything to be
simultaneous during a compression phase -- there is a single impulse Ic
generated that modifies the original velocities to produce Vc. Increment the
total impulse I+=Ic (>=0), set V=Vc (>=0).

2) Expansion phase: Generate an expansion impulse Ie such that 
Ie(i)=e(i)*Ic(i) for each of the impacters i from the compression phase. 
If Ie==0 there is no expansion to do; go to step 4. Otherwise,  
set e(i)=0 for each of the impacters; the material restitution has now been 
consumed. Increment the total impulse I+=Ie. Apply the impulse Ie to produce a 
velocity change dVe and a new velocity Ve=V+dVe, and update V=Ve. If Ve>=0 for
all proximal contacts, we are successful. In that case go to step 4 with the 
total impulse I>=0, and velocity V=Ve>=0; active constraints are those where 
V(k)=0.

3) Some contacts now have negative
vertical speeds vy (these may include both impacters and rebounders from the
compression phase). This requires a new compression phase, beginning with
these velocities and with the original impacters now having zero coefficients
of restitution. So return to step 1.

4) We have determined and applied the compression+expansion impulse I>=0 and 
have the resulting velocities V>=0. Check all contacts for which V>0 (the 
rebounders) to see if any is rebounding very slowly (<= vc). Enable those,
and calculate the impulse dI that just brings those to zero while maintaining 
other contacts. Apply that impulse to get new velocities V. If that causes 
any V(k)<0 or new V(k)<=vc, declare that a contact too and recalculate dI; 
repeat until all inactive (rebounding) V(k)>vc. Then set the final I+=dI.

5) Now calculate accelerations. If any of the active proximal contacts 
generate a zero or negative vertical reaction force they should be disabled;
otherwise we would miss the next break-free event. 
*/

#include "Simbody.h"

#include <string>
#include <iostream>
#include <exception>

using std::cout;
using std::endl;

using namespace SimTK;

const Real ReportInterval=1./30;
const Real RunTime=10;
//const Real EventIsolationTol = 1e-6;

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
        m_index(-1), m_friction(0), m_restitutionDone(false) 
    {   m_uni.setDisabledByDefault(true); }

    virtual ~MyContactElement() {}

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
    virtual void initializeForImpact(const State& state) 
    {   setRestitutionDone(false); m_Ic = m_Ie = m_I = 0; }

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


    Real getCoefRest() const {return m_coefRest;}
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
    MyFrictionElement*  m_friction; // if any

    // Runtime
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
// separate "contact master", which is one of the following:
//  - a unilateral constraint
//  - a bilateral constraint 
//  - a mobilizer
//  - a compliant force element 
// If the master is an inactive unilateral constraint, or if N=0, then no 
// friction forces are generated.
//
// For all but the compliant force element master, the normal force N is 
// acceleration-dependent and thus may be coupled to the force produced by a
// sliding friction element. This may require iteration to ensure consistency
// between the sliding friction force and its master contact's normal force.
//
// A Coulomb friction element depends on a scalar slip speed defined by the
// contact master (this might be the magnitude of a generalized speed or
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

    // TODO: is this necessary?
    // For use when the current slip speed is too small, this is the dot 
    // produce of the slip acceleration with the previous slip direction.
    // State must be realized to Acceleration stage.
    virtual Real calcSlipAccelWitness(const State&) const = 0;

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

    // Return true if the master contact element *could* be involved in an 
    // impact event (because it is touching).
    virtual bool isMasterProximal(const State&, Real posTol) const = 0;
    // Return true if the master contact element *could* be involved in contact
    // force generation (because it is touching and not separating).
    virtual bool isMasterCandidate(const State&, Real posTol, Real velTol)
        const = 0;
    // Return true if the master contact element is currently generating a
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
    explicit MyUnilateralConstraintSet(const MultibodySystem& mbs)
    :   m_mbs(mbs) {}
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

    // Return the contact and friction elements that might be involved in an
    // impact occurring in this configuration. They are the contact elements 
    // for which perr <= posTol, and friction elements whose masters can be 
    // involved in the impact. State must be realized through Position stage.
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
    // whose master is a candidate and (a) which are already sticking, or (b)
    // for which vslip <= velTol, or (c) for which vslip opposes the previous 
    // slip direction, meaning it has reversed and must have passed through 
    // zero during the last step. These are the elements that can be activated 
    // without making any changes to the configuration or velocity state 
    // variables, except for constraint projection. 
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
                candidates.m_friction.push_back(i);
            } else {
                fric.recordSlipDir(s);
                candidates.m_sliding.push_back(i);
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
        //getTriggerInfo().setRequiredLocalizationTimeWindow(EventIsolationTol);
    }

    // This is the witness function.
    Real getValue(const State& state) const {
        const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
        const MyContactElement& uni = m_unis.getContactElement(m_which);
        if (!uni.isDisabled(state)) 
            return 0; // already locked

        const Real height = uni.getPerr(state);

        if (m_stage == Stage::Position)
            return height;

        // Velocity and acceleration triggers are not needed if we're
        // above ground.
        if (height > 0) return 0;

        const Real dheight = uni.getVerr(state);

        if (m_stage == Stage::Velocity)
            return dheight;

        // Acceleration trigger is not needed if velocity is positive.
        if (dheight > 0) return 0;

        const Real ddheight = uni.getAerr(state);

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
    // in proximal, and state is updated to the new velocities and realized
    // through Velocity stage. Constraints that were stopped are enabled, those
    // that rebounded are disabled.
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

    // Examine the rebounders to see if any are rebounding with a speed at or
    // below the capture velocity. If so, enable those constraints and apply a
    // (hopefully small) negative impulse to eliminate that rebound velocity.
    // Repeat if that induces any negative velocities or any further slow
    // rebounders. This terminates will all rebounders leaving with velocities
    // greater than vCapture, or else all constraints are enabled.
    void captureSlowRebounders(Real             vCapture,
                               MyElementSubset& proximal,
                               State&           state) const;

    // This method is used at the start of compression phase to modify any
    // constraint parameters as necessary, and then enable all the proximal
    // constraints. Some or all of these will be disabled during the impact
    // analysis in compression or expansion phases. On return the state has
    // been updated and realized through Instance stage.
    void enableAllProximalConstraints(MyElementSubset&  proximal,
                                      State&            state) const;

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
        //getTriggerInfo().setRequiredLocalizationTimeWindow(EventIsolationTol);
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
            Constraint::PointInPlane(updGround(body), UnitVec3(YAxis), Zero,
                                     body, point),
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
    MySlidingFrictionForce(GeneralForceSubsystem& forces,
                           const class MyPointContactFriction& ptFriction);

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
public:
    MyPointContactFriction(MyPointContact& contact,
        Real mu_d, Real mu_s, Real mu_v, GeneralForceSubsystem& forces)
    :   MyFrictionElement(mu_d,mu_s,mu_v), m_contact(contact),
        m_noslipX(contact.updPlaneBody(), Vec3(0), UnitVec3(XAxis), 
                  contact.updPlaneBody(), contact.updBody()),
        m_noslipZ(contact.updPlaneBody(), Vec3(0), UnitVec3(ZAxis), 
                  contact.updPlaneBody(), contact.updBody()),
        m_contactPointInPlane(0), m_tIc(NaN), m_tIe(NaN), m_tI(NaN),
        m_prevN(0), m_prevSlipDir(Vec2(NaN))
    {
        assert((0 <= mu_d && mu_d <= mu_s) && (0 <= mu_v));
        contact.setFrictionElement(*this);
        m_noslipX.setDisabledByDefault(true);
        m_noslipZ.setDisabledByDefault(true);
        m_sliding = new MySlidingFrictionForce(forces, *this);
    }

    ~MyPointContactFriction() {delete m_sliding;}

    Vec2 getStictionForce(const State& s) const {
        assert(isSticking(s));
        return Vec2(-m_noslipX.getMultiplier(s), -m_noslipZ.getMultiplier(s));
    }

    // Implement pure virtuals from MyFrictionElement base class.

    bool isSticking(const State& s) const OVERRIDE_11
    {   return !m_noslipX.isDisabled(s); } // X,Z always the same

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

    Real calcSlipAccelWitness(const State& s) const OVERRIDE_11 {
        if (isSticking(s)) return 0;
        const Vec2 aNow = m_contact.getSlipAcceleration(s);
        const Vec2 vPrev = m_sliding->getPrevSlipDir(s);
        return dot(aNow, vPrev);
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
        const Real Scale = 10;
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
                geometry.push_back(DecorativeSphere(.5)
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
                    .setColor(White).setScale(.5)
                    .setTransform(loc+Vec3(0,.2,0)));
            } else {
                geometry.push_back(DecorativeText(String(i))
                    .setColor(White).setScale(.5)
                    .setTransform(loc+Vec3(0,.1,0)));
            }
        }
    }
private:
    const MyUnilateralConstraintSet& m_unis;
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
        //getTriggerInfo().setRequiredLocalizationTimeWindow(EventIsolationTol);
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
        //getTriggerInfo().setRequiredLocalizationTimeWindow(EventIsolationTol);
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
//                                  MY STOP
//==============================================================================
// Define a unilateral constraint to represent a joint stop that limits
// the allowable motion of a single generalized coordinate. You can specify
// a coefficient of restitution and whether the given limit is the upper or
// lower limit.
class MyStop : public MyContactElement {
public:
    enum Side {Lower,Upper};
    MyStop(Side side, MobilizedBody& body, int whichQ,
         Real limit, Real coefRest)
    :   MyContactElement
           (Constraint::ConstantSpeed(body, MobilizerUIndex(whichQ), Real(0)), 
            Real(side==Lower?-1:1), coefRest),
        m_body(body), m_whichq(whichQ), m_whichu(whichQ),
        m_sign(side==Lower?1.:-1.), m_limit(limit)
    {}

    String getContactType() const OVERRIDE_11 {return "Stop";}

    Real getPerr(const State& state) const OVERRIDE_11 {
        const Real q = m_body.getOneQ(state, m_whichq);
        return m_sign*(q-m_limit);
    }
    Real getVerr(const State& state) const OVERRIDE_11 {
        const Real u = m_body.getOneU(state, m_whichu);
        return m_sign*u;
    }
    Real getAerr(const State& state) const OVERRIDE_11 {
        const Real udot = m_body.getOneUDot(state, m_whichu);
        return m_sign*udot;
    }

    Vec3 whereToDisplay(const State& state) const OVERRIDE_11 {
        const Vec3& p_B = m_body.getOutboardFrame(state).p();
        return m_body.findStationLocationInGround(state,p_B);
    }

private:
    const MobilizedBody&        m_body;
    const MobilizerQIndex       m_whichq;
    const MobilizerUIndex       m_whichu;
    Real                        m_sign; // +1: lower, -1: upper
    Real                        m_limit;
};

//==============================================================================
//                                  MY ROPE
//==============================================================================
// Define a unilateral constraint to represent a "rope" that keeps the
// distance between two points at or smaller than some limit.
class MyRope : public MyContactElement {
public:
    MyRope(MobilizedBody& body1, const Vec3& pt1,
           MobilizedBody& body2, const Vec3& pt2, Real d,
           Real coefRest)
    :   MyContactElement
           (Constraint::Rod(body1, pt1, body2, pt2, d), Real(1), coefRest),
        m_body1(body1), m_point1(pt1), m_body2(body2), m_point2(pt2), m_dist(d)
    {}

    String getContactType() const OVERRIDE_11 {return "Rope";}

    Real getPerr(const State& s) const OVERRIDE_11 {
        const Vec3 p1 = m_body1.findStationLocationInGround(s,m_point1);
        const Vec3 p2 = m_body2.findStationLocationInGround(s,m_point2);
        const Vec3 p = p2-p1;
        return (square(m_dist) - dot(p,p))/2;
    }
    Real getVerr(const State& s) const OVERRIDE_11 {
        Vec3 p1, v1, p2, v2;
        m_body1.findStationLocationAndVelocityInGround(s,m_point1,p1,v1);
        m_body2.findStationLocationAndVelocityInGround(s,m_point2,p2,v2);
        const Vec3 p = p2 - p1, v = v2 - v1;
        return -dot(v, p);
    }
    Real getAerr(const State& s) const OVERRIDE_11 {
        Vec3 p1, v1, a1, p2, v2, a2;
        m_body1.findStationLocationVelocityAndAccelerationInGround
           (s,m_point1,p1,v1,a1);
        m_body2.findStationLocationVelocityAndAccelerationInGround
           (s,m_point2,p2,v2,a2);
        const Vec3 p = p2 - p1, v = v2 - v1, a = a2 - a1;
        return -(dot(a, p) + dot(v, v));
    }

    Vec3 whereToDisplay(const State& state) const OVERRIDE_11 {
        return m_body2.findStationLocationInGround(state,m_point2);
    }

private:
    const MobilizedBody&    m_body1;
    const Vec3              m_point1;
    const MobilizedBody&    m_body2;
    const Vec3              m_point2;
    const Real              m_dist;
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
                    const Vec3&          forceB,
                    Real                 onTime,
                    Real                 offTime)
    :   m_bodyB(bodyB), m_stationB(stationB), m_forceB(forceB),
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

        const Rotation& R_GB = m_bodyB.getBodyRotation(state);
        const Vec3      forceG = R_GB*m_forceB;
        m_bodyB.applyForceToBodyPoint(state, m_stationB, forceG, bodyForces);

        SimTK_DEBUG4("PUSHING @t=%g (%g,%g,%g)\n", state.getTime(),
            forceG[0], forceG[1], forceG[2]);
    }

    // No potential energy.
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {return 0;}

    void calcDecorativeGeometryAndAppend
       (const State& state, Stage stage, 
        Array_<DecorativeGeometry>& geometry) const OVERRIDE_11
    {
        if (stage != Stage::Time) return;
        if (!(m_on <= state.getTime() && state.getTime() <= m_off))
            return;
        geometry.push_back(DecorativeLine(m_stationB-m_forceB, m_stationB)
                            .setBodyId(m_bodyB.getMobilizedBodyIndex())
                            .setColor(Yellow)
                            .setLineThickness(3));
    }
private:
    const MobilizedBody& m_bodyB; 
    const Vec3           m_stationB;
    const Vec3           m_forceB;
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


        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS
    const Real CoefRest = 0.2;      // TODO: per-contact
    const Real mu_d = .4;
    const Real mu_s = .8;
    const Real mu_v = 0.05;

    MyUnilateralConstraintSet unis(mbs);

    const Vec3 CubeHalfDims(3,2,1);
    const Real CubeMass = 1;
    Body::Rigid cubeBody = 
        Body::Rigid(MassProperties(CubeMass, Vec3(0), 
                    UnitInertia::brick(CubeHalfDims)));

    // First body: cube
    MobilizedBody::Cartesian loc(Ground, MassProperties(0,Vec3(0),Inertia(0)));
    MobilizedBody::Ball cube(loc, Vec3(0),
                             cubeBody, Vec3(0));
    cube.addBodyDecoration(Transform(), DecorativeBrick(CubeHalfDims)
                                        .setColor(Red).setOpacity(.3));
    Force::Custom(forces, new MyPushForceImpl(cube,Vec3(1,1,0),2*Vec3(0,0,10),
                                               4., 6.));
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(CubeHalfDims);
        MyPointContact* contact = new MyPointContact(cube, pt, CoefRest);
        unis.addContactElement(contact);
        unis.addFrictionElement(
            new MyPointContactFriction(*contact, mu_d, mu_s, mu_v, forces));
    }

    //unis.push_back(new MyRope(Ground, Vec3(-5,10,0),
    //                          cube, Vec3(-CubeHalfDims[0],-CubeHalfDims[1],0), 
    //                          5., .5*CoefRest));
    //unis.push_back(new MyStop(MyStop::Upper,loc,1, 2.5,CoefRest));

    const Vec3 WeightEdge(-CubeHalfDims[0],-CubeHalfDims[1],0);
//#ifdef NOTDEF
    // Second body: weight
    const Vec3 ConnectEdge1(CubeHalfDims[0],0,CubeHalfDims[2]);
    MobilizedBody::Pin weight(cube, 
        Transform(Rotation(Pi/2,XAxis), ConnectEdge1),
        cubeBody, Vec3(WeightEdge));
    weight.addBodyDecoration(Transform(), DecorativeBrick(CubeHalfDims)
                                        .setColor(Gray).setOpacity(.6));
    //Force::MobilityLinearSpring(forces, weight, 0, 100, -Pi/4);
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        if (i==-1 && j==-1) continue;
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(CubeHalfDims);
        MyPointContact* contact = new MyPointContact(weight, pt, CoefRest);
        unis.addContactElement(contact);
        unis.addFrictionElement(
            new MyPointContactFriction(*contact, mu_d, mu_s, mu_v, forces));
    }
    unis.addContactElement(new MyStop(MyStop::Upper,weight,0, Pi/6,CoefRest/2));
    unis.addContactElement(new MyStop(MyStop::Lower,weight,0, -Pi/6,CoefRest/2));

//#endif
#ifdef NOTDEF
   // Third body: weight2
    const Vec3 ConnectEdge2(CubeHalfDims[0],CubeHalfDims[1],0);
    MobilizedBody::Pin weight2(cube, 
        Transform(Rotation(), ConnectEdge2),
        cubeBody, Vec3(WeightEdge));
    weight2.addBodyDecoration(Transform(), DecorativeBrick(CubeHalfDims)
                                        .setColor(Cyan).setOpacity(.6));
    Force::MobilityLinearSpring(forces, weight2, 0, 1000, Pi/4);
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        if (i==-1 && j==-1) continue;
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(CubeHalfDims);
        unis.push_back(new MyPointContact(weight2, pt, CoefRest, forces));
    }
#endif

    Visualizer viz(mbs);
    viz.setShowSimTime(true);
    viz.addDecorationGenerator(new ShowContact(unis));
    mbs.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));

    //ExplicitEulerIntegrator integ(mbs);
    //CPodesIntegrator integ(mbs,CPodes::BDF,CPodes::Newton);
    //RungeKuttaFeldbergIntegrator integ(mbs);
    Real accuracy = 1e-3;
    RungeKuttaMersonIntegrator integ(mbs);
    //RungeKutta3Integrator integ(mbs);
    //VerletIntegrator integ(mbs);
    integ.setAccuracy(accuracy);
    //integ.setAllowInterpolation(false);
    integ.setMaximumStepSize(0.5);

    StateSaver* stateSaver = new StateSaver(mbs,unis,integ,ReportInterval);
    mbs.addEventReporter(stateSaver);
    
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
  
    State s = mbs.realizeTopology(); // returns a reference to the the default state
    mbs.realizeModel(s); // define appropriate states for this System
    mbs.realize(s, Stage::Instance); // instantiate constraints if any


    // Set initial conditions so the -,-,- vertex is in the -y direction.
    const Rotation R_BC(UnitVec3(CubeHalfDims+1e-7*Vec3(1,0,0)), YAxis, Vec3(1,0,0),XAxis);
    loc.setQToFitTranslation(s, Vec3(0,10,0));
    cube.setQToFitTransform(s, Transform(~R_BC, Vec3(0)));
    cube.setUToFitAngularVelocity(s, Vec3(0,1,0));
    cube.setUToFitLinearVelocity(s, Vec3(1,0,0));

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

    Assembler(mbs).assemble(s);
    viz.report(s);
    cout << "Assembled configuration shown. Ready? ";
    getchar();

    // Now look for enabled contacts that aren't sliding; turn on stiction
    // for those.
    mbs.realize(s, Stage::Velocity);
    Array_<int> enableTheseStictions;
    for (int i=0; i < unis.getNumFrictionElements(); ++i) {
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
    
    // Simulate it.

    integ.setReturnEveryInternalStep(true);
   //integ.setAllowInterpolation(false);
    TimeStepper ts(mbs, integ);
    ts.setReportAllSignificantStates(true);
    ts.initialize(s);

    const double startReal = realTime();
    const double startCPU = cpuTime();

    Integrator::SuccessfulStepStatus status;
    do {
        status=ts.stepTo(RunTime);
        //printf("Returned at t=%g with status=%s\n", ts.getTime(),
        //    Integrator::getSuccessfulStepStatusString(status).c_str());
        if (status != Integrator::ReachedReportTime) {
            //viz.report(ts.getState());
            //stateSaver->handleEvent(ts.getState());
        }
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

    // Instant replay.
    while(true) {
        for (int i=0; i < stateSaver->getNumSavedStates(); ++i) {
            mbs.realize(stateSaver->getState(i), Stage::Velocity);
            viz.report(stateSaver->getState(i));
        }
        getchar();
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
    const Real VCapture=1e-2;

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
        m_unis.updContactElement(which).initializeForImpact(s);
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

    // Some of the rebounders may be moving so slowly that we would like
    // to be able to say they have stopped. If so, apply additional 
    // (negative) impulses necessary to stop them; enable their contact
    // constraints.
    captureSlowRebounders(VCapture, proximal, s);

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
//
// Strategy:
// (1) assume all normal & tangential constraints are on; calculate stopping
//     impulse
// (2) look for negative normal impulses; try disabling those
// (3) recapture any constraints that would be violated after disabling
// 
void ContactOn::
processCompressionPhase(MyElementSubset&    proximal,
                        State&              s) const
{
    SimTK_DEBUG("Entering processCompressionPhase() ...\n");
    Vector lambda0, lambdaTry;
    // Which constraint, whether stiction was on too.
    Array_<int> maybeDisabled, recapturing;

    const Vector u0 = s.getU(); // save presenting velocity

    // Assume at first that everyone will participate. This is necessary to
    // ensure that we get a least squares solution for the impulse involving
    // as many constraints as possible sharing the impulse.
    enableAllProximalConstraints(proximal, s);

    // First try drives all constraints to zero velocity; if that took some
    // negative impulses we'll have to remove some participants and try again.
    calcStoppingImpulse(proximal, s, lambda0);

    while (true) {
        // See if negative impulse was required to satisfy an active constraint.
        // If so, that is a candidate to be inactivated.
        maybeDisabled.clear();
        for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
            const int which = proximal.m_contact[i];
            const MyContactElement& cont = m_unis.getContactElement(which);
            if (cont.isDisabled(s))
                continue;
            Real imp = cont.getMyValueFromConstraintSpaceVector(s, lambda0);
            if (imp > 0)
                continue;
            maybeDisabled.push_back(which);

            SimTK_DEBUG2("  cont constraint %d has negative compression impulse=%g\n",
                which, imp);
        }
        if (maybeDisabled.empty())
            break;

        // Disable the candidates, then see if they rebound.
        for (unsigned d=0; d < maybeDisabled.size(); ++d) {
            const int which = maybeDisabled[d];
            const MyContactElement& cont = m_unis.getContactElement(which);
            cont.disable(s);
            if (cont.hasFrictionElement()) {
                const MyFrictionElement& fric = cont.getFrictionElement();
                fric.disableStiction(s);
            }
        }
        m_mbs.realize(s, Stage::Instance);
        calcStoppingImpulse(proximal, s, lambdaTry);
        updateVelocities(u0, lambdaTry, s);

        recapturing.clear();
        for (unsigned i=0; i<maybeDisabled.size(); ++i) {
            const int which = maybeDisabled[i];
            const MyContactElement& uni = m_unis.getContactElement(which);
            const Real newV = uni.getVerr(s);           
            SimTK_DEBUG2("  candidate uni constraint %d would have v=%g\n",
                   which, newV);
            if (newV <= 0) {
                recapturing.push_back(maybeDisabled[i]);
                SimTK_DEBUG2("  RECAPTURING uni constraint %d with v=%g\n", 
                    which, newV);
            }
        }

        // Re-enable the recaptured candidates.
        if (!recapturing.empty()) {
            for (unsigned c=0; c < recapturing.size(); ++c) {
                const int which = recapturing[c];
                m_unis.getContactElement(which).enable(s);
            }
            m_mbs.realize(s, Stage::Instance);
        }

        if (recapturing.empty()) lambda0 = lambdaTry;
        else calcStoppingImpulse(proximal, s, lambda0);

        const int numDisabled = maybeDisabled.size()-recapturing.size();
        if (numDisabled == 0) {
            SimTK_DEBUG("  None of the candidates was actually disabled.\n");
            break;
        }

        SimTK_DEBUG1("  RETRY with %d constraints disabled\n", numDisabled);
    }
    updateVelocities(u0, lambda0, s);

    // Now update the entries for each proximal constraint to reflect the
    // compression impulse and post-compression velocity.
    SimTK_DEBUG("  Compression results:\n");
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        MyContactElement& uni = m_unis.updContactElement(which);
        if (!uni.isDisabled(s))
            uni.recordImpulse(MyContactElement::Compression, s, lambda0);
        SimTK_DEBUG4("  %d %3s: Ic=%g, V=%g\n",
            which, uni.isDisabled(s) ? "off" : "ON", 
            uni.getCompressionImpulse(), uni.getVerr(s));
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
        if (uni.isDisabled(s)||uni.isRestitutionDone()||uni.getCoefRest()==0
            ||uni.getCompressionImpulse()<=0)
            continue;
        uni.setMyExpansionImpulse(s, uni.getCoefRest(), expansionImpulse);
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


//------------------------- CAPTURE SLOW REBOUNDERS ----------------------------
void ContactOn::
captureSlowRebounders(Real              vCapture,
                      MyElementSubset&  proximal,
                      State&            s) const
{
    // Capture any rebounder whose velocity is too slow.
    SimTK_DEBUG("Entering captureSlowRebounders() ...\n");

    int nCaptured=0, nPasses=0;
    while (true) {
        ++nPasses;
        SimTK_DEBUG1("  start capture pass %d:\n", nPasses);
        Array_<int> toCapture;
        for (unsigned i=0; i<proximal.m_contact.size(); ++i) {
            const int which = proximal.m_contact[i];
            MyContactElement& uni = m_unis.updContactElement(which);
            if (uni.isDisabled(s) && uni.getVerr(s) <= vCapture) {
                toCapture.push_back(which);
                ++nCaptured;
                SimTK_DEBUG2("  capturing constraint %d with v=%g\n",
                    which, uni.getVerr(s));
            }
        }

        if (toCapture.empty()) {
            if (nCaptured==0) SimTK_DEBUG("... done -- nothing captured.\n");
            else SimTK_DEBUG2("... done -- captured %d rebounders in %d passes.\n",
                nCaptured, nPasses);
            return;
        }

        // Now do the actual capturing by enabling constraints.
        for (unsigned i=0; i < toCapture.size(); ++i) {
            const int which = toCapture[i];
            MyContactElement& uni = m_unis.updContactElement(which);
            uni.enable(s);
        }

        m_mbs.realize(s, Stage::Velocity);
        Vector captureImpulse;
        calcStoppingImpulse(proximal, s, captureImpulse);
        updateVelocities(Vector(), captureImpulse, s);

        // Now update the entries for each proximal constraint to reflect the
        // capture impulse and post-capture velocity.
        for (unsigned i=0; i<proximal.m_contact.size(); ++i) {
            const int which = proximal.m_contact[i];
            MyContactElement& uni = m_unis.updContactElement(which);
            if (!uni.isDisabled(s))
                uni.recordImpulse(MyContactElement::Capture, 
                                  s, captureImpulse);
        }
    }
}



//---------------------- ENABLE ALL PROXIMAL CONSTRAINTS -----------------------
void ContactOn::
enableAllProximalConstraints(MyElementSubset&   proximal,
                             State&             state) const
{
    // TODO: don't enable stiction for now
    m_unis.enableConstraintSubset(proximal, false, state); 
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
    findCandidateElements(state, accuracy, accuracy, candidates);

    // Evaluate accelerations and reaction forces and check if 
    // any of the active contacts are generating negative ("pulling") 
    // forces; if so, inactivate them.
    findActiveCandidates(state, candidates);

    // Finally, project active constraints to the constraint manifolds.
    m_mbs.project(state, accuracy);

    // It is possible that the projection above took some marginally-sliding
    // friction and slowed it down enough to make it a stiction candidate.
    // TODO: check and restart if so.
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
// Heuristic used here:
// (1) Check stiction force magnitudes. If any exceeds mu_s*N, record the
//     impending slip direction d=-f/|f|, and estimated normal force
//     N_est = max(N,0). Switch to sliding with f=-mu_d*N_est*d.
// (2) Recalculate normal forces with sliding. Replace all estimated normal
//     forces with new values, calculating the size of the normal force change.
//     If change is large, repeat (2) until the change is within tolerance.
// (3) If there are new stiction forces exceeding their limits, return to (1).
// (4) Check for negative unilateral contact forces. Disable those contacts
//     tentatively and recalculate contact forces.
// (5) Check for negative unilateral accelerations. Enable those contacts
//     in sliding or stiction depending on velocity.
//
// ALTERNATIVE (NOT IMPLEMENTED):
// Activate all candidate constraints (no impulse needed to enable):
// - for uni constraints, perr <= ptol, verr <= vtol
// - for stiction with candidate master, |v| <= vtol.
//
// loop
// - calculate ferr for active constraints, aerr for inactive
//     (iterate sliding forces until f=mu_d*N to tol)
// - if all ferr>=0, aerr>=0, break 
// - find worst offender (TODO)
// - inactivate/activate worst offender
//     (record impending slip direction & N for stick->slide)
// end loop 
//
// Notes:
// - for stiction, ferr=mu_s*N - |f|
// - for sliding, aerr=dot(slipAccel, impendingSlipDir)
// - for active contact, ferr=lambda
// - for inactive contact, aerr=separation accel

void MyUnilateralConstraintSet::
findActiveCandidates(State& s, const MyElementSubset& candidates) const
{
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

    // Enable all candidate contact and stiction constraints.
    enableConstraintSubset(candidates, true, s);

    // Check first, disable later because we don't want to invalidate
    // the reaction forces in the state yet.
    Array_<int> stictionDisabled, contactDisabled;

    // Disable any stiction constraints that are generating too much force.
    // Note that this will be true for any stiction constraint whose master
    // unilateral contact is generating a negative normal force.
    int pass = 0; Array_<int> disabled;
    while (true) {
        ++pass; disabled.clear();
        m_mbs.realize(s, Stage::Acceleration);
        for (unsigned i=0; i < candidates.m_friction.size(); ++i) {
            const int which = candidates.m_friction[i];
            MyFrictionElement& fric = updFrictionElement(which);
            if (!fric.isSticking(s)) continue;

            const Real f = fric.calcStictionForceWitness(s);
            if (f<0) {
                SimTK_DEBUG4("  Pass %d: disable stiction %d: force capacity=%g at N=%g\n", 
                    pass, which, f, fric.getMasterNormalForce(s));
                disabled.push_back(which);
                stictionDisabled.push_back(which); // for all passes
                fric.recordImpendingSlipInfo(s);
            }
        }
        if (disabled.empty())
            break; // no stiction is greater than mu_s*N now

        // Do the necessary disabling. (Invalidates Instance stage.)
        for (unsigned tbd=0; tbd < disabled.size(); ++tbd) {
            const int which = disabled[tbd];
            const MyFrictionElement& fric = getFrictionElement(which);
            fric.disableStiction(s); // saves recorded info
        }
        // Go back for another pass.
    }

    if (!stictionDisabled.empty()) {
        // Recalculate sliding friction from normal we just obtained.
        // TODO: this should be a Newton iteration.
        SimTK_DEBUG("REVISE NORMAL\n");
        const int NIters = 6;
        for (int i=0; i < NIters; ++i) {
            s.autoUpdateDiscreteVariables();
            s.invalidateAllCacheAtOrAbove(Stage::Dynamics);
            m_mbs.realize(s, Stage::Acceleration);
        }
    }

    // Now look for contact constraints that are "sucking". We'll try disabling
    // those but will re-enable if that would lead to penetration.
    for (unsigned i=0; i < candidates.m_contact.size(); ++i) {
        const int which = candidates.m_contact[i];
        const MyContactElement& cont = getContactElement(which);
        if (cont.isDisabled(s)) continue;
        const Real f = cont.getForce(s);
        if (f<0) {
            SimTK_DEBUG2("  consider disabling contact %d because force=%g\n", 
                            which, f);
            if (cont.hasFrictionElement()) {
                MyFrictionElement& fric = cont.updFrictionElement();
                if (fric.isSticking(s)) {
                    fric.recordImpendingSlipInfo(s);
                }
            }
            contactDisabled.push_back(which);
        }
    }

    // OK, now tentatively disable the pulling contacts. If there is an 
    // associated stiction element we have to disable it too.
    for (unsigned tbd=0; tbd < contactDisabled.size(); ++tbd) {
        const int which = contactDisabled[tbd];
        const MyContactElement& cont = getContactElement(which);
        if (cont.hasFrictionElement()) {
            const MyFrictionElement& fric = cont.getFrictionElement();
            if (fric.isSticking(s)) {
                fric.disableStiction(s);
                stictionDisabled.push_back(fric.getFrictionIndex());
            }
        }
        cont.disable(s);
    }

    m_mbs.realize(s, Stage::Acceleration);
    if (!contactDisabled.empty()) {
        // Recalculate sliding friction from normal we just obtained.
        // TODO: this should be a Newton iteration.
        SimTK_DEBUG("REVISE NORMAL\n");
        const int NIters = 6;
        for (int i=0; i < NIters; ++i) {
            s.autoUpdateDiscreteVariables();
            s.invalidateAllCacheAtOrAbove(Stage::Dynamics);
            m_mbs.realize(s, Stage::Acceleration);
        }
    }

    // See which of the just-disabled constraints is violated. We'll re-enable
    // those but won't turn stiction back on if it was engaged before.
    Array_<int> violated;
    for (unsigned tbd=0; tbd < contactDisabled.size(); ++tbd) {
        const int which = contactDisabled[tbd];
        const MyContactElement& cont = getContactElement(which);
        const Real aerr = cont.getAerr(s);
        if (aerr < 0) {
            violated.push_back(which);
            SimTK_DEBUG2("  RE-ENABLE contact %d cuz aerr=%g\n", 
                            which, aerr);
        }
    }

    // Re-enable now.
    for (unsigned v=0; v < violated.size(); ++v) {
        const int which = violated[v];
        const MyContactElement& cont = getContactElement(which);
        cont.enable(s);
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

    SimTK_DEBUG2("... Done; %d contacts, %d stictions broken.\n", 
        contactDisabled.size()-violated.size(), stictionDisabled.size());

    showConstraintStatus(s, "EXIT findActiveCandidates()");
}

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
// speed |v|<=tol, or if it has reversed direction, we consider it unreliable 
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
//      shouldUpdate(v) = ~v*d_prev > 0 && |v|>tol
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
// There are two sliding-to-stiction event witness functions
//              e1(v)=dot(v, d_prev)    velocity reversal
//              e2(v)=|v|-vtol          slowdown TODO ???
// OR: e(v) = dot(v, d_prev) - vtol ?? [signed]
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
    MySlidingFrictionForceImpl(const GeneralForceSubsystem& forces,
                               const MyPointContactFriction& ptFriction)
    :   m_forces(forces), m_friction(ptFriction), 
        m_contact(ptFriction.getMyPointContact())
    {}

    bool hasSlidingForce(const State& s) const 
    {   return m_friction.isMasterActive(s) && !m_friction.isSticking(s); }

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

    // Calculate d_eff, the direction to be opposed by the sliding force.
    Vec2 getEffectiveSlipDir(const State& s) const {
        const Vec2 Vslip = m_contact.getSlipVelocity(s);
        const Vec2 prevVslipDir = getPrevSlipDir(s);
        if (shouldUpdate(Vslip, prevVslipDir, 1e-2)) { // TODO: tol?
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

        if (shouldUpdate(Vslip, prevVslipDir, 1e-2)) { // TODO: tol?
            Vec2& prevSlipUpdate = Value<Vec2>::updDowncast
               (m_forces.updDiscreteVarUpdateValue(state, m_prevSlipDirIx));
            const Real v = Vslip.norm();
            const Vec2 slipDir = Vslip / v;
            prevSlipUpdate = slipDir;
            m_forces.markDiscreteVarUpdateValueRealized(state, m_prevSlipDirIx);

            #ifndef NDEBUG
            printf("UPDATE %d: prevSlipDir=%g %g; now=%g %g; |v|=%g dot=%g vdot=%g\n",
                m_friction.getFrictionIndex(),
                prevVslipDir[0],prevVslipDir[1],slipDir[0],slipDir[1],
                v, ~slipDir*prevVslipDir, ~Vslip*prevVslipDir);
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

    const GeneralForceSubsystem&    m_forces;
    const MyPointContactFriction&   m_friction;
    const MyPointContact&           m_contact;

    mutable DiscreteVariableIndex   m_prevNix;       // previous normal force
    mutable DiscreteVariableIndex   m_prevSlipDirIx; // previous slip direction
};

// This is the force handle's constructor; it just creates the force
// implementation object.
MySlidingFrictionForce::MySlidingFrictionForce
   (GeneralForceSubsystem& forces,
    const MyPointContactFriction& friction) 
:   Force::Custom(forces, new MySlidingFrictionForceImpl(forces,friction)) 
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
