/* -------------------------------------------------------------------------- *
 *                Simbody(tm) - UnilateralPointContact Example                *
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
const Real RunTime=20;

//==============================================================================
//                           MY UNILATERAL CONSTRAINT
//==============================================================================
// This abstract class hides the details about which kind of constraint
// we're dealing with, while giving us enough to work with for deciding what's
// on and off and generating impulses.
//
// There is always a scalar associated with the constraint for making 
// decisions, although contact constraints may also have some additional
// constraint equations for stiction.
class MyUnilateralConstraint {
public:
    enum ImpulseType {Compression,Expansion,Capture};

    MyUnilateralConstraint(Constraint& uni, Real multSign, Real coefRest) 
    :   m_uni(uni), m_multSign(multSign), m_coefRest(coefRest), 
        m_restitutionDone(false) 
    {   m_uni.setDisabledByDefault(true); }

    virtual ~MyUnilateralConstraint() {}

    // These must be constructed so that a negative value means the 
    // unilateral constraint condition is violated.
    virtual Real getPerr(const State& state) const = 0;
    virtual Real getVerr(const State& state) const = 0;
    virtual Real getAerr(const State& state) const = 0;

    // This returns a point in the ground frame at which you might want to
    // say the constraint is "located", for purposes of display. This should
    // return something useful even if the constraint is currently off.
    virtual Vec3 whereToDisplay(const State& state) const = 0;

    // Returns zero if the constraint is not currently enabled.
    Real getForce(const State& s) const {
        if (isDisabled(s)) return 0;
        const Vector mult = m_uni.getMultipliersAsVector(s);
        assert(mult.size() == 1);
        return m_multSign*mult[0];
    }

    // Override these if you have auxiliary constraints but be sure to 
    // invoke superclass method too.

    virtual void enable(State& state) const {m_uni.enable(state);}
    virtual void disable(State& state) const {m_uni.disable(state);}

    virtual void setMyDesiredDeltaV(const State&    s,
                                    Vector&         desiredDeltaV) const
    {   Vector myDesiredDV(1); myDesiredDV[0] = m_multSign*getVerr(s);
        m_uni.setMyPartInConstraintSpaceVector(s, myDesiredDV, 
                                                   desiredDeltaV); }

    virtual void recordImpulse(ImpulseType type, const State& state,
                               const Vector& lambda) {
        Vector myImpulse(1);
        m_uni.getMyPartFromConstraintSpaceVector(state, lambda, myImpulse);
        const Real I = myImpulse[0];
        if (type==Compression) m_Ic = I;
        else if (type==Expansion) m_Ie = I;
        m_I += I;
    }

    // This is used by some constraints to collect position information that
    // may be used later to set instance variables when enabling the underlying
    // Simbody constraint. All constraints zero impulses here.
    virtual void initializeForImpact(const State& state) 
    {   setRestitutionDone(false); m_Ic = m_Ie = m_I = 0; }

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

    bool isDisabled(const State& state) const 
    {   return m_uni.isDisabled(state); }

    Real getCoefRest() const {return m_coefRest;}
    void setRestitutionDone(bool isDone) {m_restitutionDone=isDone;}
    bool isRestitutionDone() const {return m_restitutionDone;}

protected:
    Constraint      m_uni;
    const Real      m_multSign; // 1 or -1
    const Real      m_coefRest;

    // Runtime
    bool m_restitutionDone;
    Real m_Ic, m_Ie, m_I; // impulses
};

//==============================================================================
//                            SHOW CONTACT
//==============================================================================
// For each visualization frame, generate ephemeral geometry to show just 
// during this frame, based on the current State.
class ShowContact : public DecorationGenerator {
public:
    ShowContact(const Array_<MyUnilateralConstraint*>& unis) 
    :   m_unis(unis) {}

    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        for (unsigned i=0; i < m_unis.size(); ++i) {
            const MyUnilateralConstraint& uni = *m_unis[i];
            const Vec3 loc = uni.whereToDisplay(state);
            if (!uni.isDisabled(state)) {
                geometry.push_back(DecorativeSphere(.5)
                    .setTransform(loc)
                    .setColor(Red).setOpacity(.25));
                geometry.push_back(DecorativeText("LOCKED")
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
    const Array_<MyUnilateralConstraint*>& m_unis;
};



//==============================================================================
//                               STATE SAVER
//==============================================================================
// This reporter is called now and again to save the current state so we can
// play back a movie at the end.
class StateSaver : public PeriodicEventReporter {
public:
    StateSaver(const MultibodySystem&                   system,
               const Array_<MyUnilateralConstraint*>&   unis,
               const Integrator&                        integ,
               Real                                     reportInterval)
    :   PeriodicEventReporter(reportInterval), 
        m_system(system), m_unis(unis), m_integ(integ) 
    {   m_states.reserve(2000); }

    ~StateSaver() {}

    void clear() {m_states.clear();}
    int getNumSavedStates() const {return (int)m_states.size();}
    const State& getState(int n) const {return m_states[n];}

    void handleEvent(const State& s) const {
        const SimbodyMatterSubsystem& matter=m_system.getMatterSubsystem();
        const SpatialVec PG = matter.calcSystemMomentumAboutGroundOrigin(s);

#ifndef NDEBUG
        printf("%3d: %5g mom=%g,%g E=%g", m_integ.getNumStepsTaken(),
            s.getTime(),
            PG[0].norm(), PG[1].norm(), m_system.calcEnergy(s));
        cout << " Triggers=" << s.getEventTriggers() << endl;
        for (unsigned i=0; i < m_unis.size(); ++i) {
            const MyUnilateralConstraint& uni = *m_unis[i];
            const bool isLocked = !uni.isDisabled(s);
            printf("  Uni constraint %d is %s, h=%g dh=%g\n", i, 
                   isLocked?"LOCKED":"unlocked", uni.getPerr(s),uni.getVerr(s));
            if (isLocked) {
                m_system.realize(s, Stage::Acceleration);
                cout << "    force=" << uni.getForce(s) << endl;
            } 
        }
#endif

        m_states.push_back(s);
    }
private:
    const MultibodySystem&                  m_system;
    const Array_<MyUnilateralConstraint*>&  m_unis;
    const Integrator&                       m_integ;
    mutable Array_<State>                   m_states;
};



//==============================================================================
//                          CONTACT ON HANDLER
//==============================================================================

class ContactOn: public TriggeredEventHandler {
public:
    ContactOn(const MultibodySystem&                    system,
              const Array_<MyUnilateralConstraint*>&    unis,
              unsigned                                  which,
              Stage                                     stage) 
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
        const MyUnilateralConstraint& uni = *m_unis[m_which];
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


    // Make a list of all the unilateral constraints that could conceivably
    // receive an impulse. Any constraint that is currently enabled, or any
    // currently disabled constraint that is within posTol of contact is 
    // included.
    void findProximalConstraints(const State&       state,
                                 Real               posTol,
                                 Array_<int>&       proximal) const;



    // Given the set of proximal constraints, prevent penetration by applying
    // a nonnegative least squares impulse generating a step change in 
    // velocity. On return, the applied impulse and new velocities are recorded
    // in proximal, and state is updated to the new velocities and realized
    // through Velocity stage. Constraints that were stopped are enabled, those
    // that rebounded are disabled.
    void processCompressionPhase(Array_<int>&   proximal,
                                 State&         state) const;

    // Given a solution to the compression phase, including the compression
    // impulse, the set of impacters (enabled) and rebounders (disabled and
    // with positive rebound velocity), apply an expansion impulse based on
    // the effective coefficients of restitution of the impacters. Wherever
    // restitution is applied, the effective coefficient is reset to zero so
    // that further restitution will not be done for that contact. Returns
    // true if any expansion was done; otherwise nothing has changed.
    // Expansion may result in some negative velocities, in which case it has
    // induced further compression so another compression phase is required.
    bool processExpansionPhase(Array_<int>&     proximal,
                               State&           state) const;

    // Examine the rebounders to see if any are rebounding with a speed at or
    // below the capture velocity. If so, enable those constraints and apply a
    // (hopefully small) negative impulse to eliminate that rebound velocity.
    // Repeat if that induces any negative velocities or any further slow
    // rebounders. This terminates will all rebounders leaving with velocities
    // greater than vCapture, or else all constraints are enabled.
    void captureSlowRebounders(Real             vCapture,
                               Array_<int>&     proximal,
                               State&           state) const;

    // This method is used at the start of compression phase to modify any
    // constraint parameters as necessary, and then enable all the proximal
    // constraints. Some or all of these will be disabled during the impact
    // analysis in compression or expansion phases. On return the state has
    // been updated and realized through Instance stage.
    void enableAllProximalConstraints(Array_<int>&  proximal,
                                      State&        state) const;

    // Given only the subset of proximal constraints that are active, calculate
    // the impulse that would eliminate all their velocity errors. No change is
    // made to the set of active constraints. Some of the resulting impulses
    // may be negative.
    void calcStoppingImpulse(const Array_<int>&     proximal,
                             const State&           state,
                             Vector&                lambda0) const;

    // Given the initial generalized speeds u0, and a constraint-space impulse
    // lambda, calculate the resulting step velocity change du, modify the
    // generalized speeds in state to u0+du, and realize Velocity stage.
    void updateVelocities(const Vector& u0, 
                          const Vector& lambda, 
                          State&        state) const;


private:
    const MultibodySystem&                  m_mbs; 
    const Array_<MyUnilateralConstraint*>&  m_unis;
    const unsigned                          m_which;
    const Stage                             m_stage;
};



//==============================================================================
//                          CONTACT OFF HANDLER
//==============================================================================
// Allocate one of these for each unilateral constraint. This handler takes
// care of disabling an active constraint when its contact force crosses zero
// from positive to negative.
class ContactOff: public TriggeredEventHandler {
public:
    ContactOff(const MultibodySystem&               system,
        const Array_<MyUnilateralConstraint*>&      unis,
        unsigned                                    which) 
    :   TriggeredEventHandler(Stage::Acceleration), 
        m_mbs(system), m_unis(unis), m_which(which)
    { 
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    // This is the witness function.
    Real getValue(const State& state) const {
        const MyUnilateralConstraint& uni = *m_unis[m_which];
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

        disablePullingContacts(m_mbs,s,m_unis);

        SimTK_DEBUG("LIFTOFF DONE.\n");
        SimTK_DEBUG("----------------------------------------------------\n");
    }

    // This is also used by ContactOn at the end.
    static void disablePullingContacts
       (const MultibodySystem& mbs, State& s, 
        const Array_<MyUnilateralConstraint*>& unis); 

private:
    const MultibodySystem&                  m_mbs; 
    const Array_<MyUnilateralConstraint*>&  m_unis;
    const unsigned                          m_which; // one of the unis
};



//==============================================================================
//                             MY POINT CONTACT
//==============================================================================
// Define a unilateral constraint to represent contact of a point on a moving
// body with the ground plane. The ground normal is assumed to be +y. This
// contact constraint has "super friction" that always sticks if it contacts
// at all. Note: that can generate non-physical effects.
class MyPointContact : public MyUnilateralConstraint {
    typedef MyUnilateralConstraint Super;
public:
    MyPointContact(MobilizedBody& body, const Vec3& point,
                 Real coefRest)
    :   MyUnilateralConstraint
           (Constraint::PointInPlane(updGround(body), UnitVec3(YAxis), Zero,
                                     body, point),
             Real(-1), // multiplier sign
             coefRest),
        m_body(body), m_point(point), m_groundPoint(0),
        m_noslipX(updGround(body), Vec3(0), UnitVec3(XAxis), 
                  updGround(body), body),
        m_noslipZ(updGround(body), Vec3(0), UnitVec3(ZAxis), 
                  updGround(body), body)
    {
        m_noslipX.setDisabledByDefault(true);
        m_noslipZ.setDisabledByDefault(true);
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

    Vec3 whereToDisplay(const State& state) const OVERRIDE_11 {
        return m_body.findStationLocationInGround(state,m_point);
    }

    void recordImpulse(ImpulseType type, const State& state,
                      const Vector& lambda) OVERRIDE_11
    {
        Super::recordImpulse(type, state, lambda);

        // Record translational impulse.
        Vector myImpulseX(1), myImpulseZ(1);
        m_noslipX.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseX);
        m_noslipZ.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseZ);
        const Vec2 tI(myImpulseX[0], myImpulseZ[0]);
        if (type==Compression) m_tIc = tI;
        else if (type==Expansion) m_tIe = tI;
        m_tI += tI;
    }

    void setMyDesiredDeltaV(const State& s,
                            Vector& desiredDeltaV) const OVERRIDE_11
    {
        Super::setMyDesiredDeltaV(s, desiredDeltaV);
        const Vec3 dv = 
            m_multSign*m_body.findStationVelocityInGround(s, m_point);
        Vector myDesiredDV(1); // Nuke translational velocity also.
        myDesiredDV[0] = dv[XAxis];
        m_noslipX.setMyPartInConstraintSpaceVector(s, myDesiredDV, desiredDeltaV);
        myDesiredDV[0] = dv[ZAxis];
        m_noslipZ.setMyPartInConstraintSpaceVector(s, myDesiredDV, desiredDeltaV);
    }

    // Note that recordStartingLocation() must have been called first.
    void enable(State& state) const OVERRIDE_11 {
        Super::enable(state);
        m_noslipX.setContactPoint(state, m_groundPoint);
        m_noslipZ.setContactPoint(state, m_groundPoint);
        m_noslipX.enable(state); m_noslipZ.enable(state);
    }
    void disable(State& state) const OVERRIDE_11 {
        Super::disable(state);
        m_noslipX.disable(state); m_noslipZ.disable(state);
    }

    // Set the ground point to be the projection of the follower point
    // onto the ground plane. This will be used the next time this constraint
    // is enabled.
    void initializeForImpact(const State& s) OVERRIDE_11 {
        Super::initializeForImpact(s);
        const Vec3 p = m_body.findStationLocationInGround(s, m_point);
        m_groundPoint = p;
        m_tI = m_tIe = m_tIc = Vec2(0);
    }
private:
    MobilizedBody& updGround(MobilizedBody& body) const {
        SimbodyMatterSubsystem& matter = body.updMatterSubsystem();
        return matter.updGround();
    }

    const MobilizedBody&    m_body;
    const Vec3              m_point;
    Constraint::NoSlip1D    m_noslipX, m_noslipZ;

    // Runtime
    Vec3 m_groundPoint;
    Vec2 m_tIc; // most recent tangential compression impulse
    Vec2 m_tIe; // most recent tangential expansion impulse
    Vec2 m_tI;  // accumulated tangential impulse
};


//==============================================================================
//                                  MY STOP
//==============================================================================
// Define a unilateral constraint to represent a joint stop that limits
// the allowable motion of a single generalized coordinate. You can specify
// a coefficient of restitution and whether the given limit is the upper or
// lower limit.
class MyStop : public MyUnilateralConstraint {
public:
    enum Side {Lower,Upper};
    MyStop(Side side, MobilizedBody& body, int whichQ,
         Real limit, Real coefRest)
    :   MyUnilateralConstraint
           (Constraint::ConstantSpeed(body, MobilizerUIndex(whichQ), Real(0)), 
            Real(side==Lower?-1:1), coefRest),
        m_body(body), m_whichq(whichQ), m_whichu(whichQ),
        m_sign(side==Lower?1.:-1.), m_limit(limit)
    {}

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
class MyRope : public MyUnilateralConstraint {
public:
    MyRope(MobilizedBody& body1, const Vec3& pt1,
           MobilizedBody& body2, const Vec3& pt2, Real d,
           Real coefRest)
    :   MyUnilateralConstraint
           (Constraint::Rod(body1, pt1, body2, pt2, d), Real(1), coefRest),
        m_body1(body1), m_point1(pt1), m_body2(body2), m_point2(pt2), m_dist(d)
    {}


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
    const Real CoefRest = 0.8;
    Array_<MyUnilateralConstraint*>     unis;

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
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(CubeHalfDims);
        unis.push_back(new MyPointContact(cube, pt, CoefRest));
    }

    unis.push_back(new MyRope(Ground, Vec3(-5,10,0),
                              cube, Vec3(-CubeHalfDims[0],-CubeHalfDims[1],0), 
                              5., .5*CoefRest));
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
        unis.push_back(new MyPointContact(weight, pt, 0.5*CoefRest));
    }
    unis.push_back(new MyStop(MyStop::Upper,weight,0, Pi/9,0.1*CoefRest));
    unis.push_back(new MyStop(MyStop::Lower,weight,0, -Pi/9,0.1*CoefRest));

//#endif
//#ifdef NOTDEF
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
        unis.push_back(new MyPointContact(weight2, pt, CoefRest));
    }
//#endif

    Visualizer viz(mbs);
    viz.setShowSimTime(true);
    viz.addDecorationGenerator(new ShowContact(unis));
    mbs.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));

    //ExplicitEulerIntegrator integ(mbs);
    //CPodesIntegrator integ(mbs,CPodes::BDF,CPodes::Newton);
    //RungeKuttaFeldbergIntegrator integ(mbs);
    Real accuracy = 1e-2;
    RungeKuttaMersonIntegrator integ(mbs);
    //RungeKutta3Integrator integ(mbs);
    //VerletIntegrator integ(mbs);
    integ.setAccuracy(accuracy);
    //integ.setAllowInterpolation(false);
    integ.setMaximumStepSize(0.1);

    StateSaver* stateSaver = new StateSaver(mbs,unis,integ,ReportInterval);
    mbs.addEventReporter(stateSaver);

    for (unsigned i=0; i < unis.size(); ++i) {
        mbs.addEventHandler(new ContactOn(mbs, unis,i, Stage::Position));
        mbs.addEventHandler(new ContactOn(mbs, unis,i, Stage::Velocity));
        mbs.addEventHandler(new ContactOn(mbs, unis,i, Stage::Acceleration));
    }

    for (unsigned i=0; i < unis.size(); ++i)
        mbs.addEventHandler(new ContactOff(mbs, unis,i));
  
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

    Array_<int> enableThese;
    for (unsigned i=0; i < unis.size(); ++i) {
        const Real perr = unis[i]->getPerr(s);
        printf("uni constraint %d has perr=%g%s\n", i, perr,
            perr<=0?" (ENABLING)":"");
        if (perr <= 0)
            enableThese.push_back(i);
    }

    cout << "Initial configuration shown. Next? ";
    getchar();

    for (unsigned i=0; i < enableThese.size(); ++i)
        unis[enableThese[i]]->enable(s);

    Assembler(mbs).assemble(s);
    viz.report(s);
    cout << "Assembled configuration shown. Ready? ";
    getchar();
    
    // Simulate it.

    TimeStepper ts(mbs, integ);
    ts.initialize(s);

    const double startReal = realTime();
    const double startCPU = cpuTime();

    ts.stepTo(RunTime);

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
void ContactOn::
handleEvent(State& s, Real accuracy, bool& shouldTerminate) const 
{
    const Real VCapture=1e-2;

    Array_<int> proximal;
    findProximalConstraints(s, accuracy, proximal);

    SimTK_DEBUG4("\nIMPACT (%s) for uni constraint %d at t=%.16g; %d proximal\n", 
        m_stage.getName().c_str(), m_which, s.getTime(), proximal.size());

    bool needMoreCompression = true;
    while (needMoreCompression) {
        processCompressionPhase(proximal, s);
        needMoreCompression = false;

        if (processExpansionPhase(proximal, s)) {
            for (unsigned i=0; i<proximal.size(); ++i) {
                const MyUnilateralConstraint& uni = *m_unis[proximal[i]];
                if (uni.getVerr(s) < 0) {
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

    // Make sure all enabled position and velocity constraints 
    // are satisfied.
    m_mbs.project(s, accuracy);

    // Finally, evaluate accelerations and reaction forces and check if 
    // any of the active contacts are generating negative ("pulling") 
    // forces; if so, inactivate them.
    ContactOff::disablePullingContacts(m_mbs, s, m_unis);

#ifndef NDEBUG
    printf("END OF IMPACT for %d proximal constraints:\n",proximal.size());
    for (unsigned i=0; i < proximal.size(); ++i) {
        const int which = proximal[i];
        const MyUnilateralConstraint& uni = *m_unis[which];
        printf("  %d %3s: I=%g H=%g V=%g A=%g F=%g\n",
            which, uni.isDisabled(s) ? "off" : "ON", 
            uni.getImpulse(), uni.getPerr(s), uni.getVerr(s), 
            uni.getAerr(s), uni.getForce(s));       
    }
    printf("DONE WITH IMPACT.\n\n");
#endif
}



//------------------------ FIND PROXIMAL CONSTRAINTS ---------------------------
void ContactOn::
findProximalConstraints(const State&       s,
                        Real               posTol,
                        Array_<int>&       proximal) const
{
    m_mbs.realize(s, Stage::Position);

    proximal.clear();
    for (unsigned i=0; i<m_unis.size(); ++i) {
        MyUnilateralConstraint& uni = *m_unis[i];
        if (!uni.isDisabled(s) || uni.getPerr(s) <= posTol) 
        {
            uni.initializeForImpact(s);
            proximal.push_back(i);
        }
    }
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
processCompressionPhase(Array_<int>&    proximal,
                        State&          s) const
{
    SimTK_DEBUG("Entering processCompressionPhase() ...\n");
    Vector lambda0, lambdaTry;
    Array_<int> maybeDisabled, recapturing, turnOffStiction;

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
        for (unsigned i=0; i < proximal.size(); ++i) {
            const int which = proximal[i];
            const MyUnilateralConstraint& uni = *m_unis[which];
            if (uni.isDisabled(s))
                continue;
            Real imp = uni.getMyValueFromConstraintSpaceVector(s, lambda0);
            if (imp > 0)
                continue;
            maybeDisabled.push_back(which);
            SimTK_DEBUG2("  uni constraint %d has negative compression impulse=%g\n",
                which, imp);
        }
        if (maybeDisabled.empty())
            break;

        // Disable the candidates, then see if they rebound.
        for (unsigned d=0; d < maybeDisabled.size(); ++d)
            m_unis[maybeDisabled[d]]->disable(s);
        m_mbs.realize(s, Stage::Instance);
        calcStoppingImpulse(proximal, s, lambdaTry);
        updateVelocities(u0, lambdaTry, s);

        recapturing.clear();
        for (unsigned i=0; i<maybeDisabled.size(); ++i) {
            const int which = maybeDisabled[i];
            const MyUnilateralConstraint& uni = *m_unis[which];
            const Real newV = uni.getVerr(s);           
            SimTK_DEBUG2("  candidate uni constraint %d would have v=%g\n",
                   which, newV);
            if (newV <= 0) {
                recapturing.push_back(which);
                SimTK_DEBUG2("  RECAPTURING uni constraint %d with v=%g\n", 
                    which, newV);
            }
        }

        // Re-enable the recaptured candidates.
        if (!recapturing.empty()) {
            for (unsigned c=0; c < recapturing.size(); ++c)
                m_unis[recapturing[c]]->enable(s);
            m_mbs.realize(s, Stage::Instance);
        }

        const int numDisabled = maybeDisabled.size()-recapturing.size();
        if (numDisabled == 0) {
            SimTK_DEBUG("  None of the candidates was actually disabled.\n");
            // lambda0 is still correct
            break;
        }

        if (recapturing.empty()) lambda0 = lambdaTry;
        else calcStoppingImpulse(proximal, s, lambda0);

        SimTK_DEBUG1("  RETRY with %d constraints disabled\n", numDisabled);
    }
    updateVelocities(u0, lambda0, s);

    // Now update the entries for each proximal constraint to reflect the
    // compression impulse and post-compression velocity.
    SimTK_DEBUG("  Compression results:\n");
    for (unsigned i=0; i < proximal.size(); ++i) {
        const int which = proximal[i];
        MyUnilateralConstraint& uni = *m_unis[which];
        if (!uni.isDisabled(s))
            uni.recordImpulse(MyUnilateralConstraint::Compression, s, lambda0);
        SimTK_DEBUG4("  %d %3s: Ic=%g, V=%g\n",
            which, uni.isDisabled(s) ? "off" : "ON", 
            uni.getCompressionImpulse(), uni.getVerr(s));
    }

    SimTK_DEBUG("... compression phase done.\n");
}



//------------------------- PROCESS EXPANSION PHASE ----------------------------
bool ContactOn::
processExpansionPhase(Array_<int>&  proximal,
                      State&        s) const
{
    SimTK_DEBUG("Entering processExpansionPhase() ...\n");

    // Generate an expansion impulse if there were any active contacts that
    // still have some restitution remaining.
    Vector expansionImpulse;

    bool anyChange = false;
    for (unsigned i=0; i<proximal.size(); ++i) {
        const int which = proximal[i];
        MyUnilateralConstraint& uni = *m_unis[which];
        if (uni.isDisabled(s)||uni.isRestitutionDone()||uni.getCoefRest()==0
            ||uni.getCompressionImpulse()<=0)
            continue;
        uni.setMyExpansionImpulse(s, uni.getCoefRest(), expansionImpulse);
        uni.recordImpulse(MyUnilateralConstraint::Expansion,s,expansionImpulse);
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
    for (unsigned i=0; i < proximal.size(); ++i) {
        const int which = proximal[i];
        const MyUnilateralConstraint& uni = *m_unis[which];
        if (!uni.isDisabled(s) && uni.getVerr(s) > 0)
            toDisable.push_back(which);
    }

    // Now do the actual disabling (can't mix this with checking velocities)
    // because disabling invalidates Instance stage.
    for (unsigned i=0; i < toDisable.size(); ++i) {
        const int which = toDisable[i];
        const MyUnilateralConstraint& uni = *m_unis[which];
        uni.disable(s);
    }

    SimTK_DEBUG("  Expansion results:\n");
    m_mbs.realize(s, Stage::Velocity);
    for (unsigned i=0; i < proximal.size(); ++i) {
        const int which = proximal[i];
        const MyUnilateralConstraint& uni = *m_unis[which];
        SimTK_DEBUG4("  %d %3s: Ie=%g, V=%g\n",
            which, uni.isDisabled(s) ? "off" : "ON", 
            uni.getExpansionImpulse(), uni.getVerr(s));
    }

    SimTK_DEBUG("... expansion phase done.\n");

    return true;
}


//------------------------- CAPTURE SLOW REBOUNDERS ----------------------------
void ContactOn::
captureSlowRebounders(Real          vCapture,
                      Array_<int>&  proximal,
                      State&        s) const
{
    // Capture any rebounder whose velocity is too slow.
    SimTK_DEBUG("Entering captureSlowRebounders() ...\n");

    int nCaptured=0, nPasses=0;
    while (true) {
        ++nPasses;
        SimTK_DEBUG1("  start capture pass %d:\n", nPasses);
        Array_<int> toCapture;
        for (unsigned i=0; i < proximal.size(); ++i) {
            const int which = proximal[i];
            MyUnilateralConstraint& uni = *m_unis[which];
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
            MyUnilateralConstraint& uni = *m_unis[which];
            uni.enable(s);
        }

        m_mbs.realize(s, Stage::Velocity);
        Vector captureImpulse;
        calcStoppingImpulse(proximal, s, captureImpulse);
        updateVelocities(Vector(), captureImpulse, s);

        // Now update the entries for each proximal constraint to reflect the
        // capture impulse and post-capture velocity.
        for (unsigned i=0; i < proximal.size(); ++i) {
            const int which = proximal[i];
            MyUnilateralConstraint& uni = *m_unis[which];
            if (!uni.isDisabled(s))
                uni.recordImpulse(MyUnilateralConstraint::Capture, 
                                  s, captureImpulse);
        }
    }
}



//---------------------- ENABLE ALL PROXIMAL CONSTRAINTS -----------------------
void ContactOn::
enableAllProximalConstraints(Array_<int>&  proximal,
                             State&        state) const
{
    // Set the contact point and enable the constraints.
    for (unsigned i=0; i < proximal.size(); ++i) {
        const int which = proximal[i];
        const MyUnilateralConstraint& uni = *m_unis[which];
        if (uni.isDisabled(state))
            uni.enable(state);
    }
    m_mbs.realize(state, Stage::Instance);
}



//------------------------- CALC COMPRESSION IMPULSE ---------------------------
// Calculate the impulse that eliminates all residual velocity for the
// current set of enabled constraints.
void ContactOn::
calcStoppingImpulse(const Array_<int>&    proximal,
                    const State&          s,
                    Vector&               lambda0) const
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(s, Stage::Dynamics); // TODO: should only need Position
    Vector desiredDeltaV;  // in constraint space
    SimTK_DEBUG("  Entering calcStoppingImpulse() ...\n");
    bool gotOne = false;
    for (unsigned i=0; i < proximal.size(); ++i) {
        const int which = proximal[i];
        const MyUnilateralConstraint& uni = *m_unis[which];
        if (uni.isDisabled(s))
            continue;
        SimTK_DEBUG2("    uni constraint %d enabled, v=%g\n",
            which, uni.getVerr(s));
        uni.setMyDesiredDeltaV(s, desiredDeltaV);
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
//                               CONTACT OFF
//==============================================================================

//-------------------------- DISABLE PULLING CONTACTS --------------------------
// This method checks the active contacts to see if any of them are generating
// "pulling" forces. In that case the contact is disabled unless that would
// lead to a negative acceleration.
// This is invoked by the ContactOff handler, and as the last step of the
// ContactOn (impact) handler.
// TODO: need to search for a consistent set of active contraints.

/*static*/ void ContactOff::disablePullingContacts
   (const MultibodySystem& mbs, State& s, 
    const Array_<MyUnilateralConstraint*>& unis) 
{
    SimTK_DEBUG("Entering disablePullingContacts() ...\n");

    mbs.realize(s, Stage::Acceleration);
    // Check first, disable later because we don't want to invalidate
    // the reaction forces in the state yet.
    Array_<int> toBeDisabled;
    for (unsigned i=0; i < unis.size(); ++i) {
        const MyUnilateralConstraint& uni = *unis[i];
        if (uni.isDisabled(s)) continue;
        const Real f = uni.getForce(s);
        if (f<0) {
            SimTK_DEBUG2("  consider disabling uni %d because force=%g\n", 
                            i, f);
            toBeDisabled.push_back(i);
        }
    }

    // OK, now tentatively disable the pulling contacts.
    for (unsigned tbd=0; tbd < toBeDisabled.size(); ++tbd) {
        const int i = toBeDisabled[tbd];
        const MyUnilateralConstraint& uni = *unis[i];
        uni.disable(s);
    }

    // Now see which of the disabled constraints is violated.
    mbs.realize(s, Stage::Acceleration);
    Array_<int> violated;
    for (unsigned p=0; p < toBeDisabled.size(); ++p) {
        const int which = toBeDisabled[p];
        const MyUnilateralConstraint& uni = *unis[which];
        const Real aerr = uni.getAerr(s);
        if (aerr < 0) {
            violated.push_back(which);
            SimTK_DEBUG2("  RE-ENABLE constraint %d cuz aerr=%g\n", 
                            which, aerr);
        }
    }

    // Re-enable now.
    for (unsigned v=0; v < violated.size(); ++v) {
        const int which = violated[v];
        const MyUnilateralConstraint& uni = *unis[which];
        uni.enable(s);
    }

    // Always leave at acceleration stage.
    mbs.realize(s, Stage::Acceleration);

    SimTK_DEBUG1("... Done; %d contacts broken.\n", 
        toBeDisabled.size()-violated.size());
}
