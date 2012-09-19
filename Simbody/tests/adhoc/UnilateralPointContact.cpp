/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
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
This example shows a manual approach to simple unilateral contact between
designated points on moving bodies and a ground plane. We'll use Simbody
bilateral constraints with manual switching conditions.

For each designated contact point, we'll track the height over the ground
plane and use that as a switching ("witness") function to trigger an event
that enables the constraint. Then, we'll track the sign of the normal
reaction force and use it as a witness to disable the constraint.

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
In this example each contact consists of a constraint that prevents translation
of a point on a moving body with respect to a point on the ground plane. We
implement that here with Simbody's "Ball" constraint, so called because it
acts as though there were a ball joint between the bodies acting at the
contact point. We enable this constraint when a contact begins, defined so 
that its multipliers are the x,y,z components of the reaction force, with +y
being the ground plane normal. We monitor the reaction force y component, and
declare the contact broken if that component is negative.

How we handle impacts
---------------------
In this example, an impact is signaled by a contact point that reaches the 
ground plane with a negative vertical speed vy. This requires a step change to
the system velocities to avoid penetration. We achieve this step change by
applying a constraint-space
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

5) We are done with the impact. We now have inactive proximal contacts where
all are rebounding with velocities greater than capture velocity vc. However,
some of those may have heights that are slightly negative; this might cause
us to miss an impact if this contact reverses fast. So we calculate how much
time it would take for the rebound velocity to pull the worst-case rebounder
out of the ground, and take an explicit Euler step of that length so that
the heights of all rebounders are > 0.

6) Now calculate accelerations. If any of the active proximal contacts 
generate a zero or negative vertical reaction force they should be disabled;
otherwise we would miss the next break-free event. TODO: this can cause a 
problem if height<=0 since we'll miss the next impact if it doesn't go up
first.
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
//                            SHOW CONTACT
//==============================================================================
// For each visualization frame, generate ephemeral geometry to show just 
// during this frame, based on the current State.
class ShowContact : public DecorationGenerator {
public:
    ShowContact(const Array_<Constraint::PointInPlane>& contact) 
    :   m_contact(contact)
    {
    }
    void generateDecorations(const State& state, 
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        for (unsigned i=0; i < m_contact.size(); ++i) {
            const Constraint::PointInPlane& contact = m_contact[i];
            const MobodIndex mxP = contact.getPlaneMobilizedBodyIndex();
            const MobodIndex mxF = contact.getFollowerMobilizedBodyIndex();
            const Mobod& P = contact.getMatterSubsystem().getMobilizedBody(mxP);
            const Mobod& F = contact.getMatterSubsystem().getMobilizedBody(mxF);
            const UnitVec3& n_P  = contact.getDefaultPlaneNormal(); // +y
            const Real      h    = contact.getDefaultPlaneHeight(); // 0
            const Vec3      pt_F = contact.getDefaultFollowerPoint();
            const Vec3      pt_F_P = F.findStationLocationInAnotherBody(state, pt_F, P);
            const Real      fh = dot(pt_F_P, n_P); // follower height
            const Vec3      pt_P = pt_F_P - (fh-h)*n_P; 
            if (!contact.isDisabled(state)) {
                geometry.push_back(DecorativeSphere(.5)
                    .setTransform(pt_F)
                    .setColor(Red).setOpacity(.25)
                    .setBodyId(mxF));
                geometry.push_back(DecorativeText("LOCKED")
                    .setColor(White).setScale(.5)
                    .setTransform(pt_P+Vec3(0,.1,0))
                    .setBodyId(mxP));
            } else {
                geometry.push_back(DecorativeText(String(i))
                    .setColor(White).setScale(.5)
                    .setTransform(pt_F+Vec3(0,.1,0))
                    .setBodyId(mxF));
            }
        }
    }
private:
    const Array_<Constraint::PointInPlane>& m_contact;
};



//==============================================================================
//                               STATE SAVER
//==============================================================================
// This reporter is called now and again to save the current state so we can
// play back a movie at the end.
class StateSaver : public PeriodicEventReporter {
public:
    StateSaver(const MultibodySystem&                   system,
               const Array_<Constraint::PointInPlane>&  contact,
               const Integrator&                        integ,
               Real                                     reportInterval)
    :   PeriodicEventReporter(reportInterval), 
        m_system(system), m_contact(contact), m_integ(integ) 
    {   m_states.reserve(2000); }

    ~StateSaver() {}

    void clear() {m_states.clear();}
    int getNumSavedStates() const {return m_states.size();}
    const State& getState(int n) const {return m_states[n];}

    void handleEvent(const State& s) const {
        const SimbodyMatterSubsystem& matter=m_system.getMatterSubsystem();
        const SpatialVec PG = matter.calcSystemMomentumAboutGroundOrigin(s);

        printf("%3d: %5g mom=%g,%g E=%g", m_integ.getNumStepsTaken(),
            s.getTime(),
            PG[0].norm(), PG[1].norm(), m_system.calcEnergy(s));
        cout << " Triggers=" << s.getEventTriggers() << endl;
        for (unsigned i=0; i < m_contact.size(); ++i) {
            const Constraint::PointInPlane& contact = m_contact[i];

            const MobilizedBody& b2 = 
                matter.getMobilizedBody(contact.getFollowerMobilizedBodyIndex());
            Vec3 p2_G, v2_G;
            b2.findStationLocationAndVelocityInGround
                (s, contact.getDefaultFollowerPoint(), p2_G, v2_G);
            const bool isLocked = !contact.isDisabled(s);
            printf("  Constraint %d is %s, h=%g dh=%g\n", i, 
                   isLocked?"LOCKED":"unlocked", p2_G[YAxis], v2_G[YAxis]);
            if (isLocked) {
                m_system.realize(s, Stage::Acceleration);
                cout << "    lambda=" << contact.getMultiplier(s) << endl;
            } 
        }

        m_states.push_back(s);
    }
private:
    const MultibodySystem&                  m_system;
    const Array_<Constraint::PointInPlane>& m_contact;
    const Integrator&                       m_integ;
    mutable Array_<State,int>               m_states;
};



//==============================================================================
//                          CONTACT ON HANDLER
//==============================================================================

class CInfo;
class ContactOn: public TriggeredEventHandler {
public:
    ContactOn(const MultibodySystem&                    system,
              const Array_<Constraint::PointInPlane>&   contact,
              unsigned                                  which,
              const Array_<Real>&                       coefRest,
              const Array_<Real>&                       coefFric,
              const Array_<Constraint::NoSlip1D>&       stictionX, 
              const Array_<Constraint::NoSlip1D>&       stictionZ) 
    :   TriggeredEventHandler(Stage::Position), 
        m_mbs(system), m_contact(contact), m_which(which),
        m_coefRest(coefRest), m_coefFric(coefFric),
        m_stictionX(stictionX), m_stictionZ(stictionZ)
    { 
        // Trigger only as height goes from positive to negative.
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    Real getValue(const State& state) const {
        const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
        const Constraint::PointInPlane& contact = m_contact[m_which];
        if (!contact.isDisabled(state)) 
            return 0; // already locked
        const MobilizedBody& b2 = 
            matter.getMobilizedBody(contact.getFollowerMobilizedBodyIndex());
        const Vec3 pt2 = contact.getDefaultFollowerPoint();
        const Vec3 p_G = b2.findStationLocationInGround(state, pt2);
        const Real height = p_G[YAxis];
        //printf("    Witness@%.16g: %.16g\n", state.getTime(), height);
        return height;
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
                                 Array_<CInfo,int>& proximal) const;



    // Given the set of proximal constraints, prevent penetration by applying
    // a nonnegative least squares impulse generating a step change in 
    // velocity. On return, the applied impulse and new velocities are recorded
    // in proximal, and state is updated to the new velocities and realized
    // through Velocity stage. Constraints that were stopped are enabled, those
    // that rebounded are disabled.
    void processCompressionPhase(Array_<CInfo,int>& proximal,
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
    bool processExpansionPhase(Array_<CInfo,int>& proximal,
                               State&             state) const;

    // Examine the rebounders to see if any are rebounding with a speed at or
    // below the capture velocity. If so, enable those constraints and apply a
    // (hopefully small) negative impulse to eliminate that rebound velocity.
    // Repeat if that induces any negative velocities or any further slow
    // rebounders. This terminates will all rebounders leaving with velocities
    // greater than vCapture, or else all constraints are enabled.
    void captureSlowRebounders(Real               vCapture,
                               Array_<CInfo,int>& proximal,
                               State&             state) const;

    // All proximal, inactive constraints should have significant rebound
    // velocities. Take a small explicit Euler step to advance until all
    // such constraints have perr>0 so we won't miss the next collision.
    void satisfyPositionConditions
       (const Array_<CInfo,int>& proximal, State& state) const;

    // This is the final pass. We realize accelerations, and see if any of the
    // contact forces are negative. If so we disable those constraints.
    // TODO: need to search for a consistent set of active contraints.
    void disablePullingContacts(Array_<CInfo,int>& proximal,
                                State&             state) const;


    // This method is used at the start of compression phase to modify any
    // constraint parameters as necessary, and then enable all the proximal
    // constraints. Some or all of these will be disabled during the impact
    // analysis in compression or expansion phases. On return the state has
    // been updated and realized through Instance stage.
    void enableAllProximalConstraints(Array_<CInfo,int>& proximal,
                                      State&             state) const;

    // Given only the subset of proximal constraints that are active, calculate
    // the impulse that would eliminate all their velocity errors. No change is
    // made to the set of active constraints. Some of the resulting impulses
    // may be negative.
    void calcStoppingImpulse
       (const Array_<CInfo,int>&    proximal,
        const State&                state,
        Vector&                     lambda0) const;

    // Given the initial generalized speeds u0, and a constraint-space impulse
    // lambda, calculate the resulting step velocity change du, modify the
    // generalized speeds in state to u0+du, and realize Velocity stage.
    void updateVelocities(const Vector& u0, const Vector& lambda, 
                          State& state) const;



private:
    const MultibodySystem&                  m_mbs; 
    const Array_<Constraint::PointInPlane>& m_contact;  // no penetration in Y
    const unsigned                          m_which;
    const Array_<Real>&                     m_coefRest; // one per contact
    const Array_<Real>&                     m_coefFric; // mu_d
    const Array_<Constraint::NoSlip1D>&     m_stictionX;
    const Array_<Constraint::NoSlip1D>&     m_stictionZ;
};


//==============================================================================
//                          CONTACT OFF HANDLER
//==============================================================================
class ContactOff: public TriggeredEventHandler {
public:
    ContactOff(const MultibodySystem&           system,
        const Array_<Constraint::PointInPlane>& contact,
        unsigned                                which) 
    :   TriggeredEventHandler(Stage::Acceleration), 
        m_system(system), m_contact(contact), m_which(which)
    { 
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    Real getValue(const State& state) const {
        const Constraint::PointInPlane& contact = m_contact[m_which];
        if (contact.isDisabled(state)) return 0;
        const Real f = -contact.getMultiplier(state); // watch sign
        return f;
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const 
    {
        printf("\n------------------------------------------------------\n");
        printf("LIFTOFF trigged by constraint %d @t=%.15g\n", 
            m_which, s.getTime());
        m_system.realize(s, Stage::Acceleration);
        cout << " triggers=" << s.getEventTriggers() << "\n";
        const Vector& mults = s.getMultipliers();
        Array_<int> toBeDisabled;
        Vector myMults(3);
        for (unsigned i=0; i < m_contact.size(); ++i) {
            const Constraint::PointInPlane& contact = m_contact[i];
            if (contact.isDisabled(s)) continue;
            contact.getMyPartFromConstraintSpaceVector(s,mults,myMults);
            if (-myMults[YAxis]<0) { // watch sign
                printf("  disabling %d because force=%g", i, -myMults[YAxis]);
                toBeDisabled.push_back(i);
            }
        }
        printf("\n");

        for (unsigned tbd=0; tbd < toBeDisabled.size(); ++tbd) {
            const Constraint::PointInPlane& contact = m_contact[toBeDisabled[tbd]];
            contact.disable(s);
        }
        m_system.realize(s, Stage::Instance);

        // Leave at acceleration stage.
        m_system.realize(s, Stage::Acceleration);
        printf("LIFTOFF DONE; %d contacts broken.\n\n", toBeDisabled.size());
        printf("\n------------------------------------------------------\n");
    }

private:
    const MultibodySystem&                  m_system; 
    const Array_<Constraint::PointInPlane>& m_contact;
    const unsigned                          m_which; // one of the contacts
};

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;



static Real g = 9.8;



//==============================================================================
//                                   MAIN
//==============================================================================
int main(int argc, char** argv) {
    static const Transform GroundFrame;
    static const Rotation ZUp(UnitVec3(XAxis), XAxis, UnitVec3(YAxis), ZAxis);
    static const Vec3 TestLoc(1,0,0);

  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem             mbs;

    SimbodyMatterSubsystem      matter(mbs);
    GeneralForceSubsystem       forces(mbs);
    Force::Gravity              gravity(forces, matter, Vec3(0, -g, 0));

    MobilizedBody& Ground = matter.updGround();

    // Predefine some handy rotations.
    const Rotation Z90(Pi/2, ZAxis); // rotate +90 deg about z


        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS
    const Real CoefRest = 0.8;
    const Real CoefFric = 0.9;
    Array_<Constraint::PointInPlane> contacts;
    Array_<Real>                     coefRest; // e (Poisson's restitution)
    Array_<Real>                     coefFric; // mu_d
    Array_<Constraint::NoSlip1D>     stictionX;
    Array_<Constraint::NoSlip1D>     stictionZ;

    const Vec3 CubeHalfDims(3,2,1);
    const Real CubeMass = 1;
    Body::Rigid cubeBody = 
        Body::Rigid(MassProperties(CubeMass, Vec3(0), 
                        UnitInertia::brick(CubeHalfDims)));

    // First body: cube
    MobilizedBody::Free cube(Ground, Vec3(0),
                             cubeBody, Vec3(0));
    cube.addBodyDecoration(Transform(), DecorativeBrick(CubeHalfDims)
                                        .setColor(Red).setOpacity(.3));
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(CubeHalfDims);
        contacts.push_back
            (Constraint::PointInPlane(Ground, YAxis, Zero, cube, pt));
        coefRest.push_back(CoefRest);
        coefFric.push_back(CoefFric);
        stictionX.push_back
            (Constraint::NoSlip1D(Ground, Vec3(0), XAxis, Ground, cube));
        stictionZ.push_back
            (Constraint::NoSlip1D(Ground, Vec3(0), ZAxis, Ground, cube));

        contacts.back().setDisabledByDefault(true);
        stictionX.back().setDisabledByDefault(true);
        stictionZ.back().setDisabledByDefault(true);
    }

#ifdef NOTDEF
    // Second body: weight
    const Vec3 ConnectEdge1(CubeHalfDims[0],0,CubeHalfDims[2]);
    const Vec3 WeightEdge(-CubeHalfDims[0],-CubeHalfDims[1],0);
    MobilizedBody::Pin weight(cube, 
        Transform(Rotation(Pi/2,XAxis), ConnectEdge1),
        cubeBody, Vec3(WeightEdge));
    weight.addBodyDecoration(Transform(), DecorativeBrick(CubeHalfDims)
                                        .setColor(Gray).setOpacity(.6));
    Force::MobilityLinearSpring(forces, weight, 0, 1000, Pi/4);
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        if (i==-1 && j==-1) continue;
        constraints.push_back
            (Constraint::Ball(Ground, Vec3(0), weight, 
             Vec3(i,j,k).elementwiseMultiply(CubeHalfDims)));
        constraints.back().setDisabledByDefault(true);
        coefRest.push_back(0);
    }

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
        constraints.push_back
            (Constraint::Ball(Ground, Vec3(0), weight2, 
             Vec3(i,j,k).elementwiseMultiply(CubeHalfDims)));
        constraints.back().setDisabledByDefault(true);
        coefRest.push_back(0);
    }
#endif

    Visualizer viz(mbs);
    viz.addDecorationGenerator(new ShowContact(contacts));
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

    StateSaver* stateSaver = new StateSaver(mbs,contacts,integ,ReportInterval);
    mbs.addEventReporter(stateSaver);

    for (unsigned i=0; i < contacts.size(); ++i)
        mbs.addEventHandler(new ContactOn(mbs, contacts,i, coefRest,coefFric,
                                          stictionX, stictionZ));

    for (unsigned i=0; i < contacts.size(); ++i)
        mbs.addEventHandler(new ContactOff(mbs,contacts,i));
  
    State s = mbs.realizeTopology(); // returns a reference to the the default state
    mbs.realizeModel(s); // define appropriate states for this System
    mbs.realize(s, Stage::Instance); // instantiate constraints if any


    // Set initial conditions so the -,-,- vertex is in the -y direction.
    const Rotation R_BC(UnitVec3(CubeHalfDims+1e-7*Vec3(1,0,0)), YAxis, Vec3(1,0,0),XAxis);
    cube.setQToFitTransform(s, Transform(~R_BC, Vec3(0,10,0)));
    cube.setUToFitAngularVelocity(s, Vec3(0,1,0));
    cube.setUToFitLinearVelocity(s, Vec3(1,0,0));

    mbs.realize(s, Stage::Velocity);
    viz.report(s);

    cout << "cube X_FM=" << cube.getMobilizerTransform(s) << endl;
    cout << "cube X_GB=" << cube.getBodyTransform(s) << endl;

    mbs.realize(s, Stage::Acceleration);

    cout << "q=" << s.getQ() << endl;
    cout << "u=" << s.getU() << endl;
    cout << "qerr=" << s.getQErr() << endl;
    cout << "uerr=" << s.getUErr() << endl;
    cout << "udoterr=" << s.getUDotErr() << endl;
    cout << "mults=" << s.getMultipliers() << endl;
    cout << "qdot=" << s.getQDot() << endl;
    cout << "udot=" << s.getUDot() << endl;
    cout << "qdotdot=" << s.getQDotDot() << endl;
    viz.report(s);

    cout << "Initial configuration shown. Next? ";
    getchar();

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
    //cout << "On times: " << contactOn->getOnTimes() << endl;
    //cout << "Off times: " << contactOff->getOffTimes() << endl;

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());

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

class CInfo {
public:
    CInfo(const SimbodyMatterSubsystem&   matter,
          const State&                    state,
          const Constraint::PointInPlane& noPenetrationY,
          int                             which,
          Real                            coefRest,
          Real                            coefFric,
          const Constraint::NoSlip1D&     stictionX,
          const Constraint::NoSlip1D&     stictionZ)
    :   contact(&noPenetrationY), which(which),
        stickX(&stictionX), stickZ(&stictionZ),
        body1(&matter.getMobilizedBody(contact->getPlaneMobilizedBodyIndex())), 
        body2(&matter.getMobilizedBody(contact->getFollowerMobilizedBodyIndex())), 
        point2(contact->getDefaultFollowerPoint()),
        initp2_G(body2->findStationLocationInGround(state, point2)),
        initv2_G(body2->findStationVelocityInGround(state, point2)),
        e(coefRest), Ic(0), Ie(0), I(0), V(initv2_G), P(initp2_G)
    {
    }

    Vec3 getMyPartFromConstraintSpaceVector(const State& s, 
                                            const Vector& lambda) const {
        Vector myMults; Vec3 res(0);
        if (isActive(s)) {
            contact->getMyPartFromConstraintSpaceVector(s,lambda,myMults);
            res[YAxis] = myMults[0];
            if (isSticking(s)) {
                stickX->getMyPartFromConstraintSpaceVector(s,lambda,myMults);
                res[XAxis] = myMults[0];
                stickZ->getMyPartFromConstraintSpaceVector(s,lambda,myMults);
                res[ZAxis] = myMults[0];
            }
        }
        return res;
    }

    void setMyPartInConstraintSpaceVector(const State& s, const Vec3& myPart,
                                          Vector& lambda) const {
        Vector myMults(1);
        if (!isActive(s)) return;

        myMults[0] = myPart[YAxis];
        contact->setMyPartInConstraintSpaceVector(s,myMults,lambda);
        if (isSticking(s)) {
            myMults[0] = myPart[XAxis];
            stickX->setMyPartInConstraintSpaceVector(s,myMults,lambda);
            myMults[0] = myPart[ZAxis];
            stickZ->setMyPartInConstraintSpaceVector(s,myMults,lambda);
        }
    }


    void updateFromState(const State& state) {
        body2->findStationLocationAndVelocityInGround(state,point2,P,V);
    }

    void enableAll(State& state) const {
        contact->enable(state); stickX->enable(state); stickZ->enable(state);
    }

    void disableAll(State& state) const {
        contact->disable(state); stickX->disable(state); stickZ->disable(state);
    }

    void enableStiction(State& state) const {
        stickX->enable(state); stickZ->enable(state);
    }

    void disableStiction(State& state) const {
        stickX->disable(state); stickZ->disable(state);
    }
    bool isActive(const State& state) const 
    {   return !contact->isDisabled(state); }
    bool isSticking(const State& state) const 
    {   return !stickX->isDisabled(state); } // both are if one is
        
    Real getHeight()  const {return P[YAxis];}
    Real getDHeight() const {return V[YAxis];}
    const Constraint::PointInPlane* contact;
    const int                       which;
    const Constraint::NoSlip1D*     stickX;
    const Constraint::NoSlip1D*     stickZ;

    const MobilizedBody*    body1;
    const MobilizedBody*    body2;
    const Vec3              point2;   // contact point on body2, in B2 frame
    const Vec3              initp2_G; // location of point2 in ground
    const Vec3              initv2_G; // initial velocity of point2 in G

    // These fields change during impact processing.
    Real e; // effective coefficient of restitution (changes)
    Vec3 Ic; // most recent compression impulse
    Vec3 Ie; // most recent expansion impulse
    Vec3 I; // accumulated impulse at this contact
    Vec3 V; // current velocity (of point2 in G)
    Vec3 P; // current location of point2 in G
};



//------------------------------ HANDLE EVENT ----------------------------------
void ContactOn::
handleEvent(State& s, Real accuracy, bool& shouldTerminate) const 
{
    const Real VCapture=1e-2;

    Array_<CInfo,int> proximal;
    findProximalConstraints(s, accuracy, proximal);

    printf("\nIMPACT for constraint %d at t=%.16g; %d proximal\n", 
        m_which, s.getTime(), proximal.size());

    bool needMoreCompression = true;
    while (needMoreCompression) {
        processCompressionPhase(proximal, s);
        needMoreCompression = false;

        if (processExpansionPhase(proximal, s)) {
            for (int i=0; i<proximal.size(); ++i)
                if (proximal[i].getDHeight() < 0) {
                    needMoreCompression = true;
                    break;
                }
        }
    }

    // Some of the rebounders may be moving so slowly that we would like
    // to be able to say they have stopped. If so, apply additional 
    // (negative) impulses necessary to stop them; enable their contact
    // constraints.
    captureSlowRebounders(VCapture, proximal, s);

    // Advance state until all remaining rebounders have perr>0.
    // TODO: shouldn't need to do this.
    satisfyPositionConditions(proximal, s);

    // Make sure all enabled position and velocity constraints 
    // are satisfied.
    m_mbs.project(s, accuracy);

    // Finally, evaluate accelerations and reaction forces and check if 
    // any of the active contacts are generating negative ("pulling") 
    // forces; if so, inactivate them.
    // TODO: this causes trouble because it is ignoring the forces being
    // produced to maintain stiction, which can also induce downward 
    // acceleration.
    //disablePullingContacts(proximal, s);
    m_mbs.realize(s, Stage::Acceleration);

    printf("END OF IMPACT for %d proximal constraints:\n",proximal.size());
    const Vector& mults = s.getMultipliers();
    for (int i=0; i < proximal.size(); ++i) {
        const CInfo& ci = proximal[i];
        printf("  %d %3s: I=%g, V=%g",
            ci.which, ci.isActive(s) ? "ON" : "off", 
            ci.I[YAxis], ci.V[YAxis]);

        Vec3 myMults = ci.getMyPartFromConstraintSpaceVector(s, mults);
        const Vec3 A = ci.body2->findStationAccelerationInGround(s,ci.point2);
        printf(" F=%g, A=%g\n", -myMults[YAxis], A[YAxis]);          
    }
    printf("DONE WITH IMPACT.\n\n");
}



//------------------------ FIND PROXIMAL CONSTRAINTS ---------------------------
void ContactOn::
findProximalConstraints(const State&       s,
                        Real               posTol,
                        Array_<CInfo,int>& proximal) const
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(s, Stage::Velocity);

    proximal.clear();
    for (unsigned i=0; i<m_contact.size(); ++i) {
        CInfo ci(matter, s, m_contact[i], i, m_coefRest[i], m_coefFric[i],
                 m_stictionX[i], m_stictionZ[i]);
        if (ci.isActive(s) || ci.getHeight() <= posTol) 
        {
            proximal.push_back(ci);
        }

        //SimTK_ASSERT3_ALWAYS(height >= -posTol,
        //    "ContactOn::processImpact(): constraint %d had height %g but tol=%g\n",
        //    i, height, posTol);
    }
}



//------------------------ PROCESS COMPRESSION PHASE ---------------------------
void ContactOn::
processCompressionPhase(Array_<CInfo,int>&  proximal,
                        State&              s) const
{
    printf("Entering processCompressionPhase() ...\n");
    Vector lambda0, lambdaTry;
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
        for (int i=0; i < proximal.size(); ++i) {
            const CInfo& ci = proximal[i];
            if (!ci.isActive(s))
                continue;
            Vec3 myMults = ci.getMyPartFromConstraintSpaceVector(s, lambda0);
            if (myMults[YAxis] > 0)
                continue;
            maybeDisabled.push_back(i);
            printf("  constraint %d is candidate, because lambda_y=%g\n",
                ci.contact->getConstraintIndex(), myMults[YAxis]);
        }
        if (maybeDisabled.empty())
            break;

        // Disable the candidates, then see if they rebound.
        for (unsigned d=0; d < maybeDisabled.size(); ++d)
            proximal[maybeDisabled[d]].disableAll(s);
        m_mbs.realize(s, Stage::Instance);
        calcStoppingImpulse(proximal, s, lambdaTry);
        updateVelocities(u0, lambdaTry, s);

        recapturing.clear();
        for (unsigned i=0; i<maybeDisabled.size(); ++i) {
            const CInfo& ci = proximal[maybeDisabled[i]];
            const Vec3 newv2_G = 
                ci.body2->findStationVelocityInGround(s,ci.point2);
            printf("  candidate constraint %d would have v_y=%g\n",
                ci.contact->getConstraintIndex(), newv2_G[YAxis]);
            if (newv2_G[YAxis] <= 0) {
                recapturing.push_back(maybeDisabled[i]);
                printf("  RECAPTURING constraint %d with v_y=%g\n", 
                    ci.contact->getConstraintIndex(), newv2_G[YAxis]);
            }
        }

        // Re-enable the recaptured candidates.
        if (!recapturing.empty()) {
            for (unsigned c=0; c < recapturing.size(); ++c)
                proximal[recapturing[c]].enableAll(s);
            m_mbs.realize(s, Stage::Instance);
        }

        const int numDisabled = maybeDisabled.size()-recapturing.size();
        if (numDisabled == 0) {
            printf("  None of the candidates was actually disabled.\n");
            // lambda0 is still correct
            break;
        }

        if (recapturing.empty()) lambda0 = lambdaTry;
        else calcStoppingImpulse(proximal, s, lambda0);

        printf("  RETRY with %d constraints disabled\n", numDisabled);
    }
    updateVelocities(u0, lambda0, s);

    // Now update the entries for each proximal constraint to reflect the
    // compression impulse and post-compression velocity.
    printf("  Compression results:\n");
    for (int i=0; i < proximal.size(); ++i) {
        CInfo& ci = proximal[i];
        ci.updateFromState(s); // set current velocity V
        ci.Ic = ci.getMyPartFromConstraintSpaceVector(s, lambda0);
        ci.I += ci.Ic; // accumulate impulse
        printf("  %d %3s: Ic=%g, V=%g\n",
            ci.which, ci.isActive(s) ? "ON" : "off", 
            ci.Ic[YAxis], ci.V[YAxis]);
    }

    cout << "... compression phase done.\n";
}



//------------------------- PROCESS EXPANSION PHASE ----------------------------
bool ContactOn::
processExpansionPhase(Array_<CInfo,int>&  proximal,
                      State&              s) const
{
    printf("Entering processExpansionPhase() ...\n");

    // Generate an expansion impulse if there were any active contacts that
    // still have some restitution remaining.
    Vector expansionImpulse;

    bool anyChange = false;
    for (int i=0; i<proximal.size(); ++i) {
        CInfo& ci = proximal[i];
        if (!ci.isActive(s)) continue;
        Vec3 myLambda = ci.Ic; // compression impulse
        if (myLambda[YAxis] > 0 && ci.e != 0) {
            // NOTE: THIS IS NON-PHYSICAL and can gain energy
            myLambda[XAxis]=myLambda[ZAxis]=0; // no friction rebound
            myLambda[YAxis] *= ci.e;
            // With friction rebound energy is conserved.
            //myLambda *= ci.e;
            ci.setMyPartInConstraintSpaceVector(s, myLambda, 
                                                expansionImpulse);
            ci.e = 0; // we have now used up the material restitution
            anyChange = true;
        }
    }

    if (!anyChange) {
        printf("... no expansion impulse -- done.\n");
        return false;
    }

    // We generated an expansion impulse. Apply it and update velocities.
    updateVelocities(Vector(), expansionImpulse, s);

    // Now update the entries for each proximal constraint to reflect the
    // expansion impulse and post-expansion velocity.
    for (int i=0; i < proximal.size(); ++i) {
        CInfo& ci = proximal[i];
        ci.updateFromState(s); // set current velocity V
        ci.Ie = ci.getMyPartFromConstraintSpaceVector(s, expansionImpulse);
        ci.I += ci.Ie; // accumulate impulse

    }

    // Release any constraint that now has a positive velocity.
    printf("  Expansion results:\n");
    bool anyDisabled = false;
    for (int i=0; i < proximal.size(); ++i) {
        const CInfo& ci = proximal[i];
        if (ci.isActive(s) && ci.getDHeight() > 0) {
            ci.disableAll(s);
            anyDisabled = true;
        }
        printf("  %d %3s: Ie=%g, V=%g\n",
            ci.which, ci.isActive(s) ? "ON" : "off", 
            ci.Ie[YAxis], ci.V[YAxis]);
    }

    if (anyDisabled)
        m_mbs.realize(s, Stage::Velocity);

    cout << "... expansion phase done.\n";

    return true;
}


//------------------------- CAPTURE SLOW REBOUNDERS ----------------------------
void ContactOn::
captureSlowRebounders(Real                vCapture,
                      Array_<CInfo,int>&  proximal,
                      State&              s) const
{
    // Capture any rebounder whose velocity is too slow.
    printf("Entering captureSlowRebounders() ...\n");

    int nCaptured=0, nPasses=0;
    while (true) {
        ++nPasses;
        printf("  start capture pass %d:\n", nPasses);
        bool anyToCapture = false;
        for (int i=0; i < proximal.size(); ++i) {
            const CInfo& ci = proximal[i];
            if (!ci.isActive(s) && ci.getDHeight() <= vCapture) {
                ci.enableAll(s);
                anyToCapture = true;
                ++nCaptured;
                printf("  capturing constraint %d with v=%g\n",
                    ci.which, ci.getDHeight());
            }
        }

        if (!anyToCapture) {
            if (nCaptured==0) printf("... done -- nothing captured.\n");
            else printf("... done -- captured %d rebounders in %d passes.\n",
                nCaptured, nPasses);
            return;
        }

        m_mbs.realize(s, Stage::Velocity);
        Vector captureImpulse;
        calcStoppingImpulse(proximal, s, captureImpulse);
        updateVelocities(Vector(), captureImpulse, s);

        // Now update the entries for each proximal constraint to reflect the
        // capture impulse and post-capture velocity.
        for (int i=0; i < proximal.size(); ++i) {
            CInfo& ci = proximal[i];
            ci.updateFromState(s); // set current velocity V
            ci.I += ci.getMyPartFromConstraintSpaceVector(s, captureImpulse);
        }
    }
}



//---------------------- ENABLE ALL PROXIMAL CONSTRAINTS -----------------------
void ContactOn::
enableAllProximalConstraints(Array_<CInfo,int>&  proximal,
                             State&              state) const
{
    // Set the contact point and enable the constraints.
    for (int i=0; i < proximal.size(); ++i) {
        const CInfo& ci = proximal[i];
        if (ci.isActive(state)) continue;

        ci.stickX->setContactPoint(state, ci.P);
        ci.stickZ->setContactPoint(state, ci.P);
        ci.enableAll(state);
    }
    m_mbs.realize(state, Stage::Instance);
}



//------------------------- CALC COMPRESSION IMPULSE ---------------------------
// Calculate the impulse that eliminates all residual velocity for the
// current set of enabled constraints.
void ContactOn::
calcStoppingImpulse(const Array_<CInfo,int>&    proximal,
                    const State&                s,
                    Vector&                     lambda0) const
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(s, Stage::Dynamics); // TODO: should only need Position
    Vector desiredDeltaV;  // in constraint space
    printf("  Entering calcStoppingImpulse() ...\n");
    bool gotOne = false;
    for (int i=0; i < proximal.size(); ++i) {
        const CInfo& ci = proximal[i];
        if (!ci.isActive(s))
            continue;
        cout << "    constraint " << ci.contact->getConstraintIndex()
             << " enabled, v=" << ci.V << endl;
        ci.setMyPartInConstraintSpaceVector(s, -ci.V, desiredDeltaV);
        gotOne = true;
    }
    if (gotOne) matter.solveForConstraintImpulses(s, desiredDeltaV, lambda0);
    else lambda0.clear();

    cout << "  ... done. Stopping impulse=" << lambda0 << endl;
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



//------------------------ SATISFY POSITION CONDITIONS -------------------------
void ContactOn::
satisfyPositionConditions
    (const Array_<CInfo,int>& proximal, State& state) const
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(state, Stage::Velocity);
    const Vector& qdot = state.getQDot(); // grab before invalidated

    printf("satisfyPositionConditions(): examine %d proximals ...\n",
        proximal.size());

    // Push rebounders into positive territory so we can detect a
    // quick transition back to contact.
    Real maxDt = 0;
    for (int e=0; e < proximal.size(); ++e) {
        const CInfo& ci = proximal[e];
        if (ci.isActive(state))
            continue;

        Vec3 rp2_G, rv2_G;
        ci.body2->findStationLocationAndVelocityInGround
           (state, ci.point2, rp2_G, rv2_G);
        const Real h=rp2_G[YAxis], dh=rv2_G[YAxis];
        printf("  rebounder %d has h=%g, dh=%g\n", 
            ci.contact->getConstraintIndex(), h, dh); 
        if (h <= 0 && dh > 0) {
            const Real dt = -h/dh;
            printf("  -- needs dt=%g.\n", dt); 
            maxDt = std::max(maxDt, dt);
        } else printf("  -- no adjustment needed.\n");
    }

    state.updQ() += 2*maxDt*qdot;
    m_mbs.realize(state, Stage::Position);
    printf("... done with satisfyPositionConditions().\n");
}



//-------------------------- DISABLE PULLING CONTACTS --------------------------
void ContactOn::
disablePullingContacts(Array_<CInfo,int>& proximal,
                       State&             state) const
{
    m_mbs.realize(state, Stage::Acceleration);
    const Vector& lambda = state.getMultipliers();

    Array_<int> pulling;
    for (int i=0; i < proximal.size(); ++i) {
        const CInfo& ci = proximal[i];
        if (!ci.isActive(state)) continue;
        const Vec3 myLambda = ci.getMyPartFromConstraintSpaceVector(state, lambda);
        if (-myLambda[YAxis] < 0) { // multipliers are negative forces
            if (pulling.empty())
                printf("disablePullingContacts(): disabling");
            printf(" %d cuz f=%g %g %g", ci.which, 
                -myLambda[XAxis],-myLambda[YAxis],-myLambda[ZAxis]);
            pulling.push_back(i);
        }
    }

    if (pulling.empty()) {
        printf("disablePullingContacts(): nobody is pulling.\n");
        return;
    }

    for (unsigned p=0; p < pulling.size(); ++p) {
        const CInfo& ci = proximal[pulling[p]];
        ci.contact->disable(state);
    }
    printf(".\n");
    m_mbs.realize(state, Stage::Instance);

    // Always leave at acceleration stage.
    m_mbs.realize(state, Stage::Acceleration);
}
