/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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
 * Adhoc main program for playing with wrist-like moment arm problem.
 **/

/*
 *               q0     q1     q2                    x
 * |-----------|                 |-------|           ^
 * |  Forearm  | * aux1 * aux2 * | Hand  |           |
 * |   =====   |           |     | ===== |    y <----* z
 * |-----------|<..........O....>|-------|
 *               muscle & tendon
 *                  length l(q)
 *
 * aux1,2 are small bones and "*" is a rotational dof "|" is a sliding dof.
 * The "O" is a circular wrap surface attached to a body that is on a slider
 * connected to aux2, where the sliding dof is coupled to one of the 
 * rotational dofs and affects the length of the muscle as it passes
 * frictionlessly over the wrap circle. There are coupler constraints 
 * connecting all the dofs so that setting any q determines all the others
 * and thus the muscle length; alternatively we can turn off one of the
 * coupler constraints and use a contact constraint to remove a dof from
 * the hand, coupling all the q's indirectly.
 *
 * We would like a good 
 * instantaneous calculation of the moment arm r(q) defined as r=dl/dtheta 
 * where l(q) is muscle length and theta is the angle between a line fixed in 
 * the hand (=====) and a line fixed in the forearm. Theta is also the sum
 * of the three rotational joint angles; we're assuming the q's are 
 * scaled angles such that angle a0=q0/s0, a1=q1/s1, and a2=q2/s2 so we
 * have theta=q0/s0+q1/s1+q2/s2. In OpenSim it is common to scale the
 * generalized coordinate of the "independent" angle so that its numerical
 * value is the total angle theta.
 *
 * We can calculate r=ldot/thetadot and if all
 * velocities are zero we can also calculate r=ldotdot/thetadotdot. 
 * Unfortunately, in practice it is difficult to determine ldot and ldotdot (unless
 * the muscle is a straight line; i.e., no wrap surface), otherwise
 * this would be an easy solution: apply a muscle tension t along the
 * muscle path, measure the resulting ldotdot and thetadotdot and calculate
 * r directly. But we can get the same result with the following method:
 *
 * 1. Determine the coupling matrix C such that qdot=C * thetadot.
 * 2. Calculate f(t) the generalized forces due to an applied muscle tension t
 * 3. Calculate r=~C*f/t.
 *
 * See the document "Moment arm revisited", Sherm, 10/12/2010 for theory.
 *
 * Below we will use this method and compare it with a direct calculation from
 * the definition r=dl/dtheta done by perturbation.
 */

#include "SimTKsimbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <ctime>
using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real ReportInterval = 0.5;

// Utility for analytically calculating the two tangent points of a
// line from point p to a circle of radius r, center c. Note that this
// is limited to a coplanar point and circle.
static bool find_tangent_points(const Vec2& c, Real r, const Vec2& p,
                                Vec2& p0, Vec2& p1);
// Given chord length h in a circle of radius r, return the corresponding
// arc length, going from one end of the chord the short way around to 
// the other end.
static Real calc_arc_length(Real r, Real h);



//==============================================================================
// DEFINE THE MUSCLE
//==============================================================================

// This helper class holds cached path geometry calculations done by the 
// muscle below.
class PathInfo {
public:
    PathInfo() 
    :   p1(NaN), p2(NaN), pathLengthA(NaN), arcLength(NaN), pathLengthB(NaN) {}
    Vec2 p1, p2; // tangent points on body C, in G
    Real pathLengthA, arcLength, pathLengthB;
};

std::ostream& operator<<(std::ostream& o, const PathInfo& pi) {
    o << "PathInfo: p1=" << pi.p1 << " p2=" << pi.p2 << std::endl;
    return o;
}

// This is a model of a muscle with origin/insertion points on bodies A and B, 
// using a planar path that passes frictionlessly over one side of a wrapping 
// circle of radius r on body W. 
//
// The muscle defines a scalar discrete state variable to hold the current
// tension value t that is used to generate the muscle forces F=T(q)*t.
// T is a force transmission matrix that we don't know explicitly.
class PlanarMuscle : public Force::Custom::Implementation {
public:
    PlanarMuscle(const GeneralForceSubsystem& forces,
                 const SimbodyMatterSubsystem& matter,
                 const MobilizedBody& A, const Vec3& ptA,
                 const MobilizedBody& B, const Vec3& ptB,
                 const MobilizedBody& W, const Vec3& center,
                 Real radius, const UnitVec3& side)
    :   m_forces(forces), m_matter(matter), 
        m_A(A), m_ptA(ptA), m_B(B), m_ptB(ptB),     // origin/insertion points 
        m_W(W), m_center(center), m_radius(radius), // wrapping circle 
        m_side(side) {} // which side of the circle to wrap around (in W frame) 


    //      Specialized interface for PlanarMuscle.

    // Set the muscle tension value.
    void setTension(State& state, Real tension) const
    {   Value<Real>::updDowncast
           (m_forces.updDiscreteVariable(state,m_tensionVarIx)) = tension; }
    // Get the current muscle tension setting.
    Real getTension(const State& state) const
    {   return Value<Real>::downcast
           (m_forces.getDiscreteVariable(state,m_tensionVarIx)); }

    // Cacluate where the two insertion points are in the Ground frame.
    void findInsertionPoints(const State& state, Vec3& pA, Vec3& pB) const {
        const MobilizedBody& A = m_matter.getMobilizedBody(m_A);
        const MobilizedBody& B = m_matter.getMobilizedBody(m_B);

        pA = A.findStationLocationInGround(state, m_ptA);
        pB = B.findStationLocationInGround(state, m_ptB);
    }
    
    // These retrieve items from the PathInfo cache entry that is set
    // in realizePosition() below, so these can't be called until Position
    // stage has been realized in the supplied State.

    Real getPathLength(const State& state) const {
        const PathInfo& info = getPathInfo(state);
        return info.pathLengthA + info.arcLength + info.pathLengthB;
    }
    Real getArcLength(const State& state) const {
        const PathInfo& info = getPathInfo(state);
        return info.arcLength;
    }
    void getTangentPoints(const State& state, Vec3& pAC, Vec3& pBC) const {
        const PathInfo& info = getPathInfo(state);
        pAC = info.p1.append1(0);
        pBC = info.p2.append1(0);
    }
        
    //      Satisfy the virtual methods of Force::Custom::Implementation.
                 
    // Calculate the muscle forces and accumulate into bodyForces array.
    // (We aren't going to generate any particle or mobility forces.)
    virtual void calcForce(const State& state, 
                           Vector_<SpatialVec>& bodyForces, 
                           Vector_<Vec3>& particleForces, 
                           Vector& mobilityForces) const
    {
        const MobilizedBody& A = m_matter.getMobilizedBody(m_A);
        const MobilizedBody& B = m_matter.getMobilizedBody(m_B);
        const MobilizedBody& W = m_matter.getMobilizedBody(m_W);

        // Insertion points in Ground.
        const Vec3 iptA = A.findStationLocationInGround(state, m_ptA);
        const Vec3 iptB = B.findStationLocationInGround(state, m_ptB);
        // Wrap circle center in Ground.
        const Vec3 ctr = W.findStationLocationInGround(state, m_center);

        // Tangent points in Ground.
        Vec3 tptA, tptB;
        getTangentPoints(state, tptA, tptB);
        const Vec3 tptA_W = W.findStationAtGroundPoint(state, tptA);
        const Vec3 tptB_W = W.findStationAtGroundPoint(state, tptB);

        UnitVec3 a2w(tptA - iptA); // from insertion to tangent point
        UnitVec3 b2w(tptB - iptB); //   "
        UnitVec3 nta( ctr - tptA); // inward normal at A's tangent pt
        UnitVec3 ntb( ctr - tptB); //         "        B's   "

        const Real tension = getTension(state);

        const Vec3 fA = tension*a2w;
        const Vec3 fB = tension*b2w;

        A.applyForceToBodyPoint(state, m_ptA, fA, bodyForces);
        B.applyForceToBodyPoint(state, m_ptB, fB, bodyForces);
        W.applyForceToBodyPoint(state, m_center, -(fA+fB), bodyForces); 
    }

    // This muscle model doesn't store energy.
    virtual Real calcPotentialEnergy(const State&) const {return 0;}

    // Allocate variables and cache entries.
    virtual void realizeTopology(State& state) const {
        PlanarMuscle* mThis = const_cast<PlanarMuscle*>(this);
        mThis->m_tensionVarIx = m_forces
            .allocateDiscreteVariable(state, Stage::Dynamics, 
                                      new Value<Real>(0));

        mThis->m_pathInfoIx = m_forces
            .allocateCacheEntry(state, Stage::Position, 
                                new Value<PathInfo>());
    }

    // Calculate path geometry information.
    virtual void realizePosition(const State& state) const {
        // Get all the relevant points in Ground.
        const MobilizedBody& A = m_matter.getMobilizedBody(m_A);
        const MobilizedBody& B = m_matter.getMobilizedBody(m_B);
        const MobilizedBody& W = m_matter.getMobilizedBody(m_W);

        // Only works in x-y plane (we're dropping the z coordinate).
        const Vec2 ptA = A.findStationLocationInGround(state, m_ptA).drop1(2);
        const Vec2 ptB = B.findStationLocationInGround(state, m_ptB).drop1(2);
        const Vec2 ctr = W.findStationLocationInGround(state, m_center).drop1(2);
        const Vec2 side = W.expressVectorInGroundFrame(state, m_side).drop1(2);
        Vec2 p0,p1,a,b;
        find_tangent_points(ctr,m_radius,ptA,p0,p1);
        a = dot(p0-ctr,side) > dot(p1-ctr,side) ? p0 : p1;
        find_tangent_points(ctr,m_radius,ptB,p0,p1);
        b = dot(p0-ctr,side) > dot(p1-ctr,side) ? p0 : p1;

        PathInfo& info = updPathInfo(state);
        info.p1 = a; info.p2 = b; 
        info.pathLengthA = (a-ptA).norm();
        info.pathLengthB = (b-ptB).norm();
        info.arcLength = calc_arc_length(m_radius, (a-b).norm());
    }

private:
    const PathInfo& getPathInfo(const State& state) const
    {   return Value<PathInfo>::downcast
            (m_forces.getCacheEntry(state,m_pathInfoIx)); }
    PathInfo& updPathInfo(const State& state) const
    {   return Value<PathInfo>::updDowncast
            (m_forces.updCacheEntry(state,m_pathInfoIx)); }

    const GeneralForceSubsystem&    m_forces;
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_A, m_B, m_W;
    const Vec3                      m_ptA, m_ptB, m_center;
    const Real                      m_radius;
    const UnitVec3                  m_side;
    DiscreteVariableIndex           m_tensionVarIx;
    CacheEntryIndex                 m_pathInfoIx;
};



//==============================================================================
// DEFINE A CUSTOM MOBILIZER
//==============================================================================
// This is a pin joint but with a generalized coordinate 
//      q0=scale*theta0
// where theta0 is its actual angle. We'll use a coupler to make the 
// other pin joints have angles q1=theta0, q2=theta0 so that the 
// total angle theta0+theta1+theta2=q0.
class ScaledPin : public MobilizedBody::FunctionBased {
public:
    ScaledPin(MobilizedBody& parent, const Transform& inbFrameF, 
            const Body& body, const Transform& outbFrameM,
            Real scale)
    :   FunctionBased(parent,inbFrameF,body,outbFrameM,1,
            makeFunctions(scale), Array_<Array_<int> >(6, Array_<int>(1,0))),
        m_scale(scale)
    {}

    Real getScale() const {return m_scale;}

    void setQ(State& s, Real q) const {setOneQ(s,MobilizerQIndex(0),q);}
    void setU(State& s, Real u) const {setOneU(s,MobilizerUIndex(0),u);}
    Real getQ(const State& s) const {return getOneQ(s,MobilizerQIndex(0));}
    Real getU(const State& s) const {return getOneU(s,MobilizerQIndex(0));}

    void setAngle(State& s, Real angle) const 
    {   setOneQ(s, MobilizerQIndex(0), m_scale*angle); }
    Real getAngle(const State& s) const 
    {   return getOneQ(s, MobilizerQIndex(0)) / m_scale; }

    void setAngularRate(State& s, Real angularRate) const 
    {   setOneU(s, MobilizerQIndex(0), m_scale*angularRate); }
    Real getAngularRate(const State& s) const 
    {   return getOneU(s, MobilizerQIndex(0)) / m_scale; }

private:
    Array_<const Function*> makeFunctions(Real scale) {
        Array_<const Function*> funcs;
        funcs.push_back(new Function::Constant(0));  // x rotation
        funcs.push_back(new Function::Constant(0));  // y rotation
        funcs.push_back(new Function::Linear
                    (Vector(Vec2(1/scale,0)))); // z rotation
        funcs.push_back(new Function::Constant(0));  // x translation
        funcs.push_back(new Function::Constant(0));  // y translation
        funcs.push_back(new Function::Constant(0));  // z translation
        return funcs;
    }

    const Real m_scale;
};

//==============================================================================
// DEFINE REPORTER
//==============================================================================
// This reporter both prints to the console and drives the visualizer. We have
// to add the path lines to each frame here since the tangent end points move.
class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system, 
               Visualizer& viz,
               const PlanarMuscle& planarMuscle,
               Real reportInterval)
    :   PeriodicEventReporter(reportInterval), m_system(system), m_viz(viz), 
        m_planarMuscle(planarMuscle) {}

    ~MyReporter() {}

    // This is used by handleEvent() below but can also be called directly
    // to generate a frame from saved states for replay.
    void report(const State& state) const {
        m_system.realize(state, Stage::Position);
        m_viz.report(state);
        cout << "t=" << state.getTime() 
             << " path length=" << m_planarMuscle.getPathLength(state) << endl;
    }

    // Supply the required virtual method.
    virtual void handleEvent(const State& state) const {
        report(state);
        saveEm.push_back(state);
    }
private:
    const MultibodySystem&           m_system;
    Visualizer&                      m_viz;
    const PlanarMuscle&              m_planarMuscle;
};



//==============================================================================
// DRAW PATH LINES
//==============================================================================
// We have to add the path lines to each frame here since the tangent end points 
// move.
class DrawPathLines : public DecorationGenerator {
public:
    DrawPathLines(const MultibodySystem& system,
                  const PlanarMuscle&    muscle) 
    :   m_system(system), m_planarMuscle(muscle) {}

    virtual void generateDecorations(const State& state, 
                                     Array_<DecorativeGeometry>& geometry) 
    {
        m_system.realize(state, Stage::Position);
        Vec3 iptA, iptB, tptA, tptB;
        m_planarMuscle.findInsertionPoints(state, iptA, iptB);
        m_planarMuscle.getTangentPoints(state, tptA, tptB);
        geometry.push_back(DecorativeLine(iptA, tptA).setColor(Red));
        geometry.push_back(DecorativeLine(iptB, tptB).setColor(Red));
    }
private:
    const MultibodySystem&      m_system;
    const PlanarMuscle&         m_planarMuscle;
};



//==============================================================================
//                                 MAIN
//==============================================================================
int main() {

  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    //Force::UniformGravity   gravity(forces, matter, 1*Vec3(2, -9.8, 0));

    const Real forearmMass = 11;
    const Real handMass  = 3;
    const Real auxMass  = .001;
    const Real wrapMass = .1;
    const Vec3 forearmHDim(1,3,.25);
    const Vec3 handHDim(1,1,.25);
    const Vec3 auxHDim(.1, .1, .05);
    const Real wrapRadius = .25;  // circle in x,y
    const Real wrapHThick = .1;   // z half-depth

    const Vec3 c(2,3,5);          // constrained angle ratios q0:q1:q2
    const Real csum = sum(c);
    const Vec4 scale(csum/c[0],csum/c[1],csum/c[2],1); // scaling
    cout << "scaling=" << scale << endl;

    // Stiffnesses
    const Real k = 1000;

    // Input tension
    const Real tension = 2;

    // Global damper.
    const Real damping = 100;

    // Muscle attachment points
    const Vec3 forearmAttach(forearmHDim[0]/2,-forearmHDim[1],0);
    const Vec3 handAttach(handHDim[0], handHDim[1], 0);

    const DecorativeGeometry forearmViz = DecorativeBrick(forearmHDim).setColor(Green);
    const DecorativeGeometry handViz  = DecorativeBrick(handHDim).setColor(Orange);
    const DecorativeGeometry auxViz  = DecorativeBrick(auxHDim).setColor(Blue);
    const DecorativeGeometry wrapViz  = DecorativeCylinder(wrapRadius,wrapHThick)
        .setTransform(Rotation(Pi/2, XAxis)).setColor(Purple).setResolution(10);

    Body::Rigid forearmBody(MassProperties(forearmMass, Vec3(0), 
                                         forearmMass*UnitInertia::brick(forearmHDim)));
    forearmBody.addDecoration(Vec3(0), forearmViz);

    Body::Rigid handBody(MassProperties(handMass, Vec3(0), 
                                         handMass*UnitInertia::brick(handHDim)));
    handBody.addDecoration(Vec3(0), handViz);

    Body::Rigid auxBody(MassProperties(auxMass, Vec3(0), 
                                         auxMass*UnitInertia::brick(auxHDim)));
    auxBody.addDecoration(Vec3(0), auxViz);

    Body::Rigid wrapBody(MassProperties(wrapMass, Vec3(0),
        wrapMass*UnitInertia::cylinderAlongZ(wrapRadius, wrapHThick)));
    wrapBody.addDecoration(Vec3(0), wrapViz);

    MobilizedBody::Weld forearm(matter.Ground(), Vec3(forearmHDim[0], forearmHDim[1], 0), 
                                forearmBody, Vec3(0));


    ScaledPin aux1(forearm, Vec3(-forearmHDim[0], -forearmHDim[1]-.5, 0), 
                   auxBody, Vec3(0, auxHDim[1], 0), scale[0]);

    ScaledPin aux2(aux1, Vec3(0, -auxHDim[1]-.1, 0),
                   auxBody, Vec3(0, auxHDim[1]+.3, 0), scale[1]);

    ScaledPin hand(aux2, Vec3(0, -auxHDim[1], 0), 
                   handBody, Vec3(-handHDim[0], handHDim[1]+.5, 0), scale[2]);

    MobilizedBody::Slider wrap(aux2, Vec3(1.5,0,0),
                               wrapBody, Vec3(0));

    // Spring
    Force::MobilityLinearSpring aux1Spring(forces, aux1, 
                                    MobilizerUIndex(0), k, 0);

    // Straight-line muscle with no wrapping.
    //Force::TwoPointConstantForce muscle(forces, forearm, forearmAttach,
    //                                    hand, handAttach, -tension);

    PlanarMuscle& planarMuscle = *new PlanarMuscle(forces, matter,
                                    forearm, forearmAttach,
                                    hand, handAttach,
                                    wrap, Vec3(0), wrapRadius, UnitVec3(-1,0,0));
    Force::Custom muscle(forces, &planarMuscle);

    // Energy sucker
    Force::GlobalDamper damper(forces, matter, damping);
    damper.setDisabledByDefault(true);

    // Angles are ai=q[i]/s[i]. Then constraints are
    //    a0/c0 - a1/c1 = 0  =>  q[0]/(s0*c0) - q[1]/(s1*c1) = 0
    //    a0/c0 - a2/c2 = 0  =>  q[0]/(s0*c0) - q[2]/(s2*c2) = 0
    const Vec2 ratios1(1/(scale[0]*c[0]), -1/(scale[1]*c[1]));
    const Vec2 ratios2(1/(scale[0]*c[0]), -1/(scale[2]*c[2]));
    Array_<MobilizedBodyIndex> aux1aux2, aux1Hand, aux1Wrap; 
    aux1aux2.push_back(MobilizedBodyIndex(aux1)); 
    aux1aux2.push_back(MobilizedBodyIndex(aux2));
    aux1Hand.push_back(MobilizedBodyIndex(aux1)); 
    aux1Hand.push_back(MobilizedBodyIndex(hand)); 
    aux1Wrap.push_back(MobilizedBodyIndex(aux1)); 
    aux1Wrap.push_back(MobilizedBodyIndex(wrap)); 
    Array_<MobilizerUIndex> whichUs(2, MobilizerUIndex(0));

    Constraint::CoordinateCoupler 
    aux1toAux2(matter,
        new Function::Linear(Vector(Vec3(ratios1[0],ratios1[1],0))),
        aux1aux2, whichUs);
    aux1toAux2.setDisabledByDefault(true);

    Constraint::CoordinateCoupler 
    aux1toHand(matter,
        new Function::Linear(Vector(Vec3(ratios2[0],ratios2[1],0))),
        aux1Hand, whichUs);
    aux1toHand.setDisabledByDefault(true);

    // Coupling ratio to the wrap surface's body is hardcoded to
    // an arbitrary value here tied to q[0] rather than the angle.
    Constraint::CoordinateCoupler
    aux1toWrap(matter,
        new Function::Linear(Vector(Vec3(1, -.4, 0))),
        aux1Wrap, whichUs);

    // Enable this to make a point of the hand touch a plane fixed on ground.
    // We'll optionally use this constraint instead of one of the couplers
    // to verify that moment arm works either way.
    Constraint::PointInPlane inPlane(matter.Ground(), UnitVec3(-1,0,0), -2.25,
                                     hand, Vec3(1,-1,0));
    inPlane.setDisabledByDefault(true);

    system.realizeTopology();

    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor);

    // This draws the straight-line muscle (regardless of whether we're
    // using it).
    viz.addRubberBandLine(forearm, forearmAttach, hand, handAttach,
        DecorativeLine().setColor(Red).setLineThickness(2));

    // Add generator for muscle path lines.
    viz.addDecorationGenerator(new DrawPathLines(system,planarMuscle));
    
    MyReporter& myRep = *new MyReporter(system,viz,planarMuscle,ReportInterval);
    system.addEventReporter(&myRep);

    // Initialize the system and state.
    
    State state = system.getDefaultState();
    aux1.setAngle(state,0);
    aux2.setAngle(state,0);
    hand.setAngle(state,0);
    planarMuscle.setTension(state, 0);

    cout << "Muscle tension = " << planarMuscle.getTension(state) << endl;

    myRep.report(state);
    printf("Initial state -- hit ENTER\n");
    char ch=getchar();

    //ExplicitEulerIntegrator integ(system);
    CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
    //RungeKuttaFeldbergIntegrator integ(system);
    //RungeKuttaMersonIntegrator integ(system);
    //RungeKutta3Integrator integ(system);
    //VerletIntegrator integ(system);
    //integ.setMaximumStepSize(1e-0001);
    integ.setAccuracy(1e-6);
    TimeStepper ts(system, integ);

    damper.enable(state);

#define USE_CONTACT_CONSTRAINT
#ifdef USE_CONTACT_CONSTRAINT
    inPlane.enable(state);
#else
    aux1toAux2.enable(state);
#endif
    aux1toHand.enable(state);
    system.realize(state, Stage::Instance);

    state.updU() = 0;
    state.setTime(0);
    ts.initialize(state);
    ts.stepTo(40.0);

    state = ts.getState();
    system.realize(state, Stage::Velocity);
    const Real pathLength0 = planarMuscle.getPathLength(state);
    const Real length0 = forearm.calcStationToStationDistance(state, forearmAttach, hand, handAttach);
    const Vector q0 = state.getQ(); // equilibrium q's

    myRep.report(state);
    printf("Muscle off ke=%g -- hit ENTER\n", system.calcKineticEnergy(state));
    cout << "pathLength0=" << pathLength0 << " length0=" << length0 << " q0=" << q0 << endl;
    ch=getchar();

    //----------------------------
    // CALCULATE COUPLING MATRIX C
    //----------------------------

    state.updU() = 0;
    aux1.setU(state, 1);
    cout << "before project u=" << state.getU() << endl;
    const Vector yWeights(state.getNY(), 1);
    const Vector cWeights(state.getNMultipliers(), 1);
    Vector yErrEst;
    system.project(state, 1e-10, yWeights, cWeights, yErrEst, 
        System::ProjectOptions::VelocityOnly);
    cout << "after project u=" << state.getU() << " uerr=" << state.getUErr() << endl;
    // To compute thetadot, convert u's to angular rates and add.
    const Real thetaDot =   aux1.getAngularRate(state) 
                          + aux2.getAngularRate(state) 
                          + hand.getAngularRate(state);
    cout << "calculated thetadot=" << thetaDot << endl;
    const Vector calcc(state.getU()/thetaDot); 

    const Real calccsum = sum(calcc);
    cout << "calc c=" << calcc << " csum=" << calccsum << endl;
    state.updU() = 0;

    // None of the acceleratin- and multiplier-dependent stuff here 
    // matters for moment arm; this is just for playing around with 
    // related dynamic quantities. Feel free to ignore.

    system.realize(state, Stage::Acceleration);
    // Calculate the joint torques f0 equivalent to the equilibrium forces.
    Vector f0;
    matter.multiplyBySystemJacobianTranspose(state,
        system.getRigidBodyForces(state,Stage::Dynamics),
        f0);
    f0 += system.getMobilityForces(state,Stage::Dynamics);

    // This is how you find the scaling if you need it; we don't require
    // this knowledge to calculate moment arm since it is implicitly 
    // embedded in the coupling matrix C.
    const SpatialVec H_PB_G = aux1.getHCol(state, MobilizerUIndex(0));
    const SpatialVec H_FM   = aux1.getH_FMCol(state, MobilizerUIndex(0));
    cout << "H_PB_G=" << H_PB_G << " -> angle scale " << H_PB_G[0].norm() << endl;
    cout << "H_FM=" << H_FM << " -> angle scale " << H_FM[0].norm() << endl;

    const Vector udot0 = state.getUDot();
    const Vector lambda0 = state.getMultipliers();
    Matrix Gt;
    matter.calcGt(state, Gt);
    cout << "Gt*l0=" << Gt*lambda0 << endl;
    Vector f0l = f0 - Gt*lambda0;

    cout << "udot0=" << udot0 << endl;
    cout << "lambda0=" << lambda0 << endl;
    cout << "f0=" << f0 << endl;
    cout << "f0l=" << f0l << endl;

    const Real d2l0 = -forearm.calcStationToStationDistance2ndTimeDerivative
                        (state, forearmAttach, hand, handAttach);

    // Fire the muscle.
    planarMuscle.setTension(state, tension);
    cout << "Muscle tension = " << planarMuscle.getTension(state) << endl;

    // See note above regarding irrelevance of acceleration-related 
    // computations for moment arm.
    system.realize(state, Stage::Acceleration);

    const Vector udot1 = state.getUDot();
    const Vector lambda1 = state.getMultipliers();
    cout << "udot1=" << udot1 << endl;
    cout << "lambda1=" << lambda1 << endl;
    const Vector dudot = udot1 - udot0;
    cout << "dudot=" << dudot << endl;
    const Real d2l = -forearm.calcStationToStationDistance2ndTimeDerivative
                        (state, forearmAttach, hand, handAttach);
    const Real dd2l = d2l - d2l0;
    cout << "dd2l=" << dd2l << endl;
    Vector momentArm = dd2l * dudot.elementwiseInvert();
    cout << "dd2l ./ dudot=" << momentArm << endl;
    Vector dudotAdj(dudot); 
    dudotAdj.elementwiseDivideInPlace(Vector(scale)); // fix units so we have angles
    cout << "dudotAdj=" << dudotAdj << endl;
    cout << "dd2l / sum(dudotAdj) = " << dd2l / sum(dudotAdj) << endl;

    // Calculate the joint torques f equivalent to the muscle forces F.
    Vector f1;
    matter.multiplyBySystemJacobianTranspose(state,
        system.getRigidBodyForces(state,Stage::Dynamics),
        f1);
    f1 += system.getMobilityForces(state,Stage::Dynamics);

    matter.calcGt(state, Gt);
    cout << "Gt*l1=" << Gt*lambda1 << endl;
    Vector f1l = f1 - Gt*lambda1;

    Vector f = f1-f0, fl=f1l-f0l;
    cout << "f/tension=" << f/tension << " ~C*f/tension=" 
        << (~calcc*f)/tension << endl;
    cout << "fl/tension=" << fl/tension << " ~C*fl/tension=" 
        << (~calcc*fl)/tension << endl;

    Assembler asmb(system);
    asmb.setAccuracy(1e-10);
    Real tol = asmb.assemble(state);
    system.realize(state, Stage::Position);
    Real startAngle = aux1.getAngle(state)
                        + aux2.getAngle(state)
                        + hand.getAngle(state);
    const Real startLength = planarMuscle.getPathLength(state); 
    cout << "ASSEMBLED to tol=" << tol 
         << " startAngle=" << startAngle
         << " startLength=" << startLength << endl;

    aux1.setQ(state, aux1.getQ(state) + 1e-6); // perturb
    asmb.lockMobilizer(aux1);
    tol = asmb.assemble(state);
    system.realize(state, Stage::Position);
    Real endAngle = aux1.getAngle(state)
                        + aux2.getAngle(state)
                        + hand.getAngle(state);
    const Real endLength = planarMuscle.getPathLength(state); 
    cout << "ASSEMBLED to tol=" << tol 
         << " endAngle=" << endAngle
         << " endLength=" << endLength << endl;
    cout << "r = " << (endLength-startLength)/(endAngle-startAngle) << endl;

    // Here's another way to do it -- starting at equilibrium, apply a
    // tension to the muscle and then find the new equilibrium (actually
    // any change in configuration will work). Then calculated dl and dtheta.


    // Simulate it.
    saveEm.clear(); // Forget the saved trajectory from the equilibrium run above
    integ.setAccuracy(1e-6);
    muscle.enable(state);
    state.updQ() = q0; // back to equilibrium
    state.updU() = 0;
    state.setTime(0);
    ts.initialize(state);
    cout << "Initialized velocity=" << ts.getState().getU() << endl;
    const Real angle0 = hand.getBodyRotation(ts.getState()).convertOneAxisRotationToOneAngle( ZAxis );
    ts.stepTo(40.0);

    state = ts.getState();
    system.realize(state, Stage::Velocity);
    const Real length1 = forearm.calcStationToStationDistance(state, forearmAttach, hand, handAttach);
    const Real pathLength1 = planarMuscle.getPathLength(state); 
    
    const Vector q1 = state.getQ();

    myRep.report(state);
    printf("Tension %g; final equilibrium ke=%g -- hit ENTER\n", 
        tension, system.calcKineticEnergy(state));
    cout << "length1=" << length1 << " q1=" << q1 << endl;
    cout << "pathLength1=" << pathLength1 << endl;
    const Real dl=-(length1-length0); 
    const Real dpl=-(pathLength1-pathLength0);
    const Vector dtheta=q1-q0;
    cout << "dl=" << dl << " dpl=" << dpl << " dtheta=" << dtheta << endl;
    momentArm = dl * dtheta.elementwiseInvert();
    cout << "dl ./ dtheta=" << momentArm << endl;
    cout << "dpl ./ dtheta=" << dpl * dtheta.elementwiseInvert() << endl;
    Vector dthetaAdj(dtheta); 
    dthetaAdj.elementwiseDivideInPlace(Vector(scale)); // convert to angles
    cout << "sum of 3 joint angles=" << sum(dthetaAdj(0,3)) << " angles=" << dthetaAdj << endl;
    const Real angle1 = hand.getBodyRotation(state).convertOneAxisRotationToOneAngle( ZAxis );
    cout << "angle0=" << angle0 << " angle1=" << angle1 << endl;
    cout << "calculated dangle=" << angle1-angle0 << endl;
    cout << " dl/dangle=" << dl / (angle1-angle0) << endl;
    cout << " dpl/dangle=" << dpl / (angle1-angle0) << endl;
    system.realize(state, Stage::Acceleration);
    cout << "perr=" << state.getQErr() << endl;
    cout << "verr=" << state.getUErr() << endl;
    cout << "aerr=" << state.getUDotErr() << endl;
    ch=getchar();

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            myRep.report(saveEm[i]);
        }
        getchar();
    }

  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}

//==============================================================================
//                          PATH GEOMETRY UTILITIES
//==============================================================================
// circle_circle_intersection()
// Determine the two points where 2 circles in a common plane intersect. Returns
// false quickly if they don't intersect at all or if they are exactly coincident.
//
// Method due to Paul Bourke, 1997.
// http://local.wasp.uwa.edu.au/~pbourke/geometry/2circle/
//
// Code adapted from Tim Voght's public domain code linked to be the above
// article, dated 3/26/2005.
//
static bool circle_circle_intersection( const Vec2& c0, Real r0, 
                                        const Vec2& c1, Real r1, 
                                        Vec2& p0, Vec2& p1)
{ 
    const Vec2 c2c = c1 - c0;
    const Real d2 = c2c.normSqr();

    // Check for solvability.
    if (d2 > square(r0 + r1)) { 
      // no solution. circles do not intersect.
      p0 = p1 = NaN;
      return false; 
    } 
    if (d2 < square(r0 - r1)) {
      // no solution. one circle is contained in the other
      p0 = p1 = NaN;
      return false;
    } 
    if (d2 == 0) {
      // circles must be coincident and intersect at every point;
      // treat as unsolvable
      p0 = p1 = NaN;
      return false;
    }

    const Real d = std::sqrt(d2); // distance between the centers
    const Real ood = 1/d; // compute once

    // 'point 2' is the point where the line through the circle
    // intersection points crosses the line between the circle
    // centers.

    // Determine the distance from point 0 to point 2.
    const Real a = (r0*r0 - r1*r1 + d2) * ood / 2;

    // Determine the coordinates of point 2.
    const Vec2 p2 = c0 + c2c*a*ood;

    // Determine the distance from point 2 to the intersection points. 
    const Real h = std::sqrt(r0*r0 - a*a);
    // Now determine the offsets of the intersection points from
    // point 2.
    const Vec2 offs = Vec2(-c2c[1],c2c[0]) * h * ood;

    // Determine the absolute intersection points.
    p0 = p2 + offs;
    p1 = p2 - offs;
    return true; 
} 

// Find tangent points of line from point p to circle at c, radius r.
// By Thales' theorem via Wikipedia, the two tangent points are the points of 
// intersection of circle c with a circle centered on the midpoint
// of line cp and having cp as a diameter.
static bool find_tangent_points(const Vec2& c, Real r, const Vec2& p,
                                Vec2& p0, Vec2& p1)
{
    const Vec2 cp = p - c; // vector from c to p
    const Vec2 mid = (c+p)/2; // midpoint
    return circle_circle_intersection(c,r, mid, cp.norm()/2, p0, p1);
}

// Given a circle of radius r and a chord length h(<= d=2r), calculate the arc length
// s between them (going the shorter direction). We calculate the angle theta 
// subtended by the chord, then the arc length is r*theta. 
//      r*sin(theta/2)=h/2 => theta/2 = asin(h/d)
//      s = r*theta = d*asin(h/d)
static Real calc_arc_length(Real r, Real h) 
{
    const Real d = 2*r;
    assert(0 <= h && h <= d);
    return d==0 ? 0 : d * std::asin(h/d);
}

#define TEST
#ifdef TEST 
static void run_test(double x0, double y0, double r0, double x1, double y1, double r1) 
{   Vec2 p0, p1;
    printf("x0=%g, y0=%g, r0=%g, x1=%g, y1=%g, r1=%g :\n", x0, y0, r0, x1, y1, r1); 
    circle_circle_intersection(Vec2(x0, y0), r0, Vec2(x1, y1), r1, p0, p1); 
    printf(" p0=%g,%g, p1=%g,%g\n", p0[0],p0[1], p1[0],p1[1]); 
} 
static void testmain() {
    /* Add more! */ 
    run_test(0,0,4,4,0,4);
    run_test(0,0,4,-2,0,2);      // 1 point (-4,0 twice)
    run_test(0,0,4,-2,0,1.9999); // no points
    run_test(0,0,4,-2,0,2.0001); // two points near -4,0
    run_test(-1.0, -1.0, 1.5, 1.0, 1.0, 2.0); 
    run_test(1.0, -1.0, 1.5, -1.0, 1.0, 2.0); 
    run_test(-1.0, 1.0, 1.5, 1.0, -1.0, 2.0); 
    run_test(1.0, 1.0, 1.5, -1.0, -1.0, 2.0);  

    Vec2 p0,p1;
    find_tangent_points(Vec2(0,0),4,Vec2(-4,0), p0,p1);
    cout << "p0=" << p0 << " p1=" << p1 << endl; // 1 point (-4,0 twice)
    find_tangent_points(Vec2(0,0),4,Vec2(-4*Sqrt2,0), p0,p1);
    cout << "p0=" << p0 << " p1=" << p1 << endl; // -2 sqrt(2), +/- 2 sqrt(2)

    cout << "s=" << calc_arc_length(4, (p1-p0).norm()) << " (2pi?)\n";

} 
#endif 
