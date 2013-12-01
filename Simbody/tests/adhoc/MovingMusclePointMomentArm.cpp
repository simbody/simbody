/* -------------------------------------------------------------------------- *
 *           Simbody(tm) Adhoc test: Moving Muscle Point Moment Arm           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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


#include "Simbody.h"
#include <iostream>

using namespace SimTK;
using std::cout; using std::endl;

// Uncomment this to use an explicit rack-and-pinion mechanism rather than
// the moving muscle point.
#define USE_RACK

/*
In biomechanics muscle forces are modeled as acting on the system via 
frictionless cables anchored at an "origin point" on one bone and following a
curved path over obstacles to a final anchor at an "insertion point" on another
bone. The cables have uniform tension and apply forces to the end points and
the obstacles over which they pass.

Because it can be very expensive to calculate the actual path, obstacle 
surfaces are often replaced by simplified representations involving "via points"
(frictionless eyelets fixed to bones), or "moving muscle points" (MMPs). An MMP
is a via point that moves around on its body's surface, with its location P(q)
given in its body's local frame as a smooth kinematic function of a designated 
generalized coordinate q. The function P(q) is typically a vector-valued spline 
fit through point locations measured at sampled coordinate values taken on a 
cadaver or a more complex computational model.

We would like to ensure that the moment arm and dynamics we calculate using the
reduced model with MMPs is the same as we would have gotten with the more
complex model, assuming that P(q) is the same in both models. (Reaction forces
will unavoidably differ, but accelerations should be the same.) One way to do 
that is to calculate moment arm by perturbation r(theta)=dL/dtheta where L is
the length of the muscle path. We would like to be able to obtain the same value
for r(theta) with an instantaneous calculation r(theta)=tau_theta/s where s
is the tension in the muscle path and tau_theta is the generalized torque 
produced about theta by that tension. The two values are equivalent if all
constraints in the model are workless. For an MMP, there must be a workless
constraint that accounts for its motion. In the real system, that will be
produced by the mechanical contacts and tendons that form the joint. In the MMP 
model we don't have that mechanical system present but would like to treat it as 
though motion were caused by an equivalent "gearbox" driven by q.

In this example, we will build a simplified knee-like mechanism in which the
gearbox is explicitly present, and then attempt to get the same moment arm from
a simplified model in which only P(q) is known. The mechanism looks like this:

                  / femur
               .I/            I=insertion point fixed on femur
          M .   /             M=muscle path of interest (dots)
         .     @  q
       P---<--===-->--- rack
         .     |
            .  | tibia
              .O              O=origin point fixed on tibia (fixed to Ground)
               |
              ||| Ground

We have built a rack mechanism to replace the patella for determining where
point P is located as a function of q. Specifically, given a pitch we can find
P=(Px,Py,Pz)=(P0x + pitch*q, P0y, P0z) and dP/dq=(pitch,0,0). The muscle path
length L=|P-O|+|I-P|.
*/

/* This muscle uses a via point fixed to the rack body. */
class MuscleVP : public Force::Custom::Implementation {
public:
    MuscleVP(const SimbodyMatterSubsystem& matter,
        const MobilizedBody& A, const Vec3& origin,
        const MobilizedBody& V, const Vec3& via,
        const MobilizedBody& B, const Vec3& insertion,
        Real stiffness, Real zeroLength)
    :   m_matter(matter), m_A(A), m_ptA(origin), m_V(V), m_ptV(via),
        m_B(B), m_ptB(insertion), m_k(stiffness), m_zero(zeroLength),
        m_tension(NaN)
    {
    }

    Real calcLength(const State& state) const {
        // End points and via point in Ground.
        const Vec3 ptA = m_A.findStationLocationInGround(state, m_ptA);
        const Vec3 ptB = m_B.findStationLocationInGround(state, m_ptB);
        const Vec3 ptV = m_V.findStationLocationInGround(state, m_ptV);

        const Vec3 A2V(ptV - ptA), B2V(ptV - ptB);
        const Real Alen = A2V.norm(), Blen = B2V.norm();
        const Real len = Alen + Blen;
        return len;
    }

    Real calcTension(const State& state) const {
        const Real len = calcLength(state);
        if (len <= m_zero) return 0;
        return m_k*(len-m_zero);
    }

    // Set to NaN to enable k*x tension instead. Don't forget to invalidate
    // the state.
    void setTension(Real tension) {m_tension=tension;}
    Real getTension() const { return m_tension; }

    // Calculate the muscle forces and accumulate into bodyForces array.
    // (We aren't going to generate any particle or mobility forces.)
    void calcForce(const State& state,
        Vector_<SpatialVec>& bodyForces,
        Vector_<Vec3>& particleForces,
        Vector& mobilityForces) const OVERRIDE_11
    {
        // End points and via point in Ground.
        const Vec3 ptA = m_A.findStationLocationInGround(state, m_ptA);
        const Vec3 ptB = m_B.findStationLocationInGround(state, m_ptB);
        const Vec3 ptV = m_V.findStationLocationInGround(state, m_ptV);

        const Vec3 A2V(ptV - ptA), B2V(ptV - ptB);
        const Real Alen = A2V.norm(), Blen = B2V.norm();
        const Real len = Alen + Blen;

        const UnitVec3 uA2V(A2V / Alen, true), uB2V(B2V / Blen, true);

        const Real tension = isNaN(m_tension)
            ? ((len>m_zero) ? m_k * (len - m_zero) : Real(0))
            : m_tension;

        const Vec3 fA = tension*uA2V;
        const Vec3 fB = tension*uB2V;

        m_A.applyForceToBodyPoint(state, m_ptA, fA, bodyForces);
        m_B.applyForceToBodyPoint(state, m_ptB, fB, bodyForces);
        m_V.applyForceToBodyPoint(state, m_ptV, -(fA + fB), bodyForces);
    }

    Real calcPotentialEnergy(const State&) const OVERRIDE_11 { return 0; }
private:
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBody             m_A, m_V, m_B;
    const Vec3                      m_ptA, m_ptV, m_ptB;
    const Real                      m_k, m_zero;

    Real m_tension; // NaN to calculate as k*x
};

/* This muscle uses a differentiable function 
        P(q)_A=P0_A + (pitch, 0, 0)*q 
to determine the location of a moving muscle point on the origin body A;
there is no rack body. */
class MuscleMMP : public Force::Custom::Implementation {
public:
    MuscleMMP(const SimbodyMatterSubsystem& matter,
        const MobilizedBody& A, const Vec3& origin,
        const MobilizedBody& B, const Vec3& insertion,
        const Vec3& P0_A, Real pitch,
        Real stiffness, Real zeroLength)
        : m_matter(matter), m_A(A), m_ptA(origin), m_B(B), m_ptB(insertion), 
        m_P0(P0_A), m_pitch(pitch),
        m_k(stiffness), m_zero(zeroLength), m_tension(NaN)
    {
    }

    // Calculate P(q), in A frame.
    Vec3 calcP(const State& state) const {
        const Real q = m_B.getOneQ(state, MobilizerQIndex(0));
        return m_P0 + Vec3(m_pitch*q, 0, 0);
    }

    // Calculate dP/dq.
    Vec3 calcdPdq(const State& state) const {
        return Vec3(m_pitch, 0, 0);
    }

    void calcPathPoints(const State& state, Array_<Vec3>& pts_G) const {
        // End points and moving muscle point in Ground.
        const Vec3 ptA = m_A.findStationLocationInGround(state, m_ptA);
        const Vec3 ptP = m_A.findStationLocationInGround(state, calcP(state));
        const Vec3 ptB = m_B.findStationLocationInGround(state, m_ptB);
        pts_G.clear();
        pts_G.push_back(ptA); pts_G.push_back(ptP); pts_G.push_back(ptB);
    }


    Real calcLength(const State& state) const {
        // End points and via point in Ground.
        const Vec3 ptA = m_A.findStationLocationInGround(state, m_ptA);
        const Vec3 ptP = m_A.findStationLocationInGround(state, calcP(state));
        const Vec3 ptB = m_B.findStationLocationInGround(state, m_ptB);

        const Vec3 A2P(ptP - ptA), B2P(ptP - ptB);
        const Real Alen = A2P.norm(), Blen = B2P.norm();
        const Real len = Alen + Blen;
        return len;
    }

    Real calcTension(const State& state) const {
        const Real len = calcLength(state);
        if (len <= m_zero) return 0;
        return m_k*(len - m_zero);
    }

    // Set to NaN to enable k*x tension instead. Don't forget to invalidate
    // the state.
    void setTension(Real tension) { m_tension = tension; }
    Real getTension() const { return m_tension; }

    // Calculate the muscle forces and accumulate into bodyForces Vector, plus
    // a correction force that is added to the mobilityForces Vector.
    void calcForce(const State& state,
        Vector_<SpatialVec>& bodyForces,
        Vector_<Vec3>& particleForces,
        Vector& mobilityForces) const OVERRIDE_11
    {
        // End points and via point in Ground.
        const Vec3 ptA = m_A.findStationLocationInGround(state, m_ptA);
        const Vec3 ptB = m_B.findStationLocationInGround(state, m_ptB);
        const Vec3 ptP_A = calcP(state);
        const Vec3 ptP = m_A.findStationLocationInGround(state, ptP_A);

        const Vec3 A2P(ptP - ptA), B2P(ptP - ptB);
        const Real Alen = A2P.norm(), Blen = B2P.norm();
        const Real len = Alen + Blen;

        const UnitVec3 uA2P(A2P / Alen, true), uB2P(B2P / Blen, true);

        const Real tension = isNaN(m_tension)
            ? ((len>m_zero) ? m_k * (len - m_zero) : Real(0))
            : m_tension;

        const Vec3 fA = tension*uA2P;
        const Vec3 fB = tension*uB2P;
        const Vec3 fP = -(fA + fB);
        const Real f = ~calcdPdq(state)*fP;

        m_A.applyForceToBodyPoint(state, m_ptA, fA, bodyForces);
        m_B.applyForceToBodyPoint(state, m_ptB, fB, bodyForces);
        m_A.applyForceToBodyPoint(state, ptP_A, fP, bodyForces);
        m_B.applyOneMobilityForce(state, MobilizerQIndex(0), f, mobilityForces); 
    }

    Real calcPotentialEnergy(const State&) const OVERRIDE_11 { return 0; }
private:
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBody             m_A, m_V, m_B;
    const Vec3                      m_ptA, m_ptB, m_P0;
    const Real                      m_pitch, m_k, m_zero;

    Real    m_tension;
};

class DrawPath : public DecorationGenerator {
public:
    explicit DrawPath(const MuscleMMP& muscle)
    :   m_muscle(muscle) {}

    void generateDecorations(const State& state, 
        Array_<DecorativeGeometry>& geometry) OVERRIDE_11{
        Array_<Vec3> path;
        m_muscle.calcPathPoints(state, path);
        for (unsigned i = 1; i < path.size(); ++i) {
            geometry.push_back(DecorativeLine(path[i - 1], path[i])
                .setColor(Red).setLineThickness(3));
            if (i + 1 < path.size())
                geometry.push_back(DecorativePoint(path[i]).setColor(Cyan));
        }
    }
private:
    const MuscleMMP&        m_muscle;
};

//==============================================================================
//                                   MAIN
//==============================================================================
int main() {
  try {
    // Create the system.
    // ------------------
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    // Add gravity as a force element.
    Force::Gravity gravity(forces, matter, -YAxis, 9.8);

    // Create body specifications.
    // ---------------------------
    Vec3 tibiaSz(.07, 1., .05), femurSz(.1, 1., .05), rackSz(1, .02, .02);
    Vec3 origin_T(-tibiaSz[0] / 2, -tibiaSz[1] / 6, 0);
    Vec3 insertion_F(-femurSz[0] / 2, 0.9*femurSz[1] / 2, 0);
    Vec3 mmp_R(-rackSz[0]/2, 0, 0);
    // This is where mmp is in the tibia frame when q=0.
    Real rackHtInTibia = 0.9*tibiaSz[1]/2;
    Vec3 P0_A(-rackSz[0]/2, rackHtInTibia, 0);

    Body::Rigid tibiaBody(MassProperties(1, Vec3(0),
        UnitInertia::brick(tibiaSz / 2.)));
    tibiaBody.addDecoration(DecorativeBrick(tibiaSz / 2.).setColor(Green));
    tibiaBody.addDecoration(DecorativePoint(origin_T));

    Body::Rigid femurBody(MassProperties(1, Vec3(0),
        UnitInertia::brick(femurSz / 2.)));
    femurBody.addDecoration(DecorativeBrick(femurSz / 2.).setColor(Blue));
    femurBody.addDecoration(DecorativePoint(insertion_F));

    Body::Rigid rackBody(MassProperties(.01, Vec3(0),
        UnitInertia::brick(rackSz/2.)));
    rackBody.addDecoration(DecorativeBrick(rackSz / 2.).setColor(Red));

    // Create mobilized bodies.
    // ------------------------
    MobilizedBody::Weld tibia(
        matter.Ground(), Vec3(0),
        tibiaBody, Vec3(0, -tibiaSz[1] / 2, 0));
    MobilizedBody::Pin femur(
        tibia, Vec3(0, tibiaSz[1] / 2, 0),
        femurBody, Vec3(0, -femurSz[1] / 2, 0));

    // Visualize.
    // ----------
    Visualizer viz(system);
    const Real pitch = .3; // m/radian

#ifdef USE_RACK
    MobilizedBody::Slider rack(
        tibia, Vec3(0, rackHtInTibia, 0),
        rackBody, Vec3(0));
    // Add rack & pinion constraint.
    Array_<MobilizedBodyIndex> mobods(2);
    Array_<MobilizerQIndex> coords(2);
    mobods[0] = femur; coords[0] = MobilizerQIndex(0);
    mobods[1] = rack; coords[1] = MobilizerQIndex(0);
    Constraint::CoordinateCoupler(matter,
        new Function::Linear(Vector(Vec3(-pitch,1,0))),
        mobods, coords);
    MuscleVP* mmp = new MuscleVP(matter, tibia, origin_T, rack, mmp_R, 
                                 femur, insertion_F, 100., 2.75);

    viz.addRubberBandLine(tibia, origin_T, rack, mmp_R, 
        DecorativeLine().setColor(Red));
    viz.addRubberBandLine(femur, insertion_F, rack, mmp_R,
        DecorativeLine().setColor(Red));
#else
    MuscleMMP* mmp = new MuscleMMP(matter, tibia, origin_T, femur, insertion_F,
        P0_A, pitch, 100., 2.75);
    viz.addDecorationGenerator(new DrawPath(*mmp));
#endif

    Force::Custom muscle(forces, mmp);

    // Add some damping.
    // -----------------
    //Force::MobilityLinearDamper(forces, femur, MobilizerQIndex(0), 1.);

    // Report at the framerate (real-time).
    system.addEventReporter(new Visualizer::Reporter(viz, 1. / 30));

    // Initialize the system and state.
    // --------------------------------
    State state = system.realizeTopology();
    femur.lockAt(state, -Pi/4);
    system.projectQ(state, 1e-10);

    // Calculate moment arm by dL/dq.
    system.realize(state);
    Real q0 = femur.getAngle(state);
    Real L0 = mmp->calcLength(state);
    femur.lockAt(state, q0 + 1e-6);
    system.projectQ(state, 1e-10);
    Real q1 = femur.getAngle(state);
    Real L1 = mmp->calcLength(state);
    Real r = (L1 - L0) / (q1-q0);
    printf("q1-q0=%g, L1-L0=%g, r=%g\n",
        q1 - q0, L1 - L0, r);
    femur.unlock(state);

    printf("Assembled:\n");
    viz.report(state);
    getchar();

    //----------------------------
    // CALCULATE COUPLING MATRIX C
    //----------------------------
    state.updU() = 0;
    femur.lockAt(state, 1., Motion::Velocity);
    system.realize(state, Stage::Velocity);
    cout << "before project u=" << state.getU() << endl;
    system.projectU(state, 1e-10);
    cout << "after project u=" << state.getU() 
         << " uerr=" << state.getUErr() << endl;
    const Vector C(state.getU());
    femur.unlock(state);
    state.updU() = 0;

    Vector_<SpatialVec> bodyForces;
    Vector_<Vec3> particleForces;
    Vector mobilityForces;
    mmp->setTension(1); // override k*x spring force with unit tension
    state.invalidateAllCacheAtOrAbove(Stage::Velocity);
    system.realize(state, Stage::Velocity);
    muscle.calcForceContribution(state, bodyForces, particleForces,
                                 mobilityForces);
    mmp->setTension(NaN); // back to k*x spring
    cout << "bodyForces=" << bodyForces << endl;
    cout << "mobilityForces=" << mobilityForces << endl;

    Vector equivForces;
    matter.multiplyBySystemJacobianTranspose(state, bodyForces, equivForces);
    equivForces += mobilityForces;
    cout << "tension=1 --> equivForces=" << equivForces << endl;

    Real r2 = ~C*equivForces;
    printf("gen force r=%g\n", r2);
    getchar();

    // Simulate.
    // ---------
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5.0);

    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
};