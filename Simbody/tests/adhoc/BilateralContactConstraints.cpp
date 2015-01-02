/* -------------------------------------------------------------------------- *
 *           Simbody(tm) Adhoc Test: Bilateral Contact Constraints            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

/* This uses a variety of bilateral constraints that are intended as the
underpinnings for unilateral constraints. The most important check here is
that energy should be conserved perfectly (to integration accuracy) since
these are all non-working constraints. (Look at the total energy in the 
visualizer.) Even the friction constraints are
non-working because the underlying constraints represent rolling (a.k.a.
"stiction"); sliding is imposed elsewhere by disabling the rolling constraints
and replacing them with different conditions. */

#include "Simbody.h"
using namespace SimTK;

#include <iostream>
using std::cout; using std::endl;

class ShowEnergy : public DecorationGenerator {
public:
    explicit ShowEnergy
       (const MultibodySystem& mbs,
        const MobilizedBody::Free& brick,
        const Array_<Constraint::SphereOnPlaneContact>& balls,
        const Array_<Constraint::SphereOnSphereContact>& sphsph, 
        const Array_<Constraint::Rod>& rods) 
    :   m_mbs(mbs), m_brick(brick), m_balls(balls), m_sphsph(sphsph),
        m_rods(rods) {}

    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) override
    {
        const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
        const Real TextScale = m_mbs.getDefaultLengthScale()/10; // was .1
        m_mbs.realize(state, Stage::Dynamics);
        const Real KE=m_mbs.calcKineticEnergy(state), E=m_mbs.calcEnergy(state);
        DecorativeText energy; energy.setIsScreenText(true);
        energy.setText("Energy/KE: " + String(E, "%.6f") + String(KE, "/%.6f"));
        geometry.push_back(energy);

        //cout << "brick q=" << m_brick.getQAsVector(state) << endl;
        //cout << "brick u=" << m_brick.getUAsVector(state) << endl;

        m_mbs.realize(state, Stage::Acceleration);
        for (unsigned i=0; i < m_balls.size(); ++i) {
            const Vec3 f_GC = m_balls[i].findForceOnSphereInG(state);
            const Vec3 p_GC = m_balls[i].findContactPointInG(state);

            geometry.push_back(
                DecorativeLine(p_GC - f_GC, p_GC).setColor(Red));

            DecorativeText sep; sep.setIsScreenText(true);
            sep.setText(String(i) + ": " +
                        String(m_balls[i].findSeparation(state), "%.6f"));
            geometry.push_back(sep);
            sep.setText("  : " +
                        String(m_balls[i].getVelocityErrors(state)));
            geometry.push_back(sep);
            sep.setText("  : " +
                        String(m_balls[i].getAccelerationErrors(state)));
            geometry.push_back(sep);
        }

        for (unsigned i=0; i < m_sphsph.size(); ++i) {
            const Vec3 f_GC = m_sphsph[i].findForceOnSphereBInG(state);
            const Transform X_GC = m_sphsph[i].findContactFrameInG(state);

            geometry.push_back(
                DecorativeFrame().setTransform(X_GC).setColor(Purple));

            geometry.push_back(
                DecorativeLine(X_GC.p() - f_GC, X_GC.p()).setColor(Red));

            DecorativeText sep; sep.setIsScreenText(true);
            sep.setText(String(i) + ": " +
                        String(m_sphsph[i].findSeparation(state), "%.6f"));
            geometry.push_back(sep);
            sep.setText("  : " +
                        String(m_sphsph[i].getVelocityErrors(state)));
            geometry.push_back(sep);
            sep.setText("  : " +
                        String(m_sphsph[i].getAccelerationErrors(state)));
            geometry.push_back(sep);

            //const SimbodyMatterSubsystem& matter=m_mbs.getMatterSubsystem();
            //State s2 = state;
            //m_mbs.realize(s2,Stage::Acceleration);
            //Vector udot=s2.getUDot();
            //Vec3 aerr(0);
            //const Real du = 1e-6;
            //for (int u=0; u < matter.getNumMobilities(); ++u) {
            //    s2.updU()[u] += du; m_mbs.realize(s2,Stage::Velocity);
            //    Vec3 verrp = m_sphsph[i].getVelocityErrors(s2);
            //    s2.updU()[u] -= 2*du; m_mbs.realize(s2,Stage::Velocity);
            //    Vec3 verrm = m_sphsph[i].getVelocityErrors(s2);
            //    aerr += ((verrp-verrm) / (2*du)) * udot[u];
            //    s2.updU()[u]=state.getU()[u];
            //}
            //printf("aerr =%.15g %.15g\n", aerr[0], aerr[1]);
            //m_mbs.realize(s2, Stage::Acceleration);
        }

        for (unsigned i=0; i < m_rods.size(); ++i) {
            const Constraint::Rod& rod = m_rods[i];
            const Real t = rod.getRodTension(state);
            const UnitVec3 d = rod.findRodOrientationInG(state); // p1->p2
            const Vec3 f1_G = t*d; // force on p1
            const Vec3 f2_G = -f1_G; // force on p2
            const Vec3 p2 = rod.getPointOnBody2(state);
            const Vec3 p2_G = rod.getMobilizedBody2().
                findStationLocationInGround(state, p2);

            geometry.push_back(
                DecorativeLine(p2_G - f2_G, p2_G)
                .setColor(Red).setLineThickness(5));
        }
    }
private:
    const MultibodySystem&      m_mbs;
    const MobilizedBody::Free   m_brick;
    const Array_<Constraint::SphereOnPlaneContact>&  m_balls;
    const Array_<Constraint::SphereOnSphereContact>& m_sphsph;
    const Array_<Constraint::Rod>&                   m_rods;
};

int main() {
    // Define the system.
    MultibodySystem        system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem  forces(system);
    Force::Gravity         gravity(forces, matter, -YAxis, 9.8/10);

    //Force::GlobalDamper damp(forces, matter, 1); 

    // Describe mass and visualization properties for a generic body.
    Real mass = 2;
    Vec3 hdim(1,.5,.25);
    Body::Rigid bodyInfo(MassProperties(mass, Vec3(0), UnitInertia::brick(hdim)));
    bodyInfo.addDecoration(Transform(), 
        DecorativeBrick(hdim).setColor(Orange).setOpacity(.3));

    Real pmass = .1;
    Vec3 phdim(5,.5,2);
    Body::Rigid platformBody(MassProperties(10*mass,
        Vec3(0), UnitInertia::ellipsoid(phdim))); 
    platformBody.addDecoration(Transform(),
        DecorativeEllipsoid(phdim).setColor(Cyan).setOpacity(.1)
        .setResolution(5));

    MobilizedBody::Ball platform(matter.Ground(), Vec3(0),
                                 platformBody, phdim/2);
    //MobilizedBody platform = matter.Ground();

    // Create the moving (mobilized) bodies of the pendulum.
    //MobilizedBody::Free brick(platform, Transform(Vec3(0)),
    //                          bodyInfo,        Transform(Vec3(0)));
    MobilizedBody::Free brick(matter.Ground(), Transform(Vec3(0)),
                              bodyInfo,        Transform(Vec3(0)));


    Array_<Constraint::SphereOnPlaneContact> balls;
    Array_<Constraint::SphereOnSphereContact> sphsph;
    Array_<Constraint::Rod> rods;

    Rotation ZtoY(-Pi/2, XAxis);
    //Constraint::PointInPlaneWithStiction pt1(platform, 
    //                                     Transform(ZtoY, Vec3(0,1,0)),
    //                                     brick, hdim);
    //pt1.setPlaneDisplayHalfWidth(5);
    //Constraint::SphereOnPlaneContact ball1(platform, 
    //                                       Transform(ZtoY, Vec3(0,1,0)),
    //                                       brick, hdim, 0.5, false);
    //ball1.setPlaneDisplayHalfWidth(5);
    //balls.push_back(ball1); 

    //Constraint::SphereOnPlaneContact ball2(brick, 
    //                                       Transform(Vec3(0,0,-hdim[2])),
    //                                       platform, -phdim/2, 0.5, false);
    //ball2.setPlaneDisplayHalfWidth(5);
    //balls.push_back(ball2);

    //Constraint::SphereOnPlaneContact ball3(brick, 
    //                                       Transform(Vec3(0,0,-hdim[2])),
    //                                       platform, Vec3(-2,3,-.5), .7, true);
    //ball3.setPlaneDisplayHalfWidth(5);
    //balls.push_back(ball3);


    //MobilizedBody::Free ball(matter.Ground(), Vec3(0),
    //                         MassProperties(1,Vec3(0),UnitInertia(1,1,1)),
    //                         Vec3(0));
    //Constraint::SphereOnSphereContact ss(platform, Vec3(-2,1,-.5), .7,
    //                                     ball, Vec3(0), 1.2, true);                                         
    //Constraint::SphereOnSphereContact bb(brick, hdim, 0.5,
    //                                     ball, Vec3(0), 1.2, true);
    //sphsph.push_back(bb);

    Constraint::SphereOnSphereContact ss(brick, hdim, 0.5,
                                         platform, Vec3(-3,1,-.5), 1.2, 
                                         false);
    sphsph.push_back(ss);

    //Constraint::SphereOnSphereContact ss(platform, Vec3(-2,3,-.5), .7, 
    //                                     brick, hdim, 0.5, false);
    //Constraint::SphereOnSphereContact ss(platform, Vec3(-2,3,-.5), .7, 
    //                                     brick, hdim, 0.5, false);
    //Constraint::SphereOnSphereContact ss(brick, hdim, 0.5, 
    //                                     matter.Ground(), Vec3(-2,3,-.5), .7,true);
    Constraint::Rod rod1(brick, Vec3(0,hdim[1],hdim[2]), 
                         platform, Vec3(0,3,-.5), 1.5*1.2);
    

    // Spring to keep the brick near 000.
    //Force::TwoPointLinearSpring(forces, platform, Vec3(0),
    //                           brick, Vec3(0), 4, 1);

    // Rod to keep the brick near 000.
    //Constraint::Rod rod1(platform, Vec3(0,0,2),
    //                     brick, -hdim, 3);
    //rods.push_back(rod1);

    // Try edge/edge contact.
    Constraint::LineOnLineContact ll(platform, 
          Transform(Rotation(UnitVec3(1,1,1), XAxis, UnitVec3(-XAxis), ZAxis), 
                    Vec3(1,1,1)),
          2, // hlen
                                     brick, 
          Transform(Rotation(UnitVec3(ZAxis), XAxis, Vec3(-1,-1,0), ZAxis),
                    Vec3(-hdim[0],-hdim[1],0)),
          2, // hlen
          true);

    // Set up visualization at 30 fps.
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::SolidColor);
    viz.setShowFrameRate(true);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    // Initialize the system and acquire default state.
    State state = system.realizeTopology();
    brick.setQToFitTransform(state, Vec3(0,5,0));
    brick.setUToFitAngularVelocity(state, Vec3(10,10,10));

    //rod1.setRodLength(state, 5);

    viz.report(state); 

    printf("Initial config. Ready to assemble.\n"); getchar();
    Assembler asmb(system);
    asmb.assemble(state);

    viz.report(state);
    printf("Assembled. Ready to initialize.\n"); getchar();

    //printf("Changed ball3 from rad=%g to rad=%g\n",
    //       ball3.getSphereRadius(state), 1.5);
    //ball3.setSphereRadius(state, 1.5);
    //viz.report(state); getchar();
    //asmb.assemble(state);

    //viz.report(state);
    //printf("Re-assembled. Ready to simulate.\n"); getchar();

    // Choose integrator and simulate for 10 seconds.
    RungeKuttaMersonIntegrator integ(system);
    //RungeKutta3Integrator integ(system);
    integ.setAccuracy(1e-8);
    //integ.setConstraintTolerance(1e-3);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    viz.report(ts.getState());
    printf("Initialized. Ready to simulate.\n"); getchar();
    viz.addDecorationGenerator(new ShowEnergy(system,brick,balls,sphsph,rods));
    ts.stepTo(100.0);
    printf("# steps=%d/%d\n", 
           integ.getNumStepsTaken(), integ.getNumStepsAttempted());
}