/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
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

#include "SimTKsimbody.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

static const Transform GroundFrame;

static const Real m = 5;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters

int main(int argc, char** argv) {
try { // If anything goes wrong, an exception will be thrown.

    MultibodySystem         mbs;
    GeneralForceSubsystem    forces(mbs);
    SimbodyMatterSubsystem  pend(mbs);
    DecorationSubsystem     viz(mbs);
    Force::UniformGravity gravity(forces, pend, Vec3(0, -g, 0));

    MobilizedBody::Ball connector(pend.Ground(), 
                                    Transform(1*Vec3(0, 0, 0)),
                                  MassProperties(1, Vec3(0,0,0), Inertia(10,20,30)),
                                    Transform(1*Vec3(0, .5, 0)));

    connector.setDefaultRadius(0.05); // for the artwork

    //connector.setDefaultRotation( Rotation(Pi/4, Vec3(0,0,1) );
 
    const Real m1 = 5;
    const Real m2 = 1;
    const Real radiusRatio = std::pow(m2/m1, 1./3.);
    const Vec3 weight1Location(0, 0, -d/2); // in local frame of swinging body
    const Vec3 weight2Location(0, 0,  d/2); // in local frame of swinging body
    const Vec3 COM = (m1*weight1Location+m2*weight2Location)/(m1+m2);

    const MassProperties swingerMassProps
        (m1+m2, COM, 1*Inertia(1,1,1) + m1*UnitInertia::pointMassAt(weight1Location)
                                      + m2*UnitInertia::pointMassAt(weight2Location));
    MobilizedBody::Screw swinger(connector, 
                                    Transform( Rotation( 0*.7, Vec3(9,8,7) ),
                                              1*Vec3(0,-.5,0)),
                                 swingerMassProps,
                                    Transform( Rotation(0*1.3, Vec3(0,0,1) ),
                                              COM+0*Vec3(0,0,3)),    // inboard joint location
                                 0.3); // pitch

    // Add a blue sphere around the weight.
    viz.addBodyFixedDecoration(swinger, weight1Location, 
          DecorativeSphere(d/8).setColor(Blue).setOpacity(.2));
    viz.addBodyFixedDecoration(swinger, weight2Location, 
          DecorativeSphere(radiusRatio*d/8).setColor(Green).setOpacity(.2));
    viz.addRubberBandLine(GroundIndex, Vec3(0),
                          swinger, Vec3(0),
                          DecorativeLine().setColor(Blue).setLineThickness(10)
                                          .setRepresentation(DecorativeGeometry::DrawPoints));

    //forces.addMobilityConstantForce(swinger, 0, 10);
    Force::ConstantTorque(forces, swinger, Vec3(0,0,10));
    //forces.addConstantForce(swinger, Vec3(0), Vec3(0,10,0));
    //forces.addConstantForce(swinger, Vec3(0,0,0), Vec3(10,10,0)); // z should do nothing
    //forces.addMobilityConstantForce(swinger, 1, 10);
   // forces.addMobilityConstantForce(swinger, 2, 60-1.2);

    State s = mbs.realizeTopology(); // define appropriate states for this System


    pend.setUseEulerAngles(s, true);
    mbs.realizeModel(s);

    mbs.realize(s);

    // Create a study using the Runge Kutta Merson integrator
    RungeKuttaMersonIntegrator myStudy(mbs);
    myStudy.setAccuracy(1e-6);

    // This will pick up decorative geometry from
    // each subsystem that generates any, including of course the 
    // DecorationSubsystem, but not limited to it.
    Visualizer display(mbs);


    const Real expectedPeriod = 2*Pi*std::sqrt(d/g);
    printf("Expected period: %g seconds\n", expectedPeriod);

    const Real dt = 1./30; // output intervals
    const Real finalTime = 1*expectedPeriod;

    for (Real startAngle = 10; startAngle <= 90; startAngle += 10) {
        s = mbs.getDefaultState();
        connector.setQToFitRotation(s, Rotation(Pi/4, Vec3(1,1,1)) );

        printf("time  theta      energy           *************\n");

        swinger.setQToFitTransform(s, Transform( Rotation( BodyRotationSequence, 0*Pi/2, XAxis, 0*Pi/2, YAxis ), Vec3(0,0,0) ));
        swinger.setUToFitVelocity(s, SpatialVec(0*Vec3(1.1,1.2,1.3),Vec3(0,0,-1)));

        s.updTime() = 0;
        myStudy.initialize(s);

        cout << "MassProperties in B=" << swinger.expressMassPropertiesInAnotherBodyFrame(myStudy.getState(),swinger);
        cout << "MassProperties in G=" << swinger.expressMassPropertiesInGroundFrame(myStudy.getState());
        cout << "Spatial Inertia    =" << swinger.calcBodySpatialInertiaMatrixInGround(myStudy.getState());

        int stepNum = 0;
        for (;;) {
            // Should check for errors and other interesting status returns.
            myStudy.stepTo(myStudy.getTime() + dt);
            const State& s = myStudy.getState(); // s is now the integrator's internal state

            if ((stepNum++%10)==0) {
                // This is so we can calculate potential energy (although logically
                // one should be able to do that at Stage::Position).
                mbs.realize(s, Stage::Dynamics);

                cout << s.getTime() << ": E=" << mbs.calcEnergy(s) 
                     << " connector q=" << connector.getQ(s) 
                     << ": swinger q=" << swinger.getQ(s) << endl;

                // This is so we can look at the UDots.
                mbs.realize(s, Stage::Acceleration);

                cout << "q =" << pend.getQ(s) << endl;
                cout << "u =" << pend.getU(s) << endl;
                cout << "ud=" << pend.getUDot(s) << endl;

                cout << "Connector V=" << connector.getMobilizerVelocity(s) << endl;
                cout << "Swinger V=" << swinger.getMobilizerVelocity(s) << endl;

                const Rotation& R_MbM = swinger.getMobilizerTransform(s).R();
                Vec4 aaMb = R_MbM.convertRotationToAngleAxis();
                cout << "angle=" << aaMb[0] << endl;
                cout << "axisMb=" << aaMb.drop1(0) << endl;
                cout << "axisMb=" << ~R_MbM*aaMb.drop1(0) << endl;
            }

            display.report(s);
            if (s.getTime() >= finalTime)
                break;
        }
    }
} 
catch (const exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
}
}

