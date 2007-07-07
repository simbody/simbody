/* Portions copyright (c) 2007 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * This is a simple example of using a constraint.
 * Here we have two independent pendulums hanging from ground pins.
 * They can be connected by a spring or a distance constraint.
 */

#include "SimTKsimbody.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Pi      = (Real)SimTK_PI, 
                  Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

static const Transform GroundFrame;

static const Real m = 1;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters

class ShermsForce : public GeneralForceElements::CustomForce {
public:
    ShermsForce(const MobilizedBody& b1, const MobilizedBody& b2) : body1(b1), body2(b2) { }
    ShermsForce* clone() const {return new ShermsForce(*this);}

    void calc(const MatterSubsystem& matter, const State& state,
              Vector_<SpatialVec>& bodyForces,
              Vector_<Vec3>&       particleForces,
              Vector&              mobilityForces,
              Real&                pe) const
    {
        const Vec3& pos1 = body1.getBodyTransform(state).T();
        const Vec3& pos2 = body2.getBodyTransform(state).T();
        const Real d = (pos2-pos1).norm();
        const Real k = 1000, d0 = 1;
        const Vec3 f = k*(d-d0)*(pos2-pos1)/d;
        body1.applyBodyForce(state, SpatialVec(Vec3(0),  f), bodyForces);
        body2.applyBodyForce(state, SpatialVec(Vec3(0), -f), bodyForces);
    }
private:
    const MobilizedBody& body1;
    const MobilizedBody& body2;
};

int main(int argc, char** argv) {
  try { // If anything goes wrong, an exception will be thrown.

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem         mbs;

    SimbodyMatterSubsystem  twoPends(mbs);
    UniformGravitySubsystem gravity(mbs, Vec3(0, -g, 0));
    GeneralForceElements    forces(mbs);
    DecorationSubsystem     viz(mbs);

        // ADD BODIES AND THEIR MOBILIZERS
    MobilizedBody::Pin leftPendulum(twoPends.Ground(),
                                      Transform(Vec3(-1, 0, 0)),
                                    Body::Rigid(MassProperties(m, Vec3(0), Inertia(1))),
                                      Transform(Vec3(0, d, 0)));

    MobilizedBody::Pin rightPendulum(twoPends.Ground(),
                                     Body::Rigid(MassProperties(m, Vec3(0), Inertia(1))));

    rightPendulum.setDefaultInboardFrame(Vec3(1,0,0));
    rightPendulum.setDefaultOutboardFrame(Vec3(0,d,0));

    rightPendulum.setDefaultAngle(20*Deg2Rad);

    // Beauty is in the eye of the beholder ...
    viz.addBodyFixedDecoration(leftPendulum,  Transform(), DecorativeSphere(.1).setColor(Red));
    viz.addBodyFixedDecoration(rightPendulum, Transform(), DecorativeSphere(.1).setColor(Blue));

        // OPTIONALLY TIE TOGETHER WITH SPRING/DAMPER OR DISTANCE CONSTRAINT

    const Real distance = 2;      // nominal length for spring; length for constraint
    const Real stiffness = 100;   // only if spring is used
    const Real damping   = 10;     //          "

    char c;
    cout << "Constraint, spring, or nothing? c/s/n"; cin >> c;

    if (c == 'c')         
        Constraint::Rod(leftPendulum, Vec3(0),
                        rightPendulum, Vec3(0),
                        distance);
    else if (c == 's') {
        forces.addTwoPointLinearSpring(leftPendulum, Vec3(0),
                                       rightPendulum, Vec3(0),
                                       stiffness, distance);
        forces.addTwoPointLinearDamper(leftPendulum, Vec3(0),
                                       rightPendulum, Vec3(0),
                                       damping);
    }

    // Add visualization line (orange=spring, black=constraint)
    if (c=='c' || c=='s')
        viz.addRubberBandLine(leftPendulum, Vec3(0),
                              rightPendulum, Vec3(0),
                              DecorativeLine().setColor(c=='c' ? Black : Orange).setLineThickness(4));

    //forces.addMobilityConstantForce(rightPendulum, 0, 20);
    forces.addCustomForce(ShermsForce(leftPendulum,rightPendulum));

    State s = mbs.realizeTopology(); // returns a reference to the the default state

    mbs.realizeModel(s); // define appropriate states for this System

    // Create a study using the Runge Kutta Merson or CPODES integrator
    RungeKuttaMerson myStudy(mbs, s);
    //CPodesIntegrator myStudy(mbs, s);

    const Real dt = 0.005; // output intervals
    const Real finalTime = 5;

    leftPendulum.setAngle(s, -60*Deg2Rad);

    // TODO: this can't work unless it sets a state variable somewhere.
    // Cache entries can only be updated during a realize() operation.
    //rightPendulum.applyPinTorque(s, Stage::Instance, 2000);

    s.setTime(0);

    // visualize once before and after assembly
    VTKReporter display(mbs);
    display.report(s);
    cout << "Unassembled configuration shown. Ready to assemble? "; cin >> c;

    // Peforms assembly if constraints are violated.
    myStudy.initialize();

    display.report(s);
    cout << "Assembled configuration shown. Ready to simulate? "; cin >> c;

    for (;;) {
        mbs.realize(s);
        printf("%5g %10.4g %10.4g %10.8g h=%g\n", s.getTime(), 
            leftPendulum.getAngle(s)*Rad2Deg,
            rightPendulum.getAngle(s)*Rad2Deg,
            mbs.getEnergy(s), myStudy.getPredictedNextStep());

        cout << "Mobilizer X left =" << leftPendulum.getMobilizerTransform(s);
        cout << "Mobilizer X right=" << rightPendulum.getMobilizerTransform(s);

        cout << "Mobilizer V left =" << leftPendulum.getMobilizerVelocity(s) << endl;
        cout << "Mobilizer V right=" << rightPendulum.getMobilizerVelocity(s) << endl;

        display.report(s);
        if (s.getTime() >= finalTime)
            break;

        // TODO: should check for errors or have or teach RKM to throw. 
        myStudy.step(s.getTime() + dt);
    }

  } 
  catch (const exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
  catch (...) {
    printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

}

