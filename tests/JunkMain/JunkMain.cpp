/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
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
 * This is just a disposable outer block for tests during development of
 * Simbody. Don't include this in the nightly test suite.
 */

#include "SimTKcommon.h"
#include "Simbody.h"

#include "simbody/internal/DecorativeGeometry.h"
#include "simbody/internal/VTKReporter.h"
#include "simbody/internal/NumericalMethods.h"

#include "simbody/internal/DuMMForceFieldSubsystem.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;


static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;
static const Real EnergyUnitsPerKcal = 418.4; // exact 
static const int  Ground = 0;       // ground is always body 0
static const Transform BodyFrame;   // identity transform on any body

class MySinusoid: public GeneralForceElements::UserForce {
public:
    MySinusoid(int b, int d, const Real& amp, const Real& w, const Real& ph=0) 
      : body(b), dof(d), amplitude(amp), period(w), phase(ph)
    {
    }

    // Implementation of pure virtual.
    void calc(const MatterSubsystem& matter, 
              const State&           state,
              Vector_<SpatialVec>&   bodyForces,
              Vector_<Vec3>&         particleForces,
              Vector&                mobilityForces,
              Real&                  pe) const 
    {
        matter.addInMobilityForce(state,body,dof,
            amplitude*std::sin(2*Pi*period*state.getTime() + phase),
            mobilityForces);
    }

    // Implementation of pure virtual;
    GeneralForceElements::UserForce* clone() const { 
        return new MySinusoid(*this); 
    }
private:
    int  body, dof;
    Real amplitude, period, phase;
};


int main() {
try {
    SimbodyMatterSubsystem   bendStretchBlock;
    GeneralForceElements     forces;

    MultibodySystem mbs;
    mbs.setMatterSubsystem(bendStretchBlock);
    mbs.addForceSubsystem(forces);

    //forces.addGlobalEnergyDrain(100);


    // System modeled:
    //                
    //     <--          
    //        \           
    //      *-/-------->*
    //                    
    //
    //      | y
    //      |
    //      |-----> x
    //     /
    //    z
    // This is a mass mounted to a rotating base with a "bend stretch" joint which
    // permits rotation around mutual z axes of base & ground, followed
    // by a translation along the *block* x axis (which will have moved).
    // Base body frame is attached to ground frame by pin about z.
    //

    const Real m = 1;
    const Vec3 com = Vec3(0,0,0);
    const Inertia iner(1,1,1);

    const MassProperties mProps(m,com,iner);;
    const Transform      jointFrame;
    const Transform      baseFrame(Vec3(1,2,3));
    const Transform      groundFrame;


    int base = bendStretchBlock.addRigidBody(MassProperties(1,Vec3(0),Inertia::sphere(1)),
                        BodyFrame,
                        Ground, groundFrame,
                        Mobilizer::Pin);
    forces.addMobilityLinearSpring(base, 0, 10., 0);
    //int base = Ground;

    bool useDummy = false ;

    int mass;
    int bendBody, bendCoord;
    int stretchBody, stretchCoord;
    if (useDummy) {
        int dummy = bendStretchBlock.addRigidBody(MassProperties(0,Vec3(0),Inertia(0)),
                                                  Transform(),
                                                  base, baseFrame,
                                                  Mobilizer::Pin);

        bendBody = dummy;
        bendCoord = 0;
        mass = bendStretchBlock.addRigidBody(mProps, jointFrame,
                                dummy, Transform(),
                                Mobilizer::Sliding);
        stretchBody = mass;
        stretchCoord = 0;
    } else {
        mass = bendStretchBlock.addRigidBody(mProps, jointFrame,
                                    base, baseFrame,
                                    Mobilizer::BendStretch);
        bendBody = mass;
        bendCoord = 0;
        stretchBody = mass;
        stretchCoord = 1;
    }

    forces.addMobilityLinearSpring(bendBody, bendCoord, 1000., 0);
    forces.addMobilityLinearSpring(stretchBody, stretchCoord, 100., 1);

    VTKReporter display(mbs);
    display.addDecoration(mass, Transform(),
        DecorativeSphere(.25).setOpacity(.25));

    DecorativeLine crossBodyBond; crossBodyBond.setColor(Orange).setLineThickness(5);
    display.addRubberBandLine(base, baseFrame.T(), mass, jointFrame.T(), crossBodyBond);

    State s;
    mbs.realize(s, Stage::Topology);
    mbs.realize(s, Stage::Model);

    RungeKuttaMerson study(mbs, s);

    bendStretchBlock.setMobilizerQ(s, base, 0, 0);
    bendStretchBlock.setMobilizerU(s, base, 0, .1);
    bendStretchBlock.setMobilizerQ(s, bendBody, bendCoord, 0); // rotation
    bendStretchBlock.setMobilizerU(s, bendBody, bendCoord, 3); // rotation
    bendStretchBlock.setMobilizerQ(s, stretchBody, stretchCoord, 1); // translation
    bendStretchBlock.setMobilizerU(s, stretchBody, stretchCoord, 7);

    mbs.realize(s, Stage::Acceleration);
    const Transform&  x = bendStretchBlock.getBodyPosition(s, mass);
    cout << "Mass x=" << x.T() << " pe=" << mbs.getPotentialEnergy(s) << endl;
    const SpatialVec& v = bendStretchBlock.getBodyVelocity(s, mass);
    cout << "Mass v=" << v << " ke=" << mbs.getKineticEnergy(s) << endl;
    const SpatialVec& acc = bendStretchBlock.getBodyAcceleration(s, mass);
    cout << "Mass a=" << acc << endl;
    const SpatialVec& a=bendStretchBlock.getCoriolisAcceleration(s,mass);
    const SpatialVec& ta=bendStretchBlock.getCoriolisAcceleration(s,mass);
    const SpatialVec& g=bendStretchBlock.getGyroscopicForce(s,mass);
    const SpatialVec& c=bendStretchBlock.getCentrifugalForces(s,mass);
    cout << "coriolis a=" << a << " incl parent=" << ta << endl;
    cout << "gyro f=" << g << endl;
    cout << "centrifugal f=" << c << endl;
    cout << "ydot=" << s.getYDot() << endl;
    const SpatialMat& bai=bendStretchBlock.getArticulatedBodyInertia(s, base);
    cout << "base abi=" << bai;
    const SpatialMat& mai=bendStretchBlock.getArticulatedBodyInertia(s, mass);
    cout << "mass abi=" << mai;


    display.report(s);
    //exit(0);

    const Real h = .01;
    const int interval = 1;
    const Real tstart = 0.;
    const Real tmax = 30; //ps

    study.setAccuracy(1e-8);
    study.initialize(); 

    std::vector<State> saveEm;
    saveEm.push_back(s);
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    display.report(s);

    const Real Estart = mbs.getEnergy(s);

    s.updTime() = tstart;
    int step = 0;
    while (s.getTime() < tmax) {
        study.step(s.getTime() + h);

        cout << s.getTime();
        cout << " deltaE=" << (mbs.getEnergy(s)-Estart)/Estart
             << " (pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << ") hNext=" << study.getPredictedNextStep();
        cout << " bend=" << bendStretchBlock.getMobilizerQ(s, bendBody, bendCoord)/RadiansPerDegree
             << " stretch=" << bendStretchBlock.getMobilizerQ(s, stretchBody, stretchCoord);
        cout << endl;

        if (!(step % interval)) {
            display.report(s);
            saveEm.push_back(s);
        }
        ++step;
    }

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            display.report(saveEm[i]);
            //display.report(saveEm[i]); // half speed
        }
        getchar();
    }

}
catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
}
return 0;
}
