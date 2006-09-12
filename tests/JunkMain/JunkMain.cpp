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
#include "Simmatrix.h"
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
    //                       --------------
    //     <--             /              /|
    //        \            --------------  /
    //      *-/-------->* |              |/
    //                     -------------- 
    //
    //      | y
    //      |
    //      |-----> x
    //     /
    //    z
    // This is a rectangle block mounted on a "bend stretch" joint which
    // permits rotation around mutual z axes of block & ground, followed
    // by a translation along the *block* x axis (which will have moved).
    //
    // Block body frame is at its center. Half dimensions are
    // (hx,hy,hz)=3,1,2. Joint frame is at -hx,-hy,-hz
    //

    const Real hx=3, hy=1, hz=2;

    const Real mass = 1;
    const Vec3 com = Vec3(0);
    const InertiaMat iner(1,1,1)/* = mass*InertiaMat::brick(hx,hy,hz)*/;

    const MassProperties blockProps(mass,com,iner);
    RotationMat jf; jf.setToBodyFixed123(Vec3(0,0,/*Pi/2*/0));
    const Transform      jointFrame(jf,/*Vec3(-hx,-hy,-hz)*/Vec3(0));
    const Transform      baseFrame(jf,Vec3(0,0,0));
    const Transform      groundFrame;


    int base = bendStretchBlock.addRigidBody(MassProperties(1,Vec3(0),InertiaMat::sphere(1)),
                        baseFrame,
                        Ground, groundFrame,
                        Mobilizer::Pin);
   // forces.addMobilityLinearSpring(base, 0, 10., 0);

    bool useDummy = true;

    int block;
    int bendBody, bendCoord;
    int stretchBody, stretchCoord;
    if (useDummy) {
        int dummy = bendStretchBlock.addRigidBody(MassProperties(0,Vec3(0),InertiaMat(0)),
                                                  Transform(),
                                                  base, baseFrame,
                                                  Mobilizer::Pin);

        bendBody = dummy;
        bendCoord = 0;
        block = bendStretchBlock.addRigidBody(blockProps, jointFrame,
                                dummy, Transform(),
                                Mobilizer::Sliding);
        stretchBody = block;
        stretchCoord = 0;
    } else {
        block = bendStretchBlock.addRigidBody(blockProps, jointFrame,
                                    base, baseFrame,
                                    Mobilizer::BendStretch);
        bendBody = block;
        bendCoord = 0;
        stretchBody = block;
        stretchCoord = 1;
    }

    forces.addMobilityLinearSpring(bendBody, bendCoord, 1000., 0);
    forces.addMobilityLinearSpring(stretchBody, stretchCoord, 100., 1);

    VTKReporter display(mbs);
    //display.addDecoration(block, Transform(),
        //DecorativeBrick(Vec3(hx,hy,hz)).setOpacity(.25));
    display.addDecoration(block, Transform(),
        DecorativeSphere(.5).setOpacity(.25));

    DecorativeLine crossBodyBond; crossBodyBond.setColor(Orange).setLineThickness(5);
    display.addRubberBandLine(base, baseFrame.T(), block, jointFrame.T(), crossBodyBond);

    State s;
    mbs.realize(s, Stage::Built);
    mbs.realize(s, Stage::Modeled);

    RungeKuttaMerson study(mbs, s);

    bendStretchBlock.setMobilizerQ(s, base, 0, 0);
    bendStretchBlock.setMobilizerU(s, base, 0, 10);
    bendStretchBlock.setMobilizerQ(s, bendBody, bendCoord, 0); // rotation
    bendStretchBlock.setMobilizerQ(s, stretchBody, stretchCoord, 1); // translation

    display.report(s);

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
        cout << " bend=" << bendStretchBlock.getMobilizerQ(s, block, 0)/RadiansPerDegree
             << " stretch=" << bendStretchBlock.getMobilizerQ(s, block, 1);
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
