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
 * This is an outer block for simulating ??? in various ways with Simbody.
 * This is about testing Simbody, *not* studying ???!
 */

#include "SimTKsimbody.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
#include <vector>
using std::cout;
using std::endl;

using namespace SimTK;

static const Real Pi      = (Real)SimTK_PI;

int main() {
try
  { MultibodySystem          mbs;
    SimbodyMatterSubsystem   matter(mbs);
    UniformGravitySubsystem  gravity(mbs, Vec3(0,-9.8,0));
    GeneralForceElements     forces(mbs);
    DecorationSubsystem      artwork(mbs);

        // BUILD THE MATTER SUBSYSTEM

    MobilizedBody::Pin pend1(matter.Ground(), 
                             Body::Rigid(MassProperties(1, Vec3(-2,-2,0), Inertia(1))));

    //Constraint::Rod(matter.Ground(), Vec3(1,-1,0), pend1, Vec3(0,-2,0), 1.5);

    forces.addGlobalEnergyDrain(.1);

    Real accuracy = 1e-2;
    Real outputInterval = .10;
    Real simulationLength = 100;

    mbs.realizeTopology();
    State s = mbs.getDefaultState();

    mbs.realize(s, Stage::Model);

    //VTKReporter display(mbs, 0.1);

    RungeKuttaMerson study(mbs, s);
    //CPodesIntegrator study(mbs,s);

    const Real h = outputInterval;
    const int interval = 1;
    const Real tstart = 0.;
    const Real tmax = simulationLength; //ps

    s.updTime() = tstart;
    //display.report(s);

    study.setAccuracy(accuracy);
    study.initialize(); 

    std::vector<State> saveEm;
    saveEm.push_back(s);
    for (int i=0; i<100; ++i)
        saveEm.push_back(s);    // delay
    //display.report(s);

    const Real Estart = mbs.getEnergy(s);

    int step = 0; bool flag=false;
    while (s.getTime() <= tmax) {
        mbs.realize(s);

        cout << s.getTime();
        cout << " deltaE=" << 100*(mbs.getEnergy(s)-Estart)
                                /(std::abs(Estart)+NTraits<Real>::Tiny) 
             << "% pe=" << mbs.getPotentialEnergy(s)
             << ", ke=" << mbs.getKineticEnergy(s)
             << " hNext(fs)=" << 1000*study.getPredictedNextStep();

        cout << "\n  System COM loc=" << matter.calcSystemMassCenterLocationInGround(s);
        cout << "\n  System COM vel=" << matter.calcSystemMassCenterVelocityInGround(s);
        cout << "\n  System COM acc=" << matter.calcSystemMassCenterAccelerationInGround(s);
        cout << endl;

        cout << "     q=" << matter.getQ(s) << endl;
        cout << "     u=" << matter.getU(s) << endl;
        cout << "  udot=" << matter.getUDot(s) << endl;

        cout << endl;

        if (!(step % interval)) {
            //display.report(s);
            saveEm.push_back(s);
        }

        study.step(s.getTime() + h);
        ++step;
    }

    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            //display.report(saveEm[i]);
            //display.report(saveEm[i]); // half speed
        }
        getchar();
    }

  }
catch (const std::exception& e)
  {
    printf("EXCEPTION THROWN: %s\n", e.what());
  }
catch (...)
  {
    printf("UNKNOWN EXCEPTION THROWN\n");
  }    return 0;
}


