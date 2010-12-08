/* -------------------------------------------------------------------------- *
 *                                Simbody(tm)                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
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

/*                          Simbody ChainExample
This example demonstrates how to use the Simbody Visualizer to display and
interact with a real time simulation. It shows the use of sliders to control
"wind", uses a FrameController to track a body with the camera and show some
feedback to the user, and adds a menu to the display. A description of all user 
input received is written to the console, and some inputs are used to control 
the simulation. */

#include "SimTKsimbody.h"

#include <cstdio>
#include <iostream>

using namespace SimTK;

static const int NBodies = 100;
const Real FrameRate = 30;
const Real TimeScale = 1; // i.e., 2 -> 2X real time

// We call this "wind" but it is implemented with Force::Gravity.
static const int GravityX=1, GravityY=2, GravityZ=3, GravityMag=4; // sliders
static const Real GravityDefault=10, GravityMax=20;

// A FrameController is called by the Visualizer just prior to rendering a
// frame. Here we'll point the camera and add some geometry showing the direction
// and magnitude of gravity (which the user can change via sliders).
class MyFrameController : public Visualizer::FrameController {
public:
    MyFrameController(const SimbodyMatterSubsystem& matter,
                      MobilizedBodyIndex whichBody, // tracked with camera
                      const Force::Gravity& gravity) 
    :   m_matter(matter), m_whichBody(whichBody), m_gravity(gravity) {}

    virtual void generateControls(const Visualizer&           viz, 
                                  const State&                state,
                                  Array_<DecorativeGeometry>& geometry)
    {
        const MobilizedBody& mobod = m_matter.getMobilizedBody(m_whichBody);
        const Transform& X_GB = mobod.getBodyTransform(state);
        const UnitVec3& downDir = m_gravity.getDownDirection(state);
        const Real      gmag    = m_gravity.getMagnitude(state);

        // Point the camera at the chosen body.
        viz.pointCameraAt(X_GB.p(), Vec3(0,1,0));

        // Show gravity as a fat green line.
        geometry.push_back(DecorativeLine(Vec3(0), gmag*downDir)
            .setColor(Green).setLineThickness(4).setBodyId(0));
    }

private:
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_whichBody;
    const Force::Gravity&           m_gravity;
};

// This is a custom InputListener. We'll register it prior to the InputSilo so
// that we can intercept all input and say something about it. No input processing
// is done here other than that, and we pass on everything we receive down the
// chain to the next listener (which will be an InputSilo in this case).
class MyListener : public Visualizer::InputListener {
public:
    // Pass in the menu strings just so we can translate the index back
    // to a string to print out for testing.
    MyListener(const Array_< std::pair<std::string, int> >& menu1,
               const Array_< std::pair<std::string, int> >& menu2)
    :   m_menu1(menu1), m_menu2(menu2) {}

    ~MyListener() {}

    virtual bool keyPressed(unsigned key, unsigned modifier) {
        String mod;
        if (modifier&ControlIsDown) mod += "CTRL ";
        if (modifier&ShiftIsDown) mod += "SHIFT ";
        if (modifier&AltIsDown)  mod += "ALT ";

        const char* nm = "NoNickname";
        switch(key) {
        case KeyEsc: nm="ESC"; break;
        case KeyDelete: nm="DEL"; break;
        case KeyRightArrow: nm="Right"; break;
        case KeyLeftArrow: nm="Left"; break;
        case KeyUpArrow: nm="Up"; break;
        case KeyDownArrow: nm="Down"; break;
        case KeyEnter: nm="ENTER"; break;
        case KeyF1: nm="F1"; break;
        case KeyF12: nm="F12"; break;
        case 'a': nm="lower a"; break;
        case 'Z': nm="upper Z"; break;
        case '}': nm="right brace"; break;
        }
        if (modifier&IsSpecialKey)
            std::cout << "Listener saw special key hit: " 
                << mod << " key=" << key << " glut=" << (key & ~SpecialKeyOffset);
        else
            std::cout << "Listener saw ordinary key hit: " 
                << mod << char(key) << " (" << (int)key << ")";
        std::cout << " " << nm << std::endl;

        return false; // key passed on
    }

    virtual bool menuSelected(int menuId, int item) {
        std::cout << "Listener sees pick of menu " << menuId << " item " << item << ": ";
        Array_< std::pair<std::string, int> >& menu = menuId==1 ? m_menu1 : m_menu2;
        for (unsigned i=0; i < menu.size(); ++i)
            if (menu[i].second==item)
                std::cout << menu[i].first;
        std::cout << std::endl;
        return false; // menu click passed on
    }

    virtual bool sliderMoved(int whichSlider, Real value) {
        printf("Listener sees slider %d now at %g\n", whichSlider, value);
        return false;   // slider move passed on
    }

private:
    Array_< std::pair<std::string, int> > m_menu1, m_menu2;
};

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input. If there has been some, process it.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer& viz,
                     Visualizer::InputSilo& silo, 
                     const Force::Gravity& gravity, 
                     Real interval) 
    :   PeriodicEventHandler(interval), m_viz(viz), m_silo(silo), m_gravity(gravity) {}

    virtual void handleEvent(State& state, Real accuracy, const Vector& yWeights, 
                             const Vector& ooConstraintTols, Stage& lowestModified, 
                             bool& shouldTerminate) const 
    {
        while (m_silo.isAnyUserInput()) {
            unsigned key, modifiers;
            int whichMenu, menuItem;
            int whichSlider; Real sliderValue;

            while (m_silo.takeKeyHit(key,modifiers)) {
                if (key == Visualizer::InputListener::KeyEsc) {
                    printf("User hit ESC!!\n");
                    shouldTerminate = true;
                    m_silo.clear();
                    return;
                }
                printf("Handler sees key=%u, modifiers=%u\n",key,modifiers);
            }

            while (m_silo.takeMenuPick(whichMenu, menuItem)) {
                printf("Handler sees menu %d, pick %d\n", whichMenu, menuItem);
            }

            while (m_silo.takeSliderMove(whichSlider, sliderValue)) {
                if (whichSlider == GravityMag) {
                    m_gravity.setMagnitude(state, sliderValue);
                    continue;
                }
                Vec3 gdir = m_gravity.getDownDirection(state);
                Real remaining = std::sqrt(std::max(0., 1-square(sliderValue)));
                CoordinateAxis axis = CoordinateAxis::getCoordinateAxis(whichSlider-GravityX);
                CoordinateAxis prev = axis.getPreviousAxis();
                CoordinateAxis next = axis.getNextAxis();
                Vec2 other(gdir[prev], gdir[next]);
                if (other.norm() >= SignificantReal) other *= (remaining/other.norm());
                gdir[axis]=sliderValue; gdir[prev]=other[0]; gdir[next]=other[1];
                if (gdir.norm() < SignificantReal) gdir[next] = 1;
                m_viz.setSliderValue(GravityX+prev, gdir[prev]);
                m_viz.setSliderValue(GravityX+next, gdir[next]);
                m_gravity.setDownDirection(state, gdir);
            }
        }  
    }

private:
    Visualizer&            m_viz;
    Visualizer::InputSilo& m_silo;
    const Force::Gravity&  m_gravity;
};

int main() {
  try {
    // Create the system.

    printf("\n\n************\n");
    printf(     "ESC to quit\n");
    printf(     "************\n\n");


    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, UnitVec3(YAxis), GravityDefault);
    Force::GlobalDamper(forces, matter, 7);
    Body::Rigid pendulumBody[2]; // solid, translucent
    pendulumBody[0].setDefaultRigidBodyMassProperties(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody[0].addDecoration(Transform(), DecorativeSphere(0.49).setOpacity(1));
    pendulumBody[1].setDefaultRigidBodyMassProperties(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody[1].addDecoration(Transform(), DecorativeSphere(0.49).setOpacity(.5));
    MobilizedBody lastBody = matter.Ground();
    for (int i = 0; i < NBodies; ++i) {
        MobilizedBody::Ball pendulum(lastBody, Transform(Vec3(0)), 
            pendulumBody[i%2], Transform(Vec3(0, 1, 0))); // alternate solid, translucent
        lastBody = pendulum;
    }

    // Attach the last body back to ground.
    Constraint::Ball(matter.Ground(), Vec3(NBodies/2,0,0), lastBody, Vec3(0));

    Visualizer viz(system);

    viz.setWindowTitle("This is the so-called 'ChainExample'.");

    // Add a menu, just for fun.
    Array_< std::pair<std::string,int> > menu1, menu2;
    menu1.push_back(std::make_pair("One", 1));
    menu1.push_back(std::make_pair("Top/SubA/first", 2));
    menu1.push_back(std::make_pair("Top/SubA/second", 3));
    menu1.push_back(std::make_pair("Top/SubB/only", 4));
    menu1.push_back(std::make_pair("Two", 5));
    viz.addMenu("Test Menu", 1, menu1);

    // And another one, to check the id handling.
    menu2.push_back(std::make_pair("One", 1));
    menu2.push_back(std::make_pair("Two", 2));
    viz.addMenu("More", 2, menu2);

    MyListener*            listener = new MyListener(menu1,menu2);
    Visualizer::InputSilo* silo = new Visualizer::InputSilo();
    viz.addInputListener(listener); // order matters here
    viz.addInputListener(silo);

    // Tell the frame controller to track the middle body.
    viz.addFrameController(new MyFrameController(matter, 
        MobilizedBodyIndex(NBodies/2), gravity));

    viz.setRealTimeScale(TimeScale);
    //viz.setDesiredBufferLengthInSec(.15);
    viz.setDesiredFrameRate(FrameRate);
    //viz.setMode(Visualizer::Sampling);
    //viz.setMode(Visualizer::PassThrough);
    viz.setMode(Visualizer::RealTime);

    viz.setCameraTransform(Vec3(0,NBodies/4,2*NBodies)); 

    system.updDefaultSubsystem().addEventHandler
       (new UserInputHandler(viz,*silo, gravity, 0.1)); // check input every 100ms

    // Report visualization frames.
    Visualizer::Reporter* vr = new Visualizer::Reporter(viz, TimeScale/FrameRate);
    system.updDefaultSubsystem().addEventReporter(vr);
    
    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;
    for (int i = 0; i < state.getNQ(); ++i)
        state.updQ()[i] = random.getValue(); 

    // Use the Assembler to satisfy the loop-closing constraint.
    Assembler assembler(system);
    std::cout << "ASSEMBLING ... start configuration shown\n";
    viz.report(state);
    std::cout << "  Type something to continue:\n"; getchar();
    double asmRTstart=realTime(), asmCPUstart=cpuTime();
    assembler.addReporter(*vr);
    assembler.setSystemConstraintsWeight(1);
    Visualizer::Mode oldMode = viz.getMode();
    viz.setMode(Visualizer::PassThrough);
    assembler.assemble(state);
    viz.setMode(oldMode);
    printf("...ASSEMBLED in %gs, cpu=%gs. Final configuration shown\n",
        realTime()-asmRTstart, cpuTime()-asmCPUstart);
    viz.report(state);
    std::cout << "  Type something to continue:\n"; getchar();

    // Simulate it.

    // Add sliders to control gravity. They will display from bottom up.
    // Joy Ku thought calling this "wind direction" makes more sense.
    viz.addSlider("Wind Z", 3, -1, 1, 0); 
    viz.addSlider("Wind Y", 2, -1, 1, 1);
    viz.addSlider("Wind X", 1, -1, 1, 0);
    viz.addSlider("Wind Mag", 4, 0, GravityMax, GravityDefault);

    //RungeKutta3Integrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    //RungeKuttaFeldbergIntegrator integ(system);
    //CPodesIntegrator integ(system);
    integ.setAccuracy(1e-2);
    TimeStepper ts(system, integ);
    ts.initialize(state);

    double cpuStart = cpuTime();
    double realStart = realTime();
    //ts.stepTo(10);
    ts.stepTo(Infinity); // user must hit ESC to stop sim
    std::cout << "cpu time:  "<<cpuTime()-cpuStart<< std::endl;
    std::cout << "real time: "<<realTime()-realStart<< std::endl;
    std::cout << "steps:     "<<integ.getNumStepsTaken()<< std::endl;
    vr->getVisualizer().dumpStats(std::cout);

    std::cout << "Type something to quit: "; getchar();

  } catch (const std::exception& exc) {
      std::cout << "EXCEPTION: " << exc.what() << std::endl;
  }
}
