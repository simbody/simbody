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
gravity, uses a FrameController to track a body with the camera and show some
feedback to the user, and adds a menu to the display. A description of all user 
input received is written to the console, and some inputs are used to control 
the simulation. */

#include "SimTKsimbody.h"

using namespace SimTK;

static const int GravityX=1, GravityY=2, GravityZ=3, GravityMag=4; // sliders
static const Real GravityDefault=9.81, GravityMax=20;

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
    MyListener(const Array_< std::pair<std::string, int> >& menu)
    :   m_menu(menu) {}

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

    virtual bool menuSelected(int item) {
        std::cout << "Listener sees pick of menu item " << item << ": ";
        for (unsigned i=0; i < m_menu.size(); ++i)
            if (m_menu[i].second==item)
                std::cout << m_menu[i].first;
        std::cout << std::endl;
        return false; // menu click passed on
    }

    virtual bool sliderMoved(int whichSlider, Real value) {
        printf("Listener sees slider %d now at %g\n", whichSlider, value);
        return false;   // slider move passed on
    }

private:
    Array_< std::pair<std::string, int> > m_menu;
};

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input. If there has been some, process it.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, 
                     const Force::Gravity& gravity, 
                     Real interval) 
    :   PeriodicEventHandler(interval), m_silo(silo), m_gravity(gravity) {}

    virtual void handleEvent(State& state, Real accuracy, const Vector& yWeights, 
                             const Vector& ooConstraintTols, Stage& lowestModified, 
                             bool& shouldTerminate) const 
    {
        while (m_silo.isAnyUserInput()) {
            unsigned key, modifiers;
            int menuItem;
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

            while (m_silo.takeMenuPick(menuItem)) {
                printf("Handler sees menu pick %d\n", menuItem);
            }

            while (m_silo.takeSliderMove(whichSlider, sliderValue)) {
                switch(whichSlider) {
                case GravityMag:
                    m_gravity.setMagnitude(state, sliderValue);
                    break;
                case GravityX: case GravityY: case GravityZ: {
                    Vec3 gdir = m_gravity.getDownDirection(state);
                    gdir[whichSlider-GravityX] = sliderValue;
                    if (gdir.norm() < SignificantReal) gdir=Vec3(0,1,0);
                    m_gravity.setDownDirection(state, gdir);
                    break;
                }
                };
            }
        }  
    }

private:
    Visualizer::InputSilo& m_silo;
    const Force::Gravity&  m_gravity;
};

int main() {
  try {
    // Create the system.
    const int  NBodies = 50;
    const Real FrameRate = 30;
    const Real TimeScale = 1;


    printf("\n\n************\n");
    printf(     "ESC to quit\n");
    printf(     "************\n\n");


    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, UnitVec3(YAxis), GravityDefault);
    Force::GlobalDamper(forces, matter, 3);
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

    VisualizationReporter* vr = 
        new VisualizationReporter(system, TimeScale/FrameRate);
    system.updDefaultSubsystem().addEventReporter(vr);
    Visualizer& viz = vr->updVisualizer();

    // Add a menu, just for fun.
    Array_< std::pair<std::string,int> > items;
    items.push_back(std::make_pair("One", 1));
    items.push_back(std::make_pair("Top/SubA/first", 2));
    items.push_back(std::make_pair("Top/SubA/second", 3));
    items.push_back(std::make_pair("Top/SubB/only", 4));
    items.push_back(std::make_pair("Two", 5));
    viz.addMenu("Test Menu",items);

    // Add sliders to control gravity. They will display from bottom up.
    viz.addSlider("Gravity Z", 3, -1, 1, 0); 
    viz.addSlider("Gravity Y", 2, -1, 1, 1);
    viz.addSlider("Gravity X", 1, -1, 1, 0);
    viz.addSlider("Gravity Mag", 4, 0, GravityMax, GravityDefault);

    MyListener*            listener = new MyListener(items);
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
       (new UserInputHandler(*silo, gravity, 0.1)); // check input every 100ms
     
    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;
    for (int i = 0; i < state.getNQ(); ++i)
        state.updQ()[i] = random.getValue(); 

    Assembler(system).assemble(state);

    // Simulate it.

    //RungeKutta3Integrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    //RungeKuttaFeldbergIntegrator integ(system);
    //CPodesIntegrator integ(system);
    integ.setAccuracy(1e-3);
    TimeStepper ts(system, integ);
    ts.initialize(state);

    double cpuStart = cpuTime();
    double realStart = realTime();
    ts.stepTo(Infinity); // user must hit ESC to stop sim
    std::cout << "cpu time:  "<<cpuTime()-cpuStart<< std::endl;
    std::cout << "real time: "<<realTime()-realStart<< std::endl;
    std::cout << "steps:     "<<integ.getNumStepsTaken()<< std::endl;
    vr->getVisualizer().dumpStats(std::cout);

    std::cout << "Type something to quit: ";
    char ch; std::cin >> ch;

  } catch (const std::exception& exc) {
      std::cout << "EXCEPTION: " << exc.what() << std::endl;
  }
}
