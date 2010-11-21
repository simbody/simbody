#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h"
#include "simbody/internal/VisualizationReporter.h"
#include "simbody/internal/VisualizationEventListener.h"
//#include <sys/time.h>

using namespace SimTK;

class MyListener : public VisualizationEventListener {
public:
    MyListener(const Array_< std::pair<std::string, int> >& menu)
    :   m_menu(menu) {}

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
            std::cout << "Special key hit: " << mod << " key=" << key << " glut=" << (key & ~SpecialKeyOffset);
        else
            std::cout << "Ordinary key hit: " << mod << char(key) << " (" << (int)key << ")";
        std::cout << " " << nm << std::endl;


        return true; // key absorbed
    }

    virtual bool menuSelected(int item) {
        std::cout << "Menu item " << item << ": ";
        for (unsigned i=0; i < m_menu.size(); ++i)
            if (m_menu[i].second==item)
                std::cout << m_menu[i].first;
        std::cout << std::endl;
        return true; // menu click absorbed
    }

private:
    Array_< std::pair<std::string, int> > m_menu;
};

int main() {
  try {
    // Create the system.

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, -1*Vec3(0, -9.8, 0));
    Body::Rigid pendulumBody[2]; // solid, translucent
    pendulumBody[0].setDefaultRigidBodyMassProperties(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody[0].addDecoration(Transform(), DecorativeSphere(0.49).setOpacity(1));
    pendulumBody[1].setDefaultRigidBodyMassProperties(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody[1].addDecoration(Transform(), DecorativeSphere(0.49).setOpacity(.5));
    MobilizedBodyIndex lastBody = matter.getGround().getMobilizedBodyIndex();
    for (int i = 0; i < 50; ++i) {
        MobilizedBody::Ball pendulum(matter.updMobilizedBody(lastBody), Transform(Vec3(0)), 
            pendulumBody[i%2], Transform(Vec3(0, 1, 0))); // alternate solid, translucent
        lastBody = pendulum.getMobilizedBodyIndex();
    }
//    system.updDefaultSubsystem().addEventReporter(new VTKEventReporter(system, 0.02));
    VisualizationReporter* vr = new VisualizationReporter(system, 0.02);
    system.updDefaultSubsystem().addEventReporter(vr);

        Array_< std::pair<std::string,int> > items;
    items.push_back(std::make_pair("One", 1));
    items.push_back(std::make_pair("Top/SubA/first", 2));
    items.push_back(std::make_pair("Top/SubA/second", 3));
    items.push_back(std::make_pair("Top/SubB/only", 4));
    items.push_back(std::make_pair("Two", 5));
    vr->updVisualizer().addMenu("Test Menu",items);

    vr->updVisualizer().addEventListener(new MyListener(items));

    //vr->updVisualizer().setMode(Visualizer::RealTime);



    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;
    for (int i = 0; i < state.getNQ(); ++i)
        state.updQ()[i] = random.getValue(); 

    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-4);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    //timeval start;
    //gettimeofday(&start, NULL);
    ts.stepTo(100.0);
    //timeval end;
    //gettimeofday(&end, NULL);
    //std::cout << "time: "<<(end.tv_sec-start.tv_sec)<< std::endl;
//    std::cout << "steps: "<<integ.getNStepsTaken()<< std::endl;


    std::cout << "Type something to quit: ";
    char ch; std::cin >> ch;

  } catch (const std::exception& exc) {
      std::cout << "EXCEPTION: " << exc.what() << std::endl;
  }
}
//int main() {
//
//    // Create the system.
//
//    MultibodySystem system;
//    SimbodyMatterSubsystem matter(system);
//    GeneralForceSubsystem forces(system);
//    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
//    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
//    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));
//    MobilizedBodyIndex lastBody = matter.getGround().getMobilizedBodyIndex();
//    for (int i = 0; i < 10; ++i) {
//        MobilizedBody::Ball pendulum(matter.updMobilizedBody(lastBody), Transform(Vec3(0)), pendulumBody, Transform(Vec3(0, 1, 0)));
//        lastBody = pendulum.getMobilizedBodyIndex();
//    }
//    system.updDefaultSubsystem().addEventReporter(new VTKEventReporter(system, 0.02));
//
//    // Initialize the system and state.
//
//    system.realizeTopology();
//    State state = system.getDefaultState();
//    Random::Gaussian random;
//    for (int i = 0; i < state.getNQ(); ++i)
//        state.updQ()[i] = random.getValue();
//
//    // Simulate it.
//
//    RungeKuttaMersonIntegrator integ(system);
//    TimeStepper ts(system, integ);
//    ts.initialize(state);
//    ts.stepTo(100.0);
//}

