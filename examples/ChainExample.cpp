#include "SimTKsimbody.h"
#include "pthread.h"

#include <deque>

using namespace SimTK;

class MyFrameController : public Visualizer::FrameController {
public:
    MyFrameController(const SimbodyMatterSubsystem& matter,
                      MobilizedBodyIndex whichBody,
                      const Force::Gravity& gravity) 
    :   m_matter(matter), m_whichBody(whichBody), m_gravity(gravity) {}

    virtual void generateControls(const Visualizer&           viz, 
                                  const State&                state,
                                  Array_<DecorativeGeometry>& geometry)
    {
        const MobilizedBody& mobod = m_matter.getMobilizedBody(m_whichBody);
        const Transform& X_GB = mobod.getBodyTransform(state);
        const UnitVec3& downDir = m_gravity.getDownDirection(state);
        Vec3 cameraPos(X_GB.p()[0], X_GB.p()[1], 40);
        UnitVec3 cameraZ(0,0,1);
        viz.setCameraTransform(Transform(Rotation(cameraZ, ZAxis, Vec3(0,1,0), YAxis), cameraPos));

        geometry.push_back(DecorativeLine(Vec3(0), 10*downDir)
            .setColor(Green).setLineThickness(4).setBodyId(0));
    }

private:
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_whichBody;
    const Force::Gravity&           m_gravity;
};

class MyListener : public Visualizer::EventListener {
public:
    MyListener(const Array_< std::pair<std::string, int> >& menu)
    :   m_menu(menu) {pthread_mutex_init(&charQueueLock,0);}

    ~MyListener() {pthread_mutex_destroy(&charQueueLock);}

    // Called from the main thread.
    std::pair<unsigned,unsigned> takeNextCharIfAny() {
        std::pair<unsigned,unsigned> theChar(0,0);
        if (charQueue.size()) {
            pthread_mutex_lock(&charQueueLock);
            if (charQueue.size()) {
                theChar = charQueue.front();
                charQueue.pop_front();
            }
            pthread_mutex_unlock(&charQueueLock);
        }
        return theChar;
    }

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

        pthread_mutex_lock(&charQueueLock);
        charQueue.push_back(std::make_pair(key,modifier));
        pthread_mutex_unlock(&charQueueLock);


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
    pthread_mutex_t charQueueLock;
    std::deque<std::pair<unsigned,unsigned> > charQueue;
};

// Check for user input. If there has been some, process it.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(MyListener& listener, const Force::Gravity& gravity, Real interval) 
    :   PeriodicEventHandler(interval), m_listener(listener), m_gravity(gravity) {}

    virtual void handleEvent(State& state, Real accuracy, const Vector& yWeights, 
                             const Vector& ooConstraintTols, Stage& lowestModified, 
                             bool& shouldTerminate) const 
    {
        std::pair<unsigned,unsigned> charHit = m_listener.takeNextCharIfAny();
        if (charHit.first == 0) return;

        if (charHit.first == Visualizer::EventListener::KeyEsc) {
            printf("User hit ESC!!\n");
            shouldTerminate = true;
        } else {
            const UnitVec3& down = m_gravity.getDownDirection(state);
            bool control = (charHit.second & Visualizer::EventListener::ControlIsDown) != 0;
            switch(charHit.first) {
            case Visualizer::EventListener::KeyLeftArrow:
                    m_gravity.setDownDirection(state, down + .1*Vec3(-1,0,0));
                    lowestModified = Stage::Instance;
                    break;
            case Visualizer::EventListener::KeyRightArrow:
                    m_gravity.setDownDirection(state, down + .1*Vec3(1,0,0));
                    lowestModified = Stage::Instance;
                    break;
            case Visualizer::EventListener::KeyUpArrow:
                    m_gravity.setDownDirection(state, down + .1*Vec3(0,1,0));
                    lowestModified = Stage::Instance;
                    break;
            case Visualizer::EventListener::KeyDownArrow:
                    m_gravity.setDownDirection(state,  down + .1*Vec3(0,-1,0));
                    lowestModified = Stage::Instance;
                    break;
            case Visualizer::EventListener::KeyPageUp:
                    m_gravity.setDownDirection(state, down + .1*Vec3(0,0,-1));
                    lowestModified = Stage::Instance;
                    break;
            case Visualizer::EventListener::KeyPageDown:
                    m_gravity.setDownDirection(state, down + .1*Vec3(0,0,1));
                    lowestModified = Stage::Instance;
                    break;
            }
            if (lowestModified = Stage::Instance)
                std::cout << "New gravity down=" << m_gravity.getDownDirection(state) << std::endl;
        }
    }

private:
    MyListener& m_listener;
    const Force::Gravity& m_gravity;
};

int main() {
  try {
    // Create the system.

    const Real FrameRate = 30;
    const Real TimeScale = 1;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, UnitVec3(YAxis), 9.8);
    Force::GlobalDamper(forces, matter, 3);
    Body::Rigid pendulumBody[2]; // solid, translucent
    pendulumBody[0].setDefaultRigidBodyMassProperties(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody[0].addDecoration(Transform(), DecorativeSphere(0.49).setOpacity(1));
    pendulumBody[1].setDefaultRigidBodyMassProperties(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody[1].addDecoration(Transform(), DecorativeSphere(0.49).setOpacity(.5));
    MobilizedBody lastBody = matter.Ground();
    for (int i = 0; i < 50; ++i) {
        MobilizedBody::Ball pendulum(lastBody, Transform(Vec3(0)), 
            pendulumBody[i%2], Transform(Vec3(0, 1, 0))); // alternate solid, translucent
        lastBody = pendulum;
    }

    Constraint::Ball(matter.Ground(), Vec3(30,0,0), lastBody, Vec3(0));



    VisualizationReporter* vr = new VisualizationReporter(system, TimeScale/FrameRate);
    system.updDefaultSubsystem().addEventReporter(vr);
    Visualizer& viz = vr->updVisualizer();
   
    printf("\n\n***************************************************************\n");
    printf(    "use arrow keys and page up/down to control green gravity vector\n");
    printf(    "***************************************************************\n\n");

        Array_< std::pair<std::string,int> > items;
    items.push_back(std::make_pair("One", 1));
    items.push_back(std::make_pair("Top/SubA/first", 2));
    items.push_back(std::make_pair("Top/SubA/second", 3));
    items.push_back(std::make_pair("Top/SubB/only", 4));
    items.push_back(std::make_pair("Two", 5));
    viz.addMenu("Test Menu",items);

    MyListener& listener = *new MyListener(items);
    viz.addEventListener(&listener);

    viz.addFrameController(new MyFrameController(matter, MobilizedBodyIndex(25), gravity));

    viz.setRealTimeScale(TimeScale);
    //viz.setDesiredBufferLengthInSec(.15);
    viz.setDesiredFrameRate(FrameRate);
    //viz.setMode(Visualizer::Sampling);
    viz.setMode(Visualizer::RealTime);

    system.updDefaultSubsystem().addEventHandler(new UserInputHandler(listener, gravity, 0.1)); 
     

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
    ts.stepTo(100.0);
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
