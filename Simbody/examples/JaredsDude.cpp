// -----------------------------------------------------------------------------
// An attempt at duplicating Jared Duke's simulation for studying performance.
// -----------------------------------------------------------------------------
#include "SimTKsimbody.h"

#include <utility>
#include <map>

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
        Vec3 cameraPos(X_GB.p()[0], X_GB.p()[1], 10);
        UnitVec3 cameraZ(0,0,1);
        //viz.setCameraTransform(Transform(Rotation(cameraZ, ZAxis, Vec3(0,1,0), YAxis), cameraPos));

        geometry.push_back(DecorativeLine(Vec3(1,4,0), Vec3(1,4,0)+downDir)
            .setColor(Green).setLineThickness(3).setBodyId(0));
    }

private:
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_whichBody;
    const Force::Gravity&           m_gravity;
};

// Check for user input. If there has been some, process it.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, const Force::Gravity& gravity, Real interval) 
    :   PeriodicEventHandler(interval), m_silo(silo), m_gravity(gravity) {}

    virtual void handleEvent(State& state, Real accuracy, const Vector& yWeights, 
                             const Vector& ooConstraintTols, Stage& lowestModified, 
                             bool& shouldTerminate) const 
    {
        unsigned key, modifiers;
        if (!m_silo.takeKeyHit(key, modifiers))
            return;

        if (key == Visualizer::InputListener::KeyEsc) {
            printf("User hit ESC!!\n");
            shouldTerminate = true;
            return;
        }

        const Real KeyFactor = 0.05; 

        const UnitVec3& down = m_gravity.getDownDirection(state);
        bool control = (modifiers & Visualizer::InputListener::ControlIsDown) != 0;
        switch(key) {
        case Visualizer::InputListener::KeyLeftArrow:
                m_gravity.setDownDirection(state, down + KeyFactor*Vec3(-1,0,0));
                lowestModified = Stage::Instance;
                break;
        case Visualizer::InputListener::KeyRightArrow:
                m_gravity.setDownDirection(state, down + KeyFactor*Vec3(1,0,0));
                lowestModified = Stage::Instance;
                break;
        case Visualizer::InputListener::KeyUpArrow:
                m_gravity.setDownDirection(state, down + KeyFactor*Vec3(0,1,0));
                lowestModified = Stage::Instance;
                break;
        case Visualizer::InputListener::KeyDownArrow:
                m_gravity.setDownDirection(state,  down + KeyFactor*Vec3(0,-1,0));
                lowestModified = Stage::Instance;
                break;
        case Visualizer::InputListener::KeyPageUp:
                m_gravity.setDownDirection(state, down + KeyFactor*Vec3(0,0,-1));
                lowestModified = Stage::Instance;
                break;
        case Visualizer::InputListener::KeyPageDown:
                m_gravity.setDownDirection(state, down + KeyFactor*Vec3(0,0,1));
                lowestModified = Stage::Instance;
                break;
        }
        if (lowestModified = Stage::Instance)
            std::cout << "New gravity down=" << m_gravity.getDownDirection(state) << std::endl;

    }

private:
    Visualizer::InputSilo& m_silo;
    const Force::Gravity&  m_gravity;
};

//////////////////////////////////////////////////////////////////////////

class TwoPointMuscleDamperReflex : public Force::Custom::Implementation {
public:

  TwoPointMuscleDamperReflex(const MobilizedBody& body1,
                             const Vec3& station1,
                             const MobilizedBody& body2,
                             const Vec3& station2,
                             Real k,
                             Real d,
                             Real x0);

  virtual bool dependsOnlyOnPositions() const {
    return false;
  }

  virtual void calcForce(const State& state,
                         Vector_<SpatialVec>& bodyForces,
                         Vector_<Vec3>& particleForces,
                         Vector& mobilityForces) const;

  virtual Real calcPotentialEnergy(const State& state) const;

  void addDecorativeLine( DecorationSubsystem& viz,
                          const DecorativeLine& line,
                          Real scale1=1.0, Real scale2=1.0) const {
    viz.addRubberBandLine(mBody1, mStation1*scale1, mBody2, mStation2*scale2, line);
  }

  Vec3 getStation1( ) const { return mStation1; }
  Vec3 getStation2( ) const { return mStation2; }
  void setStation1(const Vec3& station1) { mStation1 = station1; }
  void setStation2(const Vec3& station2) { mStation2 = station2; }

  Real getK() const { return mK; }
  void setK(Real k) { mK = k; }
  Real getDamping() const { return mDamping; }
  void setDamping(Real damping) { mDamping = damping; }
  Real getX0() const { return mX0; }
  void setX0(Real x0) { mX0 = x0; }

 private:
  const MobilizedBody& mBody1;
  const MobilizedBody& mBody2;
  Vec3 mStation1, mStation2;
  Real mK, mDamping, mX0;
};

class Dude {
public:
    enum Side     {Left=0,Right,Only};
    enum BodyType {Foot=0,Shank,Thigh,Pelvis,Torso};
    enum Segment  {FootFront=0, FootToes, FootHeel,
                   ShankLower, ShankMidLow, ShankMidUp, ShankUpper,
                   ThighLower, ThighWhole,
                   PelvisFront, PelvisBack};
    enum Muscle   {Reflex11,Reflex12,Reflex21,Reflex22};
    typedef std::pair<BodyType,Side> UniqueBody;
    typedef std::pair<Muscle,Side> UniqueMuscle;

    static const int NBodyType = Torso-Foot+1;
    static const int NSegment  = PelvisBack-FootFront+1;
    static const int NMuscle   = Reflex22-Reflex11+1;

    Dude(Real scale);

    void loadDefaultState(State& state);
    
    void scaleBy(Real scale) {
        m_mass *= scale; m_length *= scale;
        m_segment *= scale;
        for (int i=0; i < NMuscle; ++i) {
            m_springW[i][2] *= scale; // last spring parameter only
            m_springR[i][2] *= scale;
        }
        std::cout << "Masses=" << m_mass << std::endl;
        std::cout << "Lengths=" << m_length << std::endl;
    }

    MultibodySystem             m_system;
    SimbodyMatterSubsystem      m_matter;
    GeneralForceSubsystem       m_forces;
    ContactTrackerSubsystem     m_tracker;
    CompliantContactSubsystem   m_contactForces;
    DecorationSubsystem         m_viz;
    Force::Gravity              m_gravity;

    Vector          m_mass, m_length;     // index by BodyType
    Vector          m_segment;            // index by Segment
    Vector_<Vec3>   m_springW, m_springR; // index by spring #

    std::map<BodyType,   Body>                          m_body;
    std::map<UniqueBody, MobilizedBody>                 m_mobod;
    std::map<UniqueMuscle,TwoPointMuscleDamperReflex*>  m_muscles;
private:
    static Real massData[NBodyType], lengthData[NBodyType];
    static Vec3 springWData[NMuscle], springRData[NMuscle];
    static Real segmentData[NSegment];
};


//////////////////////////////////////////////////////////////////////////
int main() {
  try {
    const Real FrameRate = 30;
    const Real TimeScale = 1;
    const Real scale = 10.;

    Dude dude(scale);

    MultibodySystem& system = dude.m_system;

    Visualizer viz(system);
   
    printf("\n\n***************************************************************\n");
    printf(    "use arrow keys and page up/down to control green gravity vector\n");
    printf(    "***************************************************************\n\n");

    // This menu does nothing.
    Array_< std::pair<std::string,int> > items;
    items.push_back(std::make_pair("One", 1));
    items.push_back(std::make_pair("Top/SubA/first", 2));
    items.push_back(std::make_pair("Top/SubA/second", 3));
    items.push_back(std::make_pair("Top/SubB/only", 4));
    items.push_back(std::make_pair("Two", 5));
    viz.addMenu("Test Menu", 1, items);


    // This is for per-frame camera control and single-frame geometry.
    viz.addFrameController(new MyFrameController(dude.m_matter, 
        MobilizedBodyIndex(1), dude.m_gravity));

    viz.setRealTimeScale(TimeScale);
    //viz.setDesiredBufferLengthInSec(.15);
    viz.setDesiredFrameRate(FrameRate);
    //viz.setMode(Visualizer::Sampling);
    viz.setMode(Visualizer::RealTime);

    // Use this for communication of user input from the GUI to the simulation.
    // Both the Visualizer and the simulation must know about it.
    Visualizer::InputSilo* silo = new Visualizer::InputSilo();
    viz.addInputListener(silo);

    system.addEventHandler(
        new UserInputHandler(*silo, dude.m_gravity, 0.1)); // 100ms
     
    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();
    dude.loadDefaultState(state);

    Assembler(system).assemble(state);

    // Simulate it.

    system.addEventReporter(new Visualizer::Reporter(viz, TimeScale/FrameRate));

    //RungeKutta3Integrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    //RungeKuttaFeldbergIntegrator integ(system);
    //CPodesIntegrator integ(system);
    integ.setAccuracy(.01);
    //integ.setAccuracy(1e-2);
    TimeStepper ts(system, integ);
    ts.initialize(state);

    double cpuStart = cpuTime();
    double realStart = realTime();
    ts.stepTo(Infinity);
    std::cout << "cpu time:  "<<cpuTime()-cpuStart<< std::endl;
    std::cout << "real time: "<<realTime()-realStart<< std::endl;
    std::cout << "steps:     "<<integ.getNumStepsTaken()<< std::endl;
    viz.dumpStats(std::cout);

    std::cout << "Type something to quit: ";
    char ch; std::cin >> ch;

  } catch (const std::exception& exc) {
      std::cout << "EXCEPTION: " << exc.what() << std::endl;
  }
}

Real Dude::massData[] = {.15, .15, .15, .5, .05};
Real Dude::lengthData[] = {.06, .1, .1, .1, .1};
Real Dude::segmentData[] = {lengthData[Foot]*2./3,
                            lengthData[Foot]*1./3,
                            lengthData[Foot]*1./3,
                            lengthData[Shank]*1./3,
                            lengthData[Shank]*1./3,
                            lengthData[Shank]*1./3,
                            lengthData[Shank]*1./5,
                            lengthData[Thigh]*3./10,
                            lengthData[Thigh],
                            lengthData[Pelvis]*1./2,
                            lengthData[Pelvis]*1./2};
static const Real d = 10; // more damping
Vec3 Dude::springWData[] = {Vec3(2000,d*10,.122),
                            Vec3(3000,d*10,.045),
                            Vec3(7000,d*10,.133),
                            Vec3(4000,d*10,.065)};
Vec3 Dude::springRData[] = {Vec3(12000,10,.122),
                            Vec3(8000,10,.045),
                            Vec3(20000,10,.133),
                            Vec3(18000,10,.065)};

Dude::Dude(Real scale)
:   m_matter(m_system), m_forces(m_system), m_tracker(m_system), 
    m_contactForces(m_system, m_tracker), m_viz(m_system),
    m_gravity(m_forces, m_matter, -YAxis, 9.81),
    m_mass(NBodyType,massData), m_length(NBodyType,lengthData),
    m_springW(NMuscle,springWData), m_springR(NMuscle,springRData),
    m_segment(NSegment,segmentData)
{
    const Real opacity = .3;
    const Real resolution = 3;

    scaleBy(scale);
    //Force::GlobalDamper(m_forces, m_matter, 100); // fall slowly

    const Real footHeight = m_length[Foot]  *.1; // half dimensions
    const Real footWidth  = m_length[Foot]  *.3;
    const Real boneRad    = m_length[Thigh] *.1;
    const Real contactRad = /*m_length[Foot]  *.05*/footHeight;

    const Real transitionVelocity = .1;
    m_contactForces.setTransitionVelocity(transitionVelocity);

    // Create bodies
    m_body[Pelvis] = Body::Rigid(MassProperties(m_mass[Pelvis]*1.5,Vec3(0),
        Inertia(1)));
         //Inertia::cylinderAlongY(boneRad,m_length[Pelvis]/2)*m_mass[Pelvis]*.8)));

    m_body[Torso] = Body::Rigid(MassProperties(m_mass[Torso], Vec3(0),
        m_mass[Torso]*Inertia::cylinderAlongY(boneRad,m_length[Torso]/2)));
    m_body[Thigh] = Body::Rigid(MassProperties(m_mass[Thigh], Vec3(0),
        m_mass[Thigh]*Inertia::cylinderAlongY(boneRad,m_length[Thigh]/2)));
    m_body[Shank] = Body::Rigid(MassProperties(m_mass[Shank], Vec3(0),
        m_mass[Shank]*Inertia::cylinderAlongY(boneRad,m_length[Shank]/2)));
    m_body[Foot]  = Body::Rigid(MassProperties(m_mass[Foot],  Vec3(0),
        m_mass[Foot]*Inertia::cylinderAlongY(boneRad,m_length[Foot]/2)));

    // Add DecorativeGeometry to the bodies.
    m_body[Pelvis].addDecoration(Transform(),
        DecorativeBrick(Vec3(m_length[Pelvis]/5, footHeight, footWidth))
        .setColor(Red).setOpacity(opacity));
    m_body[Torso].addDecoration(Transform(),
        DecorativeCylinder(boneRad, m_length[Torso]/2)
        .setColor(Blue).setOpacity(opacity).setResolution(resolution));
    m_body[Thigh].addDecoration(Transform(),
        DecorativeCylinder(boneRad, m_length[Thigh]/2)
        .setColor(Red).setOpacity(opacity).setResolution(resolution));
    m_body[Shank].addDecoration(Transform(),
        DecorativeCylinder(boneRad, m_length[Shank]/2)
        .setColor(Blue).setOpacity(opacity).setResolution(resolution));
    m_body[Foot].addDecoration(Transform(),
        DecorativeBrick(Vec3(footHeight, m_length[Foot]/2, footWidth))
        .setColor(Red).setOpacity(opacity));

    // Add ContactSurfaces to the feet. Since surfaces on the same body can't collide
    // anyway, the clique membership here ensures that the feet can't contact with
    // each other. 
    ContactCliqueId clique1 = ContactSurface::createNewContactClique();

    ContactMaterial material(0.02*1e7, // stiffness
                             0.9,     // dissipation
                             0.8,     // mu_static
                             0.6,     // mu_dynamic
                             1); // mu_viscous

    for (int fb=-1; fb <= 1; fb += 2)
        for (int lr=-1; lr <= 1; lr += 2) {
            const Vec3 ctr(0, fb*(m_length[Foot]/2-contactRad), lr*(footWidth-contactRad));
            m_body[Foot].addContactSurface(ctr,
                ContactSurface(ContactGeometry::Sphere(contactRad),material)
                .joinClique(clique1));
            // Visualize the contact sphere
            m_body[Foot].addDecoration(ctr,
                DecorativeSphere(contactRad).setColor(Green));
        }

    // Half space normal is -x; must rotate to make it +y.
    m_matter.Ground().updBody().addContactSurface(Rotation(-Pi/2,ZAxis),
       ContactSurface(ContactGeometry::HalfSpace(), material));

    // Now create the MobilizedBodies (bodies + joints).
    m_mobod[UniqueBody(Pelvis,Only)] =
        MobilizedBody::Free(m_matter.Ground(), m_body[Pelvis]);
   // What about Torso?

    // Add left and right legs.
    for (Side side=Left; side <= Right; side = Side(side+1)) {
        m_mobod[UniqueBody(Thigh,side)] =
            MobilizedBody::Pin( m_mobod[UniqueBody(Pelvis,Only)],
                Vec3(0,0,side==Left?-footWidth:footWidth),
                m_body[Thigh],
                Transform(Vec3(0, m_length[Thigh]/2, 0)));

        m_mobod[UniqueBody(Shank,side)] =
            MobilizedBody::Pin( m_mobod[UniqueBody(Thigh,side)],
                Transform(Vec3(0, -m_length[Thigh]/2, 0)),
                m_body[Shank],
                Transform(Vec3(0, m_length[Shank]*.5, 0.0)));

        m_mobod[UniqueBody(Foot,side)] =
            MobilizedBody::Pin( m_mobod[UniqueBody(Shank,side)],
                Transform(Vec3(0, -m_length[Shank]/2, 0)),
                m_body[Foot],
                Transform(Vec3(0, m_length[Foot]/2-m_segment[FootHeel], 0)));
    }

    DecorativeLine baseLine;
    baseLine.setColor(Red).setLineThickness(4).setOpacity(.2);

    // Add left and right leg muscles
    for (Side side=Left; side <= Right; side = Side(side+1)) {
        TwoPointMuscleDamperReflex* reflex;

        reflex = new TwoPointMuscleDamperReflex(
            m_mobod[UniqueBody(Pelvis,Only)],
            Vec3(-m_segment[PelvisBack], 0, 0),
            m_mobod[UniqueBody(Shank,side)],
            Vec3(0,m_length[Shank]/2-m_segment[ShankMidUp],0),
            m_springW[0][0], m_springW[0][1], m_springW[0][2] );
        reflex->addDecorativeLine(m_viz, baseLine, .25, 1.0),
        m_muscles[UniqueMuscle(Reflex11,side)] = reflex;
        Force::Custom(m_forces, reflex);

        reflex = new TwoPointMuscleDamperReflex(
            m_mobod[UniqueBody(Shank,side)],
            Vec3(0,-m_length[Shank]/2+m_segment[ShankLower],0),
            m_mobod[UniqueBody(Foot,side)],
            Vec3(0,-m_length[Foot]/2+m_segment[FootToes], 0),
            m_springW[1][0], m_springW[1][1], m_springW[1][2] );
        reflex->addDecorativeLine(m_viz, baseLine);
        m_muscles[UniqueMuscle(Reflex12,side)] = reflex;
        Force::Custom(m_forces, reflex);

        reflex = new TwoPointMuscleDamperReflex(
            m_mobod[UniqueBody(Thigh,side)],
            Vec3(0, -m_length[Thigh]/2+m_segment[ThighLower], 0),
            m_mobod[UniqueBody(Foot,side)],
            Vec3(0,m_length[Foot]/2, 0), /// @todo double check this
            m_springW[2][0], m_springW[2][1], m_springW[2][2] );
        reflex->addDecorativeLine(m_viz, baseLine);
        m_muscles[UniqueMuscle(Reflex21,side)] = reflex;
        Force::Custom(m_forces, reflex);

        reflex = new TwoPointMuscleDamperReflex(
            m_mobod[UniqueBody(Pelvis,Only)],
            Vec3(m_segment[PelvisFront], 0, 0),
            m_mobod[UniqueBody(Shank,side)],
            Vec3(0, m_length[Shank]/2+m_segment[ShankUpper], 0),
            m_springW[3][0], m_springW[3][1], m_springW[3][2] );
        reflex->addDecorativeLine(m_viz, baseLine, .25, 1.0),
        m_muscles[UniqueMuscle(Reflex22,side)] = reflex;
        Force::Custom(m_forces, reflex);
    }
}

void Dude::loadDefaultState(State& state) {
  const static Real hipAngle = -15*Pi/180;
  const static Real kneeAngle = -40*Pi/180;
  const static Real ankleAngle = 70*Pi/180;
  const static Real hipVelocity = .125;
  const static Real kneeVelocity = 0;

  const static Real pelvisZAngle = 20*Pi/180;

  m_mobod[UniqueBody(Pelvis,Only)].setQToFitRotation(state, Rotation(pelvisZAngle,ZAxis));
  m_mobod[UniqueBody(Pelvis,Only)].setQToFitTranslation(state, Vec3(0,2.3,0));
  m_mobod[UniqueBody(Pelvis,Only)].setOneU(state, 2, -hipVelocity);

  m_mobod[UniqueBody(Thigh,Left)].setOneQ(state, 0, -hipAngle);
  m_mobod[UniqueBody(Shank,Left)].setOneQ(state, 0, kneeAngle);
  m_mobod[UniqueBody(Foot,Left)].setOneQ(state, 0, ankleAngle);
  m_mobod[UniqueBody(Thigh,Left)].setOneU(state, 0, hipVelocity);
  m_mobod[UniqueBody(Shank,Left)].setOneU(state, 0, kneeVelocity);

  m_mobod[UniqueBody(Thigh,Right)].setOneQ(state, 0, hipAngle);
  m_mobod[UniqueBody(Shank,Right)].setOneQ(state, 0, kneeAngle);
  m_mobod[UniqueBody(Foot,Right)].setOneQ(state, 0, ankleAngle);
  m_mobod[UniqueBody(Thigh,Right)].setOneU(state, 0, hipVelocity);
  m_mobod[UniqueBody(Shank,Right)].setOneU(state, 0, -kneeVelocity);

}



//////////////////////////////////////////////////////////////////////////

TwoPointMuscleDamperReflex::TwoPointMuscleDamperReflex
   (const MobilizedBody& body1, const Vec3& station1,
    const MobilizedBody& body2, const Vec3& station2,
    Real k, Real damping, Real x0)
:   mBody1(body1), mStation1(station1), mBody2(body2), mStation2(station2),
    mK(k),
    mDamping(damping),
    mX0(x0) {
}

void TwoPointMuscleDamperReflex::calcForce(const State& state,
                                           Vector_<SpatialVec>& bodyForces,
                                           Vector_<Vec3>& particleForces,
                                           Vector& mobilityForces) const {

  const Transform& X_GB1 = mBody1.getBodyTransform(state);
  const Transform& X_GB2 = mBody2.getBodyTransform(state);

  const Vec3 s1_G = X_GB1.R() * mStation1;
  const Vec3 s2_G = X_GB2.R() * mStation2;

  const Vec3 p1_G = X_GB1.p() + s1_G; // mStation measured from ground origin
  const Vec3 p2_G = X_GB2.p() + s2_G;

  const Vec3 r_G       = p2_G - p1_G; // vector from point1 to point2
  const Real dist      = r_G.norm();  // distance between the points
  if( dist < SignificantReal ) return;
  const UnitVec3 dir(r_G);            // direction from point1 to point2

  const Real stretch    = dist - mX0;  // + -> tension, - -> compression
  const Real frcStretch = mK*stretch;  // k(x-x0)

  //////////////////////////////////////////////////////////////////////////

  const Vec3 v1_G = mBody1.findStationVelocityInGround(state, mStation1);
  const Vec3 v2_G = mBody2.findStationVelocityInGround(state, mStation2);
  const Vec3 vRel = v2_G - v1_G;                // relative velocity
  const Real frcDamp = mDamping*dot(vRel, dir); // c*v

  //////////////////////////////////////////////////////////////////////////

  const Vec3 f1_G = (frcStretch + frcDamp) * dir;

  bodyForces[mBody1.getMobilizedBodyIndex()] +=  SpatialVec(s1_G % f1_G, f1_G);
  bodyForces[mBody2.getMobilizedBodyIndex()] -=  SpatialVec(s2_G % f1_G, f1_G);

}

Real TwoPointMuscleDamperReflex::calcPotentialEnergy(const State& state) const {

  const Transform& X_GB1 = mBody1.getBodyTransform(state);
  const Transform& X_GB2 = mBody2.getBodyTransform(state);

  const Vec3 s1_G = X_GB1.R() * mStation1;
  const Vec3 s2_G = X_GB2.R() * mStation2;

  const Vec3 p1_G = X_GB1.p() + s1_G; // mStation measured from ground origin
  const Vec3 p2_G = X_GB2.p() + s2_G;

  const Vec3 r_G     = p2_G - p1_G;   // vector from point1 to point2
  const Real d       = r_G.norm();    // distance between the points
  const Real stretch = d - mX0;       // + -> tension, - -> compression

  return 0.5*mK*stretch*stretch;      // 1/2 k (x-x0)^2

}

