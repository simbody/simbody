/* -------------------------------------------------------------------------- *
 *                       Simbody(tm) - Rigid Contact 1                        *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Thomas Uchida                                                *
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

/*
This is a first try at using experimental Simbody built-in rigid contact,
where the conditional contact constraints must be specified explicitly as
part of the model.
*/

//#define NDEBUG 1

#include "Simbody.h"

#include <string>
#include <iostream>
#include <exception>

using std::cout;
using std::endl;

using namespace SimTK;

//#define USE_TIMS_PARAMS
#define ANIMATE // off to get more accurate CPU time (you can still playback)

//#define HERTZ
#define POISSON


//==============================================================================
//                            MY PUSH FORCE
//==============================================================================
// This is a force element that generates a constant force on a body for a
// specified time interval. It is useful to perturb the system to force it
// to transition from sticking to sliding, for example.
class MyPushForceImpl : public Force::Custom::Implementation {
public:
    MyPushForceImpl(const MobilizedBody& bodyB,
                    const Vec3&          stationB,
                    const Vec3&          forceG, // force direction in Ground!
                    Real                 onTime,
                    Real                 offTime)
    :   m_bodyB(bodyB), m_stationB(stationB), m_forceG(forceG),
        m_on(onTime), m_off(offTime)
    {    }

    //--------------------------------------------------------------------------
    //                       Custom force virtuals
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   override
    {
        if (!(m_on <= state.getTime() && state.getTime() <= m_off))
            return;

        m_bodyB.applyForceToBodyPoint(state, m_stationB, m_forceG, bodyForces);
    }

    // No potential energy.
    Real calcPotentialEnergy(const State& state) const override {return 0;}

    void calcDecorativeGeometryAndAppend
       (const State& state, Stage stage,
        Array_<DecorativeGeometry>& geometry) const override
    {
        const Real ScaleFactor = 0.1;
        if (stage != Stage::Time) return;
        if (!(m_on <= state.getTime() && state.getTime() <= m_off))
            return;
        const Vec3 stationG = m_bodyB.findStationLocationInGround(state, m_stationB);
        geometry.push_back(DecorativeLine(stationG-ScaleFactor*m_forceG, stationG)
                            .setColor(Yellow)
                            .setLineThickness(3));
    }
private:
    const MobilizedBody& m_bodyB;
    const Vec3           m_stationB;
    const Vec3           m_forceG;
    Real                 m_on;
    Real                 m_off;
};


//==============================================================================
//                         AUGMENTED MULTIBODY SYSTEM
//==============================================================================
/* This is a Simbody MultibodySystem with a few bells and whistles added
for visualization.
*/
const Real DefaultCaptureVelocity    = .01,
           DefaultTransitionVelocity = .01;
class AugmentedMultibodySystem : public MultibodySystem {
public:
    AugmentedMultibodySystem()
    :   m_matter(0), m_forces(0),
        m_captureVelocity(DefaultCaptureVelocity),
        m_transitionVelocity(DefaultTransitionVelocity)
    {
        m_matter = new SimbodyMatterSubsystem(*this);
        m_forces = new GeneralForceSubsystem(*this);
        //m_matter->setShowDefaultGeometry(false);
    }

    virtual ~AugmentedMultibodySystem()
    {  delete m_forces; delete m_matter; }

    virtual const MobilizedBody& getBodyToWatch() const
    {   return m_matter->getGround(); }

    virtual void addRubberBandLines(Visualizer& viz) const {}
    virtual Real getRuntime() const {return 20.;}

    virtual Vec3 getWatchOffset() const {return Vec3(0, 8, 1.5);}
    virtual void calcInitialState(State& state) const = 0;

    const SimbodyMatterSubsystem& getMatterSubsystem() const {return *m_matter;}
    SimbodyMatterSubsystem& updMatterSubsystem() {return *m_matter;}

    const GeneralForceSubsystem& getForceSubsystem() const {return *m_forces;}
    GeneralForceSubsystem& updForceSubsystem() {return *m_forces;}

    Real getCaptureVelocity() const {return m_captureVelocity;}
    Real getTransitionVelocity() const {return m_transitionVelocity;}

protected:
    Real m_captureVelocity;
    Real m_transitionVelocity;

private:
    //TODO: this shouldn't require pointers.
    SimbodyMatterSubsystem*     m_matter;
    GeneralForceSubsystem*      m_forces;

};

// Limit single-step direction change to 30 degrees.
static const Real CosMaxSlidingDirChange = std::cos(Pi/6);
static const Real MaxRollingTangVel   = 1.0e-1; //Can't roll above this velocity.


//==============================================================================
//                            SHOW CONTACT
//==============================================================================
// For each visualization frame, generate ephemeral geometry to show just
// during this frame, based on the current State.
class ShowContact : public DecorationGenerator {
public:
    explicit ShowContact(const SemiExplicitEulerTimeStepper& ts)
    :   m_ts(ts) {}

    void generateDecorations(const State&                state,
                             Array_<DecorativeGeometry>& geometry) override
    {
        const MultibodySystem& mbs = m_ts.getMultibodySystem();
        const Real TextScale = mbs.getDefaultLengthScale()/10; // was .1
        mbs.realize(state, Stage::Dynamics);
        const Real KE=mbs.calcKineticEnergy(state), E=mbs.calcEnergy(state);
        DecorativeText energy; energy.setIsScreenText(true);
        energy.setText("Energy/KE: " + String(E, "%.3f") + String(KE, "/%.3f"));
        geometry.push_back(energy);

        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
        const int nUniContacts  = matter.getNumUnilateralContacts();
        const int nLtdFrictions = matter.getNumStateLimitedFrictions();

        for (UnilateralContactIndex ux(0); ux < nUniContacts; ++ux) {
            const UnilateralContact& contact = matter.getUnilateralContact(ux);
            const Vec3 loc = contact.whereToDisplay(state);
            if (contact.isEnabled(state)) {
                MultiplierIndex mx = contact.getContactMultiplierIndex(state);
                const Vector& mults=state.updMultipliers();//TODO: can't use get
                Vec3 color = std::abs(mults[mx])<SignificantReal
                    ? Yellow : Cyan;
                geometry.push_back(DecorativeSphere(.1)
                    .setTransform(loc).setColor(color).setOpacity(.5));
                //contact.showContactForce(state, geometry);
                String text;
                if (!contact.hasFriction(state))
                    text = "-ENABLED";
                geometry.push_back(DecorativeText(String((int)ux)+text)
                    .setColor(White).setScale(TextScale)
                    .setTransform(loc+Vec3(0,.04,0)));
            } else {
                geometry.push_back(DecorativeText(String((int)ux))
                    .setColor(White).setScale(TextScale)
                    .setTransform(loc+Vec3(0,.02,0)));
            }
        }

        //for (unsigned i=0; i < m_ts.m_proximals.m_friction.size(); ++i) {
        //    const int id = m_ts.m_proximals.m_friction[i];
        //    const MyFrictionElement& felt = unis.getFrictionElement(id);
        //    const Vec3 loc = felt.whereToDisplay(state);
        //    const Frictional& fric = m_ts.m_frictional[i];
        //    String text = fric.m_wasLimited ? "slip" : "STICK";
        //    Vec3 color = fric.m_wasLimited ? Green : Orange;
        //    felt.showFrictionForce(state, geometry, color);
        //    geometry.push_back(DecorativeText(text)
        //        .setColor(color).setScale(TextScale)
        //        .setTransform(loc+Vec3(0.1,.04,0)));
        //}
    }
private:
    const SemiExplicitEulerTimeStepper& m_ts;
};

//==============================================================================
//                               TIM'S BOX
//==============================================================================
class TimsBox : public AugmentedMultibodySystem {
public:
    TimsBox();

    void calcInitialState(State& state) const override;

    const MobilizedBody& getBodyToWatch() const override
    {   return m_brick; }

    Vec3 getWatchOffset() const override
    {   return Vec3(0, 8, 8); }

private:
    Force::Gravity          m_gravity;
    Force::GlobalDamper     m_damper;
    MobilizedBody::Free     m_brick;
    MobilizedBody::Ball     m_brick2;
};


//==============================================================================
//                              BOUNCING BALLS
//==============================================================================
class BouncingBalls : public AugmentedMultibodySystem {
public:
    BouncingBalls();
    ~BouncingBalls() {delete m_contactForces; delete m_tracker;}

    void calcInitialState(State& state) const override;

    const MobilizedBody& getBodyToWatch() const override
    {   static const MobilizedBody nobod; return nobod; }
    Vec3 getWatchOffset() const override {return Vec3(0, 8, 20.);}

    const MobilizedBody::Slider& getHBall(int i) const {return m_Hballs[i];}
    const MobilizedBody::Slider& getPBall(int i) const {return m_Pballs[i];}

private:
    // Add subsystems for compliant contact. TODO: shouldn't need pointers
    ContactTrackerSubsystem*     m_tracker;
    CompliantContactSubsystem*   m_contactForces;

    static const int NBalls = 10;

    Force::Gravity          m_gravity;
    Force::GlobalDamper     m_damper;
    MobilizedBody::Slider   m_Hballs[NBalls];    // Hertz
    MobilizedBody::Slider   m_Pballs[NBalls];    // Poisson
};

//==============================================================================
//                                  PENCIL
//==============================================================================
class Pencil : public AugmentedMultibodySystem {
public:
    Pencil();
    ~Pencil() {delete m_contactForces; delete m_tracker;}

    void calcInitialState(State& state) const override;

    const MobilizedBody& getBodyToWatch() const override
    {   static const MobilizedBody none; return none; }
    Vec3 getWatchOffset() const override {return Vec3(0, 8, 30.);}
    void addRubberBandLines(Visualizer& viz) const override;
    Real getRuntime() const override {return 75.;}

    const MobilizedBody::Planar& getPencil() const {return m_pencil;}
    const MobilizedBody::Free& getBall() const {return m_ball;}

private:
    // Add subsystems for compliant contact. TODO: shouldn't need pointers
    ContactTrackerSubsystem*     m_tracker;
    CompliantContactSubsystem*   m_contactForces;

    Force::Gravity          m_gravity;
    Force::GlobalDamper     m_damper;
    MobilizedBody::Planar   m_pencil;
    MobilizedBody::Pin      m_flapper;
    MobilizedBody::Free     m_ball;

    Vec3                    m_pencilSphereCenter;
};


//==============================================================================
//                                  BLOCK
//==============================================================================
// This just a single block sitting on the ground with an off-center external
// force.
class Block : public AugmentedMultibodySystem {
public:
    Block();
    ~Block() {}

    void calcInitialState(State& state) const override;
    const MobilizedBody& getBodyToWatch() const override
    {   return m_block; }

    Vec3 getWatchOffset() const override
    {   return Vec3(0, 8, 20); }

    const MobilizedBody::Free& getBlock() const {return m_block;}

private:
    Force::Gravity          m_gravity;
    Force::GlobalDamper     m_damper;
    MobilizedBody::Free     m_block;
    MobilizedBody::Pin      m_link2;
    MobilizedBody::Pin      m_link3;
};



//==============================================================================
//                                  EDGES
//==============================================================================
// Two blocks engaged in unilateral edge-edge contact.
class Edges : public AugmentedMultibodySystem {
public:
    Edges();
    ~Edges() {}

    void calcInitialState(State& state) const override;

    const MobilizedBody& getBodyToWatch() const override
    {   static const MobilizedBody mobod; return mobod; }
    Vec3 getWatchOffset() const override
    {   return Vec3(0, 8, 30); }

    const MobilizedBody& getSwingingBlock() const {return m_swinger;}
    const MobilizedBody::Free& getGroundBlock() const {return m_grounded;}

private:
    Force::Gravity          m_gravity;
    Force::GlobalDamper     m_damper;
    MobilizedBody           m_swinger;
    MobilizedBody::Free     m_grounded;
};

//==============================================================================
//                                   MAIN
//==============================================================================
int main(int argc, char** argv) {
    #ifdef USE_TIMS_PARAMS
        const Real RunTime=16;  // Tim's time
        const Real Accuracy = 1e-4;
    #else
        const Real RunTime=20;
        const Real Accuracy = 1e-2;
    #endif

  try { // If anything goes wrong, an exception will be thrown.

    // Create the augmented multibody model.
    //TimsBox mbs;
    //BouncingBalls mbs;
    //Pencil mbs;
    //Block mbs;
    Edges mbs;

    SemiExplicitEulerTimeStepper sxe(mbs);
    sxe.setDefaultImpactCaptureVelocity(mbs.getCaptureVelocity());
    sxe.setDefaultImpactMinCORVelocity(mbs.getCaptureVelocity()); //TODO
    sxe.setDefaultFrictionTransitionVelocity(mbs.getTransitionVelocity());

    const SimbodyMatterSubsystem&    matter = mbs.getMatterSubsystem();

    Visualizer viz(mbs);
    viz.setShowSimTime(true);
    viz.setShowFrameNumber(true);
    viz.setShowFrameRate(true);
    viz.addDecorationGenerator(new ShowContact(sxe));
    mbs.addRubberBandLines(viz);

    if (!mbs.getBodyToWatch().isEmptyHandle())
        viz.addFrameController(
                new Visualizer::BodyFollower(mbs.getBodyToWatch(), Vec3(0),
                                             mbs.getWatchOffset()));

    // Simulate it.
    State s;
    mbs.calcInitialState(s);

    printf("Initial state shown. ENTER to continue.\n");
    viz.report(s);
    getchar();

    const Real ConsTol = .001;

    sxe.setRestitutionModel(SemiExplicitEulerTimeStepper::Poisson);
    //sxe.setRestitutionModel(SemiExplicitEulerTimeStepper::Newton);
    //sxe.setRestitutionModel(SemiExplicitEulerTimeStepper::NoRestitution);

    sxe.setPositionProjectionMethod(SemiExplicitEulerTimeStepper::Bilateral);
    //sxe.setPositionProjectionMethod(SemiExplicitEulerTimeStepper::Unilateral);
    //sxe.setPositionProjectionMethod(SemiExplicitEulerTimeStepper::NoPositionProjection);

    sxe.setAccuracy(Accuracy); // integration accuracy
    sxe.setConstraintTolerance(ConsTol);

    //sxe.setImpulseSolverType(SemiExplicitEulerTimeStepper::PGS);
    sxe.setImpulseSolverType(SemiExplicitEulerTimeStepper::PLUS);

    sxe.setInducedImpactModel(SemiExplicitEulerTimeStepper::Simultaneous);
    //sxe.setInducedImpactModel(SemiExplicitEulerTimeStepper::Sequential);
    //sxe.setInducedImpactModel(SemiExplicitEulerTimeStepper::Mixed);

    sxe.setMaxInducedImpactsPerStep(1000);
    //sxe.setMaxInducedImpactsPerStep(2);

    printf("RestitutionModel: %s, InducedImpactModel: %s (maxIts=%d)\n",
        sxe.getRestitutionModelName(sxe.getRestitutionModel()),
        sxe.getInducedImpactModelName(sxe.getInducedImpactModel()),
        sxe.getMaxInducedImpactsPerStep());
    printf("PositionProjectionMethod: %s, ImpulseSolverType=%s\n",
        sxe.getPositionProjectionMethodName(sxe.getPositionProjectionMethod()),
        sxe.getImpulseSolverTypeName(sxe.getImpulseSolverType()));


    sxe.initialize(s);
    mbs.resetAllCountersToZero();

    Array_<State> states; states.reserve(10000);

    int nSteps=0, nStepsWithEvent = 0;

    const double startReal = realTime();
    const double startCPU = cpuTime();

    const Real h = .001;
    const int SaveEvery = 33; // save every nth step ~= 33ms

    Vector_<SpatialVec> reactions;
    do {
        const State& sxeState = sxe.getState();
        if ((nSteps%SaveEvery)==0) {
            #ifdef ANIMATE
            viz.report(sxeState);
            //cout << "\n********** t=" << sxe.getTime();
            //cout << " u=" << sxeState.getU()[3] <<"\n\n";
            //#ifndef NDEBUG
            //printf("\nWAITING:"); getchar();
            //#endif
            #endif
            states.push_back(sxeState);
        }

        sxe.stepTo(sxeState.getTime() + h);

        //printf("Reaction forces:\n");
        //matter.calcMobilizerReactionForces(sxeState, reactions);
        //cout << reactions << endl;

        ++nSteps;
    } while (sxe.getTime() < mbs.getRuntime());
    // TODO: did you lose the last step?


    const double timeInSec = realTime()-startReal;
    const double cpuInSec = cpuTime()-startCPU;
    cout << "Done -- took " << nSteps << " steps in " <<
        timeInSec << "s for " << sxe.getTime() << "s sim (avg step="
        << (1000*sxe.getTime())/nSteps << "ms) ";
    cout << "CPUtime " << cpuInSec << endl;

    printf("Used SXETimeStepper (%s) at acc=%g consTol=%g\n",
           sxe.getRestitutionModel()==SemiExplicitEulerTimeStepper::Newton
            ? "Newton" : "Poisson",
           sxe.getAccuracyInUse(), sxe.getConstraintToleranceInUse());

    //       pgs.getPGSConvergenceTol(), pgs.getPGSMaxIters(),
    //       pgs.getPGSSOR());
    //printf("Compression: ncalls=%lld, niters=%lld (%g/call), nfail=%lld\n",
    //       pgs.getPGSNumCalls(0), pgs.getPGSNumIters(0),
    //       (double)pgs.getPGSNumIters(0)/std::max(pgs.getPGSNumCalls(0),1LL),
    //       pgs.getPGSNumFailures(0));
    //printf("Expansion: ncalls=%lld, niters=%lld (%g/call), nfail=%lld\n",
    //       pgs.getPGSNumCalls(1), pgs.getPGSNumIters(1),
    //       (double)pgs.getPGSNumIters(1)/std::max(pgs.getPGSNumCalls(1),1LL),
    //       pgs.getPGSNumFailures(1));
    //printf("Position: ncalls=%lld, niters=%lld (%g/call), nfail=%lld\n",
    //       pgs.getPGSNumCalls(2), pgs.getPGSNumIters(2),
    //       (double)pgs.getPGSNumIters(2)/std::max(pgs.getPGSNumCalls(2),1LL),
    //       pgs.getPGSNumFailures(2));

    cout << "nstates =" << states.size() << endl;

    // Instant replay.
    while(true) {
        printf("Hit ENTER for replay (%d states) ...",
                states.size());
        getchar();
        for (unsigned i=0; i < states.size(); ++i) {
            mbs.realize(states[i], Stage::Velocity);
            viz.report(states[i]);
        }
    }

  }
  catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
  catch (...) {
    printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }
}

//==============================================================================
//                               TIM'S BOX
//==============================================================================
TimsBox::TimsBox() {
    // Abbreviations.
    SimbodyMatterSubsystem&     matter = updMatterSubsystem();
    GeneralForceSubsystem&      forces = updForceSubsystem();
    MobilizedBody&              Ground = matter.updGround();

    // Build the multibody system.
    m_gravity = Force::Gravity(forces, matter, -YAxis, 9.8066);
    m_damper  = Force::GlobalDamper(forces, matter, .1);

    // Predefine some handy rotations.
    const Rotation Z90(Pi/2, ZAxis); // rotate +90 deg about z

    const Vec3 BrickHalfDims(.1, .25, .5);
    const Real BrickMass = /*10*/5;
    #ifdef USE_TIMS_PARAMS
        const Real RunTime=16;  // Tim's time
        const Real Stiffness = 2e7;
        const Real Dissipation = 1;
        const Real CoefRest = 0;
        // Painleve problem with these friction coefficients.
        //const Real mu_d = 1; /* compliant: .7*/
        //const Real mu_s = 1; /* compliant: .7*/
        const Real mu_d = .5;
        const Real mu_s = .8;
        const Real mu_v = /*0.05*/0;
        const Real CaptureVelocity = 0.01;
        const Real TransitionVelocity = 0.01;
        const Inertia brickInertia(.1,.1,.1);
        const Real Radius = .02;
    #else
        const Real RunTime=20;
        const Real Stiffness = 1e6;
        const Real CoefRest = 0.4;
        const Real TargetVelocity = 3; // speed at which to match coef rest
//        const Real Dissipation = (1-CoefRest)/TargetVelocity;
        const Real Dissipation = 0.1;
        const Real mu_d = .4;
        const Real mu_s = .8;
        const Real mu_v = 0*0.1;
        const Real CaptureVelocity = 0.01;
        const Real TransitionVelocity = 0.05;
        const Inertia brickInertia(BrickMass*UnitInertia::brick(BrickHalfDims));
        const Real Radius = BrickHalfDims[0]/3;
    #endif

    m_captureVelocity = CaptureVelocity;
    m_transitionVelocity = TransitionVelocity;

    printf("\n******************** Tim's Box ********************\n");
    printf("USING RIGID CONTACT\n");
    #ifdef USE_TIMS_PARAMS
    printf("Using Tim's parameters:\n");
    #else
    printf("Using Sherm's parameters:\n");
    #endif
    printf("  coef restitution=%g\n", CoefRest);
    printf("  mu_d=%g mu_s=%g mu_v=%g\n", mu_d, mu_s, mu_v);
    printf("  transition velocity=%g\n", TransitionVelocity);
    printf("  radius=%g\n", Radius);
    printf("  brick inertia=%g %g %g\n",
        brickInertia.getMoments()[0], brickInertia.getMoments()[1],
        brickInertia.getMoments()[2]);
    printf("******************** Tim's Box ********************\n\n");

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS

    Body::Rigid brickBody =
        Body::Rigid(MassProperties(BrickMass, Vec3(0), brickInertia));
    brickBody.addDecoration(Transform(), DecorativeBrick(BrickHalfDims)
                                   .setColor(Red).setOpacity(.3));
    m_brick = MobilizedBody::Free(Ground, Vec3(0),
                                  brickBody, Vec3(0));
    m_brick2 = MobilizedBody::Ball(m_brick, BrickHalfDims,
                                   brickBody, Vec3(-BrickHalfDims));
    //m_brick3 = MobilizedBody::Ball(brick2, BrickHalfDims,
    //                          brickBody, Vec3(-BrickHalfDims));

/*
1) t= 0.5, dt = 2 sec, pt = (0.05, 0.2, 0.4), fdir = (1,0,0), mag = 50N
2) t= 4.0, dt = 0.5 sec, pt = (0.03, 0.06, 0.09), fdir = (0.2,0.8,0), mag = 300N
3) t= 0.9, dt = 2 sec, pt = (0,0,0), fdir = (0,1,0), mag = 49.333N (half the weight of the block)
4) t= 13.0, dt = 1 sec, pt = (0 0 0), fdir = (-1,0,0), mag = 200N
*/
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0.05,0.2,0.4),
                                                    50 * Vec3(1,0,0),
                                                    0.5, 0.5+2));
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0.03, 0.06, 0.09),
                                                    300 * UnitVec3(0.2,0.8,0),
                                                    //300 * Vec3(0.2,0.8,0),
                                                    4, 4+0.5));
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0),
                                                    1.25*49.033 * Vec3(0,1,0),
                                                    9., 9.+2));
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0),
                                                    200 * Vec3(-1,0,0),
                                                    13, 13+1));

    #ifndef USE_TIMS_PARAMS
    // Extra late force.
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(.1, 0, .45),
                                                    20 * Vec3(-1,-1,.5),
                                                    15, 18));
    #endif

    const Real CornerRadius = .1;
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(BrickHalfDims);
        //PointPlaneContact* contact = new PointPlaneContact
         //  (Ground, YAxis, 0., m_brick, pt, CoefRest, mu_s, mu_d, mu_v);
        SpherePlaneContact* contact = new SpherePlaneContact
           (Ground, YAxis, 0., m_brick, pt, CornerRadius,
            CoefRest, mu_s, mu_d, mu_v);
        m_brick.addBodyDecoration(pt,
            DecorativeSphere(CornerRadius).setColor(Gray).setOpacity(.5));
        matter.adoptUnilateralContact(contact);

        if (i==-1 && j==-1 && k==-1)
            continue;
        PointPlaneContact* contact2 = new PointPlaneContact
           (Ground, YAxis, 0., m_brick2, pt, CoefRest, mu_s, mu_d, mu_v);
        matter.adoptUnilateralContact(contact2);

        //PointPlaneContact* contact3 = new PointPlaneContact
        //  (Ground, YAxis, 0., m_brick3, pt, CoefRest, mu_s, mu_d, mu_v);
        //matter.adoptUnilateralContact(contact3);
    }

    // Midline contacts
    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2) {
        const Vec3 pt = Vec3(i,j,0).elementwiseMultiply(BrickHalfDims);
        //PointPlaneContact* contact = new PointPlaneContact
        //   (Ground, YAxis, 0., m_brick, pt, CoefRest, mu_s, mu_d, mu_v);
        SpherePlaneContact* contact = new SpherePlaneContact
           (Ground, YAxis, 0., m_brick, pt, CornerRadius,
            CoefRest, mu_s, mu_d, mu_v);
        m_brick.addBodyDecoration(pt,
            DecorativeSphere(CornerRadius).setColor(Gray).setOpacity(.5));
        matter.adoptUnilateralContact(contact);

        PointPlaneContact* contact2 = new PointPlaneContact
           (Ground, YAxis, 0., m_brick2, pt, CoefRest, mu_s, mu_d, mu_v);
        matter.adoptUnilateralContact(contact2);    }
}

//---------------------------- CALC INITIAL STATE ------------------------------
void TimsBox::calcInitialState(State& s) const {
    s = realizeTopology(); // returns a reference to the the default state

    //matter.setUseEulerAngles(s, true);

    realizeModel(s); // define appropriate states for this System
    realize(s, Stage::Instance); // instantiate constraints if any


    /*
    rX_q = 0.7 rad
    rX_u = 1.0 rad/sec

    rY_q = 0.6 rad
    rY_u = 0.0 rad/sec

    rZ_q = 0.5 rad
    rZ_u = 0.2 rad/sec

    tX_q = 0.0 m
    tX_u = 10 m/s

    tY_q = 1.4 m
    tY_u = 0.0 m/s

    tZ_q = 0.0 m
    tZ_u = 0.0 m/s
    */

    #ifdef USE_TIMS_PARAMS
        m_brick.setQToFitTranslation(s, Vec3(0,10,0));
        m_brick.setUToFitLinearVelocity(s, Vec3(0,0,0));
    #else
        m_brick.setQToFitTranslation(s, Vec3(0,1.4,0));
        m_brick.setUToFitLinearVelocity(s, Vec3(10,0,0));
        const Rotation R_BC(SimTK::BodyRotationSequence,
                                    0.7, XAxis, 0.6, YAxis, 0.5, ZAxis);
        m_brick.setQToFitRotation(s, R_BC);
        m_brick.setUToFitAngularVelocity(s, Vec3(1,0,.2));
    #endif

    realize(s, Stage::Position);
    Assembler(*this).setErrorTolerance(1e-6).assemble(s);
}

//==============================================================================
//                              BOUNCING BALLS
//==============================================================================
static const Real Separation = 0*.0011;
static const Real Gravity = 9.8066;
static const Real Heavy = 10; // ratio Heavy:Light

BouncingBalls::BouncingBalls() {
    m_tracker       = new ContactTrackerSubsystem(*this);
    m_contactForces = new CompliantContactSubsystem(*this, *m_tracker);

    // Abbreviations.
    SimbodyMatterSubsystem&     matter = updMatterSubsystem();
    GeneralForceSubsystem&      forces = updForceSubsystem();
    MobilizedBody&              Ground = matter.updGround();


    // Build the multibody system.
    m_gravity = Force::Gravity(forces, matter, -YAxis, Gravity);
    //m_damper  = Force::GlobalDamper(forces, matter, .1);

    const Real BallMass = 1;
    const Real BallRadius = .25;
    const Real CoefRest = 1;
    const Real CaptureVelocity = .001;
    const Real TransitionVelocity = .001;

    // Rubber
    const Real rubber_density = 1100.;  // kg/m^3
    const Real rubber_young   = 0.01e9; // pascals (N/m)
    const Real rubber_poisson = 0.5;    // ratio
    const Real rubber_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
    const Real rubber_dissipation = 0.1;
    const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,0,0,0);
    // Nylon
    const Real nylon_density = 1100.;  // kg/m^3
    const Real nylon_young   = 10*2.5e9;  // pascals (N/m)
    const Real nylon_poisson = 0.4;    // ratio
    const Real nylon_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(nylon_young,nylon_poisson);
    const Real nylon_dissipation = 0.1;
    const ContactMaterial nylon(nylon_planestrain,nylon_dissipation,0,0,0);
    const ContactMaterial nylon_lossless
       (nylon_planestrain,0/*no dissipation*/,0,0,0);
    const ContactMaterial nylon_lossy
       (nylon_planestrain,10/*much dissipation*/,0,0,0);

    const Rotation X2Y(Pi/2, ZAxis); // rotate +90 deg about z
    const Rotation NegX2Y(-Pi/2,ZAxis); // -90

    Ground.updBody().addContactSurface(Transform(NegX2Y,Vec3(0)),
                ContactSurface(ContactGeometry::HalfSpace(),nylon));

    m_captureVelocity = CaptureVelocity;
    m_transitionVelocity = TransitionVelocity;

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS

    Body::Rigid ballBody(MassProperties(BallMass, Vec3(0),
                                        UnitInertia::sphere(BallRadius)));
    ballBody.addDecoration(Transform(), DecorativeSphere(BallRadius));

    Body::Rigid ballBody_heavy(MassProperties(Heavy*BallMass, Vec3(0),
                                        UnitInertia::sphere(BallRadius)));
    ballBody_heavy.addDecoration(Transform(), DecorativeSphere(BallRadius));

    const Vec3 HColor(Gray), PColor(Red), NColor(Orange);

#ifdef HERTZ
    m_Hballs[0] = MobilizedBody::Slider
       (Ground, Transform(X2Y,Vec3(-1,BallRadius,0)),
        ballBody, X2Y);
    m_Hballs[0].updBody().addContactSurface(Vec3(0),
            ContactSurface(ContactGeometry::Sphere(BallRadius),
            nylon
            //nylon_lossless
            ));
    m_Hballs[0].updBody().updDecoration(0).setColor(HColor);
    for (int i=1; i<NBalls; ++i) {
        m_Hballs[i] = MobilizedBody::Slider
           (m_Hballs[i-1],Transform(X2Y,Vec3(0,2*BallRadius,0)),ballBody, X2Y);
        m_Hballs[i].updBody().updDecoration(0).setColor(HColor);
        m_Hballs[i].updBody().addContactSurface(Vec3(0),
                ContactSurface(ContactGeometry::Sphere(BallRadius),
                //nylon
                nylon_lossless
                ));
    }
#endif
#ifdef POISSON
    m_Pballs[0] = MobilizedBody::Slider
       (Ground, Transform(X2Y,Vec3(0,BallRadius,0)),
        ballBody,
        //ballBody_heavy,
        X2Y);
    m_Pballs[0].updBody().updDecoration(0).setColor(PColor);
    matter.adoptUnilateralContact(new PointPlaneContact
           (Ground, YAxis, 0., m_Pballs[0], Vec3(0,-BallRadius,0),
            CoefRest,
            0,0,0)); // no friction
    for (int i=1; i<NBalls; ++i) {
        m_Pballs[i] = MobilizedBody::Slider
           (m_Pballs[i-1],Transform(X2Y,Vec3(0,2*BallRadius,0)),
           i==NBalls-1?ballBody_heavy:ballBody,
           X2Y);
        m_Pballs[i].updBody().updDecoration(0).setColor(PColor);
        //Real cor = i==NBalls/2 ? .5 : CoefRest; // middle ball different
        Real cor = CoefRest;
        matter.adoptUnilateralContact(new PointPlaneFrictionlessContact
               (m_Pballs[i-1], YAxis, BallRadius,
                m_Pballs[i], Vec3(0,-BallRadius,0), cor));
    }

#endif

}

void BouncingBalls::calcInitialState(State& s) const {
    const Real Height = 1;
    const Real Speed = -2;

    s = realizeTopology(); // returns a reference to the the default state
    realizeModel(s); // define appropriate states for this System
    realize(s, Stage::Instance); // instantiate constraints if any
    realize(s, Stage::Position);
    Assembler(*this).setErrorTolerance(1e-6).assemble(s);
    #ifdef HERTZ
        getHBall(0).setOneQ(s, MobilizerQIndex(0), Height);
        getHBall(0).setOneU(s, MobilizerUIndex(0), Speed);
        for (int i=1; i<NBalls; ++i)
            getHBall(i).setOneQ(s, MobilizerQIndex(0), Separation);
    #endif
    #ifdef POISSON
        getPBall(0).setOneQ(s, MobilizerQIndex(0), Height);
        getPBall(0).setOneU(s, MobilizerUIndex(0), Speed);
        for (int i=1; i<NBalls; ++i)
            getPBall(i).setOneQ(s, MobilizerQIndex(0), Separation);
    #endif
}

//==============================================================================
//                                  PENCIL
//==============================================================================

Pencil::Pencil() {
    m_tracker       = new ContactTrackerSubsystem(*this);
    m_contactForces = new CompliantContactSubsystem(*this, *m_tracker);

    // Abbreviations.
    SimbodyMatterSubsystem&     matter = updMatterSubsystem();
    GeneralForceSubsystem&      forces = updForceSubsystem();
    MobilizedBody&              Ground = matter.updGround();


    // Build the multibody system.
    m_gravity = Force::Gravity(forces, matter, -YAxis, 9.8066);
    //m_damper  = Force::GlobalDamper(forces, matter, .1);

    const Real PencilMass = 1;
    const Real PencilRadius = .25;
    const Real PencilHLength = 5;
    const Real CoefRest = .75;
    const Real StopCoefRest = .7;
    const Real CaptureVelocity = .001;
    const Real TransitionVelocity = .01;
    //const Real mu_d=10, mu_s=10, mu_v=0; // TODO: PAINLEVE!
    //const Real mu_d=1, mu_s=1, mu_v=0;
    const Real mu_d=.3, mu_s=.4, mu_v=0;
    //const Real mu_d=0, mu_s=0, mu_v=0;

    setDefaultLengthScale(PencilHLength);

    // Rubber
    const Real rubber_density = 1100.;  // kg/m^3
    const Real rubber_young   = 0.01e9; // pascals (N/m)
    const Real rubber_poisson = 0.5;    // ratio
    const Real rubber_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
    const Real rubber_dissipation = 0.1;
    const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,0,0,0);
    // Nylon
    const Real nylon_density = 1100.;  // kg/m^3
    const Real nylon_young   = 2.5e9;  // pascals (N/m)
    const Real nylon_poisson = 0.4;    // ratio
    const Real nylon_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(nylon_young,nylon_poisson);
    const Real nylon_dissipation = 0*0.1;
    const ContactMaterial nylon(nylon_planestrain,nylon_dissipation,0,0,0);

    const Rotation X2Y(Pi/2, ZAxis); // rotate +90 deg about z
    const Rotation NegX2Y(-Pi/2,ZAxis); // -90

    Ground.updBody().addContactSurface(Transform(NegX2Y,Vec3(0)),
                ContactSurface(ContactGeometry::HalfSpace(),nylon));

    m_captureVelocity = CaptureVelocity;
    m_transitionVelocity = TransitionVelocity;



        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS

    Body::Rigid pencilBody(MassProperties(PencilMass, Vec3(0),
           UnitInertia::cylinderAlongY(PencilRadius,PencilHLength)));
    pencilBody.addDecoration(Transform(),
                             DecorativeCylinder(PencilRadius,PencilHLength)
                             .setOpacity(.3));

    m_pencil = MobilizedBody::Planar
       (Ground, Vec3(0,PencilHLength,0), pencilBody, Vec3(0));

    const Real flapperRadius = 0.2;
    m_flapper = MobilizedBody::Pin
       (m_pencil, Vec3(0),
        MassProperties(10,Vec3(0),UnitInertia(0)), Vec3(-2,0,0));
    m_flapper.addBodyDecoration(Vec3(0), DecorativeLine(Vec3(-2,0,0),Vec3(0)));
    HardStopLower* lower =
        new HardStopLower(m_flapper, MobilizerQIndex(0), -Pi/6, StopCoefRest);
    matter.adoptUnilateralContact(lower);
    HardStopUpper* upper =
        new HardStopUpper(m_flapper, MobilizerQIndex(0), Pi/6, StopCoefRest);
    matter.adoptUnilateralContact(upper);
    SpherePlaneContact* flapperContact =
        new SpherePlaneContact(Ground, YAxis, 0.,
                                    m_flapper, Vec3(0), flapperRadius,
                                    CoefRest/2, mu_s, mu_d, mu_v);
    matter.adoptUnilateralContact(flapperContact);
    m_flapper.addBodyDecoration(Vec3(0,0,0),
                              DecorativeSphere(flapperRadius).setColor(Red));

    //PointPlaneContact* pc1 =
    //    new PointPlaneContact(Ground, YAxis, 0.,
    //                                m_pencil, Vec3(0,-PencilHLength,0),
    //                                CoefRest, mu_s, mu_d, mu_v);
    Real radius1 = 1;
    m_pencilSphereCenter = Vec3(0,-PencilHLength,0);
    SpherePlaneContact* pc1 =
        new SpherePlaneContact(Ground, YAxis, 0.,
                                    m_pencil, m_pencilSphereCenter, radius1,
                                    CoefRest/3, mu_s, mu_d, mu_v);
    m_pencil.addBodyDecoration(m_pencilSphereCenter,
        DecorativeSphere(radius1).setColor(Cyan).setOpacity(.5)
        .setResolution(3));
    m_pencil.addBodyDecoration(m_pencilSphereCenter,
        DecorativeSphere(radius1).setColor(Black).
        setRepresentation(DecorativeGeometry::DrawWireframe));
    PointPlaneContact* pc2 =
        new PointPlaneContact(Ground, YAxis, 0.,
                                    m_pencil, Vec3(0,PencilHLength,0),
                                    CoefRest, mu_s, mu_d, mu_v);
    matter.adoptUnilateralContact(pc1);
    matter.adoptUnilateralContact(pc2);

    const Real BallRad = 2, BallMass = 20;
    m_ball = MobilizedBody::Free(Ground, Vec3(0),
        MassProperties(BallMass, Vec3(0), UnitInertia::sphere(BallRad)),
        Vec3(0));
    m_ball.addBodyDecoration(Vec3(0),
                             DecorativeSphere(BallRad).setColor(Red)
                             .setResolution(3).setOpacity(.25));
    m_ball.addBodyDecoration(Vec3(0),
                             DecorativeSphere(BallRad).setColor(Black)
                             .setRepresentation(DecorativeGeometry::DrawWireframe)
                             );
    SpherePlaneContact* pc3 =
        new SpherePlaneContact(Ground, YAxis, 0.,
                               m_ball, Vec3(0), BallRad,
                               CoefRest/3, mu_s, mu_d, mu_v);
    matter.adoptUnilateralContact(pc3);

    SphereSphereContact* ss1 =
        new SphereSphereContact(m_pencil, m_pencilSphereCenter, radius1,
                                m_ball, Vec3(0), BallRad,
                                CoefRest/2, mu_s, mu_d, mu_v);
    matter.adoptUnilateralContact(ss1);


    SphereSphereContact* ss2 =
        new SphereSphereContact(m_flapper, Vec3(0), flapperRadius,
                                m_ball, Vec3(0), BallRad,
                                CoefRest/2, mu_s, mu_d, mu_v);
    matter.adoptUnilateralContact(ss2);

    Force::TwoPointLinearSpring(forces, m_pencil, m_pencilSphereCenter,
                                m_ball, Vec3(0), 7.5, 0);

    Rope* rope =
        new Rope(Ground, Vec3(1,10,-1),
                 m_ball, Vec3(0,BallRad,0),
                 6.5, 1);
    matter.adoptUnilateralContact(rope);
}

void Pencil::addRubberBandLines(Visualizer& viz) const {
    viz.addRubberBandLine(m_pencil, m_pencilSphereCenter,
                          m_ball,   Vec3(0),
                          DecorativeLine().setColor(Orange).setLineThickness(3));
}


void Pencil::calcInitialState(State& s) const {
    s = realizeTopology(); // returns a reference to the the default state
    realizeModel(s); // define appropriate states for this System
    realize(s, Stage::Instance); // instantiate constraints if any
    realize(s, Stage::Position);
    getPencil().setOneQ(s, MobilizerQIndex(0), .01);
    getPencil().setOneQ(s, MobilizerQIndex(2), 2);
    //getPencil().setOneQ(s, MobilizerQIndex(0), Pi/4);
    //getPencil().setOneQ(s, MobilizerQIndex(2), -1);
    getPencil().setOneU(s, MobilizerUIndex(1), 2);
    getPencil().setOneU(s, MobilizerUIndex(2), -2);

    getBall().setQToFitTranslation(s, Vec3(-3,3,1));
    getBall().setUToFitAngularVelocity(s, Vec3(4,5,-20));
    Assembler(*this).setErrorTolerance(1e-6).assemble(s);

}


//==============================================================================
//                                 BLOCK
//==============================================================================
Block::Block() {
    // Abbreviations.
    SimbodyMatterSubsystem&     matter = updMatterSubsystem();
    GeneralForceSubsystem&      forces = updForceSubsystem();
    MobilizedBody&              Ground = matter.updGround();

    // Build the multibody system.
    m_gravity = Force::Gravity(forces, matter, -YAxis, 50);
    //m_damper  = Force::GlobalDamper(forces, matter, .1);


    // Control gains
    const Real Kp1 = 50000; // link1-2 joint stiffness
    const Real Kp2 = 10000; // link2-3 joint stiffness
    const Real Cd1 = /*30*/100;    // link1-2 joint damping
    const Real Cd2 = 30;    // link2-3 joint damping

    // Target angles
    const Real Target1 = 0;
    const Real Target2 = Pi/4;

    const Vec3 Cube(.5,.5,.5); // half-dimensions of cube
    const Real Mass1=100, Mass2=20, Mass3=1;
    const Vec3 Centroid(.5,.5,0), Offset(0,0,-.5);
    const Vec3 COM1=Centroid, COM2=Centroid, COM3=Centroid+0*Offset;
    const Inertia Inertia1=Inertia(Vec3(1))  .shiftFromMassCenter(-COM1,Mass1);
    const Inertia Inertia2=Inertia(Vec3(.05)).shiftFromMassCenter(-COM2,Mass2);
    const Inertia Inertia3=Inertia(Vec3(.001,.001,0))
                                             .shiftFromMassCenter(-COM3,Mass3);


    const Real CoefRest = 0;
    const Real CaptureVelocity = .001;
    const Real TransitionVelocity = .01;
    const Real mu_s = 0.15, mu_d = 0.1, mu_v = 0;
    //const Real mu_s = 1, mu_d = 1, mu_v = 0;

    setDefaultLengthScale(Cube.norm());

    m_captureVelocity = CaptureVelocity;
    m_transitionVelocity = TransitionVelocity;

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS
    Ground.addBodyDecoration(Vec3(0,.05,0), DecorativeFrame(2).setColor(Green));
    DecorativeBrick drawCube(Cube); drawCube.setOpacity(0.5).setColor(Gray);

    Body::Rigid blockBody(MassProperties(Mass1, COM1, Inertia1));
    blockBody.addDecoration(Centroid,drawCube);
    blockBody.addDecoration(COM1, DecorativePoint());

    m_block = MobilizedBody::Free(Ground, Vec3(0), blockBody, Vec3(0));

    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Centroid + Vec3(i,j,k).elementwiseMultiply(Cube);
        PointPlaneContact* contact = new PointPlaneContact
           (Ground, YAxis, 0., m_block, pt, CoefRest, mu_s, mu_d, mu_v);
        matter.adoptUnilateralContact(contact);
    }

    Body::Rigid link2Info(MassProperties(Mass2, COM2, Inertia2));
    link2Info.addDecoration(Centroid, drawCube);

    Body::Rigid link3Info(MassProperties(Mass3, COM3, Inertia3));
    link3Info.addDecoration(Centroid, drawCube);
    //m_link2 = MobilizedBody::Pin(m_block, Vec3(2*Centroid),
    //                                link2Info, Vec3(0));
    //m_link3 = MobilizedBody::Pin(m_link2, Vec3(2*Centroid),
    //                                link3Info, Vec3(0));

    Force::ConstantTorque(forces, m_block, Vec3(0,0,-2000));

    // It is more stable to build the springs into the mechanism like
    // this rather than apply them discretely.
    //Force::MobilityLinearSpring
    //   (forces, m_link2, MobilizerQIndex(0), Kp1/10, Target1);
    //Force::MobilityLinearSpring
     //  (forces, m_link3, MobilizerQIndex(0), Kp2, Target2);
    //Force::MobilityLinearDamper
        //(forces, m_link2, MobilizerUIndex(0), Cd1);
   // Force::MobilityLinearDamper
     //   (forces, m_link3, MobilizerUIndex(0), Cd2);
}

void Block::calcInitialState(State& s) const {
    s = realizeTopology(); // returns a reference to the the default state
    realizeModel(s); // define appropriate states for this System
    realize(s, Stage::Instance); // instantiate constraints if any
    realize(s, Stage::Position);
    Assembler(*this).setErrorTolerance(1e-6).assemble(s);

    //m_link2.setQ(s, -0.1104);
   // m_link3.setQ(s, Pi/4);
    //m_link2.lock(s); m_link3.lock(s);

    //getBlock().setQToFitTranslation(s, Vec3(0,.2,0));
    //getBlock().setQToFitRotation(s, Rotation(Pi/20,ZAxis));
    getBlock().setUToFitLinearVelocity(s, Vec3(.2,0,0));

}



//==============================================================================
//                                 EDGES
//==============================================================================
Edges::Edges() {
    // Abbreviations.
    SimbodyMatterSubsystem&     matter = updMatterSubsystem();
    GeneralForceSubsystem&      forces = updForceSubsystem();
    MobilizedBody&              Ground = matter.updGround();

    // Build the multibody system.
    m_gravity = Force::Gravity(forces, matter, -YAxis, 50);
    //m_damper  = Force::GlobalDamper(forces, matter, .1);

    const Vec3 Grounded(2,2,2); // half-dimensions of grounded block
    const Vec3 Swinging(2,3,1); // half-dimensions of swinging block
    const Real MassGrounded=100, MassSwinging=35;

    const Real CoefRest = .7;
    const Real eCoefRest = .95;
    const Real CaptureVelocity = .001;
    const Real TransitionVelocity = .01;
    const Real mu_s = 1.5*0.15, mu_d = 1.5*0.1, mu_v = 0;
    const Real emu_s = .5, emu_d = .3, emu_v = 0;

    setDefaultLengthScale(Grounded.norm());

    m_captureVelocity = CaptureVelocity;
    m_transitionVelocity = TransitionVelocity;

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS
    Ground.addBodyDecoration(Vec3(0,.05,0), DecorativeFrame(2).setColor(Green));
    DecorativeBrick drawGrounded(Grounded);
    drawGrounded.setOpacity(0.5).setColor(Gray);

    Body::Rigid groundedBody(MassProperties(MassGrounded, Vec3(0),
                                            UnitInertia::brick(Grounded)));
    groundedBody.addDecoration(Vec3(0),drawGrounded);

    m_grounded = MobilizedBody::Free(Ground, Vec3(0), groundedBody, Vec3(0));

    DecorativeBrick drawSwinging(Swinging);
    drawSwinging.setOpacity(0.5).setColor(Cyan);

    Body::Rigid swingingBody(MassProperties(MassSwinging, Vec3(0),
                                            UnitInertia::brick(Swinging)));
    swingingBody.addDecoration(Vec3(0),drawSwinging);

    Rotation xy2xz(Pi/2,XAxis);
    m_swinger = MobilizedBody::Universal(Ground,
                                         Transform(xy2xz,Vec3(-2,7,0)),
                                         swingingBody,
                                         Transform(xy2xz,Vec3(-2,3,0.5)));

    //HardStopUpper* hsu= new HardStopUpper(m_swinger, MobilizerQIndex(2),
    //                                      .2, CoefRest);
    //HardStopLower* hsl= new HardStopLower(m_swinger, MobilizerQIndex(2),
    //                                      -.2, CoefRest);
    //matter.adoptUnilateralContact(hsu);
    //matter.adoptUnilateralContact(hsl);


    Rotation gedge(UnitVec3(ZAxis),XAxis,
                   Vec3(1,1,0),ZAxis);
    Rotation sedge1(UnitVec3(YAxis),XAxis,
                   Vec3(-1,0,1),ZAxis);
    Rotation sedge2(UnitVec3(YAxis),XAxis,
                   Vec3(-1,0,-1),ZAxis);
    Rotation sedge3(UnitVec3(YAxis),XAxis,
                   Vec3(1,0,-1),ZAxis);
    Rotation sedge4(UnitVec3(YAxis),XAxis,
                   Vec3(1,0,1),ZAxis);
    EdgeEdgeContact* ee1 = new EdgeEdgeContact
        (m_grounded, Transform(gedge,Vec3(Grounded[0],Grounded[1],0)),
                               Grounded[2],
         m_swinger, Transform(sedge1,Vec3(-Swinging[0],0,Swinging[2])),
                              Swinging[1],
                              eCoefRest, emu_s, emu_d, emu_v);
    matter.adoptUnilateralContact(ee1);

    EdgeEdgeContact* ee2 = new EdgeEdgeContact
        (m_grounded, Transform(gedge,Vec3(Grounded[0],Grounded[1],0)),
                               Grounded[2],
         m_swinger, Transform(sedge2,Vec3(-Swinging[0],0,-Swinging[2])),
                              Swinging[1],
                              eCoefRest, emu_s, emu_d, emu_v);
    matter.adoptUnilateralContact(ee2);


    //EdgeEdgeContact* ee3 = new EdgeEdgeContact
    //    (m_grounded, Transform(gedge,Vec3(Grounded[0],Grounded[1],0)),
    //                           Grounded[2],
    //     m_swinger, Transform(sedge3,Vec3(Swinging[0],0,-Swinging[2])),
    //                          Swinging[1],
    //                          eCoefRest, emu_s, emu_d, emu_v);
    //matter.adoptUnilateralContact(ee3);

    //EdgeEdgeContact* ee4 = new EdgeEdgeContact
    //    (m_grounded, Transform(gedge,Vec3(Grounded[0],Grounded[1],0)),
    //                           Grounded[2],
    //     m_swinger, Transform(sedge4,Vec3(Swinging[0],0,Swinging[2])),
    //                          Swinging[1],
    //                          eCoefRest, emu_s, emu_d, emu_v);
    //matter.adoptUnilateralContact(ee4);

    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(Grounded);
        PointPlaneContact* contact = new PointPlaneContact
           (Ground, YAxis, 0., m_grounded, pt, CoefRest, mu_s, mu_d, mu_v);
        matter.adoptUnilateralContact(contact);
    }

}

void Edges::calcInitialState(State& s) const {
    s = realizeTopology(); // returns a reference to the the default state
    realizeModel(s); // define appropriate states for this System
    realize(s, Stage::Instance); // instantiate constraints if any
    realize(s, Stage::Position);
    Assembler(*this).setErrorTolerance(1e-6).assemble(s);

    getGroundBlock().setQToFitTranslation(s, Vec3(0,2,0));
    getGroundBlock().setQToFitRotation(s, Rotation(Pi/10,YAxis));

    getSwingingBlock().setQToFitRotation(s, Rotation(BodyRotationSequence,
                                                     Pi/3, YAxis,
                                                     0*Pi/8, XAxis));

    getSwingingBlock().setUToFitAngularVelocity(s, Vec3(0,0,-2));


}


//-------------------------- SHOW CONSTRAINT STATUS ----------------------------
//void MyUnilateralConstraintSet::
//showConstraintStatus(const State& s, const String& place) const
//{
//#ifndef NDEBUG
//    printf("\n%s: uni status @t=%.15g\n", place.c_str(), s.getTime());
//    m_mbs.realize(s, Stage::Dynamics);
//    for (int i=0; i < getNumContactElements(); ++i) {
//        const MyContactElement& contact = getContactElement(i);
//        const bool isActive = !contact.isDisabled(s);
//        printf("  %6s %2d %s, h=%g dh=%g f=%g\n",
//                isActive?"ACTIVE":"off", i, contact.getContactType().c_str(),
//                contact.getPerr(s),contact.getVerr(s),
//                isActive?contact.getForce(s):Zero);
//    }
//    for (int i=0; i < getNumFrictionElements(); ++i) {
//        const MyFrictionElement& friction = getFrictionElement(i);
//        if (!friction.isMasterActive(s))
//            continue;
//        const bool isEnabled = friction.isEnabled(s);
//        printf("  %8s friction %2d\n",
//                isEnabled?"STICKING":"sliding", i);
//        friction.writeFrictionInfo(s, "    ", std::cout);
//    }
//    printf("\n");
//#endif
//}
