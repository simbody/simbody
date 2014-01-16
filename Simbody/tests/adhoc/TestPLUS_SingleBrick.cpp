/* -------------------------------------------------------------------------- *
 *          Simbody(tm) - Test PLUS impact model with a single brick          *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013-2014 Stanford University and the Authors.      *
 * Authors: Thomas Uchida                                                     *
 * Contributors: Michael Sherman                                              *
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

/* Test of Poisson-Lankarani-Uchida-Sherman (PLUS) impact model. A brick falls
under the force of gravity and spheres attached to its vertices collide with the
horizontal ground plane. Exhaustive search and successive pruning strategies are
uesd by the position projection and impact handlers to determine suitable active
sets. In this simple test, repeated impacts are used in place of an explicit
contact handler. */

#include "Simbody.h"
using namespace SimTK;
using std::cout;
using std::endl;

// Unique index types to avoid confusing the brick vertex indices (0..7) with
// indices of proximal points (of which there will be 4 or fewer).
SimTK_DEFINE_UNIQUE_INDEX_TYPE(BrickVertexIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ProximalPointIndex);


//==============================================================================
//                                  PARAMETERS
//==============================================================================
const bool SetVizModeToRealtime = false;  //Enable when saving movie.

const bool PauseBeforeClosingVizWindow      = true;
const bool PauseBeforeAndAfterProjection    = false;
const bool PauseBeforeAndAfterImpact        = false;
const bool PauseAfterEachImpactInterval     = false;
const bool PauseAfterEachActiveSetCandidate = false;

const bool PrintBasicInfo                = true;
const bool PrintSystemEnergy             = false;
const bool PrintDebugInfoProjection      = false;
const bool PrintDebugInfoImpact          = false;
const bool PrintDebugInfoActiveSetSearch = false;
const bool PrintDebugInfoStepLength      = false;

const bool ExhaustiveSearchProjection = false;
const bool ExhaustiveSearchImpact     = false;

const bool RunTestsOnePoint        = true;
const bool RunTestsTwoPoints       = true;
const bool RunTestsFourPoints      = true;
const bool RunTestsLongSimulations = true;

// Notes:
// 1. SlidingDirStepLength must be sufficiently small such that, if a sliding
//    interval advancing by this amount causes a sliding direction reversal,
//    then |v_t| < MaxStickingTangVel (and the point is actually rolling).
// 2. StepLengthConvFactor should be sufficiently small to avoid reaching
//    MaxIterStepLength, yet sufficiently large to avoid creating additional
//    intervals.

const Real TolProjectQ           = 1.0e-6;  //Accuracy used by projectQ().
const Real TolPositionFuzziness  = 1.0e-4;  //Expected position tolerance.
const Real TolVelocityFuzziness  = 1.0e-5;  //Expected velocity tolerance.
const Real TolReliableDirection  = 1.0e-4;  //Whether to trust v_t direction.
const Real TolMaxDifDirIteration = 0.05;    //Slip direction within 2.9 degrees.
const Real TolSlidingDirChange   = 0.05;    //Tolerance for MaxSlidingDirChange.
const Real MinMeaningfulImpulse  = 1.0e-6;  //Smallest acceptable impulse.
const Real MaxStickingTangVel    = 1.0e-1;  //Cannot stick above this velocity.
const Real MaxSlidingDirChange   = 0.5;     //Direction can change 28.6 degrees.
const Real SlidingDirStepLength  = 1.0e-4;  //Used to determine slip directions.
const int  MaxIterSlipDirection  = 10;      //Iteration limit for directions.
const int  MaxIterStepLength     = 20;      //Iteration limit for step length.

const Real IntegAccuracy = 1.0e-8;
const Real MaxStepSize   = 1.0e-3;
const int  DesiredFPS    = 30;
const int  DrawEveryN    = int(1.0/DesiredFPS/MaxStepSize + 0.5);
const Vec3 BrickColor    = Blue;
const Vec3 SphereColor   = Red;


//==============================================================================
//                            FREE UNILATERAL BRICK
//==============================================================================
// Establish a free brick with a unilaterally-constrained sphere attached to
// each of its vertices. Each sphere can impact the horizontal ground plane at
// Z=0. All spheres are assumed to have the same radius and material properties.
class FreeUnilateralBrick : public MobilizedBody::Free {
public:

    //--------------------------------------------------------------------------
    // Constructors
    //--------------------------------------------------------------------------
    FreeUnilateralBrick() : Mobod::Free() {}
    FreeUnilateralBrick(Mobod& parent, const Transform& X_PF,
                        const Body& bodyInfo, const Transform& X_BM,
                        const Vec3& brickHalfLengths, Real sphereRadii,
                        Real muDyn, Real vMinRebound, Real vPlasticDeform,
                        Real minCOR);

    //--------------------------------------------------------------------------
    // Getters
    //--------------------------------------------------------------------------
    Real get_muDyn() const           {return m_muDyn;}
    Real get_vMinRebound() const     {return m_vMinRebound;}
    Real get_vPlasticDeform() const  {return m_vPlasticDeform;}
    Real get_minCOR() const          {return m_minCOR;}
    int  get_numVertices() const     {return (int)m_vertices.size();}

    //--------------------------------------------------------------------------
    // Temporary point-in-plane constraints
    //--------------------------------------------------------------------------
    Real get_pipConstraintHeightInGround(const State& s,
                                         const BrickVertexIndex i) const;
    void set_pipConstraintLocation(const State& s, const BrickVertexIndex i,
                                   const Vec3& positionInGround);
    void enablePipConstraint(State& s, const BrickVertexIndex i) const;
    void disableAllPipConstraints(State& s) const;

    //--------------------------------------------------------------------------
    // Temporary ball constraints
    //--------------------------------------------------------------------------
    MultiplierIndex get_ballConstraintFirstIndex(const State& s,
        const BrickVertexIndex i) const;
    void set_ballConstraintLocation(const State& s, const BrickVertexIndex i,
                                    const Vec3& positionInGround);
    void enableBallConstraint(State& s, const BrickVertexIndex i) const;
    void disableAllBallConstraints(State& s) const;

    //--------------------------------------------------------------------------
    // Position-level information
    //--------------------------------------------------------------------------
    // Find the location of the lowest point on the ith sphere, measured from
    // the ground origin and resolved in the ground frame.
    Vec3 findLowestPointLocationInGround(const State& s,
                                         const BrickVertexIndex i) const;

    // Find the location of the lowest point on the ith sphere, measured from
    // the body's origin and resolved in the body's frame.
    Vec3 findLowestPointLocationInBodyFrame(const State& s,
                                            const BrickVertexIndex i) const;

    // Assemble an array containing the location of the lowest point on each
    // sphere, measured from the ground origin and resolved in the ground frame.
    void findAllLowestPointLocationsInGround(const State& s,
        Array_<Vec3,BrickVertexIndex>& lowestPointLocationsInG) const;

    // Check for interpenetration of the brick with the ground plane.
    bool isBrickInterpenetrating(const Array_<Vec3,BrickVertexIndex>&
                                       lowestPointLocationsInG) const;
    bool isBrickInterpenetrating(const State& s) const;

    // Check for proximity of a sphere to the ground plane.
    bool isPointProximal(const Vec3& lowestPointLocationInG) const;

    // Assemble an array containing the indices of the proximal points.
    void findProximalPointIndices(
        const Array_<Vec3,BrickVertexIndex>& lowestPointLocationsInG,
        Array_<BrickVertexIndex,ProximalPointIndex>& proximalPointIndices)
        const;

    //--------------------------------------------------------------------------
    // Velocity-level information
    //--------------------------------------------------------------------------
    // Find the velocity of the lowest point on the ith sphere resolved in the
    // ground frame.
    Vec3 findLowestPointVelocityInGround(const State& s,
                                         const BrickVertexIndex i) const;

    // Find the angle between the global X-axis and the tangential velocity
    // vector, in [-pi, pi]. Returns NaN if the magnitude of the tangential
    // velocity is too small to provide a reliable direction.
    Real findTangentialVelocityAngle(const Vec3& vel) const;

    // Check for impact of the brick with the ground plane.
    bool isBrickImpacting(const Array_<Vec3,ProximalPointIndex>&
                                proximalPointVelocities) const;

private:
    Vec3                           m_brickHalfLengths;
    Real                           m_sphereRadii;
    Real                           m_muDyn;
    Real                           m_vMinRebound;
    Real                           m_vPlasticDeform;
    Real                           m_minCOR;
    Array_<Vec3,BrickVertexIndex>  m_vertices;

    // PointInPlane and Ball constraints are created for each sphere. The
    // locations of the constraints corresponding to the proximal points are
    // adjusted by the PositionProjecter and Impacter constructors.
    Array_<Constraint::PointInPlane,BrickVertexIndex>  m_pipConstraints;
    Array_<Constraint::Ball,BrickVertexIndex>          m_ballConstraints;
};


//==============================================================================
//                              POSITION PROJECTER
//==============================================================================
class PositionProjecter {
public:

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    PositionProjecter(MultibodySystem& mbs, FreeUnilateralBrick& brick,
                      const State& s0,
                      const Array_<Vec3,BrickVertexIndex>& positionsInG);

    //--------------------------------------------------------------------------
    // Resolve position-level violations
    //--------------------------------------------------------------------------
    // Try projecting all combinations of proximal points; select the projection
    // that resolves all violations while requiring the smallest change in Q.
    void projectPositionsExhaustive(State& s) const;

    // Begin by projecting using the constraints associated with all proximal
    // points; successively prune the constraint associated with the most
    // distant proximal point until the projection is successful.
    void projectPositionsPruning(State& s) const;

private:
    MultibodySystem&                             m_mbs;
    FreeUnilateralBrick&                         m_brick;
    Array_<BrickVertexIndex,ProximalPointIndex>  m_proximalPointIndices;

    // Assemble all combinations of 1 or more indices in [0, sizeOfSet]. Creates
    // an array of 2^sizeOfSet-1 index arrays (indices into proximal positions).
    void createProximalIndexArrays(int sizeOfSet,
        Array_<Array_<ProximalPointIndex> >& arrayOfIndexArrays) const;

    // Try projecting positions using the provided combination of constraint
    // indices. Return the 2-norm distance between the original and final Q (or
    // SimTK::Infinity if projection was unsuccessful).
    Real evaluateProjection(State& sTemp, const State& sOrig,
                            const Array_<ProximalPointIndex>& indexArray) const;

    // Find and display new proximal points. For debugging only.
    void displayNewProximalPoints(State& s) const;
};


//==============================================================================
//                                   IMPACTER
//==============================================================================
class Impacter {
public:

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    Impacter(MultibodySystem& mbs, FreeUnilateralBrick& brick, const State& s0,
             const Array_<Vec3,BrickVertexIndex>& allPositionsInG,
             const Array_<BrickVertexIndex,ProximalPointIndex>&
                   proximalPointIndices);

    //--------------------------------------------------------------------------
    // Perform one complete impact
    //--------------------------------------------------------------------------
    // Evaluate all active set candidates for each impact interval, selecting
    // the most fit candidate for each interval. Returns energy dissipated.
    Real performImpactExhaustive(State& s,
        Array_<Vec3,ProximalPointIndex>& proximalVelsInG,
        Array_<bool,ProximalPointIndex>& hasRebounded) const;

    // Begin by impacting using the constraints associated with all proximal
    // points; successively prune the constraint associated with the worst
    // violation until either no violations occur or the active set is empty (in
    // which case the least objectionable active set candidate is retained).
    // Returns energy dissipated.
    Real performImpactPruning(State& s,
                              Array_<Vec3,ProximalPointIndex>& proximalVelsInG,
                              Array_<bool,ProximalPointIndex>& hasRebounded)
                              const;

private:
    enum ImpactPhase      {Compression, Restitution};
    enum TangentialState  {Observing=0, Sliding, Rolling};
    enum SolutionCategory {
        SolCat_NoViolations=0,                        //Ideal
        SolCat_ActiveConstraintDoesNothing,           //Non-minimal active set
        SolCat_RestitutionImpulsesIgnored,            //Simultaneity is lost
        SolCat_TangentialVelocityTooLargeToStick,     //Violates friction rule
        SolCat_StickingImpulseExceedsStictionLimit,   //Violates friction limit
        SolCat_GroundAppliesAttractiveImpulse,        //Violates physical law
        SolCat_NegativePostCompressionNormalVelocity, //Fatal: solution is wrong
        SolCat_NoImpulsesApplied,                     //Fatal: no progress made
        SolCat_UnableToResolveUnknownSlipDirection,   //Fatal: direction unknown
        SolCat_MinStepCausesSlipDirectionReversal,    //Fatal: direction unknown
        SolCat_NotEvaluated,
        SolCat_WorstCatWithNoSeriousViolations        //These categories can be
            = SolCat_RestitutionImpulsesIgnored,      //accepted without issue.
        SolCat_WorstCatWithSomeUsableSolution         //Use these categories
            = SolCat_GroundAppliesAttractiveImpulse   //only if no other choice.
    };

    struct ActiveSetCandidate {
        Array_<int,ProximalPointIndex>  tangentialStates;
        Vector                          systemVelocityChange;
        Vector                          localImpulses;
        SolutionCategory                solutionCategory;
        Real                            fitness;
        ProximalPointIndex              worstConstraint;
    };

    //--------------------------------------------------------------------------
    // Private methods: display
    //--------------------------------------------------------------------------
    // Display active set candidate in human-readable form.
    void printFormattedActiveSet(std::ostream& stream,
        const Array_<int,ProximalPointIndex>& tangentialStates,
        const std::string& prefix = "") const;

    // Map an enumerated solution category to a descriptive string.
    const char* getSolutionCategoryDescription(SolutionCategory solCat) const;

    // Display information about the evaluation of an active set candidate.
    void printActiveSetInfo(const ActiveSetCandidate& asc) const;

    //--------------------------------------------------------------------------
    // Private methods: helper functions
    //--------------------------------------------------------------------------
    // Add addThis to vec in base N; return false on overflow.
    bool addInBaseN(int N, Array_<int,ProximalPointIndex>& vec, int addThis)
        const;

    // Create 3^n-1 arrays of size n, where n is the number of proximal points,
    // and each element corresponds to one of the three TangentialStates.
    void initializeActiveSetCandidateArray(
        Array_<ActiveSetCandidate>& activeSetCandidates) const;

    // Initialize active set candidate for pruning search. Begin with maximal
    // active set, but disallow sticking at points whose tangential velocity
    // magnitude is too large.
    void initializeActiveSetCandidate(const State& s, ActiveSetCandidate& asc)
        const;

    // Clear data associated with active set candidate; called before each
    // active set candidate is evaluated in pruning algorithm.
    void resetActiveSetCandidate(ActiveSetCandidate& asc) const;

    // Clear data associated with active set candidates; called before each
    // impact interval begins in exhaustive search algorithm.
    void resetActiveSetCandidateArray(
        Array_<ActiveSetCandidate>& activeSetCandidates) const;

    // Return row index associated with first component of ith constraint.
    int getIndexOfFirstMultiplier(const State& s,
                                  const ProximalPointIndex i) const;

    // Find the absolute difference between two angles in [-pi, pi] radians,
    // accounting for multiples of 2*pi; returns a value in [0, pi].
    Real calcAbsDiffBetweenAngles(Real a, Real b) const;

    // Given point P and line segment AB, find the point closest to P that lies
    // on AB, which we call Q. Returns stepLength, the ratio AQ:AB. In our case,
    // P is the origin and AB is the line segment connecting the initial and
    // final tangential velocity vectors. Helps calculateIntervalStepLength().
    Real calcSlidingStepLengthToOrigin(const Vec2& A, const Vec2& B) const;

    // Determine whether active set candidate is empty.
    bool isActiveSetCandidateEmpty(const ActiveSetCandidate& asc) const;

    //--------------------------------------------------------------------------
    // Private methods: calculators
    //--------------------------------------------------------------------------
    // Calculate coefficient of restitution.
    Real calcCOR(Real vNormal) const;

    // Generate and solve a linear system of equations to determine the system
    // velocity changes and impulses; assign to ActiveSetCandidate. Resolves
    // unknown sliding directions.
    void generateAndSolveLinearSystem(const State& s0,
                                      const ImpactPhase impactPhase,
                                      const Vector& restitutionImpulses,
                                      ActiveSetCandidate& asc) const;

    // Determine category of linear system solution and calculate fitness value
    // for active set candidate (if it has not already been categorized). Assign
    // worst applicable disqualification category to ActiveSetCandidate.
    void evaluateLinearSystemSolution(const State& s,
                                      const ImpactPhase impactPhase,
                                      const Vector& restitutionImpulses,
                                      ActiveSetCandidate& asc) const;

    // Determine the step length for the selected active set. The upper bound is
    // limited by the minimum number of intervals required and the maximum
    // permissible change in the sliding direction.
    Real calculateIntervalStepLength(const State& s0,
                                     const Array_<Vec3>& currVels,
                                     const ActiveSetCandidate& asc,
                                     int intervalCtr) const;

    //--------------------------------------------------------------------------
    // Member variables
    //--------------------------------------------------------------------------
    MultibodySystem&                             m_mbs;
    FreeUnilateralBrick&                         m_brick;
    Array_<BrickVertexIndex,ProximalPointIndex>  m_proximalPointIndices;
};


//==============================================================================
//                             PAUSABLE VISUALIZER
//==============================================================================
class PausableVisualizer : public Visualizer {
public:
    PausableVisualizer(const MultibodySystem& system) : Visualizer(system) {}

    void reportAndPause(const State& s, const std::string& msg) const
    {
        report(s);
        printf("t=%0.3f  %s", s.getTime(), msg.c_str());
        char trash = getchar();
    }
};


//==============================================================================
//                                 BODY WATCHER
//==============================================================================
// Prior to rendering each frame, point the camera at the given body's origin.
// Adapted from TimsBox.cpp.
class BodyWatcher : public Visualizer::FrameController {
public:
    explicit BodyWatcher(const MobilizedBody& body) : m_body(body) {}

    void generateControls(const Visualizer&             viz,
                          const State&                  state,
                          Array_< DecorativeGeometry >& geometry) OVERRIDE_11
    {
        const Vec3 Bo = m_body.getBodyOriginLocation(state);
        const Vec3 p_GC = Bo + Vec3(0,4,2);
        const Rotation R1(-SimTK::Pi/3, XAxis);
        const Rotation R2(SimTK::Pi, ZAxis);
        viz.setCameraTransform(Transform(R1*R2, p_GC));
    }
private:
    const MobilizedBody m_body;
};


//==============================================================================
//                      Function: CREATE MULTIBODY SYSTEM
//==============================================================================
static void createMultibodySystem(MultibodySystem& mbs,
                                  FreeUnilateralBrick& brick,
                                  Real muDyn, Real minCOR)
{
    // Set up multibody system.
    SimbodyMatterSubsystem matter(mbs);
    mbs.setUpDirection(ZAxis);

    // Physical parameters.
    const Vec3 brickHalfLengths(0.2, 0.3, 0.4);
    const Real sphereRadii    = 0.1;
    const Real brickMass      = 2.0;
    const Real vMinRebound    = 1.0e-6;
    const Real vPlasticDeform = 0.1;

    // Configure brick.
    const Inertia brickInertia(brickMass*UnitInertia::brick(brickHalfLengths));
    Body::Rigid brickInfo(MassProperties(brickMass, Vec3(0), brickInertia));
    DecorativeBrick brickGeom(brickHalfLengths);
    brickGeom.setColor(BrickColor);
    brickInfo.addDecoration(brickGeom);

    // Add brick to multibody system.
    brick = FreeUnilateralBrick(mbs.updMatterSubsystem().Ground(), Vec3(0),
                                brickInfo, Vec3(0), brickHalfLengths,
                                sphereRadii, muDyn, vMinRebound, vPlasticDeform,
                                minCOR);
}


//==============================================================================
//                       Function: PRINT HORIZONTAL RULE
//==============================================================================
static void printHorizontalRule(int spacesTop, int spacesBottom,
                                char ruleCharacter, const std::string& msg = "")
{
    for (int i=0; i<spacesTop; ++i) printf("\n");
    for (int i=0; i<60; ++i) cout << ruleCharacter;
    cout << " " << msg << endl;
    for (int i=0; i<spacesBottom; ++i) printf("\n");
}


//==============================================================================
//                     Function: SIMULATE MULTIBODY SYSTEM
//==============================================================================
static void simulateMultibodySystem(const std::string& description,
                                    const Vector& initialQ,
                                    const Vector& initialU,
                                    Real simDuration, Real muDyn, Real minCOR)
{
    printHorizontalRule(3,0,'=');
    cout << description << endl;
    printHorizontalRule(0,1,'=');

    // Create the multibody system.
    MultibodySystem mbs;
    FreeUnilateralBrick brick;
    createMultibodySystem(mbs, brick, muDyn, minCOR);

    // Add gravity.
    GeneralForceSubsystem forces(mbs);
    Force::Gravity gravity(forces, mbs.updMatterSubsystem(), -ZAxis, 9.81);

    // Set up the visualizer.
    PausableVisualizer viz(mbs);
    viz.setShowSimTime(true).setShowFrameRate(true);
    if (SetVizModeToRealtime)
        viz.setMode(Visualizer::RealTime);
    viz.addFrameController(new BodyWatcher(brick));
    mbs.updMatterSubsystem().setShowDefaultGeometry(false);

    // Initialize.
    State s0 = mbs.realizeTopology();
    mbs.realize(s0, Stage::Dynamics);
    s0.setQ(initialQ);
    s0.setU(initialU);
    mbs.realize(s0, Stage::Dynamics);
    const Real initialEnergy = mbs.calcEnergy(s0);

    // Set up integrator.
    RungeKutta3Integrator integ(mbs);
    integ.setAccuracy(IntegAccuracy);
    integ.setAllowInterpolation(false);
    integ.initialize(s0);
    const int totalSteps = int(simDuration/MaxStepSize + 0.5);

    // Simulate.
    printf("Simulating for %0.1f seconds (%d steps of size %0.3f)\n",
           simDuration, totalSteps, MaxStepSize);
    viz.reportAndPause(s0, "Press <Enter> to simulate...");

    for (int stepNum=1; stepNum<totalSteps; ++stepNum) {

        //----------------------------- INTEGRATE ------------------------------
        // Advance time by MaxStepSize. Might require multiple internal steps.
        const Real tNext = stepNum*MaxStepSize;
        do {integ.stepBy(MaxStepSize);} while (integ.getTime() < tNext);

        // The state being used by the integrator.
        State& s = integ.updAdvancedState();

        //----------------- RESOLVE POSITION-LEVEL VIOLATIONS ------------------
        // Project positions to resolve interpenetration. The PositionProjecter
        // guarantees that no points are below -TolPositionFuzziness on exit.
        mbs.realize(s, Stage::Position);

        // Calculate the position of the lowest point on each sphere.
        Array_<Vec3,BrickVertexIndex> preProjectionPos;
        brick.findAllLowestPointLocationsInGround(s, preProjectionPos);

        // Can impact only if interpenetration occurred.
        bool interpenetrationDetected = false;

        if (brick.isBrickInterpenetrating(preProjectionPos)) {
            interpenetrationDetected = true;

            if (PrintSystemEnergy) {
                mbs.realize(s, Stage::Dynamics);
                cout << "  [normalized energy before project] "
                     << mbs.calcEnergy(s)/initialEnergy << endl;
            }
            if (PauseBeforeAndAfterProjection) {
                cout << "  [pos0] " << s.getQ() << endl;
                viz.reportAndPause(s, "Ready to project positions");
            }

            PositionProjecter positionProjecter(mbs,brick,s,preProjectionPos);
            if (ExhaustiveSearchProjection)
                positionProjecter.projectPositionsExhaustive(s);
            else
                positionProjecter.projectPositionsPruning(s);

            if (PrintSystemEnergy) {
                mbs.realize(s, Stage::Dynamics);
                cout << "  [normalized energy  after project] "
                     << mbs.calcEnergy(s)/initialEnergy << endl;
            }
            if (PauseBeforeAndAfterProjection) {
                cout << "  [pos1] " << s.getQ() << endl;
                viz.reportAndPause(s, "Position projection complete");
            }
        }

        //----------------- RESOLVE VELOCITY-LEVEL VIOLATIONS ------------------
        // Perform impacts to resolve negative normal velocities of proximal
        // points. The Impacter guarantees that no proximal points have normal
        // velocities less than -TolVelocityFuzziness on exit.

        // The impact handler is triggered only on actual interpenetration
        // (i.e., pos[ZAxis] < 0), not just on entering the fuzzy proximal zone
        // (i.e., pos[ZAxis] < TolPositionFuzziness). This way, other nearly-
        // penetrating points will have the chance to move into the proximal
        // zone and be included in the impact event, which will help preserve
        // symmetry and simultaneity.
        if (interpenetrationDetected) {
            mbs.realize(s, Stage::Velocity);

            // Calculate all positions after projection.
            Array_<Vec3,BrickVertexIndex> allPositionsInG;
            brick.findAllLowestPointLocationsInGround(s, allPositionsInG);

            // Find the indices of the proximal points.
            Array_<BrickVertexIndex,ProximalPointIndex> proximalPointIndices;
            brick.findProximalPointIndices(allPositionsInG,
                                           proximalPointIndices);

            // Calculate the velocity of each proximal point.
            Array_<Vec3,ProximalPointIndex>
                proximalVelsInG(proximalPointIndices.size());
            for (ProximalPointIndex i(0);
                 i<(int)proximalPointIndices.size(); ++i)
                proximalVelsInG[i] = brick.findLowestPointVelocityInGround(s,
                                         proximalPointIndices[i]);

            // An impact event occurs only if at least one negative normal
            // velocity exceeds the velocity tolerance (i.e., vel[ZAxis] <
            // -TolVelocityFuzziness).
            if (brick.isBrickImpacting(proximalVelsInG)) {

                // Record which points have already undergone a restitution
                // phase; set coefficient of restitution to zero for follow-up
                // impacts at these points.
                Array_<bool,ProximalPointIndex>
                    hasRebounded(proximalPointIndices.size());
                for (ProximalPointIndex i(0);
                     i<(int)proximalPointIndices.size(); ++i)
                    hasRebounded[i] = false;

                // Process impacts until all proximal points have non-negative
                // normal velocities.
                while (brick.isBrickImpacting(proximalVelsInG)) {

                    Real normEnergyBefore = SimTK::NaN;
                    if (PrintSystemEnergy) {
                        mbs.realize(s, Stage::Dynamics);
                        normEnergyBefore = mbs.calcEnergy(s)/initialEnergy;
                        cout << "  [normalized energy before impact ] "
                             << normEnergyBefore << endl;
                    }
                    if (PauseBeforeAndAfterImpact) {
                        cout << "  [vel0] " << s.getU() << endl;
                        viz.reportAndPause(s, "Ready to perform impact");
                    }

                    // Perform one complete impact consisting of a compression
                    // phase and (if necessary) a restitution phase.
                    Impacter impacter(mbs, brick, s, allPositionsInG,
                                      proximalPointIndices);
                    Real energyDissipated = SimTK::NaN;
                    if (ExhaustiveSearchImpact)
                        energyDissipated = impacter.performImpactExhaustive(s,
                                               proximalVelsInG, hasRebounded);
                    else
                        energyDissipated = impacter.performImpactPruning(s,
                                               proximalVelsInG, hasRebounded);

                    if (PrintSystemEnergy) {
                        mbs.realize(s, Stage::Dynamics);
                        const Real normEnergyAfter = mbs.calcEnergy(s) /
                                                     initialEnergy;
                        const Real normEnergyDiss  = energyDissipated /
                                                     initialEnergy;
                        cout << "  [normalized energy  after impact ] "
                             << normEnergyAfter << " (" << normEnergyDiss
                             << " dissipated, error = "
                             << std::abs(normEnergyBefore - normEnergyAfter
                                         - normEnergyDiss) << ")" << endl;
                    }
                    if (PauseBeforeAndAfterImpact) {
                        cout << "  [vel1] " << s.getU() << endl;
                        viz.reportAndPause(s, "Impact complete");
                    }

                } //end while processing impacts
            } //end if impacting
        } //end if able to impact

        // Output a frame to the visualizer, if necessary.
        if (stepNum%DrawEveryN == 0 || stepNum == 1 || stepNum == totalSteps-1)
            viz.report(s);

    } //end simulation loop

    printf("Simulation complete.\n");
    if (PauseBeforeClosingVizWindow) {
        printf("Press <Enter> to continue...\n");
        char trash = getchar();
    }
    viz.shutdown();
}


//==============================================================================
//                                     MAIN
//==============================================================================
int main() {
    try {

        // Perform a series of simulations using different initial conditions.

        //---------------------- (a) One point impacting -----------------------
        if (RunTestsOnePoint) {
            Vector initQ = Vector(Vec7(1,0,0,0, 0,1,0.8));
            initQ(0,4) = Vector(Quaternion(Rotation(SimTK::Pi/4, XAxis) *
                                           Rotation(SimTK::Pi/6, YAxis))
                                .asVec4());
            Vector initU = Vector(Vec6(0,0,0, 0,0,6));

            simulateMultibodySystem("Test A1: One point, no v_tangential",
                                    initQ, initU, 1.8, 0.6, 0.5);
            initU[3] = 0.5;
            simulateMultibodySystem("Test A2: One point, small v_tangential",
                                    initQ, initU, 1.8, 0.6, 0.5);
            initU[3] = -0.5;
            simulateMultibodySystem("Test A3: One point, small v_tangential",
                                    initQ, initU, 1.8, 0.6, 0.5);
            initQ[6] = 1.5;
            initU[3] = 5.0;
            initU[5] = -4.0;
            simulateMultibodySystem("Test A4: One point, large v_tangential",
                                    initQ, initU, 1.0, 0.3, 1.0);
            initU[3] = -5.0;
            simulateMultibodySystem("Test A5: One point, large v_tangential",
                                    initQ, initU, 1.0, 0.3, 1.0);
        }

        //---------------------- (b) Two points impacting ----------------------
        if (RunTestsTwoPoints) {
            Vector initQ = Vector(Vec7(1,0,0,0, 0,1,0.8));
            initQ(0,4) = Vector(Quaternion(Rotation(SimTK::Pi/4, XAxis))
                                .asVec4());
            Vector initU = Vector(Vec6(0,0,0, 0,0,6));

            simulateMultibodySystem("Test B1: Two points, no v_tangential",
                                    initQ, initU, 1.8, 0.6, 0.5);
            initU[4] = -1.0;
            simulateMultibodySystem("Test B2: Two points, small v_tangential",
                                    initQ, initU, 1.8, 0.6, 0.5);
        }

        //--------------------- (c) Four points impacting ----------------------
        if (RunTestsFourPoints) {
            Vector initQ = Vector(Vec7(1,0,0,0, 0,1,0.8));
            Vector initU = Vector(Vec6(0,0,0, 0,0,6));

            simulateMultibodySystem("Test C1: Four points, no v_tangential",
                                    initQ, initU, 1.8, 0.6, 0.5);
            initU[3] = 0.5;
            simulateMultibodySystem("Test C2: Four points, small v_tangential",
                                    initQ, initU, 1.8, 0.6, 0.5);
        }

        //------------------------ (d) Long simulation -------------------------
        if (RunTestsLongSimulations) {
            Vector initQ = Vector(Vec7(1,0,0,0, 0,1,1.5));
            initQ(0,4) = Vector(Quaternion(Rotation(SimTK::Pi/4, XAxis) *
                                           Rotation(SimTK::Pi/6, YAxis))
                                .asVec4());
            Vector initU = Vector(Vec6(0,0,0, -5,0,-4));

            simulateMultibodySystem("Test D1: low muDyn, high minCOR",
                                    initQ, initU, 5.5, 0.3, 1.0);
            simulateMultibodySystem("Test D2: low muDyn, low minCOR",
                                    initQ, initU, 3.0, 0.3, 0.5);
            simulateMultibodySystem("Test D3: high muDyn, high minCOR",
                                    initQ, initU, 3.8, 0.8, 1.0);
            simulateMultibodySystem("Test D4: high muDyn, low minCOR",
                                    initQ, initU, 3.5, 0.8, 0.5);
        }

    } catch (const std::exception& e) {
        cout << "ERROR: " << e.what() << endl;
        return 1;
    }
    return 0;
}


//==============================================================================
//                    Implementation: FREE UNILATERAL BRICK
//==============================================================================
FreeUnilateralBrick::
FreeUnilateralBrick(Mobod& parent, const Transform& X_PF, const Body& bodyInfo,
                    const Transform& X_BM, const Vec3& brickHalfLengths,
                    Real sphereRadii, Real muDyn, Real vMinRebound,
                    Real vPlasticDeform, Real minCOR)
    : Mobod::Free(parent, X_PF, bodyInfo, X_BM)
{
    // Ensure parameters are physically reasonable.
    m_brickHalfLengths = brickHalfLengths;
    m_sphereRadii      = std::max(0.0, sphereRadii);
    m_muDyn            = std::max(0.0, muDyn);
    m_vPlasticDeform   = std::max(0.0, vPlasticDeform);
    m_vMinRebound      = clamp(0.0, vMinRebound, m_vPlasticDeform);
    m_minCOR           = clamp(0.0, minCOR, 1.0);

    // Add spheres to display geometry.
    DecorativeSphere sphereGeom(m_sphereRadii);
    sphereGeom.setColor(SphereColor);

    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        Vec3 vertex = Vec3(i,j,k).elementwiseMultiply(m_brickHalfLengths);
        m_vertices.push_back(vertex);
        addBodyDecoration(Transform(vertex), sphereGeom);
    }

    // Create one PointInPlane and one Ball constraint for each sphere.
    Mobod ground = getMatterSubsystem().getGround();
    for (BrickVertexIndex i(0); i<(int)m_vertices.size(); ++i) {
        Constraint::PointInPlane pip(ground, ZAxis, 0, *this, Vec3(0));
        pip.setDisabledByDefault(true);
        m_pipConstraints.push_back(pip);

        Constraint::Ball ball(ground, Vec3(0), *this, Vec3(0));
        ball.setDisabledByDefault(true);
        m_ballConstraints.push_back(ball);
    }
}

//--------------------------- Temporary constraints ----------------------------
Real FreeUnilateralBrick::
get_pipConstraintHeightInGround(const State& s, const BrickVertexIndex i) const
{
    return getMatterSubsystem().getGround().findStationLocationInAnotherBody(s,
        m_pipConstraints[i].getDefaultFollowerPoint(), *this)[ZAxis];
}

void FreeUnilateralBrick::
set_pipConstraintLocation(const State& s, const BrickVertexIndex i,
                          const Vec3& positionInGround)
{
    m_pipConstraints[i].setDefaultFollowerPoint(
        getMatterSubsystem().getGround().
        findStationLocationInAnotherBody(s, positionInGround, *this));
}

void FreeUnilateralBrick::
enablePipConstraint(State& s, const BrickVertexIndex i) const
{   m_pipConstraints[i].enable(s); }

void FreeUnilateralBrick::disableAllPipConstraints(State& s) const
{
    for (BrickVertexIndex i(0); i<m_pipConstraints.size(); ++i)
        m_pipConstraints[i].disable(s);
}

MultiplierIndex FreeUnilateralBrick::
get_ballConstraintFirstIndex(const State& s, const BrickVertexIndex i) const
{
    MultiplierIndex px0,vx0,ax0;
    m_ballConstraints[i].getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    return px0;
}

void FreeUnilateralBrick::
set_ballConstraintLocation(const State& s, const BrickVertexIndex i,
                           const Vec3& positionInGround)
{
    m_ballConstraints[i].setDefaultPointOnBody1(positionInGround);
    m_ballConstraints[i].setDefaultPointOnBody2(
        getMatterSubsystem().getGround().
        findStationLocationInAnotherBody(s, positionInGround, *this));
}

void FreeUnilateralBrick::
enableBallConstraint(State& s, const BrickVertexIndex i) const
{   m_ballConstraints[i].enable(s); }

void FreeUnilateralBrick::disableAllBallConstraints(State& s) const
{
    for (BrickVertexIndex i(0); i<m_ballConstraints.size(); ++i)
        m_ballConstraints[i].disable(s);
}

//------------------------- Position-level information -------------------------
Vec3 FreeUnilateralBrick::
findLowestPointLocationInGround(const State& s, const BrickVertexIndex i) const
{
    return findStationLocationInGround(s, m_vertices[i])
           + Vec3(0,0,-m_sphereRadii);
}

Vec3 FreeUnilateralBrick::
findLowestPointLocationInBodyFrame(const State& s,
                                   const BrickVertexIndex i) const
{
    // This is the vector from the origin of the body to the lowest point on the
    // ith sphere, resolved in the ground frame.
    Vec3 posG = findLowestPointLocationInGround(s,i) - getBodyOriginLocation(s);

    // Transform to body-fixed frame by applying rotation only.
    return expressGroundVectorInBodyFrame(s,posG);
}

void FreeUnilateralBrick::findAllLowestPointLocationsInGround(const State& s,
     Array_<Vec3,BrickVertexIndex>& lowestPointLocationsInG) const
{
    SimTK_ASSERT(lowestPointLocationsInG.empty(), "Input array must be empty.");
    for (BrickVertexIndex i(0); i<(int)m_vertices.size(); ++i)
        lowestPointLocationsInG.push_back(findLowestPointLocationInGround(s,i));
}

bool FreeUnilateralBrick::isBrickInterpenetrating(
    const Array_<Vec3,BrickVertexIndex>& lowestPointLocationsInG) const
{
    for (BrickVertexIndex i(0); i<(int)lowestPointLocationsInG.size(); ++i)
        if (lowestPointLocationsInG[i][ZAxis] < -TolPositionFuzziness)
            return true;
    return false;
}

bool FreeUnilateralBrick::isBrickInterpenetrating(const State& s) const
{
    for (BrickVertexIndex i(0); i<(int)m_vertices.size(); ++i)
        if (findLowestPointLocationInGround(s,i)[ZAxis] < -TolPositionFuzziness)
            return true;
    return false;
}

bool FreeUnilateralBrick::
isPointProximal(const Vec3& lowestPointLocationInG) const
{   return (lowestPointLocationInG[ZAxis] < TolPositionFuzziness); }

void FreeUnilateralBrick::findProximalPointIndices(
    const Array_<Vec3,BrickVertexIndex>& lowestPointLocationsInG,
    Array_<BrickVertexIndex,ProximalPointIndex>& proximalPointIndices) const
{
    SimTK_ASSERT(proximalPointIndices.empty(), "Input array must be empty.");
    for (BrickVertexIndex i(0); i<(int)lowestPointLocationsInG.size(); ++i)
        if (isPointProximal(lowestPointLocationsInG[i]))
            proximalPointIndices.push_back(i);
}

//------------------------- Velocity-level information -------------------------
Vec3 FreeUnilateralBrick::
findLowestPointVelocityInGround(const State& s, const BrickVertexIndex i) const
{
    return findStationVelocityInGround(s,
        findLowestPointLocationInBodyFrame(s,i));
}

Real FreeUnilateralBrick::findTangentialVelocityAngle(const Vec3& vel) const
{
    if (vel.getSubVec<2>(0).norm() < TolReliableDirection)
        return SimTK::NaN;
    return atan2(vel[YAxis], vel[XAxis]);
}

bool FreeUnilateralBrick::isBrickImpacting(
    const Array_<Vec3,ProximalPointIndex>& proximalPointVelocities) const
{
    for (ProximalPointIndex i(0); i<(int)proximalPointVelocities.size(); ++i)
        if (proximalPointVelocities[i][ZAxis] < -TolVelocityFuzziness)
            return true;
    return false;
}


//==============================================================================
//                      Implementation: POSITION PROJECTER
//==============================================================================
PositionProjecter::
PositionProjecter(MultibodySystem& mbs,
                  FreeUnilateralBrick& brick,
                  const State& s0,
                  const Array_<Vec3,BrickVertexIndex>& positionsInG)
    : m_mbs(mbs), m_brick(brick)
{
    // Create an array containing the index of each proximal point.
    m_brick.findProximalPointIndices(positionsInG, m_proximalPointIndices);

    // Adjust the position of the PointInPlane constraints corresponding to the
    // proximal points.
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        m_brick.set_pipConstraintLocation(s0, m_proximalPointIndices[i],
            positionsInG[m_proximalPointIndices[i]]);

    if (PrintDebugInfoProjection) {
        printHorizontalRule(1,0,'*',"projecting positions");
        cout << "  -> " << m_proximalPointIndices.size() << " proximal point(s)"
             << endl;
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
            cout << "     [" << i << "] p="
                 << positionsInG[m_proximalPointIndices[i]] << endl;
    }
}

//--------------------- Resolve position-level violations ----------------------
void PositionProjecter::projectPositionsExhaustive(State& s) const
{
    // Assemble all combinations of 1 or more proximal points.
    Array_<Array_<ProximalPointIndex> > arrayOfIndexArrays;
    createProximalIndexArrays((int)m_proximalPointIndices.size(),
                              arrayOfIndexArrays);
    if (PrintDebugInfoProjection)
        cout << "  -> " << arrayOfIndexArrays.size() << " combination(s)"
             << endl;

    // Try projecting using every combination of constraint indices; compute
    // 2-norm distance between original and final Q.
    Real    minDistance = SimTK::Infinity;
    Vector  minQ;
    int     minIdx;
    for (int comb=0; comb<(int)arrayOfIndexArrays.size(); ++comb) {
        State sTemp = m_mbs.realizeTopology();
        sTemp.setQ(s.getQ());
        Real dist = evaluateProjection(sTemp, s, arrayOfIndexArrays[comb]);

        // Several combinations can have the same distance metric. Favor the
        // combination with the most enabled constraints.
        if (dist < minDistance || (std::abs(dist-minDistance) < TolProjectQ &&
                                   arrayOfIndexArrays[comb].size() >
                                   arrayOfIndexArrays[minIdx].size())) {
            minDistance = dist;
            minQ        = sTemp.getQ();
            minIdx      = comb;
        }

        if (PrintDebugInfoProjection)
            cout << "     [" << std::setw(2) << comb << "] "
                 << "d=" << std::setprecision(6) << std::setw(10) << dist
                 << "  " << arrayOfIndexArrays[comb] << endl;
    }

    SimTK_ASSERT_ALWAYS(minDistance < SimTK::Infinity,
        "No valid position projection found by exhaustive search.");

    // Apply the projection that resolves all violations while requiring the
    // smallest change in Q.
    s.setQ(minQ);

    if (PrintDebugInfoProjection) {
        cout << "  -> Exhaustive search selected index " << minIdx
             << ", constraints " << arrayOfIndexArrays[minIdx] << endl;
        displayNewProximalPoints(s);
        printHorizontalRule(0,1,'*');
    }
}

void PositionProjecter::projectPositionsPruning(State& s) const
{
    // Begin with indices of all proximal points.
    Array_<ProximalPointIndex> indexArray(m_proximalPointIndices.size());
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        indexArray[i] = i;

    if (PrintDebugInfoProjection)
        cout << "  -> Starting pruning search with " << indexArray.size()
             << " constraint(s)" << endl;

    // Realize topology only once.
    State sTempClean = m_mbs.realizeTopology();
    sTempClean.setQ(s.getQ());

    // Successively prune constraints until the projection is successful.
    while(true) {

        // Ensure at least one constraint will be enabled.
        SimTK_ASSERT_ALWAYS(!indexArray.empty(),
            "No valid position projection found by pruning search.");

        // Try this set of constraints.
        State sTemp(sTempClean);
        Real dist = evaluateProjection(sTemp, s, indexArray);
        if (PrintDebugInfoProjection)
            cout << "     " << indexArray << " d=" << dist << endl;

        // Exit if successful; otherwise, remove the constraint associated with
        // the most distant proximal point.
        if (dist < SimTK::Infinity) {
            s.setQ(sTemp.getQ());
            break;
        } else {
            Real maxDist = 0;
            int  maxIdx;

            for (int i=0; i<(int)indexArray.size(); ++i) {
                Real currDist = m_brick.get_pipConstraintHeightInGround(s,
                                    m_proximalPointIndices[indexArray[i]]);
                if (currDist > maxDist) {
                    maxDist = currDist;
                    maxIdx  = i;
                }
            }
            indexArray.eraseFast(&indexArray[maxIdx]);
        }
    }

    if (PrintDebugInfoProjection) {
        cout << "  -> Pruning search selected constraints " << indexArray
             << endl;
        displayNewProximalPoints(s);
        printHorizontalRule(0,1,'*');
    }
}

//------------------------------ Private methods -------------------------------
void PositionProjecter::createProximalIndexArrays(int sizeOfSet,
    Array_<Array_<ProximalPointIndex> >& arrayOfIndexArrays) const
{
    SimTK_ASSERT(arrayOfIndexArrays.empty(), "Input array must be empty.");
    const int numArrays = (1 << sizeOfSet) - 1;
    for (int ctr=1; ctr<=numArrays; ++ctr) {    //Exclude empty set.
        Array_<ProximalPointIndex> indexArray;
        for (ProximalPointIndex idx(0); idx<sizeOfSet; ++idx)
            if (ctr & (1 << idx))
                indexArray.push_back(idx);
        arrayOfIndexArrays.push_back(indexArray);
    }
}

Real PositionProjecter::evaluateProjection(State& sTemp, const State& sOrig,
                            const Array_<ProximalPointIndex>& indexArray) const
{
    Real dist = SimTK::Infinity;

    // Enable constraints.
    for (ProximalPointIndex i(0); i<(int)indexArray.size(); ++i)
        m_brick.enablePipConstraint(sTemp,
                                    m_proximalPointIndices[indexArray[i]]);

    // Try projecting. Simbody will throw an exception if the projection
    // accuracy cannot be achieved.
    try {
        m_mbs.projectQ(sTemp, TolProjectQ);
        if (!m_brick.isBrickInterpenetrating(sTemp))
            dist = (sOrig.getQ() - sTemp.getQ()).norm();
    } catch(...) {}

    return dist;
}

void PositionProjecter::displayNewProximalPoints(State& s) const
{
    // New positions.
    m_mbs.realize(s, Stage::Position);
    Array_<Vec3,BrickVertexIndex> postProjectionPos;
    m_brick.findAllLowestPointLocationsInGround(s, postProjectionPos);

    // New proximal points.
    Array_<BrickVertexIndex,ProximalPointIndex> proximalPointIndices;
    m_brick.findProximalPointIndices(postProjectionPos, proximalPointIndices);

    // Display.
    for (ProximalPointIndex i(0); i<(int)proximalPointIndices.size(); ++i)
        cout << "     [" << i << "] p="
             << postProjectionPos[proximalPointIndices[i]] << endl;
}


//==============================================================================
//                           Implementation: IMPACTER
//==============================================================================
Impacter::
Impacter(MultibodySystem& mbs, FreeUnilateralBrick& brick, const State& s0,
         const Array_<Vec3,BrickVertexIndex>& allPositionsInG,
         const Array_<BrickVertexIndex,ProximalPointIndex>&
               proximalPointIndices)
    : m_mbs(mbs), m_brick(brick), m_proximalPointIndices(proximalPointIndices)
{
    // Adjust the position of the Ball constraints corresponding to the proximal
    // points.
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        m_brick.set_ballConstraintLocation(s0, m_proximalPointIndices[i],
            allPositionsInG[m_proximalPointIndices[i]]);

    if (PrintDebugInfoImpact) {
        printHorizontalRule(1,0,'*',"starting impact");
        cout << "  -> " << m_proximalPointIndices.size() << " proximal point(s)"
             << endl;
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
            cout << "     [" << i << "] p="
                 << allPositionsInG[m_proximalPointIndices[i]] << endl;
    }
}

//------------------------ Perform one complete impact -------------------------
Real Impacter::
performImpactExhaustive(State& s,
                        Array_<Vec3,ProximalPointIndex>& proximalVelsInG,
                        Array_<bool,ProximalPointIndex>& hasRebounded) const
{
    // Calculate coefficients of restitution.
    Array_<Real> CORs(m_proximalPointIndices.size());
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i) {
        CORs [i] = hasRebounded[i] ? 0.0 : calcCOR(proximalVelsInG[i][ZAxis]);
        if (PrintDebugInfoImpact)
            cout << "  ** CORs[" << i << "] = " << CORs[i] << endl;
    }

    // Enumerate all active set candidates.
    Array_<ActiveSetCandidate> activeSetCandidates;
    initializeActiveSetCandidateArray(activeSetCandidates);
    if (PrintDebugInfoImpact)
        cout << "  -> " << activeSetCandidates.size()
             << " active set candidate(s)" << endl;

    // Interval-stepping loop.
    ImpactPhase impactPhase    = Compression;
    Vector restitutionImpulses = Vector(m_proximalPointIndices.size(), 0.0);
    Vector energyDissipated    = Vector(m_proximalPointIndices.size(), 0.0);
    int intervalCtr            = 0;
    while(true) {
        ++intervalCtr;
        if (PrintDebugInfoImpact) {
            std::stringstream sstream;
            sstream << (impactPhase == Compression ? "compression"
                                                   : "restitution")
                    << " interval " << intervalCtr;
            printHorizontalRule(1,0,'#',sstream.str());
        }

        // Clear data associated with each active set candidate.
        resetActiveSetCandidateArray(activeSetCandidates);

        // Generate and solve a linear system of equations for each active set
        // candidate. Categorize each solution and calculate fitness value.
        for (int i=0; i<(int)activeSetCandidates.size(); ++i) {
            if (PrintDebugInfoImpact) {
                cout << "  -> evaluating candidate " << i << " ";
                printFormattedActiveSet(cout,
                    activeSetCandidates[i].tangentialStates,
                    (impactPhase == Compression ? "c" : "r"));
                cout << endl;
            }

            generateAndSolveLinearSystem(s, impactPhase, restitutionImpulses,
                                         activeSetCandidates[i]);
            evaluateLinearSystemSolution(s, impactPhase, restitutionImpulses,
                                         activeSetCandidates[i]);

            if (PrintDebugInfoImpact)
                printActiveSetInfo(activeSetCandidates[i]);
            if (PauseAfterEachActiveSetCandidate)
                char trash = getchar();
        }

        if (PrintDebugInfoActiveSetSearch) {
            printHorizontalRule(2,0,'=',"exhaustive search summary");

            // Count instances of each category.
            Array_<int> solCatCtr(11,0);
            for (int i=0; i<(int)activeSetCandidates.size(); ++i)
                ++solCatCtr[activeSetCandidates[i].solutionCategory];

            for (int i=0; i<=SolCat_NotEvaluated; ++i)
                cout << "     " << std::setw(2) << solCatCtr[i] << "  "
                     << getSolutionCategoryDescription(SolutionCategory(i))
                     << endl;
            printHorizontalRule(0,1,'=');
        }

        // Select active set candidate from the best (lowest) category with the
        // best (lowest) fitness value.
        int  bestIdx     = -1;
        Real bestFitness = SimTK::Infinity;
        for (int solcat=0;
             solcat<=SolCat_WorstCatWithSomeUsableSolution; ++solcat)
        {
            for (int idx=0; idx<(int)activeSetCandidates.size(); ++idx) {
                if (activeSetCandidates[idx].solutionCategory == solcat
                    && activeSetCandidates[idx].fitness < bestFitness) {
                        bestIdx     = idx;
                        bestFitness = activeSetCandidates[idx].fitness;
                }
            }

            // Halt if a usable solution was found in this category.
            if (bestFitness < SimTK::Infinity)
                break;
        }
        SimTK_ASSERT_ALWAYS(bestFitness < SimTK::Infinity,
            "No suitable active set found by exhaustive search.");

        if (PrintBasicInfo)
            printFormattedActiveSet(cout,
                activeSetCandidates[bestIdx].tangentialStates,
                (impactPhase == Compression ? "c" : "r"));
        if (PrintDebugInfoImpact) {
            cout << "  ** selected active set candidate " << bestIdx << endl;
            printActiveSetInfo(activeSetCandidates[bestIdx]);
        }

        // Determine step length and apply impulse.
        const Real steplength = calculateIntervalStepLength(s, proximalVelsInG,
                                activeSetCandidates[bestIdx], intervalCtr);
        SimTK_ASSERT_ALWAYS(!isNaN(steplength),
            "No suitable interval step length found.");

        s.updU() += steplength*
                    activeSetCandidates[bestIdx].systemVelocityChange;
        m_mbs.realize(s, Stage::Velocity);

        if (PrintDebugInfoImpact)
            cout << "  ** steplength = " << steplength
                 << "\n     newU = " << s.getU() << endl;

        // Perform part of energy dissipation calculation now (i.e., before
        // updating velocities). Initialize powerTime0 array to zero since power
        // will not be computed for observing proximal points.
        Array_<Real,ProximalPointIndex>
            powerTime0((int)m_proximalPointIndices.size(), 0.0);
        int constraintIdx = -1;
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        {
            if (activeSetCandidates[bestIdx].tangentialStates[i] > Observing) {
                ++constraintIdx;
                const Vec2 vTang  = proximalVelsInG[i].getSubVec<2>(0);
                const Vec2 piTang = Vec2(activeSetCandidates[bestIdx].
                                         localImpulses[constraintIdx*3],
                                         activeSetCandidates[bestIdx].
                                         localImpulses[constraintIdx*3+1])
                                    * steplength;

                // Use dot product here in case tangential velocity and impulse
                // are not quite collinear.
                powerTime0[i] = dot(vTang, piTang)
                                + (proximalVelsInG[i][ZAxis]
                                   * activeSetCandidates[bestIdx].
                                     localImpulses[constraintIdx*3+2]
                                   * steplength);
            }
        }

        // Calculate the new velocity of each proximal point.
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
            proximalVelsInG[i] = m_brick.findLowestPointVelocityInGround(s,
                                     m_proximalPointIndices[i]);
        if (PrintDebugInfoImpact)
            for (ProximalPointIndex i(0);
                 i<(int)m_proximalPointIndices.size(); ++i)
                cout << "     [" << i << "] v=" << proximalVelsInG[i] << endl;

        // Calculate energy dissipated.
        constraintIdx = -1;
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        {
            if (activeSetCandidates[bestIdx].tangentialStates[i] > Observing) {
                ++constraintIdx;
                const Vec2 vTang  = proximalVelsInG[i].getSubVec<2>(0);
                const Vec2 piTang = Vec2(activeSetCandidates[bestIdx].
                                         localImpulses[constraintIdx*3],
                                         activeSetCandidates[bestIdx].
                                         localImpulses[constraintIdx*3+1])
                                    * steplength;
                const Real powerTime1 = dot(vTang, piTang)
                    + (proximalVelsInG[i][ZAxis]
                       * activeSetCandidates[bestIdx].
                         localImpulses[constraintIdx*3+2] * steplength);

                energyDissipated[i] += (powerTime0[i] + powerTime1) * 0.5;
            }
        }

        // Update the impact phase and required restitution impulses; exit if
        // complete.
        if (impactPhase == Compression) {

            int constraintIdx = -1;
            for (ProximalPointIndex i(0);
                 i<(int)m_proximalPointIndices.size(); ++i)
            {
                if (activeSetCandidates[bestIdx].tangentialStates[i]
                    > Observing)
                {
                    ++constraintIdx;
                    restitutionImpulses[i] += -activeSetCandidates[bestIdx].
                                              localImpulses[constraintIdx*3+2]
                                              * CORs[i] * steplength;
                }
            }

            if (!m_brick.isBrickImpacting(proximalVelsInG)) {
                if (PrintDebugInfoImpact)
                    cout << "  ** compression phase complete" << endl;

                // Proceed to restitution phase if any restitution impulses must
                // be applied; exit otherwise.
                Real maxRestImpulse = 0;
                for (int i=0; i<(int)restitutionImpulses.size(); ++i)
                    maxRestImpulse = std::max(maxRestImpulse,
                                              restitutionImpulses[i]);
                if (maxRestImpulse < MinMeaningfulImpulse) {
                    if (PrintDebugInfoImpact)
                        cout << "  ** no restitution impulses" << endl;
                    break;
                } else {
                    impactPhase = Restitution;
                    intervalCtr = 0;
                }
            }

        } else if (impactPhase == Restitution) {

            int constraintIdx = -1;
            for (ProximalPointIndex i(0);
                 i<(int)m_proximalPointIndices.size(); ++i)
            {
                if (activeSetCandidates[bestIdx].tangentialStates[i]
                    > Observing)
                {
                    ++constraintIdx;
                    restitutionImpulses[i] -= -activeSetCandidates[bestIdx].
                                              localImpulses[constraintIdx*3+2]
                                              * steplength;

                    if (std::abs(-activeSetCandidates[bestIdx].
                                 localImpulses[constraintIdx*3+2])
                        > MinMeaningfulImpulse)

                        hasRebounded[i] = true;
                }
            }

            Real maxRestImpulse = 0;
            for (int i=0; i<(int)restitutionImpulses.size(); ++i)
                maxRestImpulse = std::max(maxRestImpulse,
                                          restitutionImpulses[i]);
            if (maxRestImpulse < MinMeaningfulImpulse) {
                if (PrintDebugInfoImpact)
                    cout << "  ** restitution phase complete" << endl;
                break;
            }
        }

        if (PrintDebugInfoImpact) {
            cout << "     restitutionImpulses = " << restitutionImpulses
                 << "\n     hasRebounded = " << hasRebounded << endl;
        }
        if (PauseAfterEachImpactInterval)
            char trash = getchar();

    } //end interval-stepping loop
    if (PrintBasicInfo)
        cout << endl;

    // Ensure the total energy dissipated over the entire impact is positive for
    // each proximal point.
    for (int i=0; i<(int)energyDissipated.size(); ++i)
        SimTK_ASSERT(energyDissipated[i] > -SimTK::SignificantReal,
        "Negative energy dissipated.");

    return energyDissipated.sum();
}

Real Impacter::
performImpactPruning(State& s,
                     Array_<Vec3,ProximalPointIndex>& proximalVelsInG,
                     Array_<bool,ProximalPointIndex>& hasRebounded) const
{
    // Calculate coefficients of restitution.
    Array_<Real,ProximalPointIndex> CORs(m_proximalPointIndices.size());
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i) {
        CORs[i] = hasRebounded[i] ? 0.0 : calcCOR(proximalVelsInG[i][ZAxis]);
        if (PrintDebugInfoImpact)
            cout << "  ** CORs[" << i << "] = " << CORs[i] << endl;
    }

    // Interval-stepping loop.
    ImpactPhase impactPhase    = Compression;
    Vector restitutionImpulses = Vector(m_proximalPointIndices.size(), 0.0);
    Vector energyDissipated    = Vector(m_proximalPointIndices.size(), 0.0);
    int    intervalCtr         = 0;
    while(true) {
        ++intervalCtr;
        if (PrintDebugInfoImpact) {
            std::stringstream sstream;
            sstream << (impactPhase == Compression ? "compression"
                                                   : "restitution")
                    << " interval " << intervalCtr;
            printHorizontalRule(1,0,'#',sstream.str());
        }

        // Begin with maximal active set, but disallow sticking at points whose
        // tangential velocity magnitude is too large.
        ActiveSetCandidate asc;
        initializeActiveSetCandidate(s,asc);
        if (PrintDebugInfoImpact) {
            cout << "  -> Starting with active set candidate ";
            printFormattedActiveSet(cout, asc.tangentialStates, "");
            cout << endl;
        }

        // If the pruning search does not converge to an active set candidate
        // producing a solution with no violations, use the least-objectionable
        // active set candidate that was examined; throw an exception if all
        // examined active set candidates are intolerable.
        ActiveSetCandidate bestAsc;
        bestAsc.solutionCategory = SolCat_NotEvaluated;

        // Prune until active set candidate is empty.
        while (!isActiveSetCandidateEmpty(asc)) {
            resetActiveSetCandidate(asc);

            if (PrintDebugInfoImpact) {
                cout << "  -> evaluating candidate ";
                printFormattedActiveSet(cout, asc.tangentialStates,
                    (impactPhase == Compression ? "c" : "r"));
                cout << endl;
            }

            // Generate and solve a linear system of equations for this active
            // set candidate. Categorize the solution and calculate the fitness
            // value.
            generateAndSolveLinearSystem(s,impactPhase,restitutionImpulses,asc);
            evaluateLinearSystemSolution(s,impactPhase,restitutionImpulses,asc);

            if (PrintDebugInfoImpact)
                printActiveSetInfo(asc);
            if (PauseAfterEachActiveSetCandidate)
                char trash = getchar();

            // Keep track of least-objectionable active set candidate.
            if (asc.solutionCategory <= bestAsc.solutionCategory ||
                (asc.solutionCategory == bestAsc.solutionCategory
                 && asc.fitness <= bestAsc.fitness)) {
                    bestAsc = asc;
            }

            // Exit if acceptable active set candidate has been found.
            if (bestAsc.solutionCategory
                <= SolCat_WorstCatWithNoSeriousViolations)
            {
                if (PrintDebugInfoActiveSetSearch)
                    cout << "     Acceptable active set candidate found"
                         << endl;
                break;
            }

            // Exit if current active set candidate applies no impulses; cannot
            // make progress by pruning.
            if (asc.solutionCategory == SolCat_NoImpulsesApplied) {
                if (PrintDebugInfoActiveSetSearch)
                    cout << "     No impulses applied; cannot prune further"
                         << endl;
                break;
            }

            // Prune worst constraint.
            SimTK_ASSERT(asc.tangentialStates[asc.worstConstraint] > Observing,
                "Invalid worstConstraint index.");
            asc.tangentialStates[asc.worstConstraint] -= 1;

        } //end pruning loop

        // Throw an exception if no tolerable active set candidate was found.
        SimTK_ASSERT_ALWAYS(bestAsc.solutionCategory <=
            SolCat_WorstCatWithSomeUsableSolution,
            "Pruning search failed to find suitable active set candidate.");

        if (PrintBasicInfo)
            printFormattedActiveSet(cout, bestAsc.tangentialStates,
                                    (impactPhase == Compression ? "c" : "r"));
        if (PrintDebugInfoImpact) {
            cout << "  ** selected active set candidate ";
            printFormattedActiveSet(cout, bestAsc.tangentialStates);
            printActiveSetInfo(bestAsc);
        }

        // Determine step length and apply impulse.
        const Real steplength = calculateIntervalStepLength(s, proximalVelsInG,
                                                            bestAsc,
                                                            intervalCtr);
        SimTK_ASSERT_ALWAYS(!isNaN(steplength),
            "No suitable interval step length found.");

        s.updU() += steplength*bestAsc.systemVelocityChange;
        m_mbs.realize(s, Stage::Velocity);

        if (PrintDebugInfoImpact)
            cout << "  ** steplength = " << steplength
                 << "\n     newU = " << s.getU() << endl;

        // Perform part of energy dissipation calculation now (i.e., before
        // updating velocities). Initialize powerTime0 array to zero since power
        // will not be computed for observing proximal points.
        Array_<Real,ProximalPointIndex>
            powerTime0((int)m_proximalPointIndices.size(), 0.0);
        int constraintIdx = -1;
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        {
            if (bestAsc.tangentialStates[i] > Observing) {
                ++constraintIdx;
                const Vec2 vTang  = proximalVelsInG[i].getSubVec<2>(0);
                const Vec2 piTang = Vec2(bestAsc.localImpulses[constraintIdx*3],
                                    bestAsc.localImpulses[constraintIdx*3+1])
                                    * steplength;

                // Use dot product here in case tangential velocity and impulse
                // are not quite collinear.
                powerTime0[i] = dot(vTang, piTang)
                                + (proximalVelsInG[i][ZAxis]
                                   * bestAsc.localImpulses[constraintIdx*3+2]
                                   * steplength);
            }
        }

        // Calculate the new velocity of each proximal point.
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
            proximalVelsInG[i] = m_brick.findLowestPointVelocityInGround(s,
                                     m_proximalPointIndices[i]);
        if (PrintDebugInfoImpact)
            for (ProximalPointIndex i(0);
                 i<(int)m_proximalPointIndices.size(); ++i)
                cout << "     [" << i << "] v=" << proximalVelsInG[i] << endl;

        // Calculate energy dissipated.
        constraintIdx = -1;
        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        {
            if (bestAsc.tangentialStates[i] > Observing) {
                ++constraintIdx;
                const Vec2 vTang  = proximalVelsInG[i].getSubVec<2>(0);
                const Vec2 piTang = Vec2(bestAsc.localImpulses[constraintIdx*3],
                                    bestAsc.localImpulses[constraintIdx*3+1])
                                    * steplength;
                const Real powerTime1 = dot(vTang, piTang)
                    + (proximalVelsInG[i][ZAxis]
                       * bestAsc.localImpulses[constraintIdx*3+2] * steplength);

                energyDissipated[i] += (powerTime0[i] + powerTime1) * 0.5;
            }
        }

        // Update the impact phase and required restitution impulses; exit if
        // complete.
        if (impactPhase == Compression) {

            int constraintIdx = -1;
            for (ProximalPointIndex i(0);
                 i<(int)m_proximalPointIndices.size(); ++i)
            {
                if (bestAsc.tangentialStates[i] > Observing) {
                    ++constraintIdx;
                    restitutionImpulses[i] += -bestAsc.
                                              localImpulses[constraintIdx*3+2]
                                              * CORs[i] * steplength;
                }
            }

            if (!m_brick.isBrickImpacting(proximalVelsInG)) {
                if (PrintDebugInfoImpact)
                    cout << "  ** compression phase complete" << endl;

                // Proceed to restitution phase if any restitution impulses must
                // be applied; exit otherwise.
                Real maxRestImpulse = 0;
                for (int i=0; i<(int)restitutionImpulses.size(); ++i)
                    maxRestImpulse = std::max(maxRestImpulse,
                                              restitutionImpulses[i]);
                if (maxRestImpulse < MinMeaningfulImpulse) {
                    if (PrintDebugInfoImpact)
                        cout << "  ** no restitution impulses" << endl;
                    break;
                } else {
                    impactPhase = Restitution;
                    intervalCtr = 0;
                }
            }

        } else if (impactPhase == Restitution) {

            int constraintIdx = -1;
            for (ProximalPointIndex i(0);
                 i<(int)m_proximalPointIndices.size(); ++i)
            {
                if (bestAsc.tangentialStates[i] > Observing) {
                    ++constraintIdx;
                    restitutionImpulses[i] -= -bestAsc.
                                              localImpulses[constraintIdx*3+2]
                                              * steplength;

                    if (std::abs(-bestAsc.localImpulses[constraintIdx*3+2])
                        > MinMeaningfulImpulse)

                        hasRebounded[i] = true;
                }
            }

            Real maxRestImpulse = 0;
            for (int i=0; i<(int)restitutionImpulses.size(); ++i)
                maxRestImpulse = std::max(maxRestImpulse,
                                          restitutionImpulses[i]);
            if (maxRestImpulse < MinMeaningfulImpulse) {
                if (PrintDebugInfoImpact)
                    cout << "  ** restitution phase complete" << endl;
                break;
            }
        }

        if (PrintDebugInfoImpact) {
            cout << "     restitutionImpulses = " << restitutionImpulses
                 << "\n     hasRebounded = " << hasRebounded << endl;
        }
        if (PauseAfterEachImpactInterval)
            char trash = getchar();

    } //end interval-stepping loop
    if (PrintBasicInfo)
        cout << endl;

    // Ensure the total energy dissipated over the entire impact is positive for
    // each proximal point.
    for (int i=0; i<(int)energyDissipated.size(); ++i)
        SimTK_ASSERT(energyDissipated[i] > -SimTK::SignificantReal,
        "Negative energy dissipated.");

    return energyDissipated.sum();
}

//-------------------------- Private methods: display --------------------------
void Impacter::
printFormattedActiveSet(std::ostream& stream,
                        const Array_<int,ProximalPointIndex>& tangentialStates,
                        const std::string& prefix) const
{
    stream << "(" << prefix;
    for (ProximalPointIndex i(0); i<(int)tangentialStates.size(); ++i) {
        switch (tangentialStates[i]) {
          case Observing:   stream << "O";  break;
          case Rolling:     stream << "R";  break;
          case Sliding:     stream << "S";  break;
          default:          stream << "?";  break;
        }
    }
    stream << ")";
}

const char* Impacter::
getSolutionCategoryDescription(SolutionCategory solCat) const
{
    switch (solCat) {
      case SolCat_NoViolations:
          return "No violations";
      case SolCat_ActiveConstraintDoesNothing:
          return "Active constraint is doing nothing";
      case SolCat_RestitutionImpulsesIgnored:
          return "Restitution impulses were ignored";
      case SolCat_TangentialVelocityTooLargeToStick:
          return "Sticking not possible at this velocity";
      case SolCat_StickingImpulseExceedsStictionLimit:
          return "Sticking impulse exceeds stiction limit";
      case SolCat_GroundAppliesAttractiveImpulse:
          return "Ground applying attractive impulse";
      case SolCat_NegativePostCompressionNormalVelocity:
          return "Post-compression velocity is negative";
      case SolCat_NoImpulsesApplied:
          return "No impulses applied; no progress made";
      case SolCat_UnableToResolveUnknownSlipDirection:
          return "Unable to calculate unknown slip direction";
      case SolCat_MinStepCausesSlipDirectionReversal:
          return "Slip direction reverses with minimum step";
      case SolCat_NotEvaluated:
          return "Not yet evaluated";
      default:
          return "Unrecognized solution category";
    }
}

void Impacter::printActiveSetInfo(const ActiveSetCandidate& asc) const
{
    printf("\n");
    for (int i=0; i<40; ++i) printf("-");
    printf(" ");
    printFormattedActiveSet(cout, asc.tangentialStates);
    printf("\n");

    cout << "     deltaU   = " << asc.systemVelocityChange << "\n"
         << "     impulse  = " << asc.localImpulses << "\n"
         << "     category = "
         << getSolutionCategoryDescription(asc.solutionCategory) << "\n"
         << "     fitness  = " << asc.fitness << "\n"
         << "     worstIdx = " << asc.worstConstraint << endl;
}

//--------------------- Private methods: helper functions ----------------------
bool Impacter::
addInBaseN(int N, Array_<int,ProximalPointIndex>& vec, int addThis) const
{
    vec[ProximalPointIndex(0)] += addThis;
    for (ProximalPointIndex i(0); i<(int)vec.size(); ++i) {

        // Detect overflow.
        if (vec[i] >= N && i == vec.size()-1)
            return false;

        // Restore vec to base N.
        while (vec[i] >= N) {
            vec[ProximalPointIndex(i+1)] += 1;
            vec[i] -= N;
        }
    }
    return true;
}

void Impacter::initializeActiveSetCandidateArray(
    Array_<ActiveSetCandidate>& activeSetCandidates) const
{
    SimTK_ASSERT(activeSetCandidates.empty(), "Input array must be empty.");
    const int n = (int)m_proximalPointIndices.size();
    Array_<int,ProximalPointIndex> tangentialStateArray(n,0);

    // Increment until overflow in base 3 to enumerate all possibilities.
    while (addInBaseN(3, tangentialStateArray, 1)) {
        ActiveSetCandidate candidate;
        candidate.tangentialStates = tangentialStateArray;
        activeSetCandidates.push_back(candidate);
    }
}

void Impacter::initializeActiveSetCandidate(const State& s,
                                            ActiveSetCandidate& asc) const
{
    const int n = (int)m_proximalPointIndices.size();
    Array_<int,ProximalPointIndex> tangentialStateArray(n,Rolling);
    for (ProximalPointIndex i(0); i<n; ++i) {
        const Vec3 vel = m_brick.findLowestPointVelocityInGround(s,
                             m_proximalPointIndices[i]);
        const Real tangVelMag = vel.getSubVec<2>(0).norm();
        if (tangVelMag > MaxStickingTangVel)
            tangentialStateArray[i] = Sliding;
    }
    asc.tangentialStates = tangentialStateArray;
}

void Impacter::resetActiveSetCandidate(ActiveSetCandidate& asc) const
{
    asc.systemVelocityChange = Vector(1,0.0);
    asc.localImpulses        = Vector(1,0.0);
    asc.solutionCategory     = SolCat_NotEvaluated;
    asc.fitness              = SimTK::Infinity;
    asc.worstConstraint      = ProximalPointIndex(
                                   std::numeric_limits<int>::max());
}

void Impacter::resetActiveSetCandidateArray(
    Array_<ActiveSetCandidate>& activeSetCandidates) const
{
    for (int i=0; i<(int)activeSetCandidates.size(); ++i)
        resetActiveSetCandidate(activeSetCandidates[i]);
}

int Impacter::getIndexOfFirstMultiplier(const State& s,
                                        const ProximalPointIndex i) const
{   return m_brick.get_ballConstraintFirstIndex(s, m_proximalPointIndices[i]); }

Real Impacter::calcAbsDiffBetweenAngles(Real a, Real b) const
{
    // Catch unknown angles.
    if (isNaN(a) || isNaN(b))
        return SimTK::NaN;

    // Move angles from [-pi, pi] to [0, 2*pi].
    const Real twopi = 2.0*SimTK::Pi;
    if (a<0) a += twopi;
    if (b<0) b += twopi;

    // Subtract, avoiding negative answers; difference is in [0, 2*pi].
    Real absdif = std::abs(a-b);

    // Difference can be no greater than pi due to periodicity.
    return (absdif < SimTK::Pi ? absdif : twopi-absdif);
}

Real Impacter::calcSlidingStepLengthToOrigin(const Vec2& A, const Vec2& B) const
{
    // Take full step if initial tangential velocity was small (impending slip).
    if (A.norm() < MaxStickingTangVel) {
        if (PrintDebugInfoStepLength)
            cout << "     --> A.norm() < MaxStickingTangVel; returning 1.0"
                 << endl;
        return 1.0; //Q==B
    }

    const Vec2 P     = Vec2(0);
    const Vec2 AtoP  = P-A;
    const Vec2 AtoB  = B-A;
    const Real ABsqr = AtoB.normSqr();

    // Ensure line segment is of meaningful length.
    if (ABsqr < SimTK::SignificantReal) {
        if (PrintDebugInfoStepLength)
            cout << "     --> ABsqr < SimTK::SignificantReal; returning 1.0"
                 << endl;
        return 1.0; //Q==B
    }

    // Normalized distance from A to Q (Q = A + stepLength*AtoB).
    const Real stepLength = clamp(0.0, dot(AtoP,AtoB)/ABsqr, 1.0);

    if (PrintDebugInfoStepLength)
        cout << "     --> returning stepLength = " << stepLength << endl;
    return stepLength;
}

bool Impacter::isActiveSetCandidateEmpty(const ActiveSetCandidate& asc) const
{
    for (ProximalPointIndex i(0); i<(int)asc.tangentialStates.size(); ++i)
        if (asc.tangentialStates[i] > Observing)
            return false;
    return true;
}

//------------------------ Private methods: calculators ------------------------
Real Impacter::calcCOR(Real vNormal) const
{
    if (-vNormal < m_brick.get_vMinRebound())
        return 0.0;
    const Real CORline = ((m_brick.get_minCOR() - 1) /
                          m_brick.get_vPlasticDeform()) * (-vNormal) + 1;
    return std::max(CORline, m_brick.get_minCOR());
}

void Impacter::generateAndSolveLinearSystem(const State& s0,
                                            const ImpactPhase impactPhase,
                                            const Vector& restitutionImpulses,
                                            ActiveSetCandidate& asc) const
{
    // Enable constraints to initialize the Jacobian.
    State s = m_mbs.realizeTopology();
    s.setQ(s0.getQ());
    s.setU(s0.getU());
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        if (asc.tangentialStates[i] > Observing)
            m_brick.enableBallConstraint(s, m_proximalPointIndices[i]);
    m_mbs.realize(s, Stage::Velocity);

    // Begin generating linear system to solve.
    Matrix MassMatrix, GMatrix;
    m_mbs.getMatterSubsystem().calcM(s, MassMatrix);
    m_mbs.getMatterSubsystem().calcG(s, GMatrix);
    const int N = MassMatrix.nrow();
    const int M = GMatrix.nrow();

    // TODO: Should not need to form a matrix larger than MxM here.
    //       1. Solve only for impulses.
    //       2. During restitution, remove equations with no degrees of freedom
    //          (i.e., all normal impulses, and tangential impulses associated
    //          with sliding points).
    //       3. Figure out how to modify Simbody's GM\~G computation to create
    //          G'M\~G, where G' means modifying the equations associated with
    //          sliding (and normal impulses in the restitution phase).
    Matrix A = Matrix(N+M, N+M);
    A(0,0,N,N) = MassMatrix;
    A(0,N,N,M) = GMatrix.transpose();
    A(N,0,M,N) = GMatrix;
    A(N,N,M,M) = 0.0;
    Vector b = Vector(N+M, 0.0);

    // Define equations.
    Array_<Real> slidingDirections;
    for (ProximalPointIndex idx(0);
         idx<(int)m_proximalPointIndices.size(); ++idx)
    {
        if (asc.tangentialStates[idx] > Observing) {

            // Current velocity at this proximal point.
            const Vec3 currVelAtPoint = m_brick.findLowestPointVelocityInGround(
                                                s, m_proximalPointIndices[idx]);

            // Row indices into matrix A corresponding to the constraints for
            // this proximal point.
            const int row_x = N + getIndexOfFirstMultiplier(s,idx);
            const int row_y = row_x + 1;
            const int row_z = row_y + 1;

            // Tangential directions.
            if (asc.tangentialStates[idx] == Rolling) {

                // Drive both components of tangential velocity to zero.
                b[row_x] = -currVelAtPoint[XAxis];
                b[row_y] = -currVelAtPoint[YAxis];

            } else if (asc.tangentialStates[idx] == Sliding) {

                // Apply friction impulse in the direction opposing the sliding
                // direction. At this point, set muDyn=0 for all points with
                // unknown sliding directions.
                A(row_x,0,2,N)  = 0.0;
                A[row_x][row_x] = 1;    b[row_x] = 0;
                A[row_y][row_y] = 1;    b[row_y] = 0;

                // Calculate theta, the angle between the global X-axis and the
                // tangential velocity vector.
                const Real theta = m_brick.
                                   findTangentialVelocityAngle(currVelAtPoint);
                slidingDirections.push_back(theta);

                if (PrintDebugInfoImpact)
                    cout << "  ** angle of tangential velocity vector for "
                         << "proximal point " << idx << " is " << theta << endl;

                if (!isNaN(theta)) {
                    const Real impulseDir = theta + SimTK::Pi;
                    A[row_x][row_z] = -m_brick.get_muDyn() * cos(impulseDir);
                    A[row_y][row_z] = -m_brick.get_muDyn() * sin(impulseDir);
                }
            }

            // Normal direction.
            if (impactPhase == Compression) {

                // Populate with the compression equation, which drives the
                // normal velocity of the impacting point to zero.
                b[row_z] = -currVelAtPoint[ZAxis];

            } else if (impactPhase == Restitution) {

                // Populate with the restitution equation, which sets the normal
                // impulse to the impulse required in the restitution phase.
                A(row_z,0,1,N)  = 0.0;
                A[row_z][row_z] = 1;
                b[row_z]        = -restitutionImpulses[idx];
            }

        } //end if not observing
    } //end for each proximal point

    // Iterate to find sliding directions, if necessary.
    if (!slidingDirections.empty()) {
        if (PrintDebugInfoImpact) {
            cout << "  ** finding sliding directions" << endl;

            int slidingPointNum = -1;
            for (ProximalPointIndex idx(0);
                 idx<(int)m_proximalPointIndices.size(); ++idx)
                if (asc.tangentialStates[idx] == Sliding) {
                    ++slidingPointNum;

                    // Current velocity at this proximal point.
                    const Vec3 vT = m_brick.findLowestPointVelocityInGround(s,
                                        m_proximalPointIndices[idx]);
                    cout << "     v[" << idx << "] = " << vT << " (angle = "
                         << slidingDirections[slidingPointNum] << " rad)"
                         << endl;
                }
        }

        int numIter = 0;
        ProximalPointIndex worstConstraintIdx;
        while(true) {

            // Halt if maximum number of iterations is reached.
            ++numIter;
            if (numIter > MaxIterSlipDirection) {
                if (PrintDebugInfoImpact)
                    cout << "  ** maximum number of iterations reached" << endl;

                asc.solutionCategory =
                    SolCat_UnableToResolveUnknownSlipDirection;
                asc.fitness          = SimTK::Infinity;
                asc.worstConstraint  = worstConstraintIdx;
                break;
            }
            if (PrintDebugInfoImpact)
                cout << "  ** beginning iteration " << numIter << endl;

            // Solve using current directions. We have set muDyn=0 for all
            // points with unknown sliding directions.
            FactorQTZ qtzA(A);
            Vector    sol;
            qtzA.solve(b, sol);

            // Calculate new system velocities using SlidingDirStepLength, which
            // we presume is sufficiently small to avoid direction reversals
            // (except in situations where |v_t| < MaxStickingTangVel).
            Vector deltaU      = sol(0,N);
            Vector calcImpulse = sol(N,M);

            if (PrintDebugInfoImpact)
                cout << "     calculated deltaU = " << deltaU << endl
                     << "     calculated impulse = " << calcImpulse << endl;

            State sTemp(s);
            sTemp.setU(s.getU() + SlidingDirStepLength*deltaU);
            m_mbs.realize(sTemp, Stage::Velocity);

            // Update directions of all sliding points (not just those with
            // unknown sliding directions).
            Real maxAngleDif     = 0;
            int  slidingPointNum = -1;
            for (ProximalPointIndex idx(0);
                 idx<(int)m_proximalPointIndices.size(); ++idx)
            {
                if (asc.tangentialStates[idx] == Sliding) {
                    ++slidingPointNum;

                    // Determine new angle from proposed velocity at this point.
                    const Vec3 vTemp = m_brick.findLowestPointVelocityInGround(
                                           sTemp, m_proximalPointIndices[idx]);
                    Real newAngle = m_brick.findTangentialVelocityAngle(vTemp);

                    if (PrintDebugInfoImpact)
                        cout << "     v[" << idx << "] = " << vTemp
                             << " (angle = " << newAngle << " rad)" << endl;

                    // Keep track of maximum absolute difference in angle (note
                    // that oldAngle and/or newAngle might be NaN).
                    Real oldAngle = slidingDirections[slidingPointNum];
                    slidingDirections[slidingPointNum] = newAngle;

                    if (isNaN(oldAngle) || isNaN(newAngle)) {
                        maxAngleDif        = SimTK::Infinity;
                        worstConstraintIdx = idx;
                    } else {
                        Real angleDif = calcAbsDiffBetweenAngles(oldAngle,
                                                                 newAngle);
                        if (angleDif > maxAngleDif) {
                            maxAngleDif        = angleDif;
                            worstConstraintIdx = idx;
                        }
                    }

                    if (PrintDebugInfoImpact)
                        cout << "     old angle = " << oldAngle
                             << ", new angle = " << newAngle << endl;

                    // Update linear system.
                    const int row_x = N + getIndexOfFirstMultiplier(s,idx);
                    const int row_y = row_x + 1;
                    const int row_z = row_y + 1;

                    if (!isNaN(newAngle)) {
                        const Real impulseDir = newAngle + SimTK::Pi;
                        A[row_x][row_z] = -m_brick.get_muDyn()*cos(impulseDir);
                        A[row_y][row_z] = -m_brick.get_muDyn()*sin(impulseDir);
                    } else {
                        A[row_x][row_z] = 0;
                        A[row_y][row_z] = 0;
                    }

                } //end if this point is sliding
            } //end for each proximal point

            if (PrintDebugInfoImpact)
                cout << "     maximum angle change of " << maxAngleDif << endl;

            // Exit if converged.
            if (maxAngleDif < TolMaxDifDirIteration) {
                if (PrintDebugInfoImpact)
                    cout << "  ** sliding directions converged" << endl;
                break;
            }

            // Exit if a sliding direction flips. By our assumption about the
            // smallness of SlidingDirStepLength, obtaining a flipping direction
            // indicates that this point should actually be sticking.
            if (std::abs(maxAngleDif-SimTK::Pi) < TolMaxDifDirIteration) {
                if (PrintDebugInfoImpact)
                    cout << "  ** point will stick, not slide" << endl;

                asc.solutionCategory =
                    SolCat_MinStepCausesSlipDirectionReversal;
                asc.fitness          = SimTK::Infinity;
                asc.worstConstraint  = worstConstraintIdx;
                return; //Giving up.
            }

        } //end while directions are unknown
    } //end if points are sliding

    // Either no points are sliding, or have finished iterating and need one
    // more solve to reconcile friction impulses with newest directions.
    // TODO: Do not need to solve a second time if initial directions were
    //       acceptable.
    // TODO: Replace the above fixed-point iteration with a Newton iteration
    //       (provided the NaN directions are resolved first).
    FactorQTZ qtzA(A);
    Vector    sol;
    qtzA.solve(b, sol);

    // Extract system velocity changes and local impulses from sol vector.
    asc.systemVelocityChange = sol(0,N);
    asc.localImpulses        = sol(N,M);

    if (PrintDebugInfoImpact) {
        cout << "     proximal point velocities after full step:" << endl;
        State sTemp(s);
        sTemp.updU() += 1.0*asc.systemVelocityChange;
        m_mbs.realize(sTemp, Stage::Velocity);

        for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
            cout << "     [" << i << "] v="
                 << m_brick.findLowestPointVelocityInGround(sTemp,
                            m_proximalPointIndices[i]) << endl;
    }
}

void Impacter::evaluateLinearSystemSolution(const State& s,
                                            const ImpactPhase impactPhase,
                                            const Vector& restitutionImpulses,
                                            ActiveSetCandidate& asc) const
{
    // Return if already evaluated; SolCat_UnableToResolveUnknownSlipDirection
    // and SolCat_MinStepCausesSlipDirectionReversal will have been caught in
    // generateAndSolveLinearSystem.
    if (asc.solutionCategory < SolCat_NotEvaluated)
        return;

    // Gather information about active set candidate.
    const int numImpulses = (int)asc.localImpulses.size();
    SimTK_ASSERT(numImpulses%3 == 0, "Invalid number of impulses.");
    const int numConstraints = numImpulses/3;

    // Calculate proximal point velocities after taking a full step.
    State sFullStep(s);
    sFullStep.updU() += 1.0*asc.systemVelocityChange;
    m_mbs.realize(sFullStep, Stage::Velocity);

    Array_<Vec3,ProximalPointIndex> fullStepVel(m_proximalPointIndices.size());
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i)
        fullStepVel[i] = m_brick.findLowestPointVelocityInGround(sFullStep,
                                 m_proximalPointIndices[i]);

    // No impulses applied; no progress made -- avoid infinite looping.
    if (asc.localImpulses.norm() < MinMeaningfulImpulse) {
        asc.solutionCategory = SolCat_NoImpulsesApplied;
        asc.fitness          = SimTK::Infinity;
        return;
    }

    // Post-compression velocity is negative -- the linear system was generated
    // incorrectly; the solution is nonsense. Not relevant for proximal points
    // that are just observing; these would be addressed in a follow-up impact.
    if (impactPhase == Compression) {
        Real minNormVel = SimTK::Infinity;
        ProximalPointIndex minNormVelIdx;
        for (ProximalPointIndex i(0); i<(int)fullStepVel.size(); ++i) {
            if (asc.tangentialStates[i] > Observing &&
                fullStepVel[i][ZAxis] < minNormVel) {

                minNormVel    = fullStepVel[i][ZAxis];
                minNormVelIdx = i;
            }
        }
        if (minNormVel < -TolVelocityFuzziness) {
            asc.solutionCategory = SolCat_NegativePostCompressionNormalVelocity;
            asc.fitness          = -minNormVel;
            asc.worstConstraint  = minNormVelIdx;
            return;
        }
    }

    // Ground applying attractive impulse -- the normal impulse must always be
    // negative (note the sign convention).
    Real maxNormImpulse = 0;
    ProximalPointIndex maxNormImpulseIdx;
    int constraintIdx = -1;
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i) {
        if (asc.tangentialStates[i] > Observing) {
            ++constraintIdx;

            if (asc.localImpulses[constraintIdx*3+2] > maxNormImpulse) {
                maxNormImpulse    = asc.localImpulses[constraintIdx*3+2];
                maxNormImpulseIdx = i;
            }
        }
    }
    if (maxNormImpulse > MinMeaningfulImpulse) {
        asc.solutionCategory = SolCat_GroundAppliesAttractiveImpulse;
        asc.fitness          = maxNormImpulse;
        asc.worstConstraint  = maxNormImpulseIdx;
        return;
    }

    // Sticking impulse exceeds stiction limit -- should be sliding instead.
    Real maxExcessiveImpulse = 0;
    ProximalPointIndex maxExcessiveImpulseIdx;
    constraintIdx = -1;
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i) {
        if (asc.tangentialStates[i] > Observing) {
            ++constraintIdx;

            if (asc.tangentialStates[i] == Rolling) {
                const int Xidx = constraintIdx*3;
                const Real impTangMag = Vec2(asc.localImpulses[Xidx],
                                             asc.localImpulses[Xidx+1]).norm();
                const Real impNormMag = -asc.localImpulses[Xidx+2];
                const Real excessiveImpulse = impTangMag
                                              - m_brick.get_muDyn()*impNormMag;

                if (excessiveImpulse > MinMeaningfulImpulse &&
                    excessiveImpulse > maxExcessiveImpulse) {
                    maxExcessiveImpulse    = excessiveImpulse;
                    maxExcessiveImpulseIdx = i;
                }
            }
        }
    }
    if (maxExcessiveImpulse > MinMeaningfulImpulse) {
        asc.solutionCategory = SolCat_StickingImpulseExceedsStictionLimit;
        asc.fitness          = maxExcessiveImpulse;
        asc.worstConstraint  = maxExcessiveImpulseIdx;
        return;
    }

    // Sticking not possible at this velocity -- the magnitude of the initial
    // tangential velocity (i.e., at the beginning of the interval) must be
    // sufficiently small to allow sticking.
    State sCurr(s);
    m_mbs.realize(sCurr, Stage::Velocity);
    Real maxTangVelMag = 0;
    ProximalPointIndex maxTangVelMagIdx;
    for (ProximalPointIndex i(0); i<(int)m_proximalPointIndices.size(); ++i) {
        if (asc.tangentialStates[i] == Rolling) {
            const Vec3 vel = m_brick.findLowestPointVelocityInGround(sCurr,
                                 m_proximalPointIndices[i]);
            const Real tangVelMag = vel.getSubVec<2>(0).norm();

            if (tangVelMag > maxTangVelMag) {
                maxTangVelMag    = tangVelMag;
                maxTangVelMagIdx = i;
            }
        }
    }
    if (maxTangVelMag > MaxStickingTangVel) {
        asc.solutionCategory = SolCat_TangentialVelocityTooLargeToStick;
        asc.fitness          = maxTangVelMag;
        asc.worstConstraint  = maxTangVelMagIdx;
        return;
    }

    // Restitution impulses were ignored -- avoid applying restitution impulses
    // sequentially (should be applied simultaneously).
    if (impactPhase == Restitution) {
        Real ignoredImpulse = 0;
        for (int i=0; i<(int)restitutionImpulses.size(); ++i)
            ignoredImpulse += restitutionImpulses[i];
        for (int i=0; i<numConstraints; ++i)
            ignoredImpulse -= -asc.localImpulses[i*3+2];

        if (ignoredImpulse > MinMeaningfulImpulse*restitutionImpulses.size()) {
            asc.solutionCategory = SolCat_RestitutionImpulsesIgnored;
            asc.fitness          = ignoredImpulse;
            return;
        }
    }

    // Active constraint is doing nothing -- prefer to avoid active sets with
    // constraints that apply no impulses.
    for (int i=0; i<numConstraints; ++i) {
        if (Vec3(asc.localImpulses[i*3], asc.localImpulses[i*3+1],
                 asc.localImpulses[i*3+2]).norm() < MinMeaningfulImpulse) {
            asc.solutionCategory = SolCat_ActiveConstraintDoesNothing;
            asc.fitness          = asc.localImpulses.norm(); //As usual.
            return;
        }
    }

    // No violations -- ideal case.
    asc.solutionCategory = SolCat_NoViolations;
    asc.fitness          = asc.localImpulses.norm();
}

Real Impacter::calculateIntervalStepLength(const State& s0,
                                           const Array_<Vec3>& currVels,
                                           const ActiveSetCandidate& asc,
                                           int intervalCtr) const
{
    // This is the return value, which will strictly decrease as we iterate
    // through the sliding proximal points.
    Real steplength = 1.0;

    // State at the beginning of the current interval.
    State s(s0);
    m_mbs.realize(s, Stage::Velocity);

    // Loop through each sliding point and reduce the step length, if necessary.
    for (ProximalPointIndex i(0); i<(int)asc.tangentialStates.size(); ++i) {
        if (asc.tangentialStates[i] == Sliding) {

            if (PrintDebugInfoStepLength)
                cout << "  ** analyzing proximal point " << i << "..." << endl;

            // Calculate the sliding direction angle for the ith proximal point
            // (which we know is sliding) at the beginning of the interval.
            const Real ang0 = m_brick.findTangentialVelocityAngle(currVels[i]);

            // Any step length is acceptable for this point if its initial
            // tangential velocity vector lies within the "capture velocity
            // circle" defined by MaxStickingTangVel (i.e., the point is in
            // impending slip).
            if (isNaN(ang0)) {
                if (PrintDebugInfoStepLength)
                    cout << "  -- finished with proximal point " << i
                         << " (ang0 is NaN); current steplength is "
                         << steplength << endl;
                break; //Proceed to next point.
            }

            // Proposed system velocities after taking a step of size steplength
            // (which may have already been reduced from 1.0 by another point).
            State sProp(s);
            sProp.setU(s.getU() + steplength*asc.systemVelocityChange);
            m_mbs.realize(sProp, Stage::Velocity);

            // Calculate the proposed velocity and sliding direction angle of
            // this proximal point, given the current proposed step length.
            Vec3 propVel = m_brick.findLowestPointVelocityInGround(sProp,
                                   m_proximalPointIndices[i]);
            Real ang1 = m_brick.findTangentialVelocityAngle(propVel);

            // The current proposed step length is acceptable for this point if
            // its proposed tangential velocity vector lies within the "capture
            // velocity circle" (i.e., the point will transition from sliding to
            // rolling).
            if (isNaN(ang1)) {
                if (PrintDebugInfoStepLength)
                    cout << "  -- finished with proximal point " << i
                         << " (ang1 is NaN); current steplength is "
                         << steplength << endl;
                break; //Proceed to next point.
            }

            // Calculate the absolute difference between the angles of the
            // current and proposed tangential velocities.
            Real absAngDif = calcAbsDiffBetweenAngles(ang0, ang1);

            if (PrintDebugInfoStepLength)
                cout << "     currVels[" << i << "] = " << currVels[i] << " ("
                     << ang0 << " rad = " << ang0*180.0/SimTK::Pi << " deg)\n"
                     << "     propVel = " << propVel << " (" << ang1
                     << " rad = " << ang1*180.0/SimTK::Pi << " deg)\n"
                     << "     absAngDif = " << absAngDif << " rad = "
                     << absAngDif*180.0/SimTK::Pi << " deg\n" << endl;

            // Check whether the direction change is sufficiently small for this
            // proximal point.
            if (absAngDif < MaxSlidingDirChange + TolSlidingDirChange) {
                if (PrintDebugInfoStepLength)
                    cout << "  -- finished with proximal point " << i
                         << " (no iteration required); current steplength is "
                         << steplength << endl;
                break; //Proceed to next point.
            }

            // Search for a step length that results in a sliding direction
            // angle change within TolSlidingDirChange of MaxSlidingDirChange;
            // break if a transition from sliding to rolling is detected.
            // Revert to bisection method if secant fails.
            int numIter = 0;
            bool usingBisection = false;

            // Initialize variables for secant method.
            //   x-axis: multiplicative factor for steplength, in [0,1].
            //   y-axis: {absolute angle change} - MaxSlidingDirChange.
            Real x0 = 0;
            Real x1 = 1.0;
            Real y0 = -MaxSlidingDirChange;
            Real y1 = absAngDif - MaxSlidingDirChange;
            const Real y1initial = y1;

            while(true) {
                ++numIter;
                if (numIter > MaxIterStepLength) {
                    cout << "  ** maximum number of iterations reached" << endl;
                    return SimTK::NaN;
                }

                if (PrintDebugInfoStepLength)
                    cout << "     x0=" << x0 << ", y0=" << y0 << ", x1=" << x1
                         << ", y1=" << y1 << endl;

                // Use secant method to find a zero-crossing on the (steplength
                // factor) vs. (absolute angle change - MaxSlidingDirChange)
                // plot. Switch to bisection method to avoid divide-by-zero.
                Real x2 = SimTK::NaN;
                if (usingBisection) {
                    x2 = (x0 + x1) * 0.5;
                    if (PrintDebugInfoStepLength)
                        cout << "     x2=" << x2 << " (bisection)" << endl;
                } else if (std::abs(y1-y0) < SimTK::SignificantReal) {
                    numIter = 1;
                    usingBisection = true;

                    if (y0*y1 >= 0) {
                        // Reinitialize.
                        x0 = 0;
                        x1 = 1.0;
                        y0 = -MaxSlidingDirChange;
                        y1 = y1initial;
                    }
                    x2 = 0.5;

                    if (PrintDebugInfoStepLength)
                        cout << "     x2=" << x2 << " (starting bisection)"
                             << endl;
                } else {
                    x2 = x1 - y1*(x1-x0)/(y1-y0);
                    if (PrintDebugInfoStepLength)
                        cout << "     x2=" << x2 << " (secant)" << endl;
                }

                // Ensure factor remains in [0,1]; otherwise, use bisection step
                // instead.
                if (x2 <= 0 || x2 >= 1) {
                    numIter = 1;
                    usingBisection = true;

                    if (y0*y1 >= 0) {
                        // Reinitialize.
                        x0 = 0;
                        x1 = 1.0;
                        y0 = -MaxSlidingDirChange;
                        y1 = y1initial;
                    }
                    x2 = 0.5;

                    if (PrintDebugInfoStepLength)
                        cout << "     x2=" << x2 << " (starting bisection)"
                             << endl;
                }

                // Calculate the new velocity and sliding direction angle, given
                // the new proposed step length.
                if (PrintDebugInfoStepLength)
                    cout << "  -- calculating velocity and sliding direction "
                         << "angle using steplength " << steplength*x2 << endl;
                sProp.setU(s.getU() + steplength*x2*asc.systemVelocityChange);
                m_mbs.realize(sProp, Stage::Velocity);
                propVel = m_brick.findLowestPointVelocityInGround(sProp,
                                  m_proximalPointIndices[i]);
                ang1 = m_brick.findTangentialVelocityAngle(propVel);

                // Break if a transition from sliding to rolling is detected.
                if (isNaN(ang1)) {
                    steplength *= x2;
                    if (PrintDebugInfoStepLength)
                        cout << "  -- finished with proximal point " << i
                             << " (ang1 is NaN); current steplength is "
                             << steplength << endl;
                    break; //Proceed to next point.
                }

                // Calculate the absolute difference between the angles of the
                // current and new proposed tangential velocities; break if the
                // direction change is within the target range.
                absAngDif = calcAbsDiffBetweenAngles(ang0, ang1);
                if (std::abs(absAngDif - MaxSlidingDirChange)
                    < TolSlidingDirChange)
                {
                    steplength *= x2;
                    if (PrintDebugInfoStepLength)
                        cout << "  -- finished with proximal point " << i
                             << " (converged); current steplength is "
                             << steplength << endl;
                    break; //Proceed to next point.
                }

                // Update variables for next iteration of secant method.
                const Real y2 = absAngDif - MaxSlidingDirChange;
                if (usingBisection) {
                    if (y2*y0 >= 0) {
                        // Replace (x0,y0) with (x2,y2).
                        x0 = x2;
                        y0 = y2;
                    } else {
                        // Replace (x1,y1) with (x2,y2).
                        x1 = x2;
                        y1 = y2;
                    }
                } else {
                    x0 = x1;
                    x1 = x2;
                    y0 = y1;
                    y1 = y2;
                }

            } //end while for this proximal point
        } //end if sliding
    } //end for each proximal point

    if (PrintDebugInfoStepLength)
        cout << "  ** returning steplength = " << steplength << endl;

    return steplength;
}
