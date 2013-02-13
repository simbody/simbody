/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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

#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"

#include <vector>
#include <map>

using namespace SimTK;
using namespace std;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

static const int NUM_BODIES = 10;
static const Real BOND_LENGTH = 0.5;
static const int ITERATIONS = 4;
static const Real TOL = 1e-4;

bool testFitting
   (const MultibodySystem& mbs, State& state, 
    const vector<MobilizedBodyIndex>& bodyIxs, 
    const vector<vector<Vec3> >& stations, 
    const vector<vector<Vec3> >& targetLocations, 
    Real minError, Real maxError, Real endDistance) 
{    
    // Find the best fit.
    
    Real reportedError = ObservedPointFitter::findBestFit(mbs, state, bodyIxs, stations, targetLocations, TOL);
    cout << "[min,max]=" << minError << "," << maxError << " actual=" << reportedError << endl;
    bool result = (reportedError <= maxError && reportedError >= minError);
    
    // Verify that the error was calculated correctly.
    
    Real error = 0.0;
    int numStations = 0;
    mbs.realize(state, Stage::Position);
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    for (int i = 0; i < (int) bodyIxs.size(); ++i) {
        MobilizedBodyIndex id = bodyIxs[i];
        numStations += (int)stations[i].size();
        for (int j = 0; j < (int) stations[i].size(); ++j)
            error += (targetLocations[i][j]-matter.getMobilizedBody(id).getBodyTransform(state)*stations[i][j]).normSqr();
    }
    error = std::sqrt(error/numStations);
    cout << "calc wrms=" << error << endl;
    ASSERT(std::abs(1.0-error/reportedError) < 0.0001); // should match to machine precision
    
    if (endDistance >= 0) {
        // Verify that the ends are the correct distance apart.
        Real distance = (matter.getMobilizedBody(bodyIxs[0]).getBodyOriginLocation(state)-matter.getMobilizedBody(bodyIxs[bodyIxs.size()-1]).getBodyOriginLocation(state)).norm();
        cout << "required dist=" << endDistance << ", actual=" << distance << endl;
        ASSERT(std::abs(1.0-endDistance/distance) < TOL);
    }

    return result;
}


static void testObservedPointFitter(bool useConstraint) {
    int failures = 0;
    for (int iter = 0; iter < ITERATIONS; ++iter) {
        
        // Build a system consisting of a chain of bodies with occasional side chains, and
        // a variety of mobilizers.

        MultibodySystem mbs;
        SimbodyMatterSubsystem matter(mbs);
        Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)))
                                      .addDecoration(Transform(), DecorativeSphere(.1));
        MobilizedBody* lastBody = &matter.Ground();
        MobilizedBody* lastMainChainBody = &matter.Ground();
        vector<MobilizedBody*> bodies;
        Random::Uniform random(0.0, 1.0);
        random.setSeed(iter);

        for (int i = 0; i < NUM_BODIES; ++i) {
            bool mainChain = random.getValue() < 0.5;
            MobilizedBody* parent = (mainChain ? lastMainChainBody : lastBody);
            int type = (int) (random.getValue()*4);
            MobilizedBody* nextBody;
            if (type == 0) {
                MobilizedBody::Cylinder cylinder(*parent, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(cylinder.getMobilizedBodyIndex());
            }
            else if (type == 1) {
                MobilizedBody::Slider slider(*parent, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(slider.getMobilizedBodyIndex());
            }
            else if (type == 2) {
                MobilizedBody::Ball ball(*parent, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(ball.getMobilizedBodyIndex());
            }
            else {
                MobilizedBody::Pin pin(*parent, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(pin.getMobilizedBodyIndex());
            }
            bodies.push_back(nextBody);
            if (mainChain)
                lastMainChainBody = nextBody;
            lastBody = nextBody;
        }
        mbs.realizeTopology();
        State s = mbs.getDefaultState();
        matter.setUseEulerAngles(s, true);
        mbs.realizeModel(s);

        // Choose a random initial conformation.

        vector<Real> targetQ(s.getNQ(), Real(0));
        for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
            const MobilizedBody& mobod = matter.getMobilizedBody(mbx);
            for (int i = 0; i < mobod.getNumQ(s); ++i) {
                const QIndex qx0 = mobod.getFirstQIndex(s);
                s.updQ()[qx0+i] = targetQ[qx0+i] = 2.0*random.getValue();
            }
        }
        //cout << "q0=" << s.getQ() << endl;
        mbs.realize(s, Stage::Position);

        // Select some random stations on each body.

        vector<vector<Vec3> > stations(NUM_BODIES);
        vector<vector<Vec3> > targetLocations(NUM_BODIES);
        vector<MobilizedBodyIndex> bodyIxs;
        for (int i = 0; i < NUM_BODIES; ++i) {
            MobilizedBodyIndex id = bodies[i]->getMobilizedBodyIndex();
            bodyIxs.push_back(id);
            int numStations = 1 + (int) (random.getValue()*4);
            for (int j = 0; j < numStations; ++j) {
                Vec3 pos(2.0*random.getValue()-1.0, 2.0*random.getValue()-1.0, 2.0*random.getValue()-1.0);
                stations[i].push_back(pos);
                targetLocations[i].push_back(bodies[i]->getBodyTransform(s)*pos);
            }
        }

        Real distance = -1;
        if (useConstraint) {
            // Add a constraint fixing the distance between the first and last bodies.
            Real distance = (bodies[0]->getBodyOriginLocation(s)-bodies[NUM_BODIES-1]->getBodyOriginLocation(s)).norm();
            Constraint::Rod(*bodies[0], Vec3(0), *bodies[NUM_BODIES-1], Vec3(0), distance);
        }
        s = mbs.realizeTopology();
        matter.setUseEulerAngles(s, true);
        mbs.realizeModel(s);

        // Try fitting it.
        State initState = s;
        if (!testFitting(mbs, s, bodyIxs, stations, targetLocations, 0.0, 0.02, distance))
            failures++;

        //cout << "q1=" << s.getQ() << endl;

        // Now add random noise to the target locations, and see if it can still fit decently.

        Random::Gaussian gaussian(0.0, 0.15);
        for (int i = 0; i < (int) targetLocations.size(); ++i) {
            for (int j = 0; j < (int) targetLocations[i].size(); ++j) {
                targetLocations[i][j] += Vec3(gaussian.getValue(), gaussian.getValue(), gaussian.getValue());
            }
        }

        s = initState; // start from same config as before
        if (!testFitting(mbs, s, bodyIxs, stations, targetLocations, 0.1, 0.5, distance))
            failures++;

        //cout << "q2=" << s.getQ() << endl;

    }

    ASSERT(failures == 0); // It found a reasonable fit every time.
    std::cout << "Done" << std::endl;
}

static void testUnconstrained() {
    testObservedPointFitter(false);
}

static void testConstrained() {
    testObservedPointFitter(true);
}

int main() {
    SimTK_START_TEST("TestObservedPointFitter");
        SimTK_SUBTEST(testUnconstrained);
        SimTK_SUBTEST(testConstrained);
    SimTK_END_TEST();
}
