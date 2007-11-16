/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include <SimTKsimbody.h>
#include <vector>
#include <map>

using namespace SimTK;
using namespace std;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

const int NUM_BODIES = 20;
const Real BOND_LENGTH = 0.5;

void testFitting(const MultibodySystem& mbs, State& state, const vector<MobilizedBodyId>& bodyIds, const vector<vector<Vec3> >& stations, const vector<vector<Vec3> >& targetLocations, Real minError, Real maxError, Real endDistance) {
    
    // Find the best fit.
    
    Real reportedError = ObservedPointFitter::findBestFit(mbs, state, bodyIds, stations, targetLocations);
    ASSERT(reportedError <= maxError && reportedError >= minError);
    
    // Verify that the error was calculated correctly.
    
    Real error = 0.0;
    int numStations = 0;
    mbs.realize(state, Stage::Position);
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    for (int i = 0; i < (int) bodyIds.size(); ++i) {
        MobilizedBodyId id = bodyIds[i];
        numStations += stations[id].size();
        for (int j = 0; j < (int) stations[id].size(); ++j)
            error += (targetLocations[id][j]-matter.getMobilizedBody(id).getBodyTransform(state)*stations[id][j]).normSqr();
    }
    error = std::sqrt(error/numStations);
    ASSERT(std::abs(1.0-error/reportedError) < 0.0001);
    
    // Verify that the ends are the correct distance apart.
    
    Real distance = (matter.getMobilizedBody(bodyIds[0]).getBodyOriginLocation(state)-matter.getMobilizedBody(bodyIds[bodyIds.size()-1]).getBodyOriginLocation(state)).norm();
    ASSERT(std::abs(1.0-endDistance/distance) < 0.0001);
}


int main() {
    
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
    for (int i = 0; i < NUM_BODIES; ++i) {
        bool mainChain = random.getValue() < 0.5;
        MobilizedBody* parent = (mainChain ? lastMainChainBody : lastBody);
        int type = (int) (random.getValue()*3);
        MobilizedBody* nextBody;
        if (type == 0)
            nextBody = new MobilizedBody::Cylinder(*parent, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, BOND_LENGTH, 0)));
        else if (type == 1)
            nextBody = new MobilizedBody::Slider(*parent, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, BOND_LENGTH, 0)));
        else
            nextBody = new MobilizedBody::Ball(*parent, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, BOND_LENGTH, 0)));
        bodies.push_back(nextBody);
        if (mainChain)
            lastMainChainBody = nextBody;
        lastBody = nextBody;
    }
    mbs.realizeTopology();
    State& s = mbs.updDefaultState();
    mbs.realizeModel(s);
    
    // Choose a random initial conformation.
    
    vector<Real> targetQ(s.getNQ());
    for (int i = 0; i < s.getNQ(); ++i)
        s.updQ()[i] = targetQ[i] = 2.0*random.getValue();
    mbs.realize(s, Stage::Position);
    
    // Select some random stations on each body.
    
    vector<vector<Vec3> > stations(matter.getNBodies());
    vector<vector<Vec3> > targetLocations(matter.getNBodies());
    vector<MobilizedBodyId> bodyIds;
    for (int i = 0; i < NUM_BODIES; ++i) {
        MobilizedBodyId id = bodies[i]->getMobilizedBodyId();
        bodyIds.push_back(id);
        int numStations = (int) (random.getValue()*4);
        for (int j = 0; j < numStations; ++j) {
            Vec3 pos(2.0*random.getValue()-1.0, 2.0*random.getValue()-1.0, 2.0*random.getValue()-1.0);
            stations[id].push_back(pos);
            targetLocations[id].push_back(bodies[i]->getBodyTransform(s)*pos);
        }
    }
    
    // Add a constraint fixing the distance between the first and last bodies.
    
    Real distance = (bodies[0]->getBodyOriginLocation(s)-bodies[NUM_BODIES-1]->getBodyOriginLocation(s)).norm();
    Constraint::Rod(*bodies[0], Vec3(0), *bodies[NUM_BODIES-1], Vec3(0), distance);
    mbs.realizeTopology();
    
    // Try fitting it.
    
    testFitting(mbs, s, bodyIds, stations, targetLocations, 0.0, 0.02, distance);
    
    // Now add random noise to the target locations, and see if it can still fit decently.
    
    Random::Gaussian gaussian(0.0, 0.2);
    for (int i = 0; i < (int) targetLocations.size(); ++i) {
        for (int j = 0; j < (int) targetLocations[i].size(); ++j) {
            targetLocations[i][j] += Vec3(gaussian.getValue(), gaussian.getValue(), gaussian.getValue());
        }
    }
    testFitting(mbs, s, bodyIds, stations, targetLocations, 0.1, 0.4, distance);
    std::cout << "Done" << std::endl;
}
