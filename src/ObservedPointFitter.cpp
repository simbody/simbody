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

#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/ObservedPointFitter.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include <vector>
#include <map>

using namespace SimTK;
using namespace std;

/**
 * This class defines the objective function which is passed to the Optimizer.
 */

class ObservedPointFitter::OptimizerFunction : public OptimizerSystem {
public:
    OptimizerFunction(const MultibodySystem& system, const State& state, vector<MobilizedBodyIndex> bodyIxs, vector<vector<Vec3> > stations, vector<vector<Vec3> > targetLocations, vector<vector<Real> > weights) :
        OptimizerSystem(state.getNQ()), system(system), state(state), bodyIxs(bodyIxs), stations(stations), targetLocations(targetLocations), weights(weights) {
        system.realize(state, Stage::Instance);
        setNumEqualityConstraints(state.getNQErr());
    }
    int objectiveFunc(const Vector& parameters, const bool new_parameters, Real& f) const {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Position);
        f = 0.0;
        for (int i = 0; i < (int)bodyIxs.size(); ++i) {
            const MobilizedBodyIndex id = bodyIxs[i];
            const MobilizedBody& body = system.getMatterSubsystem().getMobilizedBody(id);
            for (int j = 0; j < (int)stations[i].size(); ++j) {
                f += weights[i][j]*(targetLocations[i][j]-body.getBodyTransform(state)*stations[i][j]).normSqr();
            }
        }
        return 0;
    }
    int gradientFunc(const Vector &parameters, const bool new_parameters, Vector &gradient) const  {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Position);
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        Vector_<SpatialVec> dEdR(matter.getNumBodies());
        dEdR = SpatialVec(Vec3(0), Vec3(0));
        for (int i = 0; i < (int)bodyIxs.size(); ++i) {
            const MobilizedBodyIndex id = bodyIxs[i];
            const MobilizedBody& body = matter.getMobilizedBody(id);
            for (int j = 0; j < (int)stations[i].size(); ++j) {
                Vec3 f = 2.0*weights[i][j]*(body.getBodyTransform(state)*stations[i][j]-targetLocations[i][j]);
                body.applyForceToBodyPoint(state, stations[i][j], f, dEdR);
            }
        }
        Vector dEdU;
        matter.calcInternalGradientFromSpatial(state, dEdR, dEdU);
        matter.multiplyByNInv(state, true, dEdU, gradient);
        return 0;
    }
    int constraintFunc(const Vector& parameters, const bool new_parameters, Vector& constraints) const {
        state.updQ() = parameters;
        system.realize(state, Stage::Velocity);
        constraints = state.getQErr();
        return 0;
    }
    void optimize(Vector& q, Real tolerance) {
        Optimizer opt(*this);
        opt.useNumericalJacobian(true);
        opt.setConvergenceTolerance(tolerance);
        opt.setLimitedMemoryHistory(100);
        opt.optimize(q);
    }
private:
    const MultibodySystem& system;
    const vector<MobilizedBodyIndex> bodyIxs;
    const vector<vector<Vec3> > stations;
    const vector<vector<Vec3> > targetLocations;
    const vector<vector<Real> > weights;
    mutable State state;
};

/**
 * Create a new MultibodySystem which is identical to a subset of the original MultibodySystem.  This is called once for each MobilizedBody
 * in the original system, and is used to find an initial estimate of that MobilizedBody's conformation.
 */

void ObservedPointFitter::createClonedSystem(const MultibodySystem& original, MultibodySystem& copy, const vector<MobilizedBodyIndex>& originalBodyIxs, vector<MobilizedBodyIndex>& copyBodyIxs) {
    const SimbodyMatterSubsystem& originalMatter = original.getMatterSubsystem();
    SimbodyMatterSubsystem copyMatter(copy);
    Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)))
                                  .addDecoration(Transform(), DecorativeSphere(.1));
    map<MobilizedBodyIndex, MobilizedBodyIndex> idMap;
    for (int i = 0; i < (int)originalBodyIxs.size(); ++i) {
        const MobilizedBody& originalBody = originalMatter.getMobilizedBody(originalBodyIxs[i]);
        MobilizedBody* copyBody;
        if (i == 0) {
            if (originalBody.isGround())
                copyBody = &copyMatter.Ground();
            else {
                MobilizedBody::Free free(copyMatter.Ground(), body);
                copyBody = &copyMatter.updMobilizedBody(free.getMobilizedBodyIndex());
            }
        }
        else {
            MobilizedBody& parent = copyMatter.updMobilizedBody(idMap[originalBody.getParentMobilizedBody().getMobilizedBodyIndex()]);
            copyBody = &originalBody.cloneForNewParent(parent);
        }
        copyBodyIxs.push_back(copyBody->getMobilizedBodyIndex());
        idMap[originalBodyIxs[i]] = copyBody->getMobilizedBodyIndex();
    }
    copy.realizeTopology();
    State& s = copy.updDefaultState();
    copyMatter.setUseEulerAngles(s, true);
    copy.realizeModel(s);
}

/**
 * This is invoked by findBodiesForClonedSystem().  It traces the tree upstream, adding bodies until there are sufficient stations
 * to reasonably perform a fit.
 */

void ObservedPointFitter::findUpstreamBodies(MobilizedBodyIndex currentBodyIx, const vector<int> numStations, const SimbodyMatterSubsystem& matter, vector<MobilizedBodyIndex>& bodyIxs, int requiredStations) {
    const MobilizedBody& currentBody = matter.getMobilizedBody(currentBodyIx);
    if (currentBody.isGround())
        return;
    MobilizedBodyIndex parentIx = currentBody.getParentMobilizedBody().getMobilizedBodyIndex();
    requiredStations -= numStations[parentIx];
    if (requiredStations > 0)
        findUpstreamBodies(parentIx, numStations, matter, bodyIxs, requiredStations);
    bodyIxs.push_back(parentIx);
}

/**
 * This is invoked by findBodiesForClonedSystem().  It traces the tree downstream, adding bodies until there are sufficient stations
 * to reasonably perform a fit.
 */

void ObservedPointFitter::findDownstreamBodies(MobilizedBodyIndex currentBodyIx, const vector<int> numStations, const vector<vector<MobilizedBodyIndex> > children, vector<MobilizedBodyIndex>& bodyIxs, int& requiredStations) {
    if (numStations[currentBodyIx] == 0 && children[currentBodyIx].empty())
        return; // There's no benefit from including this body.
    bodyIxs.push_back(currentBodyIx);
    requiredStations -= numStations[currentBodyIx];
    for (int i = 0; i < (int)children[currentBodyIx].size() && requiredStations > 0; ++i) {
        findDownstreamBodies(children[currentBodyIx][i], numStations, children, bodyIxs, requiredStations);
    }
}

/**
 * Find the set of bodies that will be included in the system to be created by createClonedSystem().  The goal is to have sufficient
 * stations both upstream and downstream of the MobilizedBody currently being analyzed.
 */

int ObservedPointFitter::findBodiesForClonedSystem(MobilizedBodyIndex primaryBodyIx, const vector<int> numStations, const SimbodyMatterSubsystem& matter, const vector<vector<MobilizedBodyIndex> > children, vector<MobilizedBodyIndex>& bodyIxs) {
    findUpstreamBodies(primaryBodyIx, numStations,  matter, bodyIxs, 5);
    int primaryBodyIndex = bodyIxs.size();
    int requiredStations = 5;
    findDownstreamBodies(primaryBodyIx, numStations, children, bodyIxs, requiredStations);
    return primaryBodyIndex;
}

Real ObservedPointFitter::findBestFit(const MultibodySystem& system, State& state, const vector<MobilizedBodyIndex>& bodyIxs, const vector<vector<Vec3> >& stations, const vector<vector<Vec3> >& targetLocations, Real tolerance) {
    vector<vector<Real> > weights(stations.size());
    for (int i = 0; i < (int)stations.size(); ++i)
        for (int j = 0; j < (int)stations[i].size(); ++j)
            weights[i].push_back(1.0);
    return findBestFit(system, state, bodyIxs, stations, targetLocations, weights, tolerance);
}

Real ObservedPointFitter::findBestFit(const MultibodySystem& system, State& state, const vector<MobilizedBodyIndex>& bodyIxs, const vector<vector<Vec3> >& stations, const vector<vector<Vec3> >& targetLocations, const vector<vector<Real> >& weights, Real tolerance) {
    
    // Verify the inputs.
    
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    SimTK_APIARGCHECK(bodyIxs.size() == stations.size() && stations.size() == targetLocations.size(), "ObservedPointFitter", "findBestFit", "bodyIxs, stations, and targetLocations must all be the same length");
    int numBodies = matter.getNumBodies();
    for (int i = 0; i < (int)stations.size(); ++i) {
        SimTK_APIARGCHECK(bodyIxs[i] >= 0 && bodyIxs[i] < numBodies, "ObservedPointFitter", "findBestFit", "Illegal body ID");
        SimTK_APIARGCHECK(stations[i].size() == targetLocations[i].size(), "ObservedPointFitter", "findBestFit", "Different number of stations and target locations for body");
    }
    
    // Build a list of children for each body.
    
    vector<vector<MobilizedBodyIndex> > children(matter.getNumBodies());
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        const MobilizedBody& body = matter.getMobilizedBody(MobilizedBodyIndex(i));
        if (!body.isGround())
            children[body.getParentMobilizedBody().getMobilizedBodyIndex()].push_back(body.getMobilizedBodyIndex());
    }

    // Build a mapping of body IDs to indices.
    
    vector<int> bodyIndex(matter.getNumBodies());
    for (int i = 0; i < (int) bodyIndex.size(); ++i)
        bodyIndex[i] = -1;
    for (int i = 0; i < (int)bodyIxs.size(); ++i)
        bodyIndex[bodyIxs[i]] = i;
    
    // Find the number of stations on each body with a nonzero weight.
    
    vector<int> numStations(matter.getNumBodies());
    for (int i = 0; i < (int) numStations.size(); ++i)
        numStations[i] = 0;
    for (int i = 0; i < (int)weights.size(); ++i) {
        for (int j = 0; j < (int)weights[i].size(); ++j)
            if (weights[i][j] != 0)
                numStations[bodyIxs[i]]++;
    }

    // Perform the initial estimation of Q for each mobilizer.
    
    State tempState = state;
    matter.setUseEulerAngles(tempState, true);
    system.realizeModel(tempState);
    tempState.updQ().setToZero();
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        MobilizedBodyIndex id(i);
        const MobilizedBody& body = matter.getMobilizedBody(id);
        if (body.getNumQ(tempState) == 0)
            continue; // No degrees of freedom to determine.
        if (children[id].size() == 0 && numStations[id] == 0)
            continue; // There are no stations whose positions are affected by this.
        vector<MobilizedBodyIndex> originalBodyIxs;
        int currentBodyIndex = findBodiesForClonedSystem(body.getMobilizedBodyIndex(), numStations, matter, children, originalBodyIxs);
        if (currentBodyIndex == (int)originalBodyIxs.size()-1 
            && (bodyIndex[id] == -1 || stations[bodyIndex[id]].size() == 0))
            continue; // There are no stations whose positions are affected by this.
        MultibodySystem copy;
        vector<MobilizedBodyIndex> copyBodyIxs;
        createClonedSystem(system, copy, originalBodyIxs, copyBodyIxs);
        vector<vector<Vec3> > copyStations(copy.getMatterSubsystem().getNumBodies());
        vector<vector<Vec3> > copyTargetLocations(copy.getMatterSubsystem().getNumBodies());
        vector<vector<Real> > copyWeights(copy.getMatterSubsystem().getNumBodies());
        for (int j = 0; j < (int)originalBodyIxs.size(); ++j) {
            int index = bodyIndex[originalBodyIxs[j]];
            if (index != -1) {
                copyStations[copyBodyIxs[j]] = stations[index];
                copyTargetLocations[copyBodyIxs[j]] = targetLocations[index];
                copyWeights[copyBodyIxs[j]] = weights[index];
            }
        }
        try {
            OptimizerFunction optimizer(copy, copy.getDefaultState(), copyBodyIxs, copyStations, copyTargetLocations, copyWeights);
            Vector q(copy.getDefaultState().getQ());
            optimizer.optimize(q, tolerance);
            copy.updDefaultState().updQ() = q;
            body.setQFromVector(tempState, copy.getMatterSubsystem().getMobilizedBody(copyBodyIxs[currentBodyIndex]).getQAsVector(copy.getDefaultState()));
        }
        catch (Exception::OptimizerFailed ex) {
            std::cout << "Optimization failure for body "<<i<<": "<<ex.getMessage() << std::endl;
            // Just leave this body's state variables set to 0, and rely on the final optimization to fix them.
        }
    }

    // Now do the final optimization of the whole system.

    OptimizerFunction optimizer(system, tempState, bodyIxs, stations, targetLocations, weights);
    Vector q = tempState.getQ();
    optimizer.optimize(q, tolerance);
    if (matter.getUseEulerAngles(state))
        state.updQ() = q;
    else {
        tempState.updQ() = q;
        matter.convertToQuaternions(tempState, state);
    }
    
    // Return the RMS error in the optimized system.
    
    Real error;
    optimizer.objectiveFunc(q, true, error);
    Real totalWeight = 0;
    for (int i = 0; i < (int)weights.size(); ++i)
        for (int j = 0; j < (int)weights[i].size(); ++j)
            totalWeight += weights[i][j];
    return std::sqrt(error/totalWeight);
}
