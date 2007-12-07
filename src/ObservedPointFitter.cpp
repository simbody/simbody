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

class OptimizerFunction : public OptimizerSystem {
public:
    OptimizerFunction(const MultibodySystem& system, const State& state, vector<MobilizedBodyId> bodyIds, vector<vector<Vec3> > stations, vector<vector<Vec3> > targetLocations, vector<vector<Real> > weights) :
        OptimizerSystem(state.getNQ()), system(system), state(state), bodyIds(bodyIds), stations(stations), targetLocations(targetLocations), weights(weights) {
        setNumEqualityConstraints(state.getNQErr());
    }
    int objectiveFunc(const Vector& parameters, const bool new_parameters, Real& f) const {
        state.updQ() = parameters;
        system.realize(state, Stage::Position);
        f = 0.0;
        for (int i = 0; i < (int)bodyIds.size(); ++i) {
            const MobilizedBodyId id = bodyIds[i];
            const MobilizedBody& body = system.getMatterSubsystem().getMobilizedBody(id);
            for (int j = 0; j < (int)stations[id].size(); ++j) {
                f += weights[id][j]*(targetLocations[id][j]-body.getBodyTransform(state)*stations[id][j]).normSqr();
            }
        }
        return 0;
    }
    int constraintFunc(const Vector& parameters, const bool new_parameters, Vector& constraints) const {
        state.updQ() = parameters;
        system.realize(state, Stage::Velocity);
        constraints = state.getQErr();
        return 0;
    }
    void optimize(Vector& q) {
        Optimizer opt(*this);
//        opt.setConvergenceTolerance(0.01);
        opt.useNumericalGradient(true);
        opt.useNumericalJacobian(true);
//        opt.setMaxIterations(1000);
        opt.optimize(q);
    }
private:
    const MultibodySystem& system;
    const vector<MobilizedBodyId> bodyIds;
    const vector<vector<Vec3> > stations;
    const vector<vector<Vec3> > targetLocations;
    const vector<vector<Real> > weights;
    mutable State state;
};

/**
 * Create a new MultibodySystem which is identical to a subset of the original MultibodySystem.  This is called once for each MobilizedBody
 * in the original system, and is used to find an initial estimate of that MobilizedBody's conformation.
 */

void createClonedSystem(const MultibodySystem& original, MultibodySystem& copy, const vector<MobilizedBodyId>& originalBodyIds, vector<MobilizedBodyId>& copyBodyIds) {
    const SimbodyMatterSubsystem& originalMatter = original.getMatterSubsystem();
    SimbodyMatterSubsystem copyMatter(copy);
    Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)))
                                  .addDecoration(Transform(), DecorativeSphere(.1));
    map<MobilizedBodyId, MobilizedBodyId> idMap;
    for (int i = 0; i < (int)originalBodyIds.size(); ++i) {
        const MobilizedBody& originalBody = originalMatter.getMobilizedBody(originalBodyIds[i]);
        MobilizedBody* copyBody;
        if (i == 0) {
            if (originalBody.isGround())
                copyBody = &copyMatter.Ground();
            else
                copyBody = new MobilizedBody::Free(copyMatter.Ground(), body);
        }
        else {
            MobilizedBody& parent = copyMatter.updMobilizedBody(idMap[originalBody.getParentMobilizedBody().getMobilizedBodyId()]);
            copyBody = originalBody.cloneForNewParent(parent);
        }
        copyBodyIds.push_back(copyBody->getMobilizedBodyId());
        idMap[originalBodyIds[i]] = copyBody->getMobilizedBodyId();
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

void findUpstreamBodies(MobilizedBodyId currentBodyId, const vector<int> numStations, const SimbodyMatterSubsystem& matter, vector<MobilizedBodyId>& bodyIds, int requiredStations) /*const*/ {
    const MobilizedBody& currentBody = matter.getMobilizedBody(currentBodyId);
    if (currentBody.isGround())
        return;
    MobilizedBodyId parentId = currentBody.getParentMobilizedBody().getMobilizedBodyId();
    requiredStations -= numStations[parentId];
    if (requiredStations > 0)
        findUpstreamBodies(parentId, numStations, matter, bodyIds, requiredStations);
    bodyIds.push_back(parentId);
}

/**
 * This is invoked by findBodiesForClonedSystem().  It traces the tree downstream, adding bodies until there are sufficient stations
 * to reasonably perform a fit.
 */

void findDownstreamBodies(MobilizedBodyId currentBodyId, const vector<int> numStations, const vector<vector<MobilizedBodyId> > children, vector<MobilizedBodyId>& bodyIds, int& requiredStations) /*const*/ {
    if (numStations[currentBodyId] == 0 && children[currentBodyId].empty())
        return; // There's no benefit from including this body.
    bodyIds.push_back(currentBodyId);
    requiredStations -= numStations[currentBodyId];
    for (int i = 0; i < (int)children[currentBodyId].size() && requiredStations > 0; ++i) {
        findDownstreamBodies(children[currentBodyId][i], numStations, children, bodyIds, requiredStations);
    }
}

/**
 * Find the set of bodies that will be included in the system to be created by createClonedSystem().  The goal is to have sufficient
 * stations both upstream and downstream of the MobilizedBody currently being analyzed.
 */

int findBodiesForClonedSystem(MobilizedBodyId primaryBodyId, const vector<int> numStations, const SimbodyMatterSubsystem& matter, const vector<vector<MobilizedBodyId> > children, vector<MobilizedBodyId>& bodyIds) /*const*/ {
    findUpstreamBodies(primaryBodyId, numStations,  matter, bodyIds, 5);
    int primaryBodyIndex = bodyIds.size();
    int requiredStations = 5;
    findDownstreamBodies(primaryBodyId, numStations, children, bodyIds, requiredStations);
    return primaryBodyIndex;
}

Real ObservedPointFitter::findBestFit(const MultibodySystem& system, State& state, const vector<MobilizedBodyId>& bodyIds, const vector<vector<Vec3> >& stations, const vector<vector<Vec3> >& targetLocations) {
    vector<vector<Real> > weights(stations.size());
    for (int i = 0; i < (int)stations.size(); ++i)
        for (int j = 0; j < (int)stations[i].size(); ++j)
            weights[i].push_back(1.0);
    return findBestFit(system, state, bodyIds, stations, targetLocations, weights);
}

Real ObservedPointFitter::findBestFit(const MultibodySystem& system, State& state, const vector<MobilizedBodyId>& bodyIds, const vector<vector<Vec3> >& stations, const vector<vector<Vec3> >& targetLocations, const vector<vector<Real> >& weights) {

    // Build a list of children for each body.
    
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    vector<vector<MobilizedBodyId> > children(matter.getNBodies());
    for (int i = 0; i < matter.getNBodies(); ++i) {
        const MobilizedBody& body = matter.getMobilizedBody(MobilizedBodyId(i));
        if (!body.isGround())
            children[body.getParentMobilizedBody().getMobilizedBodyId()].push_back(body.getMobilizedBodyId());
    }
    
    // Find the number of stations on each body with a nonzero weight.
    
    vector<int> numStations(stations.size());
    for (int i = 0; i < (int)weights.size(); ++i) {
        numStations[i] = 0;
        for (int j = 0; j < (int)weights[i].size(); ++j)
            if (weights[i][j] != 0)
                numStations[i]++;
    }

    // Perform the initial estimation of Q for each mobilizer.
    
    State tempState = state;
    matter.setUseEulerAngles(tempState, true);
    system.realizeModel(tempState);
    tempState.updQ().setToZero();
    for (int i = 0; i < matter.getNBodies(); ++i) {
        MobilizedBodyId id(i);
        const MobilizedBody& body = matter.getMobilizedBody(id);
        if (body.getNumQ(tempState) == 0)
            continue; // No degrees of freedom to determine.
        if (children[id].size() == 0 && stations[id].size() == 0)
            continue; // There are no stations whose positions are affected by this.
        vector<MobilizedBodyId> originalBodyIds;
        int currentBodyIndex = findBodiesForClonedSystem(body.getMobilizedBodyId(), numStations, matter, children, originalBodyIds);
        if (currentBodyIndex == (int) originalBodyIds.size()-1 && stations[id].size() == 0)
            continue; // There are no stations whose positions are affected by this.
        MultibodySystem copy;
        vector<MobilizedBodyId> copyBodyIds;
        createClonedSystem(system, copy, originalBodyIds, copyBodyIds);
        vector<vector<Vec3> > copyStations(copy.getMatterSubsystem().getNBodies());
        vector<vector<Vec3> > copyTargetLocations(copy.getMatterSubsystem().getNBodies());
        vector<vector<Real> > copyWeights(copy.getMatterSubsystem().getNBodies());
        for (int j = 0; j < (int)originalBodyIds.size(); ++j) {
            copyStations[copyBodyIds[j]] = stations[originalBodyIds[j]];
            copyTargetLocations[copyBodyIds[j]] = targetLocations[originalBodyIds[j]];
            copyWeights[copyBodyIds[j]] = weights[originalBodyIds[j]];
        }
        OptimizerFunction optimizer(copy, copy.getDefaultState(), copyBodyIds, copyStations, copyTargetLocations, copyWeights);
        Vector q(copy.getDefaultState().getQ());
        optimizer.optimize(q);
        copy.updDefaultState().updQ() = q;
        body.setQVector(tempState, copy.getMatterSubsystem().getMobilizedBody(copyBodyIds[currentBodyIndex]).getQVector(copy.getDefaultState()));
    }

    // Now do the final optimization of the whole system.

    OptimizerFunction optimizer(system, tempState, bodyIds, stations, targetLocations, weights);
    Vector q = tempState.getQ();
    optimizer.optimize(q);
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
