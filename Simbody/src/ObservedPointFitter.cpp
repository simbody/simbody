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

#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MobilizedBody_Ground.h"
#include "simbody/internal/MobilizedBody_Free.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/ObservedPointFitter.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include <map>

using namespace SimTK;

/**
 * This class defines the objective function which is passed to the Optimizer.
 */

static const bool UseWeighted = true;
static const Real MinimumShift = 1; // add to objective to get minimum away from zero

class ObservedPointFitter::OptimizerFunction : public OptimizerSystem {
public:
    OptimizerFunction(const MultibodySystem& system, const State& state, Array_<MobilizedBodyIndex> bodyIxs, Array_<Array_<Vec3> > stations, Array_<Array_<Vec3> > targetLocations, Array_<Array_<Real> > weights) :
        OptimizerSystem(state.getNQ()), system(system), state(state), bodyIxs(bodyIxs), stations(stations), targetLocations(targetLocations), weights(weights) {
        system.realize(state, Stage::Instance);
        setNumEqualityConstraints(state.getNQErr());
    }
    int objectiveFunc(const Vector& parameters, bool new_parameters, Real& f) const {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Position);
        Real wtot = 0;
        f = 0;
        for (int i = 0; i < (int)bodyIxs.size(); ++i) {
            const MobilizedBodyIndex id = bodyIxs[i];
            const MobilizedBody& body = system.getMatterSubsystem().getMobilizedBody(id);
            for (int j = 0; j < (int)stations[i].size(); ++j) {
                f += weights[i][j]*(targetLocations[i][j]-body.getBodyTransform(state)*stations[i][j]).normSqr();
                wtot += weights[i][j];
            }
        }
        if (UseWeighted && wtot > 0) f /= wtot;

        f += MinimumShift; // so minimum won't be at zero where scaling is tricky
        return 0;
    }
    int gradientFunc(const Vector &parameters, bool new_parameters, Vector &gradient) const  {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Position);
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        Vector_<SpatialVec> dEdR(matter.getNumBodies());
        dEdR = SpatialVec(Vec3(0), Vec3(0));
        Real wtot = 0;
        for (int i = 0; i < (int)bodyIxs.size(); ++i) {
            const MobilizedBodyIndex id = bodyIxs[i];
            const MobilizedBody& body = matter.getMobilizedBody(id);
            for (int j = 0; j < (int)stations[i].size(); ++j) {
                Vec3 f = 2*weights[i][j]*(body.getBodyTransform(state)*stations[i][j]-targetLocations[i][j]);
                body.applyForceToBodyPoint(state, stations[i][j], f, dEdR);
                wtot += weights[i][j];
            }
        }
        Vector dEdU;
        // Convert spatial forces dEdR to generalized forces dEdU.
        matter.multiplyBySystemJacobianTranspose(state, dEdR, dEdU);
        if (UseWeighted && wtot > 0) dEdU /= wtot;
        matter.multiplyByNInv(state, true, dEdU, gradient);
        return 0;
    }
    int constraintFunc(const Vector& parameters, bool new_parameters, Vector& constraints) const {
        state.updQ() = parameters;
        system.realize(state, Stage::Velocity);
        constraints = state.getQErr();
        return 0;
    }
    void optimize(Vector& q, Real tolerance) {
        Optimizer opt(*this
            //, LBFGSB // XXX
            //, InteriorPoint // XXX
            );
        //opt.useNumericalGradient(true); //XXX
        opt.useNumericalJacobian(true);
        opt.setConvergenceTolerance(tolerance);
        opt.setMaxIterations(3000);
        opt.setLimitedMemoryHistory(40);
        //opt.setDiagnosticsLevel(5);
        opt.optimize(q);
    }
private:
    const MultibodySystem& system;
    const Array_<MobilizedBodyIndex> bodyIxs;
    const Array_<Array_<Vec3> > stations;
    const Array_<Array_<Vec3> > targetLocations;
    const Array_<Array_<Real> > weights;
    mutable State state;
};

/**
 * Create a new MultibodySystem which is identical to a subset of the original MultibodySystem.  This is called once for each MobilizedBody
 * in the original system, and is used to find an initial estimate of that MobilizedBody's conformation.
 */

void ObservedPointFitter::
createClonedSystem(const MultibodySystem& original, MultibodySystem& copy,
                   const Array_<MobilizedBodyIndex>& originalBodyIxs,
                   Array_<MobilizedBodyIndex>& copyBodyIxs,
                   bool& hasArtificialBaseBody)
{
    const SimbodyMatterSubsystem& originalMatter = original.getMatterSubsystem();
    SimbodyMatterSubsystem copyMatter(copy);
    Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)));
    body.addDecoration(Transform(), DecorativeSphere(Real(.1)));
    std::map<MobilizedBodyIndex, MobilizedBodyIndex> idMap;
    hasArtificialBaseBody = false;
    for (int i = 0; i < (int)originalBodyIxs.size(); ++i) {
        const MobilizedBody& originalBody = originalMatter.getMobilizedBody(originalBodyIxs[i]);
        MobilizedBody* copyBody;
        if (i == 0) {
            if (originalBody.isGround())
                copyBody = &copyMatter.Ground();
            else {
                hasArtificialBaseBody = true; // not using the original joint here
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

void ObservedPointFitter::findUpstreamBodies(MobilizedBodyIndex currentBodyIx, const Array_<int> numStations, const SimbodyMatterSubsystem& matter, Array_<MobilizedBodyIndex>& bodyIxs, int requiredStations) {
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

void ObservedPointFitter::findDownstreamBodies(MobilizedBodyIndex currentBodyIx, const Array_<int> numStations, const Array_<Array_<MobilizedBodyIndex> > children, Array_<MobilizedBodyIndex>& bodyIxs, int& requiredStations) {
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

int ObservedPointFitter::findBodiesForClonedSystem(MobilizedBodyIndex primaryBodyIx, const Array_<int> numStations, const SimbodyMatterSubsystem& matter, const Array_<Array_<MobilizedBodyIndex> > children, Array_<MobilizedBodyIndex>& bodyIxs) {
    findUpstreamBodies(primaryBodyIx, numStations,  matter, bodyIxs, 5);
    int primaryBodyIndex = bodyIxs.size();
    int requiredStations = 5;
    findDownstreamBodies(primaryBodyIx, numStations, children, bodyIxs, requiredStations);
    return primaryBodyIndex;
}

Real ObservedPointFitter::findBestFit
   (const MultibodySystem& system, State& state,
    const Array_<MobilizedBodyIndex>&  bodyIxs,
    const Array_<Array_<Vec3> >&       stations,
    const Array_<Array_<Vec3> >&       targetLocations,
    Real                                        tolerance)
{
    Array_<Array_<Real> > weights(stations.size());
    for (int i = 0; i < (int)stations.size(); ++i)
        for (int j = 0; j < (int)stations[i].size(); ++j)
            weights[i].push_back(1.0);
    return findBestFit(system, state, bodyIxs, stations, targetLocations, weights, tolerance);
}

Real ObservedPointFitter::findBestFit
   (const MultibodySystem& system, State& state,
    const Array_<MobilizedBodyIndex>&  bodyIxs,
    const Array_<Array_<Vec3> >&       stations,
    const Array_<Array_<Vec3> >&       targetLocations,
    const Array_<Array_<Real> >&       weights,
    Real tolerance)
{
    // Verify the inputs.

    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    SimTK_APIARGCHECK(bodyIxs.size() == stations.size() && stations.size() == targetLocations.size(), "ObservedPointFitter", "findBestFit", "bodyIxs, stations, and targetLocations must all be the same length");
    int numBodies = matter.getNumBodies();
    for (int i = 0; i < (int)stations.size(); ++i) {
        SimTK_APIARGCHECK(bodyIxs[i] >= 0 && bodyIxs[i] < numBodies, "ObservedPointFitter", "findBestFit", "Illegal body ID");
        SimTK_APIARGCHECK(stations[i].size() == targetLocations[i].size(), "ObservedPointFitter", "findBestFit", "Different number of stations and target locations for body");
    }

    // Build a list of children for each body.

    Array_<Array_<MobilizedBodyIndex> > children(matter.getNumBodies());
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        const MobilizedBody& body = matter.getMobilizedBody(MobilizedBodyIndex(i));
        if (!body.isGround())
            children[body.getParentMobilizedBody().getMobilizedBodyIndex()].push_back(body.getMobilizedBodyIndex());
    }

    // Build a mapping of body IDs to indices.

    Array_<int> bodyIndex(matter.getNumBodies());
    for (int i = 0; i < (int) bodyIndex.size(); ++i)
        bodyIndex[i] = -1;
    for (int i = 0; i < (int)bodyIxs.size(); ++i)
        bodyIndex[bodyIxs[i]] = i;

    // Find the number of stations on each body with a nonzero weight.

    Array_<int> numStations(matter.getNumBodies());
    for (int i = 0; i < (int) numStations.size(); ++i)
        numStations[i] = 0;
    for (int i = 0; i < (int)weights.size(); ++i) {
        for (int j = 0; j < (int)weights[i].size(); ++j)
            if (weights[i][j] != 0)
                numStations[bodyIxs[i]]++;
    }

    // Perform the initial estimation of Q for each mobilizer.
    // Our first guess is the passed-in q's, with quaternions converted
    // to Euler angles if necessary. As we solve a subproblem for each
    // of the bodies in ascending order, we'll update tempState's q's
    // for that body to their solved values.
    State tempState;
    if (!matter.getUseEulerAngles(state))
        matter.convertToEulerAngles(state, tempState);
    else tempState = state;
    system.realizeModel(tempState);
    system.realize(tempState, Stage::Position);

    // This will accumulate best-guess spatial poses for the bodies as
    // they are processed. This is useful for when a body is used as
    // an artificial base body; our first guess will to be to place it
    // wherever it was the last time it was used in a subproblem.
    Array_<Transform> guessX_GB(matter.getNumBodies());
    for (MobilizedBodyIndex mbx(1); mbx < guessX_GB.size(); ++mbx)
        guessX_GB[mbx] = matter.getMobilizedBody(mbx).getBodyTransform(tempState);

    for (int i = 0; i < matter.getNumBodies(); ++i) {
        MobilizedBodyIndex id(i);
        const MobilizedBody& body = matter.getMobilizedBody(id);
        if (body.getNumQ(tempState) == 0)
            continue; // No degrees of freedom to determine.
        if (children[id].size() == 0 && numStations[id] == 0)
            continue; // There are no stations whose positions are affected by this.
        Array_<MobilizedBodyIndex> originalBodyIxs;
        int currentBodyIndex = findBodiesForClonedSystem(body.getMobilizedBodyIndex(), numStations, matter, children, originalBodyIxs);
        if (currentBodyIndex == (int)originalBodyIxs.size()-1
            && (bodyIndex[id] == -1 || stations[bodyIndex[id]].size() == 0))
            continue; // There are no stations whose positions are affected by this.
        MultibodySystem copy;
        Array_<MobilizedBodyIndex> copyBodyIxs;
        bool hasArtificialBaseBody;
        createClonedSystem(system, copy, originalBodyIxs, copyBodyIxs, hasArtificialBaseBody);
        const SimbodyMatterSubsystem& copyMatter = copy.getMatterSubsystem();
        // Construct an initial state.
        State copyState = copy.getDefaultState();
        assert(copyBodyIxs.size() == originalBodyIxs.size());
        for (int ob=0; ob < (int)originalBodyIxs.size(); ++ob) {
            const MobilizedBody& copyMobod = copyMatter.getMobilizedBody(copyBodyIxs[ob]);
            const MobilizedBody& origMobod = matter.getMobilizedBody(originalBodyIxs[ob]);
            if (ob==0 && hasArtificialBaseBody)
                copyMobod.setQToFitTransform(copyState, guessX_GB[origMobod.getMobilizedBodyIndex()]);
            else
                copyMobod.setQFromVector(copyState, origMobod.getQAsVector(tempState));
        }

        Array_<Array_<Vec3> > copyStations(copyMatter.getNumBodies());
        Array_<Array_<Vec3> > copyTargetLocations(copyMatter.getNumBodies());
        Array_<Array_<Real> > copyWeights(copyMatter.getNumBodies());
        for (int j = 0; j < (int)originalBodyIxs.size(); ++j) {
            int index = bodyIndex[originalBodyIxs[j]];
            if (index != -1) {
                copyStations[copyBodyIxs[j]] = stations[index];
                copyTargetLocations[copyBodyIxs[j]] = targetLocations[index];
                copyWeights[copyBodyIxs[j]] = weights[index];
            }
        }
        try {
            OptimizerFunction optimizer(copy, copyState, copyBodyIxs, copyStations, copyTargetLocations, copyWeights);
            Vector q(copyState.getQ());
            //std::cout << "BODY " << i << " q0=" << q << std::endl;
            optimizer.optimize(q, tolerance);
            //std::cout << "  qf=" << q << std::endl;
            copyState.updQ() = q;
            copy.realize(copyState, Stage::Position);
            // Transfer updated state back to tempState as improved initial guesses.
            // However, all but the currentBody will get overwritten later.
            for (int ob=0; ob < (int)originalBodyIxs.size(); ++ob) {
                const MobilizedBody& copyMobod = copyMatter.getMobilizedBody(copyBodyIxs[ob]);
                guessX_GB[originalBodyIxs[ob]] = copyMobod.getBodyTransform(copyState);

                if (ob==0 && hasArtificialBaseBody) continue; // leave default state
                const MobilizedBody& origMobod = matter.getMobilizedBody(originalBodyIxs[ob]);
                origMobod.setQFromVector(tempState, copyMobod.getQAsVector(copyState));
            }
            //body.setQFromVector(tempState, copyMatter.getMobilizedBody(copyBodyIxs[currentBodyIndex]).getQAsVector(copyState));
        }
        catch (const Exception::OptimizerFailed& ex) {
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
    if (UseWeighted)
        return std::sqrt(error - MinimumShift); // already weighted; this makes WRMS

    Real totalWeight = 0;
    for (int i = 0; i < (int)weights.size(); ++i)
        for (int j = 0; j < (int)weights[i].size(); ++j)
            totalWeight += weights[i][j];

    return std::sqrt((error-MinimumShift)/totalWeight);
}
