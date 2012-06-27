/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/**@file
 *
 * Implementation of SimbodyMatterSubtree and SimbodyMatterSubtreeResults.
 */

#include "SimTKcommon.h"
#include "simbody/internal/SimbodyMatterSubtree.h"

#include "MobilizedBodyImpl.h"
#include "SimbodyMatterSubsystemRep.h"
class RigidBodyNode;

#include <string>
#include <iostream>
using std::cout;
using std::endl;

namespace SimTK {

    /////////////////
    // SUBTREE REP //
    /////////////////

class SimbodyMatterSubtree::SubtreeRep {
public:
    SubtreeRep(const SimbodyMatterSubtree& handle, const SimbodyMatterSubsystem& sms) 
      : myHandle(handle), matter(&sms), stage(Stage::Empty)
    { 
    }

    // Here we don't know the matter subsystem yet.
    explicit SubtreeRep(const SimbodyMatterSubtree& handle) 
      : myHandle(handle), matter(0), stage(Stage::Empty)
    { 
    }

    void setSimbodyMatterSubsystem(const SimbodyMatterSubsystem& sms) {
        clear();
        matter = &sms; // just a reference
    }

    const SimbodyMatterSubsystem& getSimbodyMatterSubsystem() const {
        assert(matter != 0);
        return *matter;
    }

    void invalidate(Stage g) {
        if (stage >= g)
            stage = g.prev();
    }

    // Note that this retains the current handle and matter subsystem (if any).
    void clear() {
        invalidate(Stage::Topology);
        terminalBodies.clear();
        ancestor = InvalidMobilizedBodyIndex;
    }


    void setTerminalBodies(const Array_<MobilizedBodyIndex>& bids) {
        clear();
        for (int i=0; i < (int)bids.size(); ++i)
            addTerminalBody(bids[i]);
    }

    void addTerminalBody(MobilizedBodyIndex bid) {
        assert(!isTerminalBody(bid)); // can only appear once
        invalidate(Stage::Topology);
        terminalBodies.push_back(bid);
    }

    MobilizedBodyIndex getAncestorMobilizedBodyIndex() const {
        assert(stage >= Stage::Topology);
        return ancestor;
    }

    MobilizedBodyIndex getSubtreeBodyMobilizedBodyIndex(SubtreeBodyIndex b) const {
        assert(stage >= Stage::Topology);
        return allBodies[b];
    }

    SubtreeBodyIndex getParentSubtreeBodyIndex(SubtreeBodyIndex b) const {
        assert(b >= 1); // ancestor has no subtree parent
        assert(stage >= Stage::Topology);
        return parentSubtreeBodies[b];
    }

    // State must be realized to at least Stage::Model for this call to work. 
    // The supplied SubtreeResults object is allocated and properly initialized to
    // be able to hold computation results from this Subtree.
    void initializeSubtreeResults(const State&, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // This can be used as a sanity check that initializeSubtreeResults() was already called
    // in this Subtree to produce these SubtreeResults. It is by no means exhaustive but
    // will catch egregious errors.
    bool isCompatibleSubtreeResults(const SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

        // POSITION STAGE

    // State must be realized to at least Stage::Position for this to work. SubtreeResults
    // must have already been initialized to work with this Subtree. SubtreeResults stage
    // will be Stage::Position after this call. All body transforms will be the same as
    // the corresponding ones in the state, except they will be measured from the ancestor
    // frame instead of ground. Subtree q's will be identical to corresponding State q's.
    void copyPositionsFromState(const State&, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // State must be realized to Stage::Instance. subQ must be the right length for this
    // Subtree, and SubtreeResults must have been properly initialized. SubtreeResults
    // stage will be Stage::Position after this call.
    void calcPositionsFromSubtreeQ(const State&, const Vector& subQ, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // Calculates a perturbed position result starting with the subQ's and position results
    // which must already be in SubtreeResults.
    void perturbPositions(const State&, SubtreeQIndex subQIndex, Real perturbation, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;


        // VELOCITY STAGE

    // State must be realized to at least Stage::Velocity for this to work. SubtreeResults
    // must already be at Stage::Position. SubtreeResults stage
    // will be Stage::Velocity after this call. All subtree body spatial velocities will be
    // the same as in the State, except measured relative to A and expressed in A. Subtree u's
    // will be identical to corresponding State u's.
    void copyVelocitiesFromState(const State&, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // State must be realized to Stage::Instance. subU must be the right length for this
    // Subtree, and SubtreeResults must already be at Stage::Position. SubtreeResults
    // stage will be Stage::Velocity after this call.
    void calcVelocitiesFromSubtreeU(const State&, const Vector& subU, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // State must be realized to Stage::Instance and SubtreeResults must already be at
    // Stage::Position. SubtreeResults stage will be Stage::Velocity after this call, but
    // all Subtree u's and body velocities will be zero.
    void calcVelocitiesFromZeroU(const State&, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // Calculates a perturbed velocity result starting with the subU's and velocity results
    // which must already be in SubtreeResults.
    void perturbVelocities(const State&, SubtreeUIndex subUIndex, Real perturbation, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;


        // ACCELERATION STAGE

    // State must be realized to at least Stage::Acceleration for this to work. SubtreeResults
    // must already be at Stage::Velocity. SubtreeResults stage
    // will be Stage::Acceleration after this call. All subtree body spatial accelerations will be
    // the same as in the State, except measured relative to A and expressed in A. Subtree udots
    // will be identical to corresponding State udots.
    void copyAccelerationsFromState(const State&, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // State must be realized to Stage::Instance. subUDot must be the right length for this
    // Subtree, and SubtreeResults must already be at Stage::Velocity. SubtreeResults
    // stage will be Stage::Acceleration after this call.
    void calcAccelerationsFromSubtreeUDot(const State&, const Vector& subUDot, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // State must be realized to Stage::Instance and SubtreeResults must already be at
    // Stage::Velocity. SubtreeResults stage will be Stage::Acceleration after this call.
    // All Subtree udots's will be zero, body accelerations will have only their bias values
    // (coriolis accelerations from nonzero u's).
    void calcAccelerationsFromZeroUDot(const State&, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    // Calculates a perturbed velocity result starting with the subUDot's and acceleration results
    // which must already be in SubtreeResults.
    void perturbAccelerations(const State&, SubtreeUIndex subUDotIndex, Real perturbation, SimbodyMatterSubtreeResults::SubtreeResultsRep&) const;

    MobilizedBodyIndex getParentMobilizedBodyIndex(MobilizedBodyIndex childIx) const {
        return getSimbodyMatterSubsystem().getMobilizedBody(childIx)
                     .getParentMobilizedBody().getMobilizedBodyIndex();
    }

    int getLevel(MobilizedBodyIndex mbid) const {
        return getSimbodyMatterSubsystem().getMobilizedBody(mbid).getLevelInMultibodyTree();
    }

    bool hasPathToAncestor(MobilizedBodyIndex bid) const {
        assert(ancestor.isValid());
        while (bid!=ancestor && bid!=GroundIndex)
            bid = getParentMobilizedBodyIndex(bid);
        return bid == ancestor; // works if ancestor is Ground also
    }

    bool isTerminalBody(MobilizedBodyIndex bid) const {
        for (int i=0; i < (int)terminalBodies.size(); ++i)
            if (bid == terminalBodies[i])
                return true;
        return false;
    }


    void realizeTopology() {
        if (stage >= Stage::Topology)
            return;
        allBodies.clear(); parentSubtreeBodies.clear(); childSubtreeBodies.clear();

        if (terminalBodies.empty()) {
            stage = Stage::Topology;    // this is the "empty subtree"
            return;
        }

        ancestor = findAncestorBody();

        // We'll collect all the Subtree bodies in a MobilizedBodyIndex->SubtreeBodyIndex
        // map. We'll do this in two passes through the map -- the first to eliminate
        // duplicates and put the bodies in MobilizedBodyIndex order, the second to assign
        // SubtreeBodyIndexs which are just the resulting ordinals.
        // Complexity of this first pass is O(N log N) where N is the number
        // of unique bodies in the Subtree.
        typedef std::map<MobilizedBodyIndex, SubtreeBodyIndex> MapType;
        typedef MapType::value_type MapEntry;
        MapType subtreeBodyIndexMap;
        // Pre-load the map with the ancestor body and its subtree body id 0.
        subtreeBodyIndexMap.insert(MapEntry(ancestor, SubtreeBodyIndex(0)));
        for (int i=0; i < (int)terminalBodies.size(); ++i) {
            // Run down this branch adding any new bodies we encounter
            // until we hit one that is already in the map. If we get to
            // Ground without hitting the ancestor (OK if Ground *is* the
            // ancestor) then we have been given a bad terminal body which
            // should have been caught earlier.
            MobilizedBodyIndex mbid = terminalBodies[i];
            // ".second" will be true if the entry was actually inserted; otherwise
            // it was already there.
            while (subtreeBodyIndexMap.insert(MapEntry(mbid,SubtreeBodyIndex())).second) {
                assert(mbid != GroundIndex);
                mbid = getParentMobilizedBodyIndex(mbid);
            }
        }

        // Now assign the SubtreeBodyIndexs in order of MobilizedBodyIndex, and fill
        // in allBodies which serves as the SubtreeBodyIndex->MobilizedBodyIndex map,
        // and parentSubtreeBodies which maps a subtree body to its unique subtree
        // parent, and childSubtreeBodies which goes from a subtree body to the
        // list of all bodies for which it is the parent.
        // This pass is also O(N log N) because we have to look up the parent
        // mobilized body id in the map to get its assigned subtree body id.

        allBodies.resize((unsigned)subtreeBodyIndexMap.size());
        parentSubtreeBodies.resize((unsigned)subtreeBodyIndexMap.size());
        childSubtreeBodies.resize((unsigned)subtreeBodyIndexMap.size());
        allBodies[0] = ancestor;
        parentSubtreeBodies[0] = InvalidSubtreeBodyIndex;

        SubtreeBodyIndex nextFree(1); // ancestor was already assigned 0
        MapType::iterator body = subtreeBodyIndexMap.begin();
        assert(body->first == ancestor && body->second == 0);
        ++body; // skip the ancestor
        for (; body != subtreeBodyIndexMap.end(); ++body)
        {
            const MobilizedBodyIndex mbid = body->first;
            const SubtreeBodyIndex   sbid = nextFree++;
            body->second = sbid;
            allBodies[sbid] = mbid;

            // Look up the parent, which *must* be (a) present, and (b) 
            // one of the earlier Subtree bodies.
            const MobilizedBodyIndex pmbid = getParentMobilizedBodyIndex(mbid);
            MapType::const_iterator parent = subtreeBodyIndexMap.find(pmbid);
            assert(parent != subtreeBodyIndexMap.end());

            const SubtreeBodyIndex psbid = parent->second;
            assert(psbid < sbid && allBodies[psbid]==pmbid);

            parentSubtreeBodies[sbid] = psbid;
            childSubtreeBodies[psbid].push_back(sbid);
        }

        stage = Stage::Topology;
    }

    int getNumSubtreeBodies() const {
        assert(stage >= Stage::Topology);
        return (int)allBodies.size();
    }

    void realizeModel(const State& s) {
        assert(getSimbodyMatterSubsystem().getStage(s) >= Stage::Model);
        assert(stage >= Stage::Topology);
        if (stage >= Stage::Model)
            return;

        stage = Stage::Model;
    }

    void realizePosition(const State& s, const Vector& subQ) {
        assert(getSimbodyMatterSubsystem().getStage(s) >= Stage::Instance);
        assert(stage >= Stage::Model);

        stage = Stage::Position;
    }

    void realizeVelocity(const State& s, const Vector& subU) {
        assert(getSimbodyMatterSubsystem().getStage(s) >= Stage::Position);
        assert(stage >= Stage::Model);

        stage = Stage::Velocity;
    }

    void realizeAcceleration(const State& s, const Vector& subUDot) {
        assert(getSimbodyMatterSubsystem().getStage(s) >= Stage::Velocity);
        assert(stage >= Stage::Model);

        stage = Stage::Acceleration;
    }

    const SimbodyMatterSubtree& getMySubtreeOwnerHandle() const {return myHandle;}

private:
    friend class SimbodyMatterSubtree;
    const SimbodyMatterSubtree&   myHandle; // owner handle
    const SimbodyMatterSubsystem* matter;   // a reference to the full tree of which this is a subset

    Stage stage;    // initially invalid

        // TOPOLOGY STATE VARIABLES

    Array_<MobilizedBodyIndex> terminalBodies;

        // TOPOLOGY CACHE VARIABLES

    // Maps SubtreeBodyIndex to MobilizedBodyIndex
    MobilizedBodyIndex        ancestor;
    Array_<MobilizedBodyIndex> allBodies; // ancestor body is 0; ids are in increasing order

    // Maps each subtree body (by SubtreeBodyIndex) to its unique parent within the subtree
    // the base body (SubtreeBodyIndex==SubtreeBaseBodyIndex==0) returns an InvalidSubtreeBodyIndex
    // as its parent.
    Array_<SubtreeBodyIndex>   parentSubtreeBodies;

    // Maps each subtree body to its children within the subtree. Note that a subtree terminal
    // body may have children in the full matter subsystem, but which are not included in
    // the subtree.
    Array_< Array_<SubtreeBodyIndex> > childSubtreeBodies;

private:

    // This routine finds the terminal body closest to Ground in the
    // MobilizedBody tree's graph, then moves down all the other branches
    // to find a body in each branch at that same lowest level. That is,
    // we "trim" all the branches to be the same height. Then we move
    // all the branches in sync one level closer to ground until they
    // all hit the same body. That's the outmost common ancestor.
    MobilizedBodyIndex findAncestorBody() {
        assert(terminalBodies.size());

        // Copy the terminal bodies, which are the current branch tips.
        Array_<MobilizedBodyIndex> tips(&terminalBodies[0], (&terminalBodies[0])+terminalBodies.size());

        // Find the level of the lowest-level tip.
        int minTip = 0;
        int minLevel = getLevel(tips[minTip]);
        for (int i=1; i < (int)tips.size(); ++i)
            if (getLevel(tips[i]) < minLevel)
               {minTip = i; minLevel = getLevel(tips[minTip]);}

        // Trim all the other branches back to the lowest level.
        for (int i=0; i < (int)tips.size(); ++i)
            while (getLevel(tips[i]) > minLevel)
                tips[i] = getParentMobilizedBodyIndex(tips[i]);

        // All tips are at the same level now. March them in lockstep down
        // to their common ancestor or Ground.
        while (!allElementsMatch(tips))
            pruneTipsOneLevel(tips);

        return tips[0]; // all branches led here
    }

    static bool allElementsMatch(const Array_<MobilizedBodyIndex>& ids) {
        for (int i=1; i < (int)ids.size(); ++i)
            if (ids[i] != ids[0]) return false;
        return true;
    }

    void pruneTipsOneLevel(Array_<MobilizedBodyIndex>& tips) const {
        for (int i=0; i < (int)tips.size(); ++i) {
            assert(tips[i] != GroundIndex); // can't happen: they should have matched!
            tips[i] = getParentMobilizedBodyIndex(tips[i]);
        }
    }
};

    /////////////////////////
    // SUBTREE RESULTS REP //
    /////////////////////////

class SimbodyMatterSubtreeResults::SubtreeResultsRep {
public:
    explicit SubtreeResultsRep(const SimbodyMatterSubtreeResults& handle) 
      : myHandle(&handle), stage(Stage::Empty)
    { 
        clear();
    }

    void clear() {
        qSubset.clear(); uSubset.clear();
        qSeg.clear(); uSeg.clear();
        subQ.clear(); subU.clear(); subUDot.clear();
        bodyTransforms.clear();    perturbedBodyTransforms.clear();
        bodyVelocities.clear();    perturbedBodyVelocities.clear();
        bodyAccelerations.clear(); perturbedBodyAccelerations.clear();
        perturbedQ = InvalidSubtreeQIndex;
        perturbedU = perturbedUDot = InvalidSubtreeUIndex;
        stage = Stage::Empty;
    }

    // Set or change the number of Subtree bodies to be accommodated here. A
    // "change" is cheap if the number of bodies hasn't actually changed.
    void reallocateBodies(int nb) {
        assert(nb >= 1); // must be at least the ancestor body
        qSeg.resize(nb); uSeg.resize(nb);
        bodyTransforms.resize(nb);    perturbedBodyTransforms.resize(nb);
        bodyVelocities.resize(nb);    perturbedBodyVelocities.resize(nb);
        bodyAccelerations.resize(nb); perturbedBodyAccelerations.resize(nb);

        // Set the unchanging results for the ancestor body, which is treated as Ground.
        qSeg[0] = pair<SubtreeQIndex,int>(SubtreeQIndex(0), 0); // no q's
        uSeg[0] = pair<SubtreeUIndex,int>(SubtreeUIndex(0), 0); // no u's

        bodyTransforms[0]    = perturbedBodyTransforms[0]    = Transform();
        bodyVelocities[0]    = perturbedBodyVelocities[0]    = SpatialVec(Vec3(0), Vec3(0));
        bodyAccelerations[0] = perturbedBodyAccelerations[0] = SpatialVec(Vec3(0), Vec3(0));

        // Clear q and u mapping information -- that can't be known until Stage::Model.
        qSubset.clear();
        uSubset.clear();

        perturbedQ = InvalidSubtreeQIndex;
        perturbedU = perturbedUDot = InvalidSubtreeUIndex;
        stage = Stage::Topology;
    }

    // Assign the next available Subtree q and u slots to the indicated body, and
    // remember them. We are given the MatterSubsystem q's and u's associated with
    // the corresponding mobilized body so we can keep a mapping.
    void addMobilities(SubtreeBodyIndex sb, QIndex qStart, int nq, UIndex uStart, int nu) {
        assert(stage >= Stage::Topology);
        assert(nq >= 0 && nu >= 0 && nq >= nu);
        assert(1 <= sb && sb < getNumSubtreeBodies());
        stage = Stage::Topology; // back up if necessary

        qSeg[sb] = pair<SubtreeQIndex,int>(SubtreeQIndex(qSubset.size()), nq);
        uSeg[sb] = pair<SubtreeUIndex,int>(SubtreeUIndex(uSubset.size()), nu);

        for (int i=0; i<nq; ++i)
            qSubset.push_back(QIndex(qStart+i));
        for (int i=0; i<nu; ++i)
            uSubset.push_back(UIndex(uStart+i));
    }

    void packStateQIntoSubtreeQ(const Vector& stateQ, Vector& subtreeQ) const {
        assert(stage >= Stage::Model);
        assert(stateQ.size() >= getNumSubtreeQs());

        subtreeQ.resize(getNumSubtreeQs());
        for (SubtreeQIndex i(0); i<getNumSubtreeQs(); ++i)
            subtreeQ[i] = stateQ[qSubset[i]];
    }

    void packStateUIntoSubtreeU(const Vector& stateU, Vector& subtreeU) const {
        assert(stage >= Stage::Model);
        assert(stateU.size() >= getNumSubtreeUs());

        subtreeU.resize(getNumSubtreeUs());
        for (SubtreeUIndex i(0); i<getNumSubtreeUs(); ++i)
            subtreeU[i] = stateU[uSubset[i]];
    }


    void unpackSubtreeQIntoStateQ(const Vector& subtreeQ, Vector& stateQ) const {
        assert(stage >= Stage::Model);
        assert(subtreeQ.size() == getNumSubtreeQs());

        for (SubtreeQIndex i(0); i<getNumSubtreeQs(); ++i)
            stateQ[qSubset[i]] = subtreeQ[i];
    }

    // Call this when done adding mobilities.
    void realizeModel(const Vector& allQ, const Vector& allU) {
        stage = Stage::Model; // enable routines used below
        assert(allQ.size() >= getNumSubtreeQs() && allU.size() >= getNumSubtreeUs());

        subQ.resize(getNumSubtreeQs()); 
        subU.resize(getNumSubtreeUs()); 
        subUDot.resize(getNumSubtreeUs());

        packStateQIntoSubtreeQ(allQ,subQ);  // set initial values
        packStateUIntoSubtreeU(allU,subU);

        subUDot = NaN;
    }

    QIndex mapSubtreeQToSubsystemQ(SubtreeQIndex sq) const {
        assert(stage >= Stage::Model);
        assert(0 <= sq && sq < (int)qSubset.size());
        return qSubset[sq]; // range checked indexing
    }
    UIndex mapSubtreeUToSubsystemU(SubtreeUIndex su) const {
        assert(stage >= Stage::Model);
        assert(0 <= su && su < (int)uSubset.size());
        return uSubset[su];
    }

    int getNumSubtreeBodies() const {
        assert(stage >= Stage::Topology);
        return (int)bodyTransforms.size();
    }

    int getNumSubtreeQs() const {
        assert(stage >= Stage::Model);
        return (int)qSubset.size();
    }

    int getNumSubtreeUs() const {
        assert(stage >= Stage::Model);
        return (int)uSubset.size();
    }

    void setMySubtreeResultsOwnerHandle(const SimbodyMatterSubtreeResults& owner) {
        myHandle = &owner;
    }

    const SimbodyMatterSubtreeResults& getMySubtreeResultsOwnerHandle() const {
        assert(myHandle);
        return *myHandle;
    }

    Stage getStage() const {return stage;}
    void setStage(Stage g) {stage=g;}

    const Vector& getSubQ() const {
        assert(stage >= Stage::Position);
        return subQ;
    }
    Vector& updSubQ() {
        assert(stage >= Stage::Model);
        if (stage >= Stage::Position) stage=Stage::Model;
        return subQ;
    }
    const Vector& getSubU() const {
        assert(stage >= Stage::Velocity);
        return subU;
    }
    Vector& updSubU() {
        assert(stage >= Stage::Model);
        if (stage >= Stage::Velocity) stage=Stage::Position;
        return subU;
    }
    const Vector& getSubUDot() const {
        assert(stage >= Stage::Acceleration);
        return subUDot;
    }
    Vector& updSubUDot() {
        assert(stage >= Stage::Model);
        if (stage >= Stage::Acceleration) stage=Stage::Velocity;
        return subUDot;
    }

    const Transform& getSubtreeBodyTransform(SubtreeBodyIndex sb) { // X_AB
        assert(stage >= Stage::Position);
        assert(0 <= sb && sb < (int)bodyTransforms.size());
        return bodyTransforms[sb];
    }

    const SpatialVec& getSubtreeBodyVelocity(SubtreeBodyIndex sb) {
        assert(stage >= Stage::Velocity);
        assert(0 <= sb && sb < (int)bodyVelocities.size());
        return bodyVelocities[sb];
    }

    const SpatialVec& getSubtreeBodyAcceleration(SubtreeBodyIndex sb) {
        assert(stage >= Stage::Acceleration);
        assert(0 <= sb && sb < (int)bodyAccelerations.size());
        return bodyAccelerations[sb];
    }

    void setSubtreeBodyTransform(SubtreeBodyIndex sb, const Transform& X_AB) {
        assert(stage >= Stage::Model);
        assert(1 <= sb && sb < getNumSubtreeBodies()); // can't set Ancestor transform
        bodyTransforms[sb] = X_AB;
    }

    void setSubtreeBodyVelocity(SubtreeBodyIndex sb, const SpatialVec& V_AB) {
        assert(stage >= Stage::Position);
        assert(1 <= sb && sb < getNumSubtreeBodies()); // can't set Ancestor velocity
        bodyVelocities[sb] = V_AB;
    }

    void setSubtreeBodyAcceleration(SubtreeBodyIndex sb, const SpatialVec& A_AB) {
        assert(stage >= Stage::Velocity);
        assert(1 <= sb && sb < getNumSubtreeBodies()); // can't set Ancestor velocity
        bodyAccelerations[sb] = A_AB;
    }

    void findSubtreeBodyQ(SubtreeBodyIndex sb, SubtreeQIndex& q, int& nq) const {
        assert(stage >= Stage::Model);
        q  = qSeg[sb].first;
        nq = qSeg[sb].second;
    }

    void findSubtreeBodyU(SubtreeBodyIndex sb, SubtreeUIndex& u, int& nu) const {
        assert(stage >= Stage::Model);
        u  = uSeg[sb].first;
        nu = uSeg[sb].second;
    }
private:
    friend class SimbodyMatterSubtreeResults;
    const SimbodyMatterSubtreeResults* myHandle; // owner handle

    Stage stage;    // initially invalid

    // Model stage information
    Array_< QIndex > qSubset; // map from SubtreeQIndex to MatterSubsystem q
    Array_< UIndex > uSubset; // map from SubtreeUIndex to MatterSubsystem u (also udot)

    // These identify which mobilities go with which bodies.
    Array_< pair<SubtreeQIndex,int> > qSeg;  // map from SubtreeBodyIndex to qSubset offset, length
    Array_< pair<SubtreeUIndex,int> > uSeg;  // map from SubtreeBodyIndex to uSubset offset, length

    //TODO: make PIMPL
    Vector                 subQ;                        // generalized coords for Subtree bodies
    Array_<Transform> bodyTransforms;              // X_AB, index by SubtreeBodyIndex (unperturbed)

    SubtreeQIndex          perturbedQ;                  // which Q was perturbed? InvalidSubtreeQIndex if none
    Array_<Transform> perturbedBodyTransforms;     // X_AB, after perturbation

    Vector                 subU;                        // generalized speeds for Subtree bodies
    Vector_<SpatialVec>    bodyVelocities;              // V_AB, index by SubtreeBodyIndex

    SubtreeUIndex          perturbedU;                  // which u was perturbed? InvalidSubtreeUIndex if none
    Vector_<SpatialVec>    perturbedBodyVelocities;     // V_AB, after perturbation

    Vector                 subUDot;                     // generalized accelerations for Subtree bodies
    Vector_<SpatialVec>    bodyAccelerations;           // A_AB, index by SubtreeBodyIndex

    SubtreeUIndex          perturbedUDot;               // which udot was perturbed? InvalidSubtreeUIndex if none
    Vector_<SpatialVec>    perturbedBodyAccelerations;  // A_AB, after perturbation
};


    ////////////////////////////
    // SIMBODY MATTER SUBTREE //
    ////////////////////////////

// Default constructor -- we don't know the SimbodyMatterSubsystem yet.
SimbodyMatterSubtree::SimbodyMatterSubtree()
  : rep(0)
{
    rep = new SubtreeRep(*this);
}

SimbodyMatterSubtree::SimbodyMatterSubtree(const SimbodyMatterSubsystem& matter)
  : rep(0)
{
    rep = new SubtreeRep(*this, matter);
}


SimbodyMatterSubtree::SimbodyMatterSubtree(const SimbodyMatterSubsystem& matter,
                                           const Array_<MobilizedBodyIndex>& terminalBodies)
  : rep(0)
{
    rep = new SubtreeRep(*this, matter);
    rep->setTerminalBodies(terminalBodies);
}

// Copy constructor
SimbodyMatterSubtree::SimbodyMatterSubtree(const SimbodyMatterSubtree& src) 
  : rep(0)
{
    if (src.rep)
        rep = new SubtreeRep(*src.rep);
}

// Copy assignment
SimbodyMatterSubtree& 
SimbodyMatterSubtree::operator=(const SimbodyMatterSubtree& src)
{
    if (&src != this) {
        if (rep && (this == &rep->getMySubtreeOwnerHandle()))
            delete rep;
        rep = 0;
        if (src.rep)
            rep = new SubtreeRep(*src.rep);
    }
    return *this;
}

// Destructor
SimbodyMatterSubtree::~SimbodyMatterSubtree() {
    if (rep && (this == &rep->getMySubtreeOwnerHandle()))
        delete rep; 
    rep=0;
}

const SimbodyMatterSubsystem&
SimbodyMatterSubtree::getSimbodyMatterSubsystem() const {
    return getRep().getSimbodyMatterSubsystem();
}


void SimbodyMatterSubtree::setSimbodyMatterSubsystem(const SimbodyMatterSubsystem& matter) {
    return updRep().setSimbodyMatterSubsystem(matter);
}

void SimbodyMatterSubtree::clear() {
    return updRep().clear();
}


SimbodyMatterSubtree& 
SimbodyMatterSubtree::addTerminalBody(MobilizedBodyIndex i) {
    updRep().addTerminalBody(i);
    return *this;
}

void SimbodyMatterSubtree::realizeTopology() {
    updRep().realizeTopology();
}

MobilizedBodyIndex SimbodyMatterSubtree::getAncestorMobilizedBodyIndex() const {
    return getRep().ancestor;
}

const Array_<MobilizedBodyIndex>& 
SimbodyMatterSubtree::getTerminalBodies() const {
    return getRep().terminalBodies;
}

int SimbodyMatterSubtree::getNumSubtreeBodies() const {
    return (int)getRep().allBodies.size();
}

const Array_<MobilizedBodyIndex>& 
SimbodyMatterSubtree::getAllBodies() const {
    assert(getRep().stage >= Stage::Topology);
    return getRep().allBodies;
}

SubtreeBodyIndex 
SimbodyMatterSubtree::getParentSubtreeBodyIndex(SubtreeBodyIndex sbid) const {
    assert(getRep().stage >= Stage::Topology);
    return getRep().parentSubtreeBodies[sbid];
}
const Array_<SubtreeBodyIndex>& 
SimbodyMatterSubtree::getChildSubtreeBodyIndices(SubtreeBodyIndex sbid) const {
    assert(getRep().stage >= Stage::Topology);
    return getRep().childSubtreeBodies[sbid];
}

bool SimbodyMatterSubtree::
isCompatibleSubtreeResults(const SimbodyMatterSubtreeResults& sr) const {
    return getRep().isCompatibleSubtreeResults(sr.getRep());
}

void SimbodyMatterSubtree::initializeSubtreeResults(const State& s, SimbodyMatterSubtreeResults& sr) const {
    getRep().initializeSubtreeResults(s,sr.updRep());
}


void SimbodyMatterSubtree::
copyPositionsFromState(const State& s, SimbodyMatterSubtreeResults& sr) const {
    getRep().copyPositionsFromState(s,sr.updRep());
}

void SimbodyMatterSubtree::
calcPositionsFromSubtreeQ(const State& s, const Vector& subQ, SimbodyMatterSubtreeResults& sr) const {
    getRep().calcPositionsFromSubtreeQ(s,subQ,sr.updRep());
}

void SimbodyMatterSubtree::
perturbPositions(const State& s, SubtreeQIndex subQIndex, Real perturbation, SimbodyMatterSubtreeResults& sr) const {
    getRep().perturbPositions(s,subQIndex,perturbation,sr.updRep());
}

void SimbodyMatterSubtree::
copyVelocitiesFromState(const State& s, SimbodyMatterSubtreeResults& sr) const {
    getRep().copyVelocitiesFromState(s,sr.updRep());
}

void SimbodyMatterSubtree::
calcVelocitiesFromSubtreeU(const State& s, const Vector& subU, SimbodyMatterSubtreeResults& sr) const {
    getRep().calcVelocitiesFromSubtreeU(s,subU,sr.updRep());
}

void SimbodyMatterSubtree::
calcVelocitiesFromZeroU(const State& s, SimbodyMatterSubtreeResults& sr) const {
    getRep().calcVelocitiesFromZeroU(s,sr.updRep());
}

void SimbodyMatterSubtree::
perturbVelocities(const State& s, SubtreeUIndex subUIndex, Real perturbation, SimbodyMatterSubtreeResults& sr) const {
    getRep().perturbVelocities(s,subUIndex,perturbation,sr.updRep());
}

void SimbodyMatterSubtree::
copyAccelerationsFromState(const State& s, SimbodyMatterSubtreeResults& sr) const {
    getRep().copyAccelerationsFromState(s,sr.updRep());
}

void SimbodyMatterSubtree::
calcAccelerationsFromSubtreeUDot(const State& s, const Vector& subUDot, SimbodyMatterSubtreeResults& sr) const {
    getRep().calcAccelerationsFromSubtreeUDot(s,subUDot,sr.updRep());
}

void SimbodyMatterSubtree::
calcAccelerationsFromZeroUDot(const State& s, SimbodyMatterSubtreeResults& sr) const {
    getRep().calcAccelerationsFromZeroUDot(s,sr.updRep());
}

void SimbodyMatterSubtree::
perturbAccelerations(const State& s, SubtreeUIndex subUDotIndex, Real perturbation, SimbodyMatterSubtreeResults& sr) const {
    getRep().perturbAccelerations(s,subUDotIndex,perturbation,sr.updRep());
}



    ////////////////////////////////////
    // SIMBODY MATTER SUBTREE RESULTS //
    ////////////////////////////////////

SimbodyMatterSubtreeResults::SimbodyMatterSubtreeResults() : rep(0) {
    rep = new SubtreeResultsRep(*this);
}

SimbodyMatterSubtreeResults::~SimbodyMatterSubtreeResults() {
    if (rep && this == rep->myHandle)
        delete rep;
    rep = 0;
}

SimbodyMatterSubtreeResults::SimbodyMatterSubtreeResults(const SimbodyMatterSubtreeResults& src) : rep(0) {
    if (src.rep) {
        rep = new SubtreeResultsRep(*src.rep);
        rep->setMySubtreeResultsOwnerHandle(*this);
    }
}

SimbodyMatterSubtreeResults& 
SimbodyMatterSubtreeResults::operator=(const SimbodyMatterSubtreeResults& src) {
    if (&src != this) {
        if (rep && this == rep->myHandle)
            delete rep;
        rep = 0;
        if (src.rep) {
            rep = new SubtreeResultsRep(*src.rep);
            rep->setMySubtreeResultsOwnerHandle(*this);
        }
    }
    return *this;
}

int SimbodyMatterSubtreeResults::getNumSubtreeBodies() const {
    return getRep().getNumSubtreeBodies();
}

int SimbodyMatterSubtreeResults::getNumSubtreeQs() const {
    return getRep().getNumSubtreeQs();
}

int SimbodyMatterSubtreeResults::getNumSubtreeUs() const {
    return getRep().getNumSubtreeUs();
}

void SimbodyMatterSubtreeResults::reallocateBodies(int nb) {
    updRep().reallocateBodies(nb);
}

void SimbodyMatterSubtreeResults::addMobilities
   (SubtreeBodyIndex sb, QIndex qStart, int nq, UIndex uStart, int nu)
{
    updRep().addMobilities(sb, qStart, nq, uStart, nu);
}

void SimbodyMatterSubtreeResults::realizeModel(const Vector& stateQ, const Vector& stateU) {
    updRep().realizeModel(stateQ, stateU);
}

Stage SimbodyMatterSubtreeResults::getStage() const {
    return getRep().stage;
}

const Array_<QIndex>& SimbodyMatterSubtreeResults::getQSubset() const {
    assert(getRep().stage >= Stage::Model);
    return getRep().qSubset;
}

const Array_<UIndex>& SimbodyMatterSubtreeResults::getUSubset() const {
    assert(getRep().stage >= Stage::Model);
    return getRep().uSubset;
}

void SimbodyMatterSubtreeResults::findSubtreeBodyQ(SubtreeBodyIndex sbid, SubtreeQIndex& qStart, int& nq) const {
    assert(getStage() >= Stage::Model);
    const pair<SubtreeQIndex,int>& seg = getRep().qSeg[sbid];
    qStart = seg.first;
    nq     = seg.second;
}

void SimbodyMatterSubtreeResults::findSubtreeBodyU(SubtreeBodyIndex sbid, SubtreeUIndex& uStart, int& nu) const {
    assert(getStage() >= Stage::Model);
    const pair<SubtreeUIndex,int>& seg = getRep().uSeg[sbid];
    uStart = seg.first;
    nu     = seg.second;
}

const Vector& SimbodyMatterSubtreeResults::getSubtreeQ() const {
    assert(getStage() >= Stage::Position);
    return getRep().subQ;
}
const Transform& SimbodyMatterSubtreeResults::getSubtreeBodyTransform(SubtreeBodyIndex sbid) const {
    assert(getStage() >= Stage::Position);
    return getRep().bodyTransforms[sbid];
}

const Vector& SimbodyMatterSubtreeResults::getSubtreeU() const {
    assert(getStage() >= Stage::Velocity);
    return getRep().subU;
}
const SpatialVec& SimbodyMatterSubtreeResults::getSubtreeBodyVelocity(SubtreeBodyIndex sbid) const {
    assert(getStage() >= Stage::Velocity);
    return getRep().bodyVelocities[sbid];
}

const Vector& SimbodyMatterSubtreeResults::getSubtreeUDot() const {
    assert(getStage() >= Stage::Acceleration);
    return getRep().subUDot;
}
const SpatialVec& SimbodyMatterSubtreeResults::getSubtreeBodyAcceleration(SubtreeBodyIndex sbid) const {
    assert(getStage() >= Stage::Acceleration);
    return getRep().bodyAccelerations[sbid];
}

std::ostream& operator<<(std::ostream& o, const SimbodyMatterSubtree& sub) {
    o << "SUBTREE:" << endl;
    o << "  ancestor=" << sub.getAncestorMobilizedBodyIndex();

    o << "  terminalBodies=";
    for (int i=0; i < (int)sub.getTerminalBodies().size(); ++i)
        o << sub.getTerminalBodies()[i] << " ";
    o << endl;

    o << "  allBodies=";
    for (int i=0; i < (int)sub.getAllBodies().size(); ++i)
        o << sub.getAllBodies()[i] << " ";
    o << endl;

    for (SubtreeBodyIndex b(0); b < (int)sub.getAllBodies().size(); ++b) {
        o << "  parent[" << b << "]=" << sub.getParentSubtreeBodyIndex(b);

        o << "  children[" << b << "]=";
        for (int i=0; i < (int)sub.getChildSubtreeBodyIndices(b).size(); ++i)
            o << sub.getChildSubtreeBodyIndices(b)[i] << " ";
        o << endl;
    }

    return o;
}

static std::ostream& operator<<(std::ostream& o, const Array_<QIndex>& q) {
    for (int i=0; i<(int)q.size(); ++i)
        o << q[i] << " ";
    return o;
}

static std::ostream& operator<<(std::ostream& o, const Array_<UIndex>& u) {
    for (int i=0; i<(int)u.size(); ++i)
        o << u[i] << " ";
    return o;
}

std::ostream& operator<<(std::ostream& o, const SimbodyMatterSubtreeResults& sr) {
    o << "SUBTREE RESULTS (stage=" << sr.getStage() << "):" << endl;

    if (sr.getStage() >= Stage::Topology)
        o << "  " << sr.getNumSubtreeBodies() << " subtree bodies" << endl;

    if (sr.getStage() >= Stage::Model) {
        o << "  nq=" << sr.getNumSubtreeQs() << ", nu=" << sr.getNumSubtreeUs() << endl;
        o << "  QSubset: " << sr.getQSubset() << endl;
        o << "  USubset: " << sr.getUSubset() << endl;

        for (SubtreeBodyIndex sb(1); sb < sr.getNumSubtreeBodies(); ++sb) {
            SubtreeQIndex qstart; int nq;
            SubtreeUIndex ustart; int nu;
            sr.findSubtreeBodyQ(sb,qstart,nq);
            sr.findSubtreeBodyU(sb,ustart,nu);
            o << "  body " << sb << " q=" << qstart << ".." << qstart+nq-1
                << " u=" << ustart << ".." << ustart+nu-1 << endl;
        }
    }

    if (sr.getStage() >= Stage::Position) {
        o << "  POSITION RESULTS AVAILABLE:\n";
        o << "    q=" << sr.getSubtreeQ() << endl;
        for (SubtreeBodyIndex sb(0); sb < sr.getNumSubtreeBodies(); ++sb)
            o << "    X_AB" << sb << "=" << sr.getSubtreeBodyTransform(sb);
    }

    if (sr.getStage() >= Stage::Velocity) {
        o << "  VELOCITY RESULTS AVAILABLE\n";
        o << "    u=" << sr.getSubtreeU() << endl;
        for (SubtreeBodyIndex sb(0); sb < sr.getNumSubtreeBodies(); ++sb)
            o << "    V_AB" << sb << "=" << sr.getSubtreeBodyVelocity(sb) << endl;
    }

    if (sr.getStage() >= Stage::Acceleration) {
        o << "  ACCELERATION RESULTS AVAILABLE\n";
        o << "    udot=" << sr.getSubtreeUDot() << endl;
        for (SubtreeBodyIndex sb(0); sb < sr.getNumSubtreeBodies(); ++sb)
            o << "    A_AB" << sb << "=" << sr.getSubtreeBodyAcceleration(sb) << endl;
    }

    o << "END SUBTREE RESULTS." << endl;

    return o;
}


////////////////////////////////////////////////////////
// SIMBODY MATTER SUBSYSTEM :: SUBTREE :: SUBTREE REP //
////////////////////////////////////////////////////////

void SimbodyMatterSubtree::SubtreeRep::
initializeSubtreeResults(const State& s, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Model, 
        "Subtree::initializeSubtreeResults()");

    const int nSubtreeBodies = getNumSubtreeBodies();
    sr.reallocateBodies(nSubtreeBodies);

    int nSubtreeQ=0, nSubtreeU=0;
    // start at 1 because ancestor has no relevant mobilities
    for (SubtreeBodyIndex sb(1); sb < nSubtreeBodies; ++sb) {
        const MobilizedBodyIndex mb = allBodies[sb];

        QIndex qStart; int nq;
        UIndex uStart; int nu;
        matter.findMobilizerQs(s, mb, qStart, nq);
        matter.findMobilizerUs(s, mb, uStart, nu);
        nSubtreeQ += nq; nSubtreeU += nu;

        sr.addMobilities(sb, qStart, nq, uStart, nu);
    }

    sr.realizeModel(matter.getQ(s), matter.getU(s));
    assert(nSubtreeQ == sr.getNumSubtreeQs() && nSubtreeU == sr.getNumSubtreeUs());
}

bool SimbodyMatterSubtree::SubtreeRep::
isCompatibleSubtreeResults(const SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    // hardly exhaustive but better than nothing
    return sr.getStage() >= Stage::Model 
        && sr.getNumSubtreeBodies() == getNumSubtreeBodies();
}

void  SimbodyMatterSubtree::SubtreeRep::
copyPositionsFromState(const State& s, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Position, 
        "Subtree::calcPositionsFromState()");

    assert(isCompatibleSubtreeResults(sr));

    // Copy the q's; adjust the body transforms to be relative to the ancestor
    // body instead of ground.
    sr.packStateQIntoSubtreeQ(matter.getQ(s), sr.updSubQ());

    if (getAncestorMobilizedBodyIndex() == GroundIndex) {
        for (SubtreeBodyIndex sb(1); sb < getNumSubtreeBodies(); ++sb) {
            const MobilizedBodyIndex mb = getSubtreeBodyMobilizedBodyIndex(sb);
            const Transform& X_GB = matter.getBodyTransform(s,mb);
            sr.setSubtreeBodyTransform(sb, X_GB); // =X_AB
        }
    } else {
        // Ancestor A differs from Ground G so we have to adjust all the 
        // Subtree body transforms to measure from A instead of G.
        const Transform& X_GA = matter.getBodyTransform(s,getAncestorMobilizedBodyIndex());
        for (SubtreeBodyIndex sb(1); sb < getNumSubtreeBodies(); ++sb) {
            const MobilizedBodyIndex mb = getSubtreeBodyMobilizedBodyIndex(sb);
            const Transform& X_GB = matter.getBodyTransform(s,mb);
            sr.setSubtreeBodyTransform(sb, ~X_GA * X_GB); // X_AB
        }
    }

    sr.setStage(Stage::Position);
}

// Here we are given a new set of Subtree q's, and we want to calculate the resulting A-relative
// transforms for all the Subtree bodies. This requires calculating the cross-mobilizer transforms
// X_FM for each Subtree mobilizer and propagating them outwards towards the terminal bodies.
void SimbodyMatterSubtree::SubtreeRep::
calcPositionsFromSubtreeQ(const State& state, const Vector& subQ, 
                          SimbodyMatterSubtreeResults::SubtreeResultsRep& results) const
{
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(state), Stage::Instance, 
        "calcPositionsFromSubtreeQ()");

    assert(isCompatibleSubtreeResults(results));
    assert(subQ.size() == results.getNumSubtreeQs());

    results.updSubQ() = subQ; // copy in the q's

    // For high speed, find memory address for the first subQ; they are sequential after this.
    const Real* allSubQ = &results.getSubQ()[0];

    // Iterate from the ancestor outward to propagate the transforms to the terminal bodies.
    for (SubtreeBodyIndex sb(1); sb < getNumSubtreeBodies(); ++sb) {
        const SubtreeBodyIndex   sp = getParentSubtreeBodyIndex(sb);
        const MobilizedBodyIndex mb = getSubtreeBodyMobilizedBodyIndex(sb);

        SubtreeQIndex firstSubQ; int nq;
        results.findSubtreeBodyQ(sb, firstSubQ, nq);

        const Transform  X_PB = matter.calcMobilizerTransformFromQ(state, mb, nq, &allSubQ[firstSubQ]); 
        const Transform& X_AP = results.getSubtreeBodyTransform(sp);
        results.setSubtreeBodyTransform(sb, X_AP*X_PB);
    }

    results.setStage(Stage::Position);
}

void SimbodyMatterSubtree::SubtreeRep::
perturbPositions(const State& s, SubtreeQIndex subQIndex, Real perturbation, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Instance, 
        "perturbPositions()");

    assert(isCompatibleSubtreeResults(sr));
    assert(sr.getStage() >= Stage::Position);

    assert(!"not implemented yet");
}

void SimbodyMatterSubtree::SubtreeRep::
copyVelocitiesFromState(const State& s, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Velocity, 
        "calcVelocitiesFromState()");

    assert(isCompatibleSubtreeResults(sr));

    // Copy the u's; adjust the body velocities to be measured relative to,
    // and expressed in, the ancestor body instead of ground.
    sr.packStateUIntoSubtreeU(matter.getU(s), sr.updSubU());

    if (getAncestorMobilizedBodyIndex() == GroundIndex) {
        for (SubtreeBodyIndex sb(1); sb < getNumSubtreeBodies(); ++sb) {
            const MobilizedBodyIndex mb = getSubtreeBodyMobilizedBodyIndex(sb);
            const SpatialVec& V_GB = matter.getBodyVelocity(s,mb);
            sr.setSubtreeBodyVelocity(sb, V_GB); // =V_AB
        }
    } else {
        // Ancestor A differs from Ground G so we have to adjust all the 
        // Subtree body velocities to measure from A instead of G.
        const Transform&  X_GA = matter.getBodyTransform(s,getAncestorMobilizedBodyIndex());
        const SpatialVec& V_GA = matter.getBodyVelocity(s,getAncestorMobilizedBodyIndex());
        for (SubtreeBodyIndex sb(1); sb < getNumSubtreeBodies(); ++sb) {
            const MobilizedBodyIndex mb = getSubtreeBodyMobilizedBodyIndex(sb);
            const Transform&  X_GB = matter.getBodyTransform(s,mb);
            const SpatialVec& V_GB = matter.getBodyVelocity(s,mb);

            const Vec3 p_AB_G     = X_GB.p() - X_GA.p();
            const Vec3 p_AB_G_dot = V_GB[1]  - V_GA[1];        // time deriv of p taken in G

            const Vec3 w_AB_G = V_GB[0] - V_GA[0];             // relative angular velocity
            const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G; // time deriv of p in A, exp in G
            const SpatialVec V_AB = ~X_GA.R() * SpatialVec(w_AB_G, v_AB_G); // re-express in A
            sr.setSubtreeBodyVelocity(sb, V_AB);
        }
    }

    sr.setStage(Stage::Velocity);
}

void SimbodyMatterSubtree::SubtreeRep::
calcVelocitiesFromSubtreeU(const State& s, const Vector& subU, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Instance, 
        "calcVelocitiesFromSubtreeU()");

    assert(isCompatibleSubtreeResults(sr));

    assert(!"not implemented yet");
}

void SimbodyMatterSubtree::SubtreeRep::
calcVelocitiesFromZeroU(const State& s, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Instance, 
        "calcVelocitiesFromZeroU()");

    assert(isCompatibleSubtreeResults(sr));

    assert(!"not implemented yet");
}

void SimbodyMatterSubtree::SubtreeRep::
perturbVelocities(const State& s, SubtreeUIndex subUIndex, Real perturbation, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Instance, 
        "perturbVelocities()");

    assert(isCompatibleSubtreeResults(sr));
    assert(sr.getStage() >= Stage::Velocity);

    assert(!"not implemented yet");
}

void SimbodyMatterSubtree::SubtreeRep::
copyAccelerationsFromState(const State& s, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Acceleration, 
        "calcAccelerationsFromState()");

    assert(isCompatibleSubtreeResults(sr));

    // Copy the udot's; adjust the body accelerations to be measured relative to,
    // and expressed in, the ancestor body instead of ground.
    sr.packStateUIntoSubtreeU(matter.getUDot(s), sr.updSubUDot());

    if (getAncestorMobilizedBodyIndex() == GroundIndex) {
        for (SubtreeBodyIndex sb(1); sb < getNumSubtreeBodies(); ++sb) {
            const MobilizedBodyIndex mb = getSubtreeBodyMobilizedBodyIndex(sb);
            const SpatialVec& A_GB = matter.getBodyAcceleration(s,mb);
            sr.setSubtreeBodyAcceleration(sb, A_GB); // =A_AB
        }
    } else {
        // Ancestor A differs from Ground G so we have to adjust all the 
        // Subtree body accelerations to measure from A instead of G.
        const Transform&  X_GA = matter.getBodyTransform(s,getAncestorMobilizedBodyIndex());
        const SpatialVec& V_GA = matter.getBodyVelocity(s,getAncestorMobilizedBodyIndex());
        const SpatialVec& A_GA = matter.getBodyAcceleration(s,getAncestorMobilizedBodyIndex());
        for (SubtreeBodyIndex sb(1); sb < getNumSubtreeBodies(); ++sb) {
            const MobilizedBodyIndex mb = getSubtreeBodyMobilizedBodyIndex(sb);
            const Transform&  X_GB = matter.getBodyTransform(s,mb);
            const SpatialVec& V_GB = matter.getBodyVelocity(s,mb);
            const SpatialVec& A_GB = matter.getBodyAcceleration(s,mb);

            const Vec3 p_AB_G        = X_GB.p() - X_GA.p();
            const Vec3 p_AB_G_dot    = V_GB[1]  - V_GA[1];     // taken in G
            const Vec3 p_AB_G_dotdot = A_GB[1]  - A_GA[1];     // taken in G

            const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G; // taken in A, exp. in G
            const Vec3 b_AB_G = A_GB[0] - A_GA[0];             // relative angular acceleration
            const Vec3 a_AB_G = p_AB_G_dotdot - (A_GA[0] % p_AB_G + V_GA[0] % p_AB_G_dot); // taken in A, exp. in G
            const SpatialVec A_AB = ~X_GA.R() * SpatialVec(b_AB_G, a_AB_G); // re-express in A
            sr.setSubtreeBodyAcceleration(sb, A_AB);
        }
    }

    sr.setStage(Stage::Acceleration);
}

void SimbodyMatterSubtree::SubtreeRep::
calcAccelerationsFromSubtreeUDot(const State& s, const Vector& subUDot, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Instance, 
        "calcAccelerationsFromSubtreeUDot()");

    assert(isCompatibleSubtreeResults(sr));

    assert(!"not implemented yet");
}

void SimbodyMatterSubtree::SubtreeRep::
calcAccelerationsFromZeroUDot(const State& s, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Instance, 
        "calcAccelerationsFromZeroUDot()");

    assert(isCompatibleSubtreeResults(sr));

    assert(!"not implemented yet");
}

void SimbodyMatterSubtree::SubtreeRep::
perturbAccelerations(const State& s, SubtreeUIndex subUDotIndex, Real perturbation, SimbodyMatterSubtreeResults::SubtreeResultsRep& sr) const {
    const SimbodyMatterSubsystemRep& matter = getSimbodyMatterSubsystem().getRep();
    SimTK_STAGECHECK_GE_ALWAYS(matter.getStage(s), Stage::Instance, 
        "perturbAccelerations()");

    assert(isCompatibleSubtreeResults(sr));
    assert(sr.getStage() >= Stage::Acceleration);

    assert(!"not implemented yet");
}


} // namespace SimTK

