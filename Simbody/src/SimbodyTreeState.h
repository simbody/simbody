#ifndef SimTK_SIMBODY_TREE_STATE_H_
#define SimTK_SIMBODY_TREE_STATE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

/**@file
 * This file contains the classes which define the SimbodyMatterSubsystem State, 
 * that is, everything that can be changed in a SimbodyMatterSubsystem after 
 * construction.
 *
 * State variables and computation results are organized into stages:
 *    Stage::Empty         virginal state just allocated
 *    Stage::Topology      Stored in SimbodyMatterSubsystem object (construction)
 *   ---------------------------------------------------------
 *    Stage::Model         Stored in the State object
 *    Stage::Instance
 *    Stage::Time
 *    Stage::Position  
 *    Stage::Velocity  
 *    Stage::Dynamics      dynamic operators available
 *    Stage::Acceleration  response to forces in the state
 *   ---------------------------------------------------------
 *    Stage::Report        only used when outputting something
 *
 * Construction proceeds until all the bodies and constraints have been specified. 
 * After that, realizeTopology() is called. Construction-related calculations are
 * performed leading to values which are stored in the SimbodyMatterSubsystem 
 * object, NOT in the State (e.g., total number of bodies). At the same time, an
 * initial state is built, with space allocated for the state variables that will
 * be needed by the next stage (Stage::Model),and these are assigned default values. 
 * Then the stage in the SimbodyMatterSubsystem and in the initial state is set to 
 * "Topology".
 *
 * After that, values for Model stage variables can be set in the State.
 * When that's done we call realizeModel(), which evaluates the Model states
 * putting the values into state cache entries allocated for the purpose. Then
 * all remaining state variables are allocated, and set to their default values.
 * All defaults must be computable knowing only the Model stage values.
 * Then the stage is advanced to Stage::Model.
 *
 * This continues through all the stages, with realizeWhatever() expecting to 
 * receive a state evaluated to stage Whatever-1 equipped with values for stage 
 * Whatever so that it can calculate results and put them in the cache (which is 
 * allocated if necessary), and then advance to stage Whatever. 
 */

#include "simbody/internal/common.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

class SimbodyMatterSubsystemRep;
class RigidBodyNode;
template <int dof> class RigidBodyNodeSpec;

// defined below

class SBTopologyCache;
class SBModelCache;
class SBInstanceCache;
class SBTimeCache;
class SBTreePositionCache;
class SBConstrainedPositionCache;
class SBCompositeBodyInertiaCache;
class SBArticulatedBodyInertiaCache;
class SBTreeVelocityCache;
class SBConstrainedVelocityCache;
class SBDynamicsCache;
class SBTreeAccelerationCache;
class SBConstrainedAccelerationCache;

class SBModelVars;
class SBInstanceVars;
class SBTimeVars;
class SBPositionVars;
class SBVelocityVars;
class SBDynamicsVars;
class SBAccelerationVars;



// -----------------------------------------------------------------------------
//                               TOPOLOGY CACHE
// -----------------------------------------------------------------------------
// An object of this type is stored in the SimbodyMatterSubsystem after extended
// construction in complete, then copied into a slot in the State on
// realizeTopology(). It should contain enough information to size the Model
// stage, and State resource index numbers for the Model-stage state variables
// and Model-stage cache entry. Thise topology cache entry can also contain 
// whatever arbitrary data you would like to have in a State to verify that it 
// is a match for the SimbodyMatterSubsystem.
// 
// Note that this does not include all possible topological information in
// a SimbodyMatterSubsystem -- any subobjects are free to hold their own
// as long as they don't change it after realizeTopology().
class SBTopologyCache {
public:
    SBTopologyCache() {clear();}

    void clear() {
        nBodies = nParticles = nConstraints = nAncestorConstrainedBodies =
            nDOFs = maxNQs = sumSqDOFs = -1;
        modelingVarsIndex.invalidate();
        modelingCacheIndex.invalidate();
        topoInstanceVarsIndex.invalidate();
        valid = false;
    }

    // These are topological objects.
    int nBodies;
    int nParticles;
    int nConstraints;

    // This is the total number of Constrained Bodies appearing in all 
    // Constraints where the Ancestor body is not Ground, excluding the 
    // Ancestor bodies themselves even if they are also Constrained Bodies 
    // (which is common). This is used for sizing pool entries in various 
    // caches to hold precalculated Ancestor-frame data about these bodies.
    int nAncestorConstrainedBodies;

    // TODO: these should be moved to Model stage.
    int nDOFs;
    int maxNQs;
    int sumSqDOFs;

    DiscreteVariableIndex modelingVarsIndex;
    CacheEntryIndex       modelingCacheIndex;

    // These are instance variables that exist regardless of modeling
    // settings; they are instance variables corresponding to topological
    // elements of the matter subsystem (e.g. mobilized bodies and constraints).
    DiscreteVariableIndex topoInstanceVarsIndex;

    bool valid;
};
//.............................. TOPOLOGY CACHE ................................



// -----------------------------------------------------------------------------
//                                MODEL CACHE
// -----------------------------------------------------------------------------
// This cache entry contains counts of various things resulting from the 
// settings of Model-stage state variables. It also contains the resource index
// numbers for all state variable and state cache resources allocated during 
// realizeModel().
//
// Model stage is when the mobilizers settle on the meaning of the q's and u's
// they will employ, so here is where we count up the total number of q's and 
// u's and assign particular slots in those arrays to each mobilizer. We also
// determine the sizes of related "pools", including the number of q's which
// are angles (for sincos calculations), and the number of quaternions in use
// (for normalization calculations), and partition the entries in those
// pools among the mobilizers.
//
// Important things we still don't know at this stage:
//  - what constraints are enabled
//  - how motion will be driven for each coordinate

class SBModelCache {
public:
    SBModelCache() {clear();}

    // Restore this cache entry to its just-constructed condition.
    void clear() {
        totalNQInUse=totalNUInUse=totalNQuaternionsInUse=totalNAnglesInUse = -1;
        mobilizedBodyModelInfo.clear();
    }

    // Allocate the entries in this ModelCache based on information provided in
    // the TopologyCache.
    void allocate(const SBTopologyCache& tc) {
        mobilizedBodyModelInfo.resize(tc.nBodies);
    }

        // MOBILIZED BODY MODELING INFORMATION

    class PerMobilizedBodyModelInfo {
    public:
        PerMobilizedBodyModelInfo() 
          : nQInUse(-1), nUInUse(-1), hasQuaternionInUse(false), nAnglesInUse(-1) {}

        int nQInUse, nUInUse;
        QIndex firstQIndex; // Count from 0 for this SimbodyMatterSubsystem
        UIndex firstUIndex;

        // In case there is a quaternion in use by this Mobilizer. The index 
        // here can be used to find precalculated data associated with this 
        // quaternion, such as its current length.
        bool                hasQuaternionInUse;
        MobilizerQIndex     startOfQuaternion;   // 0..nQInUse-1: which local coordinate starts the quaternion if any?
        QuaternionPoolIndex quaternionPoolIndex; // assigned slot # for this MB's quat, -1 if none

        // In case there are any generalized coordinates which are angles. We 
        // require that the angular coordinates be consecutive and just store 
        // the number of angles and the coordinate index of the first one. The 
        // index can be used to find precalculated data associated with
        // angles, such as their sines and cosines.
        int             nAnglesInUse;
        MobilizerQIndex startOfAngles;  // 0..nQInUse-1: which local coordinate starts angles if any?
        AnglePoolIndex  anglePoolIndex; // start of assigned segment for this MB's angle information, if any
    };

    // Use these accessors so that you get type checking on the index types.
    int getNumMobilizedBodies() const {return (int)mobilizedBodyModelInfo.size();}
    PerMobilizedBodyModelInfo& updMobilizedBodyModelInfo(MobilizedBodyIndex mbx) 
    {   return mobilizedBodyModelInfo[mbx]; }
    const PerMobilizedBodyModelInfo& getMobilizedBodyModelInfo(MobilizedBodyIndex mbx) const
    {   return mobilizedBodyModelInfo[mbx]; }


    // These are sums over the per-MobilizedBody counts above.
    int totalNQInUse, totalNUInUse, totalNQuaternionsInUse, totalNAnglesInUse;

        // STATE ALLOCATION FOR THIS SUBSYSTEM

    // Note that a MatterSubsystem is only one of potentially many users of a 
    // System's State, so only a subset of State variables and State Cache 
    // entries belong to it. Here we record the indices we were given when 
    // we asked the State for some resources. All indices are private to this
    // Subsystem -- they'll start from zero regardless of whether there are
    // other State resource consumers.

    QIndex qIndex;  // NOTE: local, currently always zero
    UIndex uIndex;
    DiscreteVariableIndex timeVarsIndex, qVarsIndex, uVarsIndex, 
                          dynamicsVarsIndex, accelerationVarsIndex;
    CacheEntryIndex       instanceCacheIndex, timeCacheIndex, 
                          treePositionCacheIndex, constrainedPositionCacheIndex,
                          compositeBodyInertiaCacheIndex, articulatedBodyInertiaCacheIndex,
                          treeVelocityCacheIndex, constrainedVelocityCacheIndex,
                          dynamicsCacheIndex, 
                          treeAccelerationCacheIndex, constrainedAccelerationCacheIndex;

private:
    // Use accessor routines for these so that you get type checking on the index types.
    Array_<PerMobilizedBodyModelInfo> mobilizedBodyModelInfo; // MobilizedBody 0 is Ground
};

inline std::ostream& operator<<(std::ostream& o, const SBModelCache& c) { 
    o << "SBModelCache:\n";
    o << "  " << c.getNumMobilizedBodies() << " Mobilized Bodies:\n";
    for (MobilizedBodyIndex mbx(0); mbx < c.getNumMobilizedBodies(); ++mbx) {
        const SBModelCache::PerMobilizedBodyModelInfo& mInfo = c.getMobilizedBodyModelInfo(mbx);
        o << "  " << mbx << ": nq,nu="   << mInfo.nQInUse << "," << mInfo.nUInUse
                         <<  " qix,uix=" << mInfo.firstQIndex << "," << mInfo.firstUIndex << endl;
        if (mInfo.hasQuaternionInUse)
            o <<  "    firstQuat,quatPoolIx=" << mInfo.startOfQuaternion << "," << mInfo.quaternionPoolIndex << endl;
        else o << "    no quaternion in use\n";
        if (mInfo.nAnglesInUse)
             o << "    nangles,firstAngle,anglePoolIx=" << mInfo.nAnglesInUse << "," << mInfo.startOfAngles << "," << mInfo.anglePoolIndex << endl;
        else o << "    no angles in use\n";
    }
    return o; 
}
//............................... MODEL CACHE ..................................



// -----------------------------------------------------------------------------
//                               INSTANCE CACHE
// -----------------------------------------------------------------------------
// This is SimbodyMatterSubsystem information calculated during realizeInstance(), 
// possibly based on the settings of Instance-stage state variables. At this 
// point we will have determined the following information:
//  - final mass properties for all bodies (can calculate total mass)
//  - final locations for all mobilizer frames
//  - which Constraints are enabled
//  - how many and what types of constraint equations are to be included
//  - how motion is to be calculated for each mobilizer
//
// We allocate entries in the constraint error and multiplier pools among the 
// Constraints, and allocate entries in the prescribed motion and prescribed
// force pools among mobilizers whose motion is prescribed.
//
// At this point we can classify all the mobilizers based on the kind of Motion 
// they will undergo. We determine the scope of every Constraint, and classify
// them based on the kinds of mobilizers they affect.
class SBInstanceCache {
public:
        // MOBILIZED BODY INSTANCE INFORMATION

    // This is information calculated once we have seen all the Instance-stage
    // State variable values that can affect bodies, mobilizers, and motions.
    // Notes: 
    //   - all mobilities of a mobilizer must be treated identically
    //   - if any motion level is fast, then the whole mobilizer is fast
    //   - if a mobilizer is fast, so are all its outboard mobilizers
    //
    struct PerMobodInstanceInfo {
        PerMobodInstanceInfo() {clear();}

        void clear() {   
            qMethod=uMethod=udotMethod=Motion::Free;
            firstPresQ.invalidate(); firstPresU.invalidate(); 
            firstPresUDot.invalidate(); firstPresForce.invalidate();
        }

        Motion::Method      qMethod;        // how are positions calculated?
        Motion::Method      uMethod;        // how are velocities calculated?
        Motion::Method      udotMethod;     // how are accelerations calculated?

        PresQPoolIndex      firstPresQ;     // if qMethod==Prescribed
        PresUPoolIndex      firstPresU;     // if uMethod==Prescribed
        PresUDotPoolIndex   firstPresUDot;  // if udotMethod==Prescribed
        PresForcePoolIndex  firstPresForce; // if udotMethod!=Free
    };

        // CONSTRAINT INSTANCE INFORMATION

    // Store some Instance-stage information about each Constraint. Most important, we don't
    // know how many constraint equations (if any) the Constraint will generate until
    // Instance stage. In particular, a disabled Constraint won't generate any equations (it
    // will have an Info entry here, however). Also, although we know the Constrained
    // Mobilizers at Topology stage, we don't know the specific number or types of internal
    // coordinates involved until Instance stage.

    struct PerConstrainedMobilizerInstanceInfo {
        PerConstrainedMobilizerInstanceInfo() : nQInUse(0), nUInUse(0) { } // assume disabled
        // The correspondence between Constrained Mobilizers and Mobilized Bodies is
        // Topological information you can pull from the TopologyCache.
        // See the MobilizedBody for counts of its q's and u's, which define the allocated
        // number of slots for the ConstrainedMobilizer as well.
        int nQInUse, nUInUse; // same as corresponding MobilizedBody unless disabled
        ConstrainedQIndex  firstConstrainedQIndex; // these count from 0 for each Constraint
        ConstrainedUIndex  firstConstrainedUIndex;
    };
        
    class PerConstraintInstanceInfo {
    public:
        PerConstraintInstanceInfo() { }
        void clear() {
            constrainedMobilizerInstanceInfo.clear();
            constrainedQ.clear(); constrainedU.clear();
            participatingQ.clear(); participatingU.clear();
        }

        void allocateConstrainedMobilizerInstanceInfo(int nConstrainedMobilizers) {
            assert(nConstrainedMobilizers >= 0);
            constrainedMobilizerInstanceInfo.resize(nConstrainedMobilizers);
            constrainedQ.clear();   // build by appending
            constrainedU.clear();
        }

        int getNumConstrainedMobilizers() const {return (int)constrainedMobilizerInstanceInfo.size();}
        const PerConstrainedMobilizerInstanceInfo& getConstrainedMobilizerInstanceInfo(ConstrainedMobilizerIndex M) const {
            return constrainedMobilizerInstanceInfo[M];
        }
        PerConstrainedMobilizerInstanceInfo& updConstrainedMobilizerInstanceInfo(ConstrainedMobilizerIndex M) {
            return constrainedMobilizerInstanceInfo[M];
        }
        
        int getNumConstrainedQ() const {return (int)constrainedQ.size();}
        int getNumConstrainedU() const {return (int)constrainedU.size();}
        ConstrainedQIndex addConstrainedQ(QIndex qx) {
            constrainedQ.push_back(qx);
            return ConstrainedQIndex(constrainedQ.size()-1);
        }
        ConstrainedUIndex addConstrainedU(UIndex ux) {
            constrainedU.push_back(ux);
            return ConstrainedUIndex(constrainedU.size()-1);
        }
        QIndex getQIndexFromConstrainedQ(ConstrainedQIndex i) const {return constrainedQ[i];}
        UIndex getUIndexFromConstrainedU(ConstrainedUIndex i) const {return constrainedU[i];}

        int getNumParticipatingQ() const {return (int)participatingQ.size();}
        int getNumParticipatingU() const {return (int)participatingU.size();}
        ParticipatingQIndex addParticipatingQ(QIndex qx) {
            participatingQ.push_back(qx);
            return ParticipatingQIndex(participatingQ.size()-1);
        }
        ParticipatingUIndex addParticipatingU(UIndex ux) {
            participatingU.push_back(ux);
            return ParticipatingUIndex(participatingU.size()-1);
        }
        QIndex getQIndexFromParticipatingQ(ParticipatingQIndex i) const {return participatingQ[i];}
        UIndex getUIndexFromParticipatingU(ParticipatingUIndex i) const {return participatingU[i];}

        Segment holoErrSegment;    // (offset,mHolo)    for each Constraint, within subsystem qErr
        Segment nonholoErrSegment; // (offset,mNonholo) same, but for uErr slots (after holo derivs)
        Segment accOnlyErrSegment; // (offset,mAccOnly) same, but for udotErr slots (after holo/nonholo derivs)
    public:
        // Better to access using accessor methods above so you'll get type 
        // checking on the index type.
        Array_<PerConstrainedMobilizerInstanceInfo> 
            constrainedMobilizerInstanceInfo;

        // The ConstrainedBodies and ConstrainedMobilizers are set at Topology 
        // stage, but the particular generalized coordinates q and generalized 
        // speeds u which are involved can't be determined until Model stage, 
        // since the associated mobilizers have Model stage options which can 
        // affect the number and meanings of these variables. These are sorted 
        // in order of their associated ConstrainedMobilizer, not necessarily
        // in order of QIndex or UIndex. Each value appears only once.
        Array_<QIndex> constrainedQ;   // indexed by ConstrainedQIndex, maps to subsystem QIndex
        Array_<UIndex> constrainedU;   // indexed by ConstrainedUIndex, maps to subsystem UIndex

        // Participating mobilities include ALL the mobilities which may be 
        // involved in any of this Constraint's constraint equations, whether 
        // from being directly constrained or indirectly as a result of their 
        // effects on ConstrainedBodies. These are sorted in order of 
        // increasing QIndex and UIndex, and each QIndex or UIndex appears 
        // only once.
        Array_<QIndex> participatingQ; // indexed by ParticipatingQIndex, maps to subsystem QIndex
        Array_<UIndex> participatingU; // indexed by ParticipatingUIndex, maps to subsystem UIndex
    };

    // Instance variables are:
    //   body mass props; particle masses
    //   X_BM, X_PF mobilizer transforms
    //  
    // Calculations stored here derive from those states:
    //   total mass
    //   central inertia of each rigid body
    //   principal axes and corresponding principal moments of inertia of 
    //       each rigid body
    //   reference configuration X_PB when q==0 (usually that means M==F), 
    //       for each rigid body

    Real              totalMass; // sum of all rigid body and particles masses
    Array_<Inertia>   centralInertias;           // nb
    Vector_<Vec3>     principalMoments;          // nb
    Array_<Rotation>  principalAxes;             // nb
    Array_<Transform> referenceConfiguration;    // nb

    Array_<PerMobodInstanceInfo>         mobodInstanceInfo;
    Array_<PerConstraintInstanceInfo>    constraintInstanceInfo;

    // This is a sum over all the mobilizers whose q's are currently prescribed,
    // adding the number of q's (generalized coordinates) nq currently being 
    // used for each of those. An array of size totalNPresQ is allocated in the 
    // TimeCache to hold the calculated q's (which will be different from the 
    // actual q's until they are applied). Motions will also provide this many 
    // prescribed qdots and qdotdots, but we will map those to u's and udots 
    // before recording them, with nu entries being allocated in each. These 
    // nq- and nu-sized slots are allocated in order of MobilizedBodyIndex.
    int getTotalNumPresQ() const {return (int)presQ.size();}
    int getTotalNumZeroQ() const {return (int)zeroQ.size();}
    int getTotalNumFreeQ() const {return (int)freeQ.size();}
    Array_<QIndex> presQ;
    Array_<QIndex> zeroQ;
    Array_<QIndex> freeQ; // must be integrated

    // This is a sum over all the mobilizers whose u's are current prescribed, 
    // whether because of non-holonomic (velocity) prescribed motion u=u(t,q), 
    // or because the q's are prescribed via holonomic (position) prescribed 
    // motion and the u's are calculated from the qdots. We add the number u's 
    // (generalized speeds) nu currently being used for each holonomic- or 
    // nonholonomic-prescribed mobilizer. An array of this size is allocated 
    // in the PositionCache to hold the calculated u's (which will be 
    // different from the actual u's until they are applied). Nu-sized slots 
    // are allocated in order of MobilizedBodyIndex.
    int getTotalNumPresU() const {return (int)presU.size();}
    int getTotalNumZeroU() const {return (int)zeroU.size();}
    int getTotalNumFreeU() const {return (int)freeU.size();}
    Array_<UIndex> presU;
    Array_<UIndex> zeroU;
    Array_<UIndex> freeU; // must be integrated

    // This is a sum over all the mobilizers whose udots are currently prescribed, adding
    // the number of udots (mobilities) nu from each holonomic-, nonholonomic-, or 
    // acceleration-prescribed mobilizer. An array of this size is allocated
    // in the DynamicsCache, and an entry is needed in the prescrived force array
    // in the AccelerationCache as well. These nu-sized slots are allocated in 
    // order of MobilizedBodyIndex.
    int getTotalNumPresUDot() const {return (int)presUDot.size();}
    int getTotalNumZeroUDot() const {return (int)zeroUDot.size();}
    int getTotalNumFreeUDot() const {return (int)freeUDot.size();}
    Array_<UIndex> presUDot;
    Array_<UIndex> zeroUDot;
    Array_<UIndex> freeUDot; // calculated from forces

    // This is a sum over all the mobilizers whose udots are known for any reason
    // whether prescribed or not, e.g. they are Zero or Discrete (anything
    // but Free). These need slots in the array of calculated prescribed motion forces.
    int totalNPresForce;

    // Quaternion errors go in qErr also, but after all the physical contraint errors. That is,
    // they start at index totalNHolonomicConstraintEquationsInUse.
    int firstQuaternionQErrSlot;

    // These record where in the full System's State our Subsystem's qErr, uErr, and udotErr
    // entries begin. That is, this subsystem's segments can be found at
    //    qErr   (qErrIndex,    nPositionConstraintEquationsInUse + nQuaternionsInUse)
    //    uErr   (uErrIndex,    nVelocityConstraintEquationsInUse)
    //    udotErr(udotErrIndex, nAccelerationConstraintEquationsInUse)
    int qErrIndex, uErrIndex, udotErrIndex;

    // These are the sums over the per-Constraint data above. The number of
    // position constraint equations (not counting quaternion normalization constraints)
    // is the same as the number of holonomic constraints mHolo. The number of velocity
    // constraint equations is mHolo+mNonholo. The number of acceleration constraints,
    // and thus the number of multipliers, is mHolo+mNonholo+mAccOnly.
    int totalNHolonomicConstraintEquationsInUse;         // sum(mHolo)    (#position equations = mHolo)
    int totalNNonholonomicConstraintEquationsInUse;      // sum(mNonholo) (#velocity equations = mHolo+mNonholo)
    int totalNAccelerationOnlyConstraintEquationsInUse;  // sum(mAccOnly) (#acceleration eqns  = mHolo+mNonholo+mAccOnly)

public:
    void allocate(const SBTopologyCache& topo,
                  const SBModelCache&    model) 
    {
        totalMass = SimTK::NaN;
        centralInertias.resize(topo.nBodies);           // I_CB
        principalMoments.resize(topo.nBodies);          // (Ixx,Iyy,Izz)
        principalAxes.resize(topo.nBodies);             // [axx ayy azz]
        referenceConfiguration.resize(topo.nBodies);    // X0_PB

        mobodInstanceInfo.resize(topo.nBodies);
        totalNPresForce = 0;

        constraintInstanceInfo.resize(topo.nConstraints);
        firstQuaternionQErrSlot = qErrIndex = uErrIndex = udotErrIndex = -1;
        totalNHolonomicConstraintEquationsInUse        = 0;
        totalNNonholonomicConstraintEquationsInUse     = 0;
        totalNAccelerationOnlyConstraintEquationsInUse = 0;
    }

    int getNumMobilizedBodies() const {return (int)mobodInstanceInfo.size();}
    PerMobodInstanceInfo& updMobodInstanceInfo(MobilizedBodyIndex mbx)
    {   return mobodInstanceInfo[mbx]; }
    const PerMobodInstanceInfo& getMobodInstanceInfo(MobilizedBodyIndex mbx) const
    {   return mobodInstanceInfo[mbx]; }

    int getNumConstraints() const {return (int)constraintInstanceInfo.size();}
    PerConstraintInstanceInfo& updConstraintInstanceInfo(ConstraintIndex cx)
    {   return constraintInstanceInfo[cx]; }
    const PerConstraintInstanceInfo& getConstraintInstanceInfo(ConstraintIndex cx) const
    {   return constraintInstanceInfo[cx]; }
};
//.............................. INSTANCE CACHE ................................



// -----------------------------------------------------------------------------
//                                TIME CACHE
// -----------------------------------------------------------------------------
// Here we hold information that is calculated in the SimbodyMatterSubsystem's
// realizeTime() method. Currently that consists only of prescribed q's, which
// must always be defined as functions of time.

class SBTimeCache {
public:
    // This holds values from Motion prescribed position (holonomic) calculations.
    Vector presQPool;   // Index with PresQPoolIndex

public:
    void allocate(const SBTopologyCache& topo,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        presQPool.resize(instance.getTotalNumPresQ());
    }
};
//............................... TIME CACHE ...................................



// -----------------------------------------------------------------------------
//                             TREE POSITION CACHE
// -----------------------------------------------------------------------------
// Here we hold information that is calculated early in the matter subsystem's
// realizePosition() method. This includes
//
//  - mobilizer matrices H_FM, X_FM, X_PB, H_PB_G 
//  - basic kinematic information X_GB, Phi_PB_G
//  - for constrained bodies, position X_AB of each body in its ancestor A
//  - mass properties expressed in Ground (TODO: these should probably be in
//    their own cache since they aren't needed for kinematics)
//
// This cache entry can be calculated after Stage::Time and is guaranteed to have
// been calculated by the end of Stage::Position. The SimbodyMatterSubsystem's
// realizePosition() method will mark this done as soon as possible, so that
// later calculations (constraint position errors, prescribed velocities) can
// access these without a stage violation.

class SBTreePositionCache {
public:
    const Transform& getX_FM(MobilizedBodyIndex mbx) const {return bodyJointInParentJointFrame[mbx];}
    Transform&       updX_FM(MobilizedBodyIndex mbx)       {return bodyJointInParentJointFrame[mbx];}
    const Transform& getX_PB(MobilizedBodyIndex mbx) const {return bodyConfigInParent[mbx];}
    Transform&       updX_PB(MobilizedBodyIndex mbx)       {return bodyConfigInParent[mbx];}
    const Transform& getX_GB(MobilizedBodyIndex mbx) const {return bodyConfigInGround[mbx];}
    Transform&       updX_GB(MobilizedBodyIndex mbx)       {return bodyConfigInGround[mbx];}

    const Transform& getX_AB(AncestorConstrainedBodyPoolIndex cbpx) const 
    {   return constrainedBodyConfigInAncestor[cbpx]; }
    Transform&       updX_AB(AncestorConstrainedBodyPoolIndex cbpx)
    {   return constrainedBodyConfigInAncestor[cbpx]; }
public:
    // qerr cache space is provided directly by the State

    Vector sq, cq;  // nq  Sin&cos of angle q's in appropriate slots; otherwise garbage
    Vector qnorm;   // nq  Contains normalized quaternions in appropriate slots;
                    //       all else is garbage.

    // CAUTION: our definition of the H matrix is transposed from those used
    // by Jain and by Schwieters. Jain would call these H* and Schwieters
    // would call them H^T, but we call them H.
    Matrix_<Vec3> storageForH_FM; // 2 x ndof (H_FM)
    Matrix_<Vec3> storageForH;    // 2 x ndof (H_PB_G)

    Array_<Transform>    bodyJointInParentJointFrame;  // nb (X_FM)
    Array_<Transform>    bodyConfigInParent;           // nb (X_PB)
    Array_<Transform>    bodyConfigInGround;           // nb (X_GB)
    Array_<PhiMatrix>    bodyToParentShift;            // nb (phi)
    Array_<Inertia>      bodyInertiaInGround;          // nb (I_OB_G)
    Array_<SpatialInertia>    bodySpatialInertia;           // nb (Mk)
    Vector_<Vec3>             bodyCOMInGround;              // nb (p_G_CB)
    Vector_<Vec3>             bodyCOMStationInGround;       // nb (p_CB_G)


        // Ancestor Constrained Body Pool

    // For Constraints whose Ancestor body A is not Ground G, we assign pool entries
    // for each of their Constrained Bodies (call the total number 'nacb')
    // to store the above information but measured and expressed in the Ancestor frame
    // rather than Ground.
    Array_<Transform> constrainedBodyConfigInAncestor;     // nacb (X_AB)

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;   // this is the number of u's (nu)
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int nacb    = tree.nAncestorConstrainedBodies;

        // These contain uninitialized junk. Body-indexed entries get their
        // ground elements set appropriately now and forever.
        sq.resize(maxNQs);
        cq.resize(maxNQs);
        qnorm.resize(maxNQs);

        storageForH_FM.resize(2,nDofs);
        storageForH.resize(2,nDofs);

        bodyJointInParentJointFrame.resize(nBodies); 
        bodyJointInParentJointFrame[0].setToZero();

        bodyConfigInParent.resize(nBodies);          
        bodyConfigInParent[0].setToZero();

        bodyConfigInGround.resize(nBodies);          
        bodyConfigInGround[0].setToZero();

        bodyToParentShift.resize(nBodies);           
        bodyToParentShift[0].setToZero();

        bodyInertiaInGround.resize(nBodies); // TODO: ground initialization
        bodySpatialInertia.resize(nBodies);  // TODO: ground initialization

        bodyCOMInGround.resize(nBodies);             
        bodyCOMInGround[0] = 0.;

        bodyCOMStationInGround.resize(nBodies);      
        bodyCOMStationInGround[0] = 0.;

        constrainedBodyConfigInAncestor.resize(nacb);
    }
};
//.......................... TREE POSITION CACHE ...............................



// -----------------------------------------------------------------------------
//                         CONSTRAINED POSITION CACHE 
// -----------------------------------------------------------------------------
// Here we hold information that is part of the matter subsystem's 
// realizePosition() calculation but depends on the TreePositionCache having
// already been calculated. This includes:
//
//  - desired values of prescribed u's (since those are functions of at most
//    time and position)
//  - logically, position constraint errors (qerrs), although in fact that
//    array is provided as a built-in by the State
//
// This cache entry can be calculated after Stage::Time provided that
// the SBTreePositionCache entry has already been marked valid. We guarantee
// this will have been calculated by Stage::Position.

class SBConstrainedPositionCache {
public:
    // qerr cache space is provided directly by the State

    // This holds values from all the Motion prescribed velocity (nonholonomic) 
    // calculations, and those resulting from diffentiating prescribed positions.
    Vector presUPool;   // Index with PresUPoolIndex

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        presUPool.resize(instance.getTotalNumPresU());
    }
};
//........................ CONSTRAINED POSITION CACHE ..........................



// -----------------------------------------------------------------------------
//                        COMPOSITE BODY INERTIA CACHE
// -----------------------------------------------------------------------------
// Composite body inertias R are those that would be felt if all the mobilizers
// had prescribed motion (or were welded in their current configurations). These
// are convenient for inverse dynamics computations and for scaling of 
// generalized coordinates and speeds.
//
// Each spatial inertia here is expressed in the Ground frame but measured about 
// its body's origin.
//
// Note that each composite body inertia is a rigid-body spatial inertia, not
// the more complicated articulated-body spatial inertia. That means these have
// a scalar mass and well-defined mass center, and a very simple structure which
// can be exploited for speed. There are at most 10 unique elements in a rigid
// body spatial inertia matrix.
//
// Composite body inertias depend only on positions but are often not needed at 
// all. So we give them their own cache entry and expect explicit realization 
// some time after Position stage, if at all.

class SBCompositeBodyInertiaCache {
public:
    Array_<SpatialInertia> compositeBodyInertia;     // nb (R)

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies; 
        
        compositeBodyInertia.resize(nBodies); // TODO: ground initialization
    }
};
//....................... COMPOSITE BODY INERTIA CACHE .........................



// -----------------------------------------------------------------------------
//                       ARTICULATED BODY INERTIA CACHE
// -----------------------------------------------------------------------------
// These articulated body inertias take into account prescribed motion, 
// meaning that they are produced by a combination of articulated and rigid
// shift operations depending on each mobilizer's current stats as "regular"
// or "prescribed". That means that the articulated inertias here are suited
// only for "mixed" dynamics; you can't use them to calculate M^-1*f unless
// there is no prescribed motion in the system.
//
// Each articulated body inertia here is expressed in the Ground frame but 
// measured about its body's origin.
//
// Articulated body inertia matrices, though symmetric and positive
// definite, do not have the same simple structure as rigid-body (or composite-
// body) spatial inertias. For example, the apparent mass depends on direction.
// All 21 elements of this symmetric 6x6 matrix are unique, while there are only
// 10 unique elements in a rigid body spatial inertia.
//
// Note that although we use some *rigid* body shift operations here, the 
// results in general are all *articulated* body inertias, because a rigid shift 
// of an articulated body inertia is still an articulated body inertia. Only if 
// all mobilizers are prescribed will these be rigid body spatial inertias. For a 
// discussion of the properties of articulated body inertias, see Section 7.1 
// (pp. 119-123) of Roy Featherstone's excellent 2008 book, Rigid Body Dynamics 
// Algorithms. 
//
// Intermediate quantities PPlus, D, DI, and G are calculated here
// which are separately useful when dealing with "regular" mobilized bodies. 
// These quantities are not calculated for prescribed mobilizers; they will 
// remain NaN in that case. In
// particular, this means that the prescribed-mobilizer mass properties do not 
// have to be invertible, so you can have terminal massless bodies as long as 
// their motion is always prescribed.
// TODO: should D still be calculated? It doesn't require inversion.
//
// Articulated body inertias depend only on positions but are not usually needed 
// until Acceleration stage. Thus this cache entry should have dependsOn stage 
// Position, and computedBy stage Dynamics. However, it can be realized any
// time after Position.

class SBArticulatedBodyInertiaCache {
public:
    Array_<ArticulatedInertia> articulatedBodyInertia; // nb (P)

    Array_<ArticulatedInertia> pPlus;              // nb

    Vector_<Real>       storageForD;              // sum(nu[j]^2)
    Vector_<Real>       storageForDI;             // sum(nu[j]^2)
    Matrix_<Vec3>       storageForG;              // 2 X ndof

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs;   // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need     
        
        articulatedBodyInertia.resize(nBodies); // TODO: ground initialization

        pPlus.resize(nBodies); // TODO: ground initialization

        storageForD.resize(nSqDofs);
        storageForDI.resize(nSqDofs);
        storageForG.resize(2,nDofs);
    }
};
//....................... ARTICULATED BODY INERTIA CACHE .......................



// -----------------------------------------------------------------------------
//                              TREE VELOCITY CACHE
// -----------------------------------------------------------------------------
// Here we hold information that is calculated early in the matter subsystem's
// realizeVelocity() method. This includes
//
//  - mobilizer matrices HDot_FM and HDot_PB_G
//  - cross mobilizer velocities V_FM, V_PB
//  - basic kinematics V_GB giving body velocities in Ground
//  - logically, qdot, but that is provided as a built-in cache entry in State
//  - for constrained bodies, V_AB giving body velocities in their ancestor A
//  - velocity-dependent dynamics remainder terms: coriolis acceleration and
//    gyroscopic forces
//
// This cache entry can be calculated after Stage::Position and is guaranteed to 
// have been calculated by the end of Stage::Velocity. The matter subsystem's
// realizeVelocity() method will mark this done as soon as possible, so that
// later calculations (constraint velocity errors) can access these without a 
// stage violation.

class SBTreeVelocityCache {
public:
    const SpatialVec& getV_FM(MobilizedBodyIndex mbx) const {return mobilizerRelativeVelocity[mbx];}
    SpatialVec&       updV_FM(MobilizedBodyIndex mbx)       {return mobilizerRelativeVelocity[mbx];}
    const SpatialVec& getV_PB(MobilizedBodyIndex mbx) const {return bodyVelocityInParent[mbx];}
    SpatialVec&       updV_PB(MobilizedBodyIndex mbx)       {return bodyVelocityInParent[mbx];}
    const SpatialVec& getV_GB(MobilizedBodyIndex mbx) const {return bodyVelocityInGround[mbx];}
    SpatialVec&       updV_GB(MobilizedBodyIndex mbx)       {return bodyVelocityInGround[mbx];}

    const SpatialVec& getV_AB(AncestorConstrainedBodyPoolIndex cbpx) const 
    {   return constrainedBodyVelocityInAncestor[cbpx]; }
    SpatialVec&       updV_AB(AncestorConstrainedBodyPoolIndex cbpx)       
    {   return constrainedBodyVelocityInAncestor[cbpx]; }

public:
    // qdot cache space is supplied directly by the State

    Vector_<SpatialVec> mobilizerRelativeVelocity; // nb (V_FM) cross-mobilizer velocity
    Vector_<SpatialVec> bodyVelocityInParent;      // nb (V_PB) 
    Vector_<SpatialVec> bodyVelocityInGround;      // nb (V_GB)

    // CAUTION: our definition of the H matrix is transposed from those used
    // by Jain and by Schwieters.
    Matrix_<Vec3> storageForHDot_FM;  // 2 x ndof (HDot_FM)
    Matrix_<Vec3> storageForHDot;     // 2 x ndof (HDot_PB_G)

    Vector_<SpatialVec> bodyVelocityInParentDerivRemainder; // VB_PB_G=HDot_PB_G*u

    Vector_<SpatialVec> coriolisAcceleration;     // nb (a)
    Vector_<SpatialVec> totalCoriolisAcceleration;// nb (A)
    Vector_<SpatialVec> gyroscopicForces;         // nb (b)

        // Ancestor Constrained Body Pool

    // For Constraints whose Ancestor body A is not Ground G, we assign pool entries
    // for each of their Constrained Bodies (call the total number 'nacb')
    // to store the above information but measured and expressed in the Ancestor frame
    // rather than Ground.
    Array_<SpatialVec> constrainedBodyVelocityInAncestor; // nacb (V_AB)

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int nacb    = tree.nAncestorConstrainedBodies;

        mobilizerRelativeVelocity.resize(nBodies);       
        mobilizerRelativeVelocity[0] = SpatialVec(Vec3(0),Vec3(0));

        bodyVelocityInParent.resize(nBodies);       
        bodyVelocityInParent[0] = SpatialVec(Vec3(0),Vec3(0));

        bodyVelocityInGround.resize(nBodies);       
        bodyVelocityInGround[0] = SpatialVec(Vec3(0),Vec3(0));

        storageForHDot_FM.resize(2,nDofs);
        storageForHDot.resize(2,nDofs);

        bodyVelocityInParentDerivRemainder.resize(nBodies);       
        bodyVelocityInParentDerivRemainder[0] = SpatialVec(Vec3(0),Vec3(0));

        coriolisAcceleration.resize(nBodies);       
        coriolisAcceleration[0] = SpatialVec(Vec3(0),Vec3(0));

        totalCoriolisAcceleration.resize(nBodies);       
        totalCoriolisAcceleration[0] = SpatialVec(Vec3(0),Vec3(0));

        gyroscopicForces.resize(nBodies);           
        gyroscopicForces[0] = SpatialVec(Vec3(0),Vec3(0));

        constrainedBodyVelocityInAncestor.resize(nacb);
    }
};
//............................ TREE VELOCITY CACHE .............................



//------------------------------------------------------------------------------
//                         CONSTRAINED VELOCITY CACHE 
//------------------------------------------------------------------------------
// Here we hold information that is part of the matter subsystem's 
// realizeVelocity() calculation but depends on the TreeVelocityCache having
// already been calculated. This includes:
//
//  - (not prescribed udots because we delay those until Dynamics)
//  - logically, velocity constraint errors (uerrs), although in fact that
//    array is provided as a built-in by the State
//
// This cache entry can be calculated after Stage::Position provided that
// the SBTreeVelocityCache entry has already been marked valid. We guarantee
// this will have been calculated by the end of Stage::Velocity.
//
// TODO: currently there is nothing here
class SBConstrainedVelocityCache {
public:
    // uerr cache space is provided directly by the State

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        // nothing yet
    }
};
//........................ CONSTRAINED VELOCITY CACHE ..........................



// -----------------------------------------------------------------------------
//                                DYNAMICS CACHE
// -----------------------------------------------------------------------------
class SBDynamicsCache {
public:
    // This holds the values from all the Motion prescribed acceleration 
    // calculations, and those which result from diffentiating prescribed velocities,
    // or differentiating twice prescribed positions.
    Vector presUDotPool;    // Index with PresUDotPoolIndex

    // Dynamics
    Vector_<SpatialVec> centrifugalForces;        // nb (P*a+b)
    Vector_<SpatialVec> totalCentrifugalForces;   // nb (P*A+b)

    Vector_<SpatialMat> Y;                        // nb

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs; // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;    // allocate the max # q's we'll ever need     

        presUDotPool.resize(instance.getTotalNumPresUDot());
        
        centrifugalForces.resize(nBodies);           
        centrifugalForces[0] = SpatialVec(Vec3(0),Vec3(0));

        totalCentrifugalForces.resize(nBodies);           
        totalCentrifugalForces[0] = SpatialVec(Vec3(0),Vec3(0));

        Y.resize(nBodies);
        Y[0] = SpatialMat(Mat33(0));
    }
};
//............................... DYNAMICS CACHE ...............................



// -----------------------------------------------------------------------------
//                          TREE ACCELERATION CACHE
// -----------------------------------------------------------------------------
// Here we hold information that is calculated early in the matter subsystem's
// realizeAcceleration() method. This includes
//
//  - basic kinematics A_GB giving body accelerations in Ground
//  - for constrained bodies, A_AB giving body accelerations in their ancestor A
//  - prescribed motion forces
//  - logically, udot and qdotdot, but those arrays are provided as built-in 
//    cache entries in State
//
// This cache entry can be calculated after Stage::Dynamics and is guaranteed to 
// have been calculated by the end of Stage::Acceleration. The matter subsystem's
// realizeAcceleration() method will mark this done as soon as possible, so that
// later calculations (constraint acceleration errors) can access these without a 
// stage violation.

class SBTreeAccelerationCache {
public:
    const SpatialVec& getA_GB(MobilizedBodyIndex mbx) const {return bodyAccelerationInGround[mbx];}
    SpatialVec&       updA_GB(MobilizedBodyIndex mbx)       {return bodyAccelerationInGround[mbx];}

    const SpatialVec& getA_AB(AncestorConstrainedBodyPoolIndex cbpx) const 
    {   return constrainedBodyAccelerationInAncestor[cbpx]; }
    SpatialVec&       updA_AB(AncestorConstrainedBodyPoolIndex cbpx)       
    {   return constrainedBodyAccelerationInAncestor[cbpx]; }

public:
    // udot, qdotdot cache space is provided directly by the State.


    Vector_<SpatialVec> bodyAccelerationInGround; // nb (A_GB)

    // This is where the calculated prescribed motion "lambdas" go. (That is, 
    // generalized forces needed to implement prescribed generalized 
    // accelerations.) Slots here are doled out only for mobilizers that have 
    // known accelerations; there is one scalar here per mobility in those 
    // mobilizers. Look in the InstanceCache to see which slots are allocated
    // to which mobilizers.
    Vector presMotionForces;    // Index with PresForcePoolIndex

    // Temps used in calculating accelerations and prescribed forces.
    Vector              epsilon;                  // nu
    Vector_<SpatialVec> z;                        // nb
    Vector_<SpatialVec> Gepsilon;                 // nb

        // Ancestor Constrained Body Pool

    // For Constraints whose Ancestor body A is not Ground G, we assign pool entries
    // for each of their Constrained Bodies (call the total number 'nacb')
    // to store the above information but measured and expressed in the Ancestor frame
    // rather than Ground.
    Array_<SpatialVec> constrainedBodyAccelerationInAncestor; // nacb (A_AB)

public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
        // Pull out topology-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs; // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;    // allocate the max # q's we'll ever need
        const int nacb    = tree.nAncestorConstrainedBodies;

        bodyAccelerationInGround.resize(nBodies);   
        bodyAccelerationInGround[0] = SpatialVec(Vec3(0),Vec3(0));;

        presMotionForces.resize(instance.totalNPresForce);

        epsilon.resize(nDofs);
        z.resize(nBodies);
        Gepsilon.resize(nBodies); // TODO: ground initialization

        constrainedBodyAccelerationInAncestor.resize(nacb);
    }
};
//.......................... TREE ACCELERATION CACHE ...........................



//------------------------------------------------------------------------------
//                        CONSTRAINED ACCELERATION CACHE 
//------------------------------------------------------------------------------
// Here we hold information that is part of the matter subsystem's 
// realizeAcceleration() calculation but depends on the TreeAccelerationCache 
// having already been calculated. This includes:
//

//  - logically, acceleration constraint errors (udoterrs), and constraint
//    multipliers, although in fact those arrays are provided as built-ins by 
//    the State
//
// This cache entry can be calculated after Stage::Dynamics provided that
// the SBTreeAccelerationCache entry has already been marked valid. We guarantee
// this will have been calculated by the end of Stage::Acceleration.

class SBConstrainedAccelerationCache {
public:
    // udoterr and multiplier cache space is provided directly by the State.


public:
    void allocate(const SBTopologyCache& tree,
                  const SBModelCache&    model,
                  const SBInstanceCache& instance) 
    {
    }
};
//...................... CONSTRAINED ACCELERATION CACHE ........................




/* 
 * Generalized state variable collection for a SimbodyMatterSubsystem. 
 * Variables are divided into Stages, according to when their values
 * are needed during a calculation. The Stages are:
 *       (Topology: not part of the state. These are the bodies, mobilizers,
 *        and topological constraints.)
 *     Model:         choice of coordinates, knowns & unknowns, methods, etc.
 *     Instance:      setting of physical parameters, e.g. mass
 *       (Time: currently there are no time-dependent states or computations)
 *     Position:      position and orientation values q (2nd order continuous)
 *     Velocity:      generalized speeds u
 *     Dynamics:      dynamic quantities & operators available
 *     Acceleration:  applied forces and prescribed accelerations
 *     Report:        used by study for end-user reporting only; no effect on 
 *                      results
 *
 */


// -----------------------------------------------------------------------------
//                                 MODEL VARS
// -----------------------------------------------------------------------------
// This state variable is allocated during realizeTopology(). Any change made
// to it after that invalidates Stage::Model, requiring realizeModel() to be
// performed.
class SBModelVars {
public:
    bool         useEulerAngles;
    Array_<bool> prescribed;           // nb (# bodies & mobilizers, [0] always true)
public:

    // We have to allocate these without looking at any other
    // state variable or cache entries. We can only depend on topological
    // information.
    void allocate(const SBTopologyCache& tree) {
        useEulerAngles = false;
        prescribed.resize(tree.nBodies, false); 
    }

};



// -----------------------------------------------------------------------------
//                               INSTANCE VARS
// -----------------------------------------------------------------------------
// This state variable is allocated during realizeTopology(), because its 
// contents refer only to elements that form part of the fixed topology of the
// matter subsystem -- mobilized bodies, particles, and constraints that are
// specified as permanent parts of this matter subsystem.
// 
// Any change to this variable invalidates Stage::Instance (not Stage::Model), 
// requiring realize(Instance) to be performed.
//
// Note: we may at some point have instance variables whose allocation is
// deferred until realizeModel() but those would be wiped out whenever a change
// to a Model-stage variable is made (most notably useEulerAngles).
class SBInstanceVars {
public:
    Array_<MassProperties> bodyMassProperties;
    Array_<Transform>      outboardMobilizerFrames;
    Array_<Transform>      inboardMobilizerFrames;
    Vector                 particleMasses;
    Array_<bool>           disabled;             // nc (# constraints)

public:

    // We can access the tree or state variable & cache up to Modeling stage.
    void allocate(const SBTopologyCache& topology) const {
        SBInstanceVars& mutvars = *const_cast<SBInstanceVars*>(this);
        mutvars.bodyMassProperties.resize(topology.nBodies);
        mutvars.outboardMobilizerFrames.resize(topology.nBodies);
        mutvars.inboardMobilizerFrames.resize (topology.nBodies);
        mutvars.particleMasses.resize(topology.nParticles);
        mutvars.disabled.resize(topology.nConstraints, false);

        // Set default values
        for (int i = 0; i < (int)mutvars.bodyMassProperties.size(); ++i)
            mutvars.bodyMassProperties[i] = MassProperties(1,Vec3(0),Inertia(1));
        for (int i = 0; i < (int)mutvars.outboardMobilizerFrames.size(); ++i)
            mutvars.outboardMobilizerFrames[i] = Transform(); // M frame on B
        for (int i = 0; i < (int)mutvars.inboardMobilizerFrames.size(); ++i)
            mutvars.inboardMobilizerFrames[i] = Transform();  // F frame on P
        for (int i = 0; i < mutvars.particleMasses.size(); ++i)
            mutvars.particleMasses[i] = Real(1);
    }
};



// -----------------------------------------------------------------------------
//                                 TIME VARS
// -----------------------------------------------------------------------------
class SBTimeVars {
public:
    // none
public:
    void allocate(const SBTopologyCache&) const {
    }
};


// -----------------------------------------------------------------------------
//                                POSITION VARS
// -----------------------------------------------------------------------------
class SBPositionVars {
public:
    // none here -- q is supplied directly by the State
public:
    void allocate(const SBTopologyCache& tree) const {
    }
};


// -----------------------------------------------------------------------------
//                                VELOCITY VARS
// -----------------------------------------------------------------------------
class SBVelocityVars  {
public:
    // none here -- u is supplied directly by the State
public:
    void allocate(const SBTopologyCache&) const {
    }
};


// -----------------------------------------------------------------------------
//                                DYNAMICS VARS
// -----------------------------------------------------------------------------
class SBDynamicsVars {
public:
    // none here -- z is supplied directly by the State, but not
    //              used by the SimbodyMatterSubsystem anyway
public:
    void allocate(const SBTopologyCache&) const {    
    }
}; 



// -----------------------------------------------------------------------------
//                             ACCELERATION VARS
// -----------------------------------------------------------------------------
class SBAccelerationVars {
public:
    // none here
public:
    void allocate(const SBTopologyCache& topology) const {
    }
};


// These are here just so the AbstractValue's ValueHelper<> template
// will compile.
inline std::ostream& operator<<(std::ostream& o, const SBTopologyCache& c)
  { return o << "TODO: SBTopologyCache"; }
// SBModelCache output is implemented above

inline std::ostream& operator<<(std::ostream& o, const SBInstanceCache& c)
  { return o << "TODO: SBInstanceCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBTimeCache& c)
  { return o << "TODO: SBTimeCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBTreePositionCache& c)
  { return o << "TODO: SBTreePositionCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBConstrainedPositionCache& c)
  { return o << "TODO: SBConstrainedPositionCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBCompositeBodyInertiaCache& c)
  { return o << "TODO: SBCompositeBodyInertiaCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBArticulatedBodyInertiaCache& c)
  { return o << "TODO: SBArticulatedBodyInertiaCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBTreeVelocityCache& c)
  { return o << "TODO: SBTreeVelocityCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBConstrainedVelocityCache& c)
  { return o << "TODO: SBConstrainedVelocityCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBDynamicsCache& c)
  { return o << "TODO: SBDynamicsCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBTreeAccelerationCache& c)
  { return o << "TODO: SBTreeAccelerationCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBConstrainedAccelerationCache& c)
  { return o << "TODO: SBConstrainedAccelerationCache"; }

inline std::ostream& operator<<(std::ostream& o, const SBModelVars& c)
  { return o << "TODO: SBModelVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBInstanceVars& c)
  { return o << "TODO: SBInstanceVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBTimeVars& c)
  { return o << "TODO: SBTimeVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBPositionVars& c)
  { return o << "TODO: SBPositionVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBVelocityVars& c)
  { return o << "TODO: SBVelocityVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBDynamicsVars& c)
  { return o << "TODO: SBDynamicsVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBAccelerationVars& c)
  { return o << "TODO: SBAccelerationVars"; }




    /////////////////////
    // SB STATE DIGEST //
    /////////////////////

/*
 * Objects of this class are constructed for a particular State, and then used
 * briefly for a related series of computations. Depending on the stage to which the
 * State has been advanced, and the computations to be performed, some or all
 * of the pointers here will be set to refer to State and State cache data
 * for Simbody, of the types defined above.
 *
 * The idea is to do all the time consuming work of digging through the State
 * just once, then use the results repeatedly for computations which are typically
 * performed over all the nodes in the system. The low-level rigid body node computations
 * assume already-digested States.
 */
class SBStateDigest {
public:
    explicit SBStateDigest(const State& s) : state(s), modifiableState(0), stage(Stage::Empty) 
    {
    }
    SBStateDigest(const State& s, const SimbodyMatterSubsystemRep& matter, Stage g)
      : state(s), modifiableState(0), stage(Stage::Empty)
    {
        fillThroughStage(matter,g);
    }
    SBStateDigest(State& s, const SimbodyMatterSubsystemRep& matter, Stage g)
      : state(s), modifiableState(&s), stage(Stage::Empty)
    {
        fillThroughStage(matter,g);
    }

    // Stage g here is the stage we are about to compute. So we expect the referenced
    // State to have been realized to at least stage g-1.
    void fillThroughStage(const SimbodyMatterSubsystemRep& matter, Stage g);

    // The State is read only, for cache entries you have a choice.

    const State& getState() const {return state;}
    State&       updState() {
        assert(modifiableState);
        return *modifiableState;
    }
    Stage        getStage() const {return stage;}

    const SBModelVars& getModelVars() const {
        assert(stage >= Stage::Model);
        assert(mv);
        return *mv;
    }

    const SBInstanceVars& getInstanceVars() const {
        assert(stage >= Stage::Model);
        assert(iv);
        return *iv;
    }

    const SBTimeVars& getTimeVars() const {
        assert(stage >= Stage::Time);
        assert(tv);
        return *tv;
    }

    const Vector& getQ() const {
        assert(stage >= Stage::Position);
        assert(q);
        return *q;
    }

    const SBPositionVars& getPositionVars() const {
        assert(stage >= Stage::Position);
        assert(pv);
        return *pv;
    }

    const Vector& getU() const {
        assert(stage >= Stage::Velocity);
        assert(u);
        return *u;
    }

    const SBVelocityVars& getVelocityVars() const {
        assert(stage >= Stage::Velocity);
        assert(vv);
        return *vv;
    }
    const SBDynamicsVars& getDynamicsVars() const {
        assert(stage >= Stage::Dynamics);
        assert(dv);
        return *dv;
    }
    const SBAccelerationVars& getAccelerationVars() const {
        assert(stage >= Stage::Acceleration);
        assert(av);
        return *av;
    }

    // You can access the cache for update only at the stage being computed.
    // You can access the cache read-only for any stage already completed.
    // Either way you only need const access to the SBStateDigest object.

    // Model
    SBModelCache& updModelCache() const {
        assert(stage == Stage::Model);
        assert(mc);
        return *mc;
    }
    const SBModelCache& getModelCache() const {
        assert(stage > Stage::Model);
        assert(mc);
        return *mc;
    }

    // Instance
    SBInstanceCache& updInstanceCache() const {
        assert(stage == Stage::Instance);
        assert(ic);
        return *ic;
    }
    const SBInstanceCache& getInstanceCache() const {
        assert(stage > Stage::Instance);
        assert(ic);
        return *ic;
    }

    // Time
    SBTimeCache& updTimeCache() const {
        assert(stage == Stage::Time);
        assert(tc);
        return *tc;
    }
    const SBTimeCache& getTimeCache() const {
        assert(stage > Stage::Time);
        assert(tc);
        return *tc;
    }

    // Position
    Vector& updQErr() const {
        assert(stage == Stage::Position);
        assert(qErr);
        return *qErr;
    }
    const Vector& getQErr() const {
        assert(stage > Stage::Position);
        assert(qErr);
        return *qErr;
    }
    SBTreePositionCache& updTreePositionCache() const {
        assert(stage == Stage::Position);
        assert(tpc);
        return *tpc;
    }
    const SBTreePositionCache& getTreePositionCache() const {
        assert(stage > Stage::Position);
        assert(tpc);
        return *tpc;
    }
    SBConstrainedPositionCache& updConstrainedPositionCache() const {
        assert(stage == Stage::Position);
        assert(cpc);
        return *cpc;
    }
    const SBConstrainedPositionCache& getConstrainedPositionCache() const {
        assert(stage > Stage::Position);
        assert(cpc);
        return *cpc;
    }

    // Velocity
    Vector& updQDot() const {
        assert(stage == Stage::Velocity);
        assert(qdot);
        return *qdot;
    }
    const Vector& getQDot() const {
        assert(stage > Stage::Velocity);
        assert(qdot);
        return *qdot;
    }
    Vector& updUErr() const {
        assert(stage == Stage::Velocity);
        assert(uErr);
        return *uErr;
    }
    const Vector& getUErr() const {
        assert(stage > Stage::Velocity);
        assert(uErr);
        return *uErr;
    }
    SBTreeVelocityCache& updTreeVelocityCache() const {
        assert(stage == Stage::Velocity);
        assert(tvc);
        return *tvc;
    }
    const SBTreeVelocityCache& getTreeVelocityCache() const {
        assert(stage > Stage::Velocity);
        assert(tvc);
        return *tvc;
    }
    SBConstrainedVelocityCache& updConstrainedVelocityCache() const {
        assert(stage == Stage::Velocity);
        assert(cvc);
        return *cvc;
    }
    const SBConstrainedVelocityCache& getConstrainedVelocityCache() const {
        assert(stage > Stage::Velocity);
        assert(cvc);
        return *cvc;
    }

    // Dynamics
    SBDynamicsCache& updDynamicsCache() const {
        assert(stage == Stage::Dynamics);
        assert(dc);
        return *dc;
    }
    const SBDynamicsCache& getDynamicsCache() const {
        assert(stage > Stage::Dynamics);
        assert(dc);
        return *dc;
    }
    Vector& updUDot() const {
        assert(stage >= Stage::Dynamics);
        assert(udot);
        return *udot;
    }
    const Vector& getUDot() const {
        assert(stage > Stage::Dynamics);
        assert(udot);
        return *udot;
    }
    Vector& updQDotDot() const {
        assert(stage >= Stage::Dynamics);
        assert(qdotdot);
        return *qdotdot;
    }
    const Vector& getQDotDot() const {
        assert(stage > Stage::Dynamics);
        assert(qdotdot);
        return *qdotdot;
    }

    // Accelerations

    Vector& updUDotErr() const {
        assert(stage == Stage::Acceleration);
        assert(udotErr);
        return *udotErr;
    }
    const Vector& getUDotErr() const {
        assert(stage > Stage::Acceleration);
        assert(udotErr);
        return *udotErr;
    }
    SBTreeAccelerationCache& updTreeAccelerationCache() const {
        assert(stage == Stage::Acceleration);
        assert(tac);
        return *tac;
    }
    const SBTreeAccelerationCache& getTreeAccelerationCache() const {
        assert(stage > Stage::Acceleration);
        assert(tac);
        return *tac;
    }
    SBConstrainedAccelerationCache& updConstrainedAccelerationCache() const {
        assert(stage == Stage::Acceleration);
        assert(cac);
        return *cac;
    }
    const SBConstrainedAccelerationCache& getConstrainedAccelerationCache() const {
        assert(stage > Stage::Acceleration);
        assert(cac);
        return *cac;
    }

    void clear() {
        // state
        mv=0; iv=0; tv=0; 
        q=0; pv=0;
        u=0; vv=0; 
        dv=0; 
        av=0;

        // cache
        mc=0; ic=0; tc=0; 
        qErr=0; tpc=0; cpc=0;
        qdot=uErr=0; tvc=0; cvc=0; 
        dc=0; 
        udot=qdotdot=udotErr=0; tac=0; cac=0;
    }

private:
    const State&                    state;
    State*                          modifiableState;
    Stage                           stage; // the stage to be computed

    const SBModelVars*              mv;
    const SBInstanceVars*           iv;
    const SBTimeVars*               tv;

    const Vector*                   q;
    const SBPositionVars*           pv;

    const Vector*                   u;
    const SBVelocityVars*           vv;
    const SBDynamicsVars*           dv;
    const SBAccelerationVars*       av;

    SBModelCache*                   mc;
    SBInstanceCache*                ic;
    SBTimeCache*                    tc;

    Vector*                         qErr;
    SBTreePositionCache*            tpc;
    SBConstrainedPositionCache*     cpc;

    Vector*                         qdot;
    Vector*                         uErr;
    SBTreeVelocityCache*            tvc;
    SBConstrainedVelocityCache*     cvc;

    SBDynamicsCache*                dc;

    Vector*                         udot;
    Vector*                         qdotdot;
    Vector*                         udotErr;
    SBTreeAccelerationCache*        tac;
    SBConstrainedAccelerationCache* cac;
};

#endif // SimTK_SIMBODY_TREE_STATE_H_
