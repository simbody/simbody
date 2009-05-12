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
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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
 * This file contains the classes which define the SimbodyMatterSubsystem State, that is, everything
 * that can be changed in a SimbodyMatterSubsystem after construction.
 *
 * State variables and computation results are organized into stages:
 *    Stage::Empty         virginal state just allocated
 *    Stage::Topology      Stored in the SimbodyMatterSubsystem object (construction)
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
 * Construction proceeds until all the bodies and constraints have been specified. After
 * that, realizeTopology() is called. Construction-related calculations are
 * performed leading to values which are stored in the SimbodyMatterSubsystem 
 * object, NOT in the State (e.g., total number of bodies). At the same time, an
 * initial state is built, with space allocated for the state variables that will
 * be needed by the next stage (Stage::Model),and these are assigned default values. 
 * Then the stage in the SimbodyMatterSubsystem and in the initial state is set to "Topology".
 *
 * After that, values for Model stage variables can be set in the State.
 * When that's done we call realizeModel(), which evaluates the Model states
 * putting the values into state cache entries allocated for the purpose. Then
 * all remaining state variables are allocated, and set to their default values.
 * All defaults must be computable knowing only the Model stage values.
 * Then the stage is advanced to Stage::Model.
 *
 * This continues through all the stages, with realizeWhatever() expecting to receive
 * a state evaluated to stage Whatever-1 equipped with values for stage Whatever so that
 * it can calculate results and put them in the cache (which is allocated if necessary),
 * and then advance to stage Whatever. 
 */

#include "simbody/internal/common.h"

#include <cassert>
#include <vector>
#include <iostream>
using std::cout; using std::endl;


using namespace SimTK;

class SimTK::State;

class SimbodyMatterSubsystemRep;
class RigidBodyNode;
template <int dof> class RigidBodyNodeSpec;

 // defined below

class SBModelVars;
class SBInstanceVars;
class SBTimeVars;
class SBPositionVars;
class SBVelocityVars;
class SBDynamicsVars;
class SBAccelerationVars;

class SBTopologyCache;
class SBModelCache;
class SBInstanceCache;
class SBTimeCache;
class SBPositionCache;
class SBVelocityCache;
class SBDynamicsCache;
class SBAccelerationCache;


// An object of this type is stored in the SimbodyMatterSubsystem after extended
// construction in complete, then copied into a slot in the State on
// realizeTopology(). It should contain enough information to size the other
// stages, and can also contain whatever arbitrary data you would like to have
// in a State to verify that it is a match for the Subsystem.
// 
// Note that this does not include all possible topological information in
// a SimbodyMatterSubsystem -- any subobjects are free to hold their own
// as long as they don't change it after realizeTopology().
class SBTopologyCache {
public:
    SBTopologyCache() {
        clear();
    }

    void clear() {
        nBodies = nParticles = nConstraints = nAncestorConstrainedBodies =
            nDOFs = maxNQs = sumSqDOFs = -1;
        valid = false;
    }

    // These are topological objects.
    int nBodies;
    int nParticles;
    int nConstraints;

    // This is the total number of Constrained Bodies appearing in all Constraints
    // where the Ancestor body is not Ground, excluding the Ancestor bodies 
    // themselves even if they are also Constrained Bodies (which is common).
    // This is used for sizing pool entries in various caches to hold precalculated
    // Ancestor-frame data about these bodies.
    int nAncestorConstrainedBodies;

    // TODO: these should be moved to Model stage.
    int nDOFs;
    int maxNQs;
    int sumSqDOFs;

    DiscreteVariableIndex modelingVarsIndex;
    CacheEntryIndex       modelingCacheIndex;

    bool valid;
};

class SBModelCache {
public:
    // Modeling
    //   counts of various things resulting from modeling choices,
    //   constraint enabling, prescribed motion

    SBModelCache() {
        clear();
    }

    // Restore this cache entry to its just-constructed condition.
    void clear() {
        totalNQInUse = totalNUInUse = totalNQuaternionsInUse = totalNAnglesInUse = -1;

        mobilizedBodyModelInfo.clear();
    }

    void allocate(const SBTopologyCache& tc) {
        mobilizedBodyModelInfo.resize(tc.nBodies);
    }

        // MOBILIZED BODY MODELING INFORMATION

    class PerMobilizedBodyModelInfo {
    public:
        PerMobilizedBodyModelInfo() : nQInUse(-1), nUInUse(-1), hasQuaternionInUse(false), nAnglesInUse(-1) { }

        int nQInUse, nUInUse;
        QIndex firstQIndex; // These count from 0 for this SimbodyMatterSubsystem
        UIndex firstUIndex;

        // In case there is a quaternion in use by this Mobilizer. The index here can be
        // used to find precalculated data associated with this quaternion, such as its current length.
        bool                hasQuaternionInUse;
        MobilizerQIndex     startOfQuaternion;   // 0..nQInUse-1: which local coordinate starts the quaternion if any?
        QuaternionPoolIndex quaternionPoolIndex; // assigned slot # for this MB's quat, -1 if none

        // In case there are any generalized coordinates which are angles. We require that the
        // angular coordinates be consecutive and just store the number of angles and the coordinate
        // index of the first one. The index can be used to find precalculated data associated with
        // angles, such as their sines and cosines.
        int             nAnglesInUse;
        MobilizerQIndex startOfAngles;  // 0..nQInUse-1: which local coordinate starts angles if any?
        AnglePoolIndex  anglePoolIndex; // start of assigned segment for this MB's angle information, if any
    };

    // Use these accessors so that you get type checking on the index types.
    int getNumMobilizedBodies() const {return (int)mobilizedBodyModelInfo.size();}
    PerMobilizedBodyModelInfo& updMobilizedBodyModelInfo(MobilizedBodyIndex mbx) {
        return mobilizedBodyModelInfo[mbx];
    }
    const PerMobilizedBodyModelInfo& getMobilizedBodyModelInfo(MobilizedBodyIndex mbx) const {
        return mobilizedBodyModelInfo[mbx];
    }

        // STATE ALLOCATION FOR THIS SUBSYSTEM

    // Note that a MatterSubsystem is only one of potentially many users of a System's State, so only
    // a subset of State variables and State Cache entries belong to it. Here we record the indices
    // we were given when we asked the State for some resources.

    QIndex qIndex;
    UIndex uIndex;
    DiscreteVariableIndex instanceVarsIndex, timeVarsIndex, qVarsIndex, uVarsIndex, 
                          dynamicsVarsIndex, accelerationVarsIndex;
    CacheEntryIndex       instanceCacheIndex, timeCacheIndex, qCacheIndex, uCacheIndex, 
                          dynamicsCacheIndex, accelerationCacheIndex;

    // These are sums over the per-MobilizedBody counts above.
    int totalNQInUse, totalNUInUse, totalNQuaternionsInUse, totalNAnglesInUse;

private:
    // Use accessor routines for these so that you get type checking on the index types.
    std::vector<PerMobilizedBodyModelInfo> mobilizedBodyModelInfo; // MobilizedBody 0 is Ground
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

class SBInstanceCache {
public:

        // CONSTRAINT MODELING INFORMATION

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
        // Better to access using accessor methods above so you'll get type checking on the index type.
        std::vector<PerConstrainedMobilizerInstanceInfo> constrainedMobilizerInstanceInfo;

        // The ConstrainedBodies and ConstrainedMobilizers are set at Topology stage, but the
        // particular generalized coordinates q and generalized speeds u which are involved
        // can't be determined until Model stage, since the associated mobilizers have Model
        // stage options which can affect the number and meanings of these variables.
        // These are sorted in order of their associated ConstrainedMobilizer, not necessarily
        // in order of QIndex or UIndex. Each value appears only once.
        std::vector<QIndex> constrainedQ;   // indexed by ConstrainedQIndex, maps to subsystem QIndex
        std::vector<UIndex> constrainedU;   // indexed by ConstrainedUIndex, maps to subsystem UIndex

        // Participating mobilities include ALL the mobilities which may be involved in any of this
        // Constraint's constraint equations, whether from being directly constrained or indirectly
        // as a result of their effects on ConstrainedBodies. These are sorted in order of increasing
        // QIndex and UIndex, and each QIndex or UIndex appears only once.
        std::vector<QIndex> participatingQ; // indexed by ParticipatingQIndex, maps to subsystem QIndex
        std::vector<UIndex> participatingU; // indexed by ParticipatingUIndex, maps to subsystem UIndex
    };

    // Instance variables are:
    //   body mass props; particle masses
    //   X_BM, X_PF mobilizer transforms
    //  
    // Calculations stored here derive from those states:
    //   total mass
    //   central inertia of each rigid body
    //   principal axes and corresponding principal moments of inertia of each rigid body
    //   reference configuration X_PB when q==0 (usually that means M==F), for each rigid body

    Real                   totalMass; // sum of all rigid body and particles masses
    std::vector<Inertia>   centralInertias;           // nb
    Vector_<Vec3>          principalMoments;          // nb
    std::vector<Rotation>  principalAxes;             // nb
    std::vector<Transform> referenceConfiguration;    // nb
    std::vector<PerConstraintInstanceInfo>    constraintInstanceInfo;

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
    void allocate(const SBTopologyCache& topology) {
        totalMass = SimTK::NaN;
        centralInertias.resize(topology.nBodies);           // I_CB
        principalMoments.resize(topology.nBodies);          // (Ixx,Iyy,Izz)
        principalAxes.resize(topology.nBodies);             // [axx ayy azz]
        referenceConfiguration.resize(topology.nBodies);    // X0_PB
        constraintInstanceInfo.resize(topology.nConstraints);
        firstQuaternionQErrSlot = qErrIndex = uErrIndex = udotErrIndex = -1;
    }
    int getNumConstraints() const {return (int)constraintInstanceInfo.size();}
    PerConstraintInstanceInfo& updConstraintInstanceInfo(ConstraintIndex cx) {
        return constraintInstanceInfo[cx];
    }
    const PerConstraintInstanceInfo& getConstraintInstanceInfo(ConstraintIndex cx) const {
        return constraintInstanceInfo[cx];
    }
};

class SBTimeCache {
public:

    // none
public:
    void allocate(const SBTopologyCache&) {
    }
};

class SBPositionCache {
public:
    Vector sq, cq;  // nq  Sin&cos of angle q's in appropriate slots; otherwise garbage
    Vector qnorm;   // nq  Contains normalized quaternions in appropriate slots;
                    //       all else is garbage.

    // CAUTION: our definition of the H matrix is transposed from those used
    // by Jain and by Schwieters. Jain would call these H* and Schwieters
    // would call them H^T, but we call them H.
    Matrix_<Vec3> storageForH_FM; // 2 x ndof (H_FM)
    Matrix_<Vec3> storageForH;    // 2 x ndof (H_PB_G)

    std::vector<Transform>    bodyJointInParentJointFrame;  // nb (X_FM)

    std::vector<Transform>    bodyConfigInParent;           // nb (X_PB)
    std::vector<Transform>    bodyConfigInGround;           // nb (X_GB)
    std::vector<PhiMatrix>    bodyToParentShift;            // nb (phi)
    std::vector<Inertia>      bodyInertiaInGround;          // nb (I_OB_G)
    Vector_<SpatialMat>       bodySpatialInertia;           // nb (Mk)
    Vector_<Vec3>             bodyCOMInGround;              // nb (p_G_CB)
    Vector_<Vec3>             bodyCOMStationInGround;       // nb (p_CB_G)



        // Ancestor Constrained Body Pool

    // For Constraints whose Ancestor body A is not Ground G, we assign pool entries
    // for each of their Constrained Bodies (call the total number 'nacb')
    // to store the above information but measured and expressed in the Ancestor frame
    // rather than Ground.
    std::vector<Transform> constrainedBodyConfigInAncestor;     // nacb (X_AB)

public:
    void allocate(const SBTopologyCache& tree) {
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

class SBVelocityCache {
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
    std::vector<SpatialVec> constrainedBodyVelocityInAncestor; // nacb (V_AB)

public:
    void allocate(const SBTopologyCache& tree) {
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

class SBDynamicsCache {
public:
    // Dynamics
    Vector_<SpatialMat> articulatedBodyInertia;   // nb (P)

    Vector_<SpatialVec> centrifugalForces;        // nb (P*a+b)
    Vector_<SpatialVec> totalCentrifugalForces;   // nb (P*A+b)

    Vector_<SpatialMat> psi;                      // nb
    Vector_<SpatialMat> tauBar;                   // nb
    Vector_<SpatialMat> Y;                        // nb

    Vector_<Real>       storageForD;              // sum(nu[j]^2)
    Vector_<Real>       storageForDI;             // sum(nu[j]^2)
    Matrix_<Vec3>       storageForG;              // 2 X ndof

public:
    void allocate(const SBTopologyCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs;   // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need     
        
        articulatedBodyInertia.resize(nBodies); // TODO: ground initialization

        centrifugalForces.resize(nBodies);           
        centrifugalForces[0] = SpatialVec(Vec3(0),Vec3(0));

        totalCentrifugalForces.resize(nBodies);           
        totalCentrifugalForces[0] = SpatialVec(Vec3(0),Vec3(0));

        psi.resize(nBodies); // TODO: ground initialization
        tauBar.resize(nBodies); // TODO: ground initialization

        Y.resize(nBodies);
        Y[0] = SpatialMat(Mat33(0));

        storageForD.resize(nSqDofs);
        storageForDI.resize(nSqDofs);
        storageForG.resize(2,nDofs);
    }
};


class SBAccelerationCache {
public:
    // udot, qdotdot cache space is provided directly by the State
    Vector_<SpatialVec> bodyAccelerationInGround; // nb (A_GB)

    Vector              netHingeForces;           // nu (T-(~Am+R(F+C))
    Vector              nu; // greek letter, not # of u's, but there are nu of these
    Vector              epsilon;                  // nu
    Vector_<SpatialVec> z;                        // nb
    Vector_<SpatialVec> Gepsilon;                 // nb

        // Ancestor Constrained Body Pool

    // For Constraints whose Ancestor body A is not Ground G, we assign pool entries
    // for each of their Constrained Bodies (call the total number 'nacb')
    // to store the above information but measured and expressed in the Ancestor frame
    // rather than Ground.
    std::vector<SpatialVec> constrainedBodyAccelerationInAncestor; // nacb (A_AB)

public:
    void allocate(const SBTopologyCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs; // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;    // allocate the max # q's we'll ever need
        const int nacb    = tree.nAncestorConstrainedBodies;

        bodyAccelerationInGround.resize(nBodies);   
        bodyAccelerationInGround[0] = SpatialVec(Vec3(0),Vec3(0));;

        netHingeForces.resize(nDofs);

        nu.resize(nDofs);
        epsilon.resize(nDofs);
        z.resize(nBodies);

        Gepsilon.resize(nBodies); // TODO: ground initialization

        constrainedBodyAccelerationInAncestor.resize(nacb);
    }
};

/** 
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

class SBModelVars {
public:
    bool              useEulerAngles;
    std::vector<bool> prescribed;           // nb (# bodies & mobilizers, [0] always true)
public:

    // We have to allocate these without looking at any other
    // state variable or cache entries. We can only depend on topological
    // information.
    void allocate(const SBTopologyCache& tree) {
        useEulerAngles = false;
        prescribed.resize(tree.nBodies, false); 
    }

};

class SBInstanceVars {
public:
    std::vector<MassProperties> bodyMassProperties;
    std::vector<Transform>      outboardMobilizerFrames;
    std::vector<Transform>      inboardMobilizerFrames;
    std::vector<bool>           disabled;             // nc (# constraints)
    Vector                      particleMasses;

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
            mutvars.outboardMobilizerFrames[i] = Transform(); // i.e., B frame
        for (int i = 0; i < (int)mutvars.inboardMobilizerFrames.size(); ++i)
            mutvars.inboardMobilizerFrames[i] = Transform(); // i.e., P frame
        for (int i = 0; i < mutvars.particleMasses.size(); ++i)
            mutvars.particleMasses[i] = Real(1);
    }
};

class SBTimeVars {
public:
    // none
public:
    void allocate(const SBTopologyCache&) const {
    }
};

class SBPositionVars {
public:
    // none here -- q is supplied directly by the State
public:
    void allocate(const SBTopologyCache& tree) const {
    }
};

class SBVelocityVars  {
public:
    // none here -- u is supplied directly by the State
public:
    void allocate(const SBTopologyCache&) const {
    }
};

class SBDynamicsVars {
public:
    // none here -- z is supplied directly by the State, but not
    //              used by the SimbodyMatterSubsystem anyway
public:
    void allocate(const SBTopologyCache&) const {    
    }
}; 


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
inline std::ostream& operator<<(std::ostream& o, const SBPositionCache& c)
  { return o << "TODO: SBPositionCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBVelocityCache& c)
  { return o << "TODO: SBVelocityCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBDynamicsCache& c)
  { return o << "TODO: SBDynamicsCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBAccelerationCache& c)
  { return o << "TODO: SBAccelerationCache"; }

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
        assert(stage >= Stage::Instance);
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
    SBPositionCache& updPositionCache() const {
        assert(stage == Stage::Position);
        assert(pc);
        return *pc;
    }
    const SBPositionCache& getPositionCache() const {
        assert(stage > Stage::Position);
        assert(pc);
        return *pc;
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
    SBVelocityCache& updVelocityCache() const {
        assert(stage == Stage::Velocity);
        assert(vc);
        return *vc;
    }
    const SBVelocityCache& getVelocityCache() const {
        assert(stage > Stage::Velocity);
        assert(vc);
        return *vc;
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

    // Accelerations
    Vector& updUDot() const {
        assert(stage == Stage::Acceleration);
        assert(udot);
        return *udot;
    }
    const Vector& getUDot() const {
        assert(stage > Stage::Acceleration);
        assert(udot);
        return *udot;
    }
    Vector& updQDotDot() const {
        assert(stage == Stage::Acceleration);
        assert(qdotdot);
        return *qdotdot;
    }
    const Vector& getQDotDot() const {
        assert(stage > Stage::Acceleration);
        assert(qdotdot);
        return *qdotdot;
    }
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
    SBAccelerationCache& updAccelerationCache() const {
        assert(stage == Stage::Acceleration);
        assert(ac);
        return *ac;
    }
    const SBAccelerationCache& getAccelerationCache() const {
        assert(stage > Stage::Acceleration);
        assert(ac);
        return *ac;
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
        qErr=0; pc=0;
        qdot=uErr=0; vc=0; 
        dc=0; 
        udot=qdotdot=udotErr=0; ac=0;
    }

private:
    const State& state;
    State*       modifiableState;
    Stage        stage; // the stage to be computed

    const SBModelVars*          mv;
    const SBInstanceVars*       iv;
    const SBTimeVars*           tv;

    const Vector*               q;
    const SBPositionVars*       pv;

    const Vector*               u;
    const SBVelocityVars*       vv;
    const SBDynamicsVars*       dv;
    const SBAccelerationVars*   av;

    SBModelCache*               mc;
    SBInstanceCache*            ic;
    SBTimeCache*                tc;

    Vector*                     qErr;
    SBPositionCache*            pc;

    Vector*                     qdot;
    Vector*                     uErr;
    SBVelocityCache*            vc;

    SBDynamicsCache*            dc;

    Vector*                     udot;
    Vector*                     qdotdot;
    Vector*                     udotErr;
    SBAccelerationCache*        ac;
};

#endif // SimTK_SIMBODY_TREE_STATE_H_
