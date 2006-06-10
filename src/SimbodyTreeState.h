#ifndef SimTK_SIMBODY_TREE_STATE_H_
#define SimTK_SIMBODY_TREE_STATE_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * This file contains the classes which define the SimbodySubsystem State, that is, everything
 * that can be changed in a SimbodySubsystem after construction.
 *
 * State variables and computation results are organized into stages:
 *    Stage::Allocated
 *    Stage::Built             Stored in the SimbodySubsystem object (construction)
 *   ---------------------------------------------------------
 *    Stage::Modeled           Stored in the State object
 *    Stage::Parametrized
 *    Stage::Timed
 *    Stage::Configured        (positions)
 *    Stage::Moving            (velocities)
 *    Stage::Dynamics          dynamic operators available
 *    Stage::Reacting          (response to forces in the state)
 *
 * Construction proceeds until all the bodies and constraints have been specified. After
 * that, realizeConstruction() is called. Construction-related 
 * calculations are performed leading to values which are stored in the SimbodySubsystem 
 * object, NOT in the State (e.g., total number of bodies). At the same time, an
 * initial state is built, with space allocated for the state variables that will
 * be needed by the next stage (Stage::Modeled),and these are assigned default values. 
 * Then the stage in the SimbodySubsystem and in the initial state is set to "Built".
 *
 * After that, Modeling values can be set in the State. When that's done we call
 * realizeModeling(), which evaluates the Modeling states putting the values into
 * state cache entries allocated for the purpose. Then all remaining state variables
 * are allocated, and set to their default values. All defaults must be computable
 * knowing only the Modeling values. Then the stage is advanced to Modeled.
 *
 * This continues through all the stages, with realizeWhatever() expecting to receive
 * a state evaluated to stage Whatever-1 equipped with values for stage Whatever so that
 * it can calculate results and put them in the cache (which is allocated if necessary),
 * and then advance to stage Whatevered. 
 */

#include "simbody/internal/common.h"

#include <cassert>
#include <vector>

class RigidBodyTree;
class RigidBodyNode;
template <int dof> class RigidBodyNodeSpec;

namespace SimTK {

 // defined below


class SBModelingVars;
class SBParameterVars;
class SBTimeVars;
class SBConfigurationVars;
class SBMotionVars;
class SBDynamicsVars;
class SBReactionVars;

class SBConstructionCache;
class SBModelingCache;
class SBParameterCache;
class SBTimeCache;
class SBConfigurationCache;
class SBMotionCache;
class SBDynamicsCache;
class SBReactionCache;

class State;

// An object of this type is stored in the SimbodySubsystem after construction,
// then copied into a slot in the State on realizeConstruction(). It should contain
// enough information to size the other stages, and can also contain whatever
// arbitrary data you would like to have in a State to verify that it is a match
// for the Subsystem.
class SBConstructionCache {
public:
    SBConstructionCache() {
        nBodies = nConstraints = nDOFs = maxNQs = sumSqDOFs =
            nDistanceConstraints = modelingVarsIndex = modelingCacheIndex = -1;
        valid = false;
    }

    int nBodies;
    int nConstraints;

    int nDOFs;
    int maxNQs;
    int sumSqDOFs;

    int nDistanceConstraints;

    int modelingVarsIndex;
    int modelingCacheIndex;

    bool valid;
};

class SBModelingCache {
public:
    // TODO: Modeling
    //   counts of various things resulting from modeling choices,
    //   constraint enabling, prescribed motion

    SBModelingCache() {
        parameterVarsIndex = parameterCacheIndex
        = timeVarsIndex = timeCacheIndex
        = qIndex = qVarsIndex = qCacheIndex
        = uIndex = uVarsIndex = uCacheIndex
        = dynamicsVarsIndex = dynamicsCacheIndex
        = reactionVarsIndex = reactionCacheIndex
        = -1;
    }

    int parameterVarsIndex, parameterCacheIndex;
    int timeVarsIndex, timeCacheIndex;
    int qIndex; // maxNQs of these 
    int qVarsIndex, qCacheIndex;
    int uIndex; // nDOFs of these 
    int uVarsIndex, uCacheIndex;
    int dynamicsVarsIndex, dynamicsCacheIndex;
    int reactionVarsIndex, reactionCacheIndex;

public:
    void allocate(const SBConstructionCache&) {
    }
};

class SBParameterCache {
public:
    bool applyGravity;  // set after we see what value is in gravity parameter

    // TODO: Parameters
    //   body mass props; particle masses
    //   X_BJ, X_PJi transforms
    //   distance constraint distances & station positions

public:
    void allocate(const SBConstructionCache&) {
        applyGravity = false;
    }
};

class SBTimeCache {
public:

    // none
public:
    void allocate(const SBConstructionCache& tree) {

    }
};

class SBConfigurationCache {
public:
    Vector sq, cq;  // nq  Sin&cos of angle q's in appropriate slots; otherwise garbage
    Vector qnorm;   // nq  Contains normalized quaternions in appropriate slots;
                    //       all else is garbage.
    Matrix_<Vec3> storageForHt; // 2 x ndof

    Array<Transform>    bodyJointInParentJointFrame;  // nb (X_JbJ)

    Array<Transform>    bodyConfigInParent;           // nb (X_PB)
    Array<Transform>    bodyConfigInGround;           // nb (X_GB)
    Array<PhiMatrix>    bodyToParentShift;            // nb (phi)
    Array<InertiaMat>   bodyInertiaInGround;          // nb (I_OB_G)
    Vector_<SpatialMat> bodySpatialInertia;           // nb (Mk)
    Vector_<Vec3>       bodyCOMInGround;              // nb (COM_G)
    Vector_<Vec3>       bodyCOMStationInGround;       // nb (COMstation_G)

    Vector              positionConstraintErrors;     // npc

    // Distance constraint calculations. These are indexed by
    // *distance constraint* number, not *constraint* number.
    Vector_<Vec3> station_G[2];   // vec from body origin OB to station, expr. in G
    Vector_<Vec3> pos_G[2];       // vec from ground origin to station, expr. in G

    Vector_<Vec3> fromTip1ToTip2_G; // tip2.pos-tip1.pos
    Vector_<Vec3> unitDirection_G;  // fromTip1ToTip2/|fromTip1ToTip2|


public:
    void allocate(const SBConstructionCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int npc     = tree.nDistanceConstraints; // position constraints
        const int ndc     = tree.nDistanceConstraints;

        // These contain uninitialized junk. Body-indexed entries get their
        // ground elements set appropriately now and forever.
        sq.resize(maxNQs);
        cq.resize(maxNQs);
        qnorm.resize(maxNQs);
        storageForHt.resize(2,nDofs);

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

        positionConstraintErrors.resize(npc);

        station_G[0].resize(ndc); station_G[1].resize(ndc);
        pos_G[0].resize(ndc); pos_G[1].resize(ndc);

        fromTip1ToTip2_G.resize(ndc);
        unitDirection_G.resize(ndc);
    }
};

class SBMotionCache {
public:
    // qdot is supplied directly by the State
    Vector_<SpatialVec> bodyVelocityInParent;      // nb (joint velocity)
    Vector_<SpatialVec> bodyVelocityInGround;      // nb (sVel)
    Vector_<SpatialVec> mobilizerRelativeVelocity; // nb (V_JbJ)

    Vector velocityConstraintErrors;              // nvc

    // Distance constraint calculations. These are indexed by
    // *distance constraint* number, not *constraint* number.
    Vector_<Vec3> stationVel_G[2]; // vel of station relative to body origin, expr. in G
    Vector_<Vec3> vel_G[2];        // tip velocities relative to G, expr. in G
    Vector_<Vec3> relVel_G;        // spatial relative velocity tip2.velG-tip1.velG

public:
    void allocate(const SBConstructionCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int nvc     = tree.nDistanceConstraints; // velocity constraints
        const int ndc     = tree.nDistanceConstraints;

        bodyVelocityInParent.resize(nBodies);       
        bodyVelocityInParent[0] = SpatialVec(Vec3(0),Vec3(0));

        bodyVelocityInGround.resize(nBodies);       
        bodyVelocityInGround[0] = SpatialVec(Vec3(0),Vec3(0));

        mobilizerRelativeVelocity.resize(nBodies);       
        mobilizerRelativeVelocity[0] = SpatialVec(Vec3(0),Vec3(0));

        velocityConstraintErrors.resize(nvc);
        stationVel_G[0].resize(ndc); stationVel_G[1].resize(ndc);
        vel_G[0].resize(ndc); vel_G[1].resize(ndc);
        relVel_G.resize(ndc);
    }
};

class SBDynamicsCache {
public:
    // Dynamics
    Vector_<SpatialMat> articulatedBodyInertia;   // nb (P)

    Vector_<SpatialVec> coriolisAcceleration;     // nb (a)
    Vector_<SpatialVec> gyroscopicForces;         // nb (b)
    Vector_<SpatialVec> centrifugalForces;        // nb (P*a+b)

    Vector_<SpatialVec> appliedRigidBodyForces; // nb
    Vector_<Vec3>       appliedParticleForces;  // TODO
    Vector              appliedMobilityForces;  // nu
    Vector              prescribedUdot;     // nu

    Vector_<SpatialMat> psi;                      // nb
    Vector_<SpatialMat> tauBar;                   // nb
    Vector_<SpatialMat> Y;                        // nb

    Vector_<Real>       storageForD;              // sum(nu[j]^2)
    Vector_<Real>       storageForDI;             // sum(nu[j]^2)
    Matrix_<Vec3>       storageForG;              // 2 X ndof

public:
    void allocate(const SBConstructionCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs;   // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int nac     = tree.nDistanceConstraints; // acceleration constraints        
        
        articulatedBodyInertia.resize(nBodies); // TODO: ground initialization

        coriolisAcceleration.resize(nBodies);       
        coriolisAcceleration[0] = SpatialVec(Vec3(0),Vec3(0));

        gyroscopicForces.resize(nBodies);           
        gyroscopicForces[0] = SpatialVec(Vec3(0),Vec3(0));

        centrifugalForces.resize(nBodies);           
        centrifugalForces[0] = SpatialVec(Vec3(0),Vec3(0));

        appliedRigidBodyForces.resize(nBodies);
        appliedRigidBodyForces[0] = SpatialVec(Vec3(0),Vec3(0));
        appliedParticleForces.resize(0); // TODO
        appliedMobilityForces.resize(nDofs);
        prescribedUdot.resize(nDofs); // TODO

        psi.resize(nBodies); // TODO: ground initialization
        tauBar.resize(nBodies); // TODO: ground initialization

        Y.resize(nBodies);
        Y[0] = SpatialMat(Mat33(0));

        storageForD.resize(nSqDofs);
        storageForDI.resize(nSqDofs);
        storageForG.resize(2,nDofs);
    }
};


class SBReactionCache {
public:
    // udot, qdotdot are provided directly by the State
    Vector_<SpatialVec> bodyAccelerationInGround; // nb (sAcc)
    Vector              lambda;                   // nac
    Vector              netHingeForces;           // nu (T-(~Am+R(F+C))

    Vector              nu;
    Vector              epsilon;
    Vector_<SpatialVec> z;                        // nb
    Vector_<SpatialVec> Gepsilon;                 // nb

    Vector accelerationConstraintErrors;          // nac

    // Distance constraint calculations. These are indexed by
    // *distance constraint* number, not *constraint* number.
    Vector_<Vec3> acc_G[2];   // acc of tip relative to ground, expr. in G
    Vector_<Vec3> force_G[2]; // the constraint forces applied to each point

public:
    void allocate(const SBConstructionCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs;   // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int nac     = tree.nDistanceConstraints; // acceleration constraints 
        const int ndc     = tree.nDistanceConstraints;

        bodyAccelerationInGround.resize(nBodies);   
        bodyAccelerationInGround[0] = SpatialVec(Vec3(0),Vec3(0));;

        lambda.resize(nac);
        netHingeForces.resize(nDofs);

        nu.resize(nDofs);
        epsilon.resize(nDofs);
        z.resize(nBodies);

        Gepsilon.resize(nBodies); // TODO: ground initialization

        accelerationConstraintErrors.resize(nac);
        acc_G[0].resize(ndc); acc_G[1].resize(ndc);
        force_G[0].resize(ndc); force_G[1].resize(ndc);
    }
};

/** 
 * Generalized state variable collection for a SimbodyTree. 
 * Variables are divided into Stages, according to when their values
 * are needed during a calculation. The Stages that matter to the
 * MultibodyTree are:
 *       (Construction: not part of the state)
 *     Modeling:        choice of coordinates, knowns & unknowns, methods, etc.
 *     Parametrization: setting of physical parameters, e.g. mass
 *       (Time: not relevant to MultibodyTree)
 *     Configuration:   position and orientation values (2nd order continuous)
 *     Motion:          rates
 *     Dynamics;        dynamic quantities & operators available
 *     Reaction:        response to forces & prescribed accelerations in State 
 *
 */

class SBModelingVars {
public:
    bool        useEulerAngles;
    Array<bool> prescribed;           // nb  (# bodies & joints, 0 always true)
    Array<bool> enabled;              // nac (# acceleration constraints)
public:

    // We have to allocate these without looking at any other
    // state variable or cache entries. We can only depend on the tree
    // itself for information.
    void allocate(const SBConstructionCache& tree) const {
        SBModelingVars& mutvars = *const_cast<SBModelingVars*>(this);
        mutvars.useEulerAngles = false;
        mutvars.prescribed.resize(tree.nBodies); 
        mutvars.enabled.resize(tree.nConstraints);
    }

};

class SBParameterVars {
public:
    Vec3 gravity;

    // TODO: body masses, etc.
public:

    // We can access the tree or state variable & cache up to Modeling stage.
    void allocate(const SBConstructionCache&) const {
        SBParameterVars& mutvars = *const_cast<SBParameterVars*>(this);
        mutvars.gravity.setToNaN();
    }

    // Call this from Modeling stage to put some reasonable
    // defaults here.
    void initialize() {
        gravity = Vec3(0);
    }

};

class SBTimeVars {
public:
    // none
public:
    void allocate(const SBConstructionCache&) const {
    }
};

class SBConfigurationVars {
public:
    // none -- q is supplied directly by the State
public:
    void allocate(const SBConstructionCache& tree) const {
    }
};

class SBMotionVars  {
public:
    // none -- u is supplied directly by the State
public:
    void allocate(const SBConstructionCache&) const {
    }
};

class SBDynamicsVars {
public:
    // none
public:
    void allocate(const SBConstructionCache&) const {    
    }
}; 


class SBReactionVars {
public:
    // none
public:
    void allocate(const SBConstructionCache&) const {
    }
};

// These are here just so the AbstractValue's ValueHelper<> template
// will compile.
inline std::ostream& operator<<(std::ostream& o, const SBConstructionCache& c)
  { return o << "TODO: SBConstructionCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBModelingCache& c)
  { return o << "TODO: SBModelingCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBParameterCache& c)
  { return o << "TODO: SBParameterCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBTimeCache& c)
  { return o << "TODO: SBTimeCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBConfigurationCache& c)
  { return o << "TODO: SBConfigurationCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBMotionCache& c)
  { return o << "TODO: SBMotionCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBDynamicsCache& c)
  { return o << "TODO: SBDynamicsCache"; }
inline std::ostream& operator<<(std::ostream& o, const SBReactionCache& c)
  { return o << "TODO: SBReactionCache"; }

inline std::ostream& operator<<(std::ostream& o, const SBModelingVars& c)
  { return o << "TODO: SBModelingVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBParameterVars& c)
  { return o << "TODO: SBParameterVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBTimeVars& c)
  { return o << "TODO: SBTimeVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBConfigurationVars& c)
  { return o << "TODO: SBConfigurationVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBMotionVars& c)
  { return o << "TODO: SBMotionVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBDynamicsVars& c)
  { return o << "TODO: SBDynamicsVars"; }
inline std::ostream& operator<<(std::ostream& o, const SBReactionVars& c)
  { return o << "TODO: SBReactionVars"; }

}; // namespace SimTK

#endif // SimTK_SIMBODY_TREE_STATE_H_
