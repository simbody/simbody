#ifndef SIMTK_SIMBODY_TREE_STATE_H_
#define SIMTK_SIMBODY_TREE_STATE_H_

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
 * This file contains the classes which define the SimbodyTree State, that is, everything
 * that can be changed in a SimbodyTree after construction.
 *
 * State variables and computation results are organized into stages:
 *    UninitializedStage
 *    BuiltStage             Stored in the SimbodyTree object (construction)
 *   ---------------------------------------------------------
 *    ModeledStage           Stored in the State object
 *    ParametrizedStage
 *    TimedStage
 *    ConfiguredStage        (positions)
 *    MovingStage            (velocities)
 *    DynamicsStage          dynamic operators available
 *    ReactingStage          (response to forces in the state)
 *
 * Construction proceeds until all the bodies and constraints have been specified. After
 * that, realizeConstruction() is called. Construction-related 
 * calculations are performed leading to values which are stored in the SimbodyTree 
 * object, NOT in the State (e.g., total number of bodies). At the same time, an
 * initial state is built, with space allocated for the state variables that will
 * be needed by the next stage (ModeledStage),and these are assigned default values. 
 * Then the stage in the SimbodyTree and in the initial state is set to "Built".
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

#include "Simbody.h"

#include <cassert>
#include <vector>

class RigidBodyTree;
class RigidBodyNode;
template <int dof> class RigidBodyNodeSpec;

namespace simtk {

 // defined below

class SBModelingVars;
class SBParameterVars;
class SBTimeVars;
class SBConfigurationVars;
class SBMotionVars;
class SBDynamicsVars;
class SBReactionVars;

class SBModelingCache;
class SBParameterCache;
class SBTimeCache;
class SBConfigurationCache;
class SBMotionCache;
class SBDynamicsCache;
class SBReactionCache;

class SBStateRep;

    // TODO: constraint runtimes


class SBModelingCache {
public:
    // none yet
    // TODO: Modeling
    //   counts of various things resulting from modeling choices,
    //   constraint enabling, prescribed motion

public:
    void allocate(const RigidBodyTree&) {
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
    void allocate(const RigidBodyTree&, const SBStateRep&) {
        applyGravity = false;
    }
};

class SBTimeCache {
public:

    // none
public:
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {

    }
};

class SBConfigurationCache {
public:
    Vector sq, cq;  // nq  Sin&cos of angle q's in appropriate slots; otherwise garbage
    Vector qnorm;   // nq  Contains normalized quaternions in appropriate slots;
                    //       all else is garbage.
    Matrix_<Vec3> storageForHt; // 2 x ndof

    Array<TransformMat> bodyJointInParentJointFrame;  // nb (X_JbJ)

    Array<TransformMat> bodyConfigInParent;           // nb (X_PB)
    Array<TransformMat> bodyConfigInGround;           // nb (X_GB)
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
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.getNBodies();
        const int nDofs   = tree.getTotalDOF();     // this is the number of u's (nu)
        const int maxNQs  = tree.getTotalQAlloc();  // allocate the max # q's we'll ever need
        const int npc     = tree.getNDistanceConstraints(); // position constraints
        const int ndc     = tree.getNDistanceConstraints();

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
    Vector_<SpatialVec> bodyVelocityInParent;     // nb (joint velocity)
    Vector_<SpatialVec> bodyVelocityInGround;     // nb (sVel)

    Vector qdot;                                  // nq

    Vector velocityConstraintErrors;              // nvc

    // Distance constraint calculations. These are indexed by
    // *distance constraint* number, not *constraint* number.
    Vector_<Vec3> stationVel_G[2]; // vel of station relative to body origin, expr. in G
    Vector_<Vec3> vel_G[2];        // tip velocities relative to G, expr. in G
    Vector_<Vec3> relVel_G;        // spatial relative velocity tip2.velG-tip1.velG

public:
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.getNBodies();
        const int nDofs   = tree.getTotalDOF();     // this is the number of u's (nu)
        const int maxNQs  = tree.getTotalQAlloc();  // allocate the max # q's we'll ever need
        const int nvc     = tree.getNDistanceConstraints(); // velocity constraints
        const int ndc     = tree.getNDistanceConstraints();

        bodyVelocityInParent.resize(nBodies);       
        bodyVelocityInParent[0] = SpatialVec(Vec3(0),Vec3(0));

        bodyVelocityInGround.resize(nBodies);       
        bodyVelocityInGround[0] = SpatialVec(Vec3(0),Vec3(0));

        qdot.resize(maxNQs);

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

    Vector_<SpatialMat> psi;                      // nb
    Vector_<SpatialMat> tauBar;                   // nb
    Vector_<SpatialMat> Y;                        // nb

    Vector_<Real>       storageForD;              // sum(nu[j]^2)
    Vector_<Real>       storageForDI;             // sum(nu[j]^2)
    Matrix_<Vec3>       storageForG;              // 2 X ndof

public:
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.getNBodies();
        const int nDofs   = tree.getTotalDOF();     // this is the number of u's (nu)
        const int nSqDofs = tree.getTotalSqDOF();   // sum(ndof^2) for each joint
        const int maxNQs  = tree.getTotalQAlloc();  // allocate the max # q's we'll ever need
        const int nac     = tree.getNDistanceConstraints(); // acceleration constraints        
        
        articulatedBodyInertia.resize(nBodies); // TODO: ground initialization

        coriolisAcceleration.resize(nBodies);       
        coriolisAcceleration[0] = SpatialVec(Vec3(0),Vec3(0));

        gyroscopicForces.resize(nBodies);           
        gyroscopicForces[0] = SpatialVec(Vec3(0),Vec3(0));

        centrifugalForces.resize(nBodies);           
        centrifugalForces[0] = SpatialVec(Vec3(0),Vec3(0));

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
    Vector_<SpatialVec> bodyAccelerationInGround; // nb (sAcc)
    Vector udot;                                  // nu
    Vector lambda;                                // nac
    Vector netHingeForces;                        // nu (T-(~Am+R(F+C))
    Vector qdotdot;                               // nq

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
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.getNBodies();
        const int nDofs   = tree.getTotalDOF();     // this is the number of u's (nu)
        const int nSqDofs = tree.getTotalSqDOF();   // sum(ndof^2) for each joint
        const int maxNQs  = tree.getTotalQAlloc();  // allocate the max # q's we'll ever need
        const int nac     = tree.getNDistanceConstraints(); // acceleration constraints 
        const int ndc     = tree.getNDistanceConstraints();

        bodyAccelerationInGround.resize(nBodies);   
        bodyAccelerationInGround[0] = SpatialVec(Vec3(0),Vec3(0));;

        udot.resize(nDofs);
        lambda.resize(nac);
        netHingeForces.resize(nDofs);
        qdotdot.resize(maxNQs);

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
    void allocate(const RigidBodyTree& tree) const {
        assert(tree.isBuilt());

        SBModelingVars& mutvars = *const_cast<SBModelingVars*>(this);
        mutvars.useEulerAngles = false;
        mutvars.prescribed.resize(tree.getNBodies()); 
        mutvars.enabled.resize(tree.getNConstraints());
    }

};

class SBParameterVars {
public:
    Vec3 gravity;

    // TODO: body masses, etc.
public:

    // We can access the tree or state variable & cache up to Modeling stage.
    void allocate(const RigidBodyTree&, const SBStateRep&) const {
        SBParameterVars& mutvars = *const_cast<SBParameterVars*>(this);
        mutvars.gravity.setToNaN();
    }
};

class SBTimeVars {
public:
    // none
public:
    void allocate(const RigidBodyTree&, const SBStateRep&) const {
    }
};

class SBConfigurationVars {
public:
    Vector q; // [nq]
public:
    // We can access the tree or state variable & cache up to ModeledStage.
    void allocate(const RigidBodyTree& tree, const SBStateRep&) const {
        SBConfigurationVars& mutvars = *const_cast<SBConfigurationVars*>(this);
        mutvars.q.resize(tree.getTotalQAlloc()); 
        mutvars.q.setToNaN();
    }
};

class SBMotionVars  {
public:
    Vector u; // [ndof]  (== nu)
public:
    // We can access the tree or state variable & cache up to ModeledStage.
    void allocate(const RigidBodyTree& tree, const SBStateRep&) const {
        SBMotionVars& mutvars = *const_cast<SBMotionVars*>(this);
        mutvars.u.resize(tree.getTotalDOF()); 
        mutvars.u.setToNaN();
    }
};

class SBDynamicsVars {
public:
    // none
public:
    void allocate(const RigidBodyTree&, const SBStateRep&) const {
    }
};

class SBReactionVars {
public:
    Vector_<SpatialVec> appliedBodyForces;  // nb
    Vector              appliedJointForces; // nu
    Vector              prescribedUdot;     // nu
public:

    // We can access the tree or state variable & cache up to ModeledStage.
    void allocate(const RigidBodyTree& tree, const SBStateRep&) const {    
        SBReactionVars& mutvars = *const_cast<SBReactionVars*>(this);

        mutvars.appliedBodyForces.resize(tree.getNBodies());  
        mutvars.appliedBodyForces.setToNaN();

        mutvars.appliedJointForces.resize(tree.getTotalDOF());   
        mutvars.appliedJointForces.setToNaN();

        mutvars.prescribedUdot.resize(tree.getTotalDOF());       
        mutvars.prescribedUdot.setToNaN();
    }
}; 

class SBStateRep {
public:
    SBStateRep() : stage(UninitializedStage), handle(0) { }
    ~SBStateRep() {
        // in case we want to catch this in the debugger
     }

    // Default copy and assignment -- watch out for handle; the pointer
    // is copied by default and should be changed afterwards.

    void setMyHandle(SBState& s) {
        handle = &s;
    }

    SBStage getStage(const RigidBodyTree& tree) const {
        return tree.isBuilt() ? stage : UninitializedStage;
    }

    void setStage(const RigidBodyTree& tree, SBStage g) const {
        if (!tree.isBuilt() && g == UninitializedStage) {
            stage = UninitializedStage; // this is OK
            return;
        }
        assert(tree.isBuilt());
        assert(g >= BuiltStage);    // can't use this to return to uninitialized
        assert(stage >= g-1);       // can only advance one stage
        stage = g;                  // backing up any amount is OK
    }

    void initializeModelingVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= BuiltStage);
        if (getStage(tree) > BuiltStage)
            stage = BuiltStage; // back up if necessary

        SBStateRep& mutableState = *const_cast<SBStateRep*>(this);
        modelVars.allocate(tree); 
        tree.setDefaultModelingValues(*this, mutableState.modelVars); 
    }

    // Initialize the rest of the variables.
    void initializeAllVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        if (getStage(tree) > ModeledStage)
            stage = ModeledStage; // back up if necessary

        SBStateRep& mutableState = *const_cast<SBStateRep*>(this);

        // Don't initialize modeling vars!
        paramVars.allocate(tree, *this);  
        tree.setDefaultParameterValues(*this, mutableState.paramVars); 
        timeVars.allocate(tree, *this);   
        tree.setDefaultTimeValues(*this, mutableState.timeVars); 
        configVars.allocate(tree, *this);   
        tree.setDefaultConfigurationValues(*this, mutableState.configVars); 
        motionVars.allocate (tree, *this);   
        tree.setDefaultMotionValues(*this, mutableState.motionVars); 
        dynamicsVars.allocate(tree, *this);   
        tree.setDefaultDynamicsValues(*this, mutableState.dynamicsVars); 
        reactionVars.allocate(tree, *this);   
        tree.setDefaultReactionValues(*this, mutableState.reactionVars); 
    }

    // We're about to realize stage g and we want to make sure the cache entries we'll
    // need have been allocated and initialized properly. If the current stage is
    // g or greater then there is nothing to do since the cache must have been
    // allocated earlier. We expect to be able to access stage g-1 to figure out
    // the appropriate sizes and defaults here, so the stage must already be there.

    void allocateCacheIfNeeded(const RigidBodyTree& tree, SBStage g) const {
        assert(g > BuiltStage); // "built" vars & cache are in the RigidBodyTree, not the state
        assert(getStage(tree) >= g-1);
        if (getStage(tree) >= g)
            return;

        // These are uninitialized, i.e. garbage (NaN during Debug).
        switch (g) {
        case ModeledStage:      modelCache.allocate   (tree); break;
        case ParametrizedStage: paramCache.allocate   (tree, *this); break;
        case TimedStage:        timeCache.allocate    (tree, *this); break;
        case ConfiguredStage:   configCache.allocate  (tree, *this); break;
        case MovingStage:       motionCache.allocate  (tree, *this); break;
        case DynamicsStage:     dynamicsCache.allocate(tree, *this); break;
        case ReactingStage:     reactionCache.allocate(tree, *this); break;
        default: assert(false);
        }
    }

        // STAGE-CHECKED VARIABLE ACCESS

    // You can look at state variables as soon as they have been allocated.

    const SBModelingVars& getModelingVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= BuiltStage);
        return modelVars;
    }
    SBModelingVars& updModelingVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= BuiltStage);
        stage = BuiltStage; // backup if necessary
        return modelVars;
    }
    const SBParameterVars& getParameterVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        return paramVars;
    }
    SBParameterVars& updParameterVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ModeledStage);
        if (getStage(tree) >= ParametrizedStage)
            stage = SBStage(ParametrizedStage-1); // backup if necessary
        return paramVars;
    }
    const SBTimeVars& getTimeVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        return timeVars;
    }
    SBTimeVars& updTimeVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ModeledStage);
        if (getStage(tree) >= TimedStage)
            stage = SBStage(TimedStage-1); // backup if necessary
        return timeVars;
    }
    const SBConfigurationVars& getConfigurationVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        return configVars;
    }
    SBConfigurationVars& updConfigurationVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ModeledStage);
        if (getStage(tree) >= ConfiguredStage)
            stage = SBStage(ConfiguredStage-1); // backup if necessary
        return configVars;
    }
    const SBMotionVars& getMotionVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        return motionVars;
    }
    SBMotionVars& updMotionVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ModeledStage);
        if (getStage(tree) >= MovingStage)
            stage = SBStage(MovingStage-1); // backup if necessary
        return motionVars;
    }
    const SBDynamicsVars& getDynamicsVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        return dynamicsVars;
    }
    SBDynamicsVars& updDynamicsVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ModeledStage);
        if (getStage(tree) >= DynamicsStage)
            stage = SBStage(DynamicsStage-1); // backup if necessary
        return dynamicsVars;
    }
    const SBReactionVars& getReactionVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        return reactionVars;
    }
    SBReactionVars& updReactionVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ModeledStage);
        if (getStage(tree) >= ReactingStage)
            stage = SBStage(ReactingStage-1); // backup if necessary
        return reactionVars;
    }
        // STAGE-CHECKED CACHE ACCESS

    // If the "inProgress" flag is set it means that we are still calculating the
    // relevant stage and we'll trust that the caller knows that the cache entry
    // in question has already been calculated, meaning we'll only check that the
    // *previous* stage is finished.

    const SBModelingCache& getModelingCache(const RigidBodyTree& tree, 
                                            bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? ModeledStage-1 : ModeledStage));
        return modelCache;
    }
    SBModelingCache& updModelingCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage-1);
        return modelCache;
    }
    const SBParameterCache& getParameterCache(const RigidBodyTree& tree, 
                                              bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? ParametrizedStage-1 : ParametrizedStage));
        return paramCache;
    }
    SBParameterCache& updParameterCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ParametrizedStage-1);
        return paramCache;
    }
    const SBTimeCache& getTimeCache(const RigidBodyTree& tree, 
                                    bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? TimedStage-1 : TimedStage));
        return timeCache;
    }
    SBTimeCache& updTimeCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= TimedStage-1);
        return timeCache;
    }
    const SBConfigurationCache& getConfigurationCache(const RigidBodyTree& tree, 
                                                      bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? ConfiguredStage-1 : ConfiguredStage));
        return configCache;
    }
    SBConfigurationCache& updConfigurationCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ConfiguredStage-1);
        return configCache;
    }
    const SBMotionCache& getMotionCache(const RigidBodyTree& tree, 
                                        bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? MovingStage-1 : MovingStage));
        return motionCache;
    }
    SBMotionCache& updMotionCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= MovingStage-1);
        return motionCache;
    }
    const SBDynamicsCache& getDynamicsCache(const RigidBodyTree& tree, 
                                            bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? DynamicsStage-1 : DynamicsStage));
        return dynamicsCache;
    }
    SBDynamicsCache& updDynamicsCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= DynamicsStage-1);
        return dynamicsCache;
    }
    const SBReactionCache& getReactionCache(const RigidBodyTree& tree, 
                                            bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? ReactingStage-1 : ReactingStage));
        return reactionCache;
    }
    SBReactionCache& updReactionCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ReactingStage-1);
        return reactionCache;
    }
private:
    SBModelingVars      modelVars; // allocation sizes & default values only are mutable
    SBParameterVars     paramVars;
    SBTimeVars          timeVars;
    SBConfigurationVars configVars;
    SBMotionVars        motionVars;
    SBDynamicsVars      dynamicsVars;
    SBReactionVars      reactionVars;
    
    mutable SBStage              stage; // last stage completed

    mutable SBModelingCache      modelCache;
    mutable SBParameterCache     paramCache;
    mutable SBTimeCache          timeCache;
    mutable SBConfigurationCache configCache;
    mutable SBMotionCache        motionCache;
    mutable SBDynamicsCache      dynamicsCache;
    mutable SBReactionCache      reactionCache;

    SBState* handle;
    friend class RigidBodyTree;
    friend class RigidBodyNode;
    friend class RBDistanceConstraint;
    template <int dof> friend class RigidBodyNodeSpec;
};

}; // namespace simtk

#endif // SIMTK_SIMBODY_TREE_STATE_H_
