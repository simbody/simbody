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
 *    ReactingStage          (dynamics)
 *
 * Construction proceeds until all the bodies and constraints have been specified. After
 * that, realizeConstruction() is called. Construction-related 
 * calculations are performed leading to values which are stored in the SimbodyTree 
 * object, NOT in the State (e.g., total number of bodies). At the same time, an
 * initial state is built, with space allocated for the state variables that will
 * be needed by the next stage (Modeled),and these are assigned default values. 
 * Then the stage in the SimbodyTree and in the initial state is set to "Built".
 * We do not consider any of this a State change, because there can be no effect on
 * the results (that is, calling realizeConstruction() again will
 * yield the same result.
 *
 * After that, Modeling values can be set in the State. When that's done we call
 * realizeModeling(), which evaluates the Modeling states putting the values into
 * state cache entries allocated for the purpose. If necessary, we allocate the
 * state variable needed for the next stage (Parameterized), and fill in the defaults.
 * The stage is advanced to Modeled.
 *
 * This continues through all the stages, with realizeWhatever() expecting to receive
 * a state evaluated to stage Whatever-1 equipped with values for stage Whatever so that
 * it can calculate results and put them in the cache, and then advance to stage Whatevered. 
 * Then we make sure that state variables are available for stage Whatever+1, allocate and
 * set to defaults if necessary, and return.
 */

#include "simbody/Simbody.h"

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
class SBDynamicVars;

class SBModelingResults;
class SBParameterResults;
class SBTimeResults;
class SBConfigurationResults;
class SBMotionResults;
class SBDynamicResults;

class SBStateRep;

/*
 * This is the cache.
 *
class SimbodyTreeResults {
public:
    SimbodyTreeResults() : stage(UninitializedStage) { } // everything has length 0
    // default copy, copy assign, destruct

    // This allocation routine should be called after realizeModeling(). Before that
    // we don't know enough about what to put here.

    // nDofs==nu==#joint forces
    // nSqDofs=sum(nu[j]^2) for all joints j
    void allocateCache(int nBodies, int nDofs, int nSqDofs, int maxNQs,
                       int npc, int nvc, int nac) // pos, vel, acc constraints
    {
        const SpatialVec zVec(Vec3(0));
        const SpatialMat zMat(Mat33(0));
        assert(stage >= ModeledStage);
        stage = ModeledStage; // roll back if necessary

        // These contain uninitialized junk. Body-indexed entries get their
        // ground elements set appropriately now and forever.
        sq.resize(maxNQs);
        cq.resize(maxNQs);
        qnorm.resize(maxNQs);
        storageForHt.resize(2,nDofs);

        bodyJointInParentJointFrame.resize(nBodies); bodyJointInParentJointFrame[0].setToZero();
        bodyConfigInParent.resize(nBodies);          bodyConfigInParent[0].setToZero();
        bodyConfigInGround.resize(nBodies);          bodyConfigInGround[0].setToZero();
        bodyToParentShift.resize(nBodies);           bodyToParentShift[0].setToZero();
        bodyInertiaInGround.resize(nBodies); // TODO
        bodySpatialInertia.resize(nBodies);  // TODO
        bodyCOMInGround.resize(nBodies);             bodyCOMInGround[0] = 0.;
        bodyCOMStationInGround.resize(nBodies);      bodyCOMStationInGround[0] = 0.;

        positionConstraintErrors.resize(npc);

        bodyVelocityInParent.resize(nBodies);       bodyVelocityInParent[0] = zVec;
        bodyVelocityInGround.resize(nBodies);       bodyVelocityInGround[0] = zVec;
        velocityConstraintErrors.resize(nvc);
        qdot.resize(maxNQs);

        articulatedBodyInertia.resize(nBodies); // TODO
        bodyAccelerationInGround.resize(nBodies);   bodyAccelerationInGround[0] = zVec;
        coriolisAcceleration.resize(nBodies);       coriolisAcceleration[0] = zVec;
        gyroscopicForces.resize(nBodies);           gyroscopicForces[0] = zVec;
        udot.resize(nDofs);
        lambda.resize(nac);
        accelerationConstraintErrors.resize(nac);
        netHingeForces.resize(nDofs);
        qdotdot.resize(maxNQs);
        psiT.resize(nBodies); // TODO
        tau.resize(nBodies); // TODO
        Y.resize(nBodies); // TODO
        storageForDI.resize(nSqDofs);
        storageForG.resize(2,nDofs);
        nu.resize(nDofs);
        epsilon.resize(nDofs);
        z.resize(nBodies);
        Gepsilon.resize(nBodies); // TODO
    }

    SBStage stage;

    // TODO: constraint runtimes

    // TODO: Modeling
    //   counts of various things resulting from modeling choices,
    //   constraint enabling, prescribed motion

    // TODO: Parameters
    //   body mass props; particle masses
    //   X_BJ, X_PJi transforms
    //   distance constraint distances & station positions

    // Configuration
    Vector sq, cq;  // nq  Sin&cos of angle q's in appropriate slots; otherwise garbage
    Vector qnorm;   // nq  Contains normalized quaternions in appropriate slots;
                    //       all else is garbage.
    Matrix_<Vec3> storageForHt; // 2 x ndof

    std::vector<TransformMat> bodyJointInParentJointFrame;  // nb (X_JbJ)

    std::vector<TransformMat> bodyConfigInParent;           // nb (X_PB)
    std::vector<TransformMat> bodyConfigInGround;           // nb (X_GB)
    std::vector<PhiMatrix>    bodyToParentShift;            // nb (phi)
    std::vector<InertiaMat>   bodyInertiaInGround;          // nb (I_OB_G)
    Vector_<SpatialMat>       bodySpatialInertia;           // nb (Mk)
    Vector_<Vec3>             bodyCOMInGround;              // nb (COM_G)
    Vector_<Vec3>             bodyCOMStationInGround;       // nb (COMstation_G)

    Vector                    positionConstraintErrors; // npc


    // Motion
    Vector_<SpatialVec> bodyVelocityInParent;     // nb (joint velocity)
    Vector_<SpatialVec> bodyVelocityInGround;     // nb (sVel)

    Vector velocityConstraintErrors;              // nvc
    Vector qdot;                                  // nq

    // Dynamics
    Vector_<SpatialMat> articulatedBodyInertia;   // nb (P)
    Vector_<SpatialVec> bodyAccelerationInGround; // nb (sAcc)
    Vector_<SpatialVec> coriolisAcceleration;     // nb (a)
    Vector_<SpatialVec> gyroscopicForces;         // nb (b)

    Vector udot;                                  // nu
    Vector lambda;                                // nac
    Vector accelerationConstraintErrors;          // nac
    Vector netHingeForces;                        // nu (T-(~Am+R(F+C))
    Vector qdotdot;                               // nq

    // dynamic temporaries
    Vector_<SpatialMat> psiT;                     // nb
    Vector_<SpatialMat> tau;                      // nb
    Vector_<SpatialMat> Y;                        // nb

    Vector_<Real>       storageForDI;             // sum(nu[j]^2)
    Matrix_<Vec3>       storageForG;              // 2 X ndof
    Vector              nu;
    Vector              epsilon;
    Vector_<SpatialVec> z;                        // nb
    Vector_<SpatialVec> Gepsilon;                 // nb
};
*/

    // TODO: constraint runtimes


class SBModelingCache {
    friend class RigidBodyNode;

    // none yet
    // TODO: Modeling
    //   counts of various things resulting from modeling choices,
    //   constraint enabling, prescribed motion

public:
    void allocate(const RigidBodyTree& tree) {
    }
};

class SBParameterCache {
public:

    // none yet
    // TODO: Parameters
    //   body mass props; particle masses
    //   X_BJ, X_PJi transforms
    //   distance constraint distances & station positions
public:
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {

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

    // Configuration
    Vector sq, cq;  // nq  Sin&cos of angle q's in appropriate slots; otherwise garbage
    Vector qnorm;   // nq  Contains normalized quaternions in appropriate slots;
                    //       all else is garbage.
    Matrix_<Vec3> storageForHt; // 2 x ndof

    std::vector<TransformMat> bodyJointInParentJointFrame;  // nb (X_JbJ)

    std::vector<TransformMat> bodyConfigInParent;           // nb (X_PB)
    std::vector<TransformMat> bodyConfigInGround;           // nb (X_GB)
    std::vector<PhiMatrix>    bodyToParentShift;            // nb (phi)
    std::vector<InertiaMat>   bodyInertiaInGround;          // nb (I_OB_G)
    Vector_<SpatialMat>       bodySpatialInertia;           // nb (Mk)
    Vector_<Vec3>             bodyCOMInGround;              // nb (COM_G)
    Vector_<Vec3>             bodyCOMStationInGround;       // nb (COMstation_G)

    Vector                    positionConstraintErrors;     // npc

public:
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.getNBodies();
        const int nDofs   = tree.getTotalDOF();     // this is the number of u's (nu)
        const int maxNQs  = tree.getTotalQAlloc();  // allocate the max # q's we'll ever need
        const int npc     = tree.getNConstraints(); // position constraints

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
    }
};

class SBMotionCache {
public:

    Vector_<SpatialVec> bodyVelocityInParent;     // nb (joint velocity)
    Vector_<SpatialVec> bodyVelocityInGround;     // nb (sVel)

    Vector velocityConstraintErrors;              // nvc
    Vector qdot;                                  // nq
public:
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.getNBodies();
        const int nDofs   = tree.getTotalDOF();     // this is the number of u's (nu)
        const int maxNQs  = tree.getTotalQAlloc();  // allocate the max # q's we'll ever need

        const int nvc     = tree.getNConstraints(); // velocity constraints

        bodyVelocityInParent.resize(nBodies);       
        bodyVelocityInParent[0] = SpatialVec(Vec3(0),Vec3(0));

        bodyVelocityInGround.resize(nBodies);       
        bodyVelocityInGround[0] = SpatialVec(Vec3(0),Vec3(0));

        velocityConstraintErrors.resize(nvc);
        qdot.resize(maxNQs);
    }
};

class SBDynamicCache {
public:

    // Dynamics
    Vector_<SpatialMat> articulatedBodyInertia;   // nb (P)
    Vector_<SpatialVec> bodyAccelerationInGround; // nb (sAcc)
    Vector_<SpatialVec> coriolisAcceleration;     // nb (a)
    Vector_<SpatialVec> gyroscopicForces;         // nb (b)

    Vector udot;                                  // nu
    Vector lambda;                                // nac
    Vector accelerationConstraintErrors;          // nac
    Vector netHingeForces;                        // nu (T-(~Am+R(F+C))
    Vector qdotdot;                               // nq

    // dynamic temporaries
    Vector_<SpatialMat> psiT;                     // nb
    Vector_<SpatialMat> tau;                      // nb
    Vector_<SpatialMat> Y;                        // nb

    Vector_<Real>       storageForDI;             // sum(nu[j]^2)
    Matrix_<Vec3>       storageForG;              // 2 X ndof
    Vector              nu;
    Vector              epsilon;
    Vector_<SpatialVec> z;                        // nb
    Vector_<SpatialVec> Gepsilon;                 // nb
public:
    void allocate(const RigidBodyTree& tree, const SBStateRep& s) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.getNBodies();
        const int nDofs   = tree.getTotalDOF();     // this is the number of u's (nu)
        const int nSqDofs = tree.getTotalSqDOF();   // sum(ndof^2) for each joint
        const int maxNQs  = tree.getTotalQAlloc();  // allocate the max # q's we'll ever need
        const int nac     = tree.getNConstraints(); // acceleration constraints        
        
        articulatedBodyInertia.resize(nBodies); // TODO: ground initialization

        bodyAccelerationInGround.resize(nBodies);   
        bodyAccelerationInGround[0] = SpatialVec(Vec3(0),Vec3(0));;

        coriolisAcceleration.resize(nBodies);       
        coriolisAcceleration[0] = SpatialVec(Vec3(0),Vec3(0));;

        gyroscopicForces.resize(nBodies);           
        gyroscopicForces[0] = SpatialVec(Vec3(0),Vec3(0));;

        udot.resize(nDofs);
        lambda.resize(nac);
        accelerationConstraintErrors.resize(nac);
        netHingeForces.resize(nDofs);
        qdotdot.resize(maxNQs);
        psiT.resize(nBodies); // TODO: ground initialization
        tau.resize(nBodies); // TODO: ground initialization
        Y.resize(nBodies); // TODO: ground initialization
        storageForDI.resize(nSqDofs);
        storageForG.resize(2,nDofs);
        nu.resize(nDofs);
        epsilon.resize(nDofs);
        z.resize(nBodies);

        Gepsilon.resize(nBodies); // TODO: ground initialization
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
 *     Velocity:        rates
 *     Dynamics:        forces & prescribed accelerations   
 *
class SimbodyTreeVariables {
public:
    // Modeling
    bool              useEulerAngles;
    std::vector<bool> prescribed;           // nb  (# bodies & joints, 0 always true)
    std::vector<bool> enabled;              // nac (# acceleration constraints)

    // Parametrization
    // TODO: body masses, etc.

    // Configuration
    Vector q;                               // nq

    // Motion
    Vector u;                               // nu  (== ndof)

    // Dynamics
    Vector_<SpatialVec> appliedBodyForces;  // nb
    Vector              appliedJointForces; // nu
    Vector              prescribedUdot;     // nu

public:
    SimbodyTreeVariables() : useEulerAngles(false) { }

    // nDofs==nu==#joint forces

    // Call this after realizeConstruction(). These are the variables we need
    // to specify our modeling choices. We can't allocate the rest until we see
    // how we'll be modeling.
    void allocateModelingVars(int nBodies, int nConstraints) {
        useEulerAngles = false;
        prescribed.resize(nBodies); 
        prescribed.assign(nBodies,false); prescribed[0]=true; // ground
        enabled.resize(nConstraints);
        enabled.assign(nConstraints,false);
    }

    // Call this after realizeModeling(). We now know everything we need to know
    // to allocate and initialize the remaining state variables.
    void allocateAllVars(int nDofs, int maxNQs, int nac) // acc constraints
    {
        const int nBodies      = prescribed.size();
        const int nConstraints = enabled.size();
        q.resize(maxNQs);                   q.setToNaN();
        u.resize(nDofs);                    u.setToNaN();

        appliedBodyForces.resize(nBodies);  appliedBodyForces.setToNaN();
        appliedBodyForces[0] = SpatialVec(Vec3(0)); // ground

        appliedJointForces.resize(nDofs);   appliedJointForces.setToNaN();
        prescribedUdot.resize(nDofs);       prescribedUdot.setToNaN();
    }

    void setVelocitiesToZero() {
        u.setToZero();
    }

    void clearForces() {
        appliedBodyForces.setToZero();
        appliedJointForces.setToZero();
    }

    // This locks all the joints that are prescribed.
    void setPrescribedAccelerationsToZero() {
        prescribedUdot.setToZero();
    }
};
*/

class SBModelingVars {
public:
    bool              useEulerAngles;
    std::vector<bool> prescribed;           // nb  (# bodies & joints, 0 always true)
    std::vector<bool> enabled;              // nac (# acceleration constraints)
public:

    // We have to allocate and assign defaults to these without looking at any other
    // state variable or cache entries. We can only depend on the tree itself for information.
    void allocate(const RigidBodyTree& tree) const {
        SBModelingVars& mutvars = *const_cast<SBModelingVars*>(this);

        assert(tree.isBuilt());
        mutvars.useEulerAngles = false;
        mutvars.prescribed.resize(tree.getNBodies()); 
        mutvars.prescribed.assign(tree.getNBodies(),false); 
        mutvars.prescribed[0]=true; // ground
        mutvars.enabled.resize(tree.getNConstraints());
        mutvars.enabled.assign(tree.getNConstraints(),false);
    }

};

class SBParameterVars {
public:
    // TODO: body masses, etc.
public:

    // We can access the tree or state variable & cache up to Modeling stage.
    void allocate(const RigidBodyTree&, const SBStateRep&) const {
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
    // We can access the tree or state variable & cache up to Timed stage.
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
    // We can access the tree or state variable & cache up to Configured stage.
    void allocate(const RigidBodyTree& tree, const SBStateRep&) const {
        SBMotionVars& mutvars = *const_cast<SBMotionVars*>(this);
        mutvars.u.resize(tree.getTotalDOF()); 
        mutvars.u.setToNaN();
    }
};

class SBDynamicVars {
public:
    Vector_<SpatialVec> appliedBodyForces;  // nb
    Vector              appliedJointForces; // nu
    Vector              prescribedUdot;     // nu
public:

    // We can access the tree or state variable & cache up to Moving stage.
    void allocate(const RigidBodyTree& tree, const SBStateRep&) const {    
        SBDynamicVars& mutvars = *const_cast<SBDynamicVars*>(this);

        mutvars.appliedBodyForces.resize(tree.getNBodies());  
        mutvars.appliedBodyForces.setToNaN();
        mutvars.appliedBodyForces[0] = SpatialVec(Vec3(0)); // ground

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
        assert(handle==0);
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

    // Allocating vars for the next stage g is the last step in realizing stage g-1.
    // So if we find that the state is *already* at g-1 or higher we can return immediately.
    // We must be able to depend on stage g-2 and lower being available, and we have
    // a special dispensation here to look at stage g-1 variables (which we know are
    // ready), even though the stage isn't there yet.

    void allocateVarsIfNeeded(const RigidBodyTree& tree, SBStage next) const {
        assert(next > BuiltStage); // "built" vars are in the RigidBodyTree, not the state
        assert(getStage(tree) >= next-2);
        if (getStage(tree) >= next-1)
            return;

        SBStateRep& mutableState = *const_cast<SBStateRep*>(this);

        switch (next) {
        case ModeledStage:      modelVars.allocate  (tree); 
                                tree.setDefaultModelingValues     (*this, mutableState.modelVars); 
                                break;
        case ParametrizedStage: paramVars.allocate  (tree, *this);  
                                tree.setDefaultParameterValues    (*this, mutableState.paramVars); 
                                break;
        case TimedStage:        timeVars.allocate   (tree, *this);   
                                tree.setDefaultTimeValues         (*this, mutableState.timeVars); 
                                break;
        case ConfiguredStage:   configVars.allocate (tree, *this);   
                                tree.setDefaultConfigurationValues(*this, mutableState.configVars); 
                                break;
        case MovingStage:       motionVars.allocate (tree, *this);   
                                tree.setDefaultMotionValues       (*this, mutableState.motionVars); 
                                break;
        case ReactingStage:     dynamicVars.allocate(tree, *this);   
                                tree.setDefaultDynamicValues      (*this, mutableState.dynamicVars); 
                                break;
        default: assert(false);
        }
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
        case ModeledStage:      modelCache.allocate  (tree); break;
        case ParametrizedStage: paramCache.allocate  (tree, *this); break;
        case TimedStage:        timeCache.allocate   (tree, *this); break;
        case ConfiguredStage:   configCache.allocate (tree, *this); break;
        case MovingStage:       motionCache.allocate (tree, *this); break;
        case ReactingStage:     dynamicCache.allocate(tree, *this); break;
        default: assert(false);
        }
    }

        // STAGE-CHECKED VARIABLE ACCESS

    const SBModelingVars& getModelingVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ModeledStage);
        return modelVars;
    }
    SBModelingVars& updModelingVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ModeledStage-1);
        stage = SBStage(ModeledStage-1); // backup if necessary
        return modelVars;
    }
    const SBParameterVars& getParameterVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ParametrizedStage);
        return paramVars;
    }
    SBParameterVars& updParameterVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ParametrizedStage-1);
        stage = SBStage(ParametrizedStage-1); // backup if necessary
        return paramVars;
    }
    const SBTimeVars& getTimeVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= TimedStage);
        return timeVars;
    }
    SBTimeVars& updTimeVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= TimedStage-1);
        stage = SBStage(TimedStage-1); // backup if necessary
        return timeVars;
    }
    const SBConfigurationVars& getConfigurationVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ConfiguredStage);
        return configVars;
    }
    SBConfigurationVars& updConfigurationVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ConfiguredStage-1);
        stage = SBStage(ConfiguredStage-1); // backup if necessary
        return configVars;
    }
    const SBMotionVars& getMotionVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= MovingStage);
        return motionVars;
    }
    SBMotionVars& updMotionVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= MovingStage-1);
        stage = SBStage(MovingStage-1); // backup if necessary
        return motionVars;
    }
    const SBDynamicVars& getDynamicVars(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ReactingStage);
        return dynamicVars;
    }
    SBDynamicVars& updDynamicVars(const RigidBodyTree& tree) {
        assert(getStage(tree) >= ReactingStage-1);
        stage = SBStage(ReactingStage-1); // backup if necessary
        return dynamicVars;
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
    const SBDynamicCache& getDynamicCache(const RigidBodyTree& tree, 
                                          bool inProgress=false) const 
    {
        assert(getStage(tree) >= (inProgress ? ReactingStage-1 : ReactingStage));
        return dynamicCache;
    }
    SBDynamicCache& updDynamicCache(const RigidBodyTree& tree) const {
        assert(getStage(tree) >= ReactingStage-1);
        return dynamicCache;
    }
private:
    SBModelingVars      modelVars; // allocation sizes & default values only are mutable
    SBParameterVars     paramVars;
    SBTimeVars          timeVars;
    SBConfigurationVars configVars;
    SBMotionVars        motionVars;
    SBDynamicVars       dynamicVars;
    
    mutable SBStage              stage; // last stage completed

    mutable SBModelingCache      modelCache;
    mutable SBParameterCache     paramCache;
    mutable SBTimeCache          timeCache;
    mutable SBConfigurationCache configCache;
    mutable SBMotionCache        motionCache;
    mutable SBDynamicCache       dynamicCache;

    SBState* handle;
    friend class RigidBodyTree;
    friend class RigidBodyNode;
    template <int dof> friend class RigidBodyNodeSpec;
};

}; // namespace simtk

#endif // SIMTK_SIMBODY_TREE_STATE_H_
