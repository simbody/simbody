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


using SimTK::State; using SimTK::Array; using SimTK::Vector; using SimTK::Vec3;
using SimTK::Rotation; using SimTK::Transform; using SimTK::Inertia;
using SimTK::Real; using SimTK::Vector_;
using SimTK::UnitVec3; using SimTK::SpatialVec; using SimTK::SpatialRow; 
using SimTK::SpatialMat; using SimTK::Matrix_;
using SimTK::PhiMatrix; using SimTK::Mat33; using SimTK::MassProperties; using SimTK::Stage;

class State;

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
class SBTopologyCache {
public:
    SBTopologyCache() {
        clear();
    }

    void clear() {
        nBodies = nParticles = nConstraints = nDOFs = maxNQs = sumSqDOFs =
            nDistanceConstraints = nPointInPlaneConstraints = 
            modelingVarsIndex = modelingCacheIndex = -1;
        valid = false;
    }

    int nBodies;
    int nParticles;
    int nConstraints;

    int nDOFs;
    int maxNQs;
    int sumSqDOFs;

    int nDistanceConstraints;
    int nPointInPlaneConstraints;

    int modelingVarsIndex;
    int modelingCacheIndex;

    bool valid;
};

class SBModelCache {
public:
    // TODO: Modeling
    //   counts of various things resulting from modeling choices,
    //   constraint enabling, prescribed motion

    SBModelCache() {
        clear();
    }

    void clear() {
        instanceVarsIndex = instanceCacheIndex
        = timeVarsIndex = timeCacheIndex
        = qIndex = qVarsIndex = qCacheIndex
        = uIndex = uVarsIndex = uCacheIndex
        = dynamicsVarsIndex = dynamicsCacheIndex
        = accelerationVarsIndex = accelerationCacheIndex
        = nQuaternionsInUse = firstQuaternionQErrSlot 
        = qErrIndex = uErrIndex = udotErrIndex = -1;
        quaternionIndex.clear();
    }

    int instanceVarsIndex, instanceCacheIndex;
    int timeVarsIndex, timeCacheIndex;
    int qIndex; // maxNQs of these 
    int qVarsIndex, qCacheIndex;
    int uIndex; // nDOFs of these 
    int uVarsIndex, uCacheIndex;
    int dynamicsVarsIndex, dynamicsCacheIndex;
    int accelerationVarsIndex, accelerationCacheIndex;

    int nQuaternionsInUse, firstQuaternionQErrSlot;
    Array<int> quaternionIndex; // nb (-1 for bodies w/no quats)
    int qErrIndex, uErrIndex, udotErrIndex;

public:
    void allocate(const SBTopologyCache& tc) {
        quaternionIndex.resize(tc.nBodies, -1); // init to -1
    }
};

class SBInstanceCache {
public:
    // Instance variables are:
    //   body mass props; particle masses
    //   X_BM, X_PMb transforms
    //   distance constraint distances & station positions
    // Calculations stored here derive from those states:
    //   total mass
    //   central inertia of each rigid body
    //   principal axes and corresponding principal moments of inertia of each rigid body
    //   reference configuration X_PB when q==0 (usually that means M==Mb), for each rigid body

    Real             totalMass; // sum of all rigid body and particles masses
    Array<Inertia>   centralInertias;           // nb
    Vector_<Vec3>    principalMoments;          // nb
    Array<Rotation>  principalAxes;             // nb
    Array<Transform> referenceConfiguration;    // nb

public:
    void allocate(const SBTopologyCache& topology) {
        totalMass = SimTK::NaN;
        centralInertias.resize(topology.nBodies);           // I_CB
        principalMoments.resize(topology.nBodies);          // (Ixx,Iyy,Izz)
        principalAxes.resize(topology.nBodies);             // [axx ayy azz]
        referenceConfiguration.resize(topology.nBodies);    // X0_PB
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

    Matrix_<Vec3> storageForHtFM; // 2 x ndof (~H_FM)
    Matrix_<Vec3> storageForHt;    // 2 x ndof (~H_PB_G)

    Array<Transform>    bodyJointInParentJointFrame;  // nb (X_FM)

    Array<Transform>    bodyConfigInParent;           // nb (X_PB)
    Array<Transform>    bodyConfigInGround;           // nb (X_GB)
    Array<PhiMatrix>    bodyToParentShift;            // nb (phi)
    Array<Inertia>      bodyInertiaInGround;          // nb (I_OB_G)
    Vector_<SpatialMat> bodySpatialInertia;           // nb (Mk)
    Vector_<Vec3>       bodyCOMInGround;              // nb (COM_G)
    Vector_<Vec3>       bodyCOMStationInGround;       // nb (COMstation_G)

    // Distance constraint calculations. These are indexed by
    // *distance constraint* number, not *constraint* number.
    Vector_<Vec3> station_G[2];   // vec from body origin OB to station, expr. in G
    Vector_<Vec3> pos_G[2];       // vec from ground origin to station, expr. in G

    Vector_<Vec3> fromTip1ToTip2_G; // tip2.pos-tip1.pos
    Vector_<Vec3> unitDirection_G;  // fromTip1ToTip2/|fromTip1ToTip2|

    // Point-in-plane constraint calculations. These are indexed by *point-in-plane*
    // constraint number, not *constraint* number or *multiplier* number.
    Vector_<UnitVec3> pipNormal_G; // the body 1-fixed plane normal, expressed in G
    Vector_<Vec3>     pipStation_G; // the body 2-fixed vector to the station, expressed inG
    Vector_<Vec3>     pipPos_G;     // the body 2-fixed station's position in Cartesian space
    Vector_<Vec3>     pipPosInPlaneBody_G; // follower station measured from plane origin, exp. in G
    Vector            pipHeight; // height of body 2 station along body 1 normal


public:
    void allocate(const SBTopologyCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int ndistc  = tree.nDistanceConstraints;
        const int npipc   = tree.nPointInPlaneConstraints;
        const int npc     = ndistc + npipc; // all the position constraints

        // These contain uninitialized junk. Body-indexed entries get their
        // ground elements set appropriately now and forever.
        sq.resize(maxNQs);
        cq.resize(maxNQs);
        qnorm.resize(maxNQs);

        storageForHtFM.resize(2,nDofs);
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

        station_G[0].resize(ndistc); station_G[1].resize(ndistc);
        pos_G[0].resize(ndistc); pos_G[1].resize(ndistc);

        fromTip1ToTip2_G.resize(ndistc);
        unitDirection_G.resize(ndistc);

        pipNormal_G.resize(npipc);
        pipStation_G.resize(npipc);
        pipPos_G.resize(npipc);
        pipPosInPlaneBody_G.resize(npipc);
        pipHeight.resize(npipc);
    }
};

class SBVelocityCache {
public:
    // qdot is supplied directly by the State
    Vector_<SpatialVec> bodyVelocityInParent;      // nb (joint velocity)
    Vector_<SpatialVec> bodyVelocityInGround;      // nb (sVel)
    Vector_<SpatialVec> mobilizerRelativeVelocity; // nb (V_FM)

    // Distance constraint calculations. These are indexed by
    // *distance constraint* number, not *constraint* number.
    Vector_<Vec3> stationVel_G[2]; // vel of station relative to body origin, expr. in G
    Vector_<Vec3> vel_G[2];        // tip velocities relative to G, expr. in G
    Vector_<Vec3> relVel_G;        // spatial relative velocity tip2.velG-tip1.velG

    // Point-in-plane constraint calculations. These are indexed by *point-in-plane*
    // constraint number, not *constraint* number or *multiplier* number.
    Vector_<Vec3> pipNormalDot_G;  // time derivative of the body-1 fixed plane normal, taken in G
    Vector_<Vec3> pipVel_G;        // vel of station in G, expr. in G
    Vector_<Vec3> pipVelInPlaneBody_G; // vel of station relative to plane origin, expr. in G
    Vector        pipHeightDot;    // rate of change of height (ideally zero)

public:
    void allocate(const SBTopologyCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int ndistc  = tree.nDistanceConstraints;
        const int npipc   = tree.nPointInPlaneConstraints;
        const int nvc     = ndistc + npipc; // all the velocity constraints

        bodyVelocityInParent.resize(nBodies);       
        bodyVelocityInParent[0] = SpatialVec(Vec3(0),Vec3(0));

        bodyVelocityInGround.resize(nBodies);       
        bodyVelocityInGround[0] = SpatialVec(Vec3(0),Vec3(0));

        mobilizerRelativeVelocity.resize(nBodies);       
        mobilizerRelativeVelocity[0] = SpatialVec(Vec3(0),Vec3(0));

        stationVel_G[0].resize(ndistc); stationVel_G[1].resize(ndistc);
        vel_G[0].resize(ndistc); vel_G[1].resize(ndistc);
        relVel_G.resize(ndistc);

        pipNormalDot_G.resize(npipc);
        pipVel_G.resize(npipc);
        pipVelInPlaneBody_G.resize(npipc);
        pipHeightDot.resize(npipc);
    }
};

class SBDynamicsCache {
public:
    // Dynamics
    Vector_<SpatialMat> articulatedBodyInertia;   // nb (P)

    Matrix_<Vec3> storageForHtFMDot; // 2 x ndof (~H_FM_Dot)
    Matrix_<Vec3> storageForHtDot;    // 2 x ndof (~H_PB_G_Dot)

    Vector_<SpatialVec> bodyVelocityInParentDerivRemainder; // VB_PB_G=~H_PB_G_Dot*u

    Vector_<SpatialVec> coriolisAcceleration;     // nb (a)
    Vector_<SpatialVec> totalCoriolisAcceleration;// nb (A)
    Vector_<SpatialVec> gyroscopicForces;         // nb (b)
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
        const int nac     = tree.nDistanceConstraints; // acceleration constraints        
        
        articulatedBodyInertia.resize(nBodies); // TODO: ground initialization
        storageForHtFMDot.resize(2,nDofs);
        storageForHtDot.resize(2,nDofs);

        bodyVelocityInParentDerivRemainder.resize(nBodies);       
        bodyVelocityInParentDerivRemainder[0] = SpatialVec(Vec3(0),Vec3(0));

        coriolisAcceleration.resize(nBodies);       
        coriolisAcceleration[0] = SpatialVec(Vec3(0),Vec3(0));

        totalCoriolisAcceleration.resize(nBodies);       
        totalCoriolisAcceleration[0] = SpatialVec(Vec3(0),Vec3(0));

        gyroscopicForces.resize(nBodies);           
        gyroscopicForces[0] = SpatialVec(Vec3(0),Vec3(0));

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
    // udot, qdotdot are provided directly by the State
    Vector_<SpatialVec> bodyAccelerationInGround; // nb (sAcc)
    Vector              netHingeForces;           // nu (T-(~Am+R(F+C))

    Vector              nu;
    Vector              epsilon;
    Vector_<SpatialVec> z;                        // nb
    Vector_<SpatialVec> Gepsilon;                 // nb

    // Distance constraint calculations. These are indexed by
    // *distance constraint* number, not *constraint* number.
    Vector_<Vec3> acc_G[2];   // acc of tip relative to ground, expr. in G
    Vector_<Vec3> force_G[2]; // the constraint forces applied to each point

    // Point-in-plane constraint calculations. These are indexed by *point-in-plane*
    // constraint number, not *constraint* number or *multiplier* number.
    Vector_<Vec3> pipNormalDotDot_G; // 2nd time derivative of the body-1 fixed plane normal, taken in G
    Vector_<Vec3> pipAcc_G;          // acc of station in G, expr. in G
    Vector_<Vec3> pipAccInPlaneBody_G; // acc of station relative to plane origin, expr. in G
    Vector_<Vec3> pipForce_G[2];     // the constraint forces applied to each point (station on body2
                                     //   and the coincident point on body1)
    Vector        pipHeightDotDot;   // 2nd deriv of height (ideally zero)

public:
    void allocate(const SBTopologyCache& tree) {
        // Pull out construction-stage information from the tree.
        const int nBodies = tree.nBodies;
        const int nDofs   = tree.nDOFs;     // this is the number of u's (nu)
        const int nSqDofs = tree.sumSqDOFs;   // sum(ndof^2) for each joint
        const int maxNQs  = tree.maxNQs;  // allocate the max # q's we'll ever need
        const int ndistc  = tree.nDistanceConstraints;
        const int npipc   = tree.nPointInPlaneConstraints;
        const int nac     = ndistc + npipc; // all the acceleration constraints

        bodyAccelerationInGround.resize(nBodies);   
        bodyAccelerationInGround[0] = SpatialVec(Vec3(0),Vec3(0));;

        netHingeForces.resize(nDofs);

        nu.resize(nDofs);
        epsilon.resize(nDofs);
        z.resize(nBodies);

        Gepsilon.resize(nBodies); // TODO: ground initialization

        acc_G[0].resize(ndistc); acc_G[1].resize(ndistc);
        force_G[0].resize(ndistc); force_G[1].resize(ndistc);

        pipNormalDotDot_G.resize(npipc);
        pipAcc_G.resize(npipc);
        pipAccInPlaneBody_G.resize(npipc);
        pipForce_G[0].resize(npipc); pipForce_G[1].resize(npipc);
        pipHeightDotDot.resize(npipc);
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
    bool        useEulerAngles;
    Array<bool> prescribed;           // nb  (# bodies & mobilizers, [0] always true)
    Array<bool> enabled;              // nac (# acceleration constraints)
public:

    // We have to allocate these without looking at any other
    // state variable or cache entries. We can only depend on topological
    // information.
    void allocate(const SBTopologyCache& tree) const {
        SBModelVars& mutvars = *const_cast<SBModelVars*>(this);
        mutvars.useEulerAngles = false;
        mutvars.prescribed.resize(tree.nBodies); 
        mutvars.enabled.resize(tree.nConstraints);
    }

};

class SBInstanceVars {
public:
    Array<MassProperties> bodyMassProperties;
    Array<Transform>      outboardMobilizerFrames;
    Array<Transform>      inboardMobilizerFrames;
    Vector                particleMasses;

public:

    // We can access the tree or state variable & cache up to Modeling stage.
    void allocate(const SBTopologyCache& topology) const {
        SBInstanceVars& mutvars = *const_cast<SBInstanceVars*>(this);
        mutvars.bodyMassProperties.resize(topology.nBodies);
        mutvars.outboardMobilizerFrames.resize(topology.nBodies);
        mutvars.inboardMobilizerFrames.resize (topology.nBodies);
        mutvars.particleMasses.resize(topology.nParticles);

        // Set default values
        mutvars.bodyMassProperties      = MassProperties(1,Vec3(0),Inertia(1));
        mutvars.outboardMobilizerFrames = Transform(); // i.e., B frame
        mutvars.inboardMobilizerFrames  = Transform(); // i.e., P frame
        mutvars.particleMasses          = Real(1);
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
inline std::ostream& operator<<(std::ostream& o, const SBModelCache& c)
  { return o << "TODO: SBModelCache"; }
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
    explicit SBStateDigest(const State& s) : state(s) 
    {
        stage = Stage::Empty;
    }
    SBStateDigest(const State& s, const SimbodyMatterSubsystemRep& matter, Stage g)
      : state(s) 
    {
        stage = Stage::Empty;
        fillThroughStage(matter,g);
    }

    // Stage g here is the stage we are about to compute. So we expect the referenced
    // State to have been realized to at least stage g-1.
    void fillThroughStage(const SimbodyMatterSubsystemRep& matter, Stage g);

    // The State is read only, for cache entries you have a choice.

    const State& getState() const {return state;}
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

    const Real* getQ() const {
        assert(stage >= Stage::Position);
        assert(q);
        return q;
    }

    const SBPositionVars& getPositionVars() const {
        assert(stage >= Stage::Position);
        assert(pv);
        return *pv;
    }

    const Real* getU() const {
        assert(stage >= Stage::Velocity);
        assert(u);
        return u;
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
        mv=0; iv=0; tv=0; q=0; pv=0; u=0; vv=0; dv=0; av=0;
        mc=0; ic=0; tc=0;      pc=0;      vc=0; dc=0; ac=0;
    }

private:
    const State& state;
    Stage        stage; // the stage to be computed

    const SBModelVars*          mv;
    const SBInstanceVars*       iv;
    const SBTimeVars*           tv;

    const Real*                 q;
    const SBPositionVars*       pv;

    const Real*                 u;
    const SBVelocityVars*       vv;
    const SBDynamicsVars*       dv;
    const SBAccelerationVars*   av;

    SBModelCache*               mc;
    SBInstanceCache*            ic;
    SBTimeCache*                tc;
    SBPositionCache*            pc;
    SBVelocityCache*            vc;
    SBDynamicsCache*            dc;
    SBAccelerationCache*        ac;
};

#endif // SimTK_SIMBODY_TREE_STATE_H_
