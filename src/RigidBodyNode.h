#ifndef SimTK_RIGID_BODY_NODE_H_
#define SimTK_RIGID_BODY_NODE_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors: Derived from NIH IVM code written by Charles Schwieters.
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

#include "simbody/internal/common.h"
#include "SimbodyTreeState.h"

#include "internalDynamics.h"

#include <cassert>
#include <vector>

using namespace SimTK;

/**
 * This is an abstract class representing a body and its (generic) mobilizer, that is,
 * the joint connecting it to its parent. Concrete classes are derived from this one to
 * represent each specific type of mobilizer.
 *
 * RigidBodyNodes are linked into a tree structure, organized into levels as described 
 * in Schwieters' JMR paper. The root is a special 'Ground' node defined to be at 
 * level 0. The level 1 nodes (referred to
 * as 'base nodes') are those attached directly to the Ground node, level 2's attach 
 * to level 1's, etc. Every node but Ground has exactly one parent node, whose
 * level is always one less than the current node. Any node may have an arbitrary number 
 * of children, for which it is the unique parent, and all of its children have 
 * level one greater than the current node.
 *
 * Note: calling rotation matrices 'rotMat' or 'R' is a recipe for disaster.
 * We use the naming convention R_XY to mean a rotation matrix (3x3 direction
 * cosine matrix) expressing the orientation of frame Y in frame X. Given a vector
 * vY expressed in Y, you can re-express that in X via vX=R_XY*vY. To go the
 * other direction use R_YX=~R_XY (tilde '~' is SimTK for "transpose"). This convention
 * provides flawless composition of rotations as long as you "match up" the frames
 * (i.e., adjacent frame symbols must be the same). So if you have R_XY and R_ZX 
 * and you want R_YZ you can easily get it like this:
 *     R_YZ = R_YX*R_XZ = (~R_XY)*(~R_ZX) = ~(R_ZX*R_XY).
 * Also note that these are orthogonal matrices so R_XY*R_YX=I.
 *
 * Every body has a body frame B, and an inboard mobilizer frame M. For convenience, we
 * refer to the body frame of a body's unique parent as the 'P' frame. There is
 * a frame Mb on P which is where B's inboard joint attaches. When all the mobilizer
 * coordinates are 0 (=1000 for quaternions), M==Mb. The transform X_MbM tracks the
 * across-mobilizer change in configuration induced by the generalized coordinates q.
 *
 * The mobilizer frame M is fixed with respect to B, and Mb is fixed with
 * respect to P. In some cases M and B or Mb and P will be the same, but not always.
 * The constant transforms X_BM and X_PMb provide the configuration of the mobilizer
 * frames with respect to their body frames. With these definitions we can
 * easily calculate X_PB as X_PB = X_PMb*X_MbM*X_MB, where X_MbM is the q-dependent
 * cross-mobilizer transform calculated at Position stage.
 *
 * RigidBodyNodes know how to extract and deposit their own information from and
 * to the Simbody State variables and cache entries, but they don't know anything
 * about the State class, stages, etc. Instead they depend on being given appropriate
 * access by the caller, whose job it is to mine the State.
 */
class RigidBodyNode {
public:
    class VirtualBaseMethod {};    // an exception

    virtual ~RigidBodyNode() {}

    RigidBodyNode& operator=(const RigidBodyNode&);

    static RigidBodyNode* createGroundNode();

    // Register the passed-in node as a child of this one.
    void addChild(RigidBodyNode* child);
    void setParent(RigidBodyNode* p) {parent=p;}
    void setNodeNum(int n) {nodeNum=n;}
    void setLevel(int i)   {level=i;}

        // TOPOLOGICAL INFO: no State needed

    RigidBodyNode* getParent() const {return parent;}
    int            getNChildren()  const {return (int)children.size();}
    RigidBodyNode* getChild(int i) const {return (i<(int)children.size()?children[i]:0);}

    // Return this node's level, that is, how many ancestors separate it from
    // the Ground node at level 0. Level 1 nodes (directly connected to the
    // Ground node) are called 'base' nodes.
    int  getLevel() const  {return level;}

    // This is the unique "body Id" for this (body,mobilizer) node. It is used to index
    // arrays of body quantities.
    int  getNodeNum() const {return nodeNum;}

    bool isGroundNode() const { return level==0; }
    bool isBaseNode()   const { return level==1; }

    int  getUIndex() const {return uIndex;}
    int  getQIndex() const {return qIndex;}


    // Access routines for plucking the right per-body data from the pool in the State.
    const Transform&  fromB(const std::vector<Transform>& x) const {return x[nodeNum];}
    const Transform&  fromB(const Array<Transform>&       x) const {return x[nodeNum];}
    const PhiMatrix&  fromB(const std::vector<PhiMatrix>& p) const {return p[nodeNum];}
    const PhiMatrix&  fromB(const Array<PhiMatrix>&       p) const {return p[nodeNum];}
    const MassProperties& fromB(const std::vector<MassProperties>& m) const {return m[nodeNum];}
    const MassProperties& fromB(const Array<MassProperties>&       m) const {return m[nodeNum];}
    const Inertia&    fromB(const std::vector<Inertia>&   i) const {return i[nodeNum];}
    const Inertia&    fromB(const Array<Inertia>&         i) const {return i[nodeNum];}
    int               fromB(const std::vector<int>&       i) const {return i[nodeNum];}
    int               fromB(const Array<int>&             i) const {return i[nodeNum];}
    const SpatialVec& fromB(const Vector_<SpatialVec>&    v) const {return v[nodeNum];}
    const SpatialMat& fromB(const Vector_<SpatialMat>&    m) const {return m[nodeNum];}
    const Vec3&       fromB(const Vector_<Vec3>&          v) const {return v[nodeNum];}

    Transform&  toB(std::vector<Transform>& x) const {return x[nodeNum];}
    Transform&  toB(Array<Transform>&       x) const {return x[nodeNum];}
    PhiMatrix&  toB(std::vector<PhiMatrix>& p) const {return p[nodeNum];}
    PhiMatrix&  toB(Array<PhiMatrix>&       p) const {return p[nodeNum];}
    MassProperties& toB(std::vector<MassProperties>& m) const {return m[nodeNum];}
    MassProperties& toB(Array<MassProperties>&       m) const {return m[nodeNum];}
    Inertia&    toB(std::vector<Inertia>&   i) const {return i[nodeNum];}
    Inertia&    toB(Array<Inertia>&         i) const {return i[nodeNum];}
    int&        toB(std::vector<int>&       i) const {return i[nodeNum];}
    int&        toB(Array<int>&             i) const {return i[nodeNum];}
    SpatialVec& toB(Vector_<SpatialVec>&    v) const {return v[nodeNum];}
    SpatialMat& toB(Vector_<SpatialMat>&    m) const {return m[nodeNum];}
    Vec3&       toB(Vector_<Vec3>&          v) const {return v[nodeNum];}

        // MODELING INFO
    bool getUseEulerAngles(const SBModelVars& mv) const {return mv.useEulerAngles;}
    bool isPrescribed     (const SBModelVars& mv) const {return mv.prescribed[nodeNum];}

    int  getQuaternionIndex(const SBModelCache& mc) const {
        return mc.quaternionIndex[nodeNum];
    }

        // INSTANCE INFO

    // TODO: These ignore State currently since they aren't parametrizable.
    const MassProperties& getMassProperties_OB_B() const {return massProps_B;}
    const Real&           getMass          () const {return massProps_B.getMass();}
    const Vec3&           getCOM_B         () const {return massProps_B.getMassCenter();}
    const Inertia&        getInertia_OB_B  () const {return massProps_B.getInertia();}
    const Transform&      getX_BM          () const {return X_BM;}
    const Transform&      getX_PMb         () const {return X_PMb;}

    // These are calculated on construction.
    // TODO: principal axes
    const Inertia&        getInertia_CB_B  () const {return inertia_CB_B;}
    const Transform&      getX_MB          () const {return X_MB;}
    const Transform&      getRefX_PB       () const {return refX_PB;}

        // POSITION INFO

    /// Extract from the cache  X_MbM, the cross-mobilizer transformation matrix giving the configuration
    /// of this body's mobilizer frame M, measured from and expressed in the corresponding outboard
    /// mobilizer frame Mb attached to the parent. This transformation is defined to be zero (that is, Mb=M)
    /// in the reference configuration where the joint coordinates are all 0 (or 1,0,0,0 for quaternions).
    /// This is NOT a spatial (ground frame) transformation.
    const Transform& getX_MbM(const SBPositionCache& cc) const {return fromB(cc.bodyJointInParentJointFrame);}
    Transform&       updX_MbM(SBPositionCache&       cc) const {return toB  (cc.bodyJointInParentJointFrame);}

    /// Extract from the cache  X_PB, the cross-joint transformation matrix giving the configuration
    /// of this body's frame B measured from and expressed in its *parent* frame P. Thus this is NOT
    /// a spatial (ground frame) transformation.
    const Transform& getX_PB(const SBPositionCache& cc) const {return fromB(cc.bodyConfigInParent);}
    Transform&       updX_PB(SBPositionCache&       cc) const {return toB  (cc.bodyConfigInParent);}

    /// Extract from the cache X_GB, the transformation matrix giving the spatial configuration of this
    /// body's frame B measured from and expressed in ground. This consists of a rotation matrix
    /// R_GB, and a ground-frame vector r_OG_OB from ground's origin to the origin point of frame B.
    const Transform& getX_GB(const SBPositionCache& cc) const {
        return fromB(cc.bodyConfigInGround);
    }
    Transform& updX_GB(SBPositionCache& cc) const {
        return toB(cc.bodyConfigInGround);
    }

    /// Extract from the cache the body-to-parent shift matrix "phi". 
    const PhiMatrix& getPhi(const SBPositionCache& cc) const {return fromB(cc.bodyToParentShift);}
    PhiMatrix&       updPhi(SBPositionCache&       cc) const {return toB  (cc.bodyToParentShift);}

    /// Extract this body's spatial inertia matrix from the cache. This contains the mass properties
    /// measured from (and about) the body frame origin, but expressed in the *ground* frame.
    const SpatialMat& getMk(const SBPositionCache& cc) const {return fromB(cc.bodySpatialInertia);}
    SpatialMat&       updMk(SBPositionCache&       cc) const {return toB  (cc.bodySpatialInertia);}

    /// Extract from the cache the location of the body's center of mass, measured from the ground
    /// origin and expressed in ground.
    const Vec3& getCOM_G(const SBPositionCache& cc) const {return fromB(cc.bodyCOMInGround);}
    Vec3&       updCOM_G(SBPositionCache&       cc) const {return toB  (cc.bodyCOMInGround);}

    /// Extract from the cache the vector from body B's origin to its center of mass, reexpressed in Ground.
    const Vec3& getCB_G(const SBPositionCache& cc) const {return fromB(cc.bodyCOMStationInGround);}
    Vec3&       updCB_G(SBPositionCache&       cc) const {return toB  (cc.bodyCOMStationInGround);}

    /// Extract from the cache the body's inertia about the body origin OB, but reexpressed in Ground.
    const Inertia& getInertia_OB_G(const SBPositionCache& cc) const {return fromB(cc.bodyInertiaInGround);}
    Inertia&       updInertia_OB_G(SBPositionCache&       cc) const {return toB  (cc.bodyInertiaInGround);}

    /// Extract from the cache the spatial (ground-relative) location and orientation of this body's
    /// *parent's* body frame P.
    const Transform& getX_GP(const SBPositionCache& cc) const {assert(parent); return parent->getX_GB(cc);}

            // VELOCITY INFO

    const SpatialVec& getV_MbM(const SBVelocityCache& mc) const {return fromB(mc.mobilizerRelativeVelocity);}
    SpatialVec&       updV_MbM(SBVelocityCache&       mc) const {return toB  (mc.mobilizerRelativeVelocity);}


    /// Extract from the cache V_GB, the spatial velocity of this body's frame B measured in and
    /// expressed in ground. This contains the angular velocity of B in G, and the linear velocity
    /// of B's origin point OB in G, with both vectors expressed in G.
    const SpatialVec& getV_GB   (const SBVelocityCache& mc) const {return fromB(mc.bodyVelocityInGround);}
    SpatialVec&       updV_GB   (SBVelocityCache&       mc) const {return toB  (mc.bodyVelocityInGround);}

    /// Extract from the cache V_PB_G, the *spatial* velocity of this body's frame B, that is the
    /// cross-joint velocity measured with respect to the parent frame, but then expressed in the
    /// *ground* frame. This contains the angular velocity of B in P, and the linear velocity
    /// of B's origin point OB in P, with both vectors expressed in *G*.
    const SpatialVec& getV_PB_G (const SBVelocityCache& mc) const {return fromB(mc.bodyVelocityInParent);}
    SpatialVec&       updV_PB_G (SBVelocityCache&       mc) const {return toB  (mc.bodyVelocityInParent);}

    const SpatialVec& getV_GP(const SBVelocityCache& vc) const {assert(parent); return parent->getV_GB(vc);}

    const SpatialVec& getSpatialVel   (const SBVelocityCache& mc) const {return getV_GB(mc);}
    const Vec3&       getSpatialAngVel(const SBVelocityCache& mc) const {return getV_GB(mc)[0];}
    const Vec3&       getSpatialLinVel(const SBVelocityCache& mc) const {return getV_GB(mc)[1];}

        // DYNAMICS INFO
    /// Extract from the cache VD_PB_G, the *spatial* velocity derivative remainder term
    /// generated by H_PB_G_Dot*u, where H_PB_G_Dot = d/dt H_PB_G with the derivative taken
    /// in the Ground frame. This is used in calculation of coriolis acceleration.
    const SpatialVec& getVD_PB_G (const SBDynamicsCache& dc) const 
        {return fromB(dc.bodyVelocityInParentDerivRemainder);}
    SpatialVec&       updVD_PB_G (SBDynamicsCache&       dc) const 
        {return toB  (dc.bodyVelocityInParentDerivRemainder);}


    const SpatialMat& getP    (const SBDynamicsCache& dc) const {return fromB(dc.articulatedBodyInertia);}
    SpatialMat&       updP    (SBDynamicsCache&       dc) const {return toB  (dc.articulatedBodyInertia);}
 
    const SpatialVec& getCoriolisAcceleration(const SBDynamicsCache& dc) const {return fromB(dc.coriolisAcceleration);}
    SpatialVec&       updCoriolisAcceleration(SBDynamicsCache&       dc) const {return toB  (dc.coriolisAcceleration);}

    const SpatialVec& getTotalCoriolisAcceleration(const SBDynamicsCache& dc) const {return fromB(dc.totalCoriolisAcceleration);}
    SpatialVec&       updTotalCoriolisAcceleration(SBDynamicsCache&       dc) const {return toB  (dc.totalCoriolisAcceleration);}

    const SpatialVec& getGyroscopicForce(const SBDynamicsCache& dc) const {return fromB(dc.gyroscopicForces);}
    SpatialVec&       updGyroscopicForce(SBDynamicsCache&       dc) const {return toB  (dc.gyroscopicForces);}
 
    const SpatialVec& getCentrifugalForces(const SBDynamicsCache& dc) const {return fromB(dc.centrifugalForces);}
    SpatialVec&       updCentrifugalForces(SBDynamicsCache&       dc) const {return toB  (dc.centrifugalForces);}

    const SpatialVec& getTotalCentrifugalForces(const SBDynamicsCache& dc) const {return fromB(dc.totalCentrifugalForces);}
    SpatialVec&       updTotalCentrifugalForces(SBDynamicsCache&       dc) const {return toB  (dc.totalCentrifugalForces);}

    const SpatialMat& getArticulatedBodyInertia(const SBDynamicsCache& dc) const {return fromB(dc.articulatedBodyInertia);}
    SpatialMat&       updArticulatedBodyInertia(SBDynamicsCache&       dc) const {return toB  (dc.articulatedBodyInertia);}

    const SpatialVec& getZ(const SBAccelerationCache& rc) const {return fromB(rc.z);}
    SpatialVec&       updZ(SBAccelerationCache&       rc) const {return toB  (rc.z);}

    const SpatialVec& getGepsilon(const SBAccelerationCache& rc) const {return fromB(rc.Gepsilon);}
    SpatialVec&       updGepsilon(SBAccelerationCache&       rc) const {return toB  (rc.Gepsilon);}

    const SpatialMat& getPsi(const SBDynamicsCache& dc) const {return fromB(dc.psi);}
    SpatialMat&       updPsi(SBDynamicsCache&       dc) const {return toB  (dc.psi);}

    const SpatialMat& getTauBar(const SBDynamicsCache& dc) const {return fromB(dc.tauBar);}
    SpatialMat&       updTauBar(SBDynamicsCache&       dc) const {return toB  (dc.tauBar);}

    const SpatialMat& getY(const SBDynamicsCache& dc) const {return fromB(dc.Y);}
    SpatialMat&       updY(SBDynamicsCache&       dc) const {return toB  (dc.Y);}

        // ACCELERATION INFO

    const SpatialVec& getAppliedBodyForce(const SBAccelerationVars& av) const {return fromB(av.appliedRigidBodyForces);}
    SpatialVec&       updAppliedBodyForce(SBAccelerationVars& av)       const {return toB(av.appliedRigidBodyForces);}
  
    /// Extract from the cache A_GB, the spatial acceleration of this body's frame B measured in and
    /// expressed in ground. This contains the inertial angular acceleration of B in G, and the
    /// linear acceleration of B's origin point OB in G, with both vectors expressed in G.
    const SpatialVec& getA_GB (const SBAccelerationCache& ac) const {return fromB(ac.bodyAccelerationInGround);}
    SpatialVec&       updA_GB (SBAccelerationCache&       ac) const {return toB  (ac.bodyAccelerationInGround);}

    const SpatialVec& getSpatialAcc   (const SBAccelerationCache& ac) const {return getA_GB(ac);}
    const Vec3&       getSpatialAngAcc(const SBAccelerationCache& ac) const {return getA_GB(ac)[0];}
    const Vec3&       getSpatialLinAcc(const SBAccelerationCache& ac) const {return getA_GB(ac)[1];}

    virtual void realizeModel(
        const SBModelVars& mv,
        SBModelCache&      mc) const=0;

    virtual void realizeInstance(
        const SBModelVars&    mv,
        const SBInstanceVars& iv,
        SBInstanceCache&      ic) const=0;

    /// Introduce new values for generalized coordinates and calculate
    /// all the position-dependent kinematic terms, including position
    /// constraint errors. Must be called base to tip.
    virtual void realizePosition(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q,
        Vector&             qErr,
        SBPositionCache&    pc) const=0;

    /// Introduce new values for generalized speeds and calculate
    /// all the velocity-dependent kinematic terms. Assumes realizePosition()
    /// has already been called on all nodes. Must be called base to tip.
    virtual void realizeVelocity(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        SBVelocityCache&       vc,
        Vector&                qdot) const=0;

    /// Calculate base-to-tip velocity-dependent terms which will be used
    /// in Dynamics stage operators. Assumes realizeVelocity()
    /// has already been called on all nodes, as well as any Dynamics 
    /// stage calculations which must go tip-to-base (e.g. articulated
    /// body inertias).
    virtual void realizeDynamics(
        const SBModelVars&     mv,
        const SBPositionCache& pc,
        const Vector&          u,
        const SBVelocityCache& vc,
        SBDynamicsCache&       dc) const=0;

    // These are called just after new state variables are allocated,
    // in case there are any node-specific default values.
    // We can handle the "body" variables (like mass) here, but we have to forward the
    // request to the mobilizers to handle their own variables. At the Position
    // stage, for example, mobilizers which use quaternions will set the default ball 
    // joint q's to 1,0,0,0.
    // Most of these will use the default implementations here, i.e. do nothing.

    // Called after Model variables are allocated by realizeTopology()
    void setNodeDefaultModelValues(const SBTopologyCache& tc, SBModelVars& mv) const {
        // no body model variables at the moment; TODO: should forward to Body type
        setMobilizerDefaultModelValues(tc, mv);
    }

    // All the rest are called right after realizeModel() since that's when all the
    // remaining state variables are allocated.
    void setNodeDefaultInstanceValues(const SBModelVars& mv, SBInstanceVars& iv) const {
        // mass properties, inb and outb frame are handled here
        toB(iv.bodyMassProperties)      = getMassProperties_OB_B();
        toB(iv.outboardMobilizerFrames) = getX_BM();
        toB(iv.inboardMobilizerFrames)  = getX_PMb();
        setMobilizerDefaultInstanceValues(mv,iv);
    }
    void setNodeDefaultTimeValues(const SBModelVars& mv, SBTimeVars& tv)  const {
        // no body model variables at the moment; TODO: should forward to Body type
        setMobilizerDefaultTimeValues(mv, tv);
    }
    void setNodeDefaultPositionValues(const SBModelVars& mv, Vector& q) const {
        // no body model variables at the moment; TODO: should forward to Body type
        setMobilizerDefaultPositionValues(mv, q);
    }
    void setNodeDefaultVelocityValues(const SBModelVars& mv, Vector& u) const {
         // no body model variables at the moment; TODO: should forward to Body type
        setMobilizerDefaultVelocityValues(mv, u);
    }
    void setNodeDefaultDynamicsValues(const SBModelVars& mv, SBDynamicsVars& dv) const {
        // no body model variables at the moment; TODO: should forward to Body type
        setMobilizerDefaultDynamicsValues(mv, dv);
    }
    void setNodeDefaultAccelerationValues(const SBModelVars& mv, SBAccelerationVars& av) const {
         // no body model variables at the moment; TODO: should forward to Body type
        setMobilizerDefaultAccelerationValues(mv, av);
    }

    // Called after Model variables are allocated by realizeTopology()
    virtual void setMobilizerDefaultModelValues       (const SBTopologyCache&, SBModelVars&)        const {}

    // All the rest are called right after realizeModel() since that's when all the
    // remaining state variables are allocated.
    virtual void setMobilizerDefaultInstanceValues    (const SBModelVars&,     SBInstanceVars&)     const {}
    virtual void setMobilizerDefaultTimeValues        (const SBModelVars&,     SBTimeVars&)         const {}
    virtual void setMobilizerDefaultPositionValues    (const SBModelVars&,     Vector& q)           const {}
    virtual void setMobilizerDefaultVelocityValues    (const SBModelVars&,     Vector& u)           const {}
    virtual void setMobilizerDefaultDynamicsValues    (const SBModelVars&,     SBDynamicsVars&)     const {}
    virtual void setMobilizerDefaultAccelerationValues(const SBModelVars&,     SBAccelerationVars&) const {}

    // These attempt to set the mobilizer's internal configuration or velocity
    // to a specified value. This is intended to be a fast, local calculation that produces
    // an answer to machine precision *if* the mobilizer can represent the given
    // quantity exactly. The answer is returned in the appropriate slots of a caller-provided
    // "q-like" or "u-like" Vector; the implementation must not look at or change any
    // other slots.
    // It is OK for the implementation to use the current values of the coordinates or speeds,
    // and it is required to preserve any of these that are not needed to satisfy the 
    // request.
    // If the mobilizer can't satisfy the request to machine precision it should just
    // do nothing or the best it can, with the only general rule being that it shouldn't
    // make things worse. In particular, it does not need to work hard on an approximate solution.

    virtual void setMobilizerTransform
       (const SBModelVars&, const Transform& X_MbM, Vector& q) const = 0;
    virtual void setMobilizerRotation
       (const SBModelVars&, const Rotation& R_MbM, Vector& q) const = 0;
    virtual void setMobilizerTranslation
       (const SBModelVars&, const Vec3& T_MbM, Vector& q,
        bool dontChangeOrientation)                           const = 0;

    virtual void setMobilizerVelocity
       (const SBModelVars&, const Vector& q, const SpatialVec& V_MbM, Vector& u) const = 0;
    virtual void setMobilizerAngularVelocity
       (const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u)       const = 0;
    virtual void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector& q, const Vec3& v_MbM, Vector& u,
        bool dontChangeAngularVelocity)                                          const = 0;

    /// Calculate kinetic energy (from spatial quantities only).
    Real calcKineticEnergy(
        const SBPositionCache& pc,
        const SBVelocityCache& vc) const;   
  
    /// Calculate all spatial configuration quantities, assuming availability of
    /// joint-specific relative quantities.
    void calcJointIndependentKinematicsPos(
        SBPositionCache& pc) const;

    /// Calcluate all spatial velocity quantities, assuming availability of
    /// joint-specific relative quantities and all position kinematics.
    void calcJointIndependentKinematicsVel(
        const SBPositionCache& pc,
        SBVelocityCache&       vc) const;

    /// Calculate velocity-dependent quantities which will be needed for
    // computing accelerations.
    void calcJointIndependentDynamicsVel(
        const SBPositionCache& pc,
        const SBVelocityCache& vc,
        SBDynamicsCache&       dc) const;

    virtual const char* type()     const {return "unknown";}
    virtual int  getDOF()   const=0; //number of independent dofs
    virtual int  getMaxNQ() const=0; //dofs plus quaternion constraints
    virtual int  getNQ(const SBModelVars&) const=0; //actual number of q's

    // This depends on the mobilizer type and modeling options. If it returns
    // true then this node will be assigned a unique quaternion index.
    virtual bool isUsingQuaternion(const SBModelVars&) const=0;

    // This will do nothing unless the mobilizer is using a quaternion.
    virtual bool enforceQuaternionConstraints(
        const SBModelVars& mv,
        Vector&            q) const=0;

    virtual void calcArticulatedBodyInertiasInward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const=0;

    virtual void calcZ(
        const SBPositionCache&    pc,
        const SBDynamicsCache&    dc,
        const SpatialVec&         spatialForce,
        const SBAccelerationVars& av,
        SBAccelerationCache&      ac) const 
      { throw VirtualBaseMethod(); }

    virtual void calcYOutward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const                     
      { throw VirtualBaseMethod(); }

    virtual void calcAccel(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        const SBDynamicsCache& dc,
        SBAccelerationCache&   ac,
        Vector&                udot,
        Vector&                qdotdot) const 
      { throw VirtualBaseMethod(); }

    virtual void calcInternalGradientFromSpatial(
        const SBPositionCache&      pc, 
        Vector_<SpatialVec>&        zTmp,
        const Vector_<SpatialVec>&  X, 
        Vector&                     JX) const
      { throw VirtualBaseMethod(); }

    virtual void calcEquivalentJointForces(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector&                     jointForces) const
      { throw VirtualBaseMethod(); }

    virtual void calcUDotPass1Inward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               jointForces,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector_<SpatialVec>&        allGepsilon,
        Vector&                     allEpsilon) const
      { throw VirtualBaseMethod(); } 

    virtual void calcUDotPass2Outward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               epsilonTmp,
        Vector_<SpatialVec>&        allA_GB,
        Vector&                     allUDot) const
      { throw VirtualBaseMethod(); }

    virtual void calcMInverseFPass1Inward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               f,
        Vector_<SpatialVec>&        allZ,
        Vector_<SpatialVec>&        allGepsilon,
        Vector&                     allEpsilon) const
      { throw VirtualBaseMethod(); } 

    virtual void calcMInverseFPass2Outward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               epsilonTmp,
        Vector_<SpatialVec>&        allA_GB,
        Vector&                     allUDot) const
      { throw VirtualBaseMethod(); }

    virtual void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        Vector&                qdot) const
      { throw VirtualBaseMethod(); }

    virtual void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const
      { throw VirtualBaseMethod(); }

    virtual void setVelFromSVel(const SBPositionCache&, const SBVelocityCache&,
                                const SpatialVec&, Vector& u) const {throw VirtualBaseMethod();}

    // Copy the right q's from qIn to the corresponding slots in q. The number copied
    // may depend on modeling choices as supplied in the first argument.
    virtual void copyQ(
        const SBModelVars& mv, 
        const Vector&      qIn, 
        Vector&            q) const
      { throw VirtualBaseMethod(); }

    // Copy the right u's from uIn to the corresponding slots in u. Modeling choices
    // shouldn't affect the number copied since the number of u's should be identical
    // to the number of mobilities. TODO: then why the first argument? and why does this
    // need to be virtual?
    virtual void copyU(
        const SBModelVars& mv, 
        const Vector&      uIn, 
        Vector&            u) const
      { throw VirtualBaseMethod(); }

    virtual void getInternalForce(
        const SBAccelerationCache& ac, 
        Vector&                    tau) const {throw VirtualBaseMethod();}

    // Note that this requires rows of H to be packed like SpatialRow.
    virtual const SpatialRow& getHRow(const SBPositionCache&, int i) const {throw VirtualBaseMethod();}

    void nodeDump(std::ostream&) const;
    virtual void nodeSpecDump(std::ostream& o) const { o<<"NODE SPEC type="<<type()<<std::endl; }

protected:
    /// This is the constructor for the abstract base type for use by the derived
    /// concrete types in their constructors.
    RigidBodyNode(const MassProperties& mProps_B,
                  const Transform&      xform_PMb,
                  const Transform&      xform_BM)
      : qIndex(-1), uIndex(-1), uSqIndex(-1), quaternionIndex(-1),
        parent(0), children(), level(-1), nodeNum(-1),
        massProps_B(mProps_B), inertia_CB_B(mProps_B.calcCentralInertia()),
        X_BM(xform_BM), X_PMb(xform_PMb), refX_PB(xform_PMb*~xform_BM), X_MB(~xform_BM)
    {
    }

    typedef std::vector<RigidBodyNode*> RigidBodyNodeList;

    int qIndex;   // index into internal coord pos array
    int uIndex;   // index into internal coord vel,acc arrays
    int uSqIndex; // index into array of DOF^2 objects
    int quaternionIndex; // if this mobilizer has a quaternion, this is our slot

    RigidBodyNode*    parent; 
    RigidBodyNodeList children;
    int               level;        //how far from base 
    int               nodeNum;      //unique ID number in SimbodyMatterSubsystemRep

    // These are the default body properties, all supplied or calculated on
    // construction. TODO: they should be 
    // (optionally?) overrideable by Instance-level state variable entries.

    /// This is the mass, center of mass, and inertia as supplied at construction.
    /// Here the inertia is taken about the B origin OB.
    const MassProperties massProps_B;

    /// This is the supplied inertia, shifted to the center of mass. It is still
    /// a constant expressed in B, but is taken about the COM.
    const Inertia   inertia_CB_B;

    /// Orientation and location of inboard mobilizer frame M, measured
    /// and expressed in body frame B.
    const Transform X_BM; 
    const Transform X_MB; // inverse of X_BM, calculated on construction

    /// This is set when we attach this node to its parent in the tree. This is the
    /// configuration of the parent's outboard mobilizer attachment frame corresponding
    /// to body B (Mb) measured from and expressed in the parent frame P. It is 
    /// a constant in frame P. TODO: make it parameterizable.
    const Transform X_PMb;

    /// Reference configuration. This is the body frame B as measured and expressed in
    /// parent frame P *in the reference configuration*, that is, when B's mobilizer
    /// coordinates are all zero, so that Mb=M (see previous members). This is a constant
    /// after body B is attached to P in the tree: refX_PB = X_PMb * ~X_BM. 
    /// The body B frame can of course move relative to its parent, but that is not
    /// the meaning of this reference configuration transform.
    const Transform refX_PB;

    virtual void velFromCartesian() {}

    friend std::ostream& operator<<(std::ostream& s, const RigidBodyNode&);
    template<int dof> friend class RigidBodyNodeSpec;

};

#endif // SimTK_RIGID_BODY_NODE_H_
