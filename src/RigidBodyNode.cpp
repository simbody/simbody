/**@file
 * This file contains all the multibody mechanics code that involves a single body and
 * its inboard joint, that is, one node in the multibody tree.
 *
 * Most methods here expect to be called in a particular order during traversal of the
 * tree -- either base to tip or tip to base.
 */

#include "RigidBodyTree.h"
#include "RigidBodyNode.h"

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;

//////////////////////////////////////////////
// Implementation of RigidBodyNode methods. //
//////////////////////////////////////////////

void RigidBodyNode::addChild(RigidBodyNode* child) {
    children.push_back( child );
}

//
// Calc posCM, mass, Mk
//      phi, inertia
// Should be calc'd from base to tip.
// We depend on transforms X_PB and X_GB being available.
void RigidBodyNode::calcJointIndependentKinematicsPos(const State& s) const
{
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 T_PB_G = getX_GP(s).R() * getX_PB(s).T();

    // The Phi matrix conveniently performs child-to-parent shifting
    // on spatial quantities (forces); its transpose does parent-to-child
    // shifting for velocities.
    updPhi(s) = PhiMatrix(T_PB_G);

    // Calculate spatial mass properties. That means we need to transform
    // the local mass moments into the Ground frame and reconstruct the
    // spatial inertia matrix Mk.

    updInertia_OB_G(s) = getInertia_OB_B(s).changeAxes(~getX_GB(s).R());
    updCB_G(s)         = getX_GB(s).R()*getCOM_B(s);

    updCOM_G(s) = getX_GB(s).T() + getCB_G(s);

    // Calc Mk: the spatial inertia matrix about the body origin.
    // Note that this is symmetric; offDiag is *skew* symmetric so
    // that transpose(offDiag) = -offDiag.
    // Note: we need to calculate this now so that we'll be able to calculate
    // kinetic energy without going past the Motion stage.
    const Mat33 offDiag = getMass(s)*crossMat(getCB_G(s));
    updMk(s) = SpatialMat( getInertia_OB_G(s).toMat33() ,     offDiag ,
                                   -offDiag             , getMass(s)*Mat33(1) );
}

// Calculate velocity-related quantities: spatial velocity (sVel). This must be
// called base to tip: depends on parent's spatial velocity, and
// the just-calculated cross-joint spatial velocity V_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel(const State& s) const
{
    updV_GB(s) = ~getPhi(s)*parent->getV_GB(s) + getV_PB_G(s);
}

Real RigidBodyNode::calcKineticEnergy(const State& s) const {
    const Real ret = dot(getV_GB(s) , getMk(s)*getV_GB(s));
    return 0.5*ret;
}

// Calculate velocity-related quantities that are needed for building
// our dynamics operators, namely the gyroscopic force and coriolis acceleration.
// This routine expects that all spatial velocities & spatial inertias are
// already available, but does not have to be called in any particular order.
void 
RigidBodyNode::calcJointIndependentDynamicsVel(const State& s) const
{
    if (nodeNum == 0) { // ground, just in case
        updGyroscopicForce(s)      = SpatialVec(Vec3(0), Vec3(0));
        updCoriolisAcceleration(s) = SpatialVec(Vec3(0), Vec3(0));
        updCentrifugalForces(s)    = SpatialVec(Vec3(0), Vec3(0));
        return;
    }

    const Vec3& omega = getV_GB(s)[0];  // spatial angular velocity
    const Vec3& vel   = getV_GB(s)[1];  // spatial linear velocity

    updGyroscopicForce(s) = 
        SpatialVec(    omega % (getInertia_OB_G(s)*omega),    // gyroscopic moment
                   getMass(s)*(omega % (omega % getCB_G(s)))); // gyroscopic force

    // Parent velocity.
    const Vec3& pOmega = parent->getV_GB(s)[0];
    const Vec3& pVel   = parent->getV_GB(s)[1];

    // Calc a: coriolis acceleration.
    // Sherm TODO: See Schwieters & Clore Eq. [16], which uses *this*
    // body's omega (w_k) while the code below uses the *parent* pOmega
    // (w_k-1). Compare with Jain, Vaidehi, & Rodriguez 1991, Eq. 4.4,
    // esp. the following paragraph saying that the cross product
    // will be the same whether we use omega or pOmega in the 2nd term,
    // because they can only differ along H which is constant between P & B.
    // (caution: JV&R number backwards so the parent is w_k+1 there).
    // Is that also true in the first term? I.e., can we use pOmega, omega,
    // or both below and get the same answers? Anyway the code below is
    // consistent with JV&R but not obviously consistent with S&C.
    updCoriolisAcceleration(s) = 
        SpatialVec(Vec3(0), pOmega % (vel-pVel)) 
         // + crossMat(pOmega) * getV_PB_G(s); <-- IVM original
        + crossMat(omega)  * getV_PB_G(s); // JV&R paper

    updCentrifugalForces(s) =
        getP(s) * getCoriolisAcceleration(s) + getGyroscopicForce(s);
}

void RigidBodyNode::nodeDump(std::ostream& o) const {
    o << "NODE DUMP level=" << level << " type=" << type() << std::endl;
    nodeSpecDump(o);
    o << "END OF NODE type=" << type() << std::endl;
}

std::ostream& operator<<(std::ostream& o, const RigidBodyNode& node) {
    node.nodeDump(o);
    return o;
}


////////////////////////////////////////////////
// Define classes derived from RigidBodyNode. //
////////////////////////////////////////////////

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class RBGroundBody : public RigidBodyNode {
public:
    RBGroundBody() // TODO: should set mass properties to infinity
      : RigidBodyNode(MassProperties(),Transform(),Transform()) {}
    ~RBGroundBody() {}

    /*virtual*/const char* type() const { return "ground"; }
    /*virtual*/int getDOF() const { return 0; }
    /*virtual*/int getMaxNQ() const { return 0; }
    /*virtual*/int getNQ(const State&) const { return 0; }

    /*virtual*/void calcZ(const State&, const SpatialVec&) const {} 
    /*virtual*/void calcYOutward(const State&) const {}
    /*virtual*/void calcInternalForce(const State&, const SpatialVec&) const {}
    /*virtual*/void calcAccel(const State&) const {}

    /*virtual*/void realizeModeling(const State&) const {}
    /*virtual*/void realizeParameters(const State&) const {}
    /*virtual*/void realizeConfiguration(const State&) const {}
    /*virtual*/void realizeMotion(const State&) const {}
    /*virtual*/void setVelFromSVel(State&,const SpatialVec&) {}
    /*virtual*/bool enforceQuaternionConstraints(State&) const {return false;}
    
    /*virtual*/void calcArticulatedBodyInertiasInward(const State&) const {}

    /*virtual*/ void calcInternalGradientFromSpatial
        (const State&, Vector_<SpatialVec>&,
            const Vector_<SpatialVec>&, Vector&) const { }

    /*virtual*/ void calcEquivalentJointForces(const State& s,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    jointForces) const 
    { 
        allZ[0] = bodyForces[0];
        allGepsilon[0] = SpatialVec(Vec3(0), Vec3(0));
    }

    /*virtual*/void calcUDotPass1Inward(const State& s,
        const Vector&              jointForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    allEpsilon) const
    {
        allZ[0] = -bodyForces[0]; // TODO sign is weird
        allGepsilon[0] = SpatialVec(Vec3(0), Vec3(0));
    } 
    /*virtual*/void calcUDotPass2Outward(const State& s,
        const Vector&                   epsilonTmp,
        Vector_<SpatialVec>&            allA_GB,
        Vector&                         allUDot) const
    {
        allA_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    }

    /*virtual*/void setDefaultModelingValues(const State& s, 
                                             SBModelingVars& v) const
    {
        v.prescribed[0] = true; // ground's motion is prescribed to zero
    }

    /*virtual*/void getAccel(const State&, Vector&) const {}

    // /*virtual*/ const SpatialRow& getHRow(int i) const;

    void print(const State&, int) {}
};

template<int dof>
class RigidBodyNodeSpec : public RigidBodyNode {
public:
    RigidBodyNodeSpec(const MassProperties& mProps_B,
                      const Transform&   X_PJb,
                      const Transform&   X_BJ,
                      int&                  nextUSlot,
                      int&                  nextUSqSlot,
                      int&                  nextQSlot)
      : RigidBodyNode(mProps_B, X_PJb, X_BJ)
    {
        uIndex   = nextUSlot;   nextUSlot   += getDOF();
        uSqIndex = nextUSqSlot; nextUSqSlot += getDOF()*getDOF();
        qIndex   = nextQSlot;   nextQSlot   += getMaxNQ();
    }

    virtual ~RigidBodyNodeSpec() {}

    // This is the type of the joint transition matrix H.
    typedef Mat<dof,2,Row3,1,2> HType;


    // The following routines calculate joint-specific position kinematic
    // quantities. They may assume that *all* position kinematics (not just
    // joint-specific) has been done for the parent, and that the position
    // state variables q are available. Each routine may also assume that the
    // previous routines have been called, in the order below.
    // The routines are structured as operators -- they use the State but
    // do not change anything in it, including the cache. Instead they are
    // passed argument to write their results into. In practice, these
    // arguments will typically be in the State cache (see below).

    /// This mandatory routine performs expensive floating point operations sin,cos,sqrt
    /// in one place so we don't end up repeating them. sin&cos are used only
    /// for joints which have angular coordinates, and qnorm is only for joints
    /// which are using quaternions. Other joints can provide a null routine.
    /// Each of the passed-in Vectors is a "q-like" object, that is, allocated
    /// to the bodies in a manner parallel to the q state variable.
    virtual void calcJointSinCosQNorm
        (const State&, Vector& sine, Vector& cosine, Vector& qnorm) const=0;

    /// This mandatory routine calculates the across-joint transform X_JbJ generated
    /// by the current q values. This may depend on sines & cosines or normalized
    /// quaternions already being available in the State cache.
    virtual void calcAcrossJointTransform(const State&, Transform& X_JbJ) const=0;

    /// This routine is NOT joint specific, but cannot be called until the across-joint
    /// transform X_JbJ has been calculated and is available in the State cache.
    void calcBodyTransforms(const State& s, Transform& X_PB, Transform& X_GB) const 
    {
        const Transform& X_BJ  = getX_BJ(s);  // fixed
        const Transform& X_PJb = getX_PJb(s); // fixed
        const Transform& X_JbJ = getX_JbJ(s); // just calculated
        const Transform& X_GP  = getX_GP(s);  // already calculated

        X_PB = X_PJb * X_JbJ * ~X_BJ; // TODO: precalculate X_JB
        X_GB = X_GP * X_PB;
    }

    /// This mandatory routine calcluates the joint transition matrix H, giving the
    /// change of *spatial* velocity induced by the generalized speeds u for this
    /// joint. It may depend on X_PB and X_GB having been calculated already.
    virtual void calcJointTransitionMatrix(const State&, HType& H) const=0;

    /// Calculate joint-specific kinematic quantities dependent on
    /// on velocities. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for this node,
    /// that all velocity kinematics has been done for the parent, and
    /// that the velocity state variables (u) are available. The
    /// quanitites that must be computed are:
    ///   V_PB_G  relative velocity of B in P, expr. in G
    /// The code is the same for all joints, although parametrized by dof.
    void calcJointKinematicsVel(const State& s) const {
        updV_PB_G(s) = ~getH(s) * getU(s);
    }

    // These next two routines are options, but if you supply one you
    // must supply the other. (That is, ball-containing joints provide
    // both of these routines.)
    virtual void calcQDot(const State& s, const Vector& u, Vector& qdot) const {
        toQ(qdot) = fromU(u);        // default is qdot=u
    }

    virtual void calcQDotDot(const State& s, const Vector& udot, Vector& qdotdot) const {
        toQ(qdotdot) = fromU(udot);  // default is qdotdot=udot
    }

    void realizeModeling(const SimTK::State&) const {
    }

    void realizeParameters(const SimTK::State&) const {
    }

    // Set a new configuration and calculate the consequent kinematics.
    // Must call base-to-tip.
    void realizeConfiguration(const State& s, SBConfigurationCache& configCache) const {
        calcJointSinCosQNorm     (s, configCache.sq,
                                     configCache.cq,
                                     configCache.qnorm);
        calcAcrossJointTransform (s, updX_JbJ(s));
        calcBodyTransforms       (s, updX_PB(s), updX_GB(s));
        calcJointTransitionMatrix(s, updH(s));

        calcJointIndependentKinematicsPos(s, configCache);
    }

    // Set new velocities for the current configuration, and calculate
    // all the velocity-dependent terms. Must call base-to-tip.
    void realizeMotion(const State& s, const Vector& u, 
                       Vector& qdot, SBMotionCache& motionCache) const {
        calcQDot(s, u, qdot);
        calcJointKinematicsVel(s, motionCache);
        calcJointIndependentKinematicsVel(s, motionCache);
    }

    // This is a dynamics-stage calculation and must be called tip-to-base (inward).
    void calcArticulatedBodyInertiasInward(const State& s) const;

    // calcJointIndependentDynamicsVel() must be called after ArticulatedBodyInertias.

    // This dynamics-stage calculation is needed for handling constraints. It
    // must be called base-to-tip (outward);
    void calcYOutward(const State& s) const;

    // These routines give each node a chance to set appropriate defaults in a piece
    // of the state corresponding to a particular stage. Default implementations here
    // assume non-ball joint; override if necessary.
    virtual void setDefaultModelingValues (const State&, SBModelingVars&)  const {}
    virtual void setDefaultParameterValues(const State&, SBParameterVars&) const {}
    virtual void setDefaultTimeValues     (const State&, SBTimeVars&)      const {}

    virtual void setDefaultConfigurationValues(const State& s, Vector& q) const 
    {
        toQ(q) = 0.;
    }
    virtual void setDefaultMotionValues(const State&, Vector& u) const 
    {
        toU(u) = 0.;
    }
    virtual void setDefaultDynamicsValues(const State&, SBDynamicsVars&) const {}
    virtual void setDefaultReactionValues(const State&, 
                                          SBReactionVars& v) const
    {
        toB(v.appliedBodyForces) = SpatialVec(Vec3(0), Vec3(0));
        toU(v.appliedJointForces) = 0.;
        toU(v.prescribedUdot) = 0.;
    }

    // setQ and setU extract this node's values from the supplied
    // q-sized or u-sized array and put them in the corresponding
    // locations in the State. Joints which need quaternions should
    // override setQ to copy the extra q.
    virtual void setQ(State& s, const Vector& q) const {
        updQ(s) = fromQ(q);
    }

    virtual void setU(State& s, const Vector& u) const {
        updU(s) = fromU(u);
    }

    int          getDOF()            const { return dof; }
    virtual int  getMaxNQ()          const { return dof; } // maxNQ can be larger than dof
    virtual int  getNQ(const State&) const { return dof; } // DOF <= NQ <= maxNQ

    virtual void print(const State&, int) const;

    virtual void setVelFromSVel(State& s, const SpatialVec&) const;
    virtual bool enforceQuaternionConstraints(State&) const {return false;}

    const SpatialRow& getHRow(const State& s, int i) const {
        return getH(s)[i];
    }

    // Access to body-oriented state and cache entries is the same for all nodes,
    // and joint oriented access is almost the same but parametrized by dof. There is a special
    // case for quaternions because they use an extra state variable, and although we don't
    // have to we make special scalar routines available for 1-dof joints. Note that all State access
    // routines are inline, not virtual, so the cost is just an indirection and an index.

    // General joint-dependent select-my-goodies-from-the-pool routines.
    const Vec<dof>&     fromQ  (const Vector& q)   const {return Vec<dof>::getAs(&q[qIndex]);}
    Vec<dof>&           toQ    (      Vector& q)   const {return Vec<dof>::updAs(&q[qIndex]);}
    const Vec<dof>&     fromU  (const Vector& u)   const {return Vec<dof>::getAs(&u[uIndex]);}
    Vec<dof>&           toU    (      Vector& u)   const {return Vec<dof>::updAs(&u[uIndex]);}
    const Mat<dof,dof>& fromUSq(const Vector& uSq) const {return Mat<dof,dof>::getAs(&uSq[uSqIndex]);}
    Mat<dof,dof>&       toUSq  (      Vector& uSq) const {return Mat<dof,dof>::updAs(&uSq[uSqIndex]);}

    // Same, but specialized for the common case where dof=1 so everything is scalar.
    const Real& from1Q  (const Vector& q)   const {return q[qIndex];}
    Real&       to1Q    (      Vector& q)   const {return q[qIndex];}
    const Real& from1U  (const Vector& u)   const {return u[uIndex];}
    Real&       to1U    (      Vector& u)   const {return u[uIndex];}
    const Real& from1USq(const Vector& uSq) const {return uSq[uSqIndex];}
    Real&       to1USq  (      Vector& uSq) const {return uSq[uSqIndex];}

    // Same, specialized for quaternions. We're assuming that the quaternion comes first in the coordinates.
    const Vec4& fromQuat(const Vector& q) const {return Vec4::getAs(&q[qIndex]);}
    Vec4&       toQuat  (      Vector& q) const {return Vec4::updAs(&q[qIndex]);}

    // Extract a Vec3 from a Q-like or U-like object, beginning at an offset from the qIndex or uIndex.
    const Vec3& fromQVec3(const Vector& q, int offs) const {return Vec3::getAs(&q[qIndex+offs]);}
    Vec3&       toQVec3  (      Vector& q, int offs) const {return Vec3::updAs(&q[qIndex+offs]);}
    const Vec3& fromUVec3(const Vector& u, int offs) const {return Vec3::getAs(&u[uIndex+offs]);}
    Vec3&       toUVec3  (      Vector& u, int offs) const {return Vec3::updAs(&u[uIndex+offs]);}

    // Applications of the above extraction routines to particular interesting items in the State. Note
    // that you can't use these for quaternions since they extract "dof" items.

    // State variables (read only).
    const Vec<dof>&   getQ             (const State& s) const {return fromQ(s.configVars.q);}
    const Vec<dof>&   getU             (const State& s) const {return fromU(s.motionVars.u);}
    const Vec<dof>&   getAppliedJointForce(const State& s) const {return fromU(s.reactionVars.appliedJointForces);}
    const Vec<dof>&   getPrescribedUdot   (const State& s) const {return fromU(s.reactionVars.prescribedUdot);}

    // Special case state access for 1-dof joints
    const Real& get1Q             (const State& s) const {return from1Q(s.configVars.q);}
    const Real& get1U             (const State& s) const {return from1U(s.motionVars.u);}
    const Real& get1AppliedJointForce(const State& s) const {return from1U(s.reactionVars.appliedJointForces);}
    const Real& get1PrescribedUdot   (const State& s) const {return from1U(s.reactionVars.prescribedUdot);}

    // Special case for quaternions and Vec3 at offset.
    const Vec4& getQuat (const State& s)           const {return fromQuat(s.configVars.q);}
    const Vec3& getQVec3(const State& s, int offs) const {return fromQVec3(s.configVars.q, offs);}
    const Vec3& getUVec3(const State& s, int offs) const {return fromUVec3(s.motionVars.u, offs);}

    // State variables for updating; be careful. This is only appropriate for "solvers", such as
    // a method which modifies state variables to satisfy constraints.
    Vec<dof>& updQ (State& s) const {return toQ(s.configVars.q);}
    Vec<dof>& updU (State& s) const {return toU(s.motionVars.u);} 
    Real&     upd1Q(State& s) const {return to1Q(s.configVars.q);}
    Real&     upd1U(State& s) const {return to1U(s.motionVars.u);} 

    // Special case for quaternions and Vec3 at offset.
    Vec4& updQuat (State& s)           const {return toQuat(s.configVars.q);}
    Vec3& updQVec3(State& s, int offs) const {return toQVec3(s.configVars.q, offs);}
    Vec3& updUVec3(State& s, int offs) const {return toUVec3(s.motionVars.u, offs);}

    // Cache entries (cache is mutable in a const State)

        // Configuration

    // TODO: should store as H or else always reference Ht
    const Mat<dof,2,Row3,1,2>& getH(const State& s) const
      { return ~Mat<2,dof,Vec3>::getAs(&s.configCache.storageForHt(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH(const State& s) const
      { return ~Mat<2,dof,Vec3>::updAs(&s.configCache.storageForHt(0,uIndex)); }

    // These are sines and cosines of angular qs. The rest of the slots are garbage.
    const Vec<dof>&   getSinQ   (const State& s) const {return fromQ (s.configCache.sq);}
    Vec<dof>&         updSinQ   (const State& s) const {return toQ   (s.configCache.sq);}
    const Real&       get1SinQ  (const State& s) const {return from1Q(s.configCache.sq);}
    Real&             upd1SinQ  (const State& s) const {return to1Q  (s.configCache.sq);}

    const Vec<dof>&   getCosQ   (const State& s) const {return fromQ (s.configCache.cq);}
    Vec<dof>&         updCosQ   (const State& s) const {return toQ   (s.configCache.cq);}
    const Real&       get1CosQ  (const State& s) const {return from1Q(s.configCache.cq);}
    Real&             upd1CosQ  (const State& s) const {return to1Q  (s.configCache.cq);}

    // These are normalized quaternions in slots for balls. Everything else is garbage.
    const Vec4&       getQNorm  (const State& s) const {return fromQuat(s.configCache.qnorm);}
    Vec4&             updQNorm  (const State& s) const {return toQuat  (s.configCache.qnorm);}

        // Motion

    const Vec<dof>&   getQDot   (const State& s) const {return fromQ (s.motionCache.qdot);}
    Vec<dof>&         updQDot   (const State& s) const {return toQ   (s.motionCache.qdot);}
    const Real&       get1QDot  (const State& s) const {return from1Q(s.motionCache.qdot);}
    Real&             upd1QDot  (const State& s) const {return to1Q  (s.motionCache.qdot);}

        // Dynamics
    const Mat<dof,dof>& getD(const State& s) const {return fromUSq(s.dynamicsCache.storageForD);}
    Mat<dof,dof>&       updD(const State& s) const {return toUSq  (s.dynamicsCache.storageForD);}

    const Mat<dof,dof>& getDI(const State& s) const {return fromUSq(s.dynamicsCache.storageForDI);}
    Mat<dof,dof>&       updDI(const State& s) const {return toUSq  (s.dynamicsCache.storageForDI);}

    const Mat<2,dof,Vec3>& getG(const State& s) const
      { return Mat<2,dof,Vec3>::getAs(&s.dynamicsCache.storageForG(0,uIndex)); }
    Mat<2,dof,Vec3>&       updG(const State& s) const
      { return Mat<2,dof,Vec3>::updAs(&s.dynamicsCache.storageForG(0,uIndex)); }

        // Reaction
    const Vec<dof>&   getUDot   (const State& s) const {return fromU (s.reactionCache.udot);}
    Vec<dof>&         updUDot   (const State& s) const {return toU   (s.reactionCache.udot);}
    const Real&       get1UDot  (const State& s) const {return from1U(s.reactionCache.udot);}
    Real&             upd1UDot  (const State& s) const {return to1U  (s.reactionCache.udot);}

    const Vec<dof>&   getQDotDot (const State& s) const {return fromQ (&s.reactionCache.qdotdot);}
    Vec<dof>&         updQDotDot (const State& s) const {return toQ   (&s.reactionCache.qdotdot);}
    const Real&       get1QDotDot(const State& s) const {return from1Q(s.reactionCache.qdotdot);}
    Real&             upd1QDotDot(const State& s) const {return to1Q  (s.reactionCache.qdotdot);}

    const Vec<dof>&   getNetHingeForce (const State& s) const {return fromU (s.reactionCache.netHingeForces);}
    Vec<dof>&         updNetHingeForce (const State& s) const {return toU   (s.reactionCache.netHingeForces);}
    const Real&       get1NetHingeForce (const State& s) const {return from1U(s.reactionCache.netHingeForces);}
    Real&             upd1NetHingeForce(const State& s) const {return to1U  (s.reactionCache.netHingeForces);}


    const Vec<dof>&   getNu (const State& s) const {return fromU (s.reactionCache.nu);}
    Vec<dof>&         updNu (const State& s) const {return toU   (s.reactionCache.nu);}
    const Real&       get1Nu(const State& s) const {return from1U(s.reactionCache.nu);}
    Real&             upd1Nu(const State& s) const {return to1U  (s.reactionCache.nu);}

    const Vec<dof>&   getEpsilon (const State& s) const {return fromU (s.reactionCache.epsilon);}
    Vec<dof>&         updEpsilon (const State& s) const {return toU   (s.reactionCache.epsilon);}
    const Real&       get1Epsilon(const State& s) const {return from1U(s.reactionCache.epsilon);}
    Real&             upd1Epsilon(const State& s) const {return to1U  (s.reactionCache.epsilon);}

    void calcZ(const State& s, const SpatialVec& spatialForce) const;

    void calcAccel(const State& s) const;
    void calcInternalGradientFromSpatial(const State&, Vector_<SpatialVec>& zTmp,
                                         const Vector_<SpatialVec>& X, Vector& JX) const;

    void calcEquivalentJointForces(const State& s,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    jointForces) const;

    void calcUDotPass1Inward(const State& s,
        const Vector&              jointForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    allEpsilon) const; 
    void calcUDotPass2Outward(const State& s,
        const Vector&                   epsilonTmp,
        Vector_<SpatialVec>&            allA_GB,
        Vector&                         allUDot) const;

    void nodeSpecDump(std::ostream& o, const State& s) const {
        o << "uIndex=" << uIndex << " mass=" << getMass() 
            << " COM_G=" << getCOM_G(s) << std::endl;
        o << "inertia_OB_G=" << getInertia_OB_G(s) << std::endl;
        o << "H=" << getH(s) << std::endl;
        o << "SVel=" << getV_GB(s) << std::endl;
        o << "a=" << getCoriolisAcceleration(s) << std::endl;
        o << "b=" << getGyroscopicForce(s) << std::endl;
        o << "Th  =" << getQ(s) << std::endl;
        o << "dTh =" << getU(s) << std::endl;
        o << "ddTh=" << getUDot(s) << std::endl;
        o << "SAcc=" << getA_GB(s) << std::endl;
    }
};

//////////////////////////////////////////
// Derived classes for each joint type. //
//////////////////////////////////////////


/**
 * Translate (Cartesian) joint. This provides three degrees of
 * translational freedom which is suitable (e.g.) for connecting a
 * free atom to ground. The Cartesian directions are the axes of
 * the parent body's Jb frame, with J=Jb when all 3 coords are 0,
 * and the orientation of J in Jb is 0 forever.
 */
class RBNodeTranslate : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "translate"; }

    RBNodeTranslate(const MassProperties& mProps_B,
                    const Transform&   X_PJb,
                    const Transform&   X_BJ,
                    int&                  nextUSlot,
                    int&                  nextUSqSlot,
                    int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

        // Implementations of virtual methods.

    // This is required but does nothing here since we there are no rotations for this joint.
    void calcJointSinCosQNorm
        (const State&, Vector& sine, Vector& cosine, Vector& qnorm) const { }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const State& s, Transform& X_JbJ) const {
        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A Cartesian joint can't change orientation. 
        X_JbJ = Transform(RotationMat(), getQ(s));
    }

    // Calculate H.
    void calcJointTransitionMatrix(const State& s, HType& H) const {
        const Transform& X_PJb   = getX_PJb(s);      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const Transform& X_GP    = getX_GP(s); // parent orientation in ground

        // Note that H is spatial. The current spatial directions for our qs are
        // the axes of the Jb frame expressed in Ground.
        const RotationMat R_GJb = X_GP.R()*X_PJb.R();
        H[0] = SpatialRow( Row3(0), ~R_GJb.x() );
        H[1] = SpatialRow( Row3(0), ~R_GJb.y() );
        H[2] = SpatialRow( Row3(0), ~R_GJb.z() );
    }
};



/**
 * Sliding joint (1 dof translation). The translation is along the z
 * axis of the parent body's Jb frame, with J=Jb when the coordinate
 * is zero and the orientation of J in Jb frozen at 0 forever.
 */
class RBNodeSlider : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "slider"; }

    RBNodeSlider(const MassProperties& mProps_B,
                 const Transform&   X_PJb,
                 const Transform&   X_BJ,
                 int&                  nextUSlot,
                 int&                  nextUSqSlot,
                 int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }
        // Implementations of virtual methods.

    // This is required but does nothing here since we there are no rotations for this joint.
    void calcJointSinCosQNorm
        (const State&, Vector& sine, Vector& cosine, Vector& qnorm) const { }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const State& s, Transform& X_JbJ) const {
        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A sliding joint can't change orientation, and only translates along z. 
        X_JbJ = Transform(RotationMat(), Vec3(0.,0.,get1Q(s)));
    }

    // Calculate H.
    void calcJointTransitionMatrix(const State& s, HType& H) const {
        const Transform& X_PJb   = getX_PJb(s);      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const Transform& X_GP    = getX_GP(s); // parent configuration in ground

        // Note that H is spatial. The current spatial directions for our q is
        // the z axis of the Jb frame expressed in Ground.
        const Vec3 z_GJb = X_GP.R()*X_PJb.z();
        H[0] = SpatialRow( Row3(0), ~z_GJb );
    }
};

/**
 * This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
 * about a particular axis, the z axis of the parent's Jb frame, which is 
 * aligned forever with the z axis of the body's J frame. In addition, the
 * origin points of J and Jb are identical forever.
 */
class RBNodeTorsion : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "torsion"; }

    RBNodeTorsion(const MassProperties& mProps_B,
                  const Transform&   X_PJb,
                  const Transform&   X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const State& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        const Real& q = get1Q(s); // angular coordinate
        to1Q(sine)    = std::sin(q);
        to1Q(cosine)  = std::cos(q);
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const State& s, Transform& X_JbJ) const {
        const Real& q  = get1Q(s);    // angular coordinate

        // We're only updating the orientation here because a torsion joint
        // can't translate (it is defined as a rotation about the z axis).
        X_JbJ.updR().setToRotationAboutZ(q);
        X_JbJ.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(const State& s, HType& H) const {
        const Transform& X_BJ  = getX_BJ(s);  // fixed
        const Transform& X_PJb = getX_PJb(s); // fixed
        const Transform& X_GP  = getX_GP(s);  // calculated earlier
        const Transform& X_GB  = getX_GB(s);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // Calc H matrix in space-fixed coords.
        // This works because the joint z axis is the same in J & Jb
        // since that's what we rotate around.
        const Vec3 z_G = X_GP.R() * X_PJb.z();
        H[0] = SpatialRow( ~z_G, ~(z_G % T_JB_G) );
    }
};

/**
 * U-joint like joint type which allows rotation about the two axes
 * perpendicular to zDir. This is appropriate for diatoms and for allowing 
 * torsion+bond angle bending.
 */
class RBNodeRotate2 : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "rotate2"; }

    RBNodeRotate2(const MassProperties& mProps_B,
                  const Transform&   X_PJb,
                  const Transform&   X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const State& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        const Vec2& q = getQ(s); // angular coordinates
        toQ(sine)   = Vec2(std::sin(q[0]), std::sin(q[1]));
        toQ(cosine) = Vec2(std::cos(q[0]), std::cos(q[1]));
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const State& s, Transform& X_JbJ) const {
        const Vec2& q  = getQ(s); // angular coordinates

        // We're only updating the orientation here because a U-joint
        // can't translate.
        X_JbJ.updR().setToSpaceFixed12(q);
        X_JbJ.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(const State& s, HType& H) const {
        const Transform& X_BJ  = getX_BJ(s);  // fixed
        const Transform& X_PJb = getX_PJb(s); // fixed
        const Transform& X_GP  = getX_GP(s);  // calculated earlier
        const Transform& X_GB  = getX_GB(s);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // The coordinates are defined in the space-fixed (that is, Jb) frame, so
        // the orientation of Jb in ground gives the instantaneous spatial 
        // meaning of the coordinates.
        const RotationMat R_GJb = X_GP.R() * X_PJb.R();
        H[0] = SpatialRow(~R_GJb.x(), ~(R_GJb.x() % T_JB_G));
        H[1] = SpatialRow(~R_GJb.y(), ~(R_GJb.y() % T_JB_G));
    }
};

/**
 * The "diatom" joint is the equivalent of a free joint for a body with no inertia in
 * one direction, such as one composed of just two atoms. It allows unrestricted
 * translation but rotation only about directions perpendicular to the body's
 * inertialess axis.
 * The coordinate definitions are a combination of a rotate2 joint and a
 * Cartesian joint. The first 2 are rotational, the next 3 are translations.
 * However, the rotations don't affect the translations.
 */
class RBNodeTranslateRotate2 : public RigidBodyNodeSpec<5> {
public:
    virtual const char* type() { return "diatom"; }

    RBNodeTranslateRotate2(const MassProperties& mProps_B,
                           const Transform&   X_PJb,
                           const Transform&   X_BJ,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const State& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        const Vec2& q = getQ(s).getSubVec<2>(0); // angular coordinates
        toQ(sine).updSubVec<2>(0)   = Vec2(std::sin(q[0]), std::sin(q[1]));
        toQ(cosine).updSubVec<2>(0) = Vec2(std::cos(q[0]), std::cos(q[1]));
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const State& s, Transform& X_JbJ) const {
        const Vec<5>& q = getQ(s);     // joint coordinates

        X_JbJ.updR().setToSpaceFixed12(q.getSubVec<2>(0));
        X_JbJ.updT() = q.getSubVec<3>(2);
    }

    // Calculate H.
    void calcJointTransitionMatrix(const State& s, HType& H) const {
        const Transform& X_BJ  = getX_BJ(s);  // fixed
        const Transform& X_PJb = getX_PJb(s); // fixed
        const Transform& X_GP  = getX_GP(s);  // calculated earlier
        const Transform& X_GB  = getX_GB(s);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // The rotational coordinates are defined in the space-fixed 
        // (that is, Jb) frame, so the orientation of Jb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        const RotationMat R_GJb = X_GP.R() * X_PJb.R();
        H[0] = SpatialRow(~R_GJb.x(), ~(R_GJb.x() % T_JB_G));
        H[1] = SpatialRow(~R_GJb.y(), ~(R_GJb.y() % T_JB_G));
        H[2] = SpatialRow(  Row3(0) ,     ~R_GJb.x());
        H[3] = SpatialRow(  Row3(0) ,     ~R_GJb.y());
        H[4] = SpatialRow(  Row3(0) ,     ~R_GJb.z());    
    }
};

/// Ball joint. This provides three degrees of rotational freedom, i.e.,
/// unrestricted orientation of the body's J frame in the parent's Jb frame.
/// The 3 u's are the cross-joint angular velocity vector of J in Jb, and
/// udot's are the angular velocity time derivative. The q's, however, are
/// either 3 Euler angles in a 3-2-1 body-fixed sequence, or 4 quaternions.
/// In that case we calculate either 3 or 4 qdots from the u's.
class RBNodeRotate3 : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "rotate3"; }

    RBNodeRotate3(const MassProperties& mProps_B,
                  const Transform&   X_PJb,
                  const Transform&   X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const State& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        if (getUseEulerAngles(s)) {
            const Vec3& q = getQ(s); // angular coordinates
            toQ(sine)   = Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
            toQ(cosine) = Vec3(std::cos(q[0]), std::cos(q[1]), std::cos(q[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& q = getQuat(s); // unnormalized quaternion from state
            toQuat(qnorm) = q / q.norm();
        }
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const State& s, Transform& X_JbJ) const {
        X_JbJ.updT() = 0.; // This joint can't translate.
        if (getUseEulerAngles(s))
            X_JbJ.updR().setToBodyFixed123(getQ(s));
        else
            X_JbJ.updR().setToQuaternion(getQuat(s));
    }

    // Calculate H.
    void calcJointTransitionMatrix(const State& s, HType& H) const {
        const Transform& X_BJ  = getX_BJ(s);  // fixed
        const Transform& X_PJb = getX_PJb(s); // fixed
        const Transform& X_GP  = getX_GP(s);  // calculated earlier
        const Transform& X_GB  = getX_GB(s);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // The rotational coordinates are defined in the space-fixed 
        // (that is, Jb) frame, so the orientation of Jb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        const RotationMat R_GJb = X_GP.R() * X_PJb.R();
        H[0] = SpatialRow(~R_GJb.x(), ~(R_GJb.x() % T_JB_G));
        H[1] = SpatialRow(~R_GJb.y(), ~(R_GJb.y() % T_JB_G));
        H[2] = SpatialRow(~R_GJb.z(), ~(R_GJb.z() % T_JB_G));
    }

    void calcQDot(const State& s, const Vector& u, Vector& qdot) const {
        const Vec3& w_JbJ = fromU(u); // angular velocity of J in Jb 
        if (getUseEulerAngles(s)) {
            const RotationMat& R_JbJ = getX_JbJ(s).R();
            toQ(qdot) = RotationMat::convertAngVelToBodyFixed123Dot(getQ(s),
                                        ~R_JbJ*w_JbJ); // need w in *body*, not parent
        } else
            toQuat(qdot) = RotationMat::convertAngVelToQuaternionDot(getQuat(s),w_JbJ);
    }
 
    void calcQDotDot(const State& s, const Vector& udot, Vector& qdotdot) const {
        const Vec3& w_JbJ     = getU(s); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_JbJ_dot = fromU(udot);

        if (getUseEulerAngles(s)) {
            const RotationMat& R_JbJ = getX_JbJ(s).R();
            toQ(qdotdot)    = RotationMat::convertAngVelDotToBodyFixed123DotDot
                                  (getQ(s), ~R_JbJ*w_JbJ, ~R_JbJ*w_JbJ_dot);
        } else
            toQuat(qdotdot) = RotationMat::convertAngVelDotToQuaternionDotDot
                                  (getQuat(s),w_JbJ,w_JbJ_dot);
    }

    void setQ(State& s, const Vector& q) const {
        if (getUseEulerAngles(s))
            updQ(s) = fromQ(q);
        else
            updQuat(s) = fromQuat(q);
    }

    int getMaxNQ()              const {return 4;}
    int getNQ(const State& s) const {return getUseEulerAngles(s) ? 3 : 4;} 

    void setDefaultConfigurationValues(const State& s, 
                                       SBConfigurationVars& v) const 
    {
        if (getUseEulerAngles(s)) toQ(v.q) = 0.;
        else toQuat(v.q) = Vec4(1.,0.,0.,0.);
    }

    bool enforceQuaternionConstraints(State& s) const {
        if (getUseEulerAngles(s)) 
            return false;   // no change
        Vec4& quat = updQuat(s);
        quat = quat / quat.norm();
        return true;
    }

    void getInternalForce(const State& s) const {
        assert(false); // TODO: decompose cross-joint torque into 123 gimbal torques
        /* OLD BALL CODE:
        Vector& f = s.cache->netHingeForces;
        //dependency: calcR_PB must be called first
        assert( useEuler );

        const Vec<3,Vec2>& scq = getSinCosQ(s);
        const Real sPhi   = scq[0][0], cPhi   = scq[0][1];
        const Real sTheta = scq[1][0], cTheta = scq[1][1];
        const Real sPsi   = scq[2][0], cPsi   = scq[2][1];

        Vec3 torque = forceInternal;
        const Mat33 M( 0.          , 0.          , 1.    ,
                      -sPhi        , cPhi        , 0.    ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta );
        Vec3 eTorque = RigidBodyNode::DEG2RAD * M * torque;

        Vec3::updAs(&v[uIndex]) = eTorque;
        */
    }
};

/// Free joint. This provides six degrees of freedom, three rotational and
/// three translational. The rotation is like the ball joint above; the
/// translation is like the Cartesian joint above.
class RBNodeTranslateRotate3 : public RigidBodyNodeSpec<6> {
public:
    virtual const char* type() { return "full"; }

    RBNodeTranslateRotate3(const MassProperties& mProps_B,
                           const Transform&   X_PJb,
                           const Transform&   X_BJ,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<6>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const State& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        if (getUseEulerAngles(s)) {
            const Vec3& q = getQ(s).getSubVec<3>(0); // angular coordinates
            toQ(sine).updSubVec<3>(0)   = Vec3(std::sin(q[0]), std::sin(q[1]), std::sin(q[2]));
            toQ(cosine).updSubVec<3>(0) = Vec3(std::cos(q[0]), std::cos(q[1]), std::cos(q[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& q = getQuat(s); // unnormalized quaternion from state
            toQuat(qnorm) = q / q.norm();
        }
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const State& s, Transform& X_JbJ) const {
        if (getUseEulerAngles(s)) {
            X_JbJ.updR().setToBodyFixed123(getQVec3(s,0));
            X_JbJ.updT() = getQVec3(s,3);
        } else {
            X_JbJ.updR().setToQuaternion(getQuat(s));
            X_JbJ.updT() = getQVec3(s,4);
        }
    }

    // Calculate H.
    void calcJointTransitionMatrix(const State& s, HType& H) const {
        const Transform& X_BJ  = getX_BJ(s);  // fixed
        const Transform& X_PJb = getX_PJb(s); // fixed
        const Transform& X_GP  = getX_GP(s);  // calculated earlier
        const Transform& X_GB  = getX_GB(s);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // The rotational speeds are defined in the space-fixed 
        // (that is, Jb) frame, so the orientation of Jb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        const RotationMat R_GJb = X_GP.R() * X_PJb.R();
        H[0] = SpatialRow(~R_GJb.x(), ~(R_GJb.x() % T_JB_G));
        H[1] = SpatialRow(~R_GJb.y(), ~(R_GJb.y() % T_JB_G));
        H[2] = SpatialRow(~R_GJb.z(), ~(R_GJb.z() % T_JB_G));
        H[3] = SpatialRow(  Row3(0) ,     ~R_GJb.x());
        H[4] = SpatialRow(  Row3(0) ,     ~R_GJb.y());
        H[5] = SpatialRow(  Row3(0) ,     ~R_GJb.z());

        //cout << "T_JB_G=" << T_JB_G << endl;
        //cout << "H=" << H;
    }

    void calcQDot(const State& s, const Vector& u, Vector& qdot) const {
        const Vec3& w_JbJ = fromUVec3(u,0); // Angular velocity
        const Vec3& v_JbJ = fromUVec3(u,3); // Linear velocity
        if (getUseEulerAngles(s)) {
            const RotationMat& R_JbJ = getX_JbJ(s).R();
            const Vec3& theta = getQVec3(s,0); // Euler angles
            toQVec3(qdot,0) = RotationMat::convertAngVelToBodyFixed123Dot(theta,
                                            ~R_JbJ*w_JbJ); // need w in *body*, not parent
            toQVec3(qdot,3) = v_JbJ;
            //cout << "EulerAngles: " << theta << endl;
            //cout << "   w_JbJ=" << w_JbJ << "  v_JbJ=" << v_JbJ << endl;
            //cout << "   qdot=" << fromQVec3(qdot,0) << endl;
        } else {
            const Vec4& quat = getQuat(s);
            toQuat (qdot)   = RotationMat::convertAngVelToQuaternionDot(quat,w_JbJ);
            toQVec3(qdot,4) = v_JbJ;
        }
    }
 
    void calcQDotDot(const State& s, const Vector& udot, Vector& qdotdot) const {
        const Vec3& w_JbJ     = getUVec3(s,0); // angular velocity of J in Jb
        const Vec3& v_JbJ     = getUVec3(s,3); // linear velocity
        const Vec3& w_JbJ_dot = fromUVec3(udot,0);
        const Vec3& v_JbJ_dot = fromUVec3(udot,3);
        if (getUseEulerAngles(s)) {
            const RotationMat& R_JbJ = getX_JbJ(s).R();
            const Vec3& theta  = getQVec3(s,0); // Euler angles
            toQVec3(qdotdot,0) = RotationMat::convertAngVelDotToBodyFixed123DotDot
                                                (theta, ~R_JbJ*w_JbJ, ~R_JbJ*w_JbJ_dot);
            toQVec3(qdotdot,3) = v_JbJ_dot;
            //cout << "   w_JbJ_dot=" << w_JbJ_dot << "  v_JbJ_dot=" << v_JbJ_dot << endl;
            //cout << "   qdotdot=" << fromQVec3(qdotdot,0) << endl;
        } else {
            const Vec4& quat  = getQuat(s);
            toQuat(qdotdot)   = RotationMat::convertAngVelDotToQuaternionDotDot
                                                (quat,w_JbJ,w_JbJ_dot);
            toQVec3(qdotdot,4) = v_JbJ_dot;
        }
    }

    void setQ(State& s, const Vector& q) const {
        if (getUseEulerAngles(s))
            updQ(s) = fromQ(q);
        else {
            updQuat(s)    = fromQuat(q);
            updQVec3(s,4) = fromQVec3(q,4);
        }
    }

    int getMaxNQ()              const {return 7;}
    int getNQ(const State& s) const {return getUseEulerAngles(s) ? 6 : 7;} 

    void setDefaultConfigurationValues(const State& s, 
                                       SBConfigurationVars& v) const 
    {
        if (getUseEulerAngles(s)) 
            toQ(v.q) = 0.;
        else {
            toQuat(v.q) = Vec4(1.,0.,0.,0.);
            toQVec3(v.q,4) = 0.;
        }
    }

    bool enforceQuaternionConstraints(State& s) const {
        if (getUseEulerAngles(s)) 
            return false; // no change
        Vec4& quat = updQuat(s);
        quat = quat / quat.norm();
        return true;
    }

    void getInternalForce(const State& s) const {
        assert(false); // TODO: decompose cross-joint torque into 123 gimbal torques
        /* OLD BALL CODE:
        Vector& f = s.cache->netHingeForces;
        //dependency: calcR_PB must be called first
        assert( useEuler );

        const Vec<3,Vec2>& scq = getSinCosQ(s);
        const Real sPhi   = scq[0][0], cPhi   = scq[0][1];
        const Real sTheta = scq[1][0], cTheta = scq[1][1];
        const Real sPsi   = scq[2][0], cPsi   = scq[2][1];

        Vec3 torque = forceInternal;
        const Mat33 M( 0.          , 0.          , 1.    ,
                      -sPhi        , cPhi        , 0.    ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta );
        Vec3 eTorque = RigidBodyNode::DEG2RAD * M * torque;

        Vec3::updAs(&v[uIndex]) = eTorque;
        */
    }
};

////////////////////////////////////////////////
// RigidBodyNode factory based on joint type. //
////////////////////////////////////////////////

/*static*/ RigidBodyNode*
RigidBodyNode::create(
    const MassProperties& m,            // mass properties in body frame
    const Transform&   X_PJb,        // parent's attachment frame for this joint
    const Transform&   X_BJ,         // inboard joint frame J in body frame
    JointSpecification::JointType      
                          type,
    bool                  isReversed,   // child-to-parent orientation?
    int&                  nxtUSlot,
    int&                  nxtUSqSlot,
    int&                  nxtQSlot)  
{
    assert(!isReversed);

    switch(type) {
    case JointSpecification::ThisIsGround:
        return new RBGroundBody();
    case JointSpecification::Torsion:
        return new RBNodeTorsion(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case JointSpecification::Universal:        
        return new RBNodeRotate2(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case JointSpecification::Orientation:
        return new RBNodeRotate3(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case JointSpecification::Cartesian:
        return new RBNodeTranslate(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case JointSpecification::FreeLine:
        return new RBNodeTranslateRotate2(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case JointSpecification::Free:
        return new RBNodeTranslateRotate3(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case JointSpecification::Sliding:
        return new RBNodeSlider(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case JointSpecification::Cylinder:
    case JointSpecification::Planar:
    case JointSpecification::Gimbal:
    case JointSpecification::Weld:

    default: 
        assert(false);
    };

    return 0;
}

/////////////////////////////////////////////////////////////
// Implementation of RigidBodyNodeSpec base class methods. //
/////////////////////////////////////////////////////////////


//
// to be called from base to tip.
//
template<int dof> void
RigidBodyNodeSpec<dof>::setVelFromSVel(State& s, 
                                       const SpatialVec& sVel) const 
{
    updU(s) = getH(s) * (sVel - (~getPhi(s) * parent->getV_GB(s)));
}

//
// Given only position-related quantities from the State 
//      Mk  (this body's spatial inertia matrix)
//      Phi (composite body child-to-parent shift matrix)
//      H   (joint transition matrix)
// we calculate dynamic quantities 
//      P   (articulated body inertia)
//      D   (factored mass matrix LDL' diagonal part D=H*P*~H)
//      DI  (inverse of D)
//      G   (P * ~H * DI)
//   tauBar (I-G*H, a temporary not reused elsewhere)
//      Psi (articulated body child-to-parent shift matrix)
// and put them in the state cache.
// This must be called tip-to-base (inward).
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcArticulatedBodyInertiasInward(
    const State&    s) const 
{
    updP(s) = getMk(s);
    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialMat& tauBarChild = children[i]->getTauBar(s);
        const SpatialMat& PChild      = children[i]->getP(s);
        const PhiMatrix&  phiChild    = children[i]->getPhi(s);

        // TODO: this is around 450 flops but could be cut in half by
        // exploiting symmetry.
        updP(s) += phiChild * (tauBarChild * PChild) * ~phiChild;
    }

    const Mat<2,dof,Vec3> PHt = getP(s) * ~getH(s);
    updD(s)  = getH(s) * PHt;
    // this will throw an exception if the matrix is ill conditioned
    updDI(s) = getD(s).invert();
    updG(s)  = PHt * getDI(s);

    // TODO: change sign on tau to make it GH-I instead, which only requires
    // subtractions on the diagonal rather than negating all the off-diag stuff.
    // That would save 30 flops here (I know, not much).
    updTauBar(s)  = 1.; // identity matrix
    updTauBar(s) -= getG(s) * getH(s);
    updPsi(s)     = getPhi(s) * getTauBar(s);
}



// To be called base to tip.
template<int dof> void
RigidBodyNodeSpec<dof>::calcYOutward(const State& s) const {
    // TODO: this is very expensive (~1000 flops?) Could cut be at least half
    // by exploiting symmetry. Also, does Psi have special structure?
    // And does this need to be computed for every body or only those
    // which are loop "base" bodies or some such?
    updY(s) = (~getH(s) * getDI(s) * getH(s)) 
                + (~getPsi(s) * parent->getY(s) * getPsi(s));
}

//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcZ(const State& s, 
                              const SpatialVec& spatialForce) const 
{
    SpatialVec& z = updZ(s);
    z = getCentrifugalForces(s) - spatialForce;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild    = children[i]->getZ(s);
        const PhiMatrix&  phiChild  = children[i]->getPhi(s);
        const SpatialVec& GepsChild = children[i]->getGepsilon(s);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(s)  = getAppliedJointForce(s) - getH(s)*z; // TODO: pass in hinge forces
    updNu(s)       = getDI(s) * getEpsilon(s);
    updGepsilon(s) = getG(s)  * getEpsilon(s);
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were fed to calcZ (as embodied in 'nu').
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcAccel(const State& s) const {
    Vec<dof>&        udot   = updUDot(s);
    const SpatialVec alphap = ~getPhi(s) * parent->getA_GB(s); // ground A_GB is 0


    udot       = getNu(s) - (~getG(s)*alphap);
    updA_GB(s) = alphap + ~getH(s)*udot + getCoriolisAcceleration(s);  

    calcQDotDot(s, s.reactionCache.udot, s.reactionCache.qdotdot);   
}

 
//
// To be called from tip to base.
// Temps do not need to be initialized.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcUDotPass1Inward(const State& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector_<SpatialVec>&       allZ,
    Vector_<SpatialVec>&       allGepsilon,
    Vector&                    allEpsilon) const 
{
    const Vec<dof>&   myJointForce = fromU(jointForces);
    const SpatialVec& myBodyForce  = fromB(bodyForces);
    SpatialVec&       z            = toB(allZ);
    SpatialVec&       Geps         = toB(allGepsilon);
    Vec<dof>&         eps          = toU(allEpsilon);

    z = getCentrifugalForces(s) - myBodyForce;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(s);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = myJointForce - getH(s)*z;
    Geps = getG(s)  * eps;
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were reduced into epsilon (e.g., see above).
// Base to tip: temp allA_GB does not need to be initialized before
// beginning the iteration.
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcUDotPass2Outward(const State& s,
    const Vector&                   allEpsilon,
    Vector_<SpatialVec>&            allA_GB,
    Vector&                         allUDot) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = toB(allA_GB);
    Vec<dof>&       udot = toU(allUDot); // pull out this node's udot

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = parent->getNodeNum()== 0 
        ? SpatialVec(Vec3(0), Vec3(0))
        : ~getPhi(s) * allA_GB[parent->getNodeNum()];

    udot = getDI(s) * eps - (~getG(s)*A_GP);
    A_GB = A_GP + ~getH(s)*udot + getCoriolisAcceleration(s);  
}

//
// Calculate product of partial velocities J and a gradient vector on each of the
// outboard bodies. This is to be called tip to base. Requires that Phi and H are available, so this
// should only be called in Stage::Configured or higher. This does not change the cache at all.
// NOTE (sherm 060214): I reworked this from the original. This one no longer incorporates
// applied hinge gradients if there are any; just add those in at the end if you want them.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcInternalGradientFromSpatial
    (const State& s, Vector_<SpatialVec>& zTmp,
     const Vector_<SpatialVec>& X, Vector& JX) const
{
    const SpatialVec& in  = X[getNodeNum()];
    Vec<dof>&         out = Vec<dof>::updAs(&JX[getUIndex()]);
    SpatialVec&       z   = zTmp[getNodeNum()];

    z = in;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
        const PhiMatrix&  phiChild = children[i]->getPhi(s);

        z += phiChild * zChild;
    }

    out = getH(s) * z; 
}

//
// To be called from tip to base.
// Temps do not need to be initialized.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcEquivalentJointForces(const State& s,
    const Vector_<SpatialVec>& bodyForces,
    Vector_<SpatialVec>&       allZ,
    Vector_<SpatialVec>&       allGepsilon,
    Vector&                    jointForces) const 
{
    const SpatialVec& myBodyForce  = fromB(bodyForces);
    SpatialVec&       z            = toB(allZ);
    SpatialVec&       Geps         = toB(allGepsilon);
    Vec<dof>&         eps          = toU(jointForces);

    z = myBodyForce;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(s);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = getH(s) * z;
    Geps = getG(s) * eps;
}



template<int dof> void
RigidBodyNodeSpec<dof>::print(const State& s, int verbose) const {
    if (verbose&InternalDynamics::printNodePos) 
        cout << setprecision(8)
             << ": pos: " << getX_GB(s).T() << ' ' << '\n';
    if (verbose&InternalDynamics::printNodeTheta) 
        cout << setprecision(8)
             << ": theta: " 
             << getQ(s) << ' ' << getU(s)  << ' ' << getUDot(s)  << '\n';
}

