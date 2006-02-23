/**@file
 * This file contains all the multibody mechanics code that involves a single body and
 * its inboard joint, that is, one node in the multibody tree.
 *
 * Most methods here expect to be called in a particular order during traversal of the
 * tree -- either base to tip or tip to base.
 */

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
void RigidBodyNode::calcJointIndependentKinematicsPos(const SBState& s) const
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
    const Mat33 offDiag = getMass(s)*crossMat(getCB_G(s));
    updMk(s) = SpatialMat( getInertia_OB_G(s).toMat33() ,     offDiag ,
                                   -offDiag             , getMass(s)*Mat33(1) );
}

// Calculate velocity-related quantities: spatial velocity (sVel),
// gyroscopic force b, coriolis acceleration a. This must be
// called base to tip: depends on parent's sVel, V_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel(const SBState& s) const
{
    updV_GB(s) = ~getPhi(s)*parent->getV_GB(s) + getV_PB_G(s);
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
          + crossMat(pOmega) * getV_PB_G(s);
    //    + crossMat(omega)  * getV_PB_G(s);  <-- should work too?
}

Real RigidBodyNode::calcKineticEnergy(const SBState& s) const {
    const Real ret = dot(getV_GB(s) , getMk(s)*getV_GB(s));
    return 0.5*ret;
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
      : RigidBodyNode(MassProperties(),TransformMat(),TransformMat()) {}
    ~RBGroundBody() {}

    /*virtual*/const char* type() const { return "ground"; }
    /*virtual*/int getDOF() const { return 0; }
    /*virtual*/int getMaxNQ() const { return 0; }
    /*virtual*/int getNQ(const SBState&) const { return 0; }


    /*virtual*/void calcP(const SBState&) const {} 
    /*virtual*/void calcZ(const SBState&, const SpatialVec&) const {} 
    /*virtual*/void calcY(const SBState&) const {}
    /*virtual*/void calcInternalForce(const SBState&, const SpatialVec&) const {}
    /*virtual*/void calcAccel(const SBState&) const {}

    /*virtual*/void realizeModeling(const SBState&) const {}
    /*virtual*/void realizeParameters(const SBState&) const {}
    /*virtual*/void realizeConfiguration(const SBState&) const {}
    /*virtual*/void realizeMotion(const SBState&) const {}
    /*virtual*/void setVelFromSVel(SBState&,const SpatialVec&) {}
    /*virtual*/void enforceQuaternionConstraints(SBState&) const {}

    /*virtual*/void getDefaultParameters(SBState&)    const {}
    /*virtual*/void getDefaultConfiguration(SBState&) const {}
    /*virtual*/void getDefaultVelocity(SBState&)      const {}

    /*virtual*/void getAccel(const SBState&, Vector&) const {}

    /*virtual*/void getInternalForce(const SBState&, Vector&) const {}
    // /*virtual*/ const SpatialRow& getHRow(int i) const;

    void print(const SBState&, int) {}
};

template<int dof>
class RigidBodyNodeSpec : public RigidBodyNode {
public:
    RigidBodyNodeSpec(const MassProperties& mProps_B,
                      const TransformMat&   X_PJb,
                      const TransformMat&   X_BJ,
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
        (const SBState&, Vector& sine, Vector& cosine, Vector& qnorm) const=0;

    /// This mandatory routine calculates the across-joint transform X_JbJ generated
    /// by the current q values. This may depend on sines & cosines or normalized
    /// quaternions already being available in the State cache.
    virtual void calcAcrossJointTransform(const SBState&, TransformMat& X_JbJ) const=0;

    /// This routine is NOT joint specific, but cannot be called until the across-joint
    /// transform X_JbJ has been calculated and is available in the State cache.
    void calcBodyTransforms(const SBState& s, TransformMat& X_PB, TransformMat& X_GB) const 
    {
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_JbJ = getX_JbJ(s); // just calculated
        const TransformMat& X_GP  = getX_GP(s);  // already calculated

        X_PB = X_PJb * X_JbJ * ~X_BJ; // TODO: precalculate X_JB
        X_GB = X_GP * X_PB;
    }

    /// This mandatory routine calcluates the joint transition matrix H, giving the
    /// change of *spatial* velocity induced by the generalized speeds u for this
    /// joint. It may depend on X_PB and X_GB having been calculated already.
    virtual void calcJointTransitionMatrix(const SBState&, HType& H) const=0;

    /// Calculate joint-specific kinematic quantities dependent on
    /// on velocities. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for this node,
    /// that all velocity kinematics has been done for the parent, and
    /// that the velocity state variables (u) are available. The
    /// quanitites that must be computed are:
    ///   V_PB_G  relative velocity of B in P, expr. in G
    /// The code is the same for all joints, although parametrized by dof.
    void calcJointKinematicsVel(const SBState& s) const {
        updV_PB_G(s) = ~getH(s) * getU(s);
    }

    void realizeModeling(const simtk::SBState&) const {
    }

    void realizeParameters(const simtk::SBState&) const {
    }

    /// Set a new configuration and calculate the consequent kinematics.
    void realizeConfiguration(const SBState& s) const {
        calcJointSinCosQNorm     (s, s.cache->sq,s.cache->cq,s.cache->qnorm);
        calcAcrossJointTransform (s, updX_JbJ(s));
        calcBodyTransforms       (s, updX_PB(s), updX_GB(s));
        calcJointTransitionMatrix(s, updH(s));

        calcJointIndependentKinematicsPos(s);
    }

    /// Set new velocities for the current configuration, and calculate
    /// all the velocity-dependent terms.
    void realizeMotion(const SBState& s) const {
        // anything to initialize?
        calcJointKinematicsVel(s);
        calcJointIndependentKinematicsVel(s);
    }

    virtual void calcJointAccel(const SBState&) const { }

    // We are assuming that the caller is taking care of state validity.
    virtual void getDefaultParameters(SBState& s) const {
        // TODO none yet
    }
    virtual void getDefaultConfiguration(SBState& s) const {
        updQ(s) = 0.;
    }
    virtual void getDefaultVelocity(SBState& s) const {
        updU(s) = 0.;
    }

    // setQ and setU extract this node's values from the supplied
    // q-sized or u-sized array and put them in the corresponding
    // locations in the State. Joints which need quaternions should
    // override setQ to copy the extra q.
    virtual void setQ(SBState& s, const Vector& q) const {
        updQ(s) = fromQ(q);
    }

    virtual void setU(SBState& s, const Vector& u) const {
        updU(s) = fromU(u);
    }

    int          getDOF()              const { return dof; }
    virtual int  getMaxNQ()            const { return dof; } // maxNQ can be larger than dof
    virtual int  getNQ(const SBState&) const { return dof; } // DOF <= NQ <= maxNQ

    virtual void print(const SBState&, int) const;

    virtual void setVelFromSVel(SBState& s, const SpatialVec&) const;
    virtual void enforceQuaternionConstraints(SBState&) const { }

    const SpatialRow& getHRow(const SBState& s, int i) const {
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
    const Vec<dof>&   getQ             (const SBState& s) const {return fromQ(s.vars->q);}
    const Vec<dof>&   getU             (const SBState& s) const {return fromU(s.vars->u);}
    const Vec<dof>&   getJointForce    (const SBState& s) const {return fromU(s.vars->appliedJointForces);}
    const Vec<dof>&   getPrescribedUdot(const SBState& s) const {return fromU(s.vars->prescribedUdot);}

    // Special case state access for 1-dof joints
    const Real& get1Q             (const SBState& s) const {return from1Q(s.vars->q);}
    const Real& get1U             (const SBState& s) const {return from1U(s.vars->u);}
    const Real& get1JointForce    (const SBState& s) const {return from1U(s.vars->appliedJointForces);}
    const Real& get1PrescribedUdot(const SBState& s) const {return from1U(s.vars->prescribedUdot);}

    // Special case for quaternions and Vec3 at offset.
    const Vec4& getQuat (const SBState& s)           const {return fromQuat(s.vars->q);}
    const Vec3& getQVec3(const SBState& s, int offs) const {return fromQVec3(s.vars->q, offs);}
    const Vec3& getUVec3(const SBState& s, int offs) const {return fromUVec3(s.vars->u, offs);}

    // State variables for updating; be careful. This is only appropriate for "solvers", such as
    // a method which modifies state variables to satisfy constraints.
    Vec<dof>& updQ (SBState& s) const {return toQ(s.vars->q);}
    Vec<dof>& updU (SBState& s) const {return toU(s.vars->u);} 
    Real&     upd1Q(SBState& s) const {return to1Q(s.vars->q);}
    Real&     upd1U(SBState& s) const {return to1U(s.vars->u);} 

    // Special case for quaternions and Vec3 at offset.
    Vec4& updQuat (SBState& s)           const {return toQuat(s.vars->q);}
    Vec3& updQVec3(SBState& s, int offs) const {return toQVec3(s.vars->q, offs);}
    Vec3& updUVec3(SBState& s, int offs) const {return toUVec3(s.vars->u, offs);}

    // Cache entries (cache is mutable in a const State)

        // Configuration

    // TODO: should store as H or else always reference Ht
    const Mat<dof,2,Row3,1,2>& getH(const SBState& s) const
      { return ~Mat<2,dof,Vec3>::getAs(&s.cache->storageForHt(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH(const SBState& s) const
      { return ~Mat<2,dof,Vec3>::updAs(&s.cache->storageForHt(0,uIndex)); }

    // These are sines and cosines of angular qs. The rest of the slots are garbage.
    const Vec<dof>&   getSinQ   (const SBState& s) const {return fromQ (s.cache->sq);}
    Vec<dof>&         updSinQ   (const SBState& s) const {return toQ   (s.cache->sq);}
    const Real&       get1SinQ  (const SBState& s) const {return from1Q(s.cache->sq);}
    Real&             upd1SinQ  (const SBState& s) const {return to1Q  (s.cache->sq);}

    const Vec<dof>&   getCosQ   (const SBState& s) const {return fromQ (s.cache->cq);}
    Vec<dof>&         updCosQ   (const SBState& s) const {return toQ   (s.cache->cq);}
    const Real&       get1CosQ  (const SBState& s) const {return from1Q(s.cache->cq);}
    Real&             upd1CosQ  (const SBState& s) const {return to1Q  (s.cache->cq);}

    // These are normalized quaternions in slots for balls. Everything else is garbage.
    const Vec4&       getQNorm  (const SBState& s) const {return fromQuat(s.cache->qnorm);}
    Vec4&             updQNorm  (const SBState& s) const {return toQuat  (s.cache->qnorm);}

        // Motion

    const Vec<dof>&   getQdot   (const SBState& s) const {return fromQ (s.cache->qdot);}
    Vec<dof>&         updQdot   (const SBState& s) const {return toQ   (s.cache->qdot);}
    const Real&       get1Qdot  (const SBState& s) const {return from1Q(s.cache->qdot);}
    Real&             upd1Qdot  (const SBState& s) const {return to1Q  (s.cache->qdot);}

        // Dynamics
    const Vec<dof>&   getUdot   (const SBState& s) const {return fromU (s.cache->udot);}
    Vec<dof>&         updUdot   (const SBState& s) const {return toU   (s.cache->udot);}
    const Real&       get1Udot  (const SBState& s) const {return from1U(s.cache->udot);}
    Real&             upd1Udot  (const SBState& s) const {return to1U  (s.cache->udot);}

    const Vec<dof>&   getQdotDot (const SBState& s) const {return fromQ (&s.cache->qdotdot);}
    Vec<dof>&         updQdotDot (const SBState& s) const {return toQ   (&s.cache->qdotdot);}
    const Real&       get1QdotDot(const SBState& s) const {return from1Q(s.cache->qdotdot);}
    Real&             upd1QdotDot(const SBState& s) const {return to1Q  (s.cache->qdotdot);}

    const Vec<dof>&   getInternalForce (const SBState& s) const {return fromU (s.cache->netHingeForces);}
    Vec<dof>&         updInternalForce (const SBState& s) const {return toU   (s.cache->netHingeForces);}
    const Real&       get1InternalForce(const SBState& s) const {return from1U(s.cache->netHingeForces);}
    Real&             upd1InternalForce(const SBState& s) const {return to1U  (s.cache->netHingeForces);}

    const Mat<dof,dof>& getDI(const SBState& s) const {return fromUSq(s.cache->storageForDI);}
    Mat<dof,dof>&       updDI(const SBState& s) const {return toUSq  (s.cache->storageForDI);}

    const Mat<2,dof,Vec3>& getG(const SBState& s) const
      { return Mat<2,dof,Vec3>::getAs(&s.cache->storageForG(0,uIndex)); }
    Mat<2,dof,Vec3>&       updG(const SBState& s) const
      { return Mat<2,dof,Vec3>::updAs(&s.cache->storageForG(0,uIndex)); }

    const Vec<dof>&   getNu (const SBState& s) const {return fromU (s.cache->nu);}
    Vec<dof>&         updNu (const SBState& s) const {return toU   (s.cache->nu);}
    const Real&       get1Nu(const SBState& s) const {return from1U(s.cache->nu);}
    Real&             upd1Nu(const SBState& s) const {return to1U  (s.cache->nu);}

    const Vec<dof>&   getEpsilon (const SBState& s) const {return fromU (s.cache->epsilon);}
    Vec<dof>&         updEpsilon (const SBState& s) const {return toU   (s.cache->epsilon);}
    const Real&       get1Epsilon(const SBState& s) const {return from1U(s.cache->epsilon);}
    Real&             upd1Epsilon(const SBState& s) const {return to1U  (s.cache->epsilon);}

    void calcP(const SBState& s) const;
    void calcZ(const SBState& s, const SpatialVec& spatialForce) const;
    void calcY(const SBState& s) const;
    void calcAccel(const SBState& s) const;
    void calcInternalGradientFromSpatial(const SBState&, Vector_<SpatialVec>& zTmp,
                                         const Vector_<SpatialVec>& X, Vector& JX) const;

    void nodeSpecDump(std::ostream& o, const SBState& s) const {
        o << "stateOffset=" << stateOffset << " mass=" << getMass() 
            << " COM_G=" << getCOM_G(s) << std::endl;
        o << "inertia_OB_G=" << getInertia_OB_G(s) << std::endl;
        o << "H=" << getH(s) << std::endl;
        o << "SVel=" << getV_GB(s) << std::endl;
        o << "a=" << getCoriolisAcceleration(s) << std::endl;
        o << "b=" << getGyroscopicForce(s) << std::endl;
        o << "Th  =" << getQ(s) << std::endl;
        o << "dTh =" << getU(s) << std::endl;
        o << "ddTh=" << getUdot(s) << std::endl;
        o << "SAcc=" << getA_GB(s) << std::endl;
    }
protected:
    /*
    // These are the joint-specific quantities
    //      ... position level
    Vec<dof>                theta;   // internal coordinates
    Mat<dof,2,Row3,1,2>     H;       // joint transition matrix (spatial); note row-order packing
                                     //   so that transpose will be a nice column-packed matrix
    Mat<dof,dof>            DI;
    Mat<2,dof,Vec3>         G;       // structure is like ~H; that is, default column-packed matrix

    //      ... velocity level
    Vec<dof>                dTheta;  // internal coordinate time derivatives

    //      ... acceleration level
    Vec<dof>                ddTheta; // - from the eq. of motion

    Vec<dof>                nu;
    Vec<dof>                epsilon;
    Vec<dof>                forceInternal;
    */
};

/*static*/const double RigidBodyNode::DEG2RAD = std::acos(-1.) / 180.; // i.e., pi/180
//const double RigidBodyNode::DEG2RAD = 1.0;  //always use radians


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
                    const TransformMat&   X_PJb,
                    const TransformMat&   X_BJ,
                    int&                  nextUSlot,
                    int&                  nextUSqSlot,
                    int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

        // Implementations of virtual methods.

    // This is required but does nothing here since we there are no rotations for this joint.
    void calcJointSinCosQNorm
        (const SBState&, Vector& sine, Vector& cosine, Vector& qnorm) const { }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const SBState& s, TransformMat& X_JbJ) const {
        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A Cartesian joint can't change orientation. 
        X_JbJ = TransformMat(RotationMat(), getQ(s));
    }

    // Calculate H.
    void calcJointTransitionMatrix(const SBState& s, HType& H) const {
        const TransformMat& X_PJb   = getX_PJb(s);      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const TransformMat& X_GP    = getX_GP(s); // parent orientation in ground

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
                 const TransformMat&   X_PJb,
                 const TransformMat&   X_BJ,
                 int&                  nextUSlot,
                 int&                  nextUSqSlot,
                 int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }
        // Implementations of virtual methods.

    // This is required but does nothing here since we there are no rotations for this joint.
    void calcJointSinCosQNorm
        (const SBState&, Vector& sine, Vector& cosine, Vector& qnorm) const { }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const SBState& s, TransformMat& X_JbJ) const {
        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A sliding joint can't change orientation, and only translates along z. 
        X_JbJ = TransformMat(RotationMat(), Vec3(0.,0.,get1Q(s)));
    }

    // Calculate H.
    void calcJointTransitionMatrix(const SBState& s, HType& H) const {
        const TransformMat& X_PJb   = getX_PJb(s);      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const TransformMat& X_GP    = getX_GP(s); // parent configuration in ground

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
                  const TransformMat&   X_PJb,
                  const TransformMat&   X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const SBState& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        const Real& q = get1Q(s); // angular coordinate
        to1Q(sine)    = std::sin(q);
        to1Q(cosine)  = std::cos(q);
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const SBState& s, TransformMat& X_JbJ) const {
        const Real& q  = get1Q(s);    // angular coordinate

        // We're only updating the orientation here because a torsion joint
        // can't translate (it is defined as a rotation about the z axis).
        X_JbJ.updR().setToRotationAboutZ(q);
        X_JbJ.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(const SBState& s, HType& H) const {
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // calculated earlier
        const TransformMat& X_GB  = getX_GB(s);  // just calculated

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
                  const TransformMat&   X_PJb,
                  const TransformMat&   X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const SBState& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        const Vec2& q = getQ(s); // angular coordinates
        toQ(sine)   = Vec2(std::sin(q[0]), std::sin(q[1]));
        toQ(cosine) = Vec2(std::cos(q[0]), std::cos(q[1]));
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const SBState& s, TransformMat& X_JbJ) const {
        const Vec2& q  = getQ(s); // angular coordinates

        // We're only updating the orientation here because a U-joint
        // can't translate.
        X_JbJ.updR().setToSpaceFixed12(q);
        X_JbJ.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(const SBState& s, HType& H) const {
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // calculated earlier
        const TransformMat& X_GB  = getX_GB(s);  // just calculated

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
                           const TransformMat&   X_PJb,
                           const TransformMat&   X_BJ,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const SBState& s, Vector& sine, Vector& cosine, Vector& qnorm) const
    {
        const Vec2& q = getQ(s).getSubVec<2>(0); // angular coordinates
        toQ(sine).updSubVec<2>(0)   = Vec2(std::sin(q[0]), std::sin(q[1]));
        toQ(cosine).updSubVec<2>(0) = Vec2(std::cos(q[0]), std::cos(q[1]));
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(const SBState& s, TransformMat& X_JbJ) const {
        const Vec<5>& q = getQ(s);     // joint coordinates

        X_JbJ.updR().setToSpaceFixed12(q.getSubVec<2>(0));
        X_JbJ.updT() = q.getSubVec<3>(2);
    }

    // Calculate H.
    void calcJointTransitionMatrix(const SBState& s, HType& H) const {
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // calculated earlier
        const TransformMat& X_GB  = getX_GB(s);  // just calculated

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
                  const TransformMat&   X_PJb,
                  const TransformMat&   X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const SBState& s, Vector& sine, Vector& cosine, Vector& qnorm) const
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
    void calcAcrossJointTransform(const SBState& s, TransformMat& X_JbJ) const {
        X_JbJ.updT() = 0.; // This joint can't translate.
        if (getUseEulerAngles(s))
            X_JbJ.updR().setToBodyFixed321(getQ(s));
        else
            X_JbJ.updR().setToQuaternion(getQuat(s));
    }

    // Calculate H.
    void calcJointTransitionMatrix(const SBState& s, HType& H) const {
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // calculated earlier
        const TransformMat& X_GB  = getX_GB(s);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // The rotational coordinates are defined in the space-fixed 
        // (that is, Jb) frame, so the orientation of Jb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        const RotationMat R_GJb = X_GP.R() * X_PJb.R();
        H[0] = SpatialRow(~R_GJb.x(), ~(R_GJb.x() % T_JB_G));
        H[1] = SpatialRow(~R_GJb.y(), ~(R_GJb.y() % T_JB_G));
        H[2] = SpatialRow(~R_GJb.z(), ~(R_GJb.z() % T_JB_G));
    }

    void calcQdot(const SBState& s, Vector& qdot) const {
        const Vec3& w_JbJ = getU(s); // angular velocity of J in Jb 
        if (getUseEulerAngles(s))
            toQ(qdot) = RotationMat::convertAngVelToBodyFixed321Dot(getQ(s),w_JbJ);
        else
            toQuat(qdot) = RotationMat::convertAngVelToQuaternionDot(getQuat(s),w_JbJ);
    }
 
    void calcQdotdot(const SBState& s, Vector& qdotdot) const {
        const Vec3& w_JbJ     = getU(s); // angular velocity of J in Jb
        const Vec3& w_JbJ_dot = getUdot(s);
        if (getUseEulerAngles(s))
            toQ(qdotdot)    = RotationMat::convertAngVelDotToBodyFixed321DotDot
                                  (getQ(s),w_JbJ,w_JbJ_dot);
        else
            toQuat(qdotdot) = RotationMat::convertAngVelDotToQuaternionDotDot
                                  (getQuat(s),w_JbJ,w_JbJ_dot);
    }

    void setQ(SBState& s, const Vector& q) const {
        if (getUseEulerAngles(s))
            updQ(s) = fromQ(q);
        else
            updQuat(s) = fromQuat(q);
    }

    int getMaxNQ()              const {return 4;}
    int getNQ(const SBState& s) const {return getUseEulerAngles(s) ? 3 : 4;} 

    void getDefaultConfiguration(SBState& s) const {
        if (getUseEulerAngles(s)) updQ(s) = 0.;
        else updQuat(s) = Vec4(1.,0.,0.,0.);
    }
    
    void getDefaultVelocity(SBState& s) const {
        updU(s) = 0; // no funny business here
    }

    void enforceQuaternionConstraints(SBState& s) {
        if (getUseEulerAngles(s)) return;
        Vec4& quat = updQuat(s);
        quat = quat / quat.norm();
    }

    void getInternalForce(const SBState& s) const {
        assert(false); // TODO: decompose cross-joint torque into 321 gimbal torques
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
                           const TransformMat&   X_PJb,
                           const TransformMat&   X_BJ,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<6>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm
        (const SBState& s, Vector& sine, Vector& cosine, Vector& qnorm) const
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
    void calcAcrossJointTransform(const SBState& s, TransformMat& X_JbJ) const {
        if (getUseEulerAngles(s)) {
            X_JbJ.updR().setToBodyFixed321(getQVec3(s,0));
            X_JbJ.updT() = getQVec3(s,3);
        } else {
            X_JbJ.updR().setToQuaternion(getQuat(s));
            X_JbJ.updT() = getQVec3(s,4);
        }
    }

    // Calculate H.
    void calcJointTransitionMatrix(const SBState& s, HType& H) const {
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // calculated earlier
        const TransformMat& X_GB  = getX_GB(s);  // just calculated

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
    }

    void calcQdot(const SBState& s, Vector& qdot) const {
        const Vec3& w_JbJ = getUVec3(s,0); // Angular velocity
        const Vec3& v_JbJ = getUVec3(s,3); // Linear velocity
        if (getUseEulerAngles(s)) {
            const Vec3& theta = getQVec3(s,0); // Euler angles
            toQVec3(qdot,0) = RotationMat::convertAngVelToBodyFixed321Dot(theta,w_JbJ);
            toQVec3(qdot,3) = v_JbJ;
        } else {
            const Vec4& quat = getQuat(s);
            toQuat (qdot)   = RotationMat::convertAngVelToQuaternionDot(quat,w_JbJ);
            toQVec3(qdot,4) = v_JbJ;
        }
    }
 
    void calcQdotdot(const SBState& s, Vector& qdotdot) const {
        const Vec3& w_JbJ     = getUVec3(s,0); // angular velocity of J in Jb
        const Vec3& v_JbJ     = getUVec3(s,3); // linear velocity
        const Vec3& w_JbJ_dot = getUdot(s).getSubVec<3>(0);
        const Vec3& v_JbJ_dot = getUdot(s).getSubVec<3>(3);
        if (getUseEulerAngles(s)) {
            const Vec3& theta  = getQVec3(s,0); // Euler angles
            toQVec3(qdotdot,0) = RotationMat::convertAngVelDotToBodyFixed321DotDot
                                                (theta,w_JbJ,w_JbJ_dot);
            toQVec3(qdotdot,3) = v_JbJ_dot;
        } else {
            const Vec4& quat  = getQuat(s);
            toQuat(qdotdot)   = RotationMat::convertAngVelDotToQuaternionDotDot
                                                (quat,w_JbJ,w_JbJ_dot);
            toQVec3(qdotdot,4) = v_JbJ_dot;
        }
    }

    void setQ(SBState& s, const Vector& q) const {
        if (getUseEulerAngles(s))
            updQ(s) = fromQ(q);
        else {
            updQuat(s)    = fromQuat(q);
            updQVec3(s,4) = fromQVec3(q,4);
        }
    }

    int getMaxNQ()              const {return 7;}
    int getNQ(const SBState& s) const {return getUseEulerAngles(s) ? 6 : 7;} 

    void getDefaultConfiguration(SBState& s) const {
        if (getUseEulerAngles(s)) 
            updQ(s) = 0.;
        else {
            updQuat(s)    = Vec4(1,0,0,0);
            updQVec3(s,4) = 0.;
        }
    }
    
    void getDefaultVelocity(SBState& s) const {
        updU(s) = 0; // no funny business here
    }

    void enforceQuaternionConstraints(SBState& s) {
        if (getUseEulerAngles(s)) return;
        Vec4& quat = updQuat(s);
        quat = quat / quat.norm();
    }

    void getInternalForce(const SBState& s) const {
        assert(false); // TODO: decompose cross-joint torque into 321 gimbal torques
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
    const TransformMat&   X_PJb,        // parent's attachment frame for this joint
    const TransformMat&   X_BJ,         // inboard joint frame J in body frame
    Joint::JointType      type,
    bool                  isReversed,   // child-to-parent orientation?
    int&                  nxtUSlot,
    int&                  nxtUSqSlot,
    int&                  nxtQSlot)  
{
    assert(!isReversed);

    switch(type) {
    case Joint::ThisIsGround:
        return new RBGroundBody();
    case Joint::Torsion:
        return new RBNodeTorsion(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Universal:        
        return new RBNodeRotate2(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Orientation:
        return new RBNodeRotate3(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Cartesian:
        return new RBNodeTranslate(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::FreeLine:
        return new RBNodeTranslateRotate2(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Free:
        return new RBNodeTranslateRotate3(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Sliding:
        return new RBNodeSlider(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Cylinder:
    case Joint::Planar:
    case Joint::Gimbal:
    case Joint::Weld:

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
RigidBodyNodeSpec<dof>::setVelFromSVel(SBState& s, 
                                       const SpatialVec& sVel) const 
{
    updU(s) = getH(s) * (sVel - (~getPhi(s) * parent->getV_GB(s)));
}

//
// Calculate Pk and related quantities. The requires that the children
// of the node have already had their quantities calculated, i.e. this
// is a tip to base recursion.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcP(const SBState& s) const {
    //
    //how much do we need to keep around?
    // it looks like nu and G are the only ones needed for the acceleration
    // calc. The others can be freed after the parent is done with them.
    //
    SpatialMat& P = updP(s);
    P = getMk(s);
    for (int i=0 ; i<(int)children.size() ; i++) {
        // this version is readable
        // P += orthoTransform( children[i]->tau * children[i]->P ,
        //                      transpose(children[i]->phiT) );
        // this version is not
        const Mat33      lt = crossMat(children[i]->getX_GB(s).T() - getX_GB(s).T());
        const SpatialMat M  = children[i]->getTau(s) * children[i]->getP(s);
        P(0,0) += M(0,0) + lt*M(1,0) - M(0,1)*lt - lt*M(1,1)*lt;
        P(0,1) += M(0,1) + lt*M(1,1);
        P(1,0) += M(1,0) - M(1,1)*lt;
        P(1,1) += M(1,1);
    }

    const Mat<dof,dof> D = getH(s) * P * ~getH(s);
    // this will throw an exception if the matrix is ill conditioned
    updDI(s) = D.invert();
    updG(s) = P * ~getH(s) * getDI(s);

    updTau(s) = 1.; // identity matrix
    updTau(s) -= getG(s) * getH(s);
    updPsiT(s) = ~getTau(s) * ~getPhi(s);
}
 
//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcZ(const SBState& s, 
                              const SpatialVec& spatialForce) const 
{
    SpatialVec& z = updZ(s);
    z = getP(s) * getCoriolisAcceleration(s) + getGyroscopicForce(s) - spatialForce;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild    = children[i]->getZ(s);
        const PhiMatrix&  phiChild  = children[i]->getPhi(s);
        const SpatialVec& GepsChild = children[i]->getGepsilon(s);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(s)  = getInternalForce(s) - getH(s)*z; // TODO: pass in hinge forces
    updNu(s)       = getDI(s) * getEpsilon(s);
    updGepsilon(s) = getG(s)  * getEpsilon(s);
}


//
// Calculate product of partial velocities J and a gradient vector on each of the
// outboard bodies. This is to be called tip to base, with temporary zTmp initialized
// to zero before the first call. Requires that Phi and H are available, so this
// should only be called in ConfiguredStage or higher. This does not change the cache at all.
// NOTE (sherm 060214): I reworked this from the original. This one no longer incorporates
// applied hinge gradients if there are any; just add those in at the end if you want them.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcInternalGradientFromSpatial
    (const SBState& s, Vector_<SpatialVec>& zTmp,
     const Vector_<SpatialVec>& X, Vector& JX) const
{
    const SpatialVec& in  = X[getNodeNum()];
    Vec<dof>&         out = Vec<dof>::updAs(&JX[getUIndex()]);
    SpatialVec&       z   = zTmp[getNodeNum()];

    z = -in;    // sherm: why the minus sign here?

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
        const PhiMatrix&  phiChild = children[i]->getPhi(s);

        z += phiChild * zChild;
    }

    out = getH(s) * z; 
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were fed to calcZ (as embodied in 'nu').
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcAccel(const SBState& s) const {
    Vec<dof>&        udot   = updUdot(s);
    const SpatialVec alphap = ~getPhi(s) * parent->getA_GB(s); // ground A_GB is 0


    udot       = getNu(s) - (~getG(s)*alphap);
    updA_GB(s) = alphap + ~getH(s)*udot + getCoriolisAcceleration(s);  

    calcJointAccel(s);   // in case joint isn't happy with just ddTheta
}

// To be called base to tip.
template<int dof> void
RigidBodyNodeSpec<dof>::calcY(const SBState& s) const {
    updY(s) = (~getH(s) * getDI(s) * getH(s)) 
                + (getPsiT(s) * parent->getY(s) * ~getPsiT(s));
}


template<int dof> void
RigidBodyNodeSpec<dof>::print(const SBState& s, int verbose) const {
    if (verbose&InternalDynamics::printNodePos) 
        cout << setprecision(8)
             << ": pos: " << getX_GB(s).T() << ' ' << '\n';
    if (verbose&InternalDynamics::printNodeTheta) 
        cout << setprecision(8)
             << ": theta: " 
             << getQ(s) << ' ' << getU(s)  << ' ' << getUdot(s)  << '\n';
}

