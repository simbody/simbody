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
void RigidBodyNode::calcJointIndependentKinematicsPos(const SBState& s) {
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 OB_OP_G = getR_GP(s) * getX_PB(s).T();

    // The Phi matrix conveniently performs child-to-parent shifting
    // on spatial quantities (forces); its transpose does parent-to-child
    // shifting for velocities.
    updPhi(s) = PhiMatrix(OB_OP_G);

    // Get spatial configuration of this body.
    updX_GB(s) = TransformMat(getR_GP(s) * getX_PB(s).R(),
                              getOP_G(s) + OB_OP_G);

    // Calculate spatial mass properties. That means we need to transform
    // the local mass moments into the Ground frame and reconstruct the
    // spatial inertia matrix Mk.

    updInertia_OB_G(s) = getInertia_OB_B(s).changeAxes(~getR_GB(s));
    updCB_G(s)         = getR_GB(s)*getCOM_B(s);

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
RigidBodyNode::calcJointIndependentKinematicsVel(const SBState& s) {
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

    /*virtual*/void calcP(const SBState&) {} 
    /*virtual*/void calcZ(const SBState&, const SpatialVec&) {} 
    /*virtual*/void calcY(const SBState&) {}
    /*virtual*/void calcInternalForce(const SBState&, const SpatialVec&) {}
    /*virtual*/void calcAccel(const SBState&) {}

    /*virtual*/void realizeModeling(const SBState&) const {}
    /*virtual*/void realizeParameters(const SBState&) const {}
    /*virtual*/void realizeConfiguration(const SBState&) {}
    /*virtual*/void realizeVelocity(const SBState&) {}
    /*virtual*/void setVelFromSVel(SBState&,const SpatialVec&) {}
    /*virtual*/void enforceQuaternionConstraints(SBState&) {}

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
                      const TransformMat&   jointFrame,
                      int&                  nextUSlot,
                      int&                  nextUSqSlot,
                      int&                  nextQSlot)
      : RigidBodyNode(mProps_B,TransformMat(),jointFrame),
        theta(0.), dTheta(0.), ddTheta(0.), forceInternal(0.)
    {
        uIndex   = nextUSlot;   nextUSlot   += getDOF();
        uSqIndex = nextUSqSlot; nextUSqSlot += getDOF()*getDOF();
        qIndex   = nextQSlot;   nextQSlot   += getMaxNQ();
    }

    virtual ~RigidBodyNodeSpec() {}

    // This is the type of the joint transition matrix H.
    typedef Mat<dof,2,Row3,1,2> HType;

    /// Calculate joint-specific kinematic quantities dependent only
    /// on positions. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for the parent,
    /// and that the position state variables q are available. The
    /// quanitites that must be computed are:
    ///   - sines & cosines of angular coordinates
    ///   - internally normalized quaternions
    ///   - X_Jb_J cross-joint transform: configuration of body's inboard
    ///              joint frame J in corresponding parent body frame Jb           
    ///   - X_PB   the configuration of the B frame in its parent's frame
    ///   - X_GB   the configuration of the B frame in ground
    ///   - H      the joint transition matrix
    virtual void calcJointKinematicsPos
        (const SBState&, Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_Jb_J, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)=0;

    /// Calculate joint-specific kinematic quantities dependent on
    /// on velocities. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for this node,
    /// that all velocity kinematics has been done for the parent, and
    /// that the velocity state variables (u) are available. The
    /// quanitites that must be computed are:
    ///   V_PB_G  relative velocity of B in P, expr. in G
    /// The code is the same for all joints, although parametrized by dof.
    void calcJointKinematicsVel(const SBState& s) {
        updV_PB_G(s) = ~getH(s) * getU(s);
    }

    void realizeModeling(const simtk::SBState&) const {
    }

    void realizeParameters(const simtk::SBState&) const {
    }

    /// Set a new configuration and calculate the consequent kinematics.
    void realizeConfiguration(const SBState& s) {
        forceInternal = 0.;  // forget these
        setJointPos(s);
        calcJointKinematicsPos(s,
            s.cache->sq,s.cache->cq,s.cache->qnorm,
            updX_JbJ(s), updX_PB(s), updX_GB(s),
            H);
        calcJointIndependentKinematicsPos(s);
    }

    /// Set new velocities for the current configuration, and calculate
    /// all the velocity-dependent terms.
    void realizeVelocity(const SBState& s) {
        // anything to initialize?
        setJointVel(s);
        calcJointKinematicsVel(s);
        calcJointIndependentKinematicsVel(s);
    }

    // These unfortunately need to be overridden for joints using quaternions.
    virtual void setJointPos(const SBState& s) {
        theta  = getQ(s);
    }
    virtual void setJointVel(const SBState& s) {
        dTheta = getU(s);
    }
    virtual void calcJointAccel(const SBState&) { }

    // We are assuming that the caller is taking care of state validity.
    virtual void getDefaultParameters(SBState& s) const {
        // TODO none yet
    }
    virtual void getDefaultConfiguration(SBState& s) const {
        updQ(s) = theta;
    }
    virtual void getDefaultVelocity(SBState& s) const {
        updU(s) = dTheta;
    }

    virtual void getAccel(SBState& s) const {
        updUdot(s) = ddTheta;
    }
    virtual void getInternalForce(SBState& s) const {
        updInternalForce(s) = forceInternal;
    }

    int          getDOF()              const { return dof; }
    virtual int  getMaxNQ()            const { return dof; } // maxNQ can be larger than dof
    virtual int  getNQ(const SBState&) const { return dof; } // DOF <= NQ <= maxNQ

    virtual void print(const SBState&, int) const;


    virtual void setVelFromSVel(SBState& s, const SpatialVec&);
    virtual void enforceQuaternionConstraints(SBState&) {}

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

    // State variables for updating; be careful. This is only appropriate for "solvers", such as
    // a method which modifies state variables to satisfy constraints.
    Vec<dof>& updQ (SBState& s) const {return toQ(s.vars->q);}
    Vec<dof>& updU (SBState& s) const {return toU(s.vars->u);} 
    Real&     upd1Q(SBState& s) const {return to1Q(s.vars->q);}
    Real&     upd1U(SBState& s) const {return to1U(s.vars->u);} 

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

    void calcP(const SBState& s);
    void calcZ(const SBState& s, const SpatialVec& spatialForce);
    void calcY(const SBState& s);
    void calcAccel(const SBState& s);
    void calcInternalGradientFromSpatial(const SBState&, Vector_<SpatialVec>& zTmp,
                                         const Vector_<SpatialVec>& X, Vector& JX);

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
                    int&                  nextUSlot,
                    int&                  nextUSqSlot,
                    int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,TransformMat(),nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos
        (const SBState& s, 
         Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_JbJ, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)
    {
        // no angles or quaternions

        const TransformMat& refX_PB = getRefX_PB(s);    // reference config of B in P
        const TransformMat& X_BJ    = getX_BJ(s);       // fixed config of J in B
        const TransformMat& X_PJb   = getX_PJb(s);      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const TransformMat& X_GP    = getX_GP(s); // parent orientation in ground

        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A Cartesian joint can't change orientation. 
        X_JbJ = TransformMat(RotationMat(), getQ(s));

        // T_JbB is B's origin measured from Jb, expr. in Jb
        const Vec3 T_JbB = X_JbJ.T() + (~X_BJ).T(); // TODO: precalculate X_JB.
        X_PB             = TransformMat(refX_PB.R(), refX_PB.T() + X_PJb.R() * T_JbB);
        X_GB             = X_GP * X_PB;

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
                  const TransformMat&   jointFrame,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,jointFrame,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos
        (const SBState& s, 
         Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_JbJ, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)
    {
        // no angles or quaternions

        const TransformMat& refX_PB = getRefX_PB(s);    // reference config of B in P
        const TransformMat& X_BJ    = getX_BJ(s);       // fixed config of J in B
        const TransformMat& X_PJb   = getX_PJb(s);      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const TransformMat& X_GP    = getX_GP(s); // parent configuration in ground

        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A sliding joint can't change orientation, and only translates along z. 
        X_JbJ = TransformMat(RotationMat(), Vec3(0.,0.,get1Q(s)));

        // T_JbB is B's origin measured from Jb, expr. in Jb
        const Vec3 T_JbB = X_JbJ.T() + (~X_BJ).T(); // TODO: precalculate X_JB.
        X_PB             = TransformMat(refX_PB.R(), refX_PB.T() + X_PJb.R() * T_JbB);
        X_GB             = X_GP * X_PB;

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
                  const TransformMat&   jointFrame,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,jointFrame,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos
        (const SBState& s, 
         Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_JbJ, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)
    { 
        const Real&            q  = get1Q(s);    // angular coordinate
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // already calculated

        Real& sq = to1Q(sine);      // nicer names for our entries
        Real& cq = to1Q(cosine);
        
        sq = sin(q); cq = cos(q);
        // no quaternions

        // This is the rotation matrix corresponding to a rotation by q about z.
        // Any Euler sequence x-y-3 with the first two rotations zero matches.
        const Mat33 a( cq , -sq , 0. ,
                       sq ,  cq , 0. ,
                       0. ,  0. , 1. );

        // We're only updating the orientation here because a torsion joint
        // can't translate.
        X_JbJ = TransformMat(RotationMat::trustMe(a), Vec3(0.));

        // Note that both orientation and position of B in P can change
        X_PB = X_PJb * X_JbJ * ~X_BJ;
        X_GB = X_GP * X_PB;
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
                  const TransformMat&   jointFrame,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,jointFrame,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos
        (const SBState& s, 
         Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_JbJ, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)
    { 
        const Vec2& q  = getQ(s); // angular coordinates
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // already calculated

        Vec2& sq = toQ(sine);      // nicer names for our entries
        Vec2& cq = toQ(cosine);
        
        sq[0]=sin(q[0]); cq[0]=cos(q[0]); sq[1]=sin(q[1]); cq[1]=cos(q[1]);
        // no quaternions

        // This is a space-fixed rotation sequence, with a rotation first
        // about the x axis by phi==q[0], then about the original y axis
        // by psi=q[1]. TODO: should be part of RotationMat class.
        const Mat33 a //Ry(psi) * Rx(phi)
            (cq[1] , sq[1]*sq[0] , sq[1]*cq[0],
               0.  ,    cq[0]    ,   -sq[0]   ,
            -sq[1] , cq[1]*sq[0] , cq[1]*cq[0]);

        // We're only updating the orientation here because a U-joint
        // can't translate.
        X_JbJ = TransformMat(RotationMat::trustMe(a), Vec3(0.));

        // Note that both orientation and position of B in P can change
        X_PB = X_PJb * X_JbJ * ~X_BJ;
        X_GB = X_GP * X_PB;
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
                           const TransformMat&   jointFrame,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,jointFrame,nextUSlot,nextUSqSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos
        (const SBState& s, 
         Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_JbJ, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)
    { 
        const Vec<5>&       q     = getQ(s);     // joint coordinates
        const TransformMat& X_BJ  = getX_BJ(s);  // fixed
        const TransformMat& X_PJb = getX_PJb(s); // fixed
        const TransformMat& X_GP  = getX_GP(s);  // already calculated

        Vec2& sq = toQ(sine).updSubVec<2>(0);    // nicer names for our entries
        Vec2& cq = toQ(cosine).updSubVec<2>(0);
        
        sq[0]=sin(q[0]); cq[0]=cos(q[0]); sq[1]=sin(q[1]); cq[1]=cos(q[1]);
        // no quaternions

        // This is a space-fixed 1-2-0 rotation sequence, with a rotation first
        // about the x axis by phi==q[0], then about the original y axis
        // by psi=q[1]. TODO: should be part of RotationMat class.
        const Mat33 a //Ry(psi) * Rx(phi)
            (cq[1] , sq[1]*sq[0] , sq[1]*cq[0],
               0.  ,    cq[0]    ,   -sq[0]   ,
            -sq[1] , cq[1]*sq[0] , cq[1]*cq[0]);

        X_JbJ = TransformMat(RotationMat::trustMe(a), q.getSubVec<3>(2));

        X_PB = X_PJb * X_JbJ * ~X_BJ; // TODO: precalculate X_JB=~X_BJ
        X_GB = X_GP * X_PB;
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
/*
/// This class contains all the odd things required by a ball joint.
/// Any RBNode joint type which contains a ball should define a member
/// of this class and delegate to it.
class ContainedBallJoint {
public:
    ContainedBallJoint() { }

    int getBallDOF()                const {return 3;}
    int getBallMaxNQ()              const {return 4;} 
    int getBallNQ(const SBState& s) const {return getUseEulerAngles(s) ? 3 : 4;} 

    void setBallPos(int qIndex, const Vector& posv, Vec3& theta) {
        Mat43   E;  // qdot = 0.5*E*u (and u=2*~E*qdot)

        if (useEuler) 
            theta = Vec3::getAs(&posv[qIndex]);
        else {
            q = Vec4::getAs(&posv[qIndex]);
            //TODO: should normalize q here?
            E = Mat43(-q[1],-q[2],-q[3],
                       q[0], q[3],-q[2],    // TODO: signs???
                      -q[3], q[0], q[1],
                       q[2],-q[1], q[0]);
        }
    } 

    // ballIndex is the offset in q at which we can find the first ball coordinate.
    void getBallDefaultConfig(SBState& s, int ballIndex) const {
        if (getUseEulerAngles(s)) Vec3::updAs(&s.vars->q[ballIndex]) = Vec3(0,0,0);
        else                      Vec4::updAs(&s.vars->q[ballIndex]) = Vec4(1,0,0,0);
    }

    void setBallVel(int uIndex, const Vector& velv, Vec3& dTheta) { 
        Mat43 dE; // dE/dt
        dTheta = Vec3::getAs(&velv[uIndex]);
        if (!useEuler) {
            dq = 0.5*(E*dTheta);
            dE = Mat43(-dq[1],-dq[2],-dq[3],
                        dq[0], dq[3],-dq[2],    // TODO: signs???
                       -dq[3], dq[0], dq[1],
                        dq[2],-dq[1], dq[0]);
        }
    }

    void getBallVel(const Vec3& dTheta, int uIndex, Vector& velv) const {
        Vec3::updAs(&velv[uIndex]) = dTheta;
        // TODO: dq??
    }

    void calcBallAccel(const Vec3& omega, const Vec3& dOmega)
    {
        // called after calcAccel
        if (useEuler) return; // nothing to do here -- ddTheta is dOmega
        ddq = 0.5*(dE*omega + E*dOmega);
    }

    void getBallAccel(const Vec3& ddTheta, int uIndex, Vector& accv) const
    {
        Vec3::updAs(&accv[uIndex]) = ddTheta;
        // TODO: ddq ??
    }

    void calcR_PB(const SBState& s, RotationMat& R_PB) {
        if (getUseEulerAngles(s)) {
            // theta = (Phi, Theta, Psi) Euler ``3-2-1'' angles
            const Vec3&  th  = getQ(s);
            Vec<3,Vec2>& scq = updSinCosQ(s);
            for (int i=0; i<3; ++i) {
                const Real scaledTh = th[i]*RigidBodyNode::DEG2RAD; // TODO: get rid of this junk!
                scq[i] = Vec2(sin(scaledTh), cos(scaledTh));
            }

            const Real sPhi   = scq[0][0], cPhi   = scq[0][1];
            const Real sTheta = scq[1][0], cTheta = scq[1][1];
            const Real sPsi   = scq[2][0], cPsi   = scq[2][1];
            
            // (sherm 050726) This matches Kane's Body-three 3-2-1 sequence on page 423
            // of Spacecraft Dynamics.
            const Mat33 R_JiJ
                ( cPhi*cTheta , -sPhi*cPsi+cPhi*sTheta*sPsi , sPhi*sPsi+cPhi*sTheta*cPsi,
                  sPhi*cTheta ,  cPhi*cPsi+sPhi*sTheta*sPsi ,-cPhi*sPsi+sPhi*sTheta*cPsi,
                 -sTheta      ,  cTheta*sPsi                , cTheta*cPsi               );
            R_PB = RotationMat::trustMe(R_JiJ); // because P=Ji and B=J for this kind of joint
        } else {
            // (sherm 060214) Added normalization of q's here. We don't want geometry-distorting
            // rotation matrices under any circumstances.
            const Vec4 q = getQuat(s)/getQuat(s).norm();
            const Real q00=q[0]*q[0], q11=q[1]*q[1], q22=q[2]*q[2], q33=q[3]*q[3];
            const Real q01=q[0]*q[1], q02=q[0]*q[2], q03=q[0]*q[3];
            const Real q12=q[1]*q[2], q13=q[1]*q[3], q23=q[2]*q[3];

            const Mat33 R_JiJ  //rotation matrix - active-sense coordinates
                (q00+q11-q22-q33,   2*(q12-q03)  ,   2*(q13+q02),
                   2*(q12+q03)  , q00-q11+q22-q33,   2*(q23-q01),
                   2*(q13-q02)  ,   2*(q23+q01)  , q00-q11-q22+q33);
            R_PB = RotationMat::trustMe(R_JiJ); // see above
        }
    }

    // Fix up the quaternions in posv. Side effect: also cleans up q and dq locally.
    void enforceBallConstraints(int qIndex, int uIndex, Vector& posv, Vector& velv) {
        if ( !useEuler ) {
            q  = Vec4::getAs(&posv[qIndex]);
            //dq = Vec4::getAs(&velv[offset]);

            q  /= q.norm();     // Normalize Euler parameters at each time step.
            //dq -= dot(q,dq)*q; // Also fix velocity: error is prop. to position component.
            const Mat43 ENorm(-q[1],-q[2],-q[3],
                               q[0], q[3],-q[2],    // TODO: signs???
                              -q[3], q[0], q[1],
                               q[2],-q[1], q[0]);
            const Vec3 w = Vec3::getAs(&velv[uIndex]);
            dq = 0.5*(ENorm*w);

            Vec4::updAs(&posv[qIndex]) =  q;
        }
    }


    void getBallInternalForce(const SBState& s, const Vec3& forceInternal, int uIndex, Vector& v) const {
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
    }

    void setBallDerivs(const Vec3& omega) {
        assert( !useEuler );
        dq = 0.5*E*omega;
    } 
};

/// Ball joint. This provides three degrees of rotational freedom, i.e.,
/// unrestricted orientation.
/// The joint frame J is aligned with the body frame B.
class RBNodeRotate3 : public RigidBodyNodeSpec<3> {
    ContainedBallJoint ball;
public:
    virtual const char* type() { return "rotate3"; }

    RBNodeRotate3(const MassProperties& mProps_B,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,TransformMat(),nextUSlot,nextUSqSlot,nextQSlot)
    {
    }
    
    int getMaxNQ()              const {return ball.getBallMaxNQ();}
    int getNQ(const SBState& s) const {return ball.getBallNQ(s);} 

    void setJointPos(const SBState& s) {
        const Vector& posv = s.vars->q;
        ball.setBallPos(qIndex, posv, theta);
    } 

    void getDefaultConfiguration(SBState& s) const {
        ball.getBallDefaultConfig(s, qIndex+0);
    }
    
    void getDefaultVelocity(SBState& s) const {
        updU(s) = 0; // no funny business here
    }

    // setPos must have been called previously
    void setJointVel(const SBState& s) {
        const Vector& velv = s.vars->u;
        ball.setBallVel(uIndex, velv, dTheta);
    }


    void calcJointAccel() {
        ball.calcBallAccel(dTheta, ddTheta);
    }

    void getAccel(SBState& s) const {
        Vector& accv = s.cache->udot;
        ball.getBallAccel(ddTheta, uIndex, accv);
    }

    void calcJointKinematicsPos
        (const SBState& s, 
         Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_JbJ, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)
    {
        updX_PB(s).updTranslation() = refOrigin_P; // ball joint can't move B origin in P
        ball.calcR_PB(s, updX_PB(s).updRotation());

        // H matrix in space-fixed (P) coords
        updH(s)[0] = SpatialRow( ~getR_GP(s)(0), Row3(0) );
        updH(s)[1] = SpatialRow( ~getR_GP(s)(1), Row3(0) );
        updH(s)[2] = SpatialRow( ~getR_GP(s)(2), Row3(0) );
   }

    void enforceQuaternionConstraints(const SBState& s) {
        Vector& posv = s.vars->q;
        Vector& velv = s.vars->u;
        ball.enforceBallConstraints(qIndex, uIndex, posv, velv);
    }

    void getInternalForce(const SBState& s) const {
        Vector& f = s.cache->netHingeForces;
        ball.getBallInternalForce(forceInternal, uIndex, f);
    }

    void setVelFromSVel(SBState& s, const SpatialVec& sVel) {
        RigidBodyNodeSpec<3>::setVelFromSVel(s, sVel);
        ball.setBallDerivs(dTheta);
    } 
};

/// Free joint. This is a six degree of freedom joint providing unrestricted 
/// translation and rotation for a free rigid body.
/// The joint frame J is aligned with the body frame B.
class RBNodeTranslateRotate3 : public RigidBodyNodeSpec<6> {
    ContainedBallJoint ball;
public:
    virtual const char* type() { return "full"; }

    RBNodeTranslateRotate3(const MassProperties& mProps_B,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<6>(mProps_B,TransformMat(),nextUSlot,nextUSqSlot,nextQSlot)
    {
    }
    
    int getMaxNQ()              const { return ball.getBallMaxNQ() + 3; }
    int getNQ(const SBState& s) const { return ball.getBallNQ(getUseEulerAngles(s)) + 3; } 

    void setJointPos(const SBState& s) {
        const Vector& posv = s.vars->q;
        Vec3 th;
        ball.setBallPos(qIndex, posv, th);
        theta.updSubVec<3>(0) = th;
        theta.updSubVec<3>(3) = Vec3::getAs(&posv[qIndex + ball.getBallNQ(getUseEulerAngles(s))]);
    } 

    void getDefaultConfiguration(SBState& s) const {
        Vector& posv = s.vars->q;
        ball.getBallPos(theta.getSubVec<3>(0), qIndex, posv);
        Vec3::updAs(&posv[qIndex+ball.getBallNQ(getUseEulerAngles(s))]) 
            = theta.getSubVec<3>(3);
    }

    // setPos must have been called previously
    void setJointVel(const SBState& s) {
        const Vector& velv = s.vars->q;
        Vec3 dTh;
        ball.setBallVel(uIndex, velv, dTh);
        dTheta.updSubVec<3>(0) = dTh;
        dTheta.updSubVec<3>(3) = Vec3::getAs(&velv[uIndex + ball.getBallDOF()]);
    }

    void getDefaultVelocity(SBState& s) const {
        Vector& velv = s.vars->u;
        ball.getBallVel(dTheta.getSubVec<3>(0), uIndex, velv);
        Vec3::updAs(&velv[uIndex+ball.getBallDOF()]) 
            = dTheta.getSubVec<3>(3);
    }

    void calcJointAccel() {
        // get angular vel/accel in the space-fixed frame
        const Vec3 omega  =  dTheta.getSubVec<3>(0);
        const Vec3 dOmega = ddTheta.getSubVec<3>(0);
        ball.calcBallAccel(omega, dOmega);
    }

    void getAccel(SBState& s) const {
        Vector& accv = s.cache->udot;
        ball.getBallAccel(ddTheta.getSubVec<3>(0), uIndex, accv);
        Vec3::updAs(&accv[uIndex+ball.getBallDOF()]) 
            = ddTheta.getSubVec<3>(3);
    }

    void calcJointKinematicsPos
        (const SBState& s, 
         Vector& sine, Vector& cosine, Vector& qnorm,
         TransformMat& X_JbJ, TransformMat& X_PB,
         TransformMat& X_GB, HType& H)
    {
        updX_PB(s).updTranslation() = refOrigin_P + theta.getSubVec<3>(3);
        ball.calcR_PB(theta.getSubVec<3>(0), updX_PB(s).updRotation());

        // H matrix in space-fixed (P) coords
        for (int i=0; i<3; ++i) {
            updH(s)[i]   = SpatialRow( ~getR_GP(s)(i),     Row3(0) );
            updH(s)[i+3] = SpatialRow(     Row3(0),     ~getR_GP(s)(i) );
        }
    }

    void enforceQuaternionConstraints(const SBState& s) {
        Vector& posv = s.vars->q;
        Vector& velv = s.vars->u;
        ball.enforceBallConstraints(qIndex, uIndex, posv, velv);
    }

    void getInternalForce(const SBState& s) const {
        Vector& f = s.cache->netHingeForces;
        const Vec3 torque = forceInternal.getSubVec<3>(0);
        ball.getBallInternalForce(torque, uIndex, f);
        Vec3::updAs(&f[uIndex + ball.getBallDOF()]) =
            forceInternal.getSubVec<3>(3);
    }

    void setVelFromSVel(SBState& s, const SpatialVec& sVel) {
        RigidBodyNodeSpec<6>::setVelFromSVel(s, sVel);
        const Vec3 omega  = dTheta.getSubVec<3>(0);
        ball.setBallDerivs(omega);
    } 
};
*/

////////////////////////////////////////////////
// RigidBodyNode factory based on joint type. //
////////////////////////////////////////////////

/*static*/ RigidBodyNode*
RigidBodyNode::create(
    const MassProperties& m,            // mass properties in body frame
    const TransformMat&   jointFrame,   // inboard joint frame J in body frame
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
        return new RBNodeTorsion(m,jointFrame,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Universal:        
        return new RBNodeRotate2(m,jointFrame,nxtUSlot,nxtUSqSlot,nxtQSlot);
    //case Joint::Orientation:
    //    return new RBNodeRotate3(m,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Cartesian:
        return new RBNodeTranslate(m,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::FreeLine:
        return new RBNodeTranslateRotate2(m,jointFrame,nxtUSlot,nxtUSqSlot,nxtQSlot);
    //case Joint::Free:
    //    return new RBNodeTranslateRotate3(m,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Joint::Sliding:
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
RigidBodyNodeSpec<dof>::setVelFromSVel(SBState& s, const SpatialVec& sVel) {
    updU(s) = getH(s) * (sVel - (~getPhi(s) * parent->getV_GB(s)));
}

//
// Calculate Pk and related quantities. The requires that the children
// of the node have already had their quantities calculated, i.e. this
// is a tip to base recursion.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcP(const SBState& s) {
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
        const Mat33      lt = crossMat(children[i]->getOB_G(s) - getOB_G(s));
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
RigidBodyNodeSpec<dof>::calcZ(const SBState& s, const SpatialVec& spatialForce) {
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
     const Vector_<SpatialVec>& X, Vector& JX)
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
RigidBodyNodeSpec<dof>::calcAccel(const SBState& s) {
    Vec<dof>&        udot   = updUdot(s);
    const SpatialVec alphap = ~getPhi(s) * parent->getA_GB(s); // ground A_GB is 0


    udot       = getNu(s) - (~getG(s)*alphap);
    updA_GB(s) = alphap + ~getH(s)*udot + getCoriolisAcceleration(s);  

    calcJointAccel(s);   // in case joint isn't happy with just ddTheta
}

// To be called base to tip.
template<int dof> void
RigidBodyNodeSpec<dof>::calcY(const SBState& s) {
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

