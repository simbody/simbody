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

void RigidBodyNode::addChild(RigidBodyNode* child, const TransformMat& referenceFrame) {
    children.push_back( child );
    child->setParent(this);
    child->refOrigin_P = referenceFrame.getTranslation(); // ignore rotation for now, it's always identity
    child->X_GB = TransformMat(X_GB.getRotation(), 
                               X_GB.getTranslation() + child->refOrigin_P);
    child->COM_G = child->X_GB.getTranslation() + child->COMstation_G;
}

//
// Calc posCM, mass, Mk
//      phi, inertia
// Should be calc'd from base to tip.
void RigidBodyNode::calcJointIndependentKinematicsPos() {
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 OB_OP_G = getR_GP() * X_PB.getTranslation();

    // The Phi matrix conveniently performs parent-to-child shifting
    // on spatial quantities.
    phi = PhiMatrix(OB_OP_G);

    // Get spatial configuration of this body.
    X_GB = TransformMat(getR_GP() * X_PB.getRotation(),
                        getOP_G() + OB_OP_G);

    // Calculate spatial mass properties. That means we need to transform
    // the local mass moments into the Ground frame and reconstruct the
    // spatial inertia matrix Mk.

    inertia_OB_G = getInertia_OB_B().changeAxes(~getR_GB());
    COMstation_G = getR_GB()*getCOM_B();

    COM_G = X_GB.getTranslation() + COMstation_G;

    // Calc Mk: the spatial inertia matrix about the body origin.
    // Note that this is symmetric; offDiag is *skew* symmetric so
    // that transpose(offDiag) = -offDiag.
    const Mat33 offDiag = getMass()*crossMat(COMstation_G);
    Mk = SpatialMat( inertia_OB_G.toMat33() ,     offDiag ,
                           -offDiag         , getMass()*Mat33(1) );
}

// Calculate velocity-related quantities: spatial velocity (sVel),
// gyroscopic force b, coriolis acceleration a. This must be
// called base to tip: depends on parent's sVel, V_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel() {
    setSpatialVel(~phi*parent->getSpatialVel() + V_PB_G);
    const Vec3& omega   = getSpatialAngVel();

    b = SpatialVec(     omega % (inertia_OB_G*omega),         // gyroscopic moment
                   getMass()*(omega % (omega%COMstation_G))); // gyroscopic force

    const Vec3& vel    = getSpatialLinVel();
    const Vec3& pOmega = parent->getSpatialAngVel();
    const Vec3& pVel   = parent->getSpatialLinVel();

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
    a =   SpatialVec(Vec3(0), pOmega % (vel-pVel)) 
        + crossMat(pOmega) * V_PB_G;
    //  + crossMat(omega)  * V_PB_G;  <-- should work too?
}

Real RigidBodyNode::calcKineticEnergy() const {
    const Real ret = dot(sVel , Mk*sVel);
    return 0.5*ret;
}


void RigidBodyNode::nodeDump(std::ostream& o) const {
    o << "NODE DUMP level=" << level << " type=" << type() << std::endl;
    nodeSpecDump(o);
    o << "END OF NODE type=" << type() << std::endl;
}

std::ostream& operator<<(std::ostream& s, const RigidBodyNode& node) {
    node.nodeDump(s);
    return s;
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
      : RigidBodyNode(MassProperties(),Vec3(0.),TransformMat()) {}
    ~RBGroundBody() {}

    /*virtual*/const char* type() const { return "ground"; }

    /*virtual*/void calcP() {} 
    /*virtual*/void calcZ(const SpatialVec&) {} 
    /*virtual*/void calcY() {}
    /*virtual*/void calcInternalForce(const SpatialVec&) {}
    /*virtual*/void calcAccel() {}

    /*virtual*/void realizeModeling(const SBState&) const {}
    /*virtual*/void realizeParameters(const SBState&) const {}
    /*virtual*/void realizeConfiguration(const Vector&) {}
    /*virtual*/void realizeVelocity(const Vector&) {}
    /*virtual*/void setVelFromSVel(const SpatialVec&) {}
    /*virtual*/void enforceConstraints(Vector& pos, Vector& vel) {}

    /*virtual*/void getDefaultParameters(SBState&)    const {}
    /*virtual*/void getDefaultConfiguration(SBState&) const {}
    /*virtual*/void getDefaultVelocity(SBState&)      const {}

    /*virtual*/void getAccel(Vector&) const {}

    /*virtual*/void getInternalForce(Vector&) const {}
    // /*virtual*/ const SpatialRow& getHRow(int i) const;

    void print(int) {}
};

template<int dof>
class RigidBodyNodeSpec : public RigidBodyNode {
public:
    RigidBodyNodeSpec(const MassProperties& mProps_B,
                      const TransformMat&   jointFrame,
                      int&                  nextUSlot,
                      int&                  nextQSlot)
      : RigidBodyNode(mProps_B,Vec3(0.),jointFrame),
        theta(0.), dTheta(0.), ddTheta(0.), forceInternal(0.)
    {
        uIndex = nextUSlot; nextUSlot += getDOF();
        qIndex = nextQSlot; nextQSlot += getMaxNQ();
    }

    virtual ~RigidBodyNodeSpec() {}

    /// Calculate joint-specific kinematic quantities dependent only
    /// on positions. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for the parent,
    /// and that the position state variables (theta) are available. The
    /// quanitites that must be computed are:
    ///   R_PB  the orientation of the B frame in its parent's frame
    ///   OB_P  the location of B's origin in its parent (meas. from OP)
    ///   H     the joint transition matrix
    virtual void calcJointKinematicsPos()=0;

    /// Calculate joint-specific kinematic quantities dependent on
    /// on velocities. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for this node,
    /// that all velocity kinematics has been done for the parent, and
    /// that the velocity state variables (dTheta) are available. The
    /// quanitites that must be computed are:
    ///   V_PB_G  relative velocity of B in P, expr. in G
    virtual void calcJointKinematicsVel()=0;

    void realizeModeling(const simtk::SBState&) const {
    }

    void realizeParameters(const simtk::SBState&) const {
    }

    /// Set a new configuration and calculate the consequent kinematics.
    void realizeConfiguration(const Vector& posv) {
        forceInternal = 0.;  // forget these
        setJointPos(posv);
        calcJointKinematicsPos();
        calcJointIndependentKinematicsPos();
    }

    /// Set new velocities for the current configuration, and calculate
    /// all the velocity-dependent terms.
    void realizeVelocity(const Vector& velv) {
        // anything to initialize?
        setJointVel(velv);
        calcJointKinematicsVel();
        calcJointIndependentKinematicsVel();
    }

    // These unfortunately need to be overridden for joints using quaternions.
    virtual void setJointPos(const Vector& posv) {
        theta  = Vec<dof>::getAs(&posv[stateOffset]);
    }
    virtual void setJointVel(const Vector& velv) {
        dTheta = Vec<dof>::getAs(&velv[stateOffset]);
    }
    virtual void calcJointAccel() { }

    // We are assuming that the caller is taking care of state validity.
    virtual void getDefaultParameters(SBState& s) const {
        // TODO none yet
    }
    virtual void getDefaultConfiguration(SBState& s) const {
        Vec<dof>::updAs(&s.vars->q[qIndex]) = theta;
    }
    virtual void getDefaultVelocity(SBState& s) const {
        Vec<dof>::updAs(&s.vars->u[uIndex]) = dTheta;
    }

    virtual void getAccel(Vector& a) const {
        Vec<dof>::updAs(&a[stateOffset]) = ddTheta;
    }
    virtual void getInternalForce(Vector& t) const {
        Vec<dof>::updAs(&t[stateOffset]) = forceInternal;
    }

    int          getDOF()              const { return dof; }
    virtual int  getMaxNQ()            const { return dof; } // maxNQ can be larger than dof
    virtual int  getNQ(bool useAngles) const { return dof; } // DOF <= NQ <= maxNQ

    virtual void   print(int) const;


    virtual void setVelFromSVel(const SpatialVec&);
    virtual void enforceConstraints(Vector& pos, Vector& vel) {}

    const SpatialRow& getHRow(int i) const {
        return H[i];
    }

    void calcP();
    void calcZ(const SpatialVec& spatialForce);
    void calcY();
    void calcAccel();
    void calcInternalForce(const SpatialVec& spatialForce);

    void nodeSpecDump(std::ostream& o) const {
        o << "stateOffset=" << stateOffset << " mass=" << getMass() 
            << " COM_G=" << getCOM_G() << std::endl;
        o << "inertia_OB_G=" << getInertia_OB_G() << std::endl;
        o << "H=" << H << std::endl;
        o << "SVel=" << sVel << std::endl;
        o << "a=" << a << std::endl;
        o << "b=" << b << std::endl;
        o << "Th  =" << theta << std::endl;
        o << "dTh =" << dTheta << std::endl;
        o << "ddTh=" << ddTheta << std::endl;
        o << "SAcc=" << sAcc << std::endl;
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

private:
    void calcD_G(const SpatialMat& P);
};

/*static*/const double RigidBodyNode::DEG2RAD = std::acos(-1.) / 180.; // i.e., pi/180
//const double RigidBodyNode::DEG2RAD = 1.0;  //always use radians


//////////////////////////////////////////
// Derived classes for each joint type. //
//////////////////////////////////////////


/**
 * Translate (Cartesian) joint. This provides three degrees of translational freedom
 * which is suitable (e.g.) for connecting a free atom to ground.
 * The joint frame J is aligned with the body frame B.
 */
class RBNodeTranslate : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "translate"; }

    RBNodeTranslate(const MassProperties& mProps_B,
                    int&                  nextUSlot,
                    int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,TransformMat(),nextUSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos() {
        // Cartesian joint can't change orientation
        X_PB = TransformMat(RotationMat(), refOrigin_P + theta);

        // Note that this is spatial (and R_GP=R_GB for this joint)
        H[0] = SpatialRow( Row3(0), ~getR_GP()(0) );
        H[1] = SpatialRow( Row3(0), ~getR_GP()(1) );
        H[2] = SpatialRow( Row3(0), ~getR_GP()(2) );
    }

    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }
};

/**
 * This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
 * about a particular axis.
 */
class RBNodeTorsion : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "torsion"; }

    RBNodeTorsion(const MassProperties& mProps_B,
                  const TransformMat&   jointFrame,
                  int&                  nextUSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,jointFrame,nextUSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos() { 
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        const double scale=1.0;
        const double sinTau = sin( scale * theta[0] );
        const double cosTau = cos( scale * theta[0] );
        const Mat33 a( cosTau , -sinTau , 0.0 ,
                       sinTau ,  cosTau , 0.0 ,
                       0.0    ,  0.0    , 1.0 );
        const RotationMat R_JiJ = RotationMat::trustMe(a); //rotation about z-axis

        // We need R_PB=R_PJi*R_JiJ*R_JB. But R_PJi==R_BJ, so this works:
        const RotationMat& R_BJ = X_BJ.getRotation();
        X_PB = TransformMat(R_BJ * R_JiJ * ~R_BJ, refOrigin_P); // torsion can't move B origin in P
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

private:
    // Calc H matrix in space-fixed coords.
    void calcH() {
        // This only works because the joint z axis is the same in B & P
        // because that's what we rotate around.
        const RotationMat& R_BJ = X_BJ.getRotation();
        const Vec3 z = getR_GP() * R_BJ(2); // R_BJ=R_PJi
        H[0] = SpatialRow( ~z, Row3(0) );
    }  
};

/**
 * This class contains all the odd things required by a ball joint.
 * Any RBNode joint type which contains a ball should define a member
 * of this class and delegate to it.
 */
class ContainedBallJoint {
    Mat43   E;  // qdot = 0.5*E*u (and u=2*~E*qdot)
    Mat43   dE; // dE/dt
    Vec4    q; //Euler parameters for rotation relative to parent
    Vec4    dq;
    Vec4    ddq;
    double cPhi, sPhi;        //trig functions of Euler angles
    double cPsi, sPsi;        //used for minimizations
    double cTheta, sTheta;
    bool   useEuler;          //if False, use Quaternion rep.

public:
    virtual const char* type() { return "rotate3"; }

    ContainedBallJoint(bool shouldUseEuler)
      : q(1,0,0,0), dq(0), ddq(0), 
        useEuler(shouldUseEuler)
    { 
    }

    int getBallDOF() const {return 3;}
    int getBallMaxNQ() const {return 4;} 
    int getBallNQ(bool useAngles) const { 
        if (useEuler) return 3; else return 4; 
    }

    void setBallPos(int qIndex, const Vector& posv, Vec3& theta) {
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

    void getBallPos(const Vec3& theta, int qIndex, Vector& posv) const {
        if (useEuler) Vec3::updAs(&posv[qIndex]) = theta;
        else          Vec4::updAs(&posv[qIndex]) = q;
    }

    void setBallVel(int uIndex, const Vector& velv, Vec3& dTheta) { 
        dTheta = Vec3::getAs(&velv[uIndex]);
        if (useEuler) {
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

    void calcR_PB(const Vec3& theta, RotationMat& R_PB) {
        if (useEuler) {
            // theta = (Phi, Theta, Psi) Euler ``3-2-1'' angles 
            cPhi   = cos( theta[0] *RigidBodyNode::DEG2RAD );
            sPhi   = sin( theta[0] *RigidBodyNode::DEG2RAD );
            cTheta = cos( theta[1] *RigidBodyNode::DEG2RAD );
            sTheta = sin( theta[1] *RigidBodyNode::DEG2RAD );
            cPsi   = cos( theta[2] *RigidBodyNode::DEG2RAD );
            sPsi   = sin( theta[2] *RigidBodyNode::DEG2RAD );
            
            // (sherm 050726) This matches Kane's Body-three 3-2-1 sequence on page 423
            // of Spacecraft Dynamics.
            const Mat33 R_JiJ
                ( cPhi*cTheta , -sPhi*cPsi+cPhi*sTheta*sPsi , sPhi*sPsi+cPhi*sTheta*cPsi,
                  sPhi*cTheta ,  cPhi*cPsi+sPhi*sTheta*sPsi ,-cPhi*sPsi+sPhi*sTheta*cPsi,
                 -sTheta      ,  cTheta*sPsi                , cTheta*cPsi               );
            R_PB = RotationMat::trustMe(R_JiJ); // because P=Ji and B=J for this kind of joint
        } else {
            // TODO: should normalize q here??
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


    void getBallInternalForce(const Vec3& forceInternal, int uIndex, Vector& v) const {
        //dependency: calcR_PB must be called first
        assert( useEuler );

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

/**
 * Ball joint. This provides three degrees of rotational freedom, i.e.,
 * unrestricted orientation.
 * The joint frame J is aligned with the body frame B.
 */
class RBNodeRotate3 : public RigidBodyNodeSpec<3> {
    ContainedBallJoint ball;
public:
    virtual const char* type() { return "rotate3"; }

    RBNodeRotate3(const MassProperties& mProps_B,
                  int&                  nextUSlot,
                  int&                  nextQSlot,
                  bool                  useEuler)
      : RigidBodyNodeSpec<3>(mProps_B,TransformMat(),nextUSlot,nextQSlot),
        ball(useEuler)
    {
    }
    
    int getMaxNQ()           const  { return ball.getBallMaxNQ(); }
    int getNQ(bool useAngles) const { return ball.getBallNQ(useAngles); } 

    void setJointPos(const Vector& posv) {
        ball.setBallPos(qIndex, posv, theta);
    } 

    void getDefaultConfiguration(SBState& s) const {
        Vector& posv = s.vars->q;
        ball.getBallPos(theta, qIndex, posv);
    }

    // setPos must have been called previously
    void setJointVel(const Vector& velv) {
        ball.setBallVel(uIndex, velv, dTheta);
    }

    void getDefaultVelocity(SBState& s) const {
        Vector& velv = s.vars->u;
        ball.getBallVel(dTheta, uIndex, velv);
    }

    void calcJointAccel() {
        ball.calcBallAccel(dTheta, ddTheta);
    }

    void getAccel(Vector& accv) const {
        ball.getBallAccel(ddTheta, uIndex, accv);
    }

    void calcJointKinematicsPos() { 
        X_PB.updTranslation() = refOrigin_P; // ball joint can't move B origin in P
        ball.calcR_PB(theta, X_PB.updRotation());

        // H matrix in space-fixed (P) coords
        H[0] = SpatialRow( ~getR_GP()(0), Row3(0) );
        H[1] = SpatialRow( ~getR_GP()(1), Row3(0) );
        H[2] = SpatialRow( ~getR_GP()(2), Row3(0) );
   }

    // Note that dTheta = w_PB_P = ang vel of B in P, expr in P
    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

    void enforceConstraints(Vector& posv, Vector& velv) {
        ball.enforceBallConstraints(qIndex, uIndex, posv, velv);
    }

    void getInternalForce(Vector& v) const {
        ball.getBallInternalForce(forceInternal, uIndex, v);
    }

    void setVelFromSVel(const SpatialVec& sVel) {
        RigidBodyNodeSpec<3>::setVelFromSVel(sVel);
        ball.setBallDerivs(dTheta);
    } 
};

/**
 * Free joint. This is a six degree of freedom joint providing unrestricted 
 * translation and rotation for a free rigid body.
 * The joint frame J is aligned with the body frame B.
 */
class RBNodeTranslateRotate3 : public RigidBodyNodeSpec<6> {
    ContainedBallJoint ball;
public:
    virtual const char* type() { return "full"; }

    RBNodeTranslateRotate3(const MassProperties& mProps_B,
                           int&                  nextUSlot,
                           int&                  nextQSlot,
                           bool                  useEuler)
      : RigidBodyNodeSpec<6>(mProps_B,TransformMat(),nextUSlot,nextQSlot),
        ball(useEuler)
    {
    }
    
    int getMaxNQ()            const { return 3 + ball.getBallMaxNQ(); }
    int getNQ(bool useAngles) const { return 3 + ball.getBallNQ(useAngles); } 

    void setJointPos(const Vector& posv) {
        Vec3 th;
        ball.setBallPos(qIndex, posv, th);
        theta.updSubVec<3>(0) = th;
        theta.updSubVec<3>(3) = Vec3::getAs(&posv[qIndex+ball.getBallNQ()]);
    } 

    void getDefaultConfiguration(SBState& s) const {
        Vector& posv = s.vars->q;
        ball.getBallPos(theta.getSubVec<3>(0), qIndex, posv);
        Vec3::updAs(&posv[qIndex+ball.getBallNQ()]) 
            = theta.getSubVec<3>(3);
    }

    // setPos must have been called previously
    void setJointVel(const Vector& velv) {
        Vec3 dTh;
        ball.setBallVel(uIndex, velv, dTh);
        dTheta.updSubVec<3>(0) = dTh;
        dTheta.updSubVec<3>(3) = Vec3::getAs(&velv[uIndex+ball.getBallDOF()]);
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

    void getAccel(Vector& accv) const {
        ball.getBallAccel(ddTheta.getSubVec<3>(0), uIndex, accv);
        Vec3::updAs(&accv[uIndex+ball.getBallDOF()]) 
            = ddTheta.getSubVec<3>(3);
    }

    void calcJointKinematicsPos() {
        X_PB.updTranslation() = refOrigin_P + theta.getSubVec<3>(3);
        ball.calcR_PB(theta.getSubVec<3>(0), X_PB.updRotation());

        // H matrix in space-fixed (P) coords
        for (int i=0; i<3; ++i) {
            H[i]   = SpatialRow( ~getR_GP()(i),     Row3(0) );
            H[i+3] = SpatialRow(     Row3(0),     ~getR_GP()(i) );
        }
    }

    // Note that dTheta[0..2] = w_PB_P = ang vel of B in P, expr in P
    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

    void enforceConstraints(Vector& posv, Vector& velv) {
        ball.enforceBallConstraints(qIndex, uIndex, posv, velv);
    }


    void getInternalForce(Vector& v) const {
        const Vec3 torque = forceInternal.getSubVec<3>(0);
        ball.getBallInternalForce(torque, uIndex, v);
        Vec3::updAs(&v[uIndex+ball.getBallDOF()]) =
            forceInternal.getSubVec<3>(3);
    }

    void setVelFromSVel(const SpatialVec& sVel) {
        RigidBodyNodeSpec<6>::setVelFromSVel(sVel);
        const Vec3 omega  = dTheta.getSubVec<3>(0);
        ball.setBallDerivs(omega);
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
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,jointFrame,nextUSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos() { 
        X_PB.updTranslation() = refOrigin_P; // no translation with this joint
        calcR_PB(X_PB.updRotation());
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

private:
    void calcR_PB(RotationMat& R_PB) { 
        //double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0; 
        double sinPhi = sin( scale * theta[0] );
        double cosPhi = cos( scale * theta[0] );
        double sinPsi = sin( scale * theta[1] );
        double cosPsi = cos( scale * theta[1] );

        const Mat33 a //Ry(psi) * Rx(phi)
            (cosPsi , sinPsi*sinPhi , sinPsi*cosPhi,
             0      , cosPhi        , -sinPhi      ,
            -sinPsi , cosPsi*sinPhi , cosPsi*cosPhi);

        const RotationMat& R_BJ = X_BJ.getRotation();
        R_PB = R_BJ * RotationMat::trustMe(a) * ~R_BJ;
    }

    void calcH() {
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;
        const RotationMat tmpR_GB = getR_GP() * X_PB.getRotation();

        const RotationMat& R_BJ = X_BJ.getRotation();
        const Vec3 x = scale * (tmpR_GB * R_BJ(0));
        const Vec3 y = scale * (tmpR_GB * R_BJ(1));

        H[0] = SpatialRow(~x, Row3(0));
        H[1] = SpatialRow(~y, Row3(0));
    }  
};

/**
 * The "diatom" joint is the equivalent of a free joint for a body with no inertia in
 * one direction, such as one composed of just two atoms. It allows unrestricted
 * translation but rotation only about directions perpendicular to the body's
 * inertialess axis.
 */
class RBNodeTranslateRotate2 : public RigidBodyNodeSpec<5> {
public:
    virtual const char* type() { return "diatom"; }

    RBNodeTranslateRotate2(const MassProperties& mProps_B,
                           const TransformMat&   jointFrame,
                           int&                  nextUSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,jointFrame,nextUSlot,nextQSlot)
    {
    }

    void calcJointKinematicsPos() { 
        X_PB.updTranslation() = refOrigin_P + theta.getSubVec<3>(2);
        calcR_PB(X_PB.updRotation());
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

private:
    void calcR_PB(RotationMat& R_PB) { 
        // double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0;
        double sinPhi = sin( scale * theta[0] );
        double cosPhi = cos( scale * theta[0] );
        double sinPsi = sin( scale * theta[1] );
        double cosPsi = cos( scale * theta[1] );

        // space (parent)-fixed 1-2-3 sequence (rotation 3=0)
        const Mat33 R_JiJ  //Ry(psi) * Rx(phi)
            (cosPsi , sinPsi*sinPhi , sinPsi*cosPhi,
             0      , cosPhi        , -sinPhi      ,
            -sinPsi , cosPsi*sinPhi , cosPsi*cosPhi);

        // calculates R0*a*R0'  (R0=R_BJ(==R_PJi), a=R_JiJ)
        const RotationMat& R_BJ = X_BJ.getRotation();
        R_PB = R_BJ * RotationMat::trustMe(R_JiJ) * ~R_BJ; // orientation of B in parent P
    }

    void calcH() {
        //double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;
        const RotationMat& R_GP = getR_GP();
        const RotationMat tmpR_GB = R_GP * X_PB.getRotation();

        const RotationMat& R_BJ = X_BJ.getRotation();
        const Vec3 x = scale * (tmpR_GB * R_BJ(0));
        const Vec3 y = scale * (tmpR_GB * R_BJ(1));

        H[0] = SpatialRow(  ~x   ,  Row3(0));
        H[1] = SpatialRow(  ~y   ,  Row3(0));
        H[2] = SpatialRow(Row3(0), ~R_GP(0));
        H[3] = SpatialRow(Row3(0), ~R_GP(1));
        H[4] = SpatialRow(Row3(0), ~R_GP(2));
    }  
};

////////////////////////////////////////////////
// RigidBodyNode factory based on joint type. //
////////////////////////////////////////////////

/*static*/ RigidBodyNode*
RigidBodyNode::create(
    const MassProperties& m,            // mass properties in body frame
    const TransformMat&   jointFrame,   // inboard joint frame J in body frame
    Joint::JointType      type,
    bool                  isReversed,   // child-to-parent orientation?
    bool                  useEuler,
    int&                  nxtUSlot,
    int&                  nxtQSlot)  
{
    assert(!isReversed);

    switch(type) {
    case Joint::ThisIsGround:
        return new RBGroundBody();
    case Joint::Torsion:
        return new RBNodeTorsion(m,jointFrame,nxtUSlot,nxtQSlot);
    case Joint::Universal:        
        return new RBNodeRotate2(m,jointFrame,nxtUSlot,nxtQSlot);
    case Joint::Orientation:
        return new RBNodeRotate3(m,nxtUSlot,nxtQSlot,useEuler);
    case Joint::Cartesian:
        return new RBNodeTranslate(m,nxtUSlot,nxtQSlot);
    case Joint::FreeLine:
        return new RBNodeTranslateRotate2(m,jointFrame,nxtUSlot,nxtQSlot);
    case Joint::Free:
        return new RBNodeTranslateRotate3(m,nxtUSlot,nxtQSlot,useEuler);
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
RigidBodyNodeSpec<dof>::setVelFromSVel(const SpatialVec& sVel) {
    dTheta = H * (sVel - ~phi * parent->sVel);
}

template<int dof> void
RigidBodyNodeSpec<dof>::calcD_G(const SpatialMat& P) {
    const Mat<dof,dof> D = H * P * ~H;
    // this will throw an exception if the matrix is ill conditioned
    DI = D.invert();
    G = P * ~H * DI;
}

//
// Calculate Pk and related quantities. The requires that the children
// of the node have already had their quantities calculated, i.e. this
// is a tip to base recursion.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcP() {
    //
    //how much do we need to keep around?
    // it looks like nu and G are the only ones needed for the acceleration
    // calc. The others can be freed after the parent is done with them.
    //
    P = Mk;
    for (int i=0 ; i<(int)children.size() ; i++) {
        // this version is readable
        // P += orthoTransform( children[i]->tau * children[i]->P ,
        //                      transpose(children[i]->phiT) );
        // this version is not
        const Mat33      lt = crossMat(children[i]->getOB_G() - getOB_G());
        const SpatialMat M  = children[i]->tau * children[i]->P;
        P(0,0) += M(0,0) + lt*M(1,0) - M(0,1)*lt - lt*M(1,1)*lt;
        P(0,1) += M(0,1) + lt*M(1,1);
        P(1,0) += M(1,0) - M(1,1)*lt;
        P(1,1) += M(1,1);
    }

    calcD_G(P);
    tau = 1.; // identity matrix
    tau -= G * H;
    psiT = ~tau * ~phi;
}
 
//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcZ(const SpatialVec& spatialForce) {
    z = P * a + b - spatialForce;

    for (int i=0 ; i<(int)children.size() ; i++) 
        z += children[i]->phi
             * (children[i]->z + children[i]->Gepsilon);

    epsilon  = forceInternal - H*z;
    nu       = DI * epsilon;
    Gepsilon = G * epsilon;
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were fed to calcZ (as embodied in 'nu').
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcAccel() {
    // const Node* pNode = parentHinge.remNode;
    //make sure that this is phi is correct - FIX ME!
    // base alpha = 0!!!!!!!!!!!!!
    const SpatialVec alphap = ~phi * parent->sAcc;
    ddTheta = nu - ~G * alphap;
    sAcc    = alphap + ~H * ddTheta + a;  

    calcJointAccel();   // in case joint isn't happy with just ddTheta
}

// To be called base to tip.
template<int dof> void
RigidBodyNodeSpec<dof>::calcY() {
    Y = (~H * DI * H) + (psiT * parent->Y * ~psiT);
}

//
// Calculate sum of internal force and effective forces due to Cartesian
// forces.
// To be called from tip to base.
// Should be called only once after calcProps.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcInternalForce(const SpatialVec& spatialForce) {
    z = -spatialForce;

    for (int i=0 ; i<(int)children.size() ; i++) 
        z += children[i]->phi * children[i]->z;

    forceInternal += H * z; 
}

template<int dof> void
RigidBodyNodeSpec<dof>::print(int verbose) const {
    if (verbose&InternalDynamics::printNodePos) 
        cout << setprecision(8)
             << ": pos: " << X_GB.getTranslation() << ' ' << '\n';
    if (verbose&InternalDynamics::printNodeTheta) 
        cout << setprecision(8)
             << ": theta: " 
             << theta << ' ' << dTheta  << ' ' << ddTheta  << '\n';
}

