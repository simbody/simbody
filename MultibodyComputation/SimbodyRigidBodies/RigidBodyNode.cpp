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

void RigidBodyNode::addChild(RigidBodyNode* child, const Frame& referenceFrame) {
    children.push_back( child );
    child->setParent(this);
    child->refOrigin_P = referenceFrame.getOrigin(); // ignore frame for now, it's always identity
    child->R_GB = R_GB;
    child->OB_G = OB_G + child->refOrigin_P;
    child->COM_G = child->OB_G + child->COMstation_G;
}

//
// Calc posCM, mass, Mk
//      phi, inertia
// Should be calc'd from base to tip.
void RigidBodyNode::calcJointIndependentKinematicsPos() {
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 OB_OP_G = getR_GP()*OB_P;

    // The Phi matrix conveniently performs parent-to-child shifting
    // on spatial quantities.
    phi = PhiMatrix(OB_OP_G);

    // Get spatial configuration of this body.
    R_GB = getR_GP() * R_PB;
    OB_G = getOP_G() + OB_OP_G;

    // Calculate spatial mass properties. That means we need to transform
    // the local mass moments into the Ground frame and reconstruct the
    // spatial inertia matrix Mk.

    inertia_OB_G = getInertia_OB_B().changeAxes(~getR_GB());
    COMstation_G = getR_GB()*getCOM_B();

    COM_G = OB_G + COMstation_G;

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
    setSpatialVel(~phi * parent->getSpatialVel()
                  + V_PB_G);
    const Vec3& omega   = getSpatialAngVel();
    const Vec3  gMoment = omega % (inertia_OB_G * omega);
    const Vec3  gForce  = getMass() * (omega % (omega % COMstation_G));
    b = SpatialVec(gMoment, 
                   gForce);

    const Vec3& vel    = getSpatialLinVel();
    const Vec3& pOmega = parent->getSpatialAngVel();
    const Vec3& pVel   = parent->getSpatialLinVel();

    // calc a: coriolis acceleration
    a = SpatialMat(crossMat(pOmega),    Mat33(0.),
                      Mat33(0.)    , crossMat(pOmega)) 
        * V_PB_G;
    a += SpatialVec(     Vec3(0.), 
                    pOmega % (vel-pVel));
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
      : RigidBodyNode(MassProperties(),Vec3(0.),MatRotation(),Vec3(0.)) {}
    ~RBGroundBody() {}

    /*virtual*/const char* type() const { return "ground"; }

    /*virtual*/void calcP() {} 
    /*virtual*/void calcZ(const SpatialVec&) {} 
    /*virtual*/void calcY() {}
    /*virtual*/void calcInternalForce(const SpatialVec&) {}
    /*virtual*/void calcAccel() {}

    /*virtual*/void setPos(const Vector&) {}
    /*virtual*/void setVel(const Vector&) {}
    /*virtual*/void setVelFromSVel(const SpatialVec&) {}
    /*virtual*/void enforceConstraints(Vector& pos, Vector& vel) {}

    /*virtual*/void getPos(Vector&)   const {}
    /*virtual*/void getVel(Vector&)   const {}
    /*virtual*/void getAccel(Vector&) const {}

    /*virtual*/void getInternalForce(Vector&) const {}
    // /*virtual*/ const SpatialRow& getHRow(int i) const;

    void print(int) {}
};

template<int dof>
class RigidBodyNodeSpec : public RigidBodyNode {
public:
    RigidBodyNodeSpec(const MassProperties& mProps_B,
                      const Frame& jointFrame,
                      int& cnt)
      : RigidBodyNode(mProps_B,Vec3(0.),jointFrame.getAxes(),jointFrame.getOrigin()),
        theta(0.), dTheta(0.), ddTheta(0.), forceInternal(0.)
    {
        stateOffset = cnt;
        cnt += dof;
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

    /// Set a new configuration and calculate the consequent kinematics.
    void setPos(const Vector& posv) {
        forceInternal = 0.;  // forget these
        setJointPos(posv);
        calcJointKinematicsPos();
        calcJointIndependentKinematicsPos();
    }

    /// Set new velocities for the current configuration, and calculate
    /// all the velocity-dependent terms.
    void setVel(const Vector& velv) {
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

    virtual void getPos(Vector& p) const {
        Vec<dof>::updAs(&p[stateOffset]) = theta;
    }
    virtual void getVel(Vector& v) const {
        Vec<dof>::updAs(&v[stateOffset]) = dTheta;
    }   
    virtual void getAccel(Vector& a) const {
        Vec<dof>::updAs(&a[stateOffset]) = ddTheta;
    }
    virtual void getInternalForce(Vector& t) const {
        Vec<dof>::updAs(&t[stateOffset]) = forceInternal;
    }

    int          getDOF() const { return dof; }
    virtual int  getDim() const { return dof; } // dim can be larger than dof

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
    Vec<dof>            theta;   // internal coordinates
    Mat<dof,2,Row3,1,2> H;       // joint transition matrix (spatial); note row-order packing
                                 //   so that transpose will be a nice column-packed matrix
    Mat<dof,dof>        DI;
    Mat<2,dof,Vec3>     G;       // structure is like ~H; that is, default column-packed matrix

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
                    int&                  nextStateOffset)
      : RigidBodyNodeSpec<3>(mProps_B,Frame(),nextStateOffset)
    {
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P + theta;
        R_PB = MatRotation(); // Cartesian joint can't change orientation

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
                  const Frame&          jointFrame,
                  int&                  nextStateOffset)
      : RigidBodyNodeSpec<1>(mProps_B,jointFrame,nextStateOffset)
    {
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P; // torsion joint can't move B origin in P
        calcR_PB();
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

private:
    void calcR_PB() {
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0;
        double sinTau = sin( scale * theta(0) );
        double cosTau = cos( scale * theta(0) );
        const Mat33 a( cosTau , -sinTau , 0.0 ,
                       sinTau ,  cosTau , 0.0 ,
                       0.0    ,  0.0    , 1.0 );
        const MatRotation R_JiJ = MatRotation::trustMe(a); //rotation about z-axis

        // We need R_PB=R_PJi*R_JiJ*R_JB. But R_PJi==R_BJ, so this works:
        R_PB = R_BJ * R_JiJ * ~R_BJ;
    };

    // Calc H matrix in space-fixed coords.
    void calcH() {
        // This only works because the joint z axis is the same in B & P
        // because that's what we rotate around.
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
    Vec4    q; //Euler parameters for rotation relative to parent
    Vec4    dq;
    Vec4    ddq;
    double cPhi, sPhi;        //trig functions of Euler angles
    double cPsi, sPsi;        //used for minimizations
    double cTheta, sTheta;
    bool   useEuler;          //if False, use Quaternion rep.

public:
    virtual const char* type() { return "rotate3"; }

    ContainedBallJoint(int& cnt, bool shouldUseEuler)
      : q(1.0,0.0,0.0,0.0), dq(0.), ddq(0.), 
        useEuler(shouldUseEuler)
    { 
        if ( !useEuler )
            cnt++;
    }
    
    int  getBallDim() const { 
        if (useEuler) return 3; else return 4; 
    }

    void setBallPos(int stateOffset, const Vector& posv, Vec3& theta) {
        if (useEuler) theta = Vec3::getAs(&posv[stateOffset]);
        else          q     = Vec4::getAs(&posv[stateOffset]);
    } 

    void getBallPos(const Vec3& theta, int stateOffset, Vector& posv) const {
        if (useEuler) Vec3::updAs(&posv[stateOffset]) = theta;
        else          Vec4::updAs(&posv[stateOffset]) = q;
    }

    void setBallVel(int stateOffset, const Vector& velv, Vec3& dTheta) { 
        if (useEuler)
            dTheta = Vec3::getAs(&velv[stateOffset]);
        else {
            dq = Vec4::getAs(&velv[stateOffset]);
            const Mat34 M(-q(1), q(0),-q(3), q(2),
                          -q(2), q(3), q(0),-q(1),
                          -q(3),-q(2), q(1), q(0));
            dTheta = 2.0*( M * dq );
        }
    }

    void getBallVel(const Vec3& dTheta, int stateOffset, Vector& velv) const {
        if (useEuler) Vec3::updAs(&velv[stateOffset]) = dTheta;
        else          Vec4::updAs(&velv[stateOffset]) = dq;
    }

    void calcBallAccel(const Vec3& omega, const Vec3& dOmega)
    {
        // called after calcAccel
        if (useEuler) return; // nothing to do here -- ddTheta is dOmega

        const Mat43 M(-q(1),-q(2),-q(3),
                       q(0), q(3),-q(2),
                      -q(3), q(0), q(1),
                       q(2),-q(1), q(0));

        const Mat43 dM(-dq(1),-dq(2),-dq(3),
                        dq(0), dq(3),-dq(2),
                       -dq(3), dq(0), dq(1),
                        dq(2),-dq(1), dq(0));

        ddq = 0.5*(dM*omega + M*dOmega);
    }

    void getBallAccel(const Vec3& ddTheta, int stateOffset, Vector& accv) const
    {
        if (useEuler) Vec3::updAs(&accv[stateOffset]) = ddTheta;
        else          Vec4::updAs(&accv[stateOffset]) = ddq;
    }

    void calcR_PB(const Vec3& theta, MatRotation& R_PB) {
        if (useEuler) {
            // theta = (Phi, Theta, Psi) Euler ``3-2-1'' angles 
            cPhi   = cos( theta(0) *RigidBodyNode::DEG2RAD );
            sPhi   = sin( theta(0) *RigidBodyNode::DEG2RAD );
            cTheta = cos( theta(1) *RigidBodyNode::DEG2RAD );
            sTheta = sin( theta(1) *RigidBodyNode::DEG2RAD );
            cPsi   = cos( theta(2) *RigidBodyNode::DEG2RAD );
            sPsi   = sin( theta(2) *RigidBodyNode::DEG2RAD );
            
            // (sherm 050726) This matches Kane's Body-three 3-2-1 sequence on page 423
            // of Spacecraft Dynamics.
            const Mat33 R_JiJ
                ( cPhi*cTheta , -sPhi*cPsi+cPhi*sTheta*sPsi , sPhi*sPsi+cPhi*sTheta*cPsi,
                  sPhi*cTheta ,  cPhi*cPsi+sPhi*sTheta*sPsi ,-cPhi*sPsi+sPhi*sTheta*cPsi,
                 -sTheta      ,  cTheta*sPsi                , cTheta*cPsi               );
            R_PB = MatRotation::trustMe(R_JiJ); // because P=Ji and B=J for this kind of joint
        } else {
            const Real q00=q[0]*q[0], q11=q[1]*q[1], q22=q[2]*q[2], q33=q[3]*q[3];
            const Real q01=q[0]*q[1], q02=q[0]*q[2], q03=q[0]*q[3];
            const Real q12=q[1]*q[2], q13=q[1]*q[3], q23=q[2]*q[3];
            const Mat33 R_JiJ  //rotation matrix - active-sense coordinates
                (q00+q11-q22-q33,   2*(q12-q03)  ,   2*(q13+q02),
                   2*(q12+q03)  , q00-q11+q22-q33,   2*(q23-q01),
                   2*(q13-q02)  ,   2*(q23+q01)  , q00-q11-q22+q33);
            R_PB = MatRotation::trustMe(R_JiJ); // see above
        }
    }

    void enforceBallConstraints(int offset, Vector& posv, Vector& velv) {
        if ( !useEuler ) {
            q  = Vec4::getAs(&posv[offset]);
            dq = Vec4::getAs(&velv[offset]);

            q  /= q.norm();     // Normalize Euler parameters at each time step.
            dq -= dot(q,dq)*q; // Also fix velocity: error is prop. to position component.

            Vec4::updAs(&posv[offset]) =  q;
            Vec4::updAs(&velv[offset]) = dq;
        }
    }


    void getBallInternalForce(const Vec3& forceInternal, int offset, Vector& v) const {
        //dependency: calcR_PB must be called first
        assert( useEuler );

        Vec3 torque = forceInternal;
        const Mat33 M( 0.          , 0.          , 1.    ,
                      -sPhi        , cPhi        , 0.    ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta );
        Vec3 eTorque = RigidBodyNode::DEG2RAD * M * torque;

        Vec3::updAs(&v[offset]) = eTorque;
    }

    void setBallDerivs(const Vec3& omega) {
        assert( !useEuler );
        const Mat43 M(-q(1),-q(2),-q(3),
                       q(0), q(3),-q(2),
                      -q(3), q(0), q(1),
                       q(2),-q(1), q(0));
        dq = 0.5*M*omega;
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
                  int&                  nextStateOffset,
                  bool                  useEuler)
      : RigidBodyNodeSpec<3>(mProps_B,Frame(),nextStateOffset),
        ball(nextStateOffset,useEuler)
    {
    }
    
    int  getDim() const { return ball.getBallDim(); } 

    void setJointPos(const Vector& posv) {
        ball.setBallPos(stateOffset, posv, theta);
    } 

    void getPos(Vector& posv) const {
        ball.getBallPos(theta, stateOffset, posv);
    }

    // setPos must have been called previously
    void setJointVel(const Vector& velv) {
        ball.setBallVel(stateOffset, velv, dTheta);
    }

    void getVel(Vector& velv) const {
        ball.getBallVel(dTheta, stateOffset, velv);
    }

    void calcJointAccel() {
        ball.calcBallAccel(dTheta, ddTheta);
    }

    void getAccel(Vector& accv) const {
        ball.getBallAccel(ddTheta, stateOffset, accv);
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P; // ball joint can't move B origin in P
        ball.calcR_PB(theta, R_PB);
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
        ball.enforceBallConstraints(stateOffset, posv, velv);
    }

    void getInternalForce(Vector& v) const {
        ball.getBallInternalForce(forceInternal, stateOffset, v);
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
                           int&                  nextStateOffset,
                           bool                  useEuler)
      : RigidBodyNodeSpec<6>(mProps_B,Frame(),nextStateOffset),
        ball(nextStateOffset,useEuler)
    {
    }
    
    int  getDim() const { return ball.getBallDim() + 3; } 

    void setJointPos(const Vector& posv) {
        Vec3 th;
        ball.setBallPos(stateOffset, posv, th);
        theta.updSubVec<3>(0) = th;
        theta.updSubVec<3>(3) = Vec3::getAs(&posv[stateOffset+ball.getBallDim()]);
    } 

    void getPos(Vector& posv) const {
        ball.getBallPos(theta.getSubVec<3>(0), stateOffset, posv);
        Vec3::updAs(&posv[stateOffset+ball.getBallDim()]) 
            = theta.getSubVec<3>(3);
    }

    // setPos must have been called previously
    void setJointVel(const Vector& velv) {
        Vec3 dTh;
        ball.setBallVel(stateOffset, velv, dTh);
        dTheta.updSubVec<3>(0) = dTh;
        dTheta.updSubVec<3>(3) = Vec3::getAs(&velv[stateOffset+ball.getBallDim()]);
    }

    void getVel(Vector& velv) const {
        ball.getBallVel(dTheta.getSubVec<3>(0), stateOffset, velv);
        Vec3::updAs(&velv[stateOffset+ball.getBallDim()]) 
            = dTheta.getSubVec<3>(3);
    }

    void calcJointAccel() {
        // get angular vel/accel in the space-fixed frame
        const Vec3 omega  =  dTheta.getSubVec<3>(0);
        const Vec3 dOmega = ddTheta.getSubVec<3>(0);
        ball.calcBallAccel(omega, dOmega);
    }

    void getAccel(Vector& accv) const {
        ball.getBallAccel(ddTheta.getSubVec<3>(0), stateOffset, accv);
        Vec3::updAs(&accv[stateOffset+ball.getBallDim()]) 
            = ddTheta.getSubVec<3>(3);
    }

    void calcJointKinematicsPos() {
        OB_P = refOrigin_P + theta.getSubVec<3>(3);
        ball.calcR_PB(theta.getSubVec<3>(0), R_PB);

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
        ball.enforceBallConstraints(stateOffset, posv, velv);
    }


    void getInternalForce(Vector& v) const {
        const Vec3 torque = forceInternal.getSubVec<3>(0);
        ball.getBallInternalForce(torque, stateOffset, v);
        Vec3::updAs(&v[stateOffset+ball.getBallDim()]) =
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
                  const Frame&          jointFrame,
                  int&                  nextStateOffset)
      : RigidBodyNodeSpec<2>(mProps_B,jointFrame,nextStateOffset)
    {
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P; // no translation with this joint
        calcR_PB();
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

private:
    void calcR_PB() { 
        //double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0; 
        double sinPhi = sin( scale * theta(0) );
        double cosPhi = cos( scale * theta(0) );
        double sinPsi = sin( scale * theta(1) );
        double cosPsi = cos( scale * theta(1) );

        const Mat33 a //Ry(psi) * Rx(phi)
            (cosPsi , sinPsi*sinPhi , sinPsi*cosPhi,
             0      , cosPhi        , -sinPhi      ,
            -sinPsi , cosPsi*sinPhi , cosPsi*cosPhi);

        R_PB = R_BJ * MatRotation::trustMe(a) * ~R_BJ;
    }

    void calcH() {
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;
        const MatRotation tmpR_GB = getR_GP() * R_PB;

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
                           const Frame&          jointFrame,
                           int&                  nextStateOffset)
      : RigidBodyNodeSpec<5>(mProps_B,jointFrame,nextStateOffset)
    {
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P + theta.getSubVec<3>(2);
        calcR_PB();
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = ~H * dTheta;
    }

private:
    void calcR_PB() { 
        // double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0;
        double sinPhi = sin( scale * theta(0) );
        double cosPhi = cos( scale * theta(0) );
        double sinPsi = sin( scale * theta(1) );
        double cosPsi = cos( scale * theta(1) );

        // space (parent)-fixed 1-2-3 sequence (rotation 3=0)
        const Mat33 R_JiJ  //Ry(psi) * Rx(phi)
            (cosPsi , sinPsi*sinPhi , sinPsi*cosPhi,
             0      , cosPhi        , -sinPhi      ,
            -sinPsi , cosPsi*sinPhi , cosPsi*cosPhi);

        // calculates R0*a*R0'  (R0=R_BJ(==R_PJi), a=R_JiJ)
        R_PB = R_BJ * MatRotation::trustMe(R_JiJ) * ~R_BJ; // orientation of B in parent P
    }

    void calcH() {
        //double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;
        const MatRotation& R_GP = getR_GP();
        const MatRotation tmpR_GB = R_GP * R_PB;

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
    const Frame&          jointFrame,   // inboard joint frame J in body frame
    Joint::JointType      type,
    bool                  isReversed,
    bool                  useEuler,
    int&                  nxtStateOffset)   // child-to-parent orientation?
{
    assert(!isReversed);

    switch(type) {
    case Joint::ThisIsGround:
        return new RBGroundBody();
    case Joint::Torsion:
        return new RBNodeTorsion(m,jointFrame,nxtStateOffset);
    case Joint::Universal:        
        return new RBNodeRotate2(m,jointFrame,nxtStateOffset);
    case Joint::Orientation:
        return new RBNodeRotate3(m,nxtStateOffset,useEuler);
    case Joint::Cartesian:
        return new RBNodeTranslate(m,nxtStateOffset);
    case Joint::FreeLine:
        return new RBNodeTranslateRotate2(m,jointFrame,nxtStateOffset);
    case Joint::Free:
        return new RBNodeTranslateRotate3(m,nxtStateOffset,useEuler);
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
             << ": pos: " << OB_G << ' ' << '\n';
    if (verbose&InternalDynamics::printNodeTheta) 
        cout << setprecision(8)
             << ": theta: " 
             << theta << ' ' << dTheta  << ' ' << ddTheta  << '\n';
}

