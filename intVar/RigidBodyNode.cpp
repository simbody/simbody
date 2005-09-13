/**@file
 * This file contains all the multibody mechanics code that involves a single body and
 * its inboard joint, that is, one node in the multibody tree.
 *
 * Most methods here expect to be called in a particular order during traversal of the
 * tree -- either base to tip or tip to base.
 */

#include "RigidBodyNode.h"
#include "dinternal.h"
#include "dint-step.h"
#include "vec4.h"

#include "cdsMath.h"
#include "cdsVector.h"
#include "fixedVector.h"
#include "subVector.h"
#include "subMatrix.h"
#include "matrixTools.h"
#include "cdsAuto_ptr.h"

#include "cdsIomanip.h"

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
using MatrixTools::inverse;
using CDSMath::sq;
#endif /* USE_CDS_NAMESPACE */

// NOTE: the current mechanism of chosing between MatrixTools::transpose
// and that defined in phiMatrix.hh is ugly, but seems required is order
// for gcc-2.95 to compile this source.

typedef FixedMatrix<double,2,3> Mat23;
typedef FixedVector<double,5>   Vec5;

static Mat23 catRow23(const Vec3& v1, const Vec3& v2);
static Mat33 makeJointFrameFromZAxis(const Vec3& zVec);
static Mat33 makeIdentity33();
static const Mat33 ident33 = makeIdentity33(); // handy to have around
static const Mat33 zero33(0.);

//////////////////////////////////////////////
// Implementation of RigidBodyNode methods. //
//////////////////////////////////////////////

void RigidBodyNode::addChild(RigidBodyNode* child) {
    children.append( child );
    child->setParent(this);
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

    inertia_OB_G = orthoTransform(getInertia_OB_B(), getR_GB());
    COMstation_G = getR_GB()*getCOM_B();

    // Calc Mk: the spatial inertia matrix about the body origin.
    // Note that this is symmetric; offDiag is *skew* symmetric so
    // that transpose(offDiag) = -offDiag.
    const Mat33 offDiag = getMass()*crossMat(COMstation_G);
    Mk = blockMat22( inertia_OB_G , offDiag ,
                    -offDiag      , getMass()*ident33 );
}

// Calculate velocity-related quantities: spatial velocity (sVel),
// gyroscopic force b, coriolis acceleration a. This must be
// called base to tip: depends on parent's sVel, V_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel() {
    setSpatialVel(transpose(phi) * parent->getSpatialVel()
                  + V_PB_G);
    const Vec3& omega = getSpatialAngVel();
    const Vec3 gMoment = cross(omega, inertia_OB_G * omega);
    const Vec3 gForce  = getMass() * cross(omega, 
                                           cross(omega, COMstation_G));
    b = blockVec(gMoment, gForce);

    const Vec3& vel    = getSpatialLinVel();
    const Vec3& pOmega = parent->getSpatialAngVel();
    const Vec3& pVel   = parent->getSpatialLinVel();

    // calc a: coriolis acceleration
    a = blockMat22(crossMat(pOmega),   Mat33(0.0),
                      Mat33(0.0)   , crossMat(pOmega)) 
        * V_PB_G;
    a += blockVec(Vec3(0.0), cross(pOmega, vel-pVel));
}

double
RigidBodyNode::calcKineticEnergy() const {
    double ret = dot(sVel , Mk*sVel);
    return 0.5*ret;
}

ostream& 
operator<<(ostream& s, const RigidBodyNode& node) {
    s << "RIGID BODY NODE: '<<' not impl yet!";
    //for (int i=0 ; i<node.atoms.size() ; i++)
    //    s << node.atoms[i] << ", ";
    return s;
}


////////////////////////////////////////////////
// Define classes derived from RigidBodyNode. //
////////////////////////////////////////////////

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class GroundBody : public RigidBodyNode {
public:
    GroundBody(const RigidBodyNode* node)
      : RigidBodyNode(*node) {}
    ~GroundBody() {}

    /*virtual*/const char* type() { return "ground"; }

    /*virtual*/void calcP() {} 
    /*virtual*/void calcZ(const Vec6&) {} 
    /*virtual*/void calcY() {}
    /*virtual*/void calcInternalForce(const Vec6&) {}
    /*virtual*/void calcAccel() {}
    /*virtual*/void prepareVelInternal() {}
    /*virtual*/void propagateSVel(const Vec6&) {}

    /*virtual*/void setPos(const RVec&) {}
    /*virtual*/void setVel(const RVec&) {}
    /*virtual*/void setVelFromSVel(const Vec6&) {}
    /*virtual*/void enforceConstraints(RVec& pos, RVec& vel) {}

    /*virtual*/void getPos(RVec&)   const {}
    /*virtual*/void getVel(RVec&)   const {}
    /*virtual*/void getAccel(RVec&) const {}

    /*virtual*/void getInternalForce(RVec&) const {}
    // virtual RMat getH()

    void print(int) {}
};

template<int dof>
class RigidBodyNodeSpec : public RigidBodyNode {
public:

    // We're presented with a node in the reference configuration. Derive appropriate
    // quantities for the node, such as the atom stations on the body, and the reference
    // locations and orientations.
    //
    RigidBodyNodeSpec(const RigidBodyNode* node, int& cnt, const Mat33& rotBJ)
      : RigidBodyNode(*node), theta(0.), dTheta(0.), forceInternal(0.0) 
    { 
        stateOffset = cnt;
        cnt+=dof;   // leave room for this node's state variables
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
    void setPos(const RVec& posv) {
        forceInternal.set(0.);  // forget these
        setJointPos(posv);
        calcJointKinematicsPos();
        calcJointIndependentKinematicsPos();
    }

    /// Set new velocities for the current configuration, and calculate
    /// all the velocity-dependent terms.
    void setVel(const RVec& velv) {
        // anything to initialize?
        setJointVel(velv);
        calcJointKinematicsVel();
        calcJointIndependentKinematicsVel();
    }

    // These unfortunately need to be overridden for joints using quaternions.
    virtual void setJointPos(const RVec& posv) {
        theta  = ConstRSubVec(posv,stateOffset,dof).vector();
    }
    virtual void setJointVel(const RVec& velv) {
        dTheta = ConstRSubVec(velv,stateOffset,dof).vector();
    }
    virtual void getPos(RVec& p) const {
        RSubVec(p,stateOffset,dof) = theta.vector();
    }
    virtual void getVel(RVec& v) const {
        RSubVec(v,stateOffset,dof) = dTheta.vector();
    }   
    virtual void getAccel(RVec& a) const {
        RSubVec(a,stateOffset,dof) = ddTheta.vector();
    }
    virtual void getInternalForce(RVec& t) const {
        RSubVec(t,stateOffset,dof) = forceInternal.vector();
    }

    int          getDOF() const { return dof; }
    virtual int  getDim() const { return dof; } // dim can be larger than dof

    virtual void   print(int) const;


    virtual void setVelFromSVel(const Vec6&);
    virtual void enforceConstraints(RVec& pos, RVec& vel) {}

    virtual RMat getH() const { return RMat(H); }

    void calcP();
    void calcZ(const Vec6& spatialForce);
    void calcY();
    void calcAccel();
    void calcInternalForce(const Vec6& spatialForce);
    void prepareVelInternal();
    void propagateSVel(const Vec6& desiredVel);
protected:
    // These are the joint-specific quantities
    //      ... position level
    FixedVector<double,dof>     theta;   // internal coordinates
    FixedMatrix<double,dof,6>   H;       // joint transition matrix (spatial)
    FixedMatrix<double,dof,dof> DI;
    FixedMatrix<double,6,dof>   G;

    //      ... velocity level
    FixedVector<double,dof>     dTheta;  // internal coordinate time derivatives

    //      ... acceleration level
    FixedVector<double,dof>     ddTheta; // - from the eq. of motion

    FixedVector<double,dof>     nu;
    FixedVector<double,dof>     epsilon;
    FixedVector<double,dof>     forceInternal;

private:
    void calcD_G(const Mat66& P);
};

/*static*/const double RigidBodyNode::DEG2RAD = PI / 180.;
//const double RigidBodyNode::DEG2RAD = 1.0;  //always use radians


//////////////////////////////////////////
// Derived classes for each joint type. //
//////////////////////////////////////////

typedef SubVector<RVec>       RSubVec;
typedef SubVector<const RVec> ConstRSubVec;
typedef SubVector<Vec4>       RSubVec4;
typedef SubVector<Vec5>       RSubVec5;
typedef SubVector<Vec6>       RSubVec6;

/**
 * Translate (Cartesian) joint. This provides three degrees of translational freedom
 * which is suitable (e.g.) for connecting a free atom to ground.
 * The joint frame J is aligned with the body frame B.
 */
class RBNodeTranslate : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "translate"; }

    RBNodeTranslate(const RigidBodyNode* node, int& cnt)
      : RigidBodyNodeSpec<3>(node,cnt,ident33)
    { }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P + theta;
        R_PB = ident33; // Cartesian joint can't change orientation

        // Note that this is spatial (and R_GP=R_GB for this joint)
        H = blockMat12(zero33, MatrixTools::transpose(getR_GP()));
    }

    void calcJointKinematicsVel() { 
        V_PB_G = MatrixTools::transpose(H) * dTheta;
    }
};

/**
 * This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
 * about a particular axis.
 */
class RBNodeTorsion : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "torsion"; }

    RBNodeTorsion(const RigidBodyNode*   node,
                  const Vec3&            rotDir,
                  int&                   cnt)
      : RigidBodyNodeSpec<1>(node,cnt,makeJointFrameFromZAxis(rotDir))
    { 
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P; // torsion joint can't move B origin in P
        calcR_PB();
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = MatrixTools::transpose(H) * dTheta;
    }

private:
    void calcR_PB() {
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0;
        double sinTau = sin( scale * theta(0) );
        double cosTau = cos( scale * theta(0) );
        double a[] = { cosTau , -sinTau , 0.0 ,
                       sinTau ,  cosTau , 0.0 ,
                       0.0    ,  0.0    , 1.0 };
        const Mat33 R_JiJ(a); //rotation about z-axis

        // We need R_PB=R_PJi*R_JiJ*R_JB. But R_PJi==R_BJ, so this works:
        R_PB = orthoTransform( R_JiJ , R_BJ );
    };

    // Calc H matrix in space-fixed coords.
    void calcH() {
        // This only works because the joint z axis is the same in B & P
        // because that's what we rotate around.
        const Vec3 z = getR_GP() * (R_BJ * Vec3(0,0,1)); // R_BJ=R_PJi
        FixedMatrix<double,1,3> zMat(0.0);
        H = blockMat12( FixedMatrix<double,1,3>(z.getData()) , zMat);
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

    void setBallPos(int stateOffset, const RVec& posv, FixedVector<double,3>& theta) {
        if (useEuler)
            theta = ConstRSubVec(posv,stateOffset,3).vector(); 
        else //integration step using quaternion
            q     = ConstRSubVec(posv,stateOffset,4).vector();
    } 

    void getBallPos(const Vec3& theta, int stateOffset, RVec& v) {
        if (useEuler) 
            RSubVec(v,stateOffset,3) = theta.vector();
        else  //integration step
            RSubVec(v,stateOffset,4) = q.vector();
    }

    void setBallVel(int stateOffset, const RVec& velv, FixedVector<double,3>& dTheta) { 
        assert( !useEuler );

        dq = ConstRSubVec(velv,stateOffset,4).vector();
        double a2[] = {-q(1), q(0),-q(3), q(2),
                       -q(2), q(3), q(0),-q(1),
                       -q(3),-q(2), q(1), q(0)};
        FixedMatrix<double,3,4> M(a2);
        dTheta = 2.0*( M * dq );
    }

    void getBallVel(int stateOffset, RVec& velv) {
        assert( !useEuler );
        RSubVec(velv,stateOffset,4) = dq.vector();
    }

    void calcR_PB(const Vec3& theta, Mat33& R_PB) {
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
            double R_JiJ[] = 
                { cPhi*cTheta , -sPhi*cPsi+cPhi*sTheta*sPsi , sPhi*sPsi+cPhi*sTheta*cPsi,
                  sPhi*cTheta ,  cPhi*cPsi+sPhi*sTheta*sPsi ,-cPhi*sPsi+sPhi*sTheta*cPsi,
                 -sTheta      ,  cTheta*sPsi                , cTheta*cPsi               };
            R_PB = Mat33(R_JiJ); // because P=Ji and B=J for this kind of joint
        } else {
            double R_JiJ[] =  //rotation matrix - active-sense coordinates
                {sq(q(0))+sq(q(1))-
                 sq(q(2))-sq(q(3))      , 2*(q(1)*q(2)-q(0)*q(3)), 2*(q(1)*q(3)+q(0)*q(2)),
                 2*(q(1)*q(2)+q(0)*q(3)),(sq(q(0))-sq(q(1))+
                                          sq(q(2))-sq(q(3)))     , 2*(q(2)*q(3)-q(0)*q(1)),
                 2*(q(1)*q(3)-q(0)*q(2)), 2*(q(2)*q(3)+q(0)*q(1)), (sq(q(0))-sq(q(1))-
                                                                    sq(q(2))+sq(q(3)))};
            R_PB = Mat33(R_JiJ); // see above
        }
    }

    void enforceBallConstraints(int offset, RVec& posv, RVec& velv) {
        if ( !useEuler ) {
            q  = RSubVec(posv,offset,4).vector();
            dq = RSubVec(velv,offset,4).vector();

            q  /= norm(q);     // Normalize Euler parameters at each time step.
            dq -= dot(q,dq)*q; // Also fix velocity: error is prop. to position component.

            RSubVec(posv,offset,4) =  q.vector();
            RSubVec(velv,offset,4) = dq.vector();
        }
    }

    void getBallAccel(const Vec3& omega, const Vec3& dOmega, 
                      int offset, RVec& acc)
    {
        // called after calcAccel
        assert( !useEuler );

        double a1[] = {-q(1),-q(2),-q(3),
                        q(0), q(3),-q(2),
                       -q(3), q(0), q(1),
                        q(2),-q(1), q(0)};
        FixedMatrix<double,4,3> M(a1);
        double a2[] = {-dq(1),-dq(2),-dq(3),
                        dq(0), dq(3),-dq(2),
                       -dq(3), dq(0), dq(1),
                        dq(2),-dq(1), dq(0)};
        FixedMatrix<double,4,3> dM( a2 );
        ddq = 0.5*(dM*omega + M*dOmega);

        RSubVec(acc,offset,4) = ddq.vector();
    }

    void getBallInternalForce(const Vec3& forceInternal, int offset, RVec& v) {
        //dependency: calcR_PB must be called first
        assert( useEuler );

        Vec3 torque = forceInternal;
        double a[] = { 0           , 0           , 1.0   ,
                      -sPhi        , cPhi        , 0     ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta };
        Mat33 M(a);
        Vec3 eTorque = RigidBodyNode::DEG2RAD * M * torque;

        RSubVec(v,offset,3) = eTorque.vector();
    }

    void setBallDerivs(const Vec3& dTheta) {
        assert( !useEuler );
        Vec3 omega  = dTheta;
        double a[] = {-q(1),-q(2),-q(3),
                       q(0), q(3),-q(2),
                      -q(3), q(0), q(1),
                       q(2),-q(1), q(0)};
        FixedMatrix<double,4,3> M(a);
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

    RBNodeRotate3(const RigidBodyNode* node,
                  int&                 cnt,
                  bool                 useEuler)
      : RigidBodyNodeSpec<3>(node,cnt,ident33),
        ball(cnt,useEuler)
    {}
    
    int  getDim() const { return ball.getBallDim(); } 

    void setJointPos(const RVec& posv) {
        ball.setBallPos(stateOffset, posv, theta);
    } 

    void getPos(RVec& posv) {
        ball.getBallPos(theta, stateOffset, posv);
    }

    // setPos must have been called previously
    void setJointVel(const RVec& velv) {
        ball.setBallVel(stateOffset, velv, dTheta);
    }

    void getVel(RVec& velv) {
        ball.getBallVel(stateOffset, velv);
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P; // ball joint can't move B origin in P
        ball.calcR_PB(theta, R_PB);
        // H matrix in space-fixed (P) coords
        H = blockMat12(MatrixTools::transpose(getR_GP()), zero33);
    }

    // Note that dTheta = w_PB_P = ang vel of B in P, expr in P
    void calcJointKinematicsVel() { 
        V_PB_G = MatrixTools::transpose(H) * dTheta;
    }

    void enforceConstraints(RVec& posv, RVec& velv) {
        ball.enforceBallConstraints(stateOffset, posv, velv);
    }

    void getAccel(RVec& acc) {
        ball.getBallAccel(dTheta, ddTheta, stateOffset, acc);
    }

    void getInternalForce(RVec& v) {
        ball.getBallInternalForce(forceInternal, stateOffset, v);
    }

    void setVelFromSVel(const Vec6& sVel) {
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

    RBNodeTranslateRotate3(const RigidBodyNode* node,
                           int&                 cnt,
                           bool                 useEuler)
      : RigidBodyNodeSpec<6>(node,cnt,ident33),
        ball(cnt, useEuler)
    { }
    
    int  getDim() const { return ball.getBallDim() + 3; } 

    void setJointPos(const RVec& posv) {
        Vec3 th;
        ball.setBallPos(stateOffset, posv, th);
        RSubVec6(theta,0,3) = th.vector();
        RSubVec6(theta,3,3) = 
            ConstRSubVec(posv,stateOffset+ball.getBallDim(),3).vector();
    } 

    void getPos(RVec& posv) {
        ball.getBallPos(Vec3(RSubVec6(theta,0,3).vector()), stateOffset, posv);
        RSubVec(posv,stateOffset+ball.getBallDim(),3) = RSubVec6(theta,3,3).vector();
    }

    // setPos must have been called previously
    void setJointVel(const RVec& velv) {
        Vec3 dTh;
        ball.setBallVel(stateOffset, velv, dTh);
        RSubVec6(dTheta,0,3) = dTh.vector();
        RSubVec6(dTheta,3,3) = ConstRSubVec(velv,stateOffset+ball.getBallDim(),3).vector();
    }

    void getVel(RVec& velv) {
        ball.getBallVel(stateOffset, velv);
        RSubVec(velv, stateOffset+ball.getBallDim(), 3) 
            = RSubVec6(dTheta,3,3).vector();
    }

    void calcJointKinematicsPos() {
        OB_P = refOrigin_P + Vec3(RSubVec6(theta,3,3).vector());
        ball.calcR_PB(Vec3(RSubVec6(theta,0,3).vector()), R_PB);

        // H matrix in space-fixed (P) coords
        H = blockMat22( MatrixTools::transpose(getR_GP()) , zero33,
                        zero33 , MatrixTools::transpose(getR_GP()));
    }

    // Note that dTheta[0..2] = w_PB_P = ang vel of B in P, expr in P
    void calcJointKinematicsVel() { 
        V_PB_G = MatrixTools::transpose(H) * dTheta;
    }

    void enforceConstraints(RVec& posv, RVec& velv) {
        ball.enforceBallConstraints(stateOffset, posv, velv);
    }

    void getAccel(RVec& acc) {    
        // get angular vel/accel in the space-fixed frame
        const Vec3 omega  = RSubVec6( dTheta, 0, 3).vector();
        const Vec3 dOmega = RSubVec6(ddTheta, 0, 3).vector();
        ball.getBallAccel(omega, dOmega, stateOffset, acc);
        RSubVec(acc,stateOffset+ball.getBallDim(),3) = RSubVec6(ddTheta,3,3).vector();
    }

    void getInternalForce(RVec& v) {
        const Vec3 torque = RSubVec6(forceInternal, 0, 3).vector();
        ball.getBallInternalForce(torque, stateOffset, v);
        RSubVec(v,stateOffset+ball.getBallDim(),3) = RSubVec6(forceInternal,3,3).vector();

    }

    void setVelFromSVel(const Vec6& sVel) {
        RigidBodyNodeSpec<6>::setVelFromSVel(sVel);
        const Vec3 omega  = RSubVec6( dTheta, 0, 3).vector();
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

    RBNodeRotate2(const RigidBodyNode* node,
                  const Vec3&          zVec,
                  int&                 cnt)
      : RigidBodyNodeSpec<2>(node,cnt,makeJointFrameFromZAxis(zVec))
    { 
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P; // no translation with this joint
        calcR_PB();
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = MatrixTools::transpose(H) * dTheta;
    }

private:
    void calcR_PB() { 
        //double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0; 
        double sinPhi = sin( scale * theta(0) );
        double cosPhi = cos( scale * theta(0) );
        double sinPsi = sin( scale * theta(1) );
        double cosPsi = cos( scale * theta(1) );

        double a[] =  //Ry(psi) * Rx(phi)
            {cosPsi , sinPsi*sinPhi , sinPsi*cosPhi,
             0      , cosPhi        , -sinPhi      ,
            -sinPsi , cosPsi*sinPhi , cosPsi*cosPhi};

        const Mat33 R_PB = orthoTransform( Mat33(a) , R_BJ );
    }

    void calcH() {
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;
        const Mat33 tmpR_GB = getR_GP() * R_PB;

        const Vec3 x = scale * tmpR_GB * (R_BJ * Vec3(1,0,0));
        const Vec3 y = scale * tmpR_GB * (R_BJ * Vec3(0,1,0));

        const Mat23 zMat23(0.0);
        H = blockMat12(catRow23(x,y) , zMat23);
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

    RBNodeTranslateRotate2(const RigidBodyNode*  node,
                           const Vec3&           zVec,
                           int&                  cnt)
      : RigidBodyNodeSpec<5>(node,cnt,makeJointFrameFromZAxis(zVec))
    { 
    }

    void calcJointKinematicsPos() { 
        OB_P = refOrigin_P + Vec3(RSubVec5(theta,2,3).vector());
        calcR_PB();
        calcH();
    }

    void calcJointKinematicsVel() { 
        V_PB_G = MatrixTools::transpose(H) * dTheta;
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
        double R_JiJ[] =  //Ry(psi) * Rx(phi)
            {cosPsi , sinPsi*sinPhi , sinPsi*cosPhi,
             0      , cosPhi        , -sinPhi      ,
            -sinPsi , cosPsi*sinPhi , cosPsi*cosPhi};

        // calculates R0*a*R0'  (R0=R_BJ(==R_PJi), a=R_JiJ)
        const Mat33 R_PB = orthoTransform( Mat33(R_JiJ) , R_BJ ); // orientation of B in parent P
    }

    void calcH() {
        //double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;
        const Mat33 tmpR_GB = getR_GP() * R_PB;

        Vec3 x = scale * tmpR_GB * (R_BJ * Vec3(1,0,0));
        Vec3 y = scale * tmpR_GB * (R_BJ * Vec3(0,1,0));

        const Mat23 zMat23(0.0);
        H = blockMat22(catRow23(x,y) , zMat23 ,
                        zero33       , MatrixTools::transpose(getR_GP()));
    }  
};


/////////////////////////////////////////////////////////////
// Implementation of RigidBodyNodeSpec base class methods. //
/////////////////////////////////////////////////////////////


//
// to be called from base to tip.
//
template<int dof> void
RigidBodyNodeSpec<dof>::setVelFromSVel(const Vec6& sVel) {
    dTheta = H * (sVel - transpose(phi)*parent->sVel);
}

template<int dof> void
RigidBodyNodeSpec<dof>::calcD_G(const Mat66& P) {
    using InternalDynamics::Exception;
    FixedMatrix<double,dof,dof> D = orthoTransform(P,H);
    try {
        DI = inverse(D);
    }
    catch ( SingularError ) {
        cerr << "calcD_G: singular D matrix: " << D << '\n'
             << "H matrix: " << H << '\n'
             << "node level: " << level << '\n'
             << "number of children: " << children.size() << '\n'
             << endl;
        throw Exception("calcD_G: singular D matrix. Bad topology?");
    }
    G = P * MatrixTools::transpose(H) * DI;
}

//
// Calculate Pk and related quantities. The requires that the children
// of the node have already had their quantities calculated.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcP() {
    //
    //how much do we need to keep around?
    // it looks like nu and G are the only ones needed for the acceleration
    // calc. The others can be freed after the parent is done with them.
    //
    P = Mk;

    SubMatrix<Mat66> p11(P,0,0,3,3);
    SubMatrix<Mat66> p12(P,0,3,3,3);
    SubMatrix<Mat66> p21(P,3,0,3,3);
    SubMatrix<Mat66> p22(P,3,3,3,3);
    for (int i=0 ; i<children.size() ; i++) {
        // this version is readable
        // P += orthoTransform( children[i]->tau * children[i]->P ,
        //                      transpose(children[i]->phiT) );
        // this version is not
        Mat33 lt = crossMat(children[i]->getOB_G() - getOB_G());
        Mat66 M  = children[i]->tau * children[i]->P;
        SubMatrix<Mat66> m11(M,0,0,3,3);
        SubMatrix<Mat66> m12(M,0,3,3,3);
        SubMatrix<Mat66> m21(M,3,0,3,3);
        SubMatrix<Mat66> m22(M,3,3,3,3);
        p11 += m11+lt*m21-m12*lt-lt*m22*lt;
        p12 += m12+lt*m22;
        p21 += m21-m22*lt;
        p22 += m22;
    }

    calcD_G(P);
    tau.set(0.0); tau.setDiag(1.0);
    tau -= G * H;
    // calcZ();
}
 
//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcZ(const Vec6& spatialForce) {
    psiT = MatrixTools::transpose(tau) * transpose(phi);
    z = P * a + b + spatialForce;

    for (int i=0 ; i<children.size() ; i++) 
        z += children[i]->phi
             * (children[i]->z + children[i]->Gepsilon);

    epsilon  = forceInternal - H*z;
    nu       = DI * epsilon;
    Gepsilon = G * epsilon;
}

//
// Calculate acceleration in internal coordinates.
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcAccel() {
    // const Node* pNode = parentHinge.remNode;
    //make sure that this is phi is correct - FIX ME!
    // base alpha = 0!!!!!!!!!!!!!
    Vec6 alphap = transpose(phi) * parent->sAcc;
    ddTheta = nu - MatrixTools::transpose(G) * alphap;

    sAcc   = alphap + MatrixTools::transpose(H) * ddTheta + a;  
}

// To be called base to tip.
template<int dof> void
RigidBodyNodeSpec<dof>::calcY() {
    Y = orthoTransform(DI,MatrixTools::transpose(H))
        + orthoTransform(parent->Y,psiT);
}

//
// Calculate sum of internal force and effective forces due to Cartesian
// forces.
// To be called from tip to base.
// Should be called only once after calcProps.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcInternalForce(const Vec6& spatialForce) {
    z = spatialForce;

    for (int i=0 ; i<children.size() ; i++) 
        z += children[i]->phi * children[i]->z;

    forceInternal += H * z; 
}

//
// Set up various variables so that calling calcP followed by
// calling getAcc returns the time-derivative of the internal coordinates.
// This requires that sVel is previously set.
//
template<int dof> void
RigidBodyNodeSpec<dof>::prepareVelInternal() {
    dTheta.set(0.0);
    a.set(0.0);
    b.set(0.0);
    forceInternal.set(0.0);

    // Note: Schwieters was computing a new spatial mass matrix here, according
    // to Eqn. 80a in his paper. But as far as I can tell this is exactly the
    // same as the ordinary mass matrix. (sherm 050906)
    // Schwieters' comment: "set Mk = R * M * R^T"
}

//
// Set up various variables so that calling calcP followed by
// calling getAcc returns the internal velocities which best
// approximate a set of desired velocities. See Section 3.1.3 in
// Schwieters' paper: we're solving Eqn. 78. The passed-in
// velocity is a mass weighted combination of the (unachievable)
// initial atom velocities on this body. This is combined with
// outboard desired velocities, and then the result is treated
// as though it were a force. (sherm 050906)
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::propagateSVel(const Vec6& desiredVel) {
    //sAcc used as a temporary
    sVel = desiredVel; // ???
    sAcc = desiredVel;
    for (int i=0 ; i<children.size() ; i++)
        sAcc += children[i]->phi * children[i]->sAcc;
    forceInternal = H * sAcc;
}



template<int dof> void
RigidBodyNodeSpec<dof>::print(int verbose) const {
    if (verbose&InternalDynamics::printNodeForce) 
        cout << setprecision(8)
             << ": force: " << forceCartesian << '\n';
    if (verbose&InternalDynamics::printNodePos) 
        cout << setprecision(8)
             << ": pos: " << OB_G << ' ' << '\n';
    if (verbose&InternalDynamics::printNodeTheta) 
        cout << setprecision(8)
             << ": theta: " 
             << theta << ' ' << dTheta  << ' ' << ddTheta  << '\n';
}

/////////////////////////////////////
// Miscellaneous utility routines. //
/////////////////////////////////////

static Mat23
catRow23(const Vec3& v1, const Vec3& v2) {
    FixedMatrix<double,1,3> m1(v1.getData());
    FixedMatrix<double,1,3> m2(v2.getData());
    Mat23 ret = blockMat21(m1,m2);
    return ret;
}

// Calculate a rotation matrix R_BJ which defines the J
// frame by taking the B frame z axis into alignment 
// with the passed-in zDir vector. This is not unique.
// notes of 12/6/99 - CDS
static Mat33
makeJointFrameFromZAxis(const Vec3& zVec) {
    const Vec3 zDir = unitVec(zVec);

    // Calculate spherical coordinates.
    double theta = acos( zDir.z() );             // zenith (90-elevation)
    double psi   = atan2( zDir.x() , zDir.y() ); // 90-azimuth

    // This is a space fixed 1-2-3 sequence with angles
    // a1=-theta, a2=0, a3=-psi. That is, to get from B to J
    // first rotate by -theta around the B frame x axis, 
    // then rotate by -psi around the B frame z axis. (sherm)

    const double R_BJ[] = 
        { cos(psi) , cos(theta)*sin(psi) , sin(psi)*sin(theta),
         -sin(psi) , cos(theta)*cos(psi) , cos(psi)*sin(theta),
          0        , -sin(theta)         , cos(theta)         };
    return Mat33(R_BJ); // == R_PJi
}

static Mat33
makeIdentity33() {
    Mat33 ret(0.);
    ret.setDiag(1.);
    return ret;
}

