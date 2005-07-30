
#include "dint-node.h"
#include "dinternal.h"
#include "dint-step.h"
#include "dint-atom.h"
#include "vec4.h"

#include <cdsMath.h>
#include <cdsVector.h>
#include <fixedVector.h>
#include <subVector.h>
#include <subMatrix.h>
#include <matrixTools.h>
#include <cdsAuto_ptr.h>

#include <cdsIomanip.h>

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

static Mat23
catRow(const Vec3& v1, const Vec3& v2) {
    FixedMatrix<double,1,3> m1(v1.getData());
    FixedMatrix<double,1,3> m2(v2.getData());
    Mat23 ret = blockMat21(m1,m2);
    return ret;
}

// Calculate a rotation matrix R_BJ which defines the J
// frame by taking the B frame z axis into alignment 
// with the passed-in zDir vector. This is not unique.
// notes of 12/6/99 - CDS
static Mat3
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
    return Mat3(R_BJ); // == R_PJi
}

static Mat3
makeIdentityRotation() {
    Mat3 ret(0.);
    ret.setDiag(1.);
    return ret;
}

static const Mat3 R_I = makeIdentityRotation(); // handy to have around

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class GroundBody : public HingeNode {
public: 
    virtual const char* type() { return "ground"; }

    GroundBody(const HingeNode* node)
      : HingeNode(*node)
    { 
        for (int i=0 ; i<atoms.size() ; i++) atoms[i]->vel = Vec3(0.0); 
        for (l_int i=0 ; i<atoms.size() ; i++) atoms[i]->node = this; 
    }
    virtual ~GroundBody() {}

//  virtual const Vec3& posCM() const { return posCM_; }
//  virtual const double& mass() const { return mass_; }
//  virtual void velFromCartesian();

    void calcP() {} 
    void calcZ() {} 
    void calcPandZ() {} 
    void calcAccel() {}
    void calcInternalForce() {}
    void prepareVelInternal() {}
    void propagateSVel(const Vec6&) {}
    virtual void setPosVel(const RVec&, const RVec&) {}
    virtual void setVel(const RVec&) {}
    virtual void setVelFromSVel(const Vec6&) {}
    virtual void enforceConstraints(RVec& pos, RVec& vel) {}
    virtual void print(int) {}
    virtual void getPos(RVec&) {}
    virtual void getVel(RVec&) {}
    virtual void getAccel(RVec&) {}
    virtual void getInternalForce(RVec&) {}
    virtual void calcY() {}
};

template<int dof>
class HingeNodeSpec : public HingeNode {
public:

    // We're presented with a node in the reference configuration. Derive appropriate
    // quantities for the node, such as the atom stations on the body, and the reference
    // locations and orientations.
    //
    // XXX Ball and free joints require that we insert a dummy 'atom' as the origin if
    //     this is a base node and it is not bonded to ground (unless there is already
    //     a dummy atom as atom[0]). This is chosen to be
    //     a point which is the geometric center of the entire outboard tree, but
    //     attached to the current body. At the moment I'm just preserving the
    //     existing behavior; later I hope to trash it. (sherm)
    //
    HingeNodeSpec(const HingeNode* node, int& cnt, const Mat3& rotBJ,
                  bool addDummyOrigin=false)
      : HingeNode(*node), offset_(cnt),
        mass_(0), inertia(), Mk(0.), posCM_(0.),
        a(0.), forceCartesian(0.), b(0.), R_BJ(rotBJ),
        theta(0.), dTheta(0.), forceInternal(0.0) 
    { 
        cnt+=dof;   // leave room for this node's state variables

        for (int i=0 ; i<atoms.size() ; i++)
            atoms[i]->node = this; 

        // If requested, add center of mass "atom" (station) if node is not attached to
        // the origin (ground) node.
        if (addDummyOrigin && isBaseNode() &&  atoms[0]->mass != 0.0) {

            bool bound = false;
            for (int i=0 ; i<atoms[0]->bonds.size() ; i++)
                if (atoms[0]->bonds[i]->node->isGroundNode())
                    bound = true;

            if ( !bound ) {
                IVMAtom* a = new IVMAtom(-1,0.0);
                cmAtom.reset(a);
                a->pos = AtomTree::findCM( this );
                a->node = this;
                atoms.prepend(a);
            }
        }

        // drop constness just for a moment here in the constructor ...
        Vec3&              mutableRefOrigin_P    = *const_cast<Vec3*>(&refOrigin_P);
        CDSVector<Vec3,1>& mutableAtomStations_B = *const_cast<CDSVector<Vec3,1>*>(&atomStations_B);

        mutableRefOrigin_P = atoms[0]->pos - parent->getAtom(0)->pos;
        mutableAtomStations_B.resize(atoms.size()-1);
        for (int i=1 ; i<atoms.size() ; i++)
            mutableAtomStations_B(i) = atoms[i]->pos - atoms[0]->pos;
    }

    virtual ~HingeNodeSpec() {}

    // Every derived class must implement at least these two methods.

    /// Calculate the joint transition matrix H.
    virtual void calcH()=0;

    /// Calculate all kinematic quantities.
    virtual void toCartesian()=0;

    virtual const Vec3& posCM() const { return posCM_; }
    virtual const double& mass() const { return mass_; }
    virtual int offset() const {return offset_;}

    void calcPandZ() { calcP(); calcZ(); }
    void calcP();
    void calcZ();
    void calcAccel();
    void calcInternalForce();

    virtual void calcY();
    //  void calcEpsilon_nu(const Vec6& z);
    void calcD_G(const Mat6& P);
    int  getDOF() const { return dof; }
    virtual int  getDim() const { return dof; }
    virtual void print(int);
    virtual double kineticE();
    //  virtual double approxMass(int k);
    virtual double approxKE();
    void calcCartesianForce();
    void prepareVelInternal();
    void propagateSVel(const Vec6&);
    virtual void setPosVel(const RVec&,
    const RVec&);
    virtual void setVel(const RVec&);
    virtual void setVelFromSVel(const Vec6&);

    virtual void enforceConstraints(RVec& pos, RVec& vel) {}
    virtual void getPos(RVec&);
    virtual void getVel(RVec&);
    virtual void getAccel(RVec&);
    virtual void getInternalForce(RVec&);
    virtual RMat getH() { return RMat(H); }


protected:
    int           offset_;  //index into internal coord pos,vel,acc arrays

    // These are the body properties

    double        mass_;
    InertiaTensor inertia;
    Mat6          Mk;
    Vec3          posCM_;
    Vec6          a;        //coriolis acceleration
    //gyroscopic spatial force
    Vec6 forceCartesian;
    Vec6 b;

    // This serves as the owner for an extra 'dummy' IVMAtom needed by some joints.
    CDS::auto_ptr<IVMAtom> cmAtom;

    // Atom stations. These are vectors from this body's origin point to each of
    // the atoms, expressed in the body frame. These are fixed after construction;
    // they are not state dependent.
    const CDSVector<Vec3,1> atomStations_B;

    // Inboard joint frame. This is fixed forever once constructed and gives the
    // orientation of the body-fixed J frame in the body frame B. This is an 
    // identity matrix for some joint types.
    const Mat3 R_BJ;

    // Reference configuration. This is the body frame origin location, measured
    // in its parent's frame in the reference configuration. This vector is fixed
    // after construction! The body origin can of course move relative to its
    // parent, but that is not the meaning of this reference configuration vector.
    // (Note however that the body origin is also the location of the inboard joint, 
    // meaning that the origin point moves relative to the parent only due to translations.)
    // Note that by definition the orientation of the body frame is identical to P
    // in the reference configuration so we don't need to store it.
    const Vec3 refOrigin_P;

    // These are the joint properties

    FixedVector<double,dof>     theta;   //internal coordinate
    FixedVector<double,dof>     dTheta;  //internal coordinate time derivatives
    FixedVector<double,dof>     ddTheta; // - from the eq. of motion
    FixedMatrix<double,dof,6>   H;
    FixedVector<double,dof>     nu;
    FixedVector<double,dof>     epsilon;
    FixedVector<double,dof>     forceInternal;
    FixedMatrix<double,dof,dof> DI;
    FixedMatrix<double,6,dof>   G;

    void calcProps();
};

/*static*/const double HingeNode::DEG2RAD = PI / 180.;
//const double HingeNode::DEG2RAD = 1.0;  //always use radians

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
class HNodeTranslate : public HingeNodeSpec<3> {
public:
    virtual const char* type() { return "translate"; }

    HNodeTranslate(const HingeNode* node, int& cnt)
      : HingeNodeSpec<3>(node,cnt,R_I)
    { }

    void calcH() {
        using MatrixTools::transpose;
        Mat3 zMat(0.0);
        H = blockMat12(zMat,transpose(getR_GP()));
    }

    void toCartesian() { 
        R_GB = getR_GP();   // translate joint can't change orientation
        atoms[0]->pos = parent->getAtom(0)->pos
                        + getR_GP() * (refOrigin_P + theta);
        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->pos =  atoms[0]->pos + R_GB * atomStations_B(i);
        phi = PhiMatrix( atoms[0]->pos - parent->getAtom(0)->pos );
        calcH();
        sVel = transpose(phi) * parent->getSpatialVel()
                + MatrixTools::transpose(H) * dTheta;
        atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

        for (l_int i=1; i<atoms.size(); i++)
            atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                                   atoms[i]->pos - atoms[0]->pos );
        calcProps();
    }
};

/**
 * Free joint. This is a six degree of freedom joint providing unrestricted 
 * translation and rotation for a free rigid body.
 * The joint frame J is aligned with the body frame B.
 */
class HNodeTranslateRotate3 : public HingeNodeSpec<6> {
    Vec4 q; //Euler parameters for rotation relative to parent
    Vec4 dq;
    Vec4 ddq;
    double cPhi, sPhi;        //trig functions of Euler angles
    double cPsi, sPsi;        //used for minimizations
    double cTheta, sTheta;
    bool   useEuler;          //if false, use Quaternion rep.
public:
    virtual const char* type() { return "full"; }

    HNodeTranslateRotate3(const HingeNode* node,
                          int&             cnt,
                          bool             useEuler)
      : HingeNodeSpec<6>(node,cnt,R_I,true), q(1.0,0.0,0.0,0.0), dq(0.), useEuler(useEuler)
    { 
        if ( !useEuler )
            cnt++;
    }

    void calcRot() {
        Mat3 R_PB;  // rotation matrix expressing body frame in parent's frame
        if (useEuler) {
            // theta = (Phi, Theta, Psi) Euler ``3-2-1'' body-fixed angles 
            cPhi   = cos( theta(0) *DEG2RAD );
            sPhi   = sin( theta(0) *DEG2RAD );
            cTheta = cos( theta(1) *DEG2RAD );
            sTheta = sin( theta(1) *DEG2RAD );
            cPsi   = cos( theta(2) *DEG2RAD );
            sPsi   = sin( theta(2) *DEG2RAD );

            // This is the rotation matrix giving the current orientation of the J
            // frame with respect to the parent body's Ji frame, where the current
            // body B is the i'th child of P.
            const double R_JiJ[] = 
                { cPhi*cTheta , -sPhi*cPsi+cPhi*sTheta*sPsi , sPhi*sPsi+cPhi*sTheta*cPsi,
                  sPhi*cTheta ,  cPhi*cPsi+sPhi*sTheta*sPsi ,-cPhi*sPsi+sPhi*sTheta*cPsi,
                 -sTheta      ,  cTheta*sPsi                , cTheta*cPsi               };

            R_PB = Mat3(R_JiJ); // because P=Ji and B=J for this kind of joint
        } else {
            const double R_JiJ[] =  //rotation matrix - active-sense coordinates
                {sq(q(0))+sq(q(1))-
                 sq(q(2))-sq(q(3))      , 2*(q(1)*q(2)-q(0)*q(3)), 2*(q(1)*q(3)+q(0)*q(2)),
                 2*(q(1)*q(2)+q(0)*q(3)), (sq(q(0))-sq(q(1))+
                                           sq(q(2))-sq(q(3)))    , 2*(q(2)*q(3)-q(0)*q(1)),
                 2*(q(1)*q(3)-q(0)*q(2)), 2*(q(2)*q(3)+q(0)*q(1)), (sq(q(0))-sq(q(1))-
                                                                    sq(q(2))+sq(q(3)))};
            R_PB = Mat3(R_JiJ); // see above
        }
        R_GB = getR_GP() * R_PB; // the spatial rotation matrix expressing body frame B in ground
    }

    void enforceConstraints(RVec& posv, RVec& velv) {
        if (useEuler) return;
        q  = RSubVec(posv,offset_,4).vector();
        dq = RSubVec(velv,offset_,4).vector();

        q  /= norm(q);      // normalize Euler parameters at each time step
        dq -= dot(q,dq)*q; //also fix velocity: error is prop. to position component

        RSubVec(posv,offset_,4) =  q.vector();
        RSubVec(velv,offset_,4) = dq.vector();
    }

    // H matrix in body-fixed coords
    void calcH() {
        //   calcRot(); //FIX: does this need to be calculated here?
        using MatrixTools::transpose;
        Mat3 zMat(0.0);
        H = blockMat22( transpose(getR_GP()) , zMat,
                        zMat , transpose(getR_GP()));
    }

    // called after calcAccel
    void getAccel(RVec& v) {
        assert( !useEuler );

        // get angular vel/accel in the space-fixed frame
        Vec3 omega  = RSubVec6( dTheta, 0, 3).vector();
        Vec3 dOmega = RSubVec6(ddTheta, 0, 3).vector();
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

        RSubVec(v,offset_,4)   = ddq.vector();
        RSubVec(v,offset_+4,3) = RSubVec6(ddTheta,3,3).vector();
    }

    void getInternalForce(RVec& v) {
        //dependency: calcRot must be called first
        assert( useEuler );

        Vec3 torque = RSubVec6(forceInternal, 0, 3).vector();
        double a[] = { 0         , 0           , 1.0 ,
                      -sPhi      , cPhi        , 0   ,
                      cPhi*cTheta, sPhi*cTheta ,-sTheta };
        Mat3 M(a);
        Vec3 eTorque = DEG2RAD * M * torque;

        RSubVec(v,offset_,3)   = eTorque.vector();
        RSubVec(v,offset_+3,3) = RSubVec6(forceInternal,3,3).vector();
    }

    int  getDim() const { 
        if (useEuler) return 6; else return 7; 
    }

    void setPosVel(const RVec& posv, const RVec& velv) {
        if (useEuler) {
            theta = ConstRSubVec(posv,offset_,6).vector(); 
        } else { //integration step using quaternion
            q                   = ConstRSubVec(posv,offset_,4).vector();
            RSubVec6(theta,3,3) = ConstRSubVec(posv,offset_+4,3).vector();

            dq = ConstRSubVec(velv,offset_,4).vector();
            double a2[] = {-q(1), q(0),-q(3), q(2),
                           -q(2), q(3), q(0),-q(1),
                           -q(3),-q(2), q(1), q(0)};
            FixedMatrix<double,3,4> M(a2);
            RSubVec6(dTheta,0,3) = ( 2.0*( M * dq ) ).vector();
            RSubVec6(dTheta,3,3) = ConstRSubVec(velv,offset_+4,3).vector();
        }
        toCartesian();
    } 

    void setVel(const RVec& velv) { //setposvel must have been called previously
        assert( !useEuler );

        dq = ConstRSubVec(velv,offset_,4).vector();

        double a2[] = {-q(1), q(0),-q(3), q(2),
                       -q(2), q(3), q(0),-q(1),
                       -q(3),-q(2), q(1), q(0)};
        FixedMatrix<double,3,4> M(a2);
        RSubVec6(dTheta,0,3) = ( 2.0*( M * dq ) ).vector();
        RSubVec6(dTheta,3,3) = ConstRSubVec(velv,offset_+4,3).vector();

        sVel = transpose(phi) * parent->getSpatialVel()
               + MatrixTools::transpose(H) * dTheta;
        atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                                   atoms[i]->pos-atoms[0]->pos );
    } 

    void setVelFromSVel(const Vec6& sVel) {
        assert( !useEuler );

        HingeNodeSpec<6>::setVelFromSVel(sVel);
        Vec3 omega  = RSubVec6( dTheta, 0, 3).vector();
        double a[] = {-q(1),-q(2),-q(3),
                       q(0), q(3),-q(2),
                      -q(3), q(0), q(1),
                       q(2),-q(1), q(0)};
        FixedMatrix<double,4,3> M(a);
        dq = 0.5*M*omega;
    } 

    void getPos(RVec& v) {
        if (useEuler) 
            RSubVec(v,offset_,6) = theta.vector();
        else { //integration step
            RSubVec(v,offset_  ,4) = q.vector();
            RSubVec(v,offset_+4,3) = RSubVec6(theta, 3, 3).vector();
        }
    }

    void getVel(RVec& v) {
        assert( !useEuler );

        RSubVec(v,offset_,4) = dq.vector();
        RSubVec(v,offset_+4,3) = RSubVec6(dTheta, 3, 3).vector();
    }
   
    void toCartesian() { 
        atoms[0]->pos = parent->getAtom(0)->pos
                        + getR_GP() * (refOrigin_P + 
                                Vec3(RSubVec6(theta,3,3).vector()));
        calcRot();
        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->pos = atoms[0]->pos + R_GB * atomStations_B(i);

        phi = PhiMatrix( atoms[0]->pos - parent->getAtom(0)->pos );
        calcH();
        if ( !useEuler ) {
            // update vel using recursive formula
            sVel = transpose(phi) * parent->getSpatialVel()
                   + MatrixTools::transpose(H) * dTheta;
            atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

            for (l_int i=1 ; i<atoms.size() ; i++)
                atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                                       atoms[i]->pos-atoms[0]->pos );
            //   inertia.calc( atoms[0]->pos , atoms );
        }
        calcProps();
    }

    void print(int verbose) {
        if (verbose&InternalDynamics::printNodeTheta && !useEuler) 
            cout << setprecision(8) << atoms[0] << ": theta: " 
                 << q   << ' ' << RSubVec6(theta,3,3)   << ' ' 
                 << dq  << ' ' << RSubVec6(dTheta,3,3)  << ' ' 
                 << ddq << ' ' << RSubVec6(ddTheta,3,3) << '\n';
        else 
            HingeNodeSpec<6>::print(verbose);
    }
};

/**
 * Ball joint. This provides three degrees of rotational freedom, i.e.,
 * unrestricted orientation.
 * The joint frame J is aligned with the body frame B.
 */
class HNodeRotate3 : public HingeNodeSpec<3> {
    Vec4    q; //Euler parameters for rotation relative to parent
    Vec4    dq;
    Vec4    ddq;
    double cPhi, sPhi;        //trig functions of Euler angles
    double cPsi, sPsi;        //used for minimizations
    double cTheta, sTheta;
    bool   useEuler;          //if False, use Quaternion rep.

public:
    virtual const char* type() { return "rotate3"; }

    HNodeRotate3(const HingeNode* node,
                  int&            cnt,
                  bool            useEuler)
      : HingeNodeSpec<3>(node,cnt,R_I,true), 
        q(1.0,0.0,0.0,0.0), useEuler(useEuler)
    { 
        if ( !useEuler )
            cnt++;
    }

    void calcRot() {
        Mat3 R_PB;  // rotation matrix expressing body frame in parent's frame
        if (useEuler) {
            // theta = (Phi, Theta, Psi) Euler ``3-2-1'' angles 
            cPhi   = cos( theta(0) *DEG2RAD );
            sPhi   = sin( theta(0) *DEG2RAD );
            cTheta = cos( theta(1) *DEG2RAD );
            sTheta = sin( theta(1) *DEG2RAD );
            cPsi   = cos( theta(2) *DEG2RAD );
            sPsi   = sin( theta(2) *DEG2RAD );
            
            // (sherm 050726) This matches Kane's Body-three 3-2-1 sequence on page 423
            // of Spacecraft Dynamics.
            double R_JiJ[] = 
                { cPhi*cTheta , -sPhi*cPsi+cPhi*sTheta*sPsi , sPhi*sPsi+cPhi*sTheta*cPsi,
                  sPhi*cTheta ,  cPhi*cPsi+sPhi*sTheta*sPsi ,-cPhi*sPsi+sPhi*sTheta*cPsi,
                 -sTheta      ,  cTheta*sPsi                , cTheta*cPsi               };
            R_PB = Mat3(R_JiJ); // because P=Ji and B=J for this kind of joint
        } else {
            double R_JiJ[] =  //rotation matrix - active-sense coordinates
                {sq(q(0))+sq(q(1))-
                 sq(q(2))-sq(q(3))      , 2*(q(1)*q(2)-q(0)*q(3)), 2*(q(1)*q(3)+q(0)*q(2)),
                 2*(q(1)*q(2)+q(0)*q(3)),(sq(q(0))-sq(q(1))+
                                          sq(q(2))-sq(q(3)))     , 2*(q(2)*q(3)-q(0)*q(1)),
                 2*(q(1)*q(3)-q(0)*q(2)), 2*(q(2)*q(3)+q(0)*q(1)), (sq(q(0))-sq(q(1))-
                                                                    sq(q(2))+sq(q(3)))};
            R_PB = Mat3(R_JiJ); // see above
        }
        R_GB = getR_GP() * R_PB; // R_GB = R_GP*R_PB, the spatial rotation matrix
    }

    void enforceConstraints(RVec& posv, RVec& velv) {
        if ( !useEuler ) {
            q  = RSubVec(posv,offset_,4).vector();
            dq = RSubVec(velv,offset_,4).vector();

            q  /= norm(q);     // Normalize Euler parameters at each time step.
            dq -= dot(q,dq)*q; // Also fix velocity: error is prop. to position component.

            RSubVec(posv,offset_,4) =  q.vector();
            RSubVec(velv,offset_,4) = dq.vector();
        }
    }

    void calcH() { // H matrix in body-fixed coords
        //   calcRot(); //FIX: does this need to be calculated here?
        using MatrixTools::transpose;
        Mat3 zMat(0.0);
        H = blockMat12(transpose(getR_GP()) , 
                       zMat);
    }

    void getAccel(RVec& v) {
        // called after calcAccel
        assert( !useEuler );

        // get angular vel/accel in the space-fixed frame
        Vec3 omega  = dTheta;
        Vec3 dOmega = ddTheta;
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

        RSubVec(v,offset_,4)   = ddq.vector();
    }

    void getInternalForce(RVec& v) {
        //dependency: calcRot must be called first
        assert( useEuler );

        Vec3 torque = forceInternal;
        double a[] = { 0           , 0           , 1.0   ,
                      -sPhi        , cPhi        , 0     ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta };
        Mat3 M(a);
        Vec3 eTorque = DEG2RAD * M * torque;

        RSubVec(v,offset_,3)   = eTorque.vector();
    }

    int  getDim() const { 
        if (useEuler) return 3; else return 4; 
    }

    void setPosVel(const RVec& posv, const RVec& velv) {
        if (useEuler) {
            theta = ConstRSubVec(posv,offset_,3).vector(); 
        } else { //integration step using quaternion
            q  = ConstRSubVec(posv,offset_,4).vector();

            dq = ConstRSubVec(velv,offset_,4).vector();
            double a2[] = {-q(1), q(0),-q(3), q(2),
                           -q(2), q(3), q(0),-q(1),
                           -q(3),-q(2), q(1), q(0)};
            FixedMatrix<double,3,4> M(a2);
            dTheta = 2.0*( M * dq );
        }
        toCartesian();
    } 

    void setVel(const RVec& velv) { //setposvel must have been called previously
        assert( !useEuler );

        dq = ConstRSubVec(velv,offset_,4).vector();

        double a2[] = {-q(1), q(0),-q(3), q(2),
                       -q(2), q(3), q(0),-q(1),
                       -q(3),-q(2), q(1), q(0)};
        FixedMatrix<double,3,4> M(a2);
        dTheta = 2.0*( M * dq );

        sVel = transpose(phi) * parent->getSpatialVel()
                + MatrixTools::transpose(H) * dTheta;
        atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                                   atoms[i]->pos - atoms[0]->pos );
    } 

    void setVelFromSVel(const Vec6& sVel) {
        assert( !useEuler );

        HingeNodeSpec<3>::setVelFromSVel(sVel);
        Vec3 omega  = dTheta;
        double a[] = {-q(1),-q(2),-q(3),
                       q(0), q(3),-q(2),
                      -q(3), q(0), q(1),
                       q(2),-q(1), q(0)};
        FixedMatrix<double,4,3> M(a);
        dq = 0.5*M*omega;
    } 

    void getPos(RVec& v) {
        if (useEuler) 
            RSubVec(v,offset_,3) = theta.vector();
        else  //integration step
            RSubVec(v,offset_  ,4) = q.vector();
    }

    void getVel(RVec& v) {
        assert( !useEuler );
        RSubVec(v,offset_,4) = dq.vector();
    }
   
    void toCartesian() { 
        atoms[0]->pos = parent->getAtom(0)->pos
                        + getR_GP() * refOrigin_P;
        calcRot();
        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->pos = atoms[0]->pos + R_GB * atomStations_B(i);

        phi = PhiMatrix( atoms[0]->pos - parent->getAtom(0)->pos );
        calcH();
        if ( !useEuler ) {
            // update vel using recursive formula
            sVel = transpose(phi) * parent->getSpatialVel()
                    + MatrixTools::transpose(H) * dTheta;
            atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

            for (l_int i=1 ; i<atoms.size() ; i++)
                atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                                       atoms[i]->pos-atoms[0]->pos );
            //   inertia.calc( atoms[0]->pos , atoms );
        }
        calcProps();
    }

    void print(int verbose) {
        if (verbose&InternalDynamics::printNodeTheta && !useEuler) 
            cout << setprecision(8) << atoms[0] << ": theta: " 
                 << q    << ' ' 
                 << dq   << ' ' 
                 << ddq  << '\n';
        else 
            HingeNodeSpec<3>::print(verbose);
    }
};

/**
 * U-joint like joint type which allows rotation about the two axes
 * perpendicular to zDir. This is appropriate for diatoms and for allowing 
 * torsion+bond angle bending.
 */
class HNodeRotate2 : public HingeNodeSpec<2> {
public:
    virtual const char* type() { return "rotate2"; }

    HNodeRotate2(const HingeNode* node,
                 const Vec3&      zVec,
                 int&             cnt)
      : HingeNodeSpec<2>(node,cnt,makeJointFrameFromZAxis(zVec))
    { 
    }

    void calcRot() { 
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

        const Mat3 R_PB = orthoTransform( Mat3(a) , R_BJ );
        R_GB = getR_GP() * R_PB; //the spatial rotation matrix
    }

    void calcH() {
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;

        Vec3 x = scale * R_GB * (R_BJ * Vec3(1,0,0));
        Vec3 y = scale * R_GB * (R_BJ * Vec3(0,1,0));

        Mat23 zMat23(0.0);
        H = blockMat12(catRow(x,y) , zMat23);
    }

    void getInternalForce(RVec& v) {
        FixedVector<double,2> torque = forceInternal;
        //   torque *= DEG2RAD; ??
        RSubVec(v,offset_,2)   = torque.vector();
    }   

    void toCartesian() { 
        atoms[0]->pos = parent->getAtom(0)->pos + getR_GP() * refOrigin_P;
        calcRot();
        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->pos = atoms[0]->pos + R_GB * atomStations_B(i);

        // update vel using recursive formula
        phi = PhiMatrix( atoms[0]->pos - parent->getAtom(0)->pos );
        calcH();
        sVel = transpose(phi) * parent->getSpatialVel()
                + MatrixTools::transpose(H) * dTheta;
        atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

        for (l_int i=1 ; i<atoms.size() ; i++)
            atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                                   atoms[i]->pos-atoms[0]->pos );
        calcProps();
    }
};

/**
 * The "diatom" joint is the equivalent of a free joint for a body with no inertia in
 * one direction, such as one composed of just two atoms. It allows unrestricted
 * translation but rotation only about directions perpendicular to the body's
 * inertialess axis.
 */
class HNodeTranslateRotate2 : public HingeNodeSpec<5> {
public:
    virtual const char* type() { return "diatom"; }

    HNodeTranslateRotate2(const HingeNode*  node,
                          const Vec3&       zVec,
                          int&              cnt)
      : HingeNodeSpec<5>(node,cnt,makeJointFrameFromZAxis(zVec))
    { 
    }

    void calcRot() { 
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
        const Mat3 R_PB = orthoTransform( Mat3(R_JiJ) , R_BJ ); // orientation of B in parent P
        R_GB = getR_GP() * R_PB; //the spatial rotation matrix
    }

    void calcH() {
        //   calcRot(); //FIX: does this need to be calculated here?
        //double scale=InternalDynamics::minimization?DEG2RAD:1.0;
        double scale=1.0;

        Vec3 x = scale * R_GB * (R_BJ * Vec3(1,0,0));
        Vec3 y = scale * R_GB * (R_BJ * Vec3(0,1,0));

        Mat23 zMat23(0.0);
        Mat3  zMat33(0.0);
        using MatrixTools::transpose;
        H = blockMat22(catRow(x,y) , zMat23 ,
                       zMat33      , transpose(getR_GP()));
    }

    void getInternalForce(RVec& v) {
        FixedVector<double,2> torque = RSubVec5(forceInternal, 0, 2).vector();
        //torque *= DEG2RAD; ??
        RSubVec(v,offset_,2)   = torque.vector();
        RSubVec(v,offset_+2,3) = RSubVec5(forceInternal,2,3).vector();
    }   

    void toCartesian() { 
        atoms[0]->pos = parent->getAtom(0)->pos
                        + getR_GP()
                          * (refOrigin_P + Vec3(RSubVec5(theta,2,3).vector()));
        calcRot();
        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->pos = atoms[0]->pos + R_GB * atomStations_B(i);

        // update vel using recursive formula
        phi = PhiMatrix( atoms[0]->pos - parent->getAtom(0)->pos );
        calcH();
        sVel = transpose(phi) * parent->getSpatialVel()
                + MatrixTools::transpose(H) * dTheta;
        atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

        for (l_int i=1 ; i<atoms.size() ; i++)
            atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                                   atoms[i]->pos-atoms[0]->pos );
        //   inertia.calc( atoms[0]->pos , atoms );
        calcProps();
    }
};

/**
 * This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
 * about a particular axis.
 */
class HNodeTorsion : public HingeNodeSpec<1> {
public:
    virtual const char* type() { return "torsion"; }

    HNodeTorsion(const HingeNode*   node,
                 const Vec3&        rotDir,
                 int&               cnt)
      : HingeNodeSpec<1>(node,cnt,makeJointFrameFromZAxis(rotDir))
    { 
    }

    void calcRot() {
        //   double scale=InternalDynamics::minimization?DEG2RAD:1.0; ??
        double scale=1.0;
        double sinTau = sin( scale * theta(0) );
        double cosTau = cos( scale * theta(0) );
        double a[] = { cosTau , -sinTau , 0.0 ,
                       sinTau ,  cosTau , 0.0 ,
                       0.0    ,  0.0    , 1.0 };
        const Mat3 A(a); //rotation about z-axis
        const Mat3 R_PB = orthoTransform( A , R_BJ ); // orientation of B in parent's frame P
        R_GB = getR_GP() * R_PB; // spatial orientation of B in ground frame G
    };

    // Calc H matrix in space-fixed coords.
    void calcH() {
        // This only works because the joint z axis is the same in B & P
        // because that's what we rotate around.
        Vec3 z = getR_GP() * (R_BJ * Vec3(0,0,1)); // R_BJ=R_PJi
        FixedMatrix<double,1,3> zMat(0.0);
        H = blockMat12( FixedMatrix<double,1,3>(z.getData()) , zMat);
    }

//  void calcVel() { //calc vel, dTheta given the atomic positions, velocities
//   HingeNodeSpec<1>::calcVel();
//  }

    void getInternalForce(RVec& v) {
        double torque = forceInternal(0);
        //torque *= DEG2RAD; ??
        v(offset_) = torque;
    }   

#ifdef READABLE
    void toCartesian() { 
        atoms[0]->pos = parent->getAtom(0)->pos + getR_GP() * refOrigin_P;
        calcRot();
        for (int i=1 ; i<atoms.size() ; i++)
            atoms[i]->pos = atoms[0]->pos + R_GB * atomStations_B(i);

        phi = PhiMatrix( atoms[0]->pos - parent->getAtom(0)->pos );
        calcH();
        sVel = transpose(phi) * parent->getSpatialVel()
                + MatrixTools::transpose(H) * dTheta;
        atoms[0]->vel = Vec3::subVec(sVel,3,5);

        for (l_int i=1 ; i<atoms.size() ; i++)
            atoms[i]->vel = atoms[0]->vel + cross( Vec3::subVec(sVel,0,2) , 
                                                   atoms[i]->pos-atoms[0]->pos );
        //   inertia.calc( atoms[0]->pos , atoms );
        calcProps();
    }
#else /* UNREADABLE */
    void toCartesian() { 
        atoms[0]->pos = getR_GP() * refOrigin_P;
        atoms[0]->pos += parent->getAtom(0)->pos;

        calcRot();
        for (int i=1 ; i<atoms.size() ; i++) {
            atoms[i]->pos = R_GB * atomStations_B(i);
            atoms[i]->pos += atoms[0]->pos;
        }

        phi = PhiMatrix( atoms[0]->pos - parent->getAtom(0)->pos );
        calcH();
        sVel = transpose(phi) * parent->getSpatialVel()
                + MatrixTools::transpose(H) * dTheta;
        atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

        for (l_int i=1 ; i<atoms.size() ; i++) {
            atoms[i]->vel = cross( RSubVec6(sVel,0,3).vector() , 
            atoms[i]->pos-atoms[0]->pos );
            atoms[i]->vel += atoms[0]->vel;
        }
        //   inertia.calc( atoms[0]->pos , atoms );
        calcProps();
    }
#endif /* READABLE */
};

void
InertiaTensor::calc(const Vec3&     center,
                    const AtomList& aList) 
{
    set(0.0);
    Mat3 &m = *this;
    for (int cnt=0 ; cnt<aList.size() ; cnt++) {
        const IVMAtom* a = &(*aList[cnt]);
        m(0,0) += a->mass * (sq(a->pos(1)-center(1)) + sq(a->pos(2)-center(2)));
        m(1,1) += a->mass * (sq(a->pos(0)-center(0)) + sq(a->pos(2)-center(2)));
        m(2,2) += a->mass * (sq(a->pos(0)-center(0)) + sq(a->pos(1)-center(1)));
        m(0,1) -= a->mass * (a->pos(0)-center(0)) * (a->pos(1)-center(1));
        m(1,2) -= a->mass * (a->pos(1)-center(1)) * (a->pos(2)-center(2));
        m(0,2) -= a->mass * (a->pos(0)-center(0)) * (a->pos(2)-center(2));
    }
    m(1,0) = m(0,1);
    m(2,1) = m(1,2);
    m(2,0) = m(0,2);
}

//void 
//combineNodes(const HingeNode* node1,
//	       const HingeNode* node2)
//  //
//  // group the atoms of node2 with those of node 1
//  // 
//  //
//{
// cout << "combineNodes: merging node containing atoms: { ";
// cout << node1->atoms[0]->index;
// for (int i=0 ; i<node1->atoms.size() ; i++)
//   cout << " , " << node1->atoms[i]->index ;
// cout << "}\n";
// cout << "\twith that containing atoms:";
// cout << node2->atoms[0]->index;
// for (l_int i=0 ; i<node2->atoms.size() ; i++)
//   cout << " , " << node2->atoms[i]->index ;
// cout << "}\n";
//
// int groupExists=0;
// Atom* remAtom = node1->atoms[0];
// using InternalDynamics::groupList;
// for (l_int i=0 ; i<groupList.size() ; i++)
//   if ( groupList[i].contains(remAtom->index) ) {
//     groupExists=1;
//     for (int j=0 ; j<node2->atoms.size() ; j++)
//	 groupList[i].append( node2->atoms[j]->index );
//     break;
//   }
//
// if ( !groupExists ) {
//   CDSList< int > l;
//   l.append( remAtom->index );
//   for (int j=0 ; j<node2->atoms.size() ; j++)
//     l.append( node2->atoms[j]->index );
//   groupList.append( l );
// }
//} /* combineNodes */

//
// Replace the passed-in proto-HingeNode with one that has a joint
// as close as possible to the requested one. On entry, cnt indicates
// the next available state variable slot; we grab some and update cnt
// on exit.
//
HingeNode*
construct(HingeNode*                         node,
          const InternalDynamics::HingeSpec& hingeSpec,
          int&                               cnt)
{
    using InternalDynamics::Exception;

    HingeNode* newNode=0;
    const IVM* ivm = node->ivm;
    String type = hingeSpec.type;
    type.downcase();

    if ( type == "ground" ) 
        newNode = new GroundBody(node);

    if (   (type == "bendtorsion" && node->level>2)
        || (type == "bend" && !node->isGroundNode())
              && (node->atoms.size()>1 || node->children.size()>0)
              &&  node->atoms[0]->bonds.size()>0 ) 
    {
        IVMAtom* atom0 = ivm->getAtoms()[ hingeSpec.atom0 ];
        IVMAtom* atom1 = ivm->getAtoms()[ hingeSpec.atom1 ];

        //direction perp. to two bonds
        Vec3 dir = cross(atom0->pos - node->atoms[0]->pos,
                         atom1->pos - node->atoms[0]->pos);
        if ( norm(dir) > 1e-12 ) {
            newNode = new HNodeTorsion(node, dir, cnt);
        } else {
            // the two bonds are colinear - this is an error unless
            //   node has exactly two bonds
            //   and one of atom0 or atom1 has only one bond.
            if ( node->atoms[0]->bonds.size()==2
                 && (atom0->bonds.size() == 1 || atom1->bonds.size() == 1))
            {
                //in this case, the plane of the bend angle needs only to be 
                // perpendicular to the bond
                //rotation matrix which takes the z-axis
                // to (atoms[0]->pos - parentAtom->pos)
                // notes of 12/6/99 - CDS
                const Vec3 u = unitVec( atom1->pos - node->atoms[0]->pos );
                double theta = acos( u.z() );
                double psi   = atan2( u.x() , u.y() );
                double a[] = { cos(psi) , cos(theta)*sin(psi) , sin(psi)*sin(theta),
                              -sin(psi) , cos(theta)*cos(psi) , cos(psi)*sin(theta),
                               0        , -sin(theta)         , cos(theta)         };
                Vec3 x(1,0,0);
                dir = Mat3(a) * x;
                cerr << "norm: " << dot(u,dir) << endl; // should be zero
                newNode = new HNodeTorsion(node, dir, cnt);
            }
        }
    } 

    if ( type == "torsion" ) {
        if (node->isBaseNode() &&             // allow full dof for base node
            node->parentAtom->index==0)       // unless it is bonded to fixed atom
            //FIX: this seems ill-thought out
            type = "full"; 
        else if ( node->atoms.size()>1 || node->children.size()>0 )
            newNode = new HNodeTorsion(node, 
                                       node->atoms[0]->pos - node->parentAtom->pos,
                                       cnt);
    } 

    if (type == "rotate" ) {
        if ( node->atoms.size()>2 )
            newNode = new HNodeRotate3(node,cnt,
                                       ivm->minimization());
        else if ( node->atoms.size()==2 )
            newNode = new HNodeRotate2(node,
                            node->atoms[1]->pos - node->atoms[0]->pos,cnt);
    }

    if (type == "translate")
        newNode = new HNodeTranslate(node,cnt);

    if ( type == "full" || type == "unknown") { // all other cases- full freedom of movement
        if ( node->atoms.size()>2 || node->atoms.size()==2 && node->children.size() )
            newNode = new HNodeTranslateRotate3(node,cnt,
                                                ivm->minimization());
        else if ( node->atoms.size()==2 )
            newNode = new HNodeTranslateRotate2(node,
                            node->atoms[1]->pos - node->atoms[0]->pos, cnt);
        else if ( node->atoms.size()==1 )
            newNode = new HNodeTranslate(node,cnt);
    }

    if ( !newNode ) {
        cerr << "Bad Hinge type or topology not supported.\n";
        cerr << "node at level " << node->level 
             << " containing atoms: " << node->atoms << '\n';
        throw Exception("Bad Hinge type or topology not supported");
    }

    //
    // connect to other related hinge
    //
    if (!node->isGroundNode()) {
        int index = node->parent->children.getIndex(node);
        node->parent->children[index] =  newNode;
    }
    for (int i=0 ; i<node->children.size() ; i++)
        node->children[i]->parent = newNode;

    delete node; // now we're done with this.
    return newNode;
}

HingeNode::HingeNode(const IVM*     ivm,
                     IVMAtom*       hingeAtom,
                     const IVMAtom* parAtom,
                     HingeNode*     parNode) 
  : parent(parNode),
    children(0,0),
    level(hingeAtom->index==0 ? 0 : parNode->getLevel()+1),
    atoms(0,1),
    sVel(0.),
    sAcc(0.),
    phi(PhiMatrix(Vec3(0.))),
    psiT(0.),
    P(0.),
    z(0.),
    tau(0.),
    Gepsilon(0.),
    R_GB(0.),
    Y(0.), 
    ivm(ivm),
    parentAtom(parAtom)
{
    R_GB.setDiag(1.0);  // Body frame B is initially aligned with the ground frame G.
    atoms.append(hingeAtom);
    hingeAtom->node = this;
}

void
HingeNode::addChild(HingeNode* node) {
    children.append( node );
}

ostream& 
operator<<(ostream& s, const HingeNode& node) {
    for (int i=0 ; i<node.atoms.size() ; i++)
        s << node.atoms[i] << ", ";
    return s;
}

//template<int dof>
//void 
//HingeNodeSpec<dof>::velFromCartesian()
//  //
//  // calculate relative spatial velocities
//  //
//{
// // firstStep = 0;
// // const Node* pNode = parentHinge.remNode;
// // these must can be calc'd now - FIX placement
// inertia.calc( atoms[0]->pos , atoms ); 
// calcPhi( atoms[0]->pos - parent->atoms[0]->pos );
//
// // if ( level>1 ) { // base node is taken care of in AtomTree::calcVel
// Vec3 omega = Vec3::subVec( parent->vel, 0, 2 );
// if ( atoms.size()>2 ) {
//   Vec3 L(0.0); // angular momentum
//   for (int i=1 ; i<atoms.size() ; i++)
//     L += atoms[i]->mass *
//	      cross( atoms[i]->pos - atoms[0]->pos , 
//		     atoms[i]->vel - atoms[0]->vel);
//   omega = inverse(inertia) * L; // hopefully inertia is of full rank
//   //FIX: should this be included?
////   if ( level==2 && parent->atoms.size()==1 ) {
////     // correct the component of angular velocity about the bond to the
////     // parent node- with a base node consisting of one atom, omega
////     // is not entirely independent of remNode->omega
////     Vec3 u = atoms[0]->pos - parentAtom->pos;
////     u /= norm(u);
////     omega += (dot(u , Vec3::subVec( parent->vel, 0, 2 ) )-
////		 dot(u , omega                             )) * u;
////   }
// } else if ( atoms.size()==2 ) {
//   Vec3 r = atoms[1]->pos - atoms[0]->pos;
//   omega = 1.0/abs2(r) *
//	       cross( r , atoms[1]->vel - atoms[0]->vel );
////   if ( level==1 && children.size() ) {
////     // but it is specified in this case
////     Vec3 q1 = atoms[0]->pos;    	  Vec3 v1 = atoms[0]->vel;
////     Vec3 q2 = atoms[1]->pos;    	  Vec3 v2 = atoms[1]->vel;
////     Vec3 q3 = children[0]->atoms[0]->pos;Vec3 v3 = children[0]->atoms[0]->vel;
////     double a = dot( cross(q2-q1,q3-q2) ,
////		       (v3-v2) - dot(q3-q2,q2-q1)*(v2-v1)/sq(norm(q2-q1))) /
////		  sq(norm( cross(q2-q1,q3-q2)));
////     omega += a * (q2-q1);
////   } else {
//     // component along bond is not specified- it should be same as that
//     // for parent node.
//     Vec3 u = unitVec( atoms[1]->pos - atoms[0]->pos );
//     omega += dot(u , Vec3::subVec( parent->vel, 0, 2 ) ) * u;
//     //   }
// }
//
// // this velocity (velV) may have components inconsistent with the allowed 
// // internal degrees of freedom- they are projected out in the 
// // Hinge::calcVel call
//
// // parentHinge.calcVel( vel , pNode->vel );
// // calc velocity of internal coordinates
// calcH();  // build H
// 
// //
// // recalc. the cartesian velocity so it is consistent with that of a 
// // rigid body
// //
// posCM_.set( 0.0 );
// mass_ = 0; //should be calc'd beforehand - once.
// for (int i=0 ; i<atoms.size() ; i++) {
//   Atom* a = atoms[i];
//   posCM_ += a->mass * a->pos;
//   mass_ += a->mass;
// }
// posCM_ /= mass_;
//
// Vec3 velCM(0.0);
// for (l_int i=0 ; i<atoms.size() ; i++)
//   velCM += atoms[i]->mass * atoms[i]->vel;
// velCM /= mass_;
// for (l_int i=0 ; i<atoms.size() ; i++) 
//   atoms[i]->vel = velCM +
//		     cross( omega , atoms[i]->pos-posCM_ );
// vel = blockVec( omega , atoms[0]->vel );
// // the next line assumes that the columns of H are normalized: if not
// // this function must be specilized in the Node specialization
// // FixedVector<double,dof> oldDTheta = dTheta;
// dTheta = H * (vel - transpose(phi) * parent->vel);
//// if ( norm( dTheta-oldDTheta) > 1e-8 ) 
////   cout << "calcVel: dTheta is inconsistent:\n"
////	  << "\toldDTheta: " << oldDTheta << '\n'
////	  << "\t   dTheta: " << dTheta    << '\n';
// calcProps();
// 
//} /* HingeNodeSpec::velFromCartesian */

//
// Calculate acceleration in internal coordinates.
//
template<int dof> void 
HingeNodeSpec<dof>::calcAccel() {
    // const Node* pNode = parentHinge.remNode;
    //make sure that this is phi is correct - FIX ME!
    // base alpha = 0!!!!!!!!!!!!!
    Vec6 alphap = transpose(phi) * parent->sAcc;
    ddTheta = nu - MatrixTools::transpose(G) * alphap;

    sAcc   = alphap + MatrixTools::transpose(H) * ddTheta + a;  
}

//
// Calc posCM, mass, b, a, Mk
//      phi, svel, inertia
//
// Depends on atoms->(pos,vel), dTheta, parent->svel.
// Should be calc'd from base to tip.
template<int dof> void 
HingeNodeSpec<dof>::calcProps() {
    forceInternal.set( 0.0 );

    inertia.calc(atoms[0]->pos , atoms);
    phi = PhiMatrix(atoms[0]->pos - parent->atoms[0]->pos);

    posCM_.set( 0.0 );
    mass_ = 0; //should be calc'd beforehand - once.
    for (int i=0 ; i<atoms.size() ; i++) {
        IVMAtom* a = atoms[i];
        posCM_ += a->mass * a->pos;
        mass_ += a->mass;
    }
    posCM_ *= 1.0/mass_;

    // calc Mk: the mass matrix
    Mat3 uVec(0.0); uVec.setDiag(1.0);
    Mat3 offDiag = mass_*crossMat(posCM_ - atoms[0]->pos);
    Mk = blockMat22( inertia , offDiag ,
                    -offDiag , mass_*uVec );

    if ( !ivm->minimization() ) {
        sVel = transpose(phi) * parent->sVel
               + MatrixTools::transpose(H) * dTheta;
        Vec3 omega = RSubVec6(sVel, 0,3).vector();
        if ( atoms.size() > 1 ) {
            // gyroscopic spatial force = 0 if only one atom in node
            Vec3 gMoment = cross(omega, inertia * omega);
            Vec3 gForce  = mass_ * cross(omega, 
                                         cross(omega, (posCM_-atoms[0]->pos)));
            b = blockVec(gMoment, gForce);
        }

        const Vec3& rOmega = RSubVec6( parent->sVel, 0,3).vector();
        const Vec3& vel3   = RSubVec6( sVel, 3,3).vector();
        const Vec3& rVel3  = RSubVec6( parent->sVel, 3,3).vector();

        // calc a: coriolis acceleration
        a = blockMat22(crossMat(rOmega),      Mat3(0.0),
                       Mat3(0.0)       ,crossMat(rOmega)) * 
                                          MatrixTools::transpose(H) * dTheta;
        a += blockVec( Vec3(0.0), cross(rOmega, vel3 - rVel3) );
    }
}

template<int dof> void
HingeNodeSpec<dof>::calcD_G(const Mat6& P) {
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
             << "atoms:\n";
        for (int i=0 ; i<atoms.size() ; i++)
            cerr << '\t' << ivm->idAtom(atoms[i]->index)
                 << " : " << atoms[i]->pos << '\n';
        cerr << endl;
        throw Exception("calcD_G: singular D matrix. Bad topology?");
    }
    G = P * MatrixTools::transpose(H) * DI;
}

//template<int dof>
//void
//HingeNodeSpec<dof>::calcEpsilon_nu(const Vec6& z)
//{
//    epsilon = forceInternal - H*z;
//    nu = DI * epsilon;
//}

template<int dof> double
HingeNodeSpec<dof>::kineticE() {
    double ret = dot(sVel , Mk*sVel);

    // Vec3 omega = Vec3::subVec( sVel,0,2 );
    // Vec3 v     = Vec3::subVec( sVel,3,5 ); 
    //// InertiaTensor inertia;
    //// inertia.calc( posCM_ , atoms );
    // double trans = mass() * dot(v,v);
    // double rot   = dot(omega , inertia * omega);
    // double ret = trans + rot;
    // cout << ' ' << trans << ' ' << rot << ' ';

    return 0.5*ret;
}

template<int dof> void
HingeNodeSpec<dof>::setPosVel(const RVec& posv, const RVec& velv) {
    theta  = ConstRSubVec(posv,offset_,dof).vector();
    dTheta = ConstRSubVec(velv,offset_,dof).vector();
    toCartesian();
}

template<int dof> void
HingeNodeSpec<dof>::setVel(const RVec& velv) {
    dTheta = ConstRSubVec(velv,offset_,dof).vector();
    sVel   = transpose(phi) * parent->sVel
             + MatrixTools::transpose(H) * dTheta;
    atoms[0]->vel = Vec3( RSubVec6(sVel,3,3).vector() );

    for (int i=1 ; i<atoms.size() ; i++)
        atoms[i]->vel = atoms[0]->vel + cross( RSubVec6(sVel,0,3).vector() , 
                                               atoms[i]->pos-atoms[0]->pos );
    //   inertia.calc( atoms[0]->pos , atoms );
}

//
// to be called from base to tip.
//
template<int dof> void
HingeNodeSpec<dof>::setVelFromSVel(const Vec6& sVel) {
    dTheta = H * (sVel - transpose(phi)*parent->sVel);
}

template<int dof> void
HingeNodeSpec<dof>::getPos(RVec& v) {
    RSubVec(v,offset_,dof) = theta.vector();
}

template<int dof> void
HingeNodeSpec<dof>::getVel(RVec& v) {
    RSubVec(v,offset_,dof) = dTheta.vector();
}

template<int dof> void
HingeNodeSpec<dof>::print(int verbose) {
    if (verbose&InternalDynamics::printNodeForce) 
        cout << setprecision(8)
             << atoms[0] << ": force: " << forceCartesian << '\n';
    if (verbose&InternalDynamics::printNodePos) 
        cout << setprecision(8)
             << atoms[0] << ": pos: " << atoms[0]->pos << ' ' << atoms[0]->vel
             << ' ' << sAcc << '\n';
    if (verbose&InternalDynamics::printNodeTheta) 
        cout << setprecision(8)
             << atoms[0] << ": theta: " 
             << theta << ' ' << dTheta  << ' ' << ddTheta  << '\n';
}

// Called after calcAccel.
template<int dof> void
HingeNodeSpec<dof>::getAccel(RVec& v) {
    RSubVec(v,offset_,dof) = ddTheta.vector();
}

// CalcInternalForce has been called previously.
template<int dof> void
HingeNodeSpec<dof>::getInternalForce(RVec& v) {
    RSubVec(v,offset_,dof) = forceInternal.vector();
}

//
// Calculate Pk and related quantities. The requires that the children
// of the node have already had their quantities calculated.
//
template<int dof> void
HingeNodeSpec<dof>::calcP() {
    //
    //how much do we need to keep around?
    // it looks like nu and G are the only ones needed for the acceleration
    // calc. The others can be freed after the parent is done with them.
    //
    P = Mk;

    SubMatrix<Mat6> p11(P,0,0,3,3);
    SubMatrix<Mat6> p12(P,0,3,3,3);
    SubMatrix<Mat6> p21(P,3,0,3,3);
    SubMatrix<Mat6> p22(P,3,3,3,3);
    for (int i=0 ; i<children.size() ; i++) {
        // this version is readable
        // P += orthoTransform( children[i]->tau * children[i]->P ,
        //                      transpose(children[i]->phiT) );
        // this version is not
        Mat3 lt = crossMat(children[i]->getAtom(0)->pos - getAtom(0)->pos);
        Mat6 M  = children[i]->tau * children[i]->P;
        SubMatrix<Mat6> m11(M,0,0,3,3);
        SubMatrix<Mat6> m12(M,0,3,3,3);
        SubMatrix<Mat6> m21(M,3,0,3,3);
        SubMatrix<Mat6> m22(M,3,3,3,3);
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
HingeNodeSpec<dof>::calcZ() {
    psiT = MatrixTools::transpose(tau) * transpose(phi);
    calcCartesianForce();
    z = P * a + b + forceCartesian;

    for (int i=0 ; i<children.size() ; i++) 
        z += children[i]->phi
             * (children[i]->z + children[i]->Gepsilon);

    epsilon  = forceInternal - H*z;
    nu       = DI * epsilon;
    Gepsilon = G * epsilon;
}

template<int dof> void
HingeNodeSpec<dof>::calcY() {
    Y = orthoTransform(DI,MatrixTools::transpose(H))
        + orthoTransform(parent->Y,psiT);
}

template<int DOF> void
HingeNodeSpec<DOF>::calcCartesianForce() {
    Vec3 moment(0.0);
    Vec3 force(0.0);
    // notice that the sign is screwey
    for (int i=0 ; i<atoms.size() ; i++) {
        IVMAtom* a = atoms[i];
        Vec3 aForce = a->deriv;
        if ( ivm->frictionCoeff()!=0.0 )
            aForce += a->fric *(1.0-ivm->bathTemp()/ivm->currentTemp()) * a->vel;
        moment += cross( a->pos - atoms[0]->pos , aForce );
        force  += aForce;
    }
    forceCartesian = blockVec(moment, force);
}

//
// Calculate sum of internal force and effective forces due to Cartesian
// forces.
// To be called from tip to base.
// Should be called only once after calcProps.
//
template<int dof> void
HingeNodeSpec<dof>::calcInternalForce() {
    calcCartesianForce();
    z = forceCartesian;

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
HingeNodeSpec<dof>::prepareVelInternal() {
    dTheta.set(0.0);
    a.set(0.0);
    b.set(0.0);
    forceCartesian.set(0.0);
    // set Mk = R * M * R^T
    Mk.set(0.0);
    Mat3 unit(0.0); unit.setDiag(1.0);
    Vec3 q0 = atoms[0]->pos;
    for (int i=0 ; i<atoms.size() ; i++)
        Mk += atoms[i]->mass
              * blockMat22(crossMat(atoms[i]->pos-q0)*
                             crossMat(q0-atoms[i]->pos) , crossMat(atoms[i]->pos-q0),
                           crossMat(q0-atoms[i]->pos)   , unit                      );

    for (l_int i=0 ; i<atoms.size() ; i++)
        atoms[i]->deriv.set(0.0);
}

//
// Set up various variables so that calling calcP followed by
// calling getAcc returns the time-derivative of the internal coordinates.
// Note that sVel is passed in- not the member variable.
// To be called from tip to base.
//
template<int dof> void
HingeNodeSpec<dof>::propagateSVel(const Vec6& sVel) {
    //sAcc used as a temporary
    this->sVel = sVel;
    sAcc = sVel;
    for (int i=0 ; i<children.size() ; i++)
        sAcc += children[i]->phi * children[i]->sAcc;
    forceInternal = H * sAcc;
}

static double
sumMassToTip(const HingeNode* n)
{
    double ret = 0.0;
    for (int i=0 ; n->getChild(i) ; i++)
        ret += sumMassToTip( n->getChild(i) );
    for (l_int i=0 ; n->getAtom(i) ; i++)
        ret += n->getAtom(i)->mass;
    return ret;
}

static double
sumInertiaToTip(const HingeNode* n,
                const Vec3&      pos,
                const Vec3&      dir)
{
    double ret = 0.0;
    for (int i=0 ; n->getChild(i) ; i++)
        ret += sumInertiaToTip( n->getChild(i) , pos , dir);
    for (l_int i=0 ; n->getAtom(i) ; i++)
        ret += n->getAtom(i)->mass * (abs2(n->getAtom(i)->pos - pos) -
                                     sq(dot(dir,n->getAtom(i)->pos - pos)));
    return ret;
}

template<int dof> double
HingeNodeSpec<dof>::approxKE() {
    double ret=0.0;
    for (int i=0 ; i<dof ; i++) {
        double mass = 0.0;
        for (int j=0 ; j<3 ; j++)
            if ( H(i,j) != 0.0 ) {
                mass = sumInertiaToTip( this, 
                                        getAtom(0)->pos, 
                                        Vec3::subCol(MatrixTools::transpose(H),i,0,2) );
                break;
            }
        if ( mass == 0.0 ) 
            mass = sumMassToTip( this );
        ret += mass * sq(dTheta(i));
    }
    return ret;
}


