#ifndef RIGID_BODY_NODE_H_
#define RIGID_BODY_NODE_H_

#include "simbody/internal/SimbodyCommon.h"

#include "phiMatrix.h"
#include "internalDynamics.h"

#include "fixedMatrix.h"
#include "cdsMatrix.h"
#include "cdsVector.h"
#include "cdsList.h"

#include "RBMassProperties.h"

#include <cassert>

class CDSVec3;

typedef FixedVector<double,6>   CDSVec6;
typedef FixedMatrix<double,6,6> CDSMat66;
typedef CDSMatrix<double>       RMat;
typedef CDSVector<double,1>     RVec;


/**
 * This is an abstract class representing a body and its (generic) inboard joint, that is,
 * the joint connecting it to its parent. Concrete classes are derived from this one to
 * represent each specific type of joint.
 *
 * RigidBodyNodes are linked into a tree structure, organized into levels as described 
 * in Schwieters' JMR paper. The root is a special 'Ground' node defined to be at 
 * level 0 and containing only atoms fixed to ground. The level 1 nodes (referred to
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
 * other direction use R_YX=R_XY' (apostrophe is Matlab for "transpose"). This convention
 * provides flawless composition of rotations as long as you "match up" the frames
 * (i.e., adjacent frame symbols must be the same). So if you have R_XY and R_ZX 
 * and you want R_YZ you can easily get it like this:
 *     R_YZ = R_YX*R_XZ = (R_XY')*(R_ZX') = (R_ZX*R_XY)'.
 * Also note that these are orthogonal matrices so R_XY*R_YX=I.
 *
 * Every body has a body frame B. In the reference configuration that we use
 * to define the bodies, all frames B are aligned with the ground frame G. You can
 * think of this as defining the B frames by painting images of the G frame on
 * each body (translated to B's origin OB) when we first see the bodies in the reference
 * frame. Later the bodies move and take their B frames with them. For convenience, we
 * refer to the body frame of a body's unique parent as the 'P' frame. Initially frames P
 * and B are aligned (and both are aligned with Ground). Later the P and B frames
 * will differ, but only by the relative orientation and translation induced by
 * the joint connecting B to P. That is, rotation matrix R_PB expresses the
 * relative orientation between the parent and child caused by the current
 * joint coordinates. When all of that joint's coordinates are 0, B and P will
 * once again be aligned (R_PB=I) although of course in general they will not be 
 * aligned with Ground unless *all* the joint coordinates between B and G are zero.
 *
 * In addition to the B frame fixed on every body, the inboard joint has its own 
 * frame J, which is fixed with respect to B. In some cases J and B will be the
 * same, but not always. The constant rotation matrix R_BJ provides the orientation
 * of the inboard joint's frame in its body's frame. If B is the i'th child of
 * its parent P, then there is also a parent-fixed frame Ji which is the image of J
 * when the joint coordinates are zero. That is, R_JiJ is the orientation change
 * induced by the joint coordinates. Note that because of how we define the 
 * reference configuration, R_PJi = R_BJ so we don't need to store both matrices
 * explicitly. With these definitions we can easily calculate R_PB as
 *     R_PB = R_PJi*R_JiJ*R_JB = R_BJ*R_JiJ*(R_BJ)'.
 */
class RigidBodyNode {
public:
    class VirtualBaseMethod {};    // an exception

    virtual ~RigidBodyNode() {}

    RigidBodyNode& operator=(const RigidBodyNode&);

    /// Factory for producing concrete RigidBodyNodes based on joint type.
    static RigidBodyNode* create(
        const RBMassProperties& m,            // mass properties in body frame
        const RBFrame&          jointFrame,   // inboard joint frame J in body frame
        RBJointType             type,
        bool                    isReversed,   // child-to-parent orientation?
        bool                    useEuler,     // TODO: kludge (true if minimizing)
        int&                    nextStateOffset); 

    /// Register the passed-in node as a child of this one, and note in
    /// the child that this is its parent. Also set the reference frame in the child.
    void addChild(RigidBodyNode* child, const RBFrame& referenceFrame);

    RigidBodyNode*   getParent() const {return parent;}
    void             setParent(RigidBodyNode* p) { parent=p; }

    int              getNChildren()  const {return children.size();}
    RigidBodyNode*   getChild(int i) const {return (i<children.size()?children[i]:0);}

    /// Return this node's level, that is, how many ancestors separate it from
    /// the Ground node at level 0. Level 1 nodes (directly connected to the
    /// Ground node) are called 'base' nodes.
    int              getLevel() const  {return level;}
    void             setLevel(int i)   {level=i;}

    int              getNodeNum() const {return nodeNum;}
    void             setNodeNum(int n)  {nodeNum=n;}

    bool             isGroundNode() const { return level==0; }
    bool             isBaseNode()   const { return level==1; }

    int              getStateOffset() const {return stateOffset;}

    const RBMassProperties& getMassProperties() const {return massProps_B;}
    const double&    getMass()         const {return massProps_B.getMass();}
    const CDSVec3&   getCOM_B()        const {return massProps_B.getCOM();}
    const RBInertia& getInertia_OB_B() const {return massProps_B.getInertia();}
    const RBInertia& getInertia_OB_G() const {return inertia_OB_G;}

    const CDSVec3& getCOM_G()          const {return COM_G;}
    const RBInertia& getInertia_CB_B() const {return inertia_CB_B;}


    /// Return R_GB, the rotation (direction cosine) matrix giving the 
    /// spatial orientation of this body's frame B (that is, B's orientation
    /// in the ground frame G).
    const CDSMat33&  getR_GB() const {return R_GB;}

    /// Return OB_G, the spatial location of the origin of the B frame, that is, 
    /// measured from the ground origin and expressed in ground.
    const CDSVec3&   getOB_G() const {return OB_G; }

    /// Return R_GP, the rotation (direction cosine) matrix giving the
    /// orientation of this body's *parent's* body frame (which we'll call
    /// P here) in the ground frame G.
    const CDSMat33&  getR_GP() const {assert(parent); return parent->getR_GB();}

    /// Return OP_G, the spatial location of the origin of the P frame, that is, 
    /// measured from the ground origin and expressed in ground.
    const CDSVec3&   getOP_G() const {assert(parent); return parent->getOB_G();}

    void setSpatialVel(const CDSVec6& v) { sVel=v; }

    /// Return the inertial angular velocity of body frame B (i.e., angular
    /// velocity with respect to the ground frame), expressed in the ground frame.
    const CDSVec3&   getSpatialAngVel() const
        {return *reinterpret_cast<const CDSVec3*>(&sVel[0]);}

    /// Return the inertial velocity of OB (i.e., velocity with respect
    /// to the ground frame), expressed in the ground frame.
    const CDSVec3&   getSpatialLinVel() const
        {return *reinterpret_cast<const CDSVec3*>(&sVel[3]);}

    /// Return the inertial angular acceleration of body frame B (i.e., angular
    /// acceleration with respect to the ground frame), expressed in the ground frame.
    const CDSVec3&   getSpatialAngAcc() const
        {return *reinterpret_cast<const CDSVec3*>(&sAcc[0]);}

    /// Return the inertial acceleration of OB (i.e., acceleration with respect
    /// to the ground frame), expressed in the ground frame.
    const CDSVec3&   getSpatialLinAcc() const
        {return *reinterpret_cast<const CDSVec3*>(&sAcc[3]);}

    const CDSVec6&   getSpatialVel() const {return sVel;}
    const CDSVec6&   getSpatialAcc() const {return sAcc;}

    const PhiMatrix& getPhi()  const {return phi;}
    const CDSMat66&  getPsiT() const {return psiT;}
    const CDSMat66&  getY()    const {return Y;}

    /// Introduce new values for generalized coordinates and calculate
    /// all the position-dependent kinematic terms.
    virtual void setPos(const RVec&)=0;

    /// Introduce new values for generalized speeds and calculate
    /// all the velocity-dependent kinematic terms. Assumes setPos()
    /// has already been called.
    virtual void setVel(const RVec&)=0;

    double calcKineticEnergy() const;   // from spatial quantities only

    virtual const char* type()     const {return "unknown";}
    virtual int         getDOF()   const {return 0;} //number of independent dofs
    virtual int         getDim()   const {return 0;} //dofs plus quaternion constraints

    virtual void enforceConstraints(RVec& pos, RVec& vel) {throw VirtualBaseMethod();}

    virtual void calcP()                                     {throw VirtualBaseMethod();}
    virtual void calcZ    (const CDSVec6& spatialForce)      {throw VirtualBaseMethod();}
    virtual void calcY()                                     {throw VirtualBaseMethod();}
    virtual void calcAccel()                                 {throw VirtualBaseMethod();}

    virtual void calcInternalForce(const CDSVec6& spatialForce) {throw VirtualBaseMethod();}

    virtual void setVelFromSVel(const CDSVec6&) {throw VirtualBaseMethod();}


    virtual void getPos(RVec&)   const {throw VirtualBaseMethod();}
    virtual void getVel(RVec&)   const {throw VirtualBaseMethod();}
    virtual void getAccel(RVec&) const {throw VirtualBaseMethod();}

    virtual void getInternalForce(RVec&) const {throw VirtualBaseMethod();}
    virtual RMat getH() const {throw VirtualBaseMethod();}

    virtual void print(int) const { throw VirtualBaseMethod(); }

    void nodeDump(ostream&) const;
    virtual void nodeSpecDump(ostream& o) const { o<<"NODE SPEC type="<<type()<<endl; }


    static const double DEG2RAD; //using angles in degrees balances gradient

protected:
    /// This is the constructor for the abstract base type for use by the derived
    /// concrete types in their constructors.
    RigidBodyNode(const RBMassProperties& mProps_B,
                  const CDSVec3& originOfB_P, // and R_BP=I in ref config
                  const CDSMat33& rot_BJ, const CDSVec3& originOfJ_B)
      : stateOffset(-1), parent(0), children(0,0), level(-1), nodeNum(-1),
        massProps_B(mProps_B), inertia_CB_B(mProps_B.calcCentroidalInertia()),
        R_BJ(rot_BJ), OJ_B(originOfJ_B), refOrigin_P(originOfB_P)
    {
        R_PB.set(0.); R_PB.setDiag(1.); OB_P=refOrigin_P;
        V_PB_G.set(0.); sVel.set(0.); sAcc.set(0.);
        R_GB.set(0.); R_GB.setDiag(1.); OB_G.set(0.);
        COM_G.set(0.); COMstation_G = massProps_B.getCOM();
        phi = PhiMatrix(CDSVec3(0.));
        psiT.set(0.); P.set(0.); z.set(0.); tau.set(0.); Gepsilon.set(0.); Y.set(0.);
    }

    typedef CDSList<RigidBodyNode*>   RigidBodyNodeList;

    int               stateOffset;  //index into internal coord pos,vel,acc arrays
    RigidBodyNode*    parent; 
    RigidBodyNodeList children;
    int               level;        //how far from base 
    int               nodeNum;      //unique ID number in RigidBodyTree

    // These are the body properties

    // Fixed forever in the body-local frame B:
    //      ... supplied on construction
    const RBMassProperties massProps_B;
    const CDSMat33 R_BJ;          // orientation of inboard joint frame J, in B
    const CDSVec3  OJ_B;          // origin of J, measured from B origin, expr. in B

    // Reference configuration. This is the body frame origin location, measured
    // in its parent's frame in the reference configuration. This vector is fixed
    // after this node is connected to its parent. The body origin can of course move relative
    // to its parent, but that is not the meaning of this reference configuration vector.
    // (Note however that the body origin is also the location of the inboard joint, 
    // meaning that the origin point moves relative to the parent only due to translations.)
    // Note that by definition the orientation of the body frame is identical to P
    // in the reference configuration so we don't need to store it.
    CDSVec3        refOrigin_P;

    //      ... calculated on construction
    const RBInertia  inertia_CB_B;  // centroidal inertia, expr. in B

    // Calculated relative quantities (these are joint-relative quantities, 
    // but not dof dependent).
    //      ... position level
    CDSMat33 R_PB; // orientation of B in P
    CDSVec3  OB_P; // location of B origin meas & expr in P frame

    //      ... velocity level
    CDSVec6  V_PB_G; // relative velocity of B in P, but expressed in G (omega & v)

    // Calculated spatial quantities
    //      ... position level

    //      Rotation (direction cosine matrix) expressing the body frame B
    //      in the ground frame G. That is, if you have a vector vB expressed
    //      body frame and want it in ground, use vG = R_GB*vB. 
    CDSMat33 R_GB;
    CDSVec3  OB_G;          // origin of B meas & expr in G
    CDSVec3  COM_G;         // B's COM, meas & expr in G

    CDSVec3    COMstation_G;  // measured from B origin, expr. in G
    RBInertia  inertia_OB_G;  // about B's origin, expr. in G

    PhiMatrix phi;           // spatial rigid body transition matrix
    CDSMat66  Mk;            // spatial inertia matrix
    CDSMat66  P;             // articulated body spatial inertia
    CDSMat66  tau;
    CDSMat66  psiT;
    CDSMat66  Y;             // diag of Omega - for loop constraints

    //      ... velocity level
    CDSVec6   a;     // spatial coriolis acceleration
    CDSVec6   b;     // spatial gyroscopic force
    CDSVec6   sVel;  // spatial velocity

    //      ... acceleration level
    CDSVec6   z;
    CDSVec6   Gepsilon;
    CDSVec6   sAcc;              // spatial acceleration


    virtual void velFromCartesian() {}

    friend ostream& operator<<(ostream& s, const RigidBodyNode&);
    template<int dof> friend class RigidBodyNodeSpec;

private:   
    /// Calculate all spatial configuration quantities, assuming availability of
    /// joint-specific relative quantities.
    ///   R_GB
    ///   OB_G
    void calcJointIndependentKinematicsPos();

    /// Calcluate all spatial velocity quantities, assuming availability of
    /// joint-specific relative quantities and all position kinematics.
    ///   sVel  spatial velocity of B
    ///   a     spatial Coriolis acceleration
    ///   b     spatial gyroscopic force
    void calcJointIndependentKinematicsVel();
};

#endif /* RIGID_BODY_NODE_H_ */
