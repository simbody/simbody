#ifndef RIGID_BODY_NODE_H_
#define RIGID_BODY_NODE_H_

#include "phiMatrix.h"
#include "internalDynamics.h"

#include "fixedMatrix.h"
#include "cdsMatrix.h"
#include "cdsVector.h"
#include "cdsList.h"

#include "MassProperties.h"

#include <cassert>

class Vec3;

template<class CHAR> class CDSString;
typedef CDSString<char> String;


typedef FixedVector<double,6>   Vec6;
typedef FixedMatrix<double,6,6> Mat66;
typedef CDSMatrix<double>       RMat;
typedef CDSVector<double,1>     RVec;


/**
 * This is an abstract class representing a body and its (generic) inboard joint, that is,
 * the joint connecting it to its parent. Concrete classes are derived from this one to
 * represent each specific type of joint.
 *
 * RigidBodyNodes are linked into a tree structure, organized into levels as described 
 * in Schwieters' JMR paper. The root is a special 'Ground' node defined to be at 
 * level 0 and containing only atoms fixed to ground. The Level 1 nodes (referred to
 * as 'base nodes') are those attached directly to the Ground node, level 2's attach 
 * to level 1's, etc. Every node but Ground has exactly one parent node, whose
 * level is always one less than the current node. Any node may have an arbitrary number 
 * of children, for which it is the unique parent, and all of its children have 
 * level one greater than the current node.
 *
 * Every RigidBodyNode contains a list of pointers to the atoms rigidly affixed to the
 * body associated with that node. The 0th atom in that list is the one which is
 * connected to the node's parent by a joint, to a specific atom on the parent
 * node stored as 'parentAtom' here. Our 0th atom's body-fixed station serves as
 * the origin for the body frame, and thus has coordinates [0,0,0] in that frame.
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
 * each body (translated to atom 0) when we first see the bodies in the reference
 * frame. Later the bodies move and take their B frames with them. The locations
 * of all the atoms on a body are then forever fixed when measured in the B frame;
 * we call such points 'stations'. For convenience, we refer to the body frame
 * of a body's unique parent as the 'P' frame. Initially frames P and B are
 * aligned (and both are aligned with Ground). Later the P and B frames
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
 * of the inboard joint's frame in the body frame. If B is the i'th child of
 * its parent P, then there is a parent-fixed frame Ji which is the image of J
 * when the joint coordinates are zero. That is, R_JiJ is the orientation change
 * induced by the joint coordinates. Note that because of how we define the 
 * reference configuration, R_PJi = R_BJ so we don't need to store both matrices
 * explicitly. With these definitions we can easily calculate R_PB as
 *     R_PB = R_PJi*R_JiJ*R_JB = R_BJ*R_JiJ*(R_BJ)'.
 *
 * XXX Sherm's goal: migrate this to atom-free living.
 */
class RigidBodyNode {
public:
    class VirtualBaseMethod {};    // an exception

    virtual ~RigidBodyNode() {}

    RigidBodyNode(const MassProperties& mProps_B,
                  const Vec3& originOfB_P, // and R_BP=I in ref config
                  const Mat33& rot_BJ, const Vec3& originOfJ_B)
      : stateOffset(-1), parent(0), children(0,0), level(-1),
        massProps_B(mProps_B), inertia_CB_B(mProps_B.calcCentroidalInertia()),
        R_BJ(rot_BJ), OJ_B(originOfJ_B), refOrigin_P(originOfB_P)
    {
        R_GB.set(0.);       // B frame is initially aligned with Ground
        R_GB.setDiag(1.);
    }

    RigidBodyNode& operator=(const RigidBodyNode&);

    const RigidBodyNode* getParent() const {return parent;}
    RigidBodyNode*       updParent()       {return parent;}

    /// Return this node's level, that is, how many ancestors separate it from
    /// the Ground node at level 0. Level 1 nodes (directly connected to the
    /// Ground node) are called 'base' nodes.
    int              getLevel() const  {return level;}
    void             setLevel(int i)   {level=i;}

    bool             isGroundNode() const { return level==0; }
    bool             isBaseNode()   const { return level==1; }

    int              getStateOffset() const {return stateOffset;}

    const MassProperties& getMassProperties() const {return massProps_B;}
    const double&  getMass()         const {return massProps_B.getMass();}
    const Vec3&    getCOM_B()        const {return massProps_B.getCOM();}
    const Inertia& getInertia_OB_B() const {return massProps_B.getInertia();}
    const Inertia& getInertia_CB_B() const {return inertia_CB_B;}

    /// Return R_GB, the rotation (direction cosine) matrix giving the 
    /// spatial orientation of this body's frame B (that is, B's orientation
    /// in the ground frame G).
    const Mat33&     getR_GB()   const {return R_GB;}

    /// Return R_GP, the rotation (direction cosine) matrix giving the
    /// orientation of this body's *parent's* body frame (which we'll call
    /// P here) in the ground frame G.
    const Mat33&     getR_GP()   const {assert(parent); return parent->getR_GB();}

    const Vec6&      getSpatialVel()   const {return sVel;}
    const Vec6&      getSpatialAcc()   const {return sAcc;}

    const PhiMatrix& getPhi()    const {return phi;}
    const Mat66&     getPsiT()   const {return psiT;}
    const Mat66&     getY()      const {return Y;}

    virtual const char* type() { return "unknown"; }
    virtual int getDOF() const {return 0;} //number of independent dofs
    virtual int getDim() const {return 0;} //dofs plus constraints
    virtual double kineticE() { return 0; }
    virtual double approxKE() { return 0; }

    virtual void calcP() {throw VirtualBaseMethod();}
    virtual void calcZ() {throw VirtualBaseMethod();}
    virtual void calcPandZ() {throw VirtualBaseMethod();}
    virtual void calcY() {throw VirtualBaseMethod();}
    virtual void calcInternalForce()        {throw VirtualBaseMethod();}
    virtual void prepareVelInternal()       {throw VirtualBaseMethod();}
    virtual void propagateSVel(const Vec6&) {throw VirtualBaseMethod();}

    virtual void setPosVel(const RVec&, const RVec&) {throw VirtualBaseMethod();}
    virtual void setVel(const RVec&)    {throw VirtualBaseMethod();}
    virtual void setVelFromSVel(const Vec6&) {throw VirtualBaseMethod();}
    virtual void enforceConstraints(RVec& pos, RVec& vel) {throw VirtualBaseMethod();}

    virtual void getPos(RVec&) {throw VirtualBaseMethod();}
    virtual void getVel(RVec&) {throw VirtualBaseMethod();}
    virtual void getAccel(RVec&) {throw VirtualBaseMethod();}
    virtual void calcAccel() {throw VirtualBaseMethod();}
    virtual void getInternalForce(RVec&) {throw VirtualBaseMethod();}
    virtual RMat getH() {throw VirtualBaseMethod();}

    virtual void print(int) { throw VirtualBaseMethod(); }

protected:
    typedef CDSList<RigidBodyNode*>   RigidBodyNodeList;

    int               stateOffset;  //index into internal coord pos,vel,acc arrays
    RigidBodyNode*    parent; 
    RigidBodyNodeList children;
    int               level;        //how far from base 

    // These are the body properties

    // Fixed forever in the body-local frame B:
    //      ... supplied on construction
    const MassProperties massProps_B;
    const Mat33  R_BJ;          // orientation of inboard joint frame J, in B
    const Vec3   OJ_B;          // origin of J, measured from B origin, expr. in B

    // Reference configuration. This is the body frame origin location, measured
    // in its parent's frame in the reference configuration. This vector is fixed
    // after construction! The body origin can of course move relative to its
    // parent, but that is not the meaning of this reference configuration vector.
    // (Note however that the body origin is also the location of the inboard joint, 
    // meaning that the origin point moves relative to the parent only due to translations.)
    // Note that by definition the orientation of the body frame is identical to P
    // in the reference configuration so we don't need to store it.
    const Vec3   refOrigin_P;

    //      ... calculated on construction
    const Inertia inertia_CB_B;  // centroidal inertia, expr. in B

    // Calculated spatial quantities
    //      ... position level

    //      Rotation (direction cosine matrix) expressing the body frame B
    //      in the ground frame G. That is, if you have a vector vB expressed
    //      body frame and want it in ground, use vG = R_GB*vB. 
    Mat33        R_GB;
    Vec3         OB_G;          // origin of B meas & expr in G

    Vec3         COMstation_G;  // measured from B origin, expr. in G
    Mat33        inertia_OB_G;  // about B's origin, expr. in G

    PhiMatrix    phi;           // spatial rigid body transition matrix
    Mat66        Mk;            // spatial inertia matrix
    Mat66        P;             // articulated body spatial inertia
    Mat66        tau;
    Mat66        psiT;
    Mat66        Y;             // diag of Omega - for loop constraints

    //      ... velocity level
    Vec6         a;     // spatial coriolis acceleration
    Vec6         b;     // spatial gyroscopic force
    Vec6         sVel;  // spatial velocity

    //      ... acceleration level
    Vec6         forceCartesian;    // net spatial force
    Vec6         z;
    Vec6         Gepsilon;
    Vec6         sAcc;              // spatial acceleration


    static const double DEG2RAD; //using angles in degrees balances gradient

    virtual void velFromCartesian() {}

    friend ostream& operator<<(ostream& s, const RigidBodyNode&);
    template<int dof> friend class RigidBodyNodeSpec;
};

#endif /* RIGID_BODY_NODE_H_ */
