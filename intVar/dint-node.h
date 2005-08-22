#ifndef __dint_node_hh__
#define __dint_node_hh__

#include "phiMatrix.h"
#include "internalDynamics.h"

#include <fixedMatrix.h>
#include <cdsMatrix.h>
#include <cdsVector.h>
#include <cdsList.h>

#include <cassert>

class Vec3;
class IVMAtom;
template<class CHAR> class CDSString;
typedef CDSString<char> String;

typedef CDSList< IVMAtom* >    AtomList;

typedef FixedVector<double,6> Vec6;
typedef FixedMatrix<double,6> Mat66;
typedef CDSMatrix<double>     RMat;
typedef CDSVector<double,1>   RVec;

class InertiaTensor : public Mat33 {
public:
    InertiaTensor() : Mat33(0.0) {}
    //  InertiaTensor(const InertiaTensor&);
    void calc(const Vec3&     center,
              const AtomList&       );
};

class IVM;

/**
 * This is an abstract class representing a body and its (generic) inboard joint, that is,
 * the joint connecting it to its parent. Concrete classes are derived from this one to
 * represent each specific type of joint.
 *
 * HingeNodes are linked into a tree structure, organized into levels as described 
 * in Schwieters' JMR paper. The root is a special 'Ground' node defined to be at 
 * level 0 and containing only atoms fixed to ground. The Level 1 nodes (referred to
 * as 'base nodes') are those attached directly to the Ground node, level 2's attach 
 * to level 1's, etc. Every node but Ground has exactly one parent node, whose
 * level is always one less than the current node. Any node may have an arbitrary number 
 * of children, for which it is the unique parent, and all of its children have 
 * level one greater than the current node.
 *
 * Every HingeNode contains a list of pointers to the atoms rigidly affixed to the
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
 */
class HingeNode {
public:
    const IVM*    ivm;
    const IVMAtom*   parentAtom;   // atom in parent to which hinge is attached
    class VirtualBaseMethod {};    // an exception

    virtual ~HingeNode() {}
    HingeNode(const IVM*        ivm,
              IVMAtom*          hingeAtom,
              const IVMAtom*    parentAtom,
              HingeNode*        parentNode);
    HingeNode& operator=(const HingeNode&);

    HingeNode*       getParent() const {return parent;}

    /// Return R_GB, the rotation (direction cosine) matrix giving the 
    /// spatial orientation of this body's frame B (that is, B's orientation
    /// in the ground frame G).
    const Mat33&     getR_GB()   const {return R_GB;}

    /// Return R_GP, the rotation (direction cosine) matrix giving the
    /// orientation of this body's *parent's* body frame (which we'll call
    /// P here) in the ground frame G.
    const Mat33&     getR_GP()   const {assert(parent); return parent->getR_GB();}

    /// Return this node's level, that is, how many ancestors separate it from
    /// the Ground node at level 0. Level 1 nodes (directly connected to the
    /// Ground node) are called 'base' nodes.
    int              getLevel()  const {return level;}

    bool             isGroundNode() const { return level==0; }
    bool             isBaseNode()   const { return level==1; }


    const Vec6&      getSpatialVel()   const {return sVel;}
    const Vec6&      getSpatialAcc()   const {return sAcc;}

    const PhiMatrix& getPhi()    const {return phi;}
    const Mat66&     getPsiT()   const {return psiT;}
    const Mat66&     getY()      const {return Y;}

    const IVMAtom*   getAtom(int i)   const {return (i<atoms.size()?atoms[i]:0); }
    const HingeNode* getChild(int i) const {return (i<children.size()?children[i]:0);}

    void addChild(HingeNode*);

    virtual int           offset() const {return 0;}
    virtual const Vec3&   posCM()  const {throw VirtualBaseMethod();}
    virtual const double& mass()   const {throw VirtualBaseMethod();}

    virtual void calcP() {throw VirtualBaseMethod();}
    virtual void calcZ() {throw VirtualBaseMethod();}
    virtual void calcPandZ() {throw VirtualBaseMethod();}
    virtual void calcY() {throw VirtualBaseMethod();}
    virtual void calcInternalForce()        {throw VirtualBaseMethod();}
    virtual void prepareVelInternal()       {throw VirtualBaseMethod();}
    virtual void propagateSVel(const Vec6&) {throw VirtualBaseMethod();}

    virtual int getDOF() const {return 0;} //number of independent dofs
    virtual int getDim() const {return 0;} //dofs plus constraints
    virtual double kineticE() { return 0; }
    virtual double approxKE() { return 0; }
    virtual void setPosVel(const RVec&, const RVec&) {throw VirtualBaseMethod();}
    virtual void setVel(const RVec&)    {throw VirtualBaseMethod();}
    virtual void setVelFromSVel(const Vec6&) {throw VirtualBaseMethod();}
    virtual void enforceConstraints(RVec& pos, RVec& vel) {throw VirtualBaseMethod();}
    virtual const char* type() { return "unknown"; }
    virtual void print(int) { throw VirtualBaseMethod(); }
    virtual void getPos(RVec&) {throw VirtualBaseMethod();}
    virtual void getVel(RVec&) {throw VirtualBaseMethod();}
    virtual void getAccel(RVec&) {throw VirtualBaseMethod();}
    virtual void calcAccel() {throw VirtualBaseMethod();}
    virtual void getInternalForce(RVec&) {throw VirtualBaseMethod();}
    virtual RMat getH() {throw VirtualBaseMethod();}

protected:
    typedef CDSList<HingeNode*>   HingeNodeList;

    HingeNode*    parent; 
    HingeNodeList children;
    int           level;        //how far from base
    AtomList      atoms; 

    Vec6      sVel;        //spatial velocity
    Vec6      sAcc;        //spatial acceleration
    PhiMatrix phi;

    Mat66 psiT;
    Mat66 P;
    Vec6  z;
    Mat66 tau;
    Vec6  Gepsilon;

    // Rotation (direction cosine matrix) expressing the body frame B
    // in the ground frame G. That is, if you have a vector vB expressed
    // body frame and want it in ground, use vG = R_GB*vB. 
    Mat33 R_GB;

    Mat66 Y;  //diagonal components of Omega- for loop constraints

    static const double DEG2RAD; //using angles in degrees balances gradient

    virtual void velFromCartesian() {}
    void removeHinge();

    friend ostream& operator<<(ostream& s, const HingeNode&);
    friend HingeNode* construct(HingeNode*                         oldNode,
                                const InternalDynamics::HingeSpec& type,
                                int&                               cnt);
    friend void combineNodes(const HingeNode* node1,
                             const HingeNode* node2);

    //  static void groupTorsion(const HingeNode*);
    template<int dof> friend class HingeNodeSpec;
    friend class AtomTree;
    friend class AT_Build;
};

#endif /* __dint_node_hh__ */
