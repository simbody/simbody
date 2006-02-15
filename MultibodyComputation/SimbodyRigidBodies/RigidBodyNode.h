#ifndef RIGID_BODY_NODE_H_
#define RIGID_BODY_NODE_H_

#include "simbody/Simbody.h"
#include "SimbodyTree.h"
#include "SimbodyTreeState.h"

#include "internalDynamics.h"

#include <cassert>
#include <vector>

using namespace simtk;

/**
 * This is an abstract class representing a body and its (generic) inboard joint, that is,
 * the joint connecting it to its parent. Concrete classes are derived from this one to
 * represent each specific type of joint.
 *
 * RigidBodyNodes are linked into a tree structure, organized into levels as described 
 * in Schwieters' JMR paper. The root is a special 'Ground' node defined to be at 
 * level 0. The level 1 nodes (referred to
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
 * other direction use R_YX=~R_XY (tilde '~' is SimTK for "transpose"). This convention
 * provides flawless composition of rotations as long as you "match up" the frames
 * (i.e., adjacent frame symbols must be the same). So if you have R_XY and R_ZX 
 * and you want R_YZ you can easily get it like this:
 *     R_YZ = R_YX*R_XZ = (~R_XY)*(~R_ZX) = ~(R_ZX*R_XY).
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
 *     R_PB = R_PJi*R_JiJ*R_JB = R_BJ*R_JiJ*~R_BJ.
 */
class RigidBodyNode {
public:
    class VirtualBaseMethod {};    // an exception

    virtual ~RigidBodyNode() {}

    RigidBodyNode& operator=(const RigidBodyNode&);

    /// Factory for producing concrete RigidBodyNodes based on joint type.
    static RigidBodyNode* create(
        const MassProperties&   m,            // mass properties in body frame
        const TransformMat&     jointFrame,   // inboard joint frame J in body frame
        Joint::JointType        type,
        bool                    isReversed,   // child-to-parent orientation?
        bool                    useEuler,     // TODO: kludge (true if minimizing)
        int&                    nxtU,
        int&                    nxtUSq,
        int&                    nxtQ); 

    /// Register the passed-in node as a child of this one.
    void addChild(RigidBodyNode* child);

        // TOPOLOGICAL INFO: no State needed

    RigidBodyNode*   getParent() const {return parent;}
    void             setParent(RigidBodyNode* p) { parent=p; }

    int              getNChildren()  const {return (int)children.size();}
    RigidBodyNode*   getChild(int i) const {return (i<(int)children.size()?children[i]:0);}

    /// Return this node's level, that is, how many ancestors separate it from
    /// the Ground node at level 0. Level 1 nodes (directly connected to the
    /// Ground node) are called 'base' nodes.
    int              getLevel() const  {return level;}
    void             setLevel(int i)   {level=i;}

    int              getNodeNum() const {return nodeNum;}
    void             setNodeNum(int n)  {nodeNum=n;}

    bool             isGroundNode() const { return level==0; }
    bool             isBaseNode()   const { return level==1; }

    int              getUIndex() const {return uIndex;}
    int              getQIndex() const {return qIndex;}

        // MODELING INFO
    const bool getUseEulerAngles(const SBState& s) const {return s.vars->useEulerAngles;}
    const bool isPrescribed     (const SBState& s) const {return s.vars->prescribed[nodeNum];}

        // PARAMETRIZATION INFO

    // TODO: These ignore State currently since they aren't parametrizable.
    const MassProperties& getMassProperties(const SBState&) const {return massProps_B;}
    const Real&           getMass          (const SBState&) const {return massProps_B.getMass();}
    const Vec3&           getCOM_B         (const SBState&) const {return massProps_B.getCOM();}
    const InertiaMat&     getInertia_OB_B  (const SBState&) const {return massProps_B.getInertia();}
    const InertiaMat&     getInertia_CB_B  (const SBState&) const {return inertia_CB_B;}

        // CONFIGURATION INFO


    /// Extract from the cache X_GB, the transformation matrix giving the spatial configuration of this
    /// body's frame B measured from and expressed in ground. This consists of a rotation matrix
    /// R_GB, and a ground-frame vector OB_G from ground's origin to the origin point of frame B.
    const TransformMat& getX_GB(const SBState& s) const {return s.cache->bodyConfigInGround[nodeNum];}
    TransformMat&       updX_GB(const SBState& s) const {return s.cache->bodyConfigInGround[nodeNum];}

    /// Extract from the cache  X_PB, the cross-joint transformation matrix giving the configuration
    /// of this body's frame B measured from and expressed in its *parent* frame. Thus this is NOT
    /// a spatial transformation.
    const TransformMat& getX_PB(const SBState& s) const {return s.cache->bodyConfigInParent[nodeNum];}
    TransformMat&       updX_PB(const SBState& s) const {return s.cache->bodyConfigInParent[nodeNum];}

    /// Extract from the cache the body-to-parent shift matrix "phi". 
    const PhiMatrix&    getPhi(const SBState& s) const {return s.cache->bodyToParentShift[nodeNum];}
    PhiMatrix&          updPhi(const SBState& s) const {return s.cache->bodyToParentShift[nodeNum];}

    /// Extract this body's spatial inertia matrix from the cache. This contains the mass properties
    /// measured from (and about) the body frame origin, but expressed in the *ground* frame.
    const SpatialMat&   getMk(const SBState& s) const {return s.cache->bodySpatialInertia[nodeNum];}
    SpatialMat&         updMk(const SBState& s) const {return s.cache->bodySpatialInertia[nodeNum];}

    /// Extract from the cache the location of the body's center of mass, measured from the ground
    /// origin and expressed in ground.
    const Vec3& getCOM_G(const SBState& s) const {return s.cache->bodyCOMInGround[nodeNum];}
    Vec3&       updCOM_G(const SBState& s) const {return s.cache->bodyCOMInGround[nodeNum];}

    /// Extract from the cache the vector from body B's origin to its center of mass, reexpressed in Ground.
    const Vec3& getCB_G(const SBState& s) const {return s.cache->bodyCOMStationInGround[nodeNum];}
    Vec3&       updCB_G(const SBState& s) const {return s.cache->bodyCOMStationInGround[nodeNum];}

    /// Extract from the cache the body's inertia about the body origin OB, but reexpressed in Ground.
    const InertiaMat& getInertia_OB_G(const SBState& s) const {return s.cache->bodyInertiaInGround[nodeNum];}
    InertiaMat&       updInertia_OB_G(const SBState& s) const {return s.cache->bodyInertiaInGround[nodeNum];}

    /// Return R_GB, the rotation (direction cosine) matrix giving the 
    /// spatial orientation of this body's frame B (that is, B's orientation
    /// in the ground frame G).
    const RotationMat& getR_GB(const SBState& s) const {return getX_GB(s).getRotation();}

    /// Return OB_G, the spatial location of the origin of the B frame, that is, 
    /// measured from the ground origin and expressed in ground.
    const Vec3&        getOB_G(const SBState& s) const {return getX_GB(s).getTranslation(); }

    /// Return R_GP, the rotation (direction cosine) matrix giving the
    /// orientation of this body's *parent's* body frame (which we'll call
    /// P here) in the ground frame G.
    const RotationMat& getR_GP(const SBState& s) const {assert(parent); return parent->getR_GB(s);}

    /// Return OP_G, the spatial location of the origin of the P frame, that is, 
    /// measured from the ground origin and expressed in ground.
    const Vec3&        getOP_G(const SBState& s) const {assert(parent); return parent->getOB_G(s);}

            // VELOCITY INFO

    /// Extract from the cache V_GB, the spatial velocity of this body's frame B measured in and
    /// expressed in ground. This contains the angular velocity of B in G, and the linear velocity
    /// of B's origin point OB in G, with both vectors expressed in G.
    const SpatialVec& getV_GB   (const SBState& s) const {return s.cache->bodyVelocityInGround[nodeNum];}
    SpatialVec&       updV_GB   (const SBState& s) const {return s.cache->bodyVelocityInGround[nodeNum];}

    /// Extract from the cache V_PB_G, the *spatial* velocity of this body's frame B, that is the
    /// cross-joint velocity measured with respect to the parent frame, but then expressed in the
    /// *ground* frame. This contains the angular velocity of B in P, and the linear velocity
    /// of B's origin point OB in P, with both vectors expressed in *G*.
    const SpatialVec& getV_PB_G (const SBState& s) const {return s.cache->bodyVelocityInParent[nodeNum];}
    SpatialVec&       updV_PB_G (const SBState& s) const {return s.cache->bodyVelocityInParent[nodeNum];}

    const SpatialVec& getSpatialVel   (const SBState& s) const {return getV_GB(s);}
    const Vec3&       getSpatialAngVel(const SBState& s) const {return getV_GB(s)[0];}
    const Vec3&       getSpatialLinVel(const SBState& s) const {return getV_GB(s)[1];}

        // DYNAMICS INFO

    const SpatialVec& getBodyForce(const SBState& s) const {return s.vars->appliedBodyForces[nodeNum];}
 
    const SpatialVec& getCoriolisAcceleration(const SBState& s) const {return s.cache->coriolisAcceleration[nodeNum];}
    SpatialVec&       updCoriolisAcceleration(const SBState& s) const {return s.cache->coriolisAcceleration[nodeNum];}
 
    const SpatialVec& getGyroscopicForce(const SBState& s) const {return s.cache->gyroscopicForces[nodeNum];}
    SpatialVec&       updGyroscopicForce(const SBState& s) const {return s.cache->gyroscopicForces[nodeNum];}
    
    /// Extract from the cache A_GB, the spatial acceleration of this body's frame B measured in and
    /// expressed in ground. This contains the inertial angular acceleration of B in G, and the
    /// linear acceleration of B's origin point OB in G, with both vectors expressed in G.
    const SpatialVec& getA_GB (const SBState& s) const {return s.cache->bodyAccelerationInGround[nodeNum];}
    SpatialVec&       updA_GB (const SBState& s) const {return s.cache->bodyAccelerationInGround[nodeNum];}

    const SpatialVec& getSpatialAcc   (const SBState& s) const {return getA_GB(s);}
    const Vec3&       getSpatialAngAcc(const SBState& s) const {return getA_GB(s)[0];}
    const Vec3&       getSpatialLinAcc(const SBState& s) const {return getA_GB(s)[1];}

    const SpatialMat& getP    (const SBState& s) const {return s.cache->articulatedBodyInertia[nodeNum];}
    SpatialMat&       updP    (const SBState& s) const {return s.cache->articulatedBodyInertia[nodeNum];}

    const SpatialVec& getZ(const SBState& s) const {return s.cache->z[nodeNum];}
    SpatialVec&       updZ(const SBState& s) const {return s.cache->z[nodeNum];}

    const SpatialVec& getGepsilon(const SBState& s) const {return s.cache->Gepsilon[nodeNum];}
    SpatialVec&       updGepsilon(const SBState& s) const {return s.cache->Gepsilon[nodeNum];}

    const SpatialMat&  getPsiT(const SBState& s) const {return s.cache->psiT[nodeNum];}
    SpatialMat&        updPsiT(const SBState& s) const {return s.cache->psiT[nodeNum];}

    const SpatialMat&  getTau(const SBState& s) const {return s.cache->tau[nodeNum];}
    SpatialMat&        updTau(const SBState& s) const {return s.cache->tau[nodeNum];}

    const SpatialMat&  getY(const SBState& s) const {return s.cache->Y[nodeNum];}
    SpatialMat&        updY(const SBState& s) const {return s.cache->Y[nodeNum];}

    virtual void realizeModeling  (const SBState&) const=0;
    virtual void realizeParameters(const SBState&) const=0;

    /// Introduce new values for generalized coordinates and calculate
    /// all the position-dependent kinematic terms.
    virtual void realizeConfiguration(const SBState&)=0;

    /// Introduce new values for generalized speeds and calculate
    /// all the velocity-dependent kinematic terms. Assumes realizeConfiguration()
    /// has already been called.
    virtual void realizeVelocity(const SBState&)=0;

    Real calcKineticEnergy(const SBState&) const;   // from spatial quantities only

    virtual const char* type()     const {return "unknown";}
    virtual int         getDOF()   const {return 0;} //number of independent dofs
    virtual int         getMaxNQ() const {return 0;} //dofs plus quaternion constraints

    virtual void enforceQuaternionConstraints(const SBState&) {throw VirtualBaseMethod();}

    virtual void calcP(const SBState&)                                 {throw VirtualBaseMethod();}
    virtual void calcZ(const SBState&, const SpatialVec& spatialForce) {throw VirtualBaseMethod();}
    virtual void calcY(const SBState&)                                 {throw VirtualBaseMethod();}
    virtual void calcAccel(const SBState&)                             {throw VirtualBaseMethod();}

    virtual void calcInternalGradientFromSpatial(const SBState&, Vector_<SpatialVec>& zTmp,
                                                 const Vector_<SpatialVec>& X, Vector& JX)
      { throw VirtualBaseMethod(); }

    virtual void setVelFromSVel(SBState&, const SpatialVec&) {throw VirtualBaseMethod();}

    virtual void getDefaultParameters   (SBState&) const {throw VirtualBaseMethod();}
    virtual void getDefaultConfiguration(SBState&) const {throw VirtualBaseMethod();}
    virtual void getDefaultVelocity     (SBState&) const {throw VirtualBaseMethod();}

    virtual void getAccel(Vector&) const {throw VirtualBaseMethod();}

    virtual void getInternalForce(const SBState&, Vector&) const {throw VirtualBaseMethod();}

    // Note that this requires rows of H to be packed like SpatialRow.
    virtual const SpatialRow& getHRow(const SBState&, int i) const {throw VirtualBaseMethod();}

    virtual void print(int) const { throw VirtualBaseMethod(); }

    void nodeDump(std::ostream&) const;
    virtual void nodeSpecDump(std::ostream& o) const { o<<"NODE SPEC type="<<type()<<std::endl; }

    static const double DEG2RAD; //using angles in degrees balances gradient

protected:
    /// This is the constructor for the abstract base type for use by the derived
    /// concrete types in their constructors.
    RigidBodyNode(const MassProperties& mProps_B,
                  const Vec3&           originOfB_P, // and R_BP=I in ref config
                  const TransformMat&   xform_BJ)
      : uIndex(-1), qIndex(-1), uSqIndex(-1), parent(0), children(), level(-1), nodeNum(-1),
        massProps_B(mProps_B), inertia_CB_B(mProps_B.calcCentroidalInertia()),
        X_BJ(xform_BJ), refOrigin_P(originOfB_P)
    {
       // X_PB=TransformMat(RotationMat(),refOrigin_P);
       // V_PB_G=0; sVel=0; sAcc=0;
      //  X_GB=TransformMat();
       // COM_G = 0.; COMstation_G = massProps_B.getCOM();
      //  phi = PhiMatrix(Vec3(0));
       // psiT=0; P=0; z=0; tau=0; Gepsilon=0; Y=0;
    }

    typedef std::vector<RigidBodyNode*>   RigidBodyNodeList;

    int               uIndex;   // index into internal coord vel,acc arrays
    int               qIndex;   // index into internal coord pos array
    int               uSqIndex; // index into array of DOF^2 objects

    RigidBodyNode*    parent; 
    RigidBodyNodeList children;
    int               level;        //how far from base 
    int               nodeNum;      //unique ID number in RigidBodyTree

    // These are the body properties

    // Fixed forever in the body-local frame B:
    //      ... supplied on construction
    const MassProperties massProps_B;
    const TransformMat   X_BJ; // orientation and location of inboard joint frame J
                               //   measured & expressed in body frame B

    // Reference configuration. This is the body frame origin location, measured
    // in its parent's frame in the reference configuration. This vector is fixed
    // after this node is connected to its parent. The body origin can of course move relative
    // to its parent, but that is not the meaning of this reference configuration vector.
    // (Note however that the body origin is also the location of the inboard joint, 
    // meaning that the origin point moves relative to the parent only due to translations.)
    // Note that by definition the orientation of the body frame is identical to P
    // in the reference configuration so we don't need to store it.
    // TODO: this is frame Ji
    Vec3 refOrigin_P;

    //      ... calculated on construction
    const InertiaMat inertia_CB_B;  // centroidal inertia, expr. in B

    virtual void velFromCartesian() {}

    friend std::ostream& operator<<(std::ostream& s, const RigidBodyNode&);
    template<int dof> friend class RigidBodyNodeSpec;

private:   
    /// Calculate all spatial configuration quantities, assuming availability of
    /// joint-specific relative quantities.
    ///   X_GB
    void calcJointIndependentKinematicsPos(const SBState&);

    /// Calcluate all spatial velocity quantities, assuming availability of
    /// joint-specific relative quantities and all position kinematics.
    ///   sVel  spatial velocity of B
    ///   a     spatial Coriolis acceleration
    ///   b     spatial gyroscopic force
    void calcJointIndependentKinematicsVel(const SBState&);
};

#endif // RIGID_BODY_NODE_H_
