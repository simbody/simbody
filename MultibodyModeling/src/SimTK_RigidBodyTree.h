#ifndef SIMTK_RIGIDBODY_TREE_H_
#define SIMTK_RIGIDBODY_TREE_H_

#include "SimTK_Frame.h"

#include <vector>

class Vec3;
class Mat33;
class Joint;
class RigidBody;


/**
 * Built-in joint types.
 *
 * Each joint connects two frames. One is on the "Reference" body R
 * and the other is on the "Moving" body M. In a tree system, the
 * parent body is always R.
 *
 * Joints have from 0-6 coordinates. When all coordinates are zero,
 * the two defining frames R and M are aligned, and their origins 
 * OR and OM are coincident.
 *
 * Reverse joints exist for all joint types allowing a user to think
 * of the roles of R and M reversed, while we can still build a tree
 * which follows the rule that R is always the parent.
 *
 * Weld        0 Frames remain coincident forever
 *
 * Torsion     1 Mz=Rz, OM=OR forever. Coord is angle from
 *   Pin           Rx to Mx (right hand rule about z).
 *
 * Sliding     1 M=R, OMx=OMy=0 forever. Coord is (OM-OR)*Rz.
 *   Prismatic, Stretch
 *
 * Universal   2 OM=OR forever. Coordinates are body fixed 1-2
 *                 sequence about Mx=Rx, then new My
 *
 * Cylinder    2 Mz=Rz, OMx=OMy=0 forever. 1st coord is same
 *                 as Sliding, second same as Torsion.
 *
 * Planar      3 Mz=Rz, OMz=0 forever. 1st coords are translation
 *                 in x,y; 3rd is same as Torsion.
 *
 * Gimbal      3 OM=OR forever. Coordinates are 1-2-3 body fixed
 *                 Euler sequence.
 *
 * Orientation 3 OM=OR forever. Coords represent orientation
 *   Ball          of M in R as 3 Euler angles or 4 Quaternions.
 *   Spherical
 *
 * Cartesian   3 M=R forever. Coords are (OM-OR)*R, i.e. x,y,z.
 *   Translational, FreePoint
 *
 * FreeLine    5 UJoint plus Cartesian. Mz should be aligned
 *                 with the mass distribution axis.
 * 
 * Free        6 Coords are Orientation and Cartesian
 * 
 */
enum XXJointType {
    UnknownJointType    = 0,
    WeldJoint           = 1,
    TorsionJoint        = 2,
    SlidingJoint        = 3,
    UJoint              = 4,
    CylinderJoint       = 5,
    PlanarJoint         = 6,
    GimbalJoint         = 7,
    OrientationJoint    = 8,
    CartesianJoint      = 9,
    FreeLineJoint       = 10,
    FreeJoint           = 11
};




/**
 * XXRigidBodyTree is a multibody system restricted to unconstrained tree topology
 * and rigid bodies. That is, there are no loops or constraints, and body stations,
 * frames, and mass properties are constants when expressed in the body frame.
 *
 * This class contains the system topology, that is, Ground, the bodies and
 * their inboard joints. Ground is body #0, other bodies have a regular labeling
 * such that body numbers are strictly decreasing as you go from a body to its
 * parent, parent's parent, and so on down to Ground. A body and its inboard
 * joint have the same number, except that Ground does not have an inboard joint
 * so there is no joint 0.
 *
 * Bodies have mass properties, stations, vectors, and frames. One frame and
 * one station are distinguished. Those are called the body frame B and the
 * body origin OB. By definition, the body frame has measure numbers 100,010,001 and
 * body origin is a station with measure numbers 000. All other stations and 
 * vectors are measured with respect to B and OB and expressed in B.
 * 
 * Each body is defined in terms of its parent body P, which must already have
 * been defined. That is, the orientation of B is presumed to be the same
 * as the orientation of P when the inboard joint is in its reference configuration
 * (all generalized coordinates zero, where "zero" for a quaternion is 0001).
 *
 * Equations of motion are:
 *         qdot = Q(q)u
 *    M(q) udot = T(t,q,u) - C(q,u)
 * where T are the user-applied forces converted to their equivalent joint
 * force system and C are the bias forces produced by Coriolis and gyroscopic
 * effects that result from velocities. The quantity Tbar = T-C is called the
 * "bias free torque".
 *
 * Jk(q) is the kinematic Jacobian (a.k.a partial velocities) i.e. 
 * partial(V)/partial(u) where V contains the spatial velocity of each body
 * at the body origin or a specified station.
 *
 * We provide the following operators:
 *    T = Jk*F         convert spatial forces to equivalent joint forces
 *    w = M*v          mass times a user-supplied vector v (inverse dynamics)
 *    w = inv(M)*v     mass inverse times user-supplied vector v (fwd dynamics)   
 */
class XXRigidBodyTree { 
public:
    XXRigidBodyTree();
    XXRigidBodyTree(const XXRigidBodyTree&);
    XXRigidBodyTree& operator=(const XXRigidBodyTree&);
    ~XXRigidBodyTree();

    /// Add a new body connected to one of the existing bodies by
    /// a Joint. Body 0 (Ground) is predefined at construction.
    int addBody(const MassProperties&, const TransformMat& frameInB, 
                int   parent,          const TransformMat& frameInP,
                const XXJointType&);

    /// Return number of defined bodies, including ground. This is
    /// one more than the number of (tree) joints.
    int getNBodies() const;

    /// Evaluate the quantities appropriate for the indicated stage,
    /// allowing access to responses, operators and solvers for that stage.
    void realize(const State&, Stage, const RealizeOptions&) const;

    // Stage Configure

    /// Obtain the spatial location and orientation of the body frame, by
    /// returning the frame in Ground which is currently overlaid on the body frame.
    const TransformMat& getBodyFrameInGround(const State&, int body) const;

    /// Obtain the location and orientation of the body frame in its parent frame, by
    /// returning the frame in Parent which is currently overlaid on the body frame.
    const TransformMat& getBodyFrameInParent(const State&, int body) const;

    inline Vec3  xformVector2Ground (const State&, int body, const Vec3&  vB) const;
    inline Vec3  xformStation2Ground(const State&, int body, const Vec3&  sB) const;
    inline TransformMat xformFrame2Ground  (const State&, int body, const TransformMat& fB) const;

    inline Vec3  xformVector2Parent (const State&, int body, const Vec3&  vB) const;
    inline Vec3  xformStation2Parent(const State&, int body, const Vec3&  sB) const;
    inline TransformMat xformFrame2Parent  (const State&, int body, const TransformMat& fB) const;

    /// Return mass properties of subtree rooted at this body, as though
    /// all joints were welded.
    const MassProperties& getCompositeMassProperties(const State&, int body) const;

    /// Return mass properties of subtree rooted at this body, with all
    /// joints floppy and unactuated.
    const MassProperties& getArticulatedMassProperties(const State&, int body) const;

    void addInPointForce(const State& s, int body, const Vec3& station, const Vec3& forceG,
                         SpatialVector& spatialForce) const
    {
        spatialForce[0] += xformVector2Ground(s,body,station) % forceG;
        spatialForce[1] += forceG; 
    }

    void addInMoment    (const State&, int body, const Vec3& momentG,
                         SpatialVector& spatialForce) const 
      { spatialForce[0] += momentG; }

    /// Add in a spatial force for each body as produced by a uniform gravitational field.
    void addInGravity(const State&, const Vec3& gG, Vector<SpatialVector>&) const;
                            
    void convertSpatialForcesToJointForces
        (const State&, const Vector<SpatialVector>& bodyForces, Vector& jointForces);
                                

    void calcMassXVector(const State&, const Vector&, Vector&) const;
    void calcInvMassXVector(const State&, const Vector&, Vector&) const;

    /// Return partial(stationVelocity)/partial(u). We fill in only up to body's 
    /// highest inboard coordinate.
    void calcKinematicJacobian(const State&, int body, const Vec3& station, Vector<Vec3>&) const;

    // XXX calcDynamicJacobian: first term of partial(inv(M)*f)/partial(q)
    // XXX need partials of Coriolis terms too

    // Stage Move
    const Real& getKineticEnergy(const State&) const;
    
    /// These are forces produced by gyroscopic effects and Coriolis accelerations.
    void addInBiasForces(const State&, Vector<SpatialVector>& bodyForces) const;

private:
    // Nodes of the tree stored by level. There is only one
    // Node at level 0 and that is the root of this tree.
    std::vector< std::vector<MultibodyTreeNode> > level;
};

// XXRigidBodyTree inlines

inline Vec3  XXRigidBodyTree::xformVector2Ground(const State& s, int body, const Vec3& vB) const
  { return getBodyFrameInGround(s,body).xformVector2Ref(vB); }
inline Vec3  XXRigidBodyTree::xformStation2Ground(const State& s, int body, const Vec3& sB) const
  { return getBodyFrameInGround(s,body).xformStation2Ref(sB); }
inline TransformMat XXRigidBodyTree::xformFrame2Ground(const State& s, int body, const TransformMat& fB) const
  { return getBodyFrameInGround(s,body).xformFrame2Ref(fB); }

inline Vec3  XXRigidBodyTree::xformVector2Parent(const State& s, int body, const Vec3& vB) const
  { return getBodyFrameInParent(s,body).xformVector2Ref(vB); }
inline Vec3  XXRigidBodyTree::xformStation2Parent(const State& s, int body, const Vec3& sB) const
  { return getBodyFrameInParent(s,body).xformStation2Ref(sB); }
inline TransformMat XXRigidBodyTree::xformFrame2Parent(const State& s, int body, const TransformMat& fB) const
  { return getBodyFrameInParent(s,body).xformFrame2Ref(fB); }

#endif /* SIMTK_RIGIDBODY_TREE_H_ */
