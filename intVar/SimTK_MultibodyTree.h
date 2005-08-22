#ifndef SIMTK_MULTIBODY_TREE_H_
#define SIMTK_MULTIBODY_TREE_H_

#include <vector>

class Vec3;
class Body;
class Joint;
class Frame;
class Station;
class RigidBody;

/**
 * This class represents an orthogonal, right-handed coordinate frame F, 
 * measured from and expressed in a reference coordinate frame R. F consists of
 * 3 perpendicular unit vectors defining its axes as viewed from R, 
 * and a vector from R's origin point OR to F's origin point OF. Note that
 * the meaning of "R" comes from the context in which frame F is used.
 * We use the phrase "frame F is in frame R" to describe the above relationship,
 * that is, "in" means both measured in and expressed in. 
 *
 * The axis vectors are ordered 1-2-3 or x-y-z as you prefer, with
 * z = x X y, making a right-handed set. These axes are arranged as
 * columns of a 3x3 matrix Rot_RF = [ x y z ] which is a direction cosine
 * (rotation) matrix useful for conversions between frame R and F. For
 * example, given a vector vF expressed in the F frame, that same vector
 * re-expressed in R is given by vR = Rot_RF*vF. The origin point Loc_RF is 
 * stored as the vector Loc_RF=(OF-OR) and expressed in R.
 *
 * This is a "POD" (plain old data) class with a well-defined memory
 * layout on which a client of this class may depend: There are 
 * exactly 4 consecutive, packed 3-vectors in the order x,y,z,O.
 * That is, this class is equivalent to an array of 12 Reals with 
 * the order x1,x2,x3,y1,y2,y3,z1,z2,z3,O1,O2,O3. It is expressly allowed
 * to reinterpret Frame objects in any appropriate manner that depends
 * on this memory layout.
 */
class Frame {
public:
    Frame() { }
    Frame(const Mat33& axes, const Vec3& origin) {setFrame(axes,origin);}
    // default copy, assignment, destructor

    setFrame(const Mat33& axes, const Vec3& origin) {Rot_RF=axes; Loc_RF=origin;}

    // Transform various items measured and expressed in F to those same
    // items measured and expressed in F's reference frame R.
    Vec3  xformVector2Ref  (const Vec3& vF)      const { return Rot_RF*vF; }
    Vec3  xformStation2Ref (const Vec3& sF)      const { return Loc_RF + xformVector2Ref(sF); }
    Mat33 xformRotation2Ref(const Mat33& Rot_FX) const { return Rot_RF*Rot_FX; }
    Frame xformFrame2Ref(const Frame& fF) const 
      { return Frame(xformRotation2Ref(fF.Rot_RF), xformStation2Ref(fF.Loc_RF)); }

    const Mat33& getRot_RF() const { return Rot_RF; }
    Mat33&       updRot_RF()       { return Rot_RF; }

    const Vec3&  getLoc_RF() const { return Loc_RF; }
    Vec3&        updLoc_RF()       { return Loc_RF; }

    // Computation-free conversions
    const Real*  getFrameAsArray(const Frame& f) const {return reinterpret_cast<const Real*>(&f);}
    Real*        updFrameAsArray(Frame& f)             {return reinterpret_cast<Real*>(&f);}
    const Frame& getArrayAsFrame(const Real* r)  const {return *reinterpret_cast<const Frame*>(r);}
    Frame&       updArrayAsFrame(Real* r)              {return *reinterpret_cast<Frame*>(r);}

private:
    Mat33 Rot_RF;   // rotation matrix that expresses F's axes in R
    Vec3  Loc_RF;   // location of F's origin measured from R's origin, expressed in R 
};

/**
 * RigidBodyTree is a multibody system restricted to unconstrained tree topology
 * and rigid bodies. That is, there are no loops or constraints, and
 * it assumes that body stations, frames, and mass properties are constants
 * when expressed in the body frame.
 *
 * This class contains the system topology, that is, Ground, the bodies and
 * their inboard joints. Ground is body #0, other bodies have a regular labeling
 * such that body numbers are strictly decreasing as you go from a body to its
 * parent, parent's parent, and so on down to Ground. A body and its inboard
 * joint have the same number.
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
 */
class RigidBodyTree { 
public:
    RigidBodyTree();
    RigidBodyTree(const RigidBodyTree&);
    RigidBodyTree& operator=(const RigidBodyTree&);
    ~RigidBodyTree();

    /// Add a new body connected to one of the existing bodies by
    /// a Joint. Body 0 (Ground) is predefined at construction.
    int addBody(const RigidBody&, const Frame& frameInB, 
                int   parent,     const Frame& frameInP,
                const Joint&);

    /// Return number of defined bodies, including ground. This is
    /// one more than the number of joints.
    int getNBodies() const;

    /// Evaluate the quantities appropriate for the indicated stage,
    /// allowing access to responses, operators and solvers for that stage.
    void realize(const State&, Stage, const RealizeOptions&) const;

    // Stage Configure //

    /// Obtain the spatial location and orientation of the body frame, by
    /// returning the frame in Ground which is currently overlaid on the body frame.
    const Frame& getBodyFrameInGround(const State&, int body) const;

    /// Obtain the location and orientation of the body frame in its parent frame, by
    /// returning the frame in Parent which is currently overlaid on the body frame.
    const Frame& getBodyFrameInParent(const State&, int body) const;

    Vec3 xformVector2Ground(const State& s, int body, const Vec3& vB) const
      { return getBodyFrameInGround(s,body).xformVector2Ref(vB); }
    Vec3 xformStation2Ground(const State& s, int body, const Vec3& sB) const
      { return getBodyFrameInGround(s,body).xformStation2Ref(sB); }
    Frame xformFrame2Ground(const State& s, int body, const Frame& fB) const
      { return getBodyFrameInGround(s,body).xformFrame2Ref(fB); }

    Vec3 xformVector2Parent(const State& s, int body, const Vec3& vB) const
      { return getBodyFrameInParent(s,body).xformVector2Ref(vB); }
    Vec3 xformStation2Parent(const State& s, int body, const Vec3& sB) const
      { return getBodyFrameInParent(s,body).xformStation2Ref(sB); }
    Frame xformFrame2Parent(const State& s, int body, const Frame& fB) const
      { return getBodyFrameInParent(s,body).xformFrame2Ref(fB); }

    /// XXX does this interface work for deformable bodies?

    /// Return total mass of tree rooted at this body.
    const Real& getCompositeMass(int body) const;

    /// Return center of mass of tree rooted at this body, measured from body
    /// origin and expressed in body frame.
    const Vec3& getCompositeCOM(int body) const;

    /// Return inertia of tree rooted at body, as though all joints were welded.
    /// Inertia is about body origin and expressed in body frame.
    const Mat33& getCompositeInertia(int body) const;

    /// Return inertia of tree rooted at body, with all joint floppy and unactuated.
    /// Inertia is about body origin and expressed in body frame.
    const Mat33& getArticulatedInertia(int body) const;


private:
    // Nodes of the tree stored by level. There is only one
    // Node at level 0 and that is the root of this tree.
    std::vector< std::vector<MultibodyTreeNode> > level;
};

/**
 * Primitive mass elements useful for constructing bodies. These have 
 * implied reference frames depending on the kind of element, always
 * located at the center of mass. There is *no* geometry implied by
 * these objects, just pure mass.
 */
class MassElement {
private:
    Real mass;
};

/**
 * A PointMass element has only a scalar mass, with center of
 * mass 0,0,0 and an all-zero inertia matrix. Attach this to a point.
 */
class PointMass : public MassElement {
public:
private:
};

/**
 * A LineMass has mass distributed along its z axis, resulting
 * in an identical non-zero centroidal inertia about its x and y axes,
 * and zero inertia around z. Attach its center of mass to a point
 * and provide a direction for its mass distribution.
 */
class LineMass : public MassElement {
public:
private:
    Real inertia;   // about x and y
};

/**
 * A RigidBodyMass provides a full centroidal inertia. Attach this 
 * to a frame to align the inertia, and the center of mass will
 * be placed at the origin of the frame.
 */
class RigidBodyMass : public MassElement {
public:
private:
    Mat33 inertia;   // centroidal inertia
};

/**
 * A Body is fundamentally a reference frame. With this frame we can
 * express directions and measure distances from the frame origin 
 * along those directions. That can be used to define stations (point
 * locations) on the body and the
 * orientations of other reference frames attached to the same Body.
 * We can attach and orient mass elements to the resulting frames to
 * determine the Body's mass properties. Later, we can use those frames
 * to specify joints and geometric features.
 */
class Body {
public:
    Body() { }

    // construction
    void setNDirections(int);
    void setNStations();

    // parametrization
    void setMass(const Real&);
    void setCenterOfMass(const Vec3&);
    void setInertia(const Mat33&);
    void setFrame(int, const Mat33&, const Vec3&);


private:
    // Construction
    String name;

    // Frames are stored in consecutive locations in the nominalVectorPool like
    // this:  [x y z O]. Here we map frameId to the index of its x vector.
    std::vector<int>   nominalFrames;

    // State variables

    // Parametrize Stage
    Real        mass;
    Vec3        nominalCOMstation;
    SymMat33    nominalInertia;        // about origin

    // A bunch of body-fixed vectors that we always care about. These are
    // fixed after parameter stage and represent the actual vectors for
    // rigid bodies or the nominal values for deformable bodies.
    std::vector<Vec3>  nominalVectorPool;

};

class Joint {
public:
private:
    String type;
    int fixedFrame;    // in a tree, aka 'parent frame'
    int movingFrame;   //       "        'child frame'
};

class MultibodyTreeNode {
public:
private:
    Body  b;
    Joint inboardJoint;

    const MultibodyTreeNode* parent;
    std::vector<const MultibodyTreeNode*> children;
};


inline Vec3
RigidBodyTree::calcVectorInGround(const State& s, int body, const Vec3& vB) const {
    const Mat33& Rot_GB = getBodyFrameInGround(s,body).getRot_RF();
    return Rot_GB*vB;
}
inline Vec3
RigidBodyTree::calcStationInGround(const State& s, int body, const Vec3& posB) const {
    const Vec3& Loc_GB = getBodyFrameInGround(s,body).getLoc_RF();
    return Loc_GB + calcVectorInGround(s,body,posB));
}
inline Vec3
RigidBodyTree::calcVectorInParent(const State& s, int body, const Vec3& vB) const {
    const Mat33& Rot_PB = getBodyFrameInParent(s,body).getRot_RF();
    return Rot_PB*vB;
}
inline Vec3
RigidBodyTree::calcStationInParent(const State& s, int body, const Vec3& posB) const {
    const Vec3& Loc_PB = getBodyFrameInParent(s,body).getLoc_RF();
    return Loc_PB + calcVectorInParent(s,body,posB));
}
#endif /* SIMTK_MULTIBODY_TREE_H_ */
