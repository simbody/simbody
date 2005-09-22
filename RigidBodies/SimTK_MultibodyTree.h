#ifndef SIMTK_MULTIBODY_TREE_H_
#define SIMTK_MULTIBODY_TREE_H_

#include <vector>

class Vec3;
class Mat33;
class SymMat33;
class Body;
class Joint;
class Frame;
class String;

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

#endif /* SIMTK_MULTIBODY_TREE_H_ */
