#ifndef RB_MASS_PROPERTIES_H_
#define RB_MASS_PROPERTIES_H_

/** @file
 *
 * These are utility classes for dealing with mass properties, particularly
 * those messy inertias.
 */

#include "cdsVec3.h"
#include "cdsMat33.h"
#include "fixedVector.h"
#include "cdsList.h"


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
enum RBJointType {
    RBUnknownJointType    = 0,
    RBThisIsGround        = 1, // Ground's "inboard joint"
    RBWeldJoint           = 2,
    RBTorsionJoint        = 3,
    RBSlidingJoint        = 4,
    RBUJoint              = 5,
    RBCylinderJoint       = 6,
    RBPlanarJoint         = 7,
    RBGimbalJoint         = 8,
    RBOrientationJoint    = 9,
    RBCartesianJoint      = 10,
    RBFreeLineJoint       = 11,
    RBFreeJoint           = 12
};



/**
 * The physical meaning of an Inertia tensor is the distribution of
 * a rigid body's mass about a *particular* point. The measure numbers
 * of the Inertia must be expressed in a particular frame. We
 * use the notation I_pt_frame to keep out of trouble. So I_OB_B is
 * the inertia about body B's origin point OB, expressed in B, while
 * I_OB_G is the same physical quantity but expressed in Ground (the latter
 * is a component of the Spatial Inertia). 
 *
 * It is often useful to know the inertia about a body's center of mass
 * (called the "centroidal inertia"). This would be I_CB_B for body B.
 *
 * An Inertia is a symmetric matrix and is positive definite for
 * nonsingular bodies (that is, a body composed of at least three
 * noncollinear point masses).
 */
class RBInertia : public Mat33 {
public:
    /// Default is the inertia of a point mass about that point, i.e. 0.
    RBInertia() : Mat33(0.) {}

    /// Note the order of these arguments: moments of inertia first, then 
    /// products of inertia.
    RBInertia(const Real& xx, const Real& yy, const Real& zz,
            const Real& xy, const Real& xz, const Real& yz)
        { setInertia(xx,yy,zz,xy,xz,yz); }

    /// Given a point mass located at a given point p in some frame F, 
    /// construct I_OF_F, that is, the inertia of that point mass about
    /// F's origin, expressed in F. 
    ///
    /// For a collection of point masses, you can just add these together to
    /// produce a composite inertia as long as all the vectors are
    /// measured from the same point and expressed in the same frame.
    RBInertia(const Real& m, const Vec3& p) {
        Mat33& t = *this;
        const Real& x = p(0); const Real xx = x*x;
        const Real& y = p(1); const Real yy = y*y;
        const Real& z = p(2); const Real zz = z*z;

        t(0,0)          =  m*(yy + zz);
        t(1,1)          =  m*(xx + zz);
        t(2,2)          =  m*(xx + yy);
        t(0,1) = t(1,0) = -m*x*y;
        t(0,2) = t(2,0) = -m*x*z;
        t(1,2) = t(2,1) = -m*y*z;
    }

    /// We only look at the lower triangles, but fill in the whole matrix.
    RBInertia(const RBInertia& s) {
        Mat33& t = *this;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
    }

    explicit RBInertia(const Mat33& s) {
        assert(close(s(0,1),s(1,0)) 
            && close(s(0,2),s(2,0))
            && close(s(1,2),s(2,1)));
        setInertia(s(0,0),s(1,1),s(2,2),
            0.5*(s(1,0)+s(0,1)),0.5*(s(2,0)+s(0,2)),0.5*(s(2,1)+s(1,2)));
    }

    RBInertia& operator=(const RBInertia& s) {
        Mat33& t = *this;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
        return *this;
    }

    RBInertia& operator+=(const RBInertia& s) {
        Mat33& t = *this;
        t(0,0) += s(0,0); t(1,1) += s(1,1);  t(2,2) += s(2,2);
        t(0,1) = (t(1,0) += s(1,0));
        t(0,2) = (t(2,0) += s(2,0));
        t(1,2) = (t(2,1) += s(2,1));
        return *this;
    }

    RBInertia& operator-=(const RBInertia& s) {
        Mat33& t = *this;
        t(0,0) -= s(0,0); t(1,1) -= s(1,1);  t(2,2) -= s(2,2);
        t(0,1) = (t(1,0) -= s(1,0));
        t(0,2) = (t(2,0) -= s(2,0));
        t(1,2) = (t(2,1) -= s(2,1));
        return *this;
    }

    void setInertia(const Real& xx, const Real& yy, const Real& zz,
                    const Real& xy, const Real& xz, const Real& yz) {
        Mat33& t = *this;
        t(0,0) = xx; t(1,1) = yy;  t(2,2) = zz;
        t(0,1) = t(1,0) = xy;
        t(0,2) = t(2,0) = xz;
        t(1,2) = t(2,1) = yz;
    }

private:
    //TODO: the tolerance here should be a function of Real's precision
    static bool close(const Real& a, const Real& b) {
        if (fabs(a-b) < 1e-13) return true;
        if (fabs(a-b)/(fabs(a)+fabs(b)) < 0.5e-13) return true;
        return false;
    }
};

inline RBInertia operator+(const RBInertia& l, const RBInertia& r) {
    return RBInertia(l) += r;
}
inline RBInertia operator-(const RBInertia& l, const RBInertia& r) {
    return RBInertia(l) -= r;
}

/**
 * This class contains the mass, centroid, and inertia of a rigid body B.
 * The centroid is a vector from B's origin, expressed in the B frame.
 * The inertia is taken about the B origin, and expressed in B.
 */
class RBMassProperties {
public:
    RBMassProperties() { setMassProperties(0.,Vec3(0.),RBInertia()); }
    RBMassProperties(const Real& m, const Vec3& com, const RBInertia& inertia)
      { setMassProperties(m,com,inertia); }

    void setMassProperties(const Real& m, const Vec3& com, const RBInertia& inertia)
      { mass=m; comInB=com; inertia_OB_B=inertia; }

    const Real&    getMass()    const { return mass; }
    const Vec3&    getCOM()     const { return comInB; }
    const RBInertia&  getInertia() const { return inertia_OB_B; }

    RBInertia calcCentroidalInertia() const {
        return inertia_OB_B - RBInertia(mass, comInB);
    }
    RBInertia calcShiftedInertia(const Vec3& newOriginB) const {
        return calcCentroidalInertia() + RBInertia(mass, newOriginB-comInB);
    }
    RBMassProperties calcShiftedMassProps(const Vec3& newOriginB) const {
        return RBMassProperties(mass, comInB-newOriginB,
                                calcShiftedInertia(newOriginB));
    }

    bool isMassless()   const { return mass==0.; }

private:
    Real   mass;
    Vec3   comInB;         // meas. from B origin, expr. in B
    RBInertia inertia_OB_B;   // about B origin, expr. in B
};

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
class RBFrame {
public:
    RBFrame() : Rot_RF(1), Loc_RF(0) { } // identity rotation and 0 origin

    RBFrame(const Mat33& axesInR, const Vec3& originInR) {setFrame(axesInR,originInR);}
    RBFrame(const Vec3& originInR) {
        Rot_RF = 1; Loc_RF=originInR;
    }
    // default copy, assignment, destructor

    void setFrame(const Mat33& axesInR, const Vec3& originInR) 
        {Rot_RF=axesInR; Loc_RF=originInR;}

    // Transform various items measured and expressed in F to those same
    // items measured and expressed in F's reference frame R.
    Vec3  xformVector2Ref  (const Vec3& vF)      const { return Rot_RF*vF; }
    Vec3  xformStation2Ref (const Vec3& sF)      const { return Loc_RF + xformVector2Ref(sF); }
    Mat33 xformRotation2Ref(const Mat33& Rot_FX) const { return Rot_RF*Rot_FX; }
    RBFrame  xformFrame2Ref   (const RBFrame& fF) const 
      { return RBFrame(xformRotation2Ref(fF.Rot_RF), xformStation2Ref(fF.Loc_RF)); }

    const Mat33& getRot_RF() const { return Rot_RF; }
    Mat33&       updRot_RF()       { return Rot_RF; }

    const Vec3&  getLoc_RF() const { return Loc_RF; }
    Vec3&        updLoc_RF()       { return Loc_RF; }

    // Computation-free conversions
    const Real*  getFrameAsArray(const RBFrame& f)  const {return reinterpret_cast<const Real*>(&f);}
    Real*        updFrameAsArray(RBFrame& f)              {return reinterpret_cast<Real*>(&f);}
    const RBFrame&  getArrayAsFrame(const Real* r)  const {return *reinterpret_cast<const RBFrame*>(r);}
    RBFrame&        updArrayAsFrame(Real* r)              {return *reinterpret_cast<RBFrame*>(r);}

private:
    Mat33 Rot_RF;   // rotation matrix that expresses F's axes in R
    Vec3  Loc_RF;   // location of F's origin measured from R's origin, expressed in R 
};

#endif /* RB_MASS_PROPERTIES_H_ */
