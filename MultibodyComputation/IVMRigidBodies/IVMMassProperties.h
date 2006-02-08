#ifndef IVM_MASS_PROPERTIES_H_
#define IVM_MASS_PROPERTIES_H_

/** @file
 *
 * These are utility classes for dealing with mass properties, particularly
 * those messy inertias.
 */

#include "cdsVec3.h"
#include "cdsMat33.h"
#include "fixedVector.h"
#include "cdsList.h"

typedef float_type CDSReal;


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
enum IVMJointType {
    IVMUnknownJointType    = 0,
    IVMThisIsGround        = 1, // Ground's "inboard joint"
    IVMWeldJoint           = 2,
    IVMTorsionJoint        = 3,
    IVMSlidingJoint        = 4,
    IVMUJoint              = 5,
    IVMCylinderJoint       = 6,
    IVMPlanarJoint         = 7,
    IVMGimbalJoint         = 8,
    IVMOrientationJoint    = 9,
    IVMCartesianJoint      = 10,
    IVMFreeLineJoint       = 11,
    IVMFreeJoint           = 12
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
class IVMInertia : public CDSMat33 {
public:
    /// Default is the inertia of a point mass about that point, i.e. 0.
    IVMInertia() : CDSMat33(0.) {}

    /// Note the order of these arguments: moments of inertia first, then 
    /// products of inertia.
    IVMInertia(const CDSReal& xx, const CDSReal& yy, const CDSReal& zz,
            const CDSReal& xy, const CDSReal& xz, const CDSReal& yz)
        { setInertia(xx,yy,zz,xy,xz,yz); }

    /// Given a point mass located at a given point p in some frame F, 
    /// construct I_OF_F, that is, the inertia of that point mass about
    /// F's origin, expressed in F. 
    ///
    /// For a collection of point masses, you can just add these together to
    /// produce a composite inertia as long as all the vectors are
    /// measured from the same point and expressed in the same frame.
    IVMInertia(const CDSReal& m, const CDSVec3& p) {
        CDSMat33& t = *this;
        const CDSReal& x = p(0); const CDSReal xx = x*x;
        const CDSReal& y = p(1); const CDSReal yy = y*y;
        const CDSReal& z = p(2); const CDSReal zz = z*z;

        t(0,0)          =  m*(yy + zz);
        t(1,1)          =  m*(xx + zz);
        t(2,2)          =  m*(xx + yy);
        t(0,1) = t(1,0) = -m*x*y;
        t(0,2) = t(2,0) = -m*x*z;
        t(1,2) = t(2,1) = -m*y*z;
    }

    /// We only look at the lower triangles, but fill in the whole matrix.
    IVMInertia(const IVMInertia& s) {
        CDSMat33& t = *this;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
    }

    explicit IVMInertia(const CDSMat33& s) {
        assert(close(s(0,1),s(1,0)) 
            && close(s(0,2),s(2,0))
            && close(s(1,2),s(2,1)));
        setInertia(s(0,0),s(1,1),s(2,2),
            0.5*(s(1,0)+s(0,1)),0.5*(s(2,0)+s(0,2)),0.5*(s(2,1)+s(1,2)));
    }

    IVMInertia& operator=(const IVMInertia& s) {
        CDSMat33& t = *this;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
        return *this;
    }

    IVMInertia& operator+=(const IVMInertia& s) {
        CDSMat33& t = *this;
        t(0,0) += s(0,0); t(1,1) += s(1,1);  t(2,2) += s(2,2);
        t(0,1) = (t(1,0) += s(1,0));
        t(0,2) = (t(2,0) += s(2,0));
        t(1,2) = (t(2,1) += s(2,1));
        return *this;
    }

    IVMInertia& operator-=(const IVMInertia& s) {
        CDSMat33& t = *this;
        t(0,0) -= s(0,0); t(1,1) -= s(1,1);  t(2,2) -= s(2,2);
        t(0,1) = (t(1,0) -= s(1,0));
        t(0,2) = (t(2,0) -= s(2,0));
        t(1,2) = (t(2,1) -= s(2,1));
        return *this;
    }

    void setInertia(const CDSReal& xx, const CDSReal& yy, const CDSReal& zz,
                    const CDSReal& xy, const CDSReal& xz, const CDSReal& yz) {
        CDSMat33& t = *this;
        t(0,0) = xx; t(1,1) = yy;  t(2,2) = zz;
        t(0,1) = t(1,0) = xy;
        t(0,2) = t(2,0) = xz;
        t(1,2) = t(2,1) = yz;
    }

private:
    //TODO: the tolerance here should be a function of CDSReal's precision
    static bool close(const CDSReal& a, const CDSReal& b) {
        if (fabs(a-b) < 1e-13) return true;
        if (fabs(a-b)/(fabs(a)+fabs(b)) < 0.5e-13) return true;
        return false;
    }
};

inline IVMInertia operator+(const IVMInertia& l, const IVMInertia& r) {
    return IVMInertia(l) += r;
}
inline IVMInertia operator-(const IVMInertia& l, const IVMInertia& r) {
    return IVMInertia(l) -= r;
}

/**
 * This class contains the mass, centroid, and inertia of a rigid body B.
 * The centroid is a vector from B's origin, expressed in the B frame.
 * The inertia is taken about the B origin, and expressed in B.
 */
class IVMMassProperties {
public:
    IVMMassProperties() { setMassProperties(0.,CDSVec3(0.),IVMInertia()); }
    IVMMassProperties(const CDSReal& m, const CDSVec3& com, const IVMInertia& inertia)
      { setMassProperties(m,com,inertia); }

    void setMassProperties(const CDSReal& m, const CDSVec3& com, const IVMInertia& inertia)
      { mass=m; comInB=com; inertia_OB_B=inertia; }

    const CDSReal&    getMass()    const { return mass; }
    const CDSVec3&    getCOM()     const { return comInB; }
    const IVMInertia&  getInertia() const { return inertia_OB_B; }

    IVMInertia calcCentroidalInertia() const {
        return inertia_OB_B - IVMInertia(mass, comInB);
    }
    IVMInertia calcShiftedInertia(const CDSVec3& newOriginB) const {
        return calcCentroidalInertia() + IVMInertia(mass, newOriginB-comInB);
    }
    IVMMassProperties calcShiftedMassProps(const CDSVec3& newOriginB) const {
        return IVMMassProperties(mass, comInB-newOriginB,
                                calcShiftedInertia(newOriginB));
    }

    bool isMassless()   const { return mass==0.; }

private:
    CDSReal   mass;
    CDSVec3   comInB;         // meas. from B origin, expr. in B
    IVMInertia inertia_OB_B;   // about B origin, expr. in B
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
 * to reinterpret IVMFrame objects in any appropriate manner that depends
 * on this memory layout.
 */
class IVMFrame {
public:
    IVMFrame() {Rot_RF.set(0.); Rot_RF.setDiag(1.); Loc_RF.set(0.);}
    IVMFrame(const CDSMat33& axesInR, const CDSVec3& originInR) {setFrame(axesInR,originInR);}
    IVMFrame(const CDSVec3& originInR) {
        Rot_RF.set(0.); Rot_RF.setDiag(1.); Loc_RF=originInR;
    }
    // default copy, assignment, destructor

    void setFrame(const CDSMat33& axesInR, const CDSVec3& originInR) 
        {Rot_RF=axesInR; Loc_RF=originInR;}

    // Transform various items measured and expressed in F to those same
    // items measured and expressed in F's reference frame R.
    CDSVec3  xformVector2Ref  (const CDSVec3& vF)      const { return Rot_RF*vF; }
    CDSVec3  xformStation2Ref (const CDSVec3& sF)      const { return Loc_RF + xformVector2Ref(sF); }
    CDSMat33 xformRotation2Ref(const CDSMat33& Rot_FX) const { return Rot_RF*Rot_FX; }
    IVMFrame  xformFrame2Ref   (const IVMFrame& fF) const 
      { return IVMFrame(xformRotation2Ref(fF.Rot_RF), xformStation2Ref(fF.Loc_RF)); }

    const CDSMat33& getRot_RF() const { return Rot_RF; }
    CDSMat33&       updRot_RF()       { return Rot_RF; }

    const CDSVec3&  getLoc_RF() const { return Loc_RF; }
    CDSVec3&        updLoc_RF()       { return Loc_RF; }

    // Computation-free conversions
    const CDSReal*  getFrameAsArray(const IVMFrame& f)  const {return reinterpret_cast<const CDSReal*>(&f);}
    CDSReal*        updFrameAsArray(IVMFrame& f)              {return reinterpret_cast<CDSReal*>(&f);}
    const IVMFrame&  getArrayAsFrame(const CDSReal* r)  const {return *reinterpret_cast<const IVMFrame*>(r);}
    IVMFrame&        updArrayAsFrame(CDSReal* r)              {return *reinterpret_cast<IVMFrame*>(r);}

private:
    CDSMat33 Rot_RF;   // rotation matrix that expresses F's axes in R
    CDSVec3  Loc_RF;   // location of F's origin measured from R's origin, expressed in R 
};

#endif /* RB_MASS_PROPERTIES_H_ */
