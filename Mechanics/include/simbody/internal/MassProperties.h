#ifndef SIMTK_SIMBODY_MASS_PROPERTIES_H_
#define SIMTK_SIMBODY_MASS_PROPERTIES_H_

/** @file
 *
 * These are utility classes for dealing with mass properties, particularly
 * those messy inertias.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/Geometry.h"

#include <iostream>

namespace simtk {

// Spatial vectors are used for (orientation,translation) quantities.
// These include
//      spatial velocity     = (angularVelocity,linearVelocity)
//      spatial acceleration = (angularAcceleration,linearAcceleration)
//      generalized forces   = (torque,force)
// Spatial configuration has to be handled differently though since
// orientation is not a vector quantity. (We use "TransformMat" for this concept
// which includes an orientation matrix and a translation vector.)
typedef Vec<2,   Vec3>  SpatialVec;
typedef Row<2,   Row3>  SpatialRow;
typedef Mat<2,2, Mat33> SpatialMat;

/**
 * The physical meaning of an inertia is the distribution of
 * a rigid body's mass about a *particular* point. If that point is the
 * center of mass of the body, then the measured inertia is called
 * the "central inertia" of that body. To write down the inertia, we
 * need to calculate the six scalars of the inertia tensor, which
 * is a symmetric 3x3 matrix. These scalars must be expressed in 
 * an arbitrary but specified coordinate system. So a InertiaMat
 * is meaningful only in conjunction with a particular set of axes, fixed
 * to the body, whose origin is the point about which the inertia is being
 * measured, and in whose coordinate system this measurement is being
 * expressed. Note that changing the reference point results in a new
 * physical quantity, but changing the reference axes only affects the
 * measure numbers of that quantity. For any reference point, there
 * is a unique set of reference axes in which the inertia tensor is
 * diagonal; those are called the "principal axes" of the body at that
 * point, and the resulting diagonal elements are the "principal moments
 * of inertia". When we speak of an inertia being "in" a frame, we mean
 * the physical quantity measured about the frame's origin and then expressed
 * in the frame's axes.
 *
 * This low-level InertiaMat class does not attempt to keep track of *which*
 * frame it is in. It concentrates instead on providing construction and
 * operations involving inertia which can proceed using only an implicit
 * frame F. Clients of this class are responsible for keeping track of that frame.
 * In particular, in order to shift the inertia's "measured-about" point one
 * must know whether either the starting or final inertia is central,
 * because we must always shift inertias by passing through the central inertia.
 * So this class provides operations for doing the shifting, but expects
 * to be told by the client where to find the center of mass.
 *
 * Re-expressing an InertiaMat in a different coordinate system does not
 * entail a change of physical meaning in the way that shifting it to a
 * different point does. Note that because inertia is a tensor, there is
 * a "left frame" and "right frame". For our purposes, these
 * will always be the same so we'll only indicate the frame 
 * once, as in 'I_pt_frame'. This should be understood to mean
 * 'frame_I_pt_frame' and re-expressing an InertiaMat requires both a
 * left and right multiply by the rotation matrix. So I_OB_B is the
 * inertia about body B's origin point OB, expressed in B, while
 * I_OB_G is the same physical quantity but expressed in Ground (the latter
 * is a component of the Spatial InertiaMat). Conversion is done like this:
 *    I_OB_G = R_GB * I_OB_B * R_BG  (and recall that R_GB=~R_BG)
 * The central inertia would be I_CB_B for body B.
 *
 * A InertiaMat is a symmetric matrix and is positive definite for
 * nonsingular bodies (that is, a body composed of at least three
 * noncollinear point masses).
 */
class SIMTK_SIMBODY_API InertiaMat {
public:
    /// Default is a NaN-ed out mess to avoid accidents.
    InertiaMat() : I_OF_F(NTraits<Real>::getNaN()) {}

    /// Create an principal inertia matrix with identical diagonal elements.
    /// Most commonly we create InertiaMat(0) for initialization of an inertia
    /// calculation (that is also the inertia of a point mass located at
    /// the origin). This can also be used with a non-zero value as the
    /// inertia of a sphere centered at the origin.
    explicit InertiaMat(const Real& r) : I_OF_F(r) { }

    /// Create a principal inertia matrix (only non-zero on diagonal).
    InertiaMat(const Real& xx, const Real& yy, const Real& zz)
        { setInertia(xx,yy,zz,0.,0.,0.); }

    /// This is a general inertia matrix. Note the order of these
    /// arguments: moments of inertia first, then products of inertia.
    InertiaMat(const Real& xx, const Real& yy, const Real& zz,
               const Real& xy, const Real& xz, const Real& yz)
        { setInertia(xx,yy,zz,xy,xz,yz); }

    /// Given a point mass located at a given point p in some frame F, 
    /// construct I_OF_F, that is, the inertia of that point mass about
    /// F's origin, expressed in F (that is, in F's coordinate system). 
    ///
    /// For a collection of point masses, you can just add these together to
    /// produce a composite inertia as long as all the vectors are
    /// measured from the same point and expressed in the same frame.
    InertiaMat(const Vec3& p, const Real& m) {
        Mat33& t = I_OF_F;
        const Real& x = p[0]; const Real xx = x*x;
        const Real& y = p[1]; const Real yy = y*y;
        const Real& z = p[2]; const Real zz = z*z;

        t(0,0)          =  m*(yy + zz);
        t(1,1)          =  m*(xx + zz);
        t(2,2)          =  m*(xx + yy);
        t(0,1) = t(1,0) = -m*x*y;
        t(0,2) = t(2,0) = -m*x*z;
        t(1,2) = t(2,1) = -m*y*z;
    }

    /// We only look at the lower triangles, but fill in the whole matrix.
    InertiaMat(const InertiaMat& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
    }

    explicit InertiaMat(const Mat33& s) {
        assert(close(s(0,1),s(1,0)) 
            && close(s(0,2),s(2,0))
            && close(s(1,2),s(2,1)));
        setInertia(s(0,0),s(1,1),s(2,2),
            0.5*(s(1,0)+s(0,1)),0.5*(s(2,0)+s(0,2)),0.5*(s(2,1)+s(1,2)));
    }

    InertiaMat& operator=(const InertiaMat& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
        return *this;
    }

    InertiaMat& operator+=(const InertiaMat& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) += s(0,0); t(1,1) += s(1,1);  t(2,2) += s(2,2);
        t(0,1) = (t(1,0) += s(1,0));
        t(0,2) = (t(2,0) += s(2,0));
        t(1,2) = (t(2,1) += s(2,1));
        return *this;
    }

    InertiaMat& operator-=(const InertiaMat& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) -= s(0,0); t(1,1) -= s(1,1);  t(2,2) -= s(2,2);
        t(0,1) = (t(1,0) -= s(1,0));
        t(0,2) = (t(2,0) -= s(2,0));
        t(1,2) = (t(2,1) -= s(2,1));
        return *this;
    }

    InertiaMat& operator*=(const Real& r) {
        I_OF_F *= r;
        return *this;
    }
    InertiaMat& operator/=(const Real& r) {
        I_OF_F /= r;
        return *this;
    }
    void setInertia(const Real& xx, const Real& yy, const Real& zz) {
        I_OF_F = 0.; I_OF_F(0,0) = xx; I_OF_F(1,1) = yy;  I_OF_F(2,2) = zz;
    }


    void setInertia(const Real& xx, const Real& yy, const Real& zz,
                    const Real& xy, const Real& xz, const Real& yz) {
        Mat33& t = I_OF_F;
        t(0,0) = xx; t(1,1) = yy;  t(2,2) = zz;
        t(0,1) = t(1,0) = xy;
        t(0,2) = t(2,0) = xz;
        t(1,2) = t(2,1) = yz;
    }

    /// Assume that the current inertia is about the F frame's origin OF, and
    /// expressed in F. Given the vector from OF to the body center of mass CF,
    /// and the total mass of the body, we can shift the inertia to the center
    /// of mass. This produces a new InertiaMat I' whose (implicit) frame F' is
    /// aligned with F but has origin CF (an inertia like that is called a "central
    /// inertia". I' = I - Icom where Icom is the inertia of a fictitious
    /// point mass of mass mtot located at CF (measured in F) about OF.
    inline InertiaMat shiftToCOM(const Vec3& CF, const Real& mtot) const;

    /// Assuming that the current inertia I is a central inertia (that is, it is
    /// inertia about the body center of mass CF), shift it to some other point p
    /// measured from the center of mass. This produces a new inertia I' about
    /// the point p given by I' = I + Ip where Ip is the inertia of a fictitious
    /// point mass of mass mtot (the total body mass) located at p, about CF.
    inline InertiaMat shiftFromCOM(const Vec3& p, const Real& mtot) const;

    /// Re-express this inertia from frame F to frame B, given the orientation
    /// of B in F. This is a similarity transform since rotation matrices are
    /// orthogonal.
    InertiaMat changeAxes(const RotationMat& R_FB) const {
        return InertiaMat(~R_FB * I_OF_F * R_FB); // TODO can do better due to symmetry
    }

    Real trace() const {return I_OF_F(0,0) + I_OF_F(1,1) + I_OF_F(2,2);}

    // Note that we are copying into the Mat33 here so this will still
    // work if we decide to switch to a SymMat33 when they exist.
    Mat33 toMat33() const {
        return I_OF_F;
    }

    // InertiaMat factory for some common mass elements. Each defines its own
    // frame aligned (when possible) with principal moments. Each has unit
    // mass and its center of mass located at the origin. Use this with shiftFromCOM()
    // to move it somewhere else, and with xform() to express the inertia
    // in another frame. Mass enters linearly in inertia, so just multiply
    // by the actual mass to scale these properly.
    static InertiaMat point() {return InertiaMat(0.);}
    static InertiaMat sphere(const Real& r) {return InertiaMat(0.4*r*r);}

    // Cylinder is aligned along z axis, use radius and half-length.
    // If r==0 this is a thin rod; hz=0 it is a thin disk.
    static InertiaMat cylinder(const Real& r, const Real& hz) {
        const Real Ixx = 0.25*r*r + (Real(1)/Real(3))*hz*hz;
        return InertiaMat(Ixx,Ixx,0.5*r*r);
    }

    // Brick given by half-lengths in each direction. One dimension zero
    // gives inertia of a thin rectangular sheet; two zero gives inertia
    // of a thin rod in the remaining direction.
    static InertiaMat brick(const Real& hx, const Real& hy, const Real& hz) {
        const Real oo3 = Real(1)/Real(3);
        const Real hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return InertiaMat(oo3*(hy2+hz2), oo3*(hx2+hz2), oo3*(hx2+hy2));
    }

    // Ellipsoid given by half-lengths in each direction.
    static InertiaMat ellipsoid(const Real& hx, const Real& hy, const Real& hz) {
        const Real hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return InertiaMat(0.2*(hy2+hz2), 0.2*(hx2+hz2), 0.2*(hx2+hy2));
    }


private:
    //TODO: the tolerance here should be a function of Real's precision
    static bool close(const Real& a, const Real& b) {
        if (fabs(a-b) < 1e-13) return true;
        if (fabs(a-b)/(fabs(a)+fabs(b)) < 0.5e-13) return true;
        return false;
    }

private:
    // InertiaMat expressed in frame F and about F's origin OF. Note that frame F
    // is implicit here; all we actually have are the inertia scalars. This is 
    // a symmetric matrix but we keep all the elements here, and manage them
    // so that the reflected elements are *exactly* equal.
    // TODO: should use a SymMat33 type.
    Mat33 I_OF_F;
    friend Vec3 operator*(const InertiaMat& i, const Vec3& w);
};

inline InertiaMat operator+(const InertiaMat& l, const InertiaMat& r) {
    return InertiaMat(l) += r;
}
inline InertiaMat operator-(const InertiaMat& l, const InertiaMat& r) {
    return InertiaMat(l) -= r;
}
inline InertiaMat operator*(const InertiaMat& i, const Real& r) {
    return InertiaMat(i) *= r;
}
inline InertiaMat operator*(const Real& r, const InertiaMat& i) {
    return InertiaMat(i) *= r;
}
inline Vec3 operator*(const InertiaMat& i, const Vec3& w) {
    return i.I_OF_F * w;
}
inline InertiaMat operator/(const InertiaMat& i, const Real& r) {
    return InertiaMat(i) /= r;
}
inline InertiaMat InertiaMat::shiftToCOM(const Vec3& CF, const Real& mtot) const {
    return *this - InertiaMat(CF, mtot);
}
inline InertiaMat InertiaMat::shiftFromCOM(const Vec3& p, const Real& mtot) const {
    return *this + InertiaMat(p, mtot);
}

inline bool
operator==(const InertiaMat& i1, const InertiaMat& i2) {
    return i1.toMat33() == i2.toMat33();    // TODO should use underlying rep
}

SIMTK_SIMBODY_API std::ostream& 
operator<<(std::ostream& o, const InertiaMat&);


/**
 * This class contains the mass, centroid, and inertia of a rigid body B.
 * The centroid is a vector from B's origin, expressed in the B frame.
 * The inertia is taken about the B origin, and expressed in B.
 */
class SIMTK_SIMBODY_API MassProperties {
public:
    MassProperties() { setMassProperties(0.,Vec3(0.),InertiaMat()); }
    MassProperties(const Real& m, const Vec3& com, const InertiaMat& inertia)
      { setMassProperties(m,com,inertia); }

    void setMassProperties(const Real& m, const Vec3& com, const InertiaMat& inertia)
      { mass=m; comInB=com; inertia_OB_B=inertia; }

    const Real&       getMass()    const { return mass; }
    const Vec3&       getCOM()     const { return comInB; }
    const InertiaMat& getInertia() const { return inertia_OB_B; }

    InertiaMat calcCentroidalInertia() const {
        return inertia_OB_B - InertiaMat(comInB, mass);
    }
    InertiaMat calcShiftedInertia(const Vec3& newOriginB) const {
        return calcCentroidalInertia() + InertiaMat(newOriginB-comInB, mass);
    }
    MassProperties calcShiftedMassProps(const Vec3& newOriginB) const {
        return MassProperties(mass, comInB-newOriginB,
                              calcShiftedInertia(newOriginB));
    }

    bool isMassless()   const { return mass==0.; }

private:
    Real        mass;
    Vec3        comInB;         // meas. from B origin, expr. in B
    InertiaMat  inertia_OB_B;   // about B origin, expr. in B
};

} // namespace simtk

#endif /* SIMTK_SIMBODY_MASS_PROPERTIES_H_ */
