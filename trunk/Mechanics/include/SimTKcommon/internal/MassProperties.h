#ifndef SimTK_SIMMATRIX_MASS_PROPERTIES_H_
#define SimTK_SIMMATRIX_MASS_PROPERTIES_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
 *
 * These are utility classes for dealing with mass properties, particularly
 * those messy inertias.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/Orientation.h"

#include <iostream>

namespace SimTK {
typedef Mat<2,2, Mat33> SpatialMat;

// Spatial vectors are used for (orientation,translation) quantities.
// These include
//      spatial velocity     = (angularVelocity,linearVelocity)
//      spatial acceleration = (angularAcceleration,linearAcceleration)
//      generalized forces   = (torque,force)
// Spatial configuration has to be handled differently though since
// orientation is not a vector quantity. (We use "Transform" for this concept
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
 * an arbitrary but specified coordinate system. So an Inertia
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
 * This low-level Inertia class does not attempt to keep track of *which*
 * frame it is in. It concentrates instead on providing construction and
 * operations involving inertia which can proceed using only an implicit
 * frame F. Clients of this class are responsible for keeping track of that frame.
 * In particular, in order to shift the inertia's "measured-about" point one
 * must know whether either the starting or final inertia is central,
 * because we must always shift inertias by passing through the central inertia.
 * So this class provides operations for doing the shifting, but expects
 * to be told by the client where to find the center of mass.
 *
 * Re-expressing an Inertia in a different coordinate system does not
 * entail a change of physical meaning in the way that shifting it to a
 * different point does. Note that because inertia is a tensor, there is
 * a "left frame" and "right frame". For our purposes, these
 * will always be the same so we'll only indicate the frame 
 * once, as in 'I_pt_frame'. This should be understood to mean
 * 'frame_I_pt_frame' and re-expressing an Inertia requires both a
 * left and right multiply by the rotation matrix. So I_OB_B is the
 * inertia about body B's origin point OB, expressed in B, while
 * I_OB_G is the same physical quantity but expressed in Ground (the latter
 * is a component of the Spatial Inertia). Conversion is done like this:
 *    I_OB_G = R_GB * I_OB_B * R_BG  (and recall that R_GB=~R_BG)
 * The central inertia would be I_CB_B for body B.
 *
 * A Inertia is a symmetric matrix and is positive definite for
 * nonsingular bodies (that is, a body composed of at least three
 * noncollinear point masses).
 */
class Inertia {
public:
    /// Default is a NaN-ed out mess to avoid accidents.
    Inertia() : I_OF_F(NaN) {}

    /// Create a principal inertia matrix with identical diagonal elements.
    /// Most commonly we create Inertia(0) for initialization of an inertia
    /// calculation (that is also the inertia of a point mass located at
    /// the origin). This can also be used with a non-zero value as the
    /// inertia of a sphere centered at the origin.
    explicit Inertia(const Real& r) : I_OF_F(r) { }

    /// Create an inertia matrix from a vector of the *moments* of
    /// inertia (the inertia matrix diagonal) and optionally a vector of
    /// the *products* of inertia (the off-diagonals). Moments are
    /// in the order xx,yy,zz; products are xy,xz,yz.
    explicit Inertia(const Vec3& moments, const Vec3& products=Vec3(0)) {
        setInertia(moments,products);
    }
    /// Create a principal inertia matrix (only non-zero on diagonal).
    Inertia(const Real& xx, const Real& yy, const Real& zz) {
        setInertia(Vec3(xx,yy,zz)); 
    }
    /// This is a general inertia matrix. Note the order of these
    /// arguments: moments of inertia first, then products of inertia.
    Inertia(const Real& xx, const Real& yy, const Real& zz,
            const Real& xy, const Real& xz, const Real& yz) {
        setInertia(xx,yy,zz,xy,xz,yz); 
    }

    /// Given a point mass located at a given point p in some frame F, 
    /// construct I_OF_F, that is, the inertia of that point mass about
    /// F's origin, expressed in F (that is, in F's coordinate system). 
    ///
    /// For a collection of point masses, you can just add these together to
    /// produce a composite inertia as long as all the vectors are
    /// measured from the same point and expressed in the same frame.
    Inertia(const Vec3& p, const Real& m) {
        Mat33& t = I_OF_F;
        const Real& x = p[0]; const Real mx=m*x, mxx=mx*x;
        const Real& y = p[1]; const Real my=m*y, myy=my*y;
        const Real& z = p[2]; const Real mz=m*z, mzz=mz*z;

        t(0,0)          =  myy + mzz;
        t(1,1)          =  mxx + mzz;
        t(2,2)          =  mxx + myy;
        t(0,1) = t(1,0) = -mx*y;
        t(0,2) = t(2,0) = -mx*z;
        t(1,2) = t(2,1) = -my*z;
    }

    /// We only look at the lower triangle, but fill in the whole matrix.
    Inertia(const Inertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
    }

    /// Construct an Inertia from a 3x3 matrix. The matrix must be symmetric
    /// or very close to symmetric. TODO: there are other tests that should
    /// be performed here to check validity, such as the triangle inequality test.
    explicit Inertia(const Mat33& s) {
        assert(close(s(0,1),s(1,0),s.diag().norm()) 
            && close(s(0,2),s(2,0),s.diag().norm())
            && close(s(1,2),s(2,1),s.diag().norm()));
        setInertia(s(0,0),s(1,1),s(2,2),
            0.5*(s(1,0)+s(0,1)),0.5*(s(2,0)+s(0,2)),0.5*(s(2,1)+s(1,2)));
    }

    /// Copy assignment: only look at the lower triangle, but fill in the whole matrix.
    Inertia& operator=(const Inertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
        return *this;
    }

    /// Add in another inertia. Frames and reference point must be the same but
    /// we can't check. Only look at the lower triangle, but fill in the whole matrix.
    Inertia& operator+=(const Inertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) += s(0,0); t(1,1) += s(1,1);  t(2,2) += s(2,2);
        t(0,1) = (t(1,0) += s(1,0));
        t(0,2) = (t(2,0) += s(2,0));
        t(1,2) = (t(2,1) += s(2,1));
        return *this;
    }

    /// Subtract off another inertia. Frames and reference point must be the same but
    /// we can't check. Only look at the lower triangle, but fill in the whole matrix.
    Inertia& operator-=(const Inertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) -= s(0,0); t(1,1) -= s(1,1);  t(2,2) -= s(2,2);
        t(0,1) = (t(1,0) -= s(1,0));
        t(0,2) = (t(2,0) -= s(2,0));
        t(1,2) = (t(2,1) -= s(2,1));
        return *this;
    }

    /// Scale an Inertia.
    Inertia& operator*=(const Real& r) {
        I_OF_F *= r;
        return *this;
    }
    /// Scale an Inertia.
    Inertia& operator/=(const Real& r) {
        I_OF_F /= r;
        return *this;
    }
    /// Set an inertia to have only principal moments (that is, it will
    /// be diagonal). TODO: should check validity.
    Inertia& setInertia(const Real& xx, const Real& yy, const Real& zz) {
        I_OF_F = 0.; I_OF_F(0,0) = xx; I_OF_F(1,1) = yy;  I_OF_F(2,2) = zz;
        return *this;
    }

    /// Set principal moments and optionally off-diagonal terms. Behaves
    /// like an assignment statement. TODO: should check validity.
    Inertia& setInertia(const Vec3& moments, const Vec3& products=Vec3(0)) {
        Mat33& t = I_OF_F;
        t.diag() = moments;
        t(0,1) = t(1,0) = products[0];
        t(0,2) = t(2,0) = products[1];
        t(1,2) = t(2,1) = products[2];
        return *this;
    }

    /// Set this Inertia to a general inertia matrix. Note the order of these
    /// arguments: moments of inertia first, then products of inertia.
    /// Behaves like an assignment statement. TODO: should check validity.
    Inertia& setInertia(const Real& xx, const Real& yy, const Real& zz,
                    const Real& xy, const Real& xz, const Real& yz) {
        Mat33& t = I_OF_F;
        t(0,0) = xx; t(1,1) = yy;  t(2,2) = zz;
        t(0,1) = t(1,0) = xy;
        t(0,2) = t(2,0) = xz;
        t(1,2) = t(2,1) = yz;
        return *this;
    }

    /// Assume that the current inertia is about the F frame's origin OF, and
    /// expressed in F. Given the vector from OF to the body center of mass CF,
    /// and the total mass of the body, we can shift the inertia to the center
    /// of mass. This produces a new Inertia I' whose (implicit) frame F' is
    /// aligned with F but has origin CF (an inertia like that is called a "central
    /// inertia". I' = I - Icom where Icom is the inertia of a fictitious
    /// point mass of mass mtot located at CF (measured in F) about OF.
    inline Inertia shiftToMassCenter(const Vec3& CF, const Real& mtot) const;

    /// Assuming that the current inertia I is a central inertia (that is, it is
    /// inertia about the body center of mass CF), shift it to some other point p
    /// measured from the center of mass. This produces a new inertia I' about
    /// the point p given by I' = I + Ip where Ip is the inertia of a fictitious
    /// point mass of mass mtot (the total body mass) located at p, about CF.
    inline Inertia shiftFromMassCenter(const Vec3& p, const Real& mtot) const;

    /// Re-express this inertia from frame F to frame B, given the orientation
    /// of B in F. This is a similarity transform since rotation matrices are
    /// orthogonal.
    Inertia reexpress(const Rotation& R_FB) const {
        return Inertia(~R_FB * I_OF_F * R_FB); // TODO can do better due to symmetry
    }

    Real trace() const {return I_OF_F(0,0) + I_OF_F(1,1) + I_OF_F(2,2);}

    // Note that we are copying into the Mat33 here so this will still
    // work if we decide to switch to a SymMat33 when they exist.
    Mat33 toMat33() const {
        return I_OF_F;
    }
    Vec3 getMoments()  const {return I_OF_F.diag();}
    Vec3 getProducts() const {return Vec3(I_OF_F(1,0), I_OF_F(2,0), I_OF_F(2,1));}

    // Inertia factory for some common mass elements. Each defines its own
    // frame aligned (when possible) with principal moments. Each has unit
    // mass and its center of mass located at the origin. Use this with shiftFromCOM()
    // to move it somewhere else, and with xform() to express the inertia
    // in another frame. Mass enters linearly in inertia, so just multiply
    // by the actual mass to scale these properly.
    static Inertia pointMass() {return Inertia(0.);}

    // This is an exception -- the mass center will be at p.
    static Inertia pointMassAt(const Vec3& p) {
        const Real& x = p[0]; const Real xx=x*x;
        const Real& y = p[1]; const Real yy=y*y;
        const Real& z = p[2]; const Real zz=z*z;
        return Inertia(yy + zz, xx + zz, xx + yy,
                       -x*y, -x*z, -y*z);
    }

    static Inertia sphere(const Real& r) {return Inertia(0.4*r*r);}

    // Cylinder is aligned along z axis, use radius and half-length.
    // If r==0 this is a thin rod; hz=0 it is a thin disk.
    static Inertia cylinderAlongZ(const Real& r, const Real& hz) {
        const Real Ixx = 0.25*r*r + (Real(1)/Real(3))*hz*hz;
        return Inertia(Ixx,Ixx,0.5*r*r);
    }

    // Cylinder is aligned along y axis, use radius and half-length.
    // If r==0 this is a thin rod; hy=0 it is a thin disk.
    static Inertia cylinderAlongY(const Real& r, const Real& hy) {
        const Real Ixx = 0.25*r*r + (Real(1)/Real(3))*hy*hy;
        return Inertia(Ixx,0.5*r*r,Ixx);
    }

    // Cylinder is aligned along x axis, use radius and half-length.
    // If r==0 this is a thin rod; hx=0 it is a thin disk.
    static Inertia cylinderAlongX(const Real& r, const Real& hx) {
        const Real Iyy = 0.25*r*r + (Real(1)/Real(3))*hx*hx;
        return Inertia(0.5*r*r,Iyy,Iyy);
    }

    // Brick given by half-lengths in each direction. One dimension zero
    // gives inertia of a thin rectangular sheet; two zero gives inertia
    // of a thin rod in the remaining direction.
    static Inertia brick(const Real& hx, const Real& hy, const Real& hz) {
        const Real oo3 = Real(1)/Real(3);
        const Real hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return Inertia(oo3*(hy2+hz2), oo3*(hx2+hz2), oo3*(hx2+hy2));
    }

    // Ellipsoid given by half-lengths in each direction.
    static Inertia ellipsoid(const Real& hx, const Real& hy, const Real& hz) {
        const Real hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return Inertia(0.2*(hy2+hz2), 0.2*(hx2+hz2), 0.2*(hx2+hy2));
    }


private:
    // Check whether a and b are the same except for numerical error which
    // is a reasonable fraction of the overall scale, which is passed in.
    static bool close(const Real& a, const Real& b, const Real& scale) {
        const Real okErr = SignificantReal*std::abs(scale);
        const Real err = std::abs(a-b);
        return err <= okErr;
    }

private:
    // Inertia expressed in frame F and about F's origin OF. Note that frame F
    // is implicit here; all we actually have are the inertia scalars. This is 
    // a symmetric matrix but we keep all the elements here, and manage them
    // so that the reflected elements are *exactly* equal.
    // TODO: should use a SymMat33 type.
    Mat33 I_OF_F;
    friend Vec3 operator*(const Inertia& i, const Vec3& w);
};

inline Inertia operator+(const Inertia& l, const Inertia& r) {
    return Inertia(l) += r;
}
inline Inertia operator-(const Inertia& l, const Inertia& r) {
    return Inertia(l) -= r;
}
inline Inertia operator*(const Inertia& i, const Real& r) {
    return Inertia(i) *= r;
}
inline Inertia operator*(const Real& r, const Inertia& i) {
    return Inertia(i) *= r;
}
inline Vec3 operator*(const Inertia& i, const Vec3& w) {
    return i.I_OF_F * w;
}
inline Inertia operator/(const Inertia& i, const Real& r) {
    return Inertia(i) /= r;
}
inline Inertia Inertia::shiftToMassCenter(const Vec3& CF, const Real& mtot) const {
    return *this - Inertia(CF, mtot);
}
inline Inertia Inertia::shiftFromMassCenter(const Vec3& p, const Real& mtot) const {
    return *this + Inertia(p, mtot);
}

inline bool
operator==(const Inertia& i1, const Inertia& i2) {
    return i1.toMat33() == i2.toMat33();    // TODO should use underlying rep
}

SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const Inertia&);


/**
 * This class contains the mass, center of mass, and inertia of a rigid body B.
 * The center of mass is a vector from B's origin, expressed in the B frame.
 * The inertia is taken about the B origin, and expressed in B.
 */
class MassProperties {
public:
    MassProperties() { setMassProperties(0.,Vec3(0.),Inertia()); }
    MassProperties(const Real& m, const Vec3& com, const Inertia& inertia)
      { setMassProperties(m,com,inertia); }

    /// Set mass, center of mass, and inertia. Behaves like an assignment in that
    /// a reference to the modified MassProperties object is returned.
    MassProperties& setMassProperties(const Real& m, const Vec3& com, const Inertia& inertia)
      { mass=m; comInB=com; inertia_OB_B=inertia; return *this; }

    const Real&    getMass()       const { return mass; }
    const Vec3&    getMassCenter() const { return comInB; }
    const Inertia& getInertia()    const { return inertia_OB_B; }

    Inertia calcCentralInertia() const {
        return inertia_OB_B - Inertia(comInB, mass);
    }
    Inertia calcShiftedInertia(const Vec3& newOriginB) const {
        return calcCentralInertia() + Inertia(newOriginB-comInB, mass);
    }
    Inertia calcTransformedInertia(const Transform& X_BC) const {
        return calcShiftedInertia(X_BC.T()).reexpress(X_BC.R());
    }
    MassProperties calcShiftedMassProps(const Vec3& newOriginB) const {
        return MassProperties(mass, comInB-newOriginB,
                              calcShiftedInertia(newOriginB));
    }
    // Transform a body's mass properties from the (implicit) B frame
    // to a new frame C.
    MassProperties calcTransformedMassProps(const Transform& X_BC) const {
        return MassProperties(mass, ~X_BC*comInB, calcTransformedInertia(X_BC));
    }

    // Re-express these mass properties in frame C. Currently the mass properties
    // are expressed in the (implicit) frame B, so we need the Rotation matrix
    // that takes us from B to C.
    MassProperties reexpress(const Rotation& R_BC) const {
        return MassProperties(mass, ~R_BC*comInB, inertia_OB_B.reexpress(R_BC));
    }

    bool isExactlyMassless()   const { return mass==0.; }
    bool isNearlyMassless(const Real& tol=SignificantReal) const { 
        return mass <= tol; 
    }

    bool isExactlyCentral() const { return comInB==Vec3(0); }
    bool isNearlyCentral(const Real& tol=SignificantReal) const {
        return comInB.normSqr() <= tol*tol;
    }

    SpatialMat toSpatialMat() const {
        SpatialMat M;
        M(0,0) = inertia_OB_B.toMat33();
        M(0,1) = crossMat(comInB);
        M(1,0) = ~M(0,1);
        M(1,1) = mass; // a diagonal matrix
        return M;
    }

    /// Caution: this does not have the same layout in memory as
    /// a SpatialMat, although it has the same logical layout.
    Mat66 toMat66() const {
        Mat66 M;
        M.updSubMat<3,3>(0,0) = inertia_OB_B.toMat33();
        M.updSubMat<3,3>(0,3) = crossMat(comInB);
        M.updSubMat<3,3>(3,0) = ~M.getSubMat<3,3>(0,3);
        M.updSubMat<3,3>(3,3) = mass; // a diagonal matrix
        return M;
    }

private:
    Real     mass;
    Vec3     comInB;         // meas. from B origin, expr. in B
    Inertia  inertia_OB_B;   // about B origin, expr. in B
};

SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const MassProperties&);

} // namespace SimTK

#endif // SimTK_SIMMATRIX_MASS_PROPERTIES_H_
