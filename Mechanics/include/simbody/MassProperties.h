#ifndef MASS_PROPERTIES_H_
#define MASS_PROPERTIES_H_

/** @file
 *
 * These are utility classes for dealing with mass properties, particularly
 * those messy inertias.
 */

#include "simbody/SimbodyCommon.h"

#include <iostream>

namespace simtk {

/**
 * The physical meaning of an inertia is the current distribution of
 * a rigid body's mass about a *particular* point. If that point is the
 * center of mass of the body, then the measured inertia is called
 * the "central inertia" of that body. To write down the inertia, we
 * need to calculate the six scalars of the inertia tensor, which
 * is a symmetric 3x3 matrix. These scalars must be expressed in 
 * an arbitrary but specified coordinate system. So a MatInertia
 * is meaningful only in conjunction with a particular frame, fixed
 * to the body, whose origin is the point about which the inertia is being
 * measured, and in whose coordinate system this measurement is being
 * expressed.
 *
 * This low-level MatInertia class does not attempt to keep track of *which*
 * frame it is in. It concentrates instead on providing construction and
 * operations involving inertia which can proceed using only an implicit
 * frame F. Clients of this class are responsible for keeping track of that frame.
 * In particular, in order to shift the inertia's "measured-about" point one
 * must know whether either the starting or final inertia is central,
 * because we must always shift inertias by passing through the central inertia.
 * So this class provides operations for doing the shifting, but expects
 * to be told by the client where to find the center of mass.
 *
 * Re-expressing an MatInertia in a different coordinate system does not
 * entail a change of physical meaning in the way that shifting it to a
 * different point does. Note that because inertia is a tensor, there is
 * a "left frame" and "right frame". For our purposes, these
 * will always be the same so we'll only indicate the frame 
 * once, as in 'I_pt_frame'. This should be understood to mean
 * 'frame_I_pt_frame' and re-expressing an MatInertia requires both a
 * left and right multiply by the rotation matrix. So I_OB_B is the
 * inertia about body B's origin point OB, expressed in B, while
 * I_OB_G is the same physical quantity but expressed in Ground (the latter
 * is a component of the Spatial MatInertia). Conversion is done like this:
 *    I_OB_G = R_GB * I_OB_B * R_BG  (and recall that R_GB=~R_BG)
 * The central inertia would be I_CB_B for body B.
 *
 * A MatInertia is a symmetric matrix and is positive definite for
 * nonsingular bodies (that is, a body composed of at least three
 * noncollinear point masses).
 */
class MatInertia {
public:
    /// Default is a NaN-ed out mess to avoid accidents.
    MatInertia() : I_OF_F(NTraits<Real>::getNaN()) {}

    /// Create an principal inertia matrix with identical diagonal elements.
    /// Most commonly we create MatInertia(0) for initialization of an inertia
    /// calculation (that is also the inertia of a point mass located at
    /// the origin). This can also be used with a non-zero value as the
    /// inertia of a sphere centered at the origin.
    explicit MatInertia(const Real& r) : I_OF_F(r) { }

    /// Create a principal inertia matrix (only non-zero on diagonal).
    MatInertia(const Real& xx, const Real& yy, const Real& zz)
        { setInertia(xx,yy,zz,0.,0.,0.); }

    /// This is a general inertia matrix. Note the order of these
    /// arguments: moments of inertia first, then products of inertia.
    MatInertia(const Real& xx, const Real& yy, const Real& zz,
            const Real& xy, const Real& xz, const Real& yz)
        { setInertia(xx,yy,zz,xy,xz,yz); }

    /// Given a point mass located at a given point p in some frame F, 
    /// construct I_OF_F, that is, the inertia of that point mass about
    /// F's origin, expressed in F (that is, in F's coordinate system). 
    ///
    /// For a collection of point masses, you can just add these together to
    /// produce a composite inertia as long as all the vectors are
    /// measured from the same point and expressed in the same frame.
    MatInertia(const Vec3& p, const Real& m) {
        Mat33& t = I_OF_F;
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
    MatInertia(const MatInertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
    }

    explicit MatInertia(const Mat33& s) {
        assert(close(s(0,1),s(1,0)) 
            && close(s(0,2),s(2,0))
            && close(s(1,2),s(2,1)));
        setInertia(s(0,0),s(1,1),s(2,2),
            0.5*(s(1,0)+s(0,1)),0.5*(s(2,0)+s(0,2)),0.5*(s(2,1)+s(1,2)));
    }

    MatInertia& operator=(const MatInertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
        return *this;
    }

    MatInertia& operator+=(const MatInertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) += s(0,0); t(1,1) += s(1,1);  t(2,2) += s(2,2);
        t(0,1) = (t(1,0) += s(1,0));
        t(0,2) = (t(2,0) += s(2,0));
        t(1,2) = (t(2,1) += s(2,1));
        return *this;
    }

    MatInertia& operator-=(const MatInertia& src) {
        Mat33&       t = I_OF_F;
        const Mat33& s = src.I_OF_F;
        t(0,0) -= s(0,0); t(1,1) -= s(1,1);  t(2,2) -= s(2,2);
        t(0,1) = (t(1,0) -= s(1,0));
        t(0,2) = (t(2,0) -= s(2,0));
        t(1,2) = (t(2,1) -= s(2,1));
        return *this;
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
    /// of mass. This produces a new MatInertia I' whose (implicit) frame F' is
    /// aligned with F but has origin CF (an inertia like that is called a "central
    /// inertia". I' = I - Icom where Icom is the inertia of a fictitious
    /// point mass of mass mtot located at CF (measured in F) about OF.
    inline MatInertia shiftToCOM(const Vec3& CF, const Real& mtot) const;

    /// Assuming that the current inertia I is a central inertia (that is, it is
    /// inertia about the body center of mass CF), shift it to some other point p
    /// measured from the center of mass. This produces a new inertia I' about
    /// the point p given by I' = I + Ip where Ip is the inertia of a fictitious
    /// point mass of mass mtot (the total body mass) located at p, about CF.
    inline MatInertia shiftFromCOM(const Vec3& p, const Real& mtot) const;

    MatInertia shift(const Vec3& from, const Vec3& to, const Real& mtot) const {
        return shiftToCOM(from,mtot).shiftFromCOM(to,mtot);
    }

    /// Re-express this inertia from frame F to frame B, given the orientation
    /// of B in F. This is a similarity transform since rotation matrices are
    /// orthogonal.
    MatInertia rotate(const Mat33& R_FB) const {
        return MatInertia(~R_FB * I_OF_F * R_FB); // TODO can do better due to symmetry
    }

    Real trace() const {return I_OF_F(0,0) + I_OF_F(1,1) + I_OF_F(2,2);}

    // Note that we are copying into the Mat33 here so this will still
    // work if we decide to switch to a SymMat33 when they exist.
    Mat33 toMat33() const {
        return I_OF_F;
    }

    // MatInertia factory for some common mass elements. Each defines its own
    // frame aligned (when possible) with principal moments. Each has unit
    // mass and its center of mass located at the origin. Use this with shiftFromCOM()
    // to move it somewhere else, and with xform() to express the inertia
    // in another frame. Mass enters linearly in inertia, so just multiply
    // by the actual mass to scale these properly.
    static MatInertia point() {return MatInertia(0.);}
    static MatInertia sphere(const Real& r) {return MatInertia(0.4*r*r);}

    // Cylinder is aligned along z axis, use radius and half-length.
    // If r==0 this is a thin rod; hz=0 it is a thin disk.
    static MatInertia cylinder(const Real& r, const Real& hz) {
        const Real Ixx = 0.25*r*r + (Real(1)/Real(3))*hz*hz;
        return MatInertia(Ixx,Ixx,0.5*r*r);
    }

    // Brick given by half-lengths in each direction. One dimension zero
    // gives inertia of a thin rectangular sheet; two zero gives inertia
    // of a thin rod in the remaining direction.
    static MatInertia brick(const Real& hx, const Real& hy, const Real& hz) {
        const Real oo3 = Real(1)/Real(3);
        const Real hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return MatInertia(oo3*(hy2+hz2), oo3*(hx2+hz2), oo3*(hx2+hy2));
    }

    // Ellipsoid given by half-lengths in each direction.
    static MatInertia ellipsoid(const Real& hx, const Real& hy, const Real& hz) {
        const Real hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return MatInertia(0.2*(hy2+hz2), 0.2*(hx2+hz2), 0.2*(hx2+hy2));
    }


private:
    //TODO: the tolerance here should be a function of Real's precision
    static bool close(const Real& a, const Real& b) {
        if (fabs(a-b) < 1e-13) return true;
        if (fabs(a-b)/(fabs(a)+fabs(b)) < 0.5e-13) return true;
        return false;
    }

private:
    // MatInertia expressed in frame F and about F's origin OF. Note that frame F
    // is implicit here; all we actually have are the inertia scalars. This is 
    // a symmetric matrix but we keep all the elements here, and manage them
    // so that the reflected elements are *exactly* equal.
    // TODO: should use a SymMat33 type.
    Mat33 I_OF_F;                  
};

inline MatInertia operator+(const MatInertia& l, const MatInertia& r) {
    return MatInertia(l) += r;
}
inline MatInertia operator-(const MatInertia& l, const MatInertia& r) {
    return MatInertia(l) -= r;
}
inline MatInertia MatInertia::shiftToCOM(const Vec3& CF, const Real& mtot) const {
    return *this - MatInertia(CF, mtot);
}
inline MatInertia MatInertia::shiftFromCOM(const Vec3& p, const Real& mtot) const {
    return *this + MatInertia(p, mtot);
}

std::ostream& operator<<(std::ostream& o, const MatInertia&);

/**
 * This class contains the mass, centroid, and inertia of a rigid body B.
 * The centroid is a vector from B's origin, expressed in the B frame.
 * The inertia is taken about the B origin, and expressed in B.
 */
/*
class MassProperties {
public:
    MassProperties() { setMassProperties(0.,Vec3(0.),MatInertia()); }
    MassProperties(const Real& m, const Vec3& com, const MatInertia& inertia)
      { setMassProperties(m,com,inertia); }

    void setMassProperties(const Real& m, const Vec3& com, const MatInertia& inertia)
      { mass=m; comInB=com; inertia_OB_B=inertia; }

    const Real&    getMass()    const { return mass; }
    const Vec3&    getCOM()     const { return comInB; }
    const MatInertia& getInertia() const { return inertia_OB_B; }

    MatInertia calcCentroidalInertia() const {
        return inertia_OB_B - MatInertia(mass, comInB);
    }
    MatInertia calcShiftedInertia(const Vec3& newOriginB) const {
        return calcCentroidalInertia() + MatInertia(mass, newOriginB-comInB);
    }
    MassProperties calcShiftedMassProps(const Vec3& newOriginB) const {
        return MassProperties(mass, comInB-newOriginB,
                              calcShiftedInertia(newOriginB));
    }

    bool isMassless()   const { return mass==0.; }

private:
    Real mass;
    Vec3 comInB;         // meas. from B origin, expr. in B
    MatInertia inertia_OB_B;   // about B origin, expr. in B
};
*/
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
    Frame() : Rot_RF(1.), Loc_RF(0.) { }
    Frame(const Mat33& axesInR, const Vec3& originInR) {setFrame(axesInR,originInR);}
    Frame(const Vec3& originInR) : Rot_RF(1.), Loc_RF(originInR) { }
    // default copy, assignment, destructor

    void setFrame(const Mat33& axesInR, const Vec3& originInR) 
        {Rot_RF=axesInR; Loc_RF=originInR;}

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

} // namespace simtk

#endif /* MASS_PROPERTIES_H_ */
