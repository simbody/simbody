#ifndef MASS_PROPERTIES_H_
#define MASS_PROPERTIES_H_

/** @file
 *
 * These are utility classes for dealing with mass properties, particularly
 * those messy inertias.
 */

#include "vec3.h"
#include "Mat33.h"
#include "fixedVector.h"
#include "cdsList.h"

typedef float_type Real;

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
class Inertia : public Mat33 {
public:
    /// Default is the inertia of a point mass about that point, i.e. 0.
    Inertia() : Mat33(0.) {}

    /// Note the order of these arguments: moments of inertia first, then 
    /// products of inertia.
    Inertia(const Real& xx, const Real& yy, const Real& zz,
            const Real& xy, const Real& xz, const Real& yz)
        { setInertia(xx,yy,zz,xy,xz,yz); }

    /// Given a point mass located at a given point p in some frame F, 
    /// construct I_OF_F, that is, the inertia of that point mass about
    /// F's origin, expressed in F. 
    ///
    /// For a collection of point masses, you can just add these together to
    /// produce a composite inertia as long as all the vectors are
    /// measured from the same point and expressed in the same frame.
    Inertia(const double& m, const Vec3& p) {
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

    // We only look at the lower triangles, but fill in the whole matrix.
    Inertia(const Inertia& s) {
        Mat33& t = *this;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
    }

    Inertia& operator=(const Inertia& s) {
        Mat33& t = *this;
        t(0,0) = s(0,0); t(1,1) = s(1,1);  t(2,2) = s(2,2);
        t(0,1) = (t(1,0) = s(1,0));
        t(0,2) = (t(2,0) = s(2,0));
        t(1,2) = (t(2,1) = s(2,1));
        return *this;
    }

    Inertia& operator+=(const Inertia& s) {
        Mat33& t = *this;
        t(0,0) += s(0,0); t(1,1) += s(1,1);  t(2,2) += s(2,2);
        t(0,1) = (t(1,0) += s(1,0));
        t(0,2) = (t(2,0) += s(2,0));
        t(1,2) = (t(2,1) += s(2,1));
        return *this;
    }

    Inertia& operator-=(const Inertia& s) {
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
};

inline Inertia operator+(const Inertia& l, const Inertia& r) {
    Inertia t(l); t += r; return t;
}
inline Inertia operator-(const Inertia& l, const Inertia& r) {
    Inertia t(l); t -= r; return t;
}

/**
 * This class contains the mass, centroid, and inertia of a rigid body B.
 * The centroid is a vector from B's origin, expressed in the B frame.
 * The inertia is taken about the B origin, and expressed in B.
 */
class MassProperties {
public:
    MassProperties() { }
    MassProperties(const Real& m, const Vec3& com, const Inertia& inertia)
      { setMassProperties(m,com,inertia); }

    void setMassProperties(const Real& m, const Vec3& com, const Inertia& inertia)
      { mass=m; comInB=com; inertia_OB_B=inertia; }

    const Real&    getMass()    const { return mass; }
    const Vec3&    getCOM()     const { return comInB; }
    const Inertia& getInertia() const { return inertia_OB_B; }

    Inertia calcCentroidalInertia() const {
        return inertia_OB_B - Inertia(mass, comInB);
    }
    Inertia calcShiftedInertia(const Vec3& newOriginB) const {
        return calcCentroidalInertia() + Inertia(mass, newOriginB-comInB);
    }
    MassProperties calcShiftedMassProps(const Vec3& newOriginB) const {
        return MassProperties(mass, comInB-newOriginB,
                              calcShiftedInertia(newOriginB));
    }

    bool isMassless()   const { return mass==0.; }

private:
    Real    mass;
    Vec3    comInB;         // meas. from B origin, expr. in B
    Inertia inertia_OB_B;   // about B origin, expr. in B
};

#endif /* MASS_PROPERTIES_H_ */
