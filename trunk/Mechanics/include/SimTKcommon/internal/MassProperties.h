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
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
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

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/Orientation.h"

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

// These are templatized by precision (float or double).
template <class P> class Gyration_;
template <class P> class SpatialInertia_;
template <class P> class ArticulatedInertia_;

// These "no trailing underscore" typedefs use whatever the 
// compile-time precision is set to.
typedef Gyration_<Real>             Gyration;
typedef SpatialInertia_<Real>       SpatialInertia;
typedef ArticulatedInertia_<Real>   ArticulatedInertia;

class Inertia;

/**
 * A Gyration matrix is a mass covariance matrix with units of length squared.
 * This can also be considered a unit-mass inertia matrix. Gyration is measured 
 * about some point OF and expressed in some frame F, but we don't know anything
 * about that frame here.
 */
template <class P>
class SimTK_SimTKCOMMON_EXPORT Gyration_ {
    typedef P           RealP;
    typedef Vec<3,P>    Vec3P;
    typedef SymMat<3,P> SymMat33P;
    typedef Mat<3,3,P>  Mat33P;
public:
    /// Default is a NaN-ed out mess to avoid accidents, even in Release mode.
    Gyration_() : G_OF_F(NaN) {}

    // Default copy constructor, copy assignment, destructor.

    /// Create a principal gyration matrix with identical diagonal elements.
    explicit Gyration_(const RealP& r) : G_OF_F(r) { }

    /// Create a gyration matrix from a vector of the *moments* of
    /// gyration (the gyration matrix diagonal) and optionally a vector of
    /// the *products* of gyration (the off-diagonals). Moments are
    /// in the order xx,yy,zz; products are xy,xz,yz.
    explicit Gyration_(const Vec3P& moments, const Vec3P& products=Vec3P(0)) {
        setGyration(moments,products);
    }
    /// Create a principal gyration matrix (only non-zero on diagonal).
    Gyration_(const RealP& xx, const RealP& yy, const RealP& zz) {
        setGyration(Vec3P(xx,yy,zz)); 
    }
    /// This is a general gyration matrix. Note the order of these
    /// arguments: moments of gyration first, then products of gyration.
    Gyration_(const RealP& xx, const RealP& yy, const RealP& zz,
              const RealP& xy, const RealP& xz, const RealP& yz) {
        setGyration(xx,yy,zz,xy,xz,yz); 
    }

    /// Construct a Gyration from a symmetric 3x3 matrix. The diagonals must
    /// be nonnegative and satisfy the triangle inequality.
    explicit Gyration_(const SymMat33P& G) : G_OF_F(G) {
        const Vec3P& d = G.diag();
        SimTK_ERRCHK3(d >= 0, "Gyration(SymMat3)",
            "Diagonals of a Gyration matrix must be nonnegative; got %g,%g,%g.",
            (double)d[0],(double)d[1],(double)d[2]);
        SimTK_ERRCHK3(d[0]+d[1]>=d[2] && d[0]+d[2]>=d[1] && d[1]+d[2]>=d[0],
            "Gyration(SymMat3)",
            "Diagonals of a Gyration matrix must satisfy the triangle inequality; got %g,%g,%g.",
            (double)d[0],(double)d[1],(double)d[2]);
    }

    /// Construct an Inertia from a 3x3 matrix. The matrix must be symmetric
    /// or very close to symmetric. TODO: there are other tests that should
    /// be performed here to check validity, such as the triangle inequality test.
    explicit Gyration_(const Mat33P& s) {
        SimTK_ERRCHK(   close(s(0,1),s(1,0),s.diag().norm()) 
                     && close(s(0,2),s(2,0),s.diag().norm())
                     && close(s(1,2),s(2,1),s.diag().norm()), 
                     "Gyration(Mat33)", "The supplied matrix was not symmetric.");
        setGyration(s(0,0),s(1,1),s(2,2),
            0.5*(s(1,0)+s(0,1)),0.5*(s(2,0)+s(0,2)),0.5*(s(2,1)+s(1,2)));
    }

    /// Add in another gyration matrix. Frames and reference point must be the same but
    /// we can't check. (6 flops)
    Gyration_& operator+=(const Gyration_& G) {
        G_OF_F += G;
        return *this;
    }

    /// Subtract off another gyration matrix. Frames and reference point must 
    /// be the same but we can't check. (6 flops)
    Gyration_& operator-=(const Gyration_& G) {
        G_OF_F -= G;
        return *this;
    }

    /// Set a gyration matrix to have only principal moments (that is, it
    /// will be diagonal). TODO: should check validity.
    Gyration_& setGyration(const RealP& xx, const RealP& yy, const RealP& zz) {
        assert(!"check validity");
        G_OF_F = 0.; G_OF_F(0,0) = xx; G_OF_F(1,1) = yy;  G_OF_F(2,2) = zz;
        return *this;
    }

    /// Set principal moments and optionally off-diagonal terms. Behaves
    /// like an assignment statement. TODO: should check validity.
    Gyration_& setGyration(const Vec3P& moments, const Vec3P& products=Vec3P(0)) {
        G_OF_F.updDiag()  = moments;
        G_OF_F.updLower() = products;
        assert(!"check validity");
        return *this;
    }

    /// Set this Gyration to a general matrix. Note the order of these
    /// arguments: moments of gyration first, then products of gyration.
    /// Behaves like an assignment statement. TODO: should check validity.
    Gyration_& setGyration(const RealP& xx, const RealP& yy, const RealP& zz,
                           const RealP& xy, const RealP& xz, const RealP& yz) {
        setGyration(Vec3P(xx,yy,zz), Vec3P(xy,xz,yz));
        return *this;
    }

    /// Assume that the current gyration is about the F frame's origin OF, and
    /// expressed in F. Given the vector from OF to the centroid CF,
    /// we can shift the gyration to the center of mass. This produces a new 
    /// Gyration matrix G' whose (implicit) frame F' is aligned with F but has 
    /// origin CF (an gyration matrix like that is called a "central
    /// gyration". G' = G - Gcom where Gcom is the gyration of a fictitious
    /// point located at CF (measured in F) taken about OF. (17 flops)
    Gyration_ shiftToCentroid(const Vec3P& CF) const 
    {   return G_OF_F - Gyration_(CF); }

    /// Assuming that the current Gyration G is a central gyration (that is, it is
    /// gyration about the body centroid CF), shift it to some other point p
    /// measured from the centroid. This produces a new inertia G' about the
    /// point p given by G' = G + Gp where Gp is the gyration of a fictitious
    /// point located at p, taken about CF. (17 flops)
    Gyration_ shiftFromCentroid(const Vec3P& p) const
    {   return G_OF_F + Gyration_(p); }

    /// Return a new gyration matrix like this one but re-expressed in another 
    /// frame (leaving the origin point unchanged). Call this gyration matrix
    /// G_OF_F, that is, it is taken about the origin of some frame F, and 
    /// expressed in F. We want to return G_OF_B, the same gyration matrix,
    /// still taken about the origin of F, but expressed in the B frame, given
    /// by G_OF_B=R_BF*G_OF_F*R_FB where R_BF is the rotation matrix giving
    /// the orientation of frame F in B. As a pair of 3x3 multiplies, this 
    /// computation would be 90 flops, but we can take advantage of the 
    /// symmetry of G and orthogonality of R to get it down to 57 flops using
    /// a trick reported in Featherstone's 2008 book.
    /// @see reexpressInPlace()
    Gyration_ reexpress(const Rotation& R_BF) const {
        return Gyration_(R_BF.reexpressSymMat33(G_OF_F));
    }

    /// Re-express this gyration matrix in another frame, changing the object
    /// in place; see reexpress() if you want to leave this object unmolested
    /// and get a new one instead.
    Gyration_& reexpressInPlace(const Rotation& R_BF) {
        G_OF_F = R_BF.reexpressSymMat33(G_OF_F);
    }

    RealP trace() const {return G_OF_F.trace();}

    /// Obtain a reference to the underlying symmetric matrix type.
    const SymMat33P& asSymMat33() const {return G_OF_F;}
    /// Obtain the gyration moments (diagonal of the Gyration matrix) as a Vec3.
    const Vec3P& getMoments()  const {return G_OF_F.getDiag();}
    /// Obtain the gyration products (off-diagonal of the Gyration matrix)
    /// as a Vec3 with elements ordered xx, xy, yz.
    const Vec3P& getProducts() const {return G_OF_F.getLower();}

    // Gyration matrix factories for some common mass elements. Each defines its
    // own frame aligned (when possible) with principal moments. Each has unit
    // mass and its center of mass located at the origin (usually). Use this with 
    // shiftFromMassCenter() to move it somewhere else, and with reexpress() to 
    // express the Gyration matrix in another frame.

    /// Create a Gyration matrix for a point located at the origin -- that is,
    /// an all-zero matrix.
    static Gyration_ pointMassAtOrigin() {return Gyration_(0);}

    /// Create a Gyration matrix for a point located at a given location.
    static Gyration_ pointMassAt(const Vec3P& p) {return Gyration_(p);}

    /// Create a Gyration matrix for a sphere of radius \a r.
    static Gyration_ sphere(const RealP& r) {return Gyration_(0.4*r*r);}

    /// Cylinder is aligned along z axis, use radius and half-length.
    /// If r==0 this is a thin rod; hz=0 it is a thin disk.
    static Gyration_ cylinderAlongZ(const RealP& r, const RealP& hz) {
        const RealP Ixx = 0.25*r*r + (RealP(1)/RealP(3))*hz*hz;
        return Gyration_(Ixx,Ixx,0.5*r*r);
    }

    /// Cylinder is aligned along y axis, use radius and half-length.
    /// If r==0 this is a thin rod; hy=0 it is a thin disk.
    static Gyration_ cylinderAlongY(const RealP& r, const RealP& hy) {
        const RealP Ixx = 0.25*r*r + (RealP(1)/RealP(3))*hy*hy;
        return Gyration_(Ixx,0.5*r*r,Ixx);
    }

    /// Cylinder is aligned along x axis, use radius and half-length.
    /// If r==0 this is a thin rod; hx=0 it is a thin disk.
    static Gyration_ cylinderAlongX(const RealP& r, const RealP& hx) {
        const RealP Iyy = 0.25*r*r + (RealP(1)/RealP(3))*hx*hx;
        return Gyration_(0.5*r*r,Iyy,Iyy);
    }

    /// Brick given by half-lengths in each direction. One dimension zero
    /// gives inertia of a thin rectangular sheet; two zero gives inertia
    /// of a thin rod in the remaining direction.
    static Gyration_ brick(const RealP& hx, const RealP& hy, const RealP& hz) {
        const RealP oo3 = RealP(1)/RealP(3);
        const RealP hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return Gyration_(oo3*(hy2+hz2), oo3*(hx2+hz2), oo3*(hx2+hy2));
    }

    /// Ellipsoid given by half-lengths in each direction.
    static Gyration_ ellipsoid(const RealP& hx, const RealP& hy, const RealP& hz) {
        const RealP hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
        return Gyration_(0.2*(hy2+hz2), 0.2*(hx2+hz2), 0.2*(hx2+hy2));
    }

private:
    SymMat33P G_OF_F; 
};

/**
 * A spatial inertia contains the mass, center of mass point, and inertia
 * matrix for a rigid body. This is 10 independent quantities altogether; however,
 * inertia is mass-scaled making it linearly dependent on the mass. Here instead 
 * we represent inertia using a gyration matrix, which is equivalent to the 
 * inertia this body would have if it had unit mass. Then the actual inertia is 
 * given by mass*gyration. In this manner the mass, center of mass location, and 
 * gyration are completely independent so can be changed separately. That means 
 * if you double the mass, you'll also double the inertia as you would expect.
 *
 * Spatial inertia may be usefully viewed as a symmetric spatial matrix, that is, 
 * a 6x6 symmetric matrix arranged as 2x2 blocks of 3x3 matrices. Although this 
 * class represents the spatial inertia in compact form, it supports methods and
 * operators that allow it to behave as though it were a spatial matrix (except
 * much faster to work with). In spatial matrix form, the matrix has the following 
 * interpretation:
 * <pre>
 *               [  m*G   m*px ]
 *           M = [             ]
 *               [ -m*px  m*I  ]
 * </pre>
 * Here m is mass, p is the vector from the body origin to the center of mass, 
 * G is the 3x3 symmetric gyration matrix, and I is a 3x3 identity matrix.
 * "px" indicates the skew symmetric cross product matrix formed from the vector p,
 * so -px=~px.
 * 
 */
template <class P> 
class SimTK_SimTKCOMMON_EXPORT SpatialInertia_ {
    typedef P               RealP;
    typedef Vec<3,P>        Vec3P;
    typedef Gyration_<P>    GyrationP;
    typedef Mat<3,3,P>      Mat33P;
    typedef Rotation        RotationP;  // TODO: need template argument
    typedef Transform       TransformP; //   "
    typedef Inertia         InertiaP;   //   "
public:
    /// The default constructor fills everything with NaN, even in Release mode.
    SpatialInertia_() 
    :   m(nanP()), p(nanP()) {} // gyration is already NaN
    SpatialInertia_(RealP mass, const Vec3P& com, const GyrationP& gyration) 
    :   m(mass), p(com), G(gyration) {}

    // default copy constructor, copy assignment, destructor

    SpatialInertia_& setMass(RealP mass)
    {   SimTK_ERRCHK1(mass >= 0, "SpatialInertia::setMass()",
            "Negative mass %g is illegal.", (double)mass);
        m=mass; return *this; }
    SpatialInertia_& setMassCenter(const Vec3P& com)
    {   p=com; return *this;} 
    SpatialInertia_& setGyration(const GyrationP& gyration) 
    {   SimTK_ERRCHK(gyration.isValid(), "SpatialInertia::setGyration()",
        "Invalid gyration matrix.");
        G=gyration; return *this; }

    RealP            getMass()       const {return m;}
    const Vec3P&     getMassCenter() const {return p;}
    const GyrationP& getGyration()   const {return G;}

    /// Calculate the first mass moment (mass-weighted COM location)
    /// from the mass and COM vector. Cost is 3 inline flops.
    Vec3P calcMassMoment() const {return m*p;}

    /// Calculate the inertia matrix (second mass moment, mass-weighted gyration
    /// matrix) from the mass and gyration matrix. Cost is 6 inline flops.
    InertiaP calcInertia() const; //TODO

    /// Add in a compatible SpatialInertia. This is only valid if both 
    /// SpatialInertias are expressed in the same frame and measured about 
    /// the same point but there is no way for this method to check.
    /// Cost is 10 flops.
    SpatialInertia_& operator+=(const SpatialInertia_& src)
    {   m += src.m; p += src.p; G += src.G; return *this; }

    /// Subtract off a compatible SpatialInertia. This is only valid if both 
    /// SpatialInertias are expressed in the same frame and measured about 
    /// the same point but there is no way for this method to check.
    /// Cost is 10 flops.
    SpatialInertia_& operator-=(const SpatialInertia_& src)
    {   m -= src.m; p -= src.p; G -= src.G; return *this; }

    /// Return a new SpatialInertia object which is the same as this one except
    /// re-expressed in another coordinate frame. We consider this object to
    /// be expressed in some frame F and we're given a rotation matrix we
    /// can use to re-express in a new frame B. Cost is 72 flops.
    /// @see reexpressInPlace()
    SpatialInertia_ reexpress(const RotationP& R_BF) const
    {   return SpatialInertia(*this).reexpressInPlace(R_BF); }

    /// Re-express this SpatialInertia in another frame, modifying the original
    /// object. We return a reference to the object so that you can chain this
    /// operation in the manner of assignment operators. Cost is 72 flops.
    /// @see reexpress() if you want to leave this object unmolested.
    SpatialInertia_& reexpressInPlace(const RotationP& R_BF)
    {   p = R_BF*p; G.reexpressInPlace(R_BF); }

    /// Change origin from OF to OF+S.
    /// 37 flops
    SpatialInertia_ shift(const Vec3P& S) const 
    {   return SpatialInertia_(*this).shiftInPlace(S); }

    SpatialInertia_& shiftInPlace(const Vec3P& S) {
        G.shiftToCentroidInPlace(p);    // change to central gyration
        G.shiftFromCentroidInPlace(S);  // now gyration is about S
        p -= S; // was p=com-OF, now want p'=com-(OF+S)=p-S
    }

    /// Current G_OF_F, want G_OB_B. Cost is 109 flops.
    SpatialInertia_ transform(const TransformP& X_BF) const 
    {   return SpatialInertia_(*this).transformInPlace(X_BF); }

    /// Current G_OF_F, want G_OB_B. Cost is 109 flops.
    SpatialInertia_& transformInPlace(const TransformP& X_BF) {
        reexpressInPlace(X_BF.R()); // get everything in B
        shiftInPlace(X_BF.p());   // now shift to the new origin OB.
    }

private:
    RealP       m;  ///< mass of this rigid body F
    Vec3P       p;  ///< location of body's COM from OF, exp. in F
    GyrationP   G;  ///< mass distribution; inertia is mass*gyration

    static P nanP() {return NTraits<P>::getNaN();} 
};

/// Add two compatible spatial inertias. Cost is 10 inline flops.
template <class P> inline SpatialInertia_<P> 
operator+(const SpatialInertia_<P>& l, const SpatialInertia_<P>& r)
{   return SpatialInertia_<P>(l) += r; } 

/// Subtract one compatible spatial inertia from another. Cost is
/// 10 inline flops.
template <class P> inline SpatialInertia_<P> 
operator-(const SpatialInertia_<P>& l, const SpatialInertia_<P>& r)
{   return SpatialInertia_<P>(l) -= r; } 

/// This operator allows you to re-express a spatial inertia I_OF_F in 
/// assumed frame F into another frame B by writing R_BF*I_OF_F although
/// the transform is really I_OF_B = R_BF*I_OF_F*~R_BF. Note that this
/// is just a rotation of the assumed frame; the origin point is unchanged.
/// Cost is 72 flops.
template <class P> inline SpatialInertia_<P> 
operator*(const Rotation& R_BF, const SpatialInertia_<P>& I_OF_F)
{   return I_OF_F.reexpress(R_BF); } 

/// This operator allows you to efficiently transform (shift origin 
/// and re-express) a spatial inertia I_OF_F in assumed frame F 
/// into another frame B by writing X_BF*I_OF_F although the 
/// transform is really I_OB_B = X_BF*I_OF_F*~X_BF.
/// Cost is 109 flops.
template <class P> inline SpatialInertia_<P> 
operator*(const Transform& X_BF, const SpatialInertia_<P>& I_OF_F)
{   return I_OF_F.transform(X_BF); } 

/**
 * An articulated body inertia (ABI) matrix P(q) contains the spatial inertia 
 * properties that a body appears to have when it is the free base body of 
 * an articulated multibody tree in a given configuration q. Despite the 
 * complex relative motion that occurs within a multibody tree, at any given 
 * configuration q there is still a linear relationship between a spatial 
 * force F applied to a point of the base body and the resulting acceleration 
 * A of that body and that point: F = P(q)*A + c, where c is a velocity-
 * dependent inertial bias force. P is thus analogous to a rigid body 
 * spatial inertia (RBI), but for a body which has other bodies connected to it
 * by joints which are free to move.
 *
 * An ABI P is a symmetric 6x6 spatial matrix, consisting of 2x2 blocks of 3x3 
 * matrices, similar to the RBI. However, unlike the RBI which has only 10 independent
 * elements, all 21 elements of P's lower triangle are significant. For example,
 * the apparent mass of an articulated body depends on which way you push it,
 * and in general there is no well-defined center of mass. This
 * is a much more expensive matrix to manipulate than an RBI. In Simbody's formulation,
 * we only work with ABIs in the Ground frame, so there
 * is never a need to rotate or re-express them. (That is done by rotating RBIs
 * prior to using them to construct the ABIs.) Thus only shifting operations need
 * be performed when transforming ABIs from body to body.
 * Cheap rigid body shifting is done when moving an ABI 
 * within a body or across a prescribed mobilizer; otherwise we have to perform 
 * an articulated shift operation which is quite expensive.
 *
 * For a full discussion of the properties of articulated body inertias, see 
 * Section 7.1 (pp. 119-123) of Roy Featherstone's excellent 2008 book, Rigid 
 * Body Dynamics Algorithms. 
 *
 * In spatial matrix form, an ABI P may be considered to consist of the 
 * following 3x3 subblocks:
 * <pre>
 *          P =  [ J  F ]
 *               [~F  M ]
 * </pre>
 * Here M is a (symmetric) mass distribution, F is a full matrix giving the
 * first mass moment distribution, and J is a (symmetric) inertia matrix.
 */
template <class P> 
class ArticulatedInertia_ {
    typedef P               RealP;
    typedef Vec<3,P>        Vec3P;
    typedef Gyration_<P>    GyrationP;
    typedef Mat<3,3,P>      Mat33P;
    typedef SymMat<3,P>     SymMat33P;
    typedef Mat<2,2,Mat33P> SpatialMatP;
    typedef Rotation        RotationP;  // TODO: need template argument
    typedef Transform       TransformP; //   "
    typedef Inertia         InertiaP;   //   "
public:
    ArticulatedInertia_() {}
    ArticulatedInertia_(const SymMat33P& mass, const Mat33P& massMoment, const SymMat33P& inertia)
    :   M(mass), J(inertia), F(massMoment) {}

    ArticulatedInertia_& setMass      (const SymMat33P& mass)       {M=mass;       return *this;}
    ArticulatedInertia_& setMassMoment(const Mat33P&    massMoment) {F=massMoment; return *this;}
    ArticulatedInertia_& setInertia   (const SymMat33P& inertia)    {J=inertia;    return *this;}

    const SymMat33P& getMass()       const {return M;}
    const Mat33P&    getMassMoment() const {return F;}
    const SymMat33P& getInertia()    const {return J;}

    // default destructor, copy constructor, copy assignment

    ArticulatedInertia_& operator+=(const ArticulatedInertia_& src)
    {   M+=src.M; J+=src.J; F+=src.F; return *this; }
    ArticulatedInertia_& operator-=(const ArticulatedInertia_& src)
    {   M-=src.M; J-=src.J; F-=src.F; return *this; }

    /// Rigid-shift the origin of this Articulated Body Inertia P by a 
    /// shift vector s to produce a new ABI P'. The calculation is 
    /// <pre>
    /// P' =  [ J'  F' ]  =  [ 1  sx ] [ J  F ] [ 1  0 ]
    ///       [~F'  M  ]     [ 0  1  ] [~F  M ] [-sx 1 ]
    /// </pre>
    /// where sx is the cross product matrix of s. Cost is 72 flops.
    SimTK_SimTKCOMMON_EXPORT ArticulatedInertia_ shift(const Vec3P& s) const;

    /// Rigid-shift this ABI in place. 72 flops.
    /// @see shift() for details
    SimTK_SimTKCOMMON_EXPORT ArticulatedInertia_& shiftInPlace(const Vec3P& s);

    const SpatialMatP toSpatialMat() const {
        return SpatialMatP( Mat33P(J),     F,
                              ~F,       Mat33P(M) );
    }
private:
    SymMat33P M;
    SymMat33P J;
    Mat33P    F;
};

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
class SimTK_SimTKCOMMON_EXPORT Inertia {
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
class SimTK_SimTKCOMMON_EXPORT MassProperties {
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
        return calcShiftedInertia(X_BC.p()).reexpress(X_BC.R());
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
        M(0,1) = mass*crossMat(comInB);
        M(1,0) = ~M(0,1);
        M(1,1) = mass; // a diagonal matrix
        return M;
    }

    /// Caution: this does not have the same layout in memory as
    /// a SpatialMat, although it has the same logical layout.
    Mat66 toMat66() const {
        Mat66 M;
        M.updSubMat<3,3>(0,0) = inertia_OB_B.toMat33();
        M.updSubMat<3,3>(0,3) = mass*crossMat(comInB);
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
