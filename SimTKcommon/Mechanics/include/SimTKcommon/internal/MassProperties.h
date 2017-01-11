#ifndef SimTK_SIMMATRIX_MASS_PROPERTIES_H_
#define SimTK_SIMMATRIX_MASS_PROPERTIES_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
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
/** Spatial vectors are used for (rotation,translation) quantities and 
consist of a pair of Vec3 objects, arranged as a 2-vector of 3-vectors. 
Quantities represented this way include
    - spatial velocity     = (angularVelocity,linearVelocity)
    - spatial acceleration = (angularAcceleration,linearAcceleration)
    - generalized forces   = (torque,force)

Spatial configuration has to be handled differently though since
orientation is not a vector quantity. (We use Transform for this concept
which includes a Rotation matrix and a translation Vec3.) **/
typedef Vec<2,   Vec3>              SpatialVec;
/** A SpatialVec that is always single (float) precision regardless of
the compiled-in precision of Real. **/
typedef Vec<2,   Vec<3,float> >    fSpatialVec;
/** A SpatialVec that is always double precision regardless of
the compiled-in precision of Real. **/
typedef Vec<2,   Vec<3,double> >   dSpatialVec;

/** This is the type of a transposed SpatialVec; it does not usually appear
explicitly in user programs. **/
typedef Row<2,   Row3>              SpatialRow;
/** A SpatialRow that is always single (float) precision regardless of
the compiled-in precision of Real. **/
typedef Row<2,   Row<3,float> >    fSpatialRow;
/** A SpatialRow that is always double precision regardless of
the compiled-in precision of Real. **/
typedef Row<2,   Row<3,double> >   dSpatialRow;

/** Spatial matrices are used to hold 6x6 matrices that are best viewed 
as 2x2 matrices of 3x3 matrices; most commonly for spatial and articulated
body inertias and spatial shift matrices. They also arise commonly as 
intermediates in computations involving SpatialVec objects. **/
typedef Mat<2,2, Mat33>             SpatialMat;
/** A SpatialMat that is always single (float) precision regardless of
the compiled-in precision of Real. **/
typedef Mat<2,2, Mat<3,3,float> >  fSpatialMat;
/** A SpatialMat that is always double precision regardless of
the compiled-in precision of Real. **/
typedef Mat<2,2, Mat<3,3,double> > dSpatialMat;

// These are templatized by precision (float or double).
template <class P> class UnitInertia_;
template <class P> class Inertia_;
template <class P> class SpatialInertia_;
template <class P> class ArticulatedInertia_;
template <class P> class MassProperties_;

// The "no trailing underscore" typedefs use whatever the 
// compile-time precision is set to.

/** A unit inertia (gyration) tensor at default precision. **/
typedef UnitInertia_<Real>          UnitInertia;
/** A unit inertia (gyration) tensor at float precision. **/
typedef UnitInertia_<float>        fUnitInertia;
/** A unit inertia (gyration) tensor at double precision. **/
typedef UnitInertia_<double>       dUnitInertia;

/** An inertia tensor at default precision. **/
typedef Inertia_<Real>              Inertia;
/** An inertia tensor at float precision. **/
typedef Inertia_<float>            fInertia;
/** An inertia tensor at double precision. **/
typedef Inertia_<double>           dInertia;

/** Rigid body mass properties at default precision. **/
typedef MassProperties_<Real>       MassProperties;
/** Rigid body mass properties at float precision. **/
typedef MassProperties_<float>     fMassProperties;
/** Rigid body mass properties at double precision. **/
typedef MassProperties_<double>    dMassProperties;

/** A spatial (rigid body) inertia matrix at default precision. **/
typedef SpatialInertia_<Real>       SpatialInertia;
/** A spatial (rigid body) inertia matrix at float precision. **/
typedef SpatialInertia_<float>     fSpatialInertia;
/** A spatial (rigid body) inertia matrix at double precision. **/
typedef SpatialInertia_<double>    dSpatialInertia;

/** An articulated body inertia matrix at default precision. **/
typedef ArticulatedInertia_<Real>    ArticulatedInertia;
/** An articulated body inertia matrix at float precision. **/
typedef ArticulatedInertia_<float>  fArticulatedInertia;
/** An articulated body inertia matrix at double precision. **/
typedef ArticulatedInertia_<double> dArticulatedInertia;

/** For backwards compatibility only; use UnitInertia instead. **/
typedef UnitInertia  Gyration;

// -----------------------------------------------------------------------------
//                             INERTIA MATRIX
// -----------------------------------------------------------------------------
/** The physical meaning of an inertia is the distribution of a rigid body's 
mass about a \e particular point. If that point is the center of mass of the 
body, then the measured inertia is called the "central inertia" of that body. 
To write down the inertia, we need to calculate the six scalars of the inertia 
tensor, which is a symmetric 3x3 matrix. These scalars must be expressed in 
an arbitrary but specified coordinate system. So an Inertia is meaningful only 
in conjunction with a particular set of axes, fixed to the body, whose origin 
is the point about which the inertia is being measured, and in whose 
coordinate system this measurement is being expressed. Note that changing the 
reference point results in a new physical quantity, but changing the reference 
axes only affects the measure numbers of that quantity. For any reference 
point, there is a unique set of reference axes in which the inertia tensor is
diagonal; those are called the "principal axes" of the body at that point, and 
the resulting diagonal elements are the "principal moments of inertia". When 
we speak of an inertia being "in" a frame, we mean the physical quantity 
measured about the frame's origin and then expressed in the frame's axes.

This low-level Inertia class does not attempt to keep track of \e which frame 
it is in. It provides construction and operations involving inertia that can 
proceed using only an implicit frame F. Clients of this class are responsible 
for keeping track of that frame. In particular, in order to shift the 
inertia's "measured-about" point one must know whether either the starting or 
final inertia is central, because we must always shift inertias by passing 
through the central inertia. So this class provides operations for doing the 
shifting, but expects to be told by the client where to find the center of mass.

Re-expressing an Inertia in a different coordinate system does not entail a 
change of physical meaning in the way that shifting it to a different point 
does. Note that because inertia is a tensor, there is a "left frame" and 
"right frame". For our purposes, these will always be the same so we'll only 
indicate the frame once, as in 'I_pt_frame'. This should be understood to mean
'frame_I_pt_frame' and re-expressing an Inertia requires both a left and right 
multiply by the rotation matrix. So I_OB_B is the inertia about body B's 
origin point OB, expressed in B, while I_OB_G is the same physical quantity 
but expressed in Ground (the latter is a component of the Spatial Inertia
which we usually want in the Ground frame). Frame conversion is done logically 
like this:
<pre>
   I_OB_G = R_GB * I_OB_B * R_BG  (R_BG=~R_GB)
</pre>
but we can save computation time by performing this as a single operation.

The central inertia would be I_CB_B for body B.

A Inertia is a symmetric matrix and is positive definite for nonsingular bodies
(that is, a body composed of at least three noncollinear point masses).

Some attempt is made to check the validity of an Inertia matrix, at least
when running in Debug mode. Some conditions it must satisfy are:
 - must be symmetric
 - all diagonal elements must be nonnegative
 - diagonal elements must satisfy the triangle inequality (sum of any two
   is greater than or equal the other one)

<h3>Abbreviations</h3>
Typedefs exist for the most common invocations of Inertia_\<P\>:
 - \ref SimTK::Inertia "Inertia" for default Real precision (this is 
   almost always used)
 - \ref SimTK::fInertia "fInertia" for single (float) precision
 - \ref SimTK::dInertia "dInertia" for double precision
**/
template <class P>
class SimTK_SimTKCOMMON_EXPORT Inertia_ {
    typedef Rotation_<P>    RotationP;
public:
/// Default is a NaN-ed out mess to avoid accidents, even in Release mode.
/// Other than this value, an Inertia matrix should always be valid.
Inertia_() : I_OF_F(NTraits<P>::getNaN()) {}

// Default copy constructor, copy assignment, destructor.

/// Create a principal inertia matrix with identical diagonal elements,
/// like a sphere where moment=2/5 m r^2, or a cube where 
/// moment=1/6 m s^2, with m the total mass, r the sphere's radius
/// and s the length of a side of the cube. Note that many rigid
/// bodies of different shapes and masses can have the same inertia
/// matrix.
explicit Inertia_(const P& moment) : I_OF_F(moment) 
{   errChk("Inertia::Inertia(moment)"); }

/// Create an Inertia matrix for a point mass at a given location,
/// measured from the origin OF of the implicit frame F, and expressed
/// in F. Cost is 14 flops.
Inertia_(const Vec<3,P>& p, const P& mass) : I_OF_F(pointMassAt(p,mass)) {}

/// Create an inertia matrix from a vector of the \e moments of
/// inertia (the inertia matrix diagonal) and optionally a vector of
/// the \e products of inertia (the off-diagonals). Moments are
/// in the order xx,yy,zz; products are xy,xz,yz.
explicit Inertia_(const Vec<3,P>& moments, const Vec<3,P>& products=Vec<3,P>(0)) 
{   I_OF_F.updDiag()  = moments;
    I_OF_F.updLower() = products;
    errChk("Inertia::Inertia(moments,products)"); }

/// Create a principal inertia matrix (only non-zero on diagonal).
Inertia_(const P& xx, const P& yy, const P& zz) 
{   I_OF_F = SymMat<3,P>(xx,
                        0, yy,
                        0,  0, zz);
    errChk("Inertia::setInertia(xx,yy,zz)"); }

/// This is a general inertia matrix. Note the order of these
/// arguments: moments of inertia first, then products of inertia.
Inertia_(const P& xx, const P& yy, const P& zz,
            const P& xy, const P& xz, const P& yz) 
{   I_OF_F = SymMat<3,P>(xx,
                        xy, yy,
                        xz, yz, zz);
    errChk("Inertia::setInertia(xx,yy,zz,xy,xz,yz)"); }

/// Construct an Inertia from a symmetric 3x3 matrix. The diagonals must
/// be nonnegative and satisfy the triangle inequality.
explicit Inertia_(const SymMat<3,P>& inertia) : I_OF_F(inertia)
{   errChk("Inertia::Inertia(SymMat33)"); }

/// Construct an Inertia matrix from a 3x3 symmetric matrix. In Debug mode
/// we'll test that the supplied matrix is numerically close to symmetric, and
/// that it satisfies other requirements of an Inertia matrix.
explicit Inertia_(const Mat<3,3,P>& m)
{   SimTK_ERRCHK(m.isNumericallySymmetric(), 
                    "Inertia(Mat33)", "The supplied matrix was not symmetric.");
    I_OF_F = SymMat<3,P>(m);
    errChk("Inertia(Mat33)"); }


/// Set an inertia matrix to have only principal moments (that is, it
/// will be diagonal). Returns a reference to "this" like an assignment operator.
Inertia_& setInertia(const P& xx, const P& yy, const P& zz) {
    I_OF_F = P(0); I_OF_F(0,0) = xx; I_OF_F(1,1) = yy;  I_OF_F(2,2) = zz;
    errChk("Inertia::setInertia(xx,yy,zz)");
    return *this;
}

/// Set principal moments and optionally off-diagonal terms.
/// Returns a reference to "this" like an assignment operator.
Inertia_& setInertia(const Vec<3,P>& moments, const Vec<3,P>& products=Vec<3,P>(0)) {
    I_OF_F.updDiag()  = moments;
    I_OF_F.updLower() = products;
    errChk("Inertia::setInertia(moments,products)");
    return *this;
}

/// Set this Inertia to a general matrix. Note the order of these
/// arguments: moments of inertia first, then products of inertia.
/// Behaves like an assignment statement. Will throw an error message
/// in Debug mode if the supplied elements do not constitute a valid
/// Inertia matrix.
Inertia_& setInertia(const P& xx, const P& yy, const P& zz,
                        const P& xy, const P& xz, const P& yz) {
    setInertia(Vec<3,P>(xx,yy,zz), Vec<3,P>(xy,xz,yz));
    errChk("Inertia::setInertia(xx,yy,zz,xy,xz,yz)");
    return *this;
}


/// Add in another inertia matrix. Frames and reference point must be the same but
/// we can't check. (6 flops)
Inertia_& operator+=(const Inertia_& inertia) 
{   I_OF_F += inertia.I_OF_F; 
    errChk("Inertia::operator+=()");
    return *this; }

/// Subtract off another inertia matrix. Frames and reference point must 
/// be the same but we can't check. (6 flops)
Inertia_& operator-=(const Inertia_& inertia) 
{   I_OF_F -= inertia.I_OF_F; 
    errChk("Inertia::operator-=()");
    return *this; }

/// Multiply this inertia matrix by a scalar. Cost is 6 flops.
Inertia_& operator*=(const P& s) {I_OF_F *= s; return *this;}

/// Divide this inertia matrix by a scalar. Cost is about 20 flops (a divide
/// and 6 multiplies).
Inertia_& operator/=(const P& s) {I_OF_F /= s; return *this;}

/// Assume that the current inertia is about the F frame's origin OF, and
/// expressed in F. Given the vector from OF to the body center of mass CF,
/// and the mass m of the body, we can shift the inertia to the center
/// of mass. This produces a new Inertia I' whose (implicit) frame F' is
/// aligned with F but has origin CF (an inertia like that is called a "central
/// inertia". I' = I - Icom where Icom is the inertia of a fictitious
/// point mass of mass m (that is, the same as the body mass) located at CF 
/// (measured in F) about OF. Cost is 20 flops.
/// @see shiftToMassCenterInPlace(), shiftFromMassCenter()
Inertia_ shiftToMassCenter(const Vec<3,P>& CF, const P& mass) const 
{   Inertia_ I(*this); I -= pointMassAt(CF, mass);
    I.errChk("Inertia::shiftToMassCenter()");
    return I; }

/// Assume that the current inertia is about the F frame's origin OF, and
/// expressed in F. Given the vector from OF to the body center of mass CF,
/// and the mass m of the body, we can shift the inertia to the center
/// of mass. This produces a new Inertia I' whose (implicit) frame F' is
/// aligned with F but has origin CF (an inertia like that is called a "central
/// inertia". I' = I - Icom where Icom is the inertia of a fictitious
/// point mass of mass m (that is, the same as the body mass) located at CF 
/// (measured in F) about OF. Cost is 20 flops.
/// @see shiftToMassCenter() if you want to leave this object unmolested.
/// @see shiftFromMassCenterInPlace()
Inertia_& shiftToMassCenterInPlace(const Vec<3,P>& CF, const P& mass) 
{   (*this) -= pointMassAt(CF, mass);
    errChk("Inertia::shiftToMassCenterInPlace()");
    return *this; }

/// Assuming that the current inertia I is a central inertia (that is, it is
/// inertia about the body center of mass CF), shift it to some other point p
/// measured from the center of mass. This produces a new inertia I' about
/// the point p given by I' = I + Ip where Ip is the inertia of a fictitious
/// point mass of mass mtot (the total body mass) located at p, about CF.
/// Cost is 20 flops.
/// @see shiftFromMassCenterInPlace(), shiftToMassCenter()
Inertia_ shiftFromMassCenter(const Vec<3,P>& p, const P& mass) const
{   Inertia_ I(*this); I += pointMassAt(p, mass);
    I.errChk("Inertia::shiftFromMassCenter()");
    return I; }

/// Assuming that the current inertia I is a central inertia (that is, it is
/// inertia about the body center of mass CF), shift it to some other point p
/// measured from the center of mass. This produces a new inertia I' about
/// the point p given by I' = I + Ip where Ip is the inertia of a fictitious
/// point mass of mass mtot (the total body mass) located at p, about CF.
/// Cost is 20 flops.
/// @see shiftFromMassCenter() if you want to leave this object unmolested.
/// @see shiftToMassCenterInPlace()
Inertia_& shiftFromMassCenterInPlace(const Vec<3,P>& p, const P& mass)
{   (*this) += pointMassAt(p, mass);
    errChk("Inertia::shiftFromMassCenterInPlace()");
    return *this; }

/// Return a new inertia matrix like this one but re-expressed in another 
/// frame (leaving the origin point unchanged). Call this inertia matrix
/// I_OF_F, that is, it is taken about the origin of some frame F, and 
/// expressed in F. We want to return I_OF_B, the same inertia matrix,
/// still taken about the origin of F, but expressed in the B frame, given
/// by I_OF_B=R_BF*I_OF_F*R_FB where R_FB is the rotation matrix giving
/// the orientation of frame B in F. This is handled here by a special
/// method of the Rotation class which rotates a symmetric tensor
/// at a cost of 57 flops.
/// @see reexpressInPlace()
Inertia_ reexpress(const Rotation_<P>& R_FB) const 
{   return Inertia_((~R_FB).reexpressSymMat33(I_OF_F)); }

/// Rexpress using an inverse rotation to avoid having to convert it.
/// @see rexpress(Rotation) for information
Inertia_ reexpress(const InverseRotation_<P>& R_FB) const 
{   return Inertia_((~R_FB).reexpressSymMat33(I_OF_F)); }

/// Re-express this inertia matrix in another frame, changing the object
/// in place; see reexpress() if you want to leave this object unmolested
/// and get a new one instead. Cost is 57 flops.
/// @see reexpress() if you want to leave this object unmolested.
Inertia_& reexpressInPlace(const Rotation_<P>& R_FB)
{   I_OF_F = (~R_FB).reexpressSymMat33(I_OF_F); return *this; }

/// Rexpress in place using an inverse rotation to avoid having to convert it.
/// @see rexpressInPlace(Rotation) for information
Inertia_& reexpressInPlace(const InverseRotation_<P>& R_FB)
{   I_OF_F = (~R_FB).reexpressSymMat33(I_OF_F); return *this; }

P trace() const {return I_OF_F.trace();}

/// This is an implicit conversion to a const SymMat33.
operator const SymMat<3,P>&() const {return I_OF_F;}

/// Obtain a reference to the underlying symmetric matrix type.
const SymMat<3,P>& asSymMat33() const {return I_OF_F;}

/// Expand the internal packed representation into a full 3x3 symmetric
/// matrix with all elements set.
Mat<3,3,P> toMat33() const {return Mat<3,3,P>(I_OF_F);}

/// Obtain the inertia moments (diagonal of the Inertia matrix) as a Vec3
/// ordered xx, yy, zz.
const Vec<3,P>& getMoments()  const {return I_OF_F.getDiag();}
/// Obtain the inertia products (off-diagonals of the Inertia matrix)
/// as a Vec3 with elements ordered xy, xz, yz.
const Vec<3,P>& getProducts() const {return I_OF_F.getLower();}

bool isNaN()    const {return I_OF_F.isNaN();}
bool isInf()    const {return I_OF_F.isInf();}
bool isFinite() const {return I_OF_F.isFinite();}

/// Compare this inertia matrix with another one and return true if they
/// are close to within a default numerical tolerance. Cost is about
/// 30 flops.
bool isNumericallyEqual(const Inertia_<P>& other) const 
{   return I_OF_F.isNumericallyEqual(other.I_OF_F); }

/// Compare this inertia matrix with another one and return true if they
/// are close to within a specified numerical tolerance. Cost is about
/// 30 flops.
bool isNumericallyEqual(const Inertia_<P>& other, double tol) const 
{   return I_OF_F.isNumericallyEqual(other.I_OF_F, tol); }

/// %Test some conditions that must hold for a valid Inertia matrix.
/// Cost is about 25 flops.
static bool isValidInertiaMatrix(const SymMat<3,P>& m) {
    if (m.isNaN()) return false;

    const Vec<3,P>& d = m.diag();
    if (!(d >= 0)) return false;    // diagonals must be nonnegative

    const P Slop = std::max(d.sum(),P(1))
                       * NTraits<P>::getSignificant();

    if (!(d[0]+d[1]+Slop>=d[2] && d[0]+d[2]+Slop>=d[1] && d[1]+d[2]+Slop>=d[0]))
        return false;               // must satisfy triangle inequality

    // Thanks to Paul Mitiguy for this condition on products of inertia.
    const Vec<3,P>& p = m.getLower();
    if (!( d[0]+Slop>=std::abs(2*p[2])
        && d[1]+Slop>=std::abs(2*p[1])
        && d[2]+Slop>=std::abs(2*p[0])))
        return false;               // max products are limited by moments

    return true;
}

/// Create an Inertia matrix for a point located at the origin -- that is,
/// an all-zero matrix.
static Inertia_ pointMassAtOrigin() {return Inertia_(0);}

/// Create an Inertia matrix for a point of a given mass, located at 
/// a given location measured from the origin of the implicit F frame.
/// This is equivalent to m*crossMatSq(p) but is implemented elementwise
/// here for speed, giving a cost of 14 flops.
static Inertia_ pointMassAt(const Vec<3,P>& p, const P& m) {
    const Vec<3,P> mp = m*p;       // 3 flops
    const P mxx = mp[0]*p[0];
    const P myy = mp[1]*p[1];
    const P mzz = mp[2]*p[2];
    const P nmx = -mp[0];
    const P nmy = -mp[1];
    return Inertia_( myy+mzz,  mxx+mzz,  mxx+myy,
                        nmx*p[1], nmx*p[2], nmy*p[2] );
}

/// @name Unit inertia matrix factories
/// These return UnitInertia matrices (inertias of unit-mass objects) 
/// converted to Inertias. Multiply the result by the actual mass
/// to get the Inertia of an actual object of this shape. See the 
/// UnitInertia class for more information.
//@{

/// Create a UnitInertia matrix for a unit mass sphere of radius \a r centered
/// at the origin.
inline static Inertia_ sphere(const P& r);

/// Unit-mass cylinder aligned along z axis;  use radius and half-length.
/// If r==0 this is a thin rod; hz=0 it is a thin disk.
inline static Inertia_ cylinderAlongZ(const P& r, const P& hz);

/// Unit-mass cylinder aligned along y axis;  use radius and half-length.
/// If r==0 this is a thin rod; hy=0 it is a thin disk.
inline static Inertia_ cylinderAlongY(const P& r, const P& hy);

/// Unit-mass cylinder aligned along x axis; use radius and half-length.
/// If r==0 this is a thin rod; hx=0 it is a thin disk.
inline static Inertia_ cylinderAlongX(const P& r, const P& hx);

/// Unit-mass brick given by half-lengths in each direction. One dimension zero
/// gives inertia of a thin rectangular sheet; two zero gives inertia
/// of a thin rod in the remaining direction.
inline static Inertia_ brick(const P& hx, const P& hy, const P& hz);

/// Alternate interface to brick() that takes a Vec3 for the half lengths.
inline static Inertia_ brick(const Vec<3,P>& halfLengths);

/// Unit-mass ellipsoid given by half-lengths in each direction.
inline static Inertia_ ellipsoid(const P& hx, const P& hy, const P& hz);

/// Alternate interface to ellipsoid() that takes a Vec3 for the half lengths.
inline static Inertia_ ellipsoid(const Vec<3,P>& halfLengths);

//@}

protected:
// Reinterpret this Inertia matrix as a UnitInertia matrix, that is, as the
// inertia of something with unit mass. This is useful in implementing
// methods of the UnitInertia class in terms of Inertia methods. Be sure you
// know that this is a unit-mass inertia!
const UnitInertia_<P>& getAsUnitInertia() const
{   return *static_cast<const UnitInertia_<P>*>(this); }
UnitInertia_<P>& updAsUnitInertia()
{   return *static_cast<UnitInertia_<P>*>(this); }

// If error checking is enabled (only in Debug mode), this 
// method will run some tests on the current contents of this Inertia 
// matrix and throw an error message if it is not valid. This should be 
// the same set of tests as run by the isValidInertiaMatrix() method above.
void errChk(const char* methodName) const {
#ifndef NDEBUG
    SimTK_ERRCHK(!isNaN(), methodName,
        "Inertia matrix contains a NaN.");

    const Vec<3,P>& d = I_OF_F.getDiag();  // moments
    const Vec<3,P>& p = I_OF_F.getLower(); // products
    const P Ixx = d[0], Iyy = d[1], Izz = d[2];
    const P Ixy = p[0], Ixz = p[1], Iyz = p[2];

    SimTK_ERRCHK3(d >= -SignificantReal, methodName,
        "Diagonals of an Inertia matrix must be nonnegative; got %g,%g,%g.",
        (double)Ixx,(double)Iyy,(double)Izz);

    // TODO: This is looser than it should be as a workaround for distorted
    // rotation matrices that were produced by an 11,000 body chain that
    // Sam Flores encountered. 
    const P Slop = std::max(d.sum(),P(1))
                       * std::sqrt(NTraits<P>::getEps());

    SimTK_ERRCHK3(   Ixx+Iyy+Slop>=Izz 
                  && Ixx+Izz+Slop>=Iyy 
                  && Iyy+Izz+Slop>=Ixx,
        methodName,
        "Diagonals of an Inertia matrix must satisfy the triangle "
        "inequality; got %g,%g,%g.",
        (double)Ixx,(double)Iyy,(double)Izz);

    // Thanks to Paul Mitiguy for this condition on products of inertia.
    SimTK_ERRCHK(   Ixx+Slop>=std::abs(2*Iyz) 
                 && Iyy+Slop>=std::abs(2*Ixz)
                 && Izz+Slop>=std::abs(2*Ixy),
        methodName,
        "The magnitude of a product of inertia was too large to be physical.");
#endif
}

// Inertia expressed in frame F and about F's origin OF. Note that frame F
// is implicit here; all we actually have are the inertia scalars.
SymMat<3,P> I_OF_F; 
};

/// Add two compatible inertia matrices, meaning they must be taken about the
/// same point and expressed in the same frame. There is no way to verify
/// compatibility; make sure you know what you're doing. Cost is 6 flops.
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator+(const Inertia_<P>& l, const Inertia_<P>& r) 
{   return Inertia_<P>(l) += r; }

/// Subtract from one inertia matrix another one which is compatible, meaning 
/// that both must be taken about the same point and expressed in the same frame. 
/// There is no way to verify compatibility; make sure you know what you're doing. 
/// Cost is 6 flops.
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator-(const Inertia_<P>& l, const Inertia_<P>& r) 
{   return Inertia_<P>(l) -= r; }

/// Multiply an inertia matrix by a scalar. Cost is 6 flops.
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator*(const Inertia_<P>& i, const P& r) 
{   return Inertia_<P>(i) *= r; }

/// Multiply an inertia matrix by a scalar. Cost is 6 flops.
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator*(const P& r, const Inertia_<P>& i) 
{   return Inertia_<P>(i) *= r; }


/// Multiply an inertia matrix by a scalar given as an int. 
/// Cost is 6 flops.
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator*(const Inertia_<P>& i, int r) 
{   return Inertia_<P>(i) *= P(r); }

/// Multiply an inertia matrix by a scalar given as an int. 
/// Cost is 6 flops.
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator*(int r, const Inertia_<P>& i) 
{   return Inertia_<P>(i) *= P(r); }

/// Divide an inertia matrix by a scalar. Cost is about 20
/// flops (one divide and six multiplies).
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator/(const Inertia_<P>& i, const P& r) 
{   return Inertia_<P>(i) /= r; }

/// Divide an inertia matrix by a scalar provided as an int. 
/// Cost is about 20 flops (one divide and six multiplies).
/// @relates Inertia_
template <class P> inline Inertia_<P>
operator/(const Inertia_<P>& i, int r) 
{   return Inertia_<P>(i) /= P(r); }

/// Multiply an inertia matrix I on the right by a vector w giving the
/// vector result I*w.
/// @relates Inertia_
template <class P> inline Vec<3,P>
operator*(const Inertia_<P>& I, const Vec<3,P>& w) 
{   return I.asSymMat33() * w; }

/// Compare two inertia matrices for exact (bitwise) equality. This is
/// too strict for most purposes; use Inertia::isNumericallyEqual() instead
/// to test for approximate equality. Cost here is 6 flops.
/// @relates Inertia_
template <class P> inline bool
operator==(const Inertia_<P>& i1, const Inertia_<P>& i2) 
{   return i1.asSymMat33() == i2.asSymMat33(); }

/// Output a human-readable representation of an inertia matrix to the 
/// indicated stream.
/// @relates Inertia_
template <class P> inline std::ostream& 
operator<<(std::ostream& o, const Inertia_<P>& inertia)
{   return o << inertia.toMat33(); }


// -----------------------------------------------------------------------------
//                            UNIT INERTIA MATRIX
// -----------------------------------------------------------------------------
/** A UnitInertia matrix is a unit-mass inertia matrix; you can convert it to an
Inertia by multiplying it by the actual body mass. Functionality is limited
here to those few operations which ensure unit mass; most operations on a
UnitInertia matrix result in a general Inertia instead. You can use a 
UnitInertia object wherever an Inertia is expected but not vice versa.

When constructing a UnitInertia matrix, note that we cannot verify that it 
actually has unit mass because every legal Inertia matrix can be viewed as
the UnitInertia matrix for some differently-scaled object.

Unit inertia matrices are sometimes called "gyration" matrices; we will often
represent them with the symbol "G" to avoid confusion with general inertia
matrices for which the symbol "I" (or sometimes "J") is used. 

<h3>Abbreviations</h3>
Typedefs exist for the most common invocations of UnitInertia_\<P\>:
 - \ref SimTK::UnitInertia "UnitInertia" for default Real precision (this is 
   almost always used)
 - \ref SimTK::fUnitInertia "fUnitInertia" for single (float) precision
 - \ref SimTK::dUnitInertia "dUnitInertia" for double precision **/
template <class P>
class SimTK_SimTKCOMMON_EXPORT UnitInertia_ : public Inertia_<P> {
    typedef P               RealP;
    typedef Vec<3,P>        Vec3P;
    typedef SymMat<3,P>     SymMat33P;
    typedef Mat<3,3,P>      Mat33P;
    typedef Rotation_<P>    RotationP;
    typedef Inertia_<P>     InertiaP;
public:
/// Default is a NaN-ed out mess to avoid accidents, even in Release mode.
/// Other than this value, a UnitInertia_ should always be valid.
UnitInertia_() {}

// Default copy constructor, copy assignment, destructor.

/// Create a principal unit inertia matrix with identical diagonal elements.
/// This is the unit inertia matrix of a unit mass sphere of radius 
/// r = sqrt(5/2 * moment) centered on the origin.
explicit UnitInertia_(const RealP& moment) : InertiaP(moment) {}

/// Create a unit inertia matrix from a vector of the \e moments of
/// inertia (the inertia matrix diagonal) and optionally a vector of
/// the \e products of inertia (the off-diagonals). Moments are
/// in the order xx,yy,zz; products are xy,xz,yz.
explicit UnitInertia_(const Vec3P& moments, const Vec3P& products=Vec3P(0))
:   InertiaP(moments,products) {}

/// Create a principal unit inertia matrix (only non-zero on diagonal).
UnitInertia_(const RealP& xx, const RealP& yy, const RealP& zz)
:   InertiaP(xx,yy,zz) {}   

/// This is a general unit inertia matrix. Note the order of these
/// arguments: moments of inertia first, then products of inertia.
UnitInertia_(const RealP& xx, const RealP& yy, const RealP& zz,
            const RealP& xy, const RealP& xz, const RealP& yz)
:   InertiaP(xx,yy,zz,xy,xz,yz) {}

/// Construct a UnitInertia from a symmetric 3x3 matrix. The diagonals must
/// be nonnegative and satisfy the triangle inequality.
explicit UnitInertia_(const SymMat33P& m) : InertiaP(m) {}

/// Construct a UnitInertia from a 3x3 symmetric matrix. In Debug mode
/// we'll test that the supplied matrix is numerically close to symmetric, 
/// and that it satisfies other requirements of an inertia matrix.
explicit UnitInertia_(const Mat33P& m) : InertiaP(m) {}

/// Construct a UnitInertia matrix from an Inertia matrix. Note that there
/// is no way to check whether this is really a unit inertia -- \e any
/// inertia matrix may be interpreted as a unit inertia for some shape. So
/// be sure you know what you're doing before you use this constructor!
explicit UnitInertia_(const Inertia_<P>& inertia) : InertiaP(inertia) {}

/// Set a UnitInertia matrix to have only principal moments (that is, it
/// will be diagonal). Returns a reference to "this" like an assignment 
/// operator.
UnitInertia_& setUnitInertia(const RealP& xx, const RealP& yy, const RealP& zz) 
{   InertiaP::setInertia(xx,yy,zz); return *this; }

/// Set principal moments and optionally off-diagonal terms.
/// Returns a reference to "this" like an assignment operator.
UnitInertia_& setUnitInertia(const Vec3P& moments, const Vec3P& products=Vec3P(0)) 
{   InertiaP::setInertia(moments,products); return *this; }

/// Set this UnitInertia to a general matrix. Note the order of these
/// arguments: moments of inertia first, then products of inertia.
/// Behaves like an assignment statement. Will throw an error message
/// in Debug mode if the supplied elements do not constitute a valid
/// inertia matrix.
UnitInertia_& setUnitInertia(const RealP& xx, const RealP& yy, const RealP& zz,
                        const RealP& xy, const RealP& xz, const RealP& yz) 
{   InertiaP::setInertia(xx,yy,zz,xy,xz,yz); return *this; }


// No +=, -=, etc. operators because those don't result in a UnitInertia 
// matrix. The parent class ones are suppressed below.

/// Assuming that this unit inertia matrix is currently taken about some (implicit)
/// frame F's origin OF, produce a new unit inertia matrix which is the same as this one
/// except measured about the body's centroid CF. We are given the vector from OF to 
/// the centroid CF, expressed in F. This produces a new UnitInertia matrix G' whose 
/// (implicit) frame F' is aligned with F but has origin CF (an inertia matrix like 
/// that is called "central" or "centroidal"). From the parallel axis theorem for
/// inertias, G' = G - Gcom where Gcom is the inertia matrix of a fictitious, 
/// unit-mass point located at CF (measured in F) taken about OF. (17 flops)
/// @see shiftToCentroidInPlace(), shiftFromCentroid()
UnitInertia_ shiftToCentroid(const Vec3P& CF) const 
{   UnitInertia_ G(*this); 
    G.Inertia_<P>::operator-=(pointMassAt(CF));
    return G; }

/// Assuming that this unit inertia matrix is currently taken about some (implicit)
/// frame F's origin OF, modify it so that it is instead taken about the body's 
/// centroid CF. We are given the vector from OF to 
/// the centroid CF, expressed in F. This produces a new UnitInertia G' whose 
/// (implicit) frame F' is aligned with F but has origin CF (an inertia matrix like 
/// that is called "central" or "centroidal"). From the parallel axis theorem for
/// inertias, G' = G - Gcom where Gcom is the inertia matrix of a fictitious, 
/// unit-mass point located at CF (measured in F) taken about OF. A reference
/// to the modified object is returned so that you can chain this method in
/// the manner of assignment operators. Cost is 17 flops.
/// @see shiftToCentroid() if you want to leave this object unmolested.
/// @see shiftFromCentroidInPlace()
UnitInertia_& shiftToCentroidInPlace(const Vec3P& CF) 
{   InertiaP::operator-=(pointMassAt(CF));
    return *this; }

/// Assuming that the current UnitInertia G is a central inertia (that is, it is
/// inertia about the body centroid CF), create a new object that is the same
/// as this one except shifted to some other point p measured from the centroid. 
/// This produces a new inertia G' about the point p given by G' = G + Gp where 
/// Gp is the inertia of a fictitious point located at p, taken about CF. Cost
/// is 17 flops.
/// @see shiftFromCentroidInPlace(), shiftToCentroid()
UnitInertia_ shiftFromCentroid(const Vec3P& p) const
{   UnitInertia_ G(*this); 
    G.Inertia_<P>::operator+=(pointMassAt(p));
    return G; }

/// Assuming that the current UnitInertia G is a central inertia (that is, it is
/// inertia about the body centroid CF), shift it in place to some other point p
/// measured from the centroid. This changes G to a modified inertia G' taken
/// about the point p, with the parallel axis theorem for inertia giving 
/// G' = G + Gp where Gp is the inertia of a fictitious, unit-mass point located 
/// at p, taken about CF. Cost is 17 flops.
/// @see shiftFromCentroid() if you want to leave this object unmolested.
/// @see shiftToCentroidInPlace()
UnitInertia_& shiftFromCentroidInPlace(const Vec3P& p)
{   InertiaP::operator+=(pointMassAt(p));
    return *this; }

/// Return a new unit inertia matrix like this one but re-expressed in another 
/// frame (leaving the origin point unchanged). Call this inertia matrix
/// G_OF_F, that is, it is taken about the origin of some frame F, and 
/// expressed in F. We want to return G_OF_B, the same unit inertia matrix,
/// still taken about the origin of F, but expressed in the B frame, given
/// by G_OF_B=R_BF*G_OF_F*R_FB where R_FB is the rotation matrix giving
/// the orientation of frame B in F. This is handled here by a special
/// method of the Rotation class which rotates a symmetric tensor
/// at a cost of 57 flops.
/// @see reexpressInPlace()
UnitInertia_ reexpress(const Rotation_<P>& R_FB) const 
{   return UnitInertia_((~R_FB).reexpressSymMat33(this->I_OF_F)); }

/// Rexpress using an inverse rotation to avoid having to convert it.
/// @see rexpress(Rotation) for information
UnitInertia_ reexpress(const InverseRotation_<P>& R_FB) const 
{   return UnitInertia_((~R_FB).reexpressSymMat33(this->I_OF_F)); }

/// Re-express this unit inertia matrix in another frame, changing the object
/// in place; see reexpress() if you want to leave this object unmolested
/// and get a new one instead. Cost is 57 flops.
/// @see reexpress() if you want to leave this object unmolested.
UnitInertia_& reexpressInPlace(const Rotation_<P>& R_FB)
{   InertiaP::reexpressInPlace(R_FB); return *this; }

/// Rexpress using an inverse rotation to avoid having to convert it.
/// @see rexpressInPlace(Rotation) for information
UnitInertia_& reexpressInPlace(const InverseRotation_<P>& R_FB)
{   InertiaP::reexpressInPlace(R_FB); return *this; }


/// This is an implicit conversion to const SymMat33.
operator const SymMat33P&() const {return this->I_OF_F;}

/// Recast this UnitInertia matrix as a unit inertia matrix. This is just for
/// emphasis; a UnitInertia matrix is already a kind of Inertia matrix by
/// inheritance.
const Inertia_<P>& asUnitInertia() const 
{   return *static_cast<const Inertia_<P>*>(this); }

/// Set from a unit inertia matrix. Note that we can't check; every Inertia
/// matrix can be interpreted as a unit inertia for some shape.
UnitInertia_& setFromUnitInertia(const Inertia_<P>& inertia)
{   Inertia_<P>::operator=(inertia);
    return *this; }

/// %Test some conditions that must hold for a valid UnitInertia matrix.
/// Cost is about 9 flops.
/// TODO: this may not be comprehensive.
static bool isValidUnitInertiaMatrix(const SymMat33P& m) 
{   return Inertia_<P>::isValidInertiaMatrix(m); }

/// @name UnitInertia matrix factories
/// These are UnitInertia matrix factories for some common 3D solids. Each 
/// defines its own frame aligned (when possible) with principal moments. 
/// Each has unit mass and its center of mass located at the origin (usually). 
/// Use this with shiftFromCentroid() to move it somewhere else, and with 
/// reexpress() to express the UnitInertia matrix in another frame.
//@{

/// Create a UnitInertia matrix for a point located at the origin -- that is,
/// an all-zero matrix.
static UnitInertia_ pointMassAtOrigin() {return UnitInertia_(0);}

/// Create a UnitInertia matrix for a point of unit mass located at a given 
/// location measured from origin OF and expressed in F (where F is the
/// implicit frame of this UnitInertia matrix).
/// Cost is 11 flops.
static UnitInertia_ pointMassAt(const Vec3P& p) 
{   return UnitInertia_(crossMatSq(p)); }

/// Create a UnitInertia matrix for a unit mass sphere of radius \a r centered
/// at the origin.
static UnitInertia_ sphere(const RealP& r) {return UnitInertia_(RealP(0.4)*r*r);}

/// Unit-mass cylinder aligned along z axis;  use radius and half-length.
/// If r==0 this is a thin rod; hz=0 it is a thin disk.
static UnitInertia_ cylinderAlongZ(const RealP& r, const RealP& hz) {
    const RealP Ixx = (r*r)/4 + (hz*hz)/3;
    return UnitInertia_(Ixx,Ixx,(r*r)/2);
}

/// Unit-mass cylinder aligned along y axis;  use radius and half-length.
/// If r==0 this is a thin rod; hy=0 it is a thin disk.
static UnitInertia_ cylinderAlongY(const RealP& r, const RealP& hy) {
    const RealP Ixx = (r*r)/4 + (hy*hy)/3;
    return UnitInertia_(Ixx,(r*r)/2,Ixx);
}

/// Unit-mass cylinder aligned along x axis; use radius and half-length.
/// If r==0 this is a thin rod; hx=0 it is a thin disk.
static UnitInertia_ cylinderAlongX(const RealP& r, const RealP& hx) {
    const RealP Iyy = (r*r)/4 + (hx*hx)/3;
    return UnitInertia_((r*r)/2,Iyy,Iyy);
}

/// Unit-mass brick given by half-lengths in each direction. One dimension zero
/// gives inertia of a thin rectangular sheet; two zero gives inertia
/// of a thin rod in the remaining direction.
static UnitInertia_ brick(const RealP& hx, const RealP& hy, const RealP& hz) {
    const RealP oo3 = RealP(1)/RealP(3);
    const RealP hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
    return UnitInertia_(oo3*(hy2+hz2), oo3*(hx2+hz2), oo3*(hx2+hy2));
}

/// Alternate interface to brick() that takes a Vec3 for the half lengths.
static UnitInertia_ brick(const Vec3P& halfLengths)
{   return brick(halfLengths[0],halfLengths[1],halfLengths[2]); }

/// Unit-mass ellipsoid given by half-lengths in each direction.
static UnitInertia_ ellipsoid(const RealP& hx, const RealP& hy, const RealP& hz) {
    const RealP hx2=hx*hx, hy2=hy*hy, hz2=hz*hz;
    return UnitInertia_((hy2+hz2)/5, (hx2+hz2)/5, (hx2+hy2)/5);
}

/// Alternate interface to ellipsoid() that takes a Vec3 for the half lengths.
static UnitInertia_ ellipsoid(const Vec3P& halfLengths)
{   return ellipsoid(halfLengths[0],halfLengths[1],halfLengths[2]); }

//@}
private:
// Suppress Inertia_ methods which are not allowed for UnitInertia_.

// These kill all flavors of Inertia_::setInertia() and the
// Inertia_ assignment ops.
void setInertia() {}
void operator+=(int) {}
void operator-=(int) {}
void operator*=(int) {}
void operator/=(int) {}
};

// Implement Inertia methods which are pass-throughs to UnitInertia methods.

template <class P> inline Inertia_<P> Inertia_<P>::
sphere(const P& r) 
{   return UnitInertia_<P>::sphere(r); }
template <class P> inline Inertia_<P> Inertia_<P>::
cylinderAlongZ(const P& r, const P& hz)
{   return UnitInertia_<P>::cylinderAlongZ(r,hz); }
template <class P> inline Inertia_<P> Inertia_<P>::
cylinderAlongY(const P& r, const P& hy)
{   return UnitInertia_<P>::cylinderAlongY(r,hy); }
template <class P> inline Inertia_<P> Inertia_<P>::
cylinderAlongX(const P& r, const P& hx)
{   return UnitInertia_<P>::cylinderAlongX(r,hx); }
template <class P> inline Inertia_<P> Inertia_<P>::
brick(const P& hx, const P& hy, const P& hz)
{   return UnitInertia_<P>::brick(hx,hy,hz); }
template <class P> inline Inertia_<P> Inertia_<P>::
brick(const Vec<3,P>& halfLengths)
{   return UnitInertia_<P>::brick(halfLengths); }
template <class P> inline Inertia_<P> Inertia_<P>::
ellipsoid(const P& hx, const P& hy, const P& hz)
{   return UnitInertia_<P>::ellipsoid(hx,hy,hz); }
template <class P> inline Inertia_<P> Inertia_<P>::
ellipsoid(const Vec<3,P>& halfLengths)
{   return UnitInertia_<P>::ellipsoid(halfLengths); }


// -----------------------------------------------------------------------------
//                           SPATIAL INERTIA MATRIX
// -----------------------------------------------------------------------------
/** A spatial inertia contains the mass, center of mass point, and inertia
matrix for a rigid body. This is 10 independent quantities altogether; however,
inertia is mass-scaled making it linearly dependent on the mass. Here instead 
we represent inertia using a unit inertia matrix, which is equivalent to the 
inertia this body would have if it had unit mass. Then the actual inertia is 
given by mass*unitInertia. In this manner the mass, center of mass location, and 
inertia are completely independent so can be changed separately. That means 
if you double the mass, you'll also double the inertia as you would expect.

Spatial inertia may be usefully viewed as a symmetric spatial matrix, that is, 
a 6x6 symmetric matrix arranged as 2x2 blocks of 3x3 matrices. Although this 
class represents the spatial inertia in compact form, it supports methods and
operators that allow it to behave as though it were a spatial matrix (except
much faster to work with). In spatial matrix form, the matrix has the following 
interpretation:
<pre>
              [  m*G   m*px ]
          M = [             ]
              [ -m*px  m*I  ]
</pre>
Here m is mass, p is the vector from the body origin to the center of mass, 
G is the 3x3 symmetric unit inertia (gyration) matrix, and I is a 3x3 identity 
matrix. "px" indicates the skew symmetric cross product matrix formed from the 
vector p, so -px=~px. 

<h3>Abbreviations</h3>
Typedefs exist for the most common invocations of SpatialInertia_\<P\>:
 - \ref SimTK::SpatialInertia "SpatialInertia" for default Real precision (this is 
   almost always used)
 - \ref SimTK::fSpatialInertia "fSpatialInertia" for single (float) precision
 - \ref SimTK::dSpatialInertia "dSpatialInertia" for double precision **/
template <class P> 
class SimTK_SimTKCOMMON_EXPORT SpatialInertia_ {
    typedef P               RealP;
    typedef Vec<3,P>        Vec3P;
    typedef UnitInertia_<P> UnitInertiaP;
    typedef Mat<3,3,P>      Mat33P;
    typedef Vec<2, Vec3P>   SpatialVecP;
    typedef Mat<2,2,Mat33P> SpatialMatP;
    typedef Rotation_<P>    RotationP;
    typedef Transform_<P>   TransformP;
    typedef Inertia_<P>     InertiaP;
public:
/// The default constructor fills everything with NaN, even in Release mode.
SpatialInertia_() 
:   m(nanP()), p(nanP()) {} // inertia is already NaN
SpatialInertia_(RealP mass, const Vec3P& com, const UnitInertiaP& gyration) 
:   m(mass), p(com), G(gyration) {}

// default copy constructor, copy assignment, destructor

SpatialInertia_& setMass(RealP mass)
{   SimTK_ERRCHK1(mass >= 0, "SpatialInertia::setMass()",
        "Negative mass %g is illegal.", (double)mass);
    m=mass; return *this; }
SpatialInertia_& setMassCenter(const Vec3P& com)
{   p=com; return *this;} 
SpatialInertia_& setUnitInertia(const UnitInertiaP& gyration) 
{   G=gyration; return *this; }

RealP               getMass()        const {return m;}
const Vec3P&        getMassCenter()  const {return p;}
const UnitInertiaP& getUnitInertia() const {return G;}

/// Calculate the first mass moment (mass-weighted COM location)
/// from the mass and COM vector. Cost is 3 inline flops.
Vec3P calcMassMoment() const {return m*p;}

/// Calculate the inertia matrix (second mass moment, mass-weighted gyration
/// matrix) from the mass and unit inertia matrix. Cost is 6 inline flops.
InertiaP calcInertia() const {return m*G;}

/// Add in a compatible SpatialInertia. This is only valid if both 
/// SpatialInertias are expressed in the same frame and measured about 
/// the same point but there is no way for this method to check.
/// Cost is about 40 flops.
SpatialInertia_& operator+=(const SpatialInertia_& src) {
    SimTK_ERRCHK(m+src.m != 0, "SpatialInertia::operator+=()",
        "The combined mass cannot be zero.");
    const RealP mtot = m+src.m, oomtot = 1/mtot;                    // ~11 flops
    p = oomtot*(calcMassMoment() + src.calcMassMoment());           // 10 flops
    G.setFromUnitInertia(oomtot*(calcInertia()+src.calcInertia())); // 19 flops
    m = mtot; // must do this last
    return *this;
}

/// Subtract off a compatible SpatialInertia. This is only valid if both 
/// SpatialInertias are expressed in the same frame and measured about 
/// the same point but there is no way for this method to check.
/// Cost is about 40 flops.
SpatialInertia_& operator-=(const SpatialInertia_& src) {
    SimTK_ERRCHK(m != src.m, "SpatialInertia::operator-=()",
        "The combined mass cannot be zero.");
    const RealP mtot = m-src.m, oomtot = 1/mtot;                    // ~11 flops
    p = oomtot*(calcMassMoment() - src.calcMassMoment());           // 10 flops
    G.setFromUnitInertia(oomtot*(calcInertia()-src.calcInertia())); // 19 flops
    m = mtot; // must do this last
    return *this;
}

/// Multiply a SpatialInertia by a scalar. Because we keep the mass
/// factored out, this requires only a single multiply.
SpatialInertia_& operator*=(const RealP& s) {m *= s; return *this;}

/// Divide a SpatialInertia by a scalar. Because we keep the mass
/// factored out, this requires only a single divide.
SpatialInertia_& operator/=(const RealP& s) {m /= s; return *this;}

/// Multiply a SpatialInertia by a SpatialVec to produce a SpatialVec
/// result; 45 flops.
SpatialVecP operator*(const SpatialVecP& v) const
{   return m*SpatialVecP(G*v[0]+p%v[1], v[1]-p%v[0]); }

/// Return a new SpatialInertia object which is the same as this one except
/// re-expressed in another coordinate frame. We consider this object to
/// be expressed in some frame F and we're given a rotation matrix R_FB we
/// can use to re-express in a new frame B. Cost is 72 flops.
/// @see reexpressInPlace()
SpatialInertia_ reexpress(const Rotation_<P>& R_FB) const
{   return SpatialInertia_(*this).reexpressInPlace(R_FB); }

/// Rexpress using an inverse rotation to avoid having to convert it.
/// @see rexpress(Rotation) for information
SpatialInertia_ reexpress(const InverseRotation_<P>& R_FB) const
{   return SpatialInertia_(*this).reexpressInPlace(R_FB); }

/// Re-express this SpatialInertia in another frame, modifying the original
/// object. We return a reference to the object so that you can chain this
/// operation in the manner of assignment operators. Cost is 72 flops.
/// @see reexpress() if you want to leave this object unmolested.
SpatialInertia_& reexpressInPlace(const Rotation_<P>& R_FB)
{   p = (~R_FB)*p; G.reexpressInPlace(R_FB); return *this; }

/// Rexpress using an inverse rotation to avoid having to convert it.
/// @see rexpressInPlace(Rotation) for information
SpatialInertia_& reexpressInPlace(const InverseRotation_<P>& R_FB)
{   p = (~R_FB)*p; G.reexpressInPlace(R_FB); return *this; }

/// Return a new SpatialInertia object which is the same as this one except
/// the origin ("taken about" point) has changed from OF to OF+S.
/// Cost is 37 flops.
/// @see shiftInPlace()
SpatialInertia_ shift(const Vec3P& S) const 
{   return SpatialInertia_(*this).shiftInPlace(S); }

/// Change origin from OF to OF+S, modifying the original object in place.
/// Returns a reference to the modified object so that you can chain this
/// operation in the manner of assignment operators. Cost is 37 flops.
/// @see shift() if you want to leave this object unmolested.
SpatialInertia_& shiftInPlace(const Vec3P& S) {
    G.shiftToCentroidInPlace(p);      // change to central inertia
    const Vec3P pNew(p-S);            // vec from new origin to com (3 flops)
    G.shiftFromCentroidInPlace(pNew); // shift to S (pNew sign doesn't matter)
    p = pNew;                         // had p=com-OF; want p=com-(OF+S)=p-S
    return *this;
}

/// Return a new SpatialInertia object which is the same as this
/// one but measured about and expressed in a new frame. We consider
/// the current spatial inertia M to be measured (implicitly) in some
/// frame F, that is, we have M=M_OF_F. We want M_OB_B for some new
/// frame B, given the transform X_FB giving the location and orientation
/// of B in F. This combines the reexpress() and shift() operations
/// available separately. Cost is 109 flops.
/// @see transformInPlace()
SpatialInertia_ transform(const Transform_<P>& X_FB) const 
{   return SpatialInertia_(*this).transformInPlace(X_FB); }

/// Transform using an inverse transform to avoid having to convert it.
/// @see transform(Transform) for information
SpatialInertia_ transform(const InverseTransform_<P>& X_FB) const 
{   return SpatialInertia_(*this).transformInPlace(X_FB); }

/// Transform this SpatialInertia object so that it is measured about and
/// expressed in a new frame, modifying the object in place. We consider the
/// current spatial inertia M to be measured (implicitly) in some frame F, that 
/// is, we have M=M_OF_F. We want to change it to M_OB_B for some new frame B, 
/// given the transform X_FB giving the location and orientation of B in F. This 
/// combines the reexpressInPlace() and shiftInPlace() operations available 
/// separately. Returns a reference to the modified object so that you can
/// chain this operation in the manner of assignment operators. Cost is 109 flops.
/// @see transform() if you want to leave this object unmolested.
SpatialInertia_& transformInPlace(const Transform_<P>& X_FB) {
    shiftInPlace(X_FB.p());     // shift to the new origin OB.
    reexpressInPlace(X_FB.R()); // get everything in B
    return *this;
}

/// Transform using an inverse transform to avoid having to convert it.
/// @see transformInPlace(Transform) for information
SpatialInertia_& transformInPlace(const InverseTransform_<P>& X_FB) {
    shiftInPlace(X_FB.p());     // shift to the new origin OB.
    reexpressInPlace(X_FB.R()); // get everything in B
    return *this;
}

const SpatialMatP toSpatialMat() const {
    Mat33P offDiag = crossMat(m*p);
    return SpatialMatP(m*G.toMat33(), offDiag,
                       -offDiag,      Mat33P(m));
}

private:
RealP           m;  // mass of this rigid body F
Vec3P           p;  // location of body's COM from OF, expressed in F
UnitInertiaP    G;  // mass distribution; inertia is mass*gyration

static P nanP() {return NTraits<P>::getNaN();} 
};

/// Add two compatible spatial inertias. Cost is about 40 flops.
/// @relates SpatialInertia_
template <class P> inline SpatialInertia_<P> 
operator+(const SpatialInertia_<P>& l, const SpatialInertia_<P>& r)
{   return SpatialInertia_<P>(l) += r; } 

/// Subtract one compatible spatial inertia from another. Cost is
/// about 40 flops.
/// @relates SpatialInertia_
template <class P> inline SpatialInertia_<P> 
operator-(const SpatialInertia_<P>& l, const SpatialInertia_<P>& r)
{   return SpatialInertia_<P>(l) -= r; } 


// -----------------------------------------------------------------------------
//                        ARTICULATED BODY INERTIA MATRIX
// -----------------------------------------------------------------------------
/** An articulated body inertia (ABI) matrix P(q) contains the spatial inertia 
properties that a body appears to have when it is the free base body of 
an articulated multibody tree in a given configuration q. Despite the 
complex relative motion that occurs within a multibody tree, at any given 
configuration q there is still a linear relationship between a spatial 
force F applied to a point of the base body and the resulting acceleration 
A of that body and that point: F = P(q)*A + c, where c is a velocity-
dependent inertial bias force. P is thus analogous to a rigid body 
spatial inertia (RBI), but for a body which has other bodies connected to it
by joints which are free to move.

An ABI P is a symmetric 6x6 spatial matrix, consisting of 2x2 blocks of 3x3 
matrices, similar to the RBI. However, unlike the RBI which has only 10 independent
elements, all 21 elements of P's lower triangle are significant. For example,
the apparent mass of an articulated body depends on which way you push it,
and in general there is no well-defined center of mass. This
is a much more expensive matrix to manipulate than an RBI. In Simbody's formulation,
we only work with ABIs in the Ground frame, so there
is never a need to rotate or re-express them. (That is done by rotating RBIs
prior to using them to construct the ABIs.) Thus only shifting operations need
be performed when transforming ABIs from body to body.
Cheap rigid body shifting is done when moving an ABI 
within a body or across a prescribed mobilizer; otherwise we have to perform 
an articulated shift operation which is quite expensive.

For a full discussion of the properties of articulated body inertias, see 
Section 7.1 (pp. 119-123) of Roy Featherstone's excellent 2008 book, Rigid 
Body Dynamics Algorithms. 

In spatial matrix form, an ABI P may be considered to consist of the 
following 3x3 subblocks:
<pre>
         P =  [ J  F ]
              [~F  M ]
</pre>
Here M is a (symmetric) mass distribution, F is a full matrix giving the
first mass moment distribution, and J is a (symmetric) inertia matrix. 

<h3>Abbreviations</h3>
Typedefs exist for the most common invocations of ArticulatedInertia_\<P\>:
 - \ref SimTK::ArticulatedInertia "ArticulatedInertia" for default Real 
   precision (this is almost always used)
 - \ref SimTK::fArticulatedInertia "fArticulatedInertia" for single (float) 
   precision
 - \ref SimTK::dArticulatedInertia "dArticulatedInertia" for double 
   precision **/
template <class P> 
class ArticulatedInertia_ {
    typedef P               RealP;
    typedef Vec<3,P>        Vec3P;
    typedef UnitInertia_<P> UnitInertiaP;
    typedef Mat<3,3,P>      Mat33P;
    typedef SymMat<3,P>     SymMat33P;
    typedef Vec<2, Vec3P>   SpatialVecP;
    typedef Mat<2,2,Mat33P> SpatialMatP;
    typedef Rotation_<P>    RotationP;
    typedef Transform_<P>   TransformP;
    typedef Inertia_<P>     InertiaP;
public:
/// Default construction produces uninitialized junk at zero cost; be sure to 
/// fill this in before referencing it.
ArticulatedInertia_() {}
/// Construct an ArticulatedInertia from the mass, first moment, and inertia matrices it contains.
ArticulatedInertia_(const SymMat33P& mass, const Mat33P& massMoment, const SymMat33P& inertia)
:   M(mass), J(inertia), F(massMoment) {}

/// Construct an articulated body inertia (ABI) from a rigid body spatial inertia (RBI). 
/// Every RBI is also the ABI for that (unarticulated) rigid body. 12 flops.
explicit ArticulatedInertia_(const SpatialInertia_<P>& rbi)
:   M(rbi.getMass()), J(rbi.calcInertia()), F(crossMat(rbi.calcMassMoment())) {}

/// Set the mass distribution matrix M in this ArticulatedInertia (symmetric).
ArticulatedInertia_& setMass      (const SymMat33P& mass)       {M=mass;       return *this;}
/// Set the mass first moment distribution matrix F in this ArticulatedInertia (full).
ArticulatedInertia_& setMassMoment(const Mat33P&    massMoment) {F=massMoment; return *this;}
/// Set the mass second moment (inertia) matrix J in this ArticulatedInertia (symmetric).
ArticulatedInertia_& setInertia   (const SymMat33P& inertia)    {J=inertia;    return *this;}

/// Get the mass distribution matrix M from this ArticulatedInertia (symmetric).
const SymMat33P& getMass()       const {return M;}
/// Get the mass first moment distribution matrix F from this ArticulatedInertia (full).
const Mat33P&    getMassMoment() const {return F;}
/// Get the mass second moment (inertia) matrix J from this ArticulatedInertia (symmetric).
const SymMat33P& getInertia()    const {return J;}

// default destructor, copy constructor, copy assignment

/// Add in a compatible ArticulatedInertia to this one. Both inertias must be expressed
/// in the same frame and measured about the same point. 21 flops.
ArticulatedInertia_& operator+=(const ArticulatedInertia_& src)
{   M+=src.M; J+=src.J; F+=src.F; return *this; }
/// Subtract a compatible ArticulatedInertia from this one. Both inertias must be expressed
/// in the same frame and measured about the same point. 21 flops.
ArticulatedInertia_& operator-=(const ArticulatedInertia_& src)
{   M-=src.M; J-=src.J; F-=src.F; return *this; }

/// Multiply an ArticulatedIneria by a SpatialVec (66 flops).
SpatialVecP operator*(const SpatialVecP& v) const
{   return SpatialVecP(J*v[0]+F*v[1], ~F*v[0]+M*v[1]); }

/// Multiply an ArticulatedInertia by a row of SpatialVecs (66*N flops).
template <int N>
Mat<2,N,Vec3P> operator*(const Mat<2,N,Vec3P>& m) const {
    Mat<2,N,Vec3P> res;
    for (int j=0; j < N; ++j)
        res.col(j) = (*this) * m.col(j); // above this*SpatialVec operator
    return res;
}

/// Rigid-shift the origin of this Articulated Body Inertia P by a 
/// shift vector -s to produce a new ABI P'. The calculation is 
/// <pre>
/// P' =  [ J'  F' ]  =  [ 1  sx ] [ J  F ] [ 1  0 ]
///       [~F'  M  ]     [ 0  1  ] [~F  M ] [-sx 1 ]
/// </pre>
/// where sx is the cross product matrix of s (not -s). Note that for historical
/// reasons this is defined so that it shifts by -s, which is unfortunate and
/// is opposite the SpatialInertia::shift() method.
/// Cost is 72 flops.
SimTK_SimTKCOMMON_EXPORT ArticulatedInertia_ shift(const Vec3P& s) const;

/// Rigid-shift this ABI in place by -s. 72 flops.
/// @see shift() for details
SimTK_SimTKCOMMON_EXPORT ArticulatedInertia_& shiftInPlace(const Vec3P& s);

/// Convert the compactly-stored ArticulatedInertia (21 elements) into a 
/// full SpatialMat with 36 elements.
const SpatialMatP toSpatialMat() const {
    // This more-beautiful code ran into an optimization bug Microsoft 
    // introduced in Update 2 to Visual C++ 2015.
    //return SpatialMatP( Mat33P(J),     F,
    //                        ~F,       Mat33P(M) );

    // Uglier but functional.
    SpatialMatP smat;
    smat(0,0) = Mat33P(J); smat(0,1) = F;
    smat(1,0) = ~F;        smat(1,1) = Mat33P(M);
    return smat;
}
private:
SymMat33P M;
SymMat33P J;
Mat33P    F;
};

/// Add two compatible articulated inertias. Cost is 21 flops.
/// @relates ArticulatedInertia_
template <class P> inline ArticulatedInertia_<P>
operator+(const ArticulatedInertia_<P>& l, const ArticulatedInertia_<P>& r)
{   return ArticulatedInertia_<P>(l) += r; }

/// Subtract one compatible articulated inertia from another. Cost is
/// 21 flops.
/// @relates ArticulatedInertia_
template <class P> inline ArticulatedInertia_<P>
operator-(const ArticulatedInertia_<P>& l, const ArticulatedInertia_<P>& r)
{   return ArticulatedInertia_<P>(l) -= r; }


// -----------------------------------------------------------------------------
//                              MASS PROPERTIES
// -----------------------------------------------------------------------------
/** This class contains the mass, center of mass, and unit inertia matrix of a
rigid body B. The center of mass is a vector from B's origin, expressed in the
B frame. The inertia is taken about the B origin, and expressed in B. The frame
B is implicit; only the measurements are stored here. The unit inertia can be
multiplied by the mass to get the inertia matrix for the body.

<h3>Abbreviations</h3>
Typedefs exist for the most common invocations of MassProperties_\<P\>:
 - \ref SimTK::MassProperties "MassProperties" for default Real precision (this is 
   almost always used)
 - \ref SimTK::fMassProperties "fMassProperties" for single (float) precision
 - \ref SimTK::dMassProperties "dMassProperties" for double precision **/
template <class P>
class SimTK_SimTKCOMMON_EXPORT MassProperties_ {
public:
/** Create a mass properties object in which the mass, mass center, and 
inertia are meaningless; you must assign values before using this. **/
MassProperties_() { setMassProperties(0,Vec<3,P>(0),UnitInertia_<P>()); }
/** Create a mass properties object from individually supplied mass,
mass center, and inertia matrix.  The inertia matrix is divided by the
mass to produce the unit inertia. **/
MassProperties_(const P& m, const Vec<3,P>& com, const Inertia_<P>& inertia)
    { setMassProperties(m,com,inertia); }
/** Create a mass properties object from individually supplied mass,
mass center, and unit inertia (gyration) matrix. **/
MassProperties_(const P& m, const Vec<3,P>& com, const UnitInertia_<P>& gyration)
    { setMassProperties(m,com,gyration); }

/** Set mass, center of mass, and inertia. The inertia is divided by the mass to
produce the unit inertia. Behaves like an assignment in that
a reference to the modified MassProperties object is returned. **/
MassProperties_& setMassProperties(const P& m, const Vec<3,P>& com, const Inertia_<P>& inertia) {
    mass = m;
    comInB = com;
    if (m == 0) {
        SimTK_ASSERT(inertia.asSymMat33().getAsVec() == 0, "Mass is zero but inertia is nonzero");
        unitInertia_OB_B = UnitInertia_<P>(0);
    }
    else {
        unitInertia_OB_B = UnitInertia_<P>(inertia*(1/m));
    }
    return *this;
}

/** Set mass, center of mass, and unit inertia. Behaves like an assignment in
that a reference to the modified MassProperties object is returned. **/
MassProperties_& setMassProperties
   (const P& m, const Vec<3,P>& com, const UnitInertia_<P>& gyration)
{   mass=m; comInB=com; unitInertia_OB_B=gyration; return *this; }

/** Return the mass currently stored in this MassProperties object. **/
const P& getMass() const {return mass;}
/** Return the mass center currently stored in this MassProperties object;
this is expressed in an implicit frame we call "B", and measured from B's
origin, but you have to know what that frame is in order to interpret the 
returned vector. **/
const Vec<3,P>& getMassCenter() const {return comInB;}
/** Return the unit inertia currently stored in this MassProperties object;
this is expressed in an implicit frame we call "B", and measured about B's
origin, but you have to know what that frame is in order to interpret the 
returned value correctly. **/
const UnitInertia_<P>& getUnitInertia() const {return unitInertia_OB_B;}
/** Return the inertia matrix for this MassProperties object; this is equal
to the unit inertia times the mass and costs 6 flops. It is expressed in 
an implicit frame we call "B", and measured about B's origin, but you have 
to know what that frame is in order to interpret the returned value
correctly. **/
const Inertia_<P> calcInertia() const {return mass*unitInertia_OB_B;}
/** OBSOLETE -- this is just here for compatibility with 2.2; since the
UnitInertia is stored now the full inertia must be calculated at a cost of
6 flops, hence the method name should be calcInertia() and that's what
you should use unless you want getUnitInertia(). **/
const Inertia_<P> getInertia() const {return calcInertia();}

/** Return the inertia of this MassProperties object, but measured about the
mass center rather than about the (implicit) B frame origin. The result is
still expressed in B. **/
Inertia_<P> calcCentralInertia() const {
    return mass*unitInertia_OB_B - Inertia_<P>(comInB, mass);
}
/** Return the inertia of this MassProperties object, but with the "measured
about" point shifted from the (implicit) B frame origin to a new point that
is supplied in \a newOriginB which must be a vector measured from the B frame
origin and expressed in B. The result is still expressed in B. **/
Inertia_<P> calcShiftedInertia(const Vec<3,P>& newOriginB) const {
    return calcCentralInertia() + Inertia_<P>(newOriginB-comInB, mass);
}
/** Return the inertia of this MassProperties object, but transformed to
from the implicit B frame to a new frame C whose pose relative to B is
supplied. Note that this affects both the "measured about" point and the
"expressed in" frame. **/
Inertia_<P> calcTransformedInertia(const Transform_<P>& X_BC) const {
    return calcShiftedInertia(X_BC.p()).reexpress(X_BC.R());
}
/** Return a new MassProperties object that is the same as this one but with 
the origin point shifted from the (implicit) B frame origin to a new point that
is supplied in \a newOriginB which must be a vector measured from the B frame
origin and expressed in B. This affects both the mass center vector and the
inertia. The result is still expressed in B. **/
MassProperties_ calcShiftedMassProps(const Vec<3,P>& newOriginB) const {
    return MassProperties_(mass, comInB-newOriginB,
                            calcShiftedInertia(newOriginB));
}

/** %Transform these mass properties from the current frame "B" to a new
frame "C", given the pose of C in B.\ Caution: this \e shifts
the point from which the mass properties are measured from the origin of B to
the origin of C. See reexpress() to change only the measure numbers without
moving the "measured from" point. Note that the frame in which a MassProperties 
object is expressed, and the point about which the mass properties are 
measured, are implicit; we don't actually have any way to verify that 
it is in B. Make sure you are certain about the current frame before using this 
method. **/
MassProperties_ calcTransformedMassProps(const Transform_<P>& X_BC) const {
    return MassProperties_(mass, ~X_BC*comInB, calcTransformedInertia(X_BC));
}

/** Re-express these mass properties from the current frame "B" to a new
frame "C", given the orientation of C in B.\ Caution: this does not \e shift
the point from which the mass properties are measured, it just uses a different
frame to express that measurement. See calcTransformedMassProps() to perform
a shift as well. Note that the frame in which a MassProperties object is 
expressed is implicit; we don't actually have any way to verify that it is in 
B. Make sure you are certain about the current frame before using this 
method. **/
MassProperties_ reexpress(const Rotation_<P>& R_BC) const {
    return MassProperties_(mass, ~R_BC*comInB, unitInertia_OB_B.reexpress(R_BC));
}

/** Return true only if the mass stored here is \e exactly zero.\ If the mass
resulted from a computation, you should use isNearlyMassless() instead.
@see isNearlyMassless(), isExactlyCentral() **/
bool isExactlyMassless()   const { return mass==0; }
/** Return true if the mass stored here is zero to within a small tolerance.
By default we use SignificantReal (about 1e-14 in double precision) as the
tolerance but you can override that. If you are just checking to see whether
the mass was explicitly set to zero (rather than calculated) you can use
isExactlyMassless() instead. @see isExactlyMassless(), isNearlyCentral() **/
bool isNearlyMassless(const P& tol=SignificantReal) const { 
    return mass <= tol; 
}

/** Return true only if the mass center stored here is \e exactly zero.\ If 
the mass center resulted from a computation, you should use isNearlyCentral()
instead. @see isNearlyCentral(), isExactlyMassless() **/
bool isExactlyCentral() const { return comInB==Vec<3,P>(0); }
/** Return true if the mass center stored here is zero to within a small tolerance.
By default we use SignificantReal (about 1e-14 in double precision) as the
tolerance but you can override that. If you are just checking to see whether
the mass center was explicitly set to zero (rather than calculated) you can use
isExactlyCentral() instead. @see isExactlyCentral(), isNearlyMassless() **/
bool isNearlyCentral(const P& tol=SignificantReal) const {
    return comInB.normSqr() <= tol*tol;
}

/** Return true if any element of this MassProperties object is NaN. 
@see isInf(), isFinite() **/
bool isNaN() const {return SimTK::isNaN(mass) || comInB.isNaN() || unitInertia_OB_B.isNaN();}
/** Return true only if there are no NaN's in this MassProperties object, and
at least one of the elements is Infinity.\ Ground's mass properties satisfy
these conditions. 
@see isNan(), isFinite() **/
bool isInf() const {
    if (isNaN()) return false;
    return SimTK::isInf(mass) || comInB.isInf() || unitInertia_OB_B.isInf();
}
/** Return true if none of the elements of this MassProperties object are
NaN or Infinity.\ Note that Ground's mass properties are not finite. 
@see isNaN(), isInf() **/
bool isFinite() const {
    return SimTK::isFinite(mass) && comInB.isFinite() && unitInertia_OB_B.isFinite();
}

/** Convert this MassProperties object to a spatial inertia matrix and return
it as a SpatialMat, which is a 2x2 matrix of 3x3 submatrices. 
@see toMat66() **/
Mat<2,2,Mat<3,3,P>> toSpatialMat() const {
    Mat<2,2,Mat<3,3,P>> M;
    M(0,0) = mass*unitInertia_OB_B.toMat33();
    M(0,1) = mass*crossMat(comInB);
    M(1,0) = ~M(0,1);
    M(1,1) = mass; // a diagonal matrix
    return M;
}

/** Convert this MassProperties object to a spatial inertia matrix in the
form of an ordinary 6x6 matrix, \e not a SpatialMat. Logically these are
the same but the ordering of the elements in memory is different between
a Mat66 and SpatialMat. 
@see toSpatialMat() **/
Mat<6,6,P> toMat66() const {
    Mat<6,6,P> M;
    M.template updSubMat<3,3>(0,0) = mass*unitInertia_OB_B.toMat33();
    M.template updSubMat<3,3>(0,3) = mass*crossMat(comInB);
    M.template updSubMat<3,3>(3,0) = ~M.template getSubMat<3,3>(0,3);
    M.template updSubMat<3,3>(3,3) = mass; // a diagonal matrix
    return M;
}

private:
P               mass;
Vec<3,P>        comInB;         // meas. from B origin, expr. in B
UnitInertia_<P> unitInertia_OB_B;   // about B origin, expr. in B
};

/** Output a human-readable representation of a MassProperties object to
the given output stream. This is assumed to be the mass properties of some
body B. We'll show B's mass, center of mass as a vector from B's origin,
expressed in B, and \e unit inertia, abbreviated Uxx, Uyy, and so on so that 
you won't accidentally think these are mass-scaled inertias normally designated 
Ixx, Iyy, etc. Note that I=mass*U in these terms. Also note that the unit 
inertia is taken about the body frame origin, \e not about the center of 
mass. And of course the unit inertia is expressed in B. 
@relates MassProperties_ **/
template <class P> static inline std::ostream& 
operator<<(std::ostream& o, const MassProperties_<P>& mp) {
    return o << "{ mass=" << mp.getMass() 
             << "\n  com=" << mp.getMassCenter()
             << "\n  Uxx,yy,zz=" << mp.getUnitInertia().getMoments()
             << "\n  Uxy,xz,yz=" << mp.getUnitInertia().getProducts()
             << "\n}\n";
}

} // namespace SimTK

#endif // SimTK_SIMMATRIX_MASS_PROPERTIES_H_
