//-----------------------------------------------------------------------------
// File:     Rotation.h
// Class:    Rotation and InverseRotation 
// Parent:   Mat33
// Purpose:  3x3 rotation class relating two right-handed orthogonal bases
//-----------------------------------------------------------------------------
#ifndef SimTK_ROTATION_H 
#define SimTK_ROTATION_H 

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman and Paul Mitiguy                                  *
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

//-----------------------------------------------------------------------------
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/UnitVec.h"
#include "SimTKcommon/internal/Quaternion.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
namespace SimTK {

//-----------------------------------------------------------------------------
// Forward declarations
class InverseRotation;

//-----------------------------------------------------------------------------
/**
 * This class is a Mat33 plus an ironclad guarantee that the matrix represents
 * a pure rotation. Having this as a separate class allows us both to ensure
 * that we have a legitimate rotation matrix, and to take advantage of
 * knowledge about the properties of these matrices. For example, multiplication
 * by a rotation matrix preserves the length of a vector so unit vectors
 * are still unit vectors afterwards and don't need to be normalized again.
 * 
 * A rotation is an orthogonal matrix whose columns and rows
 * are directions (that is, unit vectors) which are mutually orthogonal. 
 * Furthermore, if the columns (or rows) are labeled x,y,z it always holds
 * that z = x X y (rather than -(x X y)) ensuring that this is a rotation
 * and not a reflection.
 *
 * A rotation matrix is formed from a Cartesian coordinate frame
 * simply by using the x,y,z axes of the coordinate frame as
 * the columns of the rotation matrix. That is, if you have 
 * a frame F with its axes column vectors x,y,z expressed in frame G, you can write
 * this as a matrix  rot_GF = [ x y z ] which when applied to
 * a vector with measure numbers expressed in F yields that same
 * vector re-expressed in G: v_G = rot_GF * v_F. Because a rotation
 * is orthogonal its transpose is its inverse, which is also a rotation.
 * So we write rot_FG = ~rot_GF (where "~" is SimTK for "transpose"), and
 * this matrix can be used to rotate in the other direction: 
 *      v_F = rot_FG * v_G = ~rot_GF * v_G.
 */
//-----------------------------------------------------------------------------
class Rotation : public Mat33 {
public:
    typedef Mat33 BaseMat;
    typedef UnitVec<BaseMat::RowSpacing> ColType;
    typedef UnitRow<BaseMat::ColSpacing> RowType;

    Rotation() : BaseMat(1.) { }    // default is identity
    // default copy constructor, copy assignment and destructor

    // Like copy constructor and copy assign but for inverse rotation;
    // constructor is also implicit conversion from InverseRotation to Rotation.
    inline Rotation(const InverseRotation&);
    inline Rotation& operator=(const InverseRotation&);

    // This constructor takes a Mat33 and orthogonalizes it into a (hopefully
    // nearby) rotation matrix.
    SimTK_SimTKCOMMON_EXPORT explicit Rotation(const Mat33&);

    // These static methods are like constructors with friendlier names.

    static Rotation zero() {return Rotation();}
    static Rotation NaN()  {Rotation r;r.setToNaN();return r;}

    // One-angle rotations.
    static Rotation aboutX(const Real& angInRad) {
        const Real s = std::sin(angInRad), c = std::cos(angInRad);
        return Rotation( 1 ,  0 ,  0 ,
                         0 ,  c , -s ,
                         0 ,  s ,  c );
    }
    static Rotation aboutY(const Real& angInRad) {
        const Real s = std::sin(angInRad), c = std::cos(angInRad);
        return Rotation( c ,  0 ,  s ,
                         0 ,  1 ,  0 ,
                        -s ,  0 ,  c );
    }
    static Rotation aboutZ(const Real& angInRad) {
        const Real s = std::sin(angInRad), c = std::cos(angInRad);
        return Rotation( c , -s ,  0 ,
                         s ,  c ,  0 ,
                         0 ,  0 ,  1 );
    }

    SimTK_SimTKCOMMON_EXPORT static Rotation aboutAxis
       (const Real& angInRad, const UnitVec3& axis);

    static Rotation aboutAxis(const Real& angInRad, const Vec3& axis)
      { return aboutAxis(angInRad, UnitVec3(axis)); }

    // Two-angle space-fixed rotations.
    static Rotation aboutXThenOldY(const Real& xInRad, const Real& yInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        return Rotation( cy   ,  sy*sx ,  sy*cx ,
                          0   ,   cx   ,  -sx   ,
                        -sy   ,  cy*sx ,  cy*cx );
    }
    static Rotation aboutYThenOldX(const Real& yInRad, const Real& xInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        return Rotation(  cy   ,    0   ,   sy   ,
                         sy*sx ,   cx   , -cy*sx ,
                        -sy*cx ,   sx   ,  cy*cx );
    }
    static Rotation aboutXThenOldZ(const Real& xInRad, const Real& zInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation(  cz   , -sz*cx ,  sz*sx ,
                          sz   ,  cz*cx , -cz*sx ,
                           0   ,   sx   ,   cx   );
    }
    static Rotation aboutZThenOldX(const Real& zInRad, const Real& xInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation(  cz   ,  -sz   ,    0   ,
                         sz*cx ,  cz*cx ,  -sx   ,
                         sz*sx ,  cz*sx ,   cx   );
    }
    static Rotation aboutYThenOldZ(const Real& yInRad, const Real& zInRad) {
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation( cz*cy ,  -sz   ,  cz*sy ,
                         sz*cy ,   cz   ,  sz*sy ,
                         -sy   ,    0   ,   cy   );
    }
    static Rotation aboutZThenOldY(const Real& zInRad, const Real& yInRad) {
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation( cz*cy , -sz*cy ,   sy   ,
                          sz   ,   cz   ,    0   ,
                        -cz*sy ,  sz*sy ,   cy   );
    }

    // Two-angle body fixed rotations. These are the same as the
    // reversed space-fixed ones.
    static Rotation aboutXThenNewY(const Real& xInRad, const Real& yInRad)
      { return aboutYThenOldX(yInRad, xInRad); }
    static Rotation aboutYThenNewX(const Real& yInRad, const Real& xInRad)
      { return aboutXThenOldY(xInRad, yInRad); }
    static Rotation aboutXThenNewZ(const Real& xInRad, const Real& zInRad)
      { return aboutZThenOldX(zInRad, xInRad); }
    static Rotation aboutZThenNewX(const Real& zInRad, const Real& xInRad)
      { return aboutXThenOldZ(xInRad, zInRad); }
    static Rotation aboutYThenNewZ(const Real& yInRad, const Real& zInRad)
      { return aboutZThenOldY(zInRad, yInRad); }
    static Rotation aboutZThenNewY(const Real& zInRad, const Real& yInRad)
      { return aboutYThenOldZ(yInRad, zInRad); }

    /// Create a Rotation matrix by specifying only its z axis. 
    /// The resulting x and y axes will be appropriately perpendicular
    /// but are otherwise arbitrary. This will work for any stride
    /// UnitVec because there is always an implicit conversion
    /// available to the packed form used as the argument.
    SimTK_SimTKCOMMON_EXPORT explicit Rotation(const UnitVec3& z);

    /// Create a Rotation matrix by specifying its x axis, and then
    /// a "y like" axis. We will take x seriously after normalizing, but use the y only
    /// to create z = normalize(x X y), then y = z X x. Bad things will 
    /// happen if x and y are aligned but we may not catch it.
    SimTK_SimTKCOMMON_EXPORT explicit Rotation(const Vec3& x, const Vec3& yish);


    /// Create a rotation matrix from a quaternion.
    SimTK_SimTKCOMMON_EXPORT explicit Rotation(const Quaternion&);

    /// Set this Rotation to represent a rotation of +q0 about
    /// the body frame's Z axis, followed by a rotation of +q1 about
    /// the body frame's NEW Y axis, followed by a rotation of +q3
    /// about the body frame's NEW X axis.
    /// See Kane, Spacecraft Dynamics, pg. 423, body-three: 3-2-1.
    SimTK_SimTKCOMMON_EXPORT void setToBodyFixed321(const Vec3&);

    /// Set this Rotation to represent a rotation of +q0 about
    /// the body frame's X axis, followed by a rotation of +q1 about
    /// the body frame's NEW Y axis, followed by a rotation of +q3
    /// about the body frame's NEW Z axis.
    /// See Kane, Spacecraft Dynamics, pg. 423, body-three: 1-2-3.
    SimTK_SimTKCOMMON_EXPORT void setToBodyFixed123(const Vec3&);

    /// Set this Rotation to represent the same rotation as
    /// the passed-in quaternion.
    SimTK_SimTKCOMMON_EXPORT void setToQuaternion(const Quaternion&);

    /// Convert this Rotation matrix to the equivalent quaternion.
    SimTK_SimTKCOMMON_EXPORT Quaternion convertToQuaternion() const;

    /// Convert this Rotation matrix to the equivalent 1-2-3 body fixed
    /// Euler angle sequence.
    SimTK_SimTKCOMMON_EXPORT Vec3 convertToBodyFixed123() const;

    /// Convert this Rotation matrix to the equivalent 1-2 body fixed
    /// Euler angle sequence. The result is only meaningful if the
    /// Rotation matrix is one that can be produced by such a sequence.
    SimTK_SimTKCOMMON_EXPORT Vec2 convertToBodyFixed12() const;

    /// Convert this Rotation matrix to the equivalent 1-2 space fixed
    /// Euler angle sequence. The result is only meaningful if the
    /// Rotation matrix is one that can be produced by such a sequence.
    SimTK_SimTKCOMMON_EXPORT Vec2 convertToSpaceFixed12() const;

    /// Convert this Rotation matrix to an equivalent (angle,axis)
    /// representation: (a vx vy vz), with v a unit vector, a in radians.
    SimTK_SimTKCOMMON_EXPORT Vec4 convertToAngleAxis()  const;

    /// Return true if the passed-in Rotation is identical to the current one
    /// to within a cone whose half-angle is supplied (in radians). That half-angle
    /// is called the "pointing error", and we require 0 < pointing error < pi since
    /// 0 would require perfect precision and pi (180 degrees) would return true
    /// for any pair of Rotations.
    SimTK_SimTKCOMMON_EXPORT bool isSameRotationToWithinAngle(const Rotation&, Real okPointingError) const;

    /// Return true if the passed-in Rotation is identical to the current one
    /// to within what can be expected at the precision being used.
    SimTK_SimTKCOMMON_EXPORT bool isSameRotationToMachinePrecision(const Rotation&) const;


    Rotation(const Rotation& R) : BaseMat(R) { }
    Rotation& operator=(const Rotation& R) {
        BaseMat::operator=(R.asMat33()); return *this;
    }

    /// By zero we mean "zero rotation", i.e., an identity matrix.
    Rotation& setToZero() {
        BaseMat::operator=(Real(1)); return *this;
    }

    Rotation& setToNaN() {
        BaseMat::setToNaN(); return *this;
    }

    /// Set this Rotation to represent a rotation of +q radians
    /// around the base frame's 0,0,1 axis.
    void setToRotationAboutZ(const Real& q) {
        const Real sq = std::sin(q), cq = std::cos(q);
        BaseMat::operator=(BaseMat( cq , -sq , 0. ,
                                    sq ,  cq , 0. ,
                                    0. ,  0. , 1. ));
    }

    /// Set this Rotation to represent a rotation of +q0 about
    /// the base frame's X axis, followed by a rotation of +q1 about
    /// the base frame's (unchanged) Y axis.
    void setToSpaceFixed12(const Vec2& q) {
        const Real sq0 = std::sin(q[0]), cq0 = std::cos(q[0]);
        const Real sq1 = std::sin(q[1]), cq1 = std::cos(q[1]);
        BaseMat::operator=(BaseMat( cq1 , sq1*sq0 , sq1*cq0,
                                     0. ,   cq0   ,   -sq0 ,
                                   -sq1 , cq1*sq0 , cq1*cq0));
    }


    /// Given Euler angles forming a body-fixed 3-2-1 sequence, and the relative
    /// angular velocity vector of B in the parent frame, *BUT EXPRESSED IN
    /// THE BODY FRAME*, return the Euler angle
    /// derivatives. You are dead if q[1] gets near 90 degrees!
    /// See Kane's Spacecraft Dynamics, page 428, body-three: 3-2-1.
    static Vec3 convertAngVelToBodyFixed321Dot(const Vec3& q, const Vec3& w_PB_B) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1 = Real(1)/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E( 0. ,   s2oc1  ,  c2oc1  ,
                       0. ,     c2   ,   -s2   ,
                       1. , s1*s2oc1 , s1*c2oc1 );
        return E * w_PB_B;
    }


    /// Inverse of the above routine. Returned angular velocity is B in P,
    /// expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed321DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 Einv(  -s1  ,  0. , 1. ,
                          c1*s2 ,  c2 , 0. ,
                          c1*c2 , -s2 , 0. );
        return Einv*qd;
    }
    
    /// Given Euler angles forming a body-fixed 1-2-3 sequence, and the relative
    /// angular velocity vector of B in the parent frame,  *BUT EXPRESSED IN
    /// THE BODY FRAME*, return the Euler angle
    /// derivatives. You are dead if q[1] gets near 90 degrees!
    /// See Kane's Spacecraft Dynamics, page 427, body-three: 1-2-3.
    static Vec3 convertAngVelToBodyFixed123Dot(const Vec3& q, const Vec3& w) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E(    c2oc1  , -s2oc1  , 0.,
                            s2   ,    c2   , 0.,
                       -s1*c2oc1 , s1*s2oc1, 1. );
        return E*w;
    }

    /// Inverse of the above routine. Returned angular velocity is B in P,
    /// expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed123DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 Einv( c1*c2 ,  s2 , 0. ,
                         -c1*s2 ,  c2 , 0. ,
                           s1   ,  0. , 1. );
        return Einv*qd;
    }

    // TODO: sherm: is this right? Warning: everything is measured in the
    // *PARENT* frame, but has to be expressed in the *BODY* frame.
    static Vec3 convertAngVelDotToBodyFixed321DotDot
        (const Vec3& q, const Vec3& w_PB_B, const Vec3& wdot)
    {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E( 0. ,   s2oc1  ,  c2oc1  ,
                       0. ,     c2   ,   -s2   ,
                       1. , s1*s2oc1 , s1*c2oc1 );
        const Vec3 qdot = E * w_PB_B;

        const Real t =  qdot[1]*qdot[2]*s1*ooc1;
        const Real a =  t*c2oc1; // d/dt s2oc1
        const Real b = -t*s2oc1; // d/dt c2oc1

        const Mat33 Edot( 0. ,       a           ,         b         ,
                          0. ,   -qdot[2]*s2     ,    -qdot[2]*c2    ,
                          0. , s1*a + qdot[1]*s2 , s1*b + qdot[1]*c2 );

        return E*wdot + Edot*w_PB_B;
    }

    // TODO: sherm: is this right? Warning: everything is measured in the
    // *PARENT* frame, but has to be expressed in the *BODY* frame.
    static Vec3 convertAngVelDotToBodyFixed123DotDot
        (const Vec3& q, const Vec3& w_PB_B, const Vec3& wdot)
    {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E(    c2oc1  , -s2oc1  , 0.,
                            s2   ,    c2   , 0.,
                       -s1*c2oc1 , s1*s2oc1, 1. );
        const Vec3 qdot = E * w_PB_B;

        const Real t =  qdot[1]*qdot[2]*s1*ooc1;
        const Real a =  t*c2oc1; // d/dt s2oc1
        const Real b = -t*s2oc1; // d/dt c2oc1

        const Mat33 Edot(       b           ,        -a         , 0.,
                             qdot[2]*c2     ,    -qdot[2]*s2    , 0.,
                         -s1*b - qdot[1]*c2 , s1*a + qdot[1]*s2 , 0. );

        return E*wdot + Edot*w_PB_B;
    }


    /// Given a possibly unnormalized quaternion (0th element is the scalar) and the
    /// relative angular velocity vector of B in its parent, expressed 
    /// in the *PARENT*, return the quaternion derivatives. This is never singular.
    static Vec4 convertAngVelToQuaternionDot(const Vec4& q, const Vec3& w_PB_P) {
        const Mat43 E(-q[1],-q[2],-q[3],
                       q[0], q[3],-q[2],    // TODO: signs???
                      -q[3], q[0], q[1],
                       q[2],-q[1], q[0]);
        return 0.5*(E*w_PB_P);
    }

    /// Everything is measured and expressed in the parent.
    static Vec4 convertAngVelDotToQuaternionDotDot
        (const Vec4& q, const Vec3& w_PB_P, const Vec3& wdot)
    {
        const Mat43 E(-q[1],-q[2],-q[3],
                       q[0], q[3],-q[2],    // TODO: signs???
                      -q[3], q[0], q[1],
                       q[2],-q[1], q[0]);
        const Vec4  qdot = 0.5*E*w_PB_P;
        const Mat43 Edot(-qdot[1],-qdot[2],-qdot[3],
                          qdot[0], qdot[3],-qdot[2],
                         -qdot[3], qdot[0], qdot[1],
                          qdot[2],-qdot[1], qdot[0]);

        return 0.5*(Edot*w_PB_P + E*wdot);
    }

    /// Inverse of the above routine. Returned AngVel is expressed in
    /// the *PARENT* frame: w_PB_P.
    static Vec3 convertQuaternionDotToAngVel(const Vec4& q, const Vec4& qd) {
        const Mat34 Et(-q[1], q[0],-q[3], q[2],  // TODO: signs???
                       -q[2], q[3], q[0],-q[1],
                       -q[3],-q[2], q[1], q[0]);
        return 2.*(Et*qd);
    }

    const InverseRotation& invert() const {
        return *reinterpret_cast<const InverseRotation*>(this);
    }
    InverseRotation& updInvert() {
        return *reinterpret_cast<InverseRotation*>(this);
    }

    // Override the Mat33 versions of transpose.
    const InverseRotation& transpose() const {return invert();}
    InverseRotation&       updTranspose() {return updInvert();}

    // Note that this does not have unit stride.
    const RowType& row(int i) const {
        return reinterpret_cast<const RowType&>(asMat33()[i]);
    }
    const ColType& col(int j) const {
        return reinterpret_cast<const ColType&>(asMat33()(j));
    }
    const ColType& x() const {return col(0);}
    const ColType& y() const {return col(1);}
    const ColType& z() const {return col(2);}

    inline Rotation& operator*=(const Rotation& R);
    inline Rotation& operator/=(const Rotation& R);
    inline Rotation& operator*=(const InverseRotation&);
    inline Rotation& operator/=(const InverseRotation&);

    const InverseRotation& operator~() const {return invert();}
    InverseRotation&       operator~()       {return updInvert();}

    const RowType& operator[](int i) const {return row(i);}
    const ColType& operator()(int j) const {return col(j);}

    const BaseMat& asMat33() const {
        return *static_cast<const BaseMat*>(this);
    }

    /// Less efficient version of asMat33() since it copies, but you don't
    /// have to know the internal layout.
    BaseMat toMat33() const {
        return asMat33();
    }

    static Rotation trustMe(const Mat33& m) {return Rotation(m, true);}

private:
    // We're trusting that m is a rotation.
    explicit Rotation(const Mat33& m, bool) : Mat33(m) { }

    // This is only for the most trustworthy of callers, that is, methods
    // of the Rotation class. There are a lot of ways for this NOT to
    // be a legitimate rotation matrix -- be careful!!
    // Note that these are supplied in rows.
    Rotation( const Real& xx, const Real& xy, const Real& xz,
              const Real& yx, const Real& yy, const Real& yz,
              const Real& zx, const Real& zy, const Real& zz )
      : BaseMat(xx,xy,xz, yx,yy,yz, zx,zy,zz)
    {
    }
};


///-----------------------------------------------------------------------------
///  This InverseRotation class is the inverse of a Rotation 
///  See the Rotation class for information.
///-----------------------------------------------------------------------------
class InverseRotation : public Mat33::TransposeType {
public:
    typedef Mat33::TransposeType BaseMat;
    typedef UnitVec<BaseMat::RowSpacing> ColType;
    typedef UnitRow<BaseMat::ColSpacing> RowType;

    // Should not usually construct one of these as they should only occur as expression intermediates.
    // But if you must ...
    InverseRotation() : BaseMat(1) { }

    InverseRotation( const InverseRotation& R ) : BaseMat(R) { }
    InverseRotation&  operator=( const InverseRotation& R )  { BaseMat::operator=( R.asMat33() );  return *this; }

    const Rotation&  invert() const { return *reinterpret_cast<const Rotation*>(this); }
    Rotation&  updInvert()          { return *reinterpret_cast<Rotation*>(this); }

    // Override the Mat33 versions of transpose.
    const Rotation&  transpose() const  { return invert(); }
    Rotation&        updTranspose()     { return updInvert(); }

    // Note that this does not have unit stride.
    const RowType&  row(int i) const  { return reinterpret_cast<const RowType&>(asMat33()[i]); }
    const ColType&  col(int j) const  { return reinterpret_cast<const ColType&>(asMat33()(j)); }
    const ColType&  x() const  { return col(0); }
    const ColType&  y() const  { return col(1); }
    const ColType&  z() const  { return col(2); }

    const Rotation&  operator~() const  { return invert(); }
    Rotation&        operator~()        { return updInvert(); }

    const RowType&  operator[]( int i ) const  { return row(i); }
    const ColType&  operator()( int j ) const  { return col(j); }

    const BaseMat&  asMat33() const  { return *static_cast<const BaseMat*>(this); }

    /// Less efficient version of asMat33() since it copies, but  
    /// you don't have to know the internal layout.
    BaseMat  toMat33() const  { return asMat33(); }
};


SimTK_SimTKCOMMON_EXPORT std::ostream&  operator<<( std::ostream& o, const Rotation& m );

template <int S> inline UnitVec<1>  operator*( const Rotation& R,        const UnitVec<S>& v )       { return UnitVec<1>(R.asMat33()*v.asVec3(),  true); }
template <int S> inline UnitRow<1>  operator*( const UnitRow<S>& r,      const Rotation& R   )       { return UnitRow<1>(r.asRow3(), R.asMat33(), true); }
template <int S> inline UnitVec<1>  operator*( const InverseRotation& R, const UnitVec<S>& v )       { return UnitVec<1>(R.asMat33()*v.asVec3(),  true); }
template <int S> inline UnitRow<1>  operator*( const UnitRow<S>& r,      const InverseRotation& R )  { return UnitRow<1>(r.asRow3(), R.asMat33(), true); }

inline Rotation::Rotation( const InverseRotation& R) : Mat33(R.asMat33() )  { }

inline Rotation&  Rotation::operator=(  const InverseRotation& R )  { static_cast<BaseMat&>(*this) = R.asMat33();     return *this; }
inline Rotation&  Rotation::operator*=( const Rotation& R )         { static_cast<BaseMat&>(*this) *= R.asMat33();    return *this; }
inline Rotation&  Rotation::operator/=( const Rotation& R )         { static_cast<BaseMat&>(*this) *= (~R).asMat33(); return *this; }
inline Rotation&  Rotation::operator*=( const InverseRotation& R )  { static_cast<BaseMat&>(*this) *= R.asMat33();    return *this; }
inline Rotation&  Rotation::operator/=( const InverseRotation& R )  { static_cast<BaseMat&>(*this) *= (~R).asMat33(); return *this; }

inline Rotation  operator*( const Rotation&        R1, const Rotation&        R2 )  { return Rotation(R1) *= R2; }
inline Rotation  operator*( const Rotation&        R1, const InverseRotation& R2 )  { return Rotation(R1) *= R2; }
inline Rotation  operator*( const InverseRotation& R1, const Rotation&        R2 )  { return Rotation(R1) *= R2; }
inline Rotation  operator*( const InverseRotation& R1, const InverseRotation& R2 )  { return Rotation(R1) *= R2; }

inline Rotation operator/( const Rotation&        R1, const Rotation&        R2 )  {return Rotation(R1)/=R2;}
inline Rotation operator/( const Rotation&        R1, const InverseRotation& R2 )  {return Rotation(R1)/=R2;}
inline Rotation operator/( const InverseRotation& R1, const Rotation&        R2 )  {return Rotation(R1)/=R2;}
inline Rotation operator/( const InverseRotation& R1, const InverseRotation& R2 )  {return Rotation(R1)/=R2;}


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_ROTATION_H_
//--------------------------------------------------------------------------

