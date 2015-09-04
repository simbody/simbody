#ifndef SimTK_SimTKCOMMON_ROTATION_H_ 
#define SimTK_SimTKCOMMON_ROTATION_H_ 

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-14 Stanford University and the Authors.        *
 * Authors: Paul Mitiguy, Michael Sherman                                     *
 * Contributors:                                                              *
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

//------------------------------------------------------------------------------

#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/UnitVec.h"
#include "SimTKcommon/internal/Quaternion.h"

//------------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace SimTK {


enum BodyOrSpaceType { BodyRotationSequence=0, SpaceRotationSequence=1 };

//------------------------------------------------------------------------------
// Forward declarations
template <class P> class Rotation_;
template <class P> class InverseRotation_;

typedef Rotation_<Real>             Rotation;
typedef Rotation_<float>           fRotation;
typedef Rotation_<double>          dRotation;

typedef InverseRotation_<Real>      InverseRotation;
typedef InverseRotation_<float>    fInverseRotation;
typedef InverseRotation_<double>   dInverseRotation;

//------------------------------------------------------------------------------
/** The Rotation class is a Mat33 that guarantees that the matrix can be 
interpreted as a legitimate 3x3 rotation matrix giving the relative orientation 
of two right-handed, orthogonal, unit vector bases. 

A rotation matrix, also known as a direction cosine matrix, is an orthogonal 
matrix whose columns and rows are directions (that is, unit vectors) that are 
mutually orthogonal. Furthermore, if the columns (or rows) are labeled x,y,z it 
always holds that z = x X y (rather than -(x X y)) ensuring that this is a 
right-handed rotation matrix and not a reflection. This is equivalent to saying 
that the determinant of a rotation matrix is 1, not -1.

The Rotation class takes advantage of known properties of orthogonal matrices. 
For example, multiplication by a rotation matrix preserves a vector's length so 
unit vectors are still unit vectors afterwards and don't need to be 
re-normalized.

Suppose there is a vector v_F expressed in terms of the right-handed, 
orthogonal unit vectors Fx, Fy, Fz and one would like to express v instead
as v_G, in terms of a right-handed, orthogonal unit vectors Gx, Gy, Gz. To 
calculate it, we form a rotation matrix R_GF whose columns are the F unit 
vectors re-expressed in G:
<pre>
            G F   (      |      |      )
     R_GF =  R  = ( Fx_G | Fy_G | Fz_G )
                  (      |      |      )
where
     Fx_G = ~( ~Fx*Gx, ~Fx*Gy, ~Fx*Gz ), etc.
</pre>
(~Fx*Gx means dot(Fx,Gx)). Note that we use "monogram" notation R_GF in 
code to represent the more typographically demanding superscripted notation 
for rotation matrices. Now we can re-express the vector v from frame F to 
frame G via
<pre>
     v_G = R_GF * v_F. 
</pre>
Because a rotation is orthogonal, its transpose is its inverse. Hence
R_FG = ~R_GF (where ~ is the SimTK "transpose" operator). This transpose 
matrix can be used to expressed v_G in terms of Fx, Fy, Fz as
<pre>
     v_F = R_FG * v_G  or  v_F = ~R_GF * v_G
</pre>
In either direction, correct behavior can be obtained by using the 
recommended notation and then matching up the frame labels (after
interpreting the "~" operator as reversing the labels).

The Rotation_ class is templatized by the precision P, which should be float
or double. A typedef defining type Rotation as Rotation_<Real> is always 
defined and is normally used in user programs rather than the templatized class.

\nosubgrouping
**/
//------------------------------------------------------------------------------
template <class P> // templatized by precision
class Rotation_ : public Mat<3,3,P> {
public:
typedef P               RealP; ///< These are just local abbreviations.
typedef Mat<2,2,P>      Mat22P;
typedef Mat<3,2,P>      Mat32P;
typedef Mat<3,3,P>      Mat33P;
typedef Mat<4,3,P>      Mat43P;
typedef Mat<3,4,P>      Mat34P;
typedef Vec<2,P>        Vec2P;
typedef Vec<3,P>        Vec3P;
typedef Vec<4,P>        Vec4P;
typedef UnitVec<P,1>    UnitVec3P; // stride is 1 here, length is always 3
typedef SymMat<3,P>     SymMat33P;
typedef Quaternion_<P>  QuaternionP;

/** This is the type of a column of this %Rotation matrix. It will be a Vec<3>
but will not necessarily have unit spacing. **/
typedef UnitVec<P,Mat33P::RowSpacing> ColType;
/** This is the type of a row of this %Rotation matrix. It will be a Row<3>
but will not necessarily have unit spacing. **/
typedef UnitRow<P,Mat33P::ColSpacing> RowType;

//------------------------------------------------------------------------------
/** @name Constructors, Mutators, and Assignment **/
//@{
/** Default constructor. **/
Rotation_() : Mat33P(1) {}    

/** Copy constructor. **/
Rotation_( const Rotation_& R ) : Mat33P(R)  {}
/** Like copy constructor but for inverse rotation. This allows implicit 
conversion from InverseRotation_ to Rotation_. **/
inline Rotation_( const InverseRotation_<P>& );

/** Assignment operator. **/
Rotation_&  operator=( const Rotation_& R )  
{ Mat33P::operator=( R.asMat33() );  return *this; }
/** Like copy assignment but for inverse rotation. **/
inline Rotation_& operator=( const InverseRotation_<P>& );

/** Construct Rotation_ filled with NaNs. **/
Rotation_&  setRotationToNaN()             
{ Mat33P::setToNaN();    return *this; } 

/** Construct identity Rotation_. **/
Rotation_&  setRotationToIdentityMatrix()  
{ Mat33P::operator=(RealP(1));  return *this; }

/** Constructor for right-handed rotation by an angle (in radians) about a 
coordinate axis. **/
Rotation_( RealP angle, const CoordinateAxis& axis )             
{ setRotationFromAngleAboutAxis( angle, axis ); }
/** Set this Rotation_ object to a right-handed rotation by an angle (in 
radians) about a coordinate axis. **/
Rotation_& setRotationFromAngleAboutAxis(RealP angle, const CoordinateAxis& axis)  
{ return axis.isXAxis() ? setRotationFromAngleAboutX(angle) 
      : (axis.isYAxis() ? setRotationFromAngleAboutY(angle) 
      : setRotationFromAngleAboutZ(angle) ); }

/** Constructor for right-handed rotation by an angle (in radians) about the 
X-axis. **/
Rotation_( RealP angle, const CoordinateAxis::XCoordinateAxis )  
{ setRotationFromAngleAboutX( std::cos(angle), std::sin(angle) ); }
/** Set this Rotation_ object to a right-handed rotation by an angle (in 
radians) about the X-axis. **/
Rotation_&  setRotationFromAngleAboutX( RealP angle )  
{ return setRotationFromAngleAboutX( std::cos(angle), std::sin(angle) ); }
/** Set this Rotation_ object to a right-handed rotation by an angle about the
X-axis, where the cosine and sine of the angle are specified. **/
Rotation_&  setRotationFromAngleAboutX( RealP cosAngle, RealP sinAngle )  
{ Mat33P& R = *this;  R[0][0] = 1;   R[0][1] = R[0][2] = R[1][0] = R[2][0] = 0;   
  R[1][1] = R[2][2] = cosAngle;  R[1][2] = -(R[2][1] = sinAngle);  
  return *this; }

/** Constructor for right-handed rotation by an angle (in radians) about the 
Y-axis. **/
Rotation_( RealP angle, const CoordinateAxis::YCoordinateAxis )  
{ setRotationFromAngleAboutY( std::cos(angle), std::sin(angle) ); }
/** Set this Rotation_ object to a right-handed rotation by an angle (in 
radians) about the Y-axis. **/
Rotation_&  setRotationFromAngleAboutY( RealP angle )  
{ return setRotationFromAngleAboutY( std::cos(angle), std::sin(angle) ); }
/** Set this Rotation_ object to a right-handed rotation by an angle about the
Y-axis, where the cosine and sine of the angle are specified. **/
Rotation_&  setRotationFromAngleAboutY( RealP cosAngle, RealP sinAngle )  
{ Mat33P& R = *this;  R[1][1] = 1;   R[0][1] = R[1][0] = R[1][2] = R[2][1] = 0;   
  R[0][0] = R[2][2] = cosAngle;  R[2][0] = -(R[0][2] = sinAngle);  
  return *this; }

/** Constructor for right-handed rotation by an angle (in radians) about the 
Z-axis. **/
Rotation_( RealP angle, const CoordinateAxis::ZCoordinateAxis )  
{ setRotationFromAngleAboutZ( std::cos(angle), std::sin(angle) ); }
/** Set this Rotation_ object to a right-handed rotation by an angle (in 
radians) about the Z-axis. **/
Rotation_&  setRotationFromAngleAboutZ( RealP angle )  
{ return setRotationFromAngleAboutZ( std::cos(angle), std::sin(angle) ); }
/** Set this Rotation_ object to a right-handed rotation by an angle about the
Z-axis, where the cosine and sine of the angle are specified. **/
Rotation_&  setRotationFromAngleAboutZ( RealP cosAngle, RealP sinAngle )  
{ Mat33P& R = *this;  R[2][2] = 1;   R[0][2] = R[1][2] = R[2][0] = R[2][1] = 0;   
  R[0][0] = R[1][1] = cosAngle;  R[0][1] = -(R[1][0] = sinAngle);  
  return *this; }

/** Constructor for right-handed rotation by an angle (in radians) about an 
arbitrary unit vector. **/
Rotation_( RealP angle, const UnitVec3P& unitVector ) 
{ setRotationFromAngleAboutUnitVector(angle,unitVector); }
/** Set this Rotation_ object to a right-handed rotation of an angle (in 
radians) about an arbitrary unit vector. **/
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setRotationFromAngleAboutUnitVector(RealP angle, const UnitVec3P& unitVector);

/** Constructor for right-handed rotation by an angle (in radians) about an 
arbitrary vector of arbitrary length. **/
Rotation_( RealP angle, const Vec3P& nonUnitVector )  
{ setRotationFromAngleAboutNonUnitVector(angle,nonUnitVector); }
/** Set this Rotation_ object to a right-handed rotation of an angle (in 
radians) about an arbitrary vector of arbitrary length. **/
Rotation_&  
setRotationFromAngleAboutNonUnitVector(RealP angle, const Vec3P& nonUnitVector)  
{ return setRotationFromAngleAboutUnitVector( angle, UnitVec3P(nonUnitVector) ); }

/** Constructor for two-angle, two-axes, Body-fixed or Space-fixed rotation 
sequences (angles are in radians). **/
Rotation_(BodyOrSpaceType bodyOrSpace, 
          RealP angle1, const CoordinateAxis& axis1, 
          RealP angle2, const CoordinateAxis& axis2)
{ setRotationFromTwoAnglesTwoAxes(bodyOrSpace,angle1,axis1,angle2,axis2); }
/** Set this Rotation_ object to a two-angle, two-axes, Body-fixed or 
Space-fixed rotation sequences (angles are in radians). **/
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setRotationFromTwoAnglesTwoAxes(BodyOrSpaceType bodyOrSpace, 
                                RealP angle1, const CoordinateAxis& axis1, 
                                RealP angle2, const CoordinateAxis& axis2); 

/** Constructor for three-angle Body-fixed or Space-fixed rotation sequences 
(angles are in radians). **/
Rotation_(BodyOrSpaceType bodyOrSpace, 
          RealP angle1, const CoordinateAxis& axis1, 
          RealP angle2, const CoordinateAxis& axis2, 
          RealP angle3, const CoordinateAxis& axis3 )  
{ setRotationFromThreeAnglesThreeAxes
   (bodyOrSpace,angle1,axis1,angle2,axis2,angle3,axis3); }
/** Set this Rotation_ object to a three-angle Body-fixed or Space-fixed 
rotation sequences (angles are in radians). **/
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setRotationFromThreeAnglesThreeAxes(BodyOrSpaceType bodyOrSpace, 
                                    RealP angle1, const CoordinateAxis& axis1, 
                                    RealP angle2, const CoordinateAxis& axis2, 
                                    RealP angle3, const CoordinateAxis& axis3);

/** Constructor for creating a rotation matrix from a quaternion. **/
explicit Rotation_( const QuaternionP& q )  { setRotationFromQuaternion(q); }
/** Method for creating a rotation matrix from a quaternion. **/
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setRotationFromQuaternion( const QuaternionP& q );

/** Constructs an (hopefully nearby) orthogonal rotation matrix from a 
generic Mat33P. **/
explicit Rotation_( const Mat33P& m )  { setRotationFromApproximateMat33(m); }
/** Set this Rotation_ object to an (hopefully nearby) orthogonal rotation 
matrix from a generic Mat33P. **/
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setRotationFromApproximateMat33( const Mat33P& m );

/** (Advanced) Construct a Rotation_ directly from a Mat33P (we trust that m is 
a valid Rotation_!) Things will not go well for you if it is not. **/
Rotation_(const Mat33P& m, bool) : Mat33P(m) {}
/** (Advanced) Set the Rotation_ matrix directly - but you had better know what 
you are doing! **/
Rotation_& setRotationFromMat33TrustMe(const Mat33P& m)  
{ Mat33P& R = *this; R=m;  return *this; }   
/** (Advanced) Set the Rotation_ matrix directly - but you had better know what 
you are doing! **/
Rotation_& setRotationColFromUnitVecTrustMe(int colj, const UnitVec3P& uvecj)  
{ Mat33P& R = *this; R(colj)=uvecj.asVec3(); return *this; }   
/** (Advanced) Set the Rotation_ matrix directly - but you had better know what 
you are doing! **/
Rotation_& setRotationFromUnitVecsTrustMe
   (const UnitVec3P& colA, const UnitVec3P& colB, const UnitVec3P& colC)  
{ Mat33P& R = *this; R(0)=colA.asVec3(); R(1)=colB.asVec3(); R(2)=colC.asVec3();
  return *this; }  

/** Calculate R_AB by knowing one of B's unit vectors expressed in A.
Note: The other vectors are perpendicular (but somewhat arbitrarily so). **/
Rotation_(const UnitVec3P& uvec, CoordinateAxis axis)  
{ setRotationFromOneAxis(uvec,axis); }
/** Calculate R_AB by knowing one of B's unit vectors expressed in A.
Note: The other vectors are perpendicular (but somewhat arbitrarily so). **/
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setRotationFromOneAxis(const UnitVec3P& uvec, CoordinateAxis axis);

/** Calculate R_AB by knowing one of B's unit vectors u1 (could be Bx, By, or Bz) 
expressed in A and a vector v (also expressed in A) that is approximately in the 
desired direction for a second one of B's unit vectors, u2 (!= u1). If v is not 
perpendicular to u1, no worries - we'll find a direction for u2 that is 
perpendicular to u1 and comes closest to v. The third vector u3 is +/- u1 X u2, 
as appropriate for a right-handed rotation matrix. **/
Rotation_(const UnitVec3P& uveci, const CoordinateAxis& axisi, 
          const Vec3P& vecjApprox, const CoordinateAxis& axisjApprox )  
{ setRotationFromTwoAxes(uveci,axisi,vecjApprox,axisjApprox); }
/** Calculate R_AB by knowing one of B's unit vectors u1 (could be Bx, By, or Bz) 
expressed in A and a vector v (also expressed in A) that is approximately in the 
desired direction for a second one of B's unit vectors, u2 (!= u1). If v is not 
perpendicular to u1, no worries - we'll find a direction for u2 that is 
perpendicular to u1 and comes closest to v. The third vector u3 is +/- u1 X u2, 
as appropriate for a right-handed rotation matrix. **/
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setRotationFromTwoAxes(const UnitVec3P& uveci, const CoordinateAxis& axisi, 
                       const Vec3P& vecjApprox, const CoordinateAxis& axisjApprox );

/** Set this Rotation_ to represent a rotation characterized by subsequent 
rotations of: +v[0] about the body frame's X axis, followed by a rotation of 
+v[1] about the body frame's NEW Y axis. See Kane, Spacecraft Dynamics, pg. 423, 
body-three: 1-2-3, but the last rotation is zero. **/
void setRotationToBodyFixedXY( const Vec2P& v)   
{ setRotationFromTwoAnglesTwoAxes(BodyRotationSequence,
                                  v[0],XAxis, v[1],YAxis); }

/** Set this Rotation_ to represent a rotation characterized by subsequent 
rotations of: +v[0] about the body frame's X axis, followed by a rotation of 
+v[1] about the body frame's NEW Y axis, followed by a rotation of +v[2] about 
the body frame's NEW (twice rotated) Z axis. See Kane, Spacecraft Dynamics, 
pg. 423, body-three: 1-2-3. **/
void setRotationToBodyFixedXYZ( const Vec3P& v)  
{ setRotationFromThreeAnglesThreeAxes(BodyRotationSequence, 
                                      v[0],XAxis, v[1],YAxis, v[2],ZAxis ); }
/** Given cosines and sines (in that order) of three angles, set this
%Rotation matrix to the body-fixed 1-2-3 sequence of those angles.
Cost is 18 flops. **/
void setRotationToBodyFixedXYZ(const Vec3P& c, const Vec3P& s) {
    Mat33P& R = *this;
    const RealP s0s1=s[0]*s[1], s2c0=s[2]*c[0], c0c2=c[0]*c[2], nc1= -c[1];

    R = Mat33P(     c[1]*c[2]      ,         s[2]*nc1       ,    s[1]  ,
                s2c0 + s0s1*c[2]   ,     c0c2 - s0s1*s[2]   , s[0]*nc1 ,
                s[0]*s[2] - s[1]*c0c2 ,  s[0]*c[2] + s[1]*s2c0 , c[0]*c[1] );
}
//@}

//------------------------------------------------------------------------------
/** @name Operators and Arithmetic **/
//@{
/** Transpose operator. For an orthogonal matrix like this one, 
transpose is the same thing as inversion. **/
const InverseRotation_<P>&  operator~() const  { return invert(); }
/** Transpose operator. For an orthogonal matrix like this one, 
transpose is the same thing as inversion. **/
InverseRotation_<P>&        operator~()        { return updInvert(); }

/** Transpose. For an orthogonal matrix like this one, 
transpose is the same thing as inversion. Overrides the base class 
transpose method. **/
const InverseRotation_<P>&  transpose() const  { return invert(); }
/** Transpose. For an orthogonal matrix like this one, 
transpose is the same thing as inversion. Overrides the base class 
transpose method. **/
InverseRotation_<P>&        updTranspose()     { return updInvert(); }

/** Convert from Rotation_ to InverseRotation_ (no cost). Overrides base 
class invert() method. **/
const InverseRotation_<P>& invert() const  
{ return *reinterpret_cast<const InverseRotation_<P>*>(this); }
/** Convert from Rotation_ to writable InverseRotation_ (no cost). **/
InverseRotation_<P>& updInvert()     
{ return *reinterpret_cast<InverseRotation_<P>*>(this); }

/** In-place composition of Rotation matrices. **/
inline Rotation_&  operator*=( const Rotation_<P>& R );
/** In-place composition of Rotation matrices. **/
inline Rotation_&  operator*=( const InverseRotation_<P>& );

/** In-place composition of Rotation matrices. **/
inline Rotation_&  operator/=( const Rotation_<P>& R );
/** In-place composition of Rotation matrices. **/
inline Rotation_&  operator/=( const InverseRotation_<P>& );

/** This is the fastest way to form the product qdot=N_P*w_PB for a 
body-fixed XYZ sequence where angular velocity of child in parent is
expected to be expressed in the parent. Here we assume you have
previously calculated sincos(qx), sincos(qy), and 1/cos(qy).
Cost is 10 flops, faster even than the 15 it would take if you had saved
N_P and then formed the N_P*w_PB product explicitly. **/
static Vec3P multiplyByBodyXYZ_N_P(const Vec2P& cosxy,
                                    const Vec2P& sinxy,
                                    RealP        oocosy,
                                    const Vec3P& w_PB)
{
    const RealP s0 = sinxy[0], c0 = cosxy[0];
    const RealP s1 = sinxy[1];
    const RealP w0 = w_PB[0], w1 = w_PB[1], w2 = w_PB[2];

    const RealP t = (s0*w1-c0*w2)*oocosy;
    return Vec3P( w0 + t*s1, c0*w1 + s0*w2, -t ); // qdot
}

/** This is the fastest way to form the product v_P=~N_P*q=~(~q*N_P); 
see the untransposed method multiplyByBodyXYZ_N_P() for information.
Cost is 9 flops. **/
static Vec3P multiplyByBodyXYZ_NT_P(const Vec2P& cosxy,
                                    const Vec2P& sinxy,
                                    RealP        oocosy,
                                    const Vec3P& q)
{
    const RealP s0 = sinxy[0], c0 = cosxy[0];
    const RealP s1 = sinxy[1];
    const RealP q0 = q[0], q1 = q[1], q2 = q[2];

    const RealP t = (q0*s1-q2) * oocosy;
    return Vec3P( q0, c0*q1 + t*s0, s0*q1 - t*c0 ); // v_P
}

/** Fastest way to form the product w_PB=NInv_P*qdot. This is never
singular. Cost is 9 flops. **/
static Vec3P multiplyByBodyXYZ_NInv_P(const Vec2P& cosxy,
                                      const Vec2P& sinxy,
                                      const Vec3P& qdot)
{
    const RealP s0 = sinxy[0], c0 = cosxy[0];
    const RealP s1 = sinxy[1], c1 = cosxy[1];
    const RealP q0 = qdot[0], q1 = qdot[1], q2 = qdot[2];
    const RealP c1q2 = c1*q2;

    return Vec3P( q0 + s1*q2,           // w_PB
                  c0*q1 - s0*c1q2, 
                  s0*q1 + c0*c1q2 );
}

/** Fastest way to form the product q=~NInv_P*v_P=~(~v_P*NInv_P). 
This is never singular. Cost is 10 flops. **/
static Vec3P multiplyByBodyXYZ_NInvT_P(const Vec2P& cosxy,
                                       const Vec2P& sinxy,
                                       const Vec3P& v_P)
{
    const RealP s0 = sinxy[0], c0 = cosxy[0];
    const RealP s1 = sinxy[1], c1 = cosxy[1];
    const RealP w0 = v_P[0], w1 = v_P[1], w2 = v_P[2];

    return Vec3P( w0,                           // qdot-like
                  c0*w1 + s0*w2,
                  s1*w0 - s0*c1*w1 + c0*c1*w2);
}
//@}

//------------------------------------------------------------------------------
/** @name Accessors **/
//@{
/** Return a reference to the ith row of this %Rotation matrix as 
a UnitRow3. **/
const RowType&  row( int i ) const         
{ return reinterpret_cast<const RowType&>(asMat33()[i]); }
/** Same as row(i) but nicer to look at. **/
const RowType&  operator[]( int i ) const  { return row(i); }

/** Return a reference to the jth column of this %Rotation matrix as
a UnitVec3. **/
const ColType&  col( int j ) const         
{ return reinterpret_cast<const ColType&>(asMat33()(j)); }
/** Same as col(j) but nicer to look at. **/
const ColType&  operator()( int j ) const  { return col(j); }

/** Return col(0) of this %Rotation matrix as a UnitVec3. **/
const ColType&  x() const                  { return col(0); }
/** Return col(1) of this %Rotation matrix as a UnitVec3. **/
const ColType&  y() const                  { return col(1); }
/** Return col(2) of this %Rotation matrix as a UnitVec3. **/
const ColType&  z() const                  { return col(2); }

/** Given a CoordinateAxis (XAxis,YAxis, or ZAxis) return a reference to
the corresponding column of this %Rotation matrix. The result is equivalent
to multiplying R_AB*v_B where v_B is [1,0,0],[0,1,0], or [0,0,1], which would
cost 15 flops, but requires no computation. **/
const ColType& getAxisUnitVec(CoordinateAxis axis) const 
{   return col(axis); }

/** Given a CoordinateDirection (+/-XAxis, etc.) return a unit vector in that
direction. The result is equivalent to multiplying R_AB*v_B where v_B is 
[+/-1,0,0], [0,+/-1,0], or [0,0,+/-1], which would cost 15 flops, but this
method requires at most 3 flops. **/
const UnitVec<P,1> getAxisUnitVec(CoordinateDirection dir) const {
    const ColType& axDir = getAxisUnitVec(dir.getAxis());
    return dir.getDirection() > 0 ? UnitVec<P,1>( axDir) 
                                  : UnitVec<P,1>(-axDir); // cheap 
}
//@}

//------------------------------------------------------------------------------
/** @name Calculations **/
//@{
/** Given Euler angles q forming a body-fixed X-Y-Z sequence return the block 
N_B of the system N matrix such that qdot=N_B(q)*w_PB_B where w_PB_B is the 
angular velocity of B in P EXPRESSED IN *B*!!! Note that N_B=N_P*R_PB. This 
matrix will be singular if Y (q[1]) gets near 90 degrees!

@note This version is very expensive because it has to calculate sines and 
      cosines. If you already have those, use the alternate form of this method.

Cost: about 100 flops for sin/cos plus 12 to calculate N_B.
@see Kane's Spacecraft Dynamics, page 427, body-three: 1-2-3. **/
static Mat33P calcNForBodyXYZInBodyFrame(const Vec3P& q) {
    // Note: q[0] is not referenced so we won't waste time calculating
    // its cosine and sine here.
    return calcNForBodyXYZInBodyFrame
        (Vec3P(0, std::cos(q[1]), std::cos(q[2])),
        Vec3P(0, std::sin(q[1]), std::sin(q[2])));
}

/** This faster version of calcNForBodyXYZInBodyFrame() assumes you have 
already calculated the cosine and sine of the three q's. Note that we 
only look at the cosines and sines of q[1] and q[2]; q[0] does not 
matter so you don't have to fill in the 0'th element of cq and sq.
Cost is one divide plus 6 flops, say 12 flops. **/
static Mat33P calcNForBodyXYZInBodyFrame(const Vec3P& cq, const Vec3P& sq) {
    const RealP s1 = sq[1], c1 = cq[1];
    const RealP s2 = sq[2], c2 = cq[2];
    const RealP ooc1  = 1/c1;
    const RealP s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

    return Mat33P(    c2oc1  , -s2oc1  , 0,
                        s2   ,    c2   , 0,
                    -s1*c2oc1 , s1*s2oc1, 1 );
}

/** Given Euler angles q forming a body-fixed X-Y-Z (123) sequence return the 
block N_P of the system N matrix such that qdot=N_P(q)*w_PB where w_PB is the 
angular velocity of B in P expressed in P (not the convention that Kane uses, 
where angular velocities are expressed in the outboard body B). Note that 
N_P = N_B*~R_PB. This matrix will be singular if Y (q[1]) gets near 90 degrees!

@note This version is very expensive because it has to calculate sines and 
      cosines. If you already have those, use the alternate form of this method.

Cost: about 100 flops for sin/cos plus 12 to calculate N_P. **/
static Mat33P calcNForBodyXYZInParentFrame(const Vec3P& q) {
    // Note: q[2] is not referenced so we won't waste time calculating
    // its cosine and sine here.
    return calcNForBodyXYZInParentFrame
        (Vec3P(std::cos(q[0]), std::cos(q[1]), 0),
        Vec3P(std::sin(q[0]), std::sin(q[1]), 0));
}

/** This faster version of calcNForBodyXYZInParentFrame() assumes you have 
already calculated the cosine and sine of the three q's. Note that we 
only look at the cosines and sines of q[0] and q[1]; q[2] does not 
matter so you don't have to fill in the 3rd element of cq and sq.
Cost is one divide plus 6 flops, say 12 flops.
@see Paul Mitiguy **/
static Mat33P calcNForBodyXYZInParentFrame(const Vec3P& cq, const Vec3P& sq) {
    const RealP s0 = sq[0], c0 = cq[0];
    const RealP s1 = sq[1], c1 = cq[1];
    const RealP ooc1  = 1/c1;
    const RealP s0oc1 = s0*ooc1, c0oc1 = c0*ooc1;

    return Mat33P( 1 , s1*s0oc1 , -s1*c0oc1,
                    0 ,    c0    ,    s0,
                    0 ,  -s0oc1  ,  c0oc1 );
}

/** Given Euler angles forming a body-fixed X-Y-Z (123) sequence q, and 
their time derivatives qdot, return the block of the NDot matrix such 
that qdotdot=N(q)*wdot + NDot(q,u)*w where w is the angular velocity 
of B in P EXPRESSED IN *B*!!! This matrix will be singular if Y (q[1]) 
gets near 90 degrees! See calcNForBodyXYZInBodyFrame() for the matrix 
we're differentiating here.
@note This version is very expensive because it has to calculate sines
      and cosines. If you already have those, use the alternate form
      of this method. **/
static Mat33P calcNDotForBodyXYZInBodyFrame
   (const Vec3P& q, const Vec3P& qdot) {
    // Note: q[0] is not referenced so we won't waste time calculating
    // its cosine and sine here.
    return calcNDotForBodyXYZInBodyFrame
        (Vec3P(0, std::cos(q[1]), std::cos(q[2])),
        Vec3P(0, std::sin(q[1]), std::sin(q[2])),
        qdot);
}

/** This faster version of calcNDotForBodyXYZInBodyFrame() assumes you 
have already calculated the cosine and sine of the three q's. Note 
that we only look at the cosines and sines of q[1] and q[2]; q[0] does 
not matter so you don't have to fill in the 0'th element of cq and sq.
Cost is one divide plus 21 flops, say 30 flops. **/
static Mat33P calcNDotForBodyXYZInBodyFrame
   (const Vec3P& cq, const Vec3P& sq, const Vec3P& qdot) 
{
    const RealP s1 = sq[1], c1 = cq[1];
    const RealP s2 = sq[2], c2 = cq[2];
    const RealP ooc1  = 1/c1;
    const RealP s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

    const RealP t = qdot[1]*s1*ooc1;
    const RealP a = t*s2oc1 + qdot[2]*c2oc1; // d/dt s2oc1
    const RealP b = t*c2oc1 - qdot[2]*s2oc1; // d/dt c2oc1

    return Mat33P(       b             ,        -a         , 0,
                        qdot[2]*c2       ,    -qdot[2]*s2    , 0,
                    -(s1*b + qdot[1]*c2) , s1*a + qdot[1]*s2 , 0 );
}

/** Given Euler angles forming a body-fixed X-Y-Z (123) sequence q, and 
their time derivatives qdot, return the block of the NDot matrix such 
that qdotdot=N(q)*wdot + NDot(q,u)*w where w is the angular velocity of
B in P expressed in P. This matrix will be singular if Y (q[1]) gets
near 90 degrees! See calcNForBodyXYZInParentFrame() for the matrix 
we're differentiating here.
@note This version is very expensive because it has to calculate sines
      and cosines. If you already have those, use the alternate form
      of this method. **/
static Mat33P calcNDotForBodyXYZInParentFrame
   (const Vec3P& q, const Vec3P& qdot) {
    // Note: q[2] is not referenced so we won't waste time calculating
    // its cosine and sine here.
    const RealP cy = std::cos(q[1]); // cos(y)
    return calcNDotForBodyXYZInParentFrame
        (Vec2P(std::cos(q[0]), cy), 
        Vec2P(std::sin(q[0]), std::sin(q[1])),
        1/cy, qdot);
}

/** This faster version of calcNDotForBodyXYZInParentFrame() assumes you 
have already calculated the cosine and sine of the three q's. Note that
we only look at the cosines and sines of q[0] and q[1]. Cost is 21 flops. **/
static Mat33P calcNDotForBodyXYZInParentFrame
   (const Vec2P& cq, const Vec2P& sq, RealP ooc1, const Vec3P& qdot) {
    const RealP s0 = sq[0], c0 = cq[0];
    const RealP s1 = sq[1];
    const RealP s0oc1 = s0*ooc1, c0oc1 = c0*ooc1;

    const RealP t = qdot[1]*s1*ooc1;
    const RealP a = t*s0oc1 + qdot[0]*c0oc1; // d/dt s0oc1
    const RealP b = t*c0oc1 - qdot[0]*s0oc1; // d/dt c0oc1

    return Mat33P( 0,  s1*a + qdot[1]*s0, -(s1*b + qdot[1]*c0), 
                   0,    -qdot[0]*s0    ,     qdot[0]*c0      ,
                   0,        -a         ,         b            );
}

/** Inverse of routine calcNForBodyXYZInBodyFrame(). Return the inverse 
NInv_B of the N_B block computed above, such that w_PB_B=NInv_B(q)*qdot
where w_PB_B is the angular velocity of B in P EXPRESSED IN *B*!!! 
(Kane's convention.) Note that NInv_B=~R_PB*NInv_P. This matrix is 
never singular.
@note This version is very expensive because it has to calculate sines
      and cosines. If you already have those, use the alternate form
      of this method. **/
static Mat33P calcNInvForBodyXYZInBodyFrame(const Vec3P& q) {
    // Note: q[0] is not referenced so we won't waste time calculating
    // its cosine and sine here.
    return calcNInvForBodyXYZInBodyFrame
       (Vec3P(0, std::cos(q[1]), std::cos(q[2])),
        Vec3P(0, std::sin(q[1]), std::sin(q[2])));
}

/** This faster version of calcNInvForBodyXYZInBodyFrame() assumes you have
already calculated the cosine and sine of the three q's. Note that we only look 
at the cosines and sines of q[1] and q[2]; q[0] does not matter so you don't 
have to fill in the 0'th element of cq and sq. Cost is 3 flops. **/
static Mat33P calcNInvForBodyXYZInBodyFrame
   (const Vec3P& cq, const Vec3P& sq) {
    const RealP s1 = sq[1], c1 = cq[1];
    const RealP s2 = sq[2], c2 = cq[2];

    return Mat33P( c1*c2 ,  s2 , 0 ,
                  -c1*s2 ,  c2 , 0 ,
                    s1   ,  0  , 1 );
}

/** Inverse of the above routine. Return the inverse NInv_P of the N_P 
block computed above, such that w_PB=NInv_P(q)*qdot where w_PB is the 
angular velocity of B in P (expressed in P). Note that 
NInv_P=R_PB*NInv_B. This matrix is never singular.
@note This version is very expensive because it has to calculate sines
      and cosines. If you already have those, use the alternate form
      of this method. **/
static Mat33P calcNInvForBodyXYZInParentFrame(const Vec3P& q) {
    // Note: q[0] is not referenced so we won't waste time calculating
    // its cosine and sine here.
    return calcNInvForBodyXYZInParentFrame
       (Vec3P(std::cos(q[0]), std::cos(q[1]), 0),
        Vec3P(std::sin(q[0]), std::sin(q[1]), 0));
}

/** This faster version of calcNInvForBodyXYZInParentFrame() assumes you have 
already calculated the cosine and sine of the three q's. Note that we only look 
at the cosines and sines of q[0] and q[1]; q[2] does not matter so you don't 
have to fill in the 3rd element of cq and sq. Cost is 3 flops. **/
static Mat33P calcNInvForBodyXYZInParentFrame
    (const Vec3P& cq, const Vec3P& sq) {
    const RealP s0 = sq[0], c0 = cq[0];
    const RealP s1 = sq[1], c1 = cq[1];

    return Mat33P( 1 ,  0  ,   s1   ,
                    0 ,  c0 , -s0*c1 ,
                    0 ,  s0 ,  c0*c1 );
}

/** Given a possibly unnormalized quaternion q, calculate the 4x3 matrix N which
maps angular velocity w to quaternion derivatives qdot. We expect the angular 
velocity in the parent frame, i.e. w==w_PB_P. We don't normalize, so N=|q|N' 
where N' is the normalized version. Cost is 7 flops. **/
static Mat43P calcUnnormalizedNForQuaternion(const Vec4P& q) {
    const Vec4P e = q/2;
    const RealP ne1 = -e[1], ne2 = -e[2], ne3 = -e[3];
    return Mat43P( ne1,  ne2,  ne3,
                    e[0], e[3], ne2,
                    ne3,  e[0], e[1],
                    e[2], ne1,  e[0]);
}

/** Given the time derivative qdot of a possibly unnormalized quaternion q, 
calculate the 4x3 matrix NDot which is the time derivative of the matrix N as 
described in calcUnnormalizedNForQuaternion(). Note that NDot = d/dt N = 
d/dt (|q|N') = |q|(d/dt N'), where N' is the normalized matrix, since the length
of the quaternion should be a constant. Cost is 7 flops. **/
static Mat43P calcUnnormalizedNDotForQuaternion(const Vec4P& qdot) {
    const Vec4P ed = qdot/2;
    const RealP ned1 = -ed[1], ned2 = -ed[2], ned3 = -ed[3];
    return Mat43P( ned1,  ned2,  ned3,
                    ed[0], ed[3], ned2,
                    ned3,  ed[0], ed[1],
                    ed[2], ned1,  ed[0]);
}

/** Given a (possibly unnormalized) quaternion q, calculate the 3x4 matrix
NInv (= N^-1) which maps quaternion derivatives qdot to angular velocity w, 
where the angular velocity is in the parent frame, i.e. w==w_PB_P. Note: when 
the quaternion is not normalized, this is not precisely the (pseudo)inverse of 
N. inv(N)=inv(N')/|q| but we're returning |q|*inv(N')=|q|^2*inv(N). That is, 
NInv*N =|q|^2*I, which is I if the original q was normalized. 
(Note: N*NInv != I, not even close.) Cost is 7 flops. **/
static Mat34P calcUnnormalizedNInvForQuaternion(const Vec4P& q) {
    const Vec4P e = 2*q;
    const RealP ne1 = -e[1], ne2 = -e[2], ne3 = -e[3];
    return Mat34P(ne1, e[0], ne3,  e[2],
                    ne2, e[3], e[0], ne1,
                    ne3, ne2,  e[1], e[0]);
}
//@}

//------------------------------------------------------------------------------
/** @name Conversions **/
//@{
/** Conversion from Rotation to its base class Mat33. Note: asMat33 is more 
efficient than toMat33() (no copy), but you have to know the internal 
layout. **/
const Mat33P&  asMat33() const  { return *static_cast<const Mat33P*>(this); }
/** Conversion from Rotation to its base class Mat33. Note: asMat33 is more 
efficient than toMat33() (no copy), but you have to know the internal 
layout. **/
Mat33P         toMat33() const  { return asMat33(); }

/** Perform an efficient transform of a symmetric matrix that must be 
re-expressed with a multiply from both left and right, such as an inertia 
matrix. Details: assuming this Rotation is R_AB, and given a symmetric dyadic 
matrix S_BB expressed in B, we can reexpress it in A using S_AA=R_AB*S_BB*R_BA. 
The matrix should be one that is formed as products of vectors expressed in A, 
such as inertia, gyration or covariance matrices. This can be done efficiently 
exploiting properties of R (orthogonal) and S (symmetric). 
Total cost is 57 flops. **/
SimTK_SimTKCOMMON_EXPORT SymMat33P 
reexpressSymMat33(const SymMat33P& S_BB) const;

// Converts rotation matrix to one or two or three orientation angles.
// Note:  The result is most meaningful if the Rotation_ matrix is one that can 
// be produced by such a sequence.
// Use1:  someRotation.convertOneAxisRotationToOneAngle( XAxis );
// Use2:  someRotation.convertTwoAxesRotationToTwoAngles
//                                (SpaceRotationSequence, YAxis, ZAxis );
// Use3:  someRotation.convertThreeAxesRotationToThreeAngles
//                                (SpaceRotationSequence, ZAxis, YAxis, XAxis );
// Use4:  someRotation.convertRotationToAngleAxis();   
//        Return: [angleInRadians, unitVectorX, unitVectorY, unitVectorZ].

/** Converts rotation matrix to a single orientation angle. Note:  The result is
most meaningful if the Rotation_ matrix is one that can be produced by such 
a sequence. **/
SimTK_SimTKCOMMON_EXPORT RealP  
convertOneAxisRotationToOneAngle( const CoordinateAxis& axis1 ) const;

/** Converts rotation matrix to two orientation angles. Note:  The result is 
most meaningful if the Rotation_ matrix is one that can be produced by such 
a sequence. **/
SimTK_SimTKCOMMON_EXPORT Vec2P  
convertTwoAxesRotationToTwoAngles(BodyOrSpaceType bodyOrSpace, 
                                  const CoordinateAxis& axis1, 
                                  const CoordinateAxis& axis2) const;

/** Converts rotation matrix to three orientation angles. Note:  The result is 
most meaningful if the Rotation_ matrix is one that can be produced by such 
a sequence. **/
SimTK_SimTKCOMMON_EXPORT Vec3P  
convertThreeAxesRotationToThreeAngles
   (BodyOrSpaceType bodyOrSpace, const CoordinateAxis& axis1, 
    const CoordinateAxis& axis2, const CoordinateAxis& axis3 ) const;

/** Converts rotation matrix to an equivalent quaternion in canonical form
(meaning its scalar element is nonnegative). This uses a robust,
singularity-free method due to Richard Spurrier. The cost is about 40 flops.

@par Reference
Spurrier, R.A., "Comment on 'Singularity-Free Extraction of a Quaternion 
from a Direction-Cosine Matrix'", J. Spacecraft and Rockets, 15(4):255, 
1977. 

@see Quaternion_ **/
SimTK_SimTKCOMMON_EXPORT QuaternionP convertRotationToQuaternion() const;

/** Converts rotation matrix to an equivalent angle-axis representation in
canonicalized form. The result (a,v) is returned packed into a Vec4
[a vx vy vz], with -Pi < a <= Pi and |v|=1. Cost is about 140 flops. 

If the rotation angle is zero (or very very close to zero) then the returned
unit vector is arbitrary. 

@par Theory
Euler's Rotation Theorem (1776) guarantees that any rigid body rotation is 
equivalent to a rotation by an angle about a fixed axis. This method finds
such an angle and axis. Numerically, this is a very tricky computation to
get correct in all cases. We use Spurrier's method to obtain a 
numerically-robust quaternion equivalent to this rotation matrix, then 
carefully extract and canonicalize the angle-axis form from the quaternion.

@see convertRotationToQuaternion()
**/
Vec4P convertRotationToAngleAxis() const  
{ return convertRotationToQuaternion().convertQuaternionToAngleAxis(); }

/** A convenient special case of convertTwoAxesRotationToTwoAngles(). **/
Vec2P convertRotationToBodyFixedXY() const   
{ return convertTwoAxesRotationToTwoAngles(BodyRotationSequence,XAxis,YAxis); }
/** A convenient special case of convertThreeAxesRotationToThreeAngles(). **/
Vec3P convertRotationToBodyFixedXYZ() const  
{ return convertThreeAxesRotationToThreeAngles( BodyRotationSequence, 
                                                XAxis, YAxis, ZAxis ); }

/** Given Euler angles forming a body-fixed 3-2-1 sequence, and the relative
angular velocity vector of B in the parent frame, *BUT EXPRESSED IN THE BODY 
FRAME*, return the Euler angle derivatives. You are dead if q[1] gets near 
90 degrees! See Kane's Spacecraft Dynamics, page 428, body-three: 3-2-1. **/
static Vec3P convertAngVelToBodyFixed321Dot(const Vec3P& q, const Vec3P& w_PB_B) {
    const RealP s1 = std::sin(q[1]), c1 = std::cos(q[1]);
    const RealP s2 = std::sin(q[2]), c2 = std::cos(q[2]);
    const RealP ooc1 = RealP(1)/c1;
    const RealP s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

    const Mat33P E( 0,    s2oc1  ,  c2oc1  ,
                    0,      c2   ,   -s2   ,
                    1,  s1*s2oc1 , s1*c2oc1 );
    return E * w_PB_B;
}

/** Inverse of convertAngVelToBodyFixed321Dot. Returned angular velocity is B in
P, expressed in *B*: w_PB_B. **/
static Vec3P convertBodyFixed321DotToAngVel(const Vec3P& q, const Vec3P& qd) {
    const RealP s1 = std::sin(q[1]), c1 = std::cos(q[1]);
    const RealP s2 = std::sin(q[2]), c2 = std::cos(q[2]);

    const Mat33P Einv(  -s1  ,  0  ,  1 ,
                        c1*s2 ,  c2 ,  0 ,
                        c1*c2 , -s2 ,  0 );
    return Einv*qd;
}

// TODO: sherm: is this right? Warning: everything is measured in the
// *PARENT* frame, but angular velocities and accelerations are
// expressed in the *BODY* frame.
// TODO: this is not an efficient way to do this computation.
/** Caution: needs testing. **/
static Vec3P convertAngVelDotToBodyFixed321DotDot
    (const Vec3P& q, const Vec3P& w_PB_B, const Vec3P& wdot_PB_B)
{
    const RealP s1 = std::sin(q[1]), c1 = std::cos(q[1]);
    const RealP s2 = std::sin(q[2]), c2 = std::cos(q[2]);
    const RealP ooc1  = 1/c1;
    const RealP s2oc1 = s2*ooc1, c2oc1 = c2*ooc1, s1oc1 = s1*ooc1;

    const Mat33P E( 0 ,   s2oc1  ,  c2oc1  ,
                    0 ,     c2   ,   -s2   ,
                    1 , s1*s2oc1 , s1*c2oc1 );
    const Vec3P qdot = E * w_PB_B;

    const RealP t = qdot[1]*s1oc1;
    const RealP a = t*s2oc1 + qdot[2]*c2oc1; // d/dt s2oc1
    const RealP b = t*c2oc1 - qdot[2]*s2oc1; // d/dt c2oc1

    const Mat33P Edot( 0 ,       a           ,         b         ,
                        0 ,   -qdot[2]*s2     ,    -qdot[2]*c2    ,
                        0 , s1*a + qdot[1]*s2 , s1*b + qdot[1]*c2 );

    return E*wdot_PB_B + Edot*w_PB_B;
}

/** Given Euler angles forming a body-fixed X-Y-Z (123) sequence, and the 
relative angular velocity vector w_PB_B of B in the parent frame, 
<em>BUT EXPRESSED IN THE BODY FRAME</em>, return the Euler angle 
derivatives. You are dead if q[1] gets near 90 degrees!
@note This version is very expensive because it has to calculate sines
      and cosines. If you already have those, use the alternate form
      of this method.
@see Kane's Spacecraft Dynamics, page 427, body-three: 1-2-3. **/
static Vec3P convertAngVelInBodyFrameToBodyXYZDot
    (const Vec3P& q, const Vec3P& w_PB_B) {  
    return convertAngVelInBodyFrameToBodyXYZDot
        (Vec3P(0, std::cos(q[1]), std::cos(q[2])),
        Vec3P(0, std::sin(q[1]), std::sin(q[2])),
        w_PB_B); 
}

/** This faster version of convertAngVelInBodyFrameToBodyXYZDot() assumes 
you have already calculated the cosine and sine of the three q's. Note
that we only look at the cosines and sines of q[1] and q[2]; q[0] does 
not matter so you don't have to fill in the 0'th element of cq and sq.
Cost is 27 flops. **/
//TODO: reimplement
static Vec3P convertAngVelInBodyFrameToBodyXYZDot
   (const Vec3P& cq, const Vec3P& sq, const Vec3P& w_PB_B) 
{   return calcNForBodyXYZInBodyFrame(cq,sq)*w_PB_B; }

/** Inverse of the above routine. Returned angular velocity is B in P,
expressed in *B*: w_PB_B.
@note This version is very expensive because it has to calculate sines
      and cosines. If you already have those, use the alternate form
      of this method. **/
static Vec3P convertBodyXYZDotToAngVelInBodyFrame
   (const Vec3P& q, const Vec3P& qdot) {   
        return convertBodyXYZDotToAngVelInBodyFrame
                   (Vec3P(0, std::cos(q[1]), std::cos(q[2])),
                    Vec3P(0, std::sin(q[1]), std::sin(q[2])),
                    qdot); 
}

/** This faster version of convertBodyXYZDotToAngVelInBodyFrame() assumes
you have already calculated the cosine and sine of the three q's. Note 
that we only look at the cosines and sines of q[1] and q[2]; q[0] does 
not matter so you don't have to fill in the 0'th element of cq and sq.
Cost is 18 flops. **/
// TODO: reimplement
static Vec3P convertBodyXYZDotToAngVelInBodyFrame
    (const Vec3P& cq, const Vec3P& sq, const Vec3P& qdot) 
{   return calcNInvForBodyXYZInBodyFrame(cq,sq)*qdot; }

// TODO: sherm: is this right?
/** Warning: everything is measured in the *PARENT* frame, but has to be
expressed in the *BODY* frame.
@note This version is very expensive because it has to calculate sines and
      cosines. If you already have those, use the alternate form of this method.
Caution: needs testing. **/
static Vec3P convertAngVelDotInBodyFrameToBodyXYZDotDot
   (const Vec3P& q, const Vec3P& w_PB_B, const Vec3P& wdot_PB_B)
{
    // Note: q[0] is not referenced so we won't waste time calculating
    // its cosine and sine here.
    return convertAngVelDotInBodyFrameToBodyXYZDotDot
               (Vec3P(0, std::cos(q[1]), std::cos(q[2])),
                Vec3P(0, std::sin(q[1]), std::sin(q[2])),
                w_PB_B, wdot_PB_B);
}

/** This faster version of convertAngVelDotInBodyFrameToBodyXYZDotDot() 
assumes you have already calculated the cosine and sine of the three 
q's. Note that we only look at the cosines and sines of q[1] and q[2]; 
q[0] does not matter so you don't have to fill in the 0'th element of 
cq and sq. Cost is about 93 flops. **/
// TODO: reimplement
static Vec3P convertAngVelDotInBodyFrameToBodyXYZDotDot
   (const Vec3P& cq, const Vec3P& sq, 
    const Vec3P& w_PB_B, const Vec3P& wdot_PB_B)
{
    const Mat33P N    = calcNForBodyXYZInBodyFrame(cq,sq);         // ~12 flops
    const Vec3P  qdot = N * w_PB_B;                                //  15 flops
    const Mat33P NDot = calcNDotForBodyXYZInBodyFrame(cq,sq,qdot); // ~30 flops

    return N*wdot_PB_B + NDot*w_PB_B;                              //  33 flops
}

/** Given a possibly unnormalized quaternion (0th element is the scalar) and the
relative angular velocity vector of B in its parent, expressed in the *PARENT*, 
return the quaternion derivatives. This is never singular. Cost is 27 flops. **/
static Vec4P convertAngVelToQuaternionDot(const Vec4P& q, const Vec3P& w_PB_P) {
    return calcUnnormalizedNForQuaternion(q)*w_PB_P;
}

/** Inverse of the above routine. Returned AngVel is expressed in the *PARENT* 
frame: w_PB_P. Cost is 28 flops. **/
static Vec3P convertQuaternionDotToAngVel(const Vec4P& q, const Vec4P& qdot) {
    return calcUnnormalizedNInvForQuaternion(q)*qdot;
}

/** We want to differentiate qdot=N(q)*w to get qdotdot=N*b+NDot*w where b is 
angular acceleration wdot. Note that NDot=NDot(qdot), but it is far better to 
calculate the matrix-vector product NDot(N*w)*w directly rather than calculate 
NDot separately. That gives <pre>
    NDot*w = -(w^2)/4 * q
</pre> Cost is 41 flops. **/
static Vec4P convertAngVelDotToQuaternionDotDot
    (const Vec4P& q, const Vec3P& w_PB, const Vec3P& b_PB)
{
    const Mat43P N     = calcUnnormalizedNForQuaternion(q); //  7 flops
    const Vec4P  Nb    = N*b_PB;                            // 20 flops
    const Vec4P  NDotw = RealP(-.25)*w_PB.normSqr()*q;      // 10 flops
    return Nb + NDotw;                                      //  4 flops
}

/** Calculate first time derivative qdot of body-fixed XYZ Euler angles q
given sines and cosines of the Euler angles and the angular velocity 
w_PB of child B in parent P, expressed in P. Cost is 10 flops.

Theory: calculate qdot=N_P(q)*w_PB using multiplyByBodyXYZ_N_P().
@see multiplyByBodyXYZ_N_P() **/
static Vec3P convertAngVelInParentToBodyXYZDot
    (const Vec2P& cosxy,  ///< cos(qx), cos(qy)
    const Vec2P& sinxy,  ///< sin(qx), sin(qy)
    RealP        oocosy, ///< 1/cos(qy)
    const Vec3P& w_PB)   ///< angular velocity of B in P, exp. in P
{
    return multiplyByBodyXYZ_N_P(cosxy,sinxy,oocosy,w_PB);
}

/** Calculate second time derivative qdotdot of body-fixed XYZ Euler 
angles q given sines and cosines of the Euler angles, the first 
derivative qdot and the angular acceleration b_PB of child B in 
parent P, expressed in P. Cost is 22 flops.

Theory: we have qdot=N_P*w_PB, which we differentiate in P to 
get qdotdot=N_P*b_PB + NDot_P*w_PB. Note that NDot_P=NDot_P(q,qdot) 
and w_PB=NInv_P*qdot (because N_P is invertible). We can then rewrite
qdotdot=N_P*b_PB + NDot_P*(NInv_P*qdot) which can be calculated very 
efficiently. The second term is just an acceleration remainder term
quadratic in qdot. **/
static Vec3P convertAngAccInParentToBodyXYZDotDot
    (const Vec2P& cosxy,  ///< cos(qx), cos(qy)
    const Vec2P& sinxy,  ///< sin(qx), sin(qy)
    RealP        oocosy, ///< 1/cos(qy)
    const Vec3P& qdot,   ///< previously calculated BodyXYZDot
    const Vec3P& b_PB)   ///< angular acceleration, a.k.a. wdot_PB
{
    const RealP s1 = sinxy[1], c1 = cosxy[1];
    const RealP q0 = qdot[0], q1 = qdot[1], q2 = qdot[2];

    // 10 flops
    const Vec3P Nb = multiplyByBodyXYZ_N_P(cosxy,sinxy,oocosy,b_PB);

    const RealP q1oc1 = q1*oocosy;
    const Vec3P NDotw((q0*s1-q2)*q1oc1,     //   NDot_P*w_PB
                        q0*q2*c1,            // = NDot_P*(NInv_P*qdot)
                        (q2*s1-q0)*q1oc1 );   // (9 flops)

    return Nb + NDotw; // 3 flops
}
//@}

//------------------------------------------------------------------------------
/** @name Queries **/
//@{
/** Return true if "this" Rotation is nearly identical to "R" within a specified
pointing angle error. **/
SimTK_SimTKCOMMON_EXPORT bool  
isSameRotationToWithinAngle(const Rotation_& R, RealP okPointingAngleErrorRads) 
                                                                        const;

/** Return true if "this" Rotation is nearly identical to "R" within machine
precision. **/
bool isSameRotationToWithinAngleOfMachinePrecision(const Rotation_& R) const       
{ return isSameRotationToWithinAngle( R, NTraits<P>::getSignificant() ); }

/** Returns maximum absolute difference between elements in "this" Rotation and
elements in "R". **/
RealP getMaxAbsDifferenceInRotationElements( const Rotation_& R ) const {            
    const Mat33P& A=asMat33(); const Mat33P& B=R.asMat33(); RealP maxDiff=0;  
    for( int i=0;  i<=2; i++ ) for( int j=0; j<=2; j++ ) {
        const RealP absDiff = std::abs(A[i][j] - B[i][j]);  
        if( absDiff > maxDiff ) maxDiff = absDiff; 
    }
    return maxDiff; 
} 

/** Returns true if each element of "this" Rotation is within epsilon of the
corresponding element of "R". **/
bool areAllRotationElementsSameToEpsilon(const Rotation_& R, RealP epsilon) const 
{ return getMaxAbsDifferenceInRotationElements(R) <= epsilon; }

/** Returns true if each element of "this" Rotation is within machine precision
of the corresponding element of "R". **/
bool areAllRotationElementsSameToMachinePrecision( const Rotation_& R ) const       
{ return areAllRotationElementsSameToEpsilon(R, NTraits<P>::getSignificant()); } 
//@}


private:
// This is only for the most trustworthy of callers, that is, methods of 
// the Rotation_ class.  There are a lot of ways for this NOT to be a 
// legitimate rotation matrix -- be careful!!
// Note that these are supplied in rows.
Rotation_( const RealP& xx, const RealP& xy, const RealP& xz,
            const RealP& yx, const RealP& yy, const RealP& yz,
            const RealP& zx, const RealP& zy, const RealP& zz )
:   Mat33P( xx,xy,xz, yx,yy,yz, zx,zy,zz ) {}

// These next methods are highly-efficient power-user methods. Read the 
// code to understand them.
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setTwoAngleTwoAxesBodyFixedForwardCyclicalRotation
   (RealP cosAngle1, RealP sinAngle1, const CoordinateAxis& axis1, 
    RealP cosAngle2, RealP sinAngle2, const CoordinateAxis& axis2 );
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setThreeAngleTwoAxesBodyFixedForwardCyclicalRotation
   (RealP cosAngle1, RealP sinAngle1, const CoordinateAxis& axis1, 
    RealP cosAngle2, RealP sinAngle2, const CoordinateAxis& axis2, 
    RealP cosAngle3, RealP sinAngle3 );
SimTK_SimTKCOMMON_EXPORT Rotation_&  
setThreeAngleThreeAxesBodyFixedForwardCyclicalRotation
   (RealP cosAngle1, RealP sinAngle1, const CoordinateAxis& axis1, 
    RealP cosAngle2, RealP sinAngle2, const CoordinateAxis& axis2, 
    RealP cosAngle3, RealP sinAngle3, const CoordinateAxis& axis3 );

// These next methods are highly-efficient power-user methods to convert 
// Rotation matrices to orientation angles.  Read the code to understand them.
SimTK_SimTKCOMMON_EXPORT Vec2P  
convertTwoAxesBodyFixedRotationToTwoAngles
   (const CoordinateAxis& axis1, const CoordinateAxis& axis2 ) const;
SimTK_SimTKCOMMON_EXPORT Vec3P  
convertTwoAxesBodyFixedRotationToThreeAngles
   (const CoordinateAxis& axis1, const CoordinateAxis& axis2 ) const;
SimTK_SimTKCOMMON_EXPORT Vec3P  
convertThreeAxesBodyFixedRotationToThreeAngles
   (const CoordinateAxis& axis1, const CoordinateAxis& axis2, 
    const CoordinateAxis& axis3 ) const;

//------------------------------------------------------------------------------
// These are obsolete names from a previous release, listed here so that 
// users will get a decipherable compilation error. (sherm 091101)
//------------------------------------------------------------------------------
private:
// REPLACED BY: calcNForBodyXYZInBodyFrame()
static Mat33P calcQBlockForBodyXYZInBodyFrame(const Vec3P& a)
{   return calcNForBodyXYZInBodyFrame(a); }
// REPLACED BY: calcNInvForBodyXYZInBodyFrame()
static Mat33P calcQInvBlockForBodyXYZInBodyFrame(const Vec3P& a)
{   return calcNInvForBodyXYZInBodyFrame(a); }
// REPLACED BY: calcUnnormalizedNForQuaternion()
static Mat<4,3,P> calcUnnormalizedQBlockForQuaternion(const Vec4P& q)
{   return calcUnnormalizedNForQuaternion(q); }
// REPLACED BY: calcUnnormalizedNInvForQuaternion()
static Mat<3,4,P> calcUnnormalizedQInvBlockForQuaternion(const Vec4P& q)
{   return calcUnnormalizedNInvForQuaternion(q); }
// REPLACED BY: convertAngVelInBodyFrameToBodyXYZDot
static Vec3P convertAngVelToBodyFixed123Dot(const Vec3P& q, const Vec3P& w_PB_B) 
{   return convertAngVelInBodyFrameToBodyXYZDot(q,w_PB_B); }
// REPLACED BY: convertBodyXYZDotToAngVelInBodyFrame
static Vec3P convertBodyFixed123DotToAngVel(const Vec3P& q, const Vec3P& qdot) 
{   return convertBodyXYZDotToAngVelInBodyFrame(q,qdot); }
// REPLACED BY: convertAngVelDotInBodyFrameToBodyXYZDotDot
static Vec3P convertAngVelDotToBodyFixed123DotDot
    (const Vec3P& q, const Vec3P& w_PB_B, const Vec3P& wdot_PB_B)
{   return convertAngVelDotInBodyFrameToBodyXYZDotDot(q,w_PB_B,wdot_PB_B); }

//------------------------------------------------------------------------------
// The following code is obsolete - it is here temporarily for backward 
// compatibility (Mitiguy 9/5/2007)
//------------------------------------------------------------------------------
private:
// These static methods are like constructors with friendlier names.
static Rotation_ zero() { return Rotation_(); }
static Rotation_ NaN()  { Rotation_ r;  r.setRotationToNaN();  return r; }

/// By zero we mean "zero rotation", i.e., an identity matrix.
Rotation_&  setToZero()            { return setRotationToIdentityMatrix(); }
Rotation_&  setToIdentityMatrix()  { return setRotationToIdentityMatrix(); }
Rotation_&  setToNaN()             { return setRotationToNaN(); }
static Rotation_  trustMe( const Mat33P& m )  { return Rotation_(m,true); }

// One-angle rotations.
static Rotation_ aboutX( const RealP& angleInRad ) 
{ return Rotation_( angleInRad, XAxis ); }
static Rotation_ aboutY( const RealP& angleInRad ) 
{ return Rotation_( angleInRad, YAxis ); }
static Rotation_ aboutZ( const RealP& angleInRad ) 
{ return Rotation_( angleInRad, ZAxis ); }
static Rotation_ aboutAxis( const RealP& angleInRad, const UnitVec3P& axis ) 
{ return Rotation_(angleInRad,axis); }
static Rotation_ aboutAxis( const RealP& angleInRad, const Vec3P& axis )     
{ return Rotation_(angleInRad,axis); }
void setToRotationAboutZ( const RealP& q ) { setRotationFromAngleAboutZ( q ); }

// Two-angle space-fixed rotations.
static Rotation_ aboutXThenOldY(const RealP& xInRad, const RealP& yInRad) 
{ return Rotation_( SpaceRotationSequence, xInRad, XAxis, yInRad, YAxis ); }
static Rotation_ aboutYThenOldX(const RealP& yInRad, const RealP& xInRad) 
{ return Rotation_( SpaceRotationSequence, yInRad, YAxis, xInRad, XAxis ); }
static Rotation_ aboutXThenOldZ(const RealP& xInRad, const RealP& zInRad) 
{ return Rotation_( SpaceRotationSequence, xInRad, XAxis, zInRad, ZAxis ); }
static Rotation_ aboutZThenOldX(const RealP& zInRad, const RealP& xInRad) 
{ return Rotation_( SpaceRotationSequence, zInRad, ZAxis, xInRad, XAxis ); }
static Rotation_ aboutYThenOldZ(const RealP& yInRad, const RealP& zInRad) 
{ return Rotation_( SpaceRotationSequence, yInRad, YAxis, zInRad, ZAxis ); }
static Rotation_ aboutZThenOldY(const RealP& zInRad, const RealP& yInRad) 
{ return Rotation_( SpaceRotationSequence, zInRad, ZAxis, yInRad, YAxis ); }

// Two-angle body fixed rotations (reversed space-fixed ones).
static Rotation_ aboutXThenNewY(const RealP& xInRad, const RealP& yInRad) 
{ return Rotation_( BodyRotationSequence, xInRad, XAxis, yInRad, YAxis ); }
static Rotation_ aboutYThenNewX(const RealP& yInRad, const RealP& xInRad)
{ return aboutXThenOldY(xInRad, yInRad); }
static Rotation_ aboutXThenNewZ(const RealP& xInRad, const RealP& zInRad)
{ return aboutZThenOldX(zInRad, xInRad); }
static Rotation_ aboutZThenNewX(const RealP& zInRad, const RealP& xInRad)
{ return aboutXThenOldZ(xInRad, zInRad); }
static Rotation_ aboutYThenNewZ(const RealP& yInRad, const RealP& zInRad)
{ return aboutZThenOldY(zInRad, yInRad); }
static Rotation_ aboutZThenNewY(const RealP& zInRad, const RealP& yInRad)
{ return aboutYThenOldZ(yInRad, zInRad); }

// Create a Rotation_ matrix by specifying only its z axis. 
// This will work for any stride UnitVec because there is always an implicit 
// conversion available to the packed form used as the argument.
explicit Rotation_( const UnitVec3P& uvecZ )  
{ setRotationFromOneAxis(uvecZ,ZAxis); }

// Create a Rotation_ matrix by specifying its x axis, and a "y like" axis. 
// We will take x seriously after normalizing, but use the y only to create 
// z = normalize(x X y), then y = z X x. Bad things happen if x and y are 
// aligned but we may not catch it.
Rotation_( const Vec3P& x, const Vec3P& yish )  
{ setRotationFromTwoAxes( UnitVec3P(x), XAxis, yish, YAxis ); }

// Set this Rotation_ to represent the same rotation as the passed-in quaternion.
void setToQuaternion( const QuaternionP& q )  { setRotationFromQuaternion(q); }

// Set this Rotation_ to represent a rotation of +q0 about body frame's Z axis, 
// followed by a rotation of +q1 about the body frame's NEW Y axis, 
// followed by a rotation of +q3 about the body frame's NEW X axis.
// See Kane, Spacecraft Dynamics, pg. 423, body-three: 3-2-1.
//  Similarly for BodyFixed123.
void setToBodyFixed321( const Vec3P& v)  
{ setRotationFromThreeAnglesThreeAxes(BodyRotationSequence, 
                                      v[0], ZAxis, v[1], YAxis, v[2], XAxis ); }
void setToBodyFixed123( const Vec3P& v)  
{ setRotationToBodyFixedXYZ(v); }

// Convert this Rotation_ matrix to an equivalent (angle,axis) representation: 
// Returned Vec4P is [angleInRadians, unitVectorX, unitVectorY, unitVectorZ].
Vec4P convertToAngleAxis() const  
{ return convertRotationToAngleAxis(); }

// Convert this Rotation_ matrix to equivalent quaternion representation.
QuaternionP convertToQuaternion() const  
{ return convertRotationToQuaternion(); }

// Set this Rotation_ to represent a rotation of +q0 about base frame's X axis, 
// followed by a rotation of +q1 about the base frame's (unchanged) Y axis.
void setToSpaceFixed12( const Vec2P& q ) 
{ setRotationFromTwoAnglesTwoAxes(SpaceRotationSequence,q[0],XAxis,q[1],YAxis);}

// Convert this Rotation_ matrix to the equivalent 1-2-3 body fixed Euler angle 
// sequence. Similarly, convert Rotation_ matrix to the equivalent 1-2 body  
// fixed Euler angle sequence. Similarly, convert Rotation_ matrix to the 
// equivalent 1-2 space fixed Euler angle sequence. 
Vec3P  convertToBodyFixed123() const  
{ return convertRotationToBodyFixedXYZ(); }
Vec2P  convertToBodyFixed12() const   
{ return convertRotationToBodyFixedXY(); }
Vec2P  convertToSpaceFixed12() const  
{ return convertTwoAxesRotationToTwoAngles(SpaceRotationSequence,XAxis,YAxis); }
};


//-----------------------------------------------------------------------------
/** (Advanced) This InverseRotation class is the inverse of a Rotation. See the
Rotation class for more information. **/
//-----------------------------------------------------------------------------
template <class P>
class InverseRotation_ : public Mat<3,3,P>::TransposeType {
typedef Mat<3,3,P>      Mat33P; // not the base type!
typedef Mat<2,2,P>      Mat22P;
typedef Mat<3,2,P>      Mat32P;
typedef Vec<2,P>        Vec2P;
typedef Vec<3,P>        Vec3P;
typedef Vec<4,P>        Vec4P;
typedef Quaternion_<P>  QuaternionP;
public:
/** This is the type of the underlying 3x3 matrix; note that it will have
unusual row and column spacing since we're viewing it as transposed. **/
typedef typename Mat<3,3,P>::TransposeType  BaseMat;

/** Note that the unit vectors representing the rows and columns of this
matrix do not necessarily have unit stride. **/
//@{
/** This is the type of a column of this InverseRotation. **/
typedef  UnitVec<P,BaseMat::RowSpacing>  ColType;
/** This is the type of a row of this InverseRotation. **/
typedef  UnitRow<P,BaseMat::ColSpacing>  RowType;
//@}

/** You should not ever construct one of these as they should only occur as 
expression intermediates resulting from use of the "~" operator on a Rotation.
But if you must, the default will produce an identity rotation. **/
InverseRotation_() : BaseMat(1) {}

/** An explicit implementation of the default copy constructor. **/
InverseRotation_( const InverseRotation_& R ) : BaseMat(R) {}
/** An explicit implementation of the default copy assignment operator. **/
InverseRotation_&  operator=( const InverseRotation_& R )  
{   BaseMat::operator=(R.asMat33());  return *this; }

/** Assuming this InverseRotation_ is R_AB, and given a symmetric dyadic matrix 
S_BB expressed in B, we can reexpress it in A using S_AA=R_AB*S_BB*R_BA. The 
matrix should be one that is formed as products of vectors expressed in A, such 
as inertia, unit inertia (gyration) or covariance matrices. This can be done 
efficiently exploiting properties of R and S. Cost is 57 flops.
@see Rotation::reexpressSymMat33() **/
SimTK_SimTKCOMMON_EXPORT SymMat<3,P> 
reexpressSymMat33(const SymMat<3,P>& S_BB) const;

/** We can invert an InverseRotation just by recasting it to a Rotation at 
zero cost. **/
//@{
const Rotation_<P>&  invert() const 
{return *reinterpret_cast<const Rotation_<P>*>(this);}
Rotation_<P>&        updInvert() {return *reinterpret_cast<Rotation_<P>*>(this);}
//@}

/** Transpose, and transpose operators (override BaseMat versions of transpose).
For an orthogonal matrix like this one transpose is the same as inverse. **/
//@{
const Rotation_<P>&  transpose() const  { return invert(); }
const Rotation_<P>&  operator~() const  { return invert(); }
Rotation_<P>&        updTranspose()     { return updInvert(); }
Rotation_<P>&        operator~()        { return updInvert(); }
//@}

/** Access individual rows and columns of this InverseRotation; no cost or
copying since suitably-cast references to the actual data are returned.
There are no writable versions of these methods since changing a single
row or column would violate the contract that these are always legitimate
rotation matrices. **/
//@{
const RowType&  row( int i ) const         
{ return reinterpret_cast<const RowType&>(asMat33()[i]); }
const RowType&  operator[]( int i ) const  { return row(i); }
const ColType&  col( int j ) const         
{ return reinterpret_cast<const ColType&>(asMat33()(j)); }
const ColType&  operator()( int j ) const  { return col(j); }
const ColType&  x() const                  { return col(0); }
const ColType&  y() const                  { return col(1); }
const ColType&  z() const                  { return col(2); }
//@}


/** Given a CoordinateAxis (XAxis,YAxis, or ZAxis) return a reference to
the corresponding column of this %InverseRotation matrix. The result is 
equivalent to multiplying R_AB*v_B where v_B is [1,0,0],[0,1,0], or [0,0,1], 
which would cost 15 flops, but requires no computation. **/
const ColType& getAxisUnitVec(CoordinateAxis axis) const 
{   return col(axis); }

/** Given a CoordinateDirection (+/-XAxis, etc.) return a unit vector in that
direction. The result is equivalent to multiplying R_AB*v_B where v_B is 
[+/-1,0,0], [0,+/-1,0], or [0,0,+/-1], which would cost 15 flops, but this 
method requires at most 3 flops. **/
const UnitVec<P,1> getAxisUnitVec(CoordinateDirection dir) const {
    const ColType& axDir = getAxisUnitVec(dir.getAxis());
    return dir.getDirection() > 0 ? UnitVec<P,1>( axDir) 
                                  : UnitVec<P,1>(-axDir); // cheap 
}

/** Conversion from InverseRotation_ to BaseMat. Note: asMat33 is slightly
more efficient than toMat33() (no copy), but you have to know the internal 
layout. **/
//@{
const BaseMat&  asMat33() const  { return *static_cast<const BaseMat*>(this); }
BaseMat         toMat33() const  { return asMat33(); }
//@}
};

/** Write a Rotation matrix to an output stream by writing out its underlying 
Mat33. **/
template <class P> SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream&, const Rotation_<P>&);
/** Write an InverseRotation matrix to an output stream by writing out its 
underlying Mat33. **/
template <class P> SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream&, const InverseRotation_<P>&);

/** Rotating a unit vector leaves it unit length, saving us from having to 
perform an expensive normalization. So we override the multiply operators here 
changing the return type to UnitVec or UnitRow. **/
//@{
template <class P, int S> inline UnitVec<P,1>  
operator*(const Rotation_<P>& R, const UnitVec<P,S>& v)        
{return UnitVec<P,1>(R.asMat33()* v.asVec3(),  true);}
template <class P, int S> inline UnitRow<P,1>  
operator*(const UnitRow<P,S>& r, const Rotation_<P>& R)        
{return UnitRow<P,1>(r.asRow3() * R.asMat33(), true);}
template <class P, int S> inline UnitVec<P,1>  
operator*(const InverseRotation_<P>& R, const UnitVec<P,S>& v) 
{return UnitVec<P,1>(R.asMat33()* v.asVec3(),  true);}
template <class P, int S> inline UnitRow<P,1>  
operator*(const UnitRow<P,S>& r, const InverseRotation_<P>& R) 
{return UnitRow<P,1>(r.asRow3() * R.asMat33(), true);}
//@}

// Couldn't implement these Rotation_ methods until InverseRotation_ was defined.
template <class P> inline
Rotation_<P>::Rotation_(const InverseRotation_<P>& R) 
:   Mat<3,3,P>( R.asMat33() ) {}

template <class P> inline Rotation_<P>&  
Rotation_<P>::operator=(const InverseRotation_<P>& R)  
{static_cast<Mat<3,3,P>&>(*this)  = R.asMat33();    return *this;}
template <class P> inline Rotation_<P>&  
Rotation_<P>::operator*=(const Rotation_<P>& R)        
{static_cast<Mat<3,3,P>&>(*this) *= R.asMat33();    return *this;}
template <class P> inline Rotation_<P>&  
Rotation_<P>::operator/=(const Rotation_<P>& R)        
{static_cast<Mat<3,3,P>&>(*this) *= (~R).asMat33(); return *this;}
template <class P> inline Rotation_<P>&  
Rotation_<P>::operator*=(const InverseRotation_<P>& R) 
{static_cast<Mat<3,3,P>&>(*this) *= R.asMat33();    return *this;}
template <class P> inline Rotation_<P>&  
Rotation_<P>::operator/=(const InverseRotation_<P>& R) 
{static_cast<Mat<3,3,P>&>(*this) *= (~R).asMat33(); return *this;}

/// Composition of Rotation matrices via operator*.
//@{
template <class P> inline Rotation_<P>
operator*(const Rotation_<P>&        R1, const Rotation_<P>&        R2)  
{return Rotation_<P>(R1) *= R2;}
template <class P> inline Rotation_<P>
operator*(const Rotation_<P>&        R1, const InverseRotation_<P>& R2)  
{return Rotation_<P>(R1) *= R2;}
template <class P> inline Rotation_<P>
operator*(const InverseRotation_<P>& R1, const Rotation_<P>&        R2)  
{return Rotation_<P>(R1) *= R2;}
template <class P> inline Rotation_<P>
operator*(const InverseRotation_<P>& R1, const InverseRotation_<P>& R2)  
{return Rotation_<P>(R1) *= R2;}
//@}

/// Composition of a Rotation matrix and the inverse of another Rotation via operator/, that is
/// R1/R2 == R1*(~R2).
//@{
template <class P> inline Rotation_<P>
operator/( const Rotation_<P>&        R1, const Rotation_<P>&        R2 )  
{return Rotation_<P>(R1) /= R2;}
template <class P> inline Rotation_<P>
operator/( const Rotation_<P>&        R1, const InverseRotation&     R2 )  
{return Rotation_<P>(R1) /= R2;}
template <class P> inline Rotation_<P>
operator/( const InverseRotation_<P>& R1, const Rotation_<P>&        R2 )  
{return Rotation_<P>(R1) /= R2;}
template <class P> inline Rotation_<P>
operator/( const InverseRotation_<P>& R1, const InverseRotation_<P>& R2 )  
{return Rotation_<P>(R1) /= R2;}
//@}


//------------------------------------------------------------------------------
}  // End of namespace SimTK

#endif // SimTK_SimTKCOMMON_ROTATION_H_


