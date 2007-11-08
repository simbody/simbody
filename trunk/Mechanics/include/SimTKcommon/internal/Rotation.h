//-----------------------------------------------------------------------------
// File:     Rotation.h
// Class:    Rotation and InverseRotation 
// Parent:   Mat33
// Purpose:  3x3 rotation class relating two right-handed orthogonal bases
//-----------------------------------------------------------------------------
#ifndef SIMTK_ROTATION_H 
#define SIMTK_ROTATION_H 

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Paul Mitiguy and Michael Sherman                                  *
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
#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/UnitVec.h"
#include "SimTKcommon/internal/Quaternion.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
namespace SimTK {


enum BodyOrSpaceType { BodyRotationSequence=0, SpaceRotationSequence=1 };

//-----------------------------------------------------------------------------
// Forward declarations
class InverseRotation;


//-----------------------------------------------------------------------------
/**
 * The Rotation class is a Mat33 that guarantees that the matrix is a legitimate 
 * 3x3 array associated with the relative orientation of two right-handed, 
 * orthogonal, unit vector bases.  The Rotation class takes advantage of
 * knowledge-specific information and properties of orthogonal matrices. 
 * For example, multiplication by a rotation matrix preserves a vector's length 
 * so unit vectors are still unit vectors afterwards and don't need to be re-normalized.
 * 
 * A rotation is an orthogonal matrix whose columns and rows are directions 
 * (that is, unit vectors) which are mutually orthogonal. 
 * Furthermore, if the columns (or rows) are labeled x,y,z it always holds
 * that z = x X y (rather than -(x X y)) ensuring that this is a right-handed
 * rotation matrix and not a reflection.
 *
 * Suppose there is a vector vF expressed in terms of the right-handed, orthogonal 
 * unit vectors Fx, Fy, Fz and one would like to express vF in terms of the
 * right-handed, orthogonal unit vectors Gx, Gy, Gz.  To calculate it, form
 *      vG = G_rotationMatrix_F * vF.  
 * Because a rotation is orthogonal, its transpose is its inverse. Hence
 * F_rotationMatrix_G = ~G_rotationMatrix_F  (where ~ is the SimTK "transpose").
 * This transpose matrix can be used to expressed vG in terms of Fx, Fy, Fz as
 *      vF = F_rotationMatrix_G * vG  or  vF = ~(G_rotationMatrix_F) * vG.
 *
 * Note: This entire class is exported for .dll purposes
 */
//------------------------------------------------------------------------------
class Rotation : public Mat33 {
public:
    // Default constructor and constructor-like methods
    Rotation() : Mat33(1) {}    
    Rotation&  setRotationToIdentityMatrix()  { Mat33::operator=(1);  return *this; }
    Rotation&  setRotationToNaN()             { Mat33::setToNaN();    return *this; } 

    // Default copy constructor and assignment operator
    Rotation( const Rotation& R ) : Mat33(R)  {}
    Rotation&  operator=( const Rotation& R )  { Mat33::operator=( R.asMat33() );  return *this; }

    // Constructor and related methods for right-handed rotation of an angle (in radians) about X, Y, or Z
    Rotation( Real angle, const CoordinateAxis& axis )             { setRotationFromAngleAboutAxis( angle, axis ); }
    Rotation( Real angle, const CoordinateAxis::XCoordinateAxis )  { setRotationFromAngleAboutX( std::cos(angle), std::sin(angle) ); }
    Rotation( Real angle, const CoordinateAxis::YCoordinateAxis )  { setRotationFromAngleAboutY( std::cos(angle), std::sin(angle) ); }
    Rotation( Real angle, const CoordinateAxis::ZCoordinateAxis )  { setRotationFromAngleAboutZ( std::cos(angle), std::sin(angle) ); }
    Rotation&  setRotationFromAngleAboutAxis( Real angle, const CoordinateAxis& axis )  { return axis.isXAxis() ? setRotationFromAngleAboutX(angle) : (axis.isYAxis() ? setRotationFromAngleAboutY(angle) : setRotationFromAngleAboutZ(angle) ); }
    Rotation&  setRotationFromAngleAboutX( Real angle )  { return setRotationFromAngleAboutX( std::cos(angle), std::sin(angle) ); }
    Rotation&  setRotationFromAngleAboutY( Real angle )  { return setRotationFromAngleAboutY( std::cos(angle), std::sin(angle) ); }
    Rotation&  setRotationFromAngleAboutZ( Real angle )  { return setRotationFromAngleAboutZ( std::cos(angle), std::sin(angle) ); }
    Rotation&  setRotationFromAngleAboutX( Real cosAngle, Real sinAngle )  { Mat33& R = *this;  R[0][0] = 1;   R[0][1] = R[0][2] = R[1][0] = R[2][0] = 0;   R[1][1] = R[2][2] = cosAngle;  R[1][2] = -(R[2][1] = sinAngle);  return *this; }
    Rotation&  setRotationFromAngleAboutY( Real cosAngle, Real sinAngle )  { Mat33& R = *this;  R[1][1] = 1;   R[0][1] = R[1][0] = R[1][2] = R[2][1] = 0;   R[0][0] = R[2][2] = cosAngle;  R[2][0] = -(R[0][2] = sinAngle);  return *this; }
    Rotation&  setRotationFromAngleAboutZ( Real cosAngle, Real sinAngle )  { Mat33& R = *this;  R[2][2] = 1;   R[0][2] = R[1][2] = R[2][0] = R[2][1] = 0;   R[0][0] = R[1][1] = cosAngle;  R[0][1] = -(R[1][0] = sinAngle);  return *this; }

    // Constructor and related methods for right-handed rotation of an angle (in radians) about an arbitrary vector
    Rotation( Real angle, const UnitVec3& unitVector ) { setRotationFromAngleAboutUnitVector(angle,unitVector); }
    Rotation( Real angle, const Vec3& nonUnitVector )  { setRotationFromAngleAboutNonUnitVector(angle,nonUnitVector); }
    Rotation&  setRotationFromAngleAboutNonUnitVector( Real angle, const Vec3& nonUnitVector )  { return setRotationFromAngleAboutUnitVector( angle, UnitVec3(nonUnitVector) ); }
    SimTK_SimTKCOMMON_EXPORT Rotation&  setRotationFromAngleAboutUnitVector( Real angle, const UnitVec3& unitVector );

    // Constructor and related methods for two-angle, two-axes, Body-fixed or Space-fixed rotation sequences (angles are in radians)
    // Constructor and related methods for three-angle Body-fixed or Space-fixed rotation sequences (angles are in radians)
    // Use1: Rotation someName( BodyRotationSequence,  2.1, XAxis,  3.2, YAxis );
    // Use2: Rotation someName( BodyRotationSequence,  2.1, XAxis,  3.2, YAxis,  4.3, ZAxis );
    // Use3: rotationObject.setRotationFromTwoAnglesTwoAxes(     SpaceRotationSequence,  2.1, ZAxis,  3.2, XAxis );
    // Use4: rotationObject.setRotationFromThreeAnglesThreeAxes( SpaceRotationSequence,  2.1, ZAxis,  3.2, XAxis,  4.3, ZAxis  );
    Rotation( BodyOrSpaceType bodyOrSpace, Real angle1, const CoordinateAxis& axis1, Real angle2, const CoordinateAxis& axis2 )                                            { setRotationFromTwoAnglesTwoAxes(    bodyOrSpace,angle1,axis1,angle2,axis2); }
    Rotation( BodyOrSpaceType bodyOrSpace, Real angle1, const CoordinateAxis& axis1, Real angle2, const CoordinateAxis& axis2, Real angle3, const CoordinateAxis& axis3 )  { setRotationFromThreeAnglesThreeAxes(bodyOrSpace,angle1,axis1,angle2,axis2,angle3,axis3); }
    SimTK_SimTKCOMMON_EXPORT Rotation&  setRotationFromTwoAnglesTwoAxes(     BodyOrSpaceType bodyOrSpace, Real angle1, const CoordinateAxis& axis1, Real angle2, const CoordinateAxis& axis2 ); 
    SimTK_SimTKCOMMON_EXPORT Rotation&  setRotationFromThreeAnglesThreeAxes( BodyOrSpaceType bodyOrSpace, Real angle1, const CoordinateAxis& axis1, Real angle2, const CoordinateAxis& axis2, Real angle3, const CoordinateAxis& axis3 );

    /// Set this Rotation to represent a rotation characterized by subsequent rotations of:
    /// +v[0] about the body frame's X axis,      followed by a rotation of 
    /// +v[1] about the body frame's NEW Y axis,  followed by a rotation of 
    /// +v[2] about the body frame's NEW Z axis.  See Kane, Spacecraft Dynamics, pg. 423, body-three: 1-2-3.
    void setRotationToBodyFixedXY( const Vec2& v)   { setRotationFromTwoAnglesTwoAxes(     BodyRotationSequence, v[0], XAxis, v[1], YAxis ); }
    void setRotationToBodyFixedXYZ( const Vec3& v)  { setRotationFromThreeAnglesThreeAxes( BodyRotationSequence, v[0], XAxis, v[1], YAxis, v[2], ZAxis ); }

    /// Constructor and related methods for relating a rotation matrix to a quaternion.
    explicit Rotation( const Quaternion& q )  { setRotationFromQuaternion(q); }
    SimTK_SimTKCOMMON_EXPORT Rotation&  setRotationFromQuaternion( const Quaternion& q );

    // Construct a Rotation directly from a Mat33 (we trust that m is a valid Rotation!)
    Rotation( const Mat33& m, bool ) : Mat33(m) {}

    // Constructs an (hopefully nearby) orthogonal rotation matrix from a generic Mat33.
    explicit Rotation( const Mat33& m )  { setRotationFromApproximateMat33(m); }
    SimTK_SimTKCOMMON_EXPORT Rotation&  setRotationFromApproximateMat33( const Mat33& m );

    // Calculate A.RotationMatrix.B by knowing one of B's unit vector expressed in A.
    // Note: The other vectors are perpendicular (but somewhat arbitrarily so).
    Rotation( const UnitVec3& uvec, const CoordinateAxis axis )  { setRotationFromOneAxis(uvec,axis); }
    SimTK_SimTKCOMMON_EXPORT Rotation&  setRotationFromOneAxis( const UnitVec3& uvec, const CoordinateAxis axis );

    // Calculate A.RotationMatrix.B by knowing one of B's unit vectors expressed in A and another vector that may be perpendicular 
    // If the 2nd vector is not perpendicular, no worries - we'll make it so it is perpendicular.
    Rotation( const UnitVec3& uveci, const CoordinateAxis& axisi, const Vec3& vecjApprox, const CoordinateAxis& axisjApprox )  { setRotationFromTwoAxes(uveci,axisi,vecjApprox,axisjApprox); }
    SimTK_SimTKCOMMON_EXPORT Rotation&  setRotationFromTwoAxes( const UnitVec3& uveci, const CoordinateAxis& axisi, const Vec3& vecjApprox, const CoordinateAxis& axisjApprox );

    // Converts rotation matrix to one or two or three orientation angles.
    // Note:  The result is most meaningful if the Rotation matrix is one that can be produced by such a sequence.
    // Use1:  someRotation.convertOneAxisRotationToOneAngle( XAxis );
    // Use2:  someRotation.convertTwoAxesRotationToTwoAngles(     SpaceRotationSequence, YAxis, ZAxis );
    // Use3:  someRotation.convertThreeAxesRotationToThreeAngles( SpaceRotationSequence, ZAxis, YAxis, XAxis );
    // Use4:  someRotation.convertRotationToAngleAxis();   Return: [angleInRadians, unitVectorX, unitVectorY, unitVectorZ].
    SimTK_SimTKCOMMON_EXPORT Real  convertOneAxisRotationToOneAngle( const CoordinateAxis& axis1 ) const;
    SimTK_SimTKCOMMON_EXPORT Vec2  convertTwoAxesRotationToTwoAngles(     BodyOrSpaceType bodyOrSpace, const CoordinateAxis& axis1, const CoordinateAxis& axis2 ) const;
    SimTK_SimTKCOMMON_EXPORT Vec3  convertThreeAxesRotationToThreeAngles( BodyOrSpaceType bodyOrSpace, const CoordinateAxis& axis1, const CoordinateAxis& axis2, const CoordinateAxis& axis3 ) const;
    SimTK_SimTKCOMMON_EXPORT Quaternion  convertRotationToQuaternion() const;
    Vec4  convertRotationToAngleAxis() const  { return convertRotationToQuaternion().convertQuaternionToAngleAxis(); }

    /// Special case of previous methods, e.g., convert Rotation matrix to equivalent body-fixed XYZ Euler angle sequence.
    Vec2  convertRotationToBodyFixedXY() const   { return convertTwoAxesRotationToTwoAngles( BodyRotationSequence, XAxis, YAxis ); }
    Vec3  convertRotationToBodyFixedXYZ() const  { return convertThreeAxesRotationToThreeAngles( BodyRotationSequence, XAxis, YAxis, ZAxis ); }

    /// Return true if "this" Rotation is nearly identical to "R" within a pointing angle error of less than or equal to 
    //  okPointingAngleErrorRads or machineEpsilon^{7/8), (e.g., 1E-14 rads in double precision).
    SimTK_SimTKCOMMON_EXPORT bool  isSameRotationToWithinAngle( const Rotation& R, Real okPointingAngleErrorRads ) const;
    bool  isSameRotationToWithinAngleOfMachinePrecision( const Rotation& R) const      { return isSameRotationToWithinAngle( R, SignificantReal ); }
    Real  getMaxAbsDifferenceInRotationElements( const Rotation& R ) const             { const Mat33& A = asMat33();  const Mat33& B = R.asMat33();  Real maxDiff = 0.0;  for( int i=0;  i<=2; i++ ) for( int j=0; j<=2; j++ ) { Real absDiff = std::fabs(A[i][j] - B[i][j]);  if( absDiff > maxDiff ) maxDiff = absDiff; }  return maxDiff; } 
    bool  areAllRotationElementsSameToEpsilon( const Rotation& R, Real epsilon ) const { return getMaxAbsDifferenceInRotationElements(R) <= epsilon ; }
    bool  areAllRotationElementsSameToMachinePrecision( const Rotation& R ) const      { return areAllRotationElementsSameToEpsilon( R, SignificantReal ); } 

    // Like copy constructor and copy assign but for inverse rotation.
    // Constructor allows implicit conversion from InverseRotation to Rotation.
    inline Rotation( const InverseRotation& );
    inline Rotation& operator=( const InverseRotation& );

    // Convert from Rotation to InverseRotation (no cost)
    const InverseRotation&  invert() const  { return *reinterpret_cast<const InverseRotation*>(this); }
    InverseRotation&        updInvert()     { return *reinterpret_cast<InverseRotation*>(this); }

    // Transpose, and transpose operators (override Mat33 versions of transpose).
    const InverseRotation&  transpose() const  { return invert(); }
    const InverseRotation&  operator~() const  { return invert(); }
    InverseRotation&        updTranspose()     { return updInvert(); }
    InverseRotation&        operator~()        { return updInvert(); }

    // Multiply-equals and divide-equals operators 
    inline Rotation&  operator*=( const Rotation& R );
    inline Rotation&  operator/=( const Rotation& R );
    inline Rotation&  operator*=( const InverseRotation& );
    inline Rotation&  operator/=( const InverseRotation& );

    /// Conversion from Rotation to Mat33.
    /// Note: asMat33 is more efficient than toMat33() (no copy), but you have to know the internal layout.
    const Mat33&  asMat33() const  { return *static_cast<const Mat33*>(this); }
    Mat33         toMat33() const  { return asMat33(); }

    // Note: This does not have unit stride.
    typedef  UnitVec<Mat33::RowSpacing>  ColType;
    typedef  UnitRow<Mat33::ColSpacing>  RowType;
    const RowType&  row( int i ) const         { return reinterpret_cast<const RowType&>(asMat33()[i]); }
    const ColType&  col( int j ) const         { return reinterpret_cast<const ColType&>(asMat33()(j)); }
    const ColType&  x() const                  { return col(0); }
    const ColType&  y() const                  { return col(1); }
    const ColType&  z() const                  { return col(2); }
    const RowType&  operator[]( int i ) const  { return row(i); }
    const ColType&  operator()( int j ) const  { return col(j); }

    /// Set the Rotation matrix directly - but you had better know what you are doing!
    Rotation&  setRotationFromMat33TrustMe( const Mat33& m )  { Mat33& R = *this;  R[0][0]=m[0][0];  R[0][1]=m[0][1];  R[0][2]=m[0][2];  R[1][0]=m[1][0];  R[1][1]=m[1][1];  R[1][2]=m[1][2];  R[2][0]=m[2][0];  R[2][1]=m[2][1];  R[2][2]=m[2][2];  return *this; }   
    Rotation&  setRotationColFromUnitVecTrustMe( int coli, const UnitVec3& uveci )  { Mat33& R = *this;   R[0][coli]=uveci[0];  R[1][coli]=uveci[1];  R[2][coli]=uveci[2];  return *this; }   
    Rotation&  setRotationFromUnitVecsTrustMe( const UnitVec3& colA, const UnitVec3& colB, const UnitVec3& colC )  { setRotationColFromUnitVecTrustMe(0,colA);  setRotationColFromUnitVecTrustMe(1,colB);  return setRotationColFromUnitVecTrustMe(2,colC); }   

//----------------------------------------------------------------------------------------------------
// The following code is obsolete - it is here temporarily for backward compatibility (Mitiguy 9/5/2007)
//----------------------------------------------------------------------------------------------------
private:
    // These static methods are like constructors with friendlier names.
    static Rotation zero() { return Rotation(); }
    static Rotation NaN()  { Rotation r;  r.setRotationToNaN();  return r; }

    /// By zero we mean "zero rotation", i.e., an identity matrix.
    Rotation&  setToZero()            { return setRotationToIdentityMatrix(); }
    Rotation&  setToIdentityMatrix()  { return setRotationToIdentityMatrix(); }
    Rotation&  setToNaN()             { return setRotationToNaN(); }
    static Rotation  trustMe( const Mat33& m )  { return Rotation(m,true); }

    // One-angle rotations.
    static Rotation aboutX( const Real& angleInRad ) { return Rotation( angleInRad, XAxis ); }
    static Rotation aboutY( const Real& angleInRad ) { return Rotation( angleInRad, YAxis ); }
    static Rotation aboutZ( const Real& angleInRad ) { return Rotation( angleInRad, ZAxis ); }
    static Rotation aboutAxis( const Real& angleInRad, const UnitVec3& axis ) { return Rotation(angleInRad,axis); }
    static Rotation aboutAxis( const Real& angleInRad, const Vec3& axis )     { return Rotation(angleInRad,axis); }
    void  setToRotationAboutZ( const Real& q ) { setRotationFromAngleAboutZ( q ); }

    // Two-angle space-fixed rotations.
    static Rotation aboutXThenOldY(const Real& xInRad, const Real& yInRad) { return Rotation( SpaceRotationSequence, xInRad, XAxis, yInRad, YAxis ); }
    static Rotation aboutYThenOldX(const Real& yInRad, const Real& xInRad) { return Rotation( SpaceRotationSequence, yInRad, YAxis, xInRad, XAxis ); }
    static Rotation aboutXThenOldZ(const Real& xInRad, const Real& zInRad) { return Rotation( SpaceRotationSequence, xInRad, XAxis, zInRad, ZAxis ); }
    static Rotation aboutZThenOldX(const Real& zInRad, const Real& xInRad) { return Rotation( SpaceRotationSequence, zInRad, ZAxis, xInRad, XAxis ); }
    static Rotation aboutYThenOldZ(const Real& yInRad, const Real& zInRad) { return Rotation( SpaceRotationSequence, yInRad, YAxis, zInRad, ZAxis ); }
    static Rotation aboutZThenOldY(const Real& zInRad, const Real& yInRad) { return Rotation( SpaceRotationSequence, zInRad, ZAxis, yInRad, YAxis ); }

    // Two-angle body fixed rotations (reversed space-fixed ones).
    static Rotation aboutXThenNewY(const Real& xInRad, const Real& yInRad) { return Rotation( BodyRotationSequence, xInRad, XAxis, yInRad, YAxis ); }
    static Rotation aboutYThenNewX(const Real& yInRad, const Real& xInRad) { return aboutXThenOldY(xInRad, yInRad); }
    static Rotation aboutXThenNewZ(const Real& xInRad, const Real& zInRad) { return aboutZThenOldX(zInRad, xInRad); }
    static Rotation aboutZThenNewX(const Real& zInRad, const Real& xInRad) { return aboutXThenOldZ(xInRad, zInRad); }
    static Rotation aboutYThenNewZ(const Real& yInRad, const Real& zInRad) { return aboutZThenOldY(zInRad, yInRad); }
    static Rotation aboutZThenNewY(const Real& zInRad, const Real& yInRad) { return aboutYThenOldZ(yInRad, zInRad); }

    /// Create a Rotation matrix by specifying only its z axis. 
    /// This will work for any stride UnitVec because there is always an implicit conversion available to the packed form used as the argument.
    explicit Rotation( const UnitVec3& uvecZ )  { setRotationFromOneAxis(uvecZ,ZAxis); }

    /// Create a Rotation matrix by specifying its x axis, and a "y like" axis. 
    //  We will take x seriously after normalizing, but use the y only to create z = normalize(x X y), 
    //  then y = z X x. Bad things happen if x and y are aligned but we may not catch it.
    Rotation( const Vec3& x, const Vec3& yish )  { setRotationFromTwoAxes( UnitVec3(x), XAxis, yish, YAxis ); }

    /// Set this Rotation to represent the same rotation as the passed-in quaternion.
    void setToQuaternion( const Quaternion& q )  { setRotationFromQuaternion(q); }

    /// Set this Rotation to represent a rotation of +q0 about the body frame's Z axis, 
    /// followed by a rotation of +q1 about the body frame's NEW Y axis, 
    /// followed by a rotation of +q3 about the body frame's NEW X axis.
    /// See Kane, Spacecraft Dynamics, pg. 423, body-three: 3-2-1.
    //  Similarly for BodyFixed123.
    void setToBodyFixed321( const Vec3& v)  { setRotationFromThreeAnglesThreeAxes( BodyRotationSequence, v[0], ZAxis, v[1], YAxis, v[2], XAxis ); }
    void setToBodyFixed123( const Vec3& v)  { setRotationToBodyFixedXYZ(v); }

    /// Convert this Rotation matrix to an equivalent (angle,axis) representation: 
    /// The returned Vec4 is [angleInRadians, unitVectorX, unitVectorY, unitVectorZ].
    Vec4  convertToAngleAxis() const  { return convertRotationToAngleAxis(); }

    /// Convert this Rotation matrix to equivalent quaternion representation.
    Quaternion  convertToQuaternion() const  { return convertRotationToQuaternion(); }

    /// Set this Rotation to represent a rotation of +q0 about the base frame's X axis, 
    /// followed by a rotation of +q1 about the base frame's (unchanged) Y axis.
    void setToSpaceFixed12( const Vec2& q ) { setRotationFromTwoAnglesTwoAxes( SpaceRotationSequence, q[0], XAxis, q[1], YAxis ); }

    /// Convert this Rotation matrix to the equivalent 1-2-3 body fixed Euler angle sequence.
    /// Similarly, convert Rotation matrix to the equivalent 1-2 body  fixed Euler angle sequence. 
    /// Similarly, convert Rotation matrix to the equivalent 1-2 space fixed Euler angle sequence. 
    Vec3  convertToBodyFixed123() const  { return convertRotationToBodyFixedXYZ(); }
    Vec2  convertToBodyFixed12() const   { return convertRotationToBodyFixedXY(); }
    Vec2  convertToSpaceFixed12() const  { return convertTwoAxesRotationToTwoAngles( SpaceRotationSequence, XAxis, YAxis ); }

//--------------------------------- PAUL CONTINUE FROM HERE ----------------------------------
public:
//--------------------------------------------------------------------------------------------
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

        const Mat33 E( 0,    s2oc1  ,  c2oc1  ,
                       0,      c2   ,   -s2   ,
                       1,  s1*s2oc1 , s1*c2oc1 );
        return E * w_PB_B;
    }


    /// Inverse of the above routine. Returned angular velocity is B in P,
    /// expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed321DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 Einv(  -s1  ,  0  ,  1 ,
                          c1*s2 ,  c2 ,  0 ,
                          c1*c2 , -s2 ,  0 );
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

        const Mat33 E(    c2oc1  , -s2oc1  , 0,
                            s2   ,    c2   , 0,
                       -s1*c2oc1 , s1*s2oc1, 1 );
        return E*w;
    }

    /// Inverse of the above routine. Returned angular velocity is B in P,
    /// expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed123DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 Einv( c1*c2 ,  s2 , 0 ,
                         -c1*s2 ,  c2 , 0 ,
                           s1   ,  0  , 1 );
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

        const Mat33 E( 0 ,   s2oc1  ,  c2oc1  ,
                       0 ,     c2   ,   -s2   ,
                       1 , s1*s2oc1 , s1*c2oc1 );
        const Vec3 qdot = E * w_PB_B;

        const Real t =  qdot[1]*qdot[2]*s1*ooc1;
        const Real a =  t*c2oc1; // d/dt s2oc1
        const Real b = -t*s2oc1; // d/dt c2oc1

        const Mat33 Edot( 0 ,       a           ,         b         ,
                          0 ,   -qdot[2]*s2     ,    -qdot[2]*c2    ,
                          0 , s1*a + qdot[1]*s2 , s1*b + qdot[1]*c2 );

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

        const Mat33 E(    c2oc1  , -s2oc1  , 0,
                            s2   ,    c2   , 0,
                       -s1*c2oc1 , s1*s2oc1, 1 );
        const Vec3 qdot = E * w_PB_B;

        const Real t =  qdot[1]*qdot[2]*s1*ooc1;
        const Real a =  t*c2oc1; // d/dt s2oc1
        const Real b = -t*s2oc1; // d/dt c2oc1

        const Mat33 Edot(       b           ,        -a         , 0,
                             qdot[2]*c2     ,    -qdot[2]*s2    , 0,
                         -s1*b - qdot[1]*c2 , s1*a + qdot[1]*s2 , 0 );

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

private:
    // This is only for the most trustworthy of callers, that is, methods of the Rotation class. 
    // There are a lot of ways for this NOT to be a legitimate rotation matrix -- be careful!!
    // Note that these are supplied in rows.
    Rotation( const Real& xx, const Real& xy, const Real& xz,
              const Real& yx, const Real& yy, const Real& yz,
              const Real& zx, const Real& zy, const Real& zz )
            : Mat33( xx,xy,xz, yx,yy,yz, zx,zy,zz ) {}

    // These next methods are highly-efficient power-user methods.  Read the documentation to understand them.
    SimTK_SimTKCOMMON_EXPORT Rotation&  setTwoAngleTwoAxesBodyFixedForwardCyclicalRotation(     Real cosAngle1, Real sinAngle1, const CoordinateAxis& axis1, Real cosAngle2, Real sinAngle2, const CoordinateAxis& axis2 );
    SimTK_SimTKCOMMON_EXPORT Rotation&  setThreeAngleTwoAxesBodyFixedForwardCyclicalRotation(   Real cosAngle1, Real sinAngle1, const CoordinateAxis& axis1, Real cosAngle2, Real sinAngle2, const CoordinateAxis& axis2, Real cosAngle3, Real sinAngle3 );
    SimTK_SimTKCOMMON_EXPORT Rotation&  setThreeAngleThreeAxesBodyFixedForwardCyclicalRotation( Real cosAngle1, Real sinAngle1, const CoordinateAxis& axis1, Real cosAngle2, Real sinAngle2, const CoordinateAxis& axis2, Real cosAngle3, Real sinAngle3, const CoordinateAxis& axis3 );

    // These next methods highly-efficient power-user methods convert Rotation matrices to orientation angles are highly-efficient power-user methods.  Read the documentation to understand them.
    SimTK_SimTKCOMMON_EXPORT Vec2  convertTwoAxesBodyFixedRotationToTwoAngles(     const CoordinateAxis& axis1, const CoordinateAxis& axis2 ) const;
    SimTK_SimTKCOMMON_EXPORT Vec3  convertTwoAxesBodyFixedRotationToThreeAngles(   const CoordinateAxis& axis1, const CoordinateAxis& axis2 ) const;
    SimTK_SimTKCOMMON_EXPORT Vec3  convertThreeAxesBodyFixedRotationToThreeAngles( const CoordinateAxis& axis1, const CoordinateAxis& axis2, const CoordinateAxis& axis3 ) const;

};


///-----------------------------------------------------------------------------
///  This InverseRotation class is the inverse of a Rotation 
///  See the Rotation class for information.
///-----------------------------------------------------------------------------
class InverseRotation : public Mat33::TransposeType {
public:
    // Convenient shortcut name
    typedef  Mat33::TransposeType  BaseMat;

    // Should not usually construct one of these as they should only occur as expression intermediates.
    // But if you must ...
    InverseRotation() : BaseMat(1) {}

    // Default constructor and copy constructor
    InverseRotation( const InverseRotation& R ) : BaseMat(R) {}
    InverseRotation&  operator=( const InverseRotation& R )  { BaseMat::operator=( R.asMat33() );  return *this; }

    // Convert from InverseRotation to Rotation (no cost)
    const Rotation&  invert() const { return *reinterpret_cast<const Rotation*>(this); }
    Rotation&  updInvert()          { return *reinterpret_cast<Rotation*>(this); }

    // Transpose, and transpose operators (override BaseMat versions of transpose).
    const Rotation&  transpose() const  { return invert(); }
    const Rotation&  operator~() const  { return invert(); }
    Rotation&        updTranspose()     { return updInvert(); }
    Rotation&        operator~()        { return updInvert(); }

    // Note that this does not have unit stride.
    typedef  UnitVec<BaseMat::RowSpacing>  ColType;
    typedef  UnitRow<BaseMat::ColSpacing>  RowType;
    const RowType&  row( int i ) const         { return reinterpret_cast<const RowType&>(asMat33()[i]); }
    const ColType&  col( int j ) const         { return reinterpret_cast<const ColType&>(asMat33()(j)); }
    const ColType&  x() const                  { return col(0); }
    const ColType&  y() const                  { return col(1); }
    const ColType&  z() const                  { return col(2); }
    const RowType&  operator[]( int i ) const  { return row(i); }
    const ColType&  operator()( int j ) const  { return col(j); }

    /// Conversion from InverseRotation to BaseMat.
    /// Note: asMat33 is more efficient than toMat33() (no copy), but you have to know the internal layout.
    const BaseMat&  asMat33() const  { return *static_cast<const BaseMat*>(this); }
    BaseMat         toMat33() const  { return asMat33(); }

};


SimTK_SimTKCOMMON_EXPORT std::ostream&  operator<<( std::ostream& o, const Rotation& m );

template <int S> inline UnitVec<1>  operator*( const Rotation& R,        const UnitVec<S>& v )       { return UnitVec<1>(R.asMat33()*v.asVec3(),  true); }
template <int S> inline UnitRow<1>  operator*( const UnitRow<S>& r,      const Rotation& R   )       { return UnitRow<1>(r.asRow3(), R.asMat33(), true); }
template <int S> inline UnitVec<1>  operator*( const InverseRotation& R, const UnitVec<S>& v )       { return UnitVec<1>(R.asMat33()*v.asVec3(),  true); }
template <int S> inline UnitRow<1>  operator*( const UnitRow<S>& r,      const InverseRotation& R )  { return UnitRow<1>(r.asRow3(), R.asMat33(), true); }

inline Rotation::Rotation( const InverseRotation& R) : Mat33( R.asMat33() )  {}

inline Rotation&  Rotation::operator=(  const InverseRotation& R )  { static_cast<Mat33&>(*this)  = R.asMat33();    return *this; }
inline Rotation&  Rotation::operator*=( const Rotation& R )         { static_cast<Mat33&>(*this) *= R.asMat33();    return *this; }
inline Rotation&  Rotation::operator/=( const Rotation& R )         { static_cast<Mat33&>(*this) *= (~R).asMat33(); return *this; }
inline Rotation&  Rotation::operator*=( const InverseRotation& R )  { static_cast<Mat33&>(*this) *= R.asMat33();    return *this; }
inline Rotation&  Rotation::operator/=( const InverseRotation& R )  { static_cast<Mat33&>(*this) *= (~R).asMat33(); return *this; }

inline Rotation  operator*( const Rotation&        R1, const Rotation&        R2 )  { return Rotation(R1) *= R2; }
inline Rotation  operator*( const Rotation&        R1, const InverseRotation& R2 )  { return Rotation(R1) *= R2; }
inline Rotation  operator*( const InverseRotation& R1, const Rotation&        R2 )  { return Rotation(R1) *= R2; }
inline Rotation  operator*( const InverseRotation& R1, const InverseRotation& R2 )  { return Rotation(R1) *= R2; }

inline Rotation operator/( const Rotation&        R1, const Rotation&        R2 )  {return Rotation(R1) /= R2;}
inline Rotation operator/( const Rotation&        R1, const InverseRotation& R2 )  {return Rotation(R1) /= R2;}
inline Rotation operator/( const InverseRotation& R1, const Rotation&        R2 )  {return Rotation(R1) /= R2;}
inline Rotation operator/( const InverseRotation& R1, const InverseRotation& R2 )  {return Rotation(R1) /= R2;}


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SIMTK_ROTATION_H_
//--------------------------------------------------------------------------


