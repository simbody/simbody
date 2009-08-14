//--------------------------------------------------------------------------
// File:     Rotation.cpp
// Class:    Rotation
// Parent:   Mat33
// Purpose:  3x3 rotation class relating two right-handed orthogonal bases
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

/**@file
 * Implementations of non-inline methods associated with Rotation class.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Rotation.h"

#include <cmath>
#include <utility>

//-------------------------------------------------------------------
namespace SimTK {


//----------------------------------------------------------------------
// Set Rotation for ANY two-angle ij rotation sequence (i,j = X,Y,Z)
//----------------------------------------------------------------------
Rotation&  Rotation::setRotationFromTwoAnglesTwoAxes( BodyOrSpaceType bodyOrSpace, Real angle1, const CoordinateAxis& axis1In, Real angle2, const CoordinateAxis& axis2In )  {

   // Non-const CoordinateAxes in case we have to switch axes
   CoordinateAxis axis1 = axis1In;
   CoordinateAxis axis2 = axis2In;

   // If axis2 is same as axis1, efficiently calculate with a one-angle, one-axis rotation
   if( axis1.isSameAxis(axis2) ) { return setRotationFromAngleAboutAxis( angle1+angle2, axis1 ); }

   // If using a SpaceRotationSequence, switch the order of the axes and the angles
   if( bodyOrSpace == SpaceRotationSequence )  { std::swap(angle1,angle2);  std::swap(axis1,axis2); }

   // If using a reverse cyclical, negate the signs of the angles
   if( axis1.isReverseCyclical(axis2) )  { angle1 = -angle1;   angle2 = -angle2; }

   // Calculate the sines and cosines (some hardware can do this more efficiently as one Taylor series)
   Real c1 = std::cos( angle1 ),  s1 = std::sin( angle1 );
   Real c2 = std::cos( angle2 ),  s2 = std::sin( angle2 );

   // All calculations are based on a body-fixed forward-cyclical rotation sequence
   return setTwoAngleTwoAxesBodyFixedForwardCyclicalRotation( c1,s1,axis1,  c2,s2,axis2 );
}


//----------------------------------------------------------------------
// Set Rotation for ANY three-angle ijk rotation sequence where (i,j,k = X,Y,Z)
//----------------------------------------------------------------------
Rotation&  Rotation::setRotationFromThreeAnglesThreeAxes( BodyOrSpaceType bodyOrSpace, Real angle1, const CoordinateAxis& axis1In, Real angle2, const CoordinateAxis& axis2, Real angle3, const CoordinateAxis& axis3In ) {

   // Non-const CoordinateAxes in case we have to switch axes
   CoordinateAxis axis1 = axis1In;
   CoordinateAxis axis3 = axis3In;

   // If axis2 is same as axis1 or axis3, efficiently calculate with a two-angle, two-axis rotation
   if( axis2.isSameAxis(axis1) ) return setRotationFromTwoAnglesTwoAxes( bodyOrSpace,  angle1+angle2, axis1,  angle3,axis3 );
   if( axis2.isSameAxis(axis3) ) return setRotationFromTwoAnglesTwoAxes( bodyOrSpace,  angle1,axis1,  angle2+angle3, axis3 );

   // If using a SpaceRotationSequence, switch the order of the axes and the angles
   if( bodyOrSpace == SpaceRotationSequence )  { std::swap(angle1,angle3);  std::swap(axis1,axis3); }

   // If using a reverse cyclical, negate the signs of the angles
   if( axis1.isReverseCyclical(axis2) )  { angle1 = -angle1;   angle2 = -angle2;   angle3 = -angle3; }

   // Calculate the sines and cosines (some hardware can do this more efficiently as one Taylor series)
   Real c1 = std::cos( angle1 ),  s1 = std::sin( angle1 );
   Real c2 = std::cos( angle2 ),  s2 = std::sin( angle2 );
   Real c3 = std::cos( angle3 ),  s3 = std::sin( angle3 );

   // All calculations are based on a body-fixed rotation sequence.
   // Determine whether this is a BodyXYX or BodyXYZ type of rotation sequence.
   if( axis1.isSameAxis(axis3) )  setThreeAngleTwoAxesBodyFixedForwardCyclicalRotation(   c1,s1,axis1, c2,s2,axis2, c3,s3 );
   else                           setThreeAngleThreeAxesBodyFixedForwardCyclicalRotation( c1,s1,axis1, c2,s2,axis2, c3,s3,axis3 );

   return *this;
}


//-------------------------------------------------------------------
// Calculate angle for ANY X or Y or Z rotation sequence
//-------------------------------------------------------------------
Real  Rotation::convertOneAxisRotationToOneAngle( const CoordinateAxis& axis1 ) const {

   // Get proper indices into Rotation matrix
   const CoordinateAxis axis2 = axis1.getNextAxis();
   const CoordinateAxis axis3 = axis2.getNextAxis();
   const int j = int( axis2 );
   const int k = int( axis3 );

   const Mat33& R =  asMat33();
   Real sinTheta = 0.5*( R[k][j] - R[j][k] );
   Real cosTheta = 0.5*( R[j][j] + R[k][k] );

   return std::atan2( sinTheta, cosTheta );
}


//-------------------------------------------------------------------
// Calculate angles for ANY two-angle ij rotation sequence (i,j = X,Y,Z)
//-------------------------------------------------------------------
Vec2  Rotation::convertTwoAxesRotationToTwoAngles( BodyOrSpaceType bodyOrSpace, const CoordinateAxis& axis1In, const CoordinateAxis& axis2In ) const {

   // Non-const CoordinateAxes in case we have to switch axes
   CoordinateAxis axis1 = axis1In;
   CoordinateAxis axis2 = axis2In;

   // If axis2 is same as axis1, efficiently calculate with a one-axis, one-angle method
   if( axis1.isSameAxis(axis2) ) { Real theta = 0.5*convertOneAxisRotationToOneAngle(axis1); return Vec2(theta,theta); }

   // If using a SpaceRotationSequence, switch the order of the axes (later switch the angles)
   if( bodyOrSpace == SpaceRotationSequence )  std::swap(axis1,axis2);

   // All calculations are based on a body-fixed rotation sequence
   Vec2 ans = convertTwoAxesBodyFixedRotationToTwoAngles( axis1, axis2 );

   // If using a SpaceRotationSequence, switch the angles now.
   if( bodyOrSpace == SpaceRotationSequence )  std::swap( ans[0], ans[1] );

   return ans;
}


//-------------------------------------------------------------------
// Calculate angles for ANY three-angle ijk rotation sequence where (i,j,k = X,Y,Z)
//-------------------------------------------------------------------
Vec3  Rotation::convertThreeAxesRotationToThreeAngles( BodyOrSpaceType bodyOrSpace, const CoordinateAxis& axis1In, const CoordinateAxis& axis2, const CoordinateAxis& axis3In ) const {

   // Non-const CoordinateAxes in case we have to switch axes
   CoordinateAxis axis1 = axis1In;
   CoordinateAxis axis3 = axis3In;

   // If all axes are same, efficiently calculate with a one-axis, one-angle method
   if( axis1.areAllSameAxes(axis2,axis3) ) { Real theta = 1.0/3.0*convertOneAxisRotationToOneAngle(axis1);  return Vec3(theta,theta,theta); }

   // If axis2 is same as axis1, efficiently calculate with a two-angle, two-axis rotation
   if( axis2.isSameAxis(axis1) ) {
      Vec2 xz = convertTwoAxesRotationToTwoAngles(bodyOrSpace,axis1,axis3);
      Real theta = 0.5 * xz[0];
      return Vec3( theta, theta, xz[1] );
   }

   // If axis2 is same as axis3, efficiently calculate with a two-angle, two-axis rotation
   if( axis2.isSameAxis(axis3) ) {
      Vec2 xz = convertTwoAxesRotationToTwoAngles(bodyOrSpace,axis1,axis3);
      Real theta = 0.5 * xz[1];
      return Vec3( xz[0], theta, theta);
   }

   // If using a SpaceRotationSequence, switch the order of the axes (later switch the angles)
   if( bodyOrSpace == SpaceRotationSequence )  std::swap(axis1,axis3);

   // All calculations are based on a body-fixed rotation sequence.
   // Determine whether this is a BodyXYX or BodyXYZ type of rotation sequence.
   Vec3 ans = axis1.isSameAxis(axis3) ? convertTwoAxesBodyFixedRotationToThreeAngles( axis1, axis2 ):
                                      convertThreeAxesBodyFixedRotationToThreeAngles( axis1, axis2, axis3 );

   // If using a SpaceRotationSequence, switch the angles now.
   if( bodyOrSpace == SpaceRotationSequence )  std::swap( ans[0], ans[2] );

   return ans;
}



//--------------------------------------------------------------------------
// Set Rotation ONLY for two-angle, two-axes, body-fixed, ij rotation sequence where i != j.
//--------------------------------------------------------------------------
Rotation&  Rotation::setTwoAngleTwoAxesBodyFixedForwardCyclicalRotation( Real cosAngle1, Real sinAngle1, const CoordinateAxis& axis1, Real cosAngle2, Real sinAngle2, const CoordinateAxis& axis2 ) {
   // Ensure this method has proper arguments
   assert( axis1.isDifferentAxis( axis2 ) );

   CoordinateAxis axis3 = axis1.getThirdAxis( axis2 );
   const int i = int( axis1 );
   const int j = int( axis2 );
   const int k = int( axis3 );

   Mat33& R = *this;
   R[i][i] =  cosAngle2;
   R[i][j] =  0;
   R[i][k] =  sinAngle2;
   R[j][i] =  sinAngle2 * sinAngle1;
   R[j][j] =  cosAngle1;
   R[j][k] = -sinAngle1 * cosAngle2;
   R[k][i] = -sinAngle2 * cosAngle1;
   R[k][j] =  sinAngle1;
   R[k][k] =  cosAngle1 * cosAngle2;
   return *this;
}


//--------------------------------------------------------------------------
// Set Rotation ONLY for three-angle, two-axes, body-fixed, iji rotation sequence where i != j
//--------------------------------------------------------------------------
Rotation&  Rotation::setThreeAngleTwoAxesBodyFixedForwardCyclicalRotation( Real cosAngle1, Real sinAngle1, const CoordinateAxis& axis1, Real cosAngle2, Real sinAngle2, const CoordinateAxis& axis2, Real cosAngle3, Real sinAngle3 )  {

   // Ensure this method has proper arguments
   assert( axis1.isDifferentAxis( axis2 ) );

   // Repeated calculations (for efficiency)
   Real s1c3 = sinAngle1 * cosAngle3;
   Real s3c1 = sinAngle3 * cosAngle1;
   Real s1s3 = sinAngle1 * sinAngle3;
   Real c1c3 = cosAngle1 * cosAngle3;

   CoordinateAxis axis3 = axis1.getThirdAxis( axis2 );
   const int i = int( axis1 );
   const int j = int( axis2 );
   const int k = int( axis3 );

   Mat33& R =  *this;
   R[i][i] =  cosAngle2;
   R[i][j] =  sinAngle2 * sinAngle3;
   R[i][k] =  sinAngle2 * cosAngle3;
   R[j][i] =  sinAngle1 * sinAngle2;
   R[j][j] =  c1c3 - cosAngle2 * s1s3;
   R[j][k] = -s3c1 - cosAngle2 * s1c3;
   R[k][i] = -sinAngle2 * cosAngle1;
   R[k][j] =  s1c3 + cosAngle2 * s3c1;
   R[k][k] = -s1s3 + cosAngle2 * c1c3;
   return *this;
}


//--------------------------------------------------------------------------
// Set Rotation ONLY for three-angle, three-axes, body-fixed, ijk rotation sequence where i != j != k.
//--------------------------------------------------------------------------
Rotation&  Rotation::setThreeAngleThreeAxesBodyFixedForwardCyclicalRotation( Real cosAngle1, Real sinAngle1, const CoordinateAxis& axis1, Real cosAngle2, Real sinAngle2, const CoordinateAxis& axis2, Real cosAngle3, Real sinAngle3, const CoordinateAxis& axis3 ) {

   // Ensure this method has proper arguments
   assert( axis1.areAllDifferentAxes(axis2,axis3) );

   // Repeated calculations (for efficiency)
   Real s1c3 = sinAngle1 * cosAngle3;
   Real s3c1 = sinAngle3 * cosAngle1;
   Real s1s3 = sinAngle1 * sinAngle3;
   Real c1c3 = cosAngle1 * cosAngle3;

   const int i = int(axis1);
   const int j = int(axis2);
   const int k = int(axis3);

   Mat33& R = *this;
   R[i][i] =  cosAngle2 * cosAngle3;
   R[i][j] = -sinAngle3 * cosAngle2;
   R[i][k] =  sinAngle2;
   R[j][i] =  s3c1 + sinAngle2 * s1c3;
   R[j][j] =  c1c3 - sinAngle2 * s1s3;
   R[j][k] = -sinAngle1 * cosAngle2;
   R[k][i] =  s1s3 - sinAngle2 * c1c3;
   R[k][j] =  s1c3 + sinAngle2 * s3c1;
   R[k][k] =  cosAngle1 * cosAngle2;
   return *this;
}


//-------------------------------------------------------------------
// Calculate angles ONLY for a two-angle, two-axes, body-fixed, ij rotation sequence where i != j.
//-------------------------------------------------------------------
Vec2  Rotation::convertTwoAxesBodyFixedRotationToTwoAngles( const CoordinateAxis& axis1, const CoordinateAxis& axis2 ) const {

   // Ensure this method has proper arguments
   assert( axis1.isDifferentAxis( axis2 ) );

   CoordinateAxis axis3 = axis1.getThirdAxis(axis2);
   const int i = int(axis1);
   const int j = int(axis2);
   const int k = int(axis3);
   const Mat33& R =  asMat33();

   // Can use either direct method (fast) or all matrix elements with the overhead of two additional square roots (possibly more accurate)
   const Real sinTheta1Direct = R[k][j];
   const Real signSinTheta1 = sinTheta1Direct > 0 ? 1.0 : -1.0;
   const Real sinTheta1Alternate = signSinTheta1 * std::sqrt( square(R[j][i]) + square(R[j][k]) );
   const Real sinTheta1 = 0.5*( sinTheta1Direct + sinTheta1Alternate );

   const Real cosTheta1Direct = R[j][j];
   const Real signCosTheta1 = cosTheta1Direct > 0 ? 1.0 : -1.0;
   const Real cosTheta1Alternate = signCosTheta1 * std::sqrt( square(R[k][i]) + square(R[k][k]) );
   const Real cosTheta1 = 0.5*( cosTheta1Direct + cosTheta1Alternate );

   Real theta1 = std::atan2( sinTheta1, cosTheta1 );

   // Repeat for theta2
   const Real sinTheta2Direct = R[i][k];
   const Real signSinTheta2 = sinTheta2Direct > 0 ? 1.0 : -1.0;
   const Real sinTheta2Alternate = signSinTheta2 * std::sqrt( square(R[j][i]) + square(R[k][i]) );
   const Real sinTheta2 = 0.5*( sinTheta2Direct + sinTheta2Alternate );

   const Real cosTheta2Direct = R[i][i];
   const Real signCosTheta2 = cosTheta2Direct > 0 ? 1.0 : -1.0;
   const Real cosTheta2Alternate = signCosTheta2 * std::sqrt( square(R[j][k]) + square(R[k][k]) );
   const Real cosTheta2 = 0.5*( cosTheta2Direct + cosTheta2Alternate );

   Real theta2 = std::atan2( sinTheta2, cosTheta2 );

   // If using a reverse cyclical, negate the signs of the angles
   if( axis1.isReverseCyclical(axis2) )  { theta1 = -theta1;  theta2 = -theta2; }

   // Return values have the following ranges:
   // -pi   <=  theta1  <=  +pi
   // -pi   <=  theta2  <=  +pi
   return Vec2( theta1, theta2 );
}


//-------------------------------------------------------------------
// Calculate angles ONLY for a three-angle, two-axes, body-fixed, iji rotation sequence where i != j.
//-------------------------------------------------------------------
Vec3  Rotation::convertTwoAxesBodyFixedRotationToThreeAngles( const CoordinateAxis& axis1, const CoordinateAxis& axis2 ) const {

   // Ensure this method has proper arguments
   assert( axis1.isDifferentAxis( axis2 ) );

   CoordinateAxis axis3 = axis1.getThirdAxis( axis2 );
   const int i = int( axis1 );
   const int j = int( axis2 );
   const int k = int( axis3 );
 
   // Need to know if using a forward or reverse cyclical
   Real plusMinus = 1.0,  minusPlus = -1.0;
   if( axis1.isReverseCyclical(axis2) ) { plusMinus = -1.0,  minusPlus = 1.0; }

   // Shortcut to the elements of the rotation matrix
   const Mat33& R =  asMat33();

   // Calculate theta2 using lots of information in the rotation matrix
   const Real Rsum = std::sqrt( 0.5*( square(R[i][j]) + square(R[i][k]) + square(R[j][i]) + square(R[k][i]) ) );  
   const Real theta2 =  atan2( Rsum, R[i][i] );  // Rsum = abs(sin(theta2)) is inherently positive
   Real theta1, theta3;

   // There is a "singularity" when sin(theta2) == 0
   if( Rsum > 4*Eps ) {
      theta1  =  atan2( R[j][i], minusPlus*R[k][i] );
      theta3  =  atan2( R[i][j], plusMinus*R[i][k] );
   }
   else if( R[i][i] > 0 ) {
      const Real spos = plusMinus*R[k][j] + minusPlus*R[j][k];  // 2*sin(theta1 + theta3)
      const Real cpos = R[j][j] + R[k][k];                      // 2*cos(theta1 + theta3)
      const Real theta1PlusTheta3 = atan2( spos, cpos );
      theta1 = theta1PlusTheta3;  // Arbitrary split
      theta3 = 0;                 // Arbitrary split
   }
   else {
      const Real sneg = plusMinus*R[k][j] + plusMinus*R[j][k];  // 2*sin(theta1 - theta3)
      const Real cneg = R[j][j] - R[k][k];                      // 2*cos(theta1 - theta3)
      const Real theta1MinusTheta3 = atan2( sneg, cneg );
      theta1 = theta1MinusTheta3;  // Arbitrary split
      theta3 = 0;                  // Arbitrary split
   }

   // Return values have the following ranges:
   // -pi   <=  theta1  <=  +pi
   //   0   <=  theta2  <=  +pi    (Rsum is inherently positive)
   // -pi   <=  theta3  <=  +pi
   return Vec3( theta1, theta2, theta3 );
}


//-------------------------------------------------------------------
// Calculate angles ONLY for a three-angle, three-axes, body-fixed, ijk rotation sequence where i != j and j != k.
//-------------------------------------------------------------------
Vec3  Rotation::convertThreeAxesBodyFixedRotationToThreeAngles( const CoordinateAxis& axis1, const CoordinateAxis& axis2, const CoordinateAxis& axis3 ) const {

   // Ensure this method has proper arguments
   assert( axis1.areAllDifferentAxes(axis2,axis3) );

   const int i = int(axis1);
   const int j = int(axis2);
   const int k = int(axis3);
 
   // Need to know if using a forward or reverse cyclical
   Real plusMinus = 1.0,  minusPlus = -1.0;
   if( axis1.isReverseCyclical(axis2) ) { plusMinus = -1.0,  minusPlus = 1.0; }

   // Shortcut to the elements of the rotation matrix
   const Mat33& R =  asMat33();

   // Calculate theta2 using lots of information in the rotation matrix
   Real Rsum   =  std::sqrt( 0.5*( square(R[i][i]) + square(R[i][j]) + square(R[j][k]) + square(R[k][k])) );
   Real theta2 =  atan2( plusMinus*R[i][k], Rsum ); // Rsum = abs(cos(theta2)) is inherently positive
   Real theta1, theta3;

   // There is a "singularity" when cos(theta2) == 0
   if( Rsum > 4*Eps ) {
      theta1 =  atan2( minusPlus*R[j][k], R[k][k] );
      theta3 =  atan2( minusPlus*R[i][j], R[i][i] );
   }
   else if( plusMinus*R[i][k] > 0 ) {
      const Real spos = R[j][i] + plusMinus*R[k][j];  // 2*sin(theta1 + plusMinus*theta3)
      const Real cpos = R[j][j] + minusPlus*R[k][i];  // 2*cos(theta1 + plusMinus*theta3)
      const Real theta1PlusMinusTheta3 = atan2( spos, cpos );
      theta1 = theta1PlusMinusTheta3;  // Arbitrary split
      theta3 = 0;                      // Arbitrary split
   }
   else {
      const Real sneg = plusMinus*(R[k][j] + minusPlus*R[j][i]);  // 2*sin(theta1 + minusPlus*theta3)
      const Real cneg = R[j][j] + plusMinus*R[k][i];              // 2*cos(theta1 + minusPlus*theta3)
      const Real theta1MinusPlusTheta3 = atan2( sneg, cneg );
      theta1 = theta1MinusPlusTheta3;  // Arbitrary split
      theta3 = 0;                      // Arbitrary split
   }

   // Return values have the following ranges:
   // -pi   <=  theta1  <=  +pi
   // -pi/2 <=  theta2  <=  +pi/2   (Rsum is inherently positive)
   // -pi   <=  theta3  <=  +pi
   return Vec3( theta1, theta2, theta3 );
}


//-------------------------------------------------------------------
// Calculate A.RotationMatrix.B by knowing one of B's unit vector expressed in A.
// Note: The other vectors are perpendicular (but somewhat arbitrarily so).
//-------------------------------------------------------------------
Rotation&  Rotation::setRotationFromOneAxis( const UnitVec3& uveci, const CoordinateAxis axisi )  {
   // Find a unit vector that is perpendicular to uveci
   const UnitVec3 uvecj = uveci.perp();

   // Find another unit vector that is perpendicular to both uveci and uvecj
   const UnitVec3 uveck = UnitVec3( cross( uveci, uvecj ) );

   // Determine which of the axes this corresponds to
   CoordinateAxis axisj = axisi.getNextAxis();
   CoordinateAxis axisk = axisj.getNextAxis();

   // Fill in the correct elements of the Rotation matrix
   setRotationColFromUnitVecTrustMe( int(axisi), uveci );
   setRotationColFromUnitVecTrustMe( int(axisj), uvecj );
   setRotationColFromUnitVecTrustMe( int(axisk), uveck );

   return *this;
}


//-------------------------------------------------------------------
// Calculate A.RotationMatrix.B by knowing one of B's unit vectors expressed in A and another vector that may be perpendicular
// If the 2nd vector is not perpendicular, no worries - we'll make it so it is perpendicular.
//-------------------------------------------------------------------
Rotation&  Rotation::setRotationFromTwoAxes( const UnitVec3& uveci, const CoordinateAxis& axisi, const Vec3& vecjApprox, const CoordinateAxis& axisjApprox ) {

   // Ensure vecjApprox is not a zero vector and the axes are not the same.
   Real magnitudeOfVecjApprox = vecjApprox.normSqr();
   if( magnitudeOfVecjApprox==0  ||  axisi.isSameAxis(axisjApprox) )
      return setRotationFromOneAxis( uveci, axisi );

   // Create a unit vector that is perpendicular to both uveci and vecjApprox.
   Vec3 veck = cross( uveci, vecjApprox );

   // Make sure (a x b) is not too close to the zero vector (vectors are nearly parallel)
   // Since |a x b| = |a|*|b|*sin(theta), it is easy to determine when sin(theta) = |a x b| / (|a| * |b|) is small.
   // In other words, sin(theta) = |a x b| / (|a| * |b|)  where |a| = 1 (since uveci is a unit vector).
   // I use SqrtEps (which is approx 1.0E-8) which means that the angle between the vectors is less than 1.E-08 rads = 6.0E-7 deg.
   Real magnitudeOfVeck = veck.normSqr();   // Magnitude of the cross product
   if( magnitudeOfVeck < SimTK::SqrtEps * magnitudeOfVecjApprox )
      return setRotationFromOneAxis( uveci, axisi );

   // Find a unit vector that is perpendicular to both uveci and vecjApprox
   UnitVec3 uveck = UnitVec3( veck );

   // Compute the unit vector perpendicular to both uveci and uveck
   UnitVec3 uvecj = UnitVec3( cross( uveck, uveci ) );

   // Determine which of the axes this corresponds to
   CoordinateAxis axisj = axisi.getNextAxis();
   CoordinateAxis axisk = axisj.getNextAxis();

   // If axisj is axisjApprox, all is good, otherwise switch axisj and axisk and negate the k'th axis.
   if( axisj.isDifferentAxis(axisjApprox) )
   {
      std::swap( axisj, axisk );
      uveck = -uveck;
   }

   // Fill in the correct elements of the Rotation matrix
   setRotationColFromUnitVecTrustMe( int(axisi), uveci );
   setRotationColFromUnitVecTrustMe( int(axisj), uvecj );
   setRotationColFromUnitVecTrustMe( int(axisk), uveck );

   return *this;
}


//-------------------------------------------------------------------
// Set this rotation matrix from the associated quaternion.
//-------------------------------------------------------------------
Rotation&  Rotation::setRotationFromQuaternion( const Quaternion& q )  {
    const Real q00=q[0]*q[0], q11=q[1]*q[1], q22=q[2]*q[2], q33=q[3]*q[3];
    const Real q01=q[0]*q[1], q02=q[0]*q[2], q03=q[0]*q[3];
    const Real q12=q[1]*q[2], q13=q[1]*q[3], q23=q[2]*q[3];

    Mat33::operator=( Mat33( q00+q11-q22-q33,   2.*(q12-q03)  ,   2.*(q13+q02),
                              2.*(q12+q03)  ,  q00-q11+q22-q33,   2.*(q23-q01),
                              2.*(q13-q02)  ,   2.*(q23+q01)  , q00-q11-q22+q33 )  );
    return *this;
}


//-------------------------------------------------------------------
// Convert this Rotation matrix to the equivalent quaternion.
// We use a modification of Richard Spurrier's method.
// [Spurrier, R.A., "Comment on 'Singularity-Free Extraction of a Quaternion from a
// Direction-Cosine Matrix'", J. Spacecraft and Rockets, 15(4):255, 1977.]
// Our modification avoids all but one square root and divide.
// In each of the four cases we compute 4q[m]*q where m is the "max" element, with
// m=0 if the trace is larger than any diagonal or
// m=i if the i,i element is the largest diagonal and larger than the trace.
// Then when we normalize at the end the scalar 4q[m] evaporates leaving us
// with a perfectly normalized quaternion.
//
// The returned quaternion can be interpreted as a rotation angle a about a unit
// vector v=[vx vy vz] like this:    q = [ cos(a/2) sin(a/2)*v ]
// We canonicalize the returned quaternion by insisting that
// cos(a/2) >= 0, meaning that -180 < a <= 180.
//-------------------------------------------------------------------
Quaternion  Rotation::convertRotationToQuaternion() const {
    const Mat33& R = asMat33();

    // Stores the return values [cos(theta/2), lambda1*sin(theta/2), lambda2*sin(theta/2), lambda3*sin(theta/2)]
    Vec4 q;

    // Check if the trace is larger than any diagonal
    const Real tr = R.trace();
    if( tr >= R(0,0)  &&  tr >= R(1,1)  &&  tr >= R(2,2) ) {
        q[0] = 1 + tr;
        q[1] = R(2,1) - R(1,2);
        q[2] = R(0,2) - R(2,0);
        q[3] = R(1,0) - R(0,1);

    // Check if R(0,0) is largest along the diagonal
    } else if( R(0,0) >= R(1,1)  &&  R(0,0) >= R(2,2)  ) {
        q[0] = R(2,1) - R(1,2);
        q[1] = 1 - (tr - 2*R(0,0));
        q[2] = R(0,1)+R(1,0);
        q[3] = R(0,2)+R(2,0);

    // Check if R(1,1) is largest along the diagonal
    } else if( R(1,1) >= R(2,2) ) {
        q[0] = R(0,2) - R(2,0);
        q[1] = R(0,1) + R(1,0);
        q[2] = 1 - (tr - 2*R(1,1));
        q[3] = R(1,2) + R(2,1);

    // R(2,2) is largest along the diagonal
    } else {
        q[0] = R(1,0) - R(0,1);
        q[1] = R(0,2) + R(2,0);
        q[2] = R(1,2) + R(2,1);
        q[3] = 1 - (tr - 2*R(2,2));
    }
    Real scale = q.norm();
    if( q[0] < 0 )  scale = -scale;   // canonicalize
    return Quaternion(q/scale, true); // prevent re-normalization
}



//-------------------------------------------------------------------
// Determine whether "this" Rotation matrix is nearly identical to the one passed in (R)
// within a specified "pointing angle error".  If "this" and "R" are nearly identical,
// transpose(this) * R = closeToIdentity matrix.  After finding the equivalent angle-axis
// for closeToIdentityMatrix, check the angle to see if it is less than okPointingAngleErrorRads.
//-------------------------------------------------------------------
bool  Rotation::isSameRotationToWithinAngle( const Rotation& R, Real okPointingAngleErrorRads ) const {
    // The absolute minimum pointing error is 0 if "this" and R are identical.
    // The absolute maximum pointing error should be Pi (to machine precision).
    assert( 0 <= okPointingAngleErrorRads && okPointingAngleErrorRads <= Pi);
    const Rotation closeToIdentityMatrix = ~(*this) * R;
    const Vec4 angleAxisEquivalent = closeToIdentityMatrix.convertRotationToAngleAxis();
    const Real pointingError = std::abs( angleAxisEquivalent[0] );
    return pointingError <= okPointingAngleErrorRads;
}


//-------------------------------------------------------------------
// This method is a poor man's orthogonalization from the supplied matrix to make a legitimate rotation.
// This is done by conversion to quaternions and back to Rotation and may not produce the closest rotation.
//-------------------------------------------------------------------
Rotation&  Rotation::setRotationFromApproximateMat33( const Mat33& m ) {
    // First, create a Rotation from the given Mat33
    const Rotation approximateRotation(m, true);
    const Quaternion q = approximateRotation.convertRotationToQuaternion();
    this->setRotationFromQuaternion( q );
    return *this;
}


//-------------------------------------------------------------------
// Note: It may be slightly more efficient to set this matrix directly rather than via the quaternion.
//-------------------------------------------------------------------
Rotation&  Rotation::setRotationFromAngleAboutUnitVector( Real angleInRad, const UnitVec3& unitVector ) {
    Quaternion q;
    q.setQuaternionFromAngleAxis( angleInRad, unitVector );
    return setRotationFromQuaternion( q );
}

//-------------------------------------------------------------------
// Use sneaky tricks from Featherstone to rotate a symmetric dyadic
// matrix. Consider the current Rotation matrix to be R_AB. We want
// to return S_AA=R_AB*S_BB*R_BA. Would be 90 flops for 3x3s, 75
// since we only need six elements in final result. Here we'll get
// it done in 57 flops.
// Consider S=[ a d e ]
//            [ d b f ]
//            [ e f c ]
//
// First, factor S into S=L+D+vx with v=~[-f e 0] (x means cross 
// product matrix):
//        [a-c   d   0]     [c 0 0]       [ 0  0  e]
//    L = [ d   b-c  0] D = [0 c 0]  vx = [ 0  0  f]
//        [2e   2f   0]     [0 0 c]       [-e -f  0]
// (4 flops to calculate L)
//
// A cross product matrix identity says R*vx*R'=(R*v)x, so:
//    S'=R*S*~R = R*L*~R + D + (R*v)x. 
// Let Y'=R*L, Z=Y'*~R. We only need the lower triangle of Z and a 
// 2x2 square of Y'.
//
// Don't-care's below are marked "-". Reminder: square bracket [i]
// index of a matrix means "row i", round bracket (j) means "col j".
//
//        [  -   -  0 ]
//   Y' = [ Y00 Y01 0 ]   Y = [ R[1]*L(0)  R[1]*L(1) ]  20 flops
//        [ Y10 Y11 0 ]       [ R[2]*L(0)  R[2]*L(1) ]
//
//   Z = [   Z00           -           -      ]
//       [ Y[0]*~R[0]  Y[0]*~R[1]      -      ]   15 flops (use only 2
//       [ Y[1]*~R[0]  Y[1]*~R[1]  Y[1]*~R[2] ]   elements of R's rows)
//
//   Z00 = (L00+L11)-(Z11+Z22)  3 flops ( because rotation preserves trace)
//
//        [R01*c-R00*e]            [  0       -    -  ]
//   R*v =[R11*c-R10*e]   (R*v)x = [ Rv[2]    0    -  ]
//        [R21*c-R20*e]            [-Rv[1]  Rv[0]  0  ]
// (R*v is 9 flops)
//
//        [  Z00 + c          -           -    ]
//   S' = [ Z10 + Rv[2]   Z11 + c         -    ]
//        [ Z20 - Rv[1]  Z21 + Rv[0]   Z22 + c ]
//
// which takes 6 more flops. Total 6+9Rv+18Z+20Y+4L=57.
//
// (I actually looked at the generated code in VC++ 2005 and Intel C++
//  version 11.1 and counted exactly 57 inline flops.)
//
// NOTE: there are two implementations of this routine that have
// to be kept in sync -- this one and the identical one for 
// InverseRotation right below.
//-------------------------------------------------------------------
SymMat33 Rotation::reexpressSymMat33(const SymMat33& S_BB) const {
    const Real a=S_BB(0,0), b=S_BB(1,1), c=S_BB(2,2);
    const Real d=S_BB(1,0), e=S_BB(2,0), f=S_BB(2,1);
    const Mat33& R   = this->asMat33();
    const Mat32& RR  = R.getSubMat<3,2>(0,0); // first two columns of R

    const Mat32 L( a-c ,  d,
                    d  , b-c,
                   2*e , 2*f );

    const Mat22 Y( R[1]*L(0), R[1]*L(1),
                   R[2]*L(0), R[2]*L(1) );

    const Real Z10 = Y[0]*~RR[0], Z11 = Y[0]*~RR[1],
               Z20 = Y[1]*~RR[0], Z21 = Y[1]*~RR[1], Z22= Y[1]*~RR[2];
    const Real Z00 = (L(0,0)+L(1,1)) - (Z11+Z22);

    const Vec3 Rv( R(0,1)*e-R(0,0)*f,
                   R(1,1)*e-R(1,0)*f,
                   R(2,1)*e-R(2,0)*f );

    return SymMat33( Z00 + c,
                     Z10 + Rv[2], Z11 + c,
                     Z20 - Rv[1], Z21 + Rv[0], Z22 + c );
}

// See above method for details. This method is identical except that
// the layout of the matrix used to store the rotation matrix has
// changed. Note that all the indexing here is identical to the normal
// case above; but the rotation matrix elements are drawn from different
// memory locations so that the net effect is to use the transpose of
// the original rotation from which this was created.
SymMat33 InverseRotation::reexpressSymMat33(const SymMat33& S_BB) const {
    const Real a=S_BB(0,0), b=S_BB(1,1), c=S_BB(2,2);
    const Real d=S_BB(1,0), e=S_BB(2,0), f=S_BB(2,1);
    // Note reversal of row and column spacing here (normal is 3,1).
    const Mat<3,3,Real,1,3>& R   = this->asMat33();
    const Mat<3,2,Real,1,3>& RR  = R.getSubMat<3,2>(0,0); // first two columns of R

    const Mat32 L( a-c ,  d,
                    d  , b-c,
                   2*e , 2*f );

    const Mat22 Y( R[1]*L(0), R[1]*L(1),
                   R[2]*L(0), R[2]*L(1) );

    const Real Z10 = Y[0]*~RR[0], Z11 = Y[0]*~RR[1],
               Z20 = Y[1]*~RR[0], Z21 = Y[1]*~RR[1], Z22= Y[1]*~RR[2];
    const Real Z00 = (L(0,0)+L(1,1)) - (Z11+Z22);

    const Vec3 Rv( R(0,1)*e-R(0,0)*f,
                   R(1,1)*e-R(1,0)*f,
                   R(2,1)*e-R(2,0)*f );

    return SymMat33( Z00 + c,
                     Z10 + Rv[2], Z11 + c,
                     Z20 - Rv[1], Z21 + Rv[0], Z22 + c );
}


//------------------------------------------------------------------------------
std::ostream&  operator<<( std::ostream& o, const Rotation& m )  { return o << m.asMat33(); }



//------------------------------------------------------------------------------
}  // End of namespace SimTK

