//-----------------------------------------------------------------------------
// File:     RotationTest.cpp
// Class:    None
// Parent:   None
// Purpose:  Test 3x3 Rotation class relating two right-handed orthogonal bases
//-----------------------------------------------------------------------------

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Paul Mitiguy                                                      *
 * Contributors: Michael Sherman                                              *
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
 * Tests for the classes defined in Rotation.h.
 */

//-----------------------------------------------------------------------------
#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include <iostream>
using std::cout; using std::endl;

//-----------------------------------------------------------------------------
using namespace SimTK;

bool  doRequiredTasks();
void  WriteStringToScreen( const char outputString[] )  { std::cout << outputString; }  // fputs( outputString,stdout ); if include<stdio>

//-------------------------------------------------------------------
int main() {

   // Default value is program failed
   bool programSucceeded = false;

   // It is a good programming practice to do little in the main function of a program.
   // The try-catch code in this main routine catches exceptions thrown by functions in the
   // try block, e.g., catching an exception that occurs when a NULL pointer is de-referenced.
   try {
      // Do the required programming tasks
      programSucceeded = doRequiredTasks();
   }
   // This catch statement handles certain types of exceptions
   catch( const std::exception& e ) {
      WriteStringToScreen( "\n\n Error: Programming error encountered.\n The exception thrown is: " );
      WriteStringToScreen( e.what() );
      WriteStringToScreen( "\n\n" );
   }
   // The exception-declaration statement (...) handles any type of exception,
   // including C exceptions and system/application generated exceptions.
   // This includes exceptions such as memory protection and floating-point violations.
   // An ellipsis catch handler must be the last handler for its try block.
   catch( ... ) {
      WriteStringToScreen( "\n\n Error: Programming error encountered.\n An unhandled exception was thrown.\n\n" );
   }

   if (!programSucceeded)
       WriteStringToScreen("\nError: a test failed.\n");

   // The value returned by the main function is the exit status of the program.
   // A normal program exit returns 0 (other return values usually signal an error).
   return programSucceeded == true ? 0 : 1;
}


//-------------------------------------------------------------------
// Prototypes for methods in this test
//-------------------------------------------------------------------
bool  testRotationOneAxis( const Real angle, const CoordinateAxis& axis );
bool  testRotationTwoAxes( const BodyOrSpaceType bodyOrSpace, const Real angle1, const CoordinateAxis& axis1, const Real angle2, const CoordinateAxis &axis2 );
bool  testRotationThreeAxes( const BodyOrSpaceType bodyOrSpace, const Real angle1, const CoordinateAxis& axis1, const Real angle2, const CoordinateAxis &axis2, const Real angle3, const CoordinateAxis &axis3 );
bool  testQuaternion( Real e0, Real e1, Real e2, Real e3 );

bool  testInverseRotation1Angle( Real angle, Real theta );
bool  testInverseRotation2Angle( Real angle1, Real theta1,  Real angle2, Real theta2 );
bool  testInverseRotation3AngleTwoAxes( Real angle1, Real theta1,  Real angle2, Real theta2,  Real angle3, Real theta3 );
bool  testInverseRotation3AngleThreeAxes( Real angle1, Real theta1,  Real angle2, Real theta2,  Real angle3, Real theta3 );


bool  exhaustiveTestof1AngleRotation();
bool  exhaustiveTestof2AngleRotation();
bool  exhaustiveTestof3AngleRotation();
bool  exhaustiveTestof3AngleTwoAxesRotationNearSingularity();
bool  exhaustiveTestof3AngleThreeAxesRotationNearSingularity();
bool  exhaustiveTestofQuaternions();

bool testRotationFromTwoGivenAxes( const Vec3& vi, const CoordinateAxis& ai, const Vec3& vj, const CoordinateAxis& aj);
bool testReexpressSymMat33();

//-------------------------------------------------------------------
bool  doRequiredTasks( ) {

    // Use the next Rotation to test against (simple theta-lambda rotation)
    Rotation testRotation;

    // Test default constructor
    Rotation defaultRotationConstructor;
    testRotation.setRotationToIdentityMatrix();
    bool test = defaultRotationConstructor.areAllRotationElementsSameToMachinePrecision( testRotation );

    // Test copy constructor
    testRotation.setRotationFromAngleAboutNonUnitVector( 1.0, Vec3(0.2, 0.4, 0.6) );
    Rotation rotationCopyConstructor( testRotation );
    test = test && rotationCopyConstructor.areAllRotationElementsSameToMachinePrecision( testRotation );

    // Test operator =
    Rotation rotationOperatorEqual = testRotation;
    test = test && rotationOperatorEqual.areAllRotationElementsSameToMachinePrecision( testRotation );

    // Test rotation by angle about arbitrary CoordinateAxis
    testRotation.setRotationFromAngleAboutNonUnitVector( 0.1, Vec3(1.0, 0.0, 0.0) );
    CoordinateAxis coordAxis = XAxis;
    Rotation rotationCoordAxis( 0.1, coordAxis );
    test = test && rotationCoordAxis.areAllRotationElementsSameToMachinePrecision( testRotation );

    Real testTheta = rotationCoordAxis.convertOneAxisRotationToOneAngle( coordAxis );
    test = test && fabs(0.1 - testTheta) < 10*SignificantReal;

    // Test rotation by angle about XAxis, YAxis, ZAxis
    test = test && testRotationOneAxis(  0.2, XAxis );
    test = test && testRotationOneAxis( -0.2, XAxis );
    test = test && testRotationOneAxis(  2.1, YAxis );
    test = test && testRotationOneAxis( -2.1, YAxis );
    test = test && testRotationOneAxis(  3.1, ZAxis );
    test = test && testRotationOneAxis( -3.1, ZAxis );

    // Test rotation with two angles and two axes XX, XY, XZ
    test = test && testRotationTwoAxes(  BodyRotationSequence,  0.2, XAxis, 0.3, XAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence,  0.2, XAxis, 0.3, XAxis );
    test = test && testRotationTwoAxes(  BodyRotationSequence,  1.2, XAxis,-1.3, YAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence,  1.2, XAxis,-1.3, YAxis );
    test = test && testRotationTwoAxes(  BodyRotationSequence, -3.1, XAxis, 1.2, ZAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence, -3.1, XAxis, 1.2, ZAxis );

    // Test rotation with two angles and two axes YX, YY, YZ
    test = test && testRotationTwoAxes(  BodyRotationSequence,  1.2, YAxis, 0.3, XAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence,  1.2, YAxis, 0.3, XAxis );
    test = test && testRotationTwoAxes(  BodyRotationSequence,  2.2, YAxis,-1.3, YAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence,  2.2, YAxis,-1.3, YAxis );
    test = test && testRotationTwoAxes(  BodyRotationSequence, -3.1, YAxis, 1.2, ZAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence, -3.1, YAxis, 1.2, ZAxis );

    // Test rotation with two angles and two axes ZX, ZY, ZZ
    test = test && testRotationTwoAxes(  BodyRotationSequence,  1.2, ZAxis, 0.3, XAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence,  1.2, ZAxis, 0.3, XAxis );
    test = test && testRotationTwoAxes(  BodyRotationSequence,  2.2, ZAxis,-1.3, YAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence,  2.2, ZAxis,-1.3, YAxis );
    test = test && testRotationTwoAxes(  BodyRotationSequence, -3.1, ZAxis, 1.2, ZAxis );
    test = test && testRotationTwoAxes( SpaceRotationSequence, -3.1, ZAxis, 1.2, ZAxis );

	// Test the construction of rotations from two given axes.
	const UnitVec3	vi(0.01, 0.02, .9);
	const Vec3		vj(-0.5, 0.5, 0.2);

	test = test && testRotationFromTwoGivenAxes(vi, XAxis, vj, YAxis);
	test = test && testRotationFromTwoGivenAxes(vi, YAxis, vj, XAxis);

	test = test && testRotationFromTwoGivenAxes(vi, YAxis, vj, ZAxis);
	test = test && testRotationFromTwoGivenAxes(vi, ZAxis, vj, YAxis);

	test = test && testRotationFromTwoGivenAxes(vi, ZAxis, vj, XAxis);
	test = test && testRotationFromTwoGivenAxes(vi, XAxis, vj, ZAxis);

    // Test rotation with three angles and three axes XXX, XXY, XXZ, XYX, XYY, XYZ, XZX, XZY, XZZ
    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, XAxis, 0.3, XAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, XAxis, 0.3, XAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, XAxis,-1.3, XAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, XAxis,-1.3, XAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, XAxis, 1.2, XAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, XAxis, 1.2, XAxis, 1.3, ZAxis );

    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, XAxis, 0.3, YAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, XAxis, 0.3, YAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, XAxis,-1.3, YAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, XAxis,-1.3, YAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, XAxis, 1.2, YAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, XAxis, 1.2, YAxis, 1.3, ZAxis );

    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, XAxis, 0.3, ZAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, XAxis, 0.3, ZAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, XAxis,-1.3, ZAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, XAxis,-1.3, ZAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, XAxis, 1.2, ZAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, XAxis, 1.2, ZAxis, 1.3, ZAxis );

    // Test rotation with three angles and three axes YXX, YXY, YXZ, YYX, YYY, YYZ, YZX, YZY, YZZ
    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, YAxis, 0.3, XAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, YAxis, 0.3, XAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, YAxis,-1.3, XAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, YAxis,-1.3, XAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, YAxis, 1.2, XAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, YAxis, 1.2, XAxis, 1.3, ZAxis );

    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, YAxis, 0.3, YAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, YAxis, 0.3, YAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, YAxis,-1.3, YAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, YAxis,-1.3, YAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, YAxis, 1.2, YAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, YAxis, 1.2, YAxis, 1.3, ZAxis );

    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, YAxis, 0.3, ZAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, YAxis, 0.3, ZAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, YAxis,-1.3, ZAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, YAxis,-1.3, ZAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, YAxis, 1.2, ZAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, YAxis, 1.2, ZAxis, 1.3, ZAxis );

    // Test rotation with three angles and three axes ZXX, ZXY, ZXZ, ZYX, ZYY, ZYZ, ZZX, ZZY, ZZZ
    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, ZAxis, 0.3, XAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, ZAxis, 0.3, XAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, ZAxis,-1.3, XAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, ZAxis,-1.3, XAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, ZAxis, 1.2, XAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, ZAxis, 1.2, XAxis, 1.3, ZAxis );

    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, ZAxis, 0.3, YAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, ZAxis, 0.3, YAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, ZAxis,-1.3, YAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, ZAxis,-1.3, YAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, ZAxis, 1.2, YAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, ZAxis, 1.2, YAxis, 1.3, ZAxis );

    test = test && testRotationThreeAxes(  BodyRotationSequence,  0.2, ZAxis, 0.3, ZAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  0.2, ZAxis, 0.3, ZAxis, 0.4, XAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence,  1.2, ZAxis,-1.3, ZAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence,  1.2, ZAxis,-1.3, ZAxis,-1.4, YAxis );
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.1, ZAxis, 1.2, ZAxis, 1.3, ZAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.1, ZAxis, 1.2, ZAxis, 1.3, ZAxis );

	// This used to fail
    test = test && testRotationThreeAxes(  BodyRotationSequence, -3.2288591161895095, XAxis, -3.1415926535897931, YAxis, -3.1415926535897931, XAxis );
    test = test && testRotationThreeAxes( SpaceRotationSequence, -3.2288591161895095, XAxis, -3.1415926535897931, YAxis, -3.1415926535897931, XAxis );

    // Test Rotation quaternion methods.
    test = test && testQuaternion( 0.5, 0.1, 0.2,  0.3 );
    test = test && testQuaternion(-0.5, 0.1, 0.2, -0.3 );

    // Test construction of nearby orthogonal rotation matrix from a generic Mat33.
    Rotation nearbyRotation( testRotation.asMat33() );
    test = test && nearbyRotation.areAllRotationElementsSameToMachinePrecision( testRotation );

    // Exhaustive test of 1-angle, 2-angle, and 3-angle rotations
    test = test && exhaustiveTestof1AngleRotation();
    test = test && exhaustiveTestof2AngleRotation();
    test = test && exhaustiveTestof3AngleRotation();

    // Exhaustive test of 3-angle rotations near singularity
    test = test && exhaustiveTestof3AngleTwoAxesRotationNearSingularity();
    test = test && exhaustiveTestof3AngleThreeAxesRotationNearSingularity();

    // Exhaustive test of Quaterions
    test = test && exhaustiveTestofQuaternions();

    // Test out special code for rotating symmetric matrices.
    test = test && testReexpressSymMat33();

    return test;
}


//-------------------------------------------------------------------
bool  testRotationOneAxis( const Real angle, const CoordinateAxis& axis ) {

    // Form rotation about specified axis
    Rotation rotationSpecified;
    if( axis == XAxis )  rotationSpecified.setRotationFromAngleAboutX( angle );
    if( axis == YAxis )  rotationSpecified.setRotationFromAngleAboutY( angle );
    if( axis == ZAxis )  rotationSpecified.setRotationFromAngleAboutZ( angle );

    // Form equivalent rotation by another means
    Real unitX = axis == XAxis ? 1 : 0;
    Real unitY = axis == YAxis ? 1 : 0;
    Real unitZ = axis == ZAxis ? 1 : 0;
    UnitVec3 unitVector( unitX, unitY, unitZ );
    Rotation testRotation( angle, unitVector );

    // Test to see if they are the same
    bool test = rotationSpecified.areAllRotationElementsSameToMachinePrecision( testRotation );

    // Do the inverse problem to back out the angle
    const Real theta = rotationSpecified.convertOneAxisRotationToOneAngle( axis );

    // Create a Rotation matrix with the backed-out angle and compare to the original Rotation matrix
    testRotation.setRotationFromAngleAboutAxis( theta, axis );
    test = test && rotationSpecified.areAllRotationElementsSameToMachinePrecision( testRotation );

	// Conversion should produce  angle = theta   if  angle is in proper range (-pi < angle <= pi)
	test = test && testInverseRotation1Angle( angle, theta );

    return test;
}


//-------------------------------------------------------------------
bool  testRotationTwoAxes( const BodyOrSpaceType bodyOrSpace, const Real angle1, const CoordinateAxis& axis1, const Real angle2, const CoordinateAxis &axis2 ) {

    // Form rotation about specified axes
    Rotation rotationSpecified( bodyOrSpace, angle1, axis1, angle2, axis2 );

    // Form equivalent rotation by another means
    Rotation AB( angle1, axis1 );
    Rotation BC( angle2, axis2 );
    Rotation testRotation = (bodyOrSpace == BodyRotationSequence) ? AB * BC : BC * AB;

    // Test to see if they are the same
    bool test = rotationSpecified.areAllRotationElementsSameToMachinePrecision( testRotation );

    // Do the inverse problem to back out the angles
    const Vec2 testVec = rotationSpecified.convertTwoAxesRotationToTwoAngles( bodyOrSpace, axis1, axis2 );
    const Real theta1 = testVec[0];
    const Real theta2 = testVec[1];

    // Create a Rotation matrix with the backed-out angles and compare to the original Rotation matrix
    testRotation.setRotationFromTwoAnglesTwoAxes( bodyOrSpace, theta1, axis1, theta2, axis2 );
    test = test && rotationSpecified.areAllRotationElementsSameToMachinePrecision( testRotation );

	// Conversion should produce same angles for for appropriate ranges of angle1 and angle2
    if( axis1.isSameAxis(axis2) ) 
       test = test && testInverseRotation1Angle( angle1+angle2, theta1+theta2 );
	else
       test = test && testInverseRotation2Angle( angle1,theta1, angle2,theta2 );

	return test;
}


//-------------------------------------------------------------------
bool  testRotationThreeAxes( const BodyOrSpaceType bodyOrSpace, const Real angle1, const CoordinateAxis& axis1, const Real angle2, const CoordinateAxis &axis2, const Real angle3, const CoordinateAxis &axis3 ) {

    // Form rotation about specified axes
    Rotation rotationSpecified( bodyOrSpace, angle1, axis1, angle2, axis2, angle3, axis3 );

    // Form equivalent rotation by another means
    Rotation AB( angle1, axis1 );
    Rotation BC( angle2, axis2 );
    Rotation CD( angle3, axis3 );
    Rotation testRotation = (bodyOrSpace == BodyRotationSequence) ? AB * BC * CD :  CD * BC * AB;

    // Test to see if they are the same
    bool test = rotationSpecified.areAllRotationElementsSameToMachinePrecision( testRotation );

    // Do the inverse problem to back out the angles
    const Vec3 testVec = rotationSpecified.convertThreeAxesRotationToThreeAngles( bodyOrSpace, axis1, axis2, axis3 );
    const Real theta1 = testVec[0];
    const Real theta2 = testVec[1];
    const Real theta3 = testVec[2];

    // Create a Rotation matrix with the backed-out angles and compare to the original Rotation matrix
    testRotation.setRotationFromThreeAnglesThreeAxes( bodyOrSpace, theta1, axis1, theta2, axis2, theta3, axis3 );
    test = test && rotationSpecified.areAllRotationElementsSameToMachinePrecision( testRotation );

	// Conversion should produce same angles for for appropriate ranges of angle1 and angle2
	if( axis1.areAllSameAxes(axis2,axis3) ) 
       test = test && testInverseRotation1Angle( angle1+angle2+angle3, theta1+theta2+theta3 );
    else if( axis1.isSameAxis(axis2) ) 
       test = test && testInverseRotation2Angle( angle1+angle2, theta1+theta1, angle3,theta3 );
	else if( axis2.isSameAxis(axis3) ) 
       test = test && testInverseRotation2Angle( angle1,theta1, angle2+angle3, theta2+theta3 );
	else if( axis1.isSameAxis(axis3) )
       test = test && testInverseRotation3AngleTwoAxes( angle1,theta1, angle2,theta2, angle3,theta3 );
	else 
       test = test && testInverseRotation3AngleThreeAxes( angle1,theta1, angle2,theta2, angle3,theta3 );

	return test;
}


//-------------------------------------------------------------------
bool  testInverseRotation1Angle( Real angle, Real theta ) {
    bool test = true;
	bool angleInProperRange = ( -SimTK_PI <= angle  &&  angle <= SimTK_PI );
    if( angleInProperRange )
      test = fabs( angle - theta ) < 10*SignificantReal;
	return test;
}


//-------------------------------------------------------------------
bool  testInverseRotation2Angle( Real angle1, Real theta1,  Real angle2, Real theta2 ) {
   bool test = true;
   bool angle1InProperRange = ( -SimTK_PI <= angle1  &&  angle1 <= SimTK_PI );
   bool angle2InProperRange = ( -SimTK_PI <= angle2  &&  angle2 <= SimTK_PI );
   if( angle1InProperRange && angle2InProperRange ) {
       test = test && fabs( angle1 - theta1 ) < 10*SignificantReal;
       test = test && fabs( angle2 - theta2 ) < 10*SignificantReal;
   }
   return test;
}


//-------------------------------------------------------------------
bool  testInverseRotation3AngleTwoAxes( Real angle1, Real theta1,  Real angle2, Real theta2,  Real angle3, Real theta3 ) {
    bool test = true;
	bool angle1InProperRange = (    -SimTK_PI <= angle1  &&  angle1 <= SimTK_PI );
    bool angle2InProperRange = (            0 <= angle2  &&  angle2 <= SimTK_PI );
    bool angle3InProperRange = (    -SimTK_PI <= angle3  &&  angle3 <= SimTK_PI );
    if( angle1InProperRange && angle2InProperRange && angle3InProperRange ) {
       test = test && fabs( angle1 - theta1 ) < 10*SignificantReal;
       test = test && fabs( angle2 - theta2 ) < 10*SignificantReal;
       test = test && fabs( angle3 - theta3 ) < 10*SignificantReal;

	   // Test needs to be modified if near singularity
	   const Real singularity = 0.0;
	   if( test == false && fabs(angle2-singularity) <= SignificantReal ) {
	      const Real angle1PlusAngle3 = angle1 + angle3;
	      const Real theta1PlusTheta3 = theta1 + theta3;
	      bool angleSumInProperRange = ( -SimTK_PI <= angle1PlusAngle3  &&  angle1PlusAngle3 <= SimTK_PI );
		  if( angleSumInProperRange == false ) test = true;
		  else {
              test = fabs( angle2 - theta2 ) < 10*SignificantReal;
              test = test && fabs( angle1PlusAngle3 - theta1PlusTheta3 ) < 10*SignificantReal;
		  }
	   }
    }
	return test;
}


//-------------------------------------------------------------------
bool  testInverseRotation3AngleThreeAxes( Real angle1, Real theta1,  Real angle2, Real theta2,  Real angle3, Real theta3 ) {
    bool test = true;
	bool angle1InProperRange = (    -SimTK_PI <= angle1  &&  angle1 <=     SimTK_PI );
    bool angle2InProperRange = (-0.5*SimTK_PI <= angle2  &&  angle2 <= 0.5*SimTK_PI );
    bool angle3InProperRange = (    -SimTK_PI <= angle3  &&  angle3 <=     SimTK_PI );
    if( angle1InProperRange && angle2InProperRange && angle3InProperRange ) {
       test = test && fabs( angle1 - theta1 ) < 10*SignificantReal;
       test = test && fabs( angle2 - theta2 ) < 10*SignificantReal;
       test = test && fabs( angle3 - theta3 ) < 10*SignificantReal;

	   // Test needs to be modified if near singularity
	   const Real singularity = 0.5*SimTK_PI;
	   if( test == false && fabs(angle2-singularity) <= SignificantReal ) {
	      const Real angle1PlusAngle3 = angle1 + angle3;
	      const Real theta1PlusTheta3 = theta1 + theta3;
	      bool angleSumInProperRange = ( -SimTK_PI <= angle1PlusAngle3  &&  angle1PlusAngle3 <= SimTK_PI );
		  if( angleSumInProperRange == false ) test = true;
		  else {
              test = fabs( angle2 - theta2 ) < 10*SignificantReal;
              test = test && fabs( angle1PlusAngle3 - theta1PlusTheta3 ) < 10*SignificantReal;
		  }
	   }
    }
	return test;
}


//-------------------------------------------------------------------
bool  testQuaternion( Real e0, Real e1, Real e2, Real e3 )
{
    // Construct quaternion and normalize it
    Quaternion qe0e1e2e3( e0, e1, e2, e3 );

    // Convert quaternion to a Rotation matrix
    Rotation rotationSpecified( qe0e1e2e3 );

    // Convert Rotation back to quaternion
    Quaternion qTest = rotationSpecified.convertRotationToQuaternion();

    // Convert quaternion to a Rotation matrix
    Rotation testRotation;  testRotation.setRotationFromQuaternion( qTest );

    // Test to see if they are the same
    bool test = rotationSpecified.areAllRotationElementsSameToMachinePrecision( testRotation );

    return test;
}


//-------------------------------------------------------------------
bool  exhaustiveTestof1AngleRotation( ) {
   bool test = true;

   // Range to check angles
   Real negativeStartAngle = convertDegreesToRadians( -385 );
   Real positiveStartAngle = convertDegreesToRadians(  385 );
   Real incrementAngle = convertDegreesToRadians( 0.5 );

   // Test each axis
   for( int i=0;  i<=2;  i++ ) {
      CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);
      for( Real anglei = negativeStartAngle;  anglei < positiveStartAngle;  anglei += incrementAngle )
         test = test && testRotationOneAxis(  anglei, axisi );
   }
   return test;
}


//-------------------------------------------------------------------
bool  exhaustiveTestof2AngleRotation( ) {
   bool test = true;

   // Range to check angles
   Real negativeStartAngle = convertDegreesToRadians( -200 );
   Real positiveStartAngle = convertDegreesToRadians(  200 );
   Real incrementAngle = convertDegreesToRadians( 10.0 );

   // Test each axis
   for( int i=0;  i<=2;  i++ ) {
      CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);

      for( int j=0;  j<=2;  j++ ) {
         CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);

         for( Real anglei = negativeStartAngle;  anglei < positiveStartAngle;  anglei += incrementAngle ) 
         for( Real anglej = negativeStartAngle;  anglej < positiveStartAngle;  anglej += incrementAngle ) {
               test = test && testRotationTwoAxes(  BodyRotationSequence,  anglei, axisi, anglej, axisj );
               test = test && testRotationTwoAxes( SpaceRotationSequence,  anglei, axisi, anglej, axisj );
         }
      }
   }
   return test;
}


//-------------------------------------------------------------------
bool  exhaustiveTestof3AngleRotation( ) {
   bool test = true;

   // Range to check angles
   Real negativeStartAngle = convertDegreesToRadians( -200 );
   Real positiveStartAngle = convertDegreesToRadians(  200 );
   Real incrementAngle = convertDegreesToRadians( 40.0 );

   // Test each axis
   for( int i=0;  i<=2;  i++ ) {
      CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);

      for( int j=0;  j<=2;  j++ ) {
         CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);

         for( int k=0;  k<=2;  k++ ) {
            CoordinateAxis axisk = CoordinateAxis::getCoordinateAxis(k);

            for( Real anglei = negativeStartAngle;  anglei < positiveStartAngle;  anglei += incrementAngle ) 
            for( Real anglej = negativeStartAngle;  anglej < positiveStartAngle;  anglej += incrementAngle ) 
            for( Real anglek = negativeStartAngle;  anglek < positiveStartAngle;  anglek += incrementAngle ) {
               test = test && testRotationThreeAxes(  BodyRotationSequence,  anglei, axisi,  anglej, axisj,  anglek, axisk );
               test = test && testRotationThreeAxes( SpaceRotationSequence,  anglei, axisi,  anglej, axisj,  anglek, axisk );
			}
         }
      }
   }
   return test;
}


//-------------------------------------------------------------------
bool  exhaustiveTestof3AngleTwoAxesRotationNearSingularity() {
   bool test = true;

   // Range to check angles anglei and anglek
   Real negativeStartAngle = convertDegreesToRadians( -200 );
   Real positiveStartAngle = convertDegreesToRadians(  200 );
   Real incrementStartAngle = convertDegreesToRadians( 20.0 );

   // Range to check around singularity 
   Real negativeSinglAngle = convertDegreesToRadians( 0 - 1.0E-11 );
   Real positiveSinglAngle = convertDegreesToRadians( 0 + 1.0E-11 );
   Real incrementSinglAngle = convertDegreesToRadians( 1.0E-12 );

   // Test each axis
   for( int i=0;  i<=2;  i++ ) {
      CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);

      for( int j=0;  j<=2;  j++ ) {
         CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);

         for( int k=0;  k<=2;  k++ ) {
            CoordinateAxis axisk = CoordinateAxis::getCoordinateAxis(k);

            for( Real anglei = negativeStartAngle;  anglei < positiveStartAngle;  anglei += incrementStartAngle ) 
            for( Real anglej = negativeSinglAngle;  anglej < positiveSinglAngle;  anglej += incrementSinglAngle ) 
            for( Real anglek = negativeStartAngle;  anglek < positiveStartAngle;  anglek += incrementStartAngle ) {
               test = test && testRotationThreeAxes(  BodyRotationSequence,  anglei, axisi,  anglej, axisj,  anglek, axisk );
               test = test && testRotationThreeAxes( SpaceRotationSequence,  anglei, axisi,  anglej, axisj,  anglek, axisk );
			}
         }
      }
   }
   return test;
}


//-------------------------------------------------------------------
bool  exhaustiveTestof3AngleThreeAxesRotationNearSingularity() {
   bool test = true;

   // Range to check angles anglei and anglek
   Real negativeStartAngle = convertDegreesToRadians( -200 );
   Real positiveStartAngle = convertDegreesToRadians(  200 );
   Real incrementStartAngle = convertDegreesToRadians( 20.0 );

   // Range to check around singularity 
   Real negativeSinglAngle = convertDegreesToRadians( 90 - 1.0E-11 );
   Real positiveSinglAngle = convertDegreesToRadians( 90 + 1.0E-11 );
   Real incrementSinglAngle = convertDegreesToRadians( 1.0E-12 );

   // Test each axis
   for( int i=0;  i<=2;  i++ ) {
      CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);

      for( int j=0;  j<=2;  j++ ) {
         CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);

         for( int k=0;  k<=2;  k++ ) {
            CoordinateAxis axisk = CoordinateAxis::getCoordinateAxis(k);

            for( Real anglei = negativeStartAngle;  anglei < positiveStartAngle;  anglei += incrementStartAngle ) 
            for( Real anglej = negativeSinglAngle;  anglej < positiveSinglAngle;  anglej += incrementSinglAngle ) 
            for( Real anglek = negativeStartAngle;  anglek < positiveStartAngle;  anglek += incrementStartAngle ) {
               test = test && testRotationThreeAxes(  BodyRotationSequence,  anglei, axisi,  anglej, axisj,  anglek, axisk );
               test = test && testRotationThreeAxes( SpaceRotationSequence,  anglei, axisi,  anglej, axisj,  anglek, axisk );
			}
         }
      }
   }
   return test;
}


//-------------------------------------------------------------------
bool  exhaustiveTestofQuaternions() {
   bool test = true;

   for( Real e0 = -1;  e0 <= 1;  e0 += 0.2 )
   for( Real e1 = -1;  e1 <= 1;  e1 += 0.2 )
   for( Real e2 = -1;  e2 <= 1;  e2 += 0.2 )
   for( Real e3 = -1;  e3 <= 1;  e3 += 0.2 )
     test = test && testQuaternion( e0, e1, e2, e3 );

   return test;
}

//-------------------------------------------------------------------
// Here we need to check that the resulting rotation has determinant 1 (a previous bug
// left it with determinant -1) and that the j'th axis is at least pointing in the
// direction of the given vj.
bool testRotationFromTwoGivenAxes( const Vec3& vi, const CoordinateAxis& ai, const Vec3& vj, const CoordinateAxis& aj) {
	bool test = true;

	// This makes a Rotation with vi as axis i, but axis j will only be in the general direction of vj.
	const Rotation testRotation(UnitVec3(vi), ai, vj, aj);

	test = test && std::fabs(det(testRotation) - 1) <= SignificantReal;

	test = test && dot(testRotation(ai), vi) > 0;
	test = test && dot(testRotation(aj), vj) > 0;

	return test;
}

//-------------------------------------------------------------------
// Test the tricked-out code used to rotation a symmetric dyadic
// matrix S_AA=R_AB*S_BB*R_BA.
bool testReexpressSymMat33() {
    bool test = true;

    const Rotation R_AB(Test::randRotation());
    const SymMat33 S_BB(Test::randSymMat<3>());
    const Mat33 M_BB(S_BB);

    const SymMat33 S_AA = R_AB.reexpressSymMat33(S_BB);
    const Mat33 M_AA = R_AB*M_BB*~R_AB;

    const Mat33 MS_AA(S_AA);
    test = test && (MS_AA-M_AA).norm() <= SignificantReal;

    // Now put it back with an InverseRotation.
    const SymMat33 isS_BB = (~R_AB).reexpressSymMat33(S_AA);
    test = test && (S_BB-isS_BB).norm() <= SignificantReal;

    const Rotation I; // identity
    const SymMat33 S_BB_still = I.reexpressSymMat33(S_BB);
    test = test && (S_BB_still-S_BB).norm() <= SignificantReal;

    // Test symmetric matrix multiply (doesn't belong here).
    //const SymMat33 S1(Test::randSymMat<3>()), S2(Test::randSymMat<3>());
    //const Mat33 M1(S1), M2(S2);
    //const SymMat33 S(S1*S2);
    //const Mat33 M(M1*M2);
    //cout << "S=" << S << endl;
    //cout << "M=" << M << endl;


    return test;
}
