
#include "SimTKsimbody.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>
using namespace std;
using namespace SimTK;

static const Transform GroundFrame;

int main( int numberOfCommandLineArguments, char** arrayOfCommandLineArguments ) 
{
	// Declare a multibody system (contains sub-systems) 
    MultibodySystem  mbs;

	// Create a gravity vector that is straight down
    // Create a uniform gravity sub-system
	// Add the gravity sub-system to the multibody system
	Vec3 gravityVector( 0,-9.8, 0 );
    UniformGravitySubsystem gravity( gravityVector );
    mbs.addForceSubsystem( gravity );

	// Create a matter sub-system (the apple)
    SimbodyMatterSubsystem  apple;

	// Create the mass properties of this matter sub-system (mass properties of the apple).
	// The MassProperties class holds the mass, center of mass, and inertia properties of a rigid body B.
	// The location of the center of mass is a vector from B's origin expressed in the B frame.  
	// The inertia matrix is of the body B about B's origin for the unit vectors fixed in B's frame.
	const Real  massOfApple = 2;
    const Vec3  appleCenterOfMassLocation(0, 0, 0);  // center of the apple
	const Inertia  appleInertiaMatrix( 1, 1, 1, 0, 0, 0 );
    MassProperties  appleMassProperties( massOfApple, appleCenterOfMassLocation, appleInertiaMatrix ); 

	// Define the apple and how it is connected to the ground (it is free to move relative to ground)
	const Vec3  inboardJointLocation(0,0,0);
    const BodyId appleBodyNumber = apple.addRigidBody( 
        appleMassProperties, Transform(Rotation::aboutX(NTraits<Real>::Pi/4), inboardJointLocation), 
        GroundId, Transform(Rotation::aboutX(NTraits<Real>::Pi/4)), Mobilizer::/*Free*/Cartesian() );

    // Add the matter (apple) sub-system to the system.
    mbs.setMatterSubsystem( apple );

	// Create a state for this system
    State s;
    mbs.realize(s); // define appropriate states for this System

    // Create a study using the Runge Kutta Merson integrator
    RungeKuttaMerson myStudy(mbs, s);
    //CPodesIntegrator myStudy(mbs, s);

	// Set the numerical accuracy for the integrator
	myStudy.setAccuracy( 1.0E-7 );

	// Set the initial numerical integration step and the time for the simulation to run 
    const Real dt = 0.01; 
    const Real finalTime = 4.0;

	// Open a file to record the simulation results (they are also displayed on screen)
	FILE *outputFile = fopen( "NewtonsAppleForPlottingResults.txt", "w" );
    fprintf( outputFile, "time      yPosition      mechanicalEnergy\n" );
	fprintf( stdout,     "time      yPosition      mechanicalEnergy\n" );

	// I am not sure what the next two lines do
    Transform initialApplePosition(Rotation(), Vec3(0,10,0));
    apple.setMobilizerTransform(s, appleBodyNumber, initialApplePosition);
    //apple.setMobilizerQ(s, appleBodyNumber, 5, 10.);
    s.updTime() = 0;

    myStudy.initialize();

    Matrix m(9,3); m=1;
    Matrix n(9,3); n=1;
    cout << "m+n=" << m+n << " norm(m+n)=" << (m+n).norm() << endl;

    Vec3 v(1), w(2);
    cout << "v.dot(w))=" << ~v*w << endl;

	// Run the simulation and print the results
    while( 1 )
	{
       // Query Simbody for results to be printed
       Real time = s.getTime();
	   Real kineticEnergy = mbs.getKineticEnergy(s); 
	   Real uniformGravitationalPotentialEnergy = mbs.getPotentialEnergy(s);
	   Real mechanicalEnergy = kineticEnergy + uniformGravitationalPotentialEnergy;

       // Transform nearly always has position vector from ground to point of interest
	   // Also contains a rotation matrix from ground to frame of interest
       const Transform& bodyTransform = apple.getBodyTransform( s, appleBodyNumber );

	   // .T returns translation and .R returns 3x3 rotation matrix
	   const Vec3& bodyPosition = bodyTransform.T();

	   // Gets the y position 0,
	   Real yPosition = bodyPosition[1];

	   // Print results to screen and file
	   fprintf( stdout,     "%7.3E     %10.4E     %10.8g\n", time, yPosition, mechanicalEnergy );
	   fprintf( outputFile, "%7.3E     %10.4E     %10.8g\n", time, yPosition, mechanicalEnergy );

       // Check if integration has completed
       if( time >= finalTime ) break;

	   // Increment time step
       myStudy.step( time + dt);
    }
} 
