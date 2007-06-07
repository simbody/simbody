//-----------------------------------------------------------------------------
// File:     SpinningBook.cpp
// Class:    None
// Parent:   None
// Children: None
// Purpose:  Simulates a spinning book
//-----------------------------------------------------------------------------
// The following are standard C/C++ header files.
// If a filename is enclosed inside < >  it means the header file is in the Include directory.
// If a filename is enclosed inside " "  it means the header file is in the current directory.
#include <ctype.h>      // Character Types
#include <math.h>       // Mathematical Constants
#include <stdarg.h>     // Variable Argument Lists
#include <stdio.h>      // Standard Input/Output Functions
#include <stdlib.h>     // Utility Functions
#include <string.h>     // String Operations
#include <signal.h>     // Signals (Contol-C + Unix System Calls)
#include <setjmp.h>     // Nonlocal Goto (For Control-C)
#include <time.h>       // Time and Date information
#include <assert.h>     // Verify Program Assertion
#include <errno.h>      // Error Codes (Used in Unix system())
#include <float.h>      // Floating Point Constants
#include <limits.h>     // Implementation Constants
#include <stddef.h>     // Standard Definitions
#include <exception>    // Exception handling (e.g., try, catch throw)
//-----------------------------------------------------------------------------
#include "SimTKsimbody.h"
using namespace SimTK;
using namespace std;
//-----------------------------------------------------------------------------



// Class to put a zTorque on specified body of sin(t) and an x force of a 
// specified value.
class PaulsForce : public GeneralForceElements::UserForce {
public:
    explicit PaulsForce( BodyId body, Real howHard ) : b(body), f(howHard) { }

    UserForce* clone() const { 
        return new PaulsForce(*this); 
    }

    void calc(const MatterSubsystem& matter, const State& state,
              Vector_<SpatialVec>& bodyForces,
              Vector_<Vec3>&       particleForces,
              Vector&              mobilityForces,
              Real&                pe) const 
    {
        SpatialVec& bodiesForces = bodyForces[b];
        Vec3& torqueSum = bodiesForces[0];
        Vec3& forceSum  = bodiesForces[1];
        Real time = state.getTime();
        torqueSum[2] += 0*sin(time);
        forceSum[0]  += f;

    }

private:
    BodyId b;
    Real   f;
};


//-----------------------------------------------------------------------------
// Prototypes for local functions (functions not called by code in other files)
//-----------------------------------------------------------------------------
bool  SimulateSpinningBook( void );
bool  WriteStringToFile(   const char outputString[], FILE *fptr )  { return fputs( outputString, fptr ) != 0; }
bool  WriteStringToScreen( const char outputString[] )              { return WriteStringToFile( outputString, stdout ); }
bool  WriteDoubleToFile( double x, int precision, FILE *fptr );
FILE*  FileOpenWithMessageIfCannotOpen( const char *filename, const char *attribute );
double  ConvertFromRadiansToDegrees( double angleInRadians )  { return angleInRadians * 57.295779513082320876798154814105170332405472466564321549160243861;   }
double  ConvertFromDegreesToRadians( double angleInDegrees )  { return angleInDegrees * 0.017453292519943295769236907684886127134428718885417254560971914402; }


//-----------------------------------------------------------------------------
// The executable program starts here
//-----------------------------------------------------------------------------
int  main( int numberOfCommandLineArguments, char *arrayOfCommandLineArguments[] )
{
   // Simulate the multibody system
   bool simulationSucceeded = SimulateSpinningBook();

   // Keep the screen displayed until the user presses the Enter key
   WriteStringToScreen( "\n\n Press  Enter  to terminate the program: " );
   getchar();

   // The value returned by the main function is the exit status of the program.
   // A normal program exit returns 0 (other return values usually signal an error).
   return simulationSucceeded == true ? 0 : 1;
}



//-----------------------------------------------------------------------------
bool  SimulateSpinningBook( void )
{
   // Declare a multibody system (contains one or more force and matter sub-systems)
   MultibodySystem  mbs;

   // 0. The ground's right-handed, orthogonal x,y,z unit vectors are directed with x horizontally right and y vertically upward.
   // 1. Create a gravity vector that is straight down (in the ground's frame)
   // 2. Create a uniform gravity sub-system
   // 3. Add the gravity sub-system to the multibody system
   Vec3 gravityVector( 0, 0, 0 );
   UniformGravitySubsystem gravity( gravityVector );
   mbs.addForceSubsystem( gravity );



   // Create a matter sub-system (the baseball)
   SimbodyMatterSubsystem  book;

   // Physical dimensions of a book (given in centimeters)
   Real xbookWidth     = 20.0/100;
   Real ybookHeight    = 30.0/100;
   Real zbookThickness =  5.0/100;

   // Create the mass, center of mass, and inertia properties for the book
   const Real  massOfBook = 0.4;
   const Vec3  bookCenterOfMassLocation( 0, 0, 0 );
   const Real  Ixx = massOfBook/12*( ybookHeight*ybookHeight + zbookThickness*zbookThickness );
   const Real  Iyy = massOfBook/12*( xbookWidth*xbookWidth   + zbookThickness*zbookThickness );
   const Real  Izz = massOfBook/12*( ybookHeight*ybookHeight + xbookWidth*xbookWidth) ;
   const Real  Ixy = 0,  Ixz = 0,   Iyz = 0;
   const Inertia  bookInertiaMatrix( Ixx, Iyy, Izz, Ixy, Ixz, Iyz );
   MassProperties  bookMassProperties( massOfBook, bookCenterOfMassLocation, bookInertiaMatrix );

   // Define how the book is connected to the ground
   const Rotation  inboardJointOrientationInBaseball;   // ( 1,0,0, 0,1,0, 0,0,1 );
   const Vec3  inboardJointLocationFromBaseballOrigin( 0, 0, 0 );
   const Transform  inboardJointTransformFromBaseball( inboardJointOrientationInBaseball, inboardJointLocationFromBaseballOrigin );

   // The book's motion is related to ground via "mobilizers".
   const Transform outboardFrameTransformFromGround;  // The default constructor is the identity transform
   const Transform  inboardFrameTransformFromBook;   // The default constructor is the identity transform
   Mobilizer bookToGroundMobilizer = Mobilizer::Free(); // (could also just use a ball and socket mobilizer)
   const BodyId bookBodyId = book.addRigidBody( bookMassProperties, inboardFrameTransformFromBook, GroundId, outboardFrameTransformFromGround, bookToGroundMobilizer );

   // Add the matter (book) sub-system to the system.
   mbs.setMatterSubsystem( book );

   GeneralForceElements paulsForces;
   mbs.addForceSubsystem(paulsForces);

   PaulsForce pf(bookBodyId, 27);

   // This adds a copy of the concrete user force; what happens to the original doesn't matter
   // after this call.
   paulsForces.addUserForce(pf);


   // Create a state for this system.
   // Define appropriate states for this multi-body system.
   // Set the initial time to 0.0
   State s;
   mbs.realize( s );
   s.setTime( 0.0 );

   // Set the initial values for the configuration variables (Euler parameters, then x,y,z)
   book.setMobilizerQ( s, bookBodyId, 0,  1.0 );
   book.setMobilizerQ( s, bookBodyId, 1,  0.0 );
   book.setMobilizerQ( s, bookBodyId, 2,  0.0 );
   book.setMobilizerQ( s, bookBodyId, 3,  0.0 );
   book.setMobilizerQ( s, bookBodyId, 4,  0.0 );
   book.setMobilizerQ( s, bookBodyId, 5,  0.0 );
   book.setMobilizerQ( s, bookBodyId, 6,  0.0 );

   // Set the initial values for the motion variables (angular velocity measure numbers, then vx, vy, vz)
   Real  wx = 0.2,  wy = 0.7,  wz = 0.2;
   // Real  wx = 7.0,  wy = 0.2,  wz = 0.2;
   // Real  wx = 0.2,  wy = 0.2,  wz = 7.0;
   book.setMobilizerU( s, bookBodyId, 0,  wx  );
   book.setMobilizerU( s, bookBodyId, 1,  wy  );
   book.setMobilizerU( s, bookBodyId, 2,  wz  );
   book.setMobilizerU( s, bookBodyId, 3,  0.0 );
   book.setMobilizerU( s, bookBodyId, 4,  0.0 );
   book.setMobilizerU( s, bookBodyId, 5,  0.0 );

   // Create a study using the Runge Kutta Merson integrator (alternately use the CPodesIntegrator)
   RungeKuttaMerson myStudy( mbs, s );

   // Set the numerical accuracy for the integrator
   myStudy.setAccuracy( 1.0E-7 );

   // The next statement does lots of accounting
   myStudy.initialize();

   // Open a file to record the simulation results (they are also displayed on screen)
   FILE *outputFile = FileOpenWithMessageIfCannotOpen( "SpinningBookResults.txt", "w" );
   WriteStringToFile(   "  time       wx              wy              wz             Hx             Hy              Hz              HMagnitude  kineticEnergy \n",  outputFile );
   WriteStringToScreen( "  time       wx              wy              wz             Hx             Hy              Hz              HMagnitude  kineticEnergy \n" );

   // Visualize results with VTKReporter
   //VTKReporter animationResults( mbs );

   // For visualization purposes only, create a red sphere with a radius of 0.5 meters (huge but visible)
   Vec3  halfLengthsOfBook = 0.5*Vec3( xbookWidth, ybookHeight, zbookThickness );
   DecorativeBrick redBrick = DecorativeBrick( halfLengthsOfBook );
   redBrick.setColor(Red);     // Can also specify a Vec3 with rgb which scale from 0 to 1
   redBrick.setOpacity(0.0);   // 0.0 is solid and 1.0 is transparent

   // Decorate the book with the red sphere at the book's origin
   const Rotation  bookToRedBrickOrientation; //( 1,0,0, 0,1,0, 0,0,1 );
   const Vec3      bookOriginToRedBrickOriginLocation( 0, 0, 0 );
   const Transform bookToRedBrickTransform( bookToRedBrickOrientation, bookOriginToRedBrickOriginLocation );
   //animationResults.addDecoration( bookBodyId, bookToRedBrickTransform, redBrick );

   // Set the numerical integration step and the time for the simulation to run
   const Real integrationStepDt = 0.01;
   const Real finalTime = 4.0;
   const Real finalTimeCompare = finalTime - 0.01*integrationStepDt;

   // Run the simulation and print the results
   while( 1 )
   {
      // Query for results to be printed
      Real time = s.getTime();
      Real kineticEnergy = mbs.getKineticEnergy(s);

      // Get the angular velocity of the book in ground, expressed in the ground's "x,y,z" axes.
      const Vec3 angularVelocityExpressedInGround = book.calcBodyAngularVelocityInBody( s, bookBodyId, GroundId );

      // Get the rotation matrix relating the body's x,y,z unit vectors to the ground's x,y,z unit vectors.
      Rotation bookToGroundRotationMatrixB = book.calcBodyRotationFromBody( s, bookBodyId, GroundId );

      // Get the rotation matrix from ground frame to book frame
      const Rotation groundToBookRotationMatrix = book.getBodyRotation( s, bookBodyId );

      // Create the inverse rotation matrix of book to ground
      const InverseRotation bookToGroundRotationMatrix = groundToBookRotationMatrix.invert();

      // Reexpress the angular velocity in body bases
      const Vec3 angularVelocityExpressedInBook = bookToGroundRotationMatrix * angularVelocityExpressedInGround;

      const Mat33 check = bookToGroundRotationMatrix.toMat33() - bookToGroundRotationMatrixB.toMat33();

      // Get the components of the angular velocity for plotting
      Real wx = angularVelocityExpressedInBook[0];
      Real wy = angularVelocityExpressedInBook[1];
      Real wz = angularVelocityExpressedInBook[2];

      // Form the angular momentum measure numbers and the magnitude of the angular momentum
      Real Hx = Ixx*wx;
      Real Hy = Iyy*wy;
      Real Hz = Izz*wz;
      Real Hmag = sqrt( Hx*Hx + Hy*Hy + Hz*Hz );

      // Print results to screen
      WriteDoubleToFile( time,                 2, stdout );
      WriteDoubleToFile( wx,                   7, stdout );
      WriteDoubleToFile( wy,                   7, stdout );
      WriteDoubleToFile( wz,                   7, stdout );
      WriteDoubleToFile( Hx,                   7, stdout );
      WriteDoubleToFile( Hy,                   7, stdout );
      WriteDoubleToFile( Hz,                   7, stdout );
      WriteDoubleToFile( Hmag,                 7, stdout );
      WriteDoubleToFile( kineticEnergy,        7, stdout );
      WriteStringToScreen( "\n" );
      // Print results to file
      WriteDoubleToFile( time,                 2, outputFile );
      WriteDoubleToFile( wx,                   7, outputFile );
      WriteDoubleToFile( wy,                   7, outputFile );
      WriteDoubleToFile( wz,                   7, outputFile );
      WriteDoubleToFile( Hx,                   7, outputFile );
      WriteDoubleToFile( Hy,                   7, outputFile );
      WriteDoubleToFile( Hz,                   7, outputFile );
      WriteDoubleToFile( Hmag,                 7, outputFile );
      WriteDoubleToFile( kineticEnergy,        7, outputFile );
      WriteStringToFile( "\n", outputFile );

      // Animate the results for this step
      //animationResults.report(s);

      // Check if integration has completed
      if( time >= finalTimeCompare ) break;

      // Increment time step
      myStudy.step( time + integrationStepDt );
   }

   // Simulation completed properly
   return true;
}


//-----------------------------------------------------------------------------
FILE*  FileOpenWithMessageIfCannotOpen( const char *filename, const char *attribute )
{
   // Try to open the file
   FILE *Fptr1 = fopen( filename, attribute );

   // If unable to open the file, issue a message
   if( !Fptr1 )
   {
      WriteStringToScreen( "\n\n Unable to open the file: " );
      WriteStringToScreen( filename );
      WriteStringToScreen( "\n\n" );
   }

   return Fptr1;
}


//-----------------------------------------------------------------------------
bool  WriteDoubleToFile( double x, int precision, FILE *fptr )
{
   // Ensure the precision (number of digits in the mantissa after the decimal point) makes sense.
   // Next, calculate the field width so it includes one extra space to the right of the number.
   if( precision < 0 || precision > 17 ) precision = 5;
   int fieldWidth = precision + 8;

   // Create the format specifier and print the number
   char format[20];
   sprintf( format, " %%- %d.%dE", fieldWidth, precision );
   return fprintf( fptr, format, x ) >= 0;
}





