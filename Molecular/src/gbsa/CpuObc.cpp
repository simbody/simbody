
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string.h>
#include <sstream>

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "CpuObc.h"
#include <math.h>

/**---------------------------------------------------------------------------------------

   CpuObc constructor

   obcParameters      obcParameters object
   
   --------------------------------------------------------------------------------------- */

CpuObc::CpuObc( ImplicitSolventParameters* obcParameters ) : CpuImplicitSolvent( obcParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::CpuObc";

   // ---------------------------------------------------------------------------------------

   _initializeObcDataMembers( );

   _obcParameters = static_cast<ObcParameters*> (obcParameters);

}

/**---------------------------------------------------------------------------------------

   CpuObc destructor

   --------------------------------------------------------------------------------------- */

CpuObc::~CpuObc( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::~CpuObc";

   // ---------------------------------------------------------------------------------------

   //if( _obcParameters != NULL ){
     // delete _obcParameters;
   //}

   delete[] _obcChain;
   delete[] _obcChainTemp;
}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuObc::_initializeObcDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _obcParameters = NULL;
   _obcChain      = NULL;
   _obcChainTemp  = NULL;
}

/**---------------------------------------------------------------------------------------

   Get ObcParameters reference

   @return ObcParameters reference

   --------------------------------------------------------------------------------------- */

ObcParameters* CpuObc::getObcParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getObcParameters";

   // ---------------------------------------------------------------------------------------

   return _obcParameters;
}

/**---------------------------------------------------------------------------------------

   Set ObcParameters reference

   @param ObcParameters reference

   @return SimTKOpenMMCommon::DefaultReturn;

   --------------------------------------------------------------------------------------- */

int CpuObc::setObcParameters(  ObcParameters* obcParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setObcParameters";

   // ---------------------------------------------------------------------------------------

   _obcParameters = obcParameters;
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

Real* CpuObc::getObcChain( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getObcChain";

   // ---------------------------------------------------------------------------------------

   if( _obcChain == NULL ){
      _obcChain = new Real[_obcParameters->getNumberOfAtoms()];
   }
   return _obcChain;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcParameters->getNumberOfAtoms()

   @return array

   --------------------------------------------------------------------------------------- */

Real* CpuObc::getObcChainConst( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getObcChain";

   // ---------------------------------------------------------------------------------------

   return _obcChain;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain temp work array of size=_obcParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

Real* CpuObc::getObcChainTemp( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getImplicitSolventObcChainTemp";

   // ---------------------------------------------------------------------------------------

   if( _obcChainTemp == NULL ){
      _obcChainTemp = new Real[_obcParameters->getNumberOfAtoms()];
   }
   return _obcChainTemp;
}

/**---------------------------------------------------------------------------------------

   Get Born radii based on papers:

      J. Phys. Chem. 1996 100, 19824-19839 (HCT paper)
      Proteins: Structure, Function, and Bioinformatcis 55:383-394 (2004) (OBC paper)

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii

   @return array of Born radii

   --------------------------------------------------------------------------------------- */

int CpuObc::computeBornRadii( Real** atomCoordinates, Real* bornRadii, Real* obcChain ){

   // ---------------------------------------------------------------------------------------

   static const Real zero    = (Real) 0.0;
   static const Real one     = (Real) 1.0;
   static const Real two     = (Real) 2.0;
   static const Real three   = (Real) 3.0;
   static const Real half    = (Real) 0.5;
   static const Real fourth  = (Real) 0.25;

   // static const char* methodName = "\nCpuObc::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   ObcParameters* obcParameters    = getObcParameters();

   int numberOfAtoms               = obcParameters->getNumberOfAtoms();
   Real* atomicRadii               = obcParameters->getAtomicRadii();
   const Real* scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();
   if( !obcChain ){
      obcChain                     = getObcChain();
   }

   Real dielectricOffset           = obcParameters->getDielectricOffset();
   Real alphaObc                   = obcParameters->getAlphaObc();
   Real betaObc                    = obcParameters->getBetaObc();
   Real gammaObc                   = obcParameters->getGammaObc();

   // ---------------------------------------------------------------------------------------

   // calculate Born radii

//FILE* logFile = SimTKOpenMMLog::getSimTKOpenMMLogFile( );

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
     
      Real radiusI         = atomicRadii[atomI];
      Real offsetRadiusI   = radiusI - dielectricOffset;

      Real radiusIInverse  = one/offsetRadiusI;
      Real sum             = zero;

      // HCT code

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            Real deltaX          = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            Real deltaY          = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            Real deltaZ          = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
            Real r2              = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            Real r               = SQRT( r2 );
            Real offsetRadiusJ   = atomicRadii[atomJ] - dielectricOffset; 
            Real scaledRadiusJ   = offsetRadiusJ*scaledRadiusFactor[atomJ];
            Real rScaledRadiusJ  = r + scaledRadiusJ;

            if( offsetRadiusI < rScaledRadiusJ ){
               Real rInverse = one/r;
               Real l_ij     = offsetRadiusI > FABS( r - scaledRadiusJ ) ? offsetRadiusI : FABS( r - scaledRadiusJ );
                    l_ij     = one/l_ij;

               Real u_ij     = one/rScaledRadiusJ;

               Real l_ij2    = l_ij*l_ij;
               Real u_ij2    = u_ij*u_ij;
 
               Real ratio    = LN( (u_ij/l_ij) );
               Real term     = l_ij - u_ij + fourth*r*(u_ij2 - l_ij2)  + ( half*rInverse*ratio) + (fourth*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2 - u_ij2);
               if( offsetRadiusI < (scaledRadiusJ - r) ){
                  term += two*( radiusIInverse - l_ij);
               }
               sum += term;

/*
if( atomI == 0 ){
   (void) fprintf( logFile, "\nRR %d %d r=%.4f rads[%.6f %.6f] scl=[%.3f %.3f] sum=%12.6e %12.6e %12.6e %12.6e",
                   atomI, atomJ, r, offsetRadiusI, offsetRadiusJ, scaledRadiusFactor[atomI], scaledRadiusFactor[atomJ], 0.5f*sum,
                   l_ij, u_ij, term );
} */

            }
         }
      }
 
      // OBC-specific code (Eqs. 6-8 in paper)

      sum                  *= half*offsetRadiusI;
      Real sum2             = sum*sum;
      Real sum3             = sum*sum2;
      Real tanhSum          = TANH( alphaObc*sum - betaObc*sum2 + gammaObc*sum3 );
      
      bornRadii[atomI]      = one/( one/offsetRadiusI - tanhSum/radiusI ); 
 
      obcChain[atomI]       = offsetRadiusI*( alphaObc - two*betaObc*sum + three*gammaObc*sum2 );
      obcChain[atomI]       = (one - tanhSum*tanhSum)*obcChain[atomI]/radiusI;

/*
if( atomI == 0 ){
   (void) fprintf( logFile, "\nRRQ %d sum=%12.6e tanhS=%12.6e radI=%.5f %.5f born=%12.6e obc=%12.6e",
                   atomI, sum, tanhSum, radiusI, offsetRadiusI, bornRadii[atomI], obcChain[atomI] );
} */

   }

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get Obc Born energy and forces

   @param bornRadii           Born radii -- optional; if NULL, then ObcParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   @return SimTKOpenMMCommon::DefaultReturn;

   The array bornRadii is also updated and the obcEnergy

   --------------------------------------------------------------------------------------- */

int CpuObc::computeBornEnergyForces( Real* bornRadii, Real** atomCoordinates,
                                     const Real* partialCharges, Real** forces ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::computeBornEnergyForces";

   static const Real zero    = (Real) 0.0;
   static const Real one     = (Real) 1.0;
   static const Real two     = (Real) 2.0;
   static const Real three   = (Real) 3.0;
   static const Real four    = (Real) 4.0;
   static const Real half    = (Real) 0.5;
   static const Real fourth  = (Real) 0.25;
   static const Real eighth  = (Real) 0.125;

   // ---------------------------------------------------------------------------------------

   const ObcParameters* obcParameters = getObcParameters();
   const int numberOfAtoms            = obcParameters->getNumberOfAtoms();

   if( bornRadii == NULL ){
      bornRadii   = getBornRadii();
   }

   // ---------------------------------------------------------------------------------------

   // constants

   const Real preFactor           = obcParameters->getPreFactor();
   const Real electricConstant    = obcParameters->getElectricConstant();
   const Real dielectricOffset    = obcParameters->getDielectricOffset();

   // ---------------------------------------------------------------------------------------

   // set energy/forces to zero

   Real obcEnergy                       = zero;
   const unsigned int arraySzInBytes    = sizeof( Real )*numberOfAtoms;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      memset( forces[ii], 0, 3*sizeof( Real ) );
   }

   Real* bornForces = getBornForce();
   memset( bornForces, 0, arraySzInBytes );

   // ---------------------------------------------------------------------------------------

   // N*( 8 + pow) ACE
   // compute the nonpolar solvation via ACE approximation
    
   if( includeAceApproximation() ){
      computeAceNonPolarForce( obcParameters, bornRadii, &obcEnergy, bornForces );
   }

   // ---------------------------------------------------------------------------------------

   // first main loop

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      Real partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         // 3 FLOP

         Real deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
         Real deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
         Real deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
         // 5 FLOP

         Real r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;

         // 3 FLOP

         Real alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
         Real D_ij               = r2/(four*alpha2_ij);

         // exp + 2 + sqrt FLOP 

         Real expTerm            = EXP( -D_ij );
         Real denominator2       = r2 + alpha2_ij*expTerm; 
         Real denominator        = SQRT( denominator2 ); 
         
         // 6 FLOP

         Real Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
         Real dGpol_dr           = -Gpol*( one - fourth*expTerm )/denominator2;  

         // 5 FLOP

         Real dGpol_dalpha2_ij   = -half*Gpol*expTerm*( one + D_ij )/denominator2;

         // 11 FLOP

         if( atomI != atomJ ){

             bornForces[atomJ] += dGpol_dalpha2_ij*bornRadii[atomI];

             deltaX            *= dGpol_dr;
             deltaY            *= dGpol_dr;
             deltaZ            *= dGpol_dr;

             forces[atomI][0]  -= deltaX;
             forces[atomI][1]  -= deltaY;
             forces[atomI][2]  -= deltaZ;

             forces[atomJ][0]  += deltaX;
             forces[atomJ][1]  += deltaY;
             forces[atomJ][2]  += deltaZ;

         } else {
            Gpol *= half;
         }

         // 3 FLOP

         obcEnergy         += Gpol;
         bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

      }
   }

   //obcEnergy *= getEnergyConversionFactor();

   // ---------------------------------------------------------------------------------------

   // second main loop

   // initialize Born radii & ObcChain temp arrays -- contain values
   // used in next iteration

   Real* bornRadiiTemp             = getBornRadiiTemp();
   memset( bornRadiiTemp, 0, arraySzInBytes );

   Real* obcChainTemp              = getObcChainTemp();
   memset( obcChainTemp, 0, arraySzInBytes );

   Real* obcChain                  = getObcChain();
   const Real* atomicRadii         = obcParameters->getAtomicRadii();

   const Real alphaObc             = obcParameters->getAlphaObc();
   const Real betaObc              = obcParameters->getBetaObc();
   const Real gammaObc             = obcParameters->getGammaObc();
   const Real* scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();

    // compute factor that depends only on the outer loop index

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      bornForces[atomI] *= bornRadii[atomI]*bornRadii[atomI]*obcChain[atomI];      
   }

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // radius w/ dielectric offset applied

      Real radiusI        = atomicRadii[atomI];
      Real offsetRadiusI  = radiusI - dielectricOffset;

      // used to compute Born radius for next iteration

      Real bornSum        = zero;

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            Real deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            Real deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            Real deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
    
            Real r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            Real r                  = SQRT( r2 );
 
            // radius w/ dielectric offset applied

            Real offsetRadiusJ      = atomicRadii[atomJ] - dielectricOffset;

            Real scaledRadiusJ      = offsetRadiusJ*scaledRadiusFactor[atomJ];
            Real scaledRadiusJ2     = scaledRadiusJ*scaledRadiusJ;
            Real rScaledRadiusJ     = r + scaledRadiusJ;

            // dL/dr & dU/dr are zero (this can be shown analytically)
            // removed from calculation

            if( offsetRadiusI < rScaledRadiusJ ){

               Real l_ij          = offsetRadiusI > FABS( r - scaledRadiusJ ) ? offsetRadiusI : FABS( r - scaledRadiusJ );
                    l_ij          = one/l_ij;

               Real u_ij          = one/rScaledRadiusJ;

               Real l_ij2         = l_ij*l_ij;

               Real u_ij2         = u_ij*u_ij;
 
               Real rInverse      = one/r;
               Real r2Inverse     = one/r2;

               Real t3            = eighth*(one + scaledRadiusJ2*r2Inverse)*(l_ij2 - u_ij2) + fourth*LN( u_ij/l_ij )*r2Inverse;

               Real de            = bornForces[atomI]*t3*rInverse;

               deltaX            *= de;
               deltaY            *= de;
               deltaZ            *= de;
   
               forces[atomI][0]  += deltaX;
               forces[atomI][1]  += deltaY;
               forces[atomI][2]  += deltaZ;
  
               forces[atomJ][0]  -= deltaX;
               forces[atomJ][1]  -= deltaY;
               forces[atomJ][2]  -= deltaZ;
 
               // Born radius term

               Real term          =   l_ij - u_ij  + fourth*r*(u_ij2 - l_ij2) + (half*rInverse)*LN(u_ij/l_ij)   +
                                      (fourth*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2-u_ij2);

               if( offsetRadiusI < (scaledRadiusJ - r) ){
                  term += two*( (one/offsetRadiusI) - l_ij);
               }
               bornSum += term; 
            }
         }
      }

      // OBC-specific code (Eqs. 6-8 in paper)

      bornSum             *= half*offsetRadiusI;
      Real sum2            = bornSum*bornSum;
      Real sum3            = bornSum*sum2;
      Real tanhSum         = TANH( alphaObc*bornSum - betaObc*sum2 + gammaObc*sum3 );
      
      bornRadiiTemp[atomI] = one/( one/offsetRadiusI - tanhSum/radiusI ); 
 
      obcChainTemp[atomI]  = offsetRadiusI*( alphaObc - two*betaObc*bornSum + three*gammaObc*sum2 );
      obcChainTemp[atomI]  = (one - tanhSum*tanhSum)*obcChainTemp[atomI]/radiusI;
   }

   setEnergy( obcEnergy );

   // copy new Born radii and obcChain values into permanent array

   memcpy( bornRadii, bornRadiiTemp, arraySzInBytes );
   memcpy( obcChain, obcChainTemp, arraySzInBytes );

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state 
   
   @param title               title (optional)
      
   @return string containing state
      
   --------------------------------------------------------------------------------------- */

std::string CpuObc::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << CpuImplicitSolvent::getStateString( title );

   return message.str();
}

/**---------------------------------------------------------------------------------------

   Write Born energy and forces (Simbios)

   @param atomCoordinates     atomic coordinates
   @param Real forces         forces
   @param resultsFileName     output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuObc::writeBornEnergyForces( Real** atomCoordinates,
                                   const Real* partialCharges, Real** forces,
                                   const std::string& resultsFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nCpuObc::writeBornEnergyForces";

   // ---------------------------------------------------------------------------------------

   ImplicitSolventParameters* implicitSolventParameters = getImplicitSolventParameters();
   const ObcParameters* obcParameters                   = static_cast<const ObcParameters*>(implicitSolventParameters);
   

   int numberOfAtoms              = obcParameters->getNumberOfAtoms();
   const Real* atomicRadii        = obcParameters->getAtomicRadii();
   const Real* bornRadii          = getBornRadiiConst();
   const Real* scaledRadii        = obcParameters->getScaledRadiusFactors();
   const Real* obcChain           = getObcChainConst();

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* implicitSolventResultsFile = NULL;
#ifdef WIN32
   fopen_s( &implicitSolventResultsFile, resultsFileName.c_str(), "w" );
#else
   implicitSolventResultsFile = fopen( resultsFileName.c_str(), "w" );
#endif

   // diganostics

   std::stringstream message;
   message << methodName;
   if( implicitSolventResultsFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << resultsFileName << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << "  could not open file=<" << resultsFileName << "> -- abort output.";
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;
   }

   // header

   (void) fprintf( implicitSolventResultsFile, "# %d atoms format: coords(3) bornRadii(input) q atomicRadii scaleFactors forces obcChain\n", numberOfAtoms );

   Real forceConversion  = 1.0f;
   Real lengthConversion = 1.0f;

   // output

   if( forces != NULL && atomCoordinates != NULL && partialCharges != NULL && atomicRadii != NULL ){
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
            (void) fprintf( implicitSolventResultsFile, "%.7e %.7e %.7e %.7e %.5f %.5f %.5f %.7e %.7e %.7e %.7e\n",
                            lengthConversion*atomCoordinates[ii][0],
                            lengthConversion*atomCoordinates[ii][1], 
                            lengthConversion*atomCoordinates[ii][2],
                           (bornRadii != NULL ? lengthConversion*bornRadii[ii] : 0.0f),
                            partialCharges[ii], lengthConversion*atomicRadii[ii], scaledRadii[ii],
                            forceConversion*forces[ii][0],
                            forceConversion*forces[ii][1],
                            forceConversion*forces[ii][2],
                            forceConversion*obcChain[ii]
                          );
      }
   }
   (void) fclose( implicitSolventResultsFile );

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Write  results from first loop

   @param atomCoordinates     atomic coordinates
   @param Real forces         forces
   @param outputFileName      output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuObc::writeForceLoop1( int numberOfAtoms, Real** forces, const Real* bornForce,
                             const std::string& outputFileName ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nCpuObc::writeForceLoop1";

   // ---------------------------------------------------------------------------------------

   int chunkSize;
   if( bornForce ){
      chunkSize = 3;
   } else {
      chunkSize = 4;
   }

   StringVector lineVector;
   std::stringstream header;
   lineVector.push_back( "# bornF F" );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      std::stringstream line;
      line << (atomI+1) << " ";
      SimTKOpenMMUtilities::formatRealStringStream( line, forces[atomI], chunkSize );
      if( bornForce ){
         line << " " << bornForce[atomI];
      }
      lineVector.push_back( line.str() );
   }
   return SimTKOpenMMUtilities::writeFile( lineVector, outputFileName );

}

/**---------------------------------------------------------------------------------------

   Write results

   @param numberOfAtoms       number of atoms
   @param chunkSizes          vector of chunk sizes for realRealVector
   @param realRealVector      vector of Real**
   @param realVector          vector of Real*
   @param outputFileName      output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuObc::writeForceLoop( int numberOfAtoms, const IntVector& chunkSizes,
                            const RealPtrPtrVector& realRealVector, 
                            const RealPtrVector& realVector,
                            const std::string& outputFileName ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nCpuObc::writeForceLoop";

   static const int maxChunks = 10;
   int chunks[maxChunks];

   // ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < (int) chunkSizes.size(); ii++ ){
      chunks[ii] = chunkSizes[ii];
   }
   for( int ii = (int) chunkSizes.size(); ii < maxChunks; ii++ ){
      chunks[ii] = 3;
   }

   StringVector lineVector;
   std::stringstream header;
   // lineVector.push_back( "# " );

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      std::stringstream line;
      line << (atomI+1) << " ";

      int index = 0;
      for( RealPtrPtrVectorCI ii = realRealVector.begin(); ii != realRealVector.end(); ii++ ){
         Real** forces = *ii;
         SimTKOpenMMUtilities::formatRealStringStream( line, forces[atomI], chunks[index++] );
         line << " ";
      }

      for( RealPtrVectorCI ii = realVector.begin(); ii != realVector.end(); ii++ ){
         Real* array = *ii;
         line << array[atomI] << " ";
      }

      lineVector.push_back( line.str() );
   }
   return SimTKOpenMMUtilities::writeFile( lineVector, outputFileName );

}

/**---------------------------------------------------------------------------------------

   Get Obc Born energy and forces -- used debugging

   @param bornRadii           Born radii -- optional; if NULL, then ObcParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   @return SimTKOpenMMCommon::DefaultReturn;

   The array bornRadii is also updated and the obcEnergy

   --------------------------------------------------------------------------------------- */

int CpuObc::computeBornEnergyForcesPrint( Real* bornRadii, Real** atomCoordinates,
                                         const Real* partialCharges, Real** forces ){
 
   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::computeBornEnergyForcesPrint";

   static const Real zero    = (Real) 0.0;
   static const Real one     = (Real) 1.0;
   static const Real two     = (Real) 2.0;
   static const Real three   = (Real) 3.0;
   static const Real four    = (Real) 4.0;
   static const Real half    = (Real) 0.5;
   static const Real fourth  = (Real) 0.25;
   static const Real eighth  = (Real) 0.125;

   // ---------------------------------------------------------------------------------------

   const ObcParameters* obcParameters = getObcParameters();
   const int numberOfAtoms            = obcParameters->getNumberOfAtoms();

   if( bornRadii == NULL ){
      bornRadii   = getBornRadii();
   }

FILE* logFile = SimTKOpenMMLog::getSimTKOpenMMLogFile( );

   // ---------------------------------------------------------------------------------------

   // constants

   const Real preFactor           = obcParameters->getPreFactor();
   const Real electricConstant    = obcParameters->getElectricConstant();
   const Real dielectricOffset    = obcParameters->getDielectricOffset();

   // ---------------------------------------------------------------------------------------

   // set energy/forces to zero

   Real obcEnergy                       = zero;
   const unsigned int arraySzInBytes    = sizeof( Real )*numberOfAtoms;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      memset( forces[ii], 0, 3*sizeof( Real ) );
   }

   Real* bornForces = getBornForce();
   memset( bornForces, 0, arraySzInBytes );

   // ---------------------------------------------------------------------------------------

   // N*( 8 + pow) ACE
   // compute the nonpolar solvation via ACE approximation
    
   if( includeAceApproximation() ){
      computeAceNonPolarForce( obcParameters, bornRadii, &obcEnergy, bornForces );
   }

   // ---------------------------------------------------------------------------------------

   // first main loop

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      Real partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         // 3 FLOP

         Real deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
         Real deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
         Real deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
         // 5 FLOP

         Real r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;

         // 3 FLOP

         Real alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
         Real D_ij               = r2/(four*alpha2_ij);

         // exp + 2 + sqrt FLOP 

         Real expTerm            = EXP( -D_ij );
         Real denominator2       = r2 + alpha2_ij*expTerm; 
         Real denominator        = SQRT( denominator2 ); 
         
         // 6 FLOP

         Real Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
  
         // dGpol/dr              = -1/2*(Gpol/denominator2)*(2r - r/2*exp() )
         Real dGpol_dr           = -Gpol*( one - fourth*expTerm )/denominator2;  

         // 5 FLOP

         Real dGpol_dalpha2_ij   = -half*Gpol*expTerm*( one + D_ij )/denominator2;

         // 11 FLOP

         if( atomI != atomJ ){

             bornForces[atomJ] += dGpol_dalpha2_ij*bornRadii[atomI];

             deltaX            *= dGpol_dr;
             deltaY            *= dGpol_dr;
             deltaZ            *= dGpol_dr;

             forces[atomI][0]  -= deltaX;
             forces[atomI][1]  -= deltaY;
             forces[atomI][2]  -= deltaZ;

             forces[atomJ][0]  += deltaX;
             forces[atomJ][1]  += deltaY;
             forces[atomJ][2]  += deltaZ;

         } else {
            Gpol *= half;
         }

         // 3 FLOP

         obcEnergy         += Gpol;
         bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

if( atomI == -1 || atomJ == -1 ){
//   (void) fprintf( logFile, "\nWWX %d %d F[%.6e %.6e %.6e] bF=[%.6e %.6e] Gpl[%.6e %.6e %.6e] rb[%6.4f %7.4f] rs[%6.4f %7.4f] ",
//                    atomI, atomJ,
//                    forces[atomI][0],  forces[atomI][1],  forces[atomI][2],
//                    bornForces[atomI], bornForces[atomJ],
//                    Gpol,dGpol_dr,dGpol_dalpha2_ij,
//                    bornRadii[atomI],bornRadii[atomJ],atomicRadii[atomI],atomicRadii[atomJ] );
//
   (void) fprintf( logFile, "\nWWX %d %d %.1f r2=%.4f q=%.2f bF=[%.6e %.6e] Gpl[%.6e %.6e %.6e] rb[%.5f %.5f] add[%.6e %.6e] ",
                    atomI, atomJ, preFactor, r2, partialCharges[atomJ],
                    bornForces[atomI], bornForces[atomJ],
                    Gpol,dGpol_dr,dGpol_dalpha2_ij,
                    bornRadii[atomI], bornRadii[atomJ],
                    dGpol_dalpha2_ij*bornRadii[atomJ], dGpol_dalpha2_ij*bornRadii[atomI] );
}
      }


   }

(void) fprintf( logFile, "\nWXX bF & F" );
for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
   (void) fprintf( logFile, "\nWXX %d %.6e q=%.3f F[%.6e %.6e %.6e] ",
                   atomI, partialCharges[atomI],  bornForces[atomI], forces[atomI][0],  forces[atomI][1],  forces[atomI][2] );
}

   obcEnergy *= getEnergyConversionFactor();

   int fileDebug = 1;
   if( fileDebug ){
      std::string outputFileName = "Loop1Cpu.txt";
      CpuObc::writeForceLoop1( numberOfAtoms, forces, bornForces, outputFileName );
/*
      for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
         forces[atomI][0] = forces[atomI][1] = forces[atomI][2] = (Real) 0.0;
      }
*/
   }

   // ---------------------------------------------------------------------------------------

   // second main loop

   // initialize Born radii & ObcChain temp arrays -- contain values
   // used in next iteration

   Real* bornRadiiTemp             = getBornRadiiTemp();
   memset( bornRadiiTemp, 0, arraySzInBytes );

   Real* obcChainTemp              = getObcChainTemp();
   memset( obcChainTemp, 0, arraySzInBytes );

   Real* obcChain                  = getObcChain();
   const Real* atomicRadii         = obcParameters->getAtomicRadii();

   const Real alphaObc             = obcParameters->getAlphaObc();
   const Real betaObc              = obcParameters->getBetaObc();
   const Real gammaObc             = obcParameters->getGammaObc();
   const Real* scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();

    // compute factor that depends only on the outer loop index

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      bornForces[atomI] *= bornRadii[atomI]*bornRadii[atomI]*obcChain[atomI];      
   }

Real* bornSumArray = (Real*) malloc( sizeof( Real )*numberOfAtoms );
memset( bornSumArray, 0, sizeof( Real )*numberOfAtoms );

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // radius w/ dielectric offset applied

      Real radiusI        = atomicRadii[atomI];
      Real offsetRadiusI  = radiusI - dielectricOffset;

      // used to compute Born radius for next iteration

      Real bornSum        = zero;

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            Real deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            Real deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            Real deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
    
            Real r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            Real r                  = SQRT( r2 );
 
            // radius w/ dielectric offset applied

            Real radiusJ            = atomicRadii[atomJ] - dielectricOffset;

            Real scaledRadiusJ      = radiusJ*scaledRadiusFactor[atomJ];
            Real scaledRadiusJ2     = scaledRadiusJ*scaledRadiusJ;
            Real rScaledRadiusJ     = r + scaledRadiusJ;

            // L_ij != 1 && U_ij != 1

            // dL/dr & dU/dr are zero (this can be shown analytically)
            // removed from calculation

            if( offsetRadiusI < rScaledRadiusJ ){

               Real l_ij           = offsetRadiusI > FABS( r - scaledRadiusJ ) ? offsetRadiusI : FABS( r - scaledRadiusJ );
                    l_ij           = one/l_ij;

               Real l_ij2          = l_ij*l_ij;

               Real u_ij           = one/rScaledRadiusJ;
               Real u_ij2          = u_ij*u_ij;
 
               Real rInverse       = one/r;
               Real r2Inverse      = one/r2;

               Real logRatio       = LN( u_ij/l_ij );
               Real t3             = eighth*(one + scaledRadiusJ2*r2Inverse)*(l_ij2 - u_ij2) + fourth*logRatio*r2Inverse;

               Real de             = bornForces[atomI]*t3*rInverse;

               deltaX            *= de;
               deltaY            *= de;
               deltaZ            *= de;

               forces[atomI][0]  += deltaX;
               forces[atomI][1]  += deltaY;
               forces[atomI][2]  += deltaZ;
  
               forces[atomJ][0]  -= deltaX;
               forces[atomJ][1]  -= deltaY;
               forces[atomJ][2]  -= deltaZ;
 
               // Born radius term

               Real term          =   l_ij - u_ij + fourth*r*(u_ij2 - l_ij2) + (half*rInverse)*logRatio + (fourth*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2-u_ij2);

               if( offsetRadiusI < (scaledRadiusJ - r) ){
                  term += two*( (one/offsetRadiusI) - l_ij);
               }
               bornSum += term; 

if( atomI == -1 || atomJ == -1 ){
   (void) fprintf( logFile, "\nXXY %d %d de=%.6e bo[%.6e %.6e %6e] dl %.0f t3=%.6e r=%.6e f[%.6e %.6e %.6e] f[%.6e %.6e %.6e]",
                    atomI, atomJ, de,
                    bornForces[atomI], obcChain[atomI],
                    t3, r, forces[atomI][0],  forces[atomI][1],  forces[atomI][2],
                    deltaX, deltaY, deltaZ );
}
            }
        }
      }

bornSumArray[atomI] = half*bornSum;

      // OBC-specific code (Eqs. 6-8 in paper)

      bornSum             *= half*offsetRadiusI;
      Real sum2            = bornSum*bornSum;
      Real sum3            = bornSum*sum2;
      Real tanhSum         = TANH( alphaObc*bornSum - betaObc*sum2 + gammaObc*sum3 );
      
      bornRadiiTemp[atomI] = one/( one/offsetRadiusI - tanhSum/radiusI); 
 
      obcChainTemp[atomI]  = offsetRadiusI*( alphaObc - two*betaObc*bornSum + three*gammaObc*sum2 );
      obcChainTemp[atomI]  = (one - tanhSum*tanhSum)*obcChainTemp[atomI]/radiusI;

   }

   obcEnergy *= getEnergyConversionFactor();
   setEnergy( obcEnergy );

   if( fileDebug ){

      std::string outputFileName = "Loop2Cpu.txt";

      IntVector chunkVector;
      chunkVector.push_back( 3 );

      RealPtrPtrVector realPtrPtrVector;
      realPtrPtrVector.push_back( forces );

      RealPtrVector realPtrVector;
realPtrVector.push_back( bornSumArray );
      realPtrVector.push_back( bornRadiiTemp );
      realPtrVector.push_back( obcChainTemp );

      CpuObc::writeForceLoop( numberOfAtoms, chunkVector, realPtrPtrVector, realPtrVector, outputFileName );
free( (char*) bornSumArray );
   }

   // 6 FLOP

/*
   Real forceFactor    = getForceConversionFactor();
   Real constantFactor = 1.0f/electricConstant;
   if( fabs(forceFactor - 1.0f) > 1.0e-04 ){
      constantFactor *= forceFactor;
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         forces[ii][0]  *= forceFactor;
         forces[ii][1]  *= forceFactor;
         forces[ii][2]  *= forceFactor;
      }
   } */

   // copy new Born radii and obcChain values into permanent array

//(void) fprintf( logFile, "\nBorn radii not being updated!!!!" );
   memcpy( bornRadii, bornRadiiTemp, arraySzInBytes );
   memcpy( obcChain, obcChainTemp, arraySzInBytes );

   return SimTKOpenMMCommon::DefaultReturn;

}
