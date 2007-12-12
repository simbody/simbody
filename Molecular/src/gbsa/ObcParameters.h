
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

#ifndef __ObcParameters_H__
#define __ObcParameters_H__

#include "SimTKOpenMMCommon.h"
#include "ImplicitSolventParameters.h"

// ---------------------------------------------------------------------------------------

class ObcParameters : public ImplicitSolventParameters {

   public:

       // OBC types

       enum ObcType { ObcTypeI, ObcTypeII };

       static const std::string ParameterFileName;

   private:

      // OBC constants & parameters
   
      Real _dielectricOffset;
      Real _alphaObc;
      Real _betaObc;
      Real _gammaObc;
      ObcType _obcType;

      // scaled radius factors (S_kk in HCT paper)

      int _ownScaledRadiusFactors;
      Real* _scaledRadiusFactors;

      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric (Simbios) 
      
         @param dielectricOffset         solvent dielectric
      
         @return SimTKOpenMMCommon::DefaultReturn

         --------------------------------------------------------------------------------------- */
      
      int setDielectricOffset( Real dielectricOffset );

   public:

      /**---------------------------------------------------------------------------------------
      
         ObcParameters constructor (Simbios) 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       ObcParameters( int numberOfAtoms, ObcParameters::ObcType obcType = ObcTypeII );

      /**---------------------------------------------------------------------------------------
      
         ObcParameters destructor (Simbios) 
      
         --------------------------------------------------------------------------------------- */
      
       ~ObcParameters( );

      // override of new/delete

      //static void* operator new( size_t size );
      //static void  operator delete( void *p );

      //static void* operator new[]( size_t size );
      //static void  operator delete[]( void *p );

      /**---------------------------------------------------------------------------------------
      
         Get OBC type
      
         @return OBC type
      
         --------------------------------------------------------------------------------------- */
      
      ObcParameters::ObcType getObcType( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set OBC type specific parameters
      
         @param obcType OBC type (ObcTypeI or ObcTypeII -- Eq. 7 or 8)
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int setObcTypeParameters( ObcParameters::ObcType obcType );
      
      /**---------------------------------------------------------------------------------------
      
         Get alpha OBC (Eqs. 6 & 7) in Proteins paper
      
         @return alphaObc
      
         --------------------------------------------------------------------------------------- */
      
      Real getAlphaObc( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get beta OBC (Eqs. 6 & 7) in Proteins paper
      
         @return betaObc
      
         --------------------------------------------------------------------------------------- */
      
      Real getBetaObc( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get gamma OBC (Eqs. 6 & 7) in Proteins paper
      
         @return gammaObc
      
         --------------------------------------------------------------------------------------- */
      
      Real getGammaObc( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get solvent dielectric (Simbios) 
      
         @return dielectricOffset dielectric offset
      
         --------------------------------------------------------------------------------------- */
      
      Real getDielectricOffset( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Return OBC scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const Real* getScaledRadiusFactors( void ) const;
        
      /**---------------------------------------------------------------------------------------
      
         Return OBC scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      int setScaledRadiusFactors( Real* scaledRadiusFactors );
#if RealType == 2
      int setScaledRadiusFactors( float* scaledRadiusFactors );
#endif
      int setScaledRadiusFactors( const RealVector& scaledRadiusFactors );
        
      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether scale factors arra should be deleted
      
         @param ownScaledRadiusFactors flag indicating whether scale factors 
                                       array should be deleted

         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setOwnScaleFactors( int ownScaledRadiusFactors );
      
      /**--------------------------------------------------------------------------------------- 
      
         Assign standard radii for GB/SA methods other than ACE;
         taken from Macromodel and OPLS-AA, except for hydrogens (Simbios)
      
         Logic based on logic in Tinker's ksolv.f
      
         Currently only works for standard amino acid atoms
         If invalid atom name is encountered, a message is printed to log file and the
         radius for that atom is set to 1.0f
      
         @param numberOfAtoms       number of atoms
         @param atomNames           array of atom names from GMX top data struct
         @param radii               array to store Macromodel radii for each atom
         @param log                 if set, then print error messages to log file
      
         @return SimTKOpenMMCommon::DefaultReturn always
      
         --------------------------------------------------------------------------------------- */
      
      int getMacroModelAtomicRadii( int numberOfAtoms,
                                    char*** atomNames, Real* radii, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Get AtomicRadii array w/ dielectric offset applied
      
         @return array of atom volumes
      
         --------------------------------------------------------------------------------------- */

      Real* getAtomicRadii( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii array of atomic radii
         @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                        SimTKOpenMMCommon::MdUnits 
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setAtomicRadii( Real* atomicRadii, int units = SimTKOpenMMCommon::MdUnits  );

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii vector of atomic radii
         @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                        SimTKOpenMMCommon::MdUnits 
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setAtomicRadii( const RealVector& atomicRadii, int units = SimTKOpenMMCommon::MdUnits );

      /**---------------------------------------------------------------------------------------
      
         Map Gmx atom name to Tinker atom number (Simbios)
      
         @param atomName            atom name (CA, HA, ...); upper and lower case should both work
         @param log                 if set, then print error messages to log file
      
         return Tinker atom number if atom name is valid; else return -1
      
         --------------------------------------------------------------------------------------- */
            
      int mapGmxAtomNameToTinkerAtomNumber( const char* atomName, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
            
         Get string w/ state
         
         @param title               title (optional)
            
         @return string
            
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;

      /**---------------------------------------------------------------------------------------
                  
         Return zero value if all parameters set; else return nonzero
               
         @return ready status
                  
         --------------------------------------------------------------------------------------- */
       
      int isNotReady( void ) const;
      
};
   
/**---------------------------------------------------------------------------------------
      
   Qsort/heapsort integer comparison (Simbios) 
      
   @parma a first value to compare
   @param b second value to compare

   @return -1, 0, 1
      
--------------------------------------------------------------------------------------- */

int integerComparison( const void *a, const void *b);

// ---------------------------------------------------------------------------------------

#endif // __ObcParameters_H__
