

The OBC Generalized Born model is based  on the following papers:
 
      J. Phys. Chem. 1996 100, 19824-19839 (HCT paper)
      Proteins: Structure, Function, and Bioinformatics 55:383-394 (2004) (OBC paper)

---------------------------------------------------------------------------------------

There are two main interface methods: 

	(1) setup the parameters once prior to excuting the main simulation loop 
	(2) calculate the GBSA forces given the atom coordinates

The methods and helper methods are in cpuObcInterface.cpp

---------------------------------------------------------------------------------------

The setup call is

   cpuSetObcParameters( int numberOfAtoms, Real* atomicRadii, Real* obcScaleFactors,
                        int includeAceApproximation,
                        Real soluteDielectric, Real solventDielectric, FILE* log );

The parameters include the following:

	(1) number of atoms
	(2) the solute and solvent dielectric (typically 1.0 & 78.3)
	(3) a flag indicating whether the nonpolar ACE approximation is to be included 
       in calculating the forces and energy
	(4) the OBC scale factors for each atom
	(5) the OBC atomic radii (akin to VDW radii, but different)
	(6) a log file reference (may be NULL if logging is not wanted: warnings, ... will be output to stderr
            in this case )

  cpuSetObcParameters() creates a static CpuObc object which is then used to calculate the
  forces and energy

  The default OBC model is the Type II model: Eq. 8 in the paper; for the Type I model, replace 
  the call

     ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );

   with 

     ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeI );

---------------------------------------------------------------------------------------

The OBC scale factors may be obtained via calls to either

      getObcScaleFactors( int numberOfAtoms, const int* atomicNumber, Real* scaleFactors )

   or

      getObcScaleFactorsGivenAtomMasses( int numberOfAtoms, const Real* masses, Real* scaleFactors )

Documentation for the inputs to these routines is contained in the file.

---------------------------------------------------------------------------------------

The OBC radii may be obtained via

      getGbsaRadii( int numberOfAtoms, const int* atomicNumber, 
                    const int* numberOfCovalentPartners, const int* indexOfCovalentPartner, Real* gbsaRadii );

Documentation for the inputs to these routines is contained in the file.
 
---------------------------------------------------------------------------------------

The energy/forces are calculated via the call

   cpuCalculateImplicitSolventForces( Real** atomCoordinates, const Real* partialCharges,
                                      Real** forces, Real* energy );

	The coordinate and force arrays should have the layout: atomCoordinates[atoms][3] & forces[atoms][3]

---------------------------------------------------------------------------------------
 
Notes:

	(1) The type 'Real' is set in SimTKUtilities/SimTKOpenMMRealType.h to specify 
       whether float or double is to be used

   (2) The units are Angstroms for the spatial dimension and kcal/mol.A for the forces

