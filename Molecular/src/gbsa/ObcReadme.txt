

The OBC Generalized Born model is based on the following papers:
 
      J. Phys. Chem. 1996 100, 19824-19839 (HCT paper)

      Proteins: Structure, Function, and Bioinformatics 55:383-394 (2004) (OBC paper)

---------------------------------------------------------------------------------------

The two main interface methods perform the following tasks: 

	(1) set the GBSA parameters; this method is called once prior to execution of
       the main simulation loop

	(2) calculate the GBSA forces/energy given the atom coordinates and charges
       for each simulation step

The two methods and helper methods are in the file cpuObcInterface.cpp

---------------------------------------------------------------------------------------

The setup call is

   cpuSetObcParameters( int numberOfAtoms, RealOpenMM* atomicRadii, RealOpenMM* obcScaleFactors,
                        int includeAceApproximation,
                        RealOpenMM soluteDielectric, RealOpenMM solventDielectric, FILE* log );

The parameters include the following:

	(1)   number of atoms
	(2)   the OBC atomic radii (akin to VDW radii, but optimized for GBSA calculations)
	(3)   the OBC scale factors for each atom
	(4)   a flag indicating whether the nonpolar ACE approximation is to be included 
         in calculating the forces and energy: I do not have a reference for this term
	(5-6) the solute and solvent dielectric (typically 1.0 & 78.3)
	(7)   a log file reference; the reference may be set to NULL,
         if logging is unwanted. Note: if the reference is NULL, warnings, ... will be output to stderr

  cpuSetObcParameters() creates a static CpuObc object which is then referenced when the
  GBSA forces and energy are to be calculated

  The default OBC model is the Type II model: Eq. 8 in the 2004 OBC paper; 
  for the Type I model, replace the call

     ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );

   with 

     ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeI );

---------------------------------------------------------------------------------------

The OBC scale factors may be obtained via calls to either of the following 'helper' methods

      getObcScaleFactors( int numberOfAtoms, const int* atomicNumber, RealOpenMM* scaleFactors )

   or

      getObcScaleFactorsGivenAtomMasses( int numberOfAtoms, const RealOpenMM* masses, RealOpenMM* scaleFactors )

Documentation for the inputs to these routines is contained in the file, cpuObcInterface.cpp.
The scale factors are those used in Tinker5, and are based on the element type
of each atom (H, C, N, O, ...).

---------------------------------------------------------------------------------------

The OBC radii may be obtained via the helper method:

      getGbsaRadii( int numberOfAtoms, const int* atomicNumber, 
                    const int* numberOfCovalentPartners, const int* indexOfCovalentPartner,
                    RealOpenMM* gbsaRadii );

Documentation for the inputs to this routine is contained in the file, cpuObcInterface.cpp.
The radii are based on the atom type (methyl C, carbonyl C, ...). The value for a hydrogen
depends on its covalent heavy atom partner. As for the OBC scale factors, the values 
hardwired into the method are the same as those used in Tinker5.
 
---------------------------------------------------------------------------------------

The energy/forces are calculated via the call

   cpuCalculateImplicitSolventForces( RealOpenMM** atomCoordinates, const RealOpenMM* partialCharges,
                                      RealOpenMM** forces, RealOpenMM* energy, int updateBornRadii );

	The coordinate and force arrays should have the layout: atomCoordinates[atomIndex][3] &  
                                                                    forces[atomIndex][3]

   If updateBornRadii = 0, then the Born radii from the previous iteration are used; this reduces the
   number of O(N**2) loops from 3 to 2.

   If updateBornRadii != 0, then the Born radii are calculated for the input configuration. 

   For MD runs, setting updateBornRadii = 0 is relatively safe since the radii change little 
   from one iteration to the next. However, for minimization techniques updateBornRadii 
   should probably be set to a nonzero value since the configuarations can change significantly
   between adjacent calls.

---------------------------------------------------------------------------------------
 
Notes:

	(1) The type (float or double) for floating-point numbers is set in the file SimTKOpenMMRealType.h.
       If 'RealOpenMMType' is set to 1 in the file, the calculations are performed in
       single-precision; otherwise they are performed in double-precision

   (2) The units are Angstroms for the spatial dimensions, kcal/mol for energy, and kcal/mol.A for the forces

   (3) The code in gromacsCpuObcInterface.cpp provides an example of the parameter setup, ... used for
       interfacing w/ Gromacs; the details are somewhat different from the code in cpuObcInterface.cpp

   (4) Output from Tinker for a singles ubiquitin configuration are in the directory 'Examples'.

       The xyz file gives the coordinates and atom types (Amber99 -- note the charge for atom type 268
      (line 2008 in the Tinker amber99.prm file) was changed from -0.2737 to -0.2774 to agree w/ 
      the value used in the version of Gromacs we have access to.) 

       The key file specifies that OBC should be used for the implicit solvent.

       The results are in the files 111UBQ.gbsaAce and 111UBQ.gbsaNoAce. The first file includes the
       ACE approximation, while the second does not. On the first line of each file, the number of atoms
       and GBSA energy are reported. The remaining lines (one for each atom) give the coordinates (3 entries),
       Born radius, partial charge, GBSA radius (rsolv in Tinker), the OBC scaling factor, and the forces (3). The
       last three fields give the atomic number, Tinker Amber99 atom name and type.

---------------------------------------------------------------------------------------
