#ifndef SimTK_SimTKCOMMON_CONSTANTS_H_
#define SimTK_SimTKCOMMON_CONSTANTS_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS, OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 *
 * This header is the common gathering place for numerical, *machine-independent*
 * constants in SimTK. These include unitless mathematical constants like Pi,
 * as well as physical constants and unit conversion factors. This is NOT the
 * place for computational, machine-dependent constants like NaN,
 * Infinity, number of digits in a float, etc. Those belong in templatized
 * classes in the manner of the C++ std::numeric_limits<T> class.
 *
 * These constants are provided at extremely high precision as compile-time 
 * #define macros in long double precision. By using very high precision
 * we ensure sufficient accuracy for any IEEE long double precision 
 * implementation (they can be 64, 80, or 128 bits). These constants can
 * be used as raw material for providing nicer templatized constants in
 * appropriate precisions and unit systems.
 *
 * Units: our most common unit systems are the "SI" (MKS) system, and the "MD"
 * system used for molecular dynamics. SI units are meters, kg, seconds, 
 * coulombs (ampere-s), kelvins and moles. MD units are nanometers, atomic mass units
 * (Daltons, g/mol), picoseconds, proton charge e, kelvins, and moles. Many
 * molecular dynamicists and chemists prefer kcals for energy and angstroms for
 * length. This does not constitute a consistent set of units, however, so
 * we provide for it by conversion from the MD units, which are consistent.
 * (By consistent, we mean that force units = mass-length/time^2, so f=ma!)
 */

// Unit systems
//
//            SI (MKS)           MD                 KCAL-ANGSTROM
// ---------  --------------  --------------------  ----------------------
// length     meter           nanometer             angstrom
// mass       kg              u, dalton             u, dalton
// time       second          picosecond            picosecond
// charge     coulomb         e, proton charge      e, proton charge
// temp.      kelvin          kelvin                kelvin
// substance  mole            mole                  mole
//
// energy     J (kg-m^2/s^2)  kJ/mol (u-nm^2/ps^2)  kcal/mol (418.4u-A^2/ps^2)
// force      N (kg-m/s^2)    kN/mol (u-nm/ps^2)    kcal/mol-A (418.4u-A/ps^2)
//
// We always keep angles in radians internally, which are unitless. However,
// most humans prefer degrees where 1 degree = Pi/180 radians.
//

    ////////////////////////////
    // MATHEMATICAL CONSTANTS //
    ////////////////////////////

// These are some common unitless numerical constants evaluated to 64 digits and
// written here in maximal (long double) precision. (These values were generated using
// the symbolic calculator Maple which is part of Matlab's Symbolic 
// Toolbox.) These can be cast to lower precisions when needed, and can be used
// in compile-time constant expressions like 2*SimTK_PI or 1/SimTK_SQRT2 for which the
// compiler will properly calculate a long double result with no runtime cost.
// 
// These constants are also available as type-safe, 
// already-rounded, precision-templatized values with static memory addresses
// as part of our scalar system (see NTraits<T>). You should use the
// templatized versions when possible. The templatized versions also contain
// more elaborate constants such as NaN, Infinity, and "epsilon" (machine precision)
// which can only be generated for specific types.

#define SimTK_PI     3.141592653589793238462643383279502884197169399375105820974944592L
#define SimTK_E      2.718281828459045235360287471352662497757247093699959574966967628L

// log2(e), log10(e)
#define SimTK_LOG2E  1.442695040888963407359924681001892137426645954152985934135449407L
#define SimTK_LOG10E 4.342944819032518276511289189166050822943970058036665661144537832e-1L

// sqrt(2), sqrt(3), cubeRoot(2), cubeRoot(3)
#define SimTK_SQRT2  1.414213562373095048801688724209698078569671875376948073176679738L
#define SimTK_SQRT3  1.732050807568877293527446341505872366942805253810380628055806979L
#define SimTK_CBRT2  1.259921049894873164767210607278228350570251464701507980081975112L
#define SimTK_CBRT3  1.442249570307408382321638310780109588391869253499350577546416195L

// ln(2), ln(10)
#define SimTK_LN2    6.931471805599453094172321214581765680755001343602552541206800095e-1L
#define SimTK_LN10   2.302585092994045684017991454684364207601101488628772976033327901L


    ////////////////////////
    // PHYSICAL CONSTANTS //
    ////////////////////////

// Provenance
// ----------
// These constants are from the CODATA 2002 set from the NIST Physics Laboratory web site
// physics.nist.gov/constants. (NIST SP 961 Dec,2005)
// Ref: P.J. Mohr and B.N. Taylor, Rev. Mod. Phys. 77(1) (2005).
// 
// Uncertainty
// -----------
// Uncertainties are given in the CODATA set as the one-std-deviation uncertainty in the
// last 2 digits of the given value. That means that there is about a 68% chance that
// the last two digits are as shown +/- the uncertainty.
//
// How to combine uncertainties (extracted from
// http://physics.nist.gov/cuu/Uncertainty/combination.html):
// Assume measured quantities are x1, y1 with u1=uncertainty(x1), u2=uncertainty(x2).
// We want to combine them into a new quantity y and calculate u=uncertainty(y).
// Addition rule: y    = a1*x1 + a2*x2 for factors a1,a2.
//  then          u^2  = a1*u1^2 + a2*u2^2
// Multiplication rule y = a*x1^e1*x2^e2 for factor a and exponents e1,e2.
// Let ur1=u1/|x1|, ur2=u2/|x2| be the relative uncertainties, ur is u/|y|.
//  then          ur^2 = e1^2*ur1^2 + e2^2*ur2^2, u = ur*|y|


// Avogadro's number (NA) is defined as the number of atoms in 12g of pure Carbon-12 in
// its unbound, rest state.
#define SimTK_AVOGADROS_NUMBER              6.0221415e23L       // uncertainty: 10e16

// The atomic mass unit u is defined as 1/12 of the mass of a Carbon-12 atom, unbound
// and in its rest state. This definition matched to Avogadro's number's definition
// ensures that 1 mole of particles of mass 1 u has total mass exactly 1g. This is
// synonymous with the Dalton, with units of g/mole, so 1u = 1Dalton = 1g/mole.

#define SimTK_MASS_OF_ELECTRON_IN_MD        5.4857990945e-4L    // uncertainty: 24e-14
#define SimTK_MASS_OF_PROTON_IN_MD          1.00727646688L      // uncertainty: 13e-11
#define SimTK_MASS_OF_NEUTRON_IN_MD         1.00866491560L      // uncertainty: 55e-11

// Atomic charge unit e expressed in MKS unit of Coulombs
#define SimTK_CHARGE_OF_PROTON_IN_SI        1.60217653e-19L     // uncertainty: 14e-27
#define SimTK_CHARGE_OF_PROTON_IN_MD        1.L                 // exact (duh!)

// This is the charge of 1 mole of protons, expressed in Coulombs
//    1.60217653(14)e-19 C/e * 6.0221415(10)e23 = 9.6485338(18)e+4
#define SimTK_MOLAR_CHARGE_IN_SI            9.6485338e+4L       // uncertainty: 18e-3
#define SimTK_MOLAR_CHARGE_IN_MD            SimTK_AVOGADROS_NUMBER

// Speed of light c is exact in MKS units (m/s), or in MD (nm/ps)
#define SimTK_LIGHTSPEED_IN_SI              2.99792458e+8L      // exact
#define SimTK_LIGHTSPEED_IN_MD              2.99792458e+5L      // exact

// The force between two point masses m1,m2 separated by a distance d is
//       F = -G m1*m2/d^2  ("-" indicating an attractive force)
// Newton's gravitational constant G in  N-m^2/kg^2 = m^3 kg^-1 s^-2
#define SimTK_GRAVITATIONAL_CONSTANT_IN_SI  6.6742e-11L         // uncertainty: 10e-15    

// Newton's gravitational constant G in  (kJ/mol)-nm^2/u^2 = nm^3 u^-1 ps^-2
// Conversion is (nm/m)^3 (u/kg)^-1 (ps/s)^-2
//          = 1.66053886(28)e-24L    // uncertainty: 28e-32
// This is why gravity doesn't matter in molecular systems. Don't try
// this in single precision -- you'll run out of exponent!
#define SimTK_GRAVITATIONAL_CONSTANT_IN_MD  1.10827e-34L        // uncertainty: 17e-39

// Free space magnetic permeability constant mu0
//   = 4*pi * 1e-7 exactly in N/A^2 (Newtons/Ampere^2) = kg-m/C^2
#define SimTK_MAGNETIC_PERMEABILITY_IN_SI  \
    1.256637061435917295385057353311801153678867759750042328389977837e-6L   // approx of exact

// Convert kg->g/mole, m->nm, C->e = (4*pi*1e5)*1.60217653e-19^2*6.0221415e23 (not exact)
#define SimTK_MAGNETIC_PERMEABILITY_IN_MD   1.94259179e-8L      // uncertainty: 47e-16

// Free space permittivity constant e0 = 1/(mu0*c^2) Farad/m = Coulomb^2/(N-m^2) (exact in SI units)
#define SimTK_ELECTRIC_PERMITTIVITY_IN_SI  \
    8.854187817620389850536563031710750260608370166599449808102417149e-12L  // approx of exact

// Use MD permeability and MD lightspeed: e0=1/(mu0*c^2) e^2/(kN-nm^2)
#define SimTK_ELECTRIC_PERMITTIVITY_IN_MD   5.7276575e-4L       // uncertainty: 14e-11

// Coulomb's constant kappa = 1/(4pi*e0)=1e-7*c^2 N-m^2/Coulomb^2 (exact in SI units).
// This is the constant that appears in Coulomb's law f(r)= kappa*q1*q2/r^2
#define SimTK_COULOMB_CONSTANT_IN_SI       8.9875517873681764e+9L   // exact

// Coulomb's consant in MD units uses MD e0 & c: 
// kappa = 1/(4*pi*e0)=1e5*1.60217653e-19^2*6.0221415e23*c^2 kN-nm^2/e^2 (=kJ-nm/e^2) (not exact)
#define SimTK_COULOMB_CONSTANT_IN_MD       1.38935456e+2L       // uncertainty: 33e-6

// Here is the same thing in kcal-Angstroms/e^2 (exact conversion from MD units)
#define SimTK_COULOMB_CONSTANT_IN_KCAL_ANGSTROM 3.32063711e+2L  // uncertainty: 80e-6

// This is the gas constant R in (J/mol)/K.
#define SimTK_MOLAR_GAS_CONSTANT_SI  8.314472L                  // uncertainty: 15e-6

// In this case MD and MKS differ only by kJ vs J: (kJ/mol)/K here.
#define SimTK_MOLAR_GAS_CONSTANT_MD  8.314472e-3L               // uncertainty: 15e-9

// Just need to convert kJ to kcal here
#define SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM 1.9872065e-3L    // uncertainty: 36e-10

// Units are joules/kelvin in SI; i.e. divide R by NA
#define SimTK_BOLTZMANN_CONSTANT_SI  1.3806505e-23L             // uncertainty: 24e-30

// Units are (kJ/mol)/kelvin; same as R
#define SimTK_BOLTZMANN_CONSTANT_MD  SimTK_MOLAR_GAS_CONSTANT_MD

// Units are (kcal/mol)/kelvin; same as R
#define SimTK_BOLTZMANN_CONSTANT_KCAL_ANGSTROM SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM


    /////////////////////////////
    // UNIT CONVERSION FACTORS //
    /////////////////////////////

// In each case you should *multiply* by the given constant to perform the 
// indicated conversion; divide by the constant to perform the inverse conversion.

#define SimTK_RADIAN_TO_DEGREE 5.729577951308232087679815481410517033240547246656432154916024386e+1L
#define SimTK_DEGREE_TO_RADIAN 1.745329251994329576923690768488612713442871888541725456097191440e-2L

// Kcal <-> Kjoule conversion is exact (4.184 and 1/4.184)
// Of course these work for kcal/mol <-> kj/mol too.
#define SimTK_KCAL_TO_KJOULE     4.184L // exact
#define SimTK_KJOULE_TO_KCAL     2.390057361376673040152963671128107074569789674952198852772466539e-1L

// Atomic mass unit (a.k.a. Dalton) to g. This is 1/NA (NA=avogadro's number); see below.
#define SimTK_U_TO_GRAM          1.66053886e-24L    // uncertainty: 28e-32

// Proton charge unit to Coulomb = electron volt to Joule
#define SimTK_E_TO_COULOMB       SimTK_CHARGE_OF_PROTON_IN_SI
#define SimTK_EV_TO_JOULE        SimTK_CHARGE_OF_PROTON_IN_SI

    
#endif // SimTK_SimTKCOMMON_CONSTANTS_H_
