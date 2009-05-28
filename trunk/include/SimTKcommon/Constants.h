#ifndef SimTK_SimTKCOMMON_CONSTANTS_H_
#define SimTK_SimTKCOMMON_CONSTANTS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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
 *
 * This header is the common gathering place for numerical, *machine-independent*
 * constants in SimTK. These include unitless mathematical constants like Pi,
 * as well as physical constants and unit conversion factors. This is NOT the
 * place for computational, machine-dependent constants like NaN,
 * Infinity, number of digits in a float, etc. Those belong in templatized
 * classes in the manner of the C++ std::numeric_limits<T> class.
 *
 * These constants are provided at extremely high precision as compile-time 
 * macros in long double precision. By using very high precision
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
 *
 * This file is self-contained and can be included in ANSI C programs
 * as well as C++.
 *
 * <pre>
 * Unit systems
 *
 *            SI (MKS)           MD                     KCAL-ANGSTROM
 * ---------  --------------  ------------------------  ------------------
 * length     meter           nanometer                 angstrom (Å)
 * mass       kg              amu, dalton               amu, dalton
 * time       second          picosecond                picosecond
 * charge     coulomb         e, proton charge          e, proton charge
 * temp.      kelvin          kelvin                    kelvin
 * substance  mole            mole                      mole
 *
 * velocity   m/s             km/s (nm/ps)              100m/s (Å/ps)
 *
 * energy     J (kg-m^2/s^2)  kJ/mol                    kcal/mol 
 *                              (Da-nm^2/ps^2)            (418.4 Da-Å^2/ps^2)
 * force      N (kg-m/s^2)    kJ/(mol-nm) = TN/mol      kcal/(mol-Å)
 *                              (Da-nm/ps^2) (T=10^12)    (418.4 Da-Å/ps^2)
 *
 * </pre>
 *
 * We always keep angles in radians internally, which are unitless. However,
 * most humans prefer degrees where 1 degree = Pi/180 radians so we provide
 * convenient conversions.
 */

    /**************************/
    /* MATHEMATICAL CONSTANTS */
    /**************************/

/** @name Mathematical Constants
 *
 * These are some common unitless numerical constants evaluated to 64 digits and
 * written here in maximal (long double) precision. (These values were generated using
 * the symbolic calculator Maple which is part of Matlab's Symbolic 
 * Toolbox.) These can be cast to lower precisions when needed, and can be used
 * in compile-time constant expressions like 2*SimTK_PI or 1/SimTK_SQRT2 for which the
 * compiler will properly calculate a long double result with no runtime cost.
 * 
 * These constants are also available as type-safe, 
 * already-rounded, precision-templatized values with static memory addresses
 * as part of our scalar system (see NTraits<T>). You should use the
 * templatized versions when possible. The templatized versions also contain
 * more elaborate constants such as NaN, Infinity, and "epsilon" (machine precision)
 * which can only be generated for specific types.
 */

/*@{    Mathematical Constants */

/** The ratio pi of a circle's circumference to its diameter in Euclidean geometry.
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_PI     3.141592653589793238462643383279502884197169399375105820974944592L

/** e, or exp(1).
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_E      2.718281828459045235360287471352662497757247093699959574966967628L

/** The natural (base e) logarithm of 2.
 * @par uncertainty
 *      approximation of an exact value
 * @see SimTK_E
 */
#define SimTK_LN2    6.931471805599453094172321214581765680755001343602552541206800095e-1L

/** The natural (base e) logarithm of 10.
 * @par uncertainty
 *      approximation of an exact value
 * @see SimTK_E
 */
#define SimTK_LN10   2.302585092994045684017991454684364207601101488628772976033327901L

/** log2(e).
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_LOG2E  1.442695040888963407359924681001892137426645954152985934135449407L

/** log10(e).
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_LOG10E 4.342944819032518276511289189166050822943970058036665661144537832e-1L

/** The square root of 2.
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_SQRT2  1.414213562373095048801688724209698078569671875376948073176679738L

/** The square root of 3.
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_SQRT3  1.732050807568877293527446341505872366942805253810380628055806979L

/** The cube root of 2.
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_CBRT2  1.259921049894873164767210607278228350570251464701507980081975112L

/** The cube root of 3.
 * @par uncertainty
 *      approximation of an exact value
 */
#define SimTK_CBRT3  1.442249570307408382321638310780109588391869253499350577546416195L

/*@}    end of Mathematical Constants */

    /**********************/
    /* PHYSICAL CONSTANTS */
    /**********************/

/** @name Physical Constants
 *
 * @par Provenance
 * These constants are from the CODATA 2002 set from the NIST Physics Laboratory web site
 * http://physics.nist.gov/constants. (NIST SP 961 Dec,2005)
 * Ref: P.J. Mohr and B.N. Taylor, Rev. Mod. Phys. 77(1) (2005).
 * 
 * @par Uncertainty
 * Uncertainties are given in the CODATA set as the one-std-deviation uncertainty in the
 * last 2 digits of the given value. That means that there is about a 68% chance that
 * the last two digits are as shown +/- the uncertainty.
 *
 * How to combine uncertainties (extracted from
 * http://physics.nist.gov/cuu/Uncertainty/combination.html):
 * Assume measured quantities are x1, y1 with u1=uncertainty(x1), u2=uncertainty(x2).
 * We want to combine them into a new quantity y and calculate u=uncertainty(y).
 * <pre>
 * Addition rule: y    = a1*x1 + a2*x2 for factors a1,a2.
 *  then          u^2  = a1*u1^2 + a2*u2^2
 * Multiplication rule y = a*x1^e1*x2^e2 for factor a and exponents e1,e2.
 * Let ur1=u1/|x1|, ur2=u2/|x2| be the relative uncertainties, ur is u/|y|.
 *  then          ur^2 = e1^2*ur1^2 + e2^2*ur2^2, u = ur*|y|
 * </pre>
 */

/*@{    Physical Constants */

/** 
 * Avogadro's number (NA) is defined as the number of atoms in 12g of pure Carbon-12 in
 * its unbound, rest state. The number is 1 mole (mol).
 * @par uncertainty
 *      10e16
 */
#define SimTK_AVOGADROS_NUMBER 6.0221415e23L

/**
 * Mass of a proton in MD units.
 *
 * The atomic mass unit u (or amu) is defined as 1/12 of the mass of a Carbon-12 atom,
 * unbound and in its rest state. This definition matched to Avogadro's number's definition
 * ensures that 1 mole of particles of mass 1u each has total mass exactly 1g. This is
 * synonymous with the dalton (Da), with units of g/mole, so 1u = 1Dalton = 1g/mole.
 * We will use Da for this mass unit, with kDa being a common mass measure for 
 * large biomolecules.
 *
 * @par uncertainty
 *      13e-11
 * @see SimTK_AVOGADROS_NUMBER
 */
#define SimTK_MASS_OF_PROTON_IN_MD 1.00727646688L

/**
 * Mass of a neutron in MD units.
 * @par uncertainty
 *      55e-11
 * @see SimTK_MASS_OF_PROTON_IN_MD
 */
#define SimTK_MASS_OF_NEUTRON_IN_MD 1.00866491560L

/**
 * Mass of an electron in MD units.
 * @par uncertainty
 *      24e-14
 * @see SimTK_MASS_OF_PROTON_IN_MD
 */
#define SimTK_MASS_OF_ELECTRON_IN_MD 5.4857990945e-4L

/**
 * Atomic charge unit e expressed in MKS unit of Coulombs.
 * The charge on an electron is just the negative of this value.
 * @par uncertainty
 *      14e-27
 */
#define SimTK_CHARGE_OF_PROTON_IN_SI 1.60217653e-19L

/**
 * Atomic charge unit e expressed in MD units, which uses e as its charge unit!
 * The charge on an electron is just the negative of this value.
 * @par uncertainty
 *      exact (duh!)
 */
#define SimTK_CHARGE_OF_PROTON_IN_MD 1.L

/**
 * The charge of 1 mole of protons, expressed in Coulombs.
 * <pre>
 *    1.60217653(14)e-19 C/e * 6.0221415(10)e23 = 9.6485338(18)e+4
 * </pre>
 * @par uncertainty
 *      18e-3 
 */
#define SimTK_MOLAR_CHARGE_IN_SI 9.6485338e+4L

/**
 * The charge of 1 mole of protons, expressed in MD units where the unit
 * of charge is just the charge on one proton. So in MD units this is just
 * Avogadro's number.
 * @see SimTK_AVOGADROS_NUMBER
 */
#define SimTK_MOLAR_CHARGE_IN_MD SimTK_AVOGADROS_NUMBER

/**
 * Speed of light c is exact in MKS units of m/s.
 * @par uncertainty
 *      exact
 * @see SimTK_LIGHTSPEED_IN_MD
 */
#define SimTK_LIGHTSPEED_IN_SI 2.99792458e+8L

/**
 * Speed of light c is exact in MD units of nm/ps.
 * @par uncertainty
 *      exact
 * @see SimTK_LIGHTSPEED_IN_SI
 */
#define SimTK_LIGHTSPEED_IN_MD 2.99792458e+5L

/** 
 * Newton's gravitational constant G in N-m^2/kg^2 = m^3 kg^-1 s^-2.
 * The force between two point masses m1,m2 separated by a distance d is
 *     <pre> F = -G m1*m2/d^2 </pre> 
 * (with the "-" indicating an attractive force).
 * @par uncertainty
 *      10e-15
 * @see SimTK_GRAVITATIONAL_CONSTANT_IN_MD
 */
#define SimTK_GRAVITATIONAL_CONSTANT_IN_SI 6.6742e-11L   

/**
 * Newton's gravitational constant G in (kJ/mol)-nm^2/u^2 = nm^3 u^-1 ps^-2.
 * <pre>
 * Conversion is (nm/m)^3 (u/kg)^-1 (ps/s)^-2
 *          = 1.66053886(28)e-24L     (uncertainty: 28e-32)
 * </pre>
 * This is why gravity doesn't matter in molecular systems. Don't try
 * this in single precision -- you'll run out of exponent!
 * @par uncertainty
 *      17e-39
 * @see SimTK_GRAVITATIONAL_CONSTANT_IN_SI
 */
#define SimTK_GRAVITATIONAL_CONSTANT_IN_MD 1.10827e-34L

/**
 * Free space magnetic permeability constant mu0 in SI units (exact).
 * <pre>
 *   = 4*pi * 1e-7 exactly in N/A^2 (Newtons/Ampere^2) = kg-m/C^2
 * </pre>
 * @par uncertainty
 *      approximation of an exact quantity
 * @see SimTK_ELECTRIC_PERMITTIVITY_IN_SI
 */
#define SimTK_MAGNETIC_PERMEABILITY_IN_SI  \
    1.256637061435917295385057353311801153678867759750042328389977837e-6L

/**
 * Free space magnetic permeability constant mu0 in MD units (not exact).
 * <pre>
 * Convert kg->g/mole, m->nm, C->e = (4*pi*1e5)*1.60217653e-19^2*6.0221415e23
 *      (exact in SI units, but not exact here)
 * </pre>
 * @par uncertainty
 *      47e-16
 * @see SimTK_ELECTRIC_PERMITTIVITY_IN_MD
 */
#define SimTK_MAGNETIC_PERMEABILITY_IN_MD   1.94259179e-8L

/**
 * Free space permittivity constant e0 = 1/(mu0*c^2) Farad/m = Coulomb^2/(N-m^2) (exact in SI units).
 * @par uncertainty
 *      approximation of an exact quantity
 * @see SimTK_MAGNETIC_PERMEABILITY_IN_SI
 */
#define SimTK_ELECTRIC_PERMITTIVITY_IN_SI  \
    8.854187817620389850536563031710750260608370166599449808102417149e-12L  /* approx of exact */

/**
 * Free space permittivity constant e0=1/(mu0*c^2) e^2/(kN-nm^2) using MD permeability and
 * MD lightspeed.
 * @par uncertainty
 *      14e-11
 * @see SimTK_MAGNETIC_PERMEABILITY_IN_MD
 */
#define SimTK_ELECTRIC_PERMITTIVITY_IN_MD 5.7276575e-4L

/**
 * Coulomb's constant kappa = 1/(4pi*e0)=1e-7*c^2 N-m^2/Coulomb^2 (exact in SI units).
 * This is the constant that appears in Coulomb's law f(r)= kappa*q1*q2/r^2.
 * @par uncertainty
 *      exact
 */
#define SimTK_COULOMB_CONSTANT_IN_SI 8.9875517873681764e+9L

/**
 * Coulomb's constant kappa = 1/(4*pi*e0) in MD units.
 * This is the constant that appears in Coulomb's law f(r)= kappa*q1*q2/r^2.
 * <pre>
 * Coulomb's consant in MD units uses MD e0 & c: 
 *   1/(4*pi*e0)=1e5*1.60217653e-19^2*6.0221415e23*c^2 kN-nm^2/e^2 (=kJ-nm/e^2)
 *     (exact in SI units but not exact in MD)
 * </pre>
 * @par uncertainty
 *      33e-6
 */
#define SimTK_COULOMB_CONSTANT_IN_MD 1.38935456e+2L

/**
 * Coulomb's constant kappa = 1/(4*pi*e0) in kcal-Angstroms/e^2.
 * This is the constant that appears in Coulomb's law f(r)= kappa*q1*q2/r^2.
 * This is an exact conversion from MD units (which are inexact).
 * @par uncertainty
 *      80e-6
 */
#define SimTK_COULOMB_CONSTANT_IN_KCAL_ANGSTROM 3.32063711e+2L

/**
 * This is the gas constant R in (J/mol)/K. 
 * @par uncertainty
 *      15e-6
 */
#define SimTK_MOLAR_GAS_CONSTANT_SI 8.314472L

/**
 * This is the gas constant R in (kJ/mol)/K. 
 * This is an exact conversion from SI units, differing only in the use of kJ 
 * here vs. J in SI.
 * @par uncertainty
 *      15e-9
 */
#define SimTK_MOLAR_GAS_CONSTANT_MD 8.314472e-3L

/**
 * This is the gas constant R in (kcal/mol)/K. 
 * This is an exact conversion from MD units, differing only in the use of kcal 
 * here vs. kJ in MD.
 * @par uncertainty
 *      36e-10
 */
#define SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM 1.9872065e-3L

/**
 * Boltzmann's constant in SI units of joules/kelvin; just divide R by NA.
 * @par uncertainty
 *      24e-30
 */
#define SimTK_BOLTZMANN_CONSTANT_SI  1.3806505e-23L

/**
 * Boltzmann's constant in MD units of (kJ/mol)/kelvin; same as R.
 * @see SimTK_MOLAR_GAS_CONSTANT_MD
 */
#define SimTK_BOLTZMANN_CONSTANT_MD  SimTK_MOLAR_GAS_CONSTANT_MD

/**
 * Boltzmann's constant in Kcal-Angstrom units of (kcal/mol)/kelvin; same as R.
 * @see SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM
 */
#define SimTK_BOLTZMANN_CONSTANT_KCAL_ANGSTROM SimTK_MOLAR_GAS_CONSTANT_KCAL_ANGSTROM

/*@}    end of Physical Constants */

    /***************************/
    /* UNIT CONVERSION FACTORS */
    /***************************/

/** @name Unit Conversion Factors
 *
 * In each case here, given a value in the units mentioned first in the name,
 * you should multiply by the given constant to produce the equivalent
 * quantity measured in the units that appear second in the name. You can 
 * perform the reverse conversion by dividing by the constant, or by using
 * another conversion constant with the names reversed if one is supplied here.
 */

/*@{    Unit Conversion Factors*/

/** 
 * Convert radians to degrees.
 * @par uncertainty
 *       approximation of an exact quantity
 * @see SimTK_DEGREE_TO_RADIAN
 */
#define SimTK_RADIAN_TO_DEGREE 5.729577951308232087679815481410517033240547246656432154916024386e+1L

/** 
 * Convert degrees to radians.
 * @par uncertainty
 *       approximation of an exact quantity
 * @see SimTK_RADIAN_TO_DEGREE
 */
#define SimTK_DEGREE_TO_RADIAN 1.745329251994329576923690768488612713442871888541725456097191440e-2L

/** 
 * Convert Kcal to Kjoule (also Kcal/mol to Kjoule/mol).
 * @par uncertainty
 *         exact
 * @see SimTK_KJOULE_TO_KCAL
 */
#define SimTK_KCAL_TO_KJOULE 4.184L /* exact */

/** 
 * Convert Kjoule to Kcal (also Kjoule/mol to Kcal/mol).
 * @par uncertainty
 *         approximation of an exact quantity
 * @see SimTK_KCAL_TO_KJOULE
 */
#define SimTK_KJOULE_TO_KCAL 2.390057361376673040152963671128107074569789674952198852772466539e-1L

/** 
 * Convert atomic mass unit (amu, Dalton) to g. This is 1/NA (NA=avogadro's number).
 * @par uncertainty
 *         28e-32
 * @see SimTK_AVOGADROS_NUMBER
 */
#define SimTK_DALTON_TO_GRAM     1.66053886e-24L

/** 
 * Convert proton charge units to Coulombs. This is the same as the
 * conversion from electron volts to Joules, and both are just the
 * charge of a proton in SI units.
 * @see SimTK_CHARGE_OF_PROTON_IN_SI
 * @see SimTK_EV_TO_JOULE
 */
#define SimTK_E_TO_COULOMB       SimTK_CHARGE_OF_PROTON_IN_SI

/** 
 * Convert electron volts to Joules. This is the same as the
 * conversion from proton charge units to Coulombs, and both are just the
 * charge of a proton in SI units.
 * @see SimTK_CHARGE_OF_PROTON_IN_SI
 * @see SimTK_E_TO_COULOMB
 */
#define SimTK_EV_TO_JOULE        SimTK_CHARGE_OF_PROTON_IN_SI

/*@}    end of Unit Conversion Factors */
    
#endif /* SimTK_SimTKCOMMON_CONSTANTS_H_ */
