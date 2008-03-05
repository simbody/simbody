#ifndef SimTK_SIMBODY_ELEMENT_H_
#define SimTK_SIMBODY_ELEMENT_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Christopher Bruns                                *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

namespace SimTK {

// These typedefs are a documentation tool for method signatures
typedef Real Distance; // in nanometers
typedef Real Angle; // in radians
typedef Real Mass; // in atomic mass units


// Element is a chemical element from the periodic table of
// the elements, e.g. hydrogen, carbon, gold, etc.
class Element;
class ElementRep;

// Avoid instantiating template class except in one particular library compilation unit
#ifndef DO_INSTANTIATE_ELEMENT_PIMPL_HANDLE
    extern template class SimTK_SIMBODY_EXPORT PIMPLHandle<Element, ElementRep>;
#endif

class SimTK_SIMBODY_EXPORT Element : public PIMPLHandle<SimTK_SIMBODY_EXPORT Element,ElementRep> {
public:
    typedef String Name; // e.g. "gold"
    typedef String Symbol; // e.g. "Au"
    static const int InvalidAtomicNumber = -1;

    Element();
    Element(int atomicNumber, Name name, Symbol symbol, Mass typicalMass);

    Symbol getSymbol() const;
    Name getName() const;

    static Element getByAtomicNumber(int atomicNumber);

    class Hydrogen;
    class Deuterium;
    class Helium;
    class Lithium;
    class Beryllium;
    class Boron;
    class Carbon;
    class Nitrogen;
    class Oxygen;
    class Fluorine;
    class Neon;
    class Sodium;
    class Magnesium;
    class Aluminum;
    class Silicon;
    class Phosphorus;
    class Sulfur;
    class Chlorine;
    class Argon;
    class Potassium;
    class Calcium;
    class Scandium;
    class Titanium;
    class Vanadium;
    class Chromium;
    class Manganese;
    class Iron;
    class Cobalt;
    class Nickel;
    class Copper;
    class Zinc;
    class Gallium;
    class Germanium;
    class Arsenic;
    class Selenium;
    class Bromine;
    class Krypton;
    class Rubidium;
    class Strontium;
    class Yttrium;
    class Zirconium;
    class Niobium;
    class Molybdenum;
    class Technetium;
    class Ruthenium;
    class Rhodium;
    class Palladium;
    class Silver;
    class Cadmium;
    class Indium;
    class Tin;
    class Antimony;
    class Tellurium;
    class Iodine;
    class Xenon;
    class Cesium;
    class Barium;
    class Lanthanum;
    class Cerium;
    class Praseodymium;
    class Neodymium;
    class Promethium;
    class Samarium;
    class Europium;
    class Gadolinium;
    class Terbium;
    class Dysprosium;
    class Holmium;
    class Erbium;
    class Thulium;
    class Ytterbium;
    class Lutetium;
    class Hafnium;
    class Tantalum;
    class Tungsten;
    class Rhenium;
    class Osmium;
    class Iridium;
    class Platinum;
    class Gold;
    class Mercury;
    class Thallium;
    class Lead;
    class Bismuth;
    class Polonium;
    class Astatine;
    class Radon;
    class Francium;
    class Radium;
    class Actinium;
    class Thorium;
    class Protactinium;
    class Uranium;
    class Neptunium;
    class Plutonium;
    class Americium;
    class Curium;
    class Berkelium;
    class Californium;
    class Einsteinium;
    class Fermium;
    class Mendelevium;
    class Nobelium;
    class Lawrencium;
    class Rutherfordium;
    class Dubnium;
    class Seaborgium;
    class Bohrium;
    class Hassium;
    class Meitnerium;
    class Darmstadtium;
    class Roentgenium;
    class Ununbium;
    class Ununtrium;
    class Ununquadium;
    class Ununpentium;
    class Ununhexium;

    // typedefs for multiply named elements
    typedef Sulfur Sulphur;
    typedef Aluminum Aluminium;
    typedef Cesium Caesium;
    typedef Darmstadtium Ununnilium;
    typedef Roentgenium Unununium;
};

SimTK_SIMBODY_EXPORT std::ostream& operator<<(std::ostream&, const Element&);

// TODO - populate parameters
class Element::Hydrogen : public Element {public: Hydrogen() : 
    Element(1, "hydrogen", "H", 1.007947) {} };
class Element::Deuterium : public Element {public: Deuterium() :
    Element(1, "deuterium", "D", 2.01355321270) {} };
class Element::Helium : public Element {public: Helium() :
    Element(2, "helium", "He", 4.003) {} };
class Element::Lithium : public Element {public: Lithium() :
    Element(3, "lithium", "Li", 6.9412) {} };
class Element::Beryllium : public Element {public: Beryllium() :
    Element(4, "beryllium", "Be", 9.0121823) {} };
class Element::Boron : public Element {public: Boron() :
    Element(5, "boron", "B", 10.8117) {} };
class Element::Carbon : public Element {public: Carbon() :
    Element(6, "carbon", "C", 12.01078) {} };
class Element::Nitrogen : public Element {public: Nitrogen() :
    Element(7, "nitrogen", "N", 14.00672) {} };
class Element::Oxygen : public Element {public: Oxygen() :
    Element(8, "oxygen", "O", 15.99943) {} };
class Element::Fluorine : public Element {public: Fluorine() :
    Element(9, "fluorine", "F", 18.99840325) {} };
class Element::Neon : public Element {public: Neon() :
    Element(10, "neon", "Ne", 20.17976) {} };
class Element::Sodium : public Element {public: Sodium() :
    Element(11, "sodium", "Na", 22.989769282) {} };
class Element::Magnesium : public Element {public: Magnesium() :
    Element(12, "magnesium", "Mg", 24.30506) {} };
class Element::Aluminum : public Element {public: Aluminum() :
    Element(13, "aluminum", "Al", 26.98153868) {} };
class Element::Silicon : public Element {public: Silicon() :
    Element(14, "silicon", "Si", 28.08553) {} };
class Element::Phosphorus : public Element {public: Phosphorus() :
    Element(15, "phosphorus", "P", 30.9737622) {} };
class Element::Sulfur : public Element {public: Sulfur() :
    Element(16, "sulfur", "S", 32.0655) {} };
class Element::Chlorine : public Element {public: Chlorine() :
    Element(17, "chlorine", "Cl", 35.4532) {} };
class Element::Argon : public Element {public: Argon() :
    Element(18, "argon", "Ar", 39.9481) {} };
class Element::Potassium : public Element {public: Potassium() :
    Element(19, "potassium", "K", 39.09831) {} };
class Element::Calcium : public Element {public: Calcium() :
    Element(20, "calcium", "Ca", 40.0784) {} };
class Element::Scandium : public Element {public: Scandium() :
    Element(21, "scandium", "Sc", 44.9559126) {} };
class Element::Titanium : public Element {public: Titanium() :
    Element(22, "titanium", "Ti", 47.8671) {} };
class Element::Vanadium : public Element {public: Vanadium() :
    Element(23, "vanadium", "V", 50.94151) {} };
class Element::Chromium : public Element {public: Chromium() :
    Element(24, "chromium", "Cr", 51.99616) {} };
class Element::Manganese : public Element {public: Manganese() :
    Element(25, "manganese", "Mn", 54.9380455) {} };
class Element::Iron : public Element {public: Iron() :
    Element(26, "iron", "Fe", 55.8452) {} };
class Element::Cobalt : public Element {public: Cobalt() :
    Element(27, "cobalt", "Co", 58.9331955) {} };
class Element::Nickel : public Element {public: Nickel() :
    Element(28, "nickel", "Ni", 58.69342) {} };
class Element::Copper : public Element {public: Copper() :
    Element(29, "copper", "Cu", 63.5463) {} };
class Element::Zinc : public Element {public: Zinc() :
    Element(30, "zinc", "Zn", 65.4094) {} };
class Element::Gallium : public Element {public: Gallium() :
    Element(31, "gallium", "Ga", 69.7231) {} };
class Element::Germanium : public Element {public: Germanium() :
    Element(32, "germanium", "Ge", 72.641) {} };
class Element::Arsenic : public Element {public: Arsenic() :
    Element(33, "arsenic", "As", 74.921602) {} };
class Element::Selenium : public Element {public: Selenium() :
    Element(34, "selenium", "Se", 78.963) {} };
 class Element::Bromine : public Element {public: Bromine() :
    Element(35, "bromine", "Br", 79.9041) {} };
class Element::Krypton : public Element {public: Krypton() :
    Element(36, "krypton", "Kr", 83.7982) {} };
class Element::Rubidium : public Element {public: Rubidium() :
    Element(37, "rubidium", "Rb", 85.46783) {} };
class Element::Strontium : public Element {public: Strontium() :
    Element(38, "strontium", "Sr", 87.621) {} };
class Element::Yttrium : public Element {public: Yttrium() :
    Element(39, "yttrium", "Y", 88.905852) {} };
class Element::Zirconium : public Element {public: Zirconium() :
    Element(40, "zirconium", "Zr", 91.2242) {} };
class Element::Niobium : public Element {public: Niobium() :
    Element(41, "niobium", "Nb", 92.906382) {} };
class Element::Molybdenum : public Element {public: Molybdenum() :
    Element(42, "molybdenum", "Mo", 95.942) {} };
class Element::Technetium : public Element {public: Technetium() :
    Element(43, "technetium", "Tc", 98) {} };
class Element::Ruthenium : public Element {public: Ruthenium() :
    Element(44, "ruthenium", "Ru", 101.072) {} };
class Element::Rhodium : public Element {public: Rhodium() :
    Element(45, "rhodium", "Rh", 102.905502) {} };
class Element::Palladium : public Element {public: Palladium() :
    Element(46, "palladium", "Pd", 106.421) {} };
class Element::Silver : public Element {public: Silver() :
    Element(47, "silver", "Ag", 107.86822) {} };
class Element::Cadmium : public Element {public: Cadmium() :
    Element(48, "cadmium", "Cd", 112.4118) {} };
class Element::Indium : public Element {public: Indium() :
    Element(49, "indium", "In", 114.8183) {} };
class Element::Tin : public Element {public: Tin() :
    Element(50, "tin", "Sn", 118.7107) {} };
class Element::Antimony : public Element {public: Antimony() :
    Element(51, "antimony", "Sb", 121.7601) {} };
class Element::Tellurium : public Element {public: Tellurium() :
    Element(52, "tellurium", "Te", 127.603) {} };
class Element::Iodine : public Element {public: Iodine() :
    Element(53, "iodine", "I", 126.904473) {} };
class Element::Xenon : public Element {public: Xenon() :
    Element(54, "xenon", "Xe", 131.2936) {} };
class Element::Cesium : public Element {public: Cesium() :
    Element(55, "cesium", "Cs", 132.90545192) {} };
class Element::Barium : public Element {public: Barium() :
    Element(56, "barium", "Ba", 137.3277) {} };
class Element::Lanthanum : public Element {public: Lanthanum() :
    Element(57, "lanthanum", "La", 138.905477) {} };
class Element::Cerium : public Element {public: Cerium() :
    Element(58, "cerium", "Ce", 140.1161) {} };
class Element::Praseodymium : public Element {public: Praseodymium() :
    Element(59, "praseodymium", "Pr", 140.907652) {} };
class Element::Neodymium : public Element {public: Neodymium() :
    Element(60, "neodymium", "Nd", 144.2423) {} };
class Element::Promethium : public Element {public: Promethium() :
    Element(61, "promethium", "Pm", 145) {} };
class Element::Samarium : public Element {public: Samarium() :
    Element(62, "samarium", "Sm", 150.362) {} };
class Element::Europium : public Element {public: Europium() :
    Element(63, "europium", "Eu", 151.9641) {} };
class Element::Gadolinium : public Element {public: Gadolinium() :
    Element(64, "gadolinium", "Gd", 157.253) {} };
class Element::Terbium : public Element {public: Terbium() :
    Element(65, "terbium", "Tb", 158.925352) {} };
class Element::Dysprosium : public Element {public: Dysprosium() :
    Element(66, "dysprosium", "Dy", 162.5001) {} };
class Element::Holmium : public Element {public: Holmium() :
    Element(67, "holmium", "Ho", 164.930322) {} };
class Element::Erbium : public Element {public: Erbium() :
    Element(68, "erbium", "Er", 167.2593) {} };
class Element::Thulium : public Element {public: Thulium() :
    Element(69, "thulium", "Tm", 168.934212) {} };
class Element::Ytterbium : public Element {public: Ytterbium() :
    Element(70, "ytterbium", "Yb", 173.043) {} };
class Element::Lutetium : public Element {public: Lutetium() :
    Element(71, "lutetium", "Lu", 174.9671) {} };
class Element::Hafnium : public Element {public: Hafnium() :
    Element(72, "hafnium", "Hf", 178.492) {} };
class Element::Tantalum : public Element {public: Tantalum() :
    Element(73, "tantalum", "Ta", 180.947882) {} };
class Element::Tungsten : public Element {public: Tungsten() :
    Element(74, "tungsten", "W", 183.841) {} };
class Element::Rhenium : public Element {public: Rhenium() :
    Element(75, "rhenium", "Re", 186.2071) {} };
class Element::Osmium : public Element {public: Osmium() :
    Element(76, "osmium", "Os", 190.233) {} };
class Element::Iridium : public Element {public: Iridium() :
    Element(77, "iridium", "Ir", 192.2173) {} };
class Element::Platinum : public Element {public: Platinum() :
    Element(78, "platinum", "Pt", 195.0849) {} };
class Element::Gold : public Element {public: Gold() :
    Element(79, "gold", "Au", 196.9665694) {} };
class Element::Mercury : public Element {public: Mercury() :
    Element(80, "mercury", "Hg", 200.592) {} };
class Element::Thallium : public Element {public: Thallium() :
    Element(81, "thallium", "Tl", 204.38332) {} };
class Element::Lead : public Element {public: Lead() :
    Element(82, "lead", "Pb", 207.21) {} };
class Element::Bismuth : public Element {public: Bismuth() :
    Element(83, "bismuth", "Bi", 208.980401) {} };
class Element::Polonium : public Element {public: Polonium() :
    Element(84, "polonium", "Po", 209) {} };
class Element::Astatine : public Element {public: Astatine() :
    Element(85, "astatine", "At", 210) {} };
class Element::Radon : public Element {public: Radon() :
    Element(86, "radon", "Rn", 222.018) {} };
class Element::Francium : public Element {public: Francium() :
    Element(87, "francium", "Fr", 223) {} };
class Element::Radium : public Element {public: Radium() :
    Element(88, "radium", "Ra", 226) {} };
class Element::Actinium : public Element {public: Actinium() :
    Element(89, "actinium", "Ac", 227) {} };
class Element::Thorium : public Element {public: Thorium() :
    Element(90, "thorium", "Th", 232.038062) {} };
class Element::Protactinium : public Element {public: Protactinium() :
    Element(91, "protactinium", "Pa", 231.035882) {} };
class Element::Uranium : public Element {public: Uranium() :
    Element(92, "uranium", "U", 238.028913) {} };
class Element::Neptunium : public Element {public: Neptunium() :
    Element(93, "neptunium", "Np", 237) {} };
class Element::Plutonium : public Element {public: Plutonium() :
    Element(94, "plutonium", "Pu", 244) {} };
class Element::Americium : public Element {public: Americium() :
    Element(95, "americium", "Am", 243) {} };
class Element::Curium : public Element {public: Curium() :
    Element(96, "curium", "Cm", 247) {} };
class Element::Berkelium : public Element {public: Berkelium() :
    Element(97, "berkelium", "Bk", 247) {} };
class Element::Californium : public Element {public: Californium() :
    Element(98, "californium", "Cf", 251) {} };
class Element::Einsteinium : public Element {public: Einsteinium() :
    Element(99, "einsteinium", "Es", 252) {} };
class Element::Fermium : public Element {public: Fermium() :
    Element(100, "fermium", "Fm", 257) {} };
class Element::Mendelevium : public Element {public: Mendelevium() :
    Element(101, "mendelevium", "Md", 258) {} };
class Element::Nobelium : public Element {public: Nobelium() :
    Element(102, "nobelium", "No", 259) {} };
class Element::Lawrencium : public Element {public: Lawrencium() :
    Element(103, "lawrencium", "Lr", 262) {} };
class Element::Rutherfordium : public Element {public: Rutherfordium() :
    Element(104, "rutherfordium", "Rf", 261) {} };
class Element::Dubnium : public Element {public: Dubnium() :
    Element(105, "dubnium", "Db", 262) {} };
class Element::Seaborgium : public Element {public: Seaborgium() :
    Element(106, "seaborgium", "Sg", 266) {} };
class Element::Bohrium : public Element {public: Bohrium() :
    Element(107, "bohrium", "Bh", 264) {} };
class Element::Hassium : public Element {public: Hassium() :
    Element(108, "hassium", "Hs", 269) {} };
class Element::Meitnerium : public Element {public: Meitnerium() :
    Element(109, "meitnerium", "Mt", 268) {} };
class Element::Darmstadtium : public Element {public: Darmstadtium() :
    Element(110, "darmstadtium", "Ds", 281) {} };
class Element::Roentgenium : public Element {public: Roentgenium() :
    Element(111, "roentgenium", "Rg", 272) {} };
class Element::Ununbium : public Element {public: Ununbium() :
    Element(112, "ununbium", "Uub", 285) {} };
class Element::Ununtrium : public Element {public: Ununtrium() :
    Element(113, "ununtrium", "Uut", 284) {} };
class Element::Ununquadium : public Element {public: Ununquadium() :
    Element(114, "ununquadium", "Uuq", 289) {} };
class Element::Ununpentium : public Element {public: Ununpentium() :
    Element(115, "ununpentium", "Uup", 288) {} };
class Element::Ununhexium : public Element {public: Ununhexium() :
    Element(116, "ununhexium", "Uuh", 292) {} };


} // namespace SimTK

#endif // SimTK_SIMBODY_ELEMENT_H_
