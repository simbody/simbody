/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
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

// This is the one file where we want the handle instantiated
#define DO_INSTANTIATE_ELEMENT_PIMPL_HANDLE
#include "simbody/internal/Element.h"
#undef DO_INSTANTIATE_ELEMENT_PIMPL_HANDLE

#include "ElementRep.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#include <map>

using namespace std;


namespace SimTK {

template class PIMPLImplementation<Element,ElementRep>; // explicit instantiation


///////////////
/// Element ///
///////////////

template class PIMPLHandle<Element,ElementRep>; // explicit instantiation

// TODO - modify to return references to permanent element instances
/* static */ Element Element::getByAtomicNumber(int atomicNumber) {
    switch (atomicNumber) {
        case 1: return Hydrogen();
        case 2: return Helium();
        case 3: return Lithium();
        case 4: return Beryllium();
        case 5: return Boron();
        case 6: return Carbon();
        case 7: return Nitrogen();
        case 8: return Oxygen();
        case 9: return Fluorine();
        case 10: return Neon();
        case 11: return Sodium();
        case 12: return Magnesium();
        case 13: return Aluminum();
        case 14: return Silicon();
        case 15: return Phosphorus();
        case 16: return Sulfur();
        case 17: return Chlorine();
        case 18: return Argon();
        case 19: return Potassium();
        case 20: return Calcium();
        case 21: return Scandium();
        case 22: return Titanium();
        case 23: return Vanadium();
        case 24: return Chromium();
        case 25: return Manganese();
        case 26: return Iron();
        case 27: return Cobalt();
        case 28: return Nickel();
        case 29: return Copper();
        case 30: return Zinc();
        case 31: return Gallium();
        case 32: return Germanium();
        case 33: return Arsenic();
        case 34: return Selenium();
        case 35: return Bromine();
        case 36: return Krypton();
        case 37: return Rubidium();
        case 38: return Strontium();
        case 39: return Yttrium();
        case 40: return Zirconium();
        case 41: return Niobium();
        case 42: return Molybdenum();
        case 43: return Technetium();
        case 44: return Ruthenium();
        case 45: return Rhodium();
        case 46: return Palladium();
        case 47: return Silver();
        case 48: return Cadmium();
        case 49: return Indium();
        case 50: return Tin();
        case 51: return Antimony();
        case 52: return Tellurium();
        case 53: return Iodine();
        case 54: return Xenon();
        case 55: return Cesium();
        case 56: return Barium();
        case 57: return Lanthanum();
        case 58: return Cerium();
        case 59: return Praseodymium();
        case 60: return Neodymium();
        case 61: return Promethium();
        case 62: return Samarium();
        case 63: return Europium();
        case 64: return Gadolinium();
        case 65: return Terbium();
        case 66: return Dysprosium();
        case 67: return Holmium();
        case 68: return Erbium();
        case 69: return Thulium();
        case 70: return Ytterbium();
        case 71: return Lutetium();
        case 72: return Hafnium();
        case 73: return Tantalum();
        case 74: return Tungsten();
        case 75: return Rhenium();
        case 76: return Osmium();
        case 77: return Iridium();
        case 78: return Platinum();
        case 79: return Gold();
        case 80: return Mercury();
        case 81: return Thallium();
        case 82: return Lead();
        case 83: return Bismuth();
        case 84: return Polonium();
        case 85: return Astatine();
        case 86: return Radon();
        case 87: return Francium();
        case 88: return Radium();
        case 89: return Actinium();
        case 90: return Thorium();
        case 91: return Protactinium();
        case 92: return Uranium();
        case 93: return Neptunium();
        case 94: return Plutonium();
        case 95: return Americium();
        case 96: return Curium();
        case 97: return Berkelium();
        case 98: return Californium();
        case 99: return Einsteinium();
        case 100: return Fermium();
        case 101: return Mendelevium();
        case 102: return Nobelium();
        case 103: return Lawrencium();
        case 104: return Rutherfordium();
        case 105: return Dubnium();
        case 106: return Seaborgium();
        case 107: return Bohrium();
        case 108: return Hassium();
        case 109: return Meitnerium();
        case 110: return Darmstadtium();
        case 111: return Roentgenium();
        case 112: return Ununbium();
        case 113: return Ununtrium();
        case 114: return Ununquadium();
        case 115: return Ununpentium();
        case 116: return Ununhexium();

        default: assert(false); return Hydrogen();
    }
}

// TODO - modify to return references to permanent element instances
/* static */ Element Element::getBySymbol(const SimTK::String& symbol) 
{
    if      (symbol == "H")  return Hydrogen();
    else if (symbol == "D")  return Deuterium();
    else if (symbol == "He") return Helium();
    else if (symbol == "Li") return Lithium();
    else if (symbol == "Be") return Beryllium();
    else if (symbol == "B")  return Boron();
    else if (symbol == "C")  return  Carbon();
    else if (symbol == "N")  return  Nitrogen();
    else if (symbol == "O")  return  Oxygen();
    else if (symbol == "F")  return  Fluorine();
    else if (symbol == "Ne") return Neon();
    else if (symbol == "Na") return Sodium();
    else if (symbol == "Mg") return Magnesium();
    else if (symbol == "Al") return Aluminum();
    else if (symbol == "Si") return Silicon();
    else if (symbol == "P") return  Phosphorus();
    else if (symbol == "S") return  Sulfur();
    else if (symbol == "Cl") return Chlorine();
    else if (symbol == "Ar") return Argon();
    else if (symbol == "K") return  Potassium();
    else if (symbol == "Ca") return Calcium();
    else if (symbol == "Sc") return Scandium();
    else if (symbol == "Ti") return Titanium();
    else if (symbol == "V") return  Vanadium();
    else if (symbol == "Cr") return Chromium();
    else if (symbol == "Mn") return Manganese();
    else if (symbol == "Fe") return Iron();
    else if (symbol == "Co") return Cobalt();
    else if (symbol == "Ni") return Nickel();
    else if (symbol == "Cu") return Copper();
    else if (symbol == "Zn") return Zinc();
    else if (symbol == "Ga") return Gallium();
    else if (symbol == "Ge") return Germanium();
    else if (symbol == "As") return Arsenic();
    else if (symbol == "Se") return Selenium();
    else if (symbol == "Br") return Bromine();
    else if (symbol == "Kr") return Krypton();
    else if (symbol == "Rb") return Rubidium();
    else if (symbol == "Sr") return Strontium();
    else if (symbol == "Y")  return Yttrium();
    else if (symbol == "Zr") return Zirconium();
    else if (symbol == "Nb") return Niobium();
    else if (symbol == "Mo") return Molybdenum();
    else if (symbol == "Tc") return Technetium();
    else if (symbol == "Ru") return Ruthenium();
    else if (symbol == "Rh") return Rhodium();
    else if (symbol == "Pd") return Palladium();
    else if (symbol == "Ag") return Silver();
    else if (symbol == "Cd") return Cadmium();
    else if (symbol == "In") return Indium();
    else if (symbol == "Sn") return Tin();
    else if (symbol == "Sb") return Antimony();
    else if (symbol == "Te") return Tellurium();
    else if (symbol == "I")  return Iodine();
    else if (symbol == "Xe") return Xenon();
    else if (symbol == "Cs") return Cesium();
    else if (symbol == "Ba") return Barium();
    else if (symbol == "La") return Lanthanum();
    else if (symbol == "Ce") return Cerium();
    else if (symbol == "Pr") return Praseodymium();
    else if (symbol == "Nd") return Neodymium();
    else if (symbol == "Pm") return Promethium();
    else if (symbol == "Sm") return Samarium();
    else if (symbol == "Eu") return Europium();
    else if (symbol == "Gd") return Gadolinium();
    else if (symbol == "Tb") return Terbium();
    else if (symbol == "Dy") return Dysprosium();
    else if (symbol == "Ho") return Holmium();
    else if (symbol == "Er") return Erbium();
    else if (symbol == "Tm") return Thulium();
    else if (symbol == "Yb") return Ytterbium();
    else if (symbol == "Lu") return Lutetium();
    else if (symbol == "Hf") return Hafnium();
    else if (symbol == "Ta") return Tantalum();
    else if (symbol == "W")  return Tungsten();
    else if (symbol == "Re") return Rhenium();
    else if (symbol == "Os") return Osmium();
    else if (symbol == "Ir") return Iridium();
    else if (symbol == "Pt") return Platinum();
    else if (symbol == "Au") return Gold();
    else if (symbol == "Hg") return Mercury();
    else if (symbol == "Tl") return Thallium();
    else if (symbol == "Pb") return Lead();
    else if (symbol == "Bi") return Bismuth();
    else if (symbol == "Po") return Polonium();
    else if (symbol == "At") return Astatine();
    else if (symbol == "Rn") return Radon();
    else if (symbol == "Fr") return Francium();
    else if (symbol == "Ra") return Radium();
    else if (symbol == "Ac") return Actinium();
    else if (symbol == "Th") return Thorium();
    else if (symbol == "Pa") return Protactinium();
    else if (symbol == "U")  return Uranium();
    else if (symbol == "Np") return Neptunium();
    else if (symbol == "Pu") return Plutonium();
    else if (symbol == "Am") return Americium();
    else if (symbol == "Cm") return Curium();
    else if (symbol == "Bk") return Berkelium();
    else if (symbol == "Cf") return Californium();
    else if (symbol == "Es") return Einsteinium();
    else if (symbol == "Fm") return Fermium();
    else if (symbol == "Md") return Mendelevium();
    else if (symbol == "No") return Nobelium();
    else if (symbol == "Lr") return Lawrencium();
    else if (symbol == "Rf") return Rutherfordium();
    else if (symbol == "Db") return Dubnium();
    else if (symbol == "Sg") return Seaborgium();
    else if (symbol == "Bh") return Bohrium();
    else if (symbol == "Hs") return Hassium();
    else if (symbol == "Mt") return Meitnerium();
    else if (symbol == "Ds") return Darmstadtium();
    else if (symbol == "Rg") return Roentgenium();
    else if (symbol == "Uub") return Ununbium();
    else if (symbol == "Uut") return Ununtrium();
    else if (symbol == "Uuq") return Ununquadium();
    else if (symbol == "Uup") return Ununpentium();
    else if (symbol == "Uuh") return Ununhexium();

    else {
        assert(false); 
        return Hydrogen();
    }
}

std::ostream& operator<<(std::ostream& o, const Element& e) {
    return o << e.getName();
}

Element::Element() 
    : HandleBase(new ElementRep())
{}
    
Element::Element(int atomicNumber, Name name, Symbol symbol, Mass typicalMass) 
    : HandleBase(new ElementRep(atomicNumber, name, symbol, typicalMass))
{}

Element::Symbol Element::getSymbol() const {
    return getImpl().getSymbol();
}
Element::Symbol Element::getName() const {
    return getImpl().getName();
}

//const Element& Element::Hydrogen() {
//    if (elementsByAtomicNumber.find(1) == elementsByAtomicNumber.end())
//        defineCanonicalElement(1, "hydrogen", "H", 1.007947);
//    return elementsByAtomicNumber.find(1)->second;
//}
//const Element& Element::Helium() {
//    if (elementsByAtomicNumber.find(2) == elementsByAtomicNumber.end())
//        defineCanonicalElement(2, "helium", "He", 4.003);
//    return elementsByAtomicNumber.find(2)->second;
//}
//const Element& Element::Carbon() {
//    if (elementsByAtomicNumber.find(6) == elementsByAtomicNumber.end())
//        defineCanonicalElement(6, "carbon", "C", 12.01078);
//    return elementsByAtomicNumber.find(6)->second;
//}
//const Element& Element::Nitrogen() {
//    if (elementsByAtomicNumber.find(7) == elementsByAtomicNumber.end())
//        defineCanonicalElement(7, "nitrogen", "N", 14.00672);
//    return elementsByAtomicNumber.find(7)->second;
//}
//const Element& Element::Oxygen() {
//    if (elementsByAtomicNumber.find(8) == elementsByAtomicNumber.end())
//        defineCanonicalElement(8, "oxygen", "O", 15.99943);
//    return elementsByAtomicNumber.find(8)->second;
//}
//const Element& Element::Neon() {
//    if (elementsByAtomicNumber.find(10) == elementsByAtomicNumber.end())
//        defineCanonicalElement(10, "neon", "Ne", 20.180);
//    return elementsByAtomicNumber.find(10)->second;
//}
//const Element& Element::Phosphorus() {
//    if (elementsByAtomicNumber.find(15) == elementsByAtomicNumber.end())
//        defineCanonicalElement(15, "phosphorus", "P", 30.9737622);
//    return elementsByAtomicNumber.find(15)->second;
//}
//const Element& Element::Sulfur() {
//    if (elementsByAtomicNumber.find(16) == elementsByAtomicNumber.end())
//        defineCanonicalElement(16, "sulfur", "S", 32.0655);
//    return elementsByAtomicNumber.find(16)->second;
//}
//const Element& Element::Argon() {
//    if (elementsByAtomicNumber.find(18) == elementsByAtomicNumber.end())
//        defineCanonicalElement(18, "argon", "Ar", 39.9481);
//    return elementsByAtomicNumber.find(18)->second;
//}
//const Element& Element::Selenium() {
//    if (elementsByAtomicNumber.find(34) == elementsByAtomicNumber.end())
//        defineCanonicalElement(34, "selenium", "Se", 78.963);
//    return elementsByAtomicNumber.find(34)->second;
//}
//const Element& Element::Xenon() {
//    if (elementsByAtomicNumber.find(54) == elementsByAtomicNumber.end())
//        defineCanonicalElement(54, "xenon", "Xe", 131.290);
//    return elementsByAtomicNumber.find(54)->second;
//}
//const Element& Element::Radon() {
//    if (elementsByAtomicNumber.find(86) == elementsByAtomicNumber.end())
//        defineCanonicalElement(86, "radon", "Rn", 222.018);
//    return elementsByAtomicNumber.find(86)->second;
//}

} // namespace SimTK
