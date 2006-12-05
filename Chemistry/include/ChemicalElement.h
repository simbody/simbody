/* ChemicalElement.h */

/* Portions copyright (c) 2006 Stanford University and Christopher M. Bruns.
 * Contributors: Michael Sherman
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
 
#ifndef CHEMICALELEMENT_H_
#define CHEMICALELEMENT_H_

// ChemicalElement is the client side public interface to the library implementation ChemicalElementRep
// To minimize client/library API problems, avoid the following in the public interface:
//  1) virtual functions (functions might still be actually virtual in the library implementation)
//  2) std:: classes - the implementation is public, and so can break the API
//  3) private data

// Generally one chemical element for each entry in the periodic table of elements.
// e.g. "iron" "carbon"
// One isotope, deuterium, is included, because its mass is so different from the primary isotope
class ChemicalElement
{
public:
	ChemicalElement(int number, const char* symbol, const char* name, double defaultMass);
	~ChemicalElement();

	ChemicalElement(const ChemicalElement & src);
	ChemicalElement & operator=(const ChemicalElement & src);
	bool operator==(const ChemicalElement & src) const;
	bool operator!=(const ChemicalElement & src) const;

	// Atomic number, e.g. hydrogen = 1, uranium = 92
	int number() const;

	// One or two character element symbol
	// First character is always capitalized
	// e.g. "Cl"
	const char* symbol() const;

	// Element name, NOT capitalized, as per IUPAC guidelines
	// e.g. "hydrogen"
	const char* name() const;
	
	// Atomic mass, might be overridden for isotopes
	// In atomic mass units (Daltons) (Da, g/mol, amu, u)
	double mass() const;

    static const ChemicalElement Hydrogen;
    static const ChemicalElement Deuterium;
    static const ChemicalElement Helium;
    static const ChemicalElement Lithium;
    static const ChemicalElement Beryllium;
    static const ChemicalElement Boron;
    static const ChemicalElement Carbon;
    static const ChemicalElement Nitrogen;
    static const ChemicalElement Oxygen;
    static const ChemicalElement Fluorine;
    static const ChemicalElement Neon;
    static const ChemicalElement Sodium;
    static const ChemicalElement Magnesium;
    static const ChemicalElement Aluminum;
    static const ChemicalElement Silicon;
    static const ChemicalElement Phosphorus;
    static const ChemicalElement Sulfur;
    static const ChemicalElement Chlorine;
    static const ChemicalElement Argon;
    static const ChemicalElement Potassium;
    static const ChemicalElement Calcium;
    static const ChemicalElement Scandium;
    static const ChemicalElement Titanium;
    static const ChemicalElement Vanadium;
    static const ChemicalElement Chromium;
    static const ChemicalElement Manganese;
    static const ChemicalElement Iron;
    static const ChemicalElement Cobalt;
    static const ChemicalElement Nickel;
    static const ChemicalElement Copper;
    static const ChemicalElement Zinc;
    static const ChemicalElement Gallium;
    static const ChemicalElement Germanium;
    static const ChemicalElement Arsenic;
    static const ChemicalElement Selenium;
    static const ChemicalElement Bromine;
    static const ChemicalElement Krypton;
    static const ChemicalElement Rubidium;
    static const ChemicalElement Strontium;
    static const ChemicalElement Yttrium;
    static const ChemicalElement Zirconium;
    static const ChemicalElement Niobium;
    static const ChemicalElement Molybdenum;
    static const ChemicalElement Technetium;
    static const ChemicalElement Ruthenium;
    static const ChemicalElement Rhodium;
    static const ChemicalElement Palladium;
    static const ChemicalElement Silver;
    static const ChemicalElement Cadmium;
    static const ChemicalElement Indium;
    static const ChemicalElement Tin;
    static const ChemicalElement Antimony;
    static const ChemicalElement Tellurium;
    static const ChemicalElement Iodine;
    static const ChemicalElement Xenon;
    static const ChemicalElement Cesium;
    static const ChemicalElement Barium;
    static const ChemicalElement Lanthanum;
    static const ChemicalElement Cerium;
    static const ChemicalElement Praseodymium;
    static const ChemicalElement Neodymium;
    static const ChemicalElement Promethium;
    static const ChemicalElement Samarium;
    static const ChemicalElement Europium;
    static const ChemicalElement Gadolinium;
    static const ChemicalElement Terbium;
    static const ChemicalElement Dysprosium;
    static const ChemicalElement Holmium;
    static const ChemicalElement Erbium;
    static const ChemicalElement Thulium;
    static const ChemicalElement Ytterbium;
    static const ChemicalElement Lutetium;
    static const ChemicalElement Hafnium;
    static const ChemicalElement Tantalum;
    static const ChemicalElement Tungsten;
    static const ChemicalElement Rhenium;
    static const ChemicalElement Osmium;
    static const ChemicalElement Iridium;
    static const ChemicalElement Platinum;
    static const ChemicalElement Gold;
    static const ChemicalElement Mercury;
    static const ChemicalElement Thallium;
    static const ChemicalElement Lead;
    static const ChemicalElement Bismuth;
    static const ChemicalElement Polonium;
    static const ChemicalElement Astatine;
    static const ChemicalElement Radon;
    static const ChemicalElement Francium;
    static const ChemicalElement Radium;
    static const ChemicalElement Actinium;
    static const ChemicalElement Thorium;
    static const ChemicalElement Protactinium;
    static const ChemicalElement Uranium;
    static const ChemicalElement Neptunium;
    static const ChemicalElement Plutonium;
    static const ChemicalElement Americium;
    static const ChemicalElement Curium;
    static const ChemicalElement Berkelium;
    static const ChemicalElement Californium;
    static const ChemicalElement Einsteinium;
    static const ChemicalElement Fermium;
    static const ChemicalElement Mendelevium;
    static const ChemicalElement Nobelium;
    static const ChemicalElement Lawrencium;
    static const ChemicalElement Rutherfordium;
    static const ChemicalElement Dubnium;
    static const ChemicalElement Seaborgium;
    static const ChemicalElement Bohrium;
    static const ChemicalElement Hassium;
    static const ChemicalElement Meitnerium;
    static const ChemicalElement Darmstadtium;
    static const ChemicalElement Unununium;
    static const ChemicalElement Ununbium;

private:
	class ChemicalElementRep * rep; // Handle for library implementation
};

#endif /*CHEMICALELEMENT_H_*/
