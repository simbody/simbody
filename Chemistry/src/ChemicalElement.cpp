/* ChemicalElement.cpp */

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

#include "ChemicalElement.h"
#include <string>

using namespace std;

// Library side implementation of ChemicalElement
class ChemicalElementRep
{
	friend class ChemicalElement;
private:
	ChemicalElementRep(int number, const char* symbol, const char* name, double mass, const ChemicalElement & handle)
		: number(number), symbol(symbol), name(name), defaultMass(mass), myHandle(& handle)
	{}

	int number;
	string symbol;
	string name;
	double defaultMass;
	
	const ChemicalElement * myHandle;
};

ChemicalElement::ChemicalElement(const int number, const char* symbol, const char* name, const double defaultMass)
{
	rep = new ChemicalElementRep(number, symbol, name, defaultMass, *this);
}

ChemicalElement::~ChemicalElement() {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

int ChemicalElement::number() const {
	return rep->number;
}

// One or two character element symbol
// First character is always capitalized
// e.g. "Cl"
const char * ChemicalElement::symbol() const {
	return rep->symbol.c_str();
}

// Element name, NOT capitalized, as per IUPAC
// e.g. "hydrogen"
const char * ChemicalElement::name() const {
	return rep->name.c_str();
}

double ChemicalElement::getMass() const {
	return rep->defaultMass;
}

const ChemicalElement ChemicalElement::Hydrogen    (1,  "H",  "hydrogen", 1.007947);
const ChemicalElement ChemicalElement::Deuterium   (1,  "D",  "deuterium", 2.01355321270);
const ChemicalElement ChemicalElement::Helium      (2,  "He", "helium", 4.003);
const ChemicalElement ChemicalElement::Lithium     (3,  "Li", "lithium", 6.941);
const ChemicalElement ChemicalElement::Beryllium   (4,  "Be", "beryllium", 9.012);
const ChemicalElement ChemicalElement::Boron       (5,  "B",  "boron", 10.811);
const ChemicalElement ChemicalElement::Carbon      (6,  "C",  "carbon", 12.01078);
const ChemicalElement ChemicalElement::Nitrogen    (7,  "N",  "nitrogen", 14.00672);
const ChemicalElement ChemicalElement::Oxygen      (8,  "O",  "oxygen", 15.99943);
const ChemicalElement ChemicalElement::Fluorine    (9,  "F",  "fluorine", 18.998);
const ChemicalElement ChemicalElement::Neon        (10, "Ne", "neon", 20.180);
const ChemicalElement ChemicalElement::Sodium      (11, "Na", "sodium", 22.989769282);
const ChemicalElement ChemicalElement::Magnesium   (12, "Mg", "magnesium", 24.30506);
const ChemicalElement ChemicalElement::Aluminum    (13, "Al", "aluminum", 26.982);
const ChemicalElement ChemicalElement::Silicon     (14, "Si", "silicon", 28.086);
const ChemicalElement ChemicalElement::Phosphorus  (15, "P", "phosphorus", 30.9737622);
const ChemicalElement ChemicalElement::Sulfur      (16, "S", "sulfur", 32.0655);
const ChemicalElement ChemicalElement::Chlorine    (17, "Cl", "chlorine", 35.4532);
const ChemicalElement ChemicalElement::Argon       (18, "Ar", "argon", 39.948);
const ChemicalElement ChemicalElement::Potassium   (19, "K", "potassium", 39.09831);
const ChemicalElement ChemicalElement::Calcium     (20, "Ca", "calcium", 40.0784);
const ChemicalElement ChemicalElement::Scandium    (21, "Sc", "scandium", 44.956);
const ChemicalElement ChemicalElement::Titanium    (22, "Ti", "titanium", 47.88);
const ChemicalElement ChemicalElement::Vanadium    (23, "V", "vanadium", 50.942);
const ChemicalElement ChemicalElement::Chromium    (24, "Cr", "chromium", 51.996);
const ChemicalElement ChemicalElement::Manganese   (25, "Mn", "manganese", 54.9380455);
const ChemicalElement ChemicalElement::Iron        (26, "Fe", "iron", 55.8452);
const ChemicalElement ChemicalElement::Cobalt      (27, "Co", "cobalt", 58.9331955);
const ChemicalElement ChemicalElement::Nickel      (28, "Ni", "nickel", 58.69342);
const ChemicalElement ChemicalElement::Copper      (29, "Cu", "copper", 63.5463);
const ChemicalElement ChemicalElement::Zinc        (30, "Zn", "zinc", 65.4094);
const ChemicalElement ChemicalElement::Gallium     (31, "Ga", "gallium", 69.723);
const ChemicalElement ChemicalElement::Germanium   (32, "Ge", "germanium", 72.61);
const ChemicalElement ChemicalElement::Arsenic     (33, "As", "arsenic", 74.922);
const ChemicalElement ChemicalElement::Selenium    (34, "Se", "selenium", 78.963);
const ChemicalElement ChemicalElement::Bromine     (35, "Br", "bromine", 79.9041);
const ChemicalElement ChemicalElement::Krypton     (36, "Kr", "krypton", 83.80);
const ChemicalElement ChemicalElement::Rubidium    (37, "Rb", "rubidium", 85.468);
const ChemicalElement ChemicalElement::Strontium   (38, "Sr", "strontium", 87.62);
const ChemicalElement ChemicalElement::Yttrium     (39, "Y", "yttrium", 88.906);
const ChemicalElement ChemicalElement::Zirconium   (40, "Zr", "zirconium", 91.224);
const ChemicalElement ChemicalElement::Niobium     (41, "Nb", "niobium", 92.906);
const ChemicalElement ChemicalElement::Molybdenum  (42, "Mo", "molybdenum", 95.94);
const ChemicalElement ChemicalElement::Technetium  (43, "Tc", "technetium", 97.907);
const ChemicalElement ChemicalElement::Ruthenium   (44, "Ru", "ruthenium", 101.07);
const ChemicalElement ChemicalElement::Rhodium     (45, "Rh", "rhodium", 102.906);
const ChemicalElement ChemicalElement::Palladium   (46, "Pd", "palladium", 106.42);
const ChemicalElement ChemicalElement::Silver      (47, "Ag", "silver", 107.868);
const ChemicalElement ChemicalElement::Cadmium     (48, "Cd", "cadmium", 112.411);
const ChemicalElement ChemicalElement::Indium      (49, "In", "indium", 114.82);
const ChemicalElement ChemicalElement::Tin         (50, "Sn", "tin", 118.710);
const ChemicalElement ChemicalElement::Antimony    (51, "Sb", "antimony", 121.757);
const ChemicalElement ChemicalElement::Tellurium   (52, "Te", "tellurium", 127.60);
const ChemicalElement ChemicalElement::Iodine      (53, "I", "iodine", 126.904);
const ChemicalElement ChemicalElement::Xenon       (54, "Xe", "xenon", 131.290);
const ChemicalElement ChemicalElement::Cesium      (55, "Cs", "cesium", 132.905);
const ChemicalElement ChemicalElement::Barium      (56, "Ba", "barium", 137.327);
const ChemicalElement ChemicalElement::Lanthanum   (57, "La", "lanthanum", 138.906);
const ChemicalElement ChemicalElement::Cerium      (58, "Ce", "cerium", 140.115);
const ChemicalElement ChemicalElement::Praseodymium (59, "Pr", "praseodymium", 140.908);
const ChemicalElement ChemicalElement::Neodymium   (60, "Nd", "neodymium", 144.24);
const ChemicalElement ChemicalElement::Promethium  (61, "Pm", "promethium", 144.913);
const ChemicalElement ChemicalElement::Samarium    (62, "Sm", "samarium", 150.36);
const ChemicalElement ChemicalElement::Europium    (63, "Eu", "europium", 151.965);
const ChemicalElement ChemicalElement::Gadolinium  (64, "Gd", "gadolinium", 157.25);
const ChemicalElement ChemicalElement::Terbium     (65, "Tb", "terbium", 158.925);
const ChemicalElement ChemicalElement::Dysprosium  (66, "Dy", "dysprosium", 162.50);
const ChemicalElement ChemicalElement::Holmium     (67, "Ho", "holmium", 164.930);
const ChemicalElement ChemicalElement::Erbium      (68, "Er", "erbium", 167.26);
const ChemicalElement ChemicalElement::Thulium     (69, "Tm", "thulium", 168.934);
const ChemicalElement ChemicalElement::Ytterbium   (70, "Yb", "ytterbium", 173.04);
const ChemicalElement ChemicalElement::Lutetium    (71, "Lu", "lutetium", 174.967);
const ChemicalElement ChemicalElement::Hafnium     (72, "Hf", "hafnium", 178.49);
const ChemicalElement ChemicalElement::Tantalum    (73, "Ta", "tantalum", 180.948);
const ChemicalElement ChemicalElement::Tungsten    (74, "W", "tungsten", 183.84);
const ChemicalElement ChemicalElement::Rhenium     (75, "Re", "rhenium", 186.207);
const ChemicalElement ChemicalElement::Osmium      (76, "Os", "osmium", 190.2);
const ChemicalElement ChemicalElement::Iridium     (77, "Ir", "iridium", 192.22);
const ChemicalElement ChemicalElement::Platinum    (78, "Pt", "platinum", 195.08);
const ChemicalElement ChemicalElement::Gold        (79, "Au", "gold", 196.967);
const ChemicalElement ChemicalElement::Mercury     (80, "Hg", "mercury", 200.59);
const ChemicalElement ChemicalElement::Thallium    (81, "Tl", "thallium", 204.383);
const ChemicalElement ChemicalElement::Lead        (82, "Pb", "lead", 207.2);
const ChemicalElement ChemicalElement::Bismuth     (83, "Bi", "bismuth", 208.980);
const ChemicalElement ChemicalElement::Polonium    (84, "Po", "polonium", 208.982);
const ChemicalElement ChemicalElement::Astatine    (85, "At", "astatine", 209.978);
const ChemicalElement ChemicalElement::Radon       (86, "Rn", "radon", 222.018);
const ChemicalElement ChemicalElement::Francium    (87, "Fr", "francium", 223.020);
const ChemicalElement ChemicalElement::Radium      (88, "Ra", "radium", 226.025);
const ChemicalElement ChemicalElement::Actinium    (89, "Ac", "actinium", 227.028);
const ChemicalElement ChemicalElement::Thorium     (90, "Th", "thorium", 232.038);
const ChemicalElement ChemicalElement::Protactinium (91, "Pa", "protactinium", 231.038);
const ChemicalElement ChemicalElement::Uranium     (92, "U",  "uranium", 238.028913);
const ChemicalElement ChemicalElement::Neptunium   (93, "Np", "neptunium", 237.048);
const ChemicalElement ChemicalElement::Plutonium   (94, "Pu", "plutonium", 244.064);
const ChemicalElement ChemicalElement::Americium   (95, "Am", "americium", 243.061);
const ChemicalElement ChemicalElement::Curium      (96, "Cm", "curium", 247.070);
const ChemicalElement ChemicalElement::Berkelium   (97, "Bk", "berkelium", 247.070);
const ChemicalElement ChemicalElement::Californium (98, "Cf", "californium", 251.080);
const ChemicalElement ChemicalElement::Einsteinium (99, "Es", "einsteinium", 252.083);
const ChemicalElement ChemicalElement::Fermium     (100, "Fm", "fermium", 257.095);
const ChemicalElement ChemicalElement::Mendelevium (101, "Md", "mendelevium", 258.099);
const ChemicalElement ChemicalElement::Nobelium    (102, "No", "nobelium", 259.101);
const ChemicalElement ChemicalElement::Lawrencium  (103, "Lr", "lawrencium", 260.105);
const ChemicalElement ChemicalElement::Rutherfordium (104, "Rf", "rutherfordium", 261);
const ChemicalElement ChemicalElement::Dubnium     (105, "Db", "dubnium", 262);
const ChemicalElement ChemicalElement::Seaborgium  (106, "Sg", "seaborgium", 263);
const ChemicalElement ChemicalElement::Bohrium     (107, "Bh", "bohrium", 262);
const ChemicalElement ChemicalElement::Hassium     (108, "Hs", "hassium", 265);
const ChemicalElement ChemicalElement::Meitnerium  (109, "Mt", "meitnerium", 266);
const ChemicalElement ChemicalElement::Darmstadtium (110, "Ds", "darmstadtium", 281);
const ChemicalElement ChemicalElement::Unununium   (111, "Uuu", "unununium", 272);
const ChemicalElement ChemicalElement::Ununbium    (112, "Uub", "ununbium", 285);
