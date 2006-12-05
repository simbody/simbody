/* Protein.h */

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
 
#ifndef _SIMTK_CHEM_PROTEIN_H_
#define _SIMTK_CHEM_PROTEIN_H_

#include "Molecule.h"
#include "AminoAcid.h"

class Protein : public Molecule
{
public:
	Protein(const char* sequence);
	~Protein();
	
	Protein(const Protein & src);
	Protein & operator=(const Protein & src);
	bool operator==(const Protein & src) const;
	bool operator!=(const Protein & src) const;

	const char * sequence() const;
	const AminoAcid & residue(int resNum) const;
	
private:
	class ProteinRep * rep;
};

#endif /*PROTEIN_H_*/
