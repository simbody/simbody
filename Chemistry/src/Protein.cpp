/* Protein.cpp */

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
 
#include "Protein.h"
#include "AminoAcid.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

// Library side implementation of Protein
class ProteinRep
{
	friend class Protein;
private:
	ProteinRep(const char * sequence) {
		for (const char * i = sequence; *i; i++) {
			const char olc = *i;
			residues.push_back(AminoAcid(olc));
		}
	}

	// Warning - the sequence will probably go out of scope when the protein does.
	const char * sequence() const {
		// TODO - a static string applies to all Proteins.
		static string answer;

		answer.clear();
		vector<AminoAcid>::const_iterator i = residues.begin();
		while (i != residues.end()) {
			answer += (*i).oneLetterCode();
			i++;
		}

		return answer.c_str();
	}
	
	vector<AminoAcid> residues;

	Protein * myHandle;
};


Protein::Protein(const char * sequence)
{
	rep = new ProteinRep(sequence);
	rep->myHandle = this;	
}

Protein::~Protein()
{
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

const char * Protein::sequence() const {
	return rep->sequence();
}
