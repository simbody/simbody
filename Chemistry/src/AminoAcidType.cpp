/* AminoAcidType.cpp */

#include "AminoAcidType.h"
#include <string>
#include <map>

using namespace std;

typedef map<const char, const AminoAcidType *> CharAcid;

static CharAcid olcTypes;

const AminoAcidType AminoAcidType::Alanine      ('A', "Ala", "alanine");
const AminoAcidType AminoAcidType::Cysteine     ('C', "Cys", "cysteine");
const AminoAcidType AminoAcidType::Aspartate    ('D', "Asp", "aspartate");
const AminoAcidType AminoAcidType::Glutamate    ('E', "Glu", "glutamate");
const AminoAcidType AminoAcidType::Phenylalanine('F', "Phe", "phenylalanine");
const AminoAcidType AminoAcidType::Glycine      ('G', "Gly", "glycine");
const AminoAcidType AminoAcidType::Histidine    ('H', "His", "histidine");
const AminoAcidType AminoAcidType::Isoleucine   ('I', "Ile", "isoleucine");
const AminoAcidType AminoAcidType::Lysine       ('K', "Lys", "lysine");
const AminoAcidType AminoAcidType::Leucine      ('L', "Leu", "leucine");
const AminoAcidType AminoAcidType::Methionine   ('M', "Met", "methionine");
const AminoAcidType AminoAcidType::Asparagine   ('N', "Asn", "asparagine");
const AminoAcidType AminoAcidType::Proline      ('P', "Pro", "proline");
const AminoAcidType AminoAcidType::Glutamine    ('Q', "Gln", "glutamine");
const AminoAcidType AminoAcidType::Arginine     ('R', "Arg", "arginine");
const AminoAcidType AminoAcidType::Serine       ('S', "Ser", "serine");
const AminoAcidType AminoAcidType::Threonine    ('T', "Thr", "threonine");
const AminoAcidType AminoAcidType::Valine       ('V', "Val", "valine");
const AminoAcidType AminoAcidType::Tryptophan   ('W', "Trp", "tryptophan");
const AminoAcidType AminoAcidType::Tyrosine     ('Y', "Tyr", "tyrosine");

class AminoAcidTypeRep {
	friend class AminoAcidType;

private:
	AminoAcidTypeRep(char olc, const char* tlc, const char* n, const AminoAcidType * handle)
		: oneLetterCode(olc), threeLetterCode(tlc), name(n)
	{
		// Remember this character/type association
		olcTypes[olc] = handle;
		myHandle = handle;
	}

	static const AminoAcidType * typeFor(char oneLetterCode) {
		if (olcTypes.find(oneLetterCode) == olcTypes.end()) return NULL;
		else return olcTypes[oneLetterCode];
	}

	bool operator!=(const AminoAcidTypeRep & other) const {
		if (oneLetterCode != other.oneLetterCode) return true;
		if (threeLetterCode != other.threeLetterCode) return true;
		if (name != other.name) return true;
		return false;
	}
	bool operator==(const AminoAcidTypeRep & other) const {
		return !(*this != other);
	}

	char oneLetterCode;
	string threeLetterCode;
	string name;
	const AminoAcidType* myHandle;
};

AminoAcidType::AminoAcidType(char olc, const char* tlc, const char* name)
{
	rep = new AminoAcidTypeRep(olc, tlc, name, this);
}

AminoAcidType::~AminoAcidType()
{
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}


bool AminoAcidType::operator==(const AminoAcidType & src) const {return (*(this->rep)) == (*(src.rep));}
bool AminoAcidType::operator!=(const AminoAcidType & src) const {return (*(this->rep)) != (*(src.rep));}

char AminoAcidType::oneLetterCode() const {
	return rep->oneLetterCode;
}

const AminoAcidType * AminoAcidType::typeFor(char oneLetterCode) {
	const AminoAcidType * answer = AminoAcidTypeRep::typeFor(oneLetterCode);
	return answer;
}
