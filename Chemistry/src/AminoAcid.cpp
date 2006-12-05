#include "AminoAcid.h"
#include <iostream> // NULL

class AminoAcidRep {
	friend class AminoAcid;
private:
	explicit AminoAcidRep(const AminoAcidType & aaType, const AminoAcid & handle) {
		initialize(&aaType);
		myHandle = &handle;
	}

	explicit AminoAcidRep(char oneLetterCode, const AminoAcid & handle) {
		const AminoAcidType * t = AminoAcidType::typeFor(oneLetterCode);
		initialize(t);
		myHandle = &handle;
	}

	void initialize(const AminoAcidType * aaType) {
		type = aaType;
	}

	char oneLetterCode() const {
		if (type == NULL) return 'X';
		return type->oneLetterCode();
	}
	
	AminoAcidRep* clone(AminoAcid & newHandle) const {
		AminoAcidRep* dup = new AminoAcidRep(*this);
		dup->myHandle = & newHandle;
		return dup;
	}

	bool operator!=(const AminoAcidRep & other) const {
		if (*type != *(other.type)) return true;
		return false;
	}
	bool operator==(const AminoAcidRep & other) const {
		return !(*this != other);
	}

	const AminoAcidType * type;
	const AminoAcid * myHandle;
};

AminoAcid::AminoAcid(const AminoAcidType & type)
{
	rep = new AminoAcidRep(type, *this);
	rep->myHandle = this;
}

AminoAcid::AminoAcid(char oneLetterCode)
{
	rep = new AminoAcidRep(oneLetterCode, *this);
	rep->myHandle = this;
}

AminoAcid::~AminoAcid()
{
	if (rep && (this == rep->myHandle)) delete rep;
	rep = NULL;
}

// Copy constructor
AminoAcid::AminoAcid(const AminoAcid & src) {
	if (this == &src) return;
	rep = NULL;
	*this = src;
}
AminoAcid & AminoAcid::operator=(const AminoAcid & src) {
	if (rep && (this == rep->myHandle)) delete rep;
	rep = src.rep->clone(*this);
	return *this;
}
bool AminoAcid::operator==(const AminoAcid & src) const {return (*(this->rep)) == (*(src.rep));}
bool AminoAcid::operator!=(const AminoAcid & src) const {return (*(this->rep)) != (*(src.rep));}

char AminoAcid::oneLetterCode() const {return rep->oneLetterCode();}

