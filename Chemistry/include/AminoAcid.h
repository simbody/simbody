#ifndef _SIMTK_CHEM_AMINOACID_H_
#define _SIMTK_CHEM_AMINOACID_H_

#include "AminoAcidType.h"

class AminoAcid
{
public:
	explicit AminoAcid(char oneLetterCode);
	explicit AminoAcid(const AminoAcidType & type);
	~AminoAcid();

	AminoAcid(const AminoAcid & src);
	AminoAcid & operator=(const AminoAcid & src);
	bool operator==(const AminoAcid & src) const;
	bool operator!=(const AminoAcid & src) const;
	
	char oneLetterCode() const;
private:
	class AminoAcidRep * rep;
};

#endif /*AMINOACID_H_*/
