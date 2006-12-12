#ifndef _SIMTK_CHEM_AMINOACIDTYPE_H_
#define _SIMTK_CHEM_AMINOACIDTYPE_H_

class AminoAcidType {
public:
	
	// AminoAcidType();
	AminoAcidType(char oneLetterCode, const char* threeLetterCode, const char* name);
	~AminoAcidType();

	AminoAcidType(const AminoAcidType & src);
	AminoAcidType & operator=(const AminoAcidType & src);
	bool operator==(const AminoAcidType & src) const;
	bool operator!=(const AminoAcidType & src) const;

	char oneLetterCode() const;
	
	static const AminoAcidType * typeFor(char oneLetterCode);

    static const AminoAcidType Alanine;
    static const AminoAcidType Cysteine;
    static const AminoAcidType Aspartate;
    static const AminoAcidType Glutamate;
    static const AminoAcidType Phenylalanine;
    static const AminoAcidType Glycine;
    static const AminoAcidType Histidine;
    static const AminoAcidType Isoleucine;
    static const AminoAcidType Lysine;
    static const AminoAcidType Leucine;
    static const AminoAcidType Methionine;
    static const AminoAcidType Asparagine;
    static const AminoAcidType Proline;
    static const AminoAcidType Glutamine;
    static const AminoAcidType Arginine;
    static const AminoAcidType Serine;
    static const AminoAcidType Threonine;
    static const AminoAcidType Valine;
    static const AminoAcidType Tryptophan;
    static const AminoAcidType Tyrosine;
    
private:
	class AminoAcidTypeRep * rep;
};

#endif /*AMINOACIDTYPE_H_*/
