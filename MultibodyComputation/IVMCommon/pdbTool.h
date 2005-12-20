#ifndef __pdbtool_hh__
#define __pdbtool_hh__

#include <atomSel.h>
#include <cdsList.h>
#include <cdsVector.h>
#include <cdsString.h>
#include <accessor.h>


//
// A class to read and write PDB files
//
// JJK 10/8/02
//

class PDBTool {

public:

  //
  // constructor/ (use default destructor)
  //

  PDBTool(const char*    filename=""       ,
	  const AtomSel& atomSel=AtomSel(""));

  //
  // accessors
  //

  ACCESSOR(selection,Selection,AtomSel);
  ACCESSOR(filename,Filename,String);
  ACCESSOR(makeBackup,MakeBackup,bool);

  void addRemark (const String &newVal);
  void clearRemarks ();

  CDSList< String > remarks() const;

  void useXplorNames();
  void useIupacNames();

  //
  // functionality
  //

  void read();
  void write();

  float_type aux1(const Atom&); // the occupancy field
  float_type aux2(const Atom&); // the b-factor field

  void setAux1(const Atom&,
	       const float_type&);
  void setAux2(const Atom&,
	       const float_type&);
  
  //synonyms
  float_type occupancy(const Atom& atom) { return aux1(atom); }
  float_type bfactor(const Atom& atom)   { return aux2(atom); }



protected:

  //
  // instance vbls
  //

  CDSList<String> myRemarksList_;
  int useXplorNames_;

  CDSVector< float_type > aux1_;
  CDSVector< float_type > aux2_;

};


#endif  /* __pdbtool_hh__ */
