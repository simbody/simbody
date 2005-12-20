
#ifndef __cdsRegex_hh__
#define __cdsRegex_hh__

//
// a frontend for pcre
//

#include "cdsExcept.h"
#include "cdsString.h"

class CDSRegexLetter;

class CDSRegex {
  CDSRegexLetter* letter;

public:
  CDSRegex(const char* regex);
  ~CDSRegex();
  
  int findIndex(const char* string,  //returns -1 on failure
		const int   offset);
  String match() const;              //returns most recent match
  int endIndex() const;              //index of 1st char after last match

  class exception: public CDS_NAMESPACE(exception) { 
  public:
    exception(const char* m) : CDS::exception(m) {}
  };


#ifdef TESTING
  static int test();
#endif  
};

#endif /* __cdsRegex_hh__ */

