#ifdef NOTDEF

#include "cdsRegex.h"
#include "cdsString.h"
#include "cdsVector.h"

extern "C" {
#include <pcre.h>
}

static const int OVECTOR_LEN=30;

class CDSRegexLetter {
  pcre*  re;
  pcre_extra* extra;
  int    rc;  // result count
  CDSVector<int>  ovector;
  String expr;
public:
  CDSRegexLetter(const char* regex) {
   ovector.resize(OVECTOR_LEN);
   int erroffset;
   const char* error;
   re = pcre_compile(regex,            /* the pattern */
		     0,                /* default options */
		     &error,           /* for error message */
		     &erroffset,       /* for error offset */
		     0);            /* use default character tables */
   if ( !re ) {
     cerr << "CDSRegex: error in compilation: " << error << endl;
     throw CDSRegex::exception("failure in pcre_compile");
   }
   extra = pcre_study(re, 0, &error);
   if ( error ) {
     cerr << "CDSRegex: error in study: " << error << endl;
     throw CDSRegex::exception("failure in pcre_study");
   }
  }

  ~CDSRegexLetter() { pcre_free(re); pcre_free(extra); }
  friend class CDSRegex;
};

CDSRegex::CDSRegex(const char* regex)
{
 letter = new CDSRegexLetter(regex);
}

CDSRegex::~CDSRegex()
{
 delete letter;
}
 
int
CDSRegex::findIndex(const char* string,
		    const int   offset)
{
 letter->expr = string;
 letter->rc =
   pcre_exec(letter->re,     /* result of pcre_compile() */
	     letter->extra,  /* study info */
	     letter->expr,   /* the subject string */
	     letter->expr.length(),
	     offset,         /* start at offset 0 in the subject */
	     0,              /* default options */
	     letter->ovector.pointer(),/* vector for substring information */
	     OVECTOR_LEN   );/* number of elements in the vector */
 if ( letter->rc<1 )
   return -1;
 else
   return letter->ovector[0];
} /* findIndex */

String
CDSRegex::match() const
{
 if ( letter->rc >0 )
   return subString(letter->expr,letter->ovector[0],letter->ovector[1]);
 else
   return "";
} /* match */
 
int
CDSRegex::endIndex() const
{
 if ( letter->rc >0 )
   return letter->ovector[1];
 else
   return -1;
} /* match */
 


#ifdef TESTING
#include <cdsIomanip.h>
#include "cdsString.h"

int
CDSRegex::test()
{
 int exit=0;
 cout << "testing CDSRegex...";

 String s("a string of moderate length");
 CDSRegex re("\\bof\\b");

 if ( re.findIndex(s,0) != 9 ||
      String("of") != re.match() ||
      re.endIndex() != 11 ) {
   cerr << "erroneous output: " << re.findIndex(s,0) << ' '
	<< re.match() << ' ' << re.endIndex() << endl;
   exit=1;
 }


 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

#endif /* TESTING */

#endif
