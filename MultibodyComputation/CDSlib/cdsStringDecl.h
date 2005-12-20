
//
// String template class
//

#ifndef __cdsStringDecl_h__
#define __cdsStringDecl_h__

#include "cdsList.h"

#include <cassert>

#include <cdsIostream.h>

template<class CHAR>
class CDSString_;

template<class CHAR>
class CDSStringRep {
  CHAR *sl;    
  mutable int len;
  int asize;        //number of CHAR allocated
  int blockSize;    //size for next allocation block
  int count;
  CDSStringRep() {}  //inaccessible
  CDSStringRep(const CDSStringRep&) {}  //inaccessible
  CDSStringRep<CHAR>& operator=(const CDSStringRep&) {}
  CDSStringRep(const CDSStringRep* r) : 
    len(r->len), asize(r->asize), blockSize(r->blockSize), count(1)
  { 
   assert(len>=0);
   sl = new CHAR[asize]; 
   for (int i=0 ; i<len+1 ; i++) sl[i] = r->sl[i];
  }
public:
  CDSStringRep(const int s,
	       const int bs) :
    len(s), asize(bs?((s+1)/bs+1)*bs:s+1), blockSize(bs), count(1)
  {
   sl = new CHAR[asize];
  } 
  //  CDSStringRep(const int s) : len(s) { sl = new CHAR[len+1]; count=1; }
  ~CDSStringRep() { delete [] sl; };
  friend class CDSString_<CHAR>;
};

template<class CHAR=char>
class CDSString_ {
  CDSStringRep<CHAR>* rep;
  void calc_len() const;

  inline void splitRep() 
    { if (rep->count>1) {
      rep->count--;
      calc_len();  //needed for next line
      rep= new CDSStringRep<CHAR>(rep);}
    }
public:
  CDSString_(const CHAR* s="",
		  int len=-1,
	    const int blockSize=8);
  CDSString_(const CDSString_&);
  ~CDSString_() { if (--rep->count <= 0) delete rep; }

  //  CDSString_ operator()(const char*);
  //  operator char*();
  operator const CHAR*() const {return rep->sl;}
  char* charPtr() { return rep->sl; }  // use this very carefully
  //a  operator const char*()       {return sl;}
  CHAR operator[](const unsigned int i) const {return rep->sl[i];} //check range?

  void resize(const int s);
  void blockSize(const int bs) { rep->blockSize = bs; }
  int  blockSize()             { return rep->blockSize; }

  CDSString_& operator=(const CHAR);
  CDSString_& operator=(const CDSString_&);
  CDSString_& operator=(const CHAR*);
  CDSString_& operator+=(const CDSString_&);
  CDSString_& operator+=(const CHAR*);
  CDSString_& operator+=(const CHAR);

  

  void downcase();
  void upcase();
  unsigned int length() const {if (rep->len<0) calc_len(); return rep->len;}

  CDSList<CDSString_> split(const CHAR* sep="\t ") const;

  bool contains(const CHAR*) const;
  bool matches(const CHAR* str,
		     bool  ignoreCase=0) const;

private:
  static int doGsub(      CDSString_<CHAR>& s,
		    const CHAR* s1,
		    const CHAR* s2);
public:
  int  gsub(const CHAR* s1,
	    const CHAR* s2,
	    const bool  recurse=0);

private:
  static int doGlob(const CHAR* text, 
		    const CHAR* p   );
  enum { Glob_TRUE=1 , Glob_FALSE=0 , Glob_ABORT=-1 ,
	 Glob_NEGATE='^' /* character marking inverted character class */ };
public:
  bool glob(const CHAR* str,
	    const bool  ignoreCase=0) const;
};

// comparison operators
// this could be considerably simplified if we redefine
// operatorX(const CHAR*, const CHAR*), but this seems quite unwise.
template<class CHAR>
inline bool operator==(const CDSString_<CHAR>& s1,const CDSString_<CHAR>& s2)
{return s1.matches(s2);}
template<class CHAR>
inline bool operator==(const CHAR* s1,const CDSString_<CHAR>& s2)
{return s2.matches(s1);}
template<class CHAR>
inline bool operator==(const CDSString_<CHAR>& s1,const CHAR* s2)
{return s1.matches(s2);}
template<class CHAR>
inline bool operator!=(const CDSString_<CHAR>& s1,const CDSString_<CHAR>& s2) 
{return !operator==(s1,s2);}
template<class CHAR>
inline bool operator!=(const CHAR* s1,const CDSString_<CHAR>& s2)
{return !operator==(s1,s2);}
template<class CHAR>
inline bool operator!=(const CDSString_<CHAR>& s1,const CHAR* s2)
{return !operator==(s1,s2);}

template<class CHAR>
bool operatorLess(const CHAR* s1,
		  const CHAR* s2);
template<class CHAR>
inline bool operator<(const CDSString_<CHAR>& s1,const CDSString_<CHAR>& s2)
{ return operatorLess((const CHAR*)s1,(const CHAR*)s2); }
template<class CHAR>
inline bool operator<(const CHAR* s1,const CDSString_<CHAR>& s2)
{ return operatorLess(s1,(const CHAR*)s2); }
template<class CHAR>
inline bool operator<(const CDSString_<CHAR>& s1,const CHAR* s2)
{return operatorLess((const CHAR*)s1,s2); }

template<class CHAR>
bool operatorGreater(const CHAR* s1,
		     const CHAR* s2);
template<class CHAR>
inline bool operator>(const CDSString_<CHAR>& s1,const CDSString_<CHAR>& s2)
{ return operatorGreater((const CHAR*)s1,(const CHAR*)s2); }
template<class CHAR>
inline bool operator>(const CHAR* s1,const CDSString_<CHAR>& s2)
{ return operatorGreater(s1,(const CHAR*)s2); }
template<class CHAR>
inline bool operator>(const CDSString_<CHAR>& s1,const CHAR* s2)
{return operatorGreater((const CHAR*)s1,s2); }



//other operators
template<class CHAR>
ostream& operator<<(ostream& s, const CDSString_<CHAR>& x);
template<class CHAR>
istream& operator>>(istream& s, CDSString_<CHAR>& x);
template<class CHAR>
void readline(istream&,CDSString_<CHAR>&);
template<class CHAR>
CDSString_<CHAR> operator+(const CDSString_<CHAR>&,const CDSString_<CHAR>&);
template<class CHAR>
CDSString_<CHAR> operator+(const CDSString_<CHAR>&,const CHAR*);
template<class CHAR>
CDSString_<CHAR> operator+(const CDSString_<CHAR>&,CHAR);
template<class CHAR>
CDSString_<CHAR> operator+(const CHAR* s1,const CDSString_<CHAR>& s2)
{ return CDSString_<CHAR>(s1) + s2; }
template<class CHAR>
CDSString_<CHAR> subString(const CDSString_<CHAR>&,
			  const int              beg,
				int              end=-1);

typedef CDSString_<char> CDSString;

int CDSString_test();

#endif /* __cdsStringDecl_h__ */
