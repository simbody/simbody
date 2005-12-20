#ifndef __cdsStringDefs_h__
#define __cdsStringDefs_h__

//
// issues: no nonconst operator[] 
//            (if there was, then would have to treat '\0' specially
//             with regards to length)
//

template<class CHAR> 
void
CDSString<CHAR>::resize(const int len)
  //
  // make rep large enough to hold string of length len
  //
{
 splitRep();

 if ( len+1>rep->asize ) {
   while ( len+1>rep->asize ) 
     //FIX: probably should be += asize
     rep->asize += (rep->blockSize?rep->blockSize:1);
   
   CHAR *nsl = new CHAR[rep->asize];
   strcpy(nsl,rep->sl);
   delete [] rep->sl;
   rep->sl = nsl;
 } 
 rep->len = -1;
 rep->sl[len] = '\0';

} /* resize */

 
template<class CHAR>
CDSString<CHAR>::CDSString(const CHAR* s,
                 int   slen,
               const int   bs) 
  //
  // constructor: copy from char*
  //
{
 if (slen>-1) {
   rep = new CDSStringRep<CHAR>(slen,bs);
   strncpy(rep->sl,s,slen);
   rep->sl[slen] = '\0';
   rep->len=slen;
 } else {
   if ( s ) {
     rep = new CDSStringRep<CHAR>(strlen(s),bs);
     strcpy(rep->sl,s);
   } else {
     rep = new CDSStringRep<CHAR>(0,bs);
     rep->sl[0] = '\0';
   }
 }
} /* copy constructor */

template<class CHAR>
CDSString<CHAR>::CDSString(const CDSString<CHAR>& cs) 
  //
  // copy constructor
  //
{
 rep = cs.rep;
 rep->count++;
// sl = new char[1+cs.length()];
// strcpy(sl,cs.sl);
// len=-1;
} /* copy constructor */

template<class CHAR>
CDSString<CHAR>&
CDSString<CHAR>::operator=(const CDSString<CHAR>& cs) 
  //
  // assignment operator
  //
{
  if (this == &cs) return *this; // self-assignment

  cs.rep->count++; 
  if (--rep->count <= 0) delete rep;
  
  rep = cs.rep;

//  delete sl;
//  len = cs.length();
//  sl = new char[1+len];
//  strcpy(sl,cs.sl);
  return *this;
} /* operator= */

template<class CHAR>
CDSString<CHAR>&
CDSString<CHAR>::operator=(const CHAR* s) 
  //
  // assignment operator- from CHAR*
  //
{
 int len = strlen(s);
 resize( len );
 strcpy(rep->sl,s);

 return *this;
} /* operator= (const char*) */

template<class CHAR>
CDSString<CHAR>&
CDSString<CHAR>::operator=(const CHAR s) 
  //
  // assignment operator
  //
{
 resize(1);
 rep->sl[0] = s; rep->sl[1] = '\0';

 return *this;
} /* operator= (const char*) */

template<class CHAR>
void
CDSString<CHAR>::calc_len() const
  // internal routine
{
 rep->len = strlen(rep->sl);
} /* calc_len */

template<class CHAR>
CDSString<CHAR>&
CDSString<CHAR>::operator+=(const CDSString<CHAR>& cs)
{
 int len = length() + cs.length();
 resize( len );
 strcat(rep->sl,cs.rep->sl);

 return *this;
} /* operator+= (CDSString) */
 
 
template<class CHAR>
CDSString<CHAR>&
CDSString<CHAR>::operator+=(const CHAR* s)
{
 int len = length()+strlen(s);
 resize( len );
 // cerr << "+=: " << s << endl;
 strcat(rep->sl,s);
 
 return *this;
} /* operator+= (const char*) */
 
 
template<class CHAR>
CDSString<CHAR>&
CDSString<CHAR>::operator+=(const CHAR c)
{
 int len = length()+1;
 resize( len );
 rep->sl[ len-1 ] = c;
 rep->sl[ len   ] = '\0';
 
 return *this;
} /* operator+= (const char c) */
 
template<class CHAR>
void 
CDSString<CHAR>::downcase()
{
 splitRep();

 for (int i=0 ; i<int(length()) ; i++) 
   rep->sl[i] = tolower(rep->sl[i]);
} /* downcase */

template<class CHAR>
void 
CDSString<CHAR>::upcase()
{
 splitRep();
 
 for (unsigned int i=0 ; i<length() ; i++) 
   rep->sl[i] = toupper(rep->sl[i]);
} /* upcase */

template<class CHAR>
bool
CDSString<CHAR>::contains(const CHAR* s) const
{
 int s_len = strlen(s);
 for (unsigned i=0 ; i<length() ; i++) {
   int j=0;
   if (rep->sl[i] == s[j]) {
     for (j=1 ; j<s_len ; j++)
       if (rep->sl[i+j] != s[j]) break;
     if (j==s_len)
       return 1;
   }
 }
 return 0;
} /* contains */

template<class CHAR>
int
CDSString<CHAR>::doGlob(const CHAR* text, 
            const CHAR* p   )
  //
  //  Do shell-style pattern matching for ?, \, [], and * characters.
  //  Might not be robust in face of malformed patterns; e.g., "foo[a-"
  //  could cause a segmentation violation.  It is 8bit clean.
  // 
  // from wildmat.c Revision: 1.1
  //
  //  Written by Rich $alz, mirror!rs, Wed Nov 26 19:03:17 EST 1986.
  //  Rich $alz is now <rsalz@osf.org>.

  //
  //  Match text and p, return TRUE, FALSE, or ABORT.
  //
{
 int matched;
 int reverse;

 for ( ; *p; text++, p++) {
   if (*text == '\0' && *p != '*')
     return Glob_ABORT;
   switch (*p) {
     case '\\':
       /* Literal match with following character. */
       p++;
       /* FALLTHROUGH */
     default:
       if (*text != *p)
     return Glob_FALSE;
       continue;
     case '?':
       /* Match anything. */
       continue;
     case '*':
       while (*++p == '*')
     /* Consecutive stars act just like one. */
     continue;
       if (*p == '\0')
     /* Trailing star matches everything. */
     return Glob_TRUE;
       while (*text)
     if ((matched = doGlob(text++, p)) != Glob_FALSE)
       return matched;
       return Glob_ABORT;
     case '[':
       reverse = p[1] == Glob_NEGATE ? Glob_TRUE : Glob_FALSE;
       if (reverse)
     /* Inverted character class. */
     p++;
       matched = Glob_FALSE;
       if (p[1] == ']' || p[1] == '-')
     if (*++p == *text)
       matched = Glob_TRUE;
       for (CHAR last = *p; *++p && *p != ']'; last = *p)
     /* This next line requires a good C compiler. */
     if (*p == '-' && p[1] != ']'
         ? *text <= *++p && *text >= last : *text == *p)
       matched = Glob_TRUE;
       if (matched == reverse)
     return Glob_FALSE;
       continue;
   }
 }

//    /* Do tar(1) matching rules, which ignore a trailing slash? */
//    if (*text == '/')
//      return Glob_TRUE;
 return *text == '\0';
} /* doGlob */


template<class CHAR>
bool
CDSString<CHAR>::glob(const CHAR* s         ,
              const bool  ignoreCase) const
{
// if (s == "*")   broken
//   return 1;

 CDSString<CHAR> tmpText = *this;
 CDSString<CHAR> tmpS    = s;
 if ( ignoreCase ) {
   tmpText.downcase();
   tmpS.downcase();
 }
 int ret = doGlob(tmpText, tmpS);
 return ret==Glob_TRUE;
} /* glob */


template<class CHAR>
int
CDSString<CHAR>::doGsub(      CDSString<CHAR>& s ,
            const CHAR*            s1,
            const CHAR*            s2)
{
 int cnt=0;
 int s1_len = strlen(s1);
 int s2_len = strlen(s2);
 if (s1_len<1) return cnt;

 for (unsigned int i=0 ; i<s.length() ; i++) {
   int j=0;
   for ( ; j<s1_len ; j++)
     if (s[i+j] != s1[j]) break;
   if (j==s1_len) {
     CDSString<CHAR> os = s;
     s.rep->sl[i] = '\0';
     s.calc_len();
     s += s2;
     s += ((const char*)os) + i + s1_len;
//     for (int m=0 ; m<os.length()-i-s1_len ; m++)
//    s += rep->sl[i+s1_len+m];
     i+=s2_len;
     cnt++;
   }
 }
 return cnt;
} /* doGsub */


template<class CHAR>
int
CDSString<CHAR>::gsub(const CHAR* s1,
              const CHAR* s2,
              const bool  recurse)
{
 splitRep();

 int cnt=0;
 cnt = doGsub(*this,s1,s2);
 while ( cnt && recurse ) {
   int contrib = doGsub(*this,s1,s2);
   if ( !contrib ) break;
   cnt += contrib;
 }

// for (int i=0 ; i<length() ; i++) {
//   int j=0;
//   for ( ; j<s1_len ; j++)
//     if (rep->sl[i+j] != s1[j]) break;
//   if (j==s1_len) {
//     int s2_len = strlen(s2);
//     rep->len = length() + s2_len - s1_len;
//     CHAR* new_sl = new CHAR[1+rep->len];
//     for (int k=0 ; k<i ; k++)
//     new_sl[k] = rep->sl[k];
//     for (int l=0 ; l<s2_len ; l++)
//     new_sl[l+i] = s2[l];
//     for (int m=0 ; m<length()-i-s1_len ; m++)
//     new_sl[i+s2_len+m] = rep->sl[i+s1_len+m];
//     new_sl[rep->len] = '\0';
//     delete [] rep->sl;
//     rep->sl = new_sl;
//     i+=s2_len;
//     cnt++;
//   }
// }
 
 return cnt;
} /* gsub */

template<class CHAR>
bool
CDSString<CHAR>::matches(const CHAR* str,
                   bool  ignoreCase) const
{
 if ( length() != strlen(str) )
   return 0;
 if ( ignoreCase ) {
   for (unsigned i=0 ; i<length() ; i++) 
     if ( tolower(rep->sl[i]) != tolower(str[i]) ) 
       return 0;
 } else {
   for (unsigned i=0 ; i<length() ; i++) 
     if ( rep->sl[i] != str[i] ) 
       return 0;
 }
 return 1;    
} /* matches */

template<class CHAR>
bool
operatorLess(const CHAR* s1,
         const CHAR* s2)
{
 return strcmp(s1,s2)<0;
} /* operatorLess */

template<class CHAR>
bool
operatorGreater(const CHAR* s1,
        const CHAR* s2)
{
 return strcmp(s1,s2)>0;
} /* operatorGreater */

template<class CHAR>
CDSList< CDSString<CHAR> >
CDSString<CHAR>::split(const CHAR* sep) const
{
 CDSString<CHAR> tmp = *this;
 CDSList< CDSString<CHAR> > ret;

 tmp.splitRep(); //strtok modifies the rep
 for ( CHAR* s = strtok(tmp.rep->sl,sep) ; s ; s = strtok(0,sep) )
   ret.append( s );

 return ret;
} /* split */

//bool
//operator==(const CDSString& s1,
//         const CDSString& s2)
//{
//  if ( s1.length() != s2.length() ) return 0;
//  for (int i=0 ; i<s1.length() ; i++) {
//    if ( s1[i] != s2[i] ) return 0;
//  }
//  return 1;    
//} /* operator== */
//
//bool
//operator==(const CDSString& s1,
//         const char*      str)
//{
//  if ( s1.length() != strlen(str) ) return 0;
//  for (int i=0 ; i<s1.length() ; i++) {
//    if ( s1[i] != str[i] ) return 0;
//  }
//  return 1;    
//} /* operator== */
//
//bool
//operator==(const char*      str,
//         const CDSString& s1)
//{
//  return operator==(s1,str);
//} /* operator== */

template<class CHAR>
CDSString<CHAR>
operator+(const CDSString<CHAR>& s1,
      const CDSString<CHAR>& s2)
{
 CDSString<CHAR> ret(s1);
 ret += s2;
 return ret;
} /* operator+ */

template<class CHAR>
CDSString<CHAR>
operator+(const CDSString<CHAR>& s1,
      const CHAR*            s2)
{
 CDSString<CHAR> ret(s1);
 ret += s2;
 return ret;
} /* operator+ */

template<class CHAR>
CDSString<CHAR>
operator+(const CDSString<CHAR>& s1,
        CHAR             c)
{
 CDSString<CHAR> ret(s1);
 ret += c;
 return ret;
} /* operator+ */

template<class CHAR>
ostream&
operator<<(ostream& s, const CDSString<CHAR>& x)
{
 if ( ! (s.flags()&ios::right) )
   for (int i=0 ; i<int(s.width()-x.length()) ; i++)
     s << s.fill();
 s << (const CHAR*)x;
// for (int i=0 ; i<x.length() ; i++)
//   s << x[i];
 if ( s.flags()&ios::right )
   for (int i=0 ; i<int(s.width()-x.length()) ; i++)
     s << s.fill();
 return s;
} /* operator<< */

// implementation note: the next two routines had been implemented using
// operator>> and putback, but Compaq C++ istrstream apparently doesn't make
// a copy of its string, so putback fails w/ segfault (because tries to
// write to readonly location). Hence the current implementation using
// peek and get

//
// read in up to (not including) white space
//
template<class CHAR>
istream& 
operator>>(istream&         s, 
           CDSString<CHAR>& x)
{
 FMTFLAGS_TYPE flags = s.flags();
 const bool skipWhite = (flags & ios::skipws) != 0;
 x = "";

 // if skipws not set, include all leading whitespace
 for (int c=s.peek(); c!=EOF && (c==' '||c=='\t'||c=='\n'); c=s.peek()) {
    CHAR ch = s.get();
    if (!skipWhite)
        x += ch;
 }

 s.unsetf(ios::skipws);
 for (int c=s.peek(); c!=EOF && c!=' '&&c!='\t'&&c!='\n'; c=s.peek()) 
     x += (CHAR)(s.get());
 s.setf(flags);

 return s;
}

template<class CHAR>
void 
readline(istream&         is,
         CDSString<CHAR>& x)
{
 // long flags = is.flags();
 // is.unsetf(ios::skipws);
 x = "";
 const int len=1024;
 char tmp[len];
 tmp[len-1] = '\0';
 do {
   tmp[0] = '\0';
   is.getline(tmp,len);
   if (tmp[len-1] == '\n')
     tmp[len-1] = '\0';
   x += tmp;
 } while ( strlen(tmp) == len-1 );

// for (int c = is.peek() ;
//    c!=EOF && c!='\n'  ;
//    c = is.peek()      )
//   x += (CHAR)(is.get());
 // is.setf(flags);
} /* readline */

template<class CHAR>
CDSString<CHAR> 
subString(const CDSString<CHAR>& s,
          const int              beg,
          int                    end)
{
 if ( beg>(int)s.length() || beg<0 )
   return "";
 
 if ( end<0 || end>(int)s.length() )
   end = s.length();

 CDSString<CHAR> ret(((const char*)s)+beg,end-beg);

 return ret;
} /* subString */

#endif /* __cdsStringDefs_h__ */
