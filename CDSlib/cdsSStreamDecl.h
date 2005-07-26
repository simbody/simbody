#ifndef __cdsSStreamDecl_h__
#define __cdsSStreamDecl_h__

#include <cdsIostream.h> 
#include <cdsString.h>

//  stream based on CDSString class
//
// FIX: currently only works with CHAR==char
//
//


template <class CHAR>
class CDSStringStreamBuf : public streambuf
{
public:
  CDSStringStreamBuf();
  CDSStringStreamBuf(const CDSString<CHAR>& s);
  virtual ~CDSStringStreamBuf();

  
  CDSString<CHAR> str() ;
  CDSString<CHAR> str_volatile() ;

private:
  CDSString<CHAR> str_;
  void setup(CDSString<CHAR> get,
	     CDSString<CHAR> put);

protected:                      // Overridden virtual member functions.
  virtual int underflow();
  //virtual int_type uflow();
  virtual int pbackfail(int __c = EOF);
  virtual int overflow(int __c = EOF);

  virtual streambuf* setbuf(char* __buf, int len);
//  virtual streampos seekoff(streamoff __off, ios::seek_dir __dir,
//			    int openmode= ios::in|ios::out);
//  virtual streampos seekpos(streampos __pos, int openmode=ios::in | ios::out);
};


////----------------------------------------------------------------------
//// CDSIStreamString: istream that manages a CDSStringStreamBuf.

template<class CHAR>
class CDSIStringStream : public istream
{
public:
  explicit CDSIStringStream( const CDSString<CHAR>& );
  explicit CDSIStringStream(const CHAR*);
  //CDSIStringStream(CHAR* , streamsize);
  //CDSIStringStream(const CHAR*, streamsize);
  virtual ~CDSIStringStream();
  
  //CDSStringStreamBuf<CHAR>* rdbuf() const;
  CDSString<CHAR> str() { return buf.str(); }

private:
  CDSStringStreamBuf<CHAR> buf;
};

//----------------------------------------------------------------------
// CDSOStreamString: Ostream that manages a CDSStringStreamBuf.

template<class CHAR>
class CDSOStringStream : public ostream
{
public:
  CDSOStringStream();
  //ostrstream(char*, int, ios_base::openmode = ios_base::out);
  virtual ~CDSOStringStream();

  // a copy of the strbuf CDSString is returned
  CDSString<CHAR> str() { return buf.str(); }
  //the string returned by str_volatile will change with future streambuf
  // operations - but no copy is made.
  CDSString<CHAR> str_volatile() { return buf.str_volatile(); }

private:
  CDSStringStreamBuf<CHAR> buf;
};


template<class CHAR>
class CDSStringStream : public iostream
{
public:
  CDSStringStream();
//  CDSStringStream(CDSString<CHAR>&, int, 
//		    std::ios_base::openmode = 
//		    std::ios_base::in | std::ios_base::out);
  virtual ~CDSStringStream();

  //CDSStringStreamBuf<CHAR>* rdbuf() const;
  //void freeze(bool = true);
  //int pcount() const;
  CDSString<CHAR> str() { return buf.str(); }
  CDSString<CHAR> str_volatile() { return buf.str_volatile(); }

private:
  CDSStringStreamBuf<CHAR> buf;
};

typedef CDSStringStream<char> StringStream;
typedef CDSOStringStream<char> OStringStream;
typedef CDSIStringStream<char> IStringStream;

int CDSStringStream_test();


#endif /* __cdsSStreamDecl_h__ */
