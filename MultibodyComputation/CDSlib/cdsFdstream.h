
/* The following code declares classes to read from and write to
 * file descriptore or file handles.
 *
 * See
 *      http://www.josuttis.com/cppcode
 * for details and the latest version.
 *
 * - open:
 *      - integrating BUFSIZ on some systems?
 *      - optimized reading of multiple characters
 *      - stream for reading AND writing
 *      - i18n
 *
 * (C) Copyright Nicolai M. Josuttis 2001.
 * Permission to copy, use, modify, sell and distribute this software
 * is granted provided this copyright notice appears in all copies.
 * This software is provided "as is" without express or implied
 * warranty, and with no claim as to its suitability for any purpose.
 *
 * Version: Jul 28, 2002
 * History:
 *  Sep  6, 2002: fix so can use with old stream lib.
                  now, close fd in destructor
 *  Jul 28, 2002: bugfix memcpy() => memmove()
 *                fdinbuf::underflow(): cast for return statements
 *  Aug 05, 2001: first public version
 */
#ifndef BOOST_FDSTREAM_HPP
#define BOOST_FDSTREAM_HPP

#include <cdsIostream.h>
// for memmove():
#include <string.h>


// low-level read and write functions
#ifdef _MSC_VER
# include <io.h>
#else
# include <unistd.h>
//extern "C" {
//    int write (int fd, const char* buf, int num);
//    int read (int fd, char* buf, int num);
//}
#endif


namespace CDS {


/************************************************************
 * fdostream
 * - a stream that writes on a file descriptor
 ************************************************************/


  class fdoutbuf : public streambuf {
  protected:
    int fd;    // file descriptor
  public:
    // constructor
    fdoutbuf (int _fd) : fd( _fd ) {}
    virtual ~fdoutbuf() { close(fd); }
  protected:
    // write one character
    virtual int overflow (int c) {
        if (c != EOF) {
            char z = c;
            if (write (fd, &z, 1) != 1) {
                return EOF;
            }
        }
        return c;
    }
    // write multiple characters
    virtual
    int xsputn (const char* s,
		      int num) {
        return write(fd,s,num);
    }
};

  class fdostream : public ostream {
  protected:
    fdoutbuf buf;
  public:
    fdostream(int fd) : ostream(0), buf(fd) { ios::init(&buf); }
  };


/************************************************************
 * fdistream
 * - a stream that reads on a file descriptor
 ************************************************************/

class fdinbuf : public streambuf {
  protected:
    int fd;    // file descriptor
  protected:
    /* data buffer:
     * - at most, pbSize characters in putback area plus
     * - at most, bufSize characters in ordinary read buffer
     */
    static const int pbSize = 4;        // size of putback area
    static const int bufSize = 1024;    // size of the data buffer
    char buffer[bufSize+pbSize];        // data buffer

  public:
    /* constructor
     * - initialize file descriptor
     * - initialize empty data buffer
     * - no putback area
     * => force underflow()
     */
    fdinbuf (int _fd) : fd(_fd) {
        setg (buffer+pbSize,     // beginning of putback area
              buffer+pbSize,     // read position
              buffer+pbSize);    // end position
    }
    virtual ~fdinbuf() { close(fd); }

  protected:
    // insert new characters into the buffer
    virtual int underflow () {

        // is read position before end of buffer?
        if (gptr() < egptr()) {
	  return *gptr();
        }

        /* process size of putback area
         * - use number of characters read
         * - but at most size of putback area
         */
        int numPutback;
        numPutback = gptr() - eback();
        if (numPutback > pbSize) {
            numPutback = pbSize;
        }

        /* copy up to pbSize characters previously read into
         * the putback area
         */
        memmove (buffer+(pbSize-numPutback), gptr()-numPutback,
                numPutback);

        // read at most bufSize new characters
        int num;
        num = read (fd, buffer+pbSize, bufSize);
        if (num <= 0) {
            // ERROR or EOF
            return EOF;
        }

        // reset buffer pointers
        setg (buffer+(pbSize-numPutback),   // beginning of putback area
              buffer+pbSize,                // read position
              buffer+pbSize+num);           // end of buffer

        // return next character
        return *gptr();
    }
};

  class fdistream : public istream {
  protected:
    fdinbuf buf;
  public:
    fdistream (int fd) : istream(0), buf(fd) { ios::init(&buf); }
  };

#ifdef TESTING
  int test_fdstreams();
#endif /* TESTING */

} // END namespace CDS

#endif /*BOOST_FDSTREAM_HPP*/
