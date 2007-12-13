
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __SimTKOpenMMCommon_H__
#define __SimTKOpenMMCommon_H__

// include file containing entries commonly used

// STL includes

#include <vector>
#include <map>
#include <set>
#include <string>

// generic c includes

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

// ---------------------------------------------------------------------------------------

#include "SimTKOpenMMRealType.h"
#include "SimTKOpenMMWindowLinux.h"

// ---------------------------------------------------------------------------------------

typedef std::vector<int> IntVector;
typedef IntVector::iterator IntVectorI;
typedef IntVector::const_iterator IntVectorCI;
typedef IntVector::reverse_iterator IntVectorRI;
typedef IntVector::const_reverse_iterator IntVectorCRI;

typedef std::map<int, int> IntIntMap;
typedef IntIntMap::iterator IntIntMapI;
typedef IntIntMap::const_iterator IntIntMapCI;

typedef std::set<int> IntSet;
typedef IntSet::iterator IntSetI;
typedef IntSet::const_iterator IntSetCI;

typedef std::vector<IntSet> IntSetVector;
typedef IntSetVector::iterator IntSetVectorI;
typedef IntSetVector::const_iterator IntSetVectorCI;

typedef std::multiset<int> IntMultiSet;
typedef IntMultiSet::iterator IntMultiSetI;
typedef IntMultiSet::const_iterator IntMultiSetCI;

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;

typedef std::map<std::string, std::string> StringStringMap;
typedef StringStringMap::iterator StringStringMapI;
typedef StringStringMap::const_iterator StringStringMapCI;

typedef std::vector<RealOpenMM> RealOpenMMVector;
typedef RealOpenMMVector::iterator RealOpenMMVectorI;
typedef RealOpenMMVector::const_iterator RealOpenMMVectorCI;

typedef std::vector<RealOpenMM*> RealOpenMMPtrVector;
typedef RealOpenMMPtrVector::iterator RealOpenMMPtrVectorI;
typedef RealOpenMMPtrVector::const_iterator RealOpenMMPtrVectorCI;

typedef std::vector<RealOpenMM**> RealOpenMMPtrPtrVector;
typedef RealOpenMMPtrPtrVector::iterator RealOpenMMPtrPtrVectorI;
typedef RealOpenMMPtrPtrVector::const_iterator RealOpenMMPtrPtrVectorCI;

typedef std::map<std::string, RealOpenMM> StringRealOpenMMMap;
typedef StringRealOpenMMMap::iterator StringRealOpenMMMapI;
typedef StringRealOpenMMMap::const_iterator StringRealOpenMMMapCI;

typedef std::map<std::string, int> StringIntMap;
typedef StringIntMap::iterator StringIntMapI;
typedef StringIntMap::const_iterator StringIntMapCI;

// ---------------------------------------------------------------------------------------

class SimTKOpenMMCommon {

   public:

      static const std::string NotSet;
      static const RealOpenMM BigCutoffValue;
      static const std::string Comment;
      static const std::string Tab;
      static const std::string YesU;
      static const std::string YesUl;
      static const std::string YesL;

      // subroutine returns

      static const int DefaultReturn;
      static const int ErrorReturn;

      // units

      static const int MdUnits;
      static const int KcalAngUnits;

      // specify RealOpenMM number format
  
      static const int HighStringStreamNumberWidth;
      static const int HighStringStreamNumberPrecision;

      static const RealOpenMM DegreeToRadians;

      /**---------------------------------------------------------------------------------------
      
         Get string w/ minimum desired width
      
         @param inputString   input string
         @param desiredSize   minimum desired size of output string
         @param spacer        string to append to inputString to get desired width (" ")
      
         @return DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static std::string getSpacedString( const std::string&, int desiredSize = 25,
                                          const std::string& spacer = " " );
};

// ---------------------------------------------------------------------------------------

// StringComparisonForMap: used to compare strings in STL maps w/ strings as keys

#ifndef StringComparisonForMapPtrBlcok
#define StringComparisonForMapPtrBlcok

// above define should be removed once CommonTk is removed

class StringComparisonForMapPtr : public std::binary_function< std::string*, std::string*, bool > {

   public:

      bool operator()( const std::string* str1, const std::string* str2 ) const {
         return strcmp( str1->c_str(), str2->c_str() ) < 0;
      }
};

class StringComparisonForMap : public std::binary_function< const std::string&, const std::string&, bool > {

   public:

      bool operator()( const std::string& str1, const std::string& str2 ) const {
         return strcmp( str1.c_str(), str2.c_str() ) < 0;
      }
};
#endif

// ---------------------------------------------------------------------------------------

// Used to compare integers for sorting, ...

int numericCompareI( const void* point1, const void* point2 );

// ---------------------------------------------------------------------------------------

// string-string map and iterator definitions

typedef std::map<std::string, std::string, StringComparisonForMap> StringMap;
typedef StringMap::iterator StringMapI;
typedef StringMap::const_iterator StringMapCI;

// string set and iterator definitions

typedef std::set<std::string, StringComparisonForMap> StringSet;
typedef StringSet::iterator StringSetI;
typedef StringSet::const_iterator StringSetCI;

// ---------------------------------------------------------------------------------------

#endif // __SimTKOpenMMCommon_H__
