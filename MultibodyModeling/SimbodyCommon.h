#ifndef SIMTK_SIMBODY_COMMON_H_
#define SIMTK_SIMBODY_COMMON_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * Common include file for all Simbody modules.
 */

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"

namespace simtk {

namespace Exception {

class FeatureNameNotFound : public Base {
public:
    FeatureNameNotFound(const char* fn, int ln, String method, String name) : Base(fn,ln)
    {
        setMessage(method + ": can't find any Feature named " + name);
    }
private:
};

}


} // namespace simtk

// rep helpers
// These allow a rep class to manipulate the rep pointer in its handle class.
// This avoids having to make the user-visible handle provide these routines
// or make the rep pointer public.
#define SIMTK_REP_HELPERS(HANDLECLASS,REPCLASS) \
static const REPCLASS* getRep(const HANDLECLASS& s)   {return s.rep;} \
static REPCLASS*       updRep(HANDLECLASS& s)         {return s.rep;} \
static void setRep(HANDLECLASS& s, REPCLASS* rep)     {assert(s.rep==0);s.rep=rep;} \
static void replaceRep(HANDLECLASS& s, REPCLASS* rep) {delete s.rep; s.rep=rep;}    \
static void clearRep(HANDLECLASS& s)                  {replaceRep(s,0);}

#endif // SIMTK_SIMBODY_COMMON_H_
