/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Array.h"

#include "ArrayHelperImpl.h"

SimTKimpl::ArrayHelper::ArrayHelper(const SimTKimpl::TypeManipulatorT& tmt, int n)
	: impl(new ArrayHelperImpl(tmt,n))
{
}

SimTKimpl::ArrayHelper::ArrayHelper(const SimTKimpl::TypeManipulatorT& tmt, int n, const void* init, bool repeat)
	: impl(new ArrayHelperImpl(tmt, n, init, repeat))
{
}

SimTKimpl::ArrayHelper::~ArrayHelper()
{
	delete impl;
}

// Slicing, non-owner constructors
SimTKimpl::ArrayHelper::ArrayHelper(const ArrayHelper& ah, int offset, int length)
{
    impl = new ArrayHelperImpl(ah.getImpl(), offset, length);
}

SimTKimpl::ArrayHelper::ArrayHelper(ArrayHelper& ah, int offset, int length)
{
    assert(ah.impl);
    impl = new ArrayHelperImpl(ah.updImpl(), offset, length);
}

// Copy constructor
SimTKimpl::ArrayHelper::ArrayHelper(const ArrayHelper& ah)
	: impl(ah.impl ? new ArrayHelperImpl(ah.getImpl()) : 0)
{
}

SimTKimpl::ArrayHelper&
SimTKimpl::ArrayHelper::operator=(const ArrayHelper& ah) 
{
    assert(impl);
    (*impl) = ah.getImpl();
	return *this;
}

void
SimTKimpl::ArrayHelper::reverse() {
    return updImpl().reverse();
}

const void* 
SimTKimpl::ArrayHelper::operator[](int i) const
{
	return getImpl()[i];
}

void* 
SimTKimpl::ArrayHelper::operator[](int i)
{
	return updImpl()[i];
}

int
SimTKimpl::ArrayHelper::capacity() const
{
	return impl ? getImpl().capacity() : 0;
}

int
SimTKimpl::ArrayHelper::size() const
{
	return impl ? getImpl().size() : 0;
}

void
SimTKimpl::ArrayHelper::pop_back()
{
	updImpl().pop_back();
}

void
SimTKimpl::ArrayHelper::push_back(const void* x)
{
	updImpl().push_back(x);
}

void
SimTKimpl::ArrayHelper::resize(int n, const void* x)
{
	updImpl().resize(n,x);
}

void
SimTKimpl::ArrayHelper::reserve(int n)
{
    updImpl().reserve(n);
}

void
SimTKimpl::ArrayHelper::clear()
{
    updImpl().clear();
}


