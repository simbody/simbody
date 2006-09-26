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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Array.h"

#include "ArrayHelperImpl.h"

SimTKimpl::ArrayHelper::ArrayHelper(const SimTKimpl::TypeManipulatorT& tmt, ptrdiff_t n)
	: impl(new ArrayHelperImpl(tmt,n))
{
}

SimTKimpl::ArrayHelper::ArrayHelper(const SimTKimpl::TypeManipulatorT& tmt, ptrdiff_t n, const void* init, bool repeat)
	: impl(new ArrayHelperImpl(tmt, n, init, repeat))
{
}

SimTKimpl::ArrayHelper::~ArrayHelper()
{
	delete impl;
}

// Slicing, non-owner constructors
SimTKimpl::ArrayHelper::ArrayHelper(const ArrayHelper& ah, ptrdiff_t offset, ptrdiff_t length)
{
    impl = new ArrayHelperImpl(ah.getImpl(), offset, length);
}

SimTKimpl::ArrayHelper::ArrayHelper(ArrayHelper& ah, ptrdiff_t offset, ptrdiff_t length)
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
SimTKimpl::ArrayHelper::operator[](ptrdiff_t i) const
{
	return getImpl()[i];
}

void* 
SimTKimpl::ArrayHelper::operator[](ptrdiff_t i)
{
	return updImpl()[i];
}

ptrdiff_t
SimTKimpl::ArrayHelper::capacity() const
{
	return impl ? getImpl().capacity() : 0;
}

ptrdiff_t
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
SimTKimpl::ArrayHelper::resize(ptrdiff_t n, const void* x)
{
	updImpl().resize(n,x);
}

void
SimTKimpl::ArrayHelper::reserve(ptrdiff_t n)
{
    updImpl().reserve(n);
}

void
SimTKimpl::ArrayHelper::clear()
{
    updImpl().clear();
}


