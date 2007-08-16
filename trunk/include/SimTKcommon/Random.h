#ifndef SimTK_SimTKCOMMON_RANDOM_H_
#define SimTK_SimTKCOMMON_RANDOM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon/basics.h"

class SFMTData;

namespace SimTK {

class RandomImpl;

/**
 * This is a pseudo-random number generator.  It provides methods for generating random numbers that satisfy
 * uniform or Gaussian distributions, either one at a time or in bulk.
 * 
 * This class is implemented using the SIMD-oriented Fast Mersenne Twister (SFMT) library.  It provides
 * good performance, excellent statistical properties, and a very long period.
 * 
 * The methods of this class do not provide any synchronization or other mechanism to ensure thread safety.
 * It is therefore important that a single instance of this class not be accessed from multiple threads.
 */

class SimTK_SimTKCOMMON_EXPORT Random {
private:
	RandomImpl* impl;
public:
	Random();
	explicit Random(int seed);
	~Random();
	void setSeed(int seed);
	Real getReal();
	Real getGaussian();
	int getInt(int max);
	void fillArray(Real array[], int length);
	void fillArrayGaussian(Real array[], int length);
	void fillArray(int max, int array[], int length);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_RANDOM_H_
