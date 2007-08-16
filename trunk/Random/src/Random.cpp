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
#include "SimTKcommon/Random.h"
#include "SimTKcommon/internal/SFMT.h"

namespace SimTK {

/**
 * This is the private implementation class.
 */

class RandomImpl {
private:
	SFMTData* sfmt;
	static const int bufferSize = 1024;
	uint64_t buffer[bufferSize];
	int nextIndex;
	Real nextGaussian;
	bool nextGaussianIsValid;
public:
	RandomImpl() {
		sfmt = createSFMTData();
		setSeed(0);
	}

	RandomImpl(int seed) {
		sfmt = createSFMTData();
		setSeed(seed);
	}

	~RandomImpl() {
		deleteSFMTData(sfmt);
	}

	void setSeed(int seed) {
		nextIndex = bufferSize;
		nextGaussianIsValid = false;
		init_gen_rand(seed, *sfmt);
	}

	Real getReal() {
		if (nextIndex >= bufferSize) {
			// There are no remaining values in the buffer, so we need to refill it.
			
			fill_array64(buffer, bufferSize, *sfmt);
			nextIndex = 0;
		}
		return to_res53(buffer[nextIndex++]);
	}

	Real getGaussian() {
		if (nextGaussianIsValid) {
			nextGaussianIsValid = false;
			return nextGaussian;
		}
		
		// Use the polar form of the Box-Muller transformation to generate two Gaussian random numbers.
		
		Real x, y, r2;
		do {
	        x = 2.0*getReal()-1.0;
	        y = 2.0*getReal()-1.0;
	        r2 = x*x + y*y;
		} while (r2 >= 1.0 || r2 == 0.0);
		Real multiplier = sqrt((-2.0*log(r2))/r2);
		nextGaussian = y*multiplier;
		nextGaussianIsValid = true;
		return x*multiplier;
	}
	
	int getInt(int max) {
		return (int) floor(getReal()*max);
	}

	void fillArray(Real array[], int length) {
		for (int i = 0; i < length; ++i)
			array[i] = getReal();
	}

	void fillArrayGaussian(Real array[], int length) {
		for (int i = 0; i < length; ++i)
			array[i] = getGaussian();
	}

	void fillArray(int max, int array[], int length) {
		for (int i = 0; i < length; ++i)
			array[i] = getInt(max);
	}
};

/**
 * Create a new random number generator, and initialize it with a seed of 0.
 */

Random::Random() {
	impl = new RandomImpl();
}

/**
 * Create a new random number generator, and initialize it with a specified seed value.
 */

Random::Random(int seed) {
	impl = new RandomImpl(seed);
}

Random::~Random() {
	delete impl;
}

/**
 * Reinitialize this random number generator with a new seed value.
 */

void Random::setSeed(int seed) {
	impl->setSeed(seed);
}

/**
 * Get a random number, uniformly distributed in the range [0,1).
 */

Real Random::getReal() {
	return impl->getReal();
}

/**
 * Get a random number, chosen according to a Gaussian distribution with mean 0 and standard deviation 1.
 */

Real Random::getGaussian() {
	return impl->getGaussian();
}

/**
 * Get a random integer, uniformly distributed between 0 (inclusive) and max (exclusive).
 */

int Random::getInt(int max) {
	return impl->getInt(max);
}

/**
 * Fill an array with random numbers, uniformly distributed in the range [0,1).
 */

void Random::fillArray(Real array[], int length) {
	impl->fillArray(array, length);
}

/**
 * Fill an array with random numbers, chosen according to a Gaussian distribution with mean 0 and standard deviation 1.
 */

void Random::fillArrayGaussian(Real array[], int length) {
	impl->fillArrayGaussian(array, length);
}

/**
 * Fill an array with random integers, uniformly distributed between 0 (inclusive) and max (exclusive).
 */

void Random::fillArray(int max, int array[], int length) {
	impl->fillArray(max, array, length);
}

} // namespace SimTK
