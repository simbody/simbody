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

#include "SimTKcommon.h"
#include "SimTKcommon/Random.h"

#include <assert.h>
#include <iostream>
using std::cout;
using std::endl;
using std::sqrt;
using namespace SimTK;

/**
 * Given the number of values expected and found in a set of bins, verify that the distribution is correct.
 */

void verifyDistribution(int expected[], int found[], int bins) {
	for (int i = 0; i < bins; ++i) {
		Real dev = sqrt((Real) expected[i]);
		assert(found[i] >= expected[i]-4*dev && found[i] <= expected[i]+4*dev);
	}
}

/**
 * Given a set of Reals, verify that they satisfy a uniform distribution between 0 and 1.
 */

void verifyUniformDistribution(Real value[], int length) {
	int expected[10], found[10];
	for (int i = 0; i < 10; ++i) {
		expected[i] = length/10;
		found[i] = 0;
	}
	for (int i = 0; i < length; ++i) {
		assert(value[i] >= 0.0);
		assert(value[i] < 1.0);
		found[(int) (value[i]*10)]++;
	}
	verifyDistribution(expected, found, 10);
}

/**
 * Given a set of ints, verify that they satisfy a uniform distribution between 0 and max.
 */

void verifyUniformDistribution(int max, int value[], int length) {
	int* expected = new int[max];
	int* found = new int[max];
	for (int i = 0; i < max; ++i) {
		expected[i] = length/max;
		found[i] = 0;
	}
	for (int i = 0; i < length; ++i) {
		assert(value[i] >= 0);
		assert(value[i] < max);
		found[value[i]]++;
	}
	verifyDistribution(expected, found, max);
	delete[] expected;
	delete[] found;
}

/**
 * Given a set of values, verify that they satisfy a Gaussian distribution.
 */

void verifyGaussianDistribution(Real value[], int length) {
	int expected[6], found[6];
	expected[0] = expected[5] = 0.0228*length;
	expected[1] = expected[4] = 0.1587*length-expected[0];
	expected[2] = expected[3] = 0.5*length-expected[1];
	for (int i = 0; i < 6; ++i)
		found[i] = 0;
	for (int i = 0; i < length; ++i) {
		if (value[i] < -2)
			found[0]++;
		else if (value[i] < -1)
			found[1]++;
		else if (value[i] < 0)
			found[2]++;
		else if (value[i] < 1)
			found[3]++;
		else if (value[i] < 2)
			found[4]++;
		else
			found[5]++;
	}
	verifyDistribution(expected, found, 6);
}

int main() {
    try {
    	Random rand(1);

    	// Try generating a bunch of random numbers, and make sure they are distributed uniformly between 0 and 1.
    	
    	Real value[2001];
    	value[2000] = 123.4;
    	for (int i = 0; i < 2000; ++i)
    		value[i] = rand.getReal();
    	verifyUniformDistribution(value, 2000);
    	
    	// Reset the random number generator, and make sure it produces the same values again.
    	
    	rand.setSeed(1);
    	for (int i = 0; i < 2000; ++i)
    		assert(value[i] == rand.getReal());
    	
    	// Now try asking for a whole array at a time, and verify that it still gives the same results.
    	
    	Real value2[2001];
    	value2[2000] = 567.8;
    	rand.setSeed(1);
    	rand.fillArray(value2, 2000);
    	for (int i = 0; i < 2000; ++i)
    		assert(value[i] == value2[i]);
    	
    	// Set the seed to a different value, and verify that the results are different.
    	
    	rand.setSeed(2);
    	rand.fillArray(value2, 2000);
    	for (int i = 0; i < 2000; ++i)
    		assert(value[i] != value2[i]);
    	
    	// Now try generating a bunch of Gaussian random numbers, and check the distribution.
    	
    	rand.setSeed(1);
    	for (int i = 0; i < 2000; ++i)
    		value[i] = rand.getGaussian();
    	verifyGaussianDistribution(value, 2000);
    	
    	// Try getting a whole array at a time.
    	
    	rand.setSeed(1);
    	rand.fillArrayGaussian(value2, 2000);
    	for (int i = 0; i < 2000; ++i)
    		assert(value[i] == value2[i]);
    	
    	// Try generating uniform integers.
    	
    	int value3[2001];
    	value3[2000] = -99;
    	rand.setSeed(3);
    	for (int i = 0; i < 2000; ++i)
    		value3[i] = rand.getInt(20);
    	verifyUniformDistribution(20, value3, 2000);
    	
    	// Try getting a whole array at a time.
    	
    	int value4[2001];
    	value4[2000] = -111;
    	rand.setSeed(3);
    	rand.fillArray(20, value4, 2000);
    	for (int i = 0; i < 2000; ++i)
    		assert(value3[i] == value4[i]);
    	
    	// Make sure none of the above operations has overwritten the final array element.
    	
    	assert(value[2000] == 123.4);
    	assert(value2[2000] = 567.8);
    	assert(value3[2000] = -99);
    	assert(value4[2000] = -111);
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    return 0;
}


