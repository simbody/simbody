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

namespace SimTK {

class RandomImpl;

/**
 * This class defines the interface for pseudo-random number generators.  Subclasses generate numbers according to specific
 * distributions.  Currently, there are two such subclasses: Random::Uniform and Random::Gaussian.  For example, to generate
 * a series of pseudo-random numbers uniformly distributed between 0 and 100, you would call:
 * 
 * Random::Uniform random(0.0, 100.0);
 * Real nextValue = random.getValue(); // Each time you call this, it will return a different value.
 * 
 * Although the numbers are distributed in a seemingly random way, they are nonetheless deterministic, so if you create
 * several random number generators with the same parameters, each one will return exactly the same sequence of numbers.
 * If you want to get a different sequence of numbers, you can invoke setSeed(int seed) on a Random object.  Each seed
 * value corresponds to a different sequence of numbers that is uncorrelated with all others.  When a new Random object
 * is created, it is initialized with a seed value of 0.
 * 
 * This class is implemented using the SIMD-oriented Fast Mersenne Twister (SFMT) library.  It provides
 * good performance, excellent statistical properties, and a very long period.
 * 
 * The methods of this class do not provide any synchronization or other mechanism to ensure thread safety.
 * It is therefore important that a single Random object not be accessed from multiple threads.
 */

class SimTK_SimTKCOMMON_EXPORT Random {
public:
    class Uniform;
    class Gaussian;
    Random();
    ~Random();
    void setSeed(int seed);
    Real getValue();
    void fillArray(Real array[], int length);
protected:
    RandomImpl* impl;
    RandomImpl* getImpl() const;
};

/**
 * This is a subclass of Random that generates numbers uniformly distributed within a specified range.
 */

class SimTK_SimTKCOMMON_EXPORT Random::Uniform : public Random {
public:
    Uniform();
    Uniform(Real min, Real max);
    int getIntValue();
    Real getMin() const;
    void setMin(Real min);
    Real getMax() const;
    void setMax(Real max);
};

/**
 * This is a subclass of Random that generates numbers according to a Gaussian distribution with a
 * specified mean and standard deviation.
 */

class SimTK_SimTKCOMMON_EXPORT Random::Gaussian : public Random {
public:
    Gaussian();
    Gaussian(Real mean, Real stddev);
    Real getMean() const;
    void setMean(Real mean);
    Real getStdDev() const;
    void setStdDev(Real stddev);
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_RANDOM_H_
