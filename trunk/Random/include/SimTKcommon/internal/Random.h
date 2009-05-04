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

/**
 * This class defines the interface for pseudo-random number generators.  Subclasses generate numbers according to specific
 * distributions.  Currently, there are two such subclasses: Random::Uniform and Random::Gaussian.  For example, to generate
 * a series of pseudo-random numbers uniformly distributed between 0 and 100, you would call:
 * <pre>
 *   Random::Uniform random(0.0, 100.0);
 *   Real nextValue = random.getValue(); // Each time you call this, it will return a different value.
 * </pre>
 * Although the numbers are distributed in a seemingly random way, they are nonetheless deterministic, so you can create
 * several random number generators that each returns exactly the same sequence of numbers.  The sequence is determined
 * by the seed value with which the Random object is initialized.  By default, a different seed is used for every object.
 * You can invoke setSeed(int seed) on a Random object to explicitly specify the seed to use.  Each seed
 * value corresponds to a different sequence of numbers that is uncorrelated with all others.
 * 
 * This class is implemented using the SIMD-oriented Fast Mersenne Twister (SFMT) library.  It provides
 * good performance, excellent statistical properties, and a very long period.
 * 
 * The methods of this class do not provide any synchronization or other mechanism to ensure thread safety.
 * It is therefore important that a single Random object not be accessed from multiple threads. One minor
 * concession to threads: even if you don't set the seed explicitly, each thread's Random object will
 * use a different seed so you'll get a unique series of numbers in each thread.
 */

class SimTK_SimTKCOMMON_EXPORT Random {
public:
    class Uniform;
    class Gaussian;
    class RandomImpl;
    ~Random();
    /**
     * Reinitialize this random number generator with a new seed value.
     */
    void setSeed(int seed);
    /**
     * Get the next value in the pseudo-random sequence.
     */
    Real getValue() const;
    /**
     * Fill an array with values from the pseudo-random sequence.
     */
    void fillArray(Real array[], int length) const;
protected:
    RandomImpl* impl;
    /**
     * This constructor should never be invoked directly.  Instead, create an instance of one of the subclasses.
     */
    Random();
    /**
     * Get the internal object which implements the random number generator.
     */
    RandomImpl& getImpl();
    /**
     * Get a constant reference to the internal object which implements the random number generator.
     */
    const RandomImpl& getConstImpl() const;
private:
    Random(Random& r);
    Random operator=(Random& r);
};

/**
 * This is a subclass of Random that generates numbers uniformly distributed within a specified range.
 */

class SimTK_SimTKCOMMON_EXPORT Random::Uniform : public Random {
public:
    class UniformImpl;
    /**
     * Create a new random number generator that produces values uniformly distributed between 0 (inclusive) and 1 (exclusive).
     */
    Uniform();
    /**
     * Create a new random number generator that produces values uniformly distributed between min (inclusive) and max (exclusive).
     */
    Uniform(Real min, Real max);
    /**
     * Get a random integer, uniformly distributed between 0 (inclusive) and max (exclusive).
     */
    int getIntValue();
    /**
     * Get the lower end of the range in which values are uniformly distributed.
     */
    Real getMin() const;
    /**
     * Set the lower end of the range in which values are uniformly distributed.
     */
    void setMin(Real min);
    /**
     * Get the upper end of the range in which values are uniformly distributed.
     */
    Real getMax() const;
    /**
     * Set the upper end of the range in which values are uniformly distributed.
     */
    void setMax(Real max);
protected:
    UniformImpl& getImpl();
    const UniformImpl& getConstImpl() const;
};

/**
 * This is a subclass of Random that generates numbers according to a Gaussian distribution with a
 * specified mean and standard deviation.
 */

class SimTK_SimTKCOMMON_EXPORT Random::Gaussian : public Random {
public:
    class GaussianImpl;
    /**
     * Create a new random number generator that produces values according to a Gaussian distribution with mean 0 and standard deviation 1.
     */
    Gaussian();
    /**
     * Create a new random number generator that produces values according to a Gaussian distribution with the specified mean and standard deviation.
     */
    Gaussian(Real mean, Real stddev);
    /**
     * Get the mean of the Gaussian distribution.
     */
    Real getMean() const;
    /**
     * Set the mean of the Gaussian distribution.
     */
    void setMean(Real mean);
    /**
     * Get the standard deviation of the Gaussian distribution.
     */
    Real getStdDev() const;
    /**
     * Set the standard deviation of the Gaussian distribution.
     */
    void setStdDev(Real stddev);
protected:
    GaussianImpl& getImpl();
    const GaussianImpl& getConstImpl() const;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_RANDOM_H_
