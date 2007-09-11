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

#include <cassert>
#include "SimTKcommon/basics.h"
#include "SimTKcommon/Random.h"
#include "SFMT.h"
using namespace SimTK_SFMT;

namespace SimTK {

/**
 * This is the private implementation class for Random.  It has a subclass corresponding to each subclass of Random.
 */

class Random::RandomImpl {
private:
    mutable SimTK_SFMT::SFMTData* sfmt;
    static const int bufferSize = 1024;
    mutable uint64_t buffer[bufferSize];
    mutable int nextIndex;
public:
    class UniformImpl;
    class GaussianImpl;
    
    RandomImpl() {
        sfmt = createSFMTData();
        nextIndex = bufferSize;
        init_gen_rand(0, *sfmt);
    }

    virtual ~RandomImpl() {
        deleteSFMTData(sfmt);
    }

    virtual void setSeed(int seed) {
        nextIndex = bufferSize;
        init_gen_rand(seed, *sfmt);
    }
    
    virtual Real getValue() const = 0;

    Real getNextRandom() const {
        if (nextIndex >= bufferSize) {
            // There are no remaining values in the buffer, so we need to refill it.
            
            fill_array64(buffer, bufferSize, *sfmt);
            nextIndex = 0;
        }
        return to_res53(buffer[nextIndex++]);
    }

    int getInt(int max) {
        return (int) floor(getValue()*max);
    }

    void fillArray(Real array[], int length) const {
        for (int i = 0; i < length; ++i)
            array[i] = getValue();
    }
};

/**
 * This is the private implementation class for uniform random numbers.
 */

class Random::Uniform::UniformImpl : public Random::RandomImpl {
private:
    Real min, max, range;
public:
    UniformImpl(Real min, Real max) : min(min), max(max), range(max-min) {
    }
    
    Real getValue() const {
        return min+getNextRandom()*range;
    }
    
    Real getMin() const {
        return min;
    }
    
    void setMin(Real value) {
        min = value;
        range = max-min;
    }
    
    Real getMax() const {
        return max;
    }
    
    void setMax(Real value) {
        max = value;
        range = max-min;
    }
};

/**
 * This is the private implementation class for Gaussian random numbers.
 */

class Random::Gaussian::GaussianImpl : public Random::RandomImpl {
private:
    Real mean, stddev;
    mutable Real nextGaussian;
    mutable bool nextGaussianIsValid;
public:
    GaussianImpl(Real mean, Real stddev) : mean(mean), stddev(stddev) {
        nextGaussianIsValid = false;
    }
    
    Real getValue() const {
        if (nextGaussianIsValid) {
            nextGaussianIsValid = false;
            return mean+stddev*nextGaussian;
        }
        
        // Use the polar form of the Box-Muller transformation to generate two Gaussian random numbers.
        
        Real x, y, r2;
        do {
            x = 2.0*getNextRandom()-1.0;
            y = 2.0*getNextRandom()-1.0;
            r2 = x*x + y*y;
        } while (r2 >= 1.0 || r2 == 0.0);
        Real multiplier = sqrt((-2.0*log(r2))/r2);
        nextGaussian = y*multiplier;
        nextGaussianIsValid = true;
        return mean+stddev*x*multiplier;
    }
    
    void setSeed(int seed) {
        RandomImpl::setSeed(seed);
        nextGaussianIsValid = false;
    }
    
    Real getMean() const {
        return mean;
    }
    
    void setMean(Real value) {
        mean = value;
    }
    
    Real getStdDev() const {
        return stddev;
    }
    
    void setStdDev(Real value) {
        stddev = value;
    }
};

/**
 * This constructor should never be invoked directly.  Instead, create an instance of one of the subclasses.
 */

Random::Random() : impl(0) {
}

Random::~Random() {
    delete impl;
}

/**
 * Get the internal object which implements the random number generator.
 */

Random::RandomImpl& Random::getImpl() {
    assert(impl);
    return *impl;
}

/**
 * Get a constant reference to the internal object which implements the random number generator.
 */

const Random::RandomImpl& Random::getConstImpl() const {
    assert(impl);
    return *impl;
}

/**
 * Reinitialize this random number generator with a new seed value.
 */

void Random::setSeed(int seed) {
    getImpl().setSeed(seed);
}

/**
 * Get the next value in the pseudo-random sequence.
 */

Real Random::getValue() const {
    return getConstImpl().getValue();
}

/**
 * Fill an array with values from the pseudo-random sequence.
 */

void Random::fillArray(Real array[], int length) const {
    getConstImpl().fillArray(array, length);
}

/**
 * Create a new random number generator that produces values uniformly distributed between 0 (inclusive) and 1 (exclusive).
 */

Random::Uniform::Uniform() {
    impl = new Random::Uniform::UniformImpl(0.0, 1.0);
}

/**
 * Create a new random number generator that produces values uniformly distributed between min (inclusive) and max (exclusive).
 */

Random::Uniform::Uniform(Real min, Real max) {
    impl = new Random::Uniform::UniformImpl(min, max);
}

Random::Uniform::UniformImpl& Random::Uniform::getImpl() {
    assert(impl);
    return dynamic_cast<Random::Uniform::UniformImpl&>(*impl);
}

const Random::Uniform::UniformImpl& Random::Uniform::getConstImpl() const {
    assert(impl);
    return dynamic_cast<Random::Uniform::UniformImpl&>(*impl);
}

/**
 * Get a random integer, uniformly distributed between 0 (inclusive) and max (exclusive).
 */

int Random::Uniform::getIntValue() {
    return (int) std::floor(getImpl().getValue());
}

/**
 * Get the lower end of the range in which values are uniformly distributed.
 */

Real Random::Uniform::getMin() const {
    return getConstImpl().getMin();
}

/**
 * Set the lower end of the range in which values are uniformly distributed.
 */

void Random::Uniform::setMin(Real min) {
    getImpl().setMin(min);
}

/**
 * Get the upper end of the range in which values are uniformly distributed.
 */

Real Random::Uniform::getMax() const {
    return getConstImpl().getMax();
}

/**
 * Set the upper end of the range in which values are uniformly distributed.
 */

void Random::Uniform::setMax(Real max) {
    getImpl().setMax(max);
}

/**
 * Create a new random number generator that produces values according to a Gaussian distribution with mean 0 and standard deviation 1.
 */

Random::Gaussian::Gaussian() {
    impl = new Random::Gaussian::GaussianImpl(0.0, 1.0);
}

/**
 * Create a new random number generator that produces values according to a Gaussian distribution with the specified mean and standard deviation.
 */

Random::Gaussian::Gaussian(Real mean, Real stddev) {
    impl = new Random::Gaussian::GaussianImpl(mean, stddev);
}

Random::Gaussian::GaussianImpl& Random::Gaussian::getImpl() {
    assert(impl);
    return dynamic_cast<Random::Gaussian::GaussianImpl&>(*impl);
}

const Random::Gaussian::GaussianImpl& Random::Gaussian::getConstImpl() const {
    assert(impl);
    return dynamic_cast<Random::Gaussian::GaussianImpl&>(*impl);
}

/**
 * Get the mean of the Gaussian distribution.
 */

Real Random::Gaussian::getMean() const {
    return getConstImpl().getMean();
}

/**
 * Set the mean of the Gaussian distribution.
 */

void Random::Gaussian::setMean(Real mean) {
    getImpl().setMean(mean);
}

/**
 * Get the standard deviation of the Gaussian distribution.
 */

Real Random::Gaussian::getStdDev() const {
    return getConstImpl().getStdDev();
}

/**
 * Set the standard deviation of the Gaussian distribution.
 */

void Random::Gaussian::setStdDev(Real stddev) {
    getImpl().setStdDev(stddev);
}

} // namespace SimTK
