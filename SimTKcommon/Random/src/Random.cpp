/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */


#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Random.h"
#include "SFMT.h"

#include <atomic>
#include <cassert>
#include <cmath>

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
    static std::atomic<int> nextSeed;
public:
    class UniformImpl;
    class GaussianImpl;
    
    RandomImpl() {
        sfmt = createSFMTData();
        nextIndex = bufferSize;
        init_gen_rand(++nextSeed, *sfmt);
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
        return Real(to_res53(buffer[nextIndex++]));
    }

    int getInt(int max) {
        return (int) floor(getValue()*max);
    }

    void fillArray(Real array[], int length) const {
        for (int i = 0; i < length; ++i)
            array[i] = getValue();
    }
};

std::atomic<int> Random::RandomImpl::nextSeed(0);

/**
 * This is the private implementation class for uniform random numbers.
 */

class Random::Uniform::UniformImpl : public Random::RandomImpl {
private:
    Real min, max, range;
public:
    UniformImpl(Real min, Real max) : min(min), max(max), range(max-min) {
    }
    
    Real getValue() const override {
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
    
    Real getValue() const override {
        if (nextGaussianIsValid) {
            nextGaussianIsValid = false;
            return mean+stddev*nextGaussian;
        }
        
        // Use the polar form of the Box-Muller transformation to generate two Gaussian random numbers.
        
        Real x, y, r2;
        do {
            x = 2*getNextRandom()-1;
            y = 2*getNextRandom()-1;
            r2 = x*x + y*y;
        } while (r2 >= 1.0 || r2 == 0.0);
        Real multiplier = std::sqrt((-2*std::log(r2))/r2);
        nextGaussian = y*multiplier;
        nextGaussianIsValid = true;
        return mean+stddev*x*multiplier;
    }
    
    void setSeed(int seed) override {
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

Random::Random() : impl(0) {
}

Random::~Random() {
    delete impl;
}

Random::RandomImpl& Random::getImpl() {
    assert(impl);
    return *impl;
}

const Random::RandomImpl& Random::getConstImpl() const {
    assert(impl);
    return *impl;
}

void Random::setSeed(int seed) {
    getImpl().setSeed(seed);
}

Real Random::getValue() const {
    return getConstImpl().getValue();
}


void Random::fillArray(Real array[], int length) const {
    getConstImpl().fillArray(array, length);
}

Random::Uniform::Uniform() {
    impl = new Random::Uniform::UniformImpl(0.0, 1.0);
}

Random::Uniform::Uniform(Real min, Real max) {
    impl = new Random::Uniform::UniformImpl(min, max);
}

Random::Uniform::UniformImpl& Random::Uniform::getImpl() {
    assert(impl);
    return SimTK_DYNAMIC_CAST_DEBUG<Random::Uniform::UniformImpl&>(*impl);
}

const Random::Uniform::UniformImpl& Random::Uniform::getConstImpl() const {
    assert(impl);
    return SimTK_DYNAMIC_CAST_DEBUG<Random::Uniform::UniformImpl&>(*impl);
}

int Random::Uniform::getIntValue() {
    return (int) std::floor(getImpl().getValue());
}

Real Random::Uniform::getMin() const {
    return getConstImpl().getMin();
}

void Random::Uniform::setMin(Real min) {
    getImpl().setMin(min);
}

Real Random::Uniform::getMax() const {
    return getConstImpl().getMax();
}

void Random::Uniform::setMax(Real max) {
    getImpl().setMax(max);
}

Random::Gaussian::Gaussian() {
    impl = new Random::Gaussian::GaussianImpl(0.0, 1.0);
}

Random::Gaussian::Gaussian(Real mean, Real stddev) {
    impl = new Random::Gaussian::GaussianImpl(mean, stddev);
}

Random::Gaussian::GaussianImpl& Random::Gaussian::getImpl() {
    assert(impl);
    return SimTK_DYNAMIC_CAST_DEBUG<Random::Gaussian::GaussianImpl&>(*impl);
}

const Random::Gaussian::GaussianImpl& Random::Gaussian::getConstImpl() const {
    assert(impl);
    return SimTK_DYNAMIC_CAST_DEBUG<Random::Gaussian::GaussianImpl&>(*impl);
}

Real Random::Gaussian::getMean() const {
    return getConstImpl().getMean();
}

void Random::Gaussian::setMean(Real mean) {
    getImpl().setMean(mean);
}

Real Random::Gaussian::getStdDev() const {
    return getConstImpl().getStdDev();
}

void Random::Gaussian::setStdDev(Real stddev) {
    getImpl().setStdDev(stddev);
}

} // namespace SimTK
