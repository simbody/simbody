#ifndef SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_IMPL_H_
#define SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "SimTKcommon/internal/Parallel2DExecutor.h"
#include <utility>
#include <vector>

namespace SimTK {

class ParallelExecutor;

/**
 * This is the internal implementation class for Parallel2DExecutor.
 */

class Parallel2DExecutorImpl : public PIMPLImplementation<Parallel2DExecutor, Parallel2DExecutorImpl> {
public:
    class TriangleTask;
    class SquareTask;
    Parallel2DExecutorImpl(int gridSize, int numProcessors);
    ~Parallel2DExecutorImpl();
    void addSquare(int x, int y, int pass, int level);
    void addTriangle(int x, int y, int pass, int level);
    Parallel2DExecutorImpl* clone() const;
    void execute(Parallel2DExecutor::Task& task, Parallel2DExecutor::RangeType rangeType);
    int getBinStart(int bin) const {
        return binStart[bin];
    }
private:
    const int gridSize;
    std::vector<std::vector<std::pair<int,int> > > squares;
    std::vector<int> binStart;
    ParallelExecutor* executor;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_IMPL_H_
