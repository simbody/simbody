#ifndef SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_IMPL_H_
#define SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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

#include "SimTKcommon/internal/Parallel2DExecutor.h"
#include "SimTKcommon/internal/Array.h"
#include <utility>

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
    Parallel2DExecutorImpl(int gridSize, ParallelExecutor& executor);
    ~Parallel2DExecutorImpl();
    void init(int numProcessors);
    void addSquare(int x, int y, int pass, int level);
    void addTriangle(int x, int y, int pass, int level);
    Parallel2DExecutorImpl* clone() const;
    void execute(Parallel2DExecutor::Task& task, Parallel2DExecutor::RangeType rangeType);
    ParallelExecutor& getExecutor();
    int getBinStart(int bin) const {
        return binStart[bin];
    }
private:
    const int gridSize;
    Array_<Array_<std::pair<int,int> > > squares;
    Array_<int> binStart;
    ParallelExecutor* executor;
    bool ownExecutor;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_IMPL_H_
