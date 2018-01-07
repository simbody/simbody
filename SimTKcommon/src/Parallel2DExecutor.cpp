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

#include "Parallel2DExecutorImpl.h"
#include "SimTKcommon/internal/ParallelExecutor.h"
#include <utility>

using std::pair;

namespace SimTK {

Parallel2DExecutorImpl::Parallel2DExecutorImpl(int gridSize, int numProcessors) : gridSize(gridSize), ownExecutor(true) {
    numProcessors = std::min(numProcessors, gridSize/2);
    if (numProcessors < 2)
        executor = 0;
    else
        executor = new ParallelExecutor(numProcessors);
    init(numProcessors);
}
Parallel2DExecutorImpl::Parallel2DExecutorImpl(int gridSize, ParallelExecutor& executor) : gridSize(gridSize), executor(&executor), ownExecutor(false) {
    init(executor.getNumProcessors());
}
void Parallel2DExecutorImpl::init(int numProcessors) {
    int bins;
    if (numProcessors < 2) {
        bins = 1;
    }
    else {
        
        // Determine how many levels of subdivision to use.
        
        int levels = 1;
        while (1<<levels < numProcessors)
            levels++;
        levels++;
        bins = 1<<levels;
        
        // Build the set of squares.
        
        squares.resize(bins-1);
        addTriangle(0, 0, 0, levels);
    }
    
    // Find the range of indices in each bin.
    
    binStart.resize(bins+1);
    for (int i = 0; i < bins; ++i)
        binStart[i] = (int) std::floor(0.5+i*gridSize/(double) bins);
    binStart[bins] = gridSize;
}
Parallel2DExecutorImpl::~Parallel2DExecutorImpl() {
    if (executor != 0 && ownExecutor)
        delete executor;
}
void Parallel2DExecutorImpl::addSquare(int x, int y, int pass, int level) {
    if (level == 0)
        squares[pass-1].push_back(pair<int, int>(x, y));
    else {
        addSquare(2*x+0, 2*y+1, 2*pass+1, level-1);
        addSquare(2*x+1, 2*y+2, 2*pass+1, level-1);
        addSquare(2*x+0, 2*y+2, 2*pass+2, level-1);
        addSquare(2*x+1, 2*y+1, 2*pass+2, level-1);
    }
}
void Parallel2DExecutorImpl::addTriangle(int x, int y, int pass, int level) {
    if (level > 1) {
        addSquare(2*x, 2*y, 2*pass, level-1);
        addTriangle(2*x, 2*y, 2*pass, level-1);
        addTriangle(2*x+1, 2*y+1, 2*pass, level-1);
    }
}
Parallel2DExecutorImpl* Parallel2DExecutorImpl::clone() const {
    if (ownExecutor)
        return new Parallel2DExecutorImpl(gridSize, executor->getNumProcessors());
    return new Parallel2DExecutorImpl(gridSize, *executor);
}
ParallelExecutor& Parallel2DExecutorImpl::getExecutor() {
    return *executor;
}

class Parallel2DExecutorImpl::TriangleTask : public ParallelExecutor::Task {
public:
    TriangleTask(const Parallel2DExecutorImpl& executor, Parallel2DExecutor::Task& task, Parallel2DExecutor::RangeType rangeType, int width, bool shouldInitialize, bool shouldFinish) :
        executor(executor), task(task), rangeType(rangeType), width(width), shouldInitialize(shouldInitialize), shouldFinish(shouldFinish) {
    }
    void execute(int index) override {
        int start = executor.getBinStart(width*index);
        int end = executor.getBinStart(width*(index+1));
        switch (rangeType) {
        case Parallel2DExecutor::FullMatrix:
            for (int i = start; i < end; ++i)
                for (int j = start; j < end; ++j)
                    task.execute(i, j);
            return;
        case Parallel2DExecutor::HalfMatrix:
            for (int i = start; i < end; ++i)
                for (int j = start; j < i; ++j)
                    task.execute(i, j);
            return;
        case Parallel2DExecutor::HalfPlusDiagonal:
            for (int i = start; i < end; ++i)
                for (int j = start; j <= i; ++j)
                    task.execute(i, j);
            return;
        }
    }
    void initialize() override {
        if (shouldInitialize)
            task.initialize();
    }
    void finish() override {
        if (shouldFinish)
            task.finish();
    }
private:
    const Parallel2DExecutorImpl& executor;
    Parallel2DExecutor::Task& task;
    const Parallel2DExecutor::RangeType rangeType;
    const int width;
    const bool shouldInitialize, shouldFinish;
};

class Parallel2DExecutorImpl::SquareTask : public ParallelExecutor::Task {
public:
    SquareTask(const Parallel2DExecutorImpl& executor, Parallel2DExecutor::Task& task, const Array_<pair<int,int> >& squares, Parallel2DExecutor::RangeType rangeType, bool shouldInitialize, bool shouldFinish) :
        executor(executor), task(task), squares(squares), rangeType(rangeType), shouldInitialize(shouldInitialize), shouldFinish(shouldFinish) {
    }
    void execute(int index) override {
        const pair<int,int>& square = squares[index];
        int istart = executor.getBinStart(square.second+1);
        int iend = executor.getBinStart(square.second+2);
        int jstart = executor.getBinStart(square.first);
        int jend = executor.getBinStart(square.first+1);
        switch (rangeType) {
        case Parallel2DExecutor::FullMatrix:
            for (int i = istart; i < iend; ++i)
                for (int j = jstart; j < jend; ++j) {
                    task.execute(i, j);
                    task.execute(j, i);
                }
            return;
        case Parallel2DExecutor::HalfMatrix:
        case Parallel2DExecutor::HalfPlusDiagonal:
            for (int i = istart; i < iend; ++i)
                for (int j = jstart; j < jend; ++j)
                    task.execute(i, j);
            return;
        }
    }
    void initialize() override {
        if (shouldInitialize)
            task.initialize();
    }
    void finish() override {
        if (shouldFinish)
            task.finish();
    }
private:
    const Parallel2DExecutorImpl& executor;
    Parallel2DExecutor::Task& task;
    const Array_<pair<int,int> >& squares;
    const Parallel2DExecutor::RangeType rangeType;
    bool shouldInitialize, shouldFinish;
};

void Parallel2DExecutorImpl::execute(Parallel2DExecutor::Task& task, Parallel2DExecutor::RangeType rangeType) {
    if (executor == 0) {
        task.initialize();
        TriangleTask(*this, task, rangeType, 1, false, false).execute(0);
        task.finish();
        return;
    }
    int bins = binStart.size()-1;
    
    // Execute the blocks along the diagonal.
    
    TriangleTask triangle(*this, task, rangeType, 2, true, false);
    executor->execute(triangle, bins/2);
    
    // Execute the square blocks in a series of passes.
    
    for (int i = 0; i < (int)squares.size(); ++i) {
        SquareTask square(*this, task, squares[i], rangeType, false, 
                          i == (int)squares.size()-1);
        executor->execute(square, squares[i].size());
    }
}

Parallel2DExecutor::Parallel2DExecutor(int gridSize, int numProcessors) : HandleBase(new Parallel2DExecutorImpl(gridSize, numProcessors)) {
}

Parallel2DExecutor::Parallel2DExecutor(int gridSize, ParallelExecutor& executor) : HandleBase(new Parallel2DExecutorImpl(gridSize, executor)) {
}

void Parallel2DExecutor::execute(Task& task, RangeType rangeType) {
    updImpl().execute(task, rangeType);
}

ParallelExecutor& Parallel2DExecutor::getExecutor() {
    return updImpl().getExecutor();
}

} // namespace SimTK
