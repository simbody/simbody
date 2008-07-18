/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
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

#include "Parallel2DExecutorImpl.h"
#include "SimTKcommon/internal/ParallelExecutor.h"
#include <pthread.h>
#include <utility>
#include <vector>

using std::pair;
using std::vector;

namespace SimTK {

Parallel2DExecutorImpl::Parallel2DExecutorImpl(int gridSize, int numProcessors) : gridSize(gridSize), ownExecutor(true) {
    numProcessors = std::min(numProcessors, gridSize/2);
    if (numProcessors < 2)
        executor = 0;
    else
        executor = new ParallelExecutor(numProcessors);
    init(numProcessors);
}
Parallel2DExecutorImpl::Parallel2DExecutorImpl(int gridSize, ParallelExecutor& executor) : gridSize(gridSize), ownExecutor(false), executor(&executor) {
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
        
        squares.resize(bins);
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
        squares[pass].push_back(pair<int, int>(x, y));
    else {
        addSquare(2*x+0, 2*y+1, 2*pass+1, level-1);
        addSquare(2*x+1, 2*y+2, 2*pass+1, level-1);
        addSquare(2*x+0, 2*y+2, 2*pass+2, level-1);
        addSquare(2*x+1, 2*y+1, 2*pass+2, level-1);
    }
}
void Parallel2DExecutorImpl::addTriangle(int x, int y, int pass, int level) {
    if (level > 0) {
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
    TriangleTask(const Parallel2DExecutorImpl& executor, Parallel2DExecutor::Task& task, Parallel2DExecutor::RangeType rangeType, bool shouldInitialize, bool shouldFinish) :
        executor(executor), task(task), rangeType(rangeType), shouldInitialize(shouldInitialize), shouldFinish(shouldFinish) {
    }
    void execute(int index) {
        int start = executor.getBinStart(index);
        int end = executor.getBinStart(index+1);
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
    void initialize() {
        if (shouldInitialize)
            task.initialize();
    }
    void finish() {
        if (shouldFinish)
            task.finish();
    }
private:
    const Parallel2DExecutorImpl& executor;
    Parallel2DExecutor::Task& task;
    const Parallel2DExecutor::RangeType rangeType;
    bool shouldInitialize, shouldFinish;
};

class Parallel2DExecutorImpl::SquareTask : public ParallelExecutor::Task {
public:
    SquareTask(const Parallel2DExecutorImpl& executor, Parallel2DExecutor::Task& task, const vector<pair<int,int> >& squares, Parallel2DExecutor::RangeType rangeType, bool shouldInitialize, bool shouldFinish) :
        executor(executor), task(task), squares(squares), rangeType(rangeType), shouldInitialize(shouldInitialize), shouldFinish(shouldFinish) {
    }
    void execute(int index) {
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
    void initialize() {
        if (shouldInitialize)
            task.initialize();
    }
    void finish() {
        if (shouldFinish)
            task.finish();
    }
private:
    const Parallel2DExecutorImpl& executor;
    Parallel2DExecutor::Task& task;
    const vector<pair<int,int> >& squares;
    const Parallel2DExecutor::RangeType rangeType;
    bool shouldInitialize, shouldFinish;
};

void Parallel2DExecutorImpl::execute(Parallel2DExecutor::Task& task, Parallel2DExecutor::RangeType rangeType) {
    if (executor == 0) {
        task.initialize();
        TriangleTask(*this, task, rangeType, false, false).execute(0);
        task.finish();
        return;
    }
    int bins = binStart.size()-1;
    
    // Execute the blocks along the diagonal.
    
    TriangleTask triangle(*this, task, rangeType, true, false);
    executor->execute(triangle, bins);
    
    // Execute the square blocks in a series of passes.
    
    for (int i = 0; i < (int)squares.size(); ++i) {
        SquareTask square(*this, task, squares[i], rangeType, false, i == squares.size()-1);
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
