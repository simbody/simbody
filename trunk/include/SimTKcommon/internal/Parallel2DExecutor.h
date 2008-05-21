#ifndef SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_H_
#define SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_H_

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

#include "ParallelExecutor.h"
#include "PrivateImplementation.h"

namespace SimTK {

class Parallel2DExecutorImpl;

/**
 * This class is used for performing multithreaded computations over two dimensional ranges.  That is,
 * it performs some calculation once for each pair (i, j) where i and j vary over some range.  For
 * example, it is useful for calculating pairwise forces between a set of bodies.
 * 
 * To use it, define a subclass of Parallel2DExecutor::Task that performs a computation.  Then create a
 * Parallel2DExecutor object and ask it to execute the task:
 * 
 * <pre>
 * Parallel2DExecutor executor(gridSize);
 * executor.execute(myTask, Parallel2DExecutor::FullMatrix);
 * </pre>
 * 
 * The Task's execute() method will be called once with each pair (i, j) where i and j vary between 0 and
 * gridSize-1.  You also can restrict it to only pairs with i > j or i >= j.
 * 
 * The invocations are done in parallel on multiple threads, but they are divided up in a way that avoids
 * index conflicts between simultaneous calculations.  If the task is executed with indices (i1, j1)
 * on one thread, it is guaranteed that no other thread is simultaneously executing the task with
 * either the first or second index equal to either i1 or j1.  (More precisely, if either index of one
 * invocation is equal to either index of another invocation, the two invocations are guaranteed to be
 * separated by a happens-before edge.)  This allows the task to modify data that is indexed by i and j
 * without needing to worry about concurrent modifications.
 * 
 * The threads are created in the Parallel2DExecutor's constructor and remain active until it is deleted.
 * This means that creating a Parallel2DExecutor is a somewhat expensive operation, but it may then be
 * used repeatedly for executing various calculations.  By default, the number of threads is chosen
 * to be equal to the number of available processor cores.  You can optionally specify a different number
 * of threads to create.  For example, using more threads than processors can sometimes lead to better
 * processor utilitization.
 */

class SimTK_SimTKCOMMON_EXPORT Parallel2DExecutor : public PIMPLHandle<Parallel2DExecutor, Parallel2DExecutorImpl> {
public:
    class Task;
    enum RangeType {FullMatrix, HalfMatrix, HalfPlusDiagonal};
    /**
     * Construct a Parallel2DExecutor.
     * 
     * @param gridSize   the size of the range over which i and j should vary
     * @param numThreads the number of threads to create.  By default, this is set equal to the number
     * of processors.
     */
    Parallel2DExecutor(int gridSize, int numThreads = ParallelExecutor::getNumProcessors());
    /**
     * Execute a parallel task.
     * 
     * @param task      the Task to execute
     * @param rangeType specifies what part of the range i and j should vary over.  Specify FullyMatrix to
     *                  execute the task for all values of i and j between 0 and gridSize, HalfMatrix to
     *                  restrict it to i > j, and HalfPlusDiagonal to restrict it to i >= j. 
     */
    void execute(Task& task, RangeType rangeType);
};

/**
 * Concrete subclasses of this abstract class represent tasks that can be executed by a Parallel2DExecutor.
 */

class Parallel2DExecutor::Task {
public:
    virtual ~Task() {
    }
    /**
     * This method defines the task to be performed.  When the Task is passed to a Parallel2DExecutor's execute()
     * method, this method will be called in parallel for each allowed value of i and j.
     */
    virtual void execute(int i, int j) = 0;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PARALLEL_2D_EXECUTOR_H_
