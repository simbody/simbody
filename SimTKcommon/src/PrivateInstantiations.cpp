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

#define SimTK_SIMTKCOMMON_DEFINING_POLYGONALMESH
#define SimTK_SIMTKCOMMON_DEFINING_PARALLEL_EXECUTOR
#define SimTK_SIMTKCOMMON_DEFINING_PARALLEL_2D_EXECUTOR
#define SimTK_SIMTKCOMMON_DEFINING_PARALLEL_WORK_QUEUE
#include "../Geometry/src/PolygonalMeshImpl.h"
#include "ParallelExecutorImpl.h"
#include "Parallel2DExecutorImpl.h"
#include "ParallelWorkQueueImpl.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

namespace SimTK {

template class PIMPLHandle<PolygonalMesh, PolygonalMeshImpl, true>;
template class PIMPLImplementation<PolygonalMesh, PolygonalMeshImpl>;

template class PIMPLHandle<ParallelExecutor, ParallelExecutorImpl>;
template class PIMPLImplementation<ParallelExecutor, ParallelExecutorImpl>;

template class PIMPLHandle<Parallel2DExecutor, Parallel2DExecutorImpl>;
template class PIMPLImplementation<Parallel2DExecutor, Parallel2DExecutorImpl>;

template class PIMPLHandle<ParallelWorkQueue, ParallelWorkQueueImpl>;
template class PIMPLImplementation<ParallelWorkQueue, ParallelWorkQueueImpl>;

} // namespace SimTK
