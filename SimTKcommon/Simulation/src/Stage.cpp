/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Stage.h"

using SimTK::Stage;

const Stage Stage::Empty;
const Stage Stage::Topology;
const Stage Stage::Model;
const Stage Stage::Instance;
const Stage Stage::Time;
const Stage Stage::Position;
const Stage Stage::Velocity;
const Stage Stage::Dynamics;
const Stage Stage::Acceleration;
const Stage Stage::Report;
const Stage Stage::Infinity;

const Stage Stage::LowestValid     = Empty;
const Stage Stage::HighestValid    = Infinity;
const Stage Stage::LowestRuntime   = Model;
const Stage Stage::HighestRuntime  = Report;

Stage::Stage() : Enumeration<Stage>() {
}

Stage::Stage(const Stage& thisElement, int index, const char* name) : Enumeration<Stage>(thisElement, index, name) {
}

void Stage::initValues() {
    new(&const_cast<Stage&>(Empty)) Stage(Empty, EmptyIndex, "Empty");
    new(&const_cast<Stage&>(Topology)) Stage(Topology, TopologyIndex, "Topology");
    new(&const_cast<Stage&>(Model)) Stage(Model, ModelIndex, "Model");
    new(&const_cast<Stage&>(Instance)) Stage(Instance, InstanceIndex, "Instance");
    new(&const_cast<Stage&>(Time)) Stage(Time, TimeIndex, "Time");
    new(&const_cast<Stage&>(Position)) Stage(Position, PositionIndex, "Position");
    new(&const_cast<Stage&>(Velocity)) Stage(Velocity, VelocityIndex, "Velocity");
    new(&const_cast<Stage&>(Dynamics)) Stage(Dynamics, DynamicsIndex, "Dynamics");
    new(&const_cast<Stage&>(Acceleration)) Stage(Acceleration, AccelerationIndex, "Acceleration");
    new(&const_cast<Stage&>(Report)) Stage(Report, ReportIndex, "Report");
    new(&const_cast<Stage&>(Infinity)) Stage(Infinity, InfinityIndex, "Infinity");
}
