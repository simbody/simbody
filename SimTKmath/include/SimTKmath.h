#ifndef SimTK_SIMMATH_H_
#define SimTK_SIMMATH_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton, Michael Sherman, Peter Eastman                    *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/Geo_LineSeg.h"
#include "simmath/internal/Geo_Box.h"
#include "simmath/internal/Geo_Triangle.h"
#include "simmath/internal/Geo_CubicHermiteCurve.h"
#include "simmath/internal/Geo_BicubicHermitePatch.h"
#include "simmath/internal/Geo_CubicBezierCurve.h"
#include "simmath/internal/Geo_BicubicBezierPatch.h"
#include "simmath/internal/Spline.h"
#include "simmath/internal/SplineFitter.h"
#include "simmath/internal/BicubicSurface.h"
#include "simmath/internal/Geodesic.h"
#include "simmath/internal/GeodesicIntegrator.h"
#include "simmath/internal/ContactGeometry.h"
#include "simmath/internal/OrientedBoundingBox.h"
#include "simmath/internal/Contact.h"
#include "simmath/internal/ContactTracker.h"
#include "simmath/internal/CollisionDetectionAlgorithm.h"

#include "simmath/LinearAlgebra.h"
#include "simmath/Differentiator.h"
#include "simmath/Optimizer.h"
#include "simmath/MultibodyGraphMaker.h"
#include "simmath/Integrator.h"
#include "simmath/TimeStepper.h"
#include "simmath/CPodesIntegrator.h"
#include "simmath/RungeKuttaMersonIntegrator.h"
#include "simmath/RungeKuttaFeldbergIntegrator.h"
#include "simmath/RungeKutta3Integrator.h"
#include "simmath/RungeKutta2Integrator.h"
#include "simmath/ExplicitEulerIntegrator.h"
#include "simmath/VerletIntegrator.h"
#include "simmath/SemiExplicitEulerIntegrator.h"
#include "simmath/SemiExplicitEuler2Integrator.h"

#endif // SimTK_SIMMATH_H_
