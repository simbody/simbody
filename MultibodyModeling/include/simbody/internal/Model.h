#ifndef SIMTK_SIMBODY_MODEL_H_
#define SIMTK_SIMBODY_MODEL_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * User-visible, client side declaration of Model.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Subsystem.h"
#include "simbody/internal/State.h"

namespace simtk {

class Model {
public:
    Model() : rep(0) { }
    Model(const System&);
    ~Model();
    Model(const Model&);
    Model& operator=(const Model&);

    const System& getSystem() const;
    const State& getDefaultState() const;

private:
    class ModelRep* rep;
    friend class ModelRep;
};


} // namespace simtk

#endif // SIMTK_SIMBODY_MODEL_H_
