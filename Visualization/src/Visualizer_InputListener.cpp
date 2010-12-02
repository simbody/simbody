/* -------------------------------------------------------------------------- *
 *                                Simbody(tm)                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

/** @file
This file contains the implementations of Visualizer::InputListener and any
built-in concrete InputListeners, currently just Visualizer::InputSilo. **/

#include "simbody/internal/common.h"
#include "simbody/internal/Visualizer_InputListener.h"

#include <deque>
#include <utility>
#include <pthread.h>

using namespace SimTK;

//==============================================================================
//                              INPUT SILO IMPL
//==============================================================================
/* This is the private implementation object contained in a 
Visualizer::InputSilo handle. */

class Visualizer::InputSilo::Impl {
public:
    Impl() : m_inputCount(0) {pthread_mutex_init(&m_siloLock,0);}
    ~Impl() {pthread_mutex_destroy(&m_siloLock);}

    void LOCK_silo()   {pthread_mutex_lock(&m_siloLock);}
    void UNLOCK_silo() {pthread_mutex_unlock(&m_siloLock);}

    pthread_mutex_t                             m_siloLock;
    std::deque<std::pair<unsigned,unsigned> >   m_keyHitSilo;
    std::deque<int>                             m_menuPickSilo;
    std::deque<std::pair<unsigned,Real> >       m_sliderMoveSilo;
    unsigned                                    m_inputCount;
};

//==============================================================================
//                                INPUT SILO
//==============================================================================
Visualizer::InputSilo::InputSilo() {m_impl = new InputSilo::Impl();}
Visualizer::InputSilo::~InputSilo() {delete m_impl;}

bool Visualizer::InputSilo::isAnyUserInput() const 
{   return getImpl().m_inputCount != 0; }

bool Visualizer::InputSilo::takeKeyHit(unsigned& key, unsigned& modifiers) {
    Impl& impl = updImpl(); bool gotOne;
    impl.LOCK_silo();
    if (impl.m_keyHitSilo.empty()) key=0, modifiers=0, gotOne=false;
    else {
        key       = impl.m_keyHitSilo.front().first; 
        modifiers = impl.m_keyHitSilo.front().second;
        impl.m_keyHitSilo.pop_front();
        --impl.m_inputCount;
        gotOne = true;
    }
    impl.UNLOCK_silo();
    return gotOne;
}

bool Visualizer::InputSilo::takeMenuPick(int& item) {
    Impl& impl = updImpl(); bool gotOne;
    impl.LOCK_silo();
    if (impl.m_menuPickSilo.empty()) item=0, gotOne=false;
    else {
        item = impl.m_menuPickSilo.front(); 
        impl.m_menuPickSilo.pop_front();
        --impl.m_inputCount;
        gotOne = true;
    }
    impl.UNLOCK_silo();
    return gotOne;
}

bool Visualizer::InputSilo::takeSliderMove(int& slider, Real& value) {
    Impl& impl = updImpl(); bool gotOne;
    impl.LOCK_silo();
    if (impl.m_sliderMoveSilo.empty()) slider=0, value=NaN, gotOne=false;
    else {
        slider = impl.m_sliderMoveSilo.front().first; 
        value  = impl.m_sliderMoveSilo.front().second;
        impl.m_sliderMoveSilo.pop_front();
        --impl.m_inputCount;
        gotOne = true;
    }
    impl.UNLOCK_silo();
    return gotOne;
}

void Visualizer::InputSilo::clear() {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    impl.m_inputCount = 0;
    impl.m_keyHitSilo.clear();
    impl.m_menuPickSilo.clear();
    impl.m_sliderMoveSilo.clear();
    impl.UNLOCK_silo();
}

bool Visualizer::InputSilo::keyPressed(unsigned key, unsigned modifiers) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    impl.m_keyHitSilo.push_back(std::make_pair(key,modifiers));
    ++impl.m_inputCount;
    impl.UNLOCK_silo();
    return true;
}
bool Visualizer::InputSilo::menuSelected(int item) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    impl.m_menuPickSilo.push_back(item);
    ++impl.m_inputCount;
    impl.UNLOCK_silo();
    return true;
}
bool Visualizer::InputSilo::sliderMoved(int slider, Real value) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    impl.m_sliderMoveSilo.push_back(std::make_pair(slider,value));
    ++impl.m_inputCount;
    impl.UNLOCK_silo();
    return true;
}

