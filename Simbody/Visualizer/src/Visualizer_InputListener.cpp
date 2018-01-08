/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/** @file
This file contains the implementations of Visualizer::InputListener and any
built-in concrete InputListeners, currently just Visualizer::InputSilo. **/

#include "simbody/internal/common.h"
#include "simbody/internal/Visualizer_InputListener.h"

#include <deque>
#include <utility>
#include <condition_variable>

using namespace SimTK;

//==============================================================================
//                              INPUT SILO IMPL
//==============================================================================
/* This is the private implementation object contained in a 
Visualizer::InputSilo handle. */

class Visualizer::InputSilo::Impl {
public:
    std::deque<std::pair<unsigned,unsigned> >   m_keyHitSilo;
    std::deque<std::pair<int,int> >             m_menuPickSilo;
    std::deque<std::pair<int,Real> >            m_sliderMoveSilo;
    unsigned                                    m_inputCount {0};

    mutable std::mutex m_siloMutex;
    // signal on 1st input, any silo
    mutable std::condition_variable m_someInputAvailable;
    // signal on 1st key in silo
    mutable std::condition_variable m_keyHitAvailable;
    // signal on 1st pick in silo
    mutable std::condition_variable m_menuPickAvailable;
    // signal on 1st move in silo
    mutable std::condition_variable m_sliderMoveAvailable;
};

//==============================================================================
//                                INPUT SILO
//==============================================================================
Visualizer::InputSilo::InputSilo() {m_impl = new InputSilo::Impl();}
Visualizer::InputSilo::~InputSilo() {delete m_impl;}

// Fast ... no locking required for this check.
bool Visualizer::InputSilo::isAnyUserInput() const 
{   return getImpl().m_inputCount != 0; }

// Hang until "anyInput" condition is signaled.
void Visualizer::InputSilo::waitForAnyUserInput() const {
    if (isAnyUserInput()) return;
    const Impl& impl = getImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    impl.m_someInputAvailable.wait(lock, [this] {return isAnyUserInput();});
    lock.unlock();
}

bool Visualizer::InputSilo::takeKeyHit(unsigned& key, unsigned& modifiers) {
    if (!isAnyUserInput()) return false;
    Impl& impl = updImpl(); bool gotOne;
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    if (impl.m_keyHitSilo.empty()) key=0, modifiers=0, gotOne=false;
    else {
        key       = impl.m_keyHitSilo.front().first; 
        modifiers = impl.m_keyHitSilo.front().second;
        impl.m_keyHitSilo.pop_front();
        --impl.m_inputCount;
        gotOne = true;
    }
    lock.unlock();
    return gotOne;
}

void Visualizer::InputSilo::waitForKeyHit(unsigned& key, unsigned& modifiers) {
    Impl& impl = updImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    impl.m_keyHitAvailable.wait(lock, [&] {return impl.m_keyHitSilo.size();});
    key       = impl.m_keyHitSilo.front().first; 
    modifiers = impl.m_keyHitSilo.front().second;
    impl.m_keyHitSilo.pop_front();
    --impl.m_inputCount;
    lock.unlock();
}


bool Visualizer::InputSilo::takeMenuPick(int& menuId, int& item) {
    if (!isAnyUserInput()) return false;
    Impl& impl = updImpl(); bool gotOne;
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    if (impl.m_menuPickSilo.empty()) item=0, gotOne=false;
    else {
        menuId = impl.m_menuPickSilo.front().first; 
        item   = impl.m_menuPickSilo.front().second;
        impl.m_menuPickSilo.pop_front();
        --impl.m_inputCount;
        gotOne = true;
    }
    lock.unlock();
    return gotOne;
}

void Visualizer::InputSilo::waitForMenuPick(int& menuId, int& item) {
    Impl& impl = updImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    impl.m_menuPickAvailable.wait(lock,
            [&] {return impl.m_menuPickSilo.size();});
    
    menuId = impl.m_menuPickSilo.front().first; 
    item   = impl.m_menuPickSilo.front().second;
    impl.m_menuPickSilo.pop_front();
    --impl.m_inputCount;
    lock.unlock();
}

bool Visualizer::InputSilo::takeSliderMove(int& slider, Real& value) {
    if (!isAnyUserInput()) return false;
    Impl& impl = updImpl(); bool gotOne;
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    if (impl.m_sliderMoveSilo.empty()) slider=0, value=NaN, gotOne=false;
    else {
        slider = impl.m_sliderMoveSilo.front().first; 
        value  = impl.m_sliderMoveSilo.front().second;
        impl.m_sliderMoveSilo.pop_front();
        --impl.m_inputCount;
        gotOne = true;
    }
    lock.unlock();
    return gotOne;
}

void Visualizer::InputSilo::waitForSliderMove(int& slider, Real& value) {
    Impl& impl = updImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    impl.m_sliderMoveAvailable.wait(lock, 
            [&] {return impl.m_sliderMoveSilo.size();});
    slider = impl.m_sliderMoveSilo.front().first; 
    value  = impl.m_sliderMoveSilo.front().second;
    impl.m_sliderMoveSilo.pop_front();
    --impl.m_inputCount;
    lock.unlock();
}

void Visualizer::InputSilo::clear() {
    Impl& impl = updImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    impl.m_inputCount = 0;
    impl.m_keyHitSilo.clear();
    impl.m_menuPickSilo.clear();
    impl.m_sliderMoveSilo.clear();
    lock.unlock();
}

bool Visualizer::InputSilo::keyPressed(unsigned key, unsigned modifiers) {
    Impl& impl = updImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    impl.m_keyHitSilo.push_back(std::make_pair(key,modifiers));
    if (++impl.m_inputCount == 1) {
        // in case someone was waiting for any input
        impl.m_someInputAvailable.notify_one();
    }
    if (impl.m_keyHitSilo.size() == 1) {
        // a key hit is now available
        impl.m_keyHitAvailable.notify_one();
    }
    lock.unlock();
    return true;
}
bool Visualizer::InputSilo::menuSelected(int menu, int item) {
    Impl& impl = updImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    impl.m_menuPickSilo.push_back(std::make_pair(menu,item));
    if (++impl.m_inputCount == 1) {
        // in case someone was waiting for any input
        impl.m_someInputAvailable.notify_one();
    }
    if (impl.m_menuPickSilo.size() == 1) {
        // a menu pick is now available
        impl.m_menuPickAvailable.notify_one();
    }
    lock.unlock();
    return true;
}

// We optimize here for the common case that the same slider is moving
// for a while -- in that case we just keep the most recent position.
bool Visualizer::InputSilo::sliderMoved(int slider, Real value) {
    Impl& impl = updImpl();
    std::unique_lock<std::mutex> lock(impl.m_siloMutex);
    std::deque<std::pair<int,Real> >& silo = impl.m_sliderMoveSilo;
    if (!silo.empty() && silo.front().first == slider)
        silo.front().second = value; // just replace the value; count unchanged
    else {
        silo.push_back(std::make_pair(slider,value));
        if (++impl.m_inputCount == 1) {
            // in case someone was waiting for any input
            impl.m_someInputAvailable.notify_one();
        }
        if (silo.size() == 1) {
            // a slider move is now available
            impl.m_sliderMoveAvailable.notify_one();
        }
    }
    lock.unlock();
    return true;
}

