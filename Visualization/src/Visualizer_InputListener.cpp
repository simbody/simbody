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
    Impl() : m_inputCount(0) {
        pthread_mutex_init(&m_siloLock,0);
        pthread_cond_init(&m_someInputAvailable,0);
        pthread_cond_init(&m_keyHitAvailable,0);
        pthread_cond_init(&m_menuPickAvailable,0);
        pthread_cond_init(&m_sliderMoveAvailable,0);
    }
    ~Impl() {
        pthread_cond_destroy(&m_sliderMoveAvailable);
        pthread_cond_destroy(&m_menuPickAvailable);
        pthread_cond_destroy(&m_keyHitAvailable);
        pthread_cond_destroy(&m_someInputAvailable);
        pthread_mutex_destroy(&m_siloLock);
    }

    void LOCK_silo()     const {pthread_mutex_lock(&m_siloLock);}
    void UNLOCK_silo()   const {pthread_mutex_unlock(&m_siloLock);}

    void WAIT_anyInput() const {pthread_cond_wait(&m_someInputAvailable, 
                                                  &m_siloLock);}
    void POST_anyInput() const {pthread_cond_signal(&m_someInputAvailable);}

    void WAIT_keyHit()   const {pthread_cond_wait(&m_keyHitAvailable, 
                                                  &m_siloLock);}
    void POST_keyHit()   const {pthread_cond_signal(&m_keyHitAvailable);}

    void WAIT_menuPick() const {pthread_cond_wait(&m_menuPickAvailable, 
                                                  &m_siloLock);}
    void POST_menuPick() const {pthread_cond_signal(&m_menuPickAvailable);}

    void WAIT_sliderMove() const {pthread_cond_wait(&m_sliderMoveAvailable, 
                                                    &m_siloLock);}
    void POST_sliderMove() const {pthread_cond_signal(&m_sliderMoveAvailable);}

    std::deque<std::pair<unsigned,unsigned> >   m_keyHitSilo;
    std::deque<std::pair<int,int> >             m_menuPickSilo;
    std::deque<std::pair<int,Real> >            m_sliderMoveSilo;
    unsigned                                    m_inputCount;

    mutable pthread_mutex_t m_siloLock;
    mutable pthread_cond_t  m_someInputAvailable;  // signal on 1st input, any silo
    mutable pthread_cond_t  m_keyHitAvailable;     // signal on 1st key in silo
    mutable pthread_cond_t  m_menuPickAvailable;   // signal on 1st pick in silo
    mutable pthread_cond_t  m_sliderMoveAvailable; // signal on 1st move in silo
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
    impl.LOCK_silo();
    while (!isAnyUserInput()) {impl.WAIT_anyInput();}
    impl.UNLOCK_silo();
}

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

void Visualizer::InputSilo::waitForKeyHit(unsigned& key, unsigned& modifiers) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    while (!impl.m_keyHitSilo.size()) 
        impl.WAIT_keyHit();
    key       = impl.m_keyHitSilo.front().first; 
    modifiers = impl.m_keyHitSilo.front().second;
    impl.m_keyHitSilo.pop_front();
    --impl.m_inputCount;
    impl.UNLOCK_silo();
}


bool Visualizer::InputSilo::takeMenuPick(int& menuId, int& item) {
    Impl& impl = updImpl(); bool gotOne;
    impl.LOCK_silo();
    if (impl.m_menuPickSilo.empty()) item=0, gotOne=false;
    else {
        menuId = impl.m_menuPickSilo.front().first; 
        item   = impl.m_menuPickSilo.front().second;
        impl.m_menuPickSilo.pop_front();
        --impl.m_inputCount;
        gotOne = true;
    }
    impl.UNLOCK_silo();
    return gotOne;
}

void Visualizer::InputSilo::waitForMenuPick(int& menuId, int& item) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    while (!impl.m_menuPickSilo.size()) 
        impl.WAIT_menuPick();
    menuId = impl.m_menuPickSilo.front().first; 
    item   = impl.m_menuPickSilo.front().second;
    impl.m_menuPickSilo.pop_front();
    --impl.m_inputCount;
    impl.UNLOCK_silo();
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

void Visualizer::InputSilo::waitForSliderMove(int& slider, Real& value) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    while (!impl.m_sliderMoveSilo.size()) 
        impl.WAIT_sliderMove();
    slider = impl.m_sliderMoveSilo.front().first; 
    value  = impl.m_sliderMoveSilo.front().second;
    impl.m_sliderMoveSilo.pop_front();
    --impl.m_inputCount;
    impl.UNLOCK_silo();
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
    if (++impl.m_inputCount == 1)
        impl.POST_anyInput(); // in case someone was waiting for any input
    if (impl.m_keyHitSilo.size() == 1)
        impl.POST_keyHit();   // a key hit is now available
    impl.UNLOCK_silo();
    return true;
}
bool Visualizer::InputSilo::menuSelected(int menu, int item) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    impl.m_menuPickSilo.push_back(std::make_pair(menu,item));
    if (++impl.m_inputCount == 1)
        impl.POST_anyInput(); // in case someone was waiting for any input
    if (impl.m_menuPickSilo.size() == 1)
        impl.POST_menuPick(); // a menu pick is now available
    impl.UNLOCK_silo();
    return true;
}
bool Visualizer::InputSilo::sliderMoved(int slider, Real value) {
    Impl& impl = updImpl();
    impl.LOCK_silo();
    impl.m_sliderMoveSilo.push_back(std::make_pair(slider,value));
    if (++impl.m_inputCount == 1)
        impl.POST_anyInput();   // in case someone was waiting for any input
    if (impl.m_sliderMoveSilo.size() == 1)
        impl.POST_sliderMove(); // a slider move is now available
    impl.UNLOCK_silo();
    return true;
}

