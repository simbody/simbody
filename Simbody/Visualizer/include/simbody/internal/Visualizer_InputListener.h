#ifndef SimTK_SIMBODY_VISUALIZER_INPUT_LISTENER_H_
#define SimTK_SIMBODY_VISUALIZER_INPUT_LISTENER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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
This defines the InputListener class that is internal to the Visualizer, and
the InputSilo class derived from it. **/

#include "simbody/internal/common.h"
#include "simbody/internal/Visualizer.h"

namespace SimTK {

//==============================================================================
//                             INPUT LISTENER
//==============================================================================
/** This abstract class defines methods to be called when the Visualizer
reports user activity back to the simulation process.\ Derive a concrete event 
listener whose methods take appropriate actions when event of interest occur.

You only need to implement the methods you care about. The Visualizer provides
an InputListener called InputSilo that is particularly useful for dealing with
user input that is intended to affect the progress of a running simulation.
@see Visualizer::InputSilo **/
class SimTK_SIMBODY_EXPORT Visualizer::InputListener {
public:
/** These represent modifications to the character that is passed into the 
keyPressed() method, including whether any of Shift/Control/Alt were down and 
whether a special non-ASCII code is being supplied, such as is required for an 
arrow key. These values are or'ed together to form the second argument of the 
keyPressed() method. **/
enum Modifier {
    ShiftIsDown   = 0x01, ///< Shift (left or right)
    ControlIsDown = 0x02, ///< Ctrl (left or right)
    AltIsDown     = 0x04, ///< Alt (left or right)
    IsSpecialKey  = 0xC0  ///< Special non-ASCII keycode being used
};

static const unsigned SpecialKeyOffset = 0x100; // Added to each code

/** These are the special keys that the Visualizer may report via the 
keyPressed() method. All other keys are considered "ordinary". **/
enum KeyCode {
    KeyControlC     = 3,            // some notable ASCII codes
    KeyBeep         = 7,
    KeyBackspace    = 8,            
    KeyTab          = 9,
    KeyLF           = 10,
    KeyReturn       = 13,
    KeyEnter        = KeyReturn,
    KeyEsc          = 27,
    KeyDelete       = 127,

    KeyF1  = SpecialKeyOffset + 1,  // function keys
    KeyF2  = SpecialKeyOffset + 2,
    KeyF3  = SpecialKeyOffset + 3,
    KeyF4  = SpecialKeyOffset + 4,
    KeyF5  = SpecialKeyOffset + 5,
    KeyF6  = SpecialKeyOffset + 6,
    KeyF7  = SpecialKeyOffset + 7,
    KeyF8  = SpecialKeyOffset + 8,
    KeyF9  = SpecialKeyOffset + 9,
    KeyF10 = SpecialKeyOffset + 10,
    KeyF11 = SpecialKeyOffset + 11,
    KeyF12 = SpecialKeyOffset + 12,

    KeyLeftArrow    = SpecialKeyOffset + 100,  // directional keys
    KeyUpArrow      = SpecialKeyOffset + 101,
    KeyRightArrow   = SpecialKeyOffset + 102,
    KeyDownArrow    = SpecialKeyOffset + 103,
    KeyPageUp       = SpecialKeyOffset + 104,
    KeyPageDown     = SpecialKeyOffset + 105,
    KeyHome         = SpecialKeyOffset + 106,
    KeyEnd          = SpecialKeyOffset + 107,
    KeyInsert       = SpecialKeyOffset + 108
};
    
/** Destructor is virtual; be sure to override it if you need to clean up. **/
virtual ~InputListener() {}

/** This method is called when a user hits a keyboard key in the Visualizer 
window, unless that key is being intercepted by the Visualizer for its own 
purposes. Ordinary ASCII characters 0-127 are represented by their own values; 
special keys like arrows and function keys are mapped to unique values > 255. 
You can check whether \a modifers & IsSpecialKey is true if you care; otherwise
just mix the ordinary and special codes in a case statement. You can tell if 
any or all of Shift/Control/Alt were depressed when the key was hit by checking
the \a modifiers. Note that for an ordinary capital letter you'll get the ASCII
code for the capital as well as an indication that the Shift key was down. If 
caps lock was down you'll get the capital letter but no Shift modifier.
@param[in]  key         The ASCII character or special key code.
@param[in]  modifiers   Whether Shift,Ctrl,Alt are down or key is special.
@return Return \c true if you have handled this key press and don't want any 
subsequent listeners called. **/
virtual bool keyPressed(unsigned key, unsigned modifiers) {return false;}

/** The user has clicked one of the menu items you defined; here is
the integer value you specified when you defined it.
@param[in]  menu        The id number of the menu in use.
@param[in]  item        The menu item number that the user selected.
@return Return \c true if you have handled this menu click and don't
want any subsequent listeners called. **/
virtual bool menuSelected(int menu, int item) {return false;}

/** The user has moved one of the sliders you defined; here is the integer 
value you specified when you defined it, and the new value of the slider.
@param[in]  slider      The id number of the slider that moved.
@param[in]  value       The new value represented by this slider position.
@return Return \c true if you have handled this move and don't want any 
subsequent listeners called. **/
virtual bool sliderMoved(int slider, Real value) {return false;}
};



//==============================================================================
//                              INPUT SILO
//==============================================================================
/** This pre-built InputListener is extremely useful for processing user
input that is intended to affect a running simulation. The idea is that this 
object saves up all the user input in a set of "silos", which are 
first-in-first-out (FIFO) queues. The simulation periodically checks ("polls") 
to see if there is anything in the silos that needs processing, pulling off one
user input at a time until they have all been consumed. This eliminates any 
need for tricky asynchronous handling of user input, and all thread 
synchronization issues are handled invisibly. 

You can also request to wait quietly until some input arrives, which is useful
when you can't proceed without some instruction from the user that you expect
to get through the visualizer.

When the InputSilo receives user input through one of the InputListener methods
it implements, it return \c true indicating that it has processed the input and
that no further InputListeners should be called. So if you have other 
InputListeners that you would like to have called, be sure to add them to the
Visulizer \e prior to adding an InputSilo, which is the last refuge for 
unwanted user input.

Here's how you can use this:

@code
MultibodySystem system;
// ... build system

// Set up a Visualizer to run in real time mode, and give it an
// InputSilo to gather user input.
Visualizer viz(system);
viz.setMode(Visualizer::RealTime);
InputSilo* silo = new InputSilo;
viz.addInputListener(silo);

// You create a PeriodicEventHandler to poll the input. Note that the interval
// you choose determines how responsive the simulation will be to user input,
// but it also limits the maximum step size that the integrator can take.
system.addEventHandler
    (new MyUserInputHandler(*silo, 0.1)); // check every 100ms 

// Then in MyUserInputHandler::handleEvent(...):
while (silo.isAnyUserInput()) {
    while (silo.takeKeyHit(key,modifier)) {
        // Process the key that was hit
    }
    while (silo.takeMenuPick(which, item)) {
        // Process the picked menu item
    }
    while (silo.takeSliderMove(which, value)) {
        // Process the new value for slider "which"
    }
}
@endcode

If you want to wait until some input arrives, create the InputSilo and add
it to the Visualizer as above, then in your main program (that is, not
in the Handler) use code like this:
@code
std::cout << "Hit ENTER in visualizer to continue ...\n";
unsigned key, modifiers;
do {silo->waitForKeyHit(key,modifiers);}
while (key != Visualizer::InputListener::KeyEnter);
@endcode
Similar methods are available for all the different input types, and you
can also wait on the arrival of \e any input.

<h3>Implementation</h3>

The InputSilo implementations of the InputListener methods are called from the 
Visualizer's listener thread, which is a different thread than the one that
is simultaneously running the simulation. The internal silos are double-ended
queues (deques) that allow inputs to be pushed onto one end and pulled off
the other, so that they can be consumed in FIFO order. There is a single mutex
lock associated with \e all the silos together, and the lock must be held while
anything is pushed onto or pulled off of any one of the silos.

Each of the methods for getting the input out of the silos is called from the
simulation thread, which must obtain the lock before removing anything, thus 
safely synchronizing the listener and simulation threads. 

A count is maintained of the total number of items in all the silos. It is
incremented only when the listener thread holds the lock and adds something
to a silo; it is decremented only when the simulation thread holds the lock
and pulls something from a silo. The count may be examined without locking; it
will have a value that was recently correct and can thus be used for a very
fast check on whether there is likely to be any input worth holding a lock for;
the isAnyUserInput() method returns \c true when the count is non-zero. It may
occasionally return zero in cases where there is input, but only if that input
just arrived so you can safely pick it up on the next poll. 

When possible we optimize for the case where many inputs arrive from the
same device by just keeping the most recent value. That applies to slider
and mouse moves. **/
class SimTK_SIMBODY_EXPORT Visualizer::InputSilo
:   public Visualizer::InputListener {
public:
/** Default construction is all that is needed; there are no options. **/
InputSilo();
/** Throws away any unprocessed input. **/
~InputSilo();

/** This is a very fast test that does not require locking; you don't have to 
use this but it is a good idea to do so. **/
bool isAnyUserInput() const;

/** This will wait quietly until the user has provided some input to the
visualizer.\ Any kind of input will terminate the wait; you'll have
to look to see what it was. **/
void waitForAnyUserInput() const;

/** This will return user key hits until they have all been consumed, in the 
same order they were received. The \a key and \a modifiers values are those that
were provided to our implementation of the InputListener::keyPressed() method. 
@param[out]         key         
    The key code for the key that was hit. See InputListener::KeyCode for
    interpretation.
@param[out]         modifiers    
    Status of Shift,Ctrl,Alt and "special" key code. See InputListener::Modifier
    for interpretation. 
@return \c true if a key and modifiers have been returned; \c false if the 
    character silo is now empty in which case both \a key and \a modifiers will
    be set to zero. **/
bool takeKeyHit(unsigned& key, unsigned& modifiers);

/** Same as takeKeyHit() except that if there is no key hit input available
it waits until there is, then returns the first one (which is removed from
the silo just as takeKeyHit() would do. The behavior is like calling
waitForAnyUserInput() repeatedly until takeKeyHit() returns \c true.
@see takeKeyHit(), waitForAnyUserInput() **/
void waitForKeyHit(unsigned& key, unsigned& modifiers);

/** This will return user menu picks until they have all been consumed, in the
same order they were received. The \a item value returned is the value that was
provided to our implementation of the InputListener::menuSelected() method. 
@param[out]         menu         
    The id number of the menu that was selected. This is the value that was
    assigned to this menu in the Visualizer::addMenu() call.
@param[out]         item         
    The menu item number for the entry that the user selected. This is the 
    number that was assigned at the time the menu was added via the 
    Visualizer::addMenu() method.
@return \c true if a menu item number has been returned; \c false if the menu 
    pick silo is now empty in which case \a item will be set to zero. **/
bool takeMenuPick(int& menu, int& item);

/** Same as takeMenuPick() except that if there is no menu pick input available
it waits until there is, then returns the first one (which is removed from
the silo just as takeMenuPick() would do. The behavior is like calling
waitForAnyUserInput() repeatedly until takeMenuPick() returns \c true.
@see takeMenuPick(), waitForAnyUserInput() **/
void waitForMenuPick(int& menu, int& item);

/** This will return user changes to slider positions until they have all been
consumed, in the same order they were received. The \a slider and \a value 
returns are those that were provided to our implementation of the 
InputListener::sliderMoved() method. 
@param[out]         slider         
    The id number of the slider that was moved. This is the value that was
    assigned to this slider in the Visualizer::addSlider() call.
@param[out]         value    
    This is the new value associated with the slider position to which the user
    moved it.
@return \c true if a slider move has been returned; \c false if the slider move
    silo is now empty in which case \a which will be set to zero and \a value 
    will be NaN. **/
bool takeSliderMove(int& slider, Real& value);

/** Same as takeSliderMove() except that if there is no slider move input 
available it waits until there is, then returns the first one (which is removed
from the silo just as takeSliderMove() would do. The behavior is like calling
waitForAnyUserInput() repeatedly until takeSliderMove() returns \c true.
@see takeSliderMove(), waitForAnyUserInput() **/
void waitForSliderMove(int& slider, Real& value);

/** Throw away any pending unprocessed input of all types. **/
void clear();

//------------------------------------------------------------------------------
                                 private:
// Each of these will return true to the Visualizer's listener thread, meaning 
// that the input will be absorbed and subsequent listeners (if any) will not 
// be called.
virtual bool keyPressed(unsigned key, unsigned modifiers) override;
virtual bool menuSelected(int menu, int item) override;
virtual bool sliderMoved(int slider, Real value) override;

class Impl;
const Impl& getImpl() const {assert(m_impl); return *m_impl;}
Impl&       updImpl()       {assert(m_impl); return *m_impl;}

Impl*       m_impl;   // the lone data member in this class
};

} // namespace SimTK

#endif // SimTK_SIMBODY_VISUALIZER_INPUT_LISTENER_H_
