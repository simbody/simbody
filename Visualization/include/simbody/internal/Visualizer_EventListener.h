#ifndef SimTK_SIMBODY_VISUALIZER_EVENT_LISTENER_H_
#define SimTK_SIMBODY_VISUALIZER_EVENT_LISTENER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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
This defines the EventListener class that is internal to the Visualizer. **/

#include "simbody/internal/common.h"
#include "simbody/internal/Visualizer.h"

namespace SimTK {

/** This abstract class defines methods to be called when the Visualizer
reports user activity back to the simulation process. Derive a concrete
event listener whose methods take appropriate actions when event of interest
occur. You only need to implement the methods you care about. **/
class SimTK_SIMBODY_EXPORT Visualizer::EventListener {
public:
    /** These represent modifications to the character that is passed into
    the keyPressed() method, including whether any of Shift/Control/Alt were
    down and whether a special non-ASCII code is being supplied, such as
    is required for an arrow key. These values are or'ed together to form
    the second argument of the keyPressed() method. **/
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
    
    /** Destructor is virtual; be sure to override it if you need to 
    clean up. **/
    virtual ~EventListener() {}

    /** This method is called when a user hits a keyboard key
    in the Visualizer window, unless that key is being intercepted
    by the Visualizer for its own purposes. Ordinary ASCII characters
    0-127 are represented by their own values; special keys like arrows
    and function keys are mapped to unique values > 255. You can check
    whether \a modifers & IsSpecialKey is true if you care; otherwise
    just mix the ordinary and special codes in a case statement. You can
    tell if any or all of Shift/Control/Alt were depressed when the key
    was hit by checking the \a modifiers. Note that for an ordinary 
    capital letter you'll get the ASCII code for the capital as well
    as an indication that the Shift key was down. If caps lock was down
    you'll get the capital letter but no Shift modifier. 
    @return Return true if you have handled this key press and don't
    want any subsequent listeners called. **/
    virtual bool keyPressed(unsigned key, unsigned modifiers) {return false;}

    /** The user has clicked one of the menu items you defined; here is
    the integer value you specified when you defined it.    
    @return Return true if you have handled this menu click and don't
    want any subsequent listeners called. **/
    virtual bool menuSelected(int item) {return false;}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_VISUALIZER_EVENT_LISTENER_H_
