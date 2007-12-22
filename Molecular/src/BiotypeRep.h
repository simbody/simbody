#ifndef SimTK_BIOTYPEREP_H_
#define SimTK_BIOTYPEREP_H_


/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
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

#include "SimTKcommon.h"
#include "simbody/internal/Biotype.h"

namespace SimTK {
class BiotypeRep : public PIMPLImplementation<Biotype,BiotypeRep> 
{
public:
    BiotypeRep* clone() const {return new BiotypeRep(*this);}

    BiotypeRep();

    BiotypeRep(BiotypeId b,
               TinkerBiotypeId tinkerBiotypeId, 
               const Element& e,
               int v,
               const char* r, 
               const char* a, 
               Ordinality::Residue o);

    const Element&  getElement()            const {return element;}
    int             getValence()            const {return valence;}
    BiotypeId       getId()                 const {return biotypeId;}
    TinkerBiotypeId getTinkerBiotypeIfAny() const {return tinkerBiotypeIdIfAny;}

    void setTinkerBiotypeId(TinkerBiotypeId tId) 
    {
        // If its already set, that's OK
        if (tId == tinkerBiotypeIdIfAny) return;

        assert(tinkerBiotypeIdIfAny == InvalidTinkerBiotypeId);
        tinkerBiotypeIdIfAny = tId;
        assert(tinkerBiotypeIdIfAny != InvalidTinkerBiotypeId);
    }

    const String& getAtomName() const {return atomName;}
    const String& getResidueName() const {return residueName;}
    Ordinality::Residue getOrdinality() const {return ordinality;}

    // TODO make generateSelfCode protected:
    
    // emit C++ source code that can repopulate biotype data corresponding to this biotype
    std::ostream& generateSelfCode(std::ostream& os) const 
    {
        String indent("    ");
        String indent2(indent + indent);
        String indent3(indent + indent + indent);

        String quote("\"");
        String residueString = quote + residueName + quote;
        String atomString = quote + atomName + quote;

        // Create a string with a symbolic representation of the residue ordinality
        String ordString;
        switch(ordinality) {
            case Ordinality::INITIAL:
                ordString = "Ordinality::INITIAL";
                break;
            case Ordinality::ANY:
                ordString =  "Ordinality::ANY";
                break;
            case Ordinality::FINAL:
                ordString =  "Ordinality::FINAL";
                break;
            default:
                assert(false);
                break;
        }

        os << indent << "if (! Biotype::exists(";
        os << residueString << ", " << atomString << ", " << ordString << ") )" << std::endl;

        os << indent2 << "Biotype::defineTinkerBiotype(" << std::endl;

        // Tinker biotype
        if (tinkerBiotypeIdIfAny == InvalidTinkerBiotypeId)
            os << indent3 << "InvalidTinkerBiotypeId" << std::endl;
        else
            os << indent3 << "TinkerBiotypeId(" << tinkerBiotypeIdIfAny << ")" << std::endl;

        // element
        String elementName = element.getName();
        elementName[0] = toupper(elementName[0]); // capitalize element to create class name
        os << indent3 << ", Element::" << elementName << "()" << std::endl;

        // valence
        os << indent3 << ", " << valence << std::endl;

        os << indent3 << ", " << residueString << std::endl;

        os << indent3 << ", " << atomString << std::endl;

        // ordinality
        os << indent3 << ", " << ordString << std::endl;

        // close defineTinkerBiotype method call
        os << indent3 << ");" << std::endl;

        os << std::endl;

        return os;
    }

private:
    // Moved all private data into BiotypeRep from Biotype class
    BiotypeId biotypeId;
    TinkerBiotypeId tinkerBiotypeIdIfAny;
    Element  element;
    int      valence;
    // int      formalCharge;
    String residueName;
    String atomName;
    Ordinality::Residue ordinality;
};

} // namespace SimTK

#endif // SimTK_BIOTYPEREP_H_
