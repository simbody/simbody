#ifndef SimTK_SIMBODY_BIOTYPE_H_
#define SimTK_SIMBODY_BIOTYPE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK SIMBODY                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Christopher Bruns                                *
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


#include "simbody/internal/common.h"
#include "simbody/internal/Element.h"

namespace SimTK {

namespace Ordinality {
    /// Enumeration to indicate whether a residue is at the beginning, middle, or end of a polymer chain.
    enum Residue {
        Any = 1, ///< Indicates either a residue in the middle (i.e. not at the ends) of a polymer chain, or that the ordinality is unimportant
        Initial = 2, ///<Indicates a residue at the beginning of a polymer chain
        Final = 3}; ///<Indicates a residue at the end of a polymer chain
}

// Biotype is a hook that will be used to look up molecular
// force field specific parameters for an atom type
SimTK_DEFINE_AND_EXPORT_UNIQUE_INDEX_TYPE(SimTK_SIMBODY_EXPORT,BiotypeIndex);
SimTK_DEFINE_AND_EXPORT_UNIQUE_INDEX_TYPE(SimTK_SIMBODY_EXPORT,TinkerBiotypeIndex);

class Biotype;
class BiotypeRep;

// Avoid instantiating template class except in one particular library compilation unit
#ifndef DO_INSTANTIATE_BIOTYPE_PIMPL_HANDLE
    extern template class SimTK_SIMBODY_EXPORT PIMPLHandle<Biotype, BiotypeRep>;
#endif

class SimTK_SIMBODY_EXPORT Biotype : public PIMPLHandle<SimTK_SIMBODY_EXPORT Biotype,BiotypeRep> 
{
public:

    // Noble gases
    static const Biotype& Argon();

    static const Biotype& MethaneH();
    static const Biotype& MethaneC();

    static const Biotype& EthaneH();
    static const Biotype& EthaneC();

    static const Biotype& SerineN();
    static const Biotype& SerineHN();
    static const Biotype& SerineCA();
    static const Biotype& SerineHA();

    static const Biotype& SerineC();
    static const Biotype& SerineO();
    static const Biotype& SerineCB();
    static const Biotype& SerineHB();
    static const Biotype& SerineOG();
    static const Biotype& SerineHG();

    const Element&  getElement() const;
    int             getValence() const;
    BiotypeIndex       getIndex() const;
    TinkerBiotypeIndex getTinkerBiotypeIfAny() const;

    Biotype& setTinkerBiotypeIndex(TinkerBiotypeIndex tIx);

    static void initializePopularBiotypes();

    static const Biotype& get(BiotypeIndex biotypeIndex);
    static const Biotype& get(TinkerBiotypeIndex tinkerBiotypeIndex);
    static const Biotype& get(const char* residueName, 
                              const char* atomName, 
                              Ordinality::Residue ordinality = Ordinality::Any);

    static Biotype& upd(BiotypeIndex biotypeIndex);
    static Biotype& upd(TinkerBiotypeIndex tinkerBiotypeIndex);
    static Biotype& upd(const char* residueName, 
                              const char* atomName, 
                              Ordinality::Residue ordinality = Ordinality::Any);

    static bool exists(const char* residueName, 
                       const char* atomName, 
                       Ordinality::Residue ordinality = Ordinality::Any);

    static bool exists(BiotypeIndex biotypeIndex);

    static BiotypeIndex defineBiotype(const Element& element,
                                   int valence,
                                   const char* residueName, 
                                   const char* atomName, 
                                   Ordinality::Residue ordinality = Ordinality::Any)
    {
        return defineTinkerBiotype(InvalidTinkerBiotypeIndex, 
                                   element,
                                   valence,
                                   residueName, 
                                   atomName, 
                                   ordinality);
    }

    static BiotypeIndex defineTinkerBiotype(TinkerBiotypeIndex tinkerBiotypeIndex, 
                                         const Element& element,
                                         int valence,
                                         const char* residueName, 
                                         const char* atomName, 
                                         Ordinality::Residue ordinality = Ordinality::Any); 

    static std::ostream& generateAllBiotypeCode(std::ostream& os); 

    Biotype();

    const String& getAtomName() const;
    const String& getResidueName() const;
    Ordinality::Residue getOrdinality() const;


// TODO make generateSelfCode protected:
    
    // emit C++ source code that can repopulate biotype data corresponding to this biotype
    std::ostream& generateSelfCode(std::ostream& os) const;

private:

    Biotype( BiotypeIndex biotypeIndex, 
             TinkerBiotypeIndex tinkerBiotypeIndex, 
             const Element& element,
             int valence,
             const char* residueName, 
             const char* atomName, 
             Ordinality::Residue ordinality);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_BIOTYPE_H_
