/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

// This is the one file where we want the handle instantiated
#define DO_INSTANTIATE_BIOTYPE_PIMPL_HANDLE
#include "simbody/internal/Biotype.h"
#undef DO_INSTANTIATE_BIOTYPE_PIMPL_HANDLE

#include "BiotypeRep.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#include "simbody/internal/common.h"
#include <map>
#include <algorithm>
#include <cctype>

using namespace std;

namespace SimTK {

// Store biotype data in this compilation unit, rather than in the biotype class

static void initializePopularBiotypes();

static std::map<BiotypeId, Biotype> biotypesById;
static std::map<TinkerBiotypeId, BiotypeId> biotypeIdsByTinkerId;
static BiotypeId nextUnusedBiotypeId(1);

class BiotypeKey {
public:
    BiotypeKey();
    BiotypeKey(const char* r, const char* a, Ordinality::Residue o)
        : residueName(regularizeString(r)), atomName(regularizeString(a)), ordinality(o)
    {}

    bool operator<(const BiotypeKey& other) const {
        if (residueName == other.residueName) {
            if (atomName == other.atomName)
                return (ordinality < other.ordinality);
            else
                return (atomName < other.atomName);
        }
        else {
            return (residueName < other.residueName);
        }
    }

    // convert to lower case and trim leading and trailing spaces
    static String regularizeString(const char* inputString) 
    {
        String s(inputString);

        // Trim back end first, otherwise back end indices may not be correct
        int end = s.find_last_not_of(" ");
        if ( end < (int)(s.size() - 1) )
            s.replace(++end, s.size(), "");

        // Trim front end
        int start = s.find_first_not_of(" ");
        if (start > 0)
            s.replace(0, start, "");

        // lower case
        std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) std::tolower);

        return s;
    }

    String residueName;
    String atomName;
    Ordinality::Residue ordinality;
};

static bool tinkerIdExists(TinkerBiotypeId tId) 
{
    return (biotypeIdsByTinkerId.find(tId) != biotypeIdsByTinkerId.end());
}

static bool biotypeIdExists(BiotypeId bId) 
{
    return (biotypesById.find(bId) != biotypesById.end());
}

static std::map<BiotypeKey, BiotypeId> biotypeIdsByKey;

/* static */ BiotypeId Biotype::defineTinkerBiotype(
        TinkerBiotypeId tinkerBiotypeId, 
        const Element& element,
        int valence,
        const char* residueName, 
        const char* atomName, 
        Ordinality::Residue ordinality)
{
    assert(! Biotype::exists(residueName, atomName, ordinality) );

    BiotypeId biotypeId = nextUnusedBiotypeId;
    ++nextUnusedBiotypeId;
    assert(biotypesById.find(biotypeId) == biotypesById.end());

    BiotypeKey key(residueName, atomName, ordinality);

    biotypeIdsByKey[key] = biotypeId;
    biotypesById[biotypeId] = Biotype(
        biotypeId,
        tinkerBiotypeId,
        element,
        valence,
        residueName,
        atomName,
        ordinality
        );

    if (tinkerBiotypeId != InvalidTinkerBiotypeId) {
        assert(! tinkerIdExists(tinkerBiotypeId) );
        biotypeIdsByTinkerId[tinkerBiotypeId] = biotypeId;
        assert( tinkerIdExists(tinkerBiotypeId) );
    }

    assert( Biotype::exists(residueName, atomName, ordinality) );
    assert( biotypesById.find(biotypeId) != biotypesById.end() );

    return biotypeId;
}

static bool biotypeExists(const BiotypeKey& key) {
    return biotypeIdsByKey.find(key) != biotypeIdsByKey.end();
}

/* static */ bool Biotype::exists(const char* residueName, 
                   const char* atomName, 
                   Ordinality::Residue ordinality) {
    BiotypeKey key1(residueName, atomName, ordinality);
    BiotypeKey key2(residueName, atomName, Ordinality::ANY);
    return (biotypeExists(key1) || biotypeExists(key2));
}

/* static */ bool Biotype::exists(BiotypeId biotypeId) {
    return biotypeIdExists(biotypeId);
}

////////////////
// BiotypeRep //
////////////////

BiotypeRep::BiotypeRep() {
}

BiotypeRep::BiotypeRep(BiotypeId b,
                       TinkerBiotypeId tinkerBiotypeId, 
                       const Element& e,
                       int v,
                       const char* r, 
                       const char* a, 
                       Ordinality::Residue o)
     : biotypeId(b)
     , tinkerBiotypeIdIfAny(tinkerBiotypeId)
     , element(e)
     , valence(v)
     , residueName(r)
     , atomName(a)
     , ordinality(o)
{}


///////////////
//  Biotype  //
///////////////

// These are all the constants of the Enumeration type Biotype.
///*static*/ const Biotype Biotype::Argon(Element::Argon(), );
//
///*static*/ const Biotype Biotype::MethaneH;
///*static*/ const Biotype Biotype::MethaneC;
//
///*static*/ const Biotype Biotype::EthaneH;
///*static*/ const Biotype Biotype::EthaneC;
//
///*static*/ const Biotype Biotype::SerineN;
///*static*/ const Biotype Biotype::SerineHN;
///*static*/ const Biotype Biotype::SerineCA;
///*static*/ const Biotype Biotype::SerineHA;
///*static*/ const Biotype Biotype::SerineC;
///*static*/ const Biotype Biotype::SerineO;
///*static*/ const Biotype Biotype::SerineCB;
///*static*/ const Biotype Biotype::SerineHB;
///*static*/ const Biotype Biotype::SerineOG;
///*static*/ const Biotype Biotype::SerineHG;

Biotype::Biotype() 
    : HandleBase(new BiotypeRep())
{}

Biotype::Biotype(BiotypeId b,
                 TinkerBiotypeId tinkerBiotypeId, 
                 const Element& e,
                 int v,
                 const char* r, 
                 const char* a, 
                 Ordinality::Residue o)
    : HandleBase(new BiotypeRep(b, tinkerBiotypeId, e, v, r, a, o))
{}

const Element&  Biotype::getElement() const
{
    return getImpl().getElement();
}

int             Biotype::getValence() const
{
    return getImpl().getValence();
}

BiotypeId       Biotype::getId() const
{
    return getImpl().getId();
}

TinkerBiotypeId Biotype::getTinkerBiotypeIfAny() const
{
    return getImpl().getTinkerBiotypeIfAny();
}

Biotype& Biotype::setTinkerBiotypeId(TinkerBiotypeId tId)
{
    updImpl().setTinkerBiotypeId(tId);
    return *this;
}

const String& Biotype::getAtomName() const
{
    return getImpl().getAtomName();
}

const String& Biotype::getResidueName() const
{
    return getImpl().getResidueName();
}

Ordinality::Residue Biotype::getOrdinality() const
{
    return getImpl().getOrdinality();
}

std::ostream& Biotype::generateSelfCode(std::ostream& os) const
{
    return getImpl().generateSelfCode(os);
}

/* static */ const Biotype& Biotype::get(TinkerBiotypeId tinkerBiotypeId) {
    assert(tinkerIdExists(tinkerBiotypeId));
    BiotypeId id = biotypeIdsByTinkerId.find(tinkerBiotypeId)->second;
    assert(biotypeIdExists(id));
    return biotypesById.find(id)->second;
}
 
/* static */ const Biotype& Biotype::get(const char* residueName, 
                          const char* atomName, 
                          Ordinality::Residue ordinality)
{
    assert(exists(residueName, atomName, ordinality));

    BiotypeKey key1(residueName, atomName, ordinality);
    BiotypeKey key2(residueName, atomName, Ordinality::ANY);

    BiotypeId id;
    if (biotypeExists(key1))
        id = biotypeIdsByKey.find(key1)->second;
    else 
        id = biotypeIdsByKey.find(key2)->second;

    return get(id);
}

/* static */ const Biotype& Biotype::get(BiotypeId id)
{
    assert(biotypeIdExists(id));

    return biotypesById.find(id)->second;
}

/* static */ Biotype& Biotype::upd(BiotypeId biotypeId) 
{
    assert(biotypeIdExists(biotypeId));
    return biotypesById.find(biotypeId)->second;
}

/* static */ Biotype& Biotype::upd(TinkerBiotypeId tinkerBiotypeId)
{
    assert(tinkerIdExists(tinkerBiotypeId));
    BiotypeId id = biotypeIdsByTinkerId.find(tinkerBiotypeId)->second;
    return upd(id);
}

/* static */ Biotype& Biotype::upd(const char* residueName, 
                          const char* atomName, 
                          Ordinality::Residue ordinality)
{
    assert(exists(residueName, atomName, ordinality));

    BiotypeKey key1(residueName, atomName, ordinality);
    BiotypeKey key2(residueName, atomName, Ordinality::ANY);

    BiotypeId id;
    if (biotypeExists(key1))
        id = biotypeIdsByKey.find(key1)->second;
    else 
        id = biotypeIdsByKey.find(key2)->second;

    assert(biotypeIdExists(id));

    return upd(id);
}


/* static */ std::ostream& Biotype::generateAllBiotypeCode(std::ostream& os) {
    std::map<BiotypeId, Biotype>::const_iterator biotypeI;
    for (biotypeI = biotypesById.begin(); biotypeI != biotypesById.end(); ++biotypeI) 
    {
        biotypeI->second.generateSelfCode(os);
    }

    return os;
}

 //Biotype::Biotype(const char* name,
//                 const Element& elt, int n, int chg, TinkerBiotypeId tid) 
//{
//    tinkerBiotypeIfAny = tid;
//    element            = elt;
//    valence            = n;
//    formalCharge       = chg;
//}

/* static */ const Biotype& Biotype::Argon() 
{
    initializePopularBiotypes();
    return Biotype::get("argon", "argon");
}

/* static */ const Biotype& Biotype::SerineO() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "O") )
        Biotype::defineBiotype(Element::Oxygen(), 1, "serine", "O");
    return Biotype::get("serine", "O");
}

/* static */ const Biotype& Biotype::SerineC() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "C") )
        Biotype::defineBiotype(Element::Carbon(), 3, "serine", "C");
    return Biotype::get("serine", "C");
}

/* static */ const Biotype& Biotype::SerineHN() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HN") )
        Biotype::defineBiotype(Element::Hydrogen(), 1, "serine", "HN");
    return Biotype::get("serine", "HN");
}

/* static */ const Biotype& Biotype::SerineHA() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HA") )
        Biotype::defineBiotype(Element::Hydrogen(), 1, "serine", "HA");
    return Biotype::get("serine", "HA");
}

/* static */ const Biotype& Biotype::SerineHG() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HG") )
        Biotype::defineBiotype(Element::Hydrogen(), 1, "serine", "HG");
    return Biotype::get("serine", "HG");
}

/* static */ const Biotype& Biotype::SerineOG() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "OG") )
        Biotype::defineBiotype(Element::Oxygen(), 2, "serine", "OG");
    return Biotype::get("serine", "OG");
}

/* static */ const Biotype& Biotype::SerineN() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "N") )
        Biotype::defineBiotype(Element::Nitrogen(), 3, "serine", "N");
    return Biotype::get("serine", "N");
}

/* static */ const Biotype& Biotype::SerineCA() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HA") )
        Biotype::defineBiotype(Element::Carbon(), 4, "serine", "CA");
    return Biotype::get("serine", "CA");
}

/* static */ const Biotype& Biotype::MethaneH() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("methane", "H") )
        Biotype::defineBiotype(Element::Hydrogen(), 1, "methane", "H");
    return Biotype::get("methane", "H");
}

/* static */ const Biotype& Biotype::MethaneC() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("methane", "C") )
        Biotype::defineBiotype(Element::Carbon(), 4, "methane", "C");
    return Biotype::get("methane", "C");
}

/* static */ const Biotype& Biotype::EthaneH() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("ethane", "H") )
        Biotype::defineBiotype(Element::Hydrogen(), 1, "ethane", "H");
    return Biotype::get("ethane", "H");
}

/* static */ const Biotype& Biotype::EthaneC() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("ethane", "C") )
        Biotype::defineBiotype(Element::Carbon(), 4, "ethane", "C");
    return Biotype::get("ethane", "C");
}




/* static */ void Biotype::initializePopularBiotypes() 
{
    static bool popularBiotypesAreInitialized = false;
    if (popularBiotypesAreInitialized) return;

    if (! Biotype::exists("argon", "argon", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeId
            , Element::Argon()
            , 0
            , "argon"
            , "argon"
            , Ordinality::ANY
            );

    if (! Biotype::exists("methane", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeId
            , Element::Carbon()
            , 4
            , "methane"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("methane", "H", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeId
            , Element::Hydrogen()
            , 1
            , "methane"
            , "H"
            , Ordinality::ANY
            );

    if (! Biotype::exists("ethane", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeId
            , Element::Carbon()
            , 4
            , "ethane"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("ethane", "H", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeId
            , Element::Hydrogen()
            , 1
            , "ethane"
            , "H"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glycine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1)
            , Element::Nitrogen()
            , 3
            , "Glycine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glycine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(2)
            , Element::Carbon()
            , 4
            , "Glycine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glycine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(3)
            , Element::Carbon()
            , 3
            , "Glycine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glycine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(4)
            , Element::Hydrogen()
            , 1
            , "Glycine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glycine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(5)
            , Element::Oxygen()
            , 1
            , "Glycine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glycine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(6)
            , Element::Hydrogen()
            , 1
            , "Glycine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(7)
            , Element::Nitrogen()
            , 3
            , "Alanine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(8)
            , Element::Carbon()
            , 4
            , "Alanine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(9)
            , Element::Carbon()
            , 3
            , "Alanine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(10)
            , Element::Hydrogen()
            , 1
            , "Alanine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(11)
            , Element::Oxygen()
            , 1
            , "Alanine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(12)
            , Element::Hydrogen()
            , 1
            , "Alanine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(13)
            , Element::Carbon()
            , 4
            , "Alanine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Alanine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(14)
            , Element::Hydrogen()
            , 1
            , "Alanine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(15)
            , Element::Nitrogen()
            , 3
            , "Valine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(16)
            , Element::Carbon()
            , 4
            , "Valine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(17)
            , Element::Carbon()
            , 3
            , "Valine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(18)
            , Element::Hydrogen()
            , 1
            , "Valine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(19)
            , Element::Oxygen()
            , 1
            , "Valine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(20)
            , Element::Hydrogen()
            , 1
            , "Valine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(21)
            , Element::Carbon()
            , 4
            , "Valine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(22)
            , Element::Hydrogen()
            , 1
            , "Valine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "CG1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(23)
            , Element::Carbon()
            , 4
            , "Valine"
            , "CG1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "HG1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(24)
            , Element::Hydrogen()
            , 1
            , "Valine"
            , "HG1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "CG2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(25)
            , Element::Carbon()
            , 4
            , "Valine"
            , "CG2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Valine", "HG2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(26)
            , Element::Hydrogen()
            , 1
            , "Valine"
            , "HG2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(27)
            , Element::Nitrogen()
            , 3
            , "Leucine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(28)
            , Element::Carbon()
            , 4
            , "Leucine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(29)
            , Element::Carbon()
            , 3
            , "Leucine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(30)
            , Element::Hydrogen()
            , 1
            , "Leucine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(31)
            , Element::Oxygen()
            , 1
            , "Leucine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(32)
            , Element::Hydrogen()
            , 1
            , "Leucine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(33)
            , Element::Carbon()
            , 4
            , "Leucine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(34)
            , Element::Hydrogen()
            , 1
            , "Leucine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(35)
            , Element::Carbon()
            , 4
            , "Leucine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(36)
            , Element::Hydrogen()
            , 1
            , "Leucine"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "CD1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(37)
            , Element::Carbon()
            , 4
            , "Leucine"
            , "CD1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "HD1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(38)
            , Element::Hydrogen()
            , 1
            , "Leucine"
            , "HD1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "CD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(39)
            , Element::Carbon()
            , 4
            , "Leucine"
            , "CD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Leucine", "HD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(40)
            , Element::Hydrogen()
            , 1
            , "Leucine"
            , "HD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(41)
            , Element::Nitrogen()
            , 3
            , "Isoleucine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(42)
            , Element::Carbon()
            , 4
            , "Isoleucine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(43)
            , Element::Carbon()
            , 3
            , "Isoleucine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(44)
            , Element::Hydrogen()
            , 1
            , "Isoleucine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(45)
            , Element::Oxygen()
            , 1
            , "Isoleucine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(46)
            , Element::Hydrogen()
            , 1
            , "Isoleucine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(47)
            , Element::Carbon()
            , 4
            , "Isoleucine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(48)
            , Element::Hydrogen()
            , 1
            , "Isoleucine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "CG1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(49)
            , Element::Carbon()
            , 4
            , "Isoleucine"
            , "CG1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "HG1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(50)
            , Element::Hydrogen()
            , 1
            , "Isoleucine"
            , "HG1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "CG2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(51)
            , Element::Carbon()
            , 4
            , "Isoleucine"
            , "CG2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "HG2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(52)
            , Element::Hydrogen()
            , 1
            , "Isoleucine"
            , "HG2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(53)
            , Element::Carbon()
            , 4
            , "Isoleucine"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Isoleucine", "HD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(54)
            , Element::Hydrogen()
            , 1
            , "Isoleucine"
            , "HD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(55)
            , Element::Nitrogen()
            , 3
            , "Serine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(56)
            , Element::Carbon()
            , 4
            , "Serine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(57)
            , Element::Carbon()
            , 3
            , "Serine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(58)
            , Element::Hydrogen()
            , 1
            , "Serine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(59)
            , Element::Oxygen()
            , 1
            , "Serine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(60)
            , Element::Hydrogen()
            , 1
            , "Serine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(61)
            , Element::Carbon()
            , 4
            , "Serine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(62)
            , Element::Hydrogen()
            , 1
            , "Serine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "OG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(63)
            , Element::Oxygen()
            , 2
            , "Serine"
            , "OG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Serine", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(64)
            , Element::Hydrogen()
            , 1
            , "Serine"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(65)
            , Element::Nitrogen()
            , 3
            , "Threonine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(66)
            , Element::Carbon()
            , 4
            , "Threonine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(67)
            , Element::Carbon()
            , 3
            , "Threonine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(68)
            , Element::Hydrogen()
            , 1
            , "Threonine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(69)
            , Element::Oxygen()
            , 1
            , "Threonine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(70)
            , Element::Hydrogen()
            , 1
            , "Threonine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(71)
            , Element::Carbon()
            , 4
            , "Threonine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(72)
            , Element::Hydrogen()
            , 1
            , "Threonine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "OG1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(73)
            , Element::Oxygen()
            , 2
            , "Threonine"
            , "OG1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "HG1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(74)
            , Element::Hydrogen()
            , 1
            , "Threonine"
            , "HG1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "CG2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(75)
            , Element::Carbon()
            , 4
            , "Threonine"
            , "CG2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Threonine", "HG2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(76)
            , Element::Hydrogen()
            , 1
            , "Threonine"
            , "HG2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(77)
            , Element::Nitrogen()
            , 3
            , "Cysteine (-SH)"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(78)
            , Element::Carbon()
            , 4
            , "Cysteine (-SH)"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(79)
            , Element::Carbon()
            , 3
            , "Cysteine (-SH)"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(80)
            , Element::Hydrogen()
            , 1
            , "Cysteine (-SH)"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(81)
            , Element::Oxygen()
            , 1
            , "Cysteine (-SH)"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(82)
            , Element::Hydrogen()
            , 1
            , "Cysteine (-SH)"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(83)
            , Element::Carbon()
            , 4
            , "Cysteine (-SH)"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(84)
            , Element::Hydrogen()
            , 1
            , "Cysteine (-SH)"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "SG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(85)
            , Element::Sulfur()
            , 2
            , "Cysteine (-SH)"
            , "SG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cysteine (-SH)", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(86)
            , Element::Hydrogen()
            , 1
            , "Cysteine (-SH)"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(87)
            , Element::Nitrogen()
            , 3
            , "Cystine (-SS-)"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(88)
            , Element::Carbon()
            , 4
            , "Cystine (-SS-)"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(89)
            , Element::Carbon()
            , 3
            , "Cystine (-SS-)"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(90)
            , Element::Hydrogen()
            , 1
            , "Cystine (-SS-)"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(91)
            , Element::Oxygen()
            , 1
            , "Cystine (-SS-)"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(92)
            , Element::Hydrogen()
            , 1
            , "Cystine (-SS-)"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(93)
            , Element::Carbon()
            , 4
            , "Cystine (-SS-)"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(94)
            , Element::Hydrogen()
            , 1
            , "Cystine (-SS-)"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cystine (-SS-)", "SG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(95)
            , Element::Sulfur()
            , 2
            , "Cystine (-SS-)"
            , "SG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(96)
            , Element::Nitrogen()
            , 3
            , "Proline"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(97)
            , Element::Carbon()
            , 4
            , "Proline"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(98)
            , Element::Carbon()
            , 3
            , "Proline"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(99)
            , Element::Oxygen()
            , 1
            , "Proline"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(100)
            , Element::Hydrogen()
            , 1
            , "Proline"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(101)
            , Element::Carbon()
            , 4
            , "Proline"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(102)
            , Element::Hydrogen()
            , 1
            , "Proline"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(103)
            , Element::Carbon()
            , 4
            , "Proline"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(104)
            , Element::Hydrogen()
            , 1
            , "Proline"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(105)
            , Element::Carbon()
            , 4
            , "Proline"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Proline", "HD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(106)
            , Element::Hydrogen()
            , 1
            , "Proline"
            , "HD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(107)
            , Element::Nitrogen()
            , 3
            , "Phenylalanine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(108)
            , Element::Carbon()
            , 4
            , "Phenylalanine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(109)
            , Element::Carbon()
            , 3
            , "Phenylalanine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(110)
            , Element::Hydrogen()
            , 1
            , "Phenylalanine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(111)
            , Element::Oxygen()
            , 1
            , "Phenylalanine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(112)
            , Element::Hydrogen()
            , 1
            , "Phenylalanine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(113)
            , Element::Carbon()
            , 4
            , "Phenylalanine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(114)
            , Element::Hydrogen()
            , 1
            , "Phenylalanine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(115)
            , Element::Carbon()
            , 3
            , "Phenylalanine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(116)
            , Element::Carbon()
            , 3
            , "Phenylalanine"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "HD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(117)
            , Element::Hydrogen()
            , 1
            , "Phenylalanine"
            , "HD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "CE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(118)
            , Element::Carbon()
            , 3
            , "Phenylalanine"
            , "CE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "HE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(119)
            , Element::Hydrogen()
            , 1
            , "Phenylalanine"
            , "HE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "CZ", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(120)
            , Element::Carbon()
            , 3
            , "Phenylalanine"
            , "CZ"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phenylalanine", "HZ", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(121)
            , Element::Hydrogen()
            , 1
            , "Phenylalanine"
            , "HZ"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(122)
            , Element::Nitrogen()
            , 3
            , "Tyrosine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(123)
            , Element::Carbon()
            , 4
            , "Tyrosine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(124)
            , Element::Carbon()
            , 3
            , "Tyrosine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(125)
            , Element::Hydrogen()
            , 1
            , "Tyrosine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(126)
            , Element::Oxygen()
            , 1
            , "Tyrosine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(127)
            , Element::Hydrogen()
            , 1
            , "Tyrosine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(128)
            , Element::Carbon()
            , 4
            , "Tyrosine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(129)
            , Element::Hydrogen()
            , 1
            , "Tyrosine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(130)
            , Element::Carbon()
            , 3
            , "Tyrosine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(131)
            , Element::Carbon()
            , 3
            , "Tyrosine"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "HD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(132)
            , Element::Hydrogen()
            , 1
            , "Tyrosine"
            , "HD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "CE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(133)
            , Element::Carbon()
            , 3
            , "Tyrosine"
            , "CE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "HE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(134)
            , Element::Hydrogen()
            , 1
            , "Tyrosine"
            , "HE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "CZ", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(135)
            , Element::Carbon()
            , 3
            , "Tyrosine"
            , "CZ"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "OH", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(136)
            , Element::Oxygen()
            , 2
            , "Tyrosine"
            , "OH"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tyrosine", "HH", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(137)
            , Element::Hydrogen()
            , 1
            , "Tyrosine"
            , "HH"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(138)
            , Element::Nitrogen()
            , 3
            , "Tryptophan"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(139)
            , Element::Carbon()
            , 4
            , "Tryptophan"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(140)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(141)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(142)
            , Element::Oxygen()
            , 1
            , "Tryptophan"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(143)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(144)
            , Element::Carbon()
            , 4
            , "Tryptophan"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(145)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(146)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CD1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(147)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CD1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HD1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(148)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HD1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(149)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "NE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(150)
            , Element::Nitrogen()
            , 3
            , "Tryptophan"
            , "NE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(151)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(152)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CE3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(153)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CE3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HE3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(154)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HE3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CZ2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(155)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CZ2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HZ2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(156)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HZ2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CZ3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(157)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CZ3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HZ3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(158)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HZ3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "CH2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(159)
            , Element::Carbon()
            , 3
            , "Tryptophan"
            , "CH2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Tryptophan", "HH2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(160)
            , Element::Hydrogen()
            , 1
            , "Tryptophan"
            , "HH2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(161)
            , Element::Nitrogen()
            , 3
            , "Histidine (+)"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(162)
            , Element::Carbon()
            , 4
            , "Histidine (+)"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(163)
            , Element::Carbon()
            , 3
            , "Histidine (+)"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(164)
            , Element::Hydrogen()
            , 1
            , "Histidine (+)"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(165)
            , Element::Oxygen()
            , 1
            , "Histidine (+)"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(166)
            , Element::Hydrogen()
            , 1
            , "Histidine (+)"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(167)
            , Element::Carbon()
            , 4
            , "Histidine (+)"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(168)
            , Element::Hydrogen()
            , 1
            , "Histidine (+)"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(169)
            , Element::Carbon()
            , 3
            , "Histidine (+)"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "ND1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(170)
            , Element::Nitrogen()
            , 3
            , "Histidine (+)"
            , "ND1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "HD1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(171)
            , Element::Hydrogen()
            , 1
            , "Histidine (+)"
            , "HD1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "CD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(172)
            , Element::Carbon()
            , 3
            , "Histidine (+)"
            , "CD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "HD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(173)
            , Element::Hydrogen()
            , 1
            , "Histidine (+)"
            , "HD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "CE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(174)
            , Element::Carbon()
            , 3
            , "Histidine (+)"
            , "CE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "HE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(175)
            , Element::Hydrogen()
            , 1
            , "Histidine (+)"
            , "HE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "NE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(176)
            , Element::Nitrogen()
            , 3
            , "Histidine (+)"
            , "NE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (+)", "HE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(177)
            , Element::Hydrogen()
            , 1
            , "Histidine (+)"
            , "HE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(178)
            , Element::Nitrogen()
            , 3
            , "Histidine (HD)"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(179)
            , Element::Carbon()
            , 4
            , "Histidine (HD)"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(180)
            , Element::Carbon()
            , 3
            , "Histidine (HD)"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(181)
            , Element::Hydrogen()
            , 1
            , "Histidine (HD)"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(182)
            , Element::Oxygen()
            , 1
            , "Histidine (HD)"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(183)
            , Element::Hydrogen()
            , 1
            , "Histidine (HD)"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(184)
            , Element::Carbon()
            , 4
            , "Histidine (HD)"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(185)
            , Element::Hydrogen()
            , 1
            , "Histidine (HD)"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(186)
            , Element::Carbon()
            , 3
            , "Histidine (HD)"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "ND1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(187)
            , Element::Nitrogen()
            , 3
            , "Histidine (HD)"
            , "ND1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "HD1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(188)
            , Element::Hydrogen()
            , 1
            , "Histidine (HD)"
            , "HD1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "CD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(189)
            , Element::Carbon()
            , 3
            , "Histidine (HD)"
            , "CD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "HD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(190)
            , Element::Hydrogen()
            , 1
            , "Histidine (HD)"
            , "HD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "CE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(191)
            , Element::Carbon()
            , 3
            , "Histidine (HD)"
            , "CE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "HE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(192)
            , Element::Hydrogen()
            , 1
            , "Histidine (HD)"
            , "HE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HD)", "NE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(193)
            , Element::Nitrogen()
            , 2
            , "Histidine (HD)"
            , "NE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(194)
            , Element::Nitrogen()
            , 3
            , "Histidine (HE)"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(195)
            , Element::Carbon()
            , 4
            , "Histidine (HE)"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(196)
            , Element::Carbon()
            , 3
            , "Histidine (HE)"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(197)
            , Element::Hydrogen()
            , 1
            , "Histidine (HE)"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(198)
            , Element::Oxygen()
            , 1
            , "Histidine (HE)"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(199)
            , Element::Hydrogen()
            , 1
            , "Histidine (HE)"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(200)
            , Element::Carbon()
            , 4
            , "Histidine (HE)"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(201)
            , Element::Hydrogen()
            , 1
            , "Histidine (HE)"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(202)
            , Element::Carbon()
            , 3
            , "Histidine (HE)"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "ND1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(203)
            , Element::Nitrogen()
            , 2
            , "Histidine (HE)"
            , "ND1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "CD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(204)
            , Element::Carbon()
            , 3
            , "Histidine (HE)"
            , "CD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "HD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(205)
            , Element::Hydrogen()
            , 1
            , "Histidine (HE)"
            , "HD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "CE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(206)
            , Element::Carbon()
            , 3
            , "Histidine (HE)"
            , "CE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "HE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(207)
            , Element::Hydrogen()
            , 1
            , "Histidine (HE)"
            , "HE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "NE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(208)
            , Element::Nitrogen()
            , 3
            , "Histidine (HE)"
            , "NE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Histidine (HE)", "HE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(209)
            , Element::Hydrogen()
            , 1
            , "Histidine (HE)"
            , "HE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(210)
            , Element::Nitrogen()
            , 3
            , "Aspartic Acid"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(211)
            , Element::Carbon()
            , 4
            , "Aspartic Acid"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(212)
            , Element::Carbon()
            , 3
            , "Aspartic Acid"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(213)
            , Element::Hydrogen()
            , 1
            , "Aspartic Acid"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(214)
            , Element::Oxygen()
            , 1
            , "Aspartic Acid"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(215)
            , Element::Hydrogen()
            , 1
            , "Aspartic Acid"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(216)
            , Element::Carbon()
            , 4
            , "Aspartic Acid"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(217)
            , Element::Hydrogen()
            , 1
            , "Aspartic Acid"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(218)
            , Element::Carbon()
            , 3
            , "Aspartic Acid"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Aspartic Acid", "OD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(219)
            , Element::Oxygen()
            , 1
            , "Aspartic Acid"
            , "OD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(220)
            , Element::Nitrogen()
            , 3
            , "Asparagine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(221)
            , Element::Carbon()
            , 4
            , "Asparagine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(222)
            , Element::Carbon()
            , 3
            , "Asparagine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(223)
            , Element::Hydrogen()
            , 1
            , "Asparagine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(224)
            , Element::Oxygen()
            , 1
            , "Asparagine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(225)
            , Element::Hydrogen()
            , 1
            , "Asparagine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(226)
            , Element::Carbon()
            , 4
            , "Asparagine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(227)
            , Element::Hydrogen()
            , 1
            , "Asparagine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(228)
            , Element::Carbon()
            , 3
            , "Asparagine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "OD1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(229)
            , Element::Oxygen()
            , 1
            , "Asparagine"
            , "OD1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "ND2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(230)
            , Element::Nitrogen()
            , 3
            , "Asparagine"
            , "ND2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Asparagine", "HD2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(231)
            , Element::Hydrogen()
            , 1
            , "Asparagine"
            , "HD2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(232)
            , Element::Nitrogen()
            , 3
            , "Glutamic Acid"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(233)
            , Element::Carbon()
            , 4
            , "Glutamic Acid"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(234)
            , Element::Carbon()
            , 3
            , "Glutamic Acid"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(235)
            , Element::Hydrogen()
            , 1
            , "Glutamic Acid"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(236)
            , Element::Oxygen()
            , 1
            , "Glutamic Acid"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(237)
            , Element::Hydrogen()
            , 1
            , "Glutamic Acid"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(238)
            , Element::Carbon()
            , 4
            , "Glutamic Acid"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(239)
            , Element::Hydrogen()
            , 1
            , "Glutamic Acid"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(240)
            , Element::Carbon()
            , 4
            , "Glutamic Acid"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(241)
            , Element::Hydrogen()
            , 1
            , "Glutamic Acid"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(242)
            , Element::Carbon()
            , 3
            , "Glutamic Acid"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamic Acid", "OE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(243)
            , Element::Oxygen()
            , 1
            , "Glutamic Acid"
            , "OE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(244)
            , Element::Nitrogen()
            , 3
            , "Glutamine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(245)
            , Element::Carbon()
            , 4
            , "Glutamine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(246)
            , Element::Carbon()
            , 3
            , "Glutamine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(247)
            , Element::Hydrogen()
            , 1
            , "Glutamine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(248)
            , Element::Oxygen()
            , 1
            , "Glutamine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(249)
            , Element::Hydrogen()
            , 1
            , "Glutamine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(250)
            , Element::Carbon()
            , 4
            , "Glutamine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(251)
            , Element::Hydrogen()
            , 1
            , "Glutamine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(252)
            , Element::Carbon()
            , 4
            , "Glutamine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(253)
            , Element::Hydrogen()
            , 1
            , "Glutamine"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(254)
            , Element::Carbon()
            , 3
            , "Glutamine"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "OE1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(255)
            , Element::Oxygen()
            , 1
            , "Glutamine"
            , "OE1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "NE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(256)
            , Element::Nitrogen()
            , 3
            , "Glutamine"
            , "NE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Glutamine", "HE2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(257)
            , Element::Hydrogen()
            , 1
            , "Glutamine"
            , "HE2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(258)
            , Element::Nitrogen()
            , 3
            , "Methionine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(259)
            , Element::Carbon()
            , 4
            , "Methionine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(260)
            , Element::Carbon()
            , 3
            , "Methionine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(261)
            , Element::Hydrogen()
            , 1
            , "Methionine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(262)
            , Element::Oxygen()
            , 1
            , "Methionine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(263)
            , Element::Hydrogen()
            , 1
            , "Methionine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(264)
            , Element::Carbon()
            , 4
            , "Methionine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(265)
            , Element::Hydrogen()
            , 1
            , "Methionine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(266)
            , Element::Carbon()
            , 4
            , "Methionine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(267)
            , Element::Hydrogen()
            , 1
            , "Methionine"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "SD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(268)
            , Element::Sulfur()
            , 2
            , "Methionine"
            , "SD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "CE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(269)
            , Element::Carbon()
            , 4
            , "Methionine"
            , "CE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Methionine", "HE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(270)
            , Element::Hydrogen()
            , 1
            , "Methionine"
            , "HE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(271)
            , Element::Nitrogen()
            , 3
            , "Lysine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(272)
            , Element::Carbon()
            , 4
            , "Lysine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(273)
            , Element::Carbon()
            , 3
            , "Lysine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(274)
            , Element::Hydrogen()
            , 1
            , "Lysine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(275)
            , Element::Oxygen()
            , 1
            , "Lysine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(276)
            , Element::Hydrogen()
            , 1
            , "Lysine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(277)
            , Element::Carbon()
            , 4
            , "Lysine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(278)
            , Element::Hydrogen()
            , 1
            , "Lysine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(279)
            , Element::Carbon()
            , 4
            , "Lysine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(280)
            , Element::Hydrogen()
            , 1
            , "Lysine"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(281)
            , Element::Carbon()
            , 4
            , "Lysine"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "HD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(282)
            , Element::Hydrogen()
            , 1
            , "Lysine"
            , "HD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "CE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(283)
            , Element::Carbon()
            , 4
            , "Lysine"
            , "CE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "HE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(284)
            , Element::Hydrogen()
            , 1
            , "Lysine"
            , "HE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "NZ", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(285)
            , Element::Nitrogen()
            , 4
            , "Lysine"
            , "NZ"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Lysine", "HZ", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(286)
            , Element::Hydrogen()
            , 1
            , "Lysine"
            , "HZ"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(287)
            , Element::Nitrogen()
            , 3
            , "Arginine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(288)
            , Element::Carbon()
            , 4
            , "Arginine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(289)
            , Element::Carbon()
            , 3
            , "Arginine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(290)
            , Element::Hydrogen()
            , 1
            , "Arginine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(291)
            , Element::Oxygen()
            , 1
            , "Arginine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(292)
            , Element::Hydrogen()
            , 1
            , "Arginine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(293)
            , Element::Carbon()
            , 4
            , "Arginine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(294)
            , Element::Hydrogen()
            , 1
            , "Arginine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(295)
            , Element::Carbon()
            , 4
            , "Arginine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(296)
            , Element::Hydrogen()
            , 1
            , "Arginine"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(297)
            , Element::Carbon()
            , 4
            , "Arginine"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "HD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(298)
            , Element::Hydrogen()
            , 1
            , "Arginine"
            , "HD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "NE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(299)
            , Element::Nitrogen()
            , 3
            , "Arginine"
            , "NE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "HE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(300)
            , Element::Hydrogen()
            , 1
            , "Arginine"
            , "HE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "CZ", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(301)
            , Element::Carbon()
            , 3
            , "Arginine"
            , "CZ"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "NH", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(302)
            , Element::Nitrogen()
            , 3
            , "Arginine"
            , "NH"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Arginine", "HH", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(303)
            , Element::Hydrogen()
            , 1
            , "Arginine"
            , "HH"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(304)
            , Element::Nitrogen()
            , 3
            , "Ornithine"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(305)
            , Element::Carbon()
            , 4
            , "Ornithine"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(306)
            , Element::Carbon()
            , 3
            , "Ornithine"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(307)
            , Element::Hydrogen()
            , 1
            , "Ornithine"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(308)
            , Element::Oxygen()
            , 1
            , "Ornithine"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(309)
            , Element::Hydrogen()
            , 1
            , "Ornithine"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(310)
            , Element::Carbon()
            , 4
            , "Ornithine"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(311)
            , Element::Hydrogen()
            , 1
            , "Ornithine"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(312)
            , Element::Carbon()
            , 4
            , "Ornithine"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(313)
            , Element::Hydrogen()
            , 1
            , "Ornithine"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(314)
            , Element::Carbon()
            , 4
            , "Ornithine"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "HD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(315)
            , Element::Hydrogen()
            , 1
            , "Ornithine"
            , "HD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "NE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(316)
            , Element::Nitrogen()
            , 4
            , "Ornithine"
            , "NE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Ornithine", "HE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(317)
            , Element::Hydrogen()
            , 1
            , "Ornithine"
            , "HE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(318)
            , Element::Nitrogen()
            , 3
            , "MethylAlanine (AIB)"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(319)
            , Element::Carbon()
            , 4
            , "MethylAlanine (AIB)"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(320)
            , Element::Carbon()
            , 3
            , "MethylAlanine (AIB)"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(321)
            , Element::Hydrogen()
            , 1
            , "MethylAlanine (AIB)"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(322)
            , Element::Oxygen()
            , 1
            , "MethylAlanine (AIB)"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(323)
            , Element::Carbon()
            , 4
            , "MethylAlanine (AIB)"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(324)
            , Element::Hydrogen()
            , 1
            , "MethylAlanine (AIB)"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "N", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(325)
            , Element::Nitrogen()
            , 3
            , "Pyroglutamic Acid"
            , "N"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(326)
            , Element::Carbon()
            , 4
            , "Pyroglutamic Acid"
            , "CA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "C", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(327)
            , Element::Carbon()
            , 3
            , "Pyroglutamic Acid"
            , "C"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HN", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(328)
            , Element::Hydrogen()
            , 1
            , "Pyroglutamic Acid"
            , "HN"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "O", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(329)
            , Element::Oxygen()
            , 1
            , "Pyroglutamic Acid"
            , "O"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HA", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(330)
            , Element::Hydrogen()
            , 1
            , "Pyroglutamic Acid"
            , "HA"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(331)
            , Element::Carbon()
            , 4
            , "Pyroglutamic Acid"
            , "CB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HB", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(332)
            , Element::Hydrogen()
            , 1
            , "Pyroglutamic Acid"
            , "HB"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(333)
            , Element::Carbon()
            , 4
            , "Pyroglutamic Acid"
            , "CG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HG", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(334)
            , Element::Hydrogen()
            , 1
            , "Pyroglutamic Acid"
            , "HG"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CD", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(335)
            , Element::Carbon()
            , 3
            , "Pyroglutamic Acid"
            , "CD"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Pyroglutamic Acid", "OE", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(336)
            , Element::Oxygen()
            , 1
            , "Pyroglutamic Acid"
            , "OE"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Acetyl", "CH3", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(340)
            , Element::Carbon()
            , 4
            , "Acetyl"
            , "CH3"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("Acetyl", "H", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(341)
            , Element::Hydrogen()
            , 1
            , "Acetyl"
            , "H"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("Acetyl", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(342)
            , Element::Carbon()
            , 3
            , "Acetyl"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("Acetyl", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(343)
            , Element::Oxygen()
            , 1
            , "Acetyl"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("Amide", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(344)
            , Element::Nitrogen()
            , 3
            , "Amide"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("Amide", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(345)
            , Element::Hydrogen()
            , 1
            , "Amide"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("N-MeAmide", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(346)
            , Element::Nitrogen()
            , 3
            , "N-MeAmide"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("N-MeAmide", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(347)
            , Element::Hydrogen()
            , 1
            , "N-MeAmide"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("N-MeAmide", "CH3", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(348)
            , Element::Carbon()
            , 4
            , "N-MeAmide"
            , "CH3"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("N-MeAmide", "H", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(349)
            , Element::Hydrogen()
            , 1
            , "N-MeAmide"
            , "H"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLY", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(350)
            , Element::Nitrogen()
            , 4
            , "GLY"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLY", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(351)
            , Element::Carbon()
            , 4
            , "GLY"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLY", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(352)
            , Element::Carbon()
            , 3
            , "GLY"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLY", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(353)
            , Element::Hydrogen()
            , 1
            , "GLY"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLY", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(354)
            , Element::Oxygen()
            , 1
            , "GLY"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLY", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(355)
            , Element::Hydrogen()
            , 1
            , "GLY"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ALA", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(356)
            , Element::Nitrogen()
            , 4
            , "ALA"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ALA", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(357)
            , Element::Carbon()
            , 4
            , "ALA"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ALA", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(358)
            , Element::Carbon()
            , 3
            , "ALA"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ALA", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(359)
            , Element::Hydrogen()
            , 1
            , "ALA"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ALA", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(360)
            , Element::Oxygen()
            , 1
            , "ALA"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ALA", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(361)
            , Element::Hydrogen()
            , 1
            , "ALA"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("VAL", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(362)
            , Element::Nitrogen()
            , 4
            , "VAL"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("VAL", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(363)
            , Element::Carbon()
            , 4
            , "VAL"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("VAL", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(364)
            , Element::Carbon()
            , 3
            , "VAL"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("VAL", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(365)
            , Element::Hydrogen()
            , 1
            , "VAL"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("VAL", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(366)
            , Element::Oxygen()
            , 1
            , "VAL"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("VAL", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(367)
            , Element::Hydrogen()
            , 1
            , "VAL"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LEU", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(368)
            , Element::Nitrogen()
            , 4
            , "LEU"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LEU", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(369)
            , Element::Carbon()
            , 4
            , "LEU"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LEU", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(370)
            , Element::Carbon()
            , 3
            , "LEU"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LEU", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(371)
            , Element::Hydrogen()
            , 1
            , "LEU"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LEU", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(372)
            , Element::Oxygen()
            , 1
            , "LEU"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LEU", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(373)
            , Element::Hydrogen()
            , 1
            , "LEU"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ILE", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(374)
            , Element::Nitrogen()
            , 4
            , "ILE"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ILE", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(375)
            , Element::Carbon()
            , 4
            , "ILE"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ILE", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(376)
            , Element::Carbon()
            , 3
            , "ILE"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ILE", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(377)
            , Element::Hydrogen()
            , 1
            , "ILE"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ILE", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(378)
            , Element::Oxygen()
            , 1
            , "ILE"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ILE", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(379)
            , Element::Hydrogen()
            , 1
            , "ILE"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("SER", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(380)
            , Element::Nitrogen()
            , 4
            , "SER"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("SER", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(381)
            , Element::Carbon()
            , 4
            , "SER"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("SER", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(382)
            , Element::Carbon()
            , 3
            , "SER"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("SER", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(383)
            , Element::Hydrogen()
            , 1
            , "SER"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("SER", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(384)
            , Element::Oxygen()
            , 1
            , "SER"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("SER", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(385)
            , Element::Hydrogen()
            , 1
            , "SER"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("THR", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(386)
            , Element::Nitrogen()
            , 4
            , "THR"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("THR", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(387)
            , Element::Carbon()
            , 4
            , "THR"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("THR", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(388)
            , Element::Carbon()
            , 3
            , "THR"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("THR", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(389)
            , Element::Hydrogen()
            , 1
            , "THR"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("THR", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(390)
            , Element::Oxygen()
            , 1
            , "THR"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("THR", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(391)
            , Element::Hydrogen()
            , 1
            , "THR"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SH)", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(392)
            , Element::Nitrogen()
            , 4
            , "CYS (-SH)"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SH)", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(393)
            , Element::Carbon()
            , 4
            , "CYS (-SH)"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SH)", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(394)
            , Element::Carbon()
            , 3
            , "CYS (-SH)"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SH)", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(395)
            , Element::Hydrogen()
            , 1
            , "CYS (-SH)"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SH)", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(396)
            , Element::Oxygen()
            , 1
            , "CYS (-SH)"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SH)", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(397)
            , Element::Hydrogen()
            , 1
            , "CYS (-SH)"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SS)", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(398)
            , Element::Nitrogen()
            , 4
            , "CYS (-SS)"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SS)", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(399)
            , Element::Carbon()
            , 4
            , "CYS (-SS)"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SS)", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(400)
            , Element::Carbon()
            , 3
            , "CYS (-SS)"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SS)", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(401)
            , Element::Hydrogen()
            , 1
            , "CYS (-SS)"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SS)", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(402)
            , Element::Oxygen()
            , 1
            , "CYS (-SS)"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("CYS (-SS)", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(403)
            , Element::Hydrogen()
            , 1
            , "CYS (-SS)"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(404)
            , Element::Nitrogen()
            , 4
            , "PRO"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(405)
            , Element::Carbon()
            , 4
            , "PRO"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(406)
            , Element::Carbon()
            , 3
            , "PRO"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(407)
            , Element::Hydrogen()
            , 1
            , "PRO"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(408)
            , Element::Oxygen()
            , 1
            , "PRO"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(409)
            , Element::Hydrogen()
            , 1
            , "PRO"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "CD", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(410)
            , Element::Carbon()
            , 4
            , "PRO"
            , "CD"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PRO", "HD", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(411)
            , Element::Hydrogen()
            , 1
            , "PRO"
            , "HD"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PHE", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(412)
            , Element::Nitrogen()
            , 4
            , "PHE"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PHE", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(413)
            , Element::Carbon()
            , 4
            , "PHE"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PHE", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(414)
            , Element::Carbon()
            , 3
            , "PHE"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PHE", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(415)
            , Element::Hydrogen()
            , 1
            , "PHE"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PHE", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(416)
            , Element::Oxygen()
            , 1
            , "PHE"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("PHE", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(417)
            , Element::Hydrogen()
            , 1
            , "PHE"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TYR", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(418)
            , Element::Nitrogen()
            , 4
            , "TYR"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TYR", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(419)
            , Element::Carbon()
            , 4
            , "TYR"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TYR", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(420)
            , Element::Carbon()
            , 3
            , "TYR"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TYR", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(421)
            , Element::Hydrogen()
            , 1
            , "TYR"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TYR", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(422)
            , Element::Oxygen()
            , 1
            , "TYR"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TYR", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(423)
            , Element::Hydrogen()
            , 1
            , "TYR"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TRP", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(424)
            , Element::Nitrogen()
            , 4
            , "TRP"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TRP", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(425)
            , Element::Carbon()
            , 4
            , "TRP"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TRP", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(426)
            , Element::Carbon()
            , 3
            , "TRP"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TRP", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(427)
            , Element::Hydrogen()
            , 1
            , "TRP"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TRP", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(428)
            , Element::Oxygen()
            , 1
            , "TRP"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("TRP", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(429)
            , Element::Hydrogen()
            , 1
            , "TRP"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (+)", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(430)
            , Element::Nitrogen()
            , 4
            , "HIS (+)"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (+)", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(431)
            , Element::Carbon()
            , 4
            , "HIS (+)"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (+)", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(432)
            , Element::Carbon()
            , 3
            , "HIS (+)"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (+)", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(433)
            , Element::Hydrogen()
            , 1
            , "HIS (+)"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (+)", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(434)
            , Element::Oxygen()
            , 1
            , "HIS (+)"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (+)", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(435)
            , Element::Hydrogen()
            , 1
            , "HIS (+)"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HD)", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(436)
            , Element::Nitrogen()
            , 4
            , "HIS (HD)"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HD)", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(437)
            , Element::Carbon()
            , 4
            , "HIS (HD)"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HD)", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(438)
            , Element::Carbon()
            , 3
            , "HIS (HD)"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HD)", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(439)
            , Element::Hydrogen()
            , 1
            , "HIS (HD)"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HD)", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(440)
            , Element::Oxygen()
            , 1
            , "HIS (HD)"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HD)", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(441)
            , Element::Hydrogen()
            , 1
            , "HIS (HD)"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HE)", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(442)
            , Element::Nitrogen()
            , 4
            , "HIS (HE)"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HE)", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(443)
            , Element::Carbon()
            , 4
            , "HIS (HE)"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HE)", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(444)
            , Element::Carbon()
            , 3
            , "HIS (HE)"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HE)", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(445)
            , Element::Hydrogen()
            , 1
            , "HIS (HE)"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HE)", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(446)
            , Element::Oxygen()
            , 1
            , "HIS (HE)"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("HIS (HE)", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(447)
            , Element::Hydrogen()
            , 1
            , "HIS (HE)"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASP", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(448)
            , Element::Nitrogen()
            , 4
            , "ASP"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASP", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(449)
            , Element::Carbon()
            , 4
            , "ASP"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASP", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(450)
            , Element::Carbon()
            , 3
            , "ASP"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASP", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(451)
            , Element::Hydrogen()
            , 1
            , "ASP"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASP", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(452)
            , Element::Oxygen()
            , 1
            , "ASP"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASP", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(453)
            , Element::Hydrogen()
            , 1
            , "ASP"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASN", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(454)
            , Element::Nitrogen()
            , 4
            , "ASN"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASN", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(455)
            , Element::Carbon()
            , 4
            , "ASN"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASN", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(456)
            , Element::Carbon()
            , 3
            , "ASN"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASN", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(457)
            , Element::Hydrogen()
            , 1
            , "ASN"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASN", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(458)
            , Element::Oxygen()
            , 1
            , "ASN"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ASN", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(459)
            , Element::Hydrogen()
            , 1
            , "ASN"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLU", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(460)
            , Element::Nitrogen()
            , 4
            , "GLU"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLU", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(461)
            , Element::Carbon()
            , 4
            , "GLU"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLU", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(462)
            , Element::Carbon()
            , 3
            , "GLU"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLU", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(463)
            , Element::Hydrogen()
            , 1
            , "GLU"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLU", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(464)
            , Element::Oxygen()
            , 1
            , "GLU"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLU", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(465)
            , Element::Hydrogen()
            , 1
            , "GLU"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLN", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(466)
            , Element::Nitrogen()
            , 4
            , "GLN"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLN", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(467)
            , Element::Carbon()
            , 4
            , "GLN"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLN", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(468)
            , Element::Carbon()
            , 3
            , "GLN"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLN", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(469)
            , Element::Hydrogen()
            , 1
            , "GLN"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLN", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(470)
            , Element::Oxygen()
            , 1
            , "GLN"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLN", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(471)
            , Element::Hydrogen()
            , 1
            , "GLN"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("MET", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(472)
            , Element::Nitrogen()
            , 4
            , "MET"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("MET", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(473)
            , Element::Carbon()
            , 4
            , "MET"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("MET", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(474)
            , Element::Carbon()
            , 3
            , "MET"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("MET", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(475)
            , Element::Hydrogen()
            , 1
            , "MET"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("MET", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(476)
            , Element::Oxygen()
            , 1
            , "MET"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("MET", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(477)
            , Element::Hydrogen()
            , 1
            , "MET"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LYS", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(478)
            , Element::Nitrogen()
            , 4
            , "LYS"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LYS", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(479)
            , Element::Carbon()
            , 4
            , "LYS"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LYS", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(480)
            , Element::Carbon()
            , 3
            , "LYS"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LYS", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(481)
            , Element::Hydrogen()
            , 1
            , "LYS"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LYS", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(482)
            , Element::Oxygen()
            , 1
            , "LYS"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("LYS", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(483)
            , Element::Hydrogen()
            , 1
            , "LYS"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ARG", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(484)
            , Element::Nitrogen()
            , 4
            , "ARG"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ARG", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(485)
            , Element::Carbon()
            , 4
            , "ARG"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ARG", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(486)
            , Element::Carbon()
            , 3
            , "ARG"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ARG", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(487)
            , Element::Hydrogen()
            , 1
            , "ARG"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ARG", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(488)
            , Element::Oxygen()
            , 1
            , "ARG"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ARG", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(489)
            , Element::Hydrogen()
            , 1
            , "ARG"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ORN", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(490)
            , Element::Nitrogen()
            , 4
            , "ORN"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ORN", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(491)
            , Element::Carbon()
            , 4
            , "ORN"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ORN", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(492)
            , Element::Carbon()
            , 3
            , "ORN"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ORN", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(493)
            , Element::Hydrogen()
            , 1
            , "ORN"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ORN", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(494)
            , Element::Oxygen()
            , 1
            , "ORN"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("ORN", "HA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(495)
            , Element::Hydrogen()
            , 1
            , "ORN"
            , "HA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("AIB", "N", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(496)
            , Element::Nitrogen()
            , 4
            , "AIB"
            , "N"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("AIB", "CA", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(497)
            , Element::Carbon()
            , 4
            , "AIB"
            , "CA"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("AIB", "C", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(498)
            , Element::Carbon()
            , 3
            , "AIB"
            , "C"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("AIB", "HN", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(499)
            , Element::Hydrogen()
            , 1
            , "AIB"
            , "HN"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("AIB", "O", Ordinality::INITIAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(500)
            , Element::Oxygen()
            , 1
            , "AIB"
            , "O"
            , Ordinality::INITIAL
            );

    if (! Biotype::exists("GLY", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(501)
            , Element::Nitrogen()
            , 3
            , "GLY"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLY", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(502)
            , Element::Carbon()
            , 4
            , "GLY"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLY", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(503)
            , Element::Carbon()
            , 3
            , "GLY"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLY", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(504)
            , Element::Hydrogen()
            , 1
            , "GLY"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLY", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(505)
            , Element::Oxygen()
            , 1
            , "GLY"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLY", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(506)
            , Element::Hydrogen()
            , 1
            , "GLY"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ALA", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(507)
            , Element::Nitrogen()
            , 3
            , "ALA"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ALA", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(508)
            , Element::Carbon()
            , 4
            , "ALA"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ALA", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(509)
            , Element::Carbon()
            , 3
            , "ALA"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ALA", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(510)
            , Element::Hydrogen()
            , 1
            , "ALA"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ALA", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(511)
            , Element::Oxygen()
            , 1
            , "ALA"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ALA", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(512)
            , Element::Hydrogen()
            , 1
            , "ALA"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("VAL", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(513)
            , Element::Nitrogen()
            , 3
            , "VAL"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("VAL", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(514)
            , Element::Carbon()
            , 4
            , "VAL"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("VAL", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(515)
            , Element::Carbon()
            , 3
            , "VAL"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("VAL", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(516)
            , Element::Hydrogen()
            , 1
            , "VAL"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("VAL", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(517)
            , Element::Oxygen()
            , 1
            , "VAL"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("VAL", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(518)
            , Element::Hydrogen()
            , 1
            , "VAL"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LEU", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(519)
            , Element::Nitrogen()
            , 3
            , "LEU"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LEU", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(520)
            , Element::Carbon()
            , 4
            , "LEU"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LEU", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(521)
            , Element::Carbon()
            , 3
            , "LEU"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LEU", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(522)
            , Element::Hydrogen()
            , 1
            , "LEU"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LEU", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(523)
            , Element::Oxygen()
            , 1
            , "LEU"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LEU", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(524)
            , Element::Hydrogen()
            , 1
            , "LEU"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ILE", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(525)
            , Element::Nitrogen()
            , 3
            , "ILE"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ILE", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(526)
            , Element::Carbon()
            , 4
            , "ILE"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ILE", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(527)
            , Element::Carbon()
            , 3
            , "ILE"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ILE", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(528)
            , Element::Hydrogen()
            , 1
            , "ILE"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ILE", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(529)
            , Element::Oxygen()
            , 1
            , "ILE"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ILE", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(530)
            , Element::Hydrogen()
            , 1
            , "ILE"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("SER", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(531)
            , Element::Nitrogen()
            , 3
            , "SER"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("SER", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(532)
            , Element::Carbon()
            , 4
            , "SER"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("SER", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(533)
            , Element::Carbon()
            , 3
            , "SER"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("SER", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(534)
            , Element::Hydrogen()
            , 1
            , "SER"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("SER", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(535)
            , Element::Oxygen()
            , 1
            , "SER"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("SER", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(536)
            , Element::Hydrogen()
            , 1
            , "SER"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("THR", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(537)
            , Element::Nitrogen()
            , 3
            , "THR"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("THR", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(538)
            , Element::Carbon()
            , 4
            , "THR"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("THR", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(539)
            , Element::Carbon()
            , 3
            , "THR"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("THR", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(540)
            , Element::Hydrogen()
            , 1
            , "THR"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("THR", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(541)
            , Element::Oxygen()
            , 1
            , "THR"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("THR", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(542)
            , Element::Hydrogen()
            , 1
            , "THR"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SH)", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(543)
            , Element::Nitrogen()
            , 3
            , "CYS (-SH)"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SH)", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(544)
            , Element::Carbon()
            , 4
            , "CYS (-SH)"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SH)", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(545)
            , Element::Carbon()
            , 3
            , "CYS (-SH)"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SH)", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(546)
            , Element::Hydrogen()
            , 1
            , "CYS (-SH)"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SH)", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(547)
            , Element::Oxygen()
            , 1
            , "CYS (-SH)"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SH)", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(548)
            , Element::Hydrogen()
            , 1
            , "CYS (-SH)"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SS)", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(549)
            , Element::Nitrogen()
            , 3
            , "CYS (-SS)"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SS)", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(550)
            , Element::Carbon()
            , 4
            , "CYS (-SS)"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SS)", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(551)
            , Element::Carbon()
            , 3
            , "CYS (-SS)"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SS)", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(552)
            , Element::Hydrogen()
            , 1
            , "CYS (-SS)"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SS)", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(553)
            , Element::Oxygen()
            , 1
            , "CYS (-SS)"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("CYS (-SS)", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(554)
            , Element::Hydrogen()
            , 1
            , "CYS (-SS)"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PRO", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(555)
            , Element::Nitrogen()
            , 3
            , "PRO"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PRO", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(556)
            , Element::Carbon()
            , 4
            , "PRO"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PRO", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(557)
            , Element::Carbon()
            , 3
            , "PRO"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PRO", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(558)
            , Element::Oxygen()
            , 1
            , "PRO"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PRO", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(559)
            , Element::Hydrogen()
            , 1
            , "PRO"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PHE", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(560)
            , Element::Nitrogen()
            , 3
            , "PHE"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PHE", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(561)
            , Element::Carbon()
            , 4
            , "PHE"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PHE", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(562)
            , Element::Carbon()
            , 3
            , "PHE"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PHE", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(563)
            , Element::Hydrogen()
            , 1
            , "PHE"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PHE", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(564)
            , Element::Oxygen()
            , 1
            , "PHE"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("PHE", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(565)
            , Element::Hydrogen()
            , 1
            , "PHE"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TYR", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(566)
            , Element::Nitrogen()
            , 3
            , "TYR"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TYR", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(567)
            , Element::Carbon()
            , 4
            , "TYR"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TYR", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(568)
            , Element::Carbon()
            , 3
            , "TYR"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TYR", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(569)
            , Element::Hydrogen()
            , 1
            , "TYR"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TYR", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(570)
            , Element::Oxygen()
            , 1
            , "TYR"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TYR", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(571)
            , Element::Hydrogen()
            , 1
            , "TYR"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TRP", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(572)
            , Element::Nitrogen()
            , 3
            , "TRP"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TRP", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(573)
            , Element::Carbon()
            , 4
            , "TRP"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TRP", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(574)
            , Element::Carbon()
            , 3
            , "TRP"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TRP", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(575)
            , Element::Hydrogen()
            , 1
            , "TRP"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TRP", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(576)
            , Element::Oxygen()
            , 1
            , "TRP"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("TRP", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(577)
            , Element::Hydrogen()
            , 1
            , "TRP"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (+)", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(578)
            , Element::Nitrogen()
            , 3
            , "HIS (+)"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (+)", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(579)
            , Element::Carbon()
            , 4
            , "HIS (+)"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (+)", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(580)
            , Element::Carbon()
            , 3
            , "HIS (+)"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (+)", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(581)
            , Element::Hydrogen()
            , 1
            , "HIS (+)"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (+)", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(582)
            , Element::Oxygen()
            , 1
            , "HIS (+)"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (+)", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(583)
            , Element::Hydrogen()
            , 1
            , "HIS (+)"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HD)", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(584)
            , Element::Nitrogen()
            , 3
            , "HIS (HD)"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HD)", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(585)
            , Element::Carbon()
            , 4
            , "HIS (HD)"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HD)", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(586)
            , Element::Carbon()
            , 3
            , "HIS (HD)"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HD)", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(587)
            , Element::Hydrogen()
            , 1
            , "HIS (HD)"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HD)", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(588)
            , Element::Oxygen()
            , 1
            , "HIS (HD)"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HD)", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(589)
            , Element::Hydrogen()
            , 1
            , "HIS (HD)"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HE)", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(590)
            , Element::Nitrogen()
            , 3
            , "HIS (HE)"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HE)", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(591)
            , Element::Carbon()
            , 4
            , "HIS (HE)"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HE)", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(592)
            , Element::Carbon()
            , 3
            , "HIS (HE)"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HE)", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(593)
            , Element::Hydrogen()
            , 1
            , "HIS (HE)"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HE)", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(594)
            , Element::Oxygen()
            , 1
            , "HIS (HE)"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("HIS (HE)", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(595)
            , Element::Hydrogen()
            , 1
            , "HIS (HE)"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASP", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(596)
            , Element::Nitrogen()
            , 3
            , "ASP"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASP", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(597)
            , Element::Carbon()
            , 4
            , "ASP"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASP", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(598)
            , Element::Carbon()
            , 3
            , "ASP"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASP", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(599)
            , Element::Hydrogen()
            , 1
            , "ASP"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASP", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(600)
            , Element::Oxygen()
            , 1
            , "ASP"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASP", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(601)
            , Element::Hydrogen()
            , 1
            , "ASP"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASN", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(602)
            , Element::Nitrogen()
            , 3
            , "ASN"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASN", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(603)
            , Element::Carbon()
            , 4
            , "ASN"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASN", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(604)
            , Element::Carbon()
            , 3
            , "ASN"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASN", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(605)
            , Element::Hydrogen()
            , 1
            , "ASN"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASN", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(606)
            , Element::Oxygen()
            , 1
            , "ASN"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ASN", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(607)
            , Element::Hydrogen()
            , 1
            , "ASN"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLU", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(608)
            , Element::Nitrogen()
            , 3
            , "GLU"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLU", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(609)
            , Element::Carbon()
            , 4
            , "GLU"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLU", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(610)
            , Element::Carbon()
            , 3
            , "GLU"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLU", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(611)
            , Element::Hydrogen()
            , 1
            , "GLU"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLU", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(612)
            , Element::Oxygen()
            , 1
            , "GLU"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLU", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(613)
            , Element::Hydrogen()
            , 1
            , "GLU"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLN", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(614)
            , Element::Nitrogen()
            , 3
            , "GLN"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLN", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(615)
            , Element::Carbon()
            , 4
            , "GLN"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLN", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(616)
            , Element::Carbon()
            , 3
            , "GLN"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLN", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(617)
            , Element::Hydrogen()
            , 1
            , "GLN"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLN", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(618)
            , Element::Oxygen()
            , 1
            , "GLN"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("GLN", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(619)
            , Element::Hydrogen()
            , 1
            , "GLN"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("MET", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(620)
            , Element::Nitrogen()
            , 3
            , "MET"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("MET", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(621)
            , Element::Carbon()
            , 4
            , "MET"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("MET", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(622)
            , Element::Carbon()
            , 3
            , "MET"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("MET", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(623)
            , Element::Hydrogen()
            , 1
            , "MET"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("MET", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(624)
            , Element::Oxygen()
            , 1
            , "MET"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("MET", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(625)
            , Element::Hydrogen()
            , 1
            , "MET"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LYS", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(626)
            , Element::Nitrogen()
            , 3
            , "LYS"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LYS", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(627)
            , Element::Carbon()
            , 4
            , "LYS"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LYS", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(628)
            , Element::Carbon()
            , 3
            , "LYS"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LYS", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(629)
            , Element::Hydrogen()
            , 1
            , "LYS"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LYS", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(630)
            , Element::Oxygen()
            , 1
            , "LYS"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("LYS", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(631)
            , Element::Hydrogen()
            , 1
            , "LYS"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ARG", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(632)
            , Element::Nitrogen()
            , 3
            , "ARG"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ARG", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(633)
            , Element::Carbon()
            , 4
            , "ARG"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ARG", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(634)
            , Element::Carbon()
            , 3
            , "ARG"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ARG", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(635)
            , Element::Hydrogen()
            , 1
            , "ARG"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ARG", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(636)
            , Element::Oxygen()
            , 1
            , "ARG"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ARG", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(637)
            , Element::Hydrogen()
            , 1
            , "ARG"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ORN", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(638)
            , Element::Nitrogen()
            , 3
            , "ORN"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ORN", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(639)
            , Element::Carbon()
            , 4
            , "ORN"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ORN", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(640)
            , Element::Carbon()
            , 3
            , "ORN"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ORN", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(641)
            , Element::Hydrogen()
            , 1
            , "ORN"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ORN", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(642)
            , Element::Oxygen()
            , 1
            , "ORN"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("ORN", "HA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(643)
            , Element::Hydrogen()
            , 1
            , "ORN"
            , "HA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("AIB", "N", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(644)
            , Element::Nitrogen()
            , 3
            , "AIB"
            , "N"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("AIB", "CA", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(645)
            , Element::Carbon()
            , 4
            , "AIB"
            , "CA"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("AIB", "C", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(646)
            , Element::Carbon()
            , 3
            , "AIB"
            , "C"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("AIB", "HN", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(647)
            , Element::Hydrogen()
            , 1
            , "AIB"
            , "HN"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("AIB", "OXT", Ordinality::FINAL) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(648)
            , Element::Oxygen()
            , 1
            , "AIB"
            , "OXT"
            , Ordinality::FINAL
            );

    if (! Biotype::exists("Adenosine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1001)
            , Element::Oxygen()
            , 2
            , "Adenosine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1002)
            , Element::Carbon()
            , 4
            , "Adenosine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1003)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1004)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1005)
            , Element::Carbon()
            , 4
            , "Adenosine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1006)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1007)
            , Element::Oxygen()
            , 2
            , "Adenosine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1008)
            , Element::Carbon()
            , 4
            , "Adenosine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1009)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1010)
            , Element::Carbon()
            , 4
            , "Adenosine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1011)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1012)
            , Element::Carbon()
            , 4
            , "Adenosine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1013)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "O2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1014)
            , Element::Oxygen()
            , 2
            , "Adenosine"
            , "O2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "HO*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1015)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "HO*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1016)
            , Element::Oxygen()
            , 2
            , "Adenosine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "N9", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1017)
            , Element::Nitrogen()
            , 3
            , "Adenosine"
            , "N9"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1018)
            , Element::Carbon()
            , 3
            , "Adenosine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1019)
            , Element::Carbon()
            , 3
            , "Adenosine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "N7", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1020)
            , Element::Nitrogen()
            , 2
            , "Adenosine"
            , "N7"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1021)
            , Element::Carbon()
            , 3
            , "Adenosine"
            , "C8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1022)
            , Element::Nitrogen()
            , 2
            , "Adenosine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1023)
            , Element::Carbon()
            , 3
            , "Adenosine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1024)
            , Element::Nitrogen()
            , 2
            , "Adenosine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1025)
            , Element::Carbon()
            , 3
            , "Adenosine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1026)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "N6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1027)
            , Element::Nitrogen()
            , 3
            , "Adenosine"
            , "N6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H61", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1028)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H61"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H62", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1029)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H62"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Adenosine", "H8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1030)
            , Element::Hydrogen()
            , 1
            , "Adenosine"
            , "H8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1031)
            , Element::Oxygen()
            , 2
            , "Guanosine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1032)
            , Element::Carbon()
            , 4
            , "Guanosine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1033)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1034)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1035)
            , Element::Carbon()
            , 4
            , "Guanosine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1036)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1037)
            , Element::Oxygen()
            , 2
            , "Guanosine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1038)
            , Element::Carbon()
            , 4
            , "Guanosine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1039)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1040)
            , Element::Carbon()
            , 4
            , "Guanosine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1041)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1042)
            , Element::Carbon()
            , 4
            , "Guanosine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1043)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "O2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1044)
            , Element::Oxygen()
            , 2
            , "Guanosine"
            , "O2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "HO*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1045)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "HO*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1046)
            , Element::Oxygen()
            , 2
            , "Guanosine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "N9", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1047)
            , Element::Nitrogen()
            , 3
            , "Guanosine"
            , "N9"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1048)
            , Element::Carbon()
            , 3
            , "Guanosine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1049)
            , Element::Carbon()
            , 3
            , "Guanosine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "N7", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1050)
            , Element::Nitrogen()
            , 2
            , "Guanosine"
            , "N7"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1051)
            , Element::Carbon()
            , 3
            , "Guanosine"
            , "C8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1052)
            , Element::Nitrogen()
            , 2
            , "Guanosine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1053)
            , Element::Carbon()
            , 3
            , "Guanosine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1054)
            , Element::Nitrogen()
            , 3
            , "Guanosine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1055)
            , Element::Carbon()
            , 3
            , "Guanosine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1056)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "N2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1057)
            , Element::Nitrogen()
            , 3
            , "Guanosine"
            , "N2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H21", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1058)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H21"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H22", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1059)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H22"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "O6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1060)
            , Element::Oxygen()
            , 1
            , "Guanosine"
            , "O6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Guanosine", "H8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1061)
            , Element::Hydrogen()
            , 1
            , "Guanosine"
            , "H8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1062)
            , Element::Oxygen()
            , 2
            , "Cytidine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1063)
            , Element::Carbon()
            , 4
            , "Cytidine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1064)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1065)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1066)
            , Element::Carbon()
            , 4
            , "Cytidine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1067)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1068)
            , Element::Oxygen()
            , 2
            , "Cytidine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1069)
            , Element::Carbon()
            , 4
            , "Cytidine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1070)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1071)
            , Element::Carbon()
            , 4
            , "Cytidine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1072)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1073)
            , Element::Carbon()
            , 4
            , "Cytidine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1074)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "O2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1075)
            , Element::Oxygen()
            , 2
            , "Cytidine"
            , "O2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "HO*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1076)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "HO*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1077)
            , Element::Oxygen()
            , 2
            , "Cytidine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1078)
            , Element::Nitrogen()
            , 3
            , "Cytidine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1079)
            , Element::Carbon()
            , 3
            , "Cytidine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1080)
            , Element::Nitrogen()
            , 2
            , "Cytidine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1081)
            , Element::Carbon()
            , 3
            , "Cytidine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1082)
            , Element::Carbon()
            , 3
            , "Cytidine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1083)
            , Element::Carbon()
            , 3
            , "Cytidine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "O2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1084)
            , Element::Oxygen()
            , 1
            , "Cytidine"
            , "O2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "N4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1085)
            , Element::Nitrogen()
            , 3
            , "Cytidine"
            , "N4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H41", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1086)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H41"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H42", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1087)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H42"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1088)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Cytidine", "H6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1089)
            , Element::Hydrogen()
            , 1
            , "Cytidine"
            , "H6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1090)
            , Element::Oxygen()
            , 2
            , "Uridine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1091)
            , Element::Carbon()
            , 4
            , "Uridine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1092)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1093)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1094)
            , Element::Carbon()
            , 4
            , "Uridine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1095)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1096)
            , Element::Oxygen()
            , 2
            , "Uridine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1097)
            , Element::Carbon()
            , 4
            , "Uridine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1098)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1099)
            , Element::Carbon()
            , 4
            , "Uridine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1100)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1101)
            , Element::Carbon()
            , 4
            , "Uridine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1102)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "O2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1103)
            , Element::Oxygen()
            , 2
            , "Uridine"
            , "O2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "HO*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1104)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "HO*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1105)
            , Element::Oxygen()
            , 2
            , "Uridine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1106)
            , Element::Nitrogen()
            , 3
            , "Uridine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1107)
            , Element::Carbon()
            , 3
            , "Uridine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1108)
            , Element::Nitrogen()
            , 3
            , "Uridine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1109)
            , Element::Carbon()
            , 3
            , "Uridine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1110)
            , Element::Carbon()
            , 3
            , "Uridine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1111)
            , Element::Carbon()
            , 3
            , "Uridine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "O2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1112)
            , Element::Oxygen()
            , 1
            , "Uridine"
            , "O2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1113)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "O4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1114)
            , Element::Oxygen()
            , 1
            , "Uridine"
            , "O4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1115)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Uridine", "H6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1116)
            , Element::Hydrogen()
            , 1
            , "Uridine"
            , "H6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1117)
            , Element::Oxygen()
            , 2
            , "Deoxyadenosine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1118)
            , Element::Carbon()
            , 4
            , "Deoxyadenosine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1119)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1120)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1121)
            , Element::Carbon()
            , 4
            , "Deoxyadenosine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1122)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1123)
            , Element::Oxygen()
            , 2
            , "Deoxyadenosine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1124)
            , Element::Carbon()
            , 4
            , "Deoxyadenosine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1125)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1126)
            , Element::Carbon()
            , 4
            , "Deoxyadenosine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1127)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1128)
            , Element::Carbon()
            , 4
            , "Deoxyadenosine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H2*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1129)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H2*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H2*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1130)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H2*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1131)
            , Element::Oxygen()
            , 2
            , "Deoxyadenosine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "N9", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1132)
            , Element::Nitrogen()
            , 3
            , "Deoxyadenosine"
            , "N9"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1133)
            , Element::Carbon()
            , 3
            , "Deoxyadenosine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1134)
            , Element::Carbon()
            , 3
            , "Deoxyadenosine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "N7", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1135)
            , Element::Nitrogen()
            , 2
            , "Deoxyadenosine"
            , "N7"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1136)
            , Element::Carbon()
            , 3
            , "Deoxyadenosine"
            , "C8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1137)
            , Element::Nitrogen()
            , 2
            , "Deoxyadenosine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1138)
            , Element::Carbon()
            , 3
            , "Deoxyadenosine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1139)
            , Element::Nitrogen()
            , 2
            , "Deoxyadenosine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1140)
            , Element::Carbon()
            , 3
            , "Deoxyadenosine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1141)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "N6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1142)
            , Element::Nitrogen()
            , 3
            , "Deoxyadenosine"
            , "N6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H61", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1143)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H61"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H62", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1144)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H62"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyadenosine", "H8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1145)
            , Element::Hydrogen()
            , 1
            , "Deoxyadenosine"
            , "H8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1146)
            , Element::Oxygen()
            , 2
            , "Deoxyguanosine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1147)
            , Element::Carbon()
            , 4
            , "Deoxyguanosine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1148)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1149)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1150)
            , Element::Carbon()
            , 4
            , "Deoxyguanosine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1151)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1152)
            , Element::Oxygen()
            , 2
            , "Deoxyguanosine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1153)
            , Element::Carbon()
            , 4
            , "Deoxyguanosine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1154)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1155)
            , Element::Carbon()
            , 4
            , "Deoxyguanosine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1156)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1157)
            , Element::Carbon()
            , 4
            , "Deoxyguanosine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H2*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1158)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H2*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H2*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1159)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H2*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1160)
            , Element::Oxygen()
            , 2
            , "Deoxyguanosine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "N9", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1161)
            , Element::Nitrogen()
            , 3
            , "Deoxyguanosine"
            , "N9"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1162)
            , Element::Carbon()
            , 3
            , "Deoxyguanosine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1163)
            , Element::Carbon()
            , 3
            , "Deoxyguanosine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "N7", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1164)
            , Element::Nitrogen()
            , 2
            , "Deoxyguanosine"
            , "N7"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1165)
            , Element::Carbon()
            , 3
            , "Deoxyguanosine"
            , "C8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1166)
            , Element::Nitrogen()
            , 2
            , "Deoxyguanosine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1167)
            , Element::Carbon()
            , 3
            , "Deoxyguanosine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1168)
            , Element::Nitrogen()
            , 3
            , "Deoxyguanosine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1169)
            , Element::Carbon()
            , 3
            , "Deoxyguanosine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1170)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "N2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1171)
            , Element::Nitrogen()
            , 3
            , "Deoxyguanosine"
            , "N2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H21", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1172)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H21"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H22", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1173)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H22"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "O6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1174)
            , Element::Oxygen()
            , 1
            , "Deoxyguanosine"
            , "O6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxyguanosine", "H8", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1175)
            , Element::Hydrogen()
            , 1
            , "Deoxyguanosine"
            , "H8"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1176)
            , Element::Oxygen()
            , 2
            , "Deoxycytidine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1177)
            , Element::Carbon()
            , 4
            , "Deoxycytidine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1178)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1179)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1180)
            , Element::Carbon()
            , 4
            , "Deoxycytidine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1181)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1182)
            , Element::Oxygen()
            , 2
            , "Deoxycytidine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1183)
            , Element::Carbon()
            , 4
            , "Deoxycytidine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1184)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1185)
            , Element::Carbon()
            , 4
            , "Deoxycytidine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1186)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1187)
            , Element::Carbon()
            , 4
            , "Deoxycytidine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H2*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1188)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H2*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H2*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1189)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H2*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1190)
            , Element::Oxygen()
            , 2
            , "Deoxycytidine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1191)
            , Element::Nitrogen()
            , 3
            , "Deoxycytidine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1192)
            , Element::Carbon()
            , 3
            , "Deoxycytidine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1193)
            , Element::Nitrogen()
            , 2
            , "Deoxycytidine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1194)
            , Element::Carbon()
            , 3
            , "Deoxycytidine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1195)
            , Element::Carbon()
            , 3
            , "Deoxycytidine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1196)
            , Element::Carbon()
            , 3
            , "Deoxycytidine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "O2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1197)
            , Element::Oxygen()
            , 1
            , "Deoxycytidine"
            , "O2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "N4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1198)
            , Element::Nitrogen()
            , 3
            , "Deoxycytidine"
            , "N4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H41", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1199)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H41"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H42", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1200)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H42"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1201)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxycytidine", "H6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1202)
            , Element::Hydrogen()
            , 1
            , "Deoxycytidine"
            , "H6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1203)
            , Element::Oxygen()
            , 2
            , "Deoxythymidine"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1204)
            , Element::Carbon()
            , 4
            , "Deoxythymidine"
            , "C5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H5*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1205)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H5*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H5*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1206)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H5*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1207)
            , Element::Carbon()
            , 4
            , "Deoxythymidine"
            , "C4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1208)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "O4*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1209)
            , Element::Oxygen()
            , 2
            , "Deoxythymidine"
            , "O4*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1210)
            , Element::Carbon()
            , 4
            , "Deoxythymidine"
            , "C1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H1*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1211)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H1*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1212)
            , Element::Carbon()
            , 4
            , "Deoxythymidine"
            , "C3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1213)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C2*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1214)
            , Element::Carbon()
            , 4
            , "Deoxythymidine"
            , "C2*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H2*1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1215)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H2*1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H2*2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1216)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H2*2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1217)
            , Element::Oxygen()
            , 2
            , "Deoxythymidine"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "N1", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1218)
            , Element::Nitrogen()
            , 3
            , "Deoxythymidine"
            , "N1"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1219)
            , Element::Carbon()
            , 3
            , "Deoxythymidine"
            , "C2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "N3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1220)
            , Element::Nitrogen()
            , 3
            , "Deoxythymidine"
            , "N3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1221)
            , Element::Carbon()
            , 3
            , "Deoxythymidine"
            , "C4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C5", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1222)
            , Element::Carbon()
            , 3
            , "Deoxythymidine"
            , "C5"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1223)
            , Element::Carbon()
            , 3
            , "Deoxythymidine"
            , "C6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "O2", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1224)
            , Element::Oxygen()
            , 1
            , "Deoxythymidine"
            , "O2"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H3", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1225)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H3"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "O4", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1226)
            , Element::Oxygen()
            , 1
            , "Deoxythymidine"
            , "O4"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "C7", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1227)
            , Element::Carbon()
            , 4
            , "Deoxythymidine"
            , "C7"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H7", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1228)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H7"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Deoxythymidine", "H6", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1229)
            , Element::Hydrogen()
            , 1
            , "Deoxythymidine"
            , "H6"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phosphodiester, RNA", "P", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1230)
            , Element::Phosphorus()
            , 4
            , "Phosphodiester, RNA"
            , "P"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phosphodiester, RNA", "OP", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1231)
            , Element::Oxygen()
            , 1
            , "Phosphodiester, RNA"
            , "OP"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Hydroxyl, RNA", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1232)
            , Element::Oxygen()
            , 2
            , "5'-Hydroxyl, RNA"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Hydroxyl, RNA", "H5T", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1233)
            , Element::Hydrogen()
            , 1
            , "5'-Hydroxyl, RNA"
            , "H5T"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Phosphate OS, RNA", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1234)
            , Element::Oxygen()
            , 2
            , "5'-Phosphate OS, RNA"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Phosphate P, RNA", "P", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1235)
            , Element::Phosphorus()
            , 4
            , "5'-Phosphate P, RNA"
            , "P"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Phosphate OP, RNA", "OP", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1236)
            , Element::Oxygen()
            , 1
            , "5'-Phosphate OP, RNA"
            , "OP"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Hydroxyl, RNA", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1237)
            , Element::Oxygen()
            , 2
            , "3'-Hydroxyl, RNA"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Hydroxyl, RNA", "H3T", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1238)
            , Element::Hydrogen()
            , 1
            , "3'-Hydroxyl, RNA"
            , "H3T"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Phosphate OS, RNA", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1239)
            , Element::Oxygen()
            , 2
            , "3'-Phosphate OS, RNA"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Phosphate P, RNA", "P", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1240)
            , Element::Phosphorus()
            , 4
            , "3'-Phosphate P, RNA"
            , "P"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Phosphate OP, RNA", "OP", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1241)
            , Element::Oxygen()
            , 1
            , "3'-Phosphate OP, RNA"
            , "OP"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phosphodiester, DNA", "P", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1242)
            , Element::Phosphorus()
            , 4
            , "Phosphodiester, DNA"
            , "P"
            , Ordinality::ANY
            );

    if (! Biotype::exists("Phosphodiester, DNA", "OP", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1243)
            , Element::Oxygen()
            , 1
            , "Phosphodiester, DNA"
            , "OP"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Hydroxyl, DNA", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1244)
            , Element::Oxygen()
            , 2
            , "5'-Hydroxyl, DNA"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Hydroxyl, DNA", "H5T", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1245)
            , Element::Hydrogen()
            , 1
            , "5'-Hydroxyl, DNA"
            , "H5T"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Phosphate OS, DNA", "O5*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1246)
            , Element::Oxygen()
            , 2
            , "5'-Phosphate OS, DNA"
            , "O5*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Phosphate P, DNA", "P", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1247)
            , Element::Phosphorus()
            , 4
            , "5'-Phosphate P, DNA"
            , "P"
            , Ordinality::ANY
            );

    if (! Biotype::exists("5'-Phosphate OP, DNA", "OP", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1248)
            , Element::Oxygen()
            , 1
            , "5'-Phosphate OP, DNA"
            , "OP"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Hydroxyl, DNA", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1249)
            , Element::Oxygen()
            , 2
            , "3'-Hydroxyl, DNA"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Hydroxyl, DNA", "H3T", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1250)
            , Element::Hydrogen()
            , 1
            , "3'-Hydroxyl, DNA"
            , "H3T"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Phosphate OS, DNA", "O3*", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1251)
            , Element::Oxygen()
            , 2
            , "3'-Phosphate OS, DNA"
            , "O3*"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Phosphate P, DNA", "P", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1252)
            , Element::Phosphorus()
            , 4
            , "3'-Phosphate P, DNA"
            , "P"
            , Ordinality::ANY
            );

    if (! Biotype::exists("3'-Phosphate OP, DNA", "OP", Ordinality::ANY) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeId(1253)
            , Element::Oxygen()
            , 1
            , "3'-Phosphate OP, DNA"
            , "OP"
            , Ordinality::ANY
            );

    popularBiotypesAreInitialized = true;
}


} // namespace SimTK
