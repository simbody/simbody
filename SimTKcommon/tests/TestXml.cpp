/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-16 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Chris Dembia                                                 *
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
#include "SimTKcommon/Testing.h"

#include "SimTKcommon/internal/Xml.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <type_traits>


using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::printf;

using namespace SimTK;

// This example is from Wikipedia's XML entry.
const char* xmlPainting = 
"<?xml version='1.0' encoding='UTF-8'?>\n"
"<!-- a top-level comment -->\n"
"<!-- a multiline\n comment\n third line -->\n"
"     \n" // line should be ignored
"but something like this is top level text and will need to get \n"
"moved into a new '_Root' element\n"
"<painting artist='Raphael' artist='metoo'>\n"
"  <img src=\"madonna.jpg\" alt='Foligno Madonna, by Raphael'/>\n"
"  <!-- What follows is a so-called 'caption' -->\n"
"  <caption>This is Raphael's \"Foligno\" Madonna, painted in\n"
"    <date>1511</date>-<date>1512</date>.\n"
"    <![CDATA[some non-Unicode text]]>\n"
"    <  !SOMETHING but this tag is unknown>  \n"
"  </caption>\n"
"  This part is just plain old text.\n"
"  As is this \"quoted\" thing.\n"
"</painting>\n"
"<!whazzis unknown junk>\n"
"<!-- final comment -->";

const char* xmlPlainTextFile =
"That is, the first line should be a declaration, most commonly exactly the\n"
"characters shown above, without the \"standalone\" attribute which will\n" 
"default to \"yes\" anyway. If we don't see a declaration when reading an XML\n"
"document, we'll assume we read the one above. Then the document should contain\n" 
"exactly one top-level (root) element representing the type of document &amp;\n" 
"document-level attributes.\n";

const char* xmlUnclosedComment = 
"<?xml version='1.0' encoding='UTF-8'?>\n"
"  <!-- What follows is a so-called 'caption' ->\n" // UNCLOSED!
"<!whazzis unknown junk>";

const char* xmlJustAComment = 
"<!-- this is the entire contents -->\n";

const char* xmlEmpty = "   \n \n \t "; // white space only


static void showElement(Xml::Element elt, const String& indent="") {
    cout << indent << "ELEMENT WITH TAG '" << elt.getElementTag() << "':\n";

    // Show attributes
    Xml::attribute_iterator ap = elt.attribute_begin();
    for (; ap != elt.attribute_end(); ++ap)
        cout << indent << "  ATTR '" << ap->getName() << "'='" 
             << ap->getValue() << "'\n";

    // Show all contents
    Xml::node_iterator p = elt.node_begin();
    for (; p != elt.node_end(); ++p) {
        cout << indent << p->getNodeTypeAsString() << endl;
        if (p->getNodeType() == Xml::ElementNode)
            showElement(Xml::Element::getAs(*p), indent + "  ");
    }
    cout << indent << "END OF ELEMENT.\n";
}

void testXmlFromString() {
    Xml::Document fromString;
    fromString.readFromString(xmlJustAComment);
    cout << "Just a comment: '" << fromString << "'\n";

    fromString.readFromString(xmlPlainTextFile);
    cout << "Plain text file: '" << fromString << "'\n";

    // Note that the "condense white space" setting is global, not 
    // document-specific.
    Xml::Document preserveWhite;
    Xml::Document::setXmlCondenseWhiteSpace(false);
    SimTK_TEST(!Xml::Document::isXmlWhiteSpaceCondensed());
    preserveWhite.readFromString(xmlPlainTextFile);
    cout << "Plain text file with white space preserved (raw): " 
         << preserveWhite.getRootElement().getValue() << "\n";
    cout << "... (formatted with condense=false): " 
         << preserveWhite << "\n";    
    Xml::Document::setXmlCondenseWhiteSpace(true);
    cout << "... (formatted with condense=true): " 
         << preserveWhite << "\n";    

    SimTK_TEST_MUST_THROW(fromString.readFromString(xmlEmpty));
    SimTK_TEST_MUST_THROW(fromString.readFromString(xmlUnclosedComment));

    fromString.readFromString(String(xmlPainting));
    cout << "Painting: '" << fromString << "'\n";

    cout << "Doc type is: " << fromString.getRootTag() << endl;
    cout << "  version: " << fromString.getXmlVersion() << endl;
    cout << "  encoding: " << fromString.getXmlEncoding() << endl;
    cout << "  standalone: " 
         << String(fromString.getXmlIsStandalone() ? "yes" : "no") << endl;

    cout << "All nodes in doc:\n";
    Xml::node_iterator np = fromString.node_begin();
    for (; np != fromString.node_end(); ++np)
        cout << "  " << np->getNodeTypeAsString() << endl;


    Xml::Element root = fromString.getRootElement();

    cout << "hasNode()=" << root.hasNode() << endl;
    cout << "hasNode(Comment)=" << root.hasNode(Xml::CommentNode) << endl;
    cout << "hasNode(Unknown)=" << root.hasNode(Xml::UnknownNode) << endl;

    showElement(root);

    Xml::node_iterator p = root.node_begin(Xml::NoJunkNodes);
    for (; p != root.node_end(); ++p)
        cout << p->getNodeTypeAsString() << endl;

    cout << "Caption elements:\n";
    Xml::element_iterator ep(root.element_begin("caption"));
    for (; ep != root.element_end(); ++ep)
        cout << ep->getNodeTypeAsString() << ": <" << ep->getElementTag() 
             << ">" << endl;

    cout << "All elements:\n";
    Array_<Xml::Element> all = root.getAllElements();
    for (unsigned i=0; i < all.size(); ++i)
        cout << "<" << all[i].getElementTag() << ">" << endl;


    Array_<Xml::Node> allNodes(root.node_begin(Xml::NoJunkNodes), 
                               root.node_end());
    for (unsigned i=0; i < allNodes.size(); ++i)
        cout << "Node " << allNodes[i].getNodeTypeAsString() << endl;

    String prettyString, compactString;
    fromString.writeToString(prettyString);
    cout << "String pretty: " << prettyString.size() 
         << "\n'" << prettyString << "'\n";

    fromString.writeToString(compactString, true); // compact
    cout << "String compact: " << compactString.size() 
         << "\n'" << compactString << "'\n";

    SimTK_TEST(compactString.size() < prettyString.size());

    cout << "painting.allNode=" << root.getRequiredElement("painting")
                                        .getAllNodes() << endl;
    cout << "painting.img.allAttr=" << 
        root.getRequiredElement("painting").getRequiredElement("img")
        .getAllAttributes() << endl;
    //fromString.writeToFile("TestXml.xml");

    //Xml ex("TestXml.xml");
    //cout << "Document tag: " << ex.getDocumentTag() << endl;
    //for (Xml::nodu;e_iterator xp=ex.node_begin(); xp != ex.node_end(); ++xp)
    //    cout << "Node type: " << xp->getNodeTypeAsString() << endl;

    //PolygonalMesh mesh;
    //mesh.loadVtpFile("arm_r_humerus.vtp");
    //cout << "num vertices=" << mesh.getNumVertices() << " faces="
    //    << mesh.getNumFaces() << endl;
}

void testXmlFromScratch() {
    Xml::Document scratch;
    scratch.setRootTag("MyDoc");
    cout << scratch;

    Xml::Comment c("This is a comment.");
    Xml::Unknown u("!GODONLY knows what this is!!");
    Xml::Text t("This is some\ntext on two lines with trailing blanks   ");
    Xml::Element e("elementTag");

    // We're never going to use this one so its heap space will
    // leak if we don't explicitly call clearOrphan().
    Xml::Element neverMind("neverMind");

    cout << "initially e='" << e.getValue() << "'" << endl;
    e.updValue() += "AVALUE:";
    cout << "then e='" << e.getValue() << "'" << endl;

    e.setAttributeValue("attr1", String(Vec2(9,-9)));

    cout << "attr1 is " << e.getRequiredAttributeValueAs<Vec2>("attr1") << endl;

    cout << "isOrphan? " << String(c.isOrphan()) << ":" << c;
    cout << "isOrphan? " << String(u.isOrphan()) << ":" << u;
    cout << "isOrphan? " << String(t.isOrphan()) << ":" << t;
    cout << "isOrphan? " << String(e.isOrphan()) << ":" << e;

    e.setValue("this is the only value");
    e.updValue() += " (but then I added this)";
    cout << "e value=" << e.getValue() << endl;
    cout << "e = " << e << endl;
    e.setValue("9 10 -3.2e-4");
    cout << "e value=" << e.getValueAs< Array_<float> >() << endl;
    cout << "e = " << e << endl;

    scratch.insertTopLevelNodeAfter(scratch.node_begin(), c);
    cout << "isOrphan? " << String(c.isOrphan()) << ":" << c;
    cout << scratch;

    scratch.insertTopLevelNodeBefore(scratch.node_begin(Xml::ElementNode), u);
    cout << "isOrphan? " << String(u.isOrphan()) << ":" << u;
    cout << scratch;

    scratch.insertTopLevelNodeBefore(scratch.node_begin(), 
        Xml::Comment("This should be at the top of the file, except declaration."));
    cout << scratch;

    Xml::Document scratch2;
    scratch2 = scratch; // deep copy

    scratch.eraseTopLevelNode(scratch.node_begin());
    cout << "First node gone (scratch)?\n" << scratch;
    cout << "First node still there (scratch2)?\n" << scratch2;

    Xml::Element e2("anotherElt", Vec3(.1,.2,.3));
    cout << e2;
    e.insertNodeAfter(e.element_end(), e2);
    cout << "now owns anotherElt:\n" << e;

    Xml::Element root = scratch.getRootElement();
    root.insertNodeAfter(root.node_begin(), e);
    cout << scratch;

    scratch.setIndentString("..");
    cout << scratch;

    Xml::Element ecopy = e.clone();
    ecopy.setElementTag("elementTagCopy");
    cout << "COPY of e: " << ecopy;

    //scratch.writeToFile("scratch.xml");

    e.eraseNode(e.element_begin("anotherElt"));
    cout << "in-place removal of anotherElt from e: " << scratch;
    cout << "COPY of e: " << ecopy;

    root.insertNodeAfter(root.element_begin("elementTag"), ecopy);
    cout << "After copy insert, scratch=" << scratch;

    Xml::Node extract = root.removeNode(root.element_begin("elementTagCopy"));
    cout << "Extracted copy: " << extract << endl;
    cout << "Now scratch=" << scratch << endl;
    cout << "Duplicate scratch=" << Xml::Document(scratch) << endl;

    neverMind.clearOrphan();
}

#define SHOWIT(something) \
    cout << String(#something) << String("'") << something << String("'\n")

void testStringConvert() {
    SimTK_TEST(convertStringTo<int>(" 239\n ")==239);
    SimTK_TEST(convertStringTo<String>("  lunch box\n") == "  lunch box\n");
    SimTK_TEST(convertStringTo<std::string>("  lunch box\n") == "  lunch box\n");
    SimTK_TEST(convertStringTo<int>("1234")==1234);
    SimTK_TEST(convertStringTo<unsigned>("01234")==1234);
    SimTK_TEST(convertStringTo<float>("1234.5")==1234.5);
    SimTK_TEST_MUST_THROW(convertStringTo<char*>("  lunch box\n"));
    SimTK_TEST_MUST_THROW(convertStringTo<int>(" 234 j"));
    SimTK_TEST_MUST_THROW(convertStringTo<int>("345.5"));

    SimTK_TEST(convertStringTo< std::complex<double> >("(-4,22)")
               ==std::complex<double>(-4,22));
    SimTK_TEST(convertStringTo< Vec3 >("1 2 3") == Vec3(1,2,3));
    SimTK_TEST(convertStringTo< Vec3 >("1, 2 , 3") == Vec3(1,2,3));
    SimTK_TEST(convertStringTo< Vec3 >("[ -3 , 5, 6 ] ")== Vec3(-3,5,6));
    SimTK_TEST(convertStringTo< Vec3 >("( -3  5 -6 ) ")== Vec3(-3,5,-6));
    SimTK_TEST(convertStringTo< Vec3 >(" ~ [ -3 , 5, 6 ] ")== Vec3(-3,5,6));
    SimTK_TEST(convertStringTo< Vec3 >("~( -3  5 -6 ) ")== Vec3(-3,5,-6));
    SimTK_TEST_MUST_THROW(convertStringTo< Vec3 >("( -3  5 -6 ] "));
    SimTK_TEST_MUST_THROW(convertStringTo< Vec3 >(" -3  5 -6 ] "));
    SimTK_TEST_MUST_THROW(convertStringTo< Vec3 >(" ~ -3  5 -6 "));
    typedef Vec<2,std::complex<float> > fCVec2;
    SimTK_TEST(convertStringTo<fCVec2>("[(1,2) (3,4)]")
        == fCVec2(std::complex<float>(1,2), std::complex<float>(3,4)));

    Array_<int> a = convertStringTo< Array_<int> >("1 2 3 4");

    Array_<float> af(2);
    SimTK_TEST_MUST_THROW( // because ArrayView_ is fixed size (2)
        String(" -.25, .5, 29.2e4 ").convertTo<ArrayView_<float> >(af));
    // But this should work because an Array_ can be resized.
    String(" -.25, .5, 29.2e4 ").convertTo<Array_<float> >(af);
    SimTK_TEST(af[0]==-.25 && af[1]==.5 && af[2]==292000);

}

// This class has some serialization issues to deal with. It has several 
// members of different types including a templatized one. It uses both
// member and namespace-scope functions for serialization. It defines its
// own serialization tag "TryMe" rather than letting the serialization system
// use its own tag.
// 
class TryMe {
public:
    TryMe() = default;
    //TryMe() {cout << "TryMe()" << endl;}

    TryMe(const TryMe& s) = default;
    //TryMe(const TryMe& s) : i(s.i), d(s.d), v(s.v), avs(s.avs) {
    //    cout << "TryMe copy constructor" << endl;
    //}

    TryMe& operator=(const TryMe& s) = default;
    //TryMe& operator=(const TryMe& s) {
    //    cout << "TryMe copy assignment" << endl;
    //    i = s.i; d=s.d; v=s.v; avs=s.avs;
    //    return *this;
    //}

    TryMe(TryMe&&) = default;
    //TryMe(TryMe&& s) : i(s.i), d(s.d), v(s.v), avs(s.avs) {
    //    cout << "TryMe move constructor" << endl;
    //}

    TryMe& operator=(TryMe&& s) = default;
    //TryMe& operator=(TryMe&& s) {
    //    cout << "TryMe move assignment" << endl;
    //    i = s.i; d=s.d; v=s.v; avs=s.avs;
    //    return *this;
    //}

    // A non-default constructor.
    TryMe(int ii,double dd,Vec3 vv,Array_<Vec2,short> aa)
        : i(ii),d(dd),v(vv),avs(aa) {
        cout << "TryMe member constructor" << endl;
    }

    // toXmlElement() implemented as a class member; this should be
    // preferred over a namespace-scope implementation. Note that we have our
    // own tag so the given name becomes an attribute rather than the tag.
    Xml::Element toXmlElement(const string& name) const {
        static const int version = 26; // <- member function uses 26
        Xml::Element e("TryMe");
        if (!name.empty()) e.setAttributeValue("name", name);
        e.setAttributeValue("version", String(version));
        e.setAttributeValue("free", String(false));
        e.appendNode(toXmlElementHelper(i,  "i"  , true));
        e.appendNode(toXmlElementHelper(d,  "d"  , true));
        e.appendNode(toXmlElementHelper(v,  "v"  , true));
        e.appendNode(toXmlElementHelper(avs,"avs", true));
        return e;
    }

    int                 i;
    double              d;
    Vec3                v;
    Array_<Vec2,short>  avs;
};

// Namespace-scope method suitable for a TryMe, but this shouldn't get
// called since the member function should be preferred.
Xml::Element toXmlElement(const TryMe& t, const string& name) {
    static const int version = 25; // <- free function uses 25
    Xml::Element e("TryMe");
    if (!name.empty()) e.setAttributeValue("name", name);
    e.setAttributeValue("version", String(version));
    e.setAttributeValue("free", String(true));
    e.appendNode(toXmlElement(t.i, "i"));
    e.appendNode(toXmlElement(t.d,"d"));
    e.appendNode(toXmlElement(t.v,"v"));
    e.appendNode(toXmlElement(t.avs,"avs"));
    return e;
}

// Namespace-scope method suitable for a TryMe. In this case there is no
// member function available so this should get called. Note that we're 
// expecting a particular tag; the name if any is an attribute.
void fromXmlElement(TryMe& n, Xml::Element& e, const string& reqName) {
    const int expVersion = reqName=="exp25" ? 25 : 26;
    SimTK_ERRCHK1_ALWAYS(e.getElementTag()=="TryMe", 
        "fromXmlElement<TryMe>",
        "Expected tag 'TryMe' but got '%s'.", e.getElementTag().c_str());
    if (!reqName.empty()) {
        const String& name = e.getRequiredAttributeValue("name");
        SimTK_ERRCHK2_ALWAYS(name==reqName,
        "fromXmlElement<TryMe>",
        "Expected TryMe named '%s' but got '%s'.", reqName.c_str(), name.c_str());
    }
    auto gotVersion = e.getRequiredAttributeValueAs<int>("version");
    SimTK_ERRCHK2_ALWAYS(gotVersion == expVersion, "fromXmlElement(TryMe)",
        "Expected version %d but got %d.", expVersion, gotVersion);

    auto nxt = e.element_begin();
    fromXmlElement(n.i, *nxt, "i");
    fromXmlElement(n.d, *++nxt, "d");
    fromXmlElement(n.v, *++nxt, "v");
    fromXmlElement(n.avs, *++nxt, "avs");
}

// Define an enum in a namespace and see if we can (de)serialize it properly.
namespace Nork7 {
enum class Abc {Red, Green, Blue};

// Stream extracter and inserter operators are needed for working directly
// with XML. If they are available they are used as default serialization
// methods also. That is, a default toXmlElement() and fromXmlElement() method
// will be created that uses them. 
// So in this case the explicit to/fromXmlElement() methods below are optional.
std::ostream& operator<<(std::ostream& o, Nork7::Abc color) {
    using ut = typename std::underlying_type<Nork7::Abc>::type;
    switch (color) {
    case Abc::Red:   o << "Red";   break;
    case Abc::Green: o << "Green"; break;
    case Abc::Blue:  o << "Blue";  break;
    default: o << "BAD Nork7::Abc value " << ut(color);
    }
    return o;
}

std::istream& operator>>(std::istream& i, Nork7::Abc& color) {
    std::string s;
    i >> s;
    if (s=="Red")        color = Nork7::Abc::Red;
    else if (s=="Green") color = Nork7::Abc::Green;
    else if (s=="Blue")  color = Nork7::Abc::Blue;
    else {
        SimTK_ERRCHK1_ALWAYS(!"illegal Abc",
                             "operator>>(Nork7::Abc)",
                             "Input had illegal value %s.",s.c_str());
    }
    return i;
}

// Here is an explicit, namespace-scope toXmlElement method that knows how
// output a Nork7::Abc. This is optional since the default constructed from
// the stream insertion operator would do the same thing -- try commenting it
// out and re-running the test.
Xml::Element
toXmlElement(Nork7::Abc v, const string& nm) {
    Xml::Element e(nm.empty() ? NiceTypeName<Nork7::Abc>::namestr()
                              : nm);
    switch (v) {
    case Nork7::Abc::Red:   e.setValue("Red");   break;
    case Nork7::Abc::Green: e.setValue("Green"); break;
    case Nork7::Abc::Blue:  e.setValue("Blue");  break;
    default: 
        SimTK_ERRCHK1_ALWAYS(!"illegal Abc",
                             "toXmlElement(Nork7::Abc)",
                             "Enum had illegal value %d.",(int)v);
    }
    return e;
}

// Expect a value element like <Tag>Red</Tag>. If we're given a required
// name here it is referring to the tag since we don't have our own designated 
// tag for this type.
// This is optional here since it does the same thing the default method
// would do using the provided stream extraction operator -- try commenting it
// out and re-running the test.
void fromXmlElement(Nork7::Abc& v, Xml::Element& e, const string& reqName) {
    if (!reqName.empty()) {
        const String& tag = e.getElementTag();
        SimTK_ERRCHK2_ALWAYS(tag==reqName,
          "fromXmlElement(Nork7::Abc)",
          "Expected Abc value element with tag '%s' but got '%s'.", 
          reqName.c_str(), tag.c_str());
    }
    const auto& color = e.getValue();
         if (color=="Red")   v = Nork7::Abc::Red;
    else if (color=="Green") v = Nork7::Abc::Green; 
    else if (color=="Blue")  v = Nork7::Abc::Blue;
    else {
        SimTK_ERRCHK1_ALWAYS(!"bad enum value",
                             "fromXmlElement(Nork7::abc)",
            "Unrecognized value for enumeration Nork7::Abc: '%s'",
            color.c_str());
    }
}
} // namespace Nork7


using XX = Nork7::Abc;

// Try (de)serializing Value<T> objects to and from Xml.
void testValueSerialization() {
    Value<long long> myInt(5);
    Xml::Element top("TopLevel"); // <TopLevel> ... </TopLevel>

    top.appendNode(myInt.toXmlElement(""));
    //<Value type="long long">
    //    <thing>5</thing>
    //</Value>

    // Value's toXmlElement should eventually find the one above.
    top.appendNode(Value<XX>(XX::Green).toXmlElement("color1"));
    //<Value type="Nork7::Abc" name="color1">
    //    <thing>Green</thing>
    //</Value>

    // This should also find its way to the above method.
    top.appendNode(toXmlElementHelper(Value<XX>(XX::Red), "color2", true));
    //<Value type="Nork7::Abc" name="color2">
    //    <thing>Red</thing>
    //</Value>

    cout << top << endl;
    SimTK_TEST(top.getElementTag() == "TopLevel");
    auto tels = top.getAllElements();
    SimTK_TEST(tels[0].getElementTag() == "Value");
    SimTK_TEST(tels[0].getRequiredAttributeValue("type") == "long long");
    // doesn't have a "name" attribute
    SimTK_TEST(tels[0].getRequiredElementValueAs<long long>("thing") == 5);

    SimTK_TEST(tels[1].getElementTag() == "Value");
    SimTK_TEST(tels[1].getRequiredAttributeValue("type") == "Nork7::Abc");
    SimTK_TEST(tels[1].getRequiredAttributeValue("name") == "color1");
    SimTK_TEST(tels[1].getRequiredElementValueAs<Nork7::Abc>("thing") 
               == Nork7::Abc::Green);

    SimTK_TEST(tels[2].getElementTag() == "Value");
    SimTK_TEST(tels[2].getRequiredAttributeValue("type") == "Nork7::Abc");
    SimTK_TEST(tels[2].getRequiredAttributeValue("name") == "color2");
    SimTK_TEST(tels[2].getRequiredElementValueAs<Nork7::Abc>("thing") 
               == Nork7::Abc::Red);

    // Now see if we can deserialize these into AbstractValues.
    std::unique_ptr<AbstractValue> ll = 
        AbstractValue::createFromXmlElement(tels[0], "");
    std::unique_ptr<AbstractValue> c1 = 
        AbstractValue::createFromXmlElement(tels[1], "color1");
    std::unique_ptr<AbstractValue> c2 = 
        AbstractValue::createFromXmlElement(tels[2], ""); // allow any name

    // Check the actual types.
    SimTK_TEST(dynamic_cast<Value<long long>*>(ll.get()) != nullptr)
    SimTK_TEST(dynamic_cast<Value<Nork7::Abc>*>(c1.get()) != nullptr)
    SimTK_TEST(dynamic_cast<Value<Nork7::Abc>*>(c2.get()) != nullptr)

    SimTK_TEST(ll->getValue<long long>() == 5);
    SimTK_TEST(c1->getValue<Nork7::Abc>() == Nork7::Abc::Green);
    SimTK_TEST(c2->getValue<Nork7::Abc>() == Nork7::Abc::Red);

    cout << "ll=" << *ll << endl;
    cout << "c1=" << *c1 << endl;
    cout << "c2=" << *c2 << endl;


    std::unique_ptr<AbstractValue> intValp(new Value<int>(1234));
    Xml::Element ivp = intValp->toXmlElement("int_val_name");
    SimTK_TEST(ivp.getElementTag() == "Value");
    SimTK_TEST(ivp.getRequiredAttributeValue("type") == "int");
    SimTK_TEST(ivp.getRequiredAttributeValue("name") == "int_val_name");

    std::unique_ptr<AbstractValue> intValp2 = 
        AbstractValue::createFromXmlElement(ivp, ""); // allow any name

    // Now turn it back into an element with a new name.
    Xml::Element ivp2 = intValp2->toXmlElement("int_val_name2");
    cout << "ivp2=" << ivp2 << endl;

    SimTK_TEST(ivp2.getElementTag() == "Value");
    SimTK_TEST(ivp2.getRequiredAttributeValue("type") == "int");
    SimTK_TEST(ivp2.getRequiredAttributeValue("name") == "int_val_name2");
    SimTK_TEST(ivp2.getRequiredElementValueAs<int>("thing") == 1234);


    Value<TryMe> myTry(TryMe{14,3.14,Vec3(1,2,3),
                             {Vec2(.1,.2),Vec2(3,4)}});

    // Add some more elements to the TopLevel element we built above.

    // Use the member function (version=26);
    top.appendNode(myTry.toXmlElement("myTry"));

    // This should invoke the namespace-scope method (version=25);
    top.appendNode(toXmlElement(TryMe{-3,2.7,Vec3(2,NaN,Infinity),
                                      {Vec2(NaN,3),Vec2(-Infinity,Infinity)}},
                              "exp25"));

    // This uses Value<TryMe> so should end up at the member function again.
    top.appendNode(
        Value<TryMe>(TryMe{-3,2.7,Vec3(2,NaN,Infinity),
                           {Vec2(NaN,3),Vec2(-Infinity,Infinity)}})
        .toXmlElement("exp26"));


    cout << top << endl;

    auto ee = top.getRequiredElement("TryMe");
    TryMe nw;
    fromXmlElement(nw, ee, "exp25");

    Value<TryMe> vtm;
    auto ep = top.element_begin("Value");
    while (ep->getOptionalAttributeValue("name", "NONAME") != "exp26")
        ++ep;
    vtm.fromXmlElement(*ep, "exp26");

    top.appendNode(toXmlElement(nw, "exp25Copy"));
    cout << top << endl;

    Xml::Element eltMyNo = myTry.toXmlElement("myTry2");
    auto vp = AbstractValue::createFromXmlElement(eltMyNo, "");
    cout << "vp=" << vp->toXmlElement("vp") << endl;


}



class ClassWithAbstractValue {
public:
    Real m_number;
    ClonePtr<AbstractValue> m_value;
    Xml::Element toXmlElement(const string& name) const {
        static const int version = 1;
        Xml::Element e("ClassWithAbstractValue");
        if (!name.empty()) e.setAttributeValue("name", name);
        e.setAttributeValue("version", String(version));
        e.appendNode(toXmlElementHelper(m_number,  "number" , true));
        e.appendNode(toXmlElementHelper(*m_value,  "value"  , true));
        return e;
    }
    void fromXmlElement(Xml::Element e, const string& requiredName = "") {
        // This helper knows that `m_value` is a smart pointer and will reset
        // it for us.
        fromXmlElementHelperHelper("ClassWithAbstractValue", 1, e, requiredName,
                std::make_pair(&m_number, "number"),
                std::make_pair(&m_value, "value"));
    }
};

void testFromXmlElementHelperHelper() {
    ClassWithAbstractValue obj;
    obj.m_number = 6.50;
    obj.m_value.reset(new Value<Vec3>(Vec3(6, 0, 5)));
    Xml::Element objXml = toXmlElementHelper(obj, "bar", true);
    
    ClassWithAbstractValue deserialized;
    fromXmlElementHelper(deserialized, objXml, "bar", true);

    SimTK_TEST(deserialized.m_number == 6.50);
    SimTK_TEST(deserialized.m_value->getValue<Vec3>() == Vec3(6, 0, 5));
}


class HasNoSerialization {
};

// Ensure that the error one gets when a type cannot be serialized is helpful.
void testToXmlElementException() {
    HasNoSerialization obj;
    SimTK_TEST_MUST_THROW_EXC(toXmlElementHelper(obj, "Alfred", true),
                              Exception::Cant);
    SimTK_TEST_MUST_THROW_SHOW(toXmlElementHelper(obj, "Alfred", true));
}


SimTK_DEFINE_UNIQUE_INDEX_TYPE(FooIndex);

class Outer {
public:
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Outer, LocalIndex);
};

void testXmlUniqueIndexType() {
    FooIndex fi0(73);
    Xml::Element fiXml = toXmlElementHelper(fi0, "FooIndexName", true);
    FooIndex fi1;
    fromXmlElementHelper(fi1, fiXml, "FooIndexName", true);
    SimTK_TEST(fi1 == 73);
    
    Outer::LocalIndex oli0(41);
    Xml::Element oliXml = toXmlElementHelper(oli0, "OuterInnerIndexName", true);
    Outer::LocalIndex oli1;
    fromXmlElementHelper(oli1, oliXml, "OuterInnerIndexName", true);
    SimTK_TEST(oli1 == 41);
}


int main() {
    cout << "Path of this executable: '" << Pathname::getThisExecutablePath() << "'\n";
    cout << "Executable directory: '" << Pathname::getThisExecutableDirectory() << "'\n";
    cout << "Current working directory: '" << Pathname::getCurrentWorkingDirectory() << "'\n";

    SimTK_START_TEST("TestXml");

        SimTK_SUBTEST(testStringConvert);
        SimTK_SUBTEST(testXmlFromString);
        SimTK_SUBTEST(testXmlFromScratch);
        SimTK_SUBTEST(testValueSerialization);
        SimTK_SUBTEST(testFromXmlElementHelperHelper);
        SimTK_SUBTEST(testToXmlElementException);
        SimTK_SUBTEST(testXmlUniqueIndexType);
    
    SimTK_END_TEST();
}

