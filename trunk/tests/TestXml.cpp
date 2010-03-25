/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
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

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include "SimTKcommon/internal/Xml.h"

#include <iostream>
#include <string>
#include <cstdio>
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
"     \n" // line should be ignored
"but something like this is top level text and will need to get \n"
"moved into a new 'XMLDocument' element\n"
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


static void showElement(const Xml::Element& elt, const String& indent="") {
    cout << indent << "ELEMENT WITH TAG '" << elt.getElementTag() << "':\n";

    // Show attributes
    Xml::const_attribute_iterator ap = elt.attribute_begin();
    for (; ap != elt.attribute_end(); ++ap)
        cout << indent << "  ATTR '" << ap->getName() << "'='" 
             << ap->getValue() << "'\n";

    // Show all contents
    Xml::const_node_iterator p = elt.node_begin();
    for (; p != elt.node_end(); ++p) {
        cout << indent << p->getNodeTypeAsString() << ": " << p->getValue() << endl;
        if (p->getNodeType() == Xml::ElementNode)
            showElement(Xml::Element::getAs(*p), indent + "  ");
    }
    cout << indent << "END OF ELEMENT.\n";
}

void testXmlFromString() {
    Xml fromString;
    fromString.readFromString(xmlJustAComment);
    cout << "Just a comment: '" << fromString << "'\n";

    fromString.readFromString(xmlPlainTextFile);
    cout << "Plain text file: '" << fromString << "'\n";
    cout << "'"
        << fromString.getDocumentElement().node_begin(Xml::TextNode)->getValue()
        << "'\n";

    SimTK_TEST_MUST_THROW(fromString.readFromString(xmlEmpty));
    SimTK_TEST_MUST_THROW(fromString.readFromString(xmlUnclosedComment));

    fromString.readFromString(String(xmlPainting));
    cout << "Painting: '" << fromString << "'\n";

    cout << "Doc type is: " << fromString.getDocumentTag() << endl;
    cout << "  version: " << fromString.getXmlVersion() << endl;
    cout << "  encoding: " << fromString.getXmlEncoding() << endl;
    cout << "  standalone: " 
         << String(fromString.getXmlIsStandalone() ? "yes" : "no") << endl;

    cout << "All nodes in doc:\n";
    Xml::const_node_iterator np = fromString.node_begin();
    for (; np != fromString.node_end(); ++np)
        cout << "  " << np->getNodeTypeAsString() << ": " << np->getValue() << endl;


    const Xml::Element& doc = fromString.getDocumentElement();
    Xml::Element& wdoc = fromString.updDocumentElement();

    cout << "empty()=" << doc.empty() << endl;
    cout << "empty(Comment)=" << doc.empty(Xml::CommentNode) << endl;
    cout << "empty(Unknown)=" << doc.empty(Xml::UnknownNode) << endl;

    showElement(doc);

    Xml::const_node_iterator p = doc.node_begin(Xml::NoJunkNodes);
    for (; p != doc.node_end(); ++p)
        cout << p->getNodeTypeAsString() << ": " << p->getValue() << endl;

    cout << "Caption elements:\n";
    Xml::element_iterator ep(wdoc.element_begin("caption"));
    for (; ep != wdoc.element_end(); ++ep)
        cout << ep->getNodeTypeAsString() << ": " << ep->getValue() << endl;

    cout << "All elements:\n";
    ep = wdoc.element_begin();
    for (; ep != wdoc.element_end(); ++ep)
        cout << ep->getNodeTypeAsString() << ": " << ep->getValue() << endl;

    Array_<Xml::Node> allNodes(wdoc.node_begin(Xml::NoJunkNodes), wdoc.node_end());
    for (unsigned i=0; i < allNodes.size(); ++i)
        cout << "Node " << allNodes[i].getNodeTypeAsString()
        << ": " << allNodes[i].getValue() << endl;

    String prettyString, compactString;
    fromString.writeToString(prettyString);
    cout << "String pretty: " << prettyString.size() 
         << "\n'" << prettyString << "'\n";

    fromString.writeToString(compactString, true); // compact
    cout << "String compact: " << compactString.size() 
         << "\n'" << compactString << "'\n";

    SimTK_TEST(compactString.size() < prettyString.size());

    //fromString.writeToFile("TestXml.xml");

    //Xml ex("TestXml.xml");
    //cout << "Document tag: " << ex.getDocumentTag() << endl;
    //for (Xml::node_iterator xp=ex.node_begin(); xp != ex.node_end(); ++xp)
    //    cout << "Node type: " << xp->getNodeTypeAsString() << endl;

    //PolygonalMesh mesh;
    //mesh.loadVtpFile("arm_r_humerus.vtp");
    //cout << "num vertices=" << mesh.getNumVertices() << " faces="
    //    << mesh.getNumFaces() << endl;
}

void testXmlFromScratch() {
    Xml scratch;
    scratch.setDocumentTag("MyDoc");
    cout << scratch;
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

int main() {
    cout << "Path of this executable: '" << Pathname::getThisExecutablePath() << "'\n";
    cout << "Executable directory: '" << Pathname::getThisExecutableDirectory() << "'\n";
    cout << "Current working directory: '" << Pathname::getCurrentWorkingDirectory() << "'\n";

    SimTK_START_TEST("TestXml");

        SimTK_SUBTEST(testStringConvert);
        SimTK_SUBTEST(testXmlFromScratch);
        SimTK_SUBTEST(testXmlFromString);

    SimTK_END_TEST();
}

