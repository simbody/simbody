/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
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
    Xml::setXmlCondenseWhiteSpace(false);
    SimTK_TEST(!Xml::isXmlWhiteSpaceCondensed());
    preserveWhite.readFromString(xmlPlainTextFile);
    cout << "Plain text file with white space preserved (raw): "
         << preserveWhite.getRootElement().getValue() << "\n";
    cout << "... (formatted with condense=false): "
         << preserveWhite << "\n";
    Xml::setXmlCondenseWhiteSpace(true);
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

    Xml scratch2;
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


int main() {
    cout << "Path of this executable: '" << Pathname::getThisExecutablePath() << "'\n";
    cout << "Executable directory: '" << Pathname::getThisExecutableDirectory() << "'\n";
    cout << "Current working directory: '" << Pathname::getCurrentWorkingDirectory() << "'\n";


    SimTK_START_TEST("TestXml");

        SimTK_SUBTEST(testStringConvert);
        SimTK_SUBTEST(testXmlFromString);
        SimTK_SUBTEST(testXmlFromScratch);

    SimTK_END_TEST();
}

