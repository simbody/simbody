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

// February 2024
// The ability to specify output precision was added to three methods
// in the SimTK API. In particular, an optional precision argument was added
// to the following methods:
// 
// 1) String::String(const T& t, int precision)
// 2) Element::Element(const String& tagWord, const T& value, int precision)
// 3) Element::setValueAs(const T& value, int precision)
// 
// The following subtest verifies correct execution of these methods.  
void testOutputPrecision() {

    // Input value -----
    // 24 digits are used to initialize the variable 'input' below.
    // Since the maximum number of available digits for a double is 19
    // (see std::numeric_limits<long double>::max()), this initialization
    // will result in a consistent value that utilizes all available digits
    // whether the type of 'input' is float or double.
    // A Vec<1> is used so that the templatized String constructor is evoked.
    SimTK::Vec<1> input(0.123456789012345678901234);

    // Expected conversions to String for a range of precsions -----
    // The index of the 'expected' array below is the precision.
    // Precisions past a value of 7 are not currently tested so that tests
    // will execute properly when SimTK::Real is set to float or to double.
    // Though not tested past p = 7, the following rounding behavior is what
    // is observed for a double.
    const int n{24};
    SimTK::String expected[n];
    expected[0] = "~[0.1]";
    expected[1] = "~[0.1]";
    expected[2] = "~[0.12]";
    expected[3] = "~[0.123]";
    expected[4] = "~[0.1235]"; // Rounded up
    expected[5] = "~[0.12346]"; // Rounded up
    expected[6] = "~[0.123457]"; // Rounded up
    expected[7] = "~[0.1234568]"; // Rounded up
    expected[8] = "~[0.12345679]"; // Rounded up
    expected[9] = "~[0.123456789]";
    expected[10] = "~[0.123456789]";
    expected[11] = "~[0.12345678901]";
    expected[12] = "~[0.123456789012]";
    expected[13] = "~[0.1234567890123]";
    expected[14] = "~[0.12345678901235]"; // Rounded up
    expected[15] = "~[0.123456789012346]"; // Rounded up
    expected[16] = "~[0.1234567890123457]"; // Rounded up
    expected[17] = "~[0.12345678901234568]"; // Rounded up
    expected[18] = "~[0.123456789012345677]"; // Rounded up
    expected[19] = "~[0.1234567890123456774]";
    expected[20] = "~[0.12345678901234567737]";
    expected[21] = "~[0.12345678901234567737]"; // No long changes because
    expected[22] = "~[0.12345678901234567737]"; // precision is capped at
    expected[23] = "~[0.12345678901234567737]"; // LosslessNumDigitsReal

    // Store the precision of a local std::ostringstream object before any
    // calls to String::String. This is done to verify that changes made to
    // the precision inside String::String do not affect the precision of
    // other ostringstream instances.
    std::ostringstream osLocal;
    const auto starting_precision{osLocal.precision()};

    // ================
    // String::String()
    // ================
    int p;
    for(p=0; p <= 7; ++p) {
        // Specify the precision.
        String output(input, p);
        SimTK_TEST(output == expected[p]);

        // Verify that the precision of osLocal is unchanged.
        SimTK_TEST(osLocal.precision() == starting_precision);

        // Omit the precision argument, thereby testing the default argument.
        // Make sure that String::DefaultOutputPrecision < n so that we don't
        // step out of bounds on the 'expected' array.
        if (String::DefaultOutputPrecision < n) {
            String outputDefault(input);
            SimTK_TEST(outputDefault ==
                expected[String::DefaultOutputPrecision]);
        }
    }

    // Note that nothing bad happens when p is neg or 0.
    // String::String() does not check for p < 0 or p == 0; it just
    // relies on std::ostream to handle such values.
    // If the following tests fail, it may be because the implementation
    // of std::ostream varies across operating systems or has changed
    // since these tests were first put in place.
    // -----
    // p < 0: ostream does not accept and falls back on its starting precision.
    // 'starting_precision' (see above) is used to get the default output
    // value (see below) to handle a situation in which
    // String::DefaultOuputPrecision differs from ostream.precision()
    String outputDefault(input, (int)starting_precision);
    String outputNeg(input,-2);
    SimTK_TEST(outputNeg == outputDefault);
    // -----
    // p = 0: ostream takes the min p that applies (p = 1)
    String outputOne(input, 1);
    String outputZero(input, 0); 
    SimTK_TEST(outputZero == outputOne);

    //std::ostringstream os;
    //String outputDefault =
    //    stringStreamInsertHelper<Vec<1>>(os, input, true).str();

    // Test when pecision is greater than LosslessNumDigits.
    // In String::String(), p is capped at LosslessNumDigits.
    String outputLossless(input, SimTK::LosslessNumDigitsReal);
    String outputLosslessPlus1(input, SimTK::LosslessNumDigitsReal+1);
    SimTK_TEST(outputLosslessPlus1 == outputLossless);

    // =======================
    // Test Element::Element()
    // =======================
    // Element::Element() is just passing the arguments through to a
    // String::String() call. So, just going to verify a few notable cases.
    //------ use the default precision by not passing p
    Xml::Element defaultSigFigs("default", input);
    SimTK_TEST(defaultSigFigs.getValue() ==
        expected[String::DefaultOutputPrecision]);
    //------ a precision of 0 should get lower bounded to 1
    p = 0;
    Xml::Element zeroSigFigs("zero", input, p);
    SimTK_TEST(zeroSigFigs.getValue() == expected[1]);
    //------ will not round
    p = 3;
    Xml::Element threeSigFigs("three", input, p);
    SimTK_TEST(threeSigFigs.getValue() == expected[p]);
    //------ will round
    p = 7;
    Xml::Element sevenSigFigs("seven", input, p);
    SimTK_TEST(sevenSigFigs.getValue() == expected[p]);
    //------ above the max that makes sense
    p = SimTK::LosslessNumDigitsReal;
    Xml::Element losslessSigFigs("lossless", input, p);
    p = SimTK::LosslessNumDigitsReal + 1;
    Xml::Element losslessPlusOneSigFigs("losslessPlusOne", input, p);
    SimTK_TEST(losslessPlusOneSigFigs.getValue() ==
        losslessSigFigs.getValue());

    // ==========================
    // Test Element::setValueAs()
    // ==========================
    // Element::setValueAs() is just passing the arguments through to a
    // String::String() call. So, just going to verify a few notable cases.
    // These tests use the previously constructed element nodes but first set
    // all of their values to 0.0 as way of clearing the previous value.
    //------ first set the value of existing elements to 0.0
    Vec<1> zeroValue(0.0);
    defaultSigFigs.setValueAs<Vec<1>>(zeroValue);
    zeroSigFigs.setValueAs<Vec<1>>(zeroValue);
    threeSigFigs.setValueAs<Vec<1>>(zeroValue);
    sevenSigFigs.setValueAs<Vec<1>>(zeroValue);
    losslessSigFigs.setValueAs<Vec<1>>(zeroValue);
    losslessPlusOneSigFigs.setValueAs<Vec<1>>(zeroValue);
    SimTK_TEST(defaultSigFigs.getValue() == "~[0]")
    SimTK_TEST(zeroSigFigs.getValue() == "~[0]")
    SimTK_TEST(threeSigFigs.getValue() == "~[0]")
    SimTK_TEST(sevenSigFigs.getValue() == "~[0]")
    SimTK_TEST(losslessSigFigs.getValue() == "~[0]")
    SimTK_TEST(losslessPlusOneSigFigs.getValue() == "~[0]")
    //------ use the default precision by not passing p
    defaultSigFigs.setValueAs<Vec<1>>(input);
    SimTK_TEST(defaultSigFigs.getValue() ==
        expected[String::DefaultOutputPrecision]);
    //------ a precision of 0 should get lower bounded to 1
    p = 0;
    zeroSigFigs.setValueAs<Vec<1>>(input, p);
    SimTK_TEST(zeroSigFigs.getValue() == expected[1]);
    //------ will not round
    p = 3;
    threeSigFigs.setValueAs<Vec<1>>(input, p);
    SimTK_TEST(threeSigFigs.getValue() == expected[p]);
    //------ will round
    p = 7;
    sevenSigFigs.setValueAs<Vec<1>>(input, p);
    SimTK_TEST(sevenSigFigs.getValue() == expected[p]);
    //------ above the max that makes sense
    p = SimTK::LosslessNumDigitsReal;
    losslessSigFigs.setValueAs<Vec<1>>(input, p);
    p = SimTK::LosslessNumDigitsReal + 1;
    losslessPlusOneSigFigs.setValueAs<Vec<1>>(input, p);
    SimTK_TEST(losslessPlusOneSigFigs.getValue() ==
        losslessSigFigs.getValue());
}

int main() {
    cout << "Path of this executable: '" << Pathname::getThisExecutablePath() << "'\n";
    cout << "Executable directory: '" << Pathname::getThisExecutableDirectory() << "'\n";
    cout << "Current working directory: '" << Pathname::getCurrentWorkingDirectory() << "'\n";

    SimTK_START_TEST("TestXml");

        SimTK_SUBTEST(testOutputPrecision);
        SimTK_SUBTEST(testStringConvert);
        SimTK_SUBTEST(testXmlFromString);
        SimTK_SUBTEST(testXmlFromScratch);

    SimTK_END_TEST();
}

