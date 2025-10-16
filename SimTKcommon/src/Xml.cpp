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

#include "SimTKcommon/internal/Xml.h"

#include <cstddef>

#include "SimTKcommon/internal/String.h"
#include "tinyxml2.h"

using namespace SimTK;

// Handy helper for weeding out unwanted nodes.
static bool nodeTypeIsAllowed(Xml::NodeType allowed, tinyxml2::XMLNode* found) {
    if (!found) return false;

    if (found->ToElement())
        return (allowed & Xml::ElementNode) != 0;
    else if (found->ToText())
        return (allowed & Xml::TextNode) != 0;
    else if (found->ToComment())
        return (allowed & Xml::CommentNode) != 0;
    else if (found->ToUnknown())
        return (allowed & Xml::UnknownNode) != 0;

    return false;
}

// This is an Xml namespace-scope free function.
String Xml::getNodeTypeAsString(Xml::NodeType type) {
    // Take care of special cases first.
    switch (type) {
        case NoNode:
            return "NoNode";
        case NoJunkNodes:
            return "NoJunkNodes";
        case JunkNodes:
            return "JunkNodes";
        case AnyNodes:
            return "AnyNodes";
        default:;  // fall through
    }

    // "Or" bits together in CUTE order.
    String out;
    if (type & CommentNode) out = "CommentNode";
    if (type & UnknownNode) {
        if (!out.empty()) out += "|";
        out += "UnknownNode";
    }
    if (type & TextNode) {
        if (!out.empty()) out += "|";
        out += "TextNode";
    }
    if (type & ElementNode) {
        if (!out.empty()) out += "|";
        out += "ElementNode";
    }

    return out;
}

//------------------------------------------------------------------------------
//                        XML :: DOCUMENT :: IMPL
//------------------------------------------------------------------------------
class Xml::Document::Impl {
   public:
    Impl() {}  // not canonicalized yet
    Impl(const String& pathname) { readFromFile(pathname.c_str()); }
    ~Impl() { clear(); }

    // Note that the copy must be canonicalized before use -- that's so we
    // get our root element pointing correctly into the copy rather than
    // at the source's root element.
    Impl* clone() const {
        // Impl* newImpl = new Impl();
        // newImpl->m_tixml = m_tixml;
        // root element isn't set yet
        return nullptr;
    }

    void clear() {
        m_tixml.Clear();
        m_rootElement.clear();
    }

    void setIndentString(const String& indent) {
        // m_tixml.setIndentString(indent);
    }
    const String& getIndentString() const {
        //  return m_tixml.GetIndentString();
    }

    void readFromFile(const String& pathname) {
        clear();
        m_tixml.SetValue(pathname);
        auto status = m_tixml.LoadFile(pathname);
        SimTK_ERRCHK2_ALWAYS(
            status != tinyxml2::XMLError::XML_SUCCESS, "Xml::readFromFile()",
            "Failed to load the Xml file '%s' with error '%s'.",
            pathname.c_str(), m_tixml.ErrorStr());
    }

    void writeToFile(const String& pathname) const {
        tinyxml2::XMLDocument tempDoc;
        m_tixml.DeepCopy(&tempDoc);
        auto status = tempDoc.SaveFile(pathname.c_str());
        SimTK_ERRCHK2_ALWAYS(
            status != tinyxml2::XMLError::XML_SUCCESS, "Xml::writeToFile()",
            "Failed to write to the Xml file '%s' with error '%s'.",
            pathname.c_str(), m_tixml.ErrorStr());
    }

    void readFromString(const char* xmlDocument) {
        clear();
        auto status = m_tixml.Parse(xmlDocument);
        SimTK_ERRCHK1_ALWAYS(status != tinyxml2::XMLError::XML_SUCCESS, "Xml::readFromString()",
                             "Failed to parse the Xml string with error '%s'.", m_tixml.ErrorStr());
    }

    void writeToString(String& xmlDocument, bool compact) const {
        tinyxml2::XMLPrinter printer;

        // if (compact)
        //     printer.SetStreamPrinting();
        // else
        //     printer.SetIndent(m_tixml.GetIndentChars());
        m_tixml.Print(&printer);
        xmlDocument = printer.CStr();
    }

    // Call this during construction and after a new Xml document has been
    // parsed. It guarantees: (1) there is a Declaration record and it is
    // the first node at the top level of the Xml document, and (2) there
    // is no top-level text and only one top-level element and that is the
    // "root" element whose tag name is the document type.
    void canonicalizeDocument() {
        tinyxml2::XMLDeclaration* decl = addDeclarationIfNeeded();
        tinyxml2::XMLElement* root = addRootElementIfNeeded();

        m_rootElement.setTiNodePtr(root);
    }

    const ::tinyxml2::XMLDeclaration& XMLDeclaration() const {
        const tinyxml2::XMLNode* decl = m_tixml.FirstChild();
        assert(decl && decl->Type() == tinyxml2::XMLNode::DECLARATION);
        return *decl->ToDeclaration();
    }

    tinyxml2::XMLDeclaration& XMLDeclaration() {
        tinyxml2::XMLNode* decl = m_tixml.FirstChild();
        assert(decl && decl->Type() == tinyxml2::XMLNode::DECLARATION);
        return *decl->ToDeclaration();
    }

    // The first thing in every Xml file should be the Declaration, that is,
    // the line that looks something like
    //      <?xml version="1.0" encoding="UTF-8"?>
    // That's because the encoding is needed to read the rest of the file.
    // We'll look through the file to see if there is a declaration at
    // the top level. If so, we'll move it to the first slot if necessary.
    // If not, we'll create a default one and put it first. Either way we
    // return a pointer to the declaration which will now be the first
    // node in the Xml document.
    ::tinyxml2::XMLDeclaration* addDeclarationIfNeeded() {
        tinyxml2::XMLNode* child = m_tixml.FirstChild();
        if (child && child->ToDeclaration())
            return child->ToDeclaration();  // the easy and most common case

        // Otherwise hunt for the declaration.
        for (; child; child = child->NextSibling())
            if (child->ToDeclaration()) break;

        // if (child)
        //     m_tixml.DeleteChild(child);  // it's in the wrong place
        // else
        //     child = new tinyxml2::XMLDeclaration("1.0", "UTF-8", "");

        // Insert new node as the first in the document.
        m_tixml.InsertFirstChild(child);
        return child->ToDeclaration();
    }

    // If the supplied Xml document has zero or more than one top-level element,
    // or has top-level text, we'll add <_Root> elements </_Root> to
    // encapsulate all the top-level text and element nodes (this may also
    // surround top-level comments and unknowns if they occur between the text
    // and elements. A pointer to the root element (whether original or added)
    // is returned.
    tinyxml2::XMLElement* addRootElementIfNeeded() {
        // Find the first element or text node, and remember the node that
        // preceded it since that's where we may be inserting the new root
        // element.
        tinyxml2::XMLNode* nodeBeforeFirst = 0;
        tinyxml2::XMLNode* firstEltOrText = m_tixml.FirstChild();
        while (firstEltOrText &&
               !(firstEltOrText->ToElement() || firstEltOrText->ToText())) {
            nodeBeforeFirst = firstEltOrText;
            firstEltOrText = firstEltOrText->NextSibling();
        }

        if (!firstEltOrText) {
            // No top level element or text node. We'll just append an empty
            // _Root element to the end of whatever's there.
            // tinyxml2::XMLElement* root = new tinyxml2::XMLElement("_Root");
            // m_tixml.LinkEndChild(root);
            // return root;
        }

        // There is at least one element or text node at the top level; the
        // first one is pointed to by firstEltOrText, and the node just before
        // that (if any) is pointed to by nodeBeforeFirst. Now find the last
        // top-level element or text node.
        tinyxml2::XMLNode* lastEltOrText = m_tixml.LastChild();
        while (lastEltOrText &&
               !(lastEltOrText->ToElement() || lastEltOrText->ToText()))
            lastEltOrText = lastEltOrText->PreviousSibling();

        assert(lastEltOrText);  // should have at least re-found the first one!

        // If the extremely likely case that the first and last are the same
        // node and that node is an element, then the document is already in the
        // right format and we don't have to do anything to it.
        if (firstEltOrText == lastEltOrText && firstEltOrText->ToElement())
            return firstEltOrText->ToElement();

        // Now we know there is top-level text or more than one top-level
        // element so we are going to have to surround everything between
        // first and last with a new root element.
        tinyxml2::XMLElement* root = m_tixml.NewElement("_Root");

        tinyxml2::XMLNode* nextToMove = firstEltOrText;
        // while (true) {
        //     assert(nextToMove);  // can't happen!
        //     tinyxml2::XMLNode* moveMe = nextToMove;
        //     nextToMove = moveMe->NextSibling();
        //     root->LinkEndChild(m_tixml.Child(moveMe));
        //     // Did we just move the last element or text node?
        //     if (moveMe == lastEltOrText) break;
        // }

        // Now link the new root element right where we found the first
        // element or text node.
        // if (nodeBeforeFirst)
        m_tixml.LinkEndChild(nodeBeforeFirst);
        // else
            // m_tixml.LinkEndChild(root);

        return root;
    }

    // The Tiny XML "document" is a tinyxml2::XMLNode representing the entire
    // XML object. We instead represent the document as an object of class Xml
    // which is not itself a Node.
    ::tinyxml2::XMLDocument m_tixml;

    // This references the root element withing the TinyXml document. It
    // is filled in when the document is initially canonicalized.
    Xml::Element m_rootElement;

   private:
    Impl(const Impl&);             // disable; use clone()
    Impl& operator=(const Impl&);  // "
};

//------------------------------------------------------------------------------
//                            XML :: DOCUMENT
//------------------------------------------------------------------------------

Xml::Document::Document() : impl(0) {
    impl = new Impl();
    impl->canonicalizeDocument();
}

Xml::Document::Document(const String& pathname) : impl(0) {
    impl = new Impl(pathname);
    impl->canonicalizeDocument();
}

Xml::Document::Document(const Document& source) : impl(0) {
    if (source.impl) {
        impl = source.impl->clone();
        impl->canonicalizeDocument();
    }
}

Xml::Document& Xml::Document::operator=(const Xml::Document& source) {
    if (&source != this) {
        delete impl;
        impl = 0;
        if (source.impl) {
            impl = source.impl->clone();
            impl->canonicalizeDocument();
        }
    }
    return *this;
}

Xml::Document::~Document() {
    delete impl;
    impl = 0;
}

void Xml::Document::clear() {
    updImpl().clear();
}

void Xml::Document::readFromFile(const String& pathname) {
    updImpl().readFromFile(pathname.c_str());
    updImpl().canonicalizeDocument();
}
void Xml::Document::writeToFile(const String& pathname) const {
    getImpl().writeToFile(pathname.c_str());
}
void Xml::Document::readFromString(const char* xmlDocument) {
    updImpl().readFromString(xmlDocument);
    updImpl().canonicalizeDocument();
}
void Xml::Document::readFromString(const String& xmlDocument) {
    readFromString(xmlDocument.c_str());
}
void Xml::Document::writeToString(String& xmlDocument, bool compact) const {
    getImpl().writeToString(xmlDocument, compact);
}
void Xml::Document::setIndentString(const String& indent) {
    updImpl().setIndentString(indent);
}
const String& Xml::Document::getIndentString() const {
    return getImpl().getIndentString();
}

/*static*/ void Xml::Document::setXmlCondenseWhiteSpace(bool shouldCondense) {
    // tinyxml2::XMLBase::SetCondenseWhiteSpace(shouldCondense);
}
/*static*/ bool Xml::Document::isXmlWhiteSpaceCondensed() {
    // return tinyxml2::XMLBase::IsWhiteSpaceCondensed();
}

Xml::Element Xml::Document::getRootElement() {
    assert(getImpl().m_rootElement.isValid());
    return updImpl().m_rootElement;
}
const String& Xml::Document::getRootTag() const {
    return unconst().getRootElement().getElementTag();
}
void Xml::Document::setRootTag(const String& tag) {
    getRootElement().setElementTag(tag);
}

void Xml::Document::insertTopLevelNodeAfter(const Xml::node_iterator& afterThis,
                                            Xml::Node insertThis) {
    const char* method = "Xml::insertTopLevelNodeAfter()";

    // Check that the supplied Node is OK.
    SimTK_ERRCHK_ALWAYS(insertThis.isValid(), method,
                        "The supplied Node handle was empty.");
    SimTK_ERRCHK_ALWAYS(insertThis.isOrphan(), method,
                        "The Node was not an orphan so can't be inserted.");
    SimTK_ERRCHK1_ALWAYS(
        Comment::isA(insertThis) || Unknown::isA(insertThis), method,
        "The Node had NodeType %s, but only Comment and Unknown nodes"
        " can be inserted at the topmost document level.",
        insertThis.getNodeTypeAsString().c_str());

    // If no iterator, add the Node to the end and we're done.
    if (afterThis == node_end()) {
        updImpl().m_tixml.LinkEndChild(insertThis.updTiNodePtr());
        return;
    }

    // There is an iterator, make sure it's a top-level one.
    SimTK_ERRCHK_ALWAYS(afterThis->isTopLevelNode(), method,
                        "The node_iterator did not refer to a top-level Node.");

    // updImpl().m_tixml.LinkAfterChild(afterThis->updTiNodePtr(),
    //                                  insertThis.updTiNodePtr());
}

void Xml::Document::insertTopLevelNodeBefore(
    const Xml::node_iterator& beforeThis, Xml::Node insertThis) {
    const char* method = "Xml::insertTopLevelNodeBefore()";

    // Check that the supplied Node is OK.
    SimTK_ERRCHK_ALWAYS(insertThis.isValid(), method,
                        "The supplied Node handle was empty.");
    SimTK_ERRCHK_ALWAYS(insertThis.isOrphan(), method,
                        "The Node was not an orphan so can't be inserted.");
    SimTK_ERRCHK1_ALWAYS(
        Comment::isA(insertThis) || Unknown::isA(insertThis), method,
        "The Node had NodeType %s, but only Comment and Unknown nodes"
        " can be inserted at the topmost document level.",
        insertThis.getNodeTypeAsString().c_str());

    // If no iterator, add the Node to the end and we're done.
    if (beforeThis == node_end()) {
        updImpl().m_tixml.LinkEndChild(insertThis.updTiNodePtr());
        return;
    }

    // There is an iterator, make sure it's a top-level one.
    SimTK_ERRCHK_ALWAYS(beforeThis->isTopLevelNode(), method,
                        "The node_iterator did not refer to a top-level Node.");

    // updImpl().m_tixml.LinkBeforeChild(beforeThis->updTiNodePtr(),
    //                                   insertThis.updTiNodePtr());
}

void Xml::Document::eraseTopLevelNode(const Xml::node_iterator& deleteThis) {
    const char* method = "Xml::eraseTopLevelNode()";

    // Check that the supplied iterator points to something.
    SimTK_ERRCHK_ALWAYS(
        deleteThis != node_end(), method,
        "The node_iterator is at node_end() so doesn't refer to a Node.");
    // There is an iterator, make sure it's a top-level one.
    SimTK_ERRCHK_ALWAYS(deleteThis->isTopLevelNode(), method,
                        "The node_iterator did not refer to a top-level Node.");
    SimTK_ERRCHK1_ALWAYS(
        Comment::isA(*deleteThis) || Unknown::isA(*deleteThis), method,
        "The Node had NodeType %s, but only Comment and Unknown nodes"
        " can be erased at the topmost document level.",
        deleteThis->getNodeTypeAsString().c_str());

    // updImpl().m_tixml.RemoveChild(deleteThis->updTiNodePtr());
}

Xml::Node Xml::Document::removeTopLevelNode(
    const Xml::node_iterator& removeThis) {
    const char* method = "Xml::removeTopLevelNode()";

    // Check that the supplied iterator points to something.
    SimTK_ERRCHK_ALWAYS(
        removeThis != node_end(), method,
        "The node_iterator is at node_end() so doesn't refer to a Node.");
    // There is an iterator, make sure it's a top-level one.
    SimTK_ERRCHK_ALWAYS(removeThis->isTopLevelNode(), method,
                        "The node_iterator did not refer to a top-level Node.");
    SimTK_ERRCHK1_ALWAYS(
        Comment::isA(*removeThis) || Unknown::isA(*removeThis), method,
        "The Node had NodeType %s, but only Comment and Unknown nodes"
        " can be removed at the topmost document level.",
        removeThis->getNodeTypeAsString().c_str());

    // tinyxml2::XMLNode* p =
    //     updImpl().m_tixml.DisconnectChild(removeThis->updTiNodePtr());
    // return Node(p);
    return Node();
}

String Xml::Document::getXmlVersion() const {
    // return getImpl().gettinyxml2::XMLDeclaration().Version();
}
String Xml::Document::getXmlEncoding() const {
    // return getImpl().gettinyxml2::XMLDeclaration().Encoding();
}
bool Xml::Document::getXmlIsStandalone() const {
    // return String::toLower(
    //            getImpl().gettinyxml2::XMLDeclaration().Standalone()) != "no";
}

void Xml::Document::setXmlVersion(const String& version) {
    // updImpl().updtinyxml2::XMLDeclaration().SetVersion(version.c_str());
}
void Xml::Document::setXmlEncoding(const String& encoding) {
    // updImpl().updtinyxml2::XMLDeclaration().SetEncoding(encoding.c_str());
}
// For the standalone case we set the string to "" rather than "yes" so that
// the standalone attribute won't appear in the output declaration.
void Xml::Document::setXmlIsStandalone(bool isStandalone) {
    // updImpl().updtinyxml2::XMLDeclaration().SetStandalone(isStandalone ? ""
    //                                                                    :
    //                                                                    "no");
}

// XML node_begin()
Xml::node_iterator Xml::Document::node_begin(NodeType allowed) {
    tinyxml2::XMLNode* first = updImpl().m_tixml.FirstChild();
    // while (first && !nodeTypeIsAllowed(allowed, first->Type()))
    //     first = first->NextSibling();
    // return node_iterator(first, allowed);
}

// XML node_end()
Xml::node_iterator Xml::Document::node_end() const {
    return node_iterator(0);
}

//------------------------------------------------------------------------------
//                              XML ATTRIBUTE
//------------------------------------------------------------------------------
// Xml::Attribute::Attribute(const String& name, const String& value)
//     : tiAttr(new tinyxml2::XMLAttribute(name, value)) {}

void Xml::Attribute::clear() {
    tiAttr = 0;
}

void Xml::Attribute::clearOrphan() {
    if (!tiAttr) return;
    SimTK_ERRCHK_ALWAYS(
        isOrphan(), "Xml::Attribute::clearOrphan()",
        "This Attribute is not an orphan (or it was already destructed"
        " and now contains garbage).");
    // delete tiAttr;  // not part of any document
    tiAttr = 0;
}

// Note that the criteria for orphanhood is that the referenced
// tinyxml2::XMLAttr is not in a document.
bool Xml::Attribute::isOrphan() const {
    if (!isValid()) return false;  // empty handle not considered an orphan
    // return tiAttr->GetDocument() == 0;
    return false;
}

const String& Xml::Attribute::getName() const {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::getName()",
                        "The attribute handle was empty.");
    // return getTiAttr().NameStr();
    String string = nullptr;
    return string;
}
const String& Xml::Attribute::getValue() const {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::getValue()",
                        "The attribute handle was empty.");
    // return getTiAttr().ValueStr();
    String string = nullptr;
    return string;
}

Xml::Attribute& Xml::Attribute::setName(const String& name) {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::setName()",
                        "The attribute handle was empty.");
    // updTiAttr().SetName(name);
    return *this;
}
Xml::Attribute& Xml::Attribute::setValue(const String& value) {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::setValue()",
                        "The attribute handle was empty.");
    // updTiAttr().SetValue(value);
    return *this;
}

void Xml::Attribute::writeToString(String& out) const {
    // getTiAttr().Print(0, 0, &out);
}

//------------------------------------------------------------------------------
//                         XML ATTRIBUTE ITERATOR
//------------------------------------------------------------------------------
Xml::attribute_iterator& Xml::attribute_iterator::operator++() {
    // tinyxml2::XMLAttribute* next = attr.updTiAttr().Next();
    // attr.setTiAttrPtr(next);
    return *this;
}

Xml::attribute_iterator Xml::attribute_iterator::operator++(int) {
    Attribute save(attr);
    // tinyxml2::XMLAttribute* next = attr.updTiAttr().Next();
    // attr.setTiAttrPtr(next);
    return attribute_iterator(save);
}

Xml::attribute_iterator& Xml::attribute_iterator::operator--() {
    // tinyxml2::XMLAttribute* prev = attr.updTiAttr().Previous();
    // attr.setTiAttrPtr(prev);
    return *this;
}

Xml::attribute_iterator Xml::attribute_iterator::operator--(int) {
    Attribute save(attr);
    // tinyxml2::XMLAttribute* prev = attr.updTiAttr().Previous();
    // attr.setTiAttrPtr(prev);
    return attribute_iterator(save);
}

//------------------------------------------------------------------------------
//                                 XML NODE
//------------------------------------------------------------------------------
Xml::Node Xml::Node::clone() const {
    tinyxml2::XMLNode* newNode = 0;
    // if (tiNode) newNode = tiNode->Clone();
    // return Node(newNode);
}

void Xml::Node::clear() {
    // TODO: Note that we do not clean up heap space here if this node is
    // still an orphan. To do that requires that we reference count the tiNode
    // so that we don't risk looking at deleted garbage in trying to determine
    // orphanhood.
    tiNode = 0;
}

void Xml::Node::clearOrphan() {
    if (tiNode == 0) return;
    SimTK_ERRCHK_ALWAYS(
        isOrphan(), "Xml::Node::clearOrphan()",
        "This Node is not an orphan (or it was already destructed and now"
        " contains garbage).");
    // delete tiNode;
    tiNode = 0;
}

Xml::NodeType Xml::Node::getNodeType() const {
    if (!isValid()) return NoNode;
    // switch (getTiNode().Type()) {
    //     case tinyxml2::XMLNode::COMMENT:
    //         return CommentNode;
    //     case tinyxml2::XMLNode::UNKNOWN:
    //         return UnknownNode;
    //     case tinyxml2::XMLNode::TEXT:
    //         return TextNode;
    //     case tinyxml2::XMLNode::ELEMENT:
    //         return ElementNode;
    //     default:
    //         SimTK_ASSERT1_ALWAYS(false,
    //                              "Xml::Node::getNodeType(): can't convert "
    //                              "TinyXML node type %s to any"
    //                              " SimTK::Xml node type.",
    //                              getTiNode().TypeName());
    // }
    return NoNode;
}

// Invoke the namespace-scope utility on this node's type.
String Xml::Node::getNodeTypeAsString() const {
    return Xml::getNodeTypeAsString(getNodeType());
}

const String& Xml::Node::getNodeText() const {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Node::getText()",
                        "Can't get text from an empty Node handle.");

    // return getTiNode().ValueStr();
    String string = nullptr;
    return string;
}

bool Xml::Node::isTopLevelNode() const {
    if (!isValid()) return false;
    const tinyxml2::XMLNode& n = getTiNode();
    return n.Parent() && n.Parent() == n.GetDocument();
}

// Note that the criteria for orphanhood is that the *tinyxml2::XMLNode* has no
// parent. In SimTK::Xml, none of the top-level nodes are considered to have
// a parent since the Xml document is not a node there, while in TinyXML
// the tinyxml2::XMLDocument is a tinyxml2::XMLNode, so even top-level nodes
// have a parent.
bool Xml::Node::isOrphan() const {
    if (!isValid()) return false;  // empty handle not considered an orphan
    return getTiNode().Parent() == 0;
}

void Xml::Node::writeToString(String& out, bool compact) const {
    if (!isValid()) {
        out = "<!-- EMPTY NODE -->";
        return;
    }
    // tinyxml2::XMLPrinter printer(out);
    // if (compact) printer.SetStreamPrinting();
    // getTiNode().Accept(&printer);
}

bool Xml::Node::hasParentElement() const {
    if (!isValid()) return false;
    const tinyxml2::XMLNode* parent = getTiNode().Parent();
    return parent && parent->ToElement();
}

Xml::Element Xml::Node::getParentElement() {
    SimTK_ERRCHK_ALWAYS(
        hasParentElement(), "Xml::Node::hasParentElement()",
        "This node does not have a parent element; it may be a top-level"
        " node, an orphan, or just an empty Node handle.");
    return Element(updTiNode().Parent()->ToElement());
}

//------------------------------------------------------------------------------
//                          XML NODE ITERATOR
//------------------------------------------------------------------------------
Xml::node_iterator& Xml::node_iterator::operator++() {
    tinyxml2::XMLNode* next = node.updTiNode().NextSibling();
    // while (next && !nodeTypeIsAllowed(allowed, next->Type()))
    //     next = next->NextSibling();
    // node = Node(next);
    return *this;
}

Xml::node_iterator Xml::node_iterator::operator++(int) {
    Node save(node);
    tinyxml2::XMLNode* next = node.updTiNode().NextSibling();
    // while (next && !nodeTypeIsAllowed(allowed, next->Type()))
    //     next = next->NextSibling();
    node = Node(next);
    return node_iterator(save);
}

Xml::node_iterator& Xml::node_iterator::operator--() {
    tinyxml2::XMLNode* prev = node.updTiNode().PreviousSibling();
    // while (prev && !nodeTypeIsAllowed(allowed, prev->Type()))
    //     prev = prev->PreviousSibling();
    node = Node(prev);
    return *this;
}

Xml::node_iterator Xml::node_iterator::operator--(int) {
    Node save(node);
    tinyxml2::XMLNode* prev = node.updTiNode().PreviousSibling();
    // while (prev && !nodeTypeIsAllowed(allowed, prev->Type()))
    //     prev = prev->PreviousSibling();
    node = Node(prev);
    return node_iterator(save);
}

//------------------------------------------------------------------------------
//                              XML ELEMENT
//------------------------------------------------------------------------------

// Handy helper for weeding out unwanted nodes and elements.
static bool elementIsAllowed(const String& tag,
                             const tinyxml2::XMLElement* elt) {
    if (elt == 0) return false;
    // return tag.empty() || elt->ValueStr() == tag;
    return false;
}

// Xml::Element::Element(const String& tag, const String& value)
//     : Node(new tinyxml2::XMLElement(tag)) {
//     if (value.empty()) return;
//     // We need to add a Text node.
//     updTiElement().LinkEndChild(new tinyxml2::XMLText(value));
// }

Xml::Element Xml::Element::clone() const {
    tinyxml2::XMLElement* newElt = 0;
    // if (getTiElementPtr()) newElt = getTiElementPtr()->Clone()->ToElement();
    return Element(newElt);
}

const String& Xml::Element::getElementTag() const {
    // return getTiNode().ValueStr();
}
void Xml::Element::setElementTag(const String& type) {
    updTiNode().SetValue(type);
}

bool Xml::Element::isValueElement() const {
    if (unconst().element_begin() != element_end())
        return false;  // has child elements
    node_iterator text = unconst().node_begin(TextNode);
    return text == node_end() || ++text == node_end();  // zero or one
}

const String& Xml::Element::getValue() const {
    static const String null;
    SimTK_ERRCHK1_ALWAYS(isValueElement(), "Xml::Element::getValue()",
                         "Element <%s> is not a value element.",
                         getElementTag().c_str());

    node_iterator text = unconst().node_begin(TextNode);
    return text == node_end() ? null : text->getNodeText();
}

// Must add a Text node now if this Element doesn't have one.
String& Xml::Element::updValue() {
    SimTK_ERRCHK1_ALWAYS(isValueElement(), "Xml::Element::getValue()",
                         "Element <%s> is not a value element.",
                         getElementTag().c_str());

    node_iterator text = node_begin(TextNode);
    // if (text != node_end()) return Text::getAs(*text).updText();

    // We need to add a Text node.
    // tinyxml2::XMLText* textp = new tinyxml2::XMLText("");
    // updTiElement().LinkEndChild(textp);
    // return textp->UpdValueStr();
    String string = nullptr;
    return string;
}

// If there is no Text node we'll add one; if there is just one we'll
// change its value; otherwise report an error.
void Xml::Element::setValue(const String& value) {
    SimTK_ERRCHK1_ALWAYS(isValueElement(), "Xml::Element::setValue()",
                         "Element <%s> is not a value element.",
                         getElementTag().c_str());

    node_iterator text = node_begin(TextNode);
    // if (text == node_end())
    // updTiNode().LinkEndChild(new tinyxml2::XMLText(value));
    // else
    text->updTiNode().SetValue(value);
}

bool Xml::Element::hasAttribute(const String& name) const {
    return unconst().getOptionalAttribute(name).isValid();
}

bool Xml::Element::hasElement(const String& tag) const {
    return unconst().element_begin(tag) != element_end();
}

bool Xml::Element::hasNode(NodeType allowed) const {
    return unconst().node_begin(allowed) != node_end();
}

Xml::Element Xml::Element::getRequiredElement(const String& tag) {
    element_iterator p = element_begin(tag);
    SimTK_ERRCHK2_ALWAYS(
        p != element_end(), "Xml::Element::getRequiredElement()",
        "Couldn't find required child element <%s> in element <%s>.",
        tag.c_str(), getElementTag().c_str());
    return *p;
}

Xml::Element Xml::Element::getOptionalElement(const String& tag) {
    element_iterator p = element_begin(tag);
    return p != element_end() ? *p : Element(0);
}

Xml::Attribute Xml::Element::getRequiredAttribute(const String& name) {
    Attribute attr = getOptionalAttribute(name);
    SimTK_ERRCHK2_ALWAYS(attr.isValid(), "Xml::Element::getRequiredAttribute()",
                         "Couldn't find required attribute %s in element <%s>.",
                         name.c_str(), getElementTag().c_str());
    return attr;
}

Xml::Attribute Xml::Element::getOptionalAttribute(const String& name) {
    for (attribute_iterator p = attribute_begin(); p != attribute_end(); ++p)
        if (p->getName() == name) return *p;
    return Attribute(0);
}

void Xml::Element::setAttributeValue(const String& name, const String& value) {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Element::setAttributeValue()",
                        "Can't add an attribute to an empty Element handle.");
    updTiElement().SetAttribute(name, value);
}

void Xml::Element::eraseAttribute(const String& name) {
    SimTK_ERRCHK_ALWAYS(
        isValid(), "Xml::Element::eraseAttribute()",
        "Can't erase an attribute from an empty Element handle.");
    // updTiElement().RemoveAttribute(name);
}

/*static*/ bool Xml::Element::isA(const Xml::Node& node) {
    if (!node.isValid()) return false;
    return node.getTiNode().ToElement() != 0;
}
/*static*/ const Xml::Element& Xml::Element::getAs(const Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Element::getAs()",
        "The given Node was not an Element; it is a %s. Use Element::isA()"
        " to check before calling Element::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Element&>(node);
}
/*static*/ Xml::Element& Xml::Element::getAs(Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Element::getAs()",
        "The given Node was not an Element; it is a %s. Use Element::isA()"
        " to check before calling Element::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Element&>(node);
}

void Xml::Element::insertNodeBefore(const node_iterator& beforeThis,
                                    Node node) {
    const char* method = "Xml::Element::insertNodeBefore()";
    const char* tag = getElementTag().c_str();

    SimTK_ERRCHK1_ALWAYS(
        node.isValid(), method,
        "The supplied Node handle was invalid so can't be inserted into"
        " Element <%s>.",
        tag);
    SimTK_ERRCHK1_ALWAYS(
        !node.hasParentElement(), method,
        "The supplied Node already had a parent so can't be inserted into"
        " Element <%s>.",
        tag);

    if (beforeThis == node_end()) {
        updTiNode().LinkEndChild(node.updTiNodePtr());
        return;
    }

    SimTK_ERRCHK1_ALWAYS(
        beforeThis->getParentElement() == *this, method,
        "The supplied node_iterator referred to a node that was not a"
        "child node of this Element <%s>.",
        tag);

    tinyxml2::XMLNode* p = beforeThis->updTiNodePtr();
    // updTiNode().LinkBeforeChild(p, node.updTiNodePtr());
}

void Xml::Element::insertNodeAfter(const node_iterator& afterThis, Node node) {
    const char* method = "Xml::Element::insertNodeAfter()";
    const char* tag = getElementTag().c_str();

    SimTK_ERRCHK1_ALWAYS(
        node.isValid(), method,
        "The supplied Node handle was invalid so can't be inserted into"
        " Element <%s>.",
        tag);
    SimTK_ERRCHK1_ALWAYS(
        !node.hasParentElement(), method,
        "The supplied Node already had a parent so can't be inserted into"
        " Element <%s>.",
        tag);

    if (afterThis == node_end()) {
        updTiNode().LinkEndChild(node.updTiNodePtr());
        return;
    }

    SimTK_ERRCHK1_ALWAYS(
        afterThis->getParentElement() == *this, method,
        "The supplied node_iterator referred to a node that was not a"
        "child node of this Element <%s>.",
        tag);

    tinyxml2::XMLNode* p = afterThis->updTiNodePtr();
    // updTiNode().LinkEndChild(p, node.updTiNodePtr());
}

void Xml::Element::eraseNode(const Xml::node_iterator& deleteThis) {
    const char* method = "Xml::Element::eraseNode()";

    // Check that the supplied iterator points to something.
    SimTK_ERRCHK_ALWAYS(
        deleteThis != node_end(), method,
        "The node_iterator is at node_end() so doesn't refer to a Node.");
    // There is an iterator, make sure it points to a child of this element.
    SimTK_ERRCHK_ALWAYS(
        deleteThis->hasParentElement() &&
            deleteThis->getParentElement() == *this,
        method, "The node_iterator did not refer to a child of this Element.");

    // updTiElement().RemoveChild(deleteThis->updTiNodePtr());
}

Xml::Node Xml::Element::removeNode(const Xml::node_iterator& removeThis) {
    const char* method = "Xml::Element::removeNode()";

    // Check that the supplied iterator points to something.
    SimTK_ERRCHK_ALWAYS(
        removeThis != node_end(), method,
        "The node_iterator is at node_end() so doesn't refer to a Node.");
    // There is an iterator, make sure it points to a child of this element.
    SimTK_ERRCHK_ALWAYS(
        removeThis->hasParentElement() &&
            removeThis->getParentElement() == *this,
        method, "The node_iterator did not refer to a child of this Element.");

    // tinyxml2::XMLNode* p =
    //     updTiElement().DisconnectChild(removeThis->updTiNodePtr());
    // return Xml::Node(p);
    return Xml::Node();
}

// Element node_begin()
Xml::node_iterator Xml::Element::node_begin(NodeType allowed) {
    tinyxml2::XMLNode* first = updTiNode().FirstChild();
    // while (first && !nodeTypeIsAllowed(allowed, first->Type()))
    //     first = first->NextSibling();
    return node_iterator(first, allowed);
}

// Element node_end()
Xml::node_iterator Xml::Element::node_end() const {
    return node_iterator(0);
}

// Element begin()
Xml::element_iterator Xml::Element::element_begin(const String& tag) {
    tinyxml2::XMLElement* first = updTiNode().FirstChildElement();
    while (first && !elementIsAllowed(tag, first))
        first = first->NextSiblingElement();
    return element_iterator(first, tag);
}

// Element end()
Xml::element_iterator Xml::Element::element_end() const {
    return element_iterator(0);
}

// Attribute begin()
Xml::attribute_iterator Xml::Element::attribute_begin() {
    const tinyxml2::XMLAttribute* first = getTiElement().FirstAttribute();
    // return attribute_iterator(first);
    return attribute_iterator();
}

// Attribute end()
Xml::attribute_iterator Xml::Element::attribute_end() const {
    return attribute_iterator(0);
}

//------------------------------------------------------------------------------
//                          XML ELEMENT ITERATOR
//------------------------------------------------------------------------------
Xml::element_iterator& Xml::element_iterator::operator++() {
    tinyxml2::XMLElement* next = (*this)->updTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag, next))
        next = next->NextSiblingElement();
    reassign(next);
    return *this;
}

Xml::element_iterator Xml::element_iterator::operator++(int) {
    Element save(*(*this));
    tinyxml2::XMLElement* next = (*this)->updTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag, next))
        next = next->NextSiblingElement();
    reassign(next);
    return element_iterator(save);
}

Xml::element_iterator& Xml::element_iterator::operator--() {
    tinyxml2::XMLElement* prev =
        (*this)->updTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag, prev))
        prev = prev->PreviousSiblingElement();
    reassign(prev);
    return *this;
}

Xml::element_iterator Xml::element_iterator::operator--(int) {
    Element save(*(*this));
    tinyxml2::XMLElement* prev =
        (*this)->updTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag, prev))
        prev = prev->PreviousSiblingElement();
    reassign(prev);
    return element_iterator(save);
}

//------------------------------------------------------------------------------
//                             XML TEXT NODE
//------------------------------------------------------------------------------
Xml::Text::Text(const String& text) : Node(new tinyxml2::XMLDocument(text)) {}

Xml::Text Xml::Text::clone() const {
    tinyxml2::XMLDocument* newText = 0;
    if (getTiNodePtr()) getTiNodePtr()->ShallowClone(newText)->ToText();
    // return newText->ToText();
    return Xml::Text();
}

const SimTK::String Xml::Text::getText() const {
    return SimTK::String(getTiNode().Value());
}

void Xml::Text::setText(const SimTK::String& newText) {
    updTiNode().SetValue(newText.c_str());
}

/*static*/ bool Xml::Text::isA(const Xml::Node& node) {
    if (!node.isValid()) return false;
    return node.getTiNode().ToText() != 0;
}
/*static*/ const Xml::Text& Xml::Text::getAs(const Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Text::getAs()",
        "The given Node was not a Text node; it is a %s. Use Text::isA()"
        " to check before calling Text::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Text&>(node);
}
/*static*/ Xml::Text& Xml::Text::getAs(Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Text::getAs()",
        "The given Node was not a Text node; it is a %s. Use Text::isA()"
        " to check before calling Text::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Text&>(node);
}

//------------------------------------------------------------------------------
//                           XML COMMENT NODE
//------------------------------------------------------------------------------
Xml::Comment::Comment(const String& text)
    : Node(new tinyxml2::XMLDocument(text)) {}

Xml::Comment Xml::Comment::clone() const {
    ::tinyxml2::XMLComment* newComment = nullptr;

    if (getTiNodePtr()) {
        const ::tinyxml2::XMLComment* oldComment = getTiNodePtr()->ToComment();
        ::tinyxml2::XMLDocument* doc =
            const_cast<::tinyxml2::XMLDocument*>(getTiNodePtr()->GetDocument());

        if (oldComment && doc)
            newComment = oldComment->ShallowClone(doc)->ToComment();
    }

    return Xml::Comment(newComment);
}

/*static*/ bool Xml::Comment::isA(const Xml::Node& node) {
    if (!node.isValid()) return false;
    return node.getTiNode().ToComment() != 0;
}
/*static*/ const Xml::Comment& Xml::Comment::getAs(const Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Comment::getAs()",
        "The given Node was not a Comment node; it is a %s. Use Comment::isA()"
        " to check before calling Comment::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Comment&>(node);
}
/*static*/ Xml::Comment& Xml::Comment::getAs(Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Comment::getAs()",
        "The given Node was not a Comment node; it is a %s. Use Comment::isA()"
        " to check before calling Comment::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Comment&>(node);
}

//------------------------------------------------------------------------------
//                           XML UNKNOWN NODE
//------------------------------------------------------------------------------
Xml::Unknown::Unknown(const String& contents)
    : Node(new tinyxml2::XMLDocument()) {
    updTiNode().SetValue(contents);
}

Xml::Unknown Xml::Unknown::clone() const {
    ::tinyxml2::XMLUnknown* newUnknown = nullptr;

    if (getTiNodePtr()) {
        const ::tinyxml2::XMLUnknown* oldUnknown = getTiNodePtr()->ToUnknown();
        ::tinyxml2::XMLDocument* doc = const_cast<::tinyxml2::XMLDocument*>(
            getTiNodePtr()
                ->GetDocument());  // Remove const to call ShallowClone()

        if (oldUnknown && doc)
            newUnknown = oldUnknown->ShallowClone(doc)->ToUnknown();
    }

    return Xml::Unknown(newUnknown);
}

/*static*/ bool Xml::Unknown::isA(const Xml::Node& node) {
    if (!node.isValid()) return false;
    return node.getTiNode().ToUnknown() != 0;
}
/*static*/ const Xml::Unknown& Xml::Unknown::getAs(const Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Unknown::getAs()",
        "The given Node was not an Unknown node; it is a %s. Use Unknown::isA()"
        " to check before calling Unknown::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Unknown&>(node);
}
/*static*/ Xml::Unknown& Xml::Unknown::getAs(Node& node) {
    SimTK_ERRCHK1_ALWAYS(
        isA(node), "Xml::Unknown::getAs()",
        "The given Node was not an Unknown node; it is a %s. Use Unknown::isA()"
        " to check before calling Unknown::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Unknown&>(node);
}
