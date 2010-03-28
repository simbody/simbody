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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Plugin.h"
#include "SimTKcommon/internal/Xml.h"

#include "tinyxml.h"

using namespace SimTK;

//------------------------------------------------------------------------------
//                                 XML IMPL
//------------------------------------------------------------------------------
class Xml::Impl {
public:
    Impl() {} // not canonicalized yet
    Impl(const String& pathname) {
        readFromFile(pathname.c_str());
    }
    
    void clear() {
        m_tixml.Clear();
        m_rootElement.clear();
    }

    void readFromFile(const String& pathname) {
        clear();
        m_tixml.SetValue(pathname);
        bool loadOK = m_tixml.LoadFile();
        SimTK_ERRCHK2_ALWAYS(loadOK, "Xml::readFromFile()",
            "Failed to load the Xml file '%s' with error '%s'.",
            pathname.c_str(), m_tixml.ErrorMsg().c_str());
    }

    void writeToFile(const String& pathname) const {
        bool saveOK = m_tixml.SaveFile(pathname);
        SimTK_ERRCHK2_ALWAYS(saveOK, "Xml::writeToFile()",
            "Failed to write to the Xml file '%s' with error '%s'.",
            pathname.c_str(), m_tixml.ErrorMsg().c_str());
    }

    void readFromString(const char* xmlDocument) {
        clear();
        m_tixml.Parse(xmlDocument);
        SimTK_ERRCHK1_ALWAYS(!m_tixml.Error(), "Xml::readFromString()",
            "Failed to parse the Xml string with error '%s'.",
            m_tixml.ErrorMsg().c_str());
    }

    void writeToString(String& xmlDocument, bool compact) const {
        TiXmlPrinter printer(xmlDocument);
	    if (compact) printer.SetStreamPrinting();
	    m_tixml.Accept( &printer );
    }

    // Call this during construction and after a new Xml document has been
    // parsed. It guarantees: (1) there is a Declaration record and it is
    // the first node at the top level of the Xml document, and (2) there
    // is no top-level text and only one top-level element and that is the 
    // "root" element whose tag name is the document type.
    void canonicalizeDocument() {
        TiXmlDeclaration* decl = addDeclarationIfNeeded();
        TiXmlElement*     root = addRootElementIfNeeded();

        m_rootElement.setTiNodePtr(root);
    }

    const TiXmlDeclaration& getTiXmlDeclaration() const {
        const TiXmlNode* decl = m_tixml.FirstChild();
        assert(decl && decl->Type()==TiXmlNode::DECLARATION);
        return *decl->ToDeclaration();
    }

    TiXmlDeclaration& updTiXmlDeclaration() {
        TiXmlNode* decl = m_tixml.FirstChild();
        assert(decl && decl->Type()==TiXmlNode::DECLARATION);
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
    TiXmlDeclaration* addDeclarationIfNeeded() {
        TiXmlNode* child=m_tixml.FirstChild();
        if (child && child->Type() == TiXmlNode::DECLARATION)
            return child->ToDeclaration(); // the easy and most common case

        // Otherwise hunt for the declaration.
        for (; child; child=child->NextSibling())
            if (child->Type() == TiXmlNode::DECLARATION)
                break;

        if (child)
            m_tixml.DisconnectChild(child); // it's in the wrong place
        else child = new TiXmlDeclaration("1.0", "UTF-8", "");

        // Insert new node as the first in the document.
        m_tixml.LinkBeginChild(child);
        return child->ToDeclaration();
    }

    // If the supplied Xml document has zero or more than one top-level element,
    // or has top-level text, we'll add <_Root> elements </_Root> to 
    // encapsulate all the top-level text and element nodes (this may also
    // surround top-level comments and unknowns if they occur between the text
    // and elements. A pointer to the root element (whether original or added) 
    // is returned.
    TiXmlElement* addRootElementIfNeeded() {
        // Find the first element or text node, and remember the node that
        // preceded it since that's where we may be inserting the new root
        // element.
        TiXmlNode* nodeBeforeFirst = 0;
        TiXmlNode* firstEltOrText = m_tixml.FirstChild();
        while (  firstEltOrText && 
               !(firstEltOrText->ToElement()||firstEltOrText->ToText()))
        {   nodeBeforeFirst = firstEltOrText;
            firstEltOrText = firstEltOrText->NextSibling(); }

        if (!firstEltOrText) {
            // No top level element or text node. We'll just append an empty
            // _Root element to the end of whatever's there.
            TiXmlElement* root = new TiXmlElement("_Root");
            m_tixml.LinkEndChild(root);
            return root;
        }

        // There is at least one element or text node at the top level; the 
        // first one is pointed to by firstEltOrText, and the node just before 
        // that (if any) is pointed to by nodeBeforeFirst. Now find the last 
        // top-level element or text node.
        TiXmlNode* lastEltOrText = m_tixml.LastChild();
        while (  lastEltOrText && 
               !(lastEltOrText->ToElement()||lastEltOrText->ToText()))
            lastEltOrText = lastEltOrText->PreviousSibling();

        assert(lastEltOrText); // should have at least re-found the first one!

        // If the extremely likely case that the first and last are the same 
        // node and that node is an element, then the document is already in the
        // right format and we don't have to do anything to it.
        if (firstEltOrText==lastEltOrText && firstEltOrText->ToElement())
            return firstEltOrText->ToElement();

        // Now we know there is top-level text or more than one top-level
        // element so we are going to have to surround everything between
        // first and last with a new root element.
        TiXmlElement* root = new TiXmlElement("_Root");

        TiXmlNode* nextToMove = firstEltOrText;
        while(true) {
            assert(nextToMove); // can't happen!
            TiXmlNode* moveMe = nextToMove;
            nextToMove = moveMe->NextSibling();
            root->LinkEndChild(m_tixml.DisconnectChild(moveMe));
            // Did we just move the last element or text node?
            if (moveMe == lastEltOrText) break;
        }

        // Now link the new root element right where we found the first
        // element or text node.
        if (nodeBeforeFirst) m_tixml.LinkAfterChild(nodeBeforeFirst, root);
        else                 m_tixml.LinkBeginChild(root);

        return root;
    }

    // The Tiny XML "document" is a TiXmlNode representing the entire XML 
    // object. We instead represent the document as an object of class Xml 
    // which is not itself a Node.
    TiXmlDocument   m_tixml;

    // This references the root element withing the TinyXml document. It
    // is filled in when the document is initially canonicalized.
    Xml::Element    m_rootElement;
};



//------------------------------------------------------------------------------
//                                    XML
//------------------------------------------------------------------------------

// Handy helper for weeding out unwanted nodes.
static bool nodeTypeIsAllowed(Xml::NodeType       allowed,
                              TiXmlNode::NodeType found) {
    switch(found) {
    case TiXmlNode::ELEMENT: return (allowed & Xml::ElementNode)!=0;
    case TiXmlNode::TEXT:    return (allowed & Xml::TextNode)   !=0;
    case TiXmlNode::COMMENT: return (allowed & Xml::CommentNode)!=0;
    case TiXmlNode::UNKNOWN: return (allowed & Xml::UnknownNode)!=0;
    default: return false;
    }
}

/*static*/String Xml::getNodeTypeAsString(Xml::NodeType type) {
    // Take care of special cases first.
    switch(type) {
    case NoNode: return "NoNode";
    case NoJunkNodes: return "NoJunkNodes";
    case JunkNodes: return "JunkNodes";
    case AnyNodes: return "AnyNodes";
    default: ; // fall through
    }

    // "Or" bits together in CUTE order.
    String out;
    if (type & CommentNode) out = "CommentNode";
    if (type & UnknownNode) 
    {   if (!out.empty()) out += "|"; out += "UnknownNode"; }
    if (type & TextNode) 
    {   if (!out.empty()) out += "|"; out += "TextNode"; }
    if (type & ElementNode) 
    {   if (!out.empty()) out += "|"; out += "ElementNode"; }

    return out;
}

Xml::Xml() : impl(0) {
    impl = new Impl();
    impl->canonicalizeDocument();
}

Xml::Xml(const String& pathname) : impl(0) {
    impl = new Impl(pathname);
    impl->canonicalizeDocument();
}


void Xml::readFromFile(const String& pathname) 
{   updImpl().readFromFile(pathname.c_str());
    updImpl().canonicalizeDocument(); }
void Xml::writeToFile(const String& pathname) const 
{   getImpl().writeToFile(pathname.c_str()); }
void Xml::readFromString(const char* xmlDocument) 
{   updImpl().readFromString(xmlDocument);
    updImpl().canonicalizeDocument(); }
void Xml::readFromString(const String& xmlDocument) 
{   readFromString(xmlDocument.c_str()); }
void Xml::writeToString(String& xmlDocument, bool compact) const 
{   getImpl().writeToString(xmlDocument, compact); }

Xml::Element Xml::getRootElement() 
{   assert(getImpl().m_rootElement.isValid());
    return updImpl().m_rootElement; }
const String& Xml::getRootTag() const 
{   return unconst().getRootElement().getElementTag(); }
void Xml::setRootTag(const String& tag) 
{   getRootElement().setElementTag(tag); }

void Xml::insertTopLevelNodeAfter (const Xml::node_iterator& afterThis, 
                                   Xml::Node                 insertThis) {
    const char* method = "Xml::insertTopLevelNodeAfter()";

    // Check that the supplied Node is OK.
    SimTK_ERRCHK_ALWAYS(insertThis.isValid(), method,
        "The supplied Node handle was empty.");
    SimTK_ERRCHK_ALWAYS(insertThis.isOrphan(), method,
        "The Node was not an orphan so can't be inserted.");
    SimTK_ERRCHK1_ALWAYS(Comment::isA(insertThis) || Unknown::isA(insertThis),
        method, "The Node had NodeType %s, but only Comment and Unknown nodes"
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

    updImpl().m_tixml.LinkAfterChild(afterThis->updTiNodePtr(),
                                     insertThis.updTiNodePtr());
}

void Xml::insertTopLevelNodeBefore(const Xml::node_iterator& beforeThis, 
                                   Xml::Node                 insertThis) {
    const char* method = "Xml::insertTopLevelNodeBefore()";

    // Check that the supplied Node is OK.
    SimTK_ERRCHK_ALWAYS(insertThis.isValid(), method,
        "The supplied Node handle was empty.");
    SimTK_ERRCHK_ALWAYS(insertThis.isOrphan(), method,
        "The Node was not an orphan so can't be inserted.");
    SimTK_ERRCHK1_ALWAYS(Comment::isA(insertThis) || Unknown::isA(insertThis),
        method, "The Node had NodeType %s, but only Comment and Unknown nodes"
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

    updImpl().m_tixml.LinkBeforeChild(beforeThis->updTiNodePtr(),
                                      insertThis.updTiNodePtr());
}


String Xml::getXmlVersion() const 
{   return getImpl().getTiXmlDeclaration().Version(); }
String Xml::getXmlEncoding() const 
{   return getImpl().getTiXmlDeclaration().Encoding(); }
bool Xml::getXmlIsStandalone() const 
{   return String::toLower(getImpl().getTiXmlDeclaration().Standalone())!="no"; }

void Xml::setXmlVersion(const String& version)
{   updImpl().updTiXmlDeclaration().SetVersion(version.c_str()); }
void Xml::setXmlEncoding(const String& encoding)
{   updImpl().updTiXmlDeclaration().SetEncoding(encoding.c_str()); }
// For the standalone case we set the string to "" rather than "yes" so that
// the standalone attribute won't appear in the output declaration.
void Xml::setXmlIsStandalone(bool isStandalone)
{   updImpl().updTiXmlDeclaration().SetStandalone(isStandalone ? "" : "no"); }


    // XML node_begin()
Xml::node_iterator Xml::node_begin(NodeType allowed) {
    TiXmlNode* first = updImpl().m_tixml.FirstChild();
    while (first && !nodeTypeIsAllowed(allowed, first->Type()))
        first = first->NextSibling();
    return node_iterator(first, allowed);
}

    // XML node_end()
Xml::node_iterator Xml::node_end() const 
{   return node_iterator(0); }



//------------------------------------------------------------------------------
//                              XML ATTRIBUTE
//------------------------------------------------------------------------------
Xml::Attribute::Attribute(const String& name, const String& value) 
:   tiAttr(new TiXmlAttribute(name,value)) {}

void Xml::Attribute::clear() {
    if (!tiAttr) return;
    if (!tiAttr->GetDocument())
        delete tiAttr; // not part of any document
    tiAttr = 0;
}

const String& Xml::Attribute::getName() const {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::getName()",
        "The attribute handle was empty.");
    return getTiAttr().NameStr();
}
const String& Xml::Attribute::getValue() const {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::getValue()",
        "The attribute handle was empty.");
    return getTiAttr().ValueStr();
}

Xml::Attribute& Xml::Attribute::setName(const String& name) {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::setName()",
        "The attribute handle was empty.");
    updTiAttr().SetName(name);
    return *this;
}
Xml::Attribute& Xml::Attribute::setValue(const String& value) {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Attribute::setValue()",
        "The attribute handle was empty.");
    updTiAttr().SetValue(value);
    return *this;
}

void Xml::Attribute::writeToString(String& out) const {
    getTiAttr().Print(0,0,&out);
}

//------------------------------------------------------------------------------
//                         XML ATTRIBUTE ITERATOR
//------------------------------------------------------------------------------
Xml::attribute_iterator& Xml::attribute_iterator::
operator++() {
    TiXmlAttribute* next = attr.updTiAttr().Next();
    attr.setTiAttrPtr(next);
    return *this;
}

Xml::attribute_iterator Xml::attribute_iterator::
operator++(int) {
    Attribute save(attr);
    TiXmlAttribute* next = attr.updTiAttr().Next();
    attr.setTiAttrPtr(next);
    return attribute_iterator(save);
}

Xml::attribute_iterator& Xml::attribute_iterator::
operator--() {
    TiXmlAttribute* prev = attr.updTiAttr().Previous();
    attr.setTiAttrPtr(prev);
    return *this;
}

Xml::attribute_iterator Xml::attribute_iterator::
operator--(int) {
    Attribute save(attr);
    TiXmlAttribute* prev = attr.updTiAttr().Previous();
    attr.setTiAttrPtr(prev);
    return attribute_iterator(save);
}



//------------------------------------------------------------------------------
//                                 XML NODE
//------------------------------------------------------------------------------                  
void Xml::Node::clear()
{   if (isOrphan()) delete tiNode;
    tiNode = 0; }


Xml::NodeType Xml::Node::getNodeType() const {
    if (!isValid()) return NoNode;
    switch(getTiNode().Type()) {
    case TiXmlNode::COMMENT: return CommentNode;
    case TiXmlNode::UNKNOWN: return UnknownNode;
    case TiXmlNode::TEXT:    return TextNode;
    case TiXmlNode::ELEMENT: return ElementNode;
    default: SimTK_ASSERT1_ALWAYS(false,
        "Xml::Node::getNodeType(): can't convert TinyXML node type %s to any"
        " SimTK::Xml node type.", getTiNode().TypeName());
    }
    return NoNode;
}

String Xml::Node::
getNodeTypeAsString() const {return Xml::getNodeTypeAsString(getNodeType());}

const String& Xml::Node::getNodeText() const {
    SimTK_ERRCHK_ALWAYS(isValid(), "Xml::Node::getText()",
        "Can't get text from an empty Node handle.");

    return getTiNode().ValueStr();
}

bool Xml::Node::isTopLevelNode() const 
{   if (!isValid()) return false;
    const TiXmlNode& n = getTiNode();
    return n.Parent() && n.Parent() == n.GetDocument(); }

// Note that the criteria for orphanhood is that the *TiXmlNode* has no
// parent. In SimTK::Xml, none of the top-level nodes are considered to have
// a parent since the Xml document is not a node there, while in TinyXML
// the TiXmlDocument is a TiXmlNode, so even top-level nodes have a parent.
bool Xml::Node::isOrphan() const
{   if (!isValid()) return false; // empty handle not considered an orphan
    return getTiNode().Parent() == 0; }

void Xml::Node::
writeToString(String& out, bool compact) const {
    TiXmlPrinter printer(out);
    if (compact) printer.SetStreamPrinting();
    getTiNode().Accept( &printer );
}

bool Xml::Node::hasParentElement() const 
{   if (!isValid()) return false;
    const TiXmlNode* parent = getTiNode().Parent();
    return parent && parent->ToElement(); }

Xml::Element Xml::Node::getParentElement() {
    SimTK_ERRCHK_ALWAYS(hasParentElement(),
        "Xml::Node::hasParentElement()", 
        "This node does not have a parent element; it may be a top-level"
        " node, an orphan, or just an empty Node handle.");
    return Element(updTiNode().Parent()->ToElement());
}


//------------------------------------------------------------------------------
//                          XML NODE ITERATOR
//------------------------------------------------------------------------------
Xml::node_iterator& Xml::node_iterator::
operator++() {
    TiXmlNode* next = node.updTiNode().NextSibling();
    while (next && !nodeTypeIsAllowed(allowed, next->Type()))
        next = next->NextSibling();
    node = Node(next);
    return *this;
}

Xml::node_iterator Xml::node_iterator::
operator++(int) {
    Node save(node);
    TiXmlNode* next = node.updTiNode().NextSibling();
    while (next && !nodeTypeIsAllowed(allowed, next->Type()))
        next = next->NextSibling();
    node = Node(next);
    return node_iterator(save);
}

Xml::node_iterator& Xml::node_iterator::
operator--() {
    TiXmlNode* prev = node.updTiNode().PreviousSibling();
    while (prev && !nodeTypeIsAllowed(allowed, prev->Type()))
        prev = prev->PreviousSibling();
    node = Node(prev);
    return *this;
}

Xml::node_iterator Xml::node_iterator::
operator--(int) {
    Node save(node);
    TiXmlNode* prev = node.updTiNode().PreviousSibling();
    while (prev && !nodeTypeIsAllowed(allowed, prev->Type()))
        prev = prev->PreviousSibling();
    node = Node(prev);
    return node_iterator(save);
}



//------------------------------------------------------------------------------
//                              XML ELEMENT
//------------------------------------------------------------------------------

// Handy helper for weeding out unwanted nodes and elements.
static bool elementIsAllowed(const String& tag,
                             const TiXmlElement* elt) {
    if (elt==0) return false;
    return tag.empty() || elt->ValueStr() == tag;
}

Xml::Element::Element(const String& tag, const String& value) 
:   Node(new TiXmlElement(tag)) {
    if (value.empty()) return;
    // We need to add a Text node.
    updTiElement().LinkEndChild(new TiXmlText(value));
}

const String& Xml::Element::getElementTag() const 
{   return getTiNode().ValueStr(); }
void Xml::Element::setElementTag(const String& type) 
{   updTiNode().SetValue(type); }

bool Xml::Element::isValueElement() const {
    if (unconst().element_begin() != element_end()) 
        return false; // has child elements
    node_iterator text = unconst().node_begin(TextNode);
    return text == node_end() || ++text == node_end(); // zero or one
}

const String& Xml::Element::getValue() const {
    const String null;
    SimTK_ERRCHK1_ALWAYS(isValueElement(), "Xml::Element::getValue()",
        "Element <%s> is not a value element.", getElementTag().c_str());

    node_iterator text = unconst().node_begin(TextNode);
    return text == node_end() ? null : text->getNodeText();
}

// Must add a Text node now if this Element doesn't have one.
String& Xml::Element::updValue() {
    SimTK_ERRCHK1_ALWAYS(isValueElement(), "Xml::Element::getValue()",
        "Element <%s> is not a value element.", getElementTag().c_str());

    node_iterator text = node_begin(TextNode);
    if (text != node_end()) return Text::getAs(*text).updText();

    // We need to add a Text node.
    TiXmlText* textp = new TiXmlText("");
    updTiElement().LinkEndChild(textp);
    return textp->UpdValueStr();
}

// If there is no Text node we'll add one; if there is just one we'll
// change its value; otherwise report an error.
void Xml::Element::setValue(const String& value) {
    SimTK_ERRCHK1_ALWAYS(isValueElement(), "Xml::Element::setValue()",
        "Element <%s> is not a value element.", getElementTag().c_str());

    node_iterator text = node_begin(TextNode);
    if (text == node_end()) updTiNode().LinkEndChild(new TiXmlText(value));
    else                    text->updTiNode().SetValue(value);
}

bool Xml::Element::hasAttribute(const String& name) const 
{   return unconst().find_attribute(name) != attribute_end(); }

bool Xml::Element::hasElement(const String& tag) const 
{   return unconst().element_begin(tag) != element_end(); }

bool Xml::Element::hasNode(NodeType allowed) const 
{   return unconst().node_begin(allowed) != node_end(); }

Xml::Element Xml::Element::getRequiredElement(const String& tag) {
    element_iterator p = element_begin(tag);
    SimTK_ERRCHK2_ALWAYS(p != element_end(), 
        "Xml::Element::getRequiredElement()",
        "Couldn't find required child element <%s> in element <%s>.",
        tag.c_str(), getElementTag().c_str());
    return *p;
}

Xml::Element Xml::Element::getOptionalElement(const String& tag) {
    element_iterator p = element_begin(tag);
    return p != element_end() ? *p : Element(0);
}

Xml::Attribute Xml::Element::getRequiredAttribute(const String& name) {
    attribute_iterator p = find_attribute(name);
    SimTK_ERRCHK2_ALWAYS(p != attribute_end(), 
        "Xml::Element::getRequiredAttribute()",
        "Couldn't find required attribute %s in element <%s>.",
        name.c_str(), getElementTag().c_str());
    return *p;
}


/*static*/ bool Xml::Element::isA(const Xml::Node& node) 
{   if (!node.isValid()) return false;
    return node.getTiNode().ToElement() != 0; }
/*static*/const Xml::Element& Xml::Element::getAs(const Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Element::getAs()",
        "The given Node was not an Element; it is a %s. Use Element::isA()"
        " to check before calling Element::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Element&>(node); }
/*static*/Xml::Element& Xml::Element::getAs(Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Element::getAs()",
        "The given Node was not an Element; it is a %s. Use Element::isA()"
        " to check before calling Element::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Element&>(node); }

void Xml::Element::insertNodeBefore(const node_iterator& beforeThis, Node node) {
    const char* method = "Xml::Element::insertNodeBefore()";
    const char* tag    = getElementTag().c_str();

    SimTK_ERRCHK1_ALWAYS(node.isValid(), method,
        "The supplied Node handle was invalid so can't be inserted into"
        " Element <%s>.", tag);
    SimTK_ERRCHK1_ALWAYS(!hasParentElement(), method,
        "The supplied Node already had a parent so can't be inserted into"
        " Element <%s>.", tag);

    if (beforeThis == node_end()) {
        updTiNode().LinkEndChild(node.updTiNodePtr());
        return;
    }

    SimTK_ERRCHK1_ALWAYS(beforeThis->getParentElement() == *this, method,
        "The supplied node_iterator referred to a node that was not a"
        "child node of this Element <%s>.", tag);

    TiXmlNode* p = beforeThis->updTiNodePtr();
    p->LinkBeforeChild(p, node.updTiNodePtr());
}

void Xml::Element::insertNodeAfter(const node_iterator& afterThis, Node node) {
    const char* method = "Xml::Element::insertNodeAfter()";
    const char* tag    = getElementTag().c_str();

    SimTK_ERRCHK1_ALWAYS(node.isValid(), method,
        "The supplied Node handle was invalid so can't be inserted into"
        " Element <%s>.", tag);
    SimTK_ERRCHK1_ALWAYS(!node.hasParentElement(), method,
        "The supplied Node already had a parent so can't be inserted into"
        " Element <%s>.", tag);

    if (afterThis == node_end()) {
        updTiNode().LinkEndChild(node.updTiNodePtr());
        return;
    }

    SimTK_ERRCHK1_ALWAYS(afterThis->getParentElement() == *this, method,
        "The supplied node_iterator referred to a node that was not a"
        "child node of this Element <%s>.", tag);

    TiXmlNode* p = afterThis->updTiNodePtr();
    p->LinkAfterChild(p, node.updTiNodePtr());
}



    // Element node_begin()
Xml::node_iterator Xml::Element::node_begin(NodeType allowed) {
    TiXmlNode* first = updTiNode().FirstChild();
    while (first && !nodeTypeIsAllowed(allowed, first->Type()))
        first = first->NextSibling();
    return node_iterator(first, allowed);
}

    // Element node_end()
Xml::node_iterator Xml::Element::node_end() const
{   return node_iterator(0); }

    // Element begin()
Xml::element_iterator Xml::Element::
element_begin(const String& tag) {
    TiXmlElement* first = updTiNode().FirstChildElement();
    while (first && !elementIsAllowed(tag, first))
        first = first->NextSiblingElement();
    return element_iterator(first, tag);
}

    // Element end()
Xml::element_iterator Xml::Element::element_end() const 
{   return element_iterator(0);}

    // Attribute begin()
Xml::attribute_iterator Xml::Element::
attribute_begin() {
    TiXmlAttribute* first = updTiElement().FirstAttribute();
    return attribute_iterator(first);
}

    // Attribute end()
Xml::attribute_iterator Xml::Element::
attribute_end() const
{   return attribute_iterator(0); }




//------------------------------------------------------------------------------
//                          XML ELEMENT ITERATOR
//------------------------------------------------------------------------------
Xml::element_iterator& Xml::element_iterator::
operator++() {
    TiXmlElement* next = (*this)->updTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag,next))
        next = next->NextSiblingElement();
    reassign(next);
    return *this;
}

Xml::element_iterator Xml::element_iterator::
operator++(int) {
    Element save(*(*this));
    TiXmlElement* next = (*this)->updTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag,next))
        next = next->NextSiblingElement();
    reassign(next);
    return element_iterator(save);
}

Xml::element_iterator& Xml::element_iterator::
operator--() {
    TiXmlElement* prev = (*this)->updTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag,prev))
        prev = prev->PreviousSiblingElement();
    reassign(prev);
    return *this;
}

Xml::element_iterator Xml::element_iterator::
operator--(int) {
    Element save(*(*this));
    TiXmlElement* prev = (*this)->updTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag,prev))
        prev = prev->PreviousSiblingElement();
    reassign(prev);
    return element_iterator(save);
}




//------------------------------------------------------------------------------
//                             XML TEXT NODE
//------------------------------------------------------------------------------
Xml::Text::Text(const String& text) : Node(new TiXmlText(text)) {}

const String& Xml::Text::getText() const
{   return getTiNode().ValueStr(); }
String& Xml::Text::updText()
{   return updTiNode().UpdValueStr(); }

/*static*/ bool Xml::Text::isA(const Xml::Node& node) 
{   if (!node.isValid()) return false;
    return node.getTiNode().ToText() != 0; }
/*static*/const Xml::Text& Xml::Text::getAs(const Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Text::getAs()",
        "The given Node was not a Text node; it is a %s. Use Text::isA()"
        " to check before calling Text::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Text&>(node); }
/*static*/Xml::Text& Xml::Text::getAs(Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Text::getAs()",
        "The given Node was not a Text node; it is a %s. Use Text::isA()"
        " to check before calling Text::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Text&>(node); }



//------------------------------------------------------------------------------
//                           XML COMMENT NODE
//------------------------------------------------------------------------------
Xml::Comment::Comment(const String& text) : Node(new TiXmlComment(text)) {}

/*static*/ bool Xml::Comment::isA(const Xml::Node& node) 
{   if (!node.isValid()) return false;
    return node.getTiNode().ToComment() != 0; }
/*static*/const Xml::Comment& Xml::Comment::getAs(const Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Comment::getAs()",
        "The given Node was not a Comment node; it is a %s. Use Comment::isA()"
        " to check before calling Comment::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Comment&>(node); }
/*static*/Xml::Comment& Xml::Comment::getAs(Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Comment::getAs()",
        "The given Node was not a Comment node; it is a %s. Use Comment::isA()"
        " to check before calling Comment::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Comment&>(node); }



//------------------------------------------------------------------------------
//                           XML UNKNOWN NODE
//------------------------------------------------------------------------------
Xml::Unknown::Unknown(const String& contents) : Node(new TiXmlUnknown()) 
{   updTiNode().SetValue(contents); }

/*static*/ bool Xml::Unknown::isA(const Xml::Node& node) 
{   if (!node.isValid()) return false;
    return node.getTiNode().ToUnknown() != 0; }
/*static*/const Xml::Unknown& Xml::Unknown::getAs(const Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Unknown::getAs()",
        "The given Node was not an Unknown node; it is a %s. Use Unknown::isA()"
        " to check before calling Unknown::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<const Unknown&>(node); }
/*static*/Xml::Unknown& Xml::Unknown::getAs(Node& node) 
{   SimTK_ERRCHK1_ALWAYS(isA(node), "Xml::Unknown::getAs()",
        "The given Node was not an Unknown node; it is a %s. Use Unknown::isA()"
        " to check before calling Unknown::getAs().",
        node.getNodeTypeAsString().c_str());
    return reinterpret_cast<Unknown&>(node); }

