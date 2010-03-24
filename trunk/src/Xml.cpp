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
    // or has top-level text, we'll add <XMLDocument> elements </XMLDocument> to 
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
            // XMLDocument element to the end of whatever's there.
            TiXmlElement* root = new TiXmlElement("XMLDocument");
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
        TiXmlElement* root = new TiXmlElement("XMLDocument");

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

    // Note that Tiny XML "document" is the entire XML object, while ours
    // is the document tag or what TinyXML calls the "root element". There
    // is only supposed to be one root element; if we see more than one
    // we insert a new top-level element <XML> original stuff </XML>.
    TiXmlDocument   m_tixml;

    // These contain pointers into tixml. They are filled in when the
    // document is initially canonicalized.
    Xml::Element   m_rootElement;
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

const Xml::Element& Xml::getDocumentElement() const {
    assert(getImpl().m_rootElement.isValid());
    return getImpl().m_rootElement;
}
Xml::Element& Xml::updDocumentElement() {
    assert(getImpl().m_rootElement.isValid());
    return updImpl().m_rootElement;
}

void Xml::readFromFile(const String& pathname) {
    updImpl().readFromFile(pathname.c_str());
    updImpl().canonicalizeDocument();
}
void Xml::writeToFile(const String& pathname) const {
    getImpl().writeToFile(pathname.c_str());
}

void Xml::readFromString(const char* xmlDocument) {
    updImpl().readFromString(xmlDocument);
    updImpl().canonicalizeDocument();
}
void Xml::readFromString(const String& xmlDocument) 
{   readFromString(xmlDocument.c_str()); }

void Xml::writeToString(String& xmlDocument, bool compact) const {
    getImpl().writeToString(xmlDocument, compact);
}

const String& Xml::getDocumentTag() const {
    return getDocumentElement().getElementTag();
}
void Xml::setDocumentTag(const String& tag) {
    updDocumentElement().setElementTag(tag);
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
Xml::const_node_iterator Xml::node_begin(NodeType allowed) const {
    const TiXmlNode* first = getImpl().m_tixml.FirstChild();
    while (first && !nodeTypeIsAllowed(allowed, first->Type()))
        first = first->NextSibling();
    return const_node_iterator(first, allowed);
}

Xml::node_iterator Xml::node_begin(NodeType allowed) {
    TiXmlNode* first = updImpl().m_tixml.FirstChild();
    while (first && !nodeTypeIsAllowed(allowed, first->Type()))
        first = first->NextSibling();
    return node_iterator(first, allowed);
}

    // XML node_end()
Xml::const_node_iterator Xml::node_end() const 
{   return const_node_iterator(0); }
Xml::node_iterator Xml::node_end() 
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
//                      XML CONST ATTRIBUTE ITERATOR
//------------------------------------------------------------------------------
Xml::const_attribute_iterator& Xml::const_attribute_iterator::
operator++() {
    const TiXmlAttribute* next = attr.getTiAttr().Next();
    reassign(next);
    return *this;
}

Xml::const_attribute_iterator Xml::const_attribute_iterator::
operator++(int) {
    const Attribute save(attr);
    const TiXmlAttribute* next = attr.getTiAttr().Next();
    reassign(next);
    return const_attribute_iterator(save);
}

Xml::const_attribute_iterator& Xml::const_attribute_iterator::
operator--() {
    const TiXmlAttribute* prev = attr.getTiAttr().Previous();
    reassign(prev);
    return *this;
}

Xml::const_attribute_iterator Xml::const_attribute_iterator::
operator--(int) {
    const Attribute save(attr);
    const TiXmlAttribute* prev = attr.getTiAttr().Previous();
    reassign(prev);
    return const_attribute_iterator(save);
}



//------------------------------------------------------------------------------
//                                 XML NODE
//------------------------------------------------------------------------------                  
const String& Xml::Node::getValue() const {return getTiNode().ValueStr();}
Xml::NodeType Xml::Node::getNodeType() const {
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

bool Xml::Node::empty(NodeType allowed) const 
{   return node_begin(allowed) == node_end(); }


    // NODE node_begin()
Xml::const_node_iterator Xml::Node::node_begin(NodeType allowed) const {
    const TiXmlNode* first = getTiNode().FirstChild();
    while (first && !nodeTypeIsAllowed(allowed, first->Type()))
        first = first->NextSibling();
    return const_node_iterator(first, allowed);
}

Xml::node_iterator Xml::Node::node_begin(NodeType allowed) {
    TiXmlNode* first = updTiNode().FirstChild();
    while (first && !nodeTypeIsAllowed(allowed, first->Type()))
        first = first->NextSibling();
    return node_iterator(first, allowed);
}

    // NODE node_end()
Xml::const_node_iterator Xml::Node::node_end() const 
{   return const_node_iterator(0); }
Xml::node_iterator Xml::Node::node_end() 
{   return node_iterator(0); }


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
//                          XML CONST NODE ITERATOR
//------------------------------------------------------------------------------
Xml::const_node_iterator& Xml::const_node_iterator::
operator++() {
    TiXmlNode* next = node.updTiNode().NextSibling();
    while (next && !nodeTypeIsAllowed(allowed, next->Type()))
        next = next->NextSibling();
    node = Node(next);
    return *this;
}

Xml::const_node_iterator Xml::const_node_iterator::
operator++(int) {
    const Node save(node);
    TiXmlNode* next = node.updTiNode().NextSibling();
    while (next && !nodeTypeIsAllowed(allowed, next->Type()))
        next = next->NextSibling();
    node = Node(next);
    return const_node_iterator(save);
}

Xml::const_node_iterator& Xml::const_node_iterator::
operator--() {
    TiXmlNode* prev = node.updTiNode().PreviousSibling();
    while (prev && !nodeTypeIsAllowed(allowed, prev->Type()))
        prev = prev->PreviousSibling();
    node = Node(prev);
    return *this;
}

Xml::const_node_iterator Xml::const_node_iterator::
operator--(int) {
    const Node save(node);
    TiXmlNode* prev = node.updTiNode().PreviousSibling();
    while (prev && !nodeTypeIsAllowed(allowed, prev->Type()))
        prev = prev->PreviousSibling();
    node = Node(prev);
    return const_node_iterator(save);
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

const String& Xml::Element::getElementTag() const {
    return getTiNode().ValueStr();
}

void Xml::Element::setElementTag(const String& type) {
    updTiNode().SetValue(type);
}

bool Xml::Element::isTextElement() const {
    if (element_begin() != element_end()) return false; // has child elements
    const_node_iterator text = node_begin(TextNode);
    return text == node_end() || ++text == node_end(); // zero or one
}

const String& Xml::Element::getElementText() const {
    const String null;
    SimTK_ERRCHK1_ALWAYS(isTextElement(), "Xml::Element::getElementText()",
        "Element <%s> is not a text element.", getElementTag().c_str());

    const_node_iterator text = node_begin(TextNode);
    return text == node_end() ? null : text->getValue();
}

bool Xml::Element::hasAttribute(const String& name) const {
    for (const_attribute_iterator p=attribute_begin();
            p != attribute_end(); ++p)
        if (p->getName() == name) return true;
    return false;
}

bool Xml::Element::hasElement(const String& tag) const {
    return element_begin(tag) != element_end();
}

Xml::Element Xml::Element::getRequiredElement(const String& tag) const {
    const_element_iterator p = element_begin(tag);
    SimTK_ERRCHK2_ALWAYS(p != element_end(), 
        "Xml::Element::getRequiredElement()",
        "Couldn't find required child element <%s> in element <%s>.",
        tag.c_str(), getElementTag().c_str());
    return const_cast<Element&>(*p);
}


Xml::Element Xml::Element::getOptionalElement(const String& tag) const {
    const_element_iterator p = element_begin(tag);
    return p != element_end() ? const_cast<Element&>(*p) : Element(0);
}

Xml::Attribute Xml::Element::getRequiredAttribute(const String& name) const {
    const_attribute_iterator p = find_attribute(name);
    SimTK_ERRCHK2_ALWAYS(p != attribute_end(), 
        "Xml::Element::getRequiredAttribute()",
        "Couldn't find required attribute %s in element <%s>.",
        name.c_str(), getElementTag().c_str());
    return const_cast<Attribute&>(*p);
}


/*static*/ bool Xml::Element::isA(const Xml::Node& node) {
    if (!node.isValid()) return false;
    return node.getTiNode().ToElement() != 0;
}
/*static*/const Xml::Element& Xml::Element::getAs(const Node& node) {
    assert(isA(node));
    return reinterpret_cast<const Element&>(node);
}
/*static*/Xml::Element& Xml::Element::updAs(Node& node) {
    assert(isA(node));
    return reinterpret_cast<Element&>(node);
}

    // Element begin()
Xml::element_iterator Xml::Element::
element_begin(const String& tag) {
    TiXmlElement* first = updTiNode().FirstChildElement();
    while (first && !elementIsAllowed(tag, first))
        first = first->NextSiblingElement();
    return element_iterator(first, tag);
}
Xml::const_element_iterator Xml::Element::
element_begin(const String& tag) const {
    Element* mthis = const_cast<Element*>(this);
    TiXmlElement* first = mthis->updTiNode().FirstChildElement();
    while (first && !elementIsAllowed(tag, first))
        first = first->NextSiblingElement();
    return const_element_iterator(first, tag);
}

    // Element end()
Xml::const_element_iterator Xml::Element::element_end() const 
{   return const_element_iterator(Element(0)); }
Xml::element_iterator Xml::Element::element_end() 
{   return element_iterator(Element(0));}

    // Attribute begin()
Xml::attribute_iterator Xml::Element::
attribute_begin() {
    TiXmlAttribute* first = updTiElement().FirstAttribute();
    return attribute_iterator(first);
}
Xml::const_attribute_iterator Xml::Element::
attribute_begin() const {
    const TiXmlAttribute* first = getTiElement().FirstAttribute();
    return const_attribute_iterator(first);
}

    // Attribute end()
Xml::const_attribute_iterator Xml::Element::
attribute_end() const 
{   return const_attribute_iterator(0); }
Xml::attribute_iterator Xml::Element::
attribute_end() 
{   return attribute_iterator(0); }




//------------------------------------------------------------------------------
//                          XML ELEMENT ITERATOR
//------------------------------------------------------------------------------
Xml::element_iterator& Xml::element_iterator::
operator++() {
    TiXmlElement* next = elt.updTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag,next))
        next = next->NextSiblingElement();
    elt.setTiElementPtr(next);
    return *this;
}

Xml::element_iterator Xml::element_iterator::
operator++(int) {
    Element save(elt);
    TiXmlElement* next = elt.updTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag,next))
        next = next->NextSiblingElement();
    elt.setTiElementPtr(next);
    return element_iterator(save);
}

Xml::element_iterator& Xml::element_iterator::
operator--() {
    TiXmlElement* prev = elt.updTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag,prev))
        prev = prev->PreviousSiblingElement();
    elt.setTiElementPtr(prev);
    return *this;
}

Xml::element_iterator Xml::element_iterator::
operator--(int) {
    Element save(elt);
    TiXmlElement* prev = elt.updTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag,prev))
        prev = prev->PreviousSiblingElement();
    elt.setTiElementPtr(prev);
    return element_iterator(save);
}

//------------------------------------------------------------------------------
//                       XML CONST ELEMENT ITERATOR
//------------------------------------------------------------------------------
Xml::const_element_iterator& Xml::const_element_iterator::
operator++() {
    const TiXmlElement* next = elt.getTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag,next))
        next = next->NextSiblingElement();
    reassign(next);
    return *this;
}

Xml::const_element_iterator Xml::const_element_iterator::
operator++(int) {
    const Element save(elt);
    const TiXmlElement* next = elt.getTiElement().NextSiblingElement();
    while (next && !elementIsAllowed(tag,next))
        next = next->NextSiblingElement();
    reassign(next);
    return const_element_iterator(save);
}

Xml::const_element_iterator& Xml::const_element_iterator::
operator--() {
    const TiXmlElement* prev = elt.getTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag,prev))
        prev = prev->PreviousSiblingElement();
    reassign(prev);
    return *this;
}

Xml::const_element_iterator Xml::const_element_iterator::
operator--(int) {
    const Element save(elt);
    const TiXmlElement* prev = elt.getTiElement().PreviousSiblingElement();
    while (prev && !elementIsAllowed(tag,prev))
        prev = prev->PreviousSiblingElement();
    reassign(prev);
    return const_element_iterator(save);
}
