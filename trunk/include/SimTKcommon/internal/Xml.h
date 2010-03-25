#ifndef SimTK_SimTKCOMMON_XML_H_
#define SimTK_SimTKCOMMON_XML_H_

/* -------------------------------------------------------------------------- *
 *                          SimTK Core: SimTKcommon                           *
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
#include "SimTKcommon/internal/Array.h"
#include "SimTKcommon/internal/String.h"

#include <iterator>
#include <iostream>

namespace SimTK {

// These are declared but never defined; all TinyXML code is hidden.
class TiXmlNode; 
class TiXmlElement; 
class TiXmlAttribute;


/** This class provides a minimalist capability for reading and writing XML 
documents, as files or strings. This is based with gratitude on the excellent 
open source XML parser TinyXML (http://www.grinninglizard.com/tinyxml/). Note 
that this is a <em>non-validating</em> parser, meaning it deals only with the 
XML file itself and not with a schema or other description of the XML file's 
expected contents. Instead, the structure of your code that uses this class 
encodes the expected structure and contents of the XML document.

Our in-memory model of an XML document is simplified even further than 
TinyXML's. There is much more to know about XML; you could start here: 
http://en.wikipedia.org/wiki/XML. However, everything you need to know in
order to read and write Xml documents with the SimTK::Xml class is described 
below.

<h2>Our in-memory model of an XML document</h2>

We consider an XML document to be a tree of "Nodes". There are only four
types of nodes, which you can remember with the acronym "CUTE": Comments,
Unknowns, %Text, and Elements. Only elements can contain text and other nodes.
Elements can also have "Attributes" which are name:value pairs (not nodes).

We call the top level of this tree the "document level", represented 
by an object of class Xml. There is exactly one document-level element, known
as the "document tag". The tag word associated with it conventionally 
identifies the kind of document this is. For example, XML files produced by 
VTK begin with a document tag "<VTKFile>". The document level may also have any 
number of comment and unknown nodes appearing before and after the document tag.

We go to some pain to make sure every XML document fits the above model so that
you don't have to think about anything else. For example, if the file as read in
has multiple document-level elements, or has document-level text, we will 
enclose all the element and text nodes within document start tag "<XMLDocument>"
and end tag "</XMLDocument>" thus making it fit the description above. We call
this "canonicalizing" the document.

<h2>Reading an XML document</h2>

To read an XML document, you create an Xml object and tell it to read in the
document from a file or from a string. The document will be parsed and 
canonicalized into the in-memory model described above. Then to rummage around
in the document, you ask the Xml object for its root element (the document tag),
and check the tag word to see that it is the type of document you are 
expecting. You can check the document tag's attributes, and then process its 
contents (child nodes). Iterators are provided for running through all the
attributes, all the child nodes contained in the element, or all the child
nodes of a particular type. For a child node that is an element,
you check the tag and then pass the element to some piece of code that knows
how to deal with that kind of element and its children recursively.

Here is a complete example of reading in an Xml file "example.xml", printing
the document tag and then the types of all the top level nodes, in STL
iterator style:
@code
    Xml ex("example.xml");
    cout << "Document tag: " << ex.getDocumentTag() << endl;
    for (Xml::node_iterator p=ex.node_begin(); p != ex.node_end(); ++p)
        cout << "Node type: " << p->getNodeTypeAsString() << endl;
@endcode
Exactly one of the above nodes will have type "ElementNode", and that is the
root element or document tag.

<h2>Writing an XML document</h2>

You can add, remove, and modify nodes and attributes in a document, or create a
document from scratch. Then you can write the results in a "pretty-printed" or
compact format to a file or a string. Whenever we write an XML document, we
write it in canoncial format, regardless of how it looked when we found it.

<h2>Details about XML</h2>

This section provides detailed information about the syntax of XML files as
we accept and produce them. You won't have to know these details to read and
write XML files using the SimTK::Xml class, but you may find this helpful for
when you have to look at an XML file in a text editor.

<h3>Lexical elements</h3>

(Ignore the quote characters below; those are present so I can get this text
through Doxygen.)
  - An XML document is a string of Unicode characters; all metadata is case
    sensitive.
  - The file begins with a "declaration" tag beginning with "<?xml" and ending
    with "?>"
  - Comments look like this: "<!--" \e anything "-->"
  - The characters in an XML file represent \e markup and \e content
  - Markup consists of "tags" delimited by "<" and ">", \e attributes denoted
    by \e name="value", and character escapes delimited by "&" and ";".
  - Tags come in three flavors: \e start tags like "<word>", \e end tags like
    "</word>" and <em>empty element</em> tags like "<word/>". Tag words must
    begin with a letter or an underscore and are case sensitive; "xml" is 
    reserved; don't use it.
  - Attributes are recognized only in start tags, empty element tags, and
    declaration tags. In standard XML the value must be quoted with single or 
    double quotes, but we'll supply missing quotes if there are none.
    Attribute names are case sensitive and must be unique within a tag; but if
    we see duplicates we'll just ignore all but the last.
  - There are five pre-defined escapes: "&lt;" and "&gt;" representing "<" 
    and ">", "&amp;" for ampersand, "&apos;" for apostrophe (single quote)
    and "&quot;" for double quote.
  - There are also "numeric character reference" escapes of the form "&#nnnnn;"
    (decimal) or "&#xnnnn;" (hex), with only as many digits as needed.
  - Text set off by "<![CDATA[" and "]]>" is interpreted as a raw byte stream.
  - Tags that begin "<x" where x is not a letter or underscore and isn't one
    of the above recognized patterns will be passed through uninterpreted.
  - Anything else is Unicode text. 

<h3>File structure</h3>

An XML file contains a single \e document which consists at the top level of 
  - a declaration
  - comments
  - unknowns
  - a root element

Elements can be containers and are the basis for the tree structure
of XML files. Elements can contain:
  - attributes
  - comments
  - unknowns
  - text
  - child elements, recursively

A declaration also has attributes, but there are only three: version,
encoding, and standalone ('yes' or 'no'). Unknowns are constructs found in the
file that are not recognized; they might be errors but they are likely to
be more sophisticated uses of XML that our feeble parser doesn't understand.
Unknowns are tags where the tag word doesn't begin with a letter or
underscore and isn't one of the very few other tags we recognize, like
comments.

Here is the top-level structure we expect for a well-formed XML document, and
we will impose this structure on XML documents that don't have it. This allows
us to simplify the in-memory model as discussed above.
@code
    <?xml version="1.0" encoding="UTF-8"?>
    <!-- maybe comments and unknowns -->
    <doctag attr=value ... >
        ... contents ...
    </doctag>
    <!-- maybe comments and unknowns -->
@endcode
That is, the first line should be a declaration, most commonly exactly the
characters shown above, without the "standalone" attribute which will 
default to "yes". If we don't see a declaration when reading an XML
document, we'll assume we read the one above. Then the document should contain 
exactly one top-level (root) element representing the type of document and 
document-level attributes. The tag for the root element is 
called the "document tag"; it is not literally "doctag" but some name
that makes sense for the given document. Note that the document element is
an ordinary element so "contents" can contain text and child elements (as well
as comments and unknowns).

When reading an XML document, if it has exactly one top-level element and no
top-level text, we'll assume that is the document tag. If there is 
more than one top-level element, or we find some top-level text, we'll assume 
the document tag is missing and act as though we had seen a document tag 
"<XMLDocument>" at the beginning and "</XMLDocument>" at the end so
the document tag will be "XMLDocument". Note that this means that we will
interpret even a plain text file as a well-formed XML document:
@code
    A file consisting            <?xml version="1.0" encoding="UTF-8" ?>
    of just text         ==>     <XMLDocument>
    like this.                   A file consisting of just text like this.
                                 </XMLDocument>
@endcode
The above XML document has a single top-level element and that element
contains one Text node whose value is the original text.
**/
class SimTK_SimTKCOMMON_EXPORT Xml {
public:
    // These local classes are used to describe the contents of an XML document.
    class Attribute;
    class Node;         // This is the abstract type for any node.
    class Comment;      // These are the concrete node types.
    class Unknown;      //                  "
    class Text;         //                  "
    class Element;      //                  "

    // These provide iteration over all the attributes found in a given
    // element tag.
    class attribute_iterator;
    class const_attribute_iterator;

    // These provide iteration over all the nodes at either the Xml document
    // level or over the child nodes of a node.
    class node_iterator;
    class const_node_iterator;

    // These provide iteration over all the element nodes that are children
    // of a given element, or over the child elements that have a particular
    // tag word.
    class element_iterator;
    class const_element_iterator;

    /** The NodeType enum serves as the actual type of a node and as a filter
    for allowable node types during an iteration over nodes. We consider
    Element and Text nodes to be meaningful, while Comment and Unknown 
    nodes are meaningless junk. However, you are free to extract some meaning
    from them if you know how. **/
    enum NodeType {
        NoNode      = 0x00, ///< No nodes allowed
        ElementNode = 0x01, ///< Element node type and only-Elements filter
        TextNode    = 0x02, ///< Text node type and only-Text nodes filter
        CommentNode = 0x04, ///< Comment node type and only-Comments filter
        UnknownNode = 0x08, ///< Unknown node type and only-Unknowns filter

        NoJunkNodes = ElementNode|TextNode,    ///< Filter out meaningless nodes
        JunkNodes   = CommentNode|UnknownNode, ///< Filter out meaningful nodes
        AnyNodes    = NoJunkNodes|JunkNodes    ///< Filter allowing all nodes
    };

    /** Translate a NodeType to a human-readable string. **/
    static String getNodeTypeAsString(NodeType type);

    /** Create an empty XML Document with default declaration and default
    document tag "XMLDocument". That is, if you printed out this document
    now you would see:                                          @code
            <?xml version="1.0" encoding="UTF-8"?>
            <XMLDocument/>                                      @endcode **/
    Xml();

    /** Create a new XML document and initialize it from the contents
    of the given file name. An exception will be thrown if the file doesn't 
    exist or can't be parsed. **/
    explicit Xml(const String& pathname);

    /** Restore this document to its default-constructed state. **/
    void clear();

    /** The document type is conventionally the root element's tag; if there
    wasn't a unique root element then we will have created one with document
    tag "XMLDocument", so this will always work. **/
    const String& getDocumentTag() const;
    /** The document type is the root element's tag; this changes that tag. **/
    void setDocumentTag(const String& tag);


    /** Return a const reference to the top-level element in this Xml 
    document, known as the "document tag". The tag name is considered to
    be the type of document. This is the only top-level element; all others
    are its children and descendents. **/
    const Element& getDocumentElement() const;
    /** Return a writable reference to the top-level "document tag" element. **/
    Element& updDocumentElement();

    /** Read the contents of this Xml document from the file whose pathname
    is supplied. This first clears the current document so the new one 
    completely replaces the old one. @see readFromString() **/
    void readFromFile(const String& pathname);
    /** Write the contents of this in-memory Xml document to the file whose
    pathname is supplied. The file will be created if it doesn't exist, 
    overwritten if it does exist. **/
    void writeToFile(const String& pathname) const;

    /** Read the contents of this Xml document from the supplied string. This
    first clears the current document so the new one completely replaces the 
    old one. @see readFromFile() **/
    void readFromString(const String& xmlDocument);
    /** Alternate form that reads from a null-terminated C string (char*) 
    rather than a C++ string object. This would otherwise be implicitly 
    converted to string first which would require copying. **/
    void readFromString(const char* xmlDocument);
    /** Write the contents of this in-memory Xml document to the supplied
    string. The string cleared first so will be completely overwritten.
    Normally the output is "pretty-printed" as it is for a file, but if you
    set \a compact to true the tabs and newlines will be suppressed to make
    a more compact representation. **/
    void writeToString(String& xmlDocument, bool compact = false) const;


    /** If you want to run through this document's top-level nodes (of which 
    the root element is one), these methods provide begin and end iterators.
    By default you'll see all the nodes (types Comment, Text, Element, and
    Unknown) but you can restrict the node types that you'll see. **/
    node_iterator       node_begin(NodeType allowed=AnyNodes);
    /** Const version of node_begin(). **/
    const_node_iterator node_begin(NodeType allowed=AnyNodes) const;

    /** Any node_iterator can be set node_end() regardless of the NodeTypes it 
    allows. **/
    node_iterator       node_end();
    /** Const version of node_end(). **/
    const_node_iterator node_end() const;


    /** This is the absolute path name of the file (if any) from which this
    Xml document was read in or to which it was most recently written. **/
    String getPathname() const;

    /** Returns the Xml "version" attribute as a string; this comes from
    the "declaration" line at the beginning of an Xml document; that is the
    line that begins with "<?xml" and ends with "?>". **/
    String getXmlVersion() const;
    /** Returns the Xml "encoding" attribute as a string (from the declaration
    line at the beginning of the document). **/
    String getXmlEncoding() const;
    /** Returns the Xml "standalone" attribute as a string (from the declaration
    line at the beginning of the document); default is "yes", meaning that the
    document can be parsed correctly without any other documents. **/
    bool   getXmlIsStandalone() const;

    /** Set the Xml "version" attribute; this will be written to the 
    "declaration" line which is first in any Xml document; that is the
    line that begins with "<?xml" and ends with "?>". **/
    void setXmlVersion(const String& version);
    /** Set the Xml "encoding" attribute; this does not change the encoding,
    just what gets output to the "declaration" line. **/
    void setXmlEncoding(const String& encoding);
    /** Set the Xml "standalone" attribute; this is normally true (corresponding
    to standalone="yes") and won't appear in the declaration line in that 
    case. If you set this to false then standalone="no" will appear in the
    declaration line when it is written. **/
    void setXmlIsStandalone(bool isStandalone);

private:
    class Impl;
    const Impl& getImpl() const {assert(impl); return *impl;}
    Impl&       updImpl()       {assert(impl); return *impl;}
    Impl* impl;

friend class Node;
};

/** Output a "pretty printed" textual representation of the given XML
document to an std::ostream. **/
// Do this inline so we don't have to pass the ostream through the API.
inline std::ostream& operator<<(std::ostream& o, const Xml& xmlDocument) {
    String output;
    xmlDocument.writeToString(output);
    return o << output;
}



//------------------------------------------------------------------------------
//                              XML ATTRIBUTE
//------------------------------------------------------------------------------
/** Elements can have attributes; attributes are not Nodes. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Attribute {
public:
    /** Default constructor creates a null Attribute handle. **/
    Attribute() : tiAttr(0) {}
    /** Create a new Attribute that is not connected to any Xml document. **/
    Attribute(const String& name, const String& value);
    /** Copy constructor is shallow; that is, this handle will refer to the
    same attribute as the source. **/
    Attribute(Attribute& src) : tiAttr(src.tiAttr) {} 
    /** Copy assignment is shallow; that is, this handle will refer to the
    same attribute as the source. **/
    Attribute& operator=(Attribute& src) 
    {   if (&src!=this) {clear(); tiAttr=src.tiAttr;} return *this; }
    /** Destructor will delete the referenced attribute if this handle is its
    owner, otherwise it will be left alone. **/
    ~Attribute() {clear();}
    /** Is this handle currently holding an attribute? **/
    bool isValid() const {return tiAttr!=0;}
    /** If this is a valid attribute handle, get the name of the attribute. **/
    const String& getName() const;
    /** If this is a valid attribute handle, get the value of the attribute
    as a String, not including the quotes. **/
    const String& getValue() const;
    /** If this is a valid attribute handle, change its name. 
    @return A reference to this attribute that now has the new name. **/
    Attribute& setName(const String& name);
    /** If this is a valid attribute handle, change its value to the given
    String which should not be quoted.  
    @return A reference to this attribute that now has the new value. **/
    Attribute& setValue(const String& value);

    /** Clear this attribute handle; if the handle is the attribute's owner
    (that is, it isn't part of a document) then the attribute will be deleted,
    otherwise it is left unchanged. **/
    void clear();

    /** Comparison operators return true if the same attribute is being
    referenced or both handles are empty. Note that two different attributes
    with the same properties will not test equal by this criterion. **/
    bool operator==(const Attribute& attr) const {return tiAttr==attr.tiAttr;}
    bool operator!=(const Attribute& attr) const {return tiAttr!=attr.tiAttr;}

private:
    explicit Attribute(TiXmlAttribute* attr) {tiAttr=attr;}
    const TiXmlAttribute& getTiAttr() const {assert(tiAttr);return *tiAttr;}
    TiXmlAttribute&       updTiAttr()       {assert(tiAttr);return *tiAttr;}

    // Careful; this does not clear the handle before replacing the pointer
    // so must not be used if this could be the owner handle of an attribute
    // that hasn't ever been added to a document. It is intended for use by
    // iterators, whose contained Attributes can never be owners.
    void                  setTiAttrPtr(TiXmlAttribute* attr) {tiAttr=attr;}
    const TiXmlAttribute* getTiAttrPtr() const {return tiAttr;}
    TiXmlAttribute*       updTiAttrPtr()       {return tiAttr;}

    TiXmlAttribute* tiAttr;

friend class Xml::attribute_iterator;
friend class Xml::const_attribute_iterator;
};



//------------------------------------------------------------------------------
//                          XML ATTRIBUTE ITERATOR
//------------------------------------------------------------------------------
/** This is a bidirectional iterator suitable for moving forward or backward
within a list of Attributes within an Element, for writable access. **/
class SimTK_SimTKCOMMON_EXPORT Xml::attribute_iterator 
:   public std::iterator<std::bidirectional_iterator_tag, Xml::Attribute> {
public:
    /** Default constructor creates an iterator that compares equal to 
    attribute_end(). **/
    attribute_iterator() {}
    /** Construct this iterator to point to the same attribute as does the
    supplied Attribute handle (or attribute_end() if the handle is empty). **/
    explicit attribute_iterator(Attribute& attr) : attr(attr) {}
    /** Copy constructor takes an attribute_iterator that can be const, but that
    still allows writing to the Attribute. **/
    attribute_iterator(const attribute_iterator& src) 
    :   attr(src->updTiAttrPtr()) {}
    /** Copy assignment takes an attribute_iterator that can be const, but that
    still allows writing to the Attribute. **/
    attribute_iterator& operator=(const attribute_iterator& src) 
    {   attr.setTiAttrPtr(src->updTiAttrPtr()); return *this; }

    attribute_iterator& operator++();   // prefix
    attribute_iterator operator++(int); // postfix
    attribute_iterator& operator--();   // prefix
    attribute_iterator operator--(int); // postfix
    Attribute& operator*() {return attr;}
    Attribute* operator->() {return &attr;}
    Attribute& operator*() const {return const_cast<Attribute&>(attr);}
    Attribute* operator->() const {return const_cast<Attribute*>(&attr);}
    bool operator==(const attribute_iterator& other) const 
    {   return other.attr==attr; }
    bool operator!=(const attribute_iterator& other) const 
    {   return other.attr!=attr; }
private:
    explicit attribute_iterator(TiXmlAttribute* ap) : attr(ap) {}
    Attribute  attr;
friend class Xml::Element;
friend class Xml::const_attribute_iterator;
};



//------------------------------------------------------------------------------
//                         XML CONST ATTRIBUTE ITERATOR
//------------------------------------------------------------------------------
/** This is a bidirectional iterator suitable for moving forward or backward
within a list of Attributes within an Element, for read-only access. **/
class SimTK_SimTKCOMMON_EXPORT Xml::const_attribute_iterator 
:   public std::iterator<std::bidirectional_iterator_tag, Xml::Attribute> {
public:
    const_attribute_iterator() {}
    explicit const_attribute_iterator(const Attribute& attr) 
    :   attr(const_cast<Attribute&>(attr)) {}

    /** This is an implicit conversion from writable attribute_iterator. **/
    const_attribute_iterator(const Xml::attribute_iterator& p) 
    :   attr(const_cast<Attribute&>(p.attr)) {}

    /** Copy constructor makes this iterator refer to the same attribute
    (if any) as the source iterator. **/
    const_attribute_iterator(const const_attribute_iterator& src)
    :   attr(const_cast<Attribute&>(*src)) {}

    const_attribute_iterator& operator++();   // prefix
    const_attribute_iterator operator++(int); // postfix
    const_attribute_iterator& operator--();   // prefix
    const_attribute_iterator operator--(int); // postfix
    const Attribute& operator*() const {return attr;}
    const Attribute* operator->() const {return &attr;}
    bool operator==(const const_attribute_iterator& other) const 
    {   return other.attr==attr; }
    bool operator!=(const const_attribute_iterator& other) const 
    {   return other.attr!=attr; }
private:
    explicit const_attribute_iterator(const TiXmlAttribute* ap) 
    :   attr(const_cast<TiXmlAttribute*>(ap)) {}
    void reassign(const TiXmlAttribute* ap)
    {   attr.setTiAttrPtr(const_cast<TiXmlAttribute*>(ap)); }
    Attribute  attr;
friend class Xml::Element;
};



//------------------------------------------------------------------------------
//                               XML NODE
//------------------------------------------------------------------------------
/** Abstract handle for holding any kind of node in an XML tree. The concrete
node types are: Element, Text, Comment, and Unknown. An Element may recursively
contain a list of nodes. Elements can also have Attributes, which are are
name:value pairs that are not Nodes. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Node {
public:
    Node() : tiNode(0) {}
    explicit Node(TiXmlNode* tiNode) : tiNode(tiNode) {}

    /** Get the "value" of this node; means different things for different
    types of nodes: 
      - Comment: everything between "<!--" and "-->" including spaces
      - Unknown: everything between "<" and ">" including spaces
      - Text:    the text
      - Element: the tag name (i.e. the element's type) **/
    const String& getValue() const;

    /** Get the Xml::NodeType of this node. **/
    NodeType getNodeType() const;

    /** Get the node type as a string. **/
    String getNodeTypeAsString() const;

    bool empty(NodeType allowed=AnyNodes) const;

    // For iterating through the children
    node_iterator       node_begin(NodeType allowed=AnyNodes);
    const_node_iterator node_begin(NodeType allowed=AnyNodes) const;

    /** Any node_iterator can be set node_end() regardless of the NodeTypes it 
    allows. **/
    node_iterator       node_end();
    const_node_iterator node_end() const;

    bool operator==(const Node& other) const {return other.tiNode==tiNode;}
    bool operator!=(const Node& other) const {return other.tiNode!=tiNode;}

    void clear() {tiNode=0;}
    bool isValid() const {return tiNode != 0;}

private:
    const TiXmlNode& getTiNode() const {assert(tiNode);return *tiNode;}
    TiXmlNode&       updTiNode()       {assert(tiNode);return *tiNode;}

    // Careful: these "Ptr" methods provide raw access to the contained 
    // pointer without any cleanup or error checking. In particular, 
    // setTiNodePtr() does not attempt to delete the current contents.
    void setTiNodePtr(TiXmlNode* node) {tiNode=node;}
    const TiXmlNode* getTiNodePtr() const {return tiNode;}
    TiXmlNode*       updTiNodePtr()       {return tiNode;}

    TiXmlNode* tiNode;

friend class Xml;
friend class Xml::Impl;
friend class Xml::node_iterator;
friend class Xml::const_node_iterator;
friend class Xml::Element;
};



//------------------------------------------------------------------------------
//                          XML NODE ITERATOR
//------------------------------------------------------------------------------
/** This is a bidirectional iterator suitable for moving forward or backward
within a list of Nodes, for writable access. By default we will iterate
over all nodes but you can restrict the types at construction. **/
class SimTK_SimTKCOMMON_EXPORT Xml::node_iterator 
:   public std::iterator<std::bidirectional_iterator_tag, Xml::Node> {
public:
    explicit node_iterator(NodeType allowed=AnyNodes) 
    :   allowed(allowed) {}
    explicit node_iterator(Node& node, NodeType allowed=AnyNodes) 
    :   node(node), allowed(allowed) {}
    node_iterator& operator++();   // prefix
    node_iterator operator++(int); // postfix
    node_iterator& operator--();   // prefix
    node_iterator operator--(int); // postfix
    Node& operator*() {return node;}
    Node* operator->() {return &node;}
    const Node& operator*() const {return node;}
    const Node* operator->() const {return &node;}
    bool operator==(const node_iterator& other) const {return other.node==node;}
    bool operator!=(const node_iterator& other) const {return other.node!=node;}
private:
    explicit node_iterator(TiXmlNode* tiNode, NodeType allowed=AnyNodes) 
    :   node(tiNode), allowed(allowed) {}
    Node     node;
    NodeType allowed;

friend class Xml;
friend class Xml::Node;
friend class Xml::const_node_iterator;
};



//------------------------------------------------------------------------------
//                          XML NODE CONST ITERATOR
//------------------------------------------------------------------------------
/** This is a bidirectional iterator suitable for moving forward or backward
within a list of Nodes, for const access. By default we will iterate
over all nodes but you can restrict the types at construction. **/
class SimTK_SimTKCOMMON_EXPORT Xml::const_node_iterator 
:   public std::iterator<std::bidirectional_iterator_tag, Xml::Node> {
public:
    /** This is the default constructor which leaves the node_iterator empty, 
    and you can optionally set the type(s) of Nodes which will be iterated 
    over. **/
    explicit const_node_iterator(NodeType allowed=AnyNodes) 
    :   allowed(allowed) {}
    /** Constructor an node_iterator pointing to a given Node, and optionally 
    set the type(s) of Nodes which will be iterated over. **/
    explicit const_node_iterator(const Node& node, NodeType allowed=AnyNodes) 
    :   node(node), allowed(allowed) {}

    /** This is an implicit conversion from writable node_iterator. **/
    const_node_iterator(const Xml::node_iterator& p) 
    :   node(p.node), allowed(p.allowed) {}

    const_node_iterator& operator++();   // prefix
    const_node_iterator operator++(int); // postfix
    const_node_iterator& operator--();   // prefix
    const_node_iterator operator--(int); // postfix
    const Node& operator*() const {return node;}
    const Node* operator->() const {return &node;}
    bool operator==(const const_node_iterator& other) const 
    {   return other.node==node; }
    bool operator!=(const const_node_iterator& other) const 
    {   return other.node!=node; }
private:
    explicit const_node_iterator(const TiXmlNode* tiNode, NodeType allowed=AnyNodes) 
    :   node(const_cast<TiXmlNode*>(tiNode)), allowed(allowed) {}
    Node node;
    NodeType allowed;
friend class Xml;
friend class Xml::Node;
};



//------------------------------------------------------------------------------
//                               XML ELEMENT
//------------------------------------------------------------------------------
/** An element has (1) a tag, (2) a map of name:value pairs called attributes, 
and (3) a list of nodes. The nodes can be comments, text, child elements, and
unknowns. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Element : public Xml::Node {
public:
/** Create an empty Element handle; this is suitable only for holding
references to other Elements. **/
Element() : Node() {}

/** Create an Element that uses the given tag word but is not yet part of
any XML document. Initially the Element will be empty so would print as 
"<tagWord/>", but you can add contents afterwards so that it will print as
"<tagWord>contents</tagWord>", where contents may be text and/or child
elements. **/
Element(const String& tagWord);

/** Append text to the contents of this element. If the element is 
currently empty, or if the last child node contained in the element is not
a Text node, then this will result in a new Text node with the given 
contents added to the end of the list of child nodes. Otherwise the new
text is simply appended to the last text node in the element. 
@return A handle to the Text node to which the new \a text was 
appended. **/
Text appendText(const String& text);

/** Insert text before the location indicated by a node_iterator, which
must point to a node currently in this Element or be node_end() in which
case the text is appended as with appendText(). If the indicated node
is a Text node, then then given text is prepended to that node. If not,
but the previous node is a Text node, then the given text is appended to
the previous node. Otherwise, a new Text node is created and inserted 
prior to the one indicated by the node_iterator. 
@return A handle to the Text node into which the new \a text was 
inserted. **/
Text insertText(const const_node_iterator& node, const String& text);


element_iterator            element_begin(const String& tag="");
const_element_iterator      element_begin(const String& tag="") const;

element_iterator            element_end();
const_element_iterator      element_end() const;

attribute_iterator          attribute_begin();
const_attribute_iterator    attribute_begin() const;

attribute_iterator          attribute_end();
const_attribute_iterator    attribute_end() const;

attribute_iterator find_attribute(const String& name) {
    for (attribute_iterator p = attribute_begin();
            p != attribute_end(); ++p)
        if (p->getName() == name) return p;
    return attribute_end();
}

const_attribute_iterator find_attribute(const String& name) const {
    for (const_attribute_iterator p = attribute_begin();
            p != attribute_end(); ++p)
        if (p->getName() == name) return p;
    return attribute_end();
}    

///** Return an array containing Attribute handles referencing all the
//attributes of this element. Attributes are returned in the order that
//they appear in the element tag. Attribute names within a tag are unique;
//if the source document had repeated attribute names only the last one
//to appear is retained and that's the only one we'll find here. **/
//Array_<Attribute> findAllAttributes();
//
///** Return an array containing Element handles referencing all the
//immediate child elements contained in this element, or all the child 
//elements of a particular type (that is, with a given tag word). Elements 
//are returned in the order they are seen in the document. **/
//Array_<Element> findAllElements(const String& tag="");
//
///** Return an array containing Node handles referencing all the
//immediate child nodes contained in this element, or all the child 
//nodes of a particular type or types. Nodes are returned in the order they 
//are seen in the document. **/
//Array_<Node> findAllNodes(NodeType type=AnyNodes);

/** The element tag word can be considered the "type" of the element. **/
const String& getElementTag() const;
/** Change the tag word that is used to bracket this element. **/
void setElementTag(const String& tag);
/** Return true if this element has an attribute of this name. **/
bool hasAttribute(const String& name) const;
/** Return true if this element has a child element with this tag. **/
bool hasElement(const String& tag) const;

/** Determine whether this element qualifies as a "text element", defined
as an element containing zero or one Text nodes and no child elements. 
You can treat a text element as you would an attribute -- it can be viewed
as having a single value, which is just the value of its lone Text node
(or a null string if it doesn't have any text). **/
bool isTextElement() const;

/** Get the text value of this text element. An error will be thrown if 
this is not a "text element", defined as one with zero or one Text nodes
and no child elements. A null string is returned if there is no Text 
node; otherwise, we return the value of the Text node. 
@see isTextElement() **/
const String& getElementText() const;

/** Assuming this is a "text element", convert its text value to the type
of the template argument T. It is an error if the text can not be converted,
in its entirety, to a single object of type T. (But note that type T may
be a container of some sort, like a Vector or Array.) 
@tparam T   A type that can be read from a stream using the ">>" operator.
**/
template <class T> T getElementAs() const 
{   return getElementText().convertTo<T>(); }

/** Alternate form of getElementAs() that avoids unnecessary copying and
heap allocation for reading in large container objects. **/
template <class T> void getElementAs(T& out) const 
{   return getElementText().convertTo<T>(out); }

/** Obtain a reference to a particular attribute of this element; an error
will be thrown if no such attribute is present. **/
Attribute getRequiredAttribute(const String& name) const;

/** Get the text value of an attribute (that is, its value as a string) and 
throw an error if that attribute is not present. **/
const String& getRequiredAttributeText(const String& name) const
{   return getRequiredAttribute(name).getValue(); }

/** Convert the text value of a required attribute to the type
of the template argument T. It is an error if the text can not be converted,
in its entirety, to a single object of type T. (But note that type T may
be a container of some sort, like a Vec3.) 
@tparam T   A type that can be read from a stream using the ">>" operator.
**/
template <class T> T getRequiredAttributeAs
   (const String& name) const
{   return getRequiredAttributeText(name).convertTo<T>(); }

/** Get the text value of an attribute (that is, its value as a string) if 
the attribute is present in this element, otherwise return a supplied 
default value. **/
String getOptionalAttributeText
   (const String& name, const String& def="") const
{   const_attribute_iterator p = find_attribute(name);
    return p==attribute_end() ? def : p->getValue(); }

/** Convert the text value of an optional attribute, if present, to the type
of the template argument T. It is an error if the text can not be converted,
in its entirety, to a single object of type T. (But note that type T may
be a container of some sort, like a Vec3.) If the attribute is not present,
then return a supplied default value of type T.
@tparam     T    A type that can be read from a stream with operator ">>".
@param[in]  name The name of the optional attribute.
@param[in]  def  The value of type T to return if the attribute is missing.
@return The value of attribute \a name if it is present, otherwise a copy
of the supplied default value \a def. **/
template <class T> T getOptionalAttributeAs
   (const String& name, const T& def) const
{   const_attribute_iterator p = find_attribute(name);
    return p==attribute_end() ? def : p->getValue().convertTo<T>(); }

/** Get the text value of a child text element that \e must be present in 
this element. The child is identified by its tag; if there is more than one
this refers to the first one. Then the element is expected to contain either
zero or one Text nodes; if none we'll return a null string, otherwise 
the value of the Text node. Thus an element like "<tag>stuff</tag> will
have the text value "stuff". An error will be thrown if either the element
is not found or it is not a "text element". **/
const String& getRequiredElementText(const String& tag) const
{   return getRequiredElement(tag).getElementText(); }

/** Get the text value of a child text element that \e may be present in
this element, otherwise return a default string. If the child element is 
found, it must be a "text element" as defined above. **/
String getOptionalElementText(const String& tag, const String& def="") const;

/** Convert the text value of a required child text element to the type
of the template argument T. It is an error if the element is present but is
not a text element, or if the text cannot be converted,
in its entirety, to a single object of type T. (But note that type T may
be a container of some sort, like a Vector or Array.) 
@tparam     T   A type that can be read from a stream using the ">>" operator.
@param[in]  tag The tag of the required child text element.
@return The value of the text element, converted to an object of type T. **/
template <class T> T  getRequiredElementAs(const String& tag) const
{   return getRequiredElementText(tag).convertTo<T>(); }

/** Convert the text value of an optional child text element, if present, to
the type of the template argument T. It is an error if the child element is
present but is not a text element, or if the text cannot be 
converted, in its entirety, to a single object of type T. (But note that 
type T may be a container of some sort, like a Vector or Array.) If the 
child element is not present, then return a supplied default value of type T.
@tparam     T    A type that can be read from a stream with operator ">>".
@param[in]  tag  The tag of the optional child element.
@param[in]  def  The value of type T to return if child element is missing.
@return The value of element \a tag if it is present, otherwise a copy
of the supplied default value \a def. **/
template <class T> T 
    getOptionalElementAs(const String& tag, const T& def) const
{   const Element opt(getOptionalElement(tag));
    return opt.isValid() ? opt.getElementText().convertTo<T>() : def; }

/** Get a reference to a child element that \e must be present in this 
element. The child is identified by its tag; if there is more than one
only the first one is returned. If you want to see all children with this
tag, use getAllChildElements() or use an element_iterator. **/
Element getRequiredElement(const String& tag) const;

/** Get a reference to a child element that \e may be present in this 
element; otherwise return an invalid Element handle. Test using the
Element's isValid() method. **/
Element getOptionalElement(const String& tag) const;

/** Test whether a given Node is an element node. **/
static bool isA(const Node&);
/** Recast a Node to a const Element, throwing an error if the Node is not
actually an element node. @see isA() **/
static const Element& getAs(const Node& node);
/** Recast a writable Node to a writable Element, throwing an error if the
Node is not actually an element node. @see isA() **/
static Element& updAs(Node& node);

private:
explicit Element(TiXmlElement* tiElt) 
:   Node(reinterpret_cast<TiXmlNode*>(tiElt)) {}

TiXmlElement& updTiElement() 
{   return reinterpret_cast<TiXmlElement&>(updTiNode()); } 
const TiXmlElement& getTiElement() const
{   return reinterpret_cast<const TiXmlElement&>(getTiNode()); }

// Careful: these "Ptr" methods provide raw access to the contained 
// pointer without any cleanup or error checking. In particular, 
// setTiElementPtr() does not attempt to delete the current contents.
const TiXmlElement* getTiElementPtr() const 
{   return reinterpret_cast<const TiXmlElement*>(getTiNodePtr()); }
TiXmlElement*       updTiElementPtr()
{   return reinterpret_cast<TiXmlElement*>(updTiNodePtr()); }
void                setTiElementPtr(TiXmlElement* elt)
{   setTiNodePtr(reinterpret_cast<TiXmlNode*>(elt)); }

friend class Xml::element_iterator;
friend class Xml::const_element_iterator;
};


//------------------------------------------------------------------------------
//                          XML ELEMENT ITERATOR
//------------------------------------------------------------------------------
/** This is a bidirectional iterator suitable for moving forward or backward
within a list of Nodes, for writable access. By default we will iterate
over all nodes but you can restrict the type at construction. **/
class SimTK_SimTKCOMMON_EXPORT Xml::element_iterator 
:   public std::iterator<std::bidirectional_iterator_tag, Xml::Element> {
public:
    explicit element_iterator(const String& tag="") : tag(tag) {}
    explicit element_iterator(Element& elt, const String& tag="") 
    :   elt(elt), tag(tag) {}

    /** Copy constructor takes an element_iterator that can be const, but that
    still allows writing to the Element. **/
    element_iterator(const element_iterator& src) 
    :   elt(src->updTiElementPtr()), tag(src.tag) {}
    /** Copy assignment takes an element_iterator that can be const, but that
    still allows writing to the Element. **/
    element_iterator& operator=(const element_iterator& src) 
    {   elt.setTiElementPtr(src->updTiElementPtr()); return *this; }


    element_iterator& operator++();   // prefix
    element_iterator operator++(int); // postfix
    element_iterator& operator--();   // prefix
    element_iterator operator--(int); // postfix
    Element& operator*() {return elt;}
    Element* operator->() {return &elt;}
    Element& operator*() const {return const_cast<Element&>(elt);}
    Element* operator->() const {return const_cast<Element*>(&elt);}
    bool operator==(const element_iterator& other) const {return other.elt==elt;}
    bool operator!=(const element_iterator& other) const {return other.elt!=elt;}
private:
    explicit element_iterator(TiXmlElement* tiElt, const String& tag="") 
    :   elt(tiElt), tag(tag) {}
    Element  elt;
    String   tag;
friend class Xml;
friend class Xml::Element;
friend class Xml::const_element_iterator;
};



//------------------------------------------------------------------------------
//                          XML CONST ELEMENT ITERATOR
//------------------------------------------------------------------------------
/** This is a bidirectional iterator suitable for moving forward or backward
within a list of Nodes, for const access. By default we will iterate
over all elements but you can restrict the type at construction. **/
class SimTK_SimTKCOMMON_EXPORT Xml::const_element_iterator 
:   public std::iterator<std::bidirectional_iterator_tag, Xml::Element> {
public:
    /** This is the default constructor which leaves the element_iterator empty, and
    you can optionally set the type of Element which will be iterated over. **/
    explicit const_element_iterator(const String& tag="") : tag(tag) {}
    /** Constructor an element_iterator pointing to a given Node, and optionally 
    set the type(s) of Nodes which will be iterated over. **/
    explicit const_element_iterator(const Element& elt, const String& tag="") 
    :   elt(elt), tag(tag) {}

    /** This is an implicit conversion from writable element_iterator. **/
    const_element_iterator(const Xml::element_iterator& p) 
    :   elt(p.elt), tag(p.tag) {}

    const_element_iterator& operator++();   // prefix
    const_element_iterator operator++(int); // postfix
    const_element_iterator& operator--();   // prefix
    const_element_iterator operator--(int); // postfix
    const Element& operator*() const {return elt;}
    const Element* operator->() const {return &elt;}
    bool operator==(const const_element_iterator& other) const 
    {   return other.elt==elt; }
    bool operator!=(const const_element_iterator& other) const 
    {   return other.elt!=elt; }
private:
    explicit const_element_iterator
       (const TiXmlElement* tiElt, const String& tag="") 
    :   elt(const_cast<TiXmlElement*>(tiElt)), tag(tag) {}
    void reassign(const TiXmlElement* ep)
    {   elt.setTiElementPtr(const_cast<TiXmlElement*>(ep)); }

    Element  elt;
    String   tag;

friend class Xml;
friend class Xml::Element;
};


//inline Array_<Xml::Node> Xml::Element::
//findAllNodes(Xml::NodeType type)
//{   return Array_<Xml::Node>(node_begin(type), node_end()); }
//
//inline Array_<Xml::Element> Xml::Element::
//findAllElements(const String& tag)
//{   return Array_<Element>(element_begin(tag), element_end()); }
//
//inline Array_<Xml::Attribute> Xml::Element::
//findAllAttributes()
//{   return Array_<Attribute>(attribute_begin(), attribute_end()); }

/** This is the "leaf" content of an element. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Text : public Xml::Node {
public:
    /** Create an empty Text handle, suitable only for holding references
    to other Text nodes. **/
    Text() : Node() {}
private:
};

/** A comment contains only uninterpreted text. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Comment : public Xml::Node {
public:
private:
};

/** This is something we don't understand but can carry around. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Unknown : public Xml::Node {
public:
private:
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_XML_H_


