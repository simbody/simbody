#ifndef SimTK_SimTKCOMMON_XML_H_
#define SimTK_SimTKCOMMON_XML_H_

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
 * Contributors: Peter Eastman                                                *
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
class TiXmlText;
class TiXmlComment;
class TiXmlUnknown;
    
/** This class provides a minimalist capability for reading and writing XML 
documents, as files or strings. This is based with gratitude on the excellent 
open source XML parser TinyXML (http://www.grinninglizard.com/tinyxml/). Note 
that this is a <em>non-validating</em> parser, meaning it deals only with the 
XML file itself and not with a Document Type Definition (DTD), XML Schema, or 
any other description of the XML file's expected contents. Instead, the 
structure of your code that uses this class encodes the expected structure 
and contents of the XML document.

Our in-memory model of an XML document is simplified even further than 
TinyXML's. There a lot to know about XML; you could start here: 
http://en.wikipedia.org/wiki/XML. However, everything you need to know in
order to read and write XML documents with the SimTK::Xml class is described 
below.

Much of the detailed documention is in the class Xml::Element; be sure to look
there as well as at this overview.

<h2>Our in-memory model of an XML document</h2>

We consider an XML document to be a tree of "Nodes". There are only four
types of nodes: Comments, Unknowns, %Text, and Elements. Only Elements can 
contain Text and other nodes, including recursively child Element nodes. 
Elements can also have "Attributes" which are name:value pairs (not nodes).

The XML document as a whole is represented by an object of class Xml::Document.
The Xml::Document object directly contains a short list of nodes, consisting 
only of Comments, Unknowns, and a single %Element called the "root element". The
tag word associated with the root element is called the "root tag" and 
conventionally identifies the kind of document this is. For example, XML files
produced by VTK begin with a root tag "<VTKFile>".

We go to some pain to make sure every Xml::Document fits the above model so 
that you don't have to think about anything else. For example, if the file as 
read in has multiple root-level elements, or has document-level text, we will 
enclose all the element and text nodes within document start tag "<_Root>"
and end tag "</_Root>" thus making it fit the description above. We call
this "canonicalizing" the document.

<h3>%Value Elements</h3>

Element nodes can be classified into "value elements" and "compound
elements". A value element is a "leaf" element (no child elements) that 
contains at most one Text node. For example, a document might contain value
elements like these:
@code
    <name>John Doe</name>
    <rating>7.2</rating>
    <winnings currency=euro>3429</winnings>
    <preferences/>
    <vector>1.2 -4 2e-3</vector>
@endcode
All of these have a unique value so it makes sense to talk about "the" value
of these elements (the empty "preferences" element has a null value). These 
are very common in XML documents, and the Xml::Element class makes them very 
easy to work with. For example, if Element elt is the "<vector>" element
from the example, you could retrieve its value as a Vec3 like this:
@code
    Vec3 v = elt.getValueAs<Vec3>(); 
@endcode
This would automatically throw an error if the element wasn't a value element 
or if its value didn't have the right format to convert to a Vec3.

Note that it is okay for a value element to have attributes; those are ignored
in determining the element's value. Any element that is not a value element is 
a "compound element", meaning it has either child elements and/or more than 
one %Text node.

<h2>Reading an XML document</h2>

To read an XML document, you create an Xml::Document object and tell it to 
read in the document from a file or from a string. The document will be parsed
and canonicalized into the in-memory model described above. Then to rummage 
around in the document, you ask the Xml::Document object for its root element,
and check the root tag to see that it is the type of document you are 
expecting. You can check the root element's attributes, and then process its 
contents (child nodes). Iterators are provided for running through all the
attributes, all the child nodes contained in the element, or all the child
nodes of a particular type. For a child node that is an element,
you check the tag and then pass the element to some piece of code that knows
how to deal with that kind of element and its children recursively.

Here is a complete example of reading in an Xml file "example.xml", printing
the root tag and then the types of all the document-level nodes, in STL
iterator style:
@code
    Xml::Document doc("example.xml");
    cout << "Root tag: " << ex.getRootTag() << endl;
    for (Xml::node_iterator p=doc.node_begin(); p != doc.node_end(); ++p)
        cout << "Node type: " << p->getNodeTypeAsString() << endl;
@endcode
Exactly one of the above nodes will have type "ElementNode"; that is the
root element. To print out the types of nodes contained in the root element,
you could write:
@code
    Xml::Element root = ex.getRootElement();
    for (Xml::node_iterator p=root.node_begin(); p != root.node_end(); ++p)
        cout << "Node type: " << p->getNodeTypeAsString() << endl;
@endcode
    
(Some confessions: despite appearances, "Xml" is not a namespace, it is a 
class with the other Xml classes being internal classes of Xml. An object of 
type Xml is an XML document; the name Xml::Document is a typedef synonymous 
with Xml.)

<h2>Writing an XML document</h2>

You can insert, remove, and modify nodes and attributes in a document, or 
create a document from scratch. Then you can write the results in a 
"pretty-printed" or compact format to a file or a string; for pretty-printing
you can override the default indentation string (four spaces). Whenever we 
write an XML document, we write it in canoncial format, regardless of how it 
looked when we found it.

At the document level, you can only insert Comment and Unknown nodes. Text and
Element nodes can be inserted only at the root element level and below.

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
  - %Text set off by "<![CDATA[" and "]]>" is interpreted as a raw byte stream.
  - Tags that begin "<x" where x is not a letter or underscore and isn't one
    of the above recognized patterns will be passed through uninterpreted.
  - Anything else is Unicode text. 

<h3>File structure</h3>

An XML file contains a single \e document which consists at the top level of 
  - a declaration
  - comments and unknowns
  - a root element
  - more comments and unknowns

Elements can be containers of other nodes and are thus the basis for the tree
structure of XML files. Elements can contain:
  - comments
  - unknowns
  - text
  - child elements, recursively
  - attributes

A declaration (see below) also has attributes, but there are only three: 
version, encoding, and standalone ('yes' or 'no'). Unknowns are constructs 
found in the file that are not recognized; they might be errors but they are 
likely to be more sophisticated uses of XML that our feeble parser doesn't 
understand. Unknowns are tags where the tag word doesn't begin with a letter or
underscore and isn't one of the very few other tags we recognize, like
comments. As an example, a DTD tag like this would come through as an Unknown 
node here:
@code
    <!DOCTYPE note SYSTEM "Note.dtd">
@endcode

Here is the top-level structure we expect of a well-formed XML document, and
we will impose this structure on XML documents that don't have it. This allows
us to simplify the in-memory model as discussed above.
@code
    <?xml version="1.0" encoding="UTF-8"?>
    <!-- maybe comments and unknowns -->
    <roottag attr=value ... >
        ... contents ...
    </roottag>
    <!-- maybe comments and unknowns -->
@endcode
That is, the first line should be a declaration, most commonly exactly the
characters shown above, without the "standalone" attribute which will 
default to "yes". If we don't see a declaration when reading an XML
document, we'll assume we read the one above. Then the document should contain 
exactly one root element representing the type of document and 
document-level attributes. The tag for the root element is not literally 
"roottag" but some name that makes sense for the given document. Note that the 
root element is an ordinary element so "contents" can contain text and child 
elements (as well as comments and unknowns).

When reading an XML document, if it has exactly one document-level element and 
no document-level text, we'll take the document as-is. If there is 
more than one document-level element, or we find some document-level text, 
we'll assume that the root element is missing and act as though we had seen a 
root element "<_Root>" at the beginning and "</_Root>" at the end so
the root tag will be "_Root". Note that this means that we will
interpret even a plain text file as a well-formed XML document:
@code
    A file consisting            <?xml version="1.0" encoding="UTF-8" ?>
    of just text         ==>     <_Root>
    like this.                   A file consisting of just text like this.
                                 </_Root>
@endcode
The above XML document has a single document-level element and that element
contains one Text node whose value is the original text. 

@see Xml::Element, Xml::Node **/

//------------------------------------------------------------------------------
//                                   XML
//------------------------------------------------------------------------------
class SimTK_SimTKCOMMON_EXPORT Xml {
public:

// These local classes are used to describe the contents of an XML document.
class Attribute;
class Node;         // This is the abstract handle type for any node.
class Comment;      // These are the concrete node types.
class Unknown;      //                  "
class Text;         //                  "
class Element;      //                  "

/** This typedef allows Xml::Document to be used as the type of the document
which is more conventional than using just Xml, and provides future 
compatibility should we decide to upgrade Xml::Document to a class. **/
typedef Xml Document;

// This provides iteration over all the attributes found in a given element.
class attribute_iterator;

// This provides iteration over all the nodes, or nodes of a certain type,
// at either the Xml document level or over the child nodes of an element.
class node_iterator;

// This provides iteration over all the element nodes that are children
// of a given element, or over the subset of those child elements that has
// a particular tag word.
class element_iterator;

/** The NodeType enum serves as the actual type of a node and as a filter
for allowable node types during an iteration over nodes. We consider
Element and Text nodes to be meaningful, while Comment and Unknown 
nodes are meaningless junk. However, you are free to extract some meaning
from them if you know how. In particular, DTD nodes end up as Unknown. **/
enum NodeType {
    NoNode      = 0x00, ///< Type of empty Node handle, or null filter
    ElementNode = 0x01, ///< Element node type and only-Elements filter
    TextNode    = 0x02, ///< Text node type and only-Text nodes filter
    CommentNode = 0x04, ///< Comment node type and only-Comments filter
    UnknownNode = 0x08, ///< Unknown node type and only-Unknowns filter

    NoJunkNodes = ElementNode|TextNode,    ///< Filter out meaningless nodes
    JunkNodes   = CommentNode|UnknownNode, ///< Filter out meaningful nodes
    AnyNodes    = NoJunkNodes|JunkNodes    ///< Allow all nodes
};

/** Translate a NodeType to a human-readable string. **/
static String getNodeTypeAsString(NodeType type);

/**@name                         Construction
You can start with an empty Xml::Document or initialize it from a file. **/
/*@{*/

/** Create an empty XML Document with default declaration and default root
element with tag "_Root".\ (You should invoke this as Xml::Document() instead
of just Xml().) 

If you were to print out this document now you would see:                                          
@code
    <?xml version="1.0" encoding="UTF-8"?>
    <_Root />                                      
@endcode **/
Xml();

/** Create a new XML document and initialize it from the contents of the given
file name.\ (You should invoke this as Xml::Document() instead
of just Xml().)  

An exception will be thrown if the file doesn't exist or can't be 
parsed. @see readFromFile(), readFromString() **/
explicit Xml(const String& pathname);

/** Copy constructor makes a deep copy of the entire source document; nothing
is shared between the source and the copy. **/
Xml(const Xml::Document& source);

/** Copy assignment frees all heap space associated with the current
Xml::Document and then makes a deep copy of the source document; nothing is
shared between the source and the copy. **/
Xml::Document& operator=(const Xml::Document& souce);

/** The destructor cleans up all heap space associated with this document. **/
~Xml();

/** Restore this document to its default-constructed state. **/
void clear();
/*@}*/

/**@name                    Serializing and I/O
These methods deal with conversion to and from the in-memory representation
of the XML document from and to files and strings. **/
/*@{*/
/** Read the contents of this Xml::Document from the file whose pathname
is supplied. This first clears the current document so the new one 
completely replaces the old one. @see readFromString() **/
void readFromFile(const String& pathname);
/** Write the contents of this in-memory Xml::Document to the file whose
pathname is supplied. The file will be created if it doesn't exist, 
overwritten if it does exist. The file will be "pretty-printed" using the
current indent string. @see setIndentString(), writeToString() **/
void writeToFile(const String& pathname) const;
/** Read the contents of this Xml::Document from the supplied string. This
first clears the current document so the new one completely replaces the 
old one. @see readFromFile() **/
void readFromString(const String& xmlDocument);
/** Alternate form that reads from a null-terminated C string (char*) 
rather than a C++ string object. This would otherwise be implicitly 
converted to string first which would require copying. **/
void readFromString(const char* xmlDocument);
/** Write the contents of this in-memory Xml::Document to the supplied
string. The string cleared first so will be completely overwritten.
Normally the output is "pretty-printed" as it is for a file, but if you
set \a compact to true the tabs and newlines will be suppressed to make
a more compact representation. @see setIndentString(), writeToFile() **/
void writeToString(String& xmlDocument, bool compact = false) const;
/** Set the string to be used for indentation when we produce a
"pretty-printed" serialized form of this document.\ The default is
to use four spaces for each level of indentation. 
@see writeToFile(), writeToString(), getIndentString() **/
void setIndentString(const String& indent);
/** Return the current value of the indent string.\ The default is
four spaces. @see setIndentString() **/
const String& getIndentString() const;

/** Set global mode to control whether white space is preserved or condensed 
down to a single space (affects all subsequent document reads; not document
specific). The default is to condense. **/
static void setXmlCondenseWhiteSpace(bool shouldCondense);
/** Return the current setting of the global "condense white space" option. 
Note that this option affects all Xml reads; it is not document specific. **/
static bool isXmlWhiteSpaceCondensed();
/*@}*/


/**@name                Top-level node manipulation
These methods provide access to the top-level nodes, that is, those that are
directly owned by the Xml::Document. Comment and Unknown nodes are allowed 
anywhere at the top level, but Text nodes are not allowed and there is just
one distinguished Element node, the root element. If you want to add Text
or Element nodes, add them to the root element rather than at the document
level. **/
/*@{*/

/** Return an Element handle referencing the top-level element in this 
Xml::Document, known as the "root element". The tag word of this element is 
usually the type of document. This is the \e only top-level element; all others
are its children and descendents. Once you have the root Element handle, you 
can also use any of the Element methods to manipulate it. If you need a 
node_iterator that refers to the root element (perhaps to use one of the
top-level insert methods), use node_begin() with a NodeType filter:
@code
    Xml::node_iterator rootp = Xml::node_begin(Xml::ElementNode);
@endcode
That works since there is only one element at this level. **/
Element getRootElement();

/** Shortcut for getting the tag word of the root element which is usually
the document type. This is the same as getRootElement().getElementTag(). **/
const String& getRootTag() const;
/** Shortcut for changing the tag word of the root element which is usually
the document type. This is the same as getRootElement().setElementTag(tag). **/
void setRootTag(const String& tag);

/** Insert a top-level Comment or Unknown node just \e after the location 
indicated by the node_iterator, or at the end of the list if the iterator is 
node_end(). The iterator must refer to a top-level node. The Xml::Document 
takes over ownership of the Node which must be a Comment or Unknown node and
must have been an orphan. The supplied Node handle will retain a reference 
to the node within the document and can still be used to make changes, but
will no longer by an orphan. **/
void insertTopLevelNodeAfter (const node_iterator& afterThis, 
                              Node                 insertThis);
/** Insert a top-level Comment or Unknown node just \e before the location 
indicated by the node_iterator. See insertTopLevelNodeAfter() for details. **/
void insertTopLevelNodeBefore(const node_iterator& beforeThis, 
                              Node                 insertThis);
/** Delete the indicated top-level node, which must not be the root element,
and must not be node_end(). That is, it must be a top-level Comment or
Unknown node which will be removed from the Xml::Document and deleted. The
iterator is invalid after this call; be sure not to use it again. Also, there 
must not be any handles referencing the now-deleted node. **/
void eraseTopLevelNode(const node_iterator& deleteThis);
/** Remove the indicated top-level node from the document, returning it as an
orphan rather than erasing it. The node must not be the root element,
and must not be node_end(). That is, it must be a top-level Comment or
Unknown node which will be removed from the Xml::Document and returned as
an orphan Node. The iterator is invalid after this call; be sure not to use it 
again. **/
Node removeTopLevelNode(const node_iterator& removeThis);
/*@}*/


/**@name       Iteration through top-level nodes (rarely used)
If you want to run through this document's top-level nodes (of which the
"root element" is one), these methods provide begin and end 
iterators. By default you'll see all the nodes (types Comment, Unknown, 
and the lone top-level Element) but you can restrict the node types that 
you'll see via the NodeType mask. Iteration is rarely used at this top level 
since you almost never care about about the Comment and Unknown nodes here and
you can get to the root element directly using getRootElement().
@see getRootElement() **/
/*@{*/
/** Obtain an iterator to all the top-level nodes or a subset restricted via
the \a allowed NodeType mask. **/
node_iterator       node_begin(NodeType allowed=AnyNodes);

/** This node_end() iterator indicates the end of a sequence of nodes regardless
of the NodeType restriction on the iterator being used. **/
node_iterator       node_end() const;
/*@}*/

/**@name           XML Declaration attributes (rarely used)
These methods deal with the mysterious XML "declaration" line that comes at the
beginning of every XML document; that is the line that begins with "<?xml" 
and ends with "?>". There are at most three of these attributes and they have
well-defined names that are always the same (default values shown):
  - \e version = "1.0": to what version of the XML standard does this document 
    adhere?
  - \e encoding = "UTF-8": what Unicode encoding is used to represent the 
    character in this document? Typically this is UTF-8, an 8-bit encoding in 
    which the first 128 codes match standard ASCII but where other characters 
    are represented in variable-length multibyte sequences.
  - \e standalone = "yes": can this document be correctly parsed without 
    consulting other documents?

You can examine and change these attributes with the methods in this section,
however unless you really know what you're doing you should just leave the
declaration alone; you'll get reasonable behavior automatically. **/
/*@{*/
/** Returns the Xml "version" attribute as a string (from the declaration
line at the beginning of the document). **/
String getXmlVersion() const;
/** Returns the Xml "encoding" attribute as a string (from the declaration
line at the beginning of the document). **/
String getXmlEncoding() const;
/** Returns the Xml "standalone" attribute as a bool (from the declaration
line at the beginning of the document); default is true ("yes" in a file), 
meaning that the document can be parsed correctly without any other documents. 
We won't include "standalone" in the declaration line for any Xml documents 
we generate unless the value is false ("no" in a file). **/
bool getXmlIsStandalone() const;

/** Set the Xml "version" attribute; this will be written to the 
"declaration" line which is first in any Xml document. **/
void setXmlVersion(const String& version);
/** Set the Xml "encoding" attribute; this doesn't affect the in-memory
representation but can affect how the document gets written out. **/
void setXmlEncoding(const String& encoding);
/** Set the Xml "standalone" attribute; this is normally true (corresponding
to standalone="yes") and won't appear in the declaration line in that 
case when we write it out. If you set this to false then standalone="no" 
will appear in the declaration line when it is written. **/
void setXmlIsStandalone(bool isStandalone);
/*@}*/

//------------------------------------------------------------------------------
                                  private:
friend class Node;

class Impl; // a private, local class Xml::Impl
const Impl& getImpl() const {assert(impl); return *impl;}
Impl&       updImpl()       {assert(impl); return *impl;}

Xml& unconst() const {return *const_cast<Xml*>(this);}

Impl*       impl; // This is the lone data member.
};

/** Output a "pretty printed" textual representation of the given Xml::Document
to an std::ostream, using the document's current indent string for
formatting. @see Xml::setIndentString() 
@relates Xml **/
// Do this inline so we don't have to pass the ostream through the API.
inline std::ostream& operator<<(std::ostream& o, const Xml::Document& doc) {
    String output;
    doc.writeToString(output);
    return o << output;
}



//------------------------------------------------------------------------------
//                              XML ATTRIBUTE
//------------------------------------------------------------------------------
/** Elements can have attributes, which are name="value" pairs that appear
within the element start tag in an XML document; this class represents the
in-memory representation of one of those attributes and can be used to examine
or modify the name or value. Attribute names within an element tag are unique.
**/
class SimTK_SimTKCOMMON_EXPORT Xml::Attribute {
public:
/** Default constructor creates a null Attribute handle. **/
Attribute() : tiAttr(0) {}
/** Create a new orphan Attribute, that is, an Attribute that is not 
owned by any Xml Element. **/
Attribute(const String& name, const String& value);
/** Copy constructor is shallow; that is, this handle will refer to the
same attribute as the source. Note that this handle will provide write
access to the underlying attribute, even if the source was const. **/
Attribute(const Attribute& src) : tiAttr(src.tiAttr) {} 
/** Copy assignment is shallow; the handle is first cleared and then will 
refer to the same attribute as the source. Note that this handle will 
provide write access to the underlying attribute even if the source handle
was const. @see clear() **/
Attribute& operator=(const Attribute& src) 
{   if (&src!=this) {clear(); tiAttr=src.tiAttr;} return *this; }
/** The Attribute handle destructor does not recover heap space so if you create
orphan attributes and then don't put them in a document there will be a memory 
leak unless you explicitly destruct them first with clearOrphan(). **/
~Attribute() {clear();}
/** Is this handle currently holding an attribute? **/
bool isValid() const {return tiAttr!=0;}
/** Return true if this Attribute is an orphan, meaning that it is not 
empty, but is not owned by any element or top-level document. This is 
typically an Attribute object that has just been constructed, or one that
has been cloned from another Attribute. **/
bool isOrphan() const;
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

/** This method restores the Attribute handle to its default-constructed
state but does not recover any heap space; use clearOrphan() if you know
this attribute was never put into a document. **/
void clear();
/** This method explictly frees the heap space for an orphan attribute that
was created but never inserted into a document. It is an error to call this
if the attribute is in a document. **/
void clearOrphan();

/** Serialize this attribute to the given String. The output will be as it
would appear in an XML file, i.e. name="value" or name='value' with special
characters output as suitable entities. If you don't want it that way, use
getName() and getValue() instead which return the raw strings. **/
void writeToString(String& out) const;

/** Comparison operators return true if the same attribute is being
referenced or both handles are empty. Note that two different attributes
with the same properties will not test equal by this criterion. **/
bool operator==(const Attribute& attr) const {return tiAttr==attr.tiAttr;}
bool operator!=(const Attribute& attr) const {return tiAttr!=attr.tiAttr;}

//------------------------------------------------------------------------------
                                  private:
friend class Xml::attribute_iterator;
friend class Xml::Element;

explicit Attribute(TiXmlAttribute* attr) {tiAttr=attr;}
const TiXmlAttribute& getTiAttr() const {assert(tiAttr);return *tiAttr;}
TiXmlAttribute&       updTiAttr()       {assert(tiAttr);return *tiAttr;}

// Careful; this does not clear the handle before replacing the pointer
// so should not be used if this could be the owner handle of an attribute
// that hasn't ever been added to a document. It is intended for use by
// iterators, whose contained Attributes can never be owners.
void                  setTiAttrPtr(TiXmlAttribute* attr) {tiAttr=attr;}
const TiXmlAttribute* getTiAttrPtr() const {return tiAttr;}
TiXmlAttribute*       updTiAttrPtr()       {return tiAttr;}

Attribute& unconst() const {return *const_cast<Attribute*>(this);}

TiXmlAttribute* tiAttr; // this is the lone data member
};

/** Output a textual representation of the given Attribute to an std::ostream.
This will be in the form the Attribute would appear in an XML file; that is,
name="value" or name='value' with entity substituion for odd characters, 
without surrounding blanks. @relates Xml::Attribute **/
// Do this inline so we don't have to pass the ostream through the API.
inline std::ostream& operator<<(std::ostream& o, const Xml::Attribute& attr) {
    String output;
    attr.writeToString(output);
    return o << output;
}



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
/** An iterator destructor never deletes the object to which it refers. **/
~attribute_iterator() {attr.setTiAttrPtr(0);}
/** Copy assignment takes an attribute_iterator that can be const, but that
still allows writing to the Attribute. **/
attribute_iterator& operator=(const attribute_iterator& src) 
{   attr.setTiAttrPtr(src->updTiAttrPtr()); return *this; }
/** Prefix increment operator advances the iterator to the next attribute (or
attribute_end() if it was already at the last attribute) and
returns a reference to the now-incremented iterator. **/
attribute_iterator& operator++();   // prefix
/** Postfix increment operator advances the iterator to the next attribute (or
attribute_end() if it was already at the last attribute) and
returns an iterator still referencing the previous one. **/
attribute_iterator operator++(int); // postfix
/** Prefix decrement operator moves the iterator to the previous attribute (or
attribute_end() if it was already at the first attribute) and
returns a reference to the now-decremented iterator. **/
attribute_iterator& operator--();   // prefix
/** Postfix decrement operator moves the iterator to the previous attribute (or
attribute_end() if it was already at the first attribute) and
returns an iterator still referencing the original one. **/
attribute_iterator operator--(int); // postfix

// It's the iterator that's const in these next two methods; it still points
// to a non-const object just like a char* const p.

/** Return a writable reference to the Attribute referenced by this
iterator; the handle will be invalid if the iterator was attribute_end(). **/
Attribute& operator*() const {return const_cast<Attribute&>(attr);}
/** Return a writable pointer to the Attribute referenced by this
iterator; the pointer will never be null but the handle it points to will be 
invalid if the iterator was attribute_end(). **/
Attribute* operator->() const {return const_cast<Attribute*>(&attr);}
/** Comparison return true only if both iterators refer to the same in-memory
attribute or both are at attribute_end(); iterators referencing two different 
attributes that happen to have identical properties will not test equal by 
these criteria. **/
bool operator==(const attribute_iterator& other) const 
{   return other.attr==attr; }
/** Uses same criteria as operator==(). **/
bool operator!=(const attribute_iterator& other) const 
{   return other.attr!=attr; }

//------------------------------------------------------------------------------
                                  private:
friend class Xml::Element;

explicit attribute_iterator(TiXmlAttribute* ap) : attr(ap) {}

Attribute       attr;   // the lone data member
};



//------------------------------------------------------------------------------
//                               XML NODE
//------------------------------------------------------------------------------
/** Abstract handle for holding any kind of node in an XML tree. The concrete
node handle types derived from Node are: Comment, Unknown, Text, and Element. 
Only an Element node may contain other nodes. 

A node may be classified by who owns it. There are three possibilities:
  - Top-level node: The node belongs to the top-level Xml document and does
    not have a parent node.
  - Child node: The node belongs to an element, which may be the root element
    or any lower-level element. The element that owns it is its "parent".
  - Orphan node: The node is not yet part of any Xml document and does not
    belong to an element. In that case the Node handle serves as the owner and
    the node does not have a parent node.

A Node handle may also be empty, meaning it refers to no node at all so there
is nothing to own.

Top-level nodes can only be Comment nodes, Unknown nodes, or the lone root
Element node. Child nodes and orphans can be of Element and Text type also.
Normally orphans exist only briefly during the time a new node is constructed
and the time it is adopted by some element (usually in the same constructor)
so you can ignore them for the most part. But if you must keep orphan nodes
around, be aware that they must be referenced by only one handle at a time to 
avoid ownership conflicts. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Node {
public:

/**@name                Construction and destruction
These methods are mostly used by the derived node classes; Nodes are
not generally created directly in user code. **/
/*@{*/

/** Create an empty Node handle that can be used to hold a reference to any
kind of Node. **/
Node() : tiNode(0) {}
/** Copy constructor is shallow; that is, this handle will refer to the
same node as the source. Note that this handle will provide write
access to the underlying node, even if the source was const. **/
Node(const Node& src) : tiNode(src.tiNode) {} 
/** Copy assignment is shallow; the handle is first cleared and then will 
refer to the same node as the source. Note that this handle will 
provide write access to the underlying node even if the source handle
was const. @see clear() **/
Node& operator=(const Node& src) 
{   if (&src!=this) {clear(); tiNode=src.tiNode;} return *this; }
/** The clone() method makes a deep copy of this Node and its children and
returns a new orphan Node with the same contents; ordinary assignment and
copy construction is shallow. **/
Node clone() const;
/** The Node handle destructor does not recover heap space so if you create
orphan nodes and then don't put them in a document there will be a memory 
leak unless you explicitly destruct them first with clearOrphan(). **/
~Node() {clear();}
/** This method restores the Node handle to its default-constructed
state but does not recover any heap space; use clearOrphan() if you know
this node was never put into a document. **/
void clear();
/** This method explictly frees the heap space for an orphan node that
was created but never inserted into a document. It is an error to call this
if the node is in a document. **/
void clearOrphan();
/*@}*/

/**@name                  Node classification
You can find out what concrete type of node this abstract Node handle is 
referring to (if any), who owns the node, and if it is owned by a parent
element you can get access to the parent. **/
/*@{*/

/** Get the Xml::NodeType of this node. If this Node handle is empty, the
returned NodeType will be "Xml::NoNode". **/
NodeType getNodeType() const;

/** Get the Node type as a string; an empty handle returns "NoNode". **/
String getNodeTypeAsString() const;

/** Return true if this Node handle is referencing some node, false if the
Node handle is empty. **/
bool isValid() const {return tiNode != 0;}

/** Return true if this Node is owned by the top-level Xml document, false
if the Node is owned by an Element or is an orphan, or if the Node handle
is empty. **/
bool isTopLevelNode() const;

/** Return true if this Node is an orphan, meaning that it is not empty, but
is not owned by any element or top-level document. This is typically a Node
object that has just been constructed, or one that has been cloned from another
Node. **/
bool isOrphan() const;

/** Return true if this node has a parent, i.e. it is owned by an element; 
the root element and other top-level nodes are owned by the document and thus
do not have a parent. @see getParentElement() **/
bool hasParentElement() const;

/** Return a handle referencing this node's parent if it has one, otherwise
throws an error; check first with hasParentElement() if you aren't sure. **/
Element getParentElement();
/*@}*/

/**@name               Access to node contents
Usually contents inspection is handled at the concrete node class level;
here we can only provide information for which you don't need to know what
kind of node this is. **/
/*@{*/

/** Return a text value associated with this Node (\e not including its child
nodes if any); the behavior depends on the NodeType. This is a convenience 
that saves downcasting a generic Node to a concrete type when all you want to 
do is dump out the text. It is not particularly useful for Element nodes. Here
is what you get for each type of node:
  - Comment: everything between "<!--" and "-->"
  - Unknown: everything between "<" and ">"
  - Text:    the text
  - Element: the element's tag word (\e not the element's value)
  - None:    (i.e., an empty handle) throw an error. **/
const String& getNodeText() const;

/** Serialize this node (and everything it contains) to the given String.
The output will be "pretty printed" and terminated with a newline unless you
specify \a compact = true in which case indents and newlines will be
suppressed. Pretty printing uses the containing Xml::Document's indent string,
if this Node is in a document, otherwise the default of four spaces for each
indent level is used. **/
void writeToString(String& out, bool compact=false) const;
/*@}*/

/** Comparing Nodes for equality means asking if the two Node handles are
referring to exactly the same object; two different nodes that happen to have
the same properties will not test equal by this criteria. **/
bool operator==(const Node& other) const {return other.tiNode==tiNode;}
/** Inequality test using same criteria as operator==(). **/
bool operator!=(const Node& other) const {return other.tiNode!=tiNode;}


//------------------------------------------------------------------------------
                                 protected:
/** @cond **/ // don't let Doxygen see these
explicit Node(TiXmlNode* tiNode) : tiNode(tiNode) {}

const TiXmlNode& getTiNode() const {assert(tiNode);return *tiNode;}
TiXmlNode&       updTiNode()       {assert(tiNode);return *tiNode;}

// Careful: these "Ptr" methods provide raw access to the contained 
// pointer without any cleanup or error checking. In particular, 
// setTiNodePtr() does not attempt to delete the current contents.
void setTiNodePtr(TiXmlNode* node) {tiNode=node;}
const TiXmlNode* getTiNodePtr() const {return tiNode;}
TiXmlNode*       updTiNodePtr()       {return tiNode;}
/** @endcond **/

//------------------------------------------------------------------------------
                                  private:
friend class Xml;
friend class Xml::Impl;
friend class Xml::node_iterator;
friend class Xml::Comment;
friend class Xml::Unknown;
friend class Xml::Text;
friend class Xml::Element;

Node& unconst() const {return *const_cast<Node*>(this);}

TiXmlNode*      tiNode; // the lone data member
};

/** Output a "pretty printed" textual representation of the given XML
node (and all its contents) to an std::ostream. Pretty printing uses the 
indent string from the Node's containing Xml::Document, if this Node is in a 
document, otherwise the default of four spaces for each indent level is used.
@relates Xml::Node **/
// Do this inline so we don't have to pass the ostream through the API.
inline std::ostream& operator<<(std::ostream& o, const Xml::Node& xmlNode) {
    String output;
    xmlNode.writeToString(output);
    return o << output;
}



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

/** Copy constructor takes a node_iterator that can be const, but that
still allows writing to the Node. **/
node_iterator(const node_iterator& src) 
:   node(*src), allowed(src.allowed) {}
/** An iterator destructor never deletes the object to which it refers. **/
~node_iterator() {node.setTiNodePtr(0);}
/** Copy assignment takes an node_iterator that can be const, but that
still allows writing to the Node. **/
node_iterator& operator=(const node_iterator& src) 
{   node = *src; allowed = src.allowed; return *this; }

node_iterator& operator++();   // prefix
node_iterator operator++(int); // postfix
node_iterator& operator--();   // prefix
node_iterator operator--(int); // postfix
Node& operator*() {return node;}
Node* operator->() {return &node;}
// It's the iterator that's const; it still points to a non-const object
// just like a char* const p.
Node& operator*()  const {return const_cast<Node&>(node);}
Node* operator->() const {return const_cast<Node*>(&node);}
bool operator==(const node_iterator& other) const {return other.node==node;}
bool operator!=(const node_iterator& other) const {return other.node!=node;}

//------------------------------------------------------------------------------
                                 protected:
explicit node_iterator(TiXmlNode* tiNode, NodeType allowed=AnyNodes) 
:   node(tiNode), allowed(allowed) {}
void reassign(TiXmlNode* tiNode)
{   node.setTiNodePtr(tiNode); }

//------------------------------------------------------------------------------
                                  private:
friend class Xml;
friend class Xml::Node;
friend class Xml::Element;
friend class Xml::element_iterator;

Node            node;       // data members
NodeType        allowed;
};



//------------------------------------------------------------------------------
//                          XML ELEMENT ITERATOR
//------------------------------------------------------------------------------
/** This is a bidirectional iterator suitable for moving forward or backward
within a list of Element nodes, for writable access. By default we will iterate
over all elements in a list but you can restrict the element_iterator at
construction to find only elements with a particular tag. **/
class SimTK_SimTKCOMMON_EXPORT Xml::element_iterator
:   public Xml::node_iterator {
public:

/** This is the default constructor which leaves the element_iterator empty, and
you can optionally set the type of Element which will be iterated over. **/
explicit element_iterator(const String& tag="") 
:   node_iterator(ElementNode), tag(tag) {}
/** Constructor an element_iterator pointing to a given Element, and optionally 
set the type of Element which will be iterated over. **/
inline explicit element_iterator(Element& elt, const String& tag=""); // below

/** Copy constructor takes an element_iterator that can be const, but that
still allows writing to the Element. **/
element_iterator(const element_iterator& src) 
:   node_iterator(src), tag(src.tag) {}

/** Copy assignment takes an element_iterator that can be const, but that
still allows writing to the Element. **/
element_iterator& operator=(const element_iterator& src) 
{   upcast()=src; tag = src.tag; return *this; }

element_iterator& operator++();   // prefix
element_iterator operator++(int); // postfix
element_iterator& operator--();   // prefix
element_iterator operator--(int); // postfix
inline Element& operator*() const; // below
inline Element* operator->() const; // below

bool operator==(const element_iterator& other) const 
{   return other.upcast()==upcast();}
bool operator!=(const element_iterator& other) const 
{   return other.upcast()!=upcast();}

//------------------------------------------------------------------------------
                                   private:
friend class Xml;
friend class Xml::Element;

explicit element_iterator(TiXmlElement* tiElt, const String& tag="") 
:   node_iterator((TiXmlNode*)tiElt, ElementNode), tag(tag) {}
void reassign(TiXmlElement* tiElt)
{   upcast().reassign((TiXmlNode*)tiElt); }

const node_iterator& upcast() const 
{   return *static_cast<const node_iterator*>(this); }
node_iterator& upcast() 
{   return *static_cast<node_iterator*>(this); }

String          tag;    // lone data member
};




//------------------------------------------------------------------------------
//                               XML ELEMENT
//------------------------------------------------------------------------------
/** An element has (1) a tagword, (2) a map of (name,value) pairs called 
attributes, and (3) a list of child nodes. The tag word, which begins with an 
underscore or a letter, can serve as either the type or the name of the element
depending on context. The nodes can be comments, unknowns, text, and child 
elements (recursively). It is common for "leaf" elements (elements with no 
child elements) to be supplied simply for their values, for example mass might
be provided via an element "<mass> 29.3 </mass>". We call such elements "value
elements" since they have a uniquely identifiable value similar to that of 
attributes. %Value elements have no more than one text node. They may have 
attributes, and may also have comment and unknown nodes but they cannot have 
any child elements. This class provides a special set of methods for dealing 
with value nodes very conveniently; they will fail if you attempt to use them 
on an element that is not a value element. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Element : public Xml::Node {
public:

/**@name                  Construction and destruction
As discussed elsewhere, elements come in two varieties: value elements and
compound elements. New value elements can be created easily since they are
essentially just a name,value pair. Compound elements require a series of 
method calls to create the Element node and then add child nodes to it. In
either case you may want to add attributes also. **/
/*@{*/

/** Create an empty Element handle; this is suitable only for holding
references to other Elements. **/
Element() : Node() {}

/** Create a value element that uses the given tag word but is not yet part of
any XML document, and optionally give it an inital value. Note that although
you provide the initial value as a string, you can access it as any type T to
which that string can be converted, using the getValueAs<T>() templatized
method.

If no initial value is provided, then the element will be empty so would print 
as "<tagWord />". If you provide a value (say "contents") here or add one 
later, it will print as "<tagWord>contents</tagWord>". In general you can add
child elements and other node types with subsequent method calls; that would 
change this element from a value element to a compound element. 
@see getValue(), updValue(), setValue() **/
explicit Element(const String& tagWord, const String& value="");

/** Create a new value element and set its initial value to the text 
equivalent of any type T for which a conversion construction String(T) is 
allowed (generally any type for which a stream insertion operator<<() 
exists). 
@see getValueAs<T>(), setValueAs<T>()**/
template <class T>
Element(const String& tagWord, const T& value)
{   new(this) Element(tagWord, String(value)); }

/** The clone() method makes a deep copy of this Element and its children and
returns a new orphan Element with the same contents; ordinary assignment and
copy construction are shallow. **/
Element clone() const;

/** Get the element tag word. This may represent the name or type of the 
element depending on context. **/
const String& getElementTag() const;
/** Change the tag word that is used to bracket this element. **/
void setElementTag(const String& tag);

/** Insert a node into the list of this Element's children, just before the
node pointed to by the supplied iterator (or at the end if the iterator
is node_end()). The iterator must refer to a node that is a child of this
Element. This Element takes over ownership of the node which must 
not already have a parent. **/
void insertNodeBefore(const node_iterator& pos, Node node);
/** Insert a node into the list of this Element's children, just after the
node pointed to by the supplied iterator (or at the end if the iterator
is node_end()). The iterator must refer to a node that is a child of this
Element. This Element takes over ownership of the node which must 
not already have a parent. **/
void insertNodeAfter(const node_iterator& pos, Node node);
/** Delete the indicated node, which must be a child of this element,
and must not be node_end(). The node will be removed from this element
and deleted. The iterator is invalid after this call; be sure not to use it 
again. Also, there must not be any handles referencing the now-deleted 
node. **/
void eraseNode(const node_iterator& deleteThis);
/** Remove the indicated node from this element without erasing it, returning
it as an orphan Node. The node must be a child of this element, and must not be 
node_end(). The node will be removed from this element and returned as an 
orphan. The iterator is invalid after this call; be sure not to use it 
again. **/
Node removeNode(const node_iterator& removeThis);
/*@}*/


/**@name                      Value elements
As described elsewhere, value elements are those that have no child elements
and only a single Text node, whose contents can be considered as the element's
value. Methods in this section allow you to work conveniently with value
elements, getting direct access to the value string or interpreting it as
some other type. You can easily modify the value by obtaining a writable
refence to the String object that holds it. We provide methods for working
with this element's value (if it is a value element) and with an element's
children's values (if this element is compound). **/
/*@{*/

/** Determine whether this element qualifies as a "value element", defined
as an element containing zero or one Text nodes and no child elements. 
You can treat a value element as you would an attribute -- it can be viewed
as having a single value, which is just the value of its lone Text node
(or a null string if it doesn't have any text). **/
bool isValueElement() const;

/** Get the text value of this value element. An error will be thrown if 
this is not a "value element". See the comments for this class for the
definition of a "value element". 
@note This does not return the same text as the base class method
Node::getNodeText() does in the case of an element node; that returns the
element tag word not its contents.
@see isValueElement(), setValue(), updValue() **/
const String& getValue() const;

/** Obtain a writable reference to the String containing the value of this
value element. An error will be thrown if this is not a value element. If the
element was initially empty and didn't contain a Text node, one will be added
to it here with a null-string value so that we can return a reference to it.
@see isValueElement(), getValue() **/
String& updValue();

/** Set the text value of this value element. An error will be thrown if 
this is not a value element. If the element was initially empty and didn't 
contain a Text node, one will be added to it here so that we have a place to
hold the \a value.
@see isValueElement(), setValueAs<T>() **/
void setValue(const String& value);

/** Set the value of this value element to the text equivalent of any type T
for which a conversion construction String(T) is allowed (generally any
type for which a stream insertion operator<<() exists). **/
template <class T>
void setValueAs(const T& value) 
{   setValue(String(value)); }

/** Assuming this is a "value element", convert its text value to the type
of the template argument T. It is an error if the text can not be converted,
in its entirety, to a single object of type T. (But note that type T may
be a container of some sort, like a Vector or Array.) 
@tparam T   A type that can be read from a stream using the ">>" operator.
**/
template <class T> T getValueAs() const 
{   T out; convertStringTo(getValue(),out); return out;}

/** Alternate form of getValueAs() that avoids unnecessary copying and
heap allocation for reading in large container objects. **/
template <class T> void getValueAs(T& out) const 
{   convertStringTo(getValue(),out); }

/** Get the text value of a child value element that \e must be present in 
this element. The child is identified by its tag; if there is more than one
this refers to the first one. Then the element is expected to contain either
zero or one Text nodes; if none we'll return a null string, otherwise 
the value of the Text node. Thus an element like "<tag>stuff</tag>" will
have the value "stuff". An error will be thrown if either the element
is not found or it is not a "value element". **/
const String& 
getRequiredElementValue(const String& tag) const
{   return unconst().getRequiredElement(tag).getValue(); }

/** Get the text value of a child value element that \e may be present in
this element, otherwise return a default string. If the child element is 
found, it must be a "value element" as defined above. **/
String 
getOptionalElementValue(const String& tag, const String& def="") const
{   const Element opt(unconst().getOptionalElement(tag));
    return opt.isValid() ? opt.getValue() : def; }

/** Convert the text value of a required child value element to the type
of the template argument T. It is an error if the element is present but is
not a value element, or if the text cannot be converted,
in its entirety, to a single object of type T. (But note that type T may
be a container of some sort, like a Vector or Array.) 
@tparam     T   A type that can be read from a stream using the ">>" operator.
@param[in]  tag The tag of the required child text element.
@return The value of the text element, converted to an object of type T. **/
template <class T> T  
getRequiredElementValueAs(const String& tag) const
{   T out; convertStringTo(unconst().getRequiredElementValue(tag), out); 
    return out; }

/** Convert the text value of an optional child value element, if present, to
the type of the template argument T. It is an error if the child element is
present but is not a value element, or if the text cannot be 
converted, in its entirety, to a single object of type T. (But note that 
type T may be a container of some sort, like a Vector or Array.) If the 
child element is not present, then return a supplied default value of type T.
@tparam     T    A type that can be read from a stream with operator ">>".
@param[in]  tag  The tag of the optional child element.
@param[in]  def  The value of type T to return if child element is missing.
@return The value of element \a tag if it is present, otherwise a copy
of the supplied default value \a def. **/
template <class T> T 
getOptionalElementValueAs(const String& tag, const T& def) const
{   const Element opt(unconst().getOptionalElement(tag));
    if (!opt.isValid()) return def;
    T out; convertStringTo(opt.getValue(), out); return out; }
/*@}*/


/**@name                         Attributes
You can add, modify, and remove element attributes with the methods in this
section. You can work directly with individual attributes by name, or you
can iterate through the list of attributes. **/
/*@{*/
/** Return true if this element has an attribute of this name. **/
bool hasAttribute(const String& name) const;

/** Set the value of an attribute of this element, creating a new one if
this is a new attribute name otherwise modifying an existing one. **/
void setAttributeValue(const String& name, const String& value);

/** Erase an attribute of this element if it exists, otherwise do nothing.
If you need to know if the attribute exists, use hasAttribute(). There is
no removeAttribute() that orphans an existing Attribute, but you can easily 
recreate one with the same name and value. **/
void eraseAttribute(const String& name);

/** Get the value of an attribute as a string and throw an error if that 
attribute is not present. **/
const String& 
getRequiredAttributeValue(const String& name) const
{   return unconst().getRequiredAttribute(name).getValue(); }

/** Convert the text value of a required attribute to the type
of the template argument T. It is an error if the text can not be converted,
in its entirety, to a single object of type T. (But note that type T may
be a container of some sort, like a Vec3.) 
@tparam T   A type that can be read from a stream using the ">>" operator.
**/
template <class T> T 
getRequiredAttributeValueAs(const String& name) const
{   T out; convertStringTo(getRequiredAttributeValue(name),out); return out; }

/** Get the value of an attribute as a string if the attribute is present in 
this element, otherwise return a supplied default value.
@param[in]  name The name of the optional attribute.
@param[in]  def  The string to return if the attribute is missing.
@return The value of attribute \a name if it is present, otherwise a copy
of the supplied default string \a def. **/
String 
getOptionalAttributeValue(const String& name, const String& def="") const
{   Attribute attr = unconst().getOptionalAttribute(name);
    if (!attr.isValid()) return def;
    return attr.getValue(); }

/** Convert the value of an optional attribute, if present, from a string to 
the type of the template argument T. It is an error if the text can not be 
converted, in its entirety, to a single object of type T. (But note that type 
T may be a container of some sort, like a Vec3.) If the attribute is not 
present, then return a supplied default value of type T.
@tparam     T    A type that can be read from a stream with operator ">>".
@param[in]  name The name of the optional attribute.
@param[in]  def  The value of type T to return if the attribute is missing.
@return The value of attribute \a name if it is present, otherwise a copy
of the supplied default value \a def. **/
template <class T> T 
getOptionalAttributeValueAs(const String& name, const T& def) const
{   Attribute attr = unconst().getOptionalAttribute(name);
    if (!attr.isValid()) return def;
    T out; convertStringTo(attr.getValue(), out); return out; }

/** Obtain an Attribute handle referencing a particular attribute of this 
element; an error will be thrown if no such attribute is present. **/
Attribute getRequiredAttribute(const String& name);

/** Obtain an Attribute handle referencing a particular attribute of this 
element specified by name, or an empty handle if no such attribute is 
present. **/
Attribute getOptionalAttribute(const String& name);

/** Return an array containing Attribute handles referencing all the
attributes of this element. Attributes are returned in the order that
they appear in the element tag. Attribute names within a tag are unique;
if the source document had repeated attribute names only the last one
to appear is retained and that's the only one we'll find here. This is
just a shortcut for @code
    Array_<Attribute>(attribute_begin(), attribute_end());
@endcode **/
Array_<Attribute> getAllAttributes()
{   return Array_<Attribute>(attribute_begin(), attribute_end()); }


/** For iterating through all the attributes of this element. If there are no 
attributes then the returned attribute_iterator tests equal to 
attribute_end(). **/
attribute_iterator attribute_begin();
/** This attribute_end() iterator indicates the end of a sequence of 
attributes. **/
attribute_iterator attribute_end() const;
/*@}*/

/**@name                    Compound elements
Many elements contain child nodes, including other elements. When there is
just a single child Text node and no child elements, we call the element a
"value element" and it is easiest to work with using the methods in the 
"Value elements" section. When there are child elements and/or multiple Text
nodes, the element is called a "compound element" and you need a way to 
iterate and recurse through its contents. The methods in this section support
looking through all contained nodes, nodes of specified types, element nodes,
or element nodes with a specified tags. You can obtain handles to child
Nodes or Elements and then iterate through those recursively. **/
/*@{*/

/** Return true if this element has a child element with this tag. **/
bool hasElement(const String& tag) const;
/** See if this element has any child nodes, or any child nodes of the type(s)
allowed by the NodeType filter if one is supplied. **/
bool hasNode(NodeType allowed=AnyNodes) const;

/** Get a reference to a child element that \e must be present in this 
element. The child is identified by its tag; if there is more than one
only the first one is returned. If you want to see all children with this
tag, use getAllElements() or use an element_iterator. **/
Element getRequiredElement(const String& tag);

/** Get a reference to a child element that \e may be present in this 
element; otherwise return an invalid Element handle. Test using the
Element's isValid() method. **/
Element getOptionalElement(const String& tag);

/** Return an array containing Element handles referencing all the
immediate child elements contained in this element, or all the child 
elements of a particular type (that is, with a given tag word). Elements 
are returned in the order they are seen in the document. This is just a
shortcut for @code
    Array_<Element>(element_begin(tag), element_end());
@endcode **/
Array_<Element> getAllElements(const String& tag="")
{   return Array_<Element>(element_begin(tag), element_end()); }

/** Return an array containing Node handles referencing all the
immediate child nodes contained in this element, or all the child 
nodes of a particular type or types. Nodes are returned in the order they 
are seen in the document. This is just a shortcut for @code
    Array_<Node>(node_begin(allowed), node_end());
@endcode **/
Array_<Node> getAllNodes(NodeType allowed=AnyNodes)
{   return Array_<Node>(node_begin(allowed), node_end()); }

/** For iterating through the immediate child elements of this element, or the 
child elements that have the indicated tag if one is supplied. If there are no 
children with the \a allowed tag then the returned element_iterator tests 
equal to element_end(). **/
element_iterator element_begin(const String& tag="");
/** This element_end() iterator indicates the end of any sequence of elements 
regardless of the tag restriction on the iterator being used. **/
element_iterator element_end() const;

/** For iterating through the immediate child nodes of this element, or the 
child nodes of the type(s) allowed by the NodeType filter if one is 
supplied. If there are no children of the \a allowed types then the returned
node_iterator tests equal to node_end(). **/
node_iterator node_begin(NodeType allowed=AnyNodes);
/** This node_end() iterator indicates the end of any sequence of nodes 
regardless of the NodeType restriction on the iterator being used. **/
node_iterator node_end() const;
/*@}*/

/**@name             Conversion to Element from Node
If you have a handle to a Node, such as would be returned by a node_iterator,
you can check whether that Node is an Element and if so cast it to one. **/
/*@{*/
/** Test whether a given Node is an element node. **/
static bool isA(const Node&);
/** Recast a Node to a const Element, throwing an error if the Node is not
actually an element node. @see isA() **/
static const Element& getAs(const Node& node);
/** Recast a writable Node to a writable Element, throwing an error if the
Node is not actually an element node. @see isA() **/
static Element& getAs(Node& node);
/*@}*/

//------------------------------------------------------------------------------
                                  private:
friend class Xml::Node;
friend class Xml::element_iterator;

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

Element& unconst() const {return *const_cast<Element*>(this);}


// no data members; see Node
};



// A few element_iterator inline definitions had to wait for Element to be
// defined.
inline Xml::element_iterator::element_iterator
   (Xml::Element& elt, const String& tag) 
:   Xml::node_iterator(elt, Xml::ElementNode), tag(tag) {}
inline Xml::Element& Xml::element_iterator::operator*() const 
{   return Element::getAs(*upcast());}
inline Xml::Element* Xml::element_iterator::operator->() const 
{   return &Element::getAs(*upcast());}




//------------------------------------------------------------------------------
//                               XML TEXT NODE
//------------------------------------------------------------------------------
/** This is the "leaf" content of an element. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Text : public Xml::Node {
public:
/** Create an empty Text node handle, suitable only for holding references
to other Text nodes. **/
Text() : Node() {}

/** Create a new Text node with the given text; the node is not yet owned by
any XML document. **/
explicit Text(const String& text);

/** The clone() method makes a deep copy of this Text node and returns a new 
orphan Text node with the same contents; ordinary assignment and copy 
construction are shallow. **/
Text clone() const;

/** Obtain a const reference to the String holding the value of this Text
**/
const String& getText() const;
/** Obtain a writable reference to the String holding the value of this Text
node; this can be used to alter the value. **/
String& updText();

/**@name              Conversion to Text from Node
If you have a handle to a Node, such as would be returned by a node_iterator,
you can check whether that Node is a Text node and if so cast it to one. **/
/*@{*/
/** Test whether a given Node is an Text node. **/
static bool isA(const Node&);
/** Recast a Node to a const Text node, throwing an error if the Node is not
actually a Text node. @see isA() **/
static const Text& getAs(const Node& node);
/** Recast a writable Node to a writable Text node, throwing an error if the
Node is not actually a Text node. @see isA() **/
static Text& getAs(Node& node);
/*@}*/

//------------------------------------------------------------------------------
                                   private:
// no data members; see Node

explicit Text(TiXmlText* tiText) 
:   Node(reinterpret_cast<TiXmlNode*>(tiText)) {}
};



//------------------------------------------------------------------------------
//                             XML COMMENT NODE
//------------------------------------------------------------------------------
/** A comment contains only uninterpreted text. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Comment : public Xml::Node {
public:
/** Create an empty Comment node handle, suitable only for holding references
to other Comment nodes. **/
Comment() : Node() {}

/** Create a new Comment node with the given text; the node is not yet owned by
any XML document. Don't include the comment delimiters "<!--" and "-->" in the 
text; those will be added automatically if the document is serialized to a
file or string. **/
explicit Comment(const String& text);

/** The clone() method makes a deep copy of this Comment node and returns a new 
orphan Comment node with the same contents; ordinary assignment and copy 
construction are shallow. **/
Comment clone() const;

/**@name              Conversion to Comment from Node
If you have a handle to a Node, such as would be returned by a node_iterator,
you can check whether that Node is a Comment node and if so cast it to one. **/
/*@{*/
/** Test whether a given Node is an Comment node. **/
static bool isA(const Node&);
/** Recast a Node to a const Comment, throwing an error if the Node is not
actually an Comment node. @see isA() **/
static const Comment& getAs(const Node& node);
/** Recast a writable Node to a writable Comment, throwing an error if the
Node is not actually an Comment node. @see isA() **/
static Comment& getAs(Node& node);
/*@}*/

//------------------------------------------------------------------------------
                                   private:
// no data members; see Node

explicit Comment(TiXmlComment* tiComment) 
:   Node(reinterpret_cast<TiXmlNode*>(tiComment)) {}
};



//------------------------------------------------------------------------------
//                             XML UNKNOWN NODE
//------------------------------------------------------------------------------
/** This is something we don't understand but can carry around. **/
class SimTK_SimTKCOMMON_EXPORT Xml::Unknown : public Xml::Node {
public:
/** Create an empty Unknown node handle, suitable only for holding references
to other Unknown nodes. **/
Unknown() : Node() {}

/** Create a new Unknown node with the given contents; the node is not yet 
owned by any XML document. Don't include the tag delimiters "<" and ">" in the 
contents; those will be added automatically if the document is serialized to a
file or string. That is, if you want "<!SOMETHING blah blah>", the contents
you provide should be "!SOMETHING blah blah". **/
explicit Unknown(const String& contents);

/** Create a new Unknown node and append it to the list of nodes that are
children of the given Element. The Element becomes the owner of the new
Unknown node although the handle retains a reference to it. **/
Unknown(Element& element, const String& contents)
{   new(this) Unknown(contents); 
    element.insertNodeBefore(element.node_end(), *this); }

/** The clone() method makes a deep copy of this Unknown node and returns a new 
orphan Unknown node with the same contents; ordinary assignment and copy 
construction are shallow. **/
Unknown clone() const;

/** Obtain the contents of this Unknown node. This is everything that would
be between the "<" and ">" in the XML document. **/
const String& getContents() const;
/** Change the contents of this Unknown node. This is everything that would
be between the "<" and ">" in the XML document. **/
void setContents(const String& contents);

/**@name              Conversion to Unknown from Node
If you have a handle to a Node, such as would be returned by a node_iterator,
you can check whether that Node is an Unknown node and if so cast it to one. **/
/*@{*/
/** Test whether a given Node is an Unknown node. **/
static bool isA(const Node&);
/** Recast a Node to a const Unknown, throwing an error if the Node is not
actually an Unknown node. @see isA() **/
static const Unknown& getAs(const Node& node);
/** Recast a writable Node to a writable Unknown, throwing an error if the
Node is not actually an Unknown node. @see isA() **/
static Unknown& getAs(Node& node);
/*@}*/

//------------------------------------------------------------------------------
                                   private:
// no data members; see Node

explicit Unknown(TiXmlUnknown* tiUnknown) 
:   Node(reinterpret_cast<TiXmlNode*>(tiUnknown)) {}
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_XML_H_


