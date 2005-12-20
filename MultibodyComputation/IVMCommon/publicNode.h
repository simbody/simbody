#ifndef __publicNode_hh__
#define __publicNode_hh__

#include <accessor.h>

//
// description of hinge+node info for external use
//

class PublicNode {
public:
    ACCESSOR(dim,Dim,int);
    ACCESSOR(startIndex,StartIndex,int);
    ACCESSOR(parentAtom,ParentAtom,int);

    String       type()  { return type_; }
    CDSList<int> atoms() {return atoms_;}     // note that atoms[0] is hinge atom

    void setType (const String       &v) { type_=v; }
    void setAtoms(const CDSList<int> &v) { atoms_=v;}
private:
    String       type_;
    CDSList<int> atoms_;
};

#endif /* __publicNode_hh__ */

