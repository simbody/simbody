#include "dint-atom.h"
#include "dint-node.h"
#include "dinternal.h"

ostream& 
operator<<(ostream &s,const IVMAtom* a) 
{
    if (a->node && a->node->ivm)
        s << a->node->ivm->idAtom(a->index);
    else
        s << a->index;
    return s;
}

