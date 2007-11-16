
#include "SimTKsimbody.h"
#pragma warning(disable:4661)
class Element;
class ElementRep;

//extern template class /*__declspec(dllimport)*/ SimTK::PIMPLHandle<Element, ElementRep>;


class /*__declspec(dllexport)*/ Element : public SimTK::PIMPLHandle<Element,ElementRep> {
public:
};



#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#include <iostream>

int main() 
{
    std::cout << "passed" << std::endl;

    Element e;

    return 0;
}


class ElementRep : public SimTK::PIMPLImplementation<Element,ElementRep> {
public:
    ElementRep* clone() const {return new ElementRep();}
};

template class SimTK::PIMPLHandle<Element,ElementRep>;