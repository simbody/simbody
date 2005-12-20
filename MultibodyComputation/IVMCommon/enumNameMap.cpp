
#include "enumNameMap.h"
#include <cdsSStream.h>
#include <cdsExcept.h>

const char*
toString(const EnumNameMap* map, 
	 const int          key)
{
 for (int i=0 ; String("last enum") != map[i].name ; i++)
   if ( map[i].val == key )
     return map[i].name;

 OStringStream mess; mess<< "unknown namemap key: " << key;

 cerr << "FATAL: " << mess.str() << endl;
 throw CDS::exception( mess.str() );
} /* toString */
