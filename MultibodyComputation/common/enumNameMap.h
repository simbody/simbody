
#ifndef __enumNameMap_hh__
#define __enumNameMap_hh__

typedef struct {const char* name; int val;} EnumNameMap;

// given a namemap and a key, return the associated name string
//
const char* toString(const EnumNameMap* map, 
		     const int          key);

// given a namemap and a name string, return the associated key
//
int fromString(const EnumNameMap* map, 
	       const char*        name);

#endif /* __enumNameMap_hh__ */
