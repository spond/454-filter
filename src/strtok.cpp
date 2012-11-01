
#include <cstring>

#include "common.hpp"
#include "strtok.hpp"

// strtok_t --------------------------------------------------------------------------------------------------------- //

strtok_t::strtok_t( const char * s ) : str( strdup( s ) )
{
    __CHECK_PTR( str )
    ptr = str;
}

strtok_t::~strtok_t()
{
    free( str );
}

char * strtok_t::next( const char * delim )
{
    char * pch = strpbrk( ptr, delim );

    // skip leading delims
    while ( pch && ptr == pch ) {
        ptr = pch + 1;
        pch = strpbrk( ptr, delim );
    }

    // if we advanced more than 1,
    // and aren't NULL, then set the null byte
    // advance the ptr to just after the null byte,
    // and return
    if ( pch ) {
        char * tmp = ptr;
        pch[0] = '\0';
        ptr = pch + 1;
        return tmp;
    }

    return NULL;
}
