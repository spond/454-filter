
#ifndef COMMON_H
#define COMMON_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

inline
const char * __basename( const char * path )
{
    long i;

    if ( !path )
        return NULL;

    for ( i = strlen( path ) - 1; i >= 0; --i ) {
#ifdef _WIN32

        if ( path[i] == '\\' )
#else
        if ( path[i] == '/' )
#endif
            return path + i + 1;
    }

    return path;
}

#define __CHECK_PTR( ptr ) \
{ \
    if ( !ptr ) { \
        fprintf( stderr, "\nERROR (file: %s, line: %d): memory allocation failure\n", __basename( __FILE__ ), __LINE__ ); \
        exit( 1 ); \
    } \
}

#endif // COMMON_H
