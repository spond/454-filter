
#include <cctype>

#include "str.hpp"

// str_t ------------------------------------------------------------------------------------------------------------ //

str_t::str_t() {}

str_t::~str_t() {}

str_t::str_t( const char * str )
{
    extend( str );
}

const char * str_t::c_str() const
{
    return data;
}

void str_t::extend( const char * str )
{
    extend( str, long( strlen( str ) ) );
}

void str_t::lower()
{
    long i;

    for ( i = 0; i < len; ++i )
        data[i] = tolower( data[i] );
}

void str_t::upper()
{
    long i;

    for ( i = 0; i < len; ++i )
        data[i] = toupper( data[i] );
}
