
#ifndef STR_H
#define STR_H

#include "vec.hpp"

class str_t : public vec_t<char>
{
public:
    using vec_t<char>::extend;

    str_t();
    ~str_t();
    str_t( const char * str );
    const char * c_str() const;
    void extend( const char * str );
    void lower();
    void upper();
};

#endif // STR_H
