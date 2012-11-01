
#ifndef STRTOK_H
#define STRTOK_H

class strtok_t
{
private:
    char * const str,
         * ptr;

public:
    strtok_t( const char * );
    ~strtok_t();
    char * next( const char * );
};

#endif // STRTOK_H
