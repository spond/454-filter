
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "ifile.hpp"

namespace ifile
{
    ifile_t::ifile_t( const char * path ) :
        path( path ),
        file( NULL ),
        line( 0 ),
        col( 0 ),
        end( buf ),
        ptr( buf )
    {
        if ( path ) {
            if ( !strcmp( path, "-" ) )
                file = stdin;
            else
                file = fopen( path, "rb" );
        }
    }

    ifile_t::~ifile_t()
    {
        if ( file && file != stdin ) {
            fclose( file );
            file = NULL;
        }
    }

    bool ifile_t::good() const
    {
        return file != NULL;
    }

    void ifile_t::error( const char * msg, ... ) const
    {
        va_list args;
        fprintf( stderr, "\nERROR (file: %s, line: %ld, column: %ld): ", path, line, col );
        va_start( args, msg );
        vfprintf( stderr, msg, args );
        va_end( args );
        fprintf( stderr, "\n" );
        exit( 1 );
    }

    bool ifile_t::fill()
    {
        ptr = fgets( buf, BUF_SZ, file );

        if ( ptr )
            end = buf + strlen( buf );
        else
            end = buf;

        return ptr != NULL;
    }

    char ifile_t::getc()
    {
        char chr = EOF;

        if ( ptr && ptr >= end )
            fill();

        if ( ptr ) {
            chr = ptr[0];
            if ( chr == '\n' )
                next_line();
            else if ( chr != '\r' )
                next_col();
            ++ptr;
        }

        return chr;
    }

    void ifile_t::skip_ws()
    {
        while( true ) {
            switch ( getc() ) {
            case ' ':
            case '\t':
            case '\r':
            case '\n':
                break;
            case EOF:
                return;
            default:
                --ptr;
                prev_col();
                return;
            }
        }
    }

    void ifile_t::extend_until( std::string & str, const char * delim, bool trim )
    {
        for( ; ptr != NULL ; fill() ) {
            const char * pch = strpbrk( ptr, delim );
            const size_t len = ( end > buf ) ? end - buf : 0;
            size_t truncate = 0;

            if ( pch ) {
                const size_t nchar = ( pch > ptr ) ? pch - ptr : 0;

                if ( nchar ) {
                    str.append( ptr, nchar );

                    switch ( pch[0] ) {
                    case '\n':
                        next_col( nchar - 1 );
                        next_line();
                        break;
                    case '\r':
                        next_col();
                        break;
                    }

                    ptr += nchar;
                }

                if ( trim )
                    skip_ws();

                return;
            }

            if ( len > 0 ) {
                if ( end[-1] == '\n' ) {
                    size_t nchar = 1;

                    if ( len > 1 && end[-2] == '\r' )
                        nchar = 2;

                    if ( trim )
                        truncate = nchar;

                    next_col( len - nchar );
                    next_line();
                }
                else if ( end[-1] == '\r' ) {
                    next_col( len - 1 );
                }
                else
                    next_col( len );
            }

            str.append( ptr, len - truncate );
        }
    }
}
