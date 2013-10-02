
#ifndef IFILE_H
#define IFILE_H

#include <string>
#include <vector>

#define BUF_SZ 256 

namespace ifile
{
    class ifile_t
    {
    public:
        const char * const path;

    private:
        FILE * file;
        size_t line;
        size_t col;
        char buf[BUF_SZ];
        char * end;
        char * ptr;
        std::vector<size_t> cols;

        inline
        void next_col( const size_t ncol=1 ) {
            col += ncol;
        }

        inline
        void next_line() {
            line += 1;
            cols.push_back( col );
            col = 1;
        }

        inline
        void prev_col( const size_t ncol=1 ) {
            col -= ncol;
        }

        inline
        void prev_line() {
            if ( line > 1 and cols.size() > 0 ) {
                line -= 1;
                col = cols.back();
                cols.pop_back();
            }
        }

        bool fill();

    public:
        ifile_t( const char * path=NULL );
        ~ifile_t();
        bool good() const;
        void error( const char *, ... ) const;
        char getc();
        void skip_ws();
        void extend_until( std::string &, const char *, bool trim=true );
    };
}

#endif // IFILE_H
