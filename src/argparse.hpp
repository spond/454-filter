
#ifndef ARGPARSE_H
#define ARGPARSE_H

#include "ifile.hpp"

// program name
#define PROGNAME "qfilt"

// argument defaults
#define DEFAULT_MIN_LENGTH 50
#define DEFAULT_MIN_QSCORE 20
#define DEFAULT_MODE 0
#define DEFAULT_TAG_MISMATCH 0
#define DEFAULT_FORMAT FASTA

namespace argparse
{
    enum format_t {
        FASTA,
        FASTQ
    };

    class args_t
    {
    public:
        ifile::ifile_t * fasta;
        ifile::ifile_t * fastq;
        ifile::ifile_t * qual;
        FILE * output;
        size_t min_length;
        size_t min_qscore;
        bool split; // split not truncate
        bool hpoly; // tolerate homopolymers
        bool ambig; // tolerate ambigs ('N')
        bool json; // diagnostics to JSON
        char punch;
        char tag[256];
        size_t tag_length;
        size_t tag_mismatch;
        format_t format;

        args_t( int, const char ** );
        ~args_t();
    private:
        void parse_fastq( const char * );
        void parse_fasta( const char *, const char * );
        void parse_output( const char * );
        void parse_minlength( const char * );
        void parse_minqscore( const char * );
        void parse_mode( const char * );
        void parse_split();
        void parse_hpoly();
        void parse_ambig();
        void parse_punch( const char * );
        void parse_json();
        void parse_tag( const char * );
        void parse_tagmismatch( const char * );
        void parse_format( const char * );
    };
}

#endif // ARGPARSE_H
