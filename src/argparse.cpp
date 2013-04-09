
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "argparse.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
    const char usage[] =
        "usage: " QFILT " [-h] "
        "[-o OUTPUT] "
        "[-q QSCORE] "
        "[-l LENGTH] "
        "[-m MODE] [-s] [-p] [-a] "
        "[-T PREFIX] "
        "[-t MISMATCH] "
        "[-f] "
        "( -F FASTA QUAL | -Q FASTQ )\n";

    const char help_msg[] =
        "filter sequencing data using some simple heuristics\n"
        "\n"
        "required arguments:\n"
        "  -F FASTA QUAL            FASTA and QUAL files\n"
        "  -Q FASTQ                 FASTQ file\n"
        "\n"
        "optional arguments:\n"
        "  -h, --help               show this help message and exit\n"
        "  -o OUTPUT                direct retained fragments to a file named OUTPUT (default=stdout)\n"
        "  -q QSCORE                minimum per-base quality score below which a read will be split\n"
        "                           or truncated (default=" TO_STR( DEFAULT_MIN_QSCORE ) ")\n"
        "  -l LENGTH                minimum retained fragment LENGTH (default=" TO_STR( DEFAULT_MIN_LENGTH ) ")\n"
        "  -m MODE                  MODE is a 3-bitmask (an integer in [0-7], default=" TO_STR( DEFAULT_MODE ) "):\n"
        "                           if the lowest bit is set, it is like passing -s;\n"
        "                           if the middle bit is set, it is like passing -p;\n"
        "                           and if the highest bit is set, it is like passing -a\n"
        "  -s                       when encountering a low q-score, split instead of truncate\n"
        "  -p                       tolerate low q-score homopolymeric regions\n"
        "  -a                       tolerate low q-score ambiguous nucleotides\n"
        "  -T PREFIX                if supplied, only reads with this PREFIX are retained,\n"
        "                           and the PREFIX is stripped from each contributing read\n"
        "  -t MISMATCH              if PREFIX is supplied, prefix matching tolerates at most\n"
        "                           MISMATCH mismatches (default=" TO_STR( DEFAULT_TAG_MISMATCH ) ")\n"
        "  -f FORMAT                output in FASTA or FASTQ format (default=" TO_STR( DEFAULT_FORMAT ) ")\n";

    inline
    void help()
    {
        fprintf( stderr, "%s\n%s", usage, help_msg );
        exit( 1 );
    }

    inline
    void ERROR( const char * msg, ... )
    {
        va_list args;
        fprintf( stderr, "%s" QFILT ": error: ", usage );
        va_start( args, msg );
        vfprintf( stderr, msg, args );
        va_end( args );
        fprintf( stderr, "\n" );
        exit( 1 );
    }

    args_t::args_t( int argc, const char * argv[] ) :
        fasta( NULL ),
        fastq( NULL ),
        qual( NULL ),
        output( stdout ),
        min_length( DEFAULT_MIN_LENGTH ),
        min_qscore( DEFAULT_MIN_QSCORE ),
        tag_length( 0 ),
        tag_mismatch( DEFAULT_TAG_MISMATCH ),
        format( DEFAULT_FORMAT )
    {
        int i;
        // make sure tag is an empty string
        tag[0] = '\0';
        // handle the mode separately
        parse_mode( TO_STR( DEFAULT_MODE ) );

        // skip arg[0], it's just the program name
        for ( i = 1; i < argc; ++i ) {
            const char * arg = argv[i];

            if ( arg[0] == '-' && arg[1] == '-' ) {
                if ( !strcmp( &arg[2], "help" ) ) help();
#if 0
                else if ( !strcmp( &arg[2], "fastq" ) ) parse_fastq( argv[++i] );
                else if ( !strcmp( &arg[2], "fasta" ) ) {
                    parse_fasta( argv[i + 1], argv[i + 2] );
                    i += 2;
                }
                else if ( !strcmp( &arg[2], "minlength" ) ) parse_minlength( argv[++i] );
                else if ( !strcmp( &arg[2], "minqscore" ) ) parse_minqscore( argv[++i] );
                else if ( !strcmp( &arg[2], "mode" ) ) parse_mode( argv[++i] );
                else if ( !strcmp( &arg[2], "tag" ) ) parse_tag( argv[++i] );
                else if ( !strcmp( &arg[2], "tagmismatch" ) ) parse_tagmismatch( argv[++i] );
#endif
                else
                    ERROR( "unknown argument: %s", arg );
            }
            else if ( arg[0] == '-' ) {
                if ( !strcmp( &arg[1], "h" ) ) help();
                else if ( !strcmp( &arg[1], "F" ) ) {
                    parse_fasta( argv[i + 1], argv[i + 2] );
                    i += 2;
                }
                else if ( !strcmp( &arg[1], "Q" ) ) parse_fastq( argv[++i] );
                else if ( !strcmp( &arg[1], "o" ) ) parse_output( argv[++i] );
                else if ( !strcmp( &arg[1], "l" ) ) parse_minlength( argv[++i] );
                else if ( !strcmp( &arg[1], "q" ) ) parse_minqscore( argv[++i] );
                else if ( !strcmp( &arg[1], "m" ) ) parse_mode( argv[++i] );
                else if ( !strcmp( &arg[1], "s" ) ) parse_split();
                else if ( !strcmp( &arg[1], "p" ) ) parse_hpoly();
                else if ( !strcmp( &arg[1], "a" ) ) parse_ambig();
                else if ( !strcmp( &arg[1], "T" ) ) parse_tag( argv[++i] );
                else if ( !strcmp( &arg[1], "t" ) ) parse_tagmismatch( argv[++i] );
                else if ( !strcmp( &arg[1], "f" ) ) parse_format( argv[++i] );
                else
                    ERROR( "unknown argument: %s", arg );
            }
            else
                ERROR( "unknown argument: %s", arg );
        }

        if ( !fastq && ( !fasta || !qual ) )
            ERROR( "missing required argument -F FASTA QUAL or -Q FASTQ" );
    }

    args_t::~args_t() {
        if ( fasta )
            delete fasta;
        if ( fastq )
            delete fastq;
        if ( qual )
            delete qual;
        if ( output && output != stdin )
            fclose( output );
    }

    void args_t::parse_fasta( const char * fstr, const char * qstr )
    {
        if ( fastq )
            ERROR( "-F and -Q are mutually exclusive" );

        if ( !strcmp( fstr, "-" ) && !strcmp( qstr, "-" ) )
            ERROR( "only one argument to -F FASTA and QUAL can be STDIN" );

        fasta = new ifile::ifile_t( fstr );
        qual = new ifile::ifile_t( qstr );

        if ( !fasta || !fasta->good() )
            ERROR( "failed to open the FASTA file %s", fstr );

        if ( !qual && !qual->good() )
            ERROR( "failed to open the QUAL file %s", qstr );
    }

    void args_t::parse_fastq( const char * str )
    {
        if ( fasta || qual )
            ERROR( "-Q and -F are mutually exclusive" );

        fastq = new ifile::ifile_t( str );

        if ( !fastq->good() )
            ERROR( "failed to open the FASTQ file %s", str );
    }

    void args_t::parse_output( const char * str )
    {
        if ( str && strcmp( str, "-" ) )
            output = fopen( str, "wb" );
        else
            output = stdin;

        if ( !output )
            ERROR( "failed to open the OUTPUT file %s", str );
    }

    void args_t::parse_minlength( const char * str )
    {
        min_length = atoi( str );

        if ( min_length < 1 )
            ERROR( "minimum length expected a positive integer, had: %s", str );
    }

    void args_t::parse_minqscore( const char * str )
    {
        min_qscore = atoi( str );

        if ( min_qscore < 0 )
            ERROR( "min q-score expected a non-negative integer, had: %s", str );
    }

    void args_t::parse_mode( const char * str )
    {
        int mode = atoi( str );

        if ( mode < 0 || mode > 7 )
            ERROR( "mode must be an integer in [0, 7], had: %s", str );

        split = ( mode & 1 );
        hpoly = ( mode & 2 );
        ambig = ( mode & 4 );
    }

    void args_t::parse_split()
    {
        split = true;
    }

    void args_t::parse_hpoly()
    {
        hpoly = true;
    }

    void args_t::parse_ambig()
    {
        ambig = true;
    }

    void args_t::parse_tag( const char * str )
    {
        int nvar = sscanf( str, "%256s", tag );

        if ( nvar != 1 )
            ERROR( "failed to process tag argument %s", str );

        tag_length = strlen( tag );
    }

    void args_t::parse_tagmismatch( const char * str )
    {
        tag_mismatch = atoi( str );

        if ( tag_mismatch < 0 )
            ERROR( "maximum tag mismatch expected non-negative integer, had: %s", str );
    }

    void args_t::parse_format( const char * str )
    {
        if ( !strcmp( str, "FASTA" ) )
            format = FASTA;
        else if ( !strcmp( str, "FASTQ" ) )
            format = FASTQ;
        else
            ERROR( "invalid format %s", str );
    }
}
