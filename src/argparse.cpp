
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "argparse.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

const char usage[] =
    "usage: " QFILT " [-h] "
    "[-o OUTPUT] "
    "[-q QSCORE] "
    "[-l LENGTH] "
    "[-m MODE] "
    "[-T PREFIX] "
    "[-t MISMATCH] "
    "(-F FASTQ | -Q FASTA QUAL)\n";

const char help_msg[] =
    "filter sequencing data using some simple heuristics\n"
    "\n"
    "required arguments:\n"
    "  -F FASTQ                 FASTQ file\n"
    "  -Q FASTA QUAL            FASTA and QUAL files\n"
    "\n"
    "optional arguments:\n"
    "  -h, --help               show this help message and exit\n"
    "  -o OUTPUT                direct retained fragments to a file named OUTPUT (default=stdout)\n"
    "  -q QSCORE                minimum per-base quality score below which a read will be split\n"
    "                           or truncated (default=" TO_STR( DEFAULT_MIN_QSCORE ) ")\n"
    "  -l LENGTH                minimum retained fragment LENGTH (default=" TO_STR( DEFAULT_MIN_LENGTH ) ")\n"
    "  -m MODE                  MODE is a 3-bitmask (an integer in [0-7], default=" TO_STR( DEFAULT_MODE ) "):\n"
    "                           if the lowest bit is set, a low q-score causes reads to be split,\n"
    "                           otherwise they are truncated;\n"
    "                           if the second bit is set, low q-score homopolymers are tolerated;\n"
    "                           and if the highest bit is set, low q-score 'N's are tolerated\n"
    "  -T PREFIX                if supplied, only reads with this PREFIX are retained,\n"
    "                           and the PREFIX is stripped from each contributing read\n"
    "  -t MISMATCH              if PREFIX is supplied, prefix matching tolerates at most\n"
    "                           MISMATCH mismatches (default=" TO_STR( DEFAULT_TAG_MISMATCH ) ")\n";

inline
void help()
{
    fprintf( stderr, "%s\n%s", usage, help_msg );
    exit( 1 );
}

#define ERROR( msg, args... ) \
{ \
    fprintf( stderr, "%s" QFILT ": error: " msg "\n", usage , ##args ); \
    exit( 1 ); \
}

// args_t ----------------------------------------------------------------------------------------------------------- //

args_t::args_t( int argc, const char * argv[] ) :
    fastq( NULL ),
    fasta( NULL ),
    qual( NULL ),
    output( NULL ),
    min_length( DEFAULT_MIN_LENGTH ),
    min_qscore( DEFAULT_MIN_QSCORE ),
    tag_length( 0 ),
    tag_mismatch( DEFAULT_TAG_MISMATCH )
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
            if ( !strcmp( &arg[2], "help" ) )
                help();

#if 0

            if ( !strcmp( &arg[2], "fastq" ) ) parse_fastq( argv[++i] );
            else if ( !strcmp( &arg[2], "qual" ) ) {
                parse_qual( argv[i + 1], argv[i + 2] );
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
            if ( !strcmp( &arg[1], "F" ) ) parse_fastq( argv[++i] );
            else if ( !strcmp( &arg[1], "Q" ) ) {
                parse_qual( argv[i + 1], argv[i + 2] );
                i += 2;
            }
            else if ( !strcmp( &arg[1], "o" ) ) parse_output( argv[++i] );
            else if ( !strcmp( &arg[1], "l" ) ) parse_minlength( argv[++i] );
            else if ( !strcmp( &arg[1], "q" ) ) parse_minqscore( argv[++i] );
            else if ( !strcmp( &arg[1], "m" ) ) parse_mode( argv[++i] );
            else if ( !strcmp( &arg[1], "T" ) ) parse_tag( argv[++i] );
            else if ( !strcmp( &arg[1], "t" ) ) parse_tagmismatch( argv[++i] );
            else if ( !strcmp( &arg[1], "h" ) ) help();
            else
                ERROR( "unknown argument: %s", arg );
        }
        else
            ERROR( "unknown argument: %s", arg );
    }

    if ( !fastq && ( !fasta || !qual ) )
        ERROR( "missing required argument -F FASTQ or -Q FASTA QUAL" );
}

void args_t::parse_qual( const char * fstr, const char * qstr )
{
    if ( fastq )
        ERROR( "-Q and -F are mutually exclusive" );

    fasta = fstr;
    qual = qstr;
}

void args_t::parse_fastq( const char * str )
{
    if ( fasta || qual )
        ERROR( "-F and -Q are mutually exclusive" );

    fastq = str;
}

void args_t::parse_output( const char * str )
{
    output = str;
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
