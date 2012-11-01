
#include "common.hpp"
#include "seq.hpp"
#include "strtok.hpp"

// helper classes and functions

#define MALFUNCTION() \
{ \
    fprintf( stderr, "\nERROR (file: %s, line: %d): state machine malfunction\n", basename( __FILE__ ), __LINE__ ); \
    exit( 1 ); \
}

#define PARSE_ERROR( pos, msg, args... ) \
{ \
    const char * file; \
    long line, col; \
    (pos).get( file, line, col ); \
    fprintf( stderr, "\nERROR (file: %s, line: %ld, column: %ld): " msg "\n", file, line, col , ##args ); \
    exit( 1 ); \
}

// pos_t ------------------------------------------------------------------------------------------------------------ //

pos_t::pos_t( const char * file ) : file( file ), line( 1 ), col( 1 )
{
    cols = new vec_t<unsigned long>;
}

pos_t::~pos_t()
{
    delete cols;
}

void pos_t::get( const char *& f, long & l, long & c ) const
{
    f = file;
    l = line;
    c = col;
}

// seq_t ------------------------------------------------------------------------------------------------------------ //

seq_t::seq_t() : length( 0 )
{
    id = new str_t();
    seq = new str_t();
    quals = new vec_t<long>();
}

seq_t::~seq_t()
{
    delete id;
    delete seq;
    delete quals;
}

void seq_t::clear()
{
    id->clear();
    seq->clear();
    quals->clear();
    length = 0;
}

// parser_t helper functions ---------------------------------------------------------------------------------------- //

inline
void skip_ws( FILE * file, pos_t & pos )
{
    char buf[256];

    while ( fgets( buf, 256, file ) ) {
        char * pch = strpbrk( buf, " \t\n\r" ),
               * ptr = buf;

        while ( pch && pch == ptr ) {
            // if newline, increment line,
            // otherwise increment column
            if ( pch[0] == '\n' )
                pos.next_line();
            else
                pos.next_col();

            // advance to next position
            ptr = pch + 1;
            pch = strpbrk( ptr, " \t\n\r" );
        }

        // if we moved by more than a single position,
        // rewind by what remains in buf and return
        if ( ptr ) {
            fseek( file, ptr - buf - strlen( buf ), SEEK_CUR );
            return;
        }
    }
}

long extend_until( str_t & str, const char * delim, FILE * file, pos_t & pos, bool trim )
{
    char buf[256];
    long total = 0;

    while ( fgets( buf, 256, file ) ) {
        char * const pch = strpbrk( buf, delim );
        const long len = strlen( buf );
        bool truncate = false;

        if ( pch ) {
            const long nchar = pch - buf;
            str.extend( buf, nchar );

            if ( pch[0] == '\n' ) {
                pos.next_col( nchar - 1 );
                pos.next_line();
            }
            else
                pos.next_col( nchar );

            // the +1 positions us just after the delim
            fseek( file, nchar - len + 1, SEEK_CUR );
            return total + nchar;
        }

        // if we reached a newline,
        // set truncate to trim
        // and advance pos
        if ( buf[len - 1] == '\n' ) {
            truncate = trim;
            pos.next_col( len - 1 );
            pos.next_line();
        }
        else
            pos.next_col( len );

        // extend the string by the appropriate number of chars
        str.extend( buf, len - ( truncate ? 1 : 0 ) );
        total += len - ( truncate ? 1 : 0 );

        if ( trim )
            skip_ws( file, pos );
    }

    return total;
}

// by default trim whitespace from the beginning of lines
inline
int extend_until( str_t & str, const char * delim, FILE * file, pos_t & pos )
{
    return extend_until( str, delim, file, pos, true );
}

// parser_t --------------------------------------------------------------------------------------------------------- //

void parser_t::init()
{
    qid = new str_t();
    qs = new str_t();
}

const char chr[] = ">\0@\0+";

parser_t::parser_t( const char * fastq_file ):
    fasta( NULL ),
    qual( NULL ),
    fpos( pos_t( fastq_file ) ),
    qpos( pos_t( NULL ) ),
    fstate( UNKNOWN ),
    qstate( UNKNOWN ),
    hdr( chr + 2 ),
    sep( chr + 4 )
{
    fastq = fopen( fastq_file, "rb" );

    if ( !fastq ) {
        fprintf( stderr, "\nERROR: failed to open the FASTQ file %s\n", fastq_file );
        exit( 1 );
    }

    init();
}

parser_t::parser_t( const char * fasta_file, const char * qual_file ):
    fastq( NULL ),
    fpos( pos_t( fasta_file ) ),
    qpos( pos_t( qual_file ) ),
    fstate( UNKNOWN ),
    qstate( UNKNOWN ),
    hdr( chr + 0 ),
    sep( chr + 0 )
{
    fasta = fopen( fasta_file, "rb" );

    if ( !fasta ) {
        fprintf( stderr, "\nERROR: failed to open the FASTA file %s\n", fasta_file );
        exit( 1 );
    }

    qual = fopen( qual_file, "rb" );

    if ( !qual ) {
        fprintf( stderr, "\nERROR: failed to open the QUAL file %s\n", qual_file );
        exit( 1 );
    }

    init();
}

parser_t::~parser_t()
{
    if ( fasta )
        fclose( fasta );

    if ( fastq )
        fclose( fastq );

    if ( qual )
        fclose( qual );

    delete qid;
    delete qs;
}

bool parser_t::next( seq_t & seq )
{
    FILE * file = fastq ? fastq : fasta;
    filetype_t filetype = fastq ? FASTQ : FASTA;
    pos_t * pos = &fpos;
    state_t * state = &fstate;
begin:

    do {
        switch ( *state ) {
        case UNKNOWN: {
            skip_ws( file, *pos );
            {
                const char c = fgetc( file );
                pos->next_col();

                if ( strchr( hdr, c ) )
                    *state = ID;
                else if ( c == EOF )
                    return false;
                else
                    PARSE_ERROR( *pos, "malformed file: %c", c )
                    break;
            }
        }

        case ID: {
            str_t * str = ( filetype == QUAL ) ? qid : seq.id;
            int nelem = extend_until( *str, "\n", file, *pos );

            if ( nelem < 1 )
                PARSE_ERROR( *pos, "malformed file: missing ID" )
                if ( filetype == QUAL ) {
                    qid->clear();
                    *state = QUALITY;
                }
                else
                    *state = SEQUENCE;

            break;
        }

        case SEQUENCE: {
            switch ( filetype ) {
            case FASTA:
            case FASTQ: {
                int nelem;
                nelem = extend_until( *seq.seq, sep, file, *pos );

                if ( nelem < 1 )
                    PARSE_ERROR( *pos, "malformed file: missing sequence" )
                    if ( filetype == FASTA ) {
                        if ( !feof( file ) ) {
                            const int err = fseek( file, -1, SEEK_CUR );

                            if ( err )
                                MALFUNCTION()
                                pos->prev_col();
                        }

                        *state = UNKNOWN;
                    }
                    else { // FASTQ
                        // skip the whitespace after the + separator
                        skip_ws( file, *pos );
                        *state = QUALITY;
                    }

                break;
            }

            default:
                MALFUNCTION()
            }

            break;
        }

        case QUALITY: {
            // FASTQ files permit '@' to appear as a valid quality value (31),
            // so look for a newline instead of the hdr
            int nelem = extend_until( *qs, ( filetype == QUAL ) ? hdr : "\n", file, *pos, false );

            if ( nelem < 1 )
                PARSE_ERROR( *pos, "malformed file: missing quality scores" )
                if ( filetype == QUAL ) {
                    char * buf = NULL;
                    strtok_t tok( qs->c_str() );

                    while ( ( buf = tok.next( " \t\n" ) ) )
                        seq.quals->append( atoi( buf ) );
                }
                else { // FASTQ
                    int i;

                    for ( i = 0; i < qs->length(); ++i ) {
                        // encoding: chr(phred+33)
                        seq.quals->append( long( ( *qs )[i] ) - 33 );
                    }
                }

            // clear the qual data after use
            qs->clear();

            // reset the state to UNKNOWN and just prior to the header
            if ( filetype == QUAL && !feof( file ) ) {
                fseek( file, -1, SEEK_CUR );
                pos->prev_col();
            }

            *state = UNKNOWN;
            break;
        }

        default:
            MALFUNCTION()
        }
    }
    while ( *state != UNKNOWN );

    if ( qual && file == fasta ) {
        file = qual;
        filetype = QUAL;
        pos = &qpos;
        state = &qstate;
        goto begin;
    }

    if ( seq.seq->length() != seq.quals->length() )
        PARSE_ERROR(
            *pos,
            "malformed file: sequence length (%ld) does not match the number of quality scores (%ld)",
            seq.seq->length(),
            seq.quals->length()
        )
        else
            seq.length = seq.seq->length();

    return true;
}
