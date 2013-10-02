
#include "common.hpp"
#include "seq.hpp"
#include "strtok.hpp"

// helper classes and functions

#define MALFUNCTION() \
{ \
    fprintf( stderr, "\nERROR (file: %s, line: %d): state machine malfunction\n", __basename( __FILE__ ), __LINE__ ); \
    exit( 1 ); \
}

#define IS_WHITESPACE( chr ) ( ( chr ) == ' ' || ( chr ) == '\t' || ( chr ) == '\r' || ( chr ) == '\n' )

namespace seq
{
    seq_t::seq_t() : length( 0 ) { }

    void seq_t::clear()
    {
        id.clear();
        seq.clear();
        quals.clear();
        length = 0;
    }

    const char chr[] = ">\0@\0+";

    parser_t::parser_t( ifile::ifile_t * fastq ) :
        fasta( NULL ), 
        fastq( fastq ),
        qual( NULL ), 
        fstate( UNKNOWN ),
        qstate( UNKNOWN ),
        hdr( chr + 2 ),
        sep( chr + 4 )
    {
    }

    parser_t::parser_t( ifile::ifile_t * fasta, ifile::ifile_t * qual ) :
        fasta( fasta ),
        fastq( NULL ),
        qual( qual ),
        fstate( UNKNOWN ),
        qstate( UNKNOWN ),
        hdr( chr + 0 ),
        sep( chr + 0 )
    {
    }

    bool parser_t::next( seq_t & seq )
    {
        ifile::ifile_t * file = fastq ? fastq : fasta;
        filetype_t filetype = fastq ? FASTQ : FASTA;
        state_t * state = &fstate;
    begin:

        do {
            switch ( *state ) {
            case UNKNOWN: {
                file->skip_ws();
                {
                    const char chr = file->getc();

                    if ( strchr( hdr, chr ) )
                        *state = ID;
                    else if ( chr == EOF )
                        return false;
                    else
                        file->error( "malformed file: %c", chr );
                }
                break;
            }

            case ID: {
                std::string & str = ( filetype == QUAL ) ? qid : seq.id;
                file->extend_until( str, "\r\n" );

                if ( str.length() < 1 )
                    file->error( "malformed file: missing ID" );

                if ( filetype == QUAL ) {
                    qid.clear();
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
                    file->extend_until( seq.seq, sep );

                    if ( !seq.seq.length() )
                        file->error( "malformed file: missing sequence" );

                    if ( filetype == FASTA )
                        *state = UNKNOWN;
                    else { // FASTQ
                        // skip over the + separator and any trailing whitespace
                        file->getc();
                        file->extend_until( qid, "\r\n" );
                        qid.clear();
                        *state = QUALITY;
                    }

                    break;
                }

                default:
                    MALFUNCTION();
                }

                break;
            }

            case QUALITY: {
                // FASTQ files permit '@' to appear as a valid quality value (31),
                // so look for a newline instead of the hdr
                file->extend_until( qs, ( filetype == QUAL ) ? hdr : "\r\n", false );

                if ( qs.length() < 1 )
                    file->error( "malformed file: missing quality scores" );

                if ( filetype == QUAL ) {
                    char * buf = NULL;

                    strtok_t tok( qs.c_str() );

                    while ( ( buf = tok.next( " \t\r\n" ) ) )
                        seq.quals.push_back( atoi( buf ) );
                }
                else { // FASTQ
                    size_t i;

                    for ( i = 0; i < qs.length(); ++i ) {
                        // encoding: chr(phred+33)
                        seq.quals.push_back( long( qs[i] ) - 33 );
                    }
                }

                // clear the qual data after use
                qs.clear();

                *state = UNKNOWN;
                break;
            }

            default:
                MALFUNCTION();
            }
        }
        while ( *state != UNKNOWN );

        if ( qual && file == fasta ) {
            file = qual;
            filetype = QUAL;
            state = &qstate;
            goto begin;
        }

        if ( ( qual || fastq ) && seq.seq.length() != seq.quals.size() )
            file->error(
                "malformed file: sequence length (%ld) does not match the number of quality scores (%ld)",
                seq.seq.length(),
                seq.quals.size()
            );
        else
            seq.length = seq.seq.length();

        return true;
    }
}
