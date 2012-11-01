
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "argparse.hpp"
#include "seq.hpp"
#include "vec.hpp"

#if 0
static const char * const valid_chars = "ACGTNacgtn";
static long char_lookup[256];
#endif

// vec must be sorted
void fprint_vector_stats( FILE * file, vec_t<long> & vec, const char * hdr )
{
    double sum = 0.,
           var = 0.,
           mean = 0.,
           median = 0.;
    long i,
         min = 0,
         two5 = 0,
         ninetyseven5 = 0,
         max = 0;

    for ( i = 0; i < vec.length(); ++i ) {
        sum += vec[i];
        var += vec[i] * vec[i];
    }

    // remember, i == vec.length()

    if ( vec.length() ) {
        var = ( var - ( sum * sum ) / vec.length() ) / ( vec.length() - 1 );
        mean = sum / i;
        median = ( i % 2 ) ? 1.0 * vec[i / 2] : 0.5 * ( vec[i / 2] + vec[i / 2 - 1] );
        min = vec[0];
        two5 = vec[long( 0.025 * i )];
        ninetyseven5 = vec[long( 0.975 * i )];
        max = vec[i - 1];
    }

    fprintf( file, "%s\n"
             "    mean:                %g\n"
             "    median:              %g\n"
             "    variance             %g\n"
             "    standard deviation:  %g\n"
             "    min:                 %ld\n"
             "    2.5%%:                %ld\n"
             "    97.5%%:               %ld\n"
             "    max:                 %ld\n",
             hdr,
             mean,
             median,
             var,
             sqrt( var ),
             min,
             two5,
             ninetyseven5,
             max
           );
}

// main ------------------------------------------------------------------------------------------------------------- //

int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    parser_t * parser = NULL;
    seq_t seq = seq_t();
    long ncontrib = 0;
    // if -o is unsupplied or -o - is supplied, use stdout
    FILE * output = ( args.output && strcmp( args.output, "-" ) )
                    ? fopen( args.output, "wb" )
                    : stdout;

    if ( !output ) {
        fprintf( stderr, "\nERROR: failed to open OUTPUT file %s\n", args.output );
        exit( 1 );
    }

    // initialize the parser
    if ( args.fastq )
        parser = new parser_t( args.fastq );
    else
        parser = new parser_t( args.fasta, args.qual );

    if ( !parser ) {
        fprintf( stderr, "\nERROR: failed to initialize parser\n" );
        exit( 1 );
    }

#if 0

    for ( int i = 0; i < 256; ++i )
        char_lookup[i] = -1;

    for ( int i = 0; i < valid_char_count; ++i )
        char_lookup[valid_chars[i]] = i;

#endif
    vec_t<long> read_lengths = vec_t<long>();
    vec_t<long> fragment_lengths = vec_t<long>();

    for ( ; parser->next( seq ); seq.clear() ) {
        // maxto is the maximum value of "to",
        // NOT THE UPPER BOUND
        const long maxto = seq.length - args.min_length;
        long nfragment = 0,
             to = 0;
        read_lengths.append( seq.length );

        // compare the sequence prefix to the tag,
        // if it matches by at least tag_mismatch,
        // keep the sequence, otherwise discard
        if ( args.tag_length ) {
            long mismatch = 0;

            for ( to = 0; to < args.tag_length; ++to ) {
                // tolower -> case insensitive
                if ( toupper( ( *seq.seq )[to] ) != toupper( args.tag[to] ) )
                    mismatch += 1;
            }

            if ( mismatch > args.tag_mismatch )
                continue;
        }

        // if we're splitting,
        // continue the following process until we reach the end of the sequence,
        // but only continue if there's enough left to produce a minimum-sized fragment
        while ( true ) {
            long from = 0,
                 i = 0,
                 nambigs = 0;

            // push through the sequence until the quality score meets the minimum
            while ( ( to <= maxto ) && ( ( *seq.quals )[to] < args.min_qscore ) ) {
                to += 1;
            }

            // if we don't have enough length left,
            // skip to the next sequence
            if ( to > maxto )
                break;

            // begin with positive quality score
            from = to;

            // build a read until we hit a low quality score,
            // that is, unless we're skipping Ns or retaining homopolymers
            for ( ; to < seq.length; ++to ) {
                char curr = ( *seq.seq )[to],
                     last = -1;

                if ( ( *seq.quals )[to] < args.min_qscore ) {
                    // if homopolymer (toupper -> case insensitive), continue (last == curr)
                    if ( args.hpoly && toupper( last ) == toupper( curr ) )
                        continue;
                    // if skipping Ns, continue (without assigning last)
                    else if ( args.ambig && ( curr == 'N' || curr == 'n' ) ) {
                        nambigs += 1;
                        continue;
                    }
                    // otherwise, ABORT!!!
                    else
                        break;
                }

                last = curr;
            }

            // "to" is now the upper bound

            // if our fragment isn't long enough,
            // skip to the next fragment
            if ( to - from - nambigs < args.min_length )
                continue;

            // print the read ID
            fprintf( output, ">%s", seq.id->c_str() );

            // print the fragment identifier
            if ( nfragment > 0 )
                fprintf( output, " fragment=%ld\n", nfragment + 1 );
            else {
                fprintf( output, "\n" );
                // if it's the first fragment,
                // count the contributing read
                ncontrib += 1;
            }

            // print the read sequence
            for ( i = from; i < to; i += 60 ) {
                char buf[61];
                const int nitem = ( to - i < 60 ) ? to - i : 60;
                strncpy( buf, seq.seq->c_str() + i, nitem );
                buf[nitem] = '\0';
                fprintf( output, "%s\n", buf );
            }

#if 0
            // for printing quality scores
            fprintf( output, "+\n" );

            for ( i = from; i < to; ++i ) {
                char s[] = " ";

                if ( i == from )
                    s[0] = '\0';

                fprintf( output, "%s%ld", s, ( *seq.quals )[i] );
            }

            fprintf( output, "\n" );
#endif
            fragment_lengths.append( to - from - nambigs );

            if ( !args.split )
                break;

            // only increment fragment identifier after printing
            nfragment += 1;
        }
    }

    fprintf( stderr, "run settings:\n" );

    if ( args.fasta )
        fprintf( stderr,
                 "    input fasta:         %s\n"
                 "    input qual:          %s\n",
                 args.fasta,
                 args.qual
               );
    else
        fprintf( stderr,
                 "    input fastq:         %s\n",
                 args.fastq
               );

    fprintf( stderr,
             "    min q-score:         %ld\n"
             "    min fragment length: %ld\n"
             "    run mode:            %d (%s/%s/%s)\n",
             args.min_qscore,
             args.min_length,
             ( ( args.split ? 1 : 0 ) | ( args.hpoly ? 2 : 0 ) | ( args.ambig ? 4 : 0 ) ),
             args.split ? "split" : "truncate",
             args.hpoly ? "retain homopolymers" : "don't retain homopolymers",
             args.ambig ? "skip ambigs" : "don't skip ambigs"
           );

    if ( args.tag_length )
        fprintf( stderr,
                 "    5' tag:              %s\n"
                 "    max tag mismatches:  %ld\n",
                 args.tag,
                 args.tag_mismatch
               );

    fprintf( stderr,
             "\n"
             "run summary:\n"
             "    original reads:      %ld\n"
             "    contributing reads:  %ld\n"
             "    retained fragments:  %ld\n",
             read_lengths.length(),
             ncontrib,
             fragment_lengths.length()
           );
    // print original read length and retained fragment length statistics
    read_lengths.sort();
    fragment_lengths.sort();
    fprint_vector_stats( stderr, read_lengths, "\noriginal read length distribution:" );
    fprint_vector_stats( stderr, fragment_lengths, "\nretained fragment length distribution:" );

    delete parser;

    if ( output && output != stdout )
        fclose( output );

    return 0;
}
