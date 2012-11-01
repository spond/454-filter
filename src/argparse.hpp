
#ifndef ARGPARSE_H
#define ARGPARSE_H

// program name
#define QFILT "qfilt"

// argument defaults
#define DEFAULT_MIN_LENGTH 50
#define DEFAULT_MIN_QSCORE 20
#define DEFAULT_MODE 0
#define DEFAULT_TAG_MISMATCH 0

class args_t
{
public:
    const char * fastq;
    const char * fasta;
    const char * qual;
    const char * output;
    long min_length;
    long min_qscore;
    bool split; // split not truncate
    bool hpoly; // tolerate homopolymers
    bool ambig; // tolerate ambigs ('N')
    char tag[256];
    long tag_length;
    long tag_mismatch;

    args_t( int, const char ** );
private:
    void parse_fastq( const char * );
    void parse_qual( const char *, const char * );
    void parse_output( const char * );
    void parse_minlength( const char * );
    void parse_minqscore( const char * );
    void parse_mode( const char * );
    void parse_tag( const char * );
    void parse_tagmismatch( const char * );
};

#endif // ARGPARSE_H
