
#ifndef SEQ_H
#define SEQ_H

#include <cstdio>
#include <string>
#include <vector>

#include "ifile.hpp"

namespace seq
{
    enum state_t {
        UNKNOWN,
        ID,
        SEQUENCE,
        QUALITY
    };

    enum filetype_t {
        FASTA,
        QUAL,
        FASTQ
    };

    class seq_t
    {
    public:
        std::string id;
        std::string seq;
        std::vector<size_t> quals;
        size_t length;
        seq_t();
        void clear();
    };

    class parser_t
    {
    private:
        ifile::ifile_t * fasta;
        ifile::ifile_t * fastq;
        ifile::ifile_t * qual;
        state_t fstate;
        state_t qstate;

        // these are for parsing
        const char * const hdr;
        const char * const sep;

        // these are permanent static buffer
        std::string qid;
        std::string qs;

    public:
        parser_t( ifile::ifile_t * );
        parser_t( ifile::ifile_t *, ifile::ifile_t * );
        bool next( seq_t & );
    };
}

#endif // SEQ_H
