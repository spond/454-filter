
#include <cstdio>

#include "str.h"
#include "vec.h"

#ifndef SEQ_H
#define SEQ_H

enum state_t
{
    UNKNOWN,
    ID,
    SEQUENCE,
    QUALITY
};

enum filetype_t
{
    FASTA,
    QUAL,
    FASTQ
};

class pos_t {
  private:
    const char * file;
    unsigned long line;
    unsigned long col;
    vec_t<unsigned long> * cols;

  public:
    pos_t(const char * file);
    ~pos_t();
    void get(const char * & f, long & l, long & c);

    inline
    void next_col()
    {
        col += 1;
    }
    
    inline
    void next_line()
    {
        line += 1;
        cols->append(col);
        col = 1;
    }

    inline
    void prev_col()
    {
        col -= 1;
    }

    inline
    void prev_line()
    {
        line -= 1;
        col = (*cols)[line - 1];
    }
};

class seq_t {
  public:
    str_t * id;
    str_t * seq;
    vec_t<long> * quals;
    long length;
    seq_t();
    void clear();
    ~seq_t();
};

class parser_t
{
  private:
    FILE * fasta;
    FILE * fastq;
    FILE * qual;
    pos_t fpos;
    pos_t qpos;
    state_t fstate;
    state_t qstate;

    // these are permanent static buffer
    str_t * qid;
    str_t * qs;

  protected:
    void init();

  public:
    parser_t(const char *);
    parser_t(const char *, const char *);
    ~parser_t();
    bool next(seq_t &);
};

#endif // SEQ_H