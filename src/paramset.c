#include <limits.h>

#include "paramset.h"

struct trfparamset {
    unsigned int match;
    unsigned int mismatch;
    unsigned int indel;
    unsigned int minscore;
    unsigned int maxperiod;
    unsigned int PM;
    unsigned int PI;
    int datafile;
    int maskedfile;
    int flankingsequence;
    unsigned int flankinglength;
    int HTMLoff;
    int redundoff;
    int ngs;
    int use_stdin;
    unsigned int maxwraplength;

    char inputfilename[PATH_MAX];  /* constant defined in stdlib */
    char outputprefix[PATH_MAX];
    char outputdirectory[PATH_MAX];
    char outputfilename[PATH_MAX];
    int multisequencefile;      /* flags if file has more than one sequence */
    int sequenceordinal;        /* holds seq. index starting on 1 */
    int outputcount;            /* repeats found */
    int running;
    char *endstatus;
    int percent;
} ;

struct trfparamset g_paramset;

int *paramset_datafile(void)
{
    return &g_paramset.datafile;
}

int paramset_maskedfile(void)
{
    return g_paramset.maskedfile;
}

int paramset_flankingsequence(void)
{
    return g_paramset.flankingsequence;
}

int paramset_HTMLoff(void)
{
    return g_paramset.HTMLoff;
}

int paramset_redundoff(void)
{
    return g_paramset.redundoff;
}

int paramset_ngs(void)
{
    return g_paramset.ngs;
}

