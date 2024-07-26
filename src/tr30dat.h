#ifndef TR30DAT_H
#define TR30DAT_H

#include <stdio.h>

/*
Tandem Repeats Finder
Copyright (C) 1999-2020 Gary Benson

This file is part of the Tandem Repeats Finder (TRF) program.

TRF is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

TRF is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with TRF.  If not, see <https://www.gnu.org/licenses/>.
*/

/* to disc to improve performance             */
int counterInSeq = 0;

/* uncomment only one platform target identifier */

#define UNIXCONSOLE
#define _MAX_PATH 260

/* use semantic versioning, please: https://semver.org/ */
#ifndef PACKAGE_VERSION
#define versionstring "4.10.0"
#else
#define versionstring PACKAGE_VERSION
#endif

// TODO: Replace with true and false. boa@20240726
#define TRUE 1
#define FALSE 0

#define WEIGHTCONSENSUS 0
#define REPEATCONSENSUS 0
#define MAXDISTANCECONSTANT 2000    /* should replaced by a variable, later */
#define MINDISTANCE 10
#define MINBANDRADIUS 6
#define RECENTERCRITERION 3

#define MAXBANDWIDTH 150
#define MAXTUPLESIZES 10

/* Since this is no longer a macro, use all lower case to avoid confusion. */
unsigned int maxwraplength = 0;

#define MAXPATTERNSIZECONSTANT MAXDISTANCECONSTANT // replaced by a variable 
#define DASH '-'
#define BLANK ' '

#define WITHCONSENSUS 1
#define WITHOUTCONSENSUS 0

#define MULTIPLES 3             /* We keep MULTIPLES*d distance entries */
#define FILLMULTIPLE 3


// how close the lower end of distance indices has to be to the end of one copy
#define LOWERENDMULTIPLE 3

#define ACCEPTED 1
#define NOTACCEPTED 0

// these are sizes for which we do a full array wraparound dynamic
// programming and for which we set d_range by hand 
#define SMALLDISTANCE 20        

/* minimum number of places to store a tuple match.  Usually this is the same
 * as the distance, but for small distances, we allow more tuple matches 
 * because we want see a significant number of matches */
int Min_Distance_Entries = 20;  

/* minimum size of distance window. */
/* Usually this is the same as the */
/* distance, but for small distances we */
/* allow more space because we want a */
/* significant number of matches */
int Min_Distance_Window = 20;


// index separation for tags table for linking 
// active distances for distance range addition of 
// matches 
#define TAGSEP 50

int PM;
int PI;

/* expected probability of a single character indel in the worst case tandem
 * repeat. Pindel should be tied to the indel cost parameter */
double Pindel;

int MAXDISTANCE = 500;
int MAXPATTERNSIZE = 500;
char debugbuffer[500];

/* G. Benson 1/28/2004 */

/* size of EC increased to avoid memory error when consensus length exceeds
   MAXPATTERNSIZECONSTANT after returning from get_consensus(d) */

/* Y. Hernandez 10/15/2018 */

/* If patternsize over 2000 is ever allowed, must change how this
 * variable is initialized. Must be a dynamically allocated array,
 * should use MAXPATTERNSIZE instead (but set elsewhere, after
 * user parameters have been processed).
 */
unsigned char EC[2 * (MAXPATTERNSIZECONSTANT + 1)];

int *Index;
int *ACGTcount;

unsigned char *Sequence;
int Length;

/* int S[MAXWRAPLENGTH+1][MAXPATTERNSIZE];*/
int Delta;                      /* indel penalty */
int Alpha;                      /* match bonus */
int Beta;                       /* mismatch penalty */
int AFDelta;                    /* affine gap initiation penalty */
int AFGamma;                    /* affine gap extension penalty */
int pwidth = 75;
int Reportmin, Heading;
int Classlength;
int Test;
double Rows;
double Totalcharacters;

int Wrapend;
int Maxrealrow, Maxrow, Maxcol;
int Maxscore;
int ConsClasslength;
int *Tag;                       /* list of tags for linking active distances */
int Toptag;                     /* last tag in list */

struct pairalign {
    int length;
    int score;
    char *textprime, *textsecnd;
    int *indexprime, *indexsecnd;
} AlignPair;

struct cons_data {
    char pattern[2 * (MAXPATTERNSIZECONSTANT + 1)];
    int A[2 * (MAXPATTERNSIZECONSTANT + 1)],
        C[2 * (MAXPATTERNSIZECONSTANT + 1)],
        G[2 * (MAXPATTERNSIZECONSTANT + 1)],
        T[2 * (MAXPATTERNSIZECONSTANT + 1)],
        dash[2 * (MAXPATTERNSIZECONSTANT + 1)],
        insert[2 * (MAXPATTERNSIZECONSTANT + 1)],
        letters[2 * (MAXPATTERNSIZECONSTANT + 1)], total[2 * (MAXPATTERNSIZECONSTANT + 1)];
} Consensus;

struct bestperiodlistelement {
    int indexhigh;
    int indexlow;
    int best1;
    int best2;
    int best3;
    int best4;
    int best5;
    struct bestperiodlistelement *next;
} Bestperiodlist[1];

struct distanceentry {
    int location;
    int size;
};

struct distancelist {
    int k_run_sums_criteria, waiting_time_criteria, lo_d_range, hi_d_range;
    int numentries, nummatches;
    int lowindex, highindex;
    int linked;
    int linkdown, linkup;
    struct distanceentry *entry;
} *Distance;

#define Lookratio .4

/* created 5/23/05 G. Benson */

/* this array replaces the distanceseenlist.  It stores the extent of
   alignments, both pre- or post-consensus, and is used to block an
   alignment in the same region with the same distance */

struct distanceseenarrayelement {
    int index;
    int end;
    int score;
} *Distanceseenarray;

struct distancelistelement {
    int index;
    int distance;
    int changed_from_distance;  /* use for test in
                                 * search_for_distance_match_in_distanceseenlist
                                 * 3/10/05 */
    int end;
    int score;
    int best_possible_score;    /* number of copies X length X match weight */
    int accepted;
    struct distancelistelement *next;
} Distanceseenlist[1];

/* returns max of 4 in order a,b,c,d */
static inline int max4( int a, int b, int c, int d)
{
    return ( ( a >= b )
        ? ( ( a >= c ) ? ( ( a >= d ) ? a : d ) : ( ( c >= d ) ? c : d ) ) 
        : ( ( b >= c ) ? ( ( b >= d ) ? b : d ) : ( ( c >= d ) ? c : d ) ) );
}

/* returns max of 3 in order a,b,c */
static inline int max3( int a, int b, int c )
{
    return a >= b ? a >= c ? a : c : b >= c ? b : c;
}

/* This function may be called multiple times (for different match/mismatch scores) */
int *SM = NULL;

// TODO: Replace with something saner, but note that a and b 
// differs in types here and there. char or unsigned char.
// boa@20240726
#define match( a, b ) ( SM[256 * ( ( a ) ) + ( b )] )

#define fill_align_pair( c1, c2, l, i, j ) \
    AlignPair.textprime[l]  = c1;          \
    AlignPair.textsecnd[l]  = c2;          \
    AlignPair.indexprime[l] = i;           \
    AlignPair.indexsecnd[l] = j

#define max( a, b ) ( ( ( a ) >= ( b ) ) ? ( a ) : ( b ) )
#define min( a, b ) ( ( ( a ) <= ( b ) ) ? ( a ) : ( b ) )

double Copynumber;
double WDPcount;
double OUTPUTcount;
int *Criteria_count;
int *Consensus_count;
int *Outputsize_count;
double *Cell_count;
double Try_waiting_time_count, Fail_waiting_time_count;
double Cell_total, Wasted_total;

/**********************************************************************/

/* New to 2A */

#define GLOBAL 0
#define LOCAL 1

struct paramset {
    unsigned int ps_match;
    unsigned int ps_mismatch;
    unsigned int ps_indel;
    unsigned int ps_minscore;
    unsigned int ps_maxperiod;
    unsigned int ps_PM;
    unsigned int ps_PI;
    int ps_datafile;
    int ps_maskedfile;
    int ps_flankingsequence;
    unsigned int ps_flankinglength;
    int ps_HTMLoff;
    int ps_redundoff;
    int ps_ngs;
    int ps_use_stdin;
    unsigned int ps_maxwraplength;

    char ps_inputfilename[_MAX_PATH];  /* constant defined in stdlib */
    char ps_outputprefix[_MAX_PATH];
    char ps_outputdirectory[_MAX_PATH];
    char ps_outputfilename[_MAX_PATH];
    int ps_multisequencefile;      /* flags if file has more than one sequence */
    int ps_sequenceordinal;        /* holds seq. index starting on 1 */
    int ps_outputcount;            /* repeats found */
    int ps_running;
    char *ps_endstatus;
    int ps_percent;
};

struct paramset paramset;           /* this global controls the algorithm */

/* change MAXWRAPLENGTH to MAXWRAPLENGTHCONST so MAXWRAPLENGTH can be used as an int */

/* int Bandcenter[MAXWRAPLENGTH+1]; */
int *Bandcenter = NULL;

/* version 2A changes this */

/* int S[MAXWRAPLENGTH+1][MAXDISTANCE+1];*/

/* int Up[MAXDISTANCE+1], Diag[MAXDISTANCE+1];*/

/* int S[MAXWRAPLENGTH+1][MAXBANDWIDTH+1];*/
int **S;

int Up[MAXBANDWIDTH + 1], Diag[MAXBANDWIDTH + 1];

/* version 2A adds max3 and max2 */
#define max2( a, b ) ( ( a >= b ) ? a : b )

/* returns max of 2 in order a,b */

#define max3( a, b, c ) \
    ( ( a >= b ) ? ( ( a >= c ) ? a : c ) : ( ( b >= c ) ? b : c ) )

/* returns max of 3 in order a,b,c */

int Maxrealcol;

/* new for 2Anewt */

static const int four_to_the[] = {
    1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576
};

/* number of different tuple sizes to use; preset for all distances */
int NTS;

/* int Tuplesize[NTS+1]={0,2,3,5,7};*//* what the different sizes are */

int Tuplesize[MAXTUPLESIZES + 1];
int Tuplemaxdistance[MAXTUPLESIZES + 1];

/* this is where the actual tuple codes encountered
 * at a sequence location * are stored. */
int Tuplecode[MAXTUPLESIZES + 1];   

/* points to last location of code in history list */
int *Tuplehash[MAXTUPLESIZES + 1];  

/* size of history lists */
int Historysize[MAXTUPLESIZES + 1]; 

/* next free location in history index */
int Nextfreehistoryindex[MAXTUPLESIZES + 1];

struct historyentry {
    int location, previous, code;
} *History[MAXTUPLESIZES + 1];

struct distribution_parameters {
    double exp;
    double var;
};

int Criteria_print = 0;
int Meet_criteria_print = 0;
int *Sortmultiples;

/* modified 3/25/05 G. Benson */
#define NUMBER_OF_PERIODS 5                    /* determines 5 best periods for a repeat */
#define NUMBER_OF_PERIODS_TO_TEST 3            /* only test 3 best periods for multiples test */
#define NUMBER_OF_PERIODS_INTO_SORTMULTIPLES 5 /* added 5/25/05 for compatibility with bestperiodslist */

struct MDDtype {
    int distance;
    char *direction;
};                              /* MDD[MAXWRAPLENGTH+1][MAXBANDWIDTH+1]; */

int F[MAXBANDWIDTH + 1], Fdistance[MAXBANDWIDTH + 1];

int ldong;

int *Statistics_Distance;

FILE *Fptxt;
FILE *Fpdat;

int Minsize = 1;
int Minscore;
int MaxPeriod;

int Period;

int print_flanking = 0;

#define CTRL_SUCCESS 0
#define CTRL_BADFNAME -1
#define CTRL_BADFORMAT -2
#define CTRL_NOTHINGPROCESSED -3

/* the following structure is used to pass a sequence to the algorithm */
#define MAXSEQNAMELEN 200

struct fastasequence {
    /* Changed to unsigned Feb 16, 2016 Yozen */
    unsigned int length;
    int composition[26];
    int nucleotides;
    char name[MAXSEQNAMELEN];
    char *sequence;

};

void trf_message(char *format, ...);

char* newAlignPairtext(int length);
char* newLine(int length);
int* newAlignPairindex(int length);
int* newTags(int length);
void init_bestperiodlist(void);
struct distancelist *new_distancelist(void);
void clear_distancelist(struct distancelist *objptr);
void init_links(void);
void init_index(void);
void init_distanceseenarray(void);
void init_and_fill_coin_toss_stats2000_with_4tuplesizes(void);
void newtupbo(void);
int d_range(int d);
void free_distanceseenarray(void);
void distanceentry_free(void);
void free_bestperiodlist(void);

#endif
