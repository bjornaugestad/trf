#ifndef TRF_INDEXLIST_H
#define TRF_INDEXLIST_H

struct index_list {
    int il_count;                  /* indicates order in original file */
    char il_ref[45];               /* records label for linking */
    int il_first;                  /* first index */
    int il_last;                   /* last index */
    int il_period;                 /* period size */
    float il_copies;               /* number of copies */
    int il_size;                   /* consensus size */
    int il_matches;
    int il_indels;
    int il_score;
    int il_acount;
    int il_ccount;
    int il_gcount;
    int il_tcount;
    float il_entropy;
    char *il_pattern;

    struct index_list *il_next;
};

void FreeList(struct index_list *headptr);

#endif

