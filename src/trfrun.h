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


#ifndef TRFRUN_H
#define TRFRUN_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

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

struct index_list *GlobalIndexList = NULL;
struct index_list *GlobalIndexListTail = NULL;

void FreeList(struct index_list * headptr);

#include "tr30dat.h"
#include "trfclean.h"

#ifndef _MAX_PATH
#define _MAX_PATH 260
#endif

// how much memory to allocate at first when loading a sequence
// (added by Eugene Scherba, 2010-02-16)
#define MEMBLOCK (10*1024*1024)

int LoadSequenceFromFileEugene(struct fastasequence * pseq, FILE * fp);    /* may use stdin and better file reading for handling large files over 2GB */
int LoadSequenceFromFileBenson(struct fastasequence * pseq, FILE * fp);    /* old function, uses filepos, 32bit version of this would not process a file over 2GB properly */
void TRFControlRoutine(void);
void TRF(struct fastasequence * pseq);
void PrintError(char *errortext);
void PrintProgress(char *progresstext);
void SetProgressBar(void);

#endif
