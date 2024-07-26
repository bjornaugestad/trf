#ifndef TRFCLEAN_H
#define TRFCLEAN_H

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

/****************************************************************
 *
 *   TRFCLEAN.H : Routines to remove redundancy from output files
 *
 *****************************************************************
 *   the following routine are used to remove redundancy from the
 *   alignment file as well as to produce the HTML otput for the
 *   sumary table file(s) and the optional data and masked files.
 *   The general procedure is the following :
 *   i.      Use data from .dat file to build a linked of records.
 *   ii.     Remove pattern sizes greater than MPS.
 *   iii.    Get List sorted by first index keeping track of
 *           original order.
 *   iv.     Remove redundant elements based on given criteria.
 *   v.      Sort remaining elements back to original order.
 *   vi.     In the Alignment file remove sections that correspond
 *           to records deleted in ii and iv.
 *   vii.    Break alignment file in a many files including a
 *           maximun number of alignments per file.
 *   vii.    Save all remaining elements in HTML format. Creating
 *           as many files as necesary.
 *   viii.   Delete old data file and make new one if option is
 *           set.
 *   ix.     Produce Masked file if option is set.
 *   x.      Free linked list memory.
 *   Note that the sorting of elements is made to make the
 *   different activities optimal (i.e. removing redundancy and
 *   removing from alignment file.)
 *
 *
 * Last updated Dec 14,2004
 *****************************************************************/


extern void PrintError(const char *errortext);    /* defined in trfrun.h */

/* Global strings to store non-tabulated information in html file */
char hsequence[256];
char hparameters[256];
char hlength[256];

/* max # of items in tables for extended output format*/
#define EO_MAX_TBL 120

struct index_list *GetList(const char *datafile);
struct index_list *RemoveBySize(struct index_list * headptr, int maxsize);
struct index_list *SortByIndex(struct index_list * headptr);
struct index_list *RemoveRedundancy(struct index_list * headptr);
struct index_list *SortByCount(struct index_list * headptr);
void CleanAlignments(struct index_list * headptr, const char *alignmentfile);
void BreakAlignments(struct index_list * headptr, const char *alignmentfile);
void OutputHTML(struct index_list * headptr, const char *tablefile, const char *alignmentfile);
void MakeDataFile(struct index_list * headptr, const char *datafile, int data);
void MakeMaskedFile(struct index_list * headptr, int masked, unsigned char *Sequence, const char *maskfile);

void FreeList(struct index_list * headptr);

int IntervalOverlap(struct index_list * iptr, struct index_list * jptr);
int IsRedundant(struct index_list * iptr, struct index_list * jptr);
void MakeFileName(char *newname, const char *oldname, int tag);

void TRFClean(const char *alignmentfile, const char *tablefile, int maxsize);


#endif
