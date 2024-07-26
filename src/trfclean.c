
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "trfclean.h"
#include "trfrun.h"
#include "tr30dat.h"
#include "indexlist.h"

void TRFClean(const char *alignmentfile, const char *tablefile, int maxsize)
{
    struct index_list *headptr = NULL, *currptr;
    int i;

    /* To have smaller sequences not send results */
    /* to disc to improve performance             */
    headptr = GlobalIndexList;
    headptr = RemoveBySize(headptr, maxsize);
    headptr = SortByIndex(headptr);

    if (!g_paramset.ps_redundoff) 
        headptr = RemoveRedundancy(headptr);

    headptr = SortByCount(headptr);

    if (!g_paramset.ps_HTMLoff) {
        CleanAlignments(headptr, alignmentfile);
        BreakAlignments(headptr, alignmentfile);
        OutputHTML(headptr, tablefile, alignmentfile);
    }

    /* update the global result */
    for (i = 0, currptr = headptr; currptr != NULL; i++, currptr = currptr->il_next)
        ;

    g_paramset.ps_outputcount = i;

    /* To have smaller sequences not send results */
    /* to disc to improve performance             */
    GlobalIndexList = headptr;
}

struct index_list *GetList(const char *datafile)
{
    FILE *fp;
    struct index_list *headptr, *newptr, *lastptr;
    int counter, i;
    char patbuffer[MAXDISTANCECONSTANT + 1000];

    headptr = newptr = lastptr = NULL;

    /* open file */
    fp = fopen(datafile, "r");
    if (NULL == fp)
        return NULL;

    /* get hsequence line on ninth line of data file */
    for (counter = 0; counter < 9; counter++)
        fgets(hsequence, 255, fp);

    /* get hparameters line on twelvth line of data file */
    for (counter = 0; counter < 3; counter++)
        fgets(hparameters, 255, fp);

    /* get hlength from another global variable (bad practice) */
    sprintf(hlength, "Length:  %d", Length);

    /* loop to fill out list from buffer */
    counter = 1;                /* keeps track of order they are found */
    while (1) {
        /* create new index list element */
        newptr = malloc(sizeof *newptr);
        if (newptr == NULL) {
            FreeList(headptr);
            return NULL;
        }
        newptr->il_count = counter++;

        /* get data from file */
        i = fscanf(fp, "%s %d %d %d %f %d %d %d %d %d %d %d %d %f %s",
            newptr->il_ref, &newptr->il_first, &newptr->il_last,
            &newptr->il_period, &newptr->il_copies, &newptr->il_size,
            &newptr->il_matches, &newptr->il_indels, &newptr->il_score,
            &newptr->il_acount, &newptr->il_ccount, &newptr->il_gcount, &newptr->il_tcount, &newptr->il_entropy, patbuffer);

        if (i == EOF) {
            free(newptr);
            break;
        }

        /* allocate memory to place the pattern and copy data into it */
        newptr->il_pattern = malloc(strlen(patbuffer) + 1);
        if (newptr->il_pattern == NULL) {
            free(newptr);
            FreeList(headptr);
            return NULL;
        }
        strcpy(newptr->il_pattern, patbuffer);

        if (headptr == NULL) {  /* first element */
            headptr = lastptr = newptr;
            lastptr->il_next = NULL;
        }
        else {                  /* add new element to end of list */
            lastptr->il_next = newptr;
            lastptr = newptr;
            lastptr->il_next = NULL;
        }
    }

    fclose(fp);
    return headptr;
}

/************ RemoveBySize() **********************************************/

struct index_list *RemoveBySize(struct index_list *headptr, int maxsize)
{
    struct index_list *currptr;
    struct index_list *prevptr;

    /* loop thru list removing all elements with period > maxsize */
    for (currptr = headptr; currptr != NULL;) {
        if (currptr->il_period > maxsize) {    /* remove */
            if (currptr == headptr) {
                headptr = headptr->il_next;
                free(currptr->il_pattern);
                free(currptr);
                currptr = headptr;
            }
            else {
                prevptr->il_next = currptr->il_next;
                free(currptr->il_pattern);
                free(currptr);
                currptr = prevptr->il_next;
            }

        }
        else {
            prevptr = currptr;
            currptr = currptr->il_next;
        }
    }

    return headptr;

}

struct index_list *SortByIndex(struct index_list *headptr)
{
    struct index_list *currptr;
    struct index_list *holdptr;
    struct index_list *prevptr;
    int dif;                    /* flags when changes occur in one pass */

    if (headptr == NULL)
        return headptr;         /* return if no elements */
    if (headptr->il_next == NULL)
        return headptr;         /* return if one element only */

    dif = 1;
    currptr = headptr;

    while (dif) {
        dif = 0;
        /* repeat inner loop until end is reached */
        while (currptr->il_next != NULL) {
            if (currptr->il_first > currptr->il_next->il_first) {    /* swap */
                if (currptr == headptr) {
                    holdptr = currptr->il_next->il_next;
                    headptr = currptr->il_next;
                    headptr->il_next = currptr;
                    currptr->il_next = holdptr;
                    prevptr = headptr;

                }
                else {
                    prevptr->il_next = currptr->il_next;
                    holdptr = currptr->il_next->il_next;
                    prevptr->il_next->il_next = currptr;
                    currptr->il_next = holdptr;
                    prevptr = prevptr->il_next;
                }

                dif = 1;        /* mark as changed */
            }
            else {
                prevptr = currptr;
                currptr = currptr->il_next;
            }

        }
        currptr = headptr;      /* restart from begining */
    }

    return headptr;
}

struct index_list *RemoveRedundancy(struct index_list *headptr)
{
    struct index_list *iptr;                   /* first pointer of a pair being examined */
    struct index_list *jptr;                   /* second pointer of a pair being examined */
    struct index_list *previptr;               /* points to element before iptr */
    struct index_list *prevjptr;               /* points to element before jptr */
    int overlap;                /* overlap of two intervals */
    int iinterval;              /* first pointer index interval */
    int jinterval;              /* second pointer index interval */

    if (headptr == NULL)
        return headptr;         /* return if no elements */

    if (headptr->il_next == NULL)
        return headptr;         /* return if one element only */

    iptr = headptr;             /* initialize to start at head of list */

    while (iptr != NULL) {      /* loop from start to end of list */
        //iinterval = iptr->last - iptr->first + 1;

        jptr = iptr->il_next;
        prevjptr = iptr;
        while (jptr != NULL) {  /* loop until not enough overlap */
            iinterval = iptr->il_last - iptr->il_first + 1;
            jinterval = jptr->il_last - jptr->il_first + 1;
            overlap = IntervalOverlap(iptr, jptr);
            if (overlap == 0)
                break;

            /* if neither overlap satisfies minimum of 90% break loop */
            /*if( (overlap/(double)iinterval<0.9) && */
            /*  (overlap/(double)jinterval<0.9) ) break; */

            /* overlap in relation to iinterval */
            if (!(overlap / (double)iinterval < 0.9) && IsRedundant(iptr, jptr)) {
                /* remove iptr and break from inner loop */
                if (iptr == headptr) {
                    headptr = headptr->il_next;
                    free(iptr->il_pattern);
                    free(iptr);
                    iptr = headptr;
                    jptr = iptr->il_next;
                    prevjptr = iptr;
                    continue;
                }
                else {
                    previptr->il_next = iptr->il_next;
                    free(iptr->il_pattern);
                    free(iptr);
                    iptr = previptr->il_next;
                    jptr = iptr->il_next;
                    prevjptr = iptr;
                    continue;
                }
            }

            /* overlap in relation to jinterval */
            if (!(overlap / (double)jinterval < 0.9) && IsRedundant(jptr, iptr)) {
                /* remove jptr and continue next iteration of inner loop */
                prevjptr->il_next = jptr->il_next;
                free(jptr->il_pattern);
                free(jptr);
                jptr = prevjptr->il_next;
                continue;
            }

            /* update and continue to next iteration of inner loop */
            prevjptr = jptr;
            jptr = jptr->il_next;

        }
        previptr = iptr;
        iptr = iptr->il_next;
    }
    return headptr;
}

int IntervalOverlap(struct index_list *iptr, struct index_list *jptr)
{
    int beg, end, overlap;

    beg = (iptr->il_first > jptr->il_first) ? iptr->il_first : jptr->il_first;
    end = (iptr->il_last < jptr->il_last) ? iptr->il_last : jptr->il_last;

    overlap = end - beg + 1;
    return ((overlap > 0) ? overlap : 0);
}

/* this fuction retuns 1 if iptr is redundant in reference to jptr */
int IsRedundant(struct index_list *iptr, struct index_list *jptr)
{

    if ((iptr->il_period > jptr->il_period) && (iptr->il_period % jptr->il_period == 0) && (iptr->il_score <= 1.1 * jptr->il_score))
        return 1;
    if ((iptr->il_period == jptr->il_period) && (iptr->il_score <= jptr->il_score))
        return 1;

    return 0;
}

struct index_list *SortByCount(struct index_list *headptr)
{
    struct index_list *currptr;
    struct index_list *holdptr;
    struct index_list *prevptr;
    int dif;                    /* flags when changes occur in one pass */

    if (headptr == NULL)
        return headptr;         /* return if no elements */
    if (headptr->il_next == NULL)
        return headptr;         /* return if one element only */

    dif = 1;
    currptr = headptr;

    while (dif) {
        dif = 0;
        /* repeat inner loop until end is reached */
        while (currptr->il_next != NULL) {
            if (currptr->il_count > currptr->il_next->il_count) {    /* swap */
                if (currptr == headptr) {
                    holdptr = currptr->il_next->il_next;
                    headptr = currptr->il_next;
                    headptr->il_next = currptr;
                    currptr->il_next = holdptr;
                    prevptr = headptr;

                }
                else {
                    prevptr->il_next = currptr->il_next;
                    holdptr = currptr->il_next->il_next;
                    prevptr->il_next->il_next = currptr;
                    currptr->il_next = holdptr;
                    prevptr = prevptr->il_next;
                }

                dif = 1;        /* mark as changed */
            }
            else {
                prevptr = currptr;
                currptr = currptr->il_next;
            }

        }
        currptr = headptr;      /* restart from begining */
    }

    return headptr;
}

void CleanAlignments(struct index_list *headptr, const char *alignmentfile)
{
    char string1[260];
    char string2[260];
    char string3[260];
    char tempfile[264];
    FILE *al_fp, *tmp_fp;
    struct index_list *currptr;
    int moving;                 /* in loop moving=1 to copy strings to tempfile */

    /* execute even if no repeats in list */

    /* make name of temporary file */
    strcpy(tempfile, alignmentfile);
    strcat(tempfile, ".tmp");

    /* open files */
    al_fp = fopen(alignmentfile, "r");
    if (al_fp == NULL)
        die("Unable to open alignment file for reading in CleanAlignments routine!");
    
    tmp_fp = fopen(tempfile, "w");
    if (tmp_fp == NULL)
        die("Unable to open temp file for writing in CleanAlignments routine!");

    moving = 1;                 /* starts by moving lines to temporary file */
    currptr = headptr;
    while (fgets(string1, 260, al_fp) != NULL) {
        if (string1[0] == 'F' /* Gelfand Dec 15, 2013 */  && string1[1] == 'o') {   /* Ocurrs when "Found at" is encountered */
            fgets(string2, 260, al_fp);
            fgets(string3, 260, al_fp);
            if (currptr != NULL) {
                if (strstr(string3, currptr->il_ref) != NULL) {
                    fputs(string1, tmp_fp);
                    fputs(string2, tmp_fp);
                    fputs(string3, tmp_fp);
                    currptr = currptr->il_next;
                    moving = 1;
                    continue;
                }
            }
            moving = 0;
        }

        if (string1[0] == 'D' /* Gelfand Dec 15, 2013 */  && string1[1] == 'o') {   /* Occurs when "Done." is encountered */
            fgets(string2, 260, al_fp);
            fputs(string1, tmp_fp);
            fputs(string2, tmp_fp);
            fputs("\n", tmp_fp);
            break;
        }

        if (moving)
            fputs(string1, tmp_fp);
    }

    fclose(al_fp);
    fclose(tmp_fp);

    remove(alignmentfile);
    rename(tempfile, alignmentfile);

    return;

}

void BreakAlignments(struct index_list *headptr, const char *alignmentfile)
{
    FILE *al_fp;
    FILE *out_fp;
    char outfile[260];
    int nfiles;
    struct index_list *currptr;
    int alignments;
    char headlns[30][200];      /* to hold the heading section of the file */
    char buffer[200];           /* holds one line of alignment data */
    int headcnt;
    char nextchar;
    int i, j;

    /* Find out how many alignments there are and how many files will
     * be needed */
    for (alignments = 0, currptr = headptr; currptr != NULL; currptr = currptr->il_next, alignments++);
    nfiles = alignments / EO_MAX_TBL;
    if ((alignments % EO_MAX_TBL) > 0)
        nfiles++;
    if (nfiles == 0)
        nfiles = 1;             /* make sure at least one file is generated */

    /* if only one file will be needed just rename file */
    if (nfiles == 1) {
        MakeFileName(outfile, alignmentfile, 1);
        remove(outfile);
        rename(alignmentfile, outfile);
        return;
    }

    /* get heading lines */
    al_fp = fopen(alignmentfile, "r");
    if (al_fp == NULL)
        die("Unable to open alignment file for reading in BreakAlignments routine!");

    for (i = 0, headcnt = 0; i < 30; i++) {
        /* find if next line is the "Found at" line */
        nextchar = getc(al_fp);
        ungetc(nextchar, al_fp);
        if (nextchar == 'F')
            break;
        fgets(headlns[i], 200, al_fp);
        headcnt++;
    }

    /* loop creating files */
    for (i = 1; i <= nfiles; i++) {
        /* create name of ith file */
        MakeFileName(outfile, alignmentfile, i);
        /* open the file for writing */
        out_fp = fopen(outfile, "w");
        if (out_fp == NULL)
            die("Unable to open output file for writing in BreakAlignments routine!");

        /* output heading */
        for (j = 0; j < headcnt; j++)
            fputs(headlns[j], out_fp);
        /* output n of N identifier */
        //fprintf(out_fp,"File %d of %d\n\n",i+1,nfiles);
        fprintf(out_fp, "File %d of %d\n\n", i, nfiles);

        /*move alignments to the current output file */
        for (j = 0; j < EO_MAX_TBL; j++) {
            /* copy the first line of the alignment */
            fgets(buffer, 200, al_fp);
            fputs(buffer, out_fp);
            /* copy successive lines checking for the "Found at" and
             * the "Done." lines */
            while (1) {
                nextchar = getc(al_fp);
                ungetc(nextchar, al_fp);
                if (nextchar == 'F' || nextchar == 'D')
                    break;

                // added check for NULL, Gelfand Dec 15, 2013 because of suspicion this may under
                // certain conditions enter an infinite loop.
                if (NULL == fgets(buffer, 200, al_fp)) {
                    nextchar = 'D';
                    break;
                }

                fputs(buffer, out_fp);
            }
            if (nextchar == 'F')
                continue;
            if (nextchar == 'D')
                break;
        }

        /* Output closing Lines */
        fprintf(out_fp, "\nDone.\n</PRE></BODY></HTML>\n");

        /*close file */
        fclose(out_fp);
    }

    /*close original file */
    fclose(al_fp);

    /* delete original file */
    remove(alignmentfile);

    return;
}

void MakeFileName(char *newname, const char *oldname, int tag)
{
    char newext[20];
    char oldext[10];
    int numindex;
    int count;

    /*copy oldname to newname buffer */
    strcpy(newname, oldname);

    /*find position of last number in name */
    for (count = 0, numindex = 0; newname[count] != '\0'; count++) {
        if (isdigit((int)newname[count]))
            numindex = count;
    }

    /*get old extension */
    strcpy(oldext, &newname[numindex + 1]);

    /* create new extension based on tag */
    sprintf(newext, ".%d%s", tag, oldext);

    /* copy newext in place over old newname */
    strcpy(&newname[numindex + 1], newext);

    return;
}

static void OutputHeading(FILE *fp, const char *tablefile, const char *alignmentfile)
{
    /* output fixed (old) heading */
    fprintf(fp,
        "<HTML><HEAD><TITLE>%s</TITLE><BASE TARGET=\"%s\"></HEAD><BODY bgcolor=\"#FBF8BC\"><BR><PRE>Tandem Repeats Finder Program written by:</PRE><PRE><CENTER>Gary Benson<BR>Program in Bioinformatics<BR>Boston University<BR>Version %s<BR></CENTER>",
        tablefile, alignmentfile, versionstring);
    fprintf(fp,
        "\nPlease cite:\nG. Benson,\n\"Tandem repeats finder: a program to analyze DNA sequences\"\nNucleic Acid Research(1999)\nVol. 27, No. 2, pp. 573-580.\n");

    fprintf(fp, "\n%s", hsequence);
    fprintf(fp, "%s", hparameters);
    fprintf(fp, "%s</PRE>\n", hlength);
}

void OutputHTML(struct index_list *headptr, const char *tablefile, const char *alignmentfile)
{
    FILE *fp;
    struct index_list *currptr;
    int i, j;
    int alignments;
    int nfiles;
    char outfile[260];
    char linkfile[260];
    char namebuffer[260];

    /* find out how many elements there are and how
     * many files will be needed */
    for (alignments = 0, currptr = headptr; currptr != NULL; currptr = currptr->il_next, alignments++);
    nfiles = alignments / EO_MAX_TBL;
    if ((alignments % EO_MAX_TBL) > 0)
        nfiles++;
    if (nfiles == 0)
        nfiles = 1;             /* make sure at least one file is generated */

    /* loop creating files */
    currptr = headptr;
    for (i = 1; i <= nfiles; i++) {
        /* create name of ith file */
        MakeFileName(outfile, tablefile, i);

        /* create name of file where the alignments are found */
        MakeFileName(linkfile, alignmentfile, i);

        /* open the file for writing */
        fp = fopen(outfile, "w");
        if (fp == NULL)
            die("Unable to open output file for writing in OutputHTML routine!");

        /* output heading */
        OutputHeading(fp, outfile, alignmentfile);

        /* print links to other tables */
        fprintf(fp, "\n<P><PRE>Tables:   ");
        for (j = 1; j <= nfiles; j++) {
            /* output a link for all but the current page */
            if (j != i) {
                MakeFileName(namebuffer, tablefile, j);
                fprintf(fp, "<A HREF=\"%s\" target=\"_self\">%d</A>   ", namebuffer, j);
            }
            else {
                fprintf(fp, "%d   ", j);
            }
            if (j % 16 == 0 && j < nfiles)
                fprintf(fp, "\n          ");
        }

        /* printf "Table n of N" line */
        fprintf(fp, "\n\nThis is table  %d  of  %d  ( %d repeats found )\n</PRE>", i, nfiles, alignments);

        /*print help lines */
        fprintf(fp, "<PRE>\nClick on indices to view alignment\n");

        fprintf(fp,
            "</PRE><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n");

        /* print beginning of table */
        fprintf(fp, "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");

        /* output rows of data */
        for (j = 0; j < EO_MAX_TBL && currptr != NULL; currptr = currptr->il_next, j++) {

            if (j % 22 == 0)
                fprintf(fp,
                    "<TR><TD WIDTH=140><CENTER>Indices</CENTER></TD><TD WIDTH=80><CENTER>Period<BR>Size </CENTER></TD><TD WIDTH=70><CENTER>Copy<BR>Number</CENTER></TD><TD WIDTH=70><CENTER>Consensus<BR>Size</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=60><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>A</CENTER></TD><TD WIDTH=40><CENTER>C</CENTER></TD><TD WIDTH=40><CENTER>G</CENTER></TD><TD WIDTH=40><CENTER>T</CENTER></TD><TD WIDTH=70><CENTER>Entropy<BR>(0-2)</CENTER></TD></TR>\n");
            fprintf(fp,
                "<TR><TD><CENTER><A HREF=\"%s#%s\">%d--%d</A></CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%1.1f</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%1.2f</CENTER></TD></TR>\n",
                linkfile, currptr->il_ref, currptr->il_first, currptr->il_last, currptr->il_period, currptr->il_copies, currptr->il_size,
                currptr->il_matches, currptr->il_indels, currptr->il_score, currptr->il_acount, currptr->il_ccount, currptr->il_gcount,
                currptr->il_tcount, currptr->il_entropy);
        }

        /* if empty list print at least one heading */
        if (headptr == NULL)
            fprintf(fp,
                "<TR><TD WIDTH=140><CENTER>Indices</CENTER></TD><TD WIDTH=80><CENTER>Period<BR>Size </CENTER></TD><TD WIDTH=70><CENTER>Copy<BR>Number</CENTER></TD><TD WIDTH=70><CENTER>Consensus<BR>Size</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=60><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>A</CENTER></TD><TD WIDTH=40><CENTER>C</CENTER></TD><TD WIDTH=40><CENTER>G</CENTER></TD><TD WIDTH=40><CENTER>T</CENTER></TD><TD WIDTH=70><CENTER>Entropy<BR>(0-2)</CENTER></TD></TR>\n");

        /* close table */
        fprintf(fp, "\n</TABLE>\n");

        /* if no repeats print message */
        if (headptr == NULL)
            fprintf(fp, "\nNo Repeats Found!<BR>"); /*for empty list */

        /* print links to other tables (again) */
        fprintf(fp, "\n<P><PRE>Tables:   ");
        for (j = 1; j <= nfiles; j++) {
            /* output a link for all but the current page */
            if (j != i) {
                MakeFileName(namebuffer, tablefile, j);
                fprintf(fp, "<A HREF=\"%s\" target=\"_self\">%d</A>   ", namebuffer, j);
            }
            else {
                fprintf(fp, "%d   ", j);
            }

            if (j % 16 == 0 && j < nfiles)
                fprintf(fp, "\n          ");
        }
        fprintf(fp, "\n</PRE>");

        if (i == nfiles)
            fprintf(fp, "<P>The End!\n");

        fprintf(fp, "\n</BODY></HTML>\n");
        fclose(fp);
    }
}

void MakeDataFile(struct index_list *headptr, const char *datafile, int data)
{
    FILE *fp;
    struct index_list *lpointer;
    int charcount;

    /* if data = 1 then produce new datafile overwriting old one */
    /* otherwise delete old file */
    if (data) {
        fp = fopen(datafile, "w");
        if (fp == NULL)
            die("Unable to open output file for writing in MakeDataFile routine!");

        if (g_paramset.ps_ngs != 1) {
            fprintf(fp,
                "Tandem Repeats Finder Program written by:\n\nGary Benson\nProgram in Bioinformatics\nBoston University\nVersion %s\n\n\n%s\n\n\n%s\n\n",
                versionstring, hsequence, hparameters);
        }

        for (lpointer = headptr; lpointer != NULL; lpointer = lpointer->il_next) {
            fprintf(fp, "%d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s ",
                lpointer->il_first, lpointer->il_last, lpointer->il_period,
                lpointer->il_copies, lpointer->il_size, lpointer->il_matches,
                lpointer->il_indels, lpointer->il_score, lpointer->il_acount,
                lpointer->il_ccount, lpointer->il_gcount, lpointer->il_tcount, lpointer->il_entropy, lpointer->il_pattern);

            for (charcount = lpointer->il_first; charcount <= lpointer->il_last; charcount++)
                fprintf(fp, "%c", Sequence[charcount]);

            fprintf(fp, "\n");
        }

        fclose(fp);
    }
    else
        remove(datafile);
}

void MakeMaskedFile(struct index_list *headptr, int masked, unsigned char *Sequence, const char *maskfile)
{
    int count, printcr;
    int masker;
    FILE *fp;
    struct index_list *lpointer;

    if (masked) {
        fp = fopen(maskfile, "w");
        if (fp == NULL)
            die("Unable to open output file for writing in MakeMaskedFile routine!");

        /* Ouput sequence description from global variable to file */
        fprintf(fp, ">%s", &hsequence[10]);

        for (lpointer = headptr; lpointer != NULL; lpointer = lpointer->il_next) {
            for (masker = lpointer->il_first; masker <= lpointer->il_last; masker++)
                Sequence[masker] = 'N';
        }
        printcr = 0;
        for (count = 1; Sequence[count] != '\0'; count++) {
            fputc(Sequence[count], fp);
            printcr++;
            if (printcr >= 60) {
                printcr = 0;
                fputc('\n', fp);
            }
        }
        fputc('\n', fp);
        fputc('\n', fp);
        fclose(fp);
    }

    return;
}

