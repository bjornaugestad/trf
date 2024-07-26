
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

#include <string.h>
#include <math.h>

#include "trfrun.h"
#include "tr30dat.h"

/* This routine can act on a multiple-sequence file
 * and calls TRF() routine as many times as it needs to. */
void TRFControlRoutine(void)
{
    FILE *srcfp, *outmfp, *destmfp, *destdfp = NULL;
    char source[_MAX_PATH], input[_MAX_PATH], outm[_MAX_PATH],
        prefix[_MAX_PATH], destm[_MAX_PATH], destd[_MAX_PATH],
        desth[_MAX_PATH], paramstring[_MAX_PATH], outh[_MAX_PATH];
    int a, i, loadstatus, foundsome = 0;
    unsigned int tr_count = 0;
    char line[1000];
    FILE *desthfp;
    struct fastasequence seq;

    /* save names locally so they can be replaced later */
    strcpy(source, paramset.ps_inputfilename);
    strcpy(prefix, paramset.ps_outputprefix);

    /* open input file for reading */
    if (paramset.ps_use_stdin) {
        srcfp = stdin;
    }
    else {
        srcfp = fopen(source, "rb");
        if (srcfp == NULL) {
            sprintf(line, "Error opening %s: %s\n", source, strerror(errno));
            PrintError(line);
            paramset.ps_endstatus = "Bad filename";
            paramset.ps_running = 0;
            return;
        }
    }

    /* get the first sequence */
    if (paramset.ps_ngs != 1)
        PrintProgress("Loading sequence...");

    loadstatus = LoadSequenceFromFileEugene(&seq, srcfp);
    if (loadstatus < 0) {
        PrintError("Could not load sequence. Empty file or bad format.");
        paramset.ps_endstatus = "Bad format."; /* ok for now */
        paramset.ps_running = 0;
        fclose(srcfp);
        return;
    }

    /* generate the parameter string to be used in file names */
    sprintf(paramstring, "%d.%d.%d.%d.%d.%d.%d",
        paramset.ps_match, paramset.ps_mismatch, paramset.ps_indel,
        paramset.ps_PM, paramset.ps_PI, paramset.ps_minscore, paramset.ps_maxperiod);

    sprintf(hparameters, "Parameters: %d %d %d %d %d %d %d\n",
        paramset.ps_match, paramset.ps_mismatch, paramset.ps_indel,
        paramset.ps_PM, paramset.ps_PI, paramset.ps_minscore, paramset.ps_maxperiod);

    /* based on number of sequences in file use different approach */
    if (loadstatus == 0) {      /* only one sequence in file */
        sprintf(hsequence, "Sequence: %s\n", seq.name);
        sprintf(hlength, "Length:  %d", seq.length);

        paramset.ps_multisequencefile = 0;
        paramset.ps_sequenceordinal = 1;
        /* call trf and return */
        counterInSeq = 0;
        TRF(&seq);

        if (paramset.ps_endstatus)
            return;

        if (paramset.ps_datafile) {
            /* Added by Yevgeniy Gelfand on Jan 27, 2010  */
            /* To have smaller sequences not send results */
            /* to disc to improve performance             */

            if (paramset.ps_ngs) {
                destdfp = stdout;
            }
            else {
                sprintf(destd, "%s.%s.dat", prefix, paramstring);
                destdfp = fopen(destd, "w");
                if (destdfp == NULL) {
                    PrintError("Unable to open data file for writing in TRFControlRoutine routine!");
                    exit(-1);
                }
            }

            {
                struct index_list *lpointer;
                int charcount;

                if (paramset.ps_ngs != 1) {
                    fprintf(destdfp, "Tandem Repeats Finder Program written by:\n\n");
                    fprintf(destdfp, "Gary Benson\n");
                    fprintf(destdfp, "Program in Bioinformatics\n");
                    fprintf(destdfp, "Boston University\n");
                    fprintf(destdfp, "Version %s\n", versionstring);
                }

                if (paramset.ps_ngs) {
                    /* only print if we have at least 1 record */
                    if (NULL != GlobalIndexList) {
                        fprintf(destdfp, "@%s\n", seq.name);
                    }

                }
                else {
                    fprintf(destdfp, "\n\nSequence: %s\n\n\n\nParameters: %d %d %d %d %d %d %d\n\n\n",
                        seq.name, paramset.ps_match, paramset.ps_mismatch, paramset.ps_indel, paramset.ps_PM, paramset.ps_PI,
                        paramset.ps_minscore, paramset.ps_maxperiod);
                }

                for (lpointer = GlobalIndexList; lpointer != NULL; lpointer = lpointer->il_next) {
                    fprintf(destdfp, "%d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s ",
                        lpointer->il_first, lpointer->il_last, lpointer->il_period,
                        lpointer->il_copies, lpointer->il_size, lpointer->il_matches,
                        lpointer->il_indels, lpointer->il_score, lpointer->il_acount,
                        lpointer->il_ccount, lpointer->il_gcount, lpointer->il_tcount, lpointer->il_entropy, lpointer->il_pattern);
                    for (charcount = lpointer->il_first; charcount <= lpointer->il_last; charcount++)
                        fprintf(destdfp, "%c", Sequence[charcount]);

                    /* print short flanks to .dat file */
                    if (paramset.ps_ngs) {
                        int flankstart, flankend;

                        flankstart = lpointer->il_first - 50;
                        flankstart = max(1, flankstart);
                        flankend = lpointer->il_last + 50;
                        flankend = min(Length, flankend);

                        fprintf(destdfp, " ");
                        if (lpointer->il_first == 1) {
                            fprintf(destdfp, ".");
                        }
                        else {
                            for (charcount = flankstart; charcount < lpointer->il_first; charcount++)
                                fprintf(destdfp, "%c", Sequence[charcount]);
                        }

                        fprintf(destdfp, " ");
                        if (lpointer->il_last == Length)
                            fprintf(destdfp, ".");
                        else
                            for (charcount = lpointer->il_last + 1; charcount <= flankend; charcount++)
                                fprintf(destdfp, "%c", Sequence[charcount]);
                    }

                    fprintf(destdfp, "\n");
                    ++tr_count;
                }
            }
        }

        /* masked file moved here so Sequence is not "ruined" by Ns for .dat output */
        {
            char maskstring[_MAX_PATH];

            sprintf(maskstring, "%s.%s.mask", paramset.ps_outputprefix, paramstring);
            MakeMaskedFile(GlobalIndexList, paramset.ps_maskedfile, Sequence, maskstring);
        }

        FreeList(GlobalIndexList);
        GlobalIndexList = NULL;
        GlobalIndexListTail = NULL;

        free(seq.sequence);
        fclose(srcfp);

        if (destdfp) {
            fclose(destdfp);
            destdfp = NULL;
        }

        paramset.ps_endstatus = NULL;
        paramset.ps_running = 0;
        // return CTRL_SUCCESS;
        return;
    }
    paramset.ps_multisequencefile = 1;
    paramset.ps_sequenceordinal = 1;

    /*
     *   if there are more files need to produce sumary-style
     *   output.
     */

    /* generate the parameter string to be used in file names */
    sprintf(paramstring, "%d.%d.%d.%d.%d.%d.%d",
        paramset.ps_match, paramset.ps_mismatch, paramset.ps_indel,
        paramset.ps_PM, paramset.ps_PI, paramset.ps_minscore, paramset.ps_maxperiod);

    /* open sumary table file */
    sprintf(desth, "%s.%s.summary.html", prefix, paramstring);
    if (!paramset.ps_HTMLoff) {
        desthfp = fopen(desth, "w");
        if (desthfp == NULL) {
            PrintError("Unable to open summary file for writing in TRFControlRoutine routine!");
            exit(-1);
        }
    }

    /* open masked file if requested */
    if (paramset.ps_maskedfile) {
        sprintf(destm, "%s.%s.mask", prefix, paramstring);
        destmfp = fopen(destm, "w");
        if (destmfp == NULL) {
            PrintError("Unable to open masked file for writing in TRFControlRoutine routine!");
            exit(-1);
        }
    }

    /* open datafile if requested */
    if (paramset.ps_datafile) {
        if (paramset.ps_ngs) {
            destdfp = stdout;
        }
        else {
            sprintf(destd, "%s.%s.dat", prefix, paramstring);
            destdfp = fopen(destd, "w");
            if (destdfp == NULL) {
                PrintError("Unable to open data file for writing in TRFControlRoutine routine!");
                exit(-1);
            }
        }
    }

    if (!paramset.ps_HTMLoff) {
        /* start output of sumary file */
        fprintf(desthfp, "<HTML>");
        fprintf(desthfp, "<HEAD>");
        fprintf(desthfp, "<TITLE>Output Summary</TITLE>");
        /* fprintf(desthfp,"<BASE TARGET=\"Table\">"); */
        fprintf(desthfp, "</HEAD>");
        fprintf(desthfp, "<BODY bgcolor=\"#FBF8BC\">");
        fprintf(desthfp, "<PRE>");

        fprintf(desthfp, "\nTandem Repeats Finder Program written by:<CENTER>");
        fprintf(desthfp, "\nGary Benson");
        fprintf(desthfp, "\nProgram in Bioinformatics");
        fprintf(desthfp, "\nBoston University");
        fprintf(desthfp, "\nVersion %s</CENTER>\n", versionstring);
        fprintf(desthfp,
            "\nPlease cite:\nG. Benson,\n\"Tandem repeats finder: a program to analyze DNA sequences\"\nNucleic Acid Research(1999)\nVol. 27, No. 2, pp. 573-580.\n");

        fprintf(desthfp, "\n\n<B>Multiple Sequence Summary</B>\n\n");
        fprintf(desthfp, "Only sequences containing repeats are shown!\n\n");
        fprintf(desthfp, "Click on sequence description to view repeat table.\n\n");

        /* print beginning of table */
        fprintf(desthfp, "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");
        fprintf(desthfp, "<TR><TD WIDTH=80><CENTER>Sequence\nIndex</CENTER></TD>"
            "<TD WIDTH=400><CENTER>Sequence\nDescription</CENTER></TD>"
            "<TD WIDTH=80><CENTER>Number of\nRepeats</CENTER></TD>" "</TR>\n");
    }

    // process every sequence in file
    i = 1;
    for (;;) {

        sprintf(hsequence, "Sequence: %s\n", seq.name);
        sprintf(hlength, "Length:  %d", seq.length);

        // set the prefix to be used for naming of output
        sprintf(input, "%s.s%d", prefix, i);
        strcpy(paramset.ps_inputfilename, input);
        strcpy(paramset.ps_outputprefix, input);

        /* call the tandem repeats finder routine */
        counterInSeq = 0;
        TRF(&seq);

        if (paramset.ps_datafile) {

            /* Added by Yevgeniy Gelfand on Jan 27, 2010  */
            /* To have smaller sequences not send results */
            /* to disc to improve performance             */
            {
                struct index_list *lpointer;
                int charcount;

                /* only for the first one write the header */
                if (i == 1) {
                    if (paramset.ps_ngs != 1) {
                        fprintf(destdfp, "Tandem Repeats Finder Program written by:\n\n");
                        fprintf(destdfp, "Gary Benson\n");
                        fprintf(destdfp, "Program in Bioinformatics\n");
                        fprintf(destdfp, "Boston University\n");
                        fprintf(destdfp, "Version %s\n", versionstring);
                    }
                }

                if (paramset.ps_ngs) {
                    /* only print if we have at least 1 record */
                    if (NULL != GlobalIndexList) {
                        fprintf(destdfp, "@%s\n", seq.name);
                    }
                }
                else {
                    fprintf(destdfp, "\n\nSequence: %s\n\n\n\nParameters: %d %d %d %d %d %d %d\n\n\n",
                        seq.name, paramset.ps_match, paramset.ps_mismatch, paramset.ps_indel, paramset.ps_PM, paramset.ps_PI,
                        paramset.ps_minscore, paramset.ps_maxperiod);
                }

                for (lpointer = GlobalIndexList; lpointer != NULL; lpointer = lpointer->il_next) {
                    fprintf(destdfp, "%d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s ",
                        lpointer->il_first, lpointer->il_last, lpointer->il_period,
                        lpointer->il_copies, lpointer->il_size, lpointer->il_matches,
                        lpointer->il_indels, lpointer->il_score, lpointer->il_acount,
                        lpointer->il_ccount, lpointer->il_gcount, lpointer->il_tcount, lpointer->il_entropy, lpointer->il_pattern);
                    for (charcount = lpointer->il_first; charcount <= lpointer->il_last; charcount++)
                        fprintf(destdfp, "%c", Sequence[charcount]);

                    /* print short flanks to .dat file */
                    if (paramset.ps_ngs) {
                        int flankstart, flankend;

                        flankstart = lpointer->il_first - 50;
                        flankstart = max(1, flankstart);
                        flankend = lpointer->il_last + 50;
                        flankend = min(Length, flankend);

                        fprintf(destdfp, " ");
                        if (lpointer->il_first == 1) {
                            fprintf(destdfp, ".");
                        }
                        else {
                            for (charcount = flankstart; charcount < lpointer->il_first; charcount++)
                                fprintf(destdfp, "%c", Sequence[charcount]);
                        }

                        fprintf(destdfp, " ");
                        if (lpointer->il_last == Length) {
                            fprintf(destdfp, ".");
                        }
                        else {
                            for (charcount = lpointer->il_last + 1; charcount <= flankend; charcount++)
                                fprintf(destdfp, "%c", Sequence[charcount]);
                        }
                    }

                    fprintf(destdfp, "\n");
                }
            }
        }

        if (!paramset.ps_HTMLoff) {
            /* print table rows based on repeat count */
            sprintf(outh, "%s.%s.1.html", input, paramstring);
            if (paramset.ps_outputcount > 0) {
                /* print a table raw to the summary table */
                fprintf(desthfp, "<TR><TD><CENTER>%d</CENTER></TD>"
                    "<TD><CENTER><A TARGET=\"%s\" HREF=\"%s\">%s</A>"
                    "</CENTER></TD><TD><CENTER>%d</CENTER></TD></TR>", i, outh, outh, seq.name, paramset.ps_outputcount);
                foundsome = 1;
            }
            else {
                /* remove html files if no output in it */
                remove(outh);
                sprintf(line, "%s.%s.1.txt.html", input, paramstring);
                remove(line);
            }
        }

        /* masked file moved here so Sequence is not "ruined" by Ns for .dat output */
        {
            char maskstring[_MAX_PATH];

            sprintf(maskstring, "%s.s%d.%s.mask", prefix, i, paramstring);
            MakeMaskedFile(GlobalIndexList, paramset.ps_maskedfile, Sequence, maskstring);
        }

        /* append new output to destination files */
        if (paramset.ps_maskedfile) {
            /* recreate the name of the masked sequence file */
            sprintf(outm, "%s.s%d.%s.mask", prefix, i, paramstring);
            outmfp = fopen(outm, "r");
            if (outmfp == NULL) {
                PrintError("Unable to open masked file for reading in TRFControlRoutine routine!");
                exit(-1);
            }
            /* copy until end of file */
            while (1) {
                a = getc(outmfp);
                if (a == EOF)
                    break;
                putc(a, destmfp);
            }
            fclose(outmfp);

            /* remove intermediary file */
            remove(outm);
        }

        FreeList(GlobalIndexList);
        GlobalIndexList = NULL;
        GlobalIndexListTail = NULL;

        /* free the data associated with the sequence */
        free(seq.sequence);

        /* if more sequences load and repeat */
        if (loadstatus > 0) {
            if (paramset.ps_ngs != 1)
                PrintProgress("Loading sequence file...");

            loadstatus = LoadSequenceFromFileEugene(&seq, srcfp);
            paramset.ps_sequenceordinal++;
            i++;
        }
        else {
            break;
        }
    }

    if (!paramset.ps_HTMLoff) {
        /* close table and html body */
        fprintf(desthfp, "\n</TABLE>\n");
        if (!foundsome) {
            fprintf(desthfp, "\nNo Repeats Found!<BR>");
        }
        fprintf(desthfp, "\n</BODY></HTML>\n");
    }

    /* close files */
    fclose(srcfp);
    if (paramset.ps_maskedfile)
        fclose(destmfp);

    if (paramset.ps_datafile)
        fclose(destdfp);

    if (!paramset.ps_HTMLoff)
        fclose(desthfp);

    /* set output file name to the summary table */
    strcpy(paramset.ps_outputfilename, desth);

    paramset.ps_endstatus = NULL;
    paramset.ps_running = 0;
    return;
}

/*
 *   This routine acts on single-sequence files and
 *   is used by the control routine above.
 */
void TRF(struct fastasequence *pseq)
{
    unsigned int i;             /* used at the end to free memory */
    int *stemp;
    char htmlstring[_MAX_PATH], txtstring[_MAX_PATH],
        paramstring[_MAX_PATH], datstring[_MAX_PATH], maskstring[_MAX_PATH], messagebuffer[100];

    init_bestperiodlist();

    /*  Set global print_flanking that controls the generation of flanking */
    print_flanking = paramset.ps_flankingsequence;

    /* allocate memory for file names */
    if (paramset.ps_ngs != 1)
        PrintProgress("Allocating Memory...");

    /* change made for NGS data analysis */
    /* make MAXWRAPLENGTH = 1000 for smaller for small sequences */
    maxwraplength = min(paramset.ps_maxwraplength, pseq->length);

    /* allocate memory */
    S = (int **)malloc((maxwraplength + 1) * sizeof(int *));
    if (S == NULL) {
        PrintError("Unable to allocate memory for S array");
        exit(-1);
    }

    /* Yozen Jan 26, 2016: We control the compilation and we're going to be using C99
     * or greater standard C; don't need to cast, and we can use the pointer to determine the size.
     * Also, use calloc instead of malloc+memset. */
    stemp = calloc(((size_t)(maxwraplength + 1) * (MAXBANDWIDTH + 1)), sizeof *stemp);
    if (stemp == NULL) {
        char errmsg[255];

        snprintf(errmsg, 255,
            "Unable to allocate %lu bytes for stemp array. Please set a lower value for the longest TR length. (%s:%d)\n",
            ((maxwraplength + 1) * (MAXBANDWIDTH + 1)) * sizeof(*stemp), __FILE__, __LINE__);

        PrintError(errmsg);
        exit(-1);
    }
    for (i = 0; i <= maxwraplength; i++) {
        S[i] = stemp;
        stemp += MAXBANDWIDTH + 1;
    }
    S[0][0] = 1;

    /* AlignPair holds the characters and alignments of the current */
    /* primary and secondary sequences  */
    AlignPair.textprime = newAlignPairtext(2 * maxwraplength);
    if (paramset.ps_endstatus)
        return;

    AlignPair.textsecnd = newAlignPairtext(2 * maxwraplength);
    if (paramset.ps_endstatus)
        return;

    AlignPair.indexprime = newAlignPairindex(2 * maxwraplength);
    if (paramset.ps_endstatus)
        return;
    
    AlignPair.indexsecnd = newAlignPairindex(2 * maxwraplength);
    if (paramset.ps_endstatus) 
        return;

    /* set algorithm's parameters */
    Alpha = paramset.ps_match;
    Beta = -paramset.ps_mismatch;
    Delta = -paramset.ps_indel;
    PM = paramset.ps_PM;
    PI = paramset.ps_PI;
    Minscore = paramset.ps_minscore;
    MaxPeriod = paramset.ps_maxperiod;
    MAXDISTANCE = MAXPATTERNSIZE = paramset.ps_maxperiod;
    if (MAXDISTANCE < 500) {
        MAXDISTANCE = MAXPATTERNSIZE = 500;
    }

    MAXDISTANCE = MAXPATTERNSIZE = min(MAXDISTANCE, (int)(pseq->length * .6));
    MAXDISTANCE = MAXPATTERNSIZE = max(MAXDISTANCE, 200);

    /* generate the parameter string to be used in file names */
    sprintf(paramstring, "%d.%d.%d.%d.%d.%d.%d",
        paramset.ps_match, paramset.ps_mismatch, paramset.ps_indel,
        paramset.ps_PM, paramset.ps_PI, paramset.ps_minscore, paramset.ps_maxperiod);

    Reportmin = 0;
    ldong = 0;
    Rows = 0;
    Totalcharacters = 0;
    Test = 1;

    /* print the names of the files */
    snprintf(htmlstring, sizeof htmlstring, "%s.%s.html", paramset.ps_outputprefix, paramstring);
    snprintf(txtstring, sizeof txtstring, "%s.%s.txt.html", paramset.ps_outputprefix, paramstring);
    snprintf(datstring, sizeof datstring, "%s.%s.dat", paramset.ps_outputprefix, paramstring);
    snprintf(maskstring, sizeof maskstring, "%s.%s.mask", paramset.ps_outputprefix, paramstring);

    /* start txt file */
    if (!paramset.ps_HTMLoff) {
        Fptxt = fopen(txtstring, "w");
        if (Fptxt == NULL) {
            PrintError("Unable to open alignment file for writing in TRF() routine!");
            exit(-1);
        }

        fprintf(Fptxt, "<HTML>");
        fprintf(Fptxt, "<HEAD>");
        fprintf(Fptxt, "<TITLE>%s</TITLE>", txtstring);
        fprintf(Fptxt, "</HEAD>");
        fprintf(Fptxt, "<BODY bgcolor=\"#FBF8BC\">");
        fprintf(Fptxt, "<PRE>");

        fprintf(Fptxt, "\nTandem Repeats Finder Program written by:");
        fprintf(Fptxt, "\n\n                 Gary Benson");
        fprintf(Fptxt, "\n      Program in Bioinformatics");
        fprintf(Fptxt, "\n          Boston University");
        fprintf(Fptxt, "\n\nVersion %s", versionstring);

        fprintf(Fptxt, "\n\nSequence: %s\n\nParameters: %d %d %d %d %d %d %d\n",
            pseq->name, paramset.ps_match, paramset.ps_mismatch, paramset.ps_indel,
            paramset.ps_PM, paramset.ps_PI, paramset.ps_minscore, paramset.ps_maxperiod);
    }

    if (paramset.ps_ngs != 1)
        PrintProgress("Initializing data structures...");

    Distance = new_distancelist();
    clear_distancelist(Distance);
    Tag = newTags(MAXDISTANCE / TAGSEP + 1);
    Toptag = (int)ceil(MAXDISTANCE / TAGSEP);
    init_links();

    init_index();
    /* init_distanceseenlist(); */
    /* modified 5/23/05 G. Benson */
    init_distanceseenarray();

    if (paramset.ps_ngs != 1)
        PrintProgress("Computing TR Model Statistics...");

    init_and_fill_coin_toss_stats2000_with_4tuplesizes();

    /* over allocate statistics_distance array to prevent spill in alignments
     * with execive insertion counts Jan 07, 2003 */
    Statistics_Distance = calloc(4 * MAXDISTANCE, sizeof *Statistics_Distance);
    if (Statistics_Distance == NULL) {
        PrintError("Unable to allocate memory for Statistics_Distance array");
        exit(-3);
    }

    /* set the sequence pointer. more global vars! */
    Sequence = pseq->sequence - 1;  /* start one character before */
    Length = pseq->length;

    if (!paramset.ps_HTMLoff) {
        fprintf(Fptxt, "\n\nLength: %d", Length);
        fprintf(Fptxt, "\nACGTcount: A:%3.2f, C:%3.2f, G:%3.2f, T:%3.2f\n\n",
            (double)pseq->composition['A' - 'A'] / Length,
            (double)pseq->composition['C' - 'A'] / Length,
            (double)pseq->composition['G' - 'A'] / Length, (double)pseq->composition['T' - 'A'] / Length);

        if ((pseq->length - pseq->nucleotides) > 0) {
            fprintf(Fptxt, "Warning! %d characters in sequence are not A, C, G, or T\n\n",
                (pseq->length - pseq->nucleotides));
        }
    }

    Totalcharacters += Length;
    WDPcount = 0;

    /* G. Benson 1/28/2004 */
    /* following four memory allocations increased to avoid memory error when
     * consensus length exceeds MAXDISTANCE after returning from get_consensus(d) */

    Criteria_count = calloc(2 * (MAXDISTANCE + 1), sizeof *Criteria_count);
    if (Criteria_count == NULL) {
        PrintError("Unable to allocate Criteria_count");
        exit(-4);
    }

    Consensus_count = calloc(2 * (MAXDISTANCE + 1), sizeof *Consensus_count);
    if (Consensus_count == NULL) {
        PrintError("Unable to allocate memory for Consensus_count");
        exit(-5);
    }

    Cell_count = calloc(2 * (MAXDISTANCE + 1), sizeof *Cell_count);
    if (Cell_count == NULL) {
        PrintError("Unable to allocate memory for Cell_count");
        exit(-6);
    }

    Outputsize_count = calloc(2 * (MAXDISTANCE + 1), sizeof *Outputsize_count);
    if (Outputsize_count == NULL) {
        PrintError("Unable to allocate memory for Outputsize_count");
        exit(-7);
    }

    if (paramset.ps_multisequencefile) {
        sprintf(messagebuffer, "Scanning Sequence %d...", paramset.ps_sequenceordinal);
        if (paramset.ps_ngs != 1)
            PrintProgress(messagebuffer);
    }
    else if (paramset.ps_ngs != 1) {
        PrintProgress("Scanning...");
    }

    clear_distancelist(Distance);

    newtupbo();                 /* this is the main function of the algorithm */

    Cell_total = 0;
    for (i = MAXDISTANCE; i >= 1; i--) {
        Cell_total += Cell_count[i];
        if (i <= SMALLDISTANCE)
            Wasted_total += (Criteria_count[i] + Consensus_count[i]) * i * i * 2;
        else
            Wasted_total += (Criteria_count[i] + Consensus_count[i]) * i * (2 * d_range(i) + 1) * 2;
    }

    if (!paramset.ps_HTMLoff) {
        fprintf(Fptxt, "Done.\n");
        fprintf(Fptxt, "</PRE>");
        fprintf(Fptxt, "</BODY>");
        fprintf(Fptxt, "</HTML>");
        fclose(Fptxt);
    }

    /****************************************************************
     * The following memory deallocations where not originally
     * implemented. There may still exist some marginal memory leaks
     * but due messy code structure and extensive use of macros and
     * global variables this could be impossible to fix.
     *****************************************************************/

    if (paramset.ps_ngs != 1)
        PrintProgress("Freeing Memory...");

    free(S[0]);
    free(S);
    free(Statistics_Distance);
    free(Criteria_count);
    free(Consensus_count);
    free(Cell_count);
    free(Outputsize_count);
    free(AlignPair.textprime);
    free(AlignPair.textsecnd);
    free(AlignPair.indexprime);
    free(AlignPair.indexsecnd);
    free(Tag);
    free(Index);
    free_distanceseenarray();

    /* free distance list and all its entries */
    distanceentry_free();
    free(Distance);

    for (i = 1; i <= NTS; i++) {
        free(Tuplehash[i]);
        free(History[i]);
    }

    free(Sortmultiples);

    if (paramset.ps_ngs != 1)
        PrintProgress("Resolving output...");

    TRFClean(datstring, txtstring, MaxPeriod);

    /* Set the name of the outputfilename global to name given to
     * file in routines defined in trfclean.h */
    MakeFileName(paramset.ps_outputfilename, htmlstring, 1);

    free_bestperiodlist();

    if (paramset.ps_ngs != 1)
        PrintProgress("Done.");
}

void PrintError(char *errortext)
{
    /* the code here depends on the environment defined above */
    fprintf(stderr, "Error: %s\n", errortext);
}

void PrintProgress(char *progresstext)
{
    fprintf(stderr, "\n%s", progresstext);
}

void SetProgressBar(void)
{
    static int ready = 0;

    /* if percent is minus one then un-ready the progress indicator */
    if (paramset.ps_percent < 0) {
        if (!ready)
            return;

        ready = 0;
        /* hide the progress indicator */
        putchar('\n');
        return;
    }

    /* if not ready make ready */
    if (!ready) {
        ready = 1;

        /* position the progress progress indicator */
        putchar('\n');
    }

    /* set progress indicator to percent value */
    if (paramset.ps_percent > 100)
        paramset.ps_percent = 100;

    if (!(paramset.ps_percent % 5))
        putchar('.');
}


/*
 *   loads a sequence from an open file. If the routine comes to another
 *   sequence header while reading a sequence the sequence returns 1 to
 *   indicate that more sequences are present. Otherwise the routine returns
 *   0 to indicate EOF or -1 to indicate error.
 *   The member sequence must be NULL before the routine is called.
 *   The calling function must free the allocated memory after use.
 */
int LoadSequenceFromFileEugene(struct fastasequence *pseq, FILE *fp)
{
    int i, j, c;
    int next = -1;              // whether a next sequence was encountered
    char *ptemp;
    char to_upper;

    // read the FASTA '>' symbol
    c = getc(fp);
    if ((char)c != '>' || c == EOF)
        return -1;              // invalid format

    // read name and description text
    for (i = 0; i < MAXSEQNAMELEN - 1; i++) {
        c = getc(fp);
        if (c == 10 || c == 13) {
            break;
        }
        else if (c == EOF) {
            PrintError("FASTA input terminated too early");
            return -1;
        }
        else {
            pseq->name[i] = (char)c;
        }
    }
    pseq->name[i] = '\0';

    // if line was not read completely flush the rest
    if (i == MAXSEQNAMELEN - 1) {
        c = 0;
        while (c != 13 && c != 10 && c != EOF)
            c = getc(fp);
    }

    // read sequence into buffer
    pseq->sequence = NULL;
    for (i = 0; i < 26; i++)
        pseq->composition[i] = 0;
    to_upper = 'A' - 'a';

    i = 0;
    for (j = 0; c != EOF && (char)c != '>'; j += MEMBLOCK) {
        if ((ptemp = realloc(pseq->sequence, sizeof(char) * (i + MEMBLOCK + 1))) == NULL) {
            PrintError("Insufficient memory");
            return -1;
        }
        pseq->sequence = ptemp;
        for (i = j; i < j + MEMBLOCK;) {
            c = getc(fp);       // get a character from file
            if ((char)c >= 'A' && (char)c <= 'Z') { // in upper-case range of alpha characters
                pseq->sequence[i] = (char)c;
                pseq->composition[c - (int)'A']++;
                i++;
            }
            else if ((char)c >= 'a' && (char)c <= 'z') {    // in lower-case range of alpha characters
                pseq->sequence[i] = (char)c + to_upper; // make upper-case
                pseq->composition[c - (int)'a']++;
                i++;
            }
            else if ((char)c == '>') {  // break if another sequence found
                next = 1;
                ungetc('>', fp);
                break;
            }
            else if (c == EOF) {    // break if end of file
                next = 0;
                break;
            }
        }
    }
    pseq->length = i;           // set sequence length
    if (i > 0) {
        pseq->sequence[i] = '\0';   // terminate sequence text as a string
    }

    // compute member
    pseq->nucleotides =
        pseq->composition['A' - 'A'] + pseq->composition['C' - 'A'] +
        pseq->composition['G' - 'A'] + pseq->composition['T' - 'A'];

    return next;
}

