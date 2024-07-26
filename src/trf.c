
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

/**************************************************************
 * trf.c :   Command line interface to the Tandem Repeats Finder
 *
 ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>              // ERANGE
#include <limits.h>             // LONG_MIN, LONG_MAX
#include "trfrun.h"
#include "tr30dat.h"

const char *usage = "\n\nPlease use: %s File Match Mismatch Delta PM PI Minscore MaxPeriod [options]\n"
    "\nWhere: (all weights, penalties, and scores are positive)"
    "\n  File = sequences input file"
    "\n  Match  = matching weight"
    "\n  Mismatch  = mismatching penalty"
    "\n  Delta = indel penalty"
    "\n  PM = match probability (whole number)"
    "\n  PI = indel probability (whole number)"
    "\n  Minscore = minimum alignment score to report"
    "\n  MaxPeriod = maximum period size to report. Must be between 1 and 2000, inclusive"
    "\n  [options] = one or more of the following:"
    "\n        -m        masked sequence file"
    "\n        -f        flanking sequence"
    "\n        -d        data file"
    "\n        -h        suppress html output"
    "\n        -r        no redundancy elimination"
    "\n        -l <n>    maximum TR length expected (in millions) (eg, -l 3 or -l=3 for 3 million)"
    "\n                  Human genome HG38 would need -l 6"
    "\n        -ngs      more compact .dat output on multisequence files, returns 0 on success."
    "\n                  Output is printed to the screen, not a file. You may pipe input in with"
    "\n                  this option using - for file name. Short 50 flanks are appended to .dat"
    "\n                  output."
    "\n"
    "\nSee more information on the TRF Unix Help web page: https://tandem.bu.edu/trf/trf.unix.help.html"
    "\n"
    "\nNote the sequence file should be in FASTA format:"
    "\n"
    "\n>Name of sequence"
    "\naggaaacctgccatggcctcctggtgagctgtcctcatccactgctcgctgcctctccag"
    "\natactctgacccatggatcccctgggtgcagccaagccacaatggccatggcgccgctgt"
    "\nactcccacccgccccaccctcctgatcctgctatggacatggcctttccacatccctgtg"
    "\n"
    "\n"
    "\nTandem Repeats Finder"
    "\nCopyright (C) 1999-2020 Gary Benson"
    "\n"
    "\nThis program is free software: you can redistribute it and/or modify"
    "\nit under the terms of the GNU Affero General Public License as"
    "\npublished by the Free Software Foundation, either version 3 of"
    "\nthe License, or (at your option) any later version." "\n";

char *GetNamePartAddress(char *name);
void PrintBanner(void);
static int ParseInt(const char *str, int *dest);
static int ParseUInt(const char *str, unsigned int *dest);

int main(int ac, char **av)
{
    /* Handle a lone -v argument ourselves */
    if ((ac == 2) && ((strcmp(av[1], "-v") == 0) || (strcmp(av[1], "-V") == 0))) {
        PrintBanner();
        exit(0);
    }

    /* Expects exactly 8 non-option arguments */
    if (ac < 9) {
        fprintf(stderr, usage, av[0]);
        exit(1);
    }

    /* set option defaults */
    g_paramset.ps_use_stdin = 0;
    g_paramset.ps_datafile = 0;
    g_paramset.ps_maskedfile = 0;
    g_paramset.ps_flankingsequence = 0;
    g_paramset.ps_flankinglength = 500;  /* Currently not user-configurable */
    g_paramset.ps_HTMLoff = 0;
    g_paramset.ps_redundoff = 0;
    g_paramset.ps_maxwraplength = 2000000;
    g_paramset.ps_ngs = 0;           /* this is for unix systems only */

    /* Parse command line options */
    /* Assume that since the first checks were passed, options start at argument 8
     * getopt ignores the first array element, so start at element 8 */
    char **opt_arr = &av[8];
    int remaining_opts = ac - 8;

    while (1) {
        static struct option long_options[] = {
            { "help", no_argument, 0, 'u' },    /* -u, -U */
            { "version", no_argument, 0, 'v' }, /* -v, -V */
            { "dat", no_argument, &g_paramset.ps_datafile, 1 },  /* -d, -D */
            { "mask", no_argument, &g_paramset.ps_maskedfile, 1 },   /* -m, -M */
            { "flank", no_argument, &g_paramset.ps_flankingsequence, 1 },    /* -f, -F */
            { "html-off", no_argument, &g_paramset.ps_HTMLoff, 1 },  /* -h, -H */
            { "redund-off", no_argument, &g_paramset.ps_redundoff, 1 },  /* -r, -R */
            { "ngs", no_argument, &g_paramset.ps_ngs, 1 },   /* -ngs */
            { "Ngs", no_argument, &g_paramset.ps_ngs, 1 },   /* -Ngs */
            { "NGS", no_argument, &g_paramset.ps_ngs, 1 },   /* -NGS */
            { "maxlength", required_argument, 0, 'l' }, /* -l, -L */
            { 0, 0, 0, 0 }
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        /* Accept upper and lower-case variants of options */
        int c = getopt_long_only(remaining_opts, opt_arr, "uvdmfhrUVDMFHRl:L:",
            long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 0:
                if (long_options[option_index].flag == 0)
                    break;

                /* Here you can use a function like strcasecmp to check if a
                 * particular flag long option was given, and then do something
                 * if you like (eg, print a message when testing). */

                // Example:
                /*if (strcasecmp(long_options[option_index].name, "ngs") == 0){
                 * printf("NGS option was given.\n");
                 * }
                 */

                break;

            case 'u':
            case 'U':
                printf(usage, av[0]);
                exit(0);
                break;

            case 'v':
            case 'V':
                PrintBanner();
                exit(0);
                break;

            case 'd':
            case 'D':
                g_paramset.ps_datafile = 1;
                break;

            case 'm':
            case 'M':
                g_paramset.ps_maskedfile = 1;
                break;

            case 'f':
            case 'F':
                g_paramset.ps_flankingsequence = 1;
                break;

            case 'h':
            case 'H':
                g_paramset.ps_HTMLoff = 1;
                break;

            case 'r':
            case 'R':
                g_paramset.ps_redundoff = 1;
                break;

            case 'l':
            case 'L':
                if ((atol(optarg) < 1)) {
                    fprintf(stderr, "Error: max TR length must be at least 1 million\n");
                    PrintBanner();
                    exit(2);
                }

                if (ParseUInt(av[8], &g_paramset.ps_maxwraplength) == 0) {
                    fprintf(stderr, "Error while parsing max TR length (option '-L') value\n");
                    PrintBanner();
                    exit(1);
                }
                g_paramset.ps_maxwraplength *= 1e6;

                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                break;
        }
    }

    /* get input parameters */
    strcpy(g_paramset.ps_inputfilename, av[1]);
    strcpy(g_paramset.ps_outputprefix, GetNamePartAddress(av[1]));

    /* Validate these parameters */
    if (ParseUInt(av[2], &g_paramset.ps_match) == 0) {
        g_paramset.ps_endstatus = "Error parsing match parameter." " Value must be a positive integer.";
    }
    else if (ParseUInt(av[3], &g_paramset.ps_mismatch) == 0) {
        g_paramset.ps_endstatus = "Error parsing mismatch parameter." " Value must be a positive integer.";
    }
    else if (ParseUInt(av[4], &g_paramset.ps_indel) == 0) {
        g_paramset.ps_endstatus = "Error parsing indel parameter." " Value must be a positive integer.";
    }
    else if (ParseUInt(av[5], &g_paramset.ps_PM) == 0) {
        g_paramset.ps_endstatus = "Error parsing PM parameter." " Value must be a positive integer.";
    }
    else if (ParseUInt(av[6], &g_paramset.ps_PI) == 0) {
        g_paramset.ps_endstatus = "Error parsing PI parameter." " Value must be a positive integer.";
    }
    else if (ParseUInt(av[7], &g_paramset.ps_minscore) == 0) {
        g_paramset.ps_endstatus = "Error parsing Minscore parameter." " Value must be a positive integer.";
    }
    else if ((ParseUInt(av[8], &g_paramset.ps_maxperiod) == 0) || (g_paramset.ps_maxperiod > 2000) || (g_paramset.ps_maxperiod == 0)) {
        g_paramset.ps_endstatus = "Error parsing MaxPeriod parameter." " Value must be between 1 and 2000, inclusive.";
    }

    /* Error if any validation failed */
    if (g_paramset.ps_endstatus) {
        printf("%s\n\n", g_paramset.ps_endstatus);
        printf("Please run with -h for help, or visit https://tandem.bu.edu/trf/trf.unix.help.html\n");
        exit(1);
    }

    // g_paramset.datafile must be set if HTMLoff is set
    g_paramset.ps_datafile |= g_paramset.ps_HTMLoff;

    if (g_paramset.ps_ngs == 1) {
        g_paramset.ps_datafile = 1;
    }
    else {
        PrintBanner();
    }

#if (defined(UNIXGUI)+defined(UNIXCONSOLE))>=1
    if (0 == strcmp("-", av[1])) {
        if (g_paramset.ps_ngs) {
            g_paramset.ps_use_stdin = 1;
            g_paramset.ps_HTMLoff = 1;
        }
        else {
            fprintf(stderr, "\n\nPlease use -ngs flag if piping input into TRF");
            fprintf(stderr, "\n");
            exit(-1);
        }
    }
#endif

    /* call the fuction on trfrun.h that controls execution */
    TRFControlRoutine();

    /* Check the status by looking at the output count and the endstatus
     * members of the g_paramset struct */
    if (g_paramset.ps_outputcount == 0 && !g_paramset.ps_endstatus) {
        /* If no TRs were found, print informative message to STDERR */
        printf("No TRs found. Exiting...\n");
        return 0;
    }
    else if (g_paramset.ps_endstatus) {
        printf("Error processing input: %s\n", g_paramset.ps_endstatus);
        return 1;
    }
    else {
        return 0;
    }
}

char *GetNamePartAddress(char *name)
{
    int i;
    char *pname;

    char dirsymbol = '/';

    i = strlen(name) - 1;
    pname = &name[i];
    while (i > 0 && *pname != dirsymbol) {
        pname--;
        i--;
    }
    if (*pname == dirsymbol)
        pname++;
    return pname;
}

// Use these two function for parameter validation
// TODO Move declarations to other header file
static int ParseInt(const char *str, int *dest)
{
    errno = 0;
    char *temp;
    long val = strtol(str, &temp, 0);

    // Print message and return false (0) if unable
    // to convert to a valid int
    if (temp == str || *temp != '\0' || ((val == LONG_MIN || val == LONG_MAX) && errno == ERANGE)) {
        // fprintf(stderr, "Error parsing parameter '%s' as integer value.\n",
        //         str);
        return 0;
    }

    // Set dest pointer and return true (1) if successful
    *dest = val;
    return 1;
}

static int ParseUInt(const char *str, unsigned int *dest)
{
    int temp;

    int success = ParseInt(str, &temp);

    // Return truth value of: conversion successful
    // AND the value is positive. Also set *dest if
    // those conditions are satisfied.
    if (success && (temp >= 0)) {
        *dest = temp;
        return 1;
    }
    else {
        return 0;
    }
}

void PrintBanner(void)
{
    fprintf(stderr, "\nTandem Repeats Finder, Version %s", versionstring);
    fprintf(stderr, "\nCopyright (C) Dr. Gary Benson 1999-2012. All rights reserved.\n");

    return;
}
