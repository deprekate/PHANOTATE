/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* main for shuffle
 *
 * shuffle - generate shuffled sequences
 * Mon Feb 26 16:56:08 1996
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "squid.h"

struct opt_s OPTIONS[] = {
  { "-h",     TRUE, ARG_NONE },	    /* help                                  */
  { "-n",     TRUE, ARG_INT  },     /* number of shuffled seqs per input seq */
  { "--seed", FALSE, ARG_INT },	    /* set the random number seed            */
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


char usage[]  = "Usage: shuffle [-options] <seqfile>\n\
  Generate shuffled copies of input sequences.\n\
  Available options:\n\
  -h         : help; print version and usage info\n\
  -n <num>   : make <num> shuffles per input seq (default 1)\n\
  --seed <s> : set random number seed to <s>\n\
";


int
main(int argc, char **argv)
{
  char  *seqfile;               /* name of sequence file */
  SQFILE *dbfp;			/* open sequence file */
  int    fmt;			/* format of seqfile  */
  char  *seq;			/* sequence */
  SQINFO sqinfo;                /* additional sequence info */
  char  *shuff;                 /* shuffled sequence */
  int    num;			/* number to generate */
  int    seed;			/* random number generator seed */
  int    i;

  char  *optname;               /* option name */
  char  *optarg;                /* option argument (or NULL) */
  int    optind;                /* index of next argv[] */  


  /***********************************************
   * Parse command line
   ***********************************************/

  num  = 1;
  seed = (int) time ((time_t *) NULL);

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, &optind, &optname, &optarg))
    {
      if      (strcmp(optname, "-n")     == 0) { num  = atoi(optarg); }
      else if (strcmp(optname, "--seed") == 0) { seed = atoi(optarg); }
      else if (strcmp(optname, "-h")     == 0) {
	printf("shuffle %s, %s\n%s\n", squid_version, squid_date, usage);
	exit(EXIT_SUCCESS);
      }
    }

  if (argc - optind != 1) Die("%s\n", usage); 
  seqfile = argv[optind];

  sre_srandom(seed);
		 
  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Failed to determine format of file %s", seqfile);
  if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);
  
  while (ReadSeq(dbfp, fmt, &seq, &sqinfo))
    {
      shuff = (char *) MallocOrDie ((sqinfo.len + 1) * sizeof(char));

      for (i = 0; i < num; i++)
	{
	  StrShuffle(shuff, seq);
	  WriteSeq(stdout, kPearson, shuff, &sqinfo);
	}

      free(shuff);
      FreeSequence(seq, &sqinfo);
    }

  SeqfileClose(dbfp);
  return 0;
}

