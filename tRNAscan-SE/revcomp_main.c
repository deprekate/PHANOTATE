/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* main for revcomp
 *
 * revcomp - generate reverse complement of sequences
 * SRE, Thu Aug  5 17:36:57 1993
 */

#include <stdio.h>
#include <string.h>
#include "squid.h"

#define OPTIONS "h"

char usage[]  = "Usage: revcomp [-options] <seqfile>\n\
  Reverse complement a nucleic acid sequence.\n\
  Available options:\n\
  -h    : help; print version and usage info\n";

int
main(int argc, char **argv)
{
  char  *seqfile;               /* name of sequence file */
  SQFILE *dbfp;			/* open sequence file */
  int    fmt;			/* format of seqfile  */
  char  *seq;			/* sequence */
  SQINFO sqinfo;                /* additional sequence info */
  char  *rev;			/* reverse complement */
  int    swap;

  int    optchar;		/* option character, command line */
  extern int   optind;

  /***********************************************
   * Parse command line
   ***********************************************/

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {
    case 'h': 
      printf("revcomp %s, %s\n%s\n", squid_version, squid_date, usage);
      exit(EXIT_SUCCESS);
    default:
      Die("%s\n", usage);
    }

  if (argc - optind != 1) Die("%s\n", usage); 
  seqfile = argv[optind];
		 
  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Failed to determine format of file %s", seqfile);
  if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);
  
  while (ReadSeq(dbfp, fmt, &seq, &sqinfo))
    {
      if ((rev = (char *) malloc ((sqinfo.len + 1) * sizeof(char))) == NULL)
	Die("malloc failed");

      revcomp(rev, seq);
      if (sqinfo.flags & (SQINFO_START | SQINFO_STOP))
	{
	  swap         = sqinfo.start;
	  sqinfo.start = sqinfo.stop;
	  sqinfo.stop  = swap;
	}
	/* secondary structure of reverse strand is nonsense
	 */
      if (sqinfo.flags & SQINFO_SS)
	{
	  sqinfo.flags = sqinfo.flags & ~SQINFO_SS;
	  free(sqinfo.ss);
	}

      WriteSeq(stdout, kPearson, rev, &sqinfo);

      free(rev);
      FreeSequence(seq, &sqinfo);
    }

  SeqfileClose(dbfp);
  return 0;
}
