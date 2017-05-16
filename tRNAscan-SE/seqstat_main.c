/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* seqstat_main.c
 * Wed Aug 10 15:47:14 1994
 * 
 * Look at a sequence file, determine some simple statistics.
 */

#include <stdio.h>
#include <string.h>
#include "squid.h"

#define OPTIONS "ah"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char usage[]  = "Usage: seqstat [-options] <seqfile>\n\
  Verify a sequence file; print some simple statistics and info.\n\
  Available options:\n\
  -a    : report per-sequence info, not just a summary\n\
  -h    : help; display usage and version\n";

int
main(int argc, char **argv)
{
  char     *seqfile;            /* name of sequence file     */
  SQFILE   *dbfp;		/* open sequence file        */
  int       fmt;		/* format of seqfile         */
  char     *seq;		/* sequence                  */
  SQINFO    sqinfo;             /* extra info about sequence */
  int       nseqs;
  int       small;		/* smallest length */
  int       large;		/* largest length  */
  int       total;              /* total length    */

  int    optchar;		/* option character, command line */
  extern int   optind;
  int    allreport;

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
#endif

  /***********************************************
   * Parse command line
   ***********************************************/

  allreport = FALSE;		/* default: summary only */
  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {
    case 'a': allreport = TRUE; break;
    case 'h': 
      printf("seqstat %s, %s\n%s\n", squid_version, squid_date, usage);
      exit(EXIT_SUCCESS);
    default:
      Die("%s\n", usage);
    }

  if (argc - optind != 1) Die("%s\n", usage);

  seqfile = argv[argc-1];

#ifdef MEMDEBUG
  orig_size = malloc_size(&histid1);
#endif

  /***********************************************
   * Read the file.
   ***********************************************/

  printf("seqstat %s, %s\n\n", squid_version, squid_date);

  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Can't determine format of file %s\n", seqfile);

  if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);
  
  if (allreport) {
    printf("  %-15s %-5s %s\n", "  NAME", "LEN", "DESCRIPTION");
    printf("  --------------- ----- -----------\n");
  }

  small = 9999999;
  large = 0;
  nseqs = 0;
  total = 0;
  while (ReadSeq(dbfp, fmt, &seq, &sqinfo))
    {
      if (allreport) 
	printf("* %-15s %5d %-50.50s\n", sqinfo.name, sqinfo.len, 
	       sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");

      if (sqinfo.len < small) small = sqinfo.len;
      if (sqinfo.len > large) large = sqinfo.len;
      total += sqinfo.len;
      nseqs++;
      FreeSequence(seq, &sqinfo);
    }
  if (allreport) puts("");

  printf("Format:              %s\n", SeqFormatString(fmt));
  printf("Number of sequences: %d\n", nseqs);
  printf("Total # residues:    %d\n", total);
  printf("Smallest:            %d\n", small);
  printf("Largest:             %d\n", large);
  printf("Average length:      %.1f\n", (float) total / (float) nseqs);

  SeqfileClose(dbfp);

#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "[No memory leaks]\n");
#endif
    

  return 0;
}
