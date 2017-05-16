/* score_main.c
 * Fri Feb 18 10:31:48 1994
 * 
 * main() for scoring test sequences with a model.
 *   Also, can print out alignments of model to sequence so
 *   that the pairwise assignments can be seen.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "structs.h"
#include "funcs.h"
#include "squid.h" 
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define OPTIONS "ag:ms"

static char usage[]  = "\
Usage: coves [-options] <CM file> <seqfile>\n\
where options are:\n\
    -a          : show all pairs, not just Watson-Crick\n\
    -g <gcfrac> : set expected background GC composition (default 0.5)\n\
    -m          : mountain representation of structural alignment\n\
    -s          : secondary structure string representation of \n\
                  structural alignment\n";

static char banner[] = "\
coves - scoring and structure prediction of RNA sequences\n\
        using a covariance model";

int
main(int argc, char **argv)
{ 
  char  *seq;                   /* a sequence to score            */
  SQINFO sqinfo;                /* info about seq                 */
  char  *seqfile;               /* sequence file                  */
  int    fmt;			/* format of sequence file        */
  SQFILE *dbfp;                 /* open sequence file for reading */
  char  *cmfile;                /* file containing covariance model */
  struct cm_s *cm;              /* model                          */
  double score;                 /* score of alignment             */
  struct trace_s *tr;           /* traceback of alignment         */
  struct align_s *ali;		/* alignment of seq to model      */
  char  *aseq;                  /* "aligned" sequence string      */
  char  *khstruct;              /* secondary structure string     */
  char   buffer[61];		/* output buffer for structures   */
  int    len;			/* length of aseq, khstruct       */
  int    apos;			/* position in aseq, khstruct     */
  double rfreq[ALPHASIZE];      /* expected background symbol frequencies */
  struct istate_s *icm;		/* integer log odds model         */
  int    statenum;		/* # of states in icm             */

  int     do_khstructure;	/* TRUE if we print a structure string  */
  int     do_mountain;		/* TRUE if we show a mountain structure */
  int     watsoncrick;          /* TRUE if only canonical pairs are indicated */
  double  gcfrac;		/* expected background GC fraction */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, curr_size;
#endif

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  do_khstructure = FALSE;
  do_mountain    = FALSE;
  watsoncrick    = TRUE;
  gcfrac         = 0.5;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'a': watsoncrick    = FALSE;                 break;
    case 'g': gcfrac         = (double) atof(optarg); break;
    case 'm': do_mountain    = TRUE;                  break;
    case 's': do_khstructure = TRUE;                  break;

    case 'h': 
      printf("%s\n  version %s (%s)\n%s\n", banner, RELEASE, RELEASEDATE, usage);
      exit(0);

    default: Die("unrecognized option %c\n", optc);
    }

  if (argc - optind != 2)
    Die("%s\n", usage);
  
  cmfile  = argv[argc-2];
  seqfile = argv[argc-1]; 

  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Failed to determine format of sequence database %s", seqfile);

  if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);

  if (! ReadCM(cmfile, &cm))
    Die("Failed to read model from file %s", cmfile);

  rfreq[1] = rfreq[2] = gcfrac / 2.;
  rfreq[0] = rfreq[3] = (1. - gcfrac) / 2.;
              
  if (! RearrangeCM(cm, rfreq, &icm, &statenum))
    Die("failed to make integer log odds model");

  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);
  printf("     version %s, %s\n\n", RELEASE, RELEASEDATE);

  printf("---------------------------------------------------\n");
  printf("Database to search/score:      %s\n", seqfile);
  printf("Model:                         %s\n", cmfile);
  printf("GC%% of background model:       %.0f%%\n", (gcfrac*100.));
  printf("---------------------------------------------------\n");
  puts("");

  /*********************************************** 
   * Score each sequence
   ***********************************************/

#ifdef MEMDEBUG
  orig_size = malloc_size(&histid1);
#endif

  while (ReadSeq(dbfp, fmt, &seq, &sqinfo))
    {
      char *prepseq;
      prepseq = Strdup(seq);
      PrepareSequence(prepseq);

      if (! ViterbiAlign(icm, statenum, prepseq, &score, &tr))
	Die("ViterbiAlign() failed on sequence %s", sqinfo.name);
      free(prepseq);
      
      printf("%6.2f bits : %s\n", score, sqinfo.name);
      
      if (do_khstructure || do_mountain)
	{
	  if (! Trace2ali(seq, tr, watsoncrick, &ali)) Die("Trace2ali failed");

	  if (do_khstructure)
	    {
	      if (! Align2kh(ali, &aseq, &khstruct)) Die("Align2kh failed\n");
	      
	      /* Print out the sequence and structure
	       */
	      len = strlen(aseq);
	      buffer[60] = '\0';
	      for (apos = 1; apos <= len; apos += 60)
		{
		  strncpy(buffer, aseq + apos - 1, 60);
		  printf("    %10s %s\n", sqinfo.name, buffer);
		  strncpy(buffer, khstruct + apos - 1, 60);
		  printf("    %10s %s\n", sqinfo.name, buffer);
		  puts("");
		}
	      free(aseq);
	      free(khstruct);
	    }
	  
	  if (do_mountain)
	    {
	      PrintAliLandscape(stdout, cm, ali);
	      puts("");
	    }

	  Free_align(ali);
	}

      FreeSequence(seq, &sqinfo);

#ifdef MEMDEBUG
      curr_size = malloc_size(&histid2);
      if (curr_size != orig_size)
	{
	  Warn("malloc-debug: current size %ul, starting size %ul\n", curr_size, orig_size);
	  malloc_list(2,histid1, histid2);
	}
#endif

    }

  SeqfileClose(dbfp);
  FreeCM(cm);
  free(icm);
  return 0;
}


