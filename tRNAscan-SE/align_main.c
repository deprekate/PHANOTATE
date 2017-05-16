/* align_main.c
 * SRE, Wed Jun 30 09:56:15 1993
 * 2.0 Thu Sep 30 14:23:57 1993
 *
 * main() for covea
 * Multiple sequence alignment to a covariance HMM model.
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

#define OPTIONS "aho:s:S"

static char usage[]  = "\
Usage: covea [-options] <cm file> <seqfile>\n\
where supported options are:\n\
     -a             : annotate all base pairs, not just canonical ones\n\
     -h             : print short help and version info\n\
     -o <outfile>   : write alignment to <outfile> in SELEX format\n\
     -s <scorefile> : save individual alignment scores to <scorefile>\n\
\n\
Experimental options:\n\
     -S             : use small-memory variant of alignment algorithm\n";

static char banner[] = "covea: multiple sequence alignment to a covariance model";

int
main(int argc, char **argv)
{ 
  char **rseqs;			/* raw sequences to align       */
  char **aseqs;                 /* multiple sequence alignment  */
  SQINFO *sqinfo;		/* array of info structures     */
  int    nseq;			/* number of seqs               */ 
  char  *seqfile;               /* sequence file                */
  int    format;		/* format of sequence file      */
  char  *cmfile;                /* cvhmm save file to read      */
  struct cm_s *cm;              /* model                        */
  struct trace_s **tr;          /* array of tracebacks for seqs */
  int    idx;			/* counter for sequences        */
  double score;                 /* score of indiv. alignment    */
  double tot_score;		/* sum of scores                */
  AINFO  ainfo;                 /* optional alignment info (sec structure) */
  struct istate_s *icm;         /* model, integer log odds form */
  int    statenum;		/* # of states in icm           */
  double rfreq[ALPHASIZE];	/* expected background symbol frequencies */

  char  *outfile;		/* file to write alignment to       */
  char  *scorefile;		/* file to save scores to           */
  FILE  *ofp;			/* opened outfile                   */
  FILE  *sfp;			/* opened scorefile                 */
  int    do_smallmemory;	/* use small-memory viterbi variant */
  int    watsoncrick;           /* annotate only canonical pairs    */

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/

  outfile        = NULL;
  scorefile      = NULL;
  do_smallmemory = FALSE;
  watsoncrick    = TRUE;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'a': watsoncrick = FALSE; break;
    case 'o': outfile   = optarg;  break;
    case 's': scorefile = optarg;  break;

    case 'S': do_smallmemory = TRUE; break;

    case 'h': 
      printf("%s\n  version %s (%s)\n%s\n", banner, RELEASE, RELEASEDATE, usage);
      exit(0);
    default:
      Die("unrecognized option %c\n", optc);
    }

  if (argc - optind != 2)
    Die("%s\n", usage);
  
  cmfile  = argv[argc-2];
  seqfile = argv[argc-1];

#ifdef MEMDEBUG
  orig_size = malloc_size(&histid1);
#endif


  /*********************************************** 
   * Get sequence data and model; open output ptrs
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, NULL))
    Die("Failed to determine format of file %s\n", seqfile);

  if (! ReadMultipleRseqs(seqfile, format, &rseqs, &sqinfo, &nseq))
    Die("Failed to read sequences from file %s", seqfile);

  if (! ReadCM(cmfile, &cm))
    Die("Failed to read model from file %s", cmfile);

  rfreq[0] = rfreq[1] = rfreq[2] = rfreq[3] = 0.25;
  if (! RearrangeCM(cm, rfreq, &icm, &statenum))
    Die("Failed to convert CM to integer log odds");

  if (outfile != NULL)
    if ((ofp = fopen(outfile, "w")) == NULL)
      Die("Open failed for alignment output file %s", outfile);

  if (scorefile != NULL)
    if ((sfp = fopen(scorefile, "w")) == NULL)
      Die("Open failed for score output file %s", scorefile);

  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);
  printf("     release %s, %s\n\n", RELEASE, RELEASEDATE);
  printf("---------------------------------------------------\n");
  printf("Sequence data:          %s (%d sequences)\n", seqfile, nseq);
  printf("Covariance model:       %s (%d nodes)\n", cmfile, cm->nodes);
  if (outfile != NULL) 
    printf("Alignment saved to:     %s\n", outfile);
  if (scorefile != NULL)
    printf("Indiv. scores saved to: %s\n", scorefile);
  printf("---------------------------------------------------\n");
  puts("");

  /*********************************************** 
   * Do the alignment
   ***********************************************/

  if ((tr = (struct trace_s **) malloc (nseq * sizeof(struct trace_s *))) == NULL)
    Die("Memory failure, line %d of %s", __LINE__, __FILE__);

  tot_score = 0.0;
  for (idx = 0; idx < nseq; idx++)
    {
      char *prepseq;
      prepseq = Strdup(rseqs[idx]);
      PrepareSequence(prepseq);

      if (do_smallmemory)
	{
	  if (! SmallViterbiAlign(icm, statenum, prepseq, &score, &tr[idx]))
	    Die("SmallViterbiAlign() failed on sequence %d", idx);
	}
      else if (! ViterbiAlign(icm, statenum, prepseq, &score, &tr[idx]))
	Die("ViterbiAlign() failed on sequence %d", idx);

      tot_score += score;
      if (scorefile != NULL)
	fprintf(sfp, "%-8.3f   : %s\n", score, sqinfo[idx].name);

      free(prepseq);
    }

  if (do_smallmemory)
    {
      printf("aborting... no traceback/alignment code yet for small memory variant\n");
      Free2DArray(rseqs,    nseq);
      FreeCM(cm);
      exit(0);
    }

  if (! Traces2Alignment(rseqs, sqinfo, tr, nseq, cm, watsoncrick, &aseqs, &ainfo))
    Die("Traces2Alignment() failed");

  /*********************************************** 
   * Print the alignment
   ***********************************************/

  if (outfile != NULL)
    {
      if (! WriteSELEX(ofp, aseqs, nseq, &ainfo, 60))
	Die("Write failed: can't save alignment to %s", outfile);
      fclose(ofp);
      printf("Alignment written to %s\n", outfile);
    }
  else
    {
      if (! WriteSELEX(stdout, aseqs, nseq, &ainfo, 60))
	Die("Write failed: can't print alignment");
    }

  if (scorefile != NULL) fclose(sfp);

  printf("Overall alignment score: %.2f\n", tot_score / (double) nseq);

  /*********************************************** 
   * Garbage collect and exit
   ***********************************************/

  for (idx = 0; idx < nseq; idx++)  
    {
      FreeTrace(tr[idx], NULL);
      FreeSequence(rseqs[idx],&(sqinfo[idx]));
    }
  free(tr);
  free(sqinfo);
  FreeAlignment(aseqs, nseq, &ainfo);
  FreeCM(cm);
  free(icm);

#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "No memory leaks, sir.\n");
  
#endif
  

  return 0;
}
