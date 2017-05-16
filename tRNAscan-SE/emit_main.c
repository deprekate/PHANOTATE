/* emit_main.c
 * main() for emitting sequences from a stored model
 * written as a debugging aid
 * 
 * 1.0: SRE, Tue Jun 15 09:32:43 1993
 * 2.0: SRE, Fri Sep 10 08:06:37 1993
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "structs.h"
#include "funcs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define OPTIONS "abls:L"

static char usage[]  = "\
Usage: covee [-options] <cmfile>\n\
where options are:\n\
     -a        : annotate all pairs, not just canonical ones\n\
     -b        : emit single most probable sequence\n\
     -l        : print as mountain landscape\n\
     -s <seed> : set seed for random()\n\
  EXPERIMENTAL OPTIONS:\n\
     -L        : calculate expected length distributions for states\n";


static char banner[] = "covee: emit sequences from a covariance model";

int
main(int argc, char **argv)
{
  char           *cmfile;       /* file to read model from */
  struct cm_s    *cm;           /* model                   */
  int             i;		/* counter for sequences   */
  char           *seq;          /* generated sequence      */
  char           *khseq;        /* generated structure     */
  struct align_s *ali;		/* generated "alignment"   */

  int          emitnum;		/* number of sequences to emit      */
  int          seed;		/* seed for random number generator */
  int          do_best;		/* TRUE if generating only best seq */
  int          do_landscape;	/* TRUE if printing as landscape    */
  int          watsoncrick;     /* TRUE to annotate only canonical pairs */
  int          do_lengths;	/* TRUE to do length distributions  */

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  emitnum      = 20;
  seed         = (int) time (0);
  do_best      = FALSE;
  do_landscape = FALSE;
  watsoncrick  = TRUE;
  do_lengths   = FALSE;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'a': watsoncrick  = FALSE; break;
    case 'b': do_best      = TRUE;  break;
    case 'l': do_landscape = TRUE;  break;
    case 's': seed = atoi(optarg);  break;
    case 'L': do_lengths   = TRUE;  break;

    default:
      Die("Error: unrecognized option %c\n", optc);
    }

  if (argc - optind != 1)
    Die("Wrong number of arguments.\n%s", usage);
  
  cmfile  = argv[optind];
  sre_srandom((unsigned) seed);

  /*********************************************** 
   * Read in the model
   ***********************************************/
  
  if (! ReadCM(cmfile, &cm))
    Die("Failed to read model from file %s\n", cmfile);

  /*********************************************** 
   * Generate sequences from model and print them.
   ***********************************************/

  puts(banner);
  printf("     version %s, %s\n\n", RELEASE, RELEASEDATE);

  if (do_lengths)
    {
      struct pstate_s *pcm;
      int              statenum;
      double         **lmx;
      int             *min, *max;
      int              y;
      
      MakePCM(cm, &pcm, &statenum);
      NormalizePCM(pcm, statenum);
      LengthDistribution(pcm, statenum, 200, &lmx);
      LengthBounds(lmx, statenum, 200, 1.0e-6, &min, &max);

      for (y = 0; y < statenum; y++)
	printf("%4d %4d %4d (%4d) %8.8f  %s\n", 
	       y, min[y], max[y], max[y]-min[y]+1,
	       (float) (max[y]-min[y]+1) / 200.0,
	       UstatetypeName(pcm[y].statetype));
      
      Free2DArray(lmx, statenum);
      free(min);
      free(max);
      free(pcm);
    }

  else if (do_best)
    {
      if (! EmitBestSequence(cm, watsoncrick, &ali, &khseq, &seq)) Die("EmitBestSequence() failed");
      
      if (do_landscape)
	{
	  if (! PrintAliLandscape(stdout, cm, ali))
	    Warn("PrintAliLandscape failed\n");
	}
      else
	{
	  printf("%s\n", seq);
	  printf("%s\n", khseq);
	  puts("");
	}
      
      free(ali);
      free(khseq);
      free(seq);
    }

  else
    {
      for (i = 0; i < emitnum; i++)
	{
	  if (! EmitSequence(cm, watsoncrick, &ali, &khseq, &seq))
	    Die("failed to generate a sequence from the model.");

	  if (do_landscape)
	    PrintAliLandscape(stdout, cm, ali);
	  else
	    {
	      printf("seq %2d: %s\n", i, seq);
	      printf("        %s\n", khseq);
	      puts("");
	    }
      
	  free(ali);
	  free(khseq);
	  free(seq);
	}
    }

  /*********************************************** 
   * Cleanup and exit
   ***********************************************/
  
  FreeCM(cm);
  return 0;
}
  
