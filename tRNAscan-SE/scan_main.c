/* scan_main.c
 * Mon Jan 31 11:04:57 1994
 * 
 * main() for database scanning with a covariance model.
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

#define OPTIONS "cg:ho:t:w:D:E:F"

static char usage[]  = "\
Usage: covels [-options] <CM file> <seqfile>\n\
where options are:\n\
    -c          : do complementary strand too\n\
    -g <fracGC> : set background expected GC content (0.5 default)\n\
    -h          : print short help and version info\n\
    -o <file>   : save hits in <file>\n\
    -t <thresh> : set score reporting threshold\n\
    -w <window> : set scanning window size\n\
 CRASH PROTECTION OPTIONS:\n\
    -D <filename> : save name of last sequence processed\n\
 EXPERIMENTAL OPTIONS:\n\
    -E          : set epsilon for fast search\n\
    -F          : fast heuristic search\n";


static char banner[] = "covels - scan sequences for matches to an RNA covariance model";

static char *ext_seqname;
static int   ext_seqlen;
static int   in_complement;
static FILE *ext_ofp;

static int print_hit(int i, int j, double sc);

int
main(int argc, char **argv)
{ 
  char  *seq;                   /* a sequence to score            */
  SQINFO sqinfo;                /* info about sequence            */
  char  *rev;			/* rev complement of seq          */
  char  *seqfile;               /* sequence file                  */
  int    fmt;			/* format of sequence file        */
  SQFILE *dbfp;                 /* open sequence file for reading */
  char  *cmfile;                /* file containing covariance model */
  struct cm_s *cm;              /* model                          */
  struct istate_s *icm;         /* integer log-odds search model  */
  struct pstate_s *pcm;		/* rearranged probability model   */
  int    statenum;              /* # of states in icm, pcm        */
  int   *minb, *maxb;		/* bounds for fast heuristic search */
  double rfreq[ALPHASIZE];      /* random model                   */

  int        fast;		/* TRUE if trying fast heuristic search       */
  double     epsilon;		/* how much probability we're willing to lose */
  int        do_complement;	/* TRUE if searching complementary strand too */
  double     thresh;		/* threshold score for reporting a match      */
  char      *outfile;		/* save file for scores                       */
  int        window;		/* size of search window, symbols             */
  char      *donefile;		/* crash protection: save name of last seq    */
  double     gcfrac;		/* background gc fraction                     */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, curr_size;
#endif

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  do_complement = FALSE;
  thresh        = 0.0;
  window        = 100;
  outfile       = NULL;
  fast          = FALSE;
  epsilon       = 1e-9;
  donefile      = NULL;
  gcfrac        = 0.5;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'c': do_complement = TRUE;                  break;
    case 'g': gcfrac        = (double) atof(optarg); break;
    case 'o': outfile       = optarg;                break;
    case 't': thresh        = (double) atof(optarg); break;
    case 'w': window        = atoi(optarg);          break;

    case 'D': donefile      = optarg;                break;
    case 'E': epsilon       = atof(optarg);          break;
    case 'F': fast          = TRUE;                  break;

    case 'h': 
      printf("%s\n  version %s (%s)\n%s\n", banner, RELEASE, RELEASEDATE, usage);
      exit(0);

    default: Die("unrecognized option %c\n", optc);
    }

  if (argc - optind != 2)
    Die("%s\n", usage);
  
  cmfile  = argv[argc-2];
  seqfile = argv[argc-1]; 

  /* The random model probabilities
   */
  rfreq[1] = rfreq[2] = gcfrac / 2.0;
  rfreq[0] = rfreq[3] = (1.0 - gcfrac) / 2.0;

  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Failed to determine format of sequence database %s", seqfile);

  if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);

  if (! ReadCM(cmfile, &cm))
    Die("Failed to read model from file %s", cmfile);

  if (! RearrangeCM(cm, rfreq, &icm, &statenum))
    Die("Failed to create search CM");

  /* Fast version.
   * Use lower/upper bounds on possible subsequence lengths
   * for each state. These are obtained probabilistically.
   */
  if (fast)
    {
      double **lmx;

      if (! MakePCM(cm, &pcm, &statenum))
	Die("Failed to rearrange CM for bounds calculations");
      NormalizePCM(pcm, statenum);
      LengthDistribution(pcm, statenum, window, &lmx);
      LengthBounds(lmx, statenum, window, epsilon, &minb, &maxb);
      Free2DArray(lmx, statenum);
      free(pcm);
    }

  ext_ofp = NULL;
  if (outfile != NULL)
    if ((ext_ofp = fopen(outfile, "w")) == NULL)
      Die("Failed to open score output file %s for writing", outfile);

  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);
  printf("     version %s, %s\n\n", RELEASE, RELEASEDATE);

  printf("---------------------------------------------------\n");
  printf("Database to search/score:      %s\n", seqfile);
  printf("Model:                         %s\n", cmfile);
  printf("Reporting threshold:           %.2f\n", thresh);
  printf("Maximum match size:            %d\n", window);
  printf("Complementary strand searched: %s\n", do_complement? "yes":"no");
  if (outfile != NULL)
    printf("Scores saved to file:          %s\n", outfile);  
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
      s2upper(seq);

		/* some communication to report_local through statics */
      ext_seqname   = sqinfo.name;
      ext_seqlen    = sqinfo.len;
      in_complement = FALSE;

      if (fast)
	{
	  if (! FastViterbiScan(icm, statenum, minb, maxb, seq, window, thresh, print_hit))
	    Die("Search across sequence %s failed", sqinfo.name);
	}
      else if (! ViterbiScan(icm, statenum, seq, window, thresh, print_hit))
	Die("Search across sequence %s failed", sqinfo.name);

      if (do_complement)
	{
	  in_complement = TRUE;
	  
	  if ((rev = (char *) malloc ((sqinfo.len+1) * sizeof(char))) == NULL)
	    Die("malloc failed");
	  revcomp(rev, seq);
	      
	  if (fast)
	    {
	      if (! FastViterbiScan(icm, statenum, minb, maxb, rev, window, thresh, print_hit))
		Die("Search across sequence %s failed", sqinfo.name);
	    }
	  else if (! ViterbiScan(icm, statenum, rev, window, thresh, print_hit))
	    Die("Search across complement of sequence %s failed", sqinfo.name);

	  free(rev);
	}
      
      if (donefile != NULL)
	{
	  FILE *dfp;
	  if ((dfp = fopen(donefile, "w")) == NULL)
	    Die("Failed to open file for saving name of last finished sequence");
	  fprintf(dfp, "%s\n", sqinfo.name);
	  fclose(dfp);
	}


      FreeSequence(seq, &sqinfo);
#ifdef MEMDEBUG
      curr_size = malloc_size(&histid2);

      if (curr_size != orig_size)
	{
	  Warn("memory leak: current size %ul, starting size %ul\n", curr_size, orig_size);
	  malloc_list(2,histid1, histid2);
	}
#endif
    }

  if (donefile != NULL)
    {
      FILE *dfp;
      if ((dfp = fopen(donefile, "w")) == NULL)
	Die("Failed to open file for saving name of last finished sequence");
      fprintf(dfp, "Search complete.");
      fclose(dfp);
    }


  if (fast)
    {
      free(minb);
      free(maxb);
    }
  free(icm);
  FreeCM(cm);

  SeqfileClose(dbfp);
  return 0;
}


/* Function: print_hit()
 * 
 * Purpose:  Simple function that determines the format of printing
 *           out scanning hits. It gets the start point, end point,
 *           and score of the match; any other info it wants to print
 *           must come through static external variables in this file.
 *           
 * Args:     i  - start point of match
 *           j  - end point of match
 *           sc - score of match
 *
 * Return:   (void)
 */                          
static int
print_hit(int    i,
	  int    j,
	  double sc)
{
  if (in_complement)
    printf("%6.2f  %5d  %5d   : %s\n",  sc, ext_seqlen-i+1, ext_seqlen-j+1, ext_seqname);
  else
    printf("%6.2f  %5d  %5d   : %s\n",  sc, i, j, ext_seqname);

  if (ext_ofp != NULL)
    {
      if (in_complement)
	 fprintf(ext_ofp, "%6.2f  %5d  %5d   : %s\n",  
		 sc, ext_seqlen-i+1, ext_seqlen-j+1, ext_seqname);
      else
	fprintf(ext_ofp, "%6.2f  %5d  %5d   : %s\n",  sc, i, j, ext_seqname);
      fflush(ext_ofp);
    }
  return 1;
}
