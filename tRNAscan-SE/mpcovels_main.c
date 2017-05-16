/* mpcovels_main.c
 * Sun Aug 21 12:25:52 1994
 * 
 * main() for MasPar version of covels
 * 
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
#include "maspar.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define OPTIONS "cg:ho:t:A:D:"

static char usage[]  = "\
Usage: mpcovels [-options] <CM file> <seqfile>\n\
where options are:\n\
    -c          : do complementary strand too\n\
    -g <fracGC> : set background expected GC content (0.5 default)\n\
    -h          : print short help and version info\n\
    -o <file>   : save hits in <file>\n\
    -t <thresh> : set score reporting threshold\n\
 CRASH PROTECTION OPTIONS:\n\
    -A <filename> : maintain file of names of active seqs\n\
    -D <filename> : save names of finished seqs here\n";

static char banner[] = "mpcovels - scan sequences for matches to an RNA covariance model";


extern int MPViterbiScan(struct istate_s *fe_icm, int *fe_statenum,
			 int  *fe_threshold, char *fe_buffer);

/* Because we do some back-and-forth communication between the 
 * DPU and the front end, I need to make the following variables
 * available externally to the functions that the DPU calls back to.
 */
static SQFILE  *fp;             /* open sequence file                            */
static int      format;		/* format of sequence file                       */
static char   **seq;		/* array of NYPROC sequences                     */
static SQINFO  *sqinfo;         /* array of NYPROC sequence info structures      */
static char   **sptr;           /* ptrs into NYPROC active sequences             */
static char    *buffer;         /* NYPROC blocks of processed seq to go to PE's  */
static int      moreseqs;	/* TRUE if there's still more seqs in the file   */
static int      morebases;	/* TRUE if there's still more bases in **seq     */
static int      ithresh;	/* scaled integer score threshold                */
static FILE    *ofp;		/* output file for scores                        */
static int      do_revcomp;	/* TRUE to do reverse complements too            */
static int      inrev[NYPROC];  /* TRUE if we're doing a revcomp in this row now */
static char    *activefile;	/* name of active seq file to save to, or NULL   */
static FILE    *donefp;		/* open finished seq name file, or NULL          */


int
main(int argc, char **argv)
{ 
  char    *seqfile;             /* sequence file                                 */
  char    *cmfile;              /* file containing covariance model        */
  struct cm_s *cm;              /* model                                   */
  struct istate_s *icm;		/* integer log-odds model                  */
  int      statenum;		/* number of states in icm                 */
  int      y;			/* counter for rows */
  double   rfreq[ALPHASIZE];	/* expected background symbol frequencies  */

  char    *outfile;		/* save file for scores                    */
  char    *donefile;		/* save file for finished seq names        */
  float    thresh;		/* threshold score for reporting a match   */
  double   gcfrac;		/* fraction GC expected background         */

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, curr_size;
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/

  thresh        = 0.0;
  do_revcomp    = FALSE;
  outfile       = NULL;
  activefile    = NULL;
  donefile      = NULL;
  gcfrac        = 0.5;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'c': do_revcomp = TRUE;                 break;
    case 'g': gcfrac     = (double) atof(optarg);break;      
    case 'o': outfile    = optarg;               break;
    case 't': thresh     = (float) atof(optarg); break;

    case 'A': activefile = optarg;               break;
    case 'D': donefile   = optarg;               break;

    case 'h': 
      printf("%s\n  version %s (%s)\n%s\n", 
	     banner, RELEASE, RELEASEDATE, usage);
      exit(0);

    default: Die("unrecognized option %c\n", optc);
    }

  if (argc - optind != 2)
    Die("%s\n", usage);
  
  cmfile  = argv[optind]; optind++;
  seqfile = argv[optind]; 

  /* The random model probabilities
   */
  rfreq[1] = rfreq[2] = gcfrac / 2.0;
  rfreq[0] = rfreq[3] = (1.0 - gcfrac) / 2.0;

  ofp = stdout;
  if (outfile != NULL && (ofp = fopen(outfile, "w")) == NULL)
    Die("Failed to open output file %s", outfile);

  donefp = NULL;
  if (donefile != NULL && (donefp = fopen(donefile, "w")) == NULL)
    Die("Failed to open finished sequence names file %s", donefile);

  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);
  printf("     version %s, %s\n\n", RELEASE, RELEASEDATE);

  printf("---------------------------------------------------\n");
  printf("Database to search/score:      %s\n", seqfile);
  printf("Model:                         %s\n", cmfile);
  printf("Reporting threshold:           %.2f\n", thresh);
  if (outfile != NULL)
    printf("Scores saved to file:          %s\n", outfile);  
  printf("Reverse complement too?        %s\n", do_revcomp ? "yes" : "no");
  printf("---------------------------------------------------\n");
  puts("");

  /*********************************************** 
   * Get the model, open the sequence database
   ***********************************************/

  if (! ReadCM(cmfile, &cm))	
    Die("Failed to read model from file %s", cmfile);

  if (! RearrangeCM(cm, rfreq, &icm, &statenum))
    Die("Failed to convert CM to integer log-odds form\n");  

  if (! SeqfileFormat(seqfile, &format, NULL)) 
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Sequence file %s could not be opened for reading", seqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of sequence file %s", seqfile);
    }

  if ((fp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);

  ithresh = (int) (thresh * INTPRECISION);

  /* Load the first set of sequences
   */
  if ((seq   = (char **)  malloc (sizeof(char *) * NYPROC)) == NULL ||
      (sqinfo= (SQINFO *) malloc (sizeof(SQINFO) * NYPROC)) == NULL ||
      (sptr  = (char **)  malloc (sizeof(char *) * NYPROC)) == NULL ||
      (buffer= (char *)   malloc (sizeof(char) * (NYPROC * BLOCKSIZE))) == NULL)
    Die("malloc failed");

  moreseqs  = TRUE;
  morebases = TRUE;
  for (y = 0; y < NYPROC; y++)
    {
      if (! moreseqs || ! ReadSeq(fp, format, &seq[y], &(sqinfo[y])))
	{ moreseqs = FALSE; seq[y] = NULL; }
      if (seq[y] != NULL) s2upper(seq[y]);
      sptr[y] = seq[y];
      inrev[y] = FALSE;
    }

  /* Call the parallel scanning function.
   * It will load the model up, then immediately call us
   * back and ask us to prepare the first block of sequence data.
   * We tell it where the sequence data is by passing the address of buffer.
   */
  callRequest(MPViterbiScan, 16, icm, &statenum, &ithresh, buffer);

  /* Cleanup and exit
   */
  free(seq);
  free(sqinfo);
  free(sptr);
  free(buffer);
  SeqfileClose(fp);
  fclose(ofp);
  if (donefp != NULL) fclose(donefp);
  
  return 0;
}



/* Function: NextSequenceBlock()
 * 
 * Purpose:  Prepare a sequence buffer for the DPU.
 *           The buffer consists of NYPROC blocks of size BLOCKSIZE.
 *           The values are 0,1,2,3 for A,C,G,U, or 4 for no seq.
 *           If we return 1, the DPU blockIn()'s this buffer.
 *           If we return 0, the DPU concludes that the search is complete.
 *           
 * Note:          
 *    There's some intricacy below, so be careful what you muck with.
 *    These next two functions are written carefully to make sure that
 *    sqinfo[y] is still valid for the hits that are reported to
 *    DPUReportsHit(). This has the following consequences:
 *       1) NextSequenceBlock() must make sure that there's only one
 *          sequence in a BLOCKSIZE. It must wait for the next block
 *          to start sending a new sequence.
 *       2) Moreover, it's got to make sure that the terminating '\0'
 *          was part of the block -- so the DPU can clean up and report
 *          the final hit on that sequence *before* we free it and its
 *          info. Hence, the delayed FreeSequence() call here.
 *
 *    This results in some inefficiency, since we can't start new sequences
 *    until the next block.
 */
int
NextSequenceBlock(void)
{
  int changed;			/* TRUE if active sequence list changed */
  int y, pos;
  char *rev;

  /* Load new sequences if we have to.
   */
  morebases = FALSE;
  changed   = FALSE;
  for (y = 0; y < NYPROC; y++)
    {
				/* where we have a seq, but now we're 
				 * done with it...
				 * free the old seq, load a new one */
      if (seq[y] != NULL && sptr[y] == NULL)
	{
	  if (do_revcomp && !inrev[y])
	    {
	      if ((rev = (char *) malloc (sizeof(char) * (sqinfo[y].len + 1))) == NULL)
		Die("malloc failed");
	      if (revcomp(rev, seq[y]) == NULL)
		Die("revcomp failed");
	      free(seq[y]);
	      seq[y] = sptr[y] = rev;
	      inrev[y] = TRUE;
	      changed  = TRUE;
	    }
	  else
	    {
	      if (donefp != NULL) /* save name to finished seq file */
		{ fprintf(donefp, "%s\n", sqinfo[y].name); fflush(donefp); }
	      FreeSequence(seq[y], &sqinfo[y]);
	      if (! moreseqs || ! ReadSeq(fp, format, &seq[y], &sqinfo[y]))
		{ moreseqs = FALSE; seq[y] = NULL; }
	      if (seq[y] != NULL) s2upper(seq[y]);
	      sptr[y] = seq[y];
	      inrev[y] = FALSE;
	      changed = TRUE;
	    }
	}
      if (seq[y] != NULL) morebases = TRUE;
    }
  if (morebases == FALSE) return 0;
  
  /* Construct a new block of sequence data for the DPU to blockIn().
   */
  for (y = 0; y < NYPROC; y++)
    {
      for (pos = 0; pos < BLOCKSIZE; pos++)
	{
	  if (seq[y] == NULL || sptr[y] == NULL)
	    buffer[y * BLOCKSIZE + pos] = 4;
	  else if (*(sptr[y]) == '\0')
	    {
	      buffer[y * BLOCKSIZE + pos] = 4;
	      sptr[y] = NULL;	/* signal that we're done w/ this seq */
	    }
	  else 
	    {
	      buffer[y * BLOCKSIZE + pos] = SymbolIndex(*(sptr[y]));
	      sptr[y] ++;
	    }
	}
    }
				/* crash protection: maintain list of active names */
  if (changed == TRUE && activefile != NULL)
    {
      FILE *actfp;

      if ((actfp = fopen(activefile, "w")) == NULL)
	Die("failed to open active sequence name list file %s", activefile);
      for (y = 0; y < NYPROC; y++)
	fprintf(actfp, "%s\n", sqinfo[y].name);
      fclose(actfp);
    }

  return 1;
}


/* Function: DPUReportsHit()
 * 
 * Purpose:  The DPU is reporting one or more hits (up to NYPROC).
 *           The easiest thing to do is to blockIn() from *all*
 *           NYPROC rows, and find what rows actually have real
 *           hits here. 
 */
int
DPUReportsHit(int *dpu_seeme, int *dpu_i, int *dpu_j, int *dpu_score)
{
  int seeme[NYPROC];		/* TRUE if we're reporting in this row */
  int start[NYPROC];
  int end[NYPROC];
  int score[NYPROC];
  int y;

  blockIn(dpu_seeme, seeme, 0, 0, 1, NYPROC, sizeof(int));
  blockIn(dpu_i,     start, 0, 0, 1, NYPROC, sizeof(int));
  blockIn(dpu_j,     end,   0, 0, 1, NYPROC, sizeof(int));
  blockIn(dpu_score, score, 0, 0, 1, NYPROC, sizeof(int));

  for (y = 0; y < NYPROC; y++)
    if (seeme[y])
      {
	if (inrev[y])
	  fprintf(ofp, "%6.2f  %5d  %5d   : %s %s\n",  
		  (float) score[y] / INTPRECISION, 
		  sqinfo[y].len - start[y], sqinfo[y].len - end[y], 
		  sqinfo[y].name,
		  (sqinfo[y].flags & SQINFO_DESC) ? sqinfo[y].desc : "");
	else
	  fprintf(ofp, "%6.2f  %5d  %5d   : %s %s\n",  
		  (float) score[y] / INTPRECISION, 
		  start[y] + 1, end[y] + 1, 
		  sqinfo[y].name,
		  (sqinfo[y].flags & SQINFO_DESC) ? sqinfo[y].desc : "");
	fflush(ofp);
      }
  return 1;	
}
