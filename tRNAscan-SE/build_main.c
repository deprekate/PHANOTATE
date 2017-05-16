/* build_main.c
 * SRE, Mon Sep  6 09:18:35 1993
 * 
 * coveb - construct a covariance model from aligned sequences
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

#define OPTIONS "ap:g:ho:P:"
static void mtr2rf(struct trace_s *mtr, int alen, char **ret_rf);

static char usage[]  = "\
Usage: coveb [-options] <cmfile output> <seqfile in>\n\
where options are:\n\
     -a              : annotate all pairs, not just canonical Watson-Crick\n\
     -g <gapthresh>  : -P1 or -P2 only - columns above this fractional gap \n\
                       occurrence are assigned to INS (default 0.7)\n\
     -h              : print out short help info\n\
     -o <file>       : save structure-annotated alignment to <file>\n\
     -p <priorfile>  : use prior probability info from <priorfile>\n\
 Construction plans:\n\
     (default)       : Maximum likelihood (slow)\n\
     -P1             : fast heuristic (MIXY based)\n\
     -P2             : use specified consensus structure (CS line)\n\
     -P3             : use both specified consensus and reference info\n";

static char banner[] = "coveb: construct covariance model from aligned sequences";

int
main(int argc, char **argv)
{ 
  char          **aseqs;        /* training sequences           */
  AINFO           ainfo;        /* misc. associated alignment info */
  int             nseq;		/* number of seqs               */ 
  char           *seqfile;      /* sequence file                */
  int             format;	/* format of sequence file      */
  char           *cmfile;       /* OUTPUT: saved cm             */
  FILE           *cmfp;		/* OUTPUT: fp to cmfile         */
  struct cm_s    *cm;           /* model                        */
  struct prior_s *prior;        /* prior prob. distributions    */
  int             idx;		/* index for sequences          */
  double          secinfo;      /* secondary structure info content */
  struct trace_s *mtr;		/* master traceback for alignment */
  struct trace_s *tr;		/* a traceback for indiv seq      */
  struct trmem_s *pool;         /* memory pool for traceback      */
  char          **ss;           /* secondary structures           */
  double          worstscore;
  double          bestscore;
  double          sqsum;
  double          tot_score;
  double          score;
  int             leftcount, rightcount;
  int             apos;

  enum plan_e { PLAN_ML, PLAN_MIXY, PLAN_CS, PLAN_CSRF } plan;
  char        *prifile;		/* file to get prior from                */
  FILE        *prifp;		/* open priorfile                        */
  char        *structfile;      /* file to save structure-annotation to  */
  FILE        *structfp;        /* open structfile                       */
  double       gapthresh;	/* heuristic INS assignment parameter    */
  int          watsoncrick;     /* TRUE to annotate canonical pairs only */

  int          optc;		/* for getopt() */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */


#ifdef MEMDEBUG			/* for Cahill's dbmalloc */
  unsigned long histid1, histid2, orig_size, current_size;
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/

  prifile      = NULL;          /* forces use of default prior in prior.h */
  gapthresh    = .70;	
  structfile   = NULL;
  watsoncrick  = TRUE;
  plan         = PLAN_ML;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'a': watsoncrick  = FALSE;        break;
    case 'g': gapthresh    = atof(optarg); break;
    case 'o': structfile   = optarg;       break;
    case 'p': prifile      = optarg;       break;

    case 'P':
      switch (*optarg) {
      case '1': plan = PLAN_MIXY; break;
      case '2': plan = PLAN_CS;   break;
      case '3': plan = PLAN_CSRF; break;
      default: Die("No such construction plan.\n%s", usage);
      }
      break;

    case 'h': 
      printf("%s\n  version %s (%s)\n%s", banner, RELEASE, RELEASEDATE, usage);
      exit(0);

    default:
      Die("Error: unrecognized option %c\n", optc);
    }

  if (argc - optind != 2)
    Die("Wrong number of command line arguments.\n%s\n", usage);
  
  cmfile  = argv[argc-2];
  seqfile = argv[argc-1]; 

#ifdef MEMDEBUG
  orig_size = malloc_size(&histid1);
#endif


  /*********************************************** 
   * Get sequence data and prior
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, NULL))
    Die("Can't determine format of file %s", seqfile);

				/* read the training seqs from file */
  if (! ReadAlignment(seqfile, format, &aseqs, &nseq, &ainfo))
    Die("Failed to read aligned sequence file %s", seqfile);

				/* convert seq to all upper case */
  for (idx = 0; idx < nseq; idx++)
    s2upper(aseqs[idx]);

  if (prifile == NULL)
    {
      if (! DefaultPrior(&prior))
	Die("Failed to copy prior probability distribution information");
    }
  else
    {
      if ((prifp = fopen(prifile, "r")) == NULL)
	Die("Failed to open prior probability info file %s", prifile);
      if (! ReadPrior(prifp, &prior))
	Die("Failed to read prior probability info file %s", prifile);
      fclose(prifp);
    }

  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);
  printf("     release %s, %s\n\n", RELEASE, RELEASEDATE);
  printf("---------------------------------------------------\n");
  printf("Training alignment:     %s (%d sequences)\n", seqfile, nseq);
  printf("Prior distributions:    ");
  if (prifile == NULL)
    printf("default (plus-one, Laplace)\n");
  else
    printf("from file %s\n", prifile);
  
  printf("Construction plan:      ");
  switch (plan) {
  case PLAN_ML:   printf("maximum likelihood\n");                    break;
  case PLAN_MIXY: printf("fast (MIXY-based) heuristic\n");           break;
  case PLAN_CS:   printf("specified consensus structure\n");         break;
  case PLAN_CSRF: printf("specified structure and match columns\n"); break;
  }

  if (plan == PLAN_MIXY || plan == PLAN_CS)
    printf("INS if gap freq >:      %.2f\n", gapthresh);

  printf("---------------------------------------------------\n");


  /*********************************************** 
   * Create the starting model 
   ***********************************************/

#ifdef MEMDEBUG
  printf("Checking malloc chain ");
  malloc_chain_check(0);
  printf("... done.\n");
#endif

  switch (plan) {
  case PLAN_ML: 
    if (! Maxmodelmaker(aseqs, &ainfo, nseq, -1, prior, &secinfo, &cm, &mtr))
	Die("Maxmodelmaker failed to create starting model from alignment");
    break;

  case PLAN_MIXY:
    if (! Fastmodelmaker(aseqs, &ainfo, nseq, prior, gapthresh, &secinfo, &cm, &mtr))
      Die("Fastmodelmaker failed to create starting model from alignment");
    break;

  case PLAN_CS:
    EasyModelmaker(aseqs, &ainfo, nseq, prior, gapthresh, FALSE, &cm, &mtr);
    break;

  case PLAN_CSRF:
    EasyModelmaker(aseqs, &ainfo, nseq, prior, gapthresh, TRUE, &cm, &mtr);
    break;
    
  default: Die("That must be a secret plan, pal, because I've never heard of it.");
  }

  if (! VerifyCM(cm))
    Die("Bad covariance model. Bad, bad covariance model...\n");


  /* Use master traceback to reconstruct individual traces.
   * Use individual traces to a) find individual secondary structures
   *                          b) calculate avg, high, low scores
   * Some duplication of effort here, because the modelmakers
   * already had to construct tracebacks and threw them away.
   */                         
  if ((ss = (char **) malloc (sizeof(char *) * nseq)) == NULL)
    Die("malloc failed");
  worstscore = HUGE_VAL;
  bestscore  = -HUGE_VAL;
  tot_score   = sqsum = 0.0;
  for (idx = 0; idx < nseq; idx++)
    {
      Transmogrify(mtr, aseqs[idx], &tr, &pool);

      score = TraceScore(cm, aseqs[idx], tr);
      tot_score += score;
      sqsum     += score * score;
      if (score > bestscore)  bestscore  = score;
      if (score < worstscore) worstscore = score;

      if (ainfo.sqinfo[idx].flags & SQINFO_SS) free(ainfo.sqinfo[idx].ss); 
      Trace2KHS(tr, aseqs[idx], ainfo.alen, watsoncrick, &(ss[idx]));
      MakeDealignedString(aseqs[idx], ainfo.alen, ss[idx], &(ainfo.sqinfo[idx].ss));
      ainfo.sqinfo[idx].flags |= SQINFO_SS;

      FreeTrace(tr, pool);
    }

  /* Construct a consensus structure string and reference
   * line as annotation.
   * Secondary structure strings, ss, are currently aligned to
   * the aseqs. Calculate an aligned consensus structure from
   * them. 
   */
  if (plan == PLAN_MIXY || plan == PLAN_ML)
    {
      if (ainfo.flags & AINFO_CS) free(ainfo.cs);
      if ((ainfo.cs = (char *) malloc (sizeof(char) * (ainfo.alen+1))) == NULL)
	Die("malloc failed");
      for (apos = 0; apos < ainfo.alen; apos++)
	{
	  leftcount = rightcount = 0;
	  for (idx = 0; idx < nseq; idx++)
	    if      (ss[idx][apos] == '<') rightcount++;
	    else if (ss[idx][apos] == '>') leftcount++;
	  
	  if      (rightcount >= nseq / 2) ainfo.cs[apos] = '<';
	  else if (leftcount  >= nseq / 2) ainfo.cs[apos] = '>';
	  else    ainfo.cs[apos] = '.';
	}
      ainfo.cs[ainfo.alen] = '\0';
      ainfo.flags |= AINFO_CS;
    }

  /* Construct a reference line to indicate which columns were assigned
   * as match states.
   */
  if (plan != PLAN_CSRF)
    {
      if (ainfo.flags & AINFO_RF) free(ainfo.rf);
      mtr2rf(mtr, ainfo.alen, &(ainfo.rf));
      ainfo.flags |= AINFO_RF;
    }

  if ((cmfp = fopen(cmfile, "w")) == NULL)
    Die("Failed to open %s for writing", cmfile);
  if (! WriteCM(cmfp, cm))
    Die("Failed to save the model to %s", cmfile); 
  fclose(cmfp);

  if (structfile != NULL)
    {
      if ((structfp = fopen(structfile, "w")) == NULL)
	Die("Failed to open structure annotation alignment file %s", structfile);
      if (! WriteSELEX(structfp, aseqs, nseq, &ainfo, 60))
	Die("Failed to write annotated alignment to %s", structfile);
      fclose(structfp);
      printf("Structure annotated alignment file written to %s\n", structfile);
    }

  printf("Constructed a covariance model (%d nodes)\n", cm->nodes);
  printf("Average score:                %10.2f bits\n", 
	 tot_score / (double) nseq);
  printf("Minimum score:                %10.2f bits\n", worstscore);
  printf("Maximum score:                %10.2f bits\n", bestscore);
  printf("Std. deviation:               %10.2f bits\n", 
	 sqrt((sqsum - (tot_score * tot_score / (double) nseq)) 
	      / ((double) nseq - 1.0)));

  printf("\nCM written to file %s\n", cmfile);

  FreeCM(cm);
  FreeTrace(mtr, NULL);
  free(prior);
  FreeAlignment(aseqs, nseq, &ainfo);
  Free2DArray(ss, nseq);

#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "No memory leaks, sir.\n");
#endif

  return 0;
}


/* Function: mtr2rf()
 * 
 * Purpose:  Make an #=RF line from a master traceback, to indicate
 *           which columns were used as match columns in building
 *           a model. Keep in mind that master traces use node type
 *           indices, rather than state type indices like traces
 *           are supposed to.
 */
static void
mtr2rf(struct trace_s *mtr, int alen, char **ret_rf)
{
  struct tracestack_s *dolist;
  struct trace_s      *curr;
  char *rf;

  rf = (char *) MallocOrDie(sizeof(char) * (alen+1));
  memset(rf, ' ', alen);
  rf[alen] = '\0';

  dolist = InitTracestack();
  PushTracestack(dolist, mtr->nxtl);
  while ((curr = PopTracestack(dolist)) != NULL)
    {
      if ( curr->type == MATP_NODE )
	rf[curr->emitl] = rf[curr->emitr] = '.';
      else if ( curr->type == MATL_NODE )
	rf[curr->emitl] = '.';
      else if ( curr->type == MATR_NODE )
	rf[curr->emitr] = '.';

      if (curr->nxtr) PushTracestack(dolist, curr->nxtr);
      if (curr->nxtl) PushTracestack(dolist, curr->nxtl);
    }
  FreeTracestack(dolist);
  *ret_rf = rf;
}
