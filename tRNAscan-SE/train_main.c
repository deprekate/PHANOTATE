/* train_main.c
 * 2.0: Fri Oct  1 16:02:11 1993
 * SRE, Mon Jun 28 17:52:54 1993
 * 
 * main() for covet: training of a covariance hmm from aligned seqs
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

#define OPTIONS "a:A:b:fhG:i:mp:s:X:"

static char usage[]  = "\
Usage: covet [-options] <cmfile output> <seqfile in>\n\
where options are:\n\
     -a <alignfile>  : make starting model from alignment\n\
     -A <filename>   : save alignments to filename.1, etc., for animation\n\
     -b <backupfile> : each iteration, back up curr model to <backupfile>\n\
     -f              : use flat text save formats, portable but clumsy\n\
     -G <GOP>        : gap-open prob 0 < gop < 1 for random alignment generation\n\
     -h              : print short help and version info\n\
     -i <cm file>    : take start model from <cm file>\n\
     -m              : do maximum likelihood model construction (slow!)\n\
     -p <prior file> : use prior in <file>; default is Laplace plus-one\n\
     -s <seed>       : set random() seed\n\
     -X <GEX>        : gap-extend prob 0 < gex < 1 for random alignment generation\n";

static char banner[] = "covet: training of a covariance model";

int
main(int argc, char **argv)
{ 
  char       **rseqs;		/* training sequences                        */
  SQINFO      *sqinfo;		/* array of sqinfo structures for rseqs      */
  AINFO        ainfo;           /* alignment info                            */
  int          nseq;		/* number of seqs                            */ 
  char        *seqfile;         /* training sequence file                    */
  int          format;          /* seqfile format                            */
  char        *cmfile;          /* OUTPUT: saved cvhmm                       */
  FILE        *cmfp;		/* OUTPUT: fp to cvfile                      */
  struct cm_s *cm;              /* model                                     */
  struct cm_s *newcm;           /* new model                                 */
  struct istate_s *icm;         /* model, integer log odds form              */
  int          statenum;        /* # of states in the model                  */
  struct prior_s *prior;        /* prior prob. distributions                 */
  int          idx;		/* counter for sequences                     */
  double       score;		/* score of individual alignment             */
  double       totscore;	/* summed scores over training seqs          */
  double       oldscore;	/* previous totscore for old model           */
  double       delta;	        /* fractional change in scores for iteration */
  int          iteration;	/* iteration number we're on                 */
  struct trace_s **tr;          /* tracebacks for each sequence              */
  double rfreq[ALPHASIZE];	/* expected background symbol frequencies    */

  int    max_iterations;
  double threshold;		/* fractional tolerance, test for convergence   */
  char  *aseqfile;              /* sequence alignment for model building        */
  char **aseqs;                 /* aligned sequences in aseqfile                */
  int    num_aseqs;		/* number of aseqs in aseqfile                  */
  char  *in_cmfile;		/* file containing input model                  */
  char  *bckfile;		/* backup file for saving models each iteration */
  FILE  *bckfp;		        /* open pointer to backup file                  */
  int    do_flatformat;		/* TRUE if we save in flat text format          */
  int    seed;			/* seed for random()                            */
  char  *animateroot;		/* root name for alignment animation save files */
  char   animationfile[256];	/* full animation file name (root.1, etc.)      */
  FILE  *animfp;
  double random_open;		/* insert-open probability for random alignment */
  double random_extend;		/* ins-extend probability for random alignment  */
  int    fast_version;		/* do Fastmodelmaker(), not Maxmodelmaker()     */
  char  *prifile;               /* file to obtain prior from                    */
  FILE  *prifp;			/* open prifile for reading                     */
  int    watsoncrick;           /* if TRUE, annotate only canonical pairs       */

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/

  threshold            = 0.01;		/* default: 1% */
  max_iterations       = 100;
  aseqfile             = NULL;
  in_cmfile            = NULL;
  bckfile              = NULL;
  do_flatformat        = 0;
  seed                 = (int) time (0); /* default: "random" seed */
  animateroot          = NULL;
  random_open          = 0.02;
  random_extend        = 0.39;
  fast_version         = TRUE;
  prifile              = NULL;
  watsoncrick          = TRUE;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'a': aseqfile      = optarg;       break;
    case 'A': animateroot   = optarg;       break;
    case 'b': bckfile       = optarg;       break;
    case 'f': do_flatformat = 1;            break;
    case 'G': random_open   = atof(optarg); break;
    case 'i': in_cmfile     = optarg;       break;
    case 'm': fast_version  = FALSE;        break;
    case 'p': prifile       = optarg;       break;
    case 's': seed          = atoi(optarg); break;
    case 'X': random_extend = atof(optarg); break;
    case 'h': 
      printf("%s\n  version %s (%s)\n%s\n", banner, RELEASE, RELEASEDATE, usage);
      exit(0);
    default:
      Die("Error: unrecognized option %c\n", optc);
    }

  if (argc - optind != 2)
    Die("%s\n", usage);
  
  if (aseqfile != NULL && in_cmfile != NULL)
    Die("options -i and -a are exclusive\n%s", usage);

  cmfile  = argv[argc-2];
  seqfile = argv[argc-1]; 
  sre_srandom(seed);

#ifdef MEMDEBUG
  orig_size = malloc_size(&histid1);
#endif

  /*********************************************** 
   * Get sequence data and a prior
   ***********************************************/
				/* random background model */
  rfreq[0] = rfreq[1] = rfreq[2] = rfreq[3] = 0.25;

  if (! SeqfileFormat(seqfile, &format, NULL))
    Die("Failed to determine format of sequence file %s", seqfile);

				/* read the training seqs from file */
  if (! ReadMultipleRseqs(seqfile, format, &rseqs, &sqinfo, &nseq))
    Die("Failed to read sequences from file %s", seqfile);

  for (idx = 0; idx < nseq; idx++)
    PrepareSequence(rseqs[idx]);

  if (prifile == NULL)
    {
      if (! DefaultPrior(&prior))
	Die("Failed to copy prior probability distribution information");
    }
  else
    {
      if ((prifp = fopen(prifile, "r")) == NULL)
	Die("Failed to open prior probability file %s", prifile);
      if (! ReadPrior(prifp, &prior))
	Die("Failed to read prior probabilities from %s", prifile);
      fclose(prifp);
    }
      
  /*********************************************** 
   * Create the starting model 
   ***********************************************/

  if (aseqfile != NULL)		/* A start from an alignment */
    {
      if (! SeqfileFormat(aseqfile, &format, NULL))
	Die("Failed to determine format of seed alignment file %s", aseqfile);

      if (! ReadAlignment(aseqfile, format, &aseqs, &num_aseqs, &ainfo))
	Die("Failed to read alignment from %s", aseqfile);
      
      for (idx = 0; idx < num_aseqs; idx++)
	s2upper(aseqs[idx]);
      
      if (fast_version)
	{
	  if (! Fastmodelmaker(aseqs, &ainfo, nseq, prior, 0.70, NULL, &cm, NULL))
	    Die("Fastmodelmaker failed to create starting model from alignment");
	}
      else
	{
	  if (! Maxmodelmaker(aseqs, &ainfo, num_aseqs, -1.0, prior, NULL, &cm, NULL))
	    Die("Failed to create starting model from alignment");
	}
      FreeAlignment(aseqs, num_aseqs, &ainfo);
    }

  else if (in_cmfile)		/* A start from an existing model */
    {
      if (! ReadCM(in_cmfile, &cm))
	Die("Failed to read starting model from file %s", in_cmfile);
    }

  else				/* A start from a flat model */
    {
      RandomAlignment(rseqs, sqinfo, nseq, random_open, random_extend, &aseqs, &ainfo);
      if (fast_version)
	{
	  if (! Fastmodelmaker(aseqs, &ainfo, nseq, prior, 0.70, NULL, &cm, NULL))
	    Die("Fastmodelmaker failed to create starting model from alignment");
	}
      else
	{
	  if (! Maxmodelmaker(aseqs, &ainfo, nseq, -1.0, prior, NULL, &cm, NULL))
	    Die("Failed to create starting model from alignment");
	}
      FreeAlignment(aseqs, nseq, &ainfo);
     }

  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);
  printf("     release %s, %s\n\n", RELEASE, RELEASEDATE);

  printf("---------------------------------------------------\n");
  printf("Training data:          %s (%d sequences)\n", seqfile, nseq);
  if (aseqfile != NULL)
    printf("Starting model:         from alignment in %s (%d seqs)\n", aseqfile, num_aseqs);
  else if (in_cmfile != NULL )
    printf("Starting model:         from existing model in %s\n", in_cmfile);
  else
    printf("Starting model:         random alignment\n");
  printf("Prior distributions:    %s\n", prifile == NULL ? "plus-one" : prifile);
  printf("Modelmaking strategy:   %s\n", fast_version ? "fast heuristic" : "max likelihood");
  printf("Convergence threshold:  %.4f\n", threshold);
  printf("Maximum iterations:     %d\n", max_iterations);
  if (bckfile != NULL)
    printf("Backup model file:      %s\n", bckfile);
  printf("seed for random():      %d\n", seed);
  printf("---------------------------------------------------\n");
  puts("");


  /*********************************************** 
   * Train model by expectation maximization
   ***********************************************/

  if ((tr = (struct trace_s **) malloc (nseq * sizeof(struct trace_s *))) == NULL)
    Die("Memory failure, line %d of %s", __LINE__, __FILE__);
  
  oldscore = -1.0 * HUGE_VAL;
  iteration = 0;
  while (iteration < max_iterations)
    {
      iteration++;
      printf("Iteration %4d  : model of %d nodes, ", iteration, cm->nodes);

      /* Make a search model
       */
      if (! RearrangeCM(cm, rfreq, &icm, &statenum))
	Die("Failed to make an integer log-odds model");

      /* First we align all the sequences to the model,
       * and construct a multiple sequence alignment
       */
      totscore = 0.0;
      for (idx = 0; idx < nseq; idx++)
	{
	  if (! ViterbiAlign(icm, statenum, rseqs[idx], &score, &tr[idx]))
	    Die("viterbi alignment failed on sequence %d", idx);
	  totscore += score;
	}

      /* An option for producing cool figures and animations:
       * save the alignment at each iteration, so we can animate
       * the learning process. I produced covariance matrices
       * with MIXY for each iteration, and used GNUPLOT to
       * animate the data as a series of 3D surface plots.
       */
      if (animateroot != NULL)
	{
	  sprintf(animationfile, "%s%d", animateroot, iteration);
	  if ((animfp = fopen(animationfile, "w")) == NULL)
	    Warn("Failed to open animation output file %s", animationfile); 
	  else
	    {
	      if (! Traces2Alignment(rseqs, sqinfo, tr, nseq, cm, watsoncrick, &aseqs, &ainfo))
		Warn("Traces2Alignment() failed for animation");

	      WriteSELEX(animfp, aseqs, nseq, &ainfo, 60);
	      FreeAlignment(aseqs, nseq, &ainfo);
	      fclose(animfp);
	    }
	}
	      

      /* If we've converged, stop.
       * Else, make a new model from the alignment. 
       */
      delta = (totscore - oldscore) / fabs(totscore);
      printf("score %.3f, delta %.3f\n", totscore / (double) nseq, delta);
      
      if (delta > threshold || delta < 0) 
	{
	  if (! Traces2Alignment(rseqs, sqinfo, tr, nseq, cm, watsoncrick, &aseqs, &ainfo))
	    Die("Traces2Alignment() failed");

	  for (idx = 0; idx < nseq; idx++)
	    s2upper(aseqs[idx]);

	  if (fast_version)
	    { 
	      if (! Fastmodelmaker(aseqs, &ainfo, nseq, prior, 0.70, NULL, &newcm, NULL))
		Die("Failed to create new model from alignment");
	    }
	  else
	    {
	      if (! Maxmodelmaker(aseqs, &ainfo, nseq, -1.0, prior, NULL, &newcm, NULL))
		Die("Failed to create new model from alignment");
	    }

	  FreeAlignment(aseqs, nseq, &ainfo);
	  for (idx = 0; idx < nseq; idx++)
	    FreeTrace(tr[idx], NULL);

	  oldscore = totscore;
	}
      else
	{
	  /* we've converged. Free traces and break out of iteration loop.
	   */
	  for (idx = 0; idx < nseq; idx++)
	    FreeTrace(tr[idx], NULL);
	  break;
	}

      /* switch new model for old 
       */
      FreeCM(cm);
      free(icm);
      cm = newcm;

      /* Training takes a long time, and sysadmins are nasty evil 
       * people who like to shut machines down without warning, particularly
       * during a long training run right before one gives a talk.
       * Therefore, we have an option for backing up the model every
       * iteration so we can resume after a crash or shutdown...
       */
      if (bckfile != NULL)
	{
	  if ((bckfp = fopen(bckfile, "w")) == NULL) 
	    Warn("Failed to open backup file %s\n", bckfile); 
	  else
	    { 
	      if (! WriteBinaryCM(bckfp, newcm))
		Warn("Failed to save to backup file %s\n", bckfile);
	      fclose(bckfp);
	    }
	}
    }


  /*********************************************** 
   * Save the new model and exit.
   ***********************************************/

  if ((cmfp = fopen(cmfile, "w")) == NULL)
    Die("Failed to open %s for writing", cmfile);

  if (do_flatformat && ! WriteCM(cmfp, cm))
    Die("Failed to save the model to %s", cmfile); 
  else if (! WriteBinaryCM(cmfp, cm))
    Die("Failed to save the model to %s", cmfile);
  fclose(cmfp);
  
  free(tr);
  FreeCM(cm);
  for (idx = 0; idx < nseq; idx++)
    FreeSequence(rseqs[idx], &(sqinfo[idx]));
  free(sqinfo);

  printf("New covariance model written to file %s\n", cmfile);

#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "No memory leaks, sir.\n");
#endif
    
  return 0;
}



