
/* eufindtRNA  - Eukaryotic tRNA finder
 *
 * C implementation of algorithm described by Pavesi, Conterio, 
 * Bolchi, Dieci, & Ottonello in NAR 22:1247-56 (94)
 * "Identification of new eukaryotic tRNA genes in genomic DNA
 * databases by a multistep weight matix analysis of transcriptional
 * control regions"
 *
 * To be used in tRNAscan-SE package to increase sensitivity by
 * complementing tRNAscan 1.3 first-pass scan
 *
 * by Todd MJ Lowe    4/8/96
 *
 * Uses Sean Eddy's function library for biological sequence analysis
 * (Squid v1.5g)
 *
 * v1.1: Small bug fixed 8/2000 that caused second of two consecutive tRNAs 
 * (within 40bp) to be missed if the second tRNA scored lower than the first
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "squid.h"
#include "sqfuncs.h"
#include "eufind_const.h"

char eufind_version[] = "1.1";
char eufind_date[]    = "Aug 2000";


#define OPTIONS "ho:l:X:I:rsDFi:"
char usage[] = "\n\
Usage: eufindtRNA [-options] <sequence file>\n\
Find tRNAs in eukaryotic sequences\n\n\
Available options:\n\
-h             : help - print version and usage info\n\
-o <outfile>   : save tRNAs in <outfile>\n\
-r             : relaxed mode (no terminators searched for)\n\
-s             : strict mode (discard tRNAs with missing terminators)\n\
-l <length>    : set max intron+variable loop length (default=140)\n\
-X <Score>     : manually set final score cutoff to <Score> (def=-31.8)\n\
-I <Score>     : set cutoff for intermediate score (def=-31.25)\n\
-D             : save tRNA score components (for Debugging)\n\
-F             : don't check for uppercase or DNA alphabet\n\
-i <integer>   : start nucleotide numbering at <integer> (def=1)\n\n";

int
main (int argc, char **argv)
{
  char  *seqfile;               /* file containing aligned seqs */
  char  *outfile;               /* destination file for tRNAs   */
  int       fmt;		/* format of seqfile  */
  FILE     *outfp;                /* open outfile                 */
  SQFILE   *seqfp; 
  SQINFO    sqinfo;
  int       i, errno,
    ShowScores,                 /* flag for type of info output when
				   saving tRNAs */
    RelaxedMode,                /* flag for relaxed scanning mode, do
				   not look for poly T terminator
				   signal */
    StrictMode,                 /* require poly T terminator */
    NoReformat;                 /* flag to prevent extra work of 
				/*   changing seqs to DNA & upper case */
				/*   alphabet */

  float  NoTermPenalty;         /* penalty val for tRNAs with no */
                                /* poly T terminator */

  long int sqoffset;           /* nucleotide numbering offset (set with -i param) */
  char *seq, *revseq,                  /* sequence */
    *iseq, *reviseq;           /* encoded seq & encoded reverse comp */
  int strand,                  /* 1 for orig seq, -1 for rev comp */
    seqidx;                    /* current position in seq */
  float FirstScore,            /* initial (Bbox) logodds score  */ 
    IntScore, TotScore,        /* cum tRNA logodds scores */
    IntScoreCutoff,
    TotScoreCutoff;            /* cutoff for reporting tRNAs */
  TRNA_TYPE *tRNA, *prev_tRNA,
    *swap_tRNA;                 /* current & previous tRNA info */
  int Max_AB_dist;              /* max nuc. distance searched upstream */
				/* of candidate B boxes for A boxes */
  int tRNA_ct;

  int          optchar;
  extern char *optarg;
  extern int   optind;

  /***********************************************
   * Parse command line
   ***********************************************/

  outfile    = NULL;
  TotScoreCutoff = TOT_SCORE_THRESH;
  IntScoreCutoff = INT_SCORE_THRESH;
  Max_AB_dist = MIN_AB_BOX_DIST + AB_BOX_DIST_RANGE;
  sqoffset = 0;
  ShowScores = 0;
  NoTermPenalty = MAX_PENALTY;
  RelaxedMode = 0;
  StrictMode = 0;
  NoReformat = 0;
  
  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {

    case 'o': outfile    = optarg; break;
    case 'h': 
      printf("eufindtRNA %s, %s\n%s\n", eufind_version, eufind_date, usage);
      exit(EXIT_SUCCESS);
    case 'r': RelaxedMode = 1; break;
    case 's': NoTermPenalty = 10*MAX_PENALTY; StrictMode = 1; break;
    case 'X': TotScoreCutoff = atof(optarg); break;
    case 'D': ShowScores = 1; break;
    case 'F': NoReformat = 1; break;
    case 'l': Max_AB_dist = MIN_AB_BOX_DIST + atof(optarg); break;
    case 'I': IntScoreCutoff = atof(optarg); break;
    case 'i': sqoffset = atof(optarg)-1; break;
    default:
      Die("%s\n", usage);
    }

  if (argc -optind != 1)
    Die("Wrong number of arguments specified on command line\n%s\n", usage);

  seqfile = argv[optind];

  if (outfile == NULL)
    outfp = stdout;
  else if ((outfp = fopen(outfile, "w")) == NULL)
    Die("Failed to open tRNA output file %s", outfile);

  
  if ((tRNA = (TRNA_TYPE *) malloc (sizeof(TRNA_TYPE))) == NULL)
    Die("Memory failure, couldn't allocate tRNA memory\n");
  if ((prev_tRNA = (TRNA_TYPE *) malloc (sizeof(TRNA_TYPE))) == NULL)
    Die("Memory failure, couldn't allocate tRNA memory\n");

    
  /***********************************************
   * Determine seq format & open for reading     *
   ***********************************************/

  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Can't determine format of file %s\n", seqfile);
  if ((seqfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);

  while (ReadSeq(seqfp, fmt, &seq, &sqinfo)) {

    if (ShowScores) 
      printf ("Seq: %s\n",sqinfo.name);
    
    tRNA_ct = 0;

    if (!NoReformat) {
      ToDNA(seq); 
      s2upper(seq);
    }

    /* allocate mem for integer-encoded seq (A=0,C=1,G=2,T=3) */
    if ((iseq = calloc (sqinfo.len+2, sizeof(char))) == NULL)
      Die("Memory failure, couldn't allocate sequence\n");
    
    /* integer-encode sequence */ 
    if (errno = IntEncodeSeq(iseq,seq,sqinfo.len))
      Die("Unable to encode sequence %s at base %d\n",
	  sqinfo.name,errno);
    
    /* Search both strands (0=top strand, -1=bottom strand)  */
    for (strand=0; strand >= -1; strand--) {

      Init_tRNA(prev_tRNA);     /* clear previous tRNA */
	
      /* take reverse complement of encoded seq if searching bottom strand */
      if (strand == -1) {
	if ((revseq = calloc (sqinfo.len+2, sizeof(char))) == NULL)
	  Die("Memory failure, couldn't allocate reverse sequence\n"); 
	revcomp(revseq, seq);
	free(seq);
	seq = revseq;
	if (IntEncodeSeq(iseq,seq,sqinfo.len))
	  Die("Unable to encode sequence\n");
      }
      
      
      /***********************************************
       * Find transcriptional promotor elements      *
       ***********************************************/
      
      seqidx=BBOX_START_IDX-2;
      while (GetBbox(&FirstScore,&seqidx,iseq,sqinfo.len,strand,
		     ShowScores)) {

	Init_tRNA(tRNA);
	tRNA->Bbox_st = seqidx;
	tRNA->Bbox_end = seqidx + BBOX_LEN-1;
	tRNA->BboxSc = FirstScore;

	if (((FirstScore >= SEC_LOBOUND) && 
	     (FirstScore <= SEC_HIBOUND))
	    && (GetSecABox(tRNA,seq))){

	  if (!GetBestTrxTerm(tRNA,seq,sqinfo.len,NoTermPenalty) &&
	      StrictMode) 
	    continue;       /* look for next B box */ 

	  strcpy(tRNA->acodon,"TCA");
	  /*  tRNA->Bbox_end++; */
	  tRNA->totSc = FirstScore;
	}
	
	else {           /* Searching for non-SelCys tRNA */
	  
	  
	  GetBestABox(tRNA,seq,iseq,sqinfo.len,strand,ShowScores,
		      Max_AB_dist,prev_tRNA->Abox_st);
	  IntScore = tRNA->AboxSc + tRNA->BboxSc + tRNA->ABdistSc;
	  if (IntScore < IntScoreCutoff) 
	    continue;               /* look for next B box */ 
	  
	  if (!RelaxedMode) {
	    if (!GetBestTrxTerm(tRNA,seq,sqinfo.len,NoTermPenalty) &&
		StrictMode)
	      continue;
	    TotScore = IntScore + tRNA->TermSc;
	    if (TotScore < TotScoreCutoff) 	      
	      continue;               /* look for next B box */ 

	    tRNA->totSc = TotScore;
	  }
	  else 
	    tRNA->totSc = IntScore;

	}
	
	Get_tRNA_stats(tRNA,seq,sqinfo.len,strand);
	
	if (tRNAOverlap(tRNA,prev_tRNA,strand)) {
	  
	  if  (tRNA->totSc < prev_tRNA->totSc) {
	    Init_tRNA(tRNA);  /* skip repeat tRNA */
	  }
	  else {           /* swap, but don't save yet */
	    
	    swap_tRNA = prev_tRNA;
	    prev_tRNA = tRNA;
	    Init_tRNA(swap_tRNA);
	    tRNA = swap_tRNA;
	  }      
	}
	else {             /* no overlap, save & then swap */

	  if (prev_tRNA->start > 0) {
	    prev_tRNA->idno = ++tRNA_ct;	  
	    Save_tRNA(prev_tRNA,&sqinfo,seq,strand,ShowScores,sqoffset);
	  }
	  swap_tRNA = prev_tRNA;
	  prev_tRNA = tRNA;
	  Init_tRNA(swap_tRNA);
	  tRNA = swap_tRNA;
	}

      }  /* find B box */
      
      /* save last buffered tRNA before going to other strand */

      if (prev_tRNA->start > 0) {
	prev_tRNA->idno = ++tRNA_ct;	  
	Save_tRNA(prev_tRNA,&sqinfo,seq,strand,ShowScores,sqoffset);
      }
      
    }   /* search both strands  */
        
    FreeSequence(seq, &sqinfo);
    free(iseq);
  }
  
  SeqfileClose(seqfp);
  fclose(outfp);
  return 0;
}

  








