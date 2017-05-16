/* eufindtRNA  - Eukaryotic tRNA finder
 *
 * pavesi.c - functions for finding transcriptional control regions 
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
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "squid.h"
#include "eufind_const.h"

/* #define NO_AMBIG       -use this option to eliminate conservative
 *                         calling of 'N's as best possible matches
 *                         in tRNAs -- useful for unfinished seqs with many N's
 */ 


/* log scores for each position in A Box */
/* six rows are for 1) A, 2) C, 3) G, 4) T, 5) (gap), 6) ambiguous base
   the ambiguous base value is the MIN (best score) of the ACGT rows  */

/* position 17a eliminated since always an empty pos (gap) */

float Abox_Mat[6][ABOX_LEN] = {
  {-1.268,-3.651,-0.899,-4.749,-5.442,-2.351,-3.363,-0.009,-1.977,-3.497,-5.442,
   -5.442,-5.442,-2.498,-4.749,-5.442,-0.031,-1.417,-1.180,-1.048,-4.344}, 

  {-3.651,-5.442,-4.056,-2.958,-0.480,-1.073,-0.857,-5.442,-5.442,-1.887,-2.498,
   -5.442,-5.442,-2.958,-2.224,-5.442,-5.442,-3.363,-1.417,-3.651,-0.393},

  {-0.779,-5.442,-0.598,-0.076,-3.651,-1.435,-1.614,-4.749,-0.154,-2.803,-5.442,
   0.000,0.000,-3.363,-3.651,-5.442,-3.497,-0.672,-1.012,-0.473,-3.651},

  {-1.453,-0.026,-3.651,-4.344,-1.036,-1.125,-1.073,-5.442,-5.442,-0.278,-1.399,
   -5.442,-5.442,-0.185,-0.827,-2.041,-5.442,-1.551,-2.447,-5.442,-1.253},

  {-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-0.412,
   -5.442,-5.442,-5.442,-0.868,-0.144,-5.442,-5.442,-5.442,-5.442,-5.442},
#ifdef NO_AMBIG
  {-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,
   -3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245}
#else
  {-0.779,-0.026,-0.598,-0.076,-0.480,-1.073,-0.857,-0.009,-0.154,-0.278,-1.399,
    0.000,0.000,-0.185,-0.827,-2.041,-0.031,-0.672,-1.012,-0.473,-0.393}
#endif
};

#define GAP_ROW 4             /* row in ABox mat with Gap weight */

float Bbox_Mat[6][BBOX_LEN] = {
  {-2.351,-5.442,-2.670,-5.442,-5.442,-1.472,0.000,-0.798,-2.498,-5.442,-3.497},
  {-3.245,-5.442,-5.442,-5.442,-0.004,-5.442,-5.442,-2.498,-1.435,-0.009,-0.190},
  {-0.175,-0.004,-5.442,-5.442,-5.442,-0.272,-5.442,-2.147,-5.442,-5.442,-3.651},
  {-3.651,-5.442,-0.072,0.000,-5.442,-4.749,-5.442,-1.048,-0.393,-5.442,-2.147},
  {-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442},
#ifdef NO_AMBIG
  {-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245,-3.245}
#else
  {-0.175,-0.004,-0.072,0.000,-0.004,-0.272,0.000,-0.798,-0.393,-0.009,-0.190}
#endif
};

#define ABDIST_MAT_SIZE 7
int ABDistIdx_Mat[ABDIST_MAT_SIZE] = {30,36,42,48,54,60,66};
float ABDistSc_Mat[ABDIST_MAT_SIZE] = {-0.46,-1.83,-2.35,-3.24,
				       -4.06,-3.83,-4.75};  

#define BTERM_MAT_SIZE 9
int BTermDistIdx_Mat[BTERM_MAT_SIZE] = {17,23,29,35,41,47,53,59,100};

float BTermDistSc_Mat[BTERM_MAT_SIZE] = {-0.54,-1.40,-2.80,-3.36,
					 -3.24,-5.44,-5.44,-4.06,-5.44}; 


int
Init_tRNA(TRNA_TYPE *tRNA) {
  
  strcpy(tRNA->iso_type,"???");
  strcpy(tRNA->acodon,"???");
  tRNA->start = tRNA->end = 0;
  tRNA->Abox_st = tRNA->Abox_end = tRNA->Abox_gap = 0; 
  tRNA->Bbox_st = tRNA->Bbox_end = tRNA->Term_st = tRNA->acodon_idx= 0; 
  tRNA->intron = tRNA->idno = 0;
  tRNA->totSc = tRNA->AboxSc = tRNA->BboxSc = -1000;
  tRNA->ABdistSc = tRNA->TermSc = -100;
}  

int
IntEncodeSeq (char *intseq, char *seq, int seqlen)
{
  int i;
  
  for (i=0; i<seqlen; i++) {
    switch (seq[i]) {
    case 'A': intseq[i]=0; break;
    case 'C': intseq[i]=1; break;
    case 'G': intseq[i]=2; break;
    case 'T': intseq[i]=3; break;
    default:               /* if not ACGT, assume ambiguous base */
      intseq[i]= 5;        /* ambig base will give score equal to best base
			      at that position in scoring matrix */
    }
  }
  return 0;
}


int
GetBbox (float *score, int *seqidx, char *iseq, int seqlen,
	 int strand, int verbose)
{
  int i,j, endidx;
  
  (*seqidx)++;
  endidx = seqlen - BBOX_LEN;
  for (i=*seqidx; i<endidx; i++) {
    *score = 0;
    for (j=0; j< BBOX_LEN; j++) {
      *score += Bbox_Mat[(int)(iseq[i+j])][j];
    }
    if (*score > BBOX_CUTOFF) {
      if (verbose)
	if (strand == 0) {
	  printf("Bbox at %i (End=%d), Sc= %.2f\n",i,i+BBOX_LEN+11,
		 *score);
	}
	else {
	  printf("Bbox at %i (End=%d), Sc= %.2f\n",seqlen-i+1,
		 seqlen-(i+BBOX_LEN+11)+1,*score);
	}
      *seqidx = i;
      return 1;
    }
  }
  return 0;
}

float 
Get_ABdist_weight(int ABdist)
{
  int ct;

  if (ABdist < MIN_AB_BOX_DIST)
    return MAX_PENALTY; 
  for (ct=0; ct < ABDIST_MAT_SIZE; ct++) {
    if (ABdist <= ABDistIdx_Mat[ct])
      return ABDistSc_Mat[ct];
  }
  return MAX_PENALTY;
}


int
GetSecABox(TRNA_TYPE *tRNA, char *seq)
{
  char *seqp;
  int i, startidx;

  /* Search for eukaryotic SelCys motif */

  startidx = tRNA->Bbox_st - ABOX_LEN - SEC_AB_BOX_DIST - 1;

  seqp = seq + MAX(0,startidx); 
  for (i=0; i<5; i++, seqp++) {
    if ((!strncmp(seqp,"GGTC",4) && 
	 (seqp[4] == 'T' || seqp[4] == 'C') &&
	 (seqp[5] == 'G') &&
	 (seqp[6] == 'T' || seqp[6] == 'G') &&
	 (seqp[7] == 'G') &&
	 (seqp[8] == 'G') &&
	 (seqp[9] == 'T'))) {
      
      tRNA->Abox_st = MAX(0,startidx+i-SEC_BBOX_DIST_CORR);
      strcpy(tRNA->iso_type,"SeCe");
      return 1;
    }
  }

  /* Search for Prokaryotic SelCys */

  startidx = tRNA->Bbox_st - 46 - 1;
  
  seqp = seq + MAX(0,startidx); 
  for (i=0; i<16; i++, seqp++) {
    if ((!strncmp(seqp,"GG",2) &&
	 (seqp[2] == 'A' || seqp[2] == 'T') &&
	 (seqp[3] == 'C' || seqp[3] == 'T') &&
	 (seqp[4] == 'T') &&
	 (seqp[5] == 'T') &&
	 (seqp[6] == 'C') &&
	 (seqp[7] == 'A') &&
	 (seqp[8] == 'A') &&
	 (seqp[9] == 'A') &&
	 (seqp[10] == 'A' || seqp[10] == 'T') &&
	 (seqp[11] == 'C') &&
	 (seqp[12] == 'C'))) {
      
      tRNA->Abox_st = MAX(0,startidx+i-23);
      strcpy(tRNA->iso_type,"SeCp");
      return 1;
    }
  }
  return 0;
  
}


int
GetBestABox (TRNA_TYPE *tRNA, char *seq, char *iseq, int seqlen,
	     int strand, int verbose, int Max_AB_dist, int prev_Abox_st)
{
  int i,           /* sequence position index */
    startidx, endidx,
    abox_end;
  int j,                /* matrix position index */
    offset1, offset2,   /* offset counters to keep track of gaps */
    gapct,              /* keeps track of 4 types of 2bp gaps in
			   positions 20a & 20b */
    best_gap,           /* gapct & offset1 vals for best score so far */ 
    best_offset1;

  float sc1, sc2, sc3,  /* components of total score */
    bestsc,             /* best score so far */
    abdistSc;


  startidx = MAX(MAX(0,(tRNA->Bbox_st - Max_AB_dist - ABOX_LEN)), prev_Abox_st+2);
  endidx =   MAX(0,(tRNA->Bbox_st - MIN_AB_BOX_DIST - ABOX_LEN +4));

  for (i=startidx; i < endidx; i++) {
    
    sc1=sc2=sc3=0;
    
    /* scoring Abox with weight matrix at tRNA pos 7-16 */
    for (j=0; j<=9; j++) {             
      sc1 += Abox_Mat[(int)(iseq[i+j])][j];
    }

    j=10;

    /* score gap at pos 17 by looking for conserved 'GG' at pos 18 &
       19 */

    if ((seq[i+j] == 'G')) {   /*  && (seq[i+j+1] == 'G')) { */
      sc2 = Abox_Mat[GAP_ROW][j];
      offset1 = 1;
    }
    else {
      sc2 = 0;
      offset1 = 0;
    }
    
    /* scoring Abox with weight matrix at tRNA pos 18-20 */
    for (j=10; (j+offset1) < 14; j++) {
      sc2 += Abox_Mat[(int)(iseq[i+j])][j+offset1];
    }
    
    /* score potential gap at 20a & 20b, plus rest of matrix up to
       position;  gapct, enumerates all possible 2bp gaps */

    for (gapct=0; gapct<4; gapct++) {
      j=14-offset1;
      offset2=0;
      switch (gapct) {
      case 0:
	sc3  = Abox_Mat[(int)(iseq[i+j])][j+offset1+offset2];
	j++;
	sc3 += Abox_Mat[(int)(iseq[i+j])][j+offset1+offset2];
	j++;
	break;
      case 1:
	sc3  = Abox_Mat[GAP_ROW][j+offset1+offset2];
	offset2++; 	
	sc3 += Abox_Mat[(int)(iseq[i+j])][j+offset1+offset2];
	j++;
	break;
      case 2:
	sc3 = Abox_Mat[(int)(iseq[i+j])][j+offset1+offset2];
	j++;
	sc3 += Abox_Mat[GAP_ROW][j+offset1+offset2];
	offset2++;
	break;
      case 3:
	sc3 = Abox_Mat[GAP_ROW][j+offset1+offset2];
	offset2++;
	sc3 += Abox_Mat[GAP_ROW][j+offset1+offset2];
	offset2++;
	break;
      }
      for (; (j+offset1+offset2) < ABOX_LEN; j++) {
	sc3+= Abox_Mat[(int)(iseq[i+j])][j+offset1+offset2];
      }
      
      abox_end = i+ABOX_LEN-offset1-offset2-1;
      abdistSc = Get_ABdist_weight(tRNA->Bbox_st-abox_end-1);
      
      if ((sc1 + sc2 + sc3 + abdistSc) > (tRNA->AboxSc + tRNA->ABdistSc)) {
	tRNA->Abox_st = i;
	tRNA->Abox_end = abox_end;
	tRNA->ABdistSc = abdistSc;
	tRNA->AboxSc = sc1+sc2+sc3;
	best_offset1= offset1;
	best_gap = gapct;
	if (verbose)
	  if (strand == 0) {
	  printf("Abox at %d (St=%d) A:%.2f AB(%d):%.2f I:%.2f\n",
		 i,i-5,tRNA->AboxSc, tRNA->Bbox_st-abox_end-1,tRNA->ABdistSc,
		 tRNA->AboxSc+tRNA->BboxSc+tRNA->ABdistSc);
	  }
	  else {
	    printf("Abox at %d (St=%d) A:%.2f AB(%d):%.2f I:%.2f\n",
		   seqlen-i+1,seqlen-(i-5)+1,tRNA->AboxSc, 
		   tRNA->Bbox_st-abox_end-1,tRNA->ABdistSc,
		   tRNA->AboxSc+tRNA->BboxSc+tRNA->ABdistSc);
	  }
	    
      }
    }   /* for gapct, enumerating all possible gaps */
  }   /* for i, starting pos for A box */
}


int
GetBestTrxTerm (TRNA_TYPE *tRNA, char *seq, int seqlen,
		float TermPenalty)
{
  int i,     /* current seq position */
    startidx, endidx;    /* start & end points for term search */
  int ct, BTermdist;

  float score;               /* current score */

  startidx = tRNA->Bbox_end+MIN_BTERM_DIST-1;
  endidx = MIN(startidx+MAX_TERM_SEARCH,(seqlen-4));

  for (i=startidx; i<endidx; i++) {
    
    if ((seq[i] == 'T') && (seq[i+1] == 'T') && 
	(seq[i+2] == 'T') && (seq[i+3] == 'T')) {
      
      BTermdist = i-tRNA->Bbox_end-1;
      for (ct=0; BTermDistIdx_Mat[ct] < BTermdist; ct++) {
      }
      tRNA->TermSc = BTermDistSc_Mat[ct];
      tRNA->Term_st = i;
      return 1;
    }
  }
  
  tRNA->Term_st = -1;
  
  if (endidx ==  (seqlen-4)) {
    tRNA->TermSc =  TOT_SCORE_THRESH - (INT_SCORE_THRESH);
    return 1;
  }
  else { 
    tRNA->TermSc = MAX_PENALTY;
    return 0;
  }
}


/* Uses tranlation scheme & AA lookup table from Squid library */

void
Get_IsoType (TRNA_TYPE *tRNA)
{
  int i, codon;
  char codon_seq[4];
  
  revcomp(codon_seq,tRNA->acodon);
  
  codon = 0;
  for (i = 0; i < 3; i++)
    {
      codon *= 4;
      switch (codon_seq[i]) {
      case 'A': case 'a':             break;
      case 'C': case 'c': codon += 1; break;
      case 'G': case 'g': codon += 2; break;
      case 'T': case 't': codon += 3; break;
      case 'U': case 'u': codon += 3; break;
      default: codon = 64; break;
      }
      if (codon == 64) break;
    }
  strcpy(tRNA->iso_type,stdcode3[codon]);
}

void
Get_anticodon (TRNA_TYPE *tRNA, char *seq)
{
  char *acodonp;
  int startidx, i, besti, score, bestsc;

  startidx =  tRNA->Abox_end + 7;
  acodonp = seq + startidx;
  bestsc = besti = 0;
  for (i=0; i<7; i++) {
    score = 0;
    if ((acodonp[i-2] == 'C') || (acodonp[i-2] == 'T')) 
      score++;
    if (acodonp[i-1] == 'T')  
      score++;
    if ((acodonp[i+3] == 'A') || (acodonp[i+3] == 'G')) 
      score++;
    if (score > bestsc) {
      bestsc = score;
      besti = i;
    }
  }
  strncpy(tRNA->acodon,acodonp+besti,3);
  /*  tRNA->acodon[3]='\0'; */
  tRNA->acodon_idx = startidx + besti;
}


void
Get_tRNA_stats (TRNA_TYPE *tRNA, char *seq, int seqlen, int strand)
{

  tRNA->start = MAX(1,(tRNA->Abox_st - 5));
  tRNA->end = MIN((tRNA->Bbox_end + 12),seqlen);
    
  if (strand == -1) {
    tRNA->start = seqlen - tRNA->start + 1;
    tRNA->end = seqlen - tRNA->end + 1;
  }
  
  if (!strncmp("SeC",tRNA->iso_type,3)) {
    strcpy(tRNA->acodon,"TCA");
    tRNA->Bbox_end++;
  }
  else {
    Get_anticodon(tRNA,seq);
    Get_IsoType(tRNA);
  }
}

void
Save_tRNA  (TRNA_TYPE *tRNA, SQINFO *sqinfo, char *seq, int strand,
	    int ShowScores, long int sqoffset)
{

  if (ShowScores) 
    printf("%s.%d\t%ld\t%ld\t%s\t%d\tA:%.2f B:%.2f AB:%.2f T:%.2f Tot:%.2f\n",
	   sqinfo->name,tRNA->idno,tRNA->start+sqoffset,tRNA->end+sqoffset,
	   tRNA->acodon,strand, tRNA->AboxSc,tRNA->BboxSc,tRNA->ABdistSc,
	   tRNA->TermSc,tRNA->totSc);

  else
    printf("%-10s\t%d\t%ld\t%ld\t%s\t%s\t0\t0\t%.2f\n",
	   sqinfo->name,tRNA->idno,tRNA->start+sqoffset,tRNA->end+sqoffset,
	   tRNA->iso_type,tRNA->acodon,tRNA->totSc);
    
}

int
tRNAOverlap (TRNA_TYPE *tRNA1, TRNA_TYPE *tRNA2, int strand)
{
  if (strand == 0) {
    if ((((tRNA1->start >= tRNA2->start) && 
	  (tRNA1->start <  tRNA2->end-MAX_OVLAP)) ||
	 ((tRNA1->end > tRNA2->start+MAX_OVLAP) && 
	  (tRNA1->end <= tRNA2->end))) ||

	(((tRNA2->start >= tRNA1->start) && 
	  (tRNA2->start < tRNA1->end-MAX_OVLAP)) ||
	 ((tRNA2->end > tRNA1->start+MAX_OVLAP) && 
	  (tRNA2->end <= tRNA1->end))))
      return 1;
    else
      return 0;
  }
  else {
    if ((((tRNA1->start <= tRNA2->start) && 
	  (tRNA1->start > tRNA2->end+MAX_OVLAP)) ||
	 ((tRNA1->end < tRNA2->start-MAX_OVLAP) && 
	  (tRNA1->end >= tRNA2->end))) ||

	(((tRNA2->start <= tRNA1->start) && 
	  (tRNA2->start > tRNA1->end+MAX_OVLAP)) ||
	 ((tRNA2->end < tRNA1->start-MAX_OVLAP) && 
	  (tRNA2->end >= tRNA1->end))))
      return 1;
    else
      return 0;
  }
}



