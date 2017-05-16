/*
 * Copyright, 1991, The Regents of the University of California.
 * This software was produced by the Los Alamos National Laboratory,
 * which is operated by the University of California for the United
 * States Department of Energy under contract W-7405-ENG-36.  The
 * U. S. Government is licensed to use, reproduce, and distribute
 * this software.  Permission is granted to the public to copy and
 * use this software without charge, provided that this Notice and
 * any statement of authorship are reproduced on all copies.
 * Neither the Government nor the University makes any warranty,
 * express or implied, or assumes any liability or responsibility
 * for the use of this software.
*/

#define PROGRAM  "trnascan"
#define RELEASE "1.4 (Feb 96)"
#define DERIV_VERSION  "1.3 (Oct 91)"
#define BRIEF    "'Identification of tRNA genes in genomic DNA'"
#define CITATION "Fichant and Burks, J. Mol. Biol. (1991) 220:659-671."
#define MODIF "(modified & optimized for use in tRNAscan-SE package by T. Lowe  2/96)"

/* Modified by T. Lowe  11/95  */
/* Changes:  1) Search parameters named in #define constants
             2) Print statements added to help trace progress
	         of search - VERBOSE constant must be defined
		Trace statements sent to "tscan.verb.out"
	     3) Bug that caused program to crash on any non-ACGT
	        sequence characters fixed 
	     4) fgetseq() modified to correctly read in fasta sequences
             5) compstrand() modified to increase efficiency & 
	        getseqsize() added to allow input sequences of any 
                length (memory allowing)
	     6) Numerous calls to strlen(sequence) eliminated for
	        efficiency
	     7) Calls to myindex() function eliminated - replaced
	        by in-line switch statements for efficiency
	     8) basepairings() function rewritten for efficiency
	     9) program ANSI-fied to allow compilation with gcc
            10) Fixed bug: indexing out of 'ntab' array bounds 
	        causing unpredictable side-effects
                
  
    **  Modifications result in over 200-fold speed increase **
 
    -T. Lowe  9/2000
     added "-i" option to allow alternate start of sequence nuc numbering

*/

/* #define NO_AMBIG       -use this option to eliminate conservative
 *                         calling of 'N's as base pairing matches
 *                         in tRNAs (this gives more false positives)
 */ 

static char banner[] = "trnascan: scan a sequence for tRNAs"; 

char usage[] = "\n\
Usage: trnascan [-options] <seqfile>\n\
where supported options are:\n\
-s             : use original tRNAscan 1.3 parameters (default)\n\
-r             : use relaxed search parameters (used with tRNAscan-SE)\n\
-a             : use alternate (user-set) search parameters\n\
-c             : suppress credits\n\
-o <outfile>   : write results to <outfile>\n\
-i <integer>   : start sequence numbering at <integer> (def=1)\n\
-h             : print (this) short help message\n\n";

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define TRUE 1
#define FALSE 0

#ifndef TSCANDIR
#define TSCANDIR "/usr/local/lib/trnascan"
#endif  


#define MAXLINE 1000    /* max input seq line length */
#define BUF_SIZE 100   /* extra room added onto allocated sequence */
		       /* over that determined by getseqsize()  */

#define MAX(x,y) (((x) > (y)) ? (x) : (y))

/* #define VERBOSE   */

                                  /* Original version 1.3 default parameters */
#define ST_SG_CUTOFF 5            /* general score (SG) cutoff */
#define ST_TPC_SIG_THRESH 0.40    /* TPC signal sequence matrix score */
			       /*   cutoff */
#define ST_D_SIG_THRESH 0.40      /* D signal sequence matrix score cutoff */    
#define ST_TPC_INV 2        /* Number of TPC matrix invariant bases */
			       /*   allowed NOT to be invariant */
#define ST_TPC_INCSG 5         /* Number of base pairs required in TPC */
			       /*   stem to increment general score */ 
#define ST_TPC_KEEP 4          /* Number of base pairs required in TPC */
			       /*   stem to keep trna as a candidate */
#define ST_D_INV 1          /* Number of D matrix invariant bases */
			       /*   allowed NOT to be invariant */
#define ST_LOOK_FOR_ACLOOP_SG 4   /* Minimum SG required to begin looking */
			       /*   for anticodon loop */
#define ST_ACLOOP_MIN 4            /* Minimum base pairs required in */
			       /*   anticodon loop */
#define ST_AA_INCSG 7          /* Number of base pairs in Amino acyl */
			       /*   stem needed to increment SG */ 
#define ST_AA_KEEP 6           /* Number of base pairs in Amino acyl */
			       /*   stem needed to keep tRNA candidate */

                        /* "Relaxed" parameters - used by default with */
			/* tRNAscan-SE program.  Makes tRNAscan into */
			/* rough pre-filter for Covariance tRNA prediction */
			/* program by S. Eddy  */
#define RX_SG_CUTOFF 5         
#define RX_TPC_SIG_THRESH 0.40 
#define RX_D_SIG_THRESH 0.30   
#define RX_TPC_INV 2     
#define RX_TPC_INCSG 4      
#define RX_TPC_KEEP 2        
#define RX_D_INV 2        
#define RX_LOOK_FOR_ACLOOP_SG 3 
#define RX_ACLOOP_MIN 3          
#define RX_AA_INCSG 5        
#define RX_AA_KEEP 4         

                         /* "Alternate" params - for experimenting */
			 /*    with other param values */
#define ALT_SG_CUTOFF 4         
#define ALT_TPC_SIG_THRESH 0.40 
#define ALT_D_SIG_THRESH 0.30   
#define ALT_TPC_INV 2     
#define ALT_TPC_INCSG 4      
#define ALT_TPC_KEEP 2       
#define ALT_D_INV 2        
#define ALT_LOOK_FOR_ACLOOP_SG 3 
#define ALT_ACLOOP_MIN 3          
#define ALT_AA_INCSG 6        
#define ALT_AA_KEEP 4         

#define MIN_VAR_LOOP 28     /* Minimum variable loop size, assumes min */
			    /*   intron length = 8bp */
#define MAX_INTRON_LEN 60   /* Maximum allowable intron length */
#define MIN_SEQ_LEN 70      /* Minimum length of sequence that will be */
			    /*  for tRNAs */

#define STRICT_PARAMS 1;
#define RELAXED_PARAMS 2;
#define ALT_PARAMS 3;


typedef struct pset {
  float tpc_sig_thresh, d_sig_thresh;
  int sg_cutoff, 
  tpc_inv, tpc_incsg, tpc_keep,
  d_inv,
  look_for_acloop_sg, acloop_min,  
  aa_incsg, aa_keep;
} Param_set_type;
  
Param_set_type ps;


void 
set_search_params (Param_set_type *ps,
		   int params);

/* Subroutine that accomplishes the end of the test for the presence of a
   tRNA gene */

void 
following_search(long int pos,    /* first position of the found T-Psi-C signal */ 
		 long int pos1,   /* first position of the found D signal */
		 char *ptr1,      /* pointer to the first position of the D signal */
		 char *ptr3,      /* ptr3=ptr1+2  pointer to first position of the D arm  */
		 int lpair, 
		 int nloop,       /*nloop=0 test of the sequence, */
				  /*nloop=1 test of the complementary sequence */
		 char *sequence,  
		 long int seqlen, 
		 FILE *fpo,       /*pointer to the output file */
		 FILE *fpverb, 
		 char *name, 
		 int score,       /* Value of the general score SG */
		 int match2,     /* integer testing the presence of D */
				 /* arm with 3 base-pairings */ 
		 int *ntrna, 
		 int *npred, 
		 int *match,
		 long int sqoffset   /* offset nucleotide numbering by this much (set with -i param) */
		 );


/* Subroutine to read the consensus matrix */

void 
lectval(FILE *fp,                  /* pointer to the consensus matrix file */
	float (*table_cons)[4],    /* table containing the frequency of each base
				      at each position of the signal */
	int (*table_inv)[2],  /* table containing the position and the nature of
				 the invariant bases found in the signal. Code for
				 the bases: A=0, C=1, G=2 and T=3. */
	int *lsig,            /*lsig=length of the signal */
	int *ktot,            /* ktot= number of invariant bases */
	float *maxtot         /* maxtot= sum of the maximum frequencies */
	);

/* Subroutine reading the sequence */
/* Modified to correctly read FASTA sequence files */

int fgetseq(char *name,           /* string w/name of the sequence */
	    char *sequence,       /* character string containing the sequence */
	    long int *seqlen,     /* length of sequence */
	    FILE *fpi);            /* input file pointer */
		

/* Subroutine reading & returning the sequence length */

int getseqsize(FILE *fpi);       /* input file pointer */


/* Subroutine looking for the presence of a given signal, returns 1 if a
   signal is found and 0 otherwise. It also return the table 'weight'
   containing the frequencies of the oberved bases in the windowed
   sequence and the number 'ninv' of invariant bases found in the windowed
   sequence*/

int readsignal(char *ptr,      /* pointer to the sequence */
	       int (*table_inv)[2],  /* table containing the position and nature of the
					invariant bases found in the consensus matrix */
	       int *lsig,   /* lsig= length of the signal */
	       int *ktot,   /* ktot= number of invariant bases in the consensus matrix */
	       float *weight,       /* table containing the frequencies of the observed
				       base at each position of the
				       windowed sequence tested */
	       float (*table_cons)[4],  /* table corresponding to the consensus matrix */
	       int *ninv,               /* ninv= number of invariant */
					/* bases in the windowed  sequence */
	       int threshold_inv);       /* Number of invariant bases */
					/* allowed not to b e invariant */


/* Subroutine that calculates the similarity score on the potential signal 
   previously retained by the subroutine readsignal. This subroutine
   returns 1 it the computed score is greater or equal to the defined 
   threshold and 0 otherwise. It returns also the value of the computed 
   score (score) */ 

int scoring(float *weight, /* table containing the frequencies of the observed base 
			      at each position of the potential signal */
	    int lsig,      /* length of the signal */
	    float max,     /* sum of the maximum frequencies found in */
			   /* the consensus matrix */ 
	    int ktot,      /* number of invariant bases found in the consensus matrix */
	    float *score,  /* value of the computed score on the potential signal */
	    float ThresholdValue,   /* defined threshold for the similarity score */
	    int ninv   /* number of invariant bases found in the potential signal */
	    );

/* Subroutine looking for base-pairings between two parts of the sequence.
   It returns the number of base-pairings found (ncomp) */

void 
basepairing(char *ptr,   /* pointer to the sequence */
	    int npair,   /* number of base-pairings forming a given arm */
	    int lpair,   /* number of nucleotides found between the first position of the
			    first part of the sequence involved in the stem and the last 
			    position of the second part of the
			    sequence involved in the stem */
	    int *ncomp   /* number of base-pairings observed between the two parts of
			    the sequence tested */
	    );
           

/* Subroutine that complements the sequence. It returns the complementary 
   sequence to the main */

void 
compstrand(char **sequence,   /* pointer to the sequence string */
	   long int seqlen    /* sequence length   */
	   );


/* Subroutine that codes the anticodon signal sequence by a number
   comprised between 1 and 65. It returns this number. */

void 
codage(char *anticodon,    /* anticodon signal sequence */
       int length1,        /* length of the anticodon signal, lenght1=3 */
       int *num            /* number associated to the anticodon signal sequence */
       );            
         

/* Subroutine that determines the tRNA gene family */

void 
corresaa(int num,           /* Number coding the anticodon signal sequence */
	 char *type_trna    /* tRNA gene family */
	 );
                 

/* Subroutine that prints the results of the search */

void 
printresult(FILE *fpo,     /* output file pointer */
	    FILE *fpverb,   /* character string for the name of the sequence */
	    char *name, 
	    long int pos1,    /* first position of the D signal */
	    long int pos,     /* first position of the T-Psi-C signal */
	    int lpair,        /* number of nucleotides between the */
			      /* first position  of the D arm and the last one */
	    int lpair1,     /* number of nucleotides between the first position of the
			       aminoacyl arm and the last one */
	    int lpair2,     /* number of nucleotides between the first position of the 
			       anticodon arm and the last one */
	    int nloop,      /* nloop=0, scanning of the direct strand; nloop=1 scanning of 
			       the complementary strand */
	    int *ntrna,     /* number of tRNA genes predicted in the sequence */
	    char *chaine2,  /* character string for the predicted tRNA gene sequence */
	    char *sequence,    /* character string containing the sequence tested */
	    long int length,   /* length of the sequence */
	    int *match,      /* match=1 if at least one tRNA gene has been found on the 
				direct strand and 0 otherwise */
	    int ncomp,       /* number of base-pairings in the anticodon arm of the
				predicted tRNA gene */
	    char *type_trna,     /* character string for the tRNA gene family */
	    char *anticodon,      /* character string for the anticodon signal sequence */
	    long int sqoffset   /* offset nucleotide numbering by this much (set with -i param) */
	    );

main(int argc, char **argv)
{
  /* pointers to the different files fpi=input file, fpo=output file, 
   fpcons1= T-Psi-C matrix file, fpcons2= D matrix file */

  extern Param_set_type ps;   /* search parameters */

  FILE *fpi,*fpo,*fpcons1,*fpcons2, *fpverb;

  /* lsig1= length of the T-Psi-C signal, lsig2= length of the D signal
     ktot1= number of invariant bases in the T-Psi-C signal
     ktot2= number of invariant bases in the D signal */

  int lsig1=0, lsig2=0, ktot1=0, ktot2=0;

  long int pos,pos1, begin; /* start position of the T-Psi-C signal (pos)
			       and D signal (pos1) and of the search for the
			       D signal (begin) */
  int npair=0, lpair=0; /* variables used to test the presence of a stem, 
			   see definition further in the program */
  float table_cons1[30][4]; /* table containing the frequency of each base 
			       at each position of the T-Psi-C signal */
  float table_cons2[30][4]; /* table containing the frequency of each base
			       at each position of the D signal */
  float maxtot1=0, maxtot2=0; /* sum of the maximum frequencies found in 
				 the T-Psi-C matrix (maxtot1) and in the D matrix
				 (maxtot2) */
  float weight1[30]; /* table containing the frequency of the observed base 
			at each position of the windowed sequence tested
			(T-Psi-C matrix) */
  float weight2[30]; /* table containing the frequency of the observed base
			at each position of the windowed sequence tested
			(D matrix) */
  int table_inv1[30][2]; /* table containing the position and the nature of 
			    the invariant bases found in the T-Psi-C matrix */
  int table_inv2[30][2]; /* table containing the position and the nature of
			    the invariant bases found in the D matrix */
  int ntab[5][2]; /* table containing for each potential D arm found in the
		     windowed sequence, the number of base-pairings and the 
		     number of nucleotides that separates the first position
		     of the D arm and the last one */
  char *sequence; /* character string containing the sequence */
  long int seqlen; /* keeps the length of the current sequence */
  long int sqoffset=0;  /* start numbering nucleotides by this offset (def=0) */
  char  *ptr,*ptr1,*ptr3, *ptrstart; /* pointers to the sequence, their
					definitions are given at the place 
					they are used in the program */	 
  char name[80]; /* character string containing the name of the sequence */
  char tscan_dir[120];  /* holds name of directory of consensus files */
  
  int ntrna; /* number of tRNA genes predicted in the sequence and its
		complement */
  int npred=0; /* Total number of tRNA genes predicted in the sequences of 
		  the input file */
  int nseq=0; /* Number of sequences tested (number of sequences of the
		 input file) */
  long int lseq=0; /* Number of nucleotides tested */
  int score; /* variable corresponding to the general score SG */
  int ninv=0; /* number of invariant bases oberved in the windowed sequence
		 for each signal */
  int ncomp=0; /* number of base-pairings observed in the sequence tested 
		  for the presence of one of the four arms */
  int threshold_inv; /* number of invariant bases that are allowed not to
			be invariant */
  int match; /* if at least one tRNA gene is predicted of the direct strand,
		match=1, if no tRNA gene is found on the direct strand, match=0 */
  int first_score=0,match2,h, nloop;
  
  /* ThresholdValue1= threshold value of the similarity score used to retain 
     a T-Psi-C signal
     ThresholdValue2= threshold value of the similarity score used to retain 
     a D signal 
     Changing the threshold value of the similarity score, MODIFY '0.4' in
     the following line. */
  
  float ThresholdValue1, ThresholdValue2;
  float score1=0, score2=0;
  int i,j;
  
  int params;
  int suppress_credits;      /* flag for display of credits */
  char *seqfile,   /* name of input sequence file  */
    *outfile,   /* name of output file (if not sent to stdout) */
    *verbfile;  /* name of file for verbose output */
  
  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */
  
/*********************************************** 
 * Parse command line
 ***********************************************/
  
  suppress_credits = FALSE;
  params        = STRICT_PARAMS;
  outfile          = NULL;
  verbfile         = NULL;
  
while ((optc = getopt(argc, argv, "csraho:v:i:")) != -1)
  switch (optc) {
    
  case 'c': suppress_credits = TRUE; break;
  case 's': params = STRICT_PARAMS;  break;
  case 'r': params = RELAXED_PARAMS;  break;
  case 'a': params = ALT_PARAMS;  break;
  case 'o': outfile   = optarg;  break;
  case 'v': verbfile  = optarg;  break;
  case 'i': sqoffset = atoi(optarg)-1; break;

  case 'h': 
    printf("%s\n  version %s\n%s\n", banner, RELEASE, usage);
    exit(0);
  default:
    fprintf(stderr,"unrecognized option %c\n", optc);
    exit(1);
  }

if (argc - optind != 1) {
  fprintf(stderr,"%s\n", usage);
  exit(1);
}
  
seqfile = argv[argc-1];


/* open the sequence file */

if ((fpi=fopen(seqfile,"r")) == NULL)
  {
  fprintf(stderr,"tRNAscan1.4: FATAL:  Cannot open the input sequence file %s\n",seqfile);
  exit(1);
  }

/* open the output file */

if (outfile == NULL) {
  fpo = stdout;
}
else if ((fpo=fopen(outfile,"w")) == NULL) {
  fprintf(stderr,"tRNAscan1.4: FATAL: Cannot open the output file %s\n",outfile);
  exit(1);
}
 
/* open Verbose output file */
if (verbfile != NULL) 
  if ((fpverb=fopen(verbfile,"w")) == NULL) {
    fprintf(stderr,"tRNAscan1.4: Cannot open verbose output file %s\n",verbfile);
    exit(1);
}

/* open the consensus matrix file for the T-Psi-C signal */

strcpy(tscan_dir,TSCANDIR);

if (((fpcons1=fopen("TPCsignal","r")) == NULL) &&
    ((fpcons1=fopen(strcat(tscan_dir,"/TPCsignal"),"r")) == NULL))
  {
  fprintf(stderr,"tRNAscan1.4: main cannot open %sTPCsignal consensus file\n",tscan_dir);
  exit(1);
  }

/* open the consensus matrix file for the D signal */

strcpy(tscan_dir,TSCANDIR);

if (((fpcons2=fopen("Dsignal","r")) == NULL) &&
    ((fpcons2=fopen(strcat(tscan_dir,"/Dsignal"),"r")) == NULL))
  {
  fprintf(stderr,"tRNAscan1.4: main cannot open Dsignal consensus file\n");
  exit(1);
  }


/* Set search parameters */

set_search_params(&ps,params); 

ThresholdValue1 = ps.tpc_sig_thresh;
ThresholdValue2 = ps.d_sig_thresh;

/* Credits   */

if (!suppress_credits) {
  printf("\n %s %s -- derived from version %s\n", PROGRAM, RELEASE, 
	 DERIV_VERSION);
  printf("\n Please cite: %s\n %s\n", BRIEF,CITATION); 
  printf(" %s\n\n",MODIF);
}

/* reading the two consensus file */

lectval(fpcons1,table_cons1,table_inv1,&lsig1,&ktot1,&maxtot1);
lectval(fpcons2,table_cons2,table_inv2,&lsig2,&ktot2,&maxtot2);


/* main loop for reading & analyzing one sequence at a time */

while (!feof(fpi)) {

/* find out sequence size before allocating memory & reading in */

  seqlen = getseqsize(fpi);
  sequence= (char *) calloc (seqlen+BUF_SIZE,sizeof(char));
  if (sequence == NULL ) { 
    fprintf(stderr,"tRNAscan1.4: Can't malloc for sequence\n");
    exit(-1);
  }

/* reading the name of the sequence and the sequence from the input file */	

  if (fgetseq(name,sequence,&seqlen,fpi) == 0)
    break;

  nseq++;
  lseq += seqlen;

  match=0;
  for (nloop=0; nloop <=1; nloop++)
    {
/* if the sequence is shorter than 70 bp long it is skipped */

    if (seqlen < MIN_SEQ_LEN)
    break;

    ntrna=0;

/*  search of tRNA genes starts 
    pos= first position of the T-Psi-C signal
    ptr= pointer to the first position of the T-Psi-C signal */

    for (pos=44, ptr=(sequence+43); pos<seqlen-23;pos++, ptr++)
      {
      score=0;

#ifdef VERBOSE
	  fprintf(fpverb,"--- Moving TPC window (pos=%d)\n",pos);
#endif 
	  

/*change in the number of invariant bases that must be present in the
  T-Psi-C signal, MODIFY the value of threshold_inv here */

      threshold_inv= ps.tpc_inv;

/*search for the T-Psi-C signal */

      if(readsignal(ptr,table_inv1,&lsig1,&ktot1,weight1,table_cons1,&ninv,
        threshold_inv))
        {

/*the general score SG (score) is incremented by 1, if the number of 
  invariant bases (ninv) in the found T-Psi-C signal >= 3 (ktot1=4)
  change in the threshold of the T-Psi-C signal for the increment of SG,
  MODIFY the following line */

        if (ninv >= ktot1-1)  {
	  score++;
#ifdef VERBOSE
	  fprintf(fpverb,"+ TPC invariant bp: %d.  SG++\n",ninv); }
	else  { 
	  fprintf(fpverb,"X TPC invariant bp: %d.  NO SG++\n",ninv); 	  
#endif 
	}

/* Computation of the score on the found T-Psi-C signal, signal retains if
   the computed score is >= ThresholdValue1 */

        if(scoring(weight1,lsig1,maxtot1,ktot1,&score1,ThresholdValue1,
		   ninv))
          {
#ifdef VERBOSE
	   fprintf(fpverb,"+ TPC signal over thresh: %f\n",score1);
#endif 


/* npair= number of base-pairings in the T-Psi-C stem 
   lpair= number of nucleotides between the first position of the T-Psi-C
          arm and the last position of the T-Psi-C arm */

          npair = 5;
          lpair = 16;

/* Test for the presence of the T-Psi-C stem */

          basepairing(ptr+1,npair,lpair,&ncomp);

/* If the number of base-pairings (ncomp) equal 5, SG (score) incremented 
   by 1. Change in the threshold of base-pairings in the T-Psi-C arm for
   the increment of SG, MODIFY the value '5' in the following line */

          if (ncomp >= ps.tpc_incsg) {
            score++;
#ifdef VERBOSE
 	    fprintf(fpverb,"+ TPC parings (%d) add to SG=%d\n",ncomp,score); }
	  else  {
 	    fprintf(fpverb,"X TPC parings (%d) NO add to SG=%d\n",ncomp,score); 
#endif 
	  }

/* Change in the threshold to retain a stem as a potential T-Psi-C arm,
   MODIFY the value '4' in the following line */

#ifdef VERBOSE
          if (ncomp < ps.tpc_keep)  {
	    fprintf(fpverb,"X TPC parings retain loop: %d\n",ncomp); }
	  else {
	    fprintf(fpverb,"+ TPC parings retain loop: %d\n",ncomp); 
#else
          if (ncomp >= ps.tpc_keep)  {
#endif 


/* For the same found T-Psi-C region, different potential D region can be
   found because of the possible different lengthes of the D loop. As the
   search for the tRNA gene is applied for each potential D arm, the value
   of SG computed on the T-Psi-C region is stored in 'first_score' */

            first_score=score;

/*Change in the number of invariant bases to be present in the D signal,
  MODIFY the value of threshold_inv here */

            threshold_inv= ps.d_inv;

/* Search for the presence of a D region between -120 and -37 nucleotides 
   upstream from the found T-Psi-C region. 37 nucleotides is the observed
   minimum length assuming no insertion in the D and variable loop. 120
   nucleotides allows for an intron of up to 60 nucleotides.
   Change in the length of the intron, MODIFY the value '120' in the
   definition of begin and ptrstart.
   ptr1= pointer to the first position of the D signal
   pos1= first position of the D signal 
   begin= starting position for the search of the D signal
   ptrsart= pointer to the starting position for the search of the D signal */

/* If the sequence does not have 127 nucleotides upstream of the T-Psi-C
   signal, the search for the D signal starts at position 8 */

            if (pos <=127)
              {
               begin = 8;
               ptrstart = sequence+7;
              }
            else
              {
               begin = pos-(MAX_INTRON_LEN+60);
               ptrstart = ptr-(MAX_INTRON_LEN+60);
               }

            for(pos1=begin,ptr1=ptrstart; pos1<=(pos-37);pos1++,
                ptr1++)
              {

/* Search for the D signal */

#ifdef VERBOSE
              if(!readsignal(ptr1,table_inv2,&lsig2,&ktot2,weight2,
                 table_cons2,&ninv,threshold_inv))
                {
		  fprintf(fpverb,"X D invariant bp: %d\n",ninv); }
	      else 
		{
		  fprintf(fpverb,"+ D invariant bp: %d\n",ninv); 	     
#else
              if(readsignal(ptr1,table_inv2,&lsig2,&ktot2,weight2,
                 table_cons2,&ninv,threshold_inv))
                {
#endif 


/* Number of invariant bases found in the windowed sequence equal to the
   number of invariant bases of the matrix (ktot2=3), then SG (score) is
   incremented by 1 
   Change in the increment of SG for the D signal, MODIFY the following 
   line */

		  if (ninv >= ktot2) {
		    score++;
#ifdef VERBOSE
		    fprintf(fpverb,"+ D invariant (%d) inc SG\n",ninv); }
		  else {
		    fprintf(fpverb,"X D invariant (%d) NO add to SG\n",ninv);
#endif 
		  }

/* Computation of the score on the found D signal, signal retains if the
   computed score is >= TresholdValue2 */

#ifdef VERBOSE
		if(!scoring(weight2,lsig2,maxtot2,ktot2,&score2,
                   ThresholdValue2,ninv))
                  {
		    fprintf(fpverb,"X D signal threshold: %f\n",score2); }
		else  {  
		  fprintf(fpverb,"+ D signal threshold: %f\n",score2); 
#else
		if(scoring(weight2,lsig2,maxtot2,ktot2,&score2,
                   ThresholdValue2,ninv))
                  {
#endif 


/* ptr3= pointer pointing on the first position of the D arm 
   npair= number of base-pairings in the D stem
   lpair= number of nucleotides between the first position of the D arm
   and the last one. As the D loop presents some variation in length,
   lpair can take different values:
   lpair=14 smallest length of the D loop
   lpair=18 greatest length of the D loop */

                  ptr3=ptr1+2;
                  npair = 3;
                  for (i=0; i < 5; i++)
                     for (j=0; j < 2; j++)
                        ntab[i][j]=0;

                  for (lpair=14,h=0; lpair <=18; lpair++,h++)
                    {

/* Search for the presence of the D arm */

		      basepairing(ptr3,npair,lpair,&ncomp);

/* For each potential stem found, the number of base-pairings (ncomp) and
   lpair are stored in ntab */

		      ntab[h][0]=ncomp;
		      ntab[h][1]=lpair;
                    }
       
                  match2=0;
                  for (h=0; h < 5;h++)
                    {

/* If stems with 3 base-pairings are found for the D arm, the following 
   steps of the algorithm are applied on each of these potential stems. 
   The 2 base-pairings stems are discarded */

                    if (ntab[h][0] == 3)
                      {
                      lpair=ntab[h][1];
                      match2=1;
                      following_search(pos,pos1,ptr1,ptr3,lpair,nloop,
                      sequence,seqlen,fpo,fpverb,name,score,match2,&ntrna,&npred,&match,sqoffset);
                      }
                    }

/* If no stems with 3 base-pairings have been found, then the stems with 2
   base-pairings are used in the following steps */

                  if (!match2)
                    {
                    for (h=0; h < 5;h++)
                      {
                      if (ntab[h][0] == 2)
                        {
                        lpair=ntab[h][1];
                        following_search(pos,pos1,ptr1,ptr3,lpair,nloop,
                        sequence,seqlen,fpo,fpverb,name,score,match2,&ntrna,&npred,&match,sqoffset);
                        }
                      }
                    }
                  }
                score=first_score;
                }
              }
            }
          }
        }
      }

/* The sequence is complemented and the algorithm is applied on the 
   complementary strand */

    compstrand(&sequence,seqlen);

#ifdef VERBOSE
	   if (nloop == 0) 
	     fprintf(fpverb,"\n== Trying COMPLEMENTARY strand\n");
#endif 

    }

  free(sequence);  /* free up mem to ready for next sequence */

  }
fprintf(fpo,"number of sequences= %d\n", nseq);
fprintf(fpo,"number of bases tested (one strand)=%ld\n", lseq);
lseq = 2* lseq;
fprintf(fpo,"number of bases tested (both strands)= %ld\n", lseq);
fprintf(fpo,"number of predicted tRNA=%d\n", npred);
exit(0);
}


/* Subroutine that accomplishes the end of the test for the presence of a
   tRNA gene */

void 
following_search(long int pos,    /* first position of the found T-Psi-C signal */ 
		 long int pos1,   /* first position of the found D signal */
		 char *ptr1,      /* pointer to the first position of the D signal */
		 char *ptr3,      /* ptr3=ptr1+2  pointer to first position of the D arm  */
		 int lpair, 
		 int nloop,       /*nloop=0 test of the sequence, */
				  /*nloop=1 test of the complementary sequence */
		 char *sequence,  
		 long int seqlen, 
		 FILE *fpo,       /*pointer to the output file */
		 FILE *fpverb, 
		 char *name, 
		 int score,       /* Value of the general score SG */
		 int match2,     /* integer testing the presence of D */
				 /* arm with 3 base-pairings */ 
		 int *ntrna, 
		 int *npred, 
		 int *match,
		 long int sqoffset   /* offset nucleotide numbering by this much (set with -i param) */
		 )

{
extern Param_set_type ps;

char chaine2[300]; /* character string containing the predicted tRNA gene
                   sequence */
char anticodon[4]; /* character string containing the anticodon signal
                   sequence */
char type_trna[4]; /* character string containing the tRNA gene family */
char *ptr2,*ptr4,*ptr5; /*pointers to the sequence, see their definition
                        below*/
int score2, score1; /* variables for the general score SG */
int npair1=0; /* number of base-pairings in the aminoacyl arm */
int lpair1=0; /* number of nucleotides between the first position of
              the aminoacyl arm and the last one */
int lpair2=0; /* number of nucleotides between the first position of
              the anticodon arm and the last one */
int npair2=0; /* number of base-pairings in the anticodon arm */
int pos6=0; /* number of nucleotides found between the first position of
            the anticodon arm and the first position of the T-Psi-C 
            signal. */
int pos4=0; /* number of nucleotides present in the variable loop.*/
int length1=3; /* length of the anticodon signal */
int num=0; /* variable that codes the anticodon signal sequence
           (comprised between 1 and 65). */
int match1; /* match1=1 if a tRNA without intron has been predicted 
            and 0 otherwise*/
int ncomp=0; /* number of base-pairings observed in the sequence tested 
             for the presence of a given arm */
int i;

#ifdef VERBOSE
  fprintf(fpverb,"IN following search...\n");
#endif

score1=score;

/* If match2=1, the found D arm present 3 base-pairings, the general score 
   SG (score1) is incremented by 1
   Change in the threshold of base-pairings in the D arm for the increment
   of SG, MODIFY the following line */


if (match2) {
  score1++;
#ifdef VERBOSE
  fprintf(fpverb,"+ D arm found 3 bp\n"); }
else {
  fprintf(fpverb,"X D arm found less than 3 bp\n"); 
#endif 
}
npair1=7; /* Number of base-pairings in the aminoacyl arm */
lpair1=pos-pos1+8+23; /* Number of nucleotides between the first
                      position of the aminoacyl arm and the last one */
ptr2=(ptr1-7); /* pointer to the first position of the amino acyl arm */

/* Test for the presence of the aminoacyl arm */

basepairing(ptr2,npair1,lpair1,&ncomp);

/* If the number of base-pairings (ncomp) equal 7, the general score SG 
   (score1) is incremented by 1 .
   Change in the threshold of base-pairings in the aminoacyl arm for the 
   increment of SG, MODIFY the value '7' in the following line */

if (ncomp >= ps.aa_incsg) {
  score1++;
#ifdef VERBOSE
  fprintf(fpverb,"+ AA arm found %d base pairings\n",ncomp); }
else {
  fprintf(fpverb,"X AA arm found %d base pairings\n",ncomp); 
#endif 
}

/* Change in the threshold to retain a stem as a potential aminoacyl arm
   MODIFY the value '6' in the following line */

if (ncomp >= ps.aa_keep)
  {
  i = 0;
  while (i< 300)
    {
    chaine2[i]='\0';
    i++;
  }

/* If the general score SG is >= 4, the algorithm looks for the presence of
   an anticodon stem, otherwise the algorithm is initiated again on the
   following windowed sequence */

  if (score1 >= ps.look_for_acloop_sg)
    {
#ifdef VERBOSE
      fprintf(fpverb,"Looking for anticodon stem\n");
#endif

      ptr4=ptr3+lpair+2; /* pointer to the first position of the anticodon arm */
    match1=0;
    npair2=5; /* Number of base-pairings in the anticodon arm */
    lpair2=16; /* Number of nucleotides between the first position of the
               anticodon arm and the last position */

/* Test for the presence of an anticodon arm */

    basepairing(ptr4,npair2,lpair2,&ncomp);

/* If 4 or 5 base-pairings are observed, the windowed sequence is retained
   as a potential anticodon arm */

    if (ncomp >= ps.acloop_min)
      {
      if (nloop == 0)
        *match=1;
      i=0;
      while (i< 4)
        {
        anticodon[i]='\0';
        type_trna[i]='\0';
        i++;
        }

/* Test of the presence of the residue T at the position preceding the
   anticodon signal. If the residue T is present, the general score is
   incremented by 1.*/

      if ((*(ptr4+6)) == 't')  {
        score1++;
#ifdef VERBOSE
	fprintf(fpverb,"+ Invariant T found in anticodon. SG++\n"); }
      else  {
	fprintf(fpverb,"X Invariant T NOT found in anticodon. No SG inc\n"); 
#endif
      }

/* If SG >= 5, the windowed sequence is retained as a potential tRNA gene.
   The tRNA gene predicted at that level is without intron */

#ifdef VERBOSE
      if (score1 < ps.sg_cutoff)  {
	fprintf(fpverb,"X Under SG threshold: %d\n",score1); }
      else  {
	fprintf(fpverb,"+ Over SG threshold: %d\n",score1); 
#else
      if (score1 >= ps.sg_cutoff)  {
#endif 
      

/* Identification of the tRNA gene family */

        strncpy(anticodon,ptr4+7,3);
        codage(anticodon,length1,&num);
        corresaa(num,type_trna);
        strncpy(chaine2,(ptr1-7),lpair1+2);

        (*ntrna)++;
        (*npred)++;
        match1=1; /* match1=1 indicates that a tRNA gene without
                  intron has been predicted */

/* The results are printed in the output file */

        printresult(fpo,fpverb,name,pos1,pos,lpair,lpair1,lpair2,nloop,ntrna,
        chaine2,sequence,seqlen,match,ncomp,type_trna,anticodon,sqoffset); 
      }
    }

/* If no anticodon arm has been found for a value of lpair2=16, i. e. an
   anticodon loop of 7 nucleotides (without intron), the algorithm searchs
   for the presence of an anticodon arm assuming the presence of an intron
   in the anticodon loop */
 
/* pos6 corresponds to the number of nucleotides found between the first
   position of the anticodon arm and the first position of the T-Psi-C 
   signal. */ 

#ifdef VERBOSE
    fprintf(fpverb,"+ Looking for anticodon stem WITH intron...\n");
#endif 

    pos6=pos-(pos1+lpair+3);
  
/* The search for the presence of an anticodon stem is initiated only if an
   intron of at least 8 bases long is expected to be present. The value 28
   for pos6 corresponds to the number of nucleotides found for the smallest
   length of the variable loop and for an intron length of 8 nucleotides.
   Change in the smallest length of the intron, MODIFY the value of pos6
   in the following line and also the start value of lpair2 (lpair2=24) 
   in the loop (lpair2 is the number of nucleotides found between the first
   position of the anticodon arm and the last position of the anticodon arm.
   lpair2=16 when no intron in the anticodon loop). */

    if ((!match1) && (pos6 >=MIN_VAR_LOOP))
      {
      for (lpair2=MIN_VAR_LOOP-4;lpair2<pos6;lpair2++)
        {
        score2=score1;

/* Test for the presence of the anticodon arm for the different possible
   lengthes of the intron.
   ptr4 = pointer to the first base involved in the first part of the 
          anticodon arm */

        basepairing(ptr4,npair2,lpair2,&ncomp);

/* pos4 corresponds to the number of nucleotides present in the variable
   loop. The smallest size of the variable loop in the tRNA sequences 
   is 3. Thus to retain the stem as a potential anticodon arm, pos4 has 
   to be >= 3. */

        pos4=pos-pos1-lpair-lpair2-4;

/* ptr5 = pointer to the position preceding the intron. */

        ptr5=ptr4+10;

/* The anticodon arm is retained if 4 or 5 base-pairings are formed and if
   pos4 >= 3 */

        if ((ncomp >=4) && (pos4 >= 3))
          {

/* Test for the presence of the residues A or G at the position preceding
   the intron */

          if(((*ptr5) == 'a') || ((*ptr5) == 'g'))
            {
            if (nloop == 0)
              (*match)=1;
            i=0;
            while (i<4)
              {
              anticodon[i]='\0';
              type_trna[i]='\0';
              i++;
              }

/* Test for the presence of the residue T at the position preceding the
   anticodon signal. If the residue T is present, the general score is
   incremented by 1.*/

            if ((*(ptr4+6)) == 't') {
              score2++;
#ifdef VERBOSE
	      fprintf(fpverb,"+ Found invariant T. SG= %d\n",score2); }
      else  {
	      fprintf(fpverb,"X Did NOT find invariant T. SG= %d\n",score2); 
#endif 
	    }

/* If SG >= 5, the windowed sequence is retained as a potential tRNA gene.
   The tRNA gene predicted at that level has one intron in
   the anticodon loop. */

            if (score2 >= ps.sg_cutoff)
              {

/* Identification of the tRNA gene family */

              strncpy(anticodon,ptr4+7,3);
              codage(anticodon,length1,&num);
              corresaa(num, type_trna);
              strncpy(chaine2,(ptr1-7),lpair1+2);

              (*ntrna)++;
              (*npred)++;

/* The results are printed in the output file */

              printresult(fpo,fpverb,name,pos1,pos,lpair,lpair1,lpair2,nloop,
              ntrna,chaine2,sequence,seqlen,match,ncomp,type_trna,anticodon,sqoffset); 
	    }
	  }
	}
      }
    }
    }
    }
}



/* Subroutine to read the consensus matrix */

void 
lectval(FILE *fp,                  /* pointer to the consensus matrix file */
	float (*table_cons)[4],    /* table containing the frequency of each base
				      at each position of the signal */
	int (*table_inv)[2],  /* table containing the position and the nature of
				 the invariant bases found in the signal. Code for
				 the bases: A=0, C=1, G=2 and T=3. */
	int *lsig,            /*lsig=length of the signal */
	int *ktot,            /* ktot= number of invariant bases */
	float *maxtot         /* maxtot= sum of the maximum frequencies */
	)

{
int i=0,j, k=0,l,m;
float  max=0;
int ret = 0;

for (l=0; l< 30; l++)
  for (m=0; m<4; m++)
    table_cons[l][m]=0;

for (l=0; l< 30; l++)
  for (m=0; m < 2; m++)
    table_inv[l][m]=0;

while(feof(fp) == 0)
  {
  for(j=0;j<4;j++)
    {
    ret = fscanf(fp,"%f",&table_cons[i][j]);
    }
  i++;
  }
     
*lsig=i-1;
for (i=0; i<*lsig; i++)
  {
  for (j=0; j<4; j++)
    {
    if (table_cons[i][j] == 1.0)
      {
      k++;
      table_inv[k][1]=i;
      table_inv[k][2]=j;
      }
    }
  }

for (i=0; i<*lsig; i++)
  {
  max=table_cons[i][0];
  for (j=1; j<4; j++)
    max= MAX(max,table_cons[i][j]);
  *maxtot += max;
  }

*ktot=k;
}

/* Subroutine reading the sequence */
/* Modified to correctly read FASTA sequence files */

int fgetseq(char *name,           /* string w/name of the sequence */
	    char *sequence,       /* character string containing the sequence */
	    long int *seqlen,     /* length of sequence */
	    FILE *fpi)            /* input file pointer */
		
{
char line[MAXLINE]; /* character string used to read a line */
char *ptr; /* pointer to the sequence */
long int i,j,c;
char *ptrRet;

line[0]='\0';
*sequence='\0';
*seqlen = 0;

/* Test the first character to choose between the two formats available */

if (fgets(line, MAXLINE, fpi) == NULL) {
  return 0;
}

else if (line[0] == ';') {
  
  /* File in the non-GenBank format */

  if (line[1] != ' ')
    {      
      for (i=0, ptr= &(line[1]); *ptr != ' ' && *ptr !='\n';i++)
        name[i]= *ptr++;
      name[i] = '\0';
      while ((c = getc(fpi)) == ';')
        ptrRet = fgets(line, MAXLINE, fpi);
      ungetc(c, fpi);
      ptr = sequence;
      
      while ((c= getc(fpi)) != ';' && c != EOF)
	if (isalpha(c)) {
	  *ptr = tolower(c);
	  ptr++;
	  (*seqlen)++;
	}
      
      if (c != EOF)
	ungetc(c, fpi);
      *ptr= '\0';
    }
  else if (line[1] == ' ')
    {
      /* Intelligenetics format */

      while ((c = getc(fpi)) == ';')
        ptrRet = fgets(line,MAXLINE, fpi);
      ungetc(c, fpi);
      fgets(line, MAXLINE, fpi);
      for (i=0, ptr= &(line[0]); *ptr != ' ' && *ptr !='\n';i++)
        name[i]= *ptr++;
      name[i] = '\0';
      
      ptr = sequence;
      
      while ((c= getc(fpi)) != ';' && c != EOF)
	if (isalpha(c)){
	  *ptr = tolower(c);
	  ptr++;
	  (*seqlen)++;
	}
      
      if (c != EOF)
	ungetc(c, fpi);
      *ptr= '\0';
    }
  } 
else if (line[0] == '>')
  {
    /* Fasta Format */

    for (i=1; line[i] == ' '; i++)
      ;
    for (j = 0, ptr = &line[i]; *ptr != ' '; ptr++, j++)
      name[j] = *ptr;
    name[j] = '\0';

    ptr = sequence;
    while ((c= getc(fpi)) != '>' && c != EOF)
      if (isalpha(c)) {
	*ptr = tolower(c);
	*seqlen += 1;
	ptr++;
      }
      
    if (c != EOF)
      ungetc(c, fpi);
    *ptr= '\0';
  }
else 
  {
    /* File in GenBank format */
    
    while (strncmp(line, "LOCUS", 5) != 0)
      if (fgets(line, MAXLINE, fpi) == NULL)
	exit(1);
    
    for (i = 0, ptr = &line[12]; *ptr != ' '; ptr++, i++)
      name[i] = *ptr;
    name[i] = '\0';
    
    while (strncmp(line, "ORIGIN", 6) != 0)
      if (fgets(line, MAXLINE, fpi) == NULL)
	exit(1);
	    
    ptr = sequence;
    *seqlen=0;
    ptrRet = fgets(line, MAXLINE, fpi);
    while (strncmp(line, "//", 2) != 0) {
      for (i = 0; line[i] != '\n'; i++)
	if (isalpha(line[i])) {
	  *ptr++ = tolower(line[i]);
	  (*seqlen)++;
	}
      ptrRet = fgets(line, MAXLINE, fpi);
    }
    *ptr = '\0';
    
  }
return (1);
}

/* Subroutine reading & returning the sequence length */

int getseqsize(FILE *fpi       /* input file pointer */
	       )
	
{

char line[MAXLINE]; /* character string used to read a line */
long int i,c, seqlen, fpi_save_pos;
char* ptrRet;

line[0]='\0';
seqlen = 0;

fpi_save_pos = ftell(fpi);  /* save current position in file */

/* Test the first character to choose between the two formats available */

if (fgets(line, MAXLINE, fpi) == NULL) {
  return 0;
}

else if (line[0] == ';') {
  
  /* File in the non-GenBank format */

  if (line[1] != ' ')
    {      
      while ((c = getc(fpi)) == ';')
	ptrRet = fgets(line, MAXLINE, fpi);
      ungetc(c, fpi);
      
      while ((c= getc(fpi)) != ';' && c != EOF)
	if (isalpha(c)) {
	  seqlen++;
	}
      
    }
  else if (line[1] == ' ')
    {
      /* Intelligenetics format */

      while ((c = getc(fpi)) == ';')
	ptrRet = fgets(line,MAXLINE, fpi);
      ungetc(c, fpi);
      ptrRet = fgets(line, MAXLINE, fpi);
      
      while ((c= getc(fpi)) != ';' && c != EOF)
	if (isalpha(c)){
	  seqlen++;
	}
    }
} 
else if (line[0] == '>')
  {
    /* Fasta Format */

    while ((c= getc(fpi)) != '>' && c != EOF)
      if (isalpha(c)) 
	seqlen++;
  }
else 
  {
    /* File in GenBank format */
    
    while (strncmp(line, "LOCUS", 5) != 0)
      if (fgets(line, MAXLINE, fpi) == NULL)
	exit(1);
    
    while (strncmp(line, "ORIGIN", 6) != 0)
      if (fgets(line, MAXLINE, fpi) == NULL)
	exit(1);

    ptrRet = fgets(line, MAXLINE, fpi);
    while (strncmp(line, "//", 2) != 0) {
      for (i = 0; line[i] != '\n'; i++)
	if (isalpha(line[i])) {
	  seqlen++;
	}
      ptrRet = fgets(line, MAXLINE, fpi);
    }
    
  }

fseek(fpi,fpi_save_pos,0);  /* reposition file pointer */
return seqlen;
}



/* Subroutine that returns the position in the string s where the string t
   begins or -1 if s does not contain t */

/* Calls to this function eliminated for efficiency  T. Lowe  11/95  */

myindex (char *s, char *t)
{
int i, j, k;
for (i=0; s[i] != '\0'; i++) {
  for (j=i, k=0; t[k] != '\0' && s[j] == t[k]; j++, k++)
    ;
  if (t[k] == '\0')
  return(i);
}
 return(-1);
}


/* Subroutine looking for the presence of a given signal, returns 1 if a
   signal is found and 0 otherwise. It also return the table 'weight'
   containing the frequencies of the oberved bases in the windowed
   sequence and the number 'ninv' of invariant bases found in the windowed
   sequence*/

int readsignal(char *ptr,      /* pointer to the sequence */
	       int (*table_inv)[2],  /* table containing the position and nature of the
					invariant bases found in the consensus matrix */
	       int *lsig,   /* lsig= length of the signal */
	       int *ktot,   /* ktot= number of invariant bases in the consensus matrix */
	       float *weight,       /* table containing the frequencies of the observed
				       base at each position of the
				       windowed sequence tested */
	       float (*table_cons)[4],  /* table corresponding to the consensus matrix */
	       int *ninv,               /* ninv= number of invariant */
					/* bases in the windowed  sequence */
	       int threshold_inv)       /* Number of invariant bases */
					/* allowed not to b e invariant */

  
{
  extern Param_set_type ps;

  int k = 1, i=0, j,match1;
  int l;
  int tab[30];

  (*ninv)=0;
  for (l=0; l< 30; l++) {
    weight[l]=0;
    tab[l]=0;
  }

/* If the consensus matrix contains some invariant bases, the subroutine
   calculates the number of invariant bases that are found in the windowed
   sequence (ninv) */

  if (*ktot)
    {
      for (k=1; k <= (*ktot); k++)
	{
      
/*   (original code commented out)

     temp[0]= *(ptr+table_inv[k][1]);
      j = myindex(base,temp);              */
   
	  switch (*(ptr+table_inv[k][1])) {
	  case 'a': j=0; break;
	  case 'c': j=1; break;
	  case 'g': j=2; break;
	  case 't': j=3; break;
	  default: j=-1;
	  }

/* trap for non-ATGC chars, assume match for ambiguous bases (j= -1)*/


#ifdef NO_AMBIG
	if (j == table_inv[k][2]) 
#else
	if ((j == table_inv[k][2]) || (j == -1)) 		 
#endif
	    (*ninv)++;

	}
    }
 
/* If the number of invariant bases found in the windowed sequence is < to
   the threshold allowed, the windowed sequence is discarded as a potential
   signal. */

  if ((*ninv) < (*ktot)-threshold_inv)
    return(0);

/* If the number of invariant bases is >= to the threshold, the table 
  'weight' is constructed */

  if((*ninv) >= (*ktot)-threshold_inv)
    {
      match1=1;
      while(match1 && (i<*lsig))
	{
	  
/*  (original code commented out) 

    temp[0]= *(ptr+i);
    tab[i]= myindex(base,temp);  */

	  switch (*(ptr+i)) {
	  case 'a': tab[i]=0; break;
	  case 'c': tab[i]=1; break;
	  case 'g': tab[i]=2; break;
	  case 't': tab[i]=3; break;
	  default:  tab[i]=-1;
	  }


/* trap for non-ATGC chars, assume match for ambig bases (tab[i]= -1) */

	  if (tab[i] == -1) 
#ifdef NO_AMBIG
	    weight[i] = 0;   
#else
	    weight[i] = 1;    
#endif
	  else
	    weight[i]=table_cons[i][tab[i]];
	  
	  if((ktot == 0) && (weight[i] == 0))
	    {
	      match1=0;
	      return(0);
	    }
	  else
	    {
	      i++;
	    }
	}
    }
  return(1);
}

/* Subroutine that calculates the similarity score on the potential signal 
   previously retained by the subroutine readsignal. This subroutine
   returns 1 it the computed score is greater or equal to the defined 
   threshold and 0 otherwise. It returns also the value of the computed 
   score (score) */ 

int scoring(float *weight, /* table containing the frequencies of the observed base 
			      at each position of the potential signal */
	    int lsig,      /* length of the signal */
	    float max,     /* sum of the maximum frequencies found in */
			   /* the consensus matrix */ 
	    int ktot,      /* number of invariant bases found in the */
			   /* consensus matrix */
	    float *score,  /* value of the computed score on the */
			   /* potential signal */
	    float ThresholdValue,   /* defined threshold for the */
			 	    /* similarity score */
	    int ninv       /* number of invariant bases found */
		           /* in the potential signal */
	    )
{

float tot;
int i;

/* Computation of the value of the score on the potential signal */

tot=0;
for(i=0; i< lsig; i++)
  tot += weight[i];
 
tot -= ninv;
max -= ktot;      
*score = tot / max;

/* Comparison of the computed score with the defined threshold value */

if (*score >= ThresholdValue)
  return(1);
else {
  return(0);
}
}

/* Subroutine looking for base-pairings between two parts of the sequence.
   It returns the number of base-pairings found (ncomp) */

/* rewritten to improve efficiency & eliminate use of myindex() calls */

void 
basepairing(char *ptr,   /* pointer to the sequence */
	    int npair,   /* number of base-pairings forming a given arm */
	    int lpair,   /* number of nucleotides found between the first position of the
			    first part of the sequence involved in the stem and the last 
			    position of the second part of the
			    sequence involved in the stem */
	    int *ncomp)   /* number of base-pairings observed between the two parts of
			    the sequence tested */
  
{
  int n;   /* loop counter */ 

  *ncomp=0;

#ifdef NO_AMBIG
  for(n=0; n<npair; n++)
    {
      switch (*(ptr+n)) {
	case 'a': if (*(ptr+lpair-n) == 't') (*ncomp)++;  break;
	case 'c': if (*(ptr+lpair-n) == 'g') (*ncomp)++;  break;
	case 'g': case 'r': if ((*(ptr+lpair-n) == 'c') || (*(ptr+lpair-n) == 't'))
	  (*ncomp)++;  break;
	case 't': case 'y': if ((*(ptr+lpair-n) == 'a') || (*(ptr+lpair-n) == 'g'))
	  (*ncomp)++;  break;
         }  	
    }    
#else
  for(n=0; n<npair; n++)
    {
      if ((*(ptr+lpair-n) == 'n') && (*(ptr+n) != 'n')) {
	(*ncomp)++;   }
      else  {   
	switch (*(ptr+n)) {
	case 'a': if (*(ptr+lpair-n) == 't') (*ncomp)++;  break;
	case 'c': if (*(ptr+lpair-n) == 'g') (*ncomp)++;  break;
	case 'g': case 'r': if ((*(ptr+lpair-n) == 'c') || (*(ptr+lpair-n) == 't'))
	  (*ncomp)++;  break;
	case 't': case 'y': if ((*(ptr+lpair-n) == 'a') || (*(ptr+lpair-n) == 'g'))
	  (*ncomp)++;  break;
	case 'n': if (*(ptr+lpair-n) != 'n') (*ncomp)++;  break; 
	}
      }  	
    }    
#endif

}

/* Subroutine looking for base-pairings between two parts of the sequence.
   It returns the number of base-pairings found (ncomp) */

/* replaced by new basepairing() function above  */ 

void 
old_basepairing(char *ptr, int npair, int lpair, int *ncomp)

{

char base[5]; /* character string containing the 4 bases A, C, G and T */
char temp1[2], temp2[2]; /* temporary tables */
int n,j1,j2;

/* Table 'pairing' gives the base-pairing scheme. If two nucleotides can 
   form a base-pairing, the value in the table is 1 otherwise it is 0. */

static int pairing[4][4]= {            
                          {0,0,0,1},
                          {0,0,1,0},
                          {0,1,0,1},
                          {1,0,1,0},
                          };         

base[0]='a', base[1]='c', base[2]='g', base[3]='t', base[4]='\0';
temp1[0]='\0';
temp1[1]='\0';
temp2[0]='\0';
temp2[1]='\0';

*ncomp=0;
for(n=0; n<npair; n++)
  {

    temp1[0]= *(ptr+n);  /* temp1[0] contains the base found at a given
			    position in the first part of the stem */
                        
    temp2[0]= *(ptr+lpair-n); /* temp2[0] contains the base found at the 
				 position of the second part of the stem that 
				 should form the base-pairing with temp1[0] */
    j1= myindex(base,temp1); /* Code the base found in temp1[0]. j1 can take the
				value 0 (for A), 1 (for C), 2 (for G) and 3 (for T) */
                           
    j2= myindex(base,temp2); /* same than j1 for the base in temp2[0] */

    /* Test if the bases in temp1[0] and temp2[0] can form a base-parings */

/* Assume non-ATGC base can bp with anything  */

#ifdef NO_AMBIG
    if (pairing[j1][j2] == 1) 
      *ncomp= (*ncomp)+1;
#else
    if (j1 == -1 || j2 == -1) 
      *ncomp= (*ncomp)+1;    
    else if (pairing[j1][j2] == 1) 
      *ncomp= (*ncomp)+1;
#endif
  }
}



/* Subroutine that complements the sequence. It returns the complementary 
   sequence to the main */

/* Modifed to dynamically allocate space for reverse complement & then */
/* release memory occupied by original sequence */
/* Sets *sequence pointer to new reverse complement sequence */

void 
compstrand(char **sequence,   /* pointer to the sequence string */
	   long int seqlen    /* sequence length   */
	   )
{

char *ptr, *revseq;
long int pos;


revseq= (char *) calloc (seqlen+1,sizeof(char)); 

if (revseq == NULL) 
  {
  fprintf(stderr,"tRNAscan1.4: Can't malloc for complement sequence\n");
  exit(-1);
  }


/* Complementation of the sequence */

for (pos=0, ptr= &((*sequence)[seqlen-1]); pos<seqlen; pos++, ptr--) 
  {
    switch (*ptr)
      {
      case 'a':
	revseq[pos]='t';
	break;
      case 'c':
	revseq[pos]='g';
	break;
      case 'g':
	revseq[pos]='c';
	break;
      case 't':
	revseq[pos]='a';
	break;
      default:
	revseq[pos]='n';
	break;
      }
  }

revseq[seqlen] = '\0';

free(*sequence);           /* release mem from original seq */
*sequence = revseq;        

} 


/* Subroutine that codes the anticodon signal sequence by a number
   comprised between 1 and 65. It returns this number. */

/* Modified to eliminate call to myindex()  */

void 
codage(char *anticodon,    /* anticodon signal sequence */
       int length1,        /* length of the anticodon signal, lenght1=3 */
       int *num            /* number associated to the anticodon signal sequence */
       )            
         
{
int j=0,i,iba,match=0;
iba=1;
*num=1;

/* Codage of the anticodon signal sequence. The codage is done by 
   alphabetical order, i. e., AAA is coded by 1, AAC is coded by 2 
   and so on until TTT that is coded by 64. If a base is not determined
   in the anticodon signal sequence, the anticodon is coded by 65. */

for (i=length1; i>=1;i--)
  {

    switch (anticodon[i-1]) {
    case 'a': j=0; break;
    case 'c': j=1; break;
    case 'g': j=2; break;
    case 't': j=3; break;
    default: j=-1;
    }

    if (j == -1)
      {
	match=1;
	i=1;
      }
    *num=(*num)+j*iba;
    iba=4*iba;
  }

if (match)
  *num=65;
  
}


/* Subroutine that determines the tRNA gene family */

void 
corresaa(int num,           /* Number coding the anticodon signal sequence */
	 char *type_trna    /* tRNA gene family */
	 )
                 
{

/*The table 'amino_acid' gives the correspondence between the number 
  (num) coding the anticodon signal and the amino acid that is added 
  to the protein by the tRNA. For example, if the anticodon is AAA then
  num=1 and the table will associate the amino acid Phe to that anticodon. 
  The tRNA gene in that case will be a Phe tRNA gene. */

static char *amino_acid[]= {"Phe","Val","Leu","Ile","Cys","Trp","Arg","Ser",
                           "Ser","Ala","Pro","Thr","Tyr","Asp","His","Asn",
                           "Leu","Val","Leu","Met","Trp","Gly","Arg","Arg",
                           "Ser","Ala","Pro","Thr","Sup","Glu","Gln","Lys",
                           "Phe","Val","Leu","Ile","Cys","Gly","Arg","Ser",
                           "Ser","Ala","Pro","Thr","Tyr","Asp","His","Asn",
                           "Leu","Val","Leu","Ile","Sup","Gly","Arg","Arg",
                           "Ser","Ala","Pro","Thr","Sup","Glu","Gln","Lys",
                           "Ind"};

strncpy(type_trna,amino_acid[num-1],3);
}

/* Subroutine that prints the results of the search */

void 
printresult(FILE *fpo,     /* output file pointer */
	    FILE *fpverb,   /* character string for the name of the sequence */
	    char *name, 
	    long int pos1,    /* first position of the D signal */
	    long int pos,     /* first position of the T-Psi-C signal */
	    int lpair,        /* number of nucleotides between the */
			      /* first position  of the D arm and the last one */
	    int lpair1,     /* number of nucleotides between the first position of the
			       aminoacyl arm and the last one */
	    int lpair2,     /* number of nucleotides between the first position of the 
			       anticodon arm and the last one */
	    int nloop,      /* nloop=0, scanning of the direct strand; nloop=1 scanning of 
			       the complementary strand */
	    int *ntrna,     /* number of tRNA genes predicted in the sequence */
	    char *chaine2,  /* character string for the predicted tRNA gene sequence */
	    char *sequence,    /* character string containing the sequence tested */
	    long int length,   /* length of the sequence */
	    int *match,      /* match=1 if at least one tRNA gene has been found on the 
				direct strand and 0 otherwise */
	    int ncomp,       /* number of base-pairings in the anticodon arm of the
				predicted tRNA gene */
	    char *type_trna,     /* character string for the tRNA gene family */
	    char *anticodon,      /* character string for the anticodon signal sequence */
	    long int sqoffset   /* offset nucleotide numbering by this much (set with -i param) */
	    ) 

{
long int pos2; /* first position of the predicted tRNA gene when it
               is found on the complementary strand */
long int posstart; /* first position of the intron */
long int posend; /* last position of the intron */


/* If lpair2 > 16, treatment of the tRNA gene with one intron */

/* Results for the tRNA genes predicted on the direct strand */

#ifdef VERBOSE
fprintf(fpverb,"*** tRNA found\n\n");
#endif 


if ((nloop) == 0)
  {
  if((*ntrna) == 1)
    fprintf(fpo,"sequence name= %s\n", name);

  fprintf(fpo,"start position= %ld end position= %ld\n",pos1-7+sqoffset,pos1-6+lpair1+sqoffset);
  fprintf(fpo,"potential tRNA sequence= %s\n",chaine2);
  fprintf(fpo,"D signal= %ld %ld TpsyC signal= %ld %ld\n", pos1+sqoffset,pos1+7+sqoffset, pos+sqoffset,
  pos+14+sqoffset);
  fprintf(fpo,"amino-acyl stem= %ld-%ld;%ld-%ld\n",pos1-7+sqoffset,pos1-1+sqoffset,pos1-13+lpair1+sqoffset,
  pos1-7+lpair1+sqoffset);
  fprintf(fpo,"D stem= %ld-%ld;%ld-%ld\n",pos1+2+sqoffset,pos1+4+sqoffset,pos1+lpair+sqoffset,
	  pos1+lpair+2+sqoffset);

  if(lpair2 > 16)
    {
    fprintf(fpo,"anticodon stem= %ld-%ld;%ld-%ld\n",pos1+lpair+4+sqoffset,pos1+lpair+8+sqoffset,
    pos1+lpair+lpair2+sqoffset,pos1+lpair+lpair2+4+sqoffset);
    }
  else
    {
    fprintf(fpo,"anticodon stem= %ld-%ld;%ld-%ld\n",pos1+lpair+4+sqoffset,pos1+lpair+8+sqoffset,
    pos1+lpair+16+sqoffset,pos1+lpair+20+sqoffset);
    }

  fprintf(fpo,"TpsyC stem= %ld-%ld;%ld-%ld\n",pos+1+sqoffset,pos+5+sqoffset,pos+13+sqoffset,pos+17+sqoffset);
  if (strcmp(type_trna,"Ind") != 0)
    {
    fprintf(fpo,"tRNA predict as a tRNA- %s : anticodon %s\n", type_trna,
    anticodon);
    }
  else
    {
    fprintf(fpo,"anticodon includes unknown bases\n");
    } 

  if (lpair2 > 16)
   {
    posstart=pos1+lpair+15;
    posend= pos1+lpair+lpair2-2;
    fprintf(fpo,"potential intron between positions %ld %ld\n",posstart+sqoffset,
    posend+sqoffset); 
    }
  fprintf(fpo,"number of base pairing in the anticodon stem= %d\n",ncomp);
  fprintf(fpo,"\n");
  }
else
  {

/* Results on the complementary strand */

/* If no tRNA gene is predicted on the direct strand, then for the first
   tRNA gene predicted on the complementary strand the two following lines
   are printed in the output file */  
 
  if( (!(*match)) && ((*ntrna) == 1))
    {
    fprintf(fpo,"sequence name= %s\n", name);
    fprintf(fpo,"complementary strand\n");
    }

/* If tRNA genes are predicted on the direct strand, then for the first
   tRNA gene predicted on the complemetary strand the following line is
   printed in the output file */

  else if ((*ntrna) == 1)
    {
    fprintf(fpo,"complementary strand\n");
    }

  pos2= length-pos1+8;
  fprintf(fpo,"start position= %ld end position= %ld\n",pos2+sqoffset,pos2-lpair1-1+sqoffset);
  fprintf(fpo,"potential tRNA sequence= %s\n",chaine2);
              
  fprintf(fpo,"D signal= %ld %ld TpsyC signal= %ld %ld\n",length-pos1+1+sqoffset,
  length-pos1-6+sqoffset,length-pos+1+sqoffset,length-pos-13+sqoffset); 
  fprintf(fpo,"amino-acyl stem= %ld-%ld;%ld-%ld\n",pos2+sqoffset,pos2-6+sqoffset, pos2-lpair1+6+sqoffset,
  pos2-lpair1+sqoffset);
  fprintf(fpo,"D stem= %ld-%ld;%ld-%ld\n",length-pos1-1+sqoffset,length-pos1-3+sqoffset,
  length-pos1-lpair+1+sqoffset,length-pos1-lpair-1+sqoffset);

  if (lpair2 > 16)
    {
    posstart=pos1+lpair+15;
    posend=pos1+lpair+lpair2-2;
    fprintf(fpo,"anticodon stem= %ld-%ld;%ld-%ld\n",length-pos1-lpair-3+sqoffset,
    length-pos1-lpair-7+sqoffset,length-posend-1+sqoffset,length-posend-5+sqoffset);
    }
  else
    {
    fprintf(fpo,"anticodon stem= %ld-%ld;%ld-%ld\n",length-pos1-lpair-3+sqoffset,
    length-pos1-lpair-7+sqoffset,length-pos1-lpair-lpair2+1+sqoffset,
    length-pos1-lpair-lpair2-3+sqoffset);
    }

  fprintf(fpo,"TpsyC stem= %ld-%ld;%ld-%ld\n",length-pos+sqoffset,length-pos-4+sqoffset,
  length-pos-12+sqoffset,length-pos-16+sqoffset);

  if (strcmp(type_trna,"Ind") != 0)
    {
    fprintf(fpo,"tRNA predict as a tRNA- %s : anticodon %s\n", type_trna,
    anticodon);
    }
  else
    {
    fprintf(fpo,"anticodon includes unknown bases\n"); 
    }

  if (lpair2 > 16)
    {
    posstart=pos1+lpair+15;
    posend=pos1+lpair+lpair2-2;
    fprintf(fpo,"potential intron between positions %ld %ld\n",
    length-posstart+1+sqoffset ,length-posend+1+sqoffset);
    }

  fprintf(fpo,"number of base pairing in the anticodon stem=%d\n",ncomp);
  fprintf(fpo,"\n");
  }
}


void 
set_search_params (Param_set_type *ps,
		   int params)
  {
    
    if (params == 1) {
      ps->tpc_sig_thresh = ST_TPC_SIG_THRESH;
      ps->d_sig_thresh = ST_D_SIG_THRESH;
      ps->sg_cutoff = ST_SG_CUTOFF; 
      ps->tpc_inv = ST_TPC_INV;
      ps->tpc_incsg = ST_TPC_INCSG;
      ps->tpc_keep = ST_TPC_KEEP;
      ps->d_inv = ST_D_INV;
      ps->look_for_acloop_sg = ST_LOOK_FOR_ACLOOP_SG;
      ps->acloop_min = ST_ACLOOP_MIN;  
      ps->aa_incsg = ST_AA_INCSG;
      ps->aa_keep = ST_AA_KEEP;
    }
    else if (params == 2) {
      ps->tpc_sig_thresh = RX_TPC_SIG_THRESH;
      ps->d_sig_thresh = RX_D_SIG_THRESH;
      ps->sg_cutoff = RX_SG_CUTOFF; 
      ps->tpc_inv = RX_TPC_INV;
      ps->tpc_incsg = RX_TPC_INCSG;
      ps->tpc_keep = RX_TPC_KEEP;
      ps->d_inv = RX_D_INV;
      ps->look_for_acloop_sg = RX_LOOK_FOR_ACLOOP_SG;
      ps->acloop_min = RX_ACLOOP_MIN;  
      ps->aa_incsg = RX_AA_INCSG;
      ps->aa_keep = RX_AA_KEEP;
    }
    else if (params == 3) {
      ps->tpc_sig_thresh = ALT_TPC_SIG_THRESH;
      ps->d_sig_thresh = ALT_D_SIG_THRESH;
      ps->sg_cutoff = ALT_SG_CUTOFF; 
      ps->tpc_inv = ALT_TPC_INV;
      ps->tpc_incsg = ALT_TPC_INCSG;
      ps->tpc_keep = ALT_TPC_KEEP;
      ps->d_inv = ALT_D_INV;
      ps->look_for_acloop_sg = ALT_LOOK_FOR_ACLOOP_SG;
      ps->acloop_min = ALT_ACLOOP_MIN;  
      ps->aa_incsg = ALT_AA_INCSG;
      ps->aa_keep = ALT_AA_KEEP;
    }
    else {
      fprintf (stderr,"tRNAscan1.4: FATAL: Unable to select search parameter set.\n");
      exit(1);
    }
    
  }
 











