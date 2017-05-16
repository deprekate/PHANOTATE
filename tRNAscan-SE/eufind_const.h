
/* tRNA scanning cutoffs  */

#define BBOX_CUTOFF -14.14
#define BBOX_START_IDX 45

#define SEC_LOBOUND -4.9    /* orig: -3.6 euk: -3.8  prok (hflu) -4.9 */
#define SEC_HIBOUND -2.1     /* -2.2 */
#define MAX_PENALTY -5.442   /* log(1/231) */


#define INT_SCORE_THRESH -31.25
#define TOT_SCORE_THRESH -31.8  /* -31.8 */

#define MAX_AB_BOX_DIST 140  /* not used anymore, instead
			     /* AB_BOX_DIST_RANGE used */
#define MIN_AB_BOX_DIST 24
#define AB_BOX_DIST_RANGE 116  /* check this far over MIN_AB_BOX
			       /* distance for A-B box pairs */

#define SEC_AB_BOX_DIST 26
#define SEC_BBOX_DIST_CORR 12

#define MIN_BTERM_DIST 11
#define MAX_TERM_SEARCH 133   /* Max distance to search for termination
			      signal (was 59, changed to 133 (as in
			      Pavesi paper) on  11/96
			      since was missing 4 yeast tRNAs */
#define ABOX_LEN 21
#define BBOX_LEN 11

#define MAX_OVLAP 10        /* max #bp tRNA hits are allowed to
			       overlap */

struct trna_info_s {
  char iso_type[5];
  char acodon[4];
  int start, end, 
    Abox_st, Abox_end, Abox_gap, 
    Bbox_st, Bbox_end, 
    Term_st,
    acodon_idx, intron, idno;
  float totSc, AboxSc, BboxSc, ABdistSc, TermSc;
};

typedef struct trna_info_s TRNA_TYPE;





