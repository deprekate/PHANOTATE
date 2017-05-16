/* Header file specifically for MasPar parallelized version
 */

/* These define's affect mpcovels.
 */
#define NYPROC    64	/* how many rows of processors (nyproc)             */
#define BLOCKSIZE 128	/* size of sequence blocks sent to DPU              */
#define VPENUM    3	/* # of virtual PE's per PE. N = NYPROC * VPENUM -1 */


