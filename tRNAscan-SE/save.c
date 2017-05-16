/* save.c
 * Saving models to disk, and reading them back in
 * SRE, Wed Sep  8 17:14:43 1993
 * 
 * Both binary and flat text save formats are supported.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "funcs.h"
#include "structs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* The magic number for our binary file is "cm20" + 0x80808080 */
static int v20magic = 0xe3edb2b0;

static int read_cm20   (FILE *fp, struct cm_s **ret_cm);
static int read_bincm20(FILE *fp, struct cm_s **ret_cm);

/* Function: WriteCM()
 * 
 * Purpose:  Print a flat text copy of the data from a model
 *           to a file handle.
 *           
 * Returns:  1 on success, 0 on failure.
 */
int
WriteCM(FILE *fp, struct cm_s *cm)
{
  int i, j, k;

				/* header info */
  fprintf(fp, "### cove V2\n");
  fprintf(fp, "%d \tnodes\n",  cm->nodes);

				/* over all nodes, 0..cm->nodes-1 */
  for (k = 0; k < cm->nodes; k++)
    {
      fprintf(fp, "### node %d type %d\n", k, cm->nd[k].type);
      fprintf(fp, "%d  %d\n", cm->nd[k].nxt, cm->nd[k].nxt2);
      
				/* transitions */
      for (i = 0; i < STATETYPES; i++)
	{
	  for (j = 0; j < STATETYPES; j++)
	    fprintf(fp, "%.5f ", cm->nd[k].tmx[i][j]);
	  putc('\n', fp);
	}
				/* INSL emissions */
      for (i = 0; i < ALPHASIZE; i++)
	fprintf(fp, "%.5f ", cm->nd[k].il_emit[i]);
      fputs("# INSL\n", fp);
				/* INSR emissions */
      for (i = 0; i < ALPHASIZE; i++)
	fprintf(fp, "%.5f ", cm->nd[k].ir_emit[i]);
      fputs("# INSR\n", fp);

				/* MATP emissions */
      for (i = 0; i < ALPHASIZE; i++)
	{
	  for (j = 0; j < ALPHASIZE; j++)
	    fprintf(fp, "%.5f ", cm->nd[k].mp_emit[i][j]);
	  fputs("# MATP\n", fp);
	}
				/* MATL emissions */
      for (i = 0; i < ALPHASIZE; i++)
	fprintf(fp, "%.5f ", cm->nd[k].ml_emit[i]);
      fputs("# MATL\n", fp);
				/* MATR emissions */
      for (i = 0; i < ALPHASIZE; i++)
	fprintf(fp, "%.5f ", cm->nd[k].mr_emit[i]);
      fputs("# MATR\n", fp);
    }
  return 1;
}



/* Function: WriteBinaryCM()
 * 
 * Purpose:  Save a model in binary format; somewhat nonportable
 *           but much more compressed.
 *           
 * Args:     fp - open file pointer to write to          
 *           cm - model to save
 *           
 * Return:   1 on success, 0 on failure.
 */
int
WriteBinaryCM(FILE  *fp, struct cm_s *cm)
{
  int k;			/* counter for nodes */

  /* Write the four-byte magic number. It identifies the file
   * as binary (because the high bits are set), and identifies
   * the major version of the program. It can be used to
   * identify and possibly humor byte-swapped architectures,
   * although we don't bother for now.
   */
  fwrite((void *) &v20magic, 4, 1, fp);

				/* header info */
  fwrite((void *) &(cm->nodes),  sizeof(int), 1, fp);

				/* loop over all nodes in model */
  for (k = 0; k < cm->nodes; k++)
    {
				/* type of node */
      fwrite((void *) &cm->nd[k].type,  sizeof(int), 1, fp);

				/* indices of child nodes */
      fwrite((void *) &cm->nd[k].nxt,  sizeof(int), 1, fp);
      fwrite((void *) &cm->nd[k].nxt2, sizeof(int), 1, fp);
      
				/* transitions */
      fwrite((void *) cm->nd[k].tmx, sizeof(double), STATETYPES*STATETYPES, fp);

				/* INS emissions */
      fwrite((void *) cm->nd[k].il_emit, sizeof(double), ALPHASIZE, fp);
      fwrite((void *) cm->nd[k].ir_emit, sizeof(double), ALPHASIZE, fp);
      
				/* MAT emissions */
      fwrite((void *) cm->nd[k].mp_emit, sizeof(double), ALPHASIZE * ALPHASIZE, fp);
      fwrite((void *) cm->nd[k].ml_emit, sizeof(double), ALPHASIZE, fp);
      fwrite((void *) cm->nd[k].mr_emit, sizeof(double), ALPHASIZE, fp);
    }
  return 1;
}




/* Function: ReadCM()
 * 
 * Purpose:  Read a flat text copy of the data from a model
 *           from a file.
 *           
 * Returns:  1 on success, 0 on failure.
 */
int
ReadCM(char *filename, struct cm_s **ret_cm)
{
  FILE           *fp;
  char            buffer[512];
  int             magic_number;

  /* Open file for reading
   */
  if ((fp = fopen(filename, "r")) == NULL)
    {
      Warn("Cannot open model file %s for reading", filename);
      return 0;
    }

  /* Look for "magic" header and dispatch reading
   * to the appropriate routine. First we check the leading
   * 4 bytes to see if it's a binary save file.
   */
  if (! fread((void *) &magic_number, 4, 1, fp)) 
    Die("Failed to read magic number from model file %s", filename);
  if (magic_number == v20magic)
    {
      if (! read_bincm20(fp, ret_cm))
	Die("Failed to read binary model file %s", filename);
    }
  else
    {
      rewind(fp);

      if (fgets(buffer, 512, fp) == NULL)
	{
	  Warn("ain't no data in the model file %s, pal", filename);
	  return 0;
	}
      if (strncmp(buffer, "### cove V2", 11) == 0)
	{
	  if (! read_cm20(fp, ret_cm)) return 0;
	}
      else
	{
	  Warn("File %s is not a recognized covariance model format", filename);
	  return 0;
	}
    }

  /* We're vulnerable to some roundoff error when we've read
   * files in; make sure all probabilities sum to 1.
   */
  NormalizeCM(*ret_cm);
  fclose(fp);
  return 1;
}


/* Function: read_cm20()
 * 
 * Purpose:  Read flat text model files from version 2.0 of the
 *           package. The file pointer fp is positioned on the line
 *           just after the "magic" header. Allocates, reads in,
 *           and returns the model.
 */
static int
read_cm20(FILE *fp, struct cm_s **ret_cm)
{
  struct cm_s *cm;
  int          i, j, k;
  int          nodes;
  int		   ret = 0;


				/* header info */
  ret = fscanf(fp, "%d \tnodes\n", &nodes);
  
  /* Given that header info, alloc for a model.
   */
  cm = AllocCM(nodes);
  if (cm == NULL) { Warn("Failed to allocate model"); return 0; }

				/* over all nodes, 0..nodes-1 */
  for (k = 0; k < nodes; k++)
    {
      ret = fscanf(fp, "### node %*d");
      ret = fscanf(fp, " type %d\n", &cm->nd[k].type);
      ret = fscanf(fp, "%d  %d\n", &cm->nd[k].nxt, &cm->nd[k].nxt2);
      
				/* transitions */
      for (i = 0; i < STATETYPES; i++)
	{
	  for (j = 0; j < STATETYPES; j++)
	    ret = fscanf(fp, "%lf ", &cm->nd[k].tmx[i][j]);
	  ret = fscanf(fp, "\n");
	}
      
				/* INSL emissions */
      for (i = 0; i < ALPHASIZE; i++)
	ret = fscanf(fp, "%lf ", &cm->nd[k].il_emit[i]);
      ret = fscanf(fp, "# INSL\n");
				/* INSR emissions */
      for (i = 0; i < ALPHASIZE; i++)
	ret = fscanf(fp, "%lf ", &cm->nd[k].ir_emit[i]);
      ret = fscanf(fp, "# INSR\n");

				/* MATP emissions */
      for (i = 0; i < ALPHASIZE; i++)
	{
	  for (j = 0; j < ALPHASIZE; j++)
	    ret = fscanf(fp, "%lf ", &cm->nd[k].mp_emit[i][j]);
	  ret = fscanf(fp, "# MATP\n");
	}
				/* MATL emissions */
      for (i = 0; i < ALPHASIZE; i++)
	ret = fscanf(fp, "%lf ", &cm->nd[k].ml_emit[i]);
      ret = fscanf(fp, "# MATL\n");
				/* MATR emissions */
      for (i = 0; i < ALPHASIZE; i++)
	ret = fscanf(fp, "%lf ", &cm->nd[k].mr_emit[i]);
      ret = fscanf(fp, "# MATR\n");
    }
  *ret_cm = cm;
  return 1;
}



/* Function: read_bincm20()
 * 
 * Purpose:  Read binary save files.
 *           
 * Args:     fp     - open file pointer for reading, positioned after magic number
 *           ret_cm - RETURN: model
 *                    
 * Return:   1 on success, 0 on failure                   
 */
static int
read_bincm20(FILE *fp, struct cm_s **ret_cm)
{
  struct cm_s *cm;
  int          nodes;
  int          k;            /* counter for nodes       */

  if (! fread((void *) &(nodes), sizeof(int), 1, fp)) return 0;

				/* now create space for CM. */
  cm = AllocCM(nodes);
  if (cm == NULL) return 0;
				/* everything else is nodes */
  for (k = 0; k < nodes; k++)
    {
				/* type of node */
      if (! fread((void *) &cm->nd[k].type, sizeof(int), 1, fp)) return 0;

				/* indices of child nodes */
      if (! fread((void *) &cm->nd[k].nxt,  sizeof(int), 1, fp)) return 0;
      if (! fread((void *) &cm->nd[k].nxt2, sizeof(int), 1, fp)) return 0;
      
				/* transitions */
      if (! fread((void *) cm->nd[k].tmx, sizeof(double), STATETYPES*STATETYPES, fp)) return 0;

				/* INS emissions */
      if (! fread((void *) cm->nd[k].il_emit, sizeof(double), ALPHASIZE, fp)) return 0;
      if (! fread((void *) cm->nd[k].ir_emit, sizeof(double), ALPHASIZE, fp)) return 0;
      
				/* MAT emissions */
      if (! fread((void *) cm->nd[k].mp_emit, sizeof(double), ALPHASIZE*ALPHASIZE, fp)) return 0;
      if (! fread((void *) cm->nd[k].ml_emit, sizeof(double), ALPHASIZE, fp)) return 0;
      if (! fread((void *) cm->nd[k].mr_emit, sizeof(double), ALPHASIZE, fp)) return 0;
    }
  *ret_cm = cm;
  return 1;
}
