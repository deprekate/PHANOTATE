/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* sre_math.c
 * 
 * Portability for and extensions to C math library.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static int sre_reseed = 0;	/* TRUE to reinit sre_random() */
static int sre_randseed = 666;	/* default seed for sre_random()   */

/* Function: Gaussrandom()
 * 
 * Pick a Gaussian-distributed random variable
 * with some mean and standard deviation, and
 * return it.
 * 
 * Based on RANLIB.c public domain implementation.
 * Thanks to the authors, Barry W. Brown and James Lovato,
 * University of Texas, M.D. Anderson Cancer Center, Houston TX.
 * Their implementation is from Ahrens and Dieter, "Extensions 
 * of Forsythe's method for random sampling from the normal
 * distribution", Math. Comput. 27:927-937 (1973).
 *
 * Impenetrability of the code is to be blamed on its FORTRAN/f2c lineage.
 * 
 */
float
Gaussrandom(float mean, float stddev)
{
  static float a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
  };
  static float d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
  };
  static float t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
  };
  static float h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
  };
  static long i;
  static float snorm,u,s,ustar,aa,w,y,tt;

  u = sre_random();
  s = 0.0;
  if(u > 0.5) s = 1.0;
  u += (u-s);
  u = 32.0*u;
  i = (long) (u);
  if(i == 32) i = 31;
  if(i == 0) goto S100;
  /*
   * START CENTER
   */
  ustar = u-(float)i;
  aa = *(a+i-1);
S40:
  if(ustar <= *(t+i-1)) goto S60;
  w = (ustar-*(t+i-1))**(h+i-1);
S50:
  /*
   * EXIT   (BOTH CASES)
   */
  y = aa+w;
  snorm = y;
  if(s == 1.0) snorm = -y;
  return (stddev*snorm + mean);
S60:
  /*
   * CENTER CONTINUED
   */
  u = sre_random();
  w = u*(*(a+i)-aa);
  tt = (0.5*w+aa)*w;
  goto S80;
S70:
  tt = u;
  ustar = sre_random();
S80:
  if(ustar > tt) goto S50;
  u = sre_random();
  if(ustar >= u) goto S70;
  ustar = sre_random();
  goto S40;
S100:
  /*
   * START TAIL
   */
  i = 6;
  aa = *(a+31);
  goto S120;
S110:
  aa += *(d+i-1);
  i += 1;
S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;
S140:
  w = u**(d+i-1);
  tt = (0.5*w+aa)*w;
  goto S160;
S150:
  tt = u;
S160:
  ustar = sre_random();
  if(ustar > tt) goto S50;
  u = sre_random();
  if(ustar >= u) goto S150;
  u = sre_random();
  goto S140;
}

  
/* Function: Linefit()
 * 
 * Purpose:  Given points x[0..N-1] and y[0..N-1], fit to
 *           a straight line y = a + bx.
 *           a, b, and the linear correlation coefficient r
 *           are filled in for return.
 *           
 * Args:     x     - x values of data
 *           y     - y values of data               
 *           N     - number of data points
 *           ret_a - RETURN: intercept
 *           ret_b - RETURN: slope
 *           ret_r - RETURN: correlation coefficient  
 *           
 * Return:   1 on success, 0 on failure.
 */
int          
Linefit(float *x, float *y, int N, float *ret_a, float *ret_b, float *ret_r) 
{				
  float xbar, ybar;
  float sxx, syy, sxy;
  int   i;
  
  /* Calculate averages, xbar and ybar
   */
  xbar = ybar = 0.0;
  for (i = 0; i < N; i++)
    {
      xbar += x[i];
      ybar += y[i];
    }
  xbar /= N;
  ybar /= N;

  sxx = syy = sxy = 0.0;
  for (i = 0; i < N; i++)
    {
      sxx    += (x[i] - xbar) * (x[i] - xbar);
      syy    += (y[i] - ybar) * (y[i] - xbar);
      sxy    += (x[i] - xbar) * (y[i] - ybar);
    }
  *ret_b = sxy / sxx;
  *ret_a = ybar - xbar*(*ret_b);
  *ret_r = sxy / (sqrt(sxx) * sqrt(syy));
  return 1;
}


/* Function: Gammln()
 *
 * Returns the natural log of the gamma function of x.
 * x is > 0.0.  
 *
 * Adapted from a public domain implementation in the
 * NCBI core math library. Thanks to John Spouge and
 * the NCBI. (According to the NCBI, that's Dr. John
 * "Gammas Galore" Spouge to you, pal.)
 */
float
Gammln(float x)
{
  int i;
  double xx, tx;
  double tmp, value;
  static double cof[11] = {
    4.694580336184385e+04,
    -1.560605207784446e+05,
    2.065049568014106e+05,
    -1.388934775095388e+05,
    5.031796415085709e+04,
    -9.601592329182778e+03,
    8.785855930895250e+02,
    -3.155153906098611e+01,
    2.908143421162229e-01,
    -2.319827630494973e-04,
    1.251639670050933e-10
  };
  
  /* Protect against x=0. We see this in Dirichlet code,
   * for terms alpha = 0. This is a severe hack but it is effective
   * and safe. (due to GJM)
   */ 
  if (x <= 0.0) return 999999.; 

  xx       = x - 1.0;
  tx = tmp = xx + 11.0;
  value    = 1.0;
  for (i = 10; i >= 0; i--)	/* sum least significant terms first */
    {
      value += cof[i] / tmp;
      tmp   -= 1.0;
    }
  value  = log(value);
  tx    += 0.5;
  value += 0.918938533 + (xx+0.5)*log(tx) - tx;
  return (float) value;
}


/* Vector operations for doubles and floats.
 * DNorm(), FNorm()   -- normalize a probability vector of length n.
 *                       return 0 if all values were zero.
 * DScale(), FScale() -- multiply all items in vector by scale
 * DSet(), FSet()     -- set all items in vector to value.
 */                      
int
DNorm(double *vec, int n)
{
  int    x;
  double sum;

  sum = 0.0;
  for (x = 0; x < n; x++) sum += vec[x];
  if (sum != 0.0)
    for (x = 0; x < n; x++) vec[x] /= sum;
  else
    { squid_errno = SQERR_DIVZERO; return 0; }
  return 1;
}

int
FNorm(float *vec, int n)
{
  int    x;
  float  sum;

  sum = 0.0;
  for (x = 0; x < n; x++) sum += vec[x];
  if (sum != 0.0)
    for (x = 0; x < n; x++) vec[x] /= sum;
  else
    { squid_errno = SQERR_DIVZERO; return 0; }
  return 1;
}
  
void
DScale(double *vec, int n, double scale)
{
  int x;
  for (x = 0; x < n; x++)
    vec[x] *= scale;
}
void
FScale(float *vec, int n, float scale)
{
  int x;
  for (x = 0; x < n; x++)
    vec[x] *= scale;
}

void
DSet(double *vec, int n, double value)
{
  int x; 
  for (x = 0; x < n; x++)
    vec[x] *= value;
}
void
FSet(float *vec, int n, float value)
{
  int x; 
  for (x = 0; x < n; x++)
    vec[x] *= value;
}

double 
DSum(double *vec, int n)
{
  double sum = 0.;
  int    x;
  for (x = 0; x < n; x++)
    sum += vec[x];
  return sum;
}
float 
FSum(float *vec, int n)
{
  float sum = 0.;
  int   x;
  for (x = 0; x < n; x++)
    sum += vec[x];
  return sum;
}

/* Function: sre_random()
 * 
 * Purpose:  Return a uniform deviate from 0.0 to 1.0.
 *           sre_randseed is a static variable, set
 *           by sre_srandom(). sre_reseed is a static flag
 *           raised by sre_srandom(), saying that we need
 *           to re-initialize.
 *           [0.0 <= x < 1.0]
 *           
 *           Uses a simple linear congruential generator with
 *           period 2^28. Based on discussion in Robert Sedgewick's 
 *           _Algorithms in C_, Addison-Wesley, 1990.
 *
 *           Requires that long int's have at least 32 bits.
 *
 * Reliable and portable, but slow. Benchmarks on wol,
 * using IRIX cc and IRIX C library rand() and random():
 *     sre_random():    0.8 usec/call
 *     random():        0.3 usec/call
 *     rand():          0.3 usec/call
 */
#define RANGE 268435456		/* 2^28        */
#define DIV   16384		/* sqrt(RANGE) */
#define MULT  72530821		/* my/Cathy's birthdays, x21, x even (Knuth)*/
float 
sre_random(void)
{
  static long  rnd;
  static int   firsttime = 1;
  long         high1, low1;
  long         high2, low2; 

  if (sre_reseed || firsttime) 
    {
      sre_reseed = firsttime = 0;
      if (sre_randseed <= 0) sre_randseed = 666; /* seeds of zero break me */
      high1 = sre_randseed / DIV;  low1  = sre_randseed % DIV;
      high2 = MULT / DIV;          low2  = MULT % DIV;
      rnd = (((high2*low1 + high1*low2) % DIV)*DIV + low1*low2) % RANGE;
    }
   high1 = rnd / DIV;  low1  = rnd % DIV;
   high2 = MULT / DIV; low2  = MULT % DIV;
   rnd = (((high2*low1 + high1*low2) % DIV)*DIV + low1*low2) % RANGE;

   return ((float) rnd / (float) RANGE);  
}
#undef RANGE
#undef DIV
#undef MULT


/* Function: sre_srandom()
 * 
 * Purpose:  Initialize with a random seed. Seed can be
 *           any integer.
 */
void
sre_srandom(int seed)
{
  if (seed < 0) seed = -1 * seed;
  sre_reseed   = 1;
  sre_randseed = seed;
}


/* Functions: DChoose(), FChoose()
 *
 * Purpose:   Make a random choice from a normalized distribution.
 *            DChoose() is for double-precision vectors;
 *            FChoose() is for single-precision float vectors.
 *            Returns the number of the choice.
 */
int
DChoose(double *p, int N)
{
  double roll;                  /* random fraction */
  double sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = sre_random();
  sum     = 0.0;
  for (i = 0; i < N; i++)
    {
      sum += p[i];
      if (roll < sum) return i;
    }
  return (int) (sre_random() * N);   /* bulletproof */
}
int
FChoose(float *p, int N)
{
  float roll;                   /* random fraction */
  float sum;			/* integrated prob */
  int   i;                      /* counter over the probs */

  roll    = sre_random();
  sum     = 0.0;
  for (i = 0; i < N; i++)
    {
      sum += p[i];
      if (roll < sum) return i;
    }
  return (int) (sre_random() * N);    /* bulletproof */
}

/* Functions: DLogSum(), FLogSum()
 * 
 * Calculate the sum of a log vector
 * *in normal space*, and return the log of the sum.
 */
double 
DLogSum(double *logp, int n)
{
  int x;
  double max = -1.0e30;
  double sum;
  
  for (x = 0; x < n; x++)
    if (logp[x] > max) max = logp[x];
  sum = 0.0;
  for (x = 0; x < n; x++)
    if (logp[x] > max - 50.)
      sum += exp(logp[x] - max);
  sum = log(sum) + max;
  return sum;
}
float
FLogSum(float *logp, int n)
{
  int x;
  float max = -1.0e30;
  float sum;
  
  for (x = 0; x < n; x++)
    if (logp[x] > max) max = logp[x];
  sum = 0.0;
  for (x = 0; x < n; x++)
    if (logp[x] > max - 50.)
      sum += exp(logp[x] - max);
  sum = log(sum) + max;
  return sum;
}
