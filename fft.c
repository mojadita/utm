/* $Id: fft.c,v 2.2 2007/07/13 21:55:12 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@HispaLinux.ES>
 * Date: Tue Jun  8 20:35:34 MEST 2004
 *
 * Disclaimer:
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#define IN_FFT_C

/* Standard include files */
#include <sys/types.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "fft.h"
#include "mkroots.h"

/* constants */

/* types */

/* variables */
static char FFT_C_RCSId[]="\n$Id: fft.c,v 2.2 2007/07/13 21:55:12 luis Exp $\n";

/* functions */

/* esta función invierte los bits de una palabra N de n bits (1 <= n <= 32) */
static unsigned reverse(unsigned N, int n)
{
	N = (N & 0x0000ffffU) << 16 | (N & 0xffff0000U) >> 16;
	N = (N & 0x00ff00ffU) <<  8 | (N & 0xff00ff00U) >>  8;
	N = (N & 0x0f0f0f0fU) <<  4 | (N & 0xf0f0f0f0U) >>  4;
	N = (N & 0x33333333U) <<  2 | (N & 0xccccccccU) >>  2;
	N = (N & 0x55555555U) <<  1 | (N & 0xaaaaaaaaU) >>  1;
	return N >> (32 - n);
} /* reverse */

/* suma de complejos */
static void complex_add(complex_t *res, complex_t *a, complex_t *b)
{
	res->x = a->x + b->x;
	res->y = a->y + b->y;
} /* complex_add */

/* resta de complejos */
static void complex_sub(complex_t *res, const complex_t *a, const complex_t *b)
{
	res->x = a->x - b->x;
	res->y = a->y - b->y;
} /* complex_sub */

/* multiplicación de complejos */
static void complex_mult(complex_t *res, const complex_t *a, const complex_t *b)
{
	res->x = a->x*b->x - a->y*b->y;
	res->y = a->x*b->y + a->y*b->x;
} /* complex_mult */

/* esta función aplica la función butterfly de la fft.
 * la función 
 */
static void fft_butterfly(complex_t *p, complex_t *q, const complex_t *w)
{
	complex_t P, Q, aux;

	/* (p, q) = (p + wq, p - wq) */
	complex_mult(&aux, w, q);
	complex_add(&P, p, &aux);
	complex_sub(&Q, p, &aux);

	memcpy(p, &P, sizeof(complex_t));
	memcpy(q, &Q, sizeof(complex_t));
} /* fft_butterfly */

static void fft_rec(fft_t *fft, int n, int N, complex_t *a, int I, complex_t *w)
{

	int i, j;

	if (N <= 1)
		return; /* nothing to do */

#if DEBUG
	printf("fft_rec(fft, %d, %d, a, %d): BEG\n", n, N, I);
#endif

	fft_rec(fft, n-1, N>>1, a,        I<<1, w); /* go recursively to the first half */
	fft_rec(fft, n-1, N>>1, a+(N>>1), I<<1, w); /* ... second half. */

	/* now do the butterflies of this stage */
#if DEBUG
	printf("fft_rec(fft, %d, %d, a, %d): BUTTERFLIES\n", n, N, I);
#endif
	j=0;
	for (i=0; i < (N>>1); i++) {
#if DEBUG
		printf("butterfly([%d], [%d], w[%d])\n", i, i+(N>>1), j);
#endif
		fft_butterfly(a, a+(N>>1), w+j);

		a++;
		j=(j+I)&(fft->N-1);
	} /* for */
#if DEBUG
	printf("fft_rec(fft, %d, %d, a, %d): END\n", n, N, I);
#endif
} /* fft_rec */

/* permuta los datos de entrada N debe ser igual a 2^n */
static void fft_permute(fft_t *fft, complex_t *a)
{
	int i;

	for (i = 0; i < fft->N; i++) {
		int j = reverse(i, fft->n);
		if (i < j) {
			complex_t z;

			/* intercambio a[i] <-> a[j] */
			memcpy(&z, a+i, sizeof z);
			memcpy(a+i, a+j, sizeof z);
			memcpy(a+j, &z, sizeof z);
		} /* if */
	} /* for */

} /* fft_permute */

void fft_init(fft_t *fft, int N)
{
	int i, n;

	for (i = N-1, n=0; i; i>>=1) n++;
#if DEBUG
	printf("fft_init: N=%d; n=%d;\n", N, n);
#endif

	fft->n = n;
	fft->N = N;
	fft->w = calloc(N>>1, sizeof (complex_t));
	fft->W = calloc(N>>1, sizeof (complex_t));

	/* e^(2*Pi*i/N) */
	for (i = 0; i<(N>>1); i++) {
		fft->w[i].x = cos(2.0*M_PI*i/N);
		fft->w[i].y = sin(2.0*M_PI*i/N);
		fft->W[i].x = fft->w[i].x; /* cos(2.0*M_PI*i/N); */
		fft->W[i].y = -fft->w[i].y; /* -sin(2.0*M_PI*i/N); */
	} /* for */
} /* fft_init */

void fft_direct(fft_t *fft, complex_t *a)
{
	int i;

	/* permutamos los datos de entrada */
	fft_permute(fft, a);

	/* vamos con la fft. */
	fft_rec(fft, fft->n, fft->N, a, 1, fft->W);

	/* dividimos a por N */
	for (i = 0; i < fft->N; i++) {
		a[i].x /= (double) fft->N;
		a[i].y /= (double) fft->N;
	} /* for */

} /* fft_direct */

void fft_reverse(fft_t *fft, complex_t *a)
{
	/* permutamos los datos de entrada */
	fft_permute(fft, a);

	/* vamos con la fft */
	fft_rec(fft, fft->n, fft->N, a, 1, fft->w);

} /* fft_reverse */

/* $Id: fft.c,v 2.2 2007/07/13 21:55:12 luis Exp $ */
