#include <stdio.h>
#include <math.h>

#include "fft.h"
#define NTERM	16

double f(double x)
{
	return cos(x);
} /* f */

double mod(complex_t *x)
{
	return sqrt(x->x*x->x + x->y*x->y);
} /* mod */

main()
{
	complex_t X[NTERM];
	fft_t FFT;
	int i;

	fft_init(&FFT, NTERM);
	for (i = 0; i < NTERM; i++){
		X[i].y = 0.0; X[i].x = f(2.0*M_PI*i/NTERM);
		printf("x[%d] = %0.10lg\n", i, X[i].x);
	} /* for */
	fft_direct(&FFT, X);
	for (i = 0; i < NTERM; i++) {
		printf("X[%d] = <%0.10lg, %0.10lg>, |X[%d]| = %0.10lg, arg(X[%d]) = %0.10lg\n",
			i, X[i].x, X[i].y,
			i, mod(&X[i]),
			i, atan2(X[i].y,X[i].x));
	} /* for */
} /* main */
