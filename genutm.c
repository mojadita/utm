/* $Id: genutm.c,v 1.1 1998/05/10 16:59:44 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Sun May 10 15:25:27 MET DST 1998
 */

#define IN_UTM_C

/* Standard include files */
#include <sys/types.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

/* eccentricity of earth (WGS84) */
#define E2 0.0067226700223332915
#define A  6378388.0
#define K0 0.9996

double n(double l)
{
  double sl = sin(l);
  return 1.0/sqrt(1-E2*sl*sl);
}

double m(double l)
{
  double nl = n(l);
  return (1-E2)*nl*nl*nl;
}

double simpson (double(*f)(double),double a, double b, int n)
{
  double t = a;
  double dt = (b - a) / (double) n;
  double f1 = f(a);
  double f2, f3;
  double acum = 0.0;
  int i;

  for (i = 0; i < n; i++) {
  	f2 = f(t + dt/2.0);
	f3 = f(t + dt);
	acum += (f1 + 4.0*f2 + f3)*dt/6.0;
	f1 = f3;
	t += dt;
  }
  return acum;
}

/* This function iterates with simpson until the diferences begin to
 * grow due to truncating/rounding errors.  It supposes that the
 * different iterations will give monothonic difference values. */
double simpson_e (double(*f)(double),double a, double b)
{
  double prev, next;
  double preverr, nexterr;
  int n = 1;
  prev = simpson(f, a, b, n);
  preverr = MAXDOUBLE;
  for (;;) {
    n <<= 1;
    next = simpson(f, a, b, n);
    printf ("prev[%d]=%20.17lg; next[%d]=%20.17lg\n",n>>1,prev,n, next);
    nexterr = fabs (next - prev);
    if (nexterr >= preverr)
      break;
    prev = next;
    preverr = nexterr;
  }
  return prev;
}

/* main program */
int main (int argc, char **argv)
{
	char linea [1000];
	double l, err;

	for (;;) {
	  printf ("l> ");
	  gets(linea);
	  sscanf (linea, "%lg", &l);
	  l *= PI/180.0;
	  printf ("l = %20.17g\n", l);
	  printf ("M(l) = %20.17g\n", A*m(l));
	  printf ("N(l) = %20.17g\n", A*n(l));
	  printf ("Simpson(M,0,l,1e-8/K0/A) = %20.17g\n",
	    K0*A*simpson_e(m, 0.0,l, 1e-8/K0/A));
	}
}

/* $Id: genutm.c,v 1.1 1998/05/10 16:59:44 luis Exp $ */
