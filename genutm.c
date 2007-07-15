/* $Id: genutm.c,v 2.15 2007/07/15 19:42:58 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Sun May 10 15:25:27 MET DST 1998
 * $Log: genutm.c,v $
 * Revision 2.15  2007/07/15 19:42:58  luis
 * * Modificaciones para integrar FFT en el cÃ¡lculo de las constantes.
 *
 * Revision 2.14  2007-07-13 21:55:12  luis
 * * Introduciendo la sustituciÃ³n del cÃ¡lculo de las series de Fourier por
 *   el cÃ¡lculo de los coeficientes a partir de la FFT.
 *
 * Revision 2.13  2007-07-13 20:58:28  luis
 * * Estamos incluyendo soporte para calcular los datos por fft en lugar de
 *   integrar numericamente por simpson.
 *
 * Revision 2.12  2007-07-13 19:57:54  luis
 * * Incluidas las funciones para el cÃ¡lculo de la FTT del mÃ³dulo fft.
 *   AÃºn no estÃ¡n integradas.
 *
 * Revision 2.11  2005/10/17 19:39:47  luis
 * Recuperado genutm.c
 *
 * Revision 2.9  2002/09/23 06:14:17  luis
 * Modified to support variable number of GEO_NPOT and GEO_NTERM.
 *
 * Revision 2.8  2002/09/17 19:58:27  luis
 * Added more precision.
 *
 * Revision 2.7  2002/09/06 00:12:11  luis
 * Añadidos utm_ini.h para que genutm pueda calcular por tabla los parámetros y
 * utmcalc.c para los cálculos a partir de los parámetros calculados según la
 * estructura utmparam.
 *
 * Revision 2.6  1998/08/24 13:04:47  luis
 * Changes in nomenclature of Reference Ellipsoid Names to agree with
 * WGS 1984 document.
 *
 * Revision 2.5  1998/08/06 12:07:25  luis
 * Found error in calculus of dQ2Lat to find the increment of latitude
 * from the increment in isometric latitude.
 *
 * Revision 2.4  1998/08/05 19:10:18  luis
 * Complete utm transformation(direct and inverse), but there must be an
 * error as there are some errors when going far from the central meridian,
 * giving up to 6-10meters of diference between the direct transform(exact
 * to the milimeter) and the reverse.
 *
 * Revision 2.3  1998/08/03 22:11:55  luis
 * Calculo de Ateb completo.
 *
 * Revision 2.2  1998/06/20 20:20:55  luis
 * *** empty log message ***
 *
 * Revision 2.1  1998/05/14 16:54:19  luis
 * Revisión con cálculo de coordenadas UTM a partir de geodésicas y
 * determinación del recuadro cienkilométrico y la zona.
 *
 * Revision 2.0  1998/05/12 18:53:09  luis
 * Versión directa completa con cálculo de coordenadas UTM a partir
 * de geodésicas.
 *
 * Revision 1.2  1998/05/11 18:45:00  luis
 * Cálculo directo completo.
 *
 */

#define IN_GENUTM_C

/* Standard include files */
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utm.h"
#include "utm_ini.h"

#include "fft.h"

/* Number of iterations in Simpson's numerical integration */
#ifndef Niter
#define Niter 1024
#endif

struct utmparam *sg = wgs84table;

struct utmparam *lookup(char *name)
{
	struct utmparam *p;
	for (p = wgs84table; p->name; p++)
		if(!strcmp(p->name,name))
			break;

	return p->name ? p : NULL;
} /* lookup */

/* N equatorial radius at point of latitude l given in terms of A */
double n(int dumb, double l)
{
  double sl = sin(l);
  return 1.0 / sqrt(1 - sg->e2 * sl*sl);
} /* n */

/* M meridianal radius at point of latitude l given in terms of A */
double m(int dumb, double l)
{
  double nl = n(0, l);
  return(1.0 - sg->e2) * nl*nl*nl;
} /* m */

/* Simpson's integral of function f, between a and b, subdivided in
 * n subintervals */
double simpson(double(*f)(double),double a, double b, int n)
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
  } /* for */

  return acum;
} /* simpson */

/* Auxiliary functions to calculate Fourier series. */
static double(*Fn)(int, double);
static double(*SC)(double);
static double I;
static int O;

double F_Fourier(double x)
{
	return Fn(O, x)*SC(I*x);
} /* F_Fourier */

double C_Fourier_sin(double(*f)(int, double), int ord, int i, int n)
{
  Fn = f; SC = sin; I = i; O = ord; 

  return simpson(F_Fourier, 0.0, 2.0*M_PI, n);
} /* C_Fourier_sin */

double C_Fourier_cos(double(*f)(int, double), int ord, int i, int n)
{
  double result;

  Fn = f; SC = cos; I = i; O = ord; 

  result = simpson(F_Fourier, 0.0, 2*M_PI, n);
  if(i == 0) result /= 2.0;

  return result;
} /* C_Fourier_cos */

static void myFFT(double (*f)(int n, double x), complex_t *v)
{
	static fft_t FFT;
	static int notused = TRUE;
	complex_t fft_array[2*GEO_NTERM];

	if (notused) {
		notused = FALSE;
		fft_init(&FFT, 2*GEO_NTERM);
	} /* if */

	for (i = 0; i < 2*GEO_NTERM; i++) {
		fft_array[i].x = f(0, M_PI*i/GEO_NTERM);
		fft_array[i].y = 0.0;
	} /* for */

	fft_forward(&FFT, fft_array);
} /* myFFT */

/************************************************************************
 ******************* COMIENZO DE LAS FUNCIONES CALCULADAS ***************
 ************************************************************************/

/********* TRANSFORMACIÓN GEODÉSICAS A UTM *****************/

double Beta(double x)  /* Beta */
{
  int i;
  double res;

  res = x*sg->BetaPhi;
  for (i = 1; i < GEO_NTERM; i++)
    res += sg->Beta[i] * sin(i*x);

  return res;
} /* Beta */

double preA(int ord, double x)
{
  int i;
  double res = 0.0;

  if (ord <= 1) return n(0, x)*cos(x);

  for (i = 0; i < GEO_NTERM; i++) {
  	switch (ord & 1) {
	case 0: res -= sg->A[ord-1][i] * i * sin(i*x); break;
	case 1: res += sg->A[ord-1][i] * i * cos(i*x); break;
	} /* switch */
  } /* for */

  return n(0, x)/m(0, x)*cos(x)*res/ord;
} /* preA */

double A(int ord, double x)
{
  int i;
  double res;

  res = 0.0;
  for (i = 0; i < GEO_NTERM; i++) {
    switch (ord & 1) {
	case 0: res += sg->A[ord][i] * sin(i*x); break;
	case 1: res += sg->A[ord][i] * cos(i*x); break;
	} /* switch */
  } /* for */

  return res;
} /* A */

/************** TRANSFORMACIÓN UTM -> GEODÉSICAS ****************/

/* FUNCIÓN ATEB, INVERSA DE LA FUNCIÓN BETA */
double derAteb(int, double);

double preAteb(int ord, double x)
{
  if (ord <= 1) return 1.0/m(0, x);
  return derAteb(ord-1, x) / m(0, x) / ord;
} /* preAteb */

double Ateb(int ord, double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++) {
  	switch (ord & 1) {
    case 0: res += sg->Ateb[ord][i] * sin(i*x); break;
    case 1: res += sg->Ateb[ord][i] * cos(i*x); break;
	} /* switch */
  } /* for */

  return res;
} /* Ateb */

double derAteb(int ord, double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++) {
	switch (ord & 1) {
	case 0: res += sg->Ateb[ord][i] * i * cos(i*x); break;
    case 1: res -= sg->Ateb[ord][i] * i * sin(i*x); break;
	} /* switch */
  } /* for */

  return res;
} /* derAteb */

/*******************/

double preF(int, double);
double F(int, double);
double B(int, double);
double dF(int, double);

double preF(int ord, double x)
{
  if (ord <= 1) return 1.0/n(0, x);
  return (dF(ord-1, x)*cos(x) + (ord-1)*F(ord-1, x)*sin(x))/m(0, x)/ord;
} /* preF1 */

double F(int ord, double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++) {
  	switch (ord & 1) {
	case 0: res += sg->F[ord][i] * sin(i*x); break;
    case 1: res += sg->F[ord][i] * cos(i*x); break;
	} /* switch */
  } /* for */

  return res;
} /* F1 */

double B(int ord, double x)
{
  return F(ord, x) / pow(cos(x), (double)ord);
} /* B1 */

double dF(int ord, double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++) {
    switch (ord & 1) {
	case 0: res += i * sg->F[ord][i] * cos(i*x); break;
    case 1: res -= i * sg->F[ord][i] * sin(i*x); break;
	} /* switch */
  } /* for */

  return res;
} /* dF1 */


/************ INCREMENTO DE Q -> INCREMENTO DE LATITUD **************/
double predQ2Lat(int, double);
double dQ2Lat(int, double);
double derdQ2Lat(int, double);

double predQ2Lat(int ord, double phi)
{
  if (ord <= 1) return n(0, phi)/m(0, phi)*cos(phi);
  return n(0, phi)/m(0, phi) * cos(phi) * derdQ2Lat(ord-1, phi) / ord;
} /* predQ2Lat */

double dQ2Lat(int ord, double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++) {
    switch (ord & 1) {
	case 0: res += sg->dQ2Lat[ord][i] * sin(i*phi); break;
    case 1: res += sg->dQ2Lat[ord][i] * cos(i*phi); break;
	} /* switch */
  } /* for */

  return res;
} /* dQ2Lat */

double derdQ2Lat(int ord, double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++) {
  	switch (ord & 1) {
	case 0: res += sg->dQ2Lat[ord][i] * i * cos(i*phi); break;
	case 1: res -= sg->dQ2Lat[ord][i] * i * sin(i*phi); break;
	} /* switch */
  } /* for */

  return res;
} /* derdQ2Lat */

/*******************************************************************/
/*******************************************************************/
/* main program */
int main(int argc, char **argv)
{
	char linea [1000];
	double l, L, err;
	int i, opt, ord;
	extern char *optarg;
	fft_t FFT;
	fft_init (&FFT, 2*GEO_NTERM);

  while((opt = getopt(argc, argv, "g:")) != EOF) {
    switch(opt){
    case 'g': sg=lookup(optarg);
		if(!sg) {
			fprintf(stderr, "Sistema Geodésico desconocido(%s)\n", optarg);
			exit(1);
		}
		break;
    default:
      fprintf(stderr,
       "usage: genutm [ -g geodsystem ]\n");
      exit(1);
    }
  }

  printf("divert(-1)\n");
  printf("define(ELLIPSOID, ``%s'')\n", sg->dsc);
  printf("define(NAME, ``%s'')\n", sg->name);
  printf("define(Aaxis,%0.17lG)\n", sg->a);
  printf("define(Baxis,%0.17lG)\n", sg->a*sqrt(1-sg->e2));
  printf("define(E2,%0.17lG)\n", sg->e2);
  printf("define(ak0,%0.17lG)\n", sg->a*GEO_K0);
  printf("define(GEO_NTERM,%d)\n", GEO_NTERM);
  printf("define(GEO_NPOT,%d)\n", GEO_NPOT);

  { /* M (FFT) */
  	complex_t M[2*GEO_NTERM];
	for (i = 0; i < 2*GEO_NTERM; i++) {
		M[i].x = m(0,M_PI*i/GEO_NTERM); M[i].y = 0.0;
	} /* for */
	fft_direct(&FFT, M);
	for (i = 0; i < GEO_NTERM; i++) {
		printf("define(M_%d_FFT,%0.17lG)\n", i, (i == 0) ? M[i].x : 2.0*M[i].x);
	} /* for */
  } /* bloque */
  /* M */
  for (i = 0; i < GEO_NTERM; i++) {
    sg->M[i] = (i & 1) ? 0.0 : C_Fourier_cos(m, 0, i, Niter)/M_PI;
    printf("define(M_%d,%0.17lG)\n", i, sg->M[i]);
  }

  { /* N (FFT) */
  	complex_t M[2*GEO_NTERM];
	for (i = 0; i < 2*GEO_NTERM; i++) {
		M[i].x = n(0,M_PI*i/GEO_NTERM); M[i].y = 0.0;
	} /* for */
	fft_direct(&FFT, M);
	for (i = 0; i < GEO_NTERM; i++) {
		printf("define(N_%d_FFT,%0.17lG)\n", i, (i == 0) ? M[i].x : 2.0*M[i].x);
	} /* for */
  } /* bloque */
  /* N */
  for (i = 0; i < GEO_NTERM; i++) {
    sg->N[i] = (i & 1) ? 0.0 : C_Fourier_cos(n, 0, i, Niter)/M_PI;
    printf("define(N_%d,%0.17lG)\n", i, sg->N[i]);
  }

  /* BetaPhi */
  sg->BetaPhi = sg->M[0];
  printf("define(BetaPhi,%0.17lG)\n", sg->BetaPhi);
  printf("define(BetaPhi_deg,%0.17lG)\n", sg->BetaPhi/180.0*M_PI);

  /* Beta*/
  sg->Beta[0] = 0.0;
  printf("define(Beta_0,%0.17lG)\n", 0.0);
  for (i = 1; i < GEO_NTERM; i++) {
    sg->Beta[i] = sg->M[i]/i;
    printf("define(Beta_%d,%0.17lG)\n", i, sg->Beta[i]);
  } /* for */

  /* A */
  for (i = 0; i < GEO_NTERM; i++) {
  	printf("define(A_0_%d,%0.17lG)\n", i, 0.0);
  } /* for */
  for (ord = 1; ord < GEO_NPOT; ord++) {
	for (i = 0; i < GEO_NTERM; i++) {
  		switch (ord & 1) {
		case 1: sg->A[ord][i] = i & 1
				? C_Fourier_cos(preA, ord, i, Niter)/M_PI
				: 0.0;
			break;
		case 0: sg->A[ord][i] = i & 1
				? 0.0
				: C_Fourier_sin(preA, ord, i, Niter)/M_PI;
			break;
		} /* switch */
		printf("define(A_%d_%d,%0.17lG)\n", ord, i, sg->A[ord][i]);
	} /* for */
  } /* for */

  /* Ateb */
  for (i = 0; i < GEO_NTERM; i++) {
  	printf("define(Ateb_0_%d,%0.17lG)\n", i, 0.0);
  	printf("define(Ateb_deg_0_%d,%0.17lG)\n", i, 0.0);
  } /* for */
  for (ord = 1; ord < GEO_NPOT; ord++) {
	for (i = 0; i < GEO_NTERM; i++) {
    	switch (ord & 1) {
		case 1: sg->Ateb[ord][i] = i & 1
				? 0.0
				: C_Fourier_cos(preAteb, ord, i, Niter)/M_PI;
			break;
		case 0: sg->Ateb[ord][i] = i & 1
				? 0.0
				: C_Fourier_sin(preAteb, ord, i, Niter)/M_PI;
			break;
		} /* switch */
		printf("define(Ateb_%d_%d,%0.17lG)\n", ord, i, sg->Ateb[ord][i]);
		printf("define(Ateb_deg_%d_%d,%0.17lG)\n", ord, i, sg->Ateb[ord][i]*180.0/M_PI);
	} /* for */
  } /* for */

  /* BetaPI */
  sg->BetaPI = geo_Beta(sg, M_PI)/sg->a;
  printf("define(BetaPI,%0.17lG)\n", sg->BetaPI);

  /* F */
  for (i = 0; i < GEO_NTERM; i++) {
  	printf("define(F_0_%d,%0.17lG)\n", i, 0.0);
  	printf("define(F_deg_0_%d,%0.17lG)\n", i, 0.0);
  } /* for */
  for (ord = 1; ord < GEO_NPOT; ord++) {
	for (i = 0; i < GEO_NTERM; i++) {
  		switch (ord & 1) {
		case 1: sg->F[ord][i] = i & 1
				? 0.0
				: C_Fourier_cos(preF, ord, i, Niter)/M_PI;
			break;
		case 0: sg->F[ord][i] = i & 1
				? C_Fourier_sin(preF, ord, i, Niter)/M_PI
				: 0.0;
			break;
		} /* switch */
		printf("define(F_%d_%d,%0.17lG)\n", ord, i, sg->F[ord][i]);
		printf("define(F_deg_%d_%d,%0.17lG)\n", ord, i, sg->F[ord][i]*180.0/M_PI);
	} /* for */
  } /* for */

  /* dQ2Lat */
  for (i = 0; i < GEO_NTERM; i++) {
  	printf("define(dQ2Lat_0_%d,%0.17lG)\n", i, 0.0);
  	printf("define(dQ2Lat_deg_0_%d,%0.17lG)\n", i, 0.0);
  } /* for */
  for (ord = 1; ord < GEO_NPOT; ord++) {
  	for (i = 0; i < GEO_NTERM; i++) {
		switch (ord & 1) {
		case 1: sg->dQ2Lat[ord][i] = i & 1
				? C_Fourier_cos(predQ2Lat, ord, i, Niter)/M_PI
				: 0.0;
			break;
		case 0: sg->dQ2Lat[ord][i] = i & 1
				? 0.0
				: C_Fourier_sin(predQ2Lat, ord, i, Niter)/M_PI;
			break;
		} /* switch */
		printf("define(dQ2Lat_%d_%d,%0.17lG)\n", ord, i, sg->dQ2Lat[ord][i]);
		printf("define(dQ2Lat_deg_%d_%d,%0.17lG)\n", ord, i, sg->dQ2Lat[ord][i]*180.0/M_PI);
	} /* for */
  } /* for */

  printf("divert(0)dnl\n");
  exit(0);
} /* main */

/* $Id: genutm.c,v 2.15 2007/07/15 19:42:58 luis Exp $ */
