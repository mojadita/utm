/* $Id: genutm.c,v 2.6 1998/08/24 13:04:47 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Sun May 10 15:25:27 MET DST 1998
 * $Log: genutm.c,v $
 * Revision 2.6  1998/08/24 13:04:47  luis
 * Changes in nomenclature of Reference Ellipsoid Names to agree with
 * WGS 1984 document.
 *
 * Revision 2.5  1998/08/06 12:07:25  luis
 * Found error in calculus of dQ2Lat to find the increment of latitude
 * from the increment in isometric latitude.
 *
 * Revision 2.4  1998/08/05 19:10:18  luis
 * Complete utm transformation (direct and inverse), but there must be an
 * error as there are some errors when going far from the central meridian,
 * giving up to 6-10meters of diference between the direct transform (exact
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

/* eccentricity of earth (EURO50) squared */
#define ELLIPSOID	"International 1924"  /* */
#define E2 0.0067226700223332915        /* International 1924 */
#define A  6378388.000                  /* International 1924 */
/*#define ELLIPSOID	"WGS 1984"  /* */
/*#define E2 0.00669437999014           /* WGS 1984 */
/*#define A 6378137.000                 /* WGS 1984 */

/* UTM reduction constant */
#define K0 0.9996

/* Number of iterations in Simpson's numerical integration */
#define N 1024

/* Number of terms used in fourier series */
#define NTERM 8

double e2 = E2;
double a = A;
double k0 = K0;
char *desc = ELLIPSOID;

/* N equatorial radius at point of latitude l given in terms of A */
double n(double l)
{
  double sl = sin(l);
  return 1.0/sqrt(1-e2*sl*sl);
}

/* M meridianal radius at point of latitude l given in terms of A */
double m(double l)
{
  double nl = n(l);
  return (1-e2)*nl*nl*nl;
}

/* Simpson's integral of function f, between a and b, subdivided in
 * n subintervals */
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

/* Auxiliary functions to calculate Fourier series. */
double (*F)(double);
double (*SC)(double);
double I;

double F_Fourier (double x)
{
	return F(x)*SC(I*x);
}

double C_Fourier_sin (double(*f)(double), int i, int n)
{
  F = f; SC = sin; I = i;
  return simpson (F_Fourier, 0.0, 2.0*M_PI, n);
}

double C_Fourier_cos (double(*f)(double), int i, int n)
{
  double result;
  F = f; SC = cos; I = i;
  result = simpson (F_Fourier, 0.0, 2*M_PI, n);
  if (i == 0) result /= 2.0;
  return result;
}


/************************************************************************
 ******************* COMIENZO DE LAS FUNCIONES CALCULADAS ***************
 ************************************************************************/

/********* TRANSFORMACIÓN GEODÉSICAS A UTM *****************/

double Mcos[NTERM];  /* To determine M */

double Ncos[NTERM]; /* To determine N */

double Betasin[NTERM];

double Beta(double x)  /* Beta */
{
  int i;
  double res;
  res = x*Betasin[0];
  for (i = 1; i < NTERM; i++)
    res += Betasin[i] * sin(i*x);
  return res;
}

double preA1(double x) /* N * cos(l) */
{ return n(x)*cos(x);
}

double A1cos[NTERM];

double A1(double x)
{
  int i;
  double res;
  res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A1cos[i] * cos(i*x);
  return res;
}

double preA2 (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A1cos[i] * -i * sin (i*x);
  return n(x)/m(x)*cos(x)*res/2.0;
}

double A2sin[NTERM];

double A2 (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A2sin[i] * sin (i*x);
  return res;
}

double preA3 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A2sin[i] * i * cos (i*x);
  return n(x)*cos(x)/m(x)*res/3.0;
}

double A3cos[NTERM];

double A3 (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A3cos[i] * cos (i*x);
  return res;
}

double preA4 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A3cos[i] * -i * sin (i*x);
  return n(x)*cos(x)/m(x)*res/4.0;
}

double A4sin[NTERM];

double A4 (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A4sin[i] * sin (i*x);
  return res;
}

double preA5 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A4sin[i] * i * cos (i*x);
  return n(x)*cos(x)/m(x)*res/5.0;
}

double A5cos[NTERM];

double A5 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A5cos[i] * cos (i*x);
  return res;
}

double preA6 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A5cos[i] * -i * sin (i*x);
  return n(x)*cos(x)/m(x)*res/6.0;
}

double A6sin[NTERM];

double A6 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A6sin[i] * sin (i*x);
  return res;
}

/************** TRANSFORMACIÓN UTM -> GEODÉSICAS ****************/

/* FUNCIÓN ATEB, INVERSA DE LA FUNCIÓN BETA */
double preAteb1(double x)
{
  return 1.0/m(x);
}

double Ateb1cos[NTERM];

double Ateb1(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb1cos[i] * cos (i*x);
  return res;
}

double derAteb1(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= Ateb1cos[i] * i * sin (i*x);
  return res;
}

double preAteb2(double x)
{
  return derAteb1(x) / m(x) / 2.0;
}

double Ateb2sin[NTERM];

double Ateb2(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb2sin[i] * sin (i*x);
  return res;
}

double derAteb2(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb2sin[i] * i * cos (i*x);
  return res;
}

double preAteb3(double x)
{
  return derAteb2(x) / m(x) / 3.0;
}

double Ateb3cos[NTERM];

double Ateb3(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb3cos[i] * cos (i*x);
  return res;
}

double derAteb3(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= Ateb3cos[i] * i * sin (i*x);
  return res;
}

double preAteb4(double x)
{
  return derAteb3(x) / m(x) / 4.0;
}

double Ateb4sin[NTERM];

double Ateb4(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb4sin[i] * sin (i*x);
  return res;
}

double derAteb4(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb4sin[i] * i * cos (i*x);
  return res;
}

double preAteb5(double x)
{
  return derAteb4(x) / m(x) / 5.0;
}

double Ateb5cos[NTERM];

double Ateb5(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb5cos[i] * cos (i*x);
  return res;
}

double derAteb5(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= Ateb5cos[i] * i * sin (i*x);
  return res;
}

double preAteb6(double x)
{
  return derAteb5(x) / m(x) / 6.0;
}

double Ateb6sin[NTERM];

double Ateb6(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb6sin[i] * sin (i*x);
  return res;
}

double derAteb6(double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += Ateb6sin[i] * i * cos (i*x);
  return res;
}

double preAteb7(double x)
{
  return derAteb6(x) / m(x) / 7.0;
}

double BetaPI;

double Ateb(double x)
{
	double phi0 = x / BetaPI * M_PI;
	double x0 = Beta(phi0);
	double dx = x - x0;
	return
	    phi0
	  + Ateb1(phi0) * dx
	  + Ateb2(phi0) * dx * dx
	  + Ateb3(phi0) * pow(dx, 3.0)
	  + Ateb4(phi0) * pow(dx, 4.0)
	  + Ateb5(phi0) * pow(dx, 5.0)
	  + Ateb6(phi0) * pow(dx, 6.0);
}


/*******************/

double preF1 (double x)
{
  return 1.0/n(x);
}

double F1cos[NTERM];

double F1 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F1cos[i] * cos (i*x);
  return res;
}

double B1 (double x)
{ return F1(x) / cos(x);
}

double B1_(double x)
{ return 1.0 / n(x) / cos(x);
}

double dF1 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= i*F1cos[i] * sin (i*x);
    return res;
}

double preF2 (double x)
{
  return (dF1(x)*cos(x) + F1(x)*sin(x))/m(x)/2.0;
}

double F2sin[NTERM];

double F2 (double x)
{  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F2sin[i] * sin (i*x);
  return res;
}

double B2 (double x)
{ return F2(x)/pow(cos(x), 2.0);
}

double B2_(double x)
{ return sin(x)/2.0/pow(n(x)*cos(x), 2.0);
}

double dF2 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += i*F2sin[i] * cos (i*x);
  return res;
}

double preF3 (double x)
{ return (dF2(x) * cos(x) + 2.0 * F2(x) * sin(x))/3.0/m(x);
}

double F3cos [NTERM];

double F3 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F3cos [i] * cos(i*x);
  return res;
}

double B3 (double x)
{ return F3(x)/pow(cos(x), 3.0);
}

double dF3 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= i*F3cos [i] * sin(i*x);
  return res;
}

double preF4 (double x)
{ return (dF3(x)*cos(x) + 3.0*F3(x)*sin(x))/4.0/m(x);
}

double F4sin[NTERM];

double F4 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F4sin[i] * sin (i*x);
  return res;
}

double B4 (double x)
{ return F4(x)/pow(cos(x), 4.0); }

double dF4 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += i*F4sin [i] * cos(i*x);
  return res;
}

double preF5 (double x)
{ return (dF4(x)*cos(x) + 4.0*F4(x)*sin(x))/5.0/m(x);
}

double F5cos[NTERM];

double F5 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F5cos[i] * cos (i*x);
  return res;
}

double B5 (double x)
{ return F5(x)/pow(cos(x), 5.0); }

double dF5 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= i*F5cos [i] * sin(i*x);
  return res;
}

double preF6 (double x)
{ return (dF5(x)*cos(x) + 5.0*F5(x)*sin(x))/6.0/m(x);
}

double F6sin[NTERM];

double F6 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F6sin[i] * sin (i*x);
  return res;
}

double B6 (double x)
{ return F6(x)/pow(cos(x), 6.0); }

/************ INCREMENTO DE Q -> INCREMENTO DE LATITUD **************/
double predQ2Lat1 (double phi)
{ return n(phi)/m(phi)*cos(phi);
}

double dQ2Lat1cos[NTERM];

double dQ2Lat1 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat1cos[i] * cos(i*phi);
  return res;
}

double derdQ2Lat1 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res -= dQ2Lat1cos[i] * i * sin(i*phi);
  return res;
}

double predQ2Lat2 (double phi)
{ return n(phi)/m(phi) * cos(phi) * derdQ2Lat1(phi) / 2.0;
}

double dQ2Lat2sin[NTERM];

double dQ2Lat2 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat2sin[i] * sin(i*phi);
  return res;
}

double derdQ2Lat2 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat2sin[i] * i * cos(i*phi);
  return res;
}

double predQ2Lat3 (double phi)
{ return n(phi)/m(phi) * cos(phi) * derdQ2Lat2(phi) / 3.0;
}

double dQ2Lat3cos[NTERM];

double dQ2Lat3 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat3cos[i] * cos(i*phi);
  return res;
}

double derdQ2Lat3 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res -= dQ2Lat3cos[i] * i * sin(i*phi);
  return res;
}

double predQ2Lat4 (double phi)
{ return n(phi)/m(phi) * cos(phi) * derdQ2Lat3(phi) / 4.0;
}

double dQ2Lat4sin[NTERM];

double dQ2Lat4 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat4sin[i] * sin(i*phi);
  return res;
}

double derdQ2Lat4 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat4sin[i] * i * cos(i*phi);
  return res;
}

double predQ2Lat5 (double phi)
{ return n(phi)/m(phi) * cos(phi) * derdQ2Lat4(phi) / 5.0;
}

double dQ2Lat5cos[NTERM];

double dQ2Lat5 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat5cos[i] * cos(i*phi);
  return res;
}

double derdQ2Lat5 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res -= dQ2Lat5cos[i] * i * sin(i*phi);
  return res;
}

double predQ2Lat6 (double phi)
{ return n(phi)/m(phi) * cos(phi) * derdQ2Lat5(phi) / 6.0;
}

double dQ2Lat6sin[NTERM];

double dQ2Lat6 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat6sin[i] * sin(i*phi);
  return res;
}

double derdQ2Lat6 (double phi)
{ int i;
  double res = 0.0;
  for (i=0; i < NTERM; i++)
    res += dQ2Lat6sin[i] * i * cos(i*phi);
  return res;
}

double predQ2Lat7 (double phi)
{ return n(phi)/m(phi) * cos(phi) * derdQ2Lat6(phi) / 7.0;
}

/* main program */
int main (int argc, char **argv)
{
	char linea [1000];
	double l, L, err;
	int i, opt;
	extern char *optarg;

  while ((opt = getopt(argc, argv, "e:a:k:c:")) != EOF) {
    switch (opt){
    case 'e': e2=atof(optarg); break;
    case 'a': a=atof(optarg); break;
    case 'k': k0=atof(optarg); break;
    case 'c': desc = optarg; break;
    default:
      fprintf (stderr,
       "usage: genutm [ -e eccentricity ] [ -a equatorial radius ] "
       "[ -k red. factor ]\n");
      exit(1);
    }
  }

  printf ("divert(-1)\n");
  printf ("define(ELLIPSOID, ``%s'')\n", desc);
  printf ("define(A,%0.17lG)\n", a);
  printf ("define(E2,%0.17lG)\n", e2);
  printf ("define(K0,%0.17lG)\n", k0);
  printf ("define(NTERM,%d)\n", NTERM);

  for (i = 0; i < NTERM; i++) {
    Mcos[i] = (i & 1) ? 0.0 : C_Fourier_cos(m, i, N)/M_PI;

    printf ("define(Mcos_%d,%0.17lG)\n", i, Mcos[i]);
  }

  for (i=0; i < NTERM; i++) {
    Ncos[i] = (i&1) ? 0.0 : C_Fourier_cos(n, i, N)/M_PI;
    printf ("define(Ncos_%d,%0.17lG)\n", i, Ncos[i]);
  }
  Betasin[0] = Mcos[0];
  printf ("define(BetaPhi,%0.17lG)\n", Betasin[0]);
  printf ("define(BetaPhi_deg,%0.17lG)\n", Betasin[0]/180.0*M_PI);
  for (i = 1; i < NTERM; i++) {
    Betasin[i] = Mcos[i]/i;
    printf ("define(Betasin_%d,%0.17lG)\n", i, Betasin[i]);
  }
  for (i=0;i<NTERM; i++) {
    A1cos[i] = (i&1) ? C_Fourier_cos(preA1, i, N)/M_PI : 0.0;
    printf ("define(A1cos_%d,%0.17lG)\n", i, A1cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    A2sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preA2, i, N)/M_PI;
    printf ("define(A2sin_%d,%0.17lG)\n", i, A2sin[i]);
  }
  for (i=0; i<NTERM; i++) {
    A3cos[i] = (i&1) ? C_Fourier_cos(preA3, i, N)/M_PI : 0.0;
    printf ("define(A3cos_%d,%0.17lG)\n", i, A3cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    A4sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preA4, i, N)/M_PI;
    printf ("define(A4sin_%d,%0.17lG)\n", i, A4sin[i]);
  }
  for (i=0; i<NTERM; i++) {
    A5cos[i] = (i&1) ? C_Fourier_cos(preA5, i, N)/M_PI : 0.0;
    printf ("define(A5cos_%d,%0.17lG)\n", i, A5cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    A6sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preA6, i, N)/M_PI;
    printf ("define(A6sin_%d,%0.17lG)\n", i, A6sin[i]);
  }

  for (i=0; i<NTERM; i++) {
    Ateb1cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preAteb1, i, N)/M_PI;
    printf ("define(Ateb1cos_%d,%0.17lG)\n", i, Ateb1cos[i]);
    printf ("define(Ateb1cos_deg_%d,%0.17lG)\n", i, Ateb1cos[i]*180.0/M_PI);
  }
  for (i=0; i<NTERM; i++) {
    Ateb2sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preAteb2, i, N)/M_PI;
    printf ("define(Ateb2sin_%d,%0.17lG)\n", i, Ateb2sin[i]);
    printf ("define(Ateb2sin_deg_%d,%0.17lG)\n", i, Ateb2sin[i]*180.0/M_PI);
  }
  for (i=0; i<NTERM; i++) {
    Ateb3cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preAteb3, i, N)/M_PI;
    printf ("define(Ateb3cos_%d,%0.17lG)\n", i, Ateb3cos[i]);
    printf ("define(Ateb3cos_deg_%d,%0.17lG)\n", i, Ateb3cos[i]*180.0/M_PI);
  }
  for (i=0; i<NTERM; i++) {
    Ateb4sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preAteb4, i, N)/M_PI;
    printf ("define(Ateb4sin_%d,%0.17lG)\n", i, Ateb4sin[i]);
    printf ("define(Ateb4sin_deg_%d,%0.17lG)\n", i, Ateb4sin[i]*180.0/M_PI);
  }
  for (i=0; i<NTERM; i++) {
    Ateb5cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preAteb5, i, N)/M_PI;
    printf ("define(Ateb5cos_%d,%0.17lG)\n", i, Ateb5cos[i]);
    printf ("define(Ateb5cos_deg_%d,%0.17lG)\n", i, Ateb5cos[i]*180.0/M_PI);
  }
  for (i=0; i<NTERM; i++) {
    Ateb6sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preAteb6, i, N)/M_PI;
    printf ("define(Ateb6sin_%d,%0.17lG)\n", i, Ateb6sin[i]);
    printf ("define(Ateb6sin_deg_%d,%0.17lG)\n", i, Ateb6sin[i]*180.0/M_PI);
  }
  BetaPI = Beta(M_PI);
  printf ("define(BetaPI,%0.17lG)\n", BetaPI);
  for (i=0; i < NTERM; i++) {
    F1cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF1, i, N)/M_PI;
    printf ("define(F1cos_%d,%0.17lG)\n", i, F1cos[i]);
    printf ("define(F1cos_deg_%d,%0.17lG)\n", i, F1cos[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    F2sin[i] = (i&1) ? C_Fourier_sin(preF2, i, N)/M_PI : 0.0;
    printf ("define(F2sin_%d,%0.17lG)\n", i, F2sin[i]);
    printf ("define(F2sin_deg_%d,%0.17lG)\n", i, F2sin[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    F3cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF3, i, N)/M_PI;
    printf ("define(F3cos_%d,%0.17lG)\n", i, F3cos[i]);
    printf ("define(F3cos_deg_%d,%0.17lG)\n", i, F3cos[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    F4sin[i] = (i&1) ? C_Fourier_sin(preF4, i, N)/M_PI : 0.0;
    printf ("define(F4sin_%d,%0.17lG)\n", i, F4sin[i]);
    printf ("define(F4sin_deg_%d,%0.17lG)\n", i, F4sin[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    F5cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF5, i, N)/M_PI;
    printf ("define(F5cos_%d,%0.17lG)\n", i, F5cos[i]);
    printf ("define(F5cos_deg_%d,%0.17lG)\n", i, F5cos[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    F6sin[i] = (i&1) ? C_Fourier_sin(preF6, i, N)/M_PI : 0.0;
    printf ("define(F6sin_%d,%0.17lG)\n", i, F6sin[i]);
    printf ("define(F6sin_deg_%d,%0.17lG)\n", i, F6sin[i]*180.0/M_PI);
  }

  for (i=0; i < NTERM; i++) {
    dQ2Lat1cos[i] = (i&1) ? C_Fourier_cos(predQ2Lat1, i, N)/M_PI : 0.0;
    printf ("define(dQ2Lat1cos_%d,%0.17lG)\n", i, dQ2Lat1cos[i]);
    printf ("define(dQ2Lat1cos_deg_%d,%0.17lG)\n", i, dQ2Lat1cos[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat2sin[i] = (i&1) ? 0.0 : C_Fourier_sin(predQ2Lat2, i, N)/M_PI;
    printf ("define(dQ2Lat2sin_%d,%0.17lG)\n", i, dQ2Lat2sin[i]);
    printf ("define(dQ2Lat2sin_deg_%d,%0.17lG)\n", i, dQ2Lat2sin[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat3cos[i] = (i&1) ? C_Fourier_cos(predQ2Lat3, i, N)/M_PI : 0.0;
    printf ("define(dQ2Lat3cos_%d,%0.17lG)\n", i, dQ2Lat3cos[i]);
    printf ("define(dQ2Lat3cos_deg_%d,%0.17lG)\n", i, dQ2Lat3cos[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat4sin[i] = (i&1) ? 0.0 : C_Fourier_sin(predQ2Lat4, i, N)/M_PI;
    printf ("define(dQ2Lat4sin_%d,%0.17lG)\n", i, dQ2Lat4sin[i]);
    printf ("define(dQ2Lat4sin_deg_%d,%0.17lG)\n", i, dQ2Lat4sin[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat5cos[i] = (i&1) ? C_Fourier_cos(predQ2Lat5, i, N)/M_PI : 0.0;
    printf ("define(dQ2Lat5cos_%d,%0.17lG)\n", i, dQ2Lat5cos[i]);
    printf ("define(dQ2Lat5cos_deg_%d,%0.17lG)\n", i, dQ2Lat5cos[i]*180.0/M_PI);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat6sin[i] = (i&1) ? 0.0 : C_Fourier_sin(predQ2Lat6, i, N)/M_PI;
    printf ("define(dQ2Lat6sin_%d,%0.17lG)\n", i, dQ2Lat6sin[i]);
    printf ("define(dQ2Lat6sin_deg_%d,%0.17lG)\n", i, dQ2Lat6sin[i]*180.0/M_PI);
  }

  printf ("divert(0)dnl\n");
  exit(0);
}

/* $Id: genutm.c,v 2.6 1998/08/24 13:04:47 luis Exp $ */
