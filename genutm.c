/* $Id: genutm.c,v 2.5 1998/08/06 12:07:25 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Sun May 10 15:25:27 MET DST 1998
 * $Log: genutm.c,v $
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

/* eccentricity of earth (EURO50) squared */
#define E2 0.0067226700223332915        /* Hayford 1950 */
/*#define E2 0.0                          /* Esfera */
/*#define E2 0.00669437999014           /* WGS84 */
/* semi-major axis of ellipsoid */
#define A  6378388.000                  /* Hayford 1950 */
/*#define A 6378137.000                 /* WGS84 */
/* UTM reduction constant */
#define K0 0.9996

/* Number of iterations in Simpson's numerical integration */
#define N 2000
/* Number of terms used in fourier series */
#define NTERM 12

/* N equatorial radius at point of latitude l given in terms of A */
double n(double l)
{
  double sl = sin(l);
  return 1.0/sqrt(1-E2*sl*sl);
}

/* M meridianal radius at point of latitude l given in terms of A */
double m(double l)
{
  double nl = n(l);
  return (1-E2)*nl*nl*nl;
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

double geod2utmX (double lat, double lon)
{
  double res = K0*A*(A1(lat)*lon - A3(lat)*pow(lon, 3.0) + A5(lat)*pow(lon, 5.0));
  return res + 500000;
}

double geod2utmY (double lat, double lon)
{
  double res = K0*A*(Beta(lat) - A2(lat)*pow(lon, 2.0) + A4(lat)*pow(lon, 4.0) - A6(lat)*pow(lon, 6.0));
  if (res < 0.0) res += 10000000.0;
  return res;
}

void utm2geod (double x, double y, double *lat, double *lon)
{
  double phi = Ateb(y);
  double dq;
  dq =
    - B2(phi) * pow(x, 2.0)
    + B4(phi) * pow(x, 4.0)
    - B6(phi) * pow(x, 6.0);
  printf ("utm2geod: dq -> %0.17lg\n", dq);
  *lon =
    + B1(phi) * x
    - B3(phi) * pow(x, 3.0)
    + B5(phi) * pow(x, 5.0);
  *lat = phi
    + dQ2Lat1(phi) * dq
    + dQ2Lat2(phi) * pow(dq, 2.0)
    + dQ2Lat3(phi) * pow(dq, 3.0)
    + dQ2Lat4(phi) * pow(dq, 4.0)
    + dQ2Lat5(phi) * pow(dq, 5.0)
    + dQ2Lat6(phi) * pow(dq, 6.0);
}

double hms2h (double x)
{
  double deg, min;
  x = modf (x, &deg)*100.0;
  x = modf (x, &min)*100.0;
  return deg + min / 60.0 + x / 3600.0;
}

double h2hms (double x)
{
  double deg, min;
  x = modf (x, &deg)*60.0;
  x = modf (x, &min)*60.0;
  return deg + min / 100.0 + x / 10000.0;
}

int huso (double l, double *L, char *zona)
{
  static char *t1 = "CDEFGHJKLMNPQRSTUVWXYZ";
  int h = (int)((*L + M_PI) / M_PI * 30.0) + 1;
  *L = fmod(*L + M_PI/2.0, M_PI/30.0) - M_PI/60.0;
  sprintf (zona, "%2d%c", h, t1 [(int)((l + M_PI/2.25) / M_PI * 22.5)]);
  return h;
}

char *zona (int h, double x, double y)
{
  static char *t1 = "ABCDEFGHJKLMNPQRSTUVWXYZ";
  static char res [3];
  int xx, yy;
  /* first letter */
  h--;
  xx = (int) ((x - 100000.0)/100000.0) + 8*(h % 3);
  yy = (int) (y/100000.0) % 20 + 5*(h % 2);
  sprintf (res, "%c%c", t1[xx], t1[yy]);
  return res;
}

#define OPCION_FUNCIONES 1
#define OPCION_G2U       2
#define OPCION_U2G       4
int opciones = 0;

/* main program */
int main (int argc, char **argv)
{
	char linea [1000];
	double l, L, err;
	int i, opt;

  while ((opt = getopt(argc, argv, "fgu")) != EOF) {
    switch (opt) {
    case 'f':
      opciones |= OPCION_FUNCIONES; break;
    case 'g':
      opciones |= OPCION_U2G; break;
    case 'u':
      opciones |= OPCION_G2U; break;
    default:
      fprintf (stderr, "genutm: opción incorrecta\n");
      break;
    }
  }
  if (!opciones) opciones |= OPCION_G2U;

  if (opciones & OPCION_FUNCIONES) {
    printf ("divert(-1)\n");
    printf ("define(A,%0.17lg)\n", A);
    printf ("define(E2,%0.17lg)\n", E2);
    printf ("define(K0,%0.17lg)\n", K0);
    printf ("define(NTERM,%d)\n", NTERM);
  }
  for (i = 0; i < NTERM; i++) {
    Mcos[i] = (i & 1) ? 0.0 : C_Fourier_cos(m, i, N)/M_PI;

    if (opciones & OPCION_FUNCIONES)
      printf ("define(Mcos_%d,%0.17lg)\n", i, Mcos[i]);
  }

  for (i=0; i < NTERM; i++) {
    Ncos[i] = (i&1) ? 0.0 : C_Fourier_cos(n, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Ncos_%d,%0.17lg)\n", i, Ncos[i]);
  }
  Betasin[0] = Mcos[0];
  if (opciones & OPCION_FUNCIONES)
    printf ("define(BetaPhi,%0.17lg)\n", Betasin[0]);
  for (i = 1; i < NTERM; i++) {
    Betasin[i] = Mcos[i]/i;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Betasin_%d,%0.17lg)\n", i, Betasin[i]);
  }
  for (i=0;i<NTERM; i++) {
    A1cos[i] = (i&1) ? C_Fourier_cos(preA1, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(A1cos_%d,%0.17lg)\n", i, A1cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    A2sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preA2, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(A2sin_%d,%0.17lg)\n", i, A2sin[i]);
  }
  for (i=0; i<NTERM; i++) {
    A3cos[i] = (i&1) ? C_Fourier_cos(preA3, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(A3cos_%d,%0.17lg)\n", i, A3cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    A4sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preA4, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(A4sin_%d,%0.17lg)\n", i, A4sin[i]);
  }
  for (i=0; i<NTERM; i++) {
    A5cos[i] = (i&1) ? C_Fourier_cos(preA5, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(A5cos_%d,%0.17lg)\n", i, A5cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    A6sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preA6, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(A6sin_%d,%0.17lg)\n", i, A6sin[i]);
  }

  for (i=0; i<NTERM; i++) {
    Ateb1cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preAteb1, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Ateb1cos_%d,%0.17lg)\n", i, Ateb1cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    Ateb2sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preAteb2, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Ateb2sin_%d,%0.17lg)\n", i, Ateb2sin[i]);
  }
  for (i=0; i<NTERM; i++) {
    Ateb3cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preAteb3, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Ateb3cos_%d,%0.17lg)\n", i, Ateb3cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    Ateb4sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preAteb4, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Ateb4sin_%d,%0.17lg)\n", i, Ateb4sin[i]);
  }
  for (i=0; i<NTERM; i++) {
    Ateb5cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preAteb5, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Ateb5cos_%d,%0.17lg)\n", i, Ateb5cos[i]);
  }
  for (i=0; i<NTERM; i++) {
    Ateb6sin[i] = (i&1) ? 0.0 : C_Fourier_sin(preAteb6, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(Ateb6sin_%d,%0.17lg)\n", i, Ateb6sin[i]);
  }
  BetaPI = Beta(M_PI);
  if (opciones & OPCION_FUNCIONES)
    printf ("define(BetaPI,%0.17lg)\n", BetaPI);
  for (i=0; i < NTERM; i++) {
    F1cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF1, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(F1cos_%d,%0.17lg)\n", i, F1cos[i]);
  }
  for (i=0; i < NTERM; i++) {
    F2sin[i] = (i&1) ? C_Fourier_sin(preF2, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(F2sin_%d,%0.17lg)\n", i, F2sin[i]);
  }
  for (i=0; i < NTERM; i++) {
    F3cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF3, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(F3cos_%d,%0.17lg)\n", i, F3cos[i]);
  }
  for (i=0; i < NTERM; i++) {
    F4sin[i] = (i&1) ? C_Fourier_sin(preF4, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(F4sin_%d,%0.17lg)\n", i, F4sin[i]);
  }
  for (i=0; i < NTERM; i++) {
    F5cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF5, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(F5cos_%d,%0.17lg)\n", i, F5cos[i]);
  }
  for (i=0; i < NTERM; i++) {
    F6sin[i] = (i&1) ? C_Fourier_sin(preF6, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(F6sin_%d,%0.17lg)\n", i, F6sin[i]);
  }

  for (i=0; i < NTERM; i++) {
    dQ2Lat1cos[i] = (i&1) ? C_Fourier_cos(predQ2Lat1, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(dQ2Lat1cos_%d,%0.17lg)\n", i, dQ2Lat1cos[i]);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat2sin[i] = (i&1) ? 0.0 : C_Fourier_sin(predQ2Lat2, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(dQ2Lat2sin_%d,%0.17lg)\n", i, dQ2Lat2sin[i]);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat3cos[i] = (i&1) ? C_Fourier_cos(predQ2Lat3, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(dQ2Lat3cos_%d,%0.17lg)\n", i, dQ2Lat3cos[i]);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat4sin[i] = (i&1) ? 0.0 : C_Fourier_sin(predQ2Lat4, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(dQ2Lat4sin_%d,%0.17lg)\n", i, dQ2Lat4sin[i]);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat5cos[i] = (i&1) ? C_Fourier_cos(predQ2Lat5, i, N)/M_PI : 0.0;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(dQ2Lat5cos_%d,%0.17lg)\n", i, dQ2Lat5cos[i]);
  }
  for (i=0; i < NTERM; i++) {
    dQ2Lat6sin[i] = (i&1) ? 0.0 : C_Fourier_sin(predQ2Lat6, i, N)/M_PI;
    if (opciones & OPCION_FUNCIONES)
      printf ("define(dQ2Lat6sin_%d,%0.17lg)\n", i, dQ2Lat6sin[i]);
  }

  if (opciones & OPCION_FUNCIONES) {
    printf ("divert(0)dnl\n");
    exit(0);
  }

  if (opciones & OPCION_G2U)
  for (;;) {
	char z [10];
	int h;
	double x, y;
	printf ("##(hh.mmssss hh.mmssss)> ");
  	if (!gets(linea)) break;
	l = L = 0.0;
	sscanf (linea, "%lf%lf", &l, &L);

	printf ("Lat(h.mmssss):  %0.17lg\n", l);
	printf ("Lon(h.mmssss):  %0.17lg\n", L);
	l = hms2h(l); L = hms2h(L);
	printf ("Lat(deg)        %0.17lg\n", l);
	printf ("Lon(deg)        %0.17lg\n", L);
	l *= M_PI/180.0; L *= M_PI/180.0;
	printf ("Lat(rad):       %0.17lg\n", l);
	printf ("Lon(rad):       %0.17lg\n", L);
	h = huso (l, &L, z);
	printf ("Lon(huso):      %0.17lg\n", L);
	printf ("M:              %0.17lg\n", A*m(l));
	printf ("N:              %0.17lg\n", A*n(l));
	printf ("X(m):           %0.3lf\n"
	        "Y(m):           %0.3lf\n",
	  x = geod2utmX(l, L), y = geod2utmY(l, L));
	printf ("zona:            %s%s\n", z, zona(h, x, y));
	printf ("Ateb(y):        %0.17lg\n", Ateb(y/K0/A));

  }
  if (opciones & OPCION_U2G)
  for (;;) {
	char z [10];
	int h;
	double x, y, lat, lon;
	printf ("##(x.xxx y.yyy huso)> ");
  	if (!gets(linea)) break;
	sscanf (linea, "%lf%lf%i", &x, &y, &h);
	utm2geod((x - 500000.0)/K0/A, y/K0/A, &l, &L);
	L += ((h-30) * M_PI/30.0) - M_PI/60.0;

	printf ("Huso(antes):    %d\n", h);
	printf ("Zona(antes):    %s\n", zona(h, x, y));
	printf ("Lat(rad):       %0.17lg\n", l);
	printf ("Lon(rad):       %0.17lg\n", L);
	printf ("Lat(deg)        %0.17lg\n", l * 180.0/M_PI);
	printf ("Lon(deg)        %0.17lg\n", L * 180.0/M_PI);
	printf ("Lat(h.mmssss)   %0.17lg\n", h2hms(l * 180.0/M_PI));
	printf ("Lon(h.mmssss)   %0.17lg\n", h2hms(L * 180.0/M_PI));
	h = huso (l, &L, z);
	printf ("Huso(despues):  %d\n", h);
	printf ("Zona(despues):  %s%s\n", z, zona(h, x, y));
	printf ("Lon(huso/deg):  %0.17lg\n", h2hms(L * 180.0/M_PI));
	printf ("M:              %0.17lg\n", A*m(l));
	printf ("N:              %0.17lg\n", A*n(l));
	printf ("X:              %0.3lf\n"
	        "Y:              %0.3lf\n",
	  x = geod2utmX(l, L), y = geod2utmY(l, L));
	printf ("Ateb(y):        %0.17lg\n", Ateb(y/K0/A));

  }
}

/* $Id: genutm.c,v 2.5 1998/08/06 12:07:25 luis Exp $ */
