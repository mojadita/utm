/* $Id: genutm.c,v 2.8 2002/09/17 19:58:27 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Sun May 10 15:25:27 MET DST 1998
 * $Log: genutm.c,v $
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

/* Number of iterations in Simpson's numerical integration */
#define N 3072

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
double n(double l)
{
  double sl = sin(l);
  return 1.0 / sqrt(1 - sg->e2 * sl*sl);
} /* n */

/* M meridianal radius at point of latitude l given in terms of A */
double m(double l)
{
  double nl = n(l);
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
double(*F)(double);
double(*SC)(double);
double I;

double F_Fourier(double x)
{
	return F(x)*SC(I*x);
} /* F_Fourier */

double C_Fourier_sin(double(*f)(double), int i, int n)
{
  F = f; SC = sin; I = i;

  return simpson(F_Fourier, 0.0, 2.0*M_PI, n);
} /* C_Fourier_sin */

double C_Fourier_cos(double(*f)(double), int i, int n)
{
  double result;

  F = f; SC = cos; I = i;
  result = simpson(F_Fourier, 0.0, 2*M_PI, n);
  if(i == 0) result /= 2.0;

  return result;
} /* C_Fourier_cos */


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
    res += sg->Betasin[i] * sin((double)i*x);

  return res;
} /* Beta */

double preA1(double x) /* N * cos(l) */
{
  return n(x)*cos(x);
} /* preA1 */


double A1(double x)
{
  int i;
  double res;

  res = 0.0;
  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A1cos[i] * cos(i*x);

  return res;
} /* A1 */

double preA2(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->A1cos[i] * i * sin(i*x);

  return n(x)/m(x)*cos(x)*res/2.0;
} /* preA2 */

double A2(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A2sin[i] * sin(i*x);

  return res;
} /* A2 */

double preA3(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A2sin[i] * i * cos(i*x);

  return n(x)*cos(x)/m(x)*res/3.0;
} /* preA3 */

double A3(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A3cos[i] * cos(i*x);

  return res;
} /* A3 */

double preA4(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->A3cos[i] * i * sin(i*x);

  return n(x)*cos(x)/m(x)*res/4.0;
} /* preA4 */

double A4(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A4sin[i] * sin(i*x);

  return res;
} /* A4 */

double preA5(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A4sin[i] * i * cos(i*x);

  return n(x)*cos(x)/m(x)*res/5.0;
} /* preA5 */

double A5(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A5cos[i] * cos(i*x);

  return res;
} /* A5 */

double preA6(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->A5cos[i] * i * sin(i*x);

  return n(x)*cos(x)/m(x)*res/6.0;
} /* preA6 */

double A6(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A6sin[i] * sin(i*x);

  return res;
} /* A6 */

double preA7(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A6sin[i] * i * cos(i*x);

  return n(x)*cos(x)/m(x)*res/7.0;
} /* preA7 */

double A7(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A7cos[i] * cos(i*x);

  return res;
} /* A7 */

double preA8(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->A7cos[i] * i * sin(i*x);

  return n(x)*cos(x)/m(x)*res/8.0;
} /* preA8 */

double A8(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->A8sin[i] * sin(i*x);

  return res;
} /* A8 */

/************** TRANSFORMACIÓN UTM -> GEODÉSICAS ****************/

/* FUNCIÓN ATEB, INVERSA DE LA FUNCIÓN BETA */
double preAteb1(double x)
{
  return 1.0/m(x);
} /* preAteb1 */

double Ateb1(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb1cos[i] * cos(i*x);

  return res;
} /* Ateb1 */

double derAteb1(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->Ateb1cos[i] * i * sin(i*x);

  return res;
} /* derAteb1 */

double preAteb2(double x)
{
  return derAteb1(x) / m(x) / 2.0;
} /* preAteb2 */

double Ateb2(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb2sin[i] * sin(i*x);

  return res;
} /* Ateb2 */

double derAteb2(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb2sin[i] * i * cos(i*x);

  return res;
} /* derAteb2 */

double preAteb3(double x)
{
  return derAteb2(x) / m(x) / 3.0;
} /* preAteb3 */

double Ateb3(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb3cos[i] * cos(i*x);

  return res;
} /* Ateb3 */

double derAteb3(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->Ateb3cos[i] * i * sin(i*x);

  return res;
} /* derAteb3 */

double preAteb4(double x)
{
  return derAteb3(x) / m(x) / 4.0;
} /* preAteb4 */

double Ateb4(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb4sin[i] * sin(i*x);

  return res;
} /* Ateb4 */

double derAteb4(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb4sin[i] * i * cos(i*x);

  return res;
} /* derAteb4 */

double preAteb5(double x)
{
  return derAteb4(x) / m(x) / 5.0;
} /* preAteb5 */

double Ateb5(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb5cos[i] * cos(i*x);

  return res;
} /* Ateb5 */

double derAteb5(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->Ateb5cos[i] * i * sin(i*x);

  return res;
} /* derAteb5 */

double preAteb6(double x)
{
  return derAteb5(x) / m(x) / 6.0;
} /* preAteb6 */

double Ateb6(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb6sin[i] * sin(i*x);

  return res;
} /* Ateb6 */

double derAteb6(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb6sin[i] * i * cos(i*x);

  return res;
} /* derAteb6 */

double preAteb7(double x)
{
  return derAteb6(x) / m(x) / 7.0;
} /* preAteb7 */

double Ateb7(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb7cos[i] * cos(i*x);

  return res;
} /* Ateb7 */

double derAteb7(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= sg->Ateb7cos[i] * i * sin(i*x);

  return res;
} /* derAteb7 */

double preAteb8(double x)
{
  return derAteb7(x) / m(x) / 8.0;
} /* preAteb8 */

double Ateb8(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->Ateb8sin[i] * sin(i*x);

  return res;
} /* Ateb8 */

double Ateb(double x)
{
	double phi0 = x / sg->BetaPI * M_PI;
	double x0 = Beta(phi0);
	double dx = x - x0;
	return
	    phi0
	  + Ateb1(phi0) * dx
	  + Ateb2(phi0) * dx * dx
	  + Ateb3(phi0) * pow(dx, 3.0)
	  + Ateb4(phi0) * pow(dx, 4.0)
	  + Ateb5(phi0) * pow(dx, 5.0)
	  + Ateb6(phi0) * pow(dx, 6.0)
	  + Ateb7(phi0) * pow(dx, 7.0)
	  + Ateb8(phi0) * pow(dx, 8.0);
} /* Ateb */

/*******************/

double preF1(double x)
{
  return 1.0/n(x);
} /* preF1 */

double F1(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F1cos[i] * cos(i*x);

  return res;
} /* F1 */

double B1(double x)
{
  return F1(x) / cos(x);
} /* B1 */

double dF1(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= i * sg->F1cos[i] * sin(i*x);

  return res;
} /* dF1 */

double preF2(double x)
{
  return(dF1(x)*cos(x) + F1(x)*sin(x))/m(x)/2.0;
} /* preF2 */

double F2(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F2sin[i] * sin(i*x);

  return res;
} /* F2 */

double B2(double x)
{
  return F2(x)/pow(cos(x), 2.0);
} /* B2 */

double dF2(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += i * sg->F2sin[i] * cos(i*x);

  return res;
} /* dF2 */

double preF3(double x)
{
  return(dF2(x) * cos(x) + 2.0 * F2(x) * sin(x))/3.0/m(x);
} /* preF3 */

double F3(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F3cos[i] * cos(i*x);

  return res;
} /* F3 */

double B3(double x)
{
  return F3(x)/pow(cos(x), 3.0);
} /* B3 */

double dF3(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= i * sg->F3cos[i] * sin(i*x);

  return res;
} /* dF3 */

double preF4(double x)
{
  return(dF3(x)*cos(x) + 3.0*F3(x)*sin(x))/4.0/m(x);
} /* preF4 */

double F4(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F4sin[i] * sin(i*x);

  return res;
} /* F4 */

double B4(double x)
{
  return F4(x)/pow(cos(x), 4.0);
} /* B4 */

double dF4(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += i * sg->F4sin [i] * cos(i*x);

  return res;
} /* dF4 */

double preF5(double x)
{
  return(dF4(x)*cos(x) + 4.0*F4(x)*sin(x))/5.0/m(x);
} /* preF5 */

double F5(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F5cos[i] * cos(i*x);

  return res;
} /* F5 */

double B5(double x)
{
  return F5(x)/pow(cos(x), 5.0);
} /* B5 */

double dF5(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= i * sg->F5cos [i] * sin(i*x);

  return res;
} /* dF5 */

double preF6(double x)
{
  return(dF5(x)*cos(x) + 5.0*F5(x)*sin(x))/6.0/m(x);
} /* preF6 */

double F6(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F6sin[i] * sin(i*x);

  return res;
} /* F6 */

double B6(double x)
{
  return F6(x)/pow(cos(x), 6.0);
} /* B6 */

double dF6(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += i * sg->F6sin [i] * cos(i*x);

  return res;
} /* dF6 */

double preF7(double x)
{
  return(dF6(x)*cos(x) + 6.0*F6(x)*sin(x))/7.0/m(x);
} /* preF7 */

double F7(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F7cos[i] * cos(i*x);

  return res;
} /* F7 */

double B7(double x)
{
  return F7(x)/pow(cos(x), 7.0);
} /* B7 */

double dF7(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res -= i * sg->F7cos [i] * sin(i*x);

  return res;
} /* dF7 */

double preF8(double x)
{
  return(dF7(x)*cos(x) + 7.0*F7(x)*sin(x))/8.0/m(x);
} /* preF8 */

double F8(double x)
{
  int i;
  double res = 0.0;

  for (i = 0; i < GEO_NTERM; i++)
    res += sg->F8sin[i] * sin(i*x);

  return res;
} /* F8 */

double B8(double x)
{
  return F8(x)/pow(cos(x), 8.0);
} /* B8 */

/************ INCREMENTO DE Q -> INCREMENTO DE LATITUD **************/

double predQ2Lat1(double phi)
{
  return n(phi)/m(phi)*cos(phi);
} /* predQ2Lat1 */

double dQ2Lat1(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat1cos[i] * cos(i*phi);

  return res;
} /* dQ2Lat1 */

double derdQ2Lat1(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res -= sg->dQ2Lat1cos[i] * i * sin(i*phi);

  return res;
} /* derdQ2Lat1 */

double predQ2Lat2(double phi)
{
  return n(phi)/m(phi) * cos(phi) * derdQ2Lat1(phi) / 2.0;
} /* predQ2Lat2 */

double dQ2Lat2(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat2sin[i] * sin(i*phi);

  return res;
} /* dQ2Lat2 */

double derdQ2Lat2(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat2sin[i] * i * cos(i*phi);

  return res;
} /* derdQ2Lat2 */

double predQ2Lat3(double phi)
{
  return n(phi)/m(phi) * cos(phi) * derdQ2Lat2(phi) / 3.0;
} /* predQ2Lat3 */

double dQ2Lat3(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat3cos[i] * cos(i*phi);

  return res;
} /* dQ2Lat3 */

double derdQ2Lat3(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res -= sg->dQ2Lat3cos[i] * i * sin(i*phi);

  return res;
} /* derdQ2Lat3 */

double predQ2Lat4(double phi)
{
  return n(phi)/m(phi) * cos(phi) * derdQ2Lat3(phi) / 4.0;
} /* predQ2Lat4 */

double dQ2Lat4(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat4sin[i] * sin(i*phi);

  return res;
} /* dQ2Lat4 */

double derdQ2Lat4(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat4sin[i] * i * cos(i*phi);

  return res;
} /* derdQ2Lat4 */

double predQ2Lat5(double phi)
{
  return n(phi)/m(phi) * cos(phi) * derdQ2Lat4(phi) / 5.0;
} /* predQ2Lat5 */

double dQ2Lat5(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat5cos[i] * cos(i*phi);

  return res;
} /* dQ2Lat5 */

double derdQ2Lat5(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res -= sg->dQ2Lat5cos[i] * i * sin(i*phi);

  return res;
} /* derdQ2Lat5 */

double predQ2Lat6(double phi)
{
  return n(phi)/m(phi) * cos(phi) * derdQ2Lat5(phi) / 6.0;
} /* predQ2Lat6 */

double dQ2Lat6(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat6sin[i] * sin(i*phi);

  return res;
} /* dQ2Lat6 */

double derdQ2Lat6(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat6sin[i] * i * cos(i*phi);

  return res;
} /* derdQ2Lat6 */

double predQ2Lat7(double phi)
{
  return n(phi)/m(phi) * cos(phi) * derdQ2Lat6(phi) / 7.0;
} /* predQ2Lat7 */

double dQ2Lat7(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat7cos[i] * cos(i*phi);

  return res;
} /* dQ2Lat7 */

double derdQ2Lat7(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res -= sg->dQ2Lat7cos[i] * i * sin(i*phi);

  return res;
} /* derdQ2Lat7 */

double predQ2Lat8(double phi)
{
  return n(phi)/m(phi) * cos(phi) * derdQ2Lat7(phi) / 8.0;
} /* predQ2Lat8 */

double dQ2Lat8(double phi)
{
  int i;
  double res = 0.0;

  for (i=0; i < GEO_NTERM; i++)
    res += sg->dQ2Lat8sin[i] * sin(i*phi);

  return res;
} /* dQ2Lat8 */

/*******************************************************************/
/*******************************************************************/
/* main program */
int main(int argc, char **argv)
{
	char linea [1000];
	double l, L, err;
	int i, opt;
	extern char *optarg;

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
  printf("define(A,%0.17lG)\n", sg->A);
  printf("define(B,%0.17lG)\n", sg->A*sqrt(1-sg->e2));
  printf("define(E2,%0.17lG)\n", sg->e2);
  printf("define(AK0,%0.17lG)\n", sg->A*GEO_K0);
  printf("define(NTERM,%d)\n", GEO_NTERM);

  /* Mcos */
  for (i = 0; i < GEO_NTERM; i++) {
    sg->Mcos[i] = (i & 1) ? 0.0 : C_Fourier_cos(m, i, N)/M_PI;
    printf("define(Mcos_%d,%0.17lG)\n", i, sg->Mcos[i]);
  }

  /* Ncos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->Ncos[i] = (i & 1) ? 0.0 : C_Fourier_cos(n, i, N)/M_PI;
    printf("define(Ncos_%d,%0.17lG)\n", i, sg->Ncos[i]);
  }

  /* BetaPhi */
  sg->BetaPhi = sg->Mcos[0];
  printf("define(BetaPhi,%0.17lG)\n", sg->BetaPhi);
  printf("define(BetaPhi_deg,%0.17lG)\n", sg->BetaPhi/180.0*M_PI);

  /* Betasin */
  sg->Betasin[0] = 0.0;
  printf("define(Betasin_0,%0.17lG)\n", 0.0);
  for (i = 1; i < GEO_NTERM; i++) {
    sg->Betasin[i] = sg->Mcos[i]/i;
    printf("define(Betasin_%d,%0.17lG)\n", i, sg->Betasin[i]);
  }

  /* A1cos */
  for (i = 0; i < GEO_NTERM; i++) {
    sg->A1cos[i] = (i & 1) ? C_Fourier_cos(preA1, i, N)/M_PI : 0.0;
    printf("define(A1cos_%d,%0.17lG)\n", i, sg->A1cos[i]);
  }

  /* A2sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->A2sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preA2, i, N)/M_PI;
    printf("define(A2sin_%d,%0.17lG)\n", i, sg->A2sin[i]);
  }

  /* A3cos */
  for (i=0; i<GEO_NTERM; i++) {
    sg->A3cos[i] = (i & 1) ? C_Fourier_cos(preA3, i, N)/M_PI : 0.0;
    printf("define(A3cos_%d,%0.17lG)\n", i, sg->A3cos[i]);
  }

  /* A4sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->A4sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preA4, i, N)/M_PI;
    printf("define(A4sin_%d,%0.17lG)\n", i, sg->A4sin[i]);
  }

  /* A5cos */
  for (i=0; i<GEO_NTERM; i++) {
    sg->A5cos[i] = (i & 1) ? C_Fourier_cos(preA5, i, N)/M_PI : 0.0;
    printf("define(A5cos_%d,%0.17lG)\n", i, sg->A5cos[i]);
  }

  /* A6sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->A6sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preA6, i, N)/M_PI;
    printf("define(A6sin_%d,%0.17lG)\n", i, sg->A6sin[i]);
  }

  /* A7cos */
  for (i=0; i<GEO_NTERM; i++) {
    sg->A7cos[i] = (i & 1) ? C_Fourier_cos(preA7, i, N)/M_PI : 0.0;
    printf("define(A7cos_%d,%0.17lG)\n", i, sg->A7cos[i]);
  }

  /* A8sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->A8sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preA8, i, N)/M_PI;
    printf("define(A8sin_%d,%0.17lG)\n", i, sg->A8sin[i]);
  }

  /* Ateb1cos */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb1cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preAteb1, i, N)/M_PI;
    printf("define(Ateb1cos_%d,%0.17lG)\n", i, sg->Ateb1cos[i]);
    printf("define(Ateb1cos_deg_%d,%0.17lG)\n", i, sg->Ateb1cos[i]*180.0/M_PI);
  }

  /* Ateb2sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb2sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preAteb2, i, N)/M_PI;
    printf("define(Ateb2sin_%d,%0.17lG)\n", i, sg->Ateb2sin[i]);
    printf("define(Ateb2sin_deg_%d,%0.17lG)\n", i, sg->Ateb2sin[i]*180.0/M_PI);
  }

  /* Ateb3cos */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb3cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preAteb3, i, N)/M_PI;
    printf("define(Ateb3cos_%d,%0.17lG)\n", i, sg->Ateb3cos[i]);
    printf("define(Ateb3cos_deg_%d,%0.17lG)\n", i, sg->Ateb3cos[i]*180.0/M_PI);
  }

  /* Ateb4sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb4sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preAteb4, i, N)/M_PI;
    printf("define(Ateb4sin_%d,%0.17lG)\n", i, sg->Ateb4sin[i]);
    printf("define(Ateb4sin_deg_%d,%0.17lG)\n", i, sg->Ateb4sin[i]*180.0/M_PI);
  }

  /* Ateb5cos */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb5cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preAteb5, i, N)/M_PI;
    printf("define(Ateb5cos_%d,%0.17lG)\n", i, sg->Ateb5cos[i]);
    printf("define(Ateb5cos_deg_%d,%0.17lG)\n", i, sg->Ateb5cos[i]*180.0/M_PI);
  }

  /* Ateb6sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb6sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preAteb6, i, N)/M_PI;
    printf("define(Ateb6sin_%d,%0.17lG)\n", i, sg->Ateb6sin[i]);
    printf("define(Ateb6sin_deg_%d,%0.17lG)\n", i, sg->Ateb6sin[i]*180.0/M_PI);
  }

  /* Ateb7cos */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb7cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preAteb7, i, N)/M_PI;
    printf("define(Ateb7cos_%d,%0.17lG)\n", i, sg->Ateb7cos[i]);
    printf("define(Ateb7cos_deg_%d,%0.17lG)\n", i, sg->Ateb7cos[i]*180.0/M_PI);
  }

  /* Ateb8sin */
  for (i=0; i<GEO_NTERM; i++) {
    sg->Ateb8sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(preAteb8, i, N)/M_PI;
    printf("define(Ateb8sin_%d,%0.17lG)\n", i, sg->Ateb8sin[i]);
    printf("define(Ateb8sin_deg_%d,%0.17lG)\n", i, sg->Ateb8sin[i]*180.0/M_PI);
  }

  /* BetaPI */
  sg->BetaPI = Beta(M_PI);
  printf("define(BetaPI,%0.17lG)\n", sg->BetaPI);

  /* F1cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F1cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preF1, i, N)/M_PI;
    printf("define(F1cos_%d,%0.17lG)\n", i, sg->F1cos[i]);
    printf("define(F1cos_deg_%d,%0.17lG)\n", i, sg->F1cos[i]*180.0/M_PI);
  }

  /* F2sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F2sin[i] = (i & 1) ? C_Fourier_sin(preF2, i, N)/M_PI : 0.0;
    printf("define(F2sin_%d,%0.17lG)\n", i, sg->F2sin[i]);
    printf("define(F2sin_deg_%d,%0.17lG)\n", i, sg->F2sin[i]*180.0/M_PI);
  }

  /* F3cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F3cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preF3, i, N)/M_PI;
    printf("define(F3cos_%d,%0.17lG)\n", i, sg->F3cos[i]);
    printf("define(F3cos_deg_%d,%0.17lG)\n", i, sg->F3cos[i]*180.0/M_PI);
  }

  /* F4sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F4sin[i] = (i & 1) ? C_Fourier_sin(preF4, i, N)/M_PI : 0.0;
    printf("define(F4sin_%d,%0.17lG)\n", i, sg->F4sin[i]);
    printf("define(F4sin_deg_%d,%0.17lG)\n", i, sg->F4sin[i]*180.0/M_PI);
  }

  /* F5cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F5cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preF5, i, N)/M_PI;
    printf("define(F5cos_%d,%0.17lG)\n", i, sg->F5cos[i]);
    printf("define(F5cos_deg_%d,%0.17lG)\n", i, sg->F5cos[i]*180.0/M_PI);
  }

  /* F6sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F6sin[i] = (i & 1) ? C_Fourier_sin(preF6, i, N)/M_PI : 0.0;
    printf("define(F6sin_%d,%0.17lG)\n", i, sg->F6sin[i]);
    printf("define(F6sin_deg_%d,%0.17lG)\n", i, sg->F6sin[i]*180.0/M_PI);
  }

  /* F7cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F7cos[i] = (i & 1) ? 0.0 : C_Fourier_cos(preF7, i, N)/M_PI;
    printf("define(F7cos_%d,%0.17lG)\n", i, sg->F7cos[i]);
    printf("define(F7cos_deg_%d,%0.17lG)\n", i, sg->F7cos[i]*180.0/M_PI);
  }

  /* F8sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->F8sin[i] = (i & 1) ? C_Fourier_sin(preF8, i, N)/M_PI : 0.0;
    printf("define(F8sin_%d,%0.17lG)\n", i, sg->F8sin[i]);
    printf("define(F8sin_deg_%d,%0.17lG)\n", i, sg->F8sin[i]*180.0/M_PI);
  }

  /* dQ2Lat1cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat1cos[i] = (i & 1) ? C_Fourier_cos(predQ2Lat1, i, N)/M_PI : 0.0;
    printf("define(dQ2Lat1cos_%d,%0.17lG)\n", i, sg->dQ2Lat1cos[i]);
    printf("define(dQ2Lat1cos_deg_%d,%0.17lG)\n", i, sg->dQ2Lat1cos[i]*180.0/M_PI);
  }

  /* dQ2Lat2sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat2sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(predQ2Lat2, i, N)/M_PI;
    printf("define(dQ2Lat2sin_%d,%0.17lG)\n", i, sg->dQ2Lat2sin[i]);
    printf("define(dQ2Lat2sin_deg_%d,%0.17lG)\n", i, sg->dQ2Lat2sin[i]*180.0/M_PI);
  }

  /* dQ2Lat3cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat3cos[i] = (i & 1) ? C_Fourier_cos(predQ2Lat3, i, N)/M_PI : 0.0;
    printf("define(dQ2Lat3cos_%d,%0.17lG)\n", i, sg->dQ2Lat3cos[i]);
    printf("define(dQ2Lat3cos_deg_%d,%0.17lG)\n", i, sg->dQ2Lat3cos[i]*180.0/M_PI);
  }

  /* dQ2Lat4sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat4sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(predQ2Lat4, i, N)/M_PI;
    printf("define(dQ2Lat4sin_%d,%0.17lG)\n", i, sg->dQ2Lat4sin[i]);
    printf("define(dQ2Lat4sin_deg_%d,%0.17lG)\n", i, sg->dQ2Lat4sin[i]*180.0/M_PI);
  }

  /* dQ2Lat5cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat5cos[i] = (i & 1) ? C_Fourier_cos(predQ2Lat5, i, N)/M_PI : 0.0;
    printf("define(dQ2Lat5cos_%d,%0.17lG)\n", i, sg->dQ2Lat5cos[i]);
    printf("define(dQ2Lat5cos_deg_%d,%0.17lG)\n", i, sg->dQ2Lat5cos[i]*180.0/M_PI);
  }

  /* dQ2Lat6sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat6sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(predQ2Lat6, i, N)/M_PI;
    printf("define(dQ2Lat6sin_%d,%0.17lG)\n", i, sg->dQ2Lat6sin[i]);
    printf("define(dQ2Lat6sin_deg_%d,%0.17lG)\n", i, sg->dQ2Lat6sin[i]*180.0/M_PI);
  }

  /* dQ2Lat7cos */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat7cos[i] = (i & 1) ? C_Fourier_cos(predQ2Lat7, i, N)/M_PI : 0.0;
    printf("define(dQ2Lat7cos_%d,%0.17lG)\n", i, sg->dQ2Lat7cos[i]);
    printf("define(dQ2Lat7cos_deg_%d,%0.17lG)\n", i, sg->dQ2Lat7cos[i]*180.0/M_PI);
  }

  /* dQ2Lat8sin */
  for (i=0; i < GEO_NTERM; i++) {
    sg->dQ2Lat8sin[i] = (i & 1) ? 0.0 : C_Fourier_sin(predQ2Lat8, i, N)/M_PI;
    printf("define(dQ2Lat8sin_%d,%0.17lG)\n", i, sg->dQ2Lat8sin[i]);
    printf("define(dQ2Lat8sin_deg_%d,%0.17lG)\n", i, sg->dQ2Lat8sin[i]*180.0/M_PI);
  }

  printf("divert(0)dnl\n");
  exit(0);
} /* main */

/* $Id: genutm.c,v 2.8 2002/09/17 19:58:27 luis Exp $ */
