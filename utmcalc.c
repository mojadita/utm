/* $Id: utmcalc.c,v 2.3 2002/09/23 06:14:17 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Thu Sep  5 22:44:30 MEST 2002
 * $Log: utmcalc.c,v $
 * Revision 2.3  2002/09/23 06:14:17  luis
 * Modified to support variable number of GEO_NPOT and GEO_NTERM.
 *
 * Revision 2.2  2002/09/17 19:58:27  luis
 * Added more precision.
 *
 * Revision 2.1  2002/09/06 00:12:11  luis
 * Añadidos utm_ini.h para que genutm pueda calcular por tabla los parámetros y
 * utmcalc.c para los cálculos a partir de los parámetros calculados según la
 * estructura utmparam.
 *
 */

#define IN_UTM_C

/* Standard include files */
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include "utm.h"

/*******************************************/
/*       STATIC PRIVATE FUNCTIONS:         */
/*******************************************/

/* dot product between a and b */
static double dot(double a[], double b[], int n)
{
  int i;
  double res = 0.0;

  for (i = 0; i < n; i++)
    res += a[i] * b[i];

  return res;
} /* dot */

/* this function constructs a vector of components to be used later in several
 * calculus of fourier series. v gets instantiated to sin or cos.
 */
static void vfourier(double v[], int n, double x, double(*f)(double))
{
  int i;
  double aux = 0.0;

  for (i = 0; i < n; i++) {
    v[i] = f(aux);
	aux += x;
  }
} /* vfourier */

/* this function gets a vector of powers of x to be used later in several
 * calculus of taylor series */
static void vpot(double v[], int n, double x)
{
  int i;
  double aux = 1.0;

  for (i = 0; i < n; i++) {
    v[i] = aux;
	aux *= x;
  }
} /* vpot */

/* this function gets the fourier approximation of t[] components on x
 * from sin() or cos() (pointer passed as a parameter) */
static double fourier_s(double (*f)(double), double t[], double x)
{
  double v[GEO_NTERM];

  vfourier(v, GEO_NTERM, x, f);

  return dot(v, t, GEO_NTERM);
} /* fourier_s */


/*******************************************/
/*           PUBLIC FUNCTIONS              */
/*******************************************/

/* PUBLIC: Equatorial Radius */
double geo_N(struct utmparam *gs, double l)
{
  return gs->a * fourier_s(cos,gs->N, l);
} /* geo_N */

/* PUBLIC: M meridional radius at point of latitude l given in terms of A */
double geo_M(struct utmparam *gs, double l)
{
  return gs->a * fourier_s(cos,gs->M, l);
} /* geo_M */

/* PUBLIC: Beta, distance from point to equator */
double geo_Beta(struct utmparam *gs, double x)  /* Beta */
{
  return gs->a * (x*gs->BetaPhi + fourier_s(sin, gs->Beta, x));
} /* geo_Beta */

/* PUBLIC: Geodetic to UTM conversion routine */
void geo_geod2utm (struct utmparam *gs, double lat, double lon,
	double *x, double *y)
{
  double vs[GEO_NTERM], vc[GEO_NTERM], vp[GEO_NPOT];
  int i;

  vpot(vp, GEO_NPOT, lon);

  if (x) {
  	vfourier(vc, GEO_NTERM, lat, cos);
	*x = 0.0;
  }
  if (y) {
#define BETA(x,vs) ((x)*gs->BetaPhi + dot((vs),gs->Beta,GEO_NTERM))
  	vfourier(vs, GEO_NTERM, lat, sin);
	*y = BETA(lat,vs);
  }
  for (i = 1; i < GEO_NPOT; i++) {
  	switch (i & 3) {
	case 1: if (x) *x += dot(vc, gs->A[i], GEO_NTERM) * vp[i]; break;
	case 2: if (y) *y -= dot(vs, gs->A[i], GEO_NTERM) * vp[i]; break;
	case 3: if (x) *x -= dot(vc, gs->A[i], GEO_NTERM) * vp[i]; break;
	case 0: if (y) *y += dot(vs, gs->A[i], GEO_NTERM) * vp[i]; break;
	} /* switch */
  } /* for */
  if (x) {
  	*x *= gs->ak0;
	*x += GEO_FALSE_EASTING;
  }
  if (y) {
  	*y *= gs->ak0;
	*y += GEO_FALSE_NORTHING;
  }
} /* geo_geod2utm */

/* PUBLIC: K MODULUS/CONVERGENCE CALCULUS */
void geo_K_conv (struct utmparam *gs, double lat, double lon,
	double *kres, double *conv)
{
  double vs[GEO_NTERM], vc[GEO_NTERM], vp[GEO_NPOT];
  double resx, resy;
  int i;

  vpot(vp, GEO_NPOT, lon);
  vfourier(vc, GEO_NTERM, lat, cos);
  vfourier(vs, GEO_NTERM, lat, sin);

  resx = resy = 0.0;
  for (i = 0; i < GEO_NPOT; i++) {
  	switch (i & 3) {
	case 0: resy += dot(vc, gs->A[i+1], GEO_NTERM) * vp[i] * (double)(i+1); break;
	case 1: resx += dot(vs, gs->A[i+1], GEO_NTERM) * vp[i] * (double)(i+1); break;
	case 2: resy -= dot(vc, gs->A[i+1], GEO_NTERM) * vp[i] * (double)(i+1); break;
	case 3: resx -= dot(vs, gs->A[i+1], GEO_NTERM) * vp[i] * (double)(i+1); break;
	} /* switch */
  }

  if (kres)
    *kres = sqrt(resx*resx+resy*resy)
		/ dot(vc, gs->N, GEO_NTERM)
		/ vc[1]
		* GEO_K0;

  if (conv)
    *conv = atan2(resx, resy);
} /* geo_K_conv */

/* PUBLIC: Inverse of the function geo_Beta */
double geo_Ateb(struct utmparam *gs, double y)
{
	double phi0;
	double y0, dy;
	double vs[GEO_NTERM], vc[GEO_NTERM], vp[GEO_NPOT], ateb[GEO_NPOT];
	int i;

	y /= gs->a;

	phi0 = y / gs->BetaPI * M_PI;
	vfourier(vs, GEO_NTERM, phi0, sin);
	vfourier(vc, GEO_NTERM, phi0, cos);
	y0 = BETA(phi0,vs);
	dy = y - y0;
	vpot(vp, GEO_NPOT, dy);

	ateb[0] = phi0;
	for (i = 1; i < GEO_NPOT; i++) {
		ateb[i] = dot(gs->Ateb[i], (i & 1 ? vc : vs), GEO_NTERM);
	} /* for */

	return dot(ateb, vp, GEO_NPOT);
} /* geo_Ateb */ 


/* PUBLIC: UTM to geodetic conversion routine */
void geo_utm2geod (struct utmparam *gs, double x, double y,
	double *lat, double *lon)
{
  double vs[GEO_NTERM], vc[GEO_NTERM], vp[GEO_NPOT];
  double phi;
  double dq;
  int i;

  phi = geo_Ateb(gs, y/GEO_K0);

  x -= GEO_FALSE_EASTING;
  y -= GEO_FALSE_NORTHING;
  x /= gs->ak0;
  y /= gs->ak0;

  vpot(vp, GEO_NPOT, x/cos(phi));

  if (lon) {
    vfourier(vc, GEO_NTERM, phi, cos);
	*lon = 0.0;
  }
  if (lat) {
    vfourier(vs, GEO_NTERM, phi, sin);
	dq = 0.0;
  }

  for (i = 1; i < GEO_NPOT; i++) {
  	switch (i & 3) {
	case 1: if (lon) *lon += dot(vc, gs->F[i], GEO_NTERM) * vp[i]; break;
	case 2: if (lat)  dq  -= dot(vs, gs->F[i], GEO_NTERM) * vp[i]; break;
	case 3: if (lon) *lon -= dot(vc, gs->F[i], GEO_NTERM) * vp[i]; break;
	case 0: if (lat)  dq  += dot(vs, gs->F[i], GEO_NTERM) * vp[i]; break;
	} /* switch */
  } /* for */

  if (lat) {
    double aux[GEO_NPOT];

    vpot(vp, GEO_NPOT, dq);
    aux[0] = phi;
	for (i = 1; i < GEO_NPOT; i++) {
		aux[i] = dot(gs->dQ2Lat[i], (i & 1 ? vc : vs), GEO_NTERM);
	} /* for */

    *lat = dot(aux, vp, GEO_NPOT);
  } /* if */
} /* geo_utm2geod */

/* $Id: utmcalc.c,v 2.3 2002/09/23 06:14:17 luis Exp $ */
