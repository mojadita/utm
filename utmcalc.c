/* $Id: utmcalc.c,v 2.2 2002/09/17 19:58:27 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Thu Sep  5 22:44:30 MEST 2002
 * $Log: utmcalc.c,v $
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
  return gs->A * fourier_s(cos,gs->Ncos, l);
} /* geo_N */

/* PUBLIC: M meridional radius at point of latitude l given in terms of A */
double geo_M(struct utmparam *gs, double l)
{
  return gs->A * fourier_s(cos,gs->Mcos, l);
} /* geo_M */

/* PUBLIC: Beta, distance from point to equator */
double geo_Beta(struct utmparam *gs, double x)  /* Beta */
{
  return gs->A * (x*gs->BetaPhi + fourier_s(sin, gs->Betasin, x));
} /* geo_Beta */

/* PUBLIC: Geodetic to UTM conversion routine */
void geo_geod2utm (struct utmparam *gs, double lat, double lon,
	double *x, double *y)
{
  double vs[GEO_NTERM], vc[GEO_NTERM], vp[GEO_NPOT];

  vpot(vp, GEO_NPOT, lon);

  if (x) {
    vfourier(vc, GEO_NTERM, lat, cos);
    *x = gs->AK0 * (
        dot(vc, gs->A1cos, GEO_NTERM) * vp[1]
      - dot(vc, gs->A3cos, GEO_NTERM) * vp[3]
      + dot(vc, gs->A5cos, GEO_NTERM) * vp[5]
      - dot(vc, gs->A7cos, GEO_NTERM) * vp[7]
	  ) + 500000.0;
  }

  if (y) {
#define BETA(x,vs) ((x)*gs->BetaPhi + dot((vs),gs->Betasin,GEO_NTERM))
    vfourier(vs, GEO_NTERM, lat, sin);
    *y = gs->AK0 * (
        BETA(lat,vs)
      - dot(vs, gs->A2sin, GEO_NTERM) * vp[2]
      + dot(vs, gs->A4sin, GEO_NTERM) * vp[4]
      - dot(vs, gs->A6sin, GEO_NTERM) * vp[6]
	  + dot(vs, gs->A8sin, GEO_NTERM) * vp[8]
	  );
  }
} /* geo_geod2utm */

/* PUBLIC: K MODULUS/CONVERGENCE CALCULUS */
void geo_K_conv (struct utmparam *gs, double lat, double lon,
	double *kres, double *conv)
{
  double vs[GEO_NTERM], vc[GEO_NTERM], vp[GEO_NPOT];
  double resx, resy;

  vpot(vp, GEO_NPOT, lon);
  vfourier(vc, GEO_NTERM, lat, cos);
  vfourier(vs, GEO_NTERM, lat, sin);

  resy =
      dot(vc, gs->A1cos, GEO_NTERM) /* ... * vp[0] * 1.0 */
    - dot(vc, gs->A3cos, GEO_NTERM) * vp[2] * 3.0
    + dot(vc, gs->A5cos, GEO_NTERM) * vp[4] * 5.0
    - dot(vc, gs->A7cos, GEO_NTERM) * vp[6] * 7.0
	;

  resx =
      dot(vs, gs->A2sin, GEO_NTERM) * vp[1] * 2.0
    - dot(vs, gs->A4sin, GEO_NTERM) * vp[3] * 4.0
    + dot(vs, gs->A6sin, GEO_NTERM) * vp[5] * 6.0
    - dot(vs, gs->A8sin, GEO_NTERM) * vp[7] * 8.0
	;

  if (kres)
    *kres = sqrt(resx*resx+resy*resy)
		/ dot(vc, gs->Ncos, GEO_NTERM)
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

	y /= gs->A;

	phi0 = y / gs->BetaPI * M_PI;
	vfourier(vs, GEO_NTERM, phi0, sin);
	vfourier(vc, GEO_NTERM, phi0, cos);
	y0 = BETA(phi0,vs);
	dy = y - y0;
	vpot(vp, GEO_NPOT, dy);

	ateb[0] = phi0;
	ateb[1] = dot(gs->Ateb1cos,vc,GEO_NTERM);
	ateb[2] = dot(gs->Ateb2sin,vs,GEO_NTERM);
	ateb[3] = dot(gs->Ateb3cos,vc,GEO_NTERM);
	ateb[4] = dot(gs->Ateb4sin,vs,GEO_NTERM);
	ateb[5] = dot(gs->Ateb5cos,vc,GEO_NTERM);
	ateb[6] = dot(gs->Ateb6sin,vs,GEO_NTERM);
	ateb[7] = dot(gs->Ateb7cos,vc,GEO_NTERM);
	ateb[8] = dot(gs->Ateb8sin,vs,GEO_NTERM);

	return dot(ateb, vp, GEO_NPOT);
} /* geo_Ateb */ 


/* PUBLIC: UTM to geodetic conversion routine */
void geo_utm2geod (struct utmparam *gs, double x, double y,
	double *lat, double *lon)
{
  double vs[GEO_NTERM], vc[GEO_NTERM], vp[GEO_NPOT];
  double phi;
  double dq;

  phi = geo_Ateb(gs, y/GEO_K0);

  x -= 500000.0;
  x /= gs->AK0;
  y /= gs->AK0;

  vpot(vp, GEO_NPOT, x/cos(phi));

  if (lon) {
    vfourier(vc, GEO_NTERM, phi, cos);
    *lon =
      + dot(vc, gs->F1cos, GEO_NTERM) * vp[1]
      - dot(vc, gs->F3cos, GEO_NTERM) * vp[3]
      + dot(vc, gs->F5cos, GEO_NTERM) * vp[5]
      - dot(vc, gs->F7cos, GEO_NTERM) * vp[7]
	  ;
  }

  if (lat) {
    double aux[GEO_NPOT];

    vfourier(vs, GEO_NTERM, phi, sin);
    dq =
      - dot(vs, gs->F2sin, GEO_NTERM) * vp[2]
      + dot(vs, gs->F4sin, GEO_NTERM) * vp[4]
      - dot(vs, gs->F6sin, GEO_NTERM) * vp[6]
      + dot(vs, gs->F8sin, GEO_NTERM) * vp[8]
	  ;

    vpot(vp, GEO_NPOT, dq);
    aux[0] = phi;
    aux[1] = dot(gs->dQ2Lat1cos, vc, GEO_NTERM);
    aux[2] = dot(gs->dQ2Lat2sin, vs, GEO_NTERM);
    aux[3] = dot(gs->dQ2Lat3cos, vc, GEO_NTERM);
    aux[4] = dot(gs->dQ2Lat4sin, vs, GEO_NTERM);
    aux[5] = dot(gs->dQ2Lat5cos, vc, GEO_NTERM);
    aux[6] = dot(gs->dQ2Lat6sin, vs, GEO_NTERM);
    aux[7] = dot(gs->dQ2Lat7cos, vc, GEO_NTERM);
    aux[8] = dot(gs->dQ2Lat8sin, vs, GEO_NTERM);

    *lat = dot(aux, vp, GEO_NPOT);
  }
} /* geo_utm2geod */

/* $Id: utmcalc.c,v 2.2 2002/09/17 19:58:27 luis Exp $ */
