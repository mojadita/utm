divert(-1)
changecom(,)
define(`lq',changequote([,])[changequote([,])`changequote(`,')]changequote(`,'))
define(`rq',changequote([,])[changequote([,])'changequote(`,')]changequote(`,'))
define(`id', ``'$1_$2`'')
define(`for',`ifelse(eval((`$2') <= (`$3')),1,`pushdef(`$1',`$2')`'$4`'dnl
popdef(`$1')for(`$1',eval((`$2')+1),`$3', `$4')')')
define(`defineTable',`static double $1 [`NTERM'] = {
for(`I',0,NTERM-1,`    id(`$1',I),   /* #I */
')dnl
}; /* $1 */
')
ifndef(`PREFIX',define(PREFIX,`'))
divert(0)dnl
/* $Id: utm.c.m4,v 2.3 1998/08/17 19:02:20 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Mon Aug 10 15:54:07 MET DST 1998
 * $Log: utm.c.m4,v $
 * Revision 2.3  1998/08/17 19:02:20  luis
 * Elimination of spanish comments in the code.
 *
 * Revision 2.2  1998/08/17 18:58:12  luis
 * Inclusion of Automatically generated message into the source code.
 *
 * Revision 2.1  1998/08/17 18:54:55  luis
 * *** empty log message ***
 *
 * Revision 2.0  1998/08/17 18:50:13  luis
 * Inclusion of Linear Deformation Modulus lq`'K`'rq and meridian convergence `for'
 * ease of distance and azimuth calculus.
 *
 * Revision 1.1  1998/08/17 13:43:24  luis
 * Initial revision
 *
 * CAUTION: THIS FILE HAS BEEN GENERATED AUTOMATICALLY FROM CONSTANTS
 * CALCULUS AND M4 lq`'TEMPLATE`'rq FILE.
 * EDIT ONLY IF YOU ARE UNABLE TO GENERATE IT AUTOMATICALLY.
 */

#`define' IN_UTM_C

/* Standard `include' files */
#`include' <sys/types.h>
#`include' <stdlib.h>
#`include' <math.h>

/* eccentricity of earth squared */
#`define' `E2' E2
double PREFIX`'e2 = `E2';
/* semi-major axis of ellipsoid */
#`define' `A'  A 
double PREFIX`'a = `A';
/* UTM reduction constant */
#`define' `K0' K0
double PREFIX`'k0 = `K0';
#define `ELLIPSOID' "ELLIPSOID"
char *PREFIX`'dsc = `ELLIPSOID';


/* Number of terms used in fourier series */
#`define' `NTERM' NTERM
#`define' NPOT 7

static double dot(double a[], double b[], int n)
{ int i;
  double res = 0.0;
  `for' (i=0; i < n; i++)
    res += a[i]*b[i];
  return res;
}

static void vfourier(double v[], double x, double(*f)(double))
{ int i;
  double aux = 0.0;

  `for' (i=0; i < `NTERM'; i++) {
    v[i] = f(aux); aux += x;
  }
}

static void vpot(double v[], double x)
{ int i;
  double aux = 1.0;
  `for' (i=0; i < `NTERM'; i++) {
    v[i] = aux; aux *= x;
  }
}

static double ffourier(double (*f)(double), double t[], double x)
{ double v[`NTERM'];

  vfourier(v, x, f);
  return dot(v, t, `NTERM');
}

/* N equatorial radius at point of latitude l given in terms of `A' */
defineTable(`Ncos')
double PREFIX`'n(double l)
{ return ffourier(cos,Ncos, l);
}

/* M meridianal radius at point of latitude l given in terms of `A' */
defineTable(`Mcos')
double PREFIX`'m(double l)
{ return ffourier(cos,Mcos, l);
}

/* `Beta', distance from point to equator */
static double PREFIX`BetaPhi' = BetaPhi;
define(`Betasin_0',0)dnl
defineTable(`Betasin')
double PREFIX`'Beta(double x)  /* Beta */
{
  return x*`BetaPhi' + ffourier(sin, `Betasin', x);
}
#`define' BETA(x,vs) ((x)*`BetaPhi' + dot((vs),`Betasin',`NTERM'))


/* `A' tables, to calculate geod->utm transformation */
defineTable(`A1cos')
defineTable(`A2sin')
defineTable(`A3cos')
defineTable(`A4sin')
defineTable(`A5cos')
defineTable(`A6sin')

void PREFIX`'geod2utm (double lat, double lon, double *x, double *y)
{
  double vs[`NTERM'], vc[`NTERM'], vp[NPOT];

  vpot(vp, lon);

  if (x) {
    vfourier(vc, lat, cos);
    *x =
        dot(vc, `A1cos', `NTERM') * vp[1]
      - dot(vc, `A3cos', `NTERM') * vp[3]
      + dot(vc, `A5cos', `NTERM') * vp[5];
  }

  if (y) {
    vfourier(vs, lat, sin);
    *y =
        BETA(lat,vs)
      - dot(vs, `A2sin', `NTERM') * vp[2]
      + dot(vs, `A4sin', `NTERM') * vp[4]
      - dot(vs, `A6sin', `NTERM') * vp[6];
  }
}

/************** K MODULUS/CONVERGENCE CALCULUS ******************/

void PREFIX`'K_conv (double lat, double lon, double *kres, double *conv)
{ double vs[`NTERM'], vc[`NTERM'], vp[NPOT];
  double resx, resy;

  vpot(vp, lon);
  vfourier(vc, lat, cos);
  vfourier(vs, lat, sin);
  resy =
      dot(vc, `A1cos', `NTERM')
    - dot(vc, `A3cos', `NTERM') * vp[2] * 3.0
    + dot(vc, `A5cos', `NTERM') * vp[4] * 5.0;
  resx =
      dot(vs, `A2sin', `NTERM') * vp[1] * 2.0
    - dot(vs, `A4sin', `NTERM') * vp[3] * 4.0
    + dot(vs, `A6sin', `NTERM') * vp[5] * 6.0;
  if (kres)
    *kres = sqrt(resx*resx+resy*resy) / dot(vc, Ncos, `NTERM') / vc[1];
  if (conv)
    *conv = atan2(resx, resy);
}

/************** GEODETIC -> UTM CALCULATION ****************/

/* FUNCIÓN ATEB, INVERSA DE LA FUNCIÓN BETA */
defineTable(`Ateb1cos')
defineTable(`Ateb2sin')
defineTable(`Ateb3cos')
defineTable(`Ateb4sin')
defineTable(`Ateb5cos')
defineTable(`Ateb6sin')

static double `BetaPI' = BetaPI;

double PREFIX`'Ateb(double y)
{
	double phi0 = y / `BetaPI' * M_PI;
	double y0, dy;
	double vs[`NTERM'], vc[`NTERM'], vp[NPOT], ateb[NPOT];

	vfourier(vs, phi0, sin);
	vfourier(vc, phi0, cos);
	y0 = BETA(phi0,vs);
	dy = y - y0;
	vpot(vp, dy);

	ateb[0] = phi0;
	ateb[1] = dot(Ateb1cos,vc,`NTERM');
	ateb[2] = dot(Ateb2sin,vs,`NTERM');
	ateb[3] = dot(Ateb3cos,vc,`NTERM');
	ateb[4] = dot(Ateb4sin,vs,`NTERM');
	ateb[5] = dot(Ateb5cos,vc,`NTERM');
	ateb[6] = dot(Ateb6sin,vs,`NTERM');

	return dot(ateb, vp, NPOT);
}


/********** INVERSE CALCULUS *********/

defineTable(`F1cos')
defineTable(`F2sin')
defineTable(`F3cos')
defineTable(`F4sin')
defineTable(`F5cos')
defineTable(`F6sin')

defineTable(`dQ2Lat1cos')
defineTable(`dQ2Lat2sin')
defineTable(`dQ2Lat3cos')
defineTable(`dQ2Lat4sin')
defineTable(`dQ2Lat5cos')
defineTable(`dQ2Lat6sin')

void PREFIX`'utm2geod (double x, double y, double *lat, double *lon)
{
  double vs[`NTERM'], vc[`NTERM'], vp[NPOT];
  double phi = PREFIX`'Ateb(y);
  double dq;

  vpot(vp, x/cos(phi));

  if (lon) {
    vfourier(vc, phi, cos);
    *lon =
      + dot(vc, `F1cos', `NTERM') * vp[1]
      - dot(vc, `F3cos', `NTERM') * vp[3]
      + dot(vc, `F5cos', `NTERM') * vp[5];
  }
  if (lat) {
    double aux[NPOT];

    vfourier(vs, phi, sin);
    dq =
      - dot(vs, `F2sin', `NTERM') * vp[2]
      + dot(vs, `F4sin', `NTERM') * vp[4]
      - dot(vs, `F6sin', `NTERM') * vp[6];
    vpot(vp, dq);
    aux[0] = phi;
    aux[1] = dot(`dQ2Lat1cos', vc, `NTERM');
    aux[2] = dot(`dQ2Lat2sin', vs, `NTERM');
    aux[3] = dot(`dQ2Lat3cos', vc, `NTERM');
    aux[4] = dot(`dQ2Lat4sin', vs, `NTERM');
    aux[5] = dot(`dQ2Lat5cos', vc, `NTERM');
    aux[6] = dot(`dQ2Lat6sin', vs, `NTERM');

    *lat = dot(aux,vp, NPOT);
  }
}

/* $Id: utm.c.m4,v 2.3 1998/08/17 19:02:20 luis Exp $ */
