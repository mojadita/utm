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
define(PREFIX,`')
divert(0)dnl
/* $Id: utm.c.m4,v 1.1 1998/08/17 13:43:24 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Mon Aug 10 15:54:07 MET DST 1998
 * $Log: utm.c.m4,v $
 * Revision 1.1  1998/08/17 13:43:24  luis
 * Initial revision
 *
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

static double f(double (*g)(double), double t[], double x)
{ int i;
  double res = 0.0;
  `for' (i=0; i < `NTERM'; i++)
    res += t[i] * g((double) i * x);
  return res;
}

/* N equatorial radius at point of latitude l given in terms of `A' */
defineTable(`Ncos')
double PREFIX`'n(double l)
{
  return f(cos,Ncos, l);
}

/* M meridianal radius at point of latitude l given in terms of `A' */
defineTable(`Mcos')
double PREFIX`'m(double l)
{
  return f(cos,Mcos, l);
}

/* `Beta', distance from point to equator */
double `BetaPhi' = BetaPhi;
define(`Betasin_0',0)dnl
defineTable(`Betasin')
double PREFIX`'Beta(double x)  /* Beta */
{
  return x*`BetaPhi' + f(sin,Betasin,x);
}

/* `A' functions, to calculate geod->utm transformation */
defineTable(`A1cos')
static double A1(double x)
{
  return f(cos, A1cos, x);
}

defineTable(`A2sin')
static double A2(double x)
{
  return f(sin, A2sin, x);
}

defineTable(`A3cos')
static double A3(double x)
{
  return f(cos, A3cos, x);
}

defineTable(`A4sin')
static double A4(double x)
{
  return f(sin, A4sin, x);
}

defineTable(`A5cos')
static double A5(double x)
{
  return f(cos, A5cos, x);
}

defineTable(`A6sin')
static double A6 (double x)
{
  return f(sin, A6sin, x);
}


/************** TRANSFORMACIÓN UTM -> GEODÉSICAS ****************/

/* FUNCIÓN ATEB, INVERSA DE LA FUNCIÓN BETA */
defineTable(`Ateb1cos')
static double Ateb1(double x)
{
  return f(cos, `Ateb1cos', x);
}

defineTable(`Ateb2sin')
static double Ateb2(double x)
{
  return f(sin, `Ateb2sin', x);
}

defineTable(`Ateb3cos')
static double Ateb3(double x)
{
  return f(cos, `Ateb3cos', x);
}

defineTable(`Ateb4sin')
static double Ateb4(double x)
{
  return f(sin, `Ateb4sin', x);
}

defineTable(`Ateb5cos')
static double Ateb5(double x)
{
  return f(cos, `Ateb5cos', x);
}

defineTable(`Ateb6sin')
static double Ateb6(double x)
{
  return f(sin, `Ateb6sin', x);
}

double PREFIX`BetaPI' = BetaPI;

double PREFIX`'Ateb(double x)
{
	double phi0 = x / PREFIX`BetaPI' * M_PI;
	double x0 = PREFIX`'Beta(phi0);
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


/********** INVERSE TRANSFORMATION *********/

defineTable(`F1cos')
static double B1 (double x)
{
  return f(cos, `F1cos', x) / cos(x);
}

defineTable(`F2sin')
static double B2 (double x)
{
  return f(sin, `F2sin', x) / pow(cos(x), 2.0);
}

defineTable(`F3cos')
static double B3 (double x)
{
  return f(cos, `F3cos', x) / pow(cos(x), 3.0);
}

defineTable(`F4sin')
static double B4 (double x)
{
  return f(sin, `F4sin', x) / pow(cos(x), 4.0);
}

defineTable(`F5cos')
static double B5 (double x)
{
  return f(cos, `F5cos', x) / pow(cos(x), 5.0);
}

defineTable(`F6sin')
static double B6 (double x)
{
  return f(sin, `F6sin', x) / pow(cos(x), 6.0);
}

/************ INCREMENT IN Q -> INCREMENT IN LATITUD **************/
defineTable(`dQ2Lat1cos')
static double dQ2Lat1 (double phi)
{
  return f(cos, `dQ2Lat1cos', phi);
}

defineTable(`dQ2Lat2sin')
static double dQ2Lat2 (double phi)
{
  return f(sin, `dQ2Lat2sin', phi);
}

defineTable(`dQ2Lat3cos')
static double dQ2Lat3 (double phi)
{
  return f(cos, `dQ2Lat3cos', phi);
}

defineTable(`dQ2Lat4sin')
static double dQ2Lat4 (double phi)
{
  return f(sin, `dQ2Lat4sin', phi);
}

defineTable(`dQ2Lat5cos')
static double dQ2Lat5 (double phi)
{
  return f(cos, `dQ2Lat5cos', phi);
}

defineTable(`dQ2Lat6sin')
static double dQ2Lat6 (double phi)
{
  return f(sin, `dQ2Lat6sin', phi);
}

void PREFIX`'geod2utm (double lat, double lon, double *x, double *y)
{
  if (x) {
    *x =
        A1(lat) * lon
      - A3(lat) * pow(lon, 3.0)
      + A5(lat) * pow(lon, 5.0);
  }
  if (y) {
    *y =
        Beta(lat)
      - A2(lat)*pow(lon, 2.0)
      + A4(lat)*pow(lon, 4.0)
      - A6(lat)*pow(lon, 6.0);
  }
}

void PREFIX`'utm2geod (double x, double y, double *lat, double *lon)
{
  double phi = PREFIX`'Ateb(y);
  double dq;
  if (lon) {
    *lon =
      + B1(phi) * x
      - B3(phi) * pow(x, 3.0)
      + B5(phi) * pow(x, 5.0);
  }
  if (lat) {
    dq =
      - B2(phi) * pow(x, 2.0)
      + B4(phi) * pow(x, 4.0)
      - B6(phi) * pow(x, 6.0);
    *lat = phi
      + dQ2Lat1(phi) * dq
      + dQ2Lat2(phi) * pow(dq, 2.0)
      + dQ2Lat3(phi) * pow(dq, 3.0)
      + dQ2Lat4(phi) * pow(dq, 4.0)
      + dQ2Lat5(phi) * pow(dq, 5.0)
      + dQ2Lat6(phi) * pow(dq, 6.0);
  }
}


/* $Id: utm.c.m4,v 1.1 1998/08/17 13:43:24 luis Exp $ */
