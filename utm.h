/* $Id: utm.h,v 2.5 2002/10/17 09:32:59 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Thu Aug 29 21:21:22 MEST 2002
 * $Log: utm.h,v $
 * Revision 2.5  2002/10/17 09:32:59  luis
 * Borrado por error genutm.c
 *
 * Revision 2.4  2002/09/23 06:14:17  luis
 * Modified to support variable number of GEO_NPOT and GEO_NTERM.
 *
 * Revision 2.3  2002/09/17 19:58:27  luis
 * Added more precision.
 *
 * Revision 2.2  2002/09/06 00:12:11  luis
 * Añadidos utm_ini.h para que genutm pueda calcular por tabla los parámetros y
 * utmcalc.c para los cálculos a partir de los parámetros calculados según la
 * estructura utmparam.
 *
 * Revision 2.1  2002/09/02 08:05:19  luis
 * Added utm.h so we can use different GS in the calculus.
 *
 */

#ifndef _UTM_H
#define _UTM_H

#ifndef M_PI
#define M_PI 3.141592653589793238462643383
#endif

/* Number of terms used in fourier series */
#ifndef GEO_NTERM
#define GEO_NTERM 9
#endif

#ifndef GEO_NPOT
#define GEO_NPOT 7
#endif

#ifndef GEO_K0
#define GEO_K0	0.9996	/* UTM reduction factor */
#endif

#ifndef GEO_FALSE_EASTING
#define GEO_FALSE_EASTING 500000.0
#endif

#ifndef GEO_FALSE_NORTHING
#define GEO_FALSE_NORTHING	0.0
#endif

struct utmparam {
	char *name; /* Name */
	char *dsc; /* Description */
	double e2; /* eccentricity, squared */
	double a; /* Semimajor axis */
	double b; /* Minor Semiaxis */
	double ak0; /* A semiaxis * K0 (==0.9996) */
	double N[GEO_NTERM]; /* N calculus */
	double M[GEO_NTERM]; /* M calculus */
	double BetaPhi; /* Beta calculus */
	double Beta[GEO_NTERM];
	double A[GEO_NPOT][GEO_NTERM]; /* geod2utm calculus */
	double BetaPI; /* Ateb calculus */
	double Ateb[GEO_NPOT][GEO_NTERM];
	double F[GEO_NPOT][GEO_NTERM]; /* utm2geod calculus */
	double dQ2Lat[GEO_NPOT][GEO_NTERM];
}; /* struct utmparam */

struct geodsys {
	char *name;
	char *dsc;
	struct utmparam *ellip;
	double dX, dY, dZ;
};

double geo_N(struct utmparam *gs, double l);
double geo_M(struct utmparam *gs, double l);
double geo_Beta(struct utmparam *gs, double x);
void geo_geod2utm (struct utmparam *gs, double lat, double lon,
	double *x, double *y);
void geo_K_conv (struct utmparam *gs, double lat, double lon,
	double *kres, double *conv);
double geo_Ateb(struct utmparam *gs, double y);
void geo_utm2geod (struct utmparam *gs, double x, double y,
	double *lat, double *lon);

#endif /* _UTM_H */
/* $Id: utm.h,v 2.5 2002/10/17 09:32:59 luis Exp $ */
