/* $Id: utm.h,v 2.2 2002/09/06 00:12:11 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Thu Aug 29 21:21:22 MEST 2002
 * $Log: utm.h,v $
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

/* Number of terms used in fourier series */
#define GEO_NTERM 10
#define GEO_NPOT 7
#define GEO_K0	0.9996	/* UTM reduction factor */

struct utmparam {
	char *name; /* Name */
	char *dsc; /* Description */
	double e2; /* eccentricity, squared */
	double A; /* Semimajor axis */
	double X0, Y0, Z0; /* offsets of center of ellipsoid referred to WGS84 */
	double B; /* Minor Semiaxis */
	double AK0; /* A semiaxis * K0 (==0.9996) */
	double Ncos[GEO_NTERM]; /* N calculus */
	double Mcos[GEO_NTERM]; /* M calculus */
	double BetaPhi; /* Beta calculus */
	double Betasin[GEO_NTERM];
	double A1cos[GEO_NTERM]; /* geod2utm calculus */
	double A2sin[GEO_NTERM];
	double A3cos[GEO_NTERM];
	double A4sin[GEO_NTERM];
	double A5cos[GEO_NTERM];
	double A6sin[GEO_NTERM];
	double BetaPI; /* Ateb calculus */
	double Ateb1cos[GEO_NTERM];
	double Ateb2sin[GEO_NTERM];
	double Ateb3cos[GEO_NTERM];
	double Ateb4sin[GEO_NTERM];
	double Ateb5cos[GEO_NTERM];
	double Ateb6sin[GEO_NTERM];
	double F1cos[GEO_NTERM]; /* utm2geod calculus */
	double F2sin[GEO_NTERM];
	double F3cos[GEO_NTERM];
	double F4sin[GEO_NTERM];
	double F5cos[GEO_NTERM];
	double F6sin[GEO_NTERM];
	double dQ2Lat1cos[GEO_NTERM];
	double dQ2Lat2sin[GEO_NTERM];
	double dQ2Lat3cos[GEO_NTERM];
	double dQ2Lat4sin[GEO_NTERM];
	double dQ2Lat5cos[GEO_NTERM];
	double dQ2Lat6sin[GEO_NTERM];
}; /* struct utmparam */

#endif /* _UTM_H */
/* $Id: utm.h,v 2.2 2002/09/06 00:12:11 luis Exp $ */
