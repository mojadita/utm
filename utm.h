/* $Id: utm.h,v 2.1 2002/09/02 08:05:19 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Thu Aug 29 21:21:22 MEST 2002
 * $Log: utm.h,v $
 * Revision 2.1  2002/09/02 08:05:19  luis
 * Added utm.h so we can use different GS in the calculus.
 *
 */

#ifndef _UTM_H
#define _UTM_H

/* Number of terms used in fourier series */
#define NTERM 8
#define NPOT 7

struct utmparam {
	char *name; /* Name */
	char *dsc; /* Description */
	double e2; /* eccentricity, squared */
	double A; /* Major Semiaxis, equatorial radius */
	double x0, y0, z0; /* offsets of center of ellipsoid referred to WGS84 */
	double B; /* Minor Semiaxis */
	double k0; /* UTM reduction constant. */
	double *Ncos; /* N calculus */
	double *Mcos; /* M calculus */
	double BetaPhi; /* Beta calculus */
	double *Betasin;
	double *A1cos; /* geod2utm calculus */
	double *A2sin;
	double *A3cos;
	double *A4sin;
	double *A5cos;
	double *A6sin;
	double BetaPI; /* Ateb calculus */
	double *Ateb1cos;
	double *Ateb2sin;
	double *Ateb3cos;
	double *Ateb4sin;
	double *Ateb5cos;
	double *Ateb6sin;
	double *F1cos; /* utm2geod calculus */
	double *F2sin;
	double *F3cos;
	double *F4sin;
	double *F5cos;
	double *F6sin;
	double *dQ2Lat1cos;
	double *dQ2Lat2sin;
	double *dQ2Lat3cos;
	double *dQ2Lat4sin;
	double *dQ2Lat5cos;
	double *dQ2Lat6sin;
}; /* struct utmparam */

#endif /* _UTM_H */
/* $Id: utm.h,v 2.1 2002/09/02 08:05:19 luis Exp $ */
