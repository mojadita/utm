/* $Id: utm_ini.h,v 2.1 2002/09/06 00:12:11 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Fri Sep  6 02:07:35 MEST 2002
 * $Log: utm_ini.h,v $
 * Revision 2.1  2002/09/06 00:12:11  luis
 * Añadidos utm_ini.h para que genutm pueda calcular por tabla los parámetros y
 * utmcalc.c para los cálculos a partir de los parámetros calculados según la
 * estructura utmparam.
 *
 */

#include "utm.h"

struct utmparam wgs84table[] = {
	{ "AA", "Airy, 1.830", 0.006670540001, 6377563.396,  /* ... */ },
	{ "AM", "Modified Airy", 0.006670540001, 6377340.189,  /* ... */ },
	{ "AN", "Australian National", 0.006694541854, 6378160.0,  /* ... */ },
	{ "BN", "Bessel, 1.841, Namibia", 0.006674372231, 6377483.865,  /* ... */ },
	{ "BR", "Bessel, 1.841, Ethiopia, Indonesia, Japan, and Korea", 0.006674372231, 6377397.155,  /* ... */ },
	{ "CC", "Clarke, 1.866", 0.006768657997, 6378206.4,  /* ... */ },
	{ "CD", "Clarke, 1.880", 0.006803511283, 6378249.145,  /* ... */ },
	{ "EA", "Everest, India 1830", 0.006637846631, 6377276.345,  /* ... */ },
	{ "EB", "Everest, Brunei and E. Malaysia (Sabah and Sarawak)", 0.006637846631, 6377298.556,  /* ... */ },
	{ "EC", "Everest, India 1956", 0.006637846631, 6377301.243,  /* ... */ },
	{ "ED", "Everest, W. Malaysia 1969", 0.006637846631, 6377295.664,  /* ... */ },
	{ "EE", "Everest, W. Malaysia and Singapore 1948", 0.006637846631, 6377304.063,  /* ... */ },
	{ "EF", "Everest, Pakistan", 0.006637846631, 6377309.613,  /* ... */ },
	{ "FA", "Modified Fischer, 1.960", 0.006693421622, 6378155.0,  /* ... */ },
	{ "HE", "Helmert 1906", 0.006693421622, 6378200.0,  /* ... */ },
	{ "HO", "Hough, 1.960", 0.006722670022, 6378270.0,  /* ... */ },
	{ "ID", "Indonesian, 1.974", 0.00669460908, 6378160.0,  /* ... */ },
	{ "IN", "European/Hayford/International, 1.924", 0.0067226700223332915, 6378388.000,  /* ... */ },
	{ "KA", "Krassovsky, 1.940", 0.006693421622, 6378245.0,  /* ... */ },
	{ "RF", "Geodetic Reference System 1980", 0.006694380023, 6378137,  /* ... */ },
	{ "SA", "South American, 1.969", 0.006694541854, 6378160.0,  /* ... */ },
	{ "WD", "World Geodetic, 1.972", 0.006694317778, 6378135.0,  /* ... */ },
	{ "WE", "World Geodetic, 1.984", 0.00669437999014, 6378137.000,  /* ... */ },
	{ NULL /* ... */ },
}; /* wgs84table [] */

/* $Id: utm_ini.h,v 2.1 2002/09/06 00:12:11 luis Exp $ */
