/* $Id: main.c,v 1.9 2007/07/13 19:57:54 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * 		Luis Colorado <lc@luiscoloradosistemas.com>
 * 		Luis Colorado <luis.colorado@HispaLinux.ES>
 * Date: Mon Aug 10 20:15:25 MET DST 1998
 */

#define IN_MAIN_C

/* Standard include files */
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdarg.h>
#include <math.h>
#include "utm.h"

#define FMT_DIST	"%17.5lf"
#define FMT_NUMB	"%17.8lf"
#define FMT_ANGL	"%17.7lf"

#ifndef DEBUG
#define DEBUG 0
#endif

/* functions */

int menu(char *arg1, ...)
{
	va_list p;
	int n;
	char buffer[100];
	int opt;
	char *str;

	do {
		va_start(p, arg1);
		str = arg1;
		n = 1;
		do {
			printf("%2d.- %s\n", n, str);
			str = va_arg(p, char *);
			n++;
		} while (str);
		va_end(p);
		printf("Pulse una opción (0 para salir) > ");
		if (!fgets(buffer, sizeof buffer, stdin))
			return 0;
	} while ((sscanf(buffer, "%d", &opt) != 1) || (opt < 0) || (opt > n));
	return opt;
} /* menu */


extern struct utmparam
	geo_IN, geo_WE, geo_WD, geo_SA, geo_FA, geo_AM,
	geo_KA, geo_HO, geo_ID, geo_AA, geo_AN, geo_BR,
	geo_BN, geo_CC, geo_CD, geo_EA, geo_EB, geo_EC,
	geo_ED, geo_EE, geo_EF, geo_RF, geo_HE, geo_ST; 

struct utmparam *tabla[] = {
	&geo_IN, &geo_WE, &geo_WD, &geo_SA, &geo_FA, &geo_AM,
	&geo_KA, &geo_HO, &geo_ID, &geo_AA, &geo_AN, &geo_BR,
	&geo_BN, &geo_CC, &geo_CD, &geo_EA, &geo_EB, &geo_EC,
	&geo_ED, &geo_EE, &geo_EF, &geo_RF, &geo_HE, &geo_ST,
	NULL
};

struct utmparam *lookup(char *name)
{
	struct utmparam **p;

	for (p = tabla; *p; p++)
		if (!strcmp((*p)->name,name))
			break;

	return *p;
} /* lookup */

double hms2h (double x)
{
  double deg, min, sec;

#if DEBUG
  printf("hms2h: x = %lg\n", x);
#endif
  x = modf(x, &deg)*100.0;
  if (x > 80.0) x *= 0.60;
#if DEBUG
  printf("hms2h: x = %lg, deg = %lg\n", x, deg);
#endif
  x = modf(x, &min)*100.0;
  if (x > 80.0) x *= 0.60;
#if DEBUG
  printf("hms2h: x = %lg, min = %lg\n", x, min);
#endif

  x = deg + min/60.0 + x/3600.0;
#if DEBUG
  printf("hms2h: --> x = %lg\n", x);
#endif
  return x;
} /* hms2h */

double h2hms (double x)
{
  double deg, min;

#if DEBUG
  printf("h2hms: x = %lg\n", x);
#endif
  x = modf(x, &deg)*60.0;
#if DEBUG
  printf("h2hms: x = %lg, deg = %lg\n", x, deg);
#endif
  x = modf(x, &min)*60.0;
#if DEBUG
  printf("h2hms: x = %lg, min = %lg\n", x, min);
#endif

  x = deg + min/100.0 + x/10000.0;
#if DEBUG
  printf("h2hms: --> x = %lg\n", x);
#endif
  return x;
} /* h2hms */

int huso (double l, double *L, char *zona)
{
  static char *t1 = "CDEFGHJKLMNPQRSTUVWXYZ";
  int h = (int)((*L + M_PI) / M_PI * 30.0);

  *L = *L + M_PI - h*M_PI/30.0 - M_PI/60.0;
  h++;
  sprintf (zona, "%02d %c", h, t1 [(int)((l + M_PI/2.25) / M_PI * 22.5)]);

  return h;
} /* huso */

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
} /* zona */


/* main program */
int main (int argc, char **argv)
{
	char linea [1000];
	double l, L, err, k, d;
	int i, opt;
	struct utmparam *sg = &geo_IN;

  while ((opt = getopt(argc, argv, "g:l")) != EOF) {
    switch (opt) {
    case 'g':
		sg = lookup(optarg);
		if (!sg) {
			fprintf(stderr, "utm: sistema geodésico incorrecto: %s\n", optarg);
			exit(1);
		}
      break;
	case 'l': {
		struct utmparam **i;
		for (i = tabla; *i; i++) {
			printf("[%s]: %s\n", i[0]->name, i[0]->dsc);
		}
		exit(0);
	}
    default:
      fprintf (stderr, "utm: opción incorrecta\n");
      break;
    } /* switch */
  } /* while */

  printf("PROGRAMA PARA CONVERSIÓN DE DATOS UTM/GEODESICAS\n");
  printf("(C) 2002 Luis.Colorado@hispalinux.es\n");
  printf("Sistema Geodésico: [%s] %s\n", sg->name, sg->dsc);
  printf ("A:              "FMT_DIST"\n", sg->a);
  printf ("E2:             "FMT_NUMB"\n", sg->e2);

  while (opt = menu(
		"GEODESICAS -> UTM",
  		"UTM -> GEODESICAS",
		"GEODESICAS -> UTM (dentro del huso)",
  		"UTM -> GEODESICAS (dentro del huso)",
		NULL))
  {
  	switch(opt) {
	case 1: { /* geodesicas -> utm */
		char z [10];
		int h;
		double x, y;
	
		printf ("Geodésicas->UTM (%s)\n"
			"##(hh.mmssss hh.mmssss)> ", sg->dsc);
	  	if (!fgets(linea, sizeof linea, stdin)) break;

		l = L = 0.0;
		sscanf (linea, "%lf%lf", &l, &L);
	
		printf ("Lat(h.mmssss):  "FMT_ANGL"\n", l);
		printf ("Lon(h.mmssss):  "FMT_ANGL"\n", L);
		l = hms2h(l); L = hms2h(L);
		printf ("Lat(deg)     :  "FMT_ANGL"\n", l);
		printf ("Lon(deg)     :  "FMT_ANGL"\n", L);
		l *= M_PI/180.0; L *= M_PI/180.0;
		h = huso (l, &L, z);
		printf ("Lon(huso/deg):  "FMT_ANGL"\n", L*180.0/M_PI);
		printf ("Beta         :  "FMT_DIST"\n", geo_Beta(sg, l));
		printf ("M            :  "FMT_DIST"\n", geo_M(sg, l));
		printf ("N            :  "FMT_DIST"\n", geo_N(sg, l));

		geo_K_conv(sg, l, L, &k, &d);
		printf ("K            :  "FMT_NUMB"\n", k);
		printf ("Converg.     :  "FMT_ANGL"\n", 180.0/M_PI*d);
		printf ("Conv(R.)     :  "FMT_DIST"\n", geo_N(sg, l)/tan(l));

		geo_geod2utm(sg, l, L, &x, &y);
		if (y < 0.0) y += 1.0e7;
		printf ("X(m)         :  "FMT_DIST"\n"
		        "Y(m)         :  "FMT_DIST"\n", x, y);
		printf ("MGRS         :  %s%s %05d %05d\n", z, zona(h, x, y),
			(int)(x) % 100000, (int)(y) % 100000);
		break;
	}
	case 2: { /* utm->geodesicas */
		char z [10];
		int h;
		double x, y, lat, lon;

		printf ("UTM->Geodésicas (%s)\n"
			"##(x.xxx y.yyy huso)> ", sg->dsc);
	  	if (!fgets(linea, sizeof linea, stdin))
			break;
		sscanf (linea, "%lf%lf%i", &x, &y, &h);

		geo_utm2geod(sg, x, y, &l, &L);
		L += ((h-30) * M_PI/30.0) - M_PI/60.0;
	
		printf ("Lat(h.mmssss):  "FMT_ANGL"\n", h2hms(l * 180.0/M_PI));
		printf ("Lon(h.mmssss):  "FMT_ANGL"\n", h2hms(L * 180.0/M_PI));
		printf ("Lat(deg)     :  "FMT_ANGL"\n", l * 180.0/M_PI);
		printf ("Lon(deg)     :  "FMT_ANGL"\n", L * 180.0/M_PI);

		h = huso (l, &L, z);
		printf ("Huso         :  %d\n", h);

		printf ("Beta         :  "FMT_DIST"\n", geo_Beta(sg, l));
		printf ("M            :  "FMT_DIST"\n", geo_M(sg, l));
		printf ("N            :  "FMT_DIST"\n", geo_N(sg, l));

		geo_K_conv(sg, l, L, &k, &d);
		printf ("K            :  "FMT_NUMB"\n", k);
		printf ("Converg.     :  "FMT_ANGL"\n", 180.0/M_PI*d);
		printf ("Conv(R.)     :  "FMT_DIST"\n", geo_N(sg, l)/tan(l));

		geo_geod2utm(sg, l, L, &x, &y);
		if (y < 0.0) y += 1.0e7;
		printf ("X(m)         :  "FMT_DIST"\n"
		        "Y(m)         :  "FMT_DIST"\n", x, y);
		printf ("MGRS         :  %s%s %05d %05d\n", z, zona(h, x, y),
			(int)(x) % 100000, (int)(y) % 100000);

		break;
	}
	case 3: { /* geodesicas -> utm  (dentro del huso) */
		char z [10];
		int h;
		double x, y;
	
		printf ("Geodésicas->UTM (%s)\n"
			"##(hh.mmssss hh.mmssss)> ", sg->dsc);
	  	if (!fgets(linea, sizeof linea, stdin)) break;

		l = L = 0.0;
		sscanf (linea, "%lf%lf", &l, &L);
	
		printf ("Lat(h.mmssss):  "FMT_ANGL"\n", l);
		printf ("Lon(h.mmssss):  "FMT_ANGL"\n", L);
		l = hms2h(l); L = hms2h(L);
		printf ("Lat(deg)     :  "FMT_ANGL"\n", l);
		printf ("Lon(deg)     :  "FMT_ANGL"\n", L);
		l *= M_PI/180.0; L *= M_PI/180.0;
		printf ("Beta         :  "FMT_DIST"\n", geo_Beta(sg, l));
		printf ("M            :  "FMT_DIST"\n", geo_M(sg, l));
		printf ("N            :  "FMT_DIST"\n", geo_N(sg, l));

		geo_K_conv(sg, l, L, &k, &d);
		printf ("K            :  "FMT_NUMB"\n", k);
		printf ("Converg.     :  "FMT_ANGL"\n", 180.0/M_PI*d);
		printf ("Conv(R.)     :  "FMT_DIST"\n", geo_N(sg, l)/tan(l));

		geo_geod2utm(sg, l, L, &x, &y);
		if (y < 0.0) y += 1.0e7;
		printf ("X(m)         :  "FMT_DIST"\n"
		        "Y(m)         :  "FMT_DIST"\n", x, y);
		break;
	}
	case 4: { /* utm->geodesicas (dentro del huso) */
		char z [10];
		int h;
		double x, y, lat, lon;

		printf ("UTM->Geodésicas (%s)\n"
			"##(x.xxx y.yyy)> ", sg->dsc);
	  	if (!fgets(linea, sizeof linea, stdin))
			break;
		sscanf (linea, "%lf%lf", &x, &y);

		geo_utm2geod(sg, x, y, &l, &L);
	
		printf ("Lat(h.mmssss):  "FMT_ANGL"\n", h2hms(l * 180.0/M_PI));
		printf ("Lon(h.mmssss):  "FMT_ANGL"\n", h2hms(L * 180.0/M_PI));
		printf ("Lat(deg)     :  "FMT_ANGL"\n", l * 180.0/M_PI);
		printf ("Lon(deg)     :  "FMT_ANGL"\n", L * 180.0/M_PI);

		printf ("Beta         :  "FMT_DIST"\n", geo_Beta(sg, l));
		printf ("M            :  "FMT_DIST"\n", geo_M(sg, l));
		printf ("N            :  "FMT_DIST"\n", geo_N(sg, l));

		geo_K_conv(sg, l, L, &k, &d);
		printf ("K            :  "FMT_NUMB"\n", k);
		printf ("Converg.     :  "FMT_ANGL"\n", 180.0/M_PI*d);
		printf ("Conv(R.)     :  "FMT_DIST"\n", geo_N(sg, l)/tan(l));

		geo_geod2utm(sg, l, L, &x, &y);
		if (y < 0.0) y += 1.0e7;
		printf ("X(m)         :  "FMT_DIST"\n"
		        "Y(m)         :  "FMT_DIST"\n", x, y);
		break;
	}
	} /* switch */
  } /* while */

} /* main */

/* $Id: main.c,v 1.9 2007/07/13 19:57:54 luis Exp $ */
