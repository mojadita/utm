/* $Id: main.c,v 1.1 1998/08/24 14:20:02 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Mon Aug 10 20:15:25 MET DST 1998
 */

#define IN_MAIN_C

/* Standard include files */
#include <sys/types.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

/* constants */

/* types */

/* prototypes */
double Beta(double), m(double), n(double), Ateb(double);

/* variables */
extern double a, k0, e2;
extern char *dsc;

/* functions */

double hms2h (double x)
{
  double deg, min;
  x = modf (x, &deg)*100.0;
  x = modf (x, &min)*100.0;
  return deg + min / 60.0 + x / 3600.0;
}

double h2hms (double x)
{
  double deg, min;
  x = modf (x, &deg)*60.0;
  x = modf (x, &min)*60.0;
  return deg + min / 100.0 + x / 10000.0;
}

int huso (double l, double *L, char *zona)
{
  static char *t1 = "CDEFGHJKLMNPQRSTUVWXYZ";
  int h = (int)((*L + M_PI) / M_PI * 30.0) + 1;
  *L = fmod(*L + M_PI/2.0, M_PI/30.0) - M_PI/60.0;
  sprintf (zona, "%2d%c", h, t1 [(int)((l + M_PI/2.25) / M_PI * 22.5)]);
  return h;
}

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
}

#include <stdio.h>

#define OPCION_FUNCIONES 1
#define OPCION_G2U       2
#define OPCION_U2G       4
int opciones = 0;

/* main program */
int main (int argc, char **argv)
{
	char linea [1000];
	double l, L, err, k, d;
	int i, opt;

  while ((opt = getopt(argc, argv, "gu")) != EOF) {
    switch (opt) {
    case 'g':
      opciones |= OPCION_U2G; break;
    case 'u':
      opciones |= OPCION_G2U; break;
    default:
      fprintf (stderr, "utm: opción incorrecta\n");
      break;
    }
  }
  if (!opciones) opciones |= OPCION_G2U;

  if (opciones & OPCION_G2U)
  for (;;) {
	char z [10];
	int h;
	double x, y;
	printf ("Geodésicas->UTM (%s)\n"
		"##(hh.mmssss hh.mmssss)> ", dsc);
  	if (!gets(linea)) break;
	l = L = 0.0;
	sscanf (linea, "%lf%lf", &l, &L);

	printf ("A:              %0.17lg\n", a);
	printf ("E2:             %0.17lg\n", e2);
	printf ("K0:             %0.17lg\n", k0);
	printf ("Lat(h.mmssss):  %0.17lg\n", l);
	printf ("Lon(h.mmssss):  %0.17lg\n", L);
	l = hms2h(l); L = hms2h(L);
	printf ("Lat(deg)        %0.17lg\n", l);
	printf ("Lon(deg)        %0.17lg\n", L);
	l *= M_PI/180.0; L *= M_PI/180.0;
	printf ("Lat(rad):       %0.17lg\n", l);
	printf ("Lon(rad):       %0.17lg\n", L);
	h = huso (l, &L, z);
	printf ("Lon(huso):      %0.17lg\n", L);
	printf ("Beta:           %0.17lg\n", a*Beta(l));
	printf ("M:              %0.17lg\n", a*m(l));
	printf ("N:              %0.17lg\n", a*n(l));
	K_conv(l, L, &k, &d);
	printf ("K:              %0.17lg\n", k*k0);
	printf ("Converg:        %0.17lg\n", 180.0/M_PI*d);
	geod2utm(l, L, &x, &y);
	x *= a*k0; y *= a*k0; x += 5.0e5;
	printf ("X(m):           %0.3lf\n"
	        "Y(m):           %0.3lf\n", x, y);
	printf ("zona:           %s%s\n", z, zona(h, x, y));
	printf ("Ateb(y):        %0.17lg\n", 180.0/M_PI*Ateb(y/k0/a));

  }
  if (opciones & OPCION_U2G)
  for (;;) {
	char z [10];
	int h;
	double x, y, lat, lon;
	printf ("UTM->Geodésicas (%s)\n"
		"##(x.xxx y.yyy huso)> ", dsc);
  	if (!gets(linea)) break;
	sscanf (linea, "%lf%lf%i", &x, &y, &h);
	utm2geod((x - 500000.0)/k0/a, y/k0/a, &l, &L);
	L += ((h-30) * M_PI/30.0) - M_PI/60.0;

	printf ("Huso(antes):    %d\n", h);
	printf ("Zona(antes):    %s\n", zona(h, x, y));
	printf ("Lat(rad):       %0.17lg\n", l);
	printf ("Lon(rad):       %0.17lg\n", L);
	printf ("Lat(deg)        %0.17lg\n", l * 180.0/M_PI);
	printf ("Lon(deg)        %0.17lg\n", L * 180.0/M_PI);
	printf ("Lat(h.mmssss)   %0.17lg\n", h2hms(l * 180.0/M_PI));
	printf ("Lon(h.mmssss)   %0.17lg\n", h2hms(L * 180.0/M_PI));
	h = huso (l, &L, z);
	printf ("Huso(despues):  %d\n", h);
	printf ("Zona(despues):  %s%s\n", z, zona(h, x, y));
	printf ("Lon(huso/deg):  %0.17lg\n", h2hms(L * 180.0/M_PI));
	printf ("M:              %0.17lg\n", a*m(l));
	printf ("N:              %0.17lg\n", a*n(l));
	K_conv(l, L, &k, &d);
	printf ("K:              %0.17lg\n", k*k0);
	printf ("Converg:        %0.17lg\n", 180.0/M_PI*d);
	geod2utm(l, L, &x, &y);
	x *= a*k0; y *= a*k0; x += 500000.0;
	printf ("X:              %0.3lf\n"
	        "Y:              %0.3lf\n", x, y);
	printf ("Ateb(y):        %0.17lg\n", 180.0/M_PI*Ateb(y/k0/a));

  }
}

/* $Id: main.c,v 1.1 1998/08/24 14:20:02 luis Exp $ */
