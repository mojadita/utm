/* $Id: genutm.c,v 2.2 1998/06/20 20:20:55 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Sun May 10 15:25:27 MET DST 1998
 * $Log: genutm.c,v $
 * Revision 2.2  1998/06/20 20:20:55  luis
 * *** empty log message ***
 *
 * Revision 2.1  1998/05/14 16:54:19  luis
 * Revisión con cálculo de coordenadas UTM a partir de geodésicas y
 * determinación del recuadro cienkilométrico y la zona.
 *
 * Revision 2.0  1998/05/12 18:53:09  luis
 * Versión directa completa con cálculo de coordenadas UTM a partir
 * de geodésicas.
 *
 * Revision 1.2  1998/05/11 18:45:00  luis
 * Cálculo directo completo.
 *
 */

#define IN_UTM_C

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

/* eccentricity of earth (EURO50) */
#define E2 0.0067226700223332915
#define A  6378388.0
#define K0 0.9996

#define N 1000
#define NTERM 12

double n(double l)
{
  double sl = sin(l);
  return 1.0/sqrt(1-E2*sl*sl);
}

double m(double l)
{
  double nl = n(l);
  return (1-E2)*nl*nl*nl;
}

double simpson (double(*f)(double),double a, double b, int n)
{
  double t = a;
  double dt = (b - a) / (double) n;
  double f1 = f(a);
  double f2, f3;
  double acum = 0.0;
  int i;

  for (i = 0; i < n; i++) {
  	f2 = f(t + dt/2.0);
	f3 = f(t + dt);
	acum += (f1 + 4.0*f2 + f3)*dt/6.0;
	f1 = f3;
	t += dt;
  }
  return acum;
}

double (*F)(double);
double (*SC)(double);
double I;

double F_Fourier (double x)
{
	return F(x)*SC(I*x);
}

double C_Fourier_sin (double(*f)(double), int i, int n)
{
  F = f; SC = sin; I = i;
  return simpson (F_Fourier, 0.0, 2.0*M_PI, n);
}

double C_Fourier_cos (double(*f)(double), int i, int n)
{
  F = f; SC = cos; I = i;
  if (i) return simpson (F_Fourier, 0.0, 2*M_PI, n);
  else return simpson (F_Fourier, 0.0, 2*M_PI, n)/2.0;
}


double Mcos[NTERM];  /* To determine M */

double A0sin[NTERM];

double A0(double x)  /* Beta */
{
  int i;
  double res;
  res = x*A0sin[0];
  for (i = 1; i < NTERM; i++)
    res += A0sin[i] * sin(i*x);
  return res;
}

double dA0Ncos_M(double x) /* N * cos(l) */
{ return n(x)*cos(x);
}

double A1cos[NTERM];

double A1(double x)
{
  int i;
  double res;
  res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A1cos[i] * cos(i*x);
  return res;
}

double dA1Ncos_M (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A1cos[i] * -i * sin (i*x);
  return n(x)/m(x)*cos(x)*res/2.0;
}

double A2sin[NTERM];

double A2 (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A2sin[i] * sin (i*x);
  return res;
}

double dA2Ncos_M (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A2sin[i] * i * cos (i*x);
  return n(x)*cos(x)/m(x)*res/3.0;
}

double A3cos[NTERM];

double A3 (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A3cos[i] * cos (i*x);
  return res;
}

double dA3Ncos_M (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A3cos[i] * -i * sin (i*x);
  return n(x)*cos(x)/m(x)*res/4.0;
}

double A4sin[NTERM];

double A4 (double x)
{
  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A4sin[i] * sin (i*x);
  return res;
}

double dA4Ncos_M (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A4sin[i] * i * cos (i*x);
  return n(x)*cos(x)/m(x)*res/5.0;
}

double A5cos[NTERM];

double A5 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A5cos[i] * cos (i*x);
  return res;
}

double dA5Ncos_M (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A5cos[i] * -i * sin (i*x);
  return n(x)*cos(x)/m(x)*res/6.0;
}

double A6sin[NTERM];

double A6 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += A6sin[i] * sin (i*x);
  return res;
}

double preF1 (double x)
{
  return 1.0/n(x);
}

double F1cos[NTERM];

double F1 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F1cos[i] * cos (i*x);
  return res;
}

double B1 (double x)
{ return F1(x) / cos(x);
}

double B1_(double x)
{ return 1.0 / n(x) / cos(x);
}

double dF1 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= i*F1cos[i] * sin (i*x);
    return res;
}

double preF2 (double x)
{
  return (dF1(x)*cos(x) + F1(x)*sin(x))/m(x)/2.0;
}

double F2sin[NTERM];

double F2 (double x)
{  int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F2sin[i] * sin (i*x);
  return res;
}

double B2 (double x)
{ return F2(x)/pow(cos(x), 2.0);
}

double B2_(double x)
{ return sin(x)/2.0/pow(n(x)*cos(x), 2.0);
}

double dF2 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += i*F2sin[i] * cos (i*x);
  return res;
}

double preF3 (double x)
{ return (dF2(x) * cos(x) + 2.0 * F2(x) * sin(x))/3.0/m(x);
}

double F3cos [NTERM];

double F3 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F3cos [i] * cos(i*x);
  return res;
}

double B3 (double x)
{ return F3(x)/pow(cos(x), 3.0);
}

double dF3 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= i*F3cos [i] * sin(i*x);
  return res;
}

double preF4 (double x)
{ return (dF3(x)*cos(x) + 3.0*F3(x)*sin(x))/4.0/m(x);
}

double F4sin[NTERM];

double F4 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F4sin[i] * sin (i*x);
  return res;
}

double B4 (double x)
{ return F4(x)/pow(cos(x), 4.0); }

double dF4 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += i*F4sin [i] * cos(i*x);
  return res;
}

double preF5 (double x)
{ return (dF4(x)*cos(x) + 4.0*F4(x)*sin(x))/5.0/m(x);
}

double F5cos[NTERM];

double F5 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F5cos[i] * cos (i*x);
  return res;
}

double B5 (double x)
{ return F5(x)/pow(cos(x), 5.0); }

double dF5 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res -= i*F5cos [i] * sin(i*x);
  return res;
}

double preF6 (double x)
{ return (dF5(x)*cos(x) + 5.0*F5(x)*sin(x))/6.0/m(x);
}

double F6sin[NTERM];

double F6 (double x)
{ int i;
  double res = 0.0;
  for (i = 0; i < NTERM; i++)
    res += F6sin[i] * sin (i*x);
  return res;
}

double B6 (double x)
{ return F6(x)/pow(cos(x), 6.0); }

double geod2utmX (double lat, double lon)
{
  double res = K0*A*(A1(lat)*lon - A3(lat)*pow(lon, 3.0) + A5(lat)*pow(lon, 5.0));
  return res + 500000;
}

double geod2utmY (double lat, double lon)
{
  double res = K0*A*(A0(lat) - A2(lat)*pow(lon, 2.0) + A4(lat)*pow(lon, 4.0) - A6(lat)*pow(lon, 6.0));
  if (res < 0.0) res += 10000000.0;
  return res;
}

double hms2h (double x)
{
  double deg, min;
  x = modf (x, &deg)*100.0;
  x = modf (x, &min)*100.0;
  return deg + min / 60.0 + x / 3600.0;
}

int huso (double l, double L, char *zona)
{
  static char *t1 = "CDEFGHJKLMNPQRSTUVWXYZ";
  int h = (int)((L + 180.0) / 6.0) + 1;
  sprintf (zona, "%2d%c", h, t1 [(int)((l + 80.0) / 8.0)]);
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


/* main program */
int main (int argc, char **argv)
{
	char linea [1000];
	double l, L, err;
	int i;

  printf ("Determinación de beta como serie de funciones... Primero determinamos\n");
  printf ("M como serie de funciones cos(i*phi).\n");
  for (i = 0; i < NTERM; i++) {
    Mcos[i] = (i & 1) ? 0.0 : C_Fourier_cos(m, i, N)/M_PI;
    if (!(i&1)) printf ("Mcos[%d] = %-20.17lg;\n", 
      i, Mcos[i]);
  }
  printf ("Ahora vamos a calcular los coeficientes en sin(i*phi) para la\n");
  printf ("función beta(l), integrando los de M(phi)\n");
  A0sin[0] = Mcos[0];
  printf ("Bphi = %-20.17lg --> %-20.17lg(deg.)\n",
    A0sin[0], M_PI/180.0*A0sin[0]);
  for (i = 1; i < NTERM; i++) {
    A0sin[i] = 1.0/i*Mcos[i];
    if (!(i&1)) printf ("A0sin[%d] = %-20.17lg\n", i, A0sin[i]);
  }
  printf ("Ahora calculamos A1[..] funciones en cos(i*phi)\n");
  for (i=0;i<NTERM; i++) {
    A1cos[i] = (i&1) ? C_Fourier_cos(dA0Ncos_M, i, N)/M_PI : 0.0;
    if (i&1) printf ("A1cos[%d] = %-20.17lg\n", i, A1cos[i]);
  }
  printf ("...A2[..]\n");
  A2sin[0] = 0.0;
  for (i=1; i<NTERM; i++) {
    A2sin[i] = (i&1) ? 0.0 : C_Fourier_sin(dA1Ncos_M, i, N)/M_PI;
    if (!(i&1)) printf ("A2sin[%d] = %-20.17lg\n", i, A2sin[i]);
  }
  printf ("...A3[..]\n");
  for (i=0; i<NTERM; i++) {
    A3cos[i] = (i&1) ? C_Fourier_cos(dA2Ncos_M, i, N)/M_PI : 0.0;
    if (i&1) printf ("A3cos[%d] = %-20.17lg\n", i, A3cos[i]);
  }
  printf ("...A4[..]\n");
  for (i=0; i<NTERM; i++) {
    A4sin[i] = (i&1) ? 0.0 : C_Fourier_sin(dA3Ncos_M, i, N)/M_PI;
    if (!(i&1)) printf ("A4sin[%d] = %-20.17lg\n", i, A4sin[i]);
  }
  printf ("...A5[..]\n");
  for (i=0; i<NTERM; i++) {
    A5cos[i] = (i&1) ? C_Fourier_cos(dA4Ncos_M, i, N)/M_PI : 0.0;
    if (i&1) printf ("A5cos[%d] = %-20.17lg\n", i, A5cos[i]);
  }
  printf ("...A6[..]\n");
  for (i=0; i<NTERM; i++) {
    A6sin[i] = (i&1) ? 0.0 : C_Fourier_sin(dA5Ncos_M, i, N)/M_PI;
    if (!(i&1)) printf ("A6sin[%d] = %-20.17lg\n", i, A6sin[i]);
  }
  printf ("...F1[..] (B1 = F1 / cos(phi))\n");
  for (i=0; i < NTERM; i++) {
    F1cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF1, i, N)/M_PI;
    if (!(i&1)) printf ("F1cos[%d] = %-20.17lg\n", i, F1cos[i]);
  }
  printf ("...F2[..] (B2 = F2 / cos(phi)^2.0)\n");
  for (i=0; i < NTERM; i++) {
    F2sin[i] = (i&1) ? C_Fourier_sin(preF2, i, N)/M_PI : 0.0;
    if (i&1) printf ("F2sin[%d] = %-20.17lg\n", i, F2sin[i]);
  }
  printf ("...F3[..] (B3 = F3 / cos(phi)^3.0)\n");
  for (i=0; i < NTERM; i++) {
    F3cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF3, i, N)/M_PI;
    if (!(i&1)) printf ("F3cos[%d] = %-20.17lg\n", i, F3cos[i]);
  }
  printf ("...F4[..] (B4 = F4 / cos(phi)^4.0)\n");
  for (i=0; i < NTERM; i++) {
    F4sin[i] = (i&1) ? C_Fourier_sin(preF4, i, N)/M_PI : 0.0;
    if (i&1) printf ("F4sin[%d] = %-20.17lg\n", i, F4sin[i]);
  }
  printf ("...F5[..] (B5 = F5 / cos(phi)^5.0)\n");
  for (i=0; i < NTERM; i++) {
    F5cos[i] = (i&1) ? 0.0 : C_Fourier_cos(preF5, i, N)/M_PI;
    if (!(i&1)) printf ("F5cos[%d] = %-20.17lg\n", i, F5cos[i]);
  }
  printf ("...F6[..] (B6 = F6 / cos(phi)^6.0)\n");
  for (i=0; i < NTERM; i++) {
    F6sin[i] = (i&1) ? C_Fourier_sin(preF6, i, N)/M_PI : 0.0;
    if (i&1) printf ("F6sin[%d] = %-20.17lg\n", i, F6sin[i]);
  }

  for (;;) {
	char z [10];
	int h;
	double x, y;
	printf ("##(hh.mmssss hh.mmssss)> ");
  	if (!gets(linea)) break;
	l = L = 0.0;
	sscanf (linea, "%lf%lf", &l, &L);
	l = hms2h(l); L = hms2h(L);
	h = huso (l, L, z);
	L = fmod (L + 180.0, 6.0) - 3.0;

	printf ("Huso:           %d\n", h);
	printf ("Zona:           %s\n", z);
	printf ("Lat(deg)        %-20.17lg\n", l);
	printf ("Lon(deg)        %-20.17lg\n", L);
	l *= M_PI/180.0; L *= M_PI/180.0;
	printf ("Lat:            %-20.17lg\n", l);
	printf ("Lon:            %-20.17lg\n", L);
	printf ("M:              %-20.17lg\n", A*m(l));
	printf ("N:              %-20.17lg\n", A*n(l));
#define KK 4.84813681108e-2
	printf ("TRANSFORMACIÓN DIRECTA:\n");
	printf (" A0(Beta): (I)   %20.3lf(%lg)\n", K0*A*A0(l), A0(l));
	printf (" A2:       (II)  %20.3lf(%lg)\n", -pow(KK,2.0)*K0*A*A2(l), A2(l));
	printf (" A4:       (III) %20.3lf(%lg)\n", pow(KK,4.0)*K0*A*A4(l)), A4(l);
	printf (" A6:       (A6)  %20.3lf(%lg)\n", -pow(KK,6.0)*K0*A*A6(l), A6(l));
	printf (" A1(IV):   (IV)  %20.3lf(%lg)\n", pow(KK,1.0)*K0*A*A1(l), A1(l));
	printf (" A3:       (V)   %20.3lf(%lg)\n", -pow(KK,3.0)*K0*A*A3(l), A3(l));
	printf (" A5:       (B5)  %20.3lf(%lg)\n", pow(KK,5.0)*K0*A*A5(l), A5(l));
	printf ("TRANSFORMACIÓN INVERSA:\n");
	printf (" B1:       (IX)  %20.17lf\n", B1(l));
	printf (" B1':      (ix)  %20.17lf\n", B1_(l));
	printf (" B2:       (VII) %20.17lf\n", B2(l));
	printf (" B2':      (vii) %20.17lf\n", B2_(l));
	printf (" B3:       (X)   %20.17lf\n", B3(l));
	printf (" B4:       (VIII)%20.17lf\n", B4(l));
	printf (" B5:       (E5)  %20.17lf\n", B5(l));
	printf (" B6:       (D6)  %20.17lf\n", B6(l));
	printf ("\n");
	printf ("X: %20.3lf\nY: %20.3lf\n", x = geod2utmX(l, L), y = geod2utmY(l, L));
	printf ("zona:            %s%s\n", z, zona(h, x, y));

  }
}

/* $Id: genutm.c,v 2.2 1998/06/20 20:20:55 luis Exp $ */
