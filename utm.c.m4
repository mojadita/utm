divert(-1)
changecom(,)
define(`lq',changequote([,])[changequote([,])`changequote(`,')]changequote(`,'))
define(`rq',changequote([,])[changequote([,])'changequote(`,')]changequote(`,'))
define(`id', ``'$1_$2`'')
define(`for',`ifelse(eval((`$2') <= (`$3')),1,`pushdef(`$1',`$2')`'$4`'dnl
popdef(`$1')for(`$1',eval((`$2')+1),`$3', `$4')')')
define(`defineTable',`{
for(`I',0,NTERM-1,`      /* `#'I */ id(`$1',I),
')dnl
    } /* $1 */dnl
')

ifdef(`X0',,`define(`X0',`0.0')')dnl
ifdef(`Y0',,`define(`Y0',`0.0')')dnl
ifdef(`Z0',,`define(`Z0',`0.0')')dnl

divert(0)dnl
`/* $Id: utm.c.m4,v 2.7 2002/09/17 20:18:44 luis Exp $'
` * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>'
` * Date: Mon Aug 10 15:54:07 MET DST 1998'
 * $Log: utm.c.m4,v $
 * Revision 2.7  2002/09/17 20:18:44  luis
 * There was an unescaped `for' in a cvs comment.
 *
 * Revision 2.6  2002/09/17 20:16:47  luis
 * Preparing to use [GEO_NPOT][GEO_NTERM] matrices, instead of vectors.
 *
 * Revision 2.5  2002/09/17 19:58:27  luis
 * Added more precision.
 *
 * Revision 2.4  2002/09/06 00:12:11  luis
 * `Añadidos' utm_ini.h para que genutm pueda calcular por tabla los parámetros y
 * utmcalc.c para los cálculos a partir de los parámetros calculados según la
 * estructura utmparam.
 *
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
 * Inclusion of Linear Deformation Modulus K and meridian convergence `for'
 * ease of distance and azimuth calculus.
 *
 * Revision 1.1  1998/08/17 13:43:24  luis
 * Initial revision
 *
 * CAUTION: THIS FILE HAS BEEN GENERATED AUTOMATICALLY FROM CONSTANTS
 * CALCULUS AND M4 TEMPLATE FILE.
 * EDIT ONLY IF YOU ARE UNABLE TO GENERATE IT AUTOMATICALLY.
 */
`'
`#define IN_UTM_C'
`'
`#include "utm.h"'
`'
`struct utmparam geo_'NAME` = {'
`    /* name */ "'NAME`",'
`    /* dsc */ "'ELLIPSOID`",'
`	/* e2 */ 'E2`,'
`    /* A */ 'A`,'
`    /* B */ 'B`,'
`    /* AK0 */ 'AK0`,'
`    /* Ncos */ 'defineTable(Ncos)`,'
`    /* Mcos */ 'defineTable(Mcos)`,'
`    /* BetaPhi */ 'BetaPhi`,'
`    /* Betasin */ 'defineTable(Betasin)`,'
`    /* A1cos */ 'defineTable(A1cos)`,'
`    /* A2sin */ 'defineTable(A2sin)`,'
`    /* A3cos */ 'defineTable(A3cos)`,'
`    /* A4sin */ 'defineTable(A4sin)`,'
`    /* A5cos */ 'defineTable(A5cos)`,'
`    /* A6sin */ 'defineTable(A6sin)`,'
`    /* A7cos */ 'defineTable(A7cos)`,'
`    /* A8sin */ 'defineTable(A8sin)`,'
`    /* BetaPI */ 'BetaPI`,'
`    /* Ateb1cos */ 'defineTable(Ateb1cos)`,'
`    /* Ateb2sin */ 'defineTable(Ateb2sin)`,'
`    /* Ateb3cos */ 'defineTable(Ateb3cos)`,'
`    /* Ateb4sin */ 'defineTable(Ateb4sin)`,'
`    /* Ateb5cos */ 'defineTable(Ateb5cos)`,'
`    /* Ateb6sin */ 'defineTable(Ateb6sin)`,'
`    /* Ateb7cos */ 'defineTable(Ateb7cos)`,'
`    /* Ateb8sin */ 'defineTable(Ateb8sin)`,'
`    /* F1cos */ 'defineTable(F1cos)`,'
`    /* F2sin */ 'defineTable(F2sin)`,'
`    /* F3cos */ 'defineTable(F3cos)`,'
`    /* F4sin */ 'defineTable(F4sin)`,'
`    /* F5cos */ 'defineTable(F5cos)`,'
`    /* F6sin */ 'defineTable(F6sin)`,'
`    /* F7cos */ 'defineTable(F7cos)`,'
`    /* F8sin */ 'defineTable(F8sin)`,'
`    /* dQ2Lat1cos */ 'defineTable(dQ2Lat1cos)`,'
`    /* dQ2Lat2sin */ 'defineTable(dQ2Lat2sin)`,'
`    /* dQ2Lat3cos */ 'defineTable(dQ2Lat3cos)`,'
`    /* dQ2Lat4sin */ 'defineTable(dQ2Lat4sin)`,'
`    /* dQ2Lat5cos */ 'defineTable(dQ2Lat5cos)`,'
`    /* dQ2Lat6sin */ 'defineTable(dQ2Lat6sin)`,'
`    /* dQ2Lat7cos */ 'defineTable(dQ2Lat7cos)`,'
`    /* dQ2Lat8sin */ 'defineTable(dQ2Lat8sin)`'
`}; /* s_'NAME` */'
`'
`/* $Id: utm.c.m4,v 2.7 2002/09/17 20:18:44 luis Exp $ */'
