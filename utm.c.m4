divert(-1)
changecom(,)
define(`lq',changequote([,])[changequote([,])`changequote(`,')]changequote(`,'))
define(`rq',changequote([,])[changequote([,])'changequote(`,')]changequote(`,'))
define(`id', ``'$1_$2`'')
define(`for',`ifelse(eval((`$2') <= (`$3')),1,`pushdef(`$1',`$2')`'$4`'dnl
popdef(`$1')for(`$1',eval((`$2')+1),`$3', `$4')')')
define(`defineTable',`	/* `$1' */ {
for(`I',0,GEO_NTERM-1,`      /* `#'I */ id(`$1',I),
')dnl
	}, /* $1 */
')
define(`defineTable2',`  /* `$1' */ {
for(`I',0,GEO_NPOT-1,`defineTable(id(`$1', I))')
  }, /* $1 */

')

divert(0)dnl
`/* $Id: utm.c.m4,v 2.8 2002/09/23 06:14:17 luis Exp $'
` * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>'
` * Date: Mon Aug 10 15:54:07 MET DST 1998'
 * $Log: utm.c.m4,v $
 * Revision 2.8  2002/09/23 06:14:17  luis
 * Modified to support variable number of GEO_NPOT and GEO_NTERM.
 *
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
`	 /* e2 */ 'E2`,'
`    /* a */ 'Aaxis`,'
`    /* b */ 'Baxis`,'
`    /* ak0 */ 'ak0`,'
defineTable(N)
defineTable(M)
`	 /* BetaPhi */ 'BetaPhi`,'
defineTable(Beta)
defineTable2(A)
`    /* BetaPI */ 'BetaPI`,'
defineTable2(Ateb)
defineTable2(F)
defineTable2(dQ2Lat)
`}; /* s_'NAME` */'
`'
`/* $Id: utm.c.m4,v 2.8 2002/09/23 06:14:17 luis Exp $ */'
