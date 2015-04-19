divert(-1)
changecom(,)
define(`lq',changequote([,])[changequote([,])`changequote(`,')]changequote(`,'))
define(`rq',changequote([,])[changequote([,])'changequote(`,')]changequote(`,'))
define(`id', ``'$1_$2`'')
define(`for',`ifelse(eval((`$2') <= (`$3')),1,`pushdef(`$1',`$2')`'$4`'dnl
popdef(`$1')for(`$1',eval((`$2')+1),`$3', `$4')')')

define(`defineTable',`	[for(`I',0,GEO_NTERM-1,` id(`$1',I) ') ]
')
define(`defineTable2',`
@ `$1'
[for(`I',0,GEO_NPOT-1,` defineTable(id(`$1',I))')]
')
	
ifndef(`PREFIX',define(PREFIX,`'))

divert(0)dnl
@ $Id: utmtables_hp48.m4,v 1.3 2002/09/23 06:14:17 luis Exp $
@ Autor: Luis Colorado <luis.colorado@@slug.ctv.es>
@ Date: Mon Aug 24 16:29:12 MET DST 1998
  {
    `@NAME@'         "NAME"
    `@ELLIPSOID@'    "ELLIPSOID"
    `@A@'            Aaxis
    `@K0@'           0.9996
    `@E2@'           E2
    `@BetaPI@'       BetaPI
    `@BetaPhi@'      BetaPhi_deg
    `@Nv@'           defineTable(`N')
    `@Mv@'           defineTable(`M')
    `@Betav@'        defineTable(`Beta')
    `@Av@'           defineTable2(`A')
    `@Atebv@'        defineTable2(`Ateb_deg')
    `@Fv@'           defineTable2(`F')
    `@dQ2Latv@'      defineTable2(`dQ2Lat_deg')
  }
@ $Id: utmtables_hp48.m4,v 1.3 2002/09/23 06:14:17 luis Exp $
