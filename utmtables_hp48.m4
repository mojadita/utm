divert(-1)
changecom(,)
define(`lq',changequote([,])[changequote([,])`changequote(`,')]changequote(`,'))
define(`rq',changequote([,])[changequote([,])'changequote(`,')]changequote(`,'))
define(`id', ``'$1_$2`'')
define(`for',`ifelse(eval((`$2') <= (`$3')),1,`pushdef(`$1',`$2')`'$4`'dnl
popdef(`$1')for(`$1',eval((`$2')+1),`$3', `$4')')')
define(`defineTable',`$1
    [for(`I',0,NTERM-1,` id(`$1'ifelse(`$2',,,`$2'),I)') ]')
ifndef(`PREFIX',define(PREFIX,`'))
divert(0)dnl
lq
@ $Id: utmtables_hp48.m4,v 1.2 1998/08/24 16:33:36 luis Exp $
@ Autor: Luis Colorado <luis.colorado@@slug.ctv.es>
@ Date: Mon Aug 24 16:29:12 MET DST 1998
  DIR
    `A' A
    `K0' K0
    `E2' E2
    `BetaPI' BetaPI
    `BetaPhi' BetaPhi_deg
    defineTable(`Ncos')
    defineTable(`Mcos')
    define(`Betasin_0',0)dnl
    defineTable(`Betasin')
    defineTable(`A1cos')
    defineTable(`A2sin')
    defineTable(`A3cos')
    defineTable(`A4sin')
    defineTable(`A5cos')
    defineTable(`A6sin')
    defineTable(`Ateb1cos',`_deg')
    defineTable(`Ateb2sin',`_deg')
    defineTable(`Ateb3cos',`_deg')
    defineTable(`Ateb4sin',`_deg')
    defineTable(`Ateb5cos',`_deg')
    defineTable(`Ateb6sin',`_deg')
    defineTable(`F1cos')
    defineTable(`F2sin')
    defineTable(`F3cos')
    defineTable(`F4sin')
    defineTable(`F5cos')
    defineTable(`F6sin')
    defineTable(`dQ2Lat1cos',`_deg')
    defineTable(`dQ2Lat2sin',`_deg')
    defineTable(`dQ2Lat3cos',`_deg')
    defineTable(`dQ2Lat4sin',`_deg')
    defineTable(`dQ2Lat5cos',`_deg')
    defineTable(`dQ2Lat6sin',`_deg')
  END
@ $Id: utmtables_hp48.m4,v 1.2 1998/08/24 16:33:36 luis Exp $
rq
