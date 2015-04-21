divert(-1)
changecom(,)
define(`lq',changequote([,])[changequote([,])`changequote(`,')]changequote(`,'))
define(`rq',changequote([,])[changequote([,])'changequote(`,')]changequote(`,'))
define(`id', ``'$1_$2`'')
define(`for',`ifelse(eval((`$2') <= (`$3')),1,`pushdef(`$1',`$2')`'$4`'dnl
popdef(`$1')for(`$1',eval((`$2')+1),`$3', `$4')')')
ifndef(`PREFIX',define(PREFIX,`'))
divert(0)dnl
`%%HP: T(3)A(D)F(.);
@ $Id: utm.hp48.m4,v 1.5 2002/09/23 06:14:17 luis Exp $
@ Author: Luis Colorado <luis.colorado@@slug.ctv.es>
@ Date: Mon Aug 24 16:26:38 MET DST 1998
DIR
  U\->G
    \<< U\-> U\->Gint \>>

  G\->U
    \<< G\->Uint \->U \>>

  Beta
    \<< Betaint A * \>>

  Ateb
    \<< A / Atebint \>>

  N
    \<< Nint A * \>>

  M
    \<< Mint A * \>>

  K
    \<< \-> l L
      \<< l L dA
        l Nint / l COS / K0 *
        DUP ABS "K" \->TAG
        SWAP ARG "\Gd" \->TAG
      \>>
    \>>

  D
    \<< \-> P1 P2
      \<< P1 P2 - ABS
        P2 P1 DUP2 + 2 / U\->G K DROP
	INV 4 * SWAP U\->G K DROP INV +
	SWAP U\->G K DROP INV + 6 / *
      \>>
    \>>

  \->U
    \<< A * K0 * OFFSET + \>>

  U\->
    \<< OFFSET - A / K0 / \>>

  U\->Gint
    \<< C\->R 0 SWAP R\->C \-> Z
      \<< Atebint \-> L
        \<< L Vsin L Vcos \-> VS VC
          \<<
	    0
'for(`I', 1, GEO_NTERM, `dnl
        ifelse(eval(I & 1),1,
            `F`'defn(`I')cos VC',
            `F`'defn(`I')sin VS') DOT
')dnl
`	    \->V
	    Z L COS / Vpot DOT
	    C\->R R\->D \-> dQ LON
            \<<
	          L
'for(`I', 1, GEO_NTERM, `dnl
              ifelse(eval(I & 1),1,
                    `dQ2Lat`'defn(`I')cos VC',
                    `dQ2Lat`'defn(`I')sin VS') DOT
')dnl
`	      \->V
	      dQ Vpot DOT "Lat" \->TAG
	      LON "Lon" \->TAG
            \>>
          \>>
        \>>
      \>>
    \>>
	
  G\->Uint
    \<< D\->R \-> l L
      \<< l Vsin l Vcos 0 L R\->C Vpot \-> VS VC VP
        \<<
	      l Betaint
'for(`I', 1, GEO_NTERM, `dnl
          ifelse(eval(I & 1),1,
                `VC A'defn(`I')`cos',
                `VS A'defn(`I')`sin') DOT
')dnl
`          \->V VP DOT
        \>>
      \>>
    \>>

  Betaint
    \<< \-> L
      \<< L BetaPhi * L Vsin Betav DOT + \>>
    \>>

  Atebint
    \<< \-> D
      \<< D BetaPI / 180 * \-> L
        \<< L Vsin L Vcos \-> VS VC
          \<< L
'for(`I', 1, GEO_NTERM,`dnl
            ifelse(eval(I & 1),1,
              `Ateb'defn(`I')`cos VC',
              `Ateb'defn(`I')`sin VS') DOT
')dnl
`           \->V D L Betaint - Vpot DOT
          \>>
        \>>
      \>>
    \>>

  Mint
    \<< \-> L
      \<< L Vcos Mcos DOT \>>
    \>>

  Nint
    \<< \-> L
      \<< L Vcos Ncos DOT \>>
    \>>

  dA
    \<< D\->R \-> l L
      \<< l Vsin l Vcos 0 L R\->C Vpot
        \-> VS VC VP
        \<<
'for(`I', 1, GEO_NTERM,`dnl
          ifelse(eval(I & 1),1,
              `A'defn(`I')`cos VC',
              `A'defn(`I')`sin VS') DOT`'ifelse(I,1,,` defn(`I') *')
')dnl
`         0 \->V VP DOT
        \>>
      \>>
    \>>

  \->V
    \<< 'GEO_NTERM` ROW\-> \>>

  CST { U\->G G\->U Beta Ateb N M
        D K A E2 K0 { "" \<< \>> }
	UTM
        EURO50
        WGS84
        WGS72
        HOUGH60
        FISCHERMOD
        AIRYMOD
        KRASSOVSKY40
        HOUGH60
        INDONESIAN74
        AIRY1830
        AUSTRALIAN
        BESSEL1841
        CLARKE1866
        CLARKE1880
        EVEREST1956
        GRS80
        HELMERT1906
    }

  Vpot
    \<< \-> X
      \<< 1 1 'GEO_NPOT`
        START DUP X * NEXT 'GEO_NPOT` ROW\->
      \>>
    \>>

  Vsin
    \<< \-> L
      \<<
        0 'eval(GEO_NTERM `- 1')` FOR I
          'rq`SIN(I*L)'rq` EVAL
        NEXT
	'GEO_NTERM` ROW\->
      \>>
    \>>

  Vcos
    \<< \-> L
      \<<
        0 'eval(GEO_NTERM `- 1')` FOR I
	      'rq`COS(I*L)'rq` EVAL
        NEXT
	'GEO_NTERM` ROW\->
      \>>
    \>>

  OFFSET (0,5E5)
'
  EURO50
include(`IN.hp48')
  WGS84
include(`WE.hp48')
  WGS72
include(`WD.hp48')
  HOUGH60
include(`SA.hp48')
  FISCHERMOD
include(`FA.hp48')
  AIRYMOD
include(`AM.hp48')
  KRASSOVSKY40
include(`KA.hp48')
  HOUGH60
include(`HO.hp48')
  INDONESIAN74
include(`ID.hp48')
  AIRY1830
include(`AA.hp48')
  AUSTRALIAN
include(`AN.hp48')
  BESSEL1841
include(`BR.hp48')
  CLARKE1866
include(`CC.hp48')
  CLARKE1880
include(`CD.hp48')
  EVEREST1956
include(`EC.hp48')
  GRS80
include(`RF.hp48')
  HELMERT1906
include(`HE.hp48')

END
@ $Id: utm.hp48.m4,v 1.5 2002/09/23 06:14:17 luis Exp $'
