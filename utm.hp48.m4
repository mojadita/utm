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
	    F1cos VC DOT
	    F2sin VS DOT
	    F3cos VC DOT
	    F4sin VS DOT
	    F5cos VC DOT
	    F6sin VS DOT
	    \->V
	    Z L COS / Vpot DOT
	    C\->R R\->D \-> dQ LON
            \<<
	      L
	      dQ2Lat1cos VC DOT
	      dQ2Lat2sin VS DOT
	      dQ2Lat3cos VC DOT
	      dQ2Lat4sin VS DOT
	      dQ2Lat5cos VC DOT
	      dQ2Lat6sin VS DOT
	      \->V
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
	  VC A1cos DOT
	  VS A2sin DOT
	  VC A3cos DOT
	  VS A4sin DOT
	  VC A5cos DOT
	  VS A6sin DOT \->V VP DOT
        \>>
      \>>
    \>>

  Betaint
    \<< \-> L
      \<< L BetaPhi * L Vsin Betasin DOT + \>>
    \>>

  Atebint
    \<< \-> D
      \<< D BetaPI / 180 * \-> L
        \<< L Vsin L Vcos \-> VS VC
          \<< L
	    Ateb1cos VC DOT
	    Ateb2sin VS DOT
	    Ateb3cos VC DOT
	    Ateb4sin VS DOT
	    Ateb5cos VC DOT
	    Ateb6sin VS DOT
	    \->V D L Betaint - Vpot DOT
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
      \<< l Vsin l Vcos 0 L R\->C Vpot \-> VS VC VP
        \<<
	  A1cos VC DOT
	  A2sin VS DOT 2 *
	  A3cos VC DOT 3 *
	  A4sin VS DOT 4 *
	  A5cos VC DOT 5 *
	  A6sin VS DOT 6 *
	  0 \->V VP DOT
        \>>
      \>>
    \>>

  \->V
    \<< 7 ROW\-> \>>

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
      \<< 1 1 6
        START DUP X * NEXT 7 ROW\->
      \>>
    \>>

  Vsin
    \<< \-> L
      \<<
        0 7 FOR I
          'rq`SIN(I*L)'rq` EVAL
        NEXT
	8 ROW\->
      \>>
    \>>

  Vcos
    \<< \-> L
      \<<
        0 7 FOR I
	  'rq`COS(I*L)'rq` EVAL
        NEXT
	8 ROW\->
      \>>
    \>>

  OFFSET (0,5E5)
'
  EURO50 include(`IN.hp48')dnl
  WGS84 include(`WE.hp48')dnl
  WGS72 include(`WD.hp48')dnl
  HOUGH60 include(`SA.hp48')dnl
  FISCHERMOD include(`FA.hp48')dnl
  AIRYMOD include(`AM.hp48')dnl
  KRASSOVSKY40 include(`KA.hp48')dnl
  HOUGH60 include(`HO.hp48')dnl
  INDONESIAN74 include(`ID.hp48')dnl
  AIRY1830 include(`AA.hp48')dnl
  AUSTRALIAN include(`AN.hp48')dnl
  BESSEL1841 include(`BR.hp48')dnl
  CLARKE1866 include(`CC.hp48')dnl
  CLARKE1880 include(`CD.hp48')dnl
  EVEREST1956 include(`EC.hp48')dnl
  GRS80 include(`RF.hp48')dnl
  HELMERT1906 include(`HE.hp48')dnl

END
@ $Id: utm.hp48.m4,v 1.5 2002/09/23 06:14:17 luis Exp $
