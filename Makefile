# $Id: Makefile,v 1.5 2002/09/06 00:12:11 luis Exp $
# Author: Luis Colorado <luis.colorado@slug.ctv.es>
# Date: Mon Aug 24 16:53:10 MET DST 1998

CFLAGS=-O2 -g
LDFLAGS=-lm

.SUFFIXES: .m4 .hp48

progs = IN WE WD SA FA AM KA HO ID AA AN BR \
	BN CC CD EA EB EC ED EE EF RF HE 
modulos_m4 = $(progs:=.m4)
modulos_c = $(progs:=.c)
modulos_o = $(progs:=.o)
modulos_hp48 = $(progs:=.hp48)

all: $(progs) utm.hp48

AA_objs = AA.o main.o
AA: $(AA_objs)
	$(CC) $(LDFLAGS) -o AA $(AA_objs)

AN_objs = AN.o main.o
AN: $(AN_objs)
	$(CC) $(LDFLAGS) -o AN $(AN_objs)

BR_objs = BR.o main.o
BR: $(BR_objs)
	$(CC) $(LDFLAGS) -o BR $(BR_objs)

BN_objs = BN.o main.o
BN: $(BN_objs)
	$(CC) $(LDFLAGS) -o BN $(BN_objs)

CC_objs = CC.o main.o
CC: $(CC_objs)
	$(CC) $(LDFLAGS) -o CC $(CC_objs)

CD_objs = CD.o main.o
CD: $(CD_objs)
	$(CC) $(LDFLAGS) -o CD $(CD_objs)

EA_objs = EA.o main.o
EA: $(EA_objs)
	$(CC) $(LDFLAGS) -o EA $(EA_objs)

EB_objs = EB.o main.o
EB: $(EB_objs)
	$(CC) $(LDFLAGS) -o EB $(EB_objs)

EC_objs = EC.o main.o
EC: $(EC_objs)
	$(CC) $(LDFLAGS) -o EC $(EC_objs)

ED_objs = ED.o main.o
ED: $(ED_objs)
	$(CC) $(LDFLAGS) -o ED $(ED_objs)

EE_objs = EE.o main.o
EE: $(EE_objs)
	$(CC) $(LDFLAGS) -o EE $(EE_objs)

EF_objs = EF.o main.o
EF: $(EF_objs)
	$(CC) $(LDFLAGS) -o EF $(EF_objs)

RF_objs = RF.o main.o
RF: $(RF_objs)
	$(CC) $(LDFLAGS) -o RF $(RF_objs)

HE_objs = HE.o main.o
HE: $(HE_objs)
	$(CC) $(LDFLAGS) -o HE $(HE_objs)

IN_objs = IN.o main.o
IN: $(IN_objs)
	$(CC) $(LDFLAGS) -o IN $(IN_objs)

WE_objs = WE.o main.o
WE: $(WE_objs)
	$(CC) $(LDFLAGS) -o WE $(WE_objs)

WD_objs = WD.o main.o
WD: $(WD_objs)
	$(CC) $(LDFLAGS) -o WD $(WD_objs)

SA_objs = SA.o main.o
SA: $(SA_objs)
	$(CC) $(LDFLAGS) -o SA $(SA_objs)

FA_objs = FA.o main.o
FA: $(FA_objs)
	$(CC) $(LDFLAGS) -o FA $(FA_objs)

AM_objs = AM.o main.o
AM: $(AM_objs)
	$(CC) $(LDFLAGS) -o AM $(AM_objs)

KA_objs = KA.o main.o
KA: $(KA_objs)
	$(CC) $(LDFLAGS) -o KA $(KA_objs)

HO_objs = HO.o main.o
HO: $(HO_objs)
	$(CC) $(LDFLAGS) -o HO $(HO_objs)

ID_objs = ID.o main.o
ID: $(ID_objs)
	$(CC) $(LDFLAGS) -o ID $(ID_objs)

.m4.o:
	m4 $< utm.c.m4 >$*.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o

.m4.hp48:
	m4 $*.m4 utmtables_hp48.m4 >$*.hp48

HE.m4: genutm
	genutm -a 6378200.0 -e 0.006693421622 -k 0.9996 \
	  -c "Helmert 1906" -n HE >HE.m4

RF.m4: genutm
	genutm -a 6378137 -e 0.006694380023 -k 0.9996 \
	  -c "Geodetic Reference System 1980" -n RF >RF.m4

EF.m4: genutm
	genutm -a 6377309.613 -e 0.006637846631 -k 0.9996 \
	  -c "Everest, Pakistan" -n EF >EF.m4

EE.m4: genutm
	genutm -a 6377304.063 -e 0.006637846631 -k 0.9996 \
	  -c "Everest, W. Malaysia and Singapore 1948" -n EE >EE.m4

ED.m4: genutm
	genutm -a 6377295.664 -e 0.006637846631 -k 0.9996 \
	  -c "Everest, W. Malaysia 1969" -n ED >ED.m4

EC.m4: genutm
	genutm -a 6377301.243 -e 0.006637846631 -k 0.9996 \
	  -c "Everest, India 1956" -n EC >EC.m4

EB.m4: genutm
	genutm -a 6377298.556 -e 0.006637846631 -k 0.9996 \
	  -c "Everest, Brunei and E. Malaysia (Sabah and Sarawak)" -n EB >EB.m4

EA.m4: genutm
	genutm -a 6377276.345 -e 0.006637846631 -k 0.9996 \
	  -c "Everest, India 1830" -n EA >EA.m4

CD.m4: genutm
	genutm -a 6378249.145 -e 0.006803511283 -k 0.9996 \
	  -c "Clarke, 1.880" -n CD >CD.m4

CC.m4: genutm
	genutm -a 6378206.4 -e 0.006768657997 -k 0.9996 \
	  -c "Clarke, 1.866" -n CC >CC.m4

BN.m4: genutm
	genutm -a 6377483.865 -e 0.006674372231 -k 0.9996 \
	  -c "Bessel, 1.841, Namibia" -n BN >BN.m4

BR.m4: genutm
	genutm -a 6377397.155 -e 0.006674372231 -k 0.9996 \
	  -c "Bessel, 1.841, Ethiopia, Indonesia, Japan, and Korea" -n BR >BR.m4

AN.m4: genutm
	genutm -a 6378160.0 -e 0.006694541854 -k 0.9996 \
	  -c "Australian National" -n AN >AN.m4

AA.m4: genutm
	genutm -a 6377563.396 -e 0.006670540001 -k 0.9996 \
	  -c "Airy, 1.830" -n AA >AA.m4

IN.m4: genutm
	genutm -a 6378388.000 -e 0.0067226700223332915 -k 0.9996 \
	  -c "European/Hayford/International, 1.924" -n IN >IN.m4

WE.m4: genutm
	genutm -a 6378137.000 -e 0.00669437999014 -k 0.9996 \
	  -c "World Geodetic, 1.984" -n WE >WE.m4

WD.m4: genutm
	genutm -a 6378135.0 -e 0.006694317778 -k 0.9996 \
	  -c "World Geodetic, 1.972" -n WD >WD.m4

SA.m4: genutm
	genutm -a 6378160.0 -e 0.006694541854 -k 0.9996 \
	  -c "South American, 1.969" -n SA >SA.m4

FA.m4: genutm
	genutm -a 6378155.0 -e 0.006693421622 -k 0.9996 \
	  -c "Modified Fischer, 1.960" -n FA >FA.m4

AM.m4: genutm
	genutm -a 6377340.189 -e 0.006670540001 -k 0.9996 \
	  -c "Modified Airy" -n AM >AM.m4

KA.m4: genutm
	genutm -a 6378245.0 -e 0.006693421622 -k 0.9996 \
	  -c "Krassovsky, 1.940" -n KA >KA.m4

HO.m4: genutm
	genutm -a 6378270.0 -e 0.006722670022 -k 0.9996 \
	  -c "Hough, 1.960" -n HO >HO.m4

ID.m4: genutm
	genutm -a 6378160.0 -e 0.00669460908 -k 0.9996 \
	  -c "Indonesian, 1.974" -n ID >ID.m4

utm.hp48: utm.hp48.m4 $(modulos_hp48)
	m4 utm.hp48.m4 >utm.hp48

$(modulos_hp48): utmtables_hp48.m4

$(modulos_o): utm.c.m4

genutm_objs = genutm.o
genutm: $(genutm_objs)
	$(CC) $(LDFLAGS) -o genutm $(genutm_objs)

clean:
	rm -f $(modulos_o) $(modulos_c) $(modulos_m4) $(modulos_hp48) $(progs)
	rm -f *.o genutm utm.hp48
