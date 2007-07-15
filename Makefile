# $Id: Makefile,v 1.15 2007/07/15 23:39:05 luis Exp $
# Author: Luis Colorado <luis.colorado@slug.ctv.es>
# Date: Mon Aug 24 16:53:10 MET DST 1998

OPTIONS=-DGEO_NTERM=16 -DGEO_NPOT=15 -DUSE_FFT #-DNiter=4096
CFLAGS=-O2 -g
LIBS=-lm
M4=m4
CC=gcc

.SUFFIXES: .m4 .c .o .hp48

progs = genutm utm pru #media

modulos = AA AM AN BN BR CC CD EA EB EC ED EE EF FA HE HO ID IN KA RF SA ST WD WE

modulos_c = $(modulos:=.c)
modulos_o = $(modulos:=.o)
modulos_m4 = $(modulos:=.m4)
modulos_hp48 = $(modulos:=.hp48)

all: $(progs) Makefile $(modulos_hp48)
modulos_c: $(modulos_c) Makefile
modulos_m4: $(modulos_m4) Makefile
modulos_o: $(modulos_o) Makefile
modulos_hp48: $(modulos_hp48) Makefile

.m4.c:
	$(M4) $(M4FLAGS) $(OPTIONS) $< utm.c.m4 >$*.c

.m4.hp48:
	$(M4) $(M4FLAGS) $(OPTIONS) $< utmtables_hp48.m4 >$*.hp48

.c.o:
	$(CC) $(CFLAGS) $(OPTIONS) -c $<

utm_objs = $(modulos_o) utmcalc.o main.o
utm: $(utm_objs) Makefile
	$(CC) $(CFLAGS) $(LDFLAGS) $(utm_objs) $(LIBS) -o utm
$(utm_objs): utm.h Makefile

media_objs = $(modulos_o) media.o utmcalc.o
media: $(media_objs) Makefile
	$(CC) $(CFLAGS) $(LDFLAGS) $(media_objs) $(LIBS) -o media

genutm_objs = genutm.o utmcalc.o fft.o
genutm: $(genutm_objs) Makefile
	$(CC) $(LDFLAGS) $(genutm_objs) $(LIBS) -o genutm
genutm.o: utm.h utm_ini.h
fft.o: fft.h mkroots.h

pru_objs = pru.o
pru: $(pru_objs) Makefile
	$(CC) $(LDFLAGS) $(pru_objs) -o pru
pru.o: utm.h Makefile

clean:
	rm -f $(modulos_o) $(modulos_c) $(modulos_m4) $(modulos_hp48) $(progs)
	rm -f *.o genutm utm.hp48

$(modulos_m4): genutm  Makefile
	genutm -g $* >$*.m4

$(modulos_c): utm.c.m4 Makefile

$(modulos_o): utm.h Makefile
$(modulos_hp48): utmtables_hp48.m4 Makefile

# HE.m4: genutm
# 	genutm -a 6378200.0 -e 0.006693421622 -k 0.9996 \
# 	  -c "Helmert 1906" -n HE >HE.m4
# 
# RF.m4: genutm
# 	genutm -a 6378137 -e 0.006694380023 -k 0.9996 \
# 	  -c "Geodetic Reference System 1980" -n RF >RF.m4
# 
# EF.m4: genutm
# 	genutm -a 6377309.613 -e 0.006637846631 -k 0.9996 \
# 	  -c "Everest, Pakistan" -n EF >EF.m4
# 
# EE.m4: genutm
# 	genutm -a 6377304.063 -e 0.006637846631 -k 0.9996 \
# 	  -c "Everest, W. Malaysia and Singapore 1948" -n EE >EE.m4
# 
# ED.m4: genutm
# 	genutm -a 6377295.664 -e 0.006637846631 -k 0.9996 \
# 	  -c "Everest, W. Malaysia 1969" -n ED >ED.m4
# 
# EC.m4: genutm
# 	genutm -a 6377301.243 -e 0.006637846631 -k 0.9996 \
# 	  -c "Everest, India 1956" -n EC >EC.m4
# 
# EB.m4: genutm
# 	genutm -a 6377298.556 -e 0.006637846631 -k 0.9996 \
# 	  -c "Everest, Brunei and E. Malaysia (Sabah and Sarawak)" -n EB >EB.m4
# 
# EA.m4: genutm
# 	genutm -a 6377276.345 -e 0.006637846631 -k 0.9996 \
# 	  -c "Everest, India 1830" -n EA >EA.m4
# 
# CD.m4: genutm
# 	genutm -a 6378249.145 -e 0.006803511283 -k 0.9996 \
# 	  -c "Clarke, 1.880" -n CD >CD.m4
# 
# CC.m4: genutm
# 	genutm -a 6378206.4 -e 0.006768657997 -k 0.9996 \
# 	  -c "Clarke, 1.866" -n CC >CC.m4
# 
# BN.m4: genutm
# 	genutm -a 6377483.865 -e 0.006674372231 -k 0.9996 \
# 	  -c "Bessel, 1.841, Namibia" -n BN >BN.m4
# 
# BR.m4: genutm
# 	genutm -a 6377397.155 -e 0.006674372231 -k 0.9996 \
# 	  -c "Bessel, 1.841, Ethiopia, Indonesia, Japan, and Korea" -n BR >BR.m4
# 
# AN.m4: genutm
# 	genutm -a 6378160.0 -e 0.006694541854 -k 0.9996 \
# 	  -c "Australian National" -n AN >AN.m4
# 
# AA.m4: genutm
# 	genutm -a 6377563.396 -e 0.006670540001 -k 0.9996 \
# 	  -c "Airy, 1.830" -n AA >AA.m4
# 
# IN.m4: genutm
# 	genutm -a 6378388.000 -e 0.0067226700223332915 -k 0.9996 \
# 	  -c "European/Hayford/International, 1.924" -n IN >IN.m4
# 
# WE.m4: genutm
# 	genutm -a 6378137.000 -e 0.00669437999014 -k 0.9996 \
# 	  -c "World Geodetic, 1.984" -n WE >WE.m4
# 
# WD.m4: genutm
# 	genutm -a 6378135.0 -e 0.006694317778 -k 0.9996 \
# 	  -c "World Geodetic, 1.972" -n WD >WD.m4
# 
# SA.m4: genutm
# 	genutm -a 6378160.0 -e 0.006694541854 -k 0.9996 \
# 	  -c "South American, 1.969" -n SA >SA.m4
# 
# FA.m4: genutm
# 	genutm -a 6378155.0 -e 0.006693421622 -k 0.9996 \
# 	  -c "Modified Fischer, 1.960" -n FA >FA.m4
# 
# AM.m4: genutm
# 	genutm -a 6377340.189 -e 0.006670540001 -k 0.9996 \
# 	  -c "Modified Airy" -n AM >AM.m4
# 
# KA.m4: genutm
# 	genutm -a 6378245.0 -e 0.006693421622 -k 0.9996 \
# 	  -c "Krassovsky, 1.940" -n KA >KA.m4
# 
# HO.m4: genutm
# 	genutm -a 6378270.0 -e 0.006722670022 -k 0.9996 \
# 	  -c "Hough, 1.960" -n HO >HO.m4
# 
# ID.m4: genutm
# 	genutm -a 6378160.0 -e 0.00669460908 -k 0.9996 \
# 	  -c "Indonesian, 1.974" -n ID >ID.m4

