CFLAGS=-O2
LDFLAGS=-lm

.SUFFIXES: .m4 .hp48

modulos_o = IN.o WE.o WD.o SA.o FA.o AM.o KA.o HO.o ID.o
modulos_m4 = $(modulos_o:o=m4)
modulos_c = $(modulos_o:o=c)
progs = $(modulos_o:.o=)
modulos_hp48 = $(modulos_o:o=hp48)

all: $(progs) utm.hp48

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
	m4 $*.m4 utm.c.m4 >$*.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o

IN.m4: genutm
	genutm -a 6378388.000 -e 0.0067226700223332915 -k 0.9996 \
	  -c "European/Hayford/International, 1.924" >IN.m4

WE.m4: genutm
	genutm -a 6378137.000 -e 0.00669437999014 -k 0.9996 \
	  -c "World Geodetic, 1.984" >WE.m4

WD.m4: genutm
	genutm -a 6378135.0 -e 0.006694317778 -k 0.9996 \
	  -c "World Geodetic, 1.972" >WD.m4

SA.m4: genutm
	genutm -a 6378160.0 -e 0.006694541854 -k 0.9996 \
	  -c "South American, 1.969" >SA.m4

FA.m4: genutm
	genutm -a 6378155.0 -e 0.006693421622 -k 0.9996 \
	  -c "Modified Fischer, 1.960" >FA.m4

AM.m4: genutm
	genutm -a 6377340.189 -e 0.006670540001 -k 0.9996 \
	  -c "Modified Airy" >AM.m4

KA.m4: genutm
	genutm -a 6378245.0 -e 0.006693421622 -k 0.9996 \
	  -c "Krassovsky, 1.940" >KA.m4

HO.m4: genutm
	genutm -a 6378270.0 -e 0.006722670022 -k 0.9996 \
	  -c "Hough, 1.960" >HO.m4

ID.m4: genutm
	genutm -a 6378160.0 -e 0.00669460908 -k 0.9996 \
	  -c "Indonesian, 1.974" >ID.m4

utm.hp48: utm.hp48.m4 $(modulos_hp48)
	m4 utm.hp48.m4 >utm.hp48

$(modulos_hp48): utmtables_hp48.m4

$(modulos_o): utm.c.m4

.m4.hp48:
	m4 $*.m4 utmtables_hp48.m4 >$*.hp48
.m4.o:
	m4 $*.m4 $*.m4 utm.c.m4 >$*.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o

genutm_objs = genutm.o
genutm: $(genutm_objs)
	$(CC) $(LDFLAGS) -o genutm $(genutm_objs)

clean:
	rm -f $(modulos_o) $(modulos_c) $(modulos_m4) $(modulos_hp48) $(progs)
	rm -f *.o genutm utm.hp48
