CFLAGS=-O2
LDFLAGS=-lm

.SUFFIXES: .m4 

all: hayford24 wgs84 wgs72 southam69 modfischer60 modairy krassovsky40 houg60 indo74 utm.hp48

hayford24_objs = hayford24.o main.o
hayford24: $(hayford24_objs)
	$(CC) $(LDFLAGS) -o hayford24 $(hayford24_objs)

wgs84_objs = wgs84.o main.o
wgs84: $(wgs84_objs)
	$(CC) $(LDFLAGS) -o wgs84 $(wgs84_objs)

wgs72_objs = wgs72.o main.o
wgs72: $(wgs72_objs)
	$(CC) $(LDFLAGS) -o wgs72 $(wgs72_objs)

southam69_objs = southam69.o main.o
southam69: $(southam69_objs)
	$(CC) $(LDFLAGS) -o southam69 $(southam69_objs)

modfischer60_objs = modfischer60.o main.o
modfischer60: $(modfischer60_objs)
	$(CC) $(LDFLAGS) -o modfischer60 $(modfischer60_objs)

modairy_objs = modairy.o main.o
modairy: $(modairy_objs)
	$(CC) $(LDFLAGS) -o modairy $(modairy_objs)

krassovsky40_objs = krassovsky40.o main.o
krassovsky40: $(krassovsky40_objs)
	$(CC) $(LDFLAGS) -o krassovsky40 $(krassovsky40_objs)

houg60_objs = houg60.o main.o
houg60: $(houg60_objs)
	$(CC) $(LDFLAGS) -o houg60 $(houg60_objs)

indo74_objs = indo74.o main.o
indo74: $(indo74_objs)
	$(CC) $(LDFLAGS) -o indo74 $(indo74_objs)

.m4.o:
	m4 $*.m4 utm.c.m4 >$*.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o

hayford24.m4: genutm
	genutm -a 6378388.000 -e 0.0067226700223332915 -k 0.9996 \
	  -c "European/Hayford/International, 1.924" >hayford24.m4

wgs84.m4: genutm
	genutm -a 6378137.000 -e 0.00669437999014 -k 0.9996 \
	  -c "World Geodetic, 1.984" >wgs84.m4

wgs72.m4: genutm
	genutm -a 6378135.0 -e 0.006694317778 -k 0.9996 \
	  -c "World Geodetic, 1.972" >wgs72.m4

southam69.m4: genutm
	genutm -a 6378160.0 -e 0.006694541854 -k 0.9996 \
	  -c "South American, 1.969" >southam69.m4

modfischer60.m4: genutm
	genutm -a 6378155.0 -e 0.006693421622 -k 0.9996 \
	  -c "Modified Fischer, 1.960" >modfischer60.m4

modairy.m4: genutm
	genutm -a 6377340.189 -e 0.006670540001 -k 0.9996 \
	  -c "Modified Airy" >modairy.m4

krassovsky40.m4: genutm
	genutm -a 6378245.0 -e 0.006693421622 -k 0.9996 \
	  -c "Krassovsky, 1.940" >krassovsky40.m4

houg60.m4: genutm
	genutm -a 6378270.0 -e 0.006722670022 -k 0.9996 \
	  -c "Hough, 1.960" >houg60.m4

indo74.m4: genutm
	genutm -a 6378160.0 -e 0.00669460908 -k 0.9996 \
	  -c "Indonesian, 1.974" >indo74.m4

utm.hp48: utm.hp48.m4 hayford24.hp48 wgs84.hp48
	m4 utm.hp48.m4 >utm.hp48

hayford24.hp48: hayford24.m4 utmtables.hp48.m4
	m4 hayford24.m4 utmtables.hp48.m4 >hayford24.hp48

wgs84.hp48: wgs84.m4 utmtables.hp48.m4
	m4 wgs84.m4 utmtables.hp48.m4 >wgs84.hp48

genutm_objs = genutm.o
genutm: $(genutm_objs)
	$(CC) $(LDFLAGS) -o genutm $(genutm_objs)
