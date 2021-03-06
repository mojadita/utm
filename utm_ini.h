/* $Id: utm_ini.h,v 2.4 2007/07/01 23:50:48 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@SLUG.CTV.ES>
 * Date: Fri Sep  6 02:07:35 MEST 2002
 * $Log: utm_ini.h,v $
 * Revision 2.4  2007/07/01 23:50:48  luis
 * * Eliminadas las trazas que se habían colocado en las funciones hms2h() y
 *   h2hms().  Ahora dependen de la definición de DEBUG.
 *
 * Revision 2.3  2007-07-01 22:31:20  luis
 * * Corregido un error debido al redondeo en hms2h() que hacía que se calculara
 *   mal las transformaciones de geodésicas a utm.
 *
 * Revision 2.2  2002-09-17 19:58:27  luis
 * Added more precision.
 *
 * Revision 2.1  2002/09/06 00:12:11  luis
 * A�adidos utm_ini.h para que genutm pueda calcular por tabla los par�metros y
 * utmcalc.c para los c�lculos a partir de los par�metros calculados seg�n la
 * estructura utmparam.
 *
 */

#include "utm.h"

struct utmparam wgs84table[] = {
	{ "WE", "World Geodetic, 1.984", 0.00669437999014, 6378137.000,  /* ... */ },
	{ "AA", "Airy, 1.830", 0.006670540001, 6377563.396,  /* ... */ },
	{ "AM", "Modified Airy", 0.006670540001, 6377340.189,  /* ... */ },
	{ "AN", "Australian National", 0.006694541854, 6378160.0,  /* ... */ },
	{ "BN", "Bessel, 1.841, Namibia", 0.006674372231, 6377483.865,  /* ... */ },
	{ "BR", "Bessel, 1.841, Ethiopia, Indonesia, Japan, and Korea", 0.006674372231, 6377397.155,  /* ... */ },
	{ "CC", "Clarke, 1.866", 0.006768657997, 6378206.4,  /* ... */ },
	{ "CD", "Clarke, 1.880", 0.006803511283, 6378249.145,  /* ... */ },
	{ "EA", "Everest, India 1830", 0.006637846631, 6377276.345,  /* ... */ },
	{ "EB", "Everest, Brunei and E. Malaysia (Sabah and Sarawak)", 0.006637846631, 6377298.556,  /* ... */ },
	{ "EC", "Everest, India 1956", 0.006637846631, 6377301.243,  /* ... */ },
	{ "ED", "Everest, W. Malaysia 1969", 0.006637846631, 6377295.664,  /* ... */ },
	{ "EE", "Everest, W. Malaysia and Singapore 1948", 0.006637846631, 6377304.063,  /* ... */ },
	{ "EF", "Everest, Pakistan", 0.006637846631, 6377309.613,  /* ... */ },
	{ "FA", "Modified Fischer, 1.960", 0.006693421622, 6378155.0,  /* ... */ },
	{ "HE", "Helmert 1906", 0.006693421622, 6378200.0,  /* ... */ },
	{ "HO", "Hough, 1.960", 0.006722670022, 6378270.0,  /* ... */ },
	{ "ID", "Indonesian, 1.974", 0.00669460908, 6378160.0,  /* ... */ },
	{ "IN", "European/Hayford/International, 1.924", 0.0067226700223332915, 6378388.000,  /* ... */ },
	{ "KA", "Krassovsky, 1.940", 0.006693421622, 6378245.0,  /* ... */ },
	{ "RF", "Geodetic Reference System 1980", 0.006694380023, 6378137,  /* ... */ },
	{ "SA", "South American, 1.969", 0.006694541854, 6378160.0,  /* ... */ },
	{ "ST", "Struve, Espa�a, anterior a 1970", 0.00677436, 6378298.3, /* ... */ },
	{ "WD", "World Geodetic, 1.972", 0.006694317778, 6378135.0,  /* ... */ },
	{ NULL /* ... */ },
}; /* wgs84table [] */

#if 0
struct geosys geo_sys[] {
	{ "ADI-M", "Adindan/Mean Solution (Ethiopia and Sudan)", &geo_CD, -166.0, -15.0, +204.0 },
	{ "ADI-E", "Adindan/Burkina Faso", &geo_CD, -118.0, -14.0, +218.0 },
	{ "ADI-F", "Adindan/Cameroon", &geo_CD, -134.0, -2.0, +210.0 },
	{ "ADI-A", "Adindan/Ethiopia", &geo_CD, -165.0, -11.0, +206.0 },
	{ "ADI-C", "Adindan/Mali", &geo_CD, -123.0, -20.0, +220.0 },
	{ "ADI-D", "Adindan/Senegal", &geo_CD, -128.0, -18.0, +224.0 },
	{ "ADI-B", "Adindan/Sudan", &geo_CD, -161.0 -14.0, +205.0 },
	{ "AFG", "Afgooye/Somalia", &geo_KA, -43.0, -163.0, +45.0 },
	{ "ARF-A", "Arc 1950/Botswana", &geo_CD, -138.0, -105.0, -289.0 },
	{ "ARF-B", "Arc 1950/Lesotho", &geo_CD, -125.0, -108.0, -295.0 },
	{ "ARF-C", "Arc 1950/Malawi", &geo_CD, -161.0, -73.0, -317.0 },
	{ "ARF-D", "Arc 1950/Swaziland", &geo_CD, -134.0, -105.0, -295.0 },
	{ "ARF-E", "Arc 1950/Zaire", &geo_CD, -169.0, -19.0, -278.0 },
	{ "ARF-F", "Arc 1950/Zambia", &geo_CD, -147.0, -74.0, -283.0 },
	{ "ARF-G", "Arc 1950/Zimbabwe", &geo_CD, -142.0, -96.0, -293.0 },
	{ "ARF-H", "Arc 1950/Burundi", &geo_CD, -153.0, -5.0, -292.0 },
	{ "ARF-M", "Arc 1950/Mean", &geo_CD, -143.0, -90.0, -294.0 },
	{ "", "Ain el Abd 1970", &geo_IN, },
	{ "", "American Samoa 1962", &geo_CC, },
	{ "", "Anna 1 Astro 1965", &geo_AN, },
	{ "", "Antigua Island Astro 1943", &geo_CD, },
	{ "", "Arc 1960", &geo_CD, },
	{ "", "Ascension Island 1958", &geo_IN, },
	{ "", "Astro Beacon E 1945", &geo_IN, },
	{ "", "Astro DOS 71/4", &geo_IN, },
	{ "", "Astro Tern Island (FRIG) 1961", &geo_IN, },
	{ "", "Astronomical Station 1952", &geo_IN, },
	{ "", "Australian Geodetic 1966", &geo_AN, },
	{ "", "Australian Geodetic 1984", &geo_AN, },
	{ "", "Ayabelle Lighthouse", &geo_CD, },
	{ "", "Bellevue (IGN)", &geo_IN, },
	{ "", "Bermuda 1957", &geo_CC, },
	{ "", "Bissau", &geo_IN, },
	{ "", "Bogota Observatory", &geo_IN, },
	{ "", "Campo Inchauspe", &geo_IN, },
	{ "", "Canton Astro 1966", &geo_IN, },
	{ "", "Cape", &geo_CD, },
	{ "", "Cape Canaveral", &geo_CC, },
	{ "", "Carthage", &geo_CD, },
	{ "", "Chatham Island Astro 1971", &geo_IN, },
	{ "", "Chua Astro", &geo_IN, },
	{ "", "Co-Ordinate System 1937 of Estonia", &geo_BR, },
	{ "", "Corrego Alegre", &geo_IN, },
	{ "", "Dabola", &geo_CD, },
	{ "", "Deception Island", &geo_CD, },
	{ "", "Djakarta (Batavia)", &geo_BR, },
	{ "", "DOS 1968", &geo_IN, },
	{ "", "Easter Island 1967", &geo_IN, },
	{ "", "European 1950", &geo_IN, },
	{ "", "European 1979", &geo_IN, },
	{ "", "Fort Thomas 1955", &geo_CD, },
	{ "", "Gan 1970", &geo_IN, },
	{ "", "Geodetic Datum 1949", &geo_IN, },
	{ "", "Graciosa Base SW 1948", &geo_IN, },
	{ "", "Guam 1963", &geo_CC, },
	{ "", "GUX 1 Astro", &geo_IN, },
	{ "", "Hjorsey 1955", &geo_IN, },
	{ "", "Hong Kong 1963", &geo_IN, },
	{ "", "Hu-Tzu-Shan", &geo_IN, },
	{ "", "Indian", &geo_EA, },
	{ "", "Indian 1954", &geo_EA, },
	{ "", "Indian 1960", &geo_EA, },
	{ "", "Indian 1975", &geo_EA, },
	{ "", "Indonesian 1974", &geo_ID, },
	{ "", "Ireland 1965", &geo_AM, },
	{ "", "ISTS 061 Astro 1968", &geo_IN, },
	{ "", "ISTS 073 Astro 1969", &geo_IN, },
	{ "", "Johnston Island 1961", &geo_IN, },
	{ "", "Kandawala", &geo_EA, },
	{ "", "Kerguelen Island 1949", &geo_IN, },
	{ "", "Kertau 1948", &geo_EE, },
	{ "", "Korean Geodetic System 1995", &geo_WE, },
	{ "", "Kusaie Astro 1951", &geo_IN, },
	{ "", "L. C. 5 Astro 1961", &geo_CC, },
	{ "", "Leigon", &geo_CD, },
	{ "", "Liberia 1964", &geo_CD, },
	{ "", "Luzon", &geo_CC, },
	{ "", "Mahe 1971", &geo_CD, },
	{ "", "Massawa", &geo_BR, },
	{ "", "Merchich", &geo_CD, },
	{ "", "Midway Astro 1961", &geo_IN, },
	{ "", "Minna", &geo_CD, },
	{ "", "Montserrat Island Astro 1958", &geo_CD, },
	{ "", "M'Poraloko", &geo_CD, },
	{ "", "Nahrwan", &geo_CD, },
	{ "", "Naparima, BWI", &geo_IN, },
	{ "", "North American 1927", &geo_CC, },
	{ "", "North American 1983", &geo_RF, },
	{ "", "North Sahara 1959", &geo_CD, },
	{ "", "Observatorio Meteorologico 1939", &geo_IN, },
	{ "", "Old Egyptian 1907", &geo_HE, },
	{ "", "Old Hawaiian", &geo_CC, },
	{ "", "Old Hawaiian", &geo_IN, },
	{ "", "Oman", &geo_CD, },
	{ "", "Ordnance Survey of Great Britain 1936", &geo_AA, },
	{ "", "Pico de las Nieves", &geo_IN, },
	{ "", "Pitcairn Astro 1967", &geo_IN, },
	{ "", "Point 58", &geo_CD, },
	{ "", "Pointe Noire 1948", &geo_CD, },
	{ "", "Porto Santo 1936", &geo_IN, },
	{ "", "Provisional South American 1956", &geo_IN, },
	{ "", "Provisional South Chilean 1963", &geo_IN, },
	{ "", "Puerto Rico", &geo_CC, },
	{ "", "Qatar National", &geo_IN, },
	{ "", "Qornoq", &geo_IN, },
	{ "", "Reunion", &geo_IN, },
	{ "", "Rome 1940", &geo_IN, },
	{ "", "S-42 (Pulkovo 1942)", &geo_KA, },
	{ "", "Santo (DOS) 1965", &geo_IN, },
	{ "", "Sao Braz", &geo_IN, },
	{ "", "Sapper Hill 1943", &geo_IN, },
	{ "", "Schwarzeck", &geo_BN, },
	{ "", "Selvagem Grande 1938", &geo_IN, },
	{ "", "Sierra Leone 1960", &geo_CD, },
	{ "", "S-JTSK", &geo_BR, },
	{ "", "South American 1969", &geo_SA, },
	{ "", "South American Geocentric Reference System (SIRGAS)", &geo_RF, },
	{ "", "South Asia", &geo_FA, },
	{ "", "Timbalai 1948", &geo_EB, },
	{ "", "Tokyo", &geo_BR, },
	{ "", "Tristan Astro 1968", &geo_IN, },
	{ "", "Viti Levu 1916", &geo_CD, },
	{ "", "Voirol 1960", &geo_CD, },
	{ "", "Wake-Eniwetok 1960", &geo_HO, },
	{ "", "Wake Island Astro 1952", &geo_IN, },
	{ "", "Zanderij", &geo_IN, },
}; /* geo_sys */
#endif

/* $Id: utm_ini.h,v 2.4 2007/07/01 23:50:48 luis Exp $ */
