/* $Id: mkroots.h,v 2.2 2007/07/13 21:55:12 luis Exp $
 * Author: Luis.Colorado@HispaLinux.ES
 */

static complex_t e_0[1] = {{ 1.0, 0.0 }};

static complex_t e_1[1] = {
	{ -1.00000000000000000, 0.00000000000000012}
}; /* static complex_t e_1[1] */

static complex_t e_2[2] = {
	{ 0.00000000000000006, 1.00000000000000000},
	{ -0.00000000000000018, -1.00000000000000000}
}; /* static complex_t e_2[2] */

static complex_t e_3[4] = {
	{ 0.70710678118654757, 0.70710678118654746},
	{ -0.70710678118654746, 0.70710678118654757},
	{ -0.70710678118654768, -0.70710678118654746},
	{ 0.70710678118654735, -0.70710678118654768}
}; /* static complex_t e_3[4] */

static complex_t e_4[8] = {
	{ 0.92387953251128674, 0.38268343236508978},
	{ 0.38268343236508984, 0.92387953251128674},
	{ -0.38268343236508973, 0.92387953251128674},
	{ -0.92387953251128674, 0.38268343236508989},
	{ -0.92387953251128685, -0.38268343236508967},
	{ -0.38268343236509034, -0.92387953251128652},
	{ 0.38268343236509000, -0.92387953251128663},
	{ 0.92387953251128652, -0.38268343236509039}
}; /* static complex_t e_4[8] */

static complex_t e_5[16] = {
	{ 0.98078528040323043, 0.19509032201612825},
	{ 0.83146961230254524, 0.55557023301960218},
	{ 0.55557023301960229, 0.83146961230254524},
	{ 0.19509032201612833, 0.98078528040323043},
	{ -0.19509032201612819, 0.98078528040323043},
	{ -0.55557023301960196, 0.83146961230254546},
	{ -0.83146961230254535, 0.55557023301960218},
	{ -0.98078528040323043, 0.19509032201612861},
	{ -0.98078528040323043, -0.19509032201612836},
	{ -0.83146961230254546, -0.55557023301960196},
	{ -0.55557023301960218, -0.83146961230254524},
	{ -0.19509032201612866, -0.98078528040323032},
	{ 0.19509032201612830, -0.98078528040323043},
	{ 0.55557023301960184, -0.83146961230254546},
	{ 0.83146961230254524, -0.55557023301960218},
	{ 0.98078528040323032, -0.19509032201612872}
}; /* static complex_t e_5[16] */

static complex_t e_6[32] = {
	{ 0.99518472667219693, 0.09801714032956060},
	{ 0.95694033573220882, 0.29028467725446233},
	{ 0.88192126434835505, 0.47139673682599764},
	{ 0.77301045336273699, 0.63439328416364549},
	{ 0.63439328416364549, 0.77301045336273699},
	{ 0.47139673682599781, 0.88192126434835494},
	{ 0.29028467725446233, 0.95694033573220894},
	{ 0.09801714032956077, 0.99518472667219682},
	{ -0.09801714032956065, 0.99518472667219693},
	{ -0.29028467725446216, 0.95694033573220894},
	{ -0.47139673682599770, 0.88192126434835505},
	{ -0.63439328416364538, 0.77301045336273710},
	{ -0.77301045336273699, 0.63439328416364549},
	{ -0.88192126434835494, 0.47139673682599786},
	{ -0.95694033573220882, 0.29028467725446239},
	{ -0.99518472667219682, 0.09801714032956083},
	{ -0.99518472667219693, -0.09801714032956059},
	{ -0.95694033573220894, -0.29028467725446211},
	{ -0.88192126434835505, -0.47139673682599764},
	{ -0.77301045336273710, -0.63439328416364527},
	{ -0.63439328416364593, -0.77301045336273666},
	{ -0.47139673682599786, -0.88192126434835494},
	{ -0.29028467725446244, -0.95694033573220882},
	{ -0.09801714032956045, -0.99518472667219693},
	{ 0.09801714032956009, -0.99518472667219693},
	{ 0.29028467725446205, -0.95694033573220894},
	{ 0.47139673682599759, -0.88192126434835505},
	{ 0.63439328416364560, -0.77301045336273688},
	{ 0.77301045336273666, -0.63439328416364593},
	{ 0.88192126434835483, -0.47139673682599792},
	{ 0.95694033573220882, -0.29028467725446250},
	{ 0.99518472667219693, -0.09801714032956051}
}; /* static complex_t e_6[32] */

static complex_t e_7[64] = {
	{ 0.99879545620517241, 0.04906767432741801},
	{ 0.98917650996478101, 0.14673047445536175},
	{ 0.97003125319454397, 0.24298017990326387},
	{ 0.94154406518302081, 0.33688985339222005},
	{ 0.90398929312344334, 0.42755509343028208},
	{ 0.85772861000027212, 0.51410274419322166},
	{ 0.80320753148064494, 0.59569930449243336},
	{ 0.74095112535495911, 0.67155895484701833},
	{ 0.67155895484701833, 0.74095112535495911},
	{ 0.59569930449243347, 0.80320753148064483},
	{ 0.51410274419322166, 0.85772861000027212},
	{ 0.42755509343028220, 0.90398929312344334},
	{ 0.33688985339222005, 0.94154406518302081},
	{ 0.24298017990326398, 0.97003125319454397},
	{ 0.14673047445536175, 0.98917650996478101},
	{ 0.04906767432741813, 0.99879545620517241},
	{ -0.04906767432741801, 0.99879545620517241},
	{ -0.14673047445536164, 0.98917650996478101},
	{ -0.24298017990326387, 0.97003125319454397},
	{ -0.33688985339221994, 0.94154406518302081},
	{ -0.42755509343028186, 0.90398929312344345},
	{ -0.51410274419322166, 0.85772861000027212},
	{ -0.59569930449243336, 0.80320753148064494},
	{ -0.67155895484701844, 0.74095112535495899},
	{ -0.74095112535495888, 0.67155895484701855},
	{ -0.80320753148064483, 0.59569930449243347},
	{ -0.85772861000027201, 0.51410274419322177},
	{ -0.90398929312344334, 0.42755509343028203},
	{ -0.94154406518302070, 0.33688985339222033},
	{ -0.97003125319454397, 0.24298017990326407},
	{ -0.98917650996478101, 0.14673047445536180},
	{ -0.99879545620517241, 0.04906767432741797},
	{ -0.99879545620517241, -0.04906767432741772},
	{ -0.98917650996478101, -0.14673047445536158},
	{ -0.97003125319454397, -0.24298017990326382},
	{ -0.94154406518302081, -0.33688985339222011},
	{ -0.90398929312344345, -0.42755509343028181},
	{ -0.85772861000027212, -0.51410274419322155},
	{ -0.80320753148064494, -0.59569930449243325},
	{ -0.74095112535495911, -0.67155895484701844},
	{ -0.67155895484701866, -0.74095112535495888},
	{ -0.59569930449243313, -0.80320753148064505},
	{ -0.51410274419322177, -0.85772861000027201},
	{ -0.42755509343028247, -0.90398929312344312},
	{ -0.33688985339221994, -0.94154406518302081},
	{ -0.24298017990326412, -0.97003125319454397},
	{ -0.14673047445536230, -0.98917650996478090},
	{ -0.04906767432741803, -0.99879545620517241},
	{ 0.04906767432741766, -0.99879545620517241},
	{ 0.14673047445536194, -0.98917650996478090},
	{ 0.24298017990326376, -0.97003125319454397},
	{ 0.33688985339221961, -0.94154406518302092},
	{ 0.42755509343028214, -0.90398929312344334},
	{ 0.51410274419322155, -0.85772861000027223},
	{ 0.59569930449243291, -0.80320753148064528},
	{ 0.67155895484701833, -0.74095112535495911},
	{ 0.74095112535495888, -0.67155895484701866},
	{ 0.80320753148064505, -0.59569930449243325},
	{ 0.85772861000027201, -0.51410274419322188},
	{ 0.90398929312344312, -0.42755509343028253},
	{ 0.94154406518302081, -0.33688985339222000},
	{ 0.97003125319454397, -0.24298017990326418},
	{ 0.98917650996478090, -0.14673047445536239},
	{ 0.99879545620517241, -0.04906767432741809}
}; /* static complex_t e_7[64] */

static complex_t *W_0[1] = {
	&e_0[0]
}; /* static complex_t *W_0[1] */

static complex_t *W_1[2] = {
	&e_0[0],
	&e_1[0]
}; /* static complex_t *W_1[2] */

static complex_t *W_2[4] = {
	&e_0[0],
	&e_2[0],
	&e_1[0],
	&e_2[1]
}; /* static complex_t *W_2[4] */

static complex_t *W_3[8] = {
	&e_0[0],
	&e_3[0],
	&e_2[0],
	&e_3[1],
	&e_1[0],
	&e_3[2],
	&e_2[1],
	&e_3[3]
}; /* static complex_t *W_3[8] */

static complex_t *W_4[16] = {
	&e_0[0],
	&e_4[0],
	&e_3[0],
	&e_4[1],
	&e_2[0],
	&e_4[2],
	&e_3[1],
	&e_4[3],
	&e_1[0],
	&e_4[4],
	&e_3[2],
	&e_4[5],
	&e_2[1],
	&e_4[6],
	&e_3[3],
	&e_4[7]
}; /* static complex_t *W_4[16] */

static complex_t *W_5[32] = {
	&e_0[0],
	&e_5[0],
	&e_4[0],
	&e_5[1],
	&e_3[0],
	&e_5[2],
	&e_4[1],
	&e_5[3],
	&e_2[0],
	&e_5[4],
	&e_4[2],
	&e_5[5],
	&e_3[1],
	&e_5[6],
	&e_4[3],
	&e_5[7],
	&e_1[0],
	&e_5[8],
	&e_4[4],
	&e_5[9],
	&e_3[2],
	&e_5[10],
	&e_4[5],
	&e_5[11],
	&e_2[1],
	&e_5[12],
	&e_4[6],
	&e_5[13],
	&e_3[3],
	&e_5[14],
	&e_4[7],
	&e_5[15]
}; /* static complex_t *W_5[32] */

static complex_t *W_6[64] = {
	&e_0[0],
	&e_6[0],
	&e_5[0],
	&e_6[1],
	&e_4[0],
	&e_6[2],
	&e_5[1],
	&e_6[3],
	&e_3[0],
	&e_6[4],
	&e_5[2],
	&e_6[5],
	&e_4[1],
	&e_6[6],
	&e_5[3],
	&e_6[7],
	&e_2[0],
	&e_6[8],
	&e_5[4],
	&e_6[9],
	&e_4[2],
	&e_6[10],
	&e_5[5],
	&e_6[11],
	&e_3[1],
	&e_6[12],
	&e_5[6],
	&e_6[13],
	&e_4[3],
	&e_6[14],
	&e_5[7],
	&e_6[15],
	&e_1[0],
	&e_6[16],
	&e_5[8],
	&e_6[17],
	&e_4[4],
	&e_6[18],
	&e_5[9],
	&e_6[19],
	&e_3[2],
	&e_6[20],
	&e_5[10],
	&e_6[21],
	&e_4[5],
	&e_6[22],
	&e_5[11],
	&e_6[23],
	&e_2[1],
	&e_6[24],
	&e_5[12],
	&e_6[25],
	&e_4[6],
	&e_6[26],
	&e_5[13],
	&e_6[27],
	&e_3[3],
	&e_6[28],
	&e_5[14],
	&e_6[29],
	&e_4[7],
	&e_6[30],
	&e_5[15],
	&e_6[31]
}; /* static complex_t *W_6[64] */

static complex_t *W_7[128] = {
	&e_0[0],
	&e_7[0],
	&e_6[0],
	&e_7[1],
	&e_5[0],
	&e_7[2],
	&e_6[1],
	&e_7[3],
	&e_4[0],
	&e_7[4],
	&e_6[2],
	&e_7[5],
	&e_5[1],
	&e_7[6],
	&e_6[3],
	&e_7[7],
	&e_3[0],
	&e_7[8],
	&e_6[4],
	&e_7[9],
	&e_5[2],
	&e_7[10],
	&e_6[5],
	&e_7[11],
	&e_4[1],
	&e_7[12],
	&e_6[6],
	&e_7[13],
	&e_5[3],
	&e_7[14],
	&e_6[7],
	&e_7[15],
	&e_2[0],
	&e_7[16],
	&e_6[8],
	&e_7[17],
	&e_5[4],
	&e_7[18],
	&e_6[9],
	&e_7[19],
	&e_4[2],
	&e_7[20],
	&e_6[10],
	&e_7[21],
	&e_5[5],
	&e_7[22],
	&e_6[11],
	&e_7[23],
	&e_3[1],
	&e_7[24],
	&e_6[12],
	&e_7[25],
	&e_5[6],
	&e_7[26],
	&e_6[13],
	&e_7[27],
	&e_4[3],
	&e_7[28],
	&e_6[14],
	&e_7[29],
	&e_5[7],
	&e_7[30],
	&e_6[15],
	&e_7[31],
	&e_1[0],
	&e_7[32],
	&e_6[16],
	&e_7[33],
	&e_5[8],
	&e_7[34],
	&e_6[17],
	&e_7[35],
	&e_4[4],
	&e_7[36],
	&e_6[18],
	&e_7[37],
	&e_5[9],
	&e_7[38],
	&e_6[19],
	&e_7[39],
	&e_3[2],
	&e_7[40],
	&e_6[20],
	&e_7[41],
	&e_5[10],
	&e_7[42],
	&e_6[21],
	&e_7[43],
	&e_4[5],
	&e_7[44],
	&e_6[22],
	&e_7[45],
	&e_5[11],
	&e_7[46],
	&e_6[23],
	&e_7[47],
	&e_2[1],
	&e_7[48],
	&e_6[24],
	&e_7[49],
	&e_5[12],
	&e_7[50],
	&e_6[25],
	&e_7[51],
	&e_4[6],
	&e_7[52],
	&e_6[26],
	&e_7[53],
	&e_5[13],
	&e_7[54],
	&e_6[27],
	&e_7[55],
	&e_3[3],
	&e_7[56],
	&e_6[28],
	&e_7[57],
	&e_5[14],
	&e_7[58],
	&e_6[29],
	&e_7[59],
	&e_4[7],
	&e_7[60],
	&e_6[30],
	&e_7[61],
	&e_5[15],
	&e_7[62],
	&e_6[31],
	&e_7[63]
}; /* static complex_t *W_7[128] */

complex_t **FFT_W [32] = {
	W_0,
	W_1,
	W_2,
	W_3,
	W_4,
	W_5,
	W_6,
	W_7
}; /* complex_t **FFT_ROOTS [32] */


/* $Id: mkroots.h,v 2.2 2007/07/13 21:55:12 luis Exp $ */