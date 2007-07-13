/* $Id: fft.h,v 2.1 2007/07/13 19:57:54 luis Exp $
 * Author: Luis Colorado <Luis.Colorado@HispaLinux.ES>
 * Date: Tue Jun  8 20:35:59 MEST 2004
 *
 * Disclaimer:
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Do not include anything BEFORE the line below, as it would not be
 * protected against double inclusion from other files
 */
#ifndef FFT_H
#define FFT_H

static char FFT_H_RCSId[] = "\n$Id: fft.h,v 2.1 2007/07/13 19:57:54 luis Exp $\n";

/* constants */

/* types */
typedef struct complex {
	double x;
	double y;
} complex_t;

typedef struct fft {
	int n; /* log2(N) */
	int N;
	complex_t *w;
	complex_t *W;
} fft_t;

/* prototypes */

void fft_init(fft_t *fft, int N);
void fft_direct(fft_t *fft, complex_t *a);
void fft_reverse(fft_t *fft, complex_t *a);

#endif /* FFT_H */
/* Do not include anything AFTER the line above, as it would not be
 * protected against double inclusion from other files.
 */

/* $Id: fft.h,v 2.1 2007/07/13 19:57:54 luis Exp $ */
