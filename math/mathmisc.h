/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef AOS_LIB_MATHMISC_H
#define AOS_LIB_MATHMISC_H
#include "../sys/sys.h"
/**
   \file mathmisc.h
   A few math routines
*/
long double factorial(long n1, long n2);


#define sinc(x) (fabs(x)<1.e-5?1:sin(M_PI*x)/(M_PI*x))

void invsq(long n, double *restrict A);
#define mysqrt(A) (A<0?-sqrt(-A):sqrt(A))
long *invperm(long *p, long np);
void maxminlong(const long *restrict p, long N,
		long *restrict max, long *restrict min);
long nextpow2(long n);
long nextfftsize(long n);
unsigned long mylog2(unsigned long n);
typedef double(*golden_section_fun)(void *param, double x);
double golden_section_search(golden_section_fun f, void *param, 
			     double x1, double x4, double tau);
void readspintdata(file_t *fp, uint32_t magic, spint *out, long len);
spint *readspint(file_t *fp, long* nx, long* ny);
void readvec(void *p, uint32_t magic_p, uint32_t magic_file, size_t size, size_t nmemb, const file_t *fp);

#endif
