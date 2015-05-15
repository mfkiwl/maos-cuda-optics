/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_LIB_SP_H
#define AOS_LIB_SP_H
#ifndef AOS_LIB_MATH_H
#error "Don't include this file directly"
#endif
#include "random.h"

#define dspfree(A)      {dspfree_do(A); A=NULL;}
#define dspcellfree(A)  {dspcellfree_do(A); A=NULL;}
#define PDSPCELL(M,P)  dsp* (*restrict P)[(M)->nx]=(dsp*(*)[(M)->nx])(M)->p

#define sspfree(A)      {sspfree_do(A); A=NULL;}
#define sspcellfree(A)  {sspcellfree_do(A); A=NULL;}
#define PSSPCELL(M,P)   ssp* (*restrict P)[(M)->nx]=(ssp*(*)[(M)->nx])(M)->p

#define cspfree(A)     {cspfree_do(A); A=NULL;}
#define cspcellfree(A) {cspcellfree_do(A); A=NULL;}
#define PCSPCELL(M,P)  csp* (*restrict P)[(M)->nx]=(csp*(*)[(M)->nx])(M)->p

#define zspfree(A)     {zspfree_do(A); A=NULL;}
#define zspcellfree(A) {zspcellfree_do(A); A=NULL;}
#define PZSPCELL(M,P)  zsp* (*restrict P)[(M)->nx]=(zsp*(*)[(M)->nx])(M)->p

#define AOS_SP_DEF(X,T,R,RI) \
X(sp)* X(spnew)(long nx, long ny, long nzmax) CHECK_UNUSED_RESULT;\
X(sp) *X(spref)(X(sp) *A) CHECK_UNUSED_RESULT;\
X(sp) *X(spdup)(const X(sp) *A) CHECK_UNUSED_RESULT;\
X(sp) *X(spnew2)(const X(sp) *A) CHECK_UNUSED_RESULT;\
X(sp) *X(spnewrandu)(int nx, int ny, const T mean, R fill,rand_t *rstat) CHECK_UNUSED_RESULT;\
void X(spsetnzmax)(X(sp) *sp, long nzmax);\
void X(spfree_do)(X(sp) *sp);\
void X(spdisp)(const X(sp) *sp);\
int X(spcheck)(const X(sp) *sp);\
void X(spscale)(X(sp) *A, const T beta);\
void X(spcellscale)(X(spcell) *A, const T beta);\
X(sp)* X(spnewdiag)(long N, T *vec, T alpha) CHECK_UNUSED_RESULT;\
X(mat) *X(spdiag)(const X(sp) *A) CHECK_UNUSED_RESULT;\
void X(spmuldiag)(X(sp) *restrict A, const T* w, T alpha);\
void X(spmulvec)(T *restrict y, const X(sp) *A, const T * restrict x, char trans, T alpha); \
void X(spmulcreal)(T *restrict y, const X(sp) *A, const RI * restrict x, T alpha); \
void X(mulsp)(X(mat) **yout, const X(mat) *x, const X(sp) *A, const T alpha); \
void X(spmm)(X(mat) **yout, const X(sp) *A, const X(mat) *x, char trans, const T alpha); \
T X(spwdinn)(const X(mat) *y, const X(sp) *A, const X(mat) *x) CHECK_UNUSED_RESULT;\
T X(spcellwdinn)(const X(cell) *y, const X(spcell) *A, const X(cell) *x) CHECK_UNUSED_RESULT;\
void X(cellmulsp)(X(cell) **C0, const X(cell) *A, const X(spcell) *B, R alpha); \
void X(cellmm)(X(cell) **C0, const void *A_, const X(cell) *B, const char trans[2], const R alpha); \
X(sp) *X(2sp)(X(mat)*A);					\
void X(spfull)(X(mat) **out0, const X(sp) *A, const T f);\
void X(sptfull)(X(mat) **out0, const X(sp) *A, const T f);\
void X(spcellfull)(X(cell) **out0, const X(spcell) *A, const T f);\
void X(sptcellfull)(X(cell) **out0, const X(spcell) *A, const T f);\
X(sp) *X(spadd2)(X(sp) *A,X(sp)*B,T a,T b) CHECK_UNUSED_RESULT;\
void X(spadd)(X(sp) **A0, const X(sp) *B);\
void X(spcelladd)(X(spcell) **A0, const X(spcell) *B);\
void X(spaddI)(X(sp) **A0, R alpha);\
void X(spcelladdI)(X(spcell) *A0, R alpha);\
X(sp) *X(sptrans)(const X(sp) *A) CHECK_UNUSED_RESULT;\
X(sp) *X(spmulsp)(const X(sp) *A, const X(sp) *B) CHECK_UNUSED_RESULT;\
X(sp) *X(sptmulsp)(const X(sp) *A, const X(sp) *B) CHECK_UNUSED_RESULT;\
void X(spmulsp2)(X(sp) **C0, const X(sp) *A, const X(sp) *B, const T scale);\
X(spcell) *X(spcellmulspcell)(const X(spcell) *A, const X(spcell) *B, const T scale) CHECK_UNUSED_RESULT;\
X(spcell) *X(spcelltrans)(const X(spcell) *spc) CHECK_UNUSED_RESULT;\
void X(spcellfree_do)(X(spcell) *spc);\
X(sp) *X(spcat)(const X(sp) *A, const X(sp) *B, int type) CHECK_UNUSED_RESULT;\
X(sp) *X(spcell2sp)(const X(spcell) *A) CHECK_UNUSED_RESULT;\
X(mat) *X(spsum)(const X(sp) *A, int col) CHECK_UNUSED_RESULT;\
X(mat) *X(spsumabs)(const X(sp) *A, int col) CHECK_UNUSED_RESULT;\
void X(spclean)(X(sp) *A);\
void X(spdroptol)(X(sp) *A, R thres);		\
void X(spcelldroptol)(X(spcell) *A, R thres); \
void X(spsort)(X(sp) *A);\
void X(spcellsort)(X(spcell) *A);\
void X(spsym)(X(sp) **A);\
void X(spcellsym)(X(spcell) **A);\
X(sp) *X(spconvolvop)(X(mat) *A) CHECK_UNUSED_RESULT;\
X(sp) *X(spperm)(X(sp) *A, int reverse, long *pcol, long *prow) CHECK_UNUSED_RESULT;\
X(sp) *X(spinvbdiag)(const X(sp) *A, long bs) CHECK_UNUSED_RESULT;\
X(cell) *X(spblockextract)(const X(sp) *A, long bs) CHECK_UNUSED_RESULT;
#endif