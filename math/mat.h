/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_LIB_MAT_H
#define AOS_LIB_MAT_H
#ifndef AOS_LIB_MATH_H
#error "Don't include this file directly"
#endif

#define AOS_MAT_DEF(X,T,R)						\
    X(mat) *X(new_ref)(long nx, long ny, T *p) CHECK_UNUSED_RESULT;	\
    X(mat) *X(new_data)(long nx, long ny, T *p) CHECK_UNUSED_RESULT;	\
    X(mat) *X(new)(long nx, long ny) CHECK_UNUSED_RESULT;		\
    X(mat) *X(mat_cast)(const TwoDim *A) CHECK_UNUSED_RESULT;		\
    void X(init)(X(mat)**A, long nx, long ny) ;				\
    void X(free_do)(X(mat) *A, int keepdata);				\
    X(mat) *X(ref)(const X(mat) *in) CHECK_UNUSED_RESULT;		\
    X(mat) *X(ref_reshape)(const X(mat) *in, long nx, long ny) CHECK_UNUSED_RESULT; \
    X(mat) *X(refcols)(const X(mat) *in, long icol, long ncol) CHECK_UNUSED_RESULT; \
    void X(resize)(X(mat) *A, long nx, long ny);			\
    X(mat) *X(sub)(const X(mat) *in, long sx, long nx, long sy, long ny) CHECK_UNUSED_RESULT; \
    X(mat) *X(cat)(const X(mat) *in1, const X(mat) *in2, int dim) CHECK_UNUSED_RESULT; \
    X(mat) *X(dup)(const X(mat) *in) CHECK_UNUSED_RESULT;		\
    void X(zero)(X(mat)*A);						\
    void X(zerocol)(X(mat)*A, int icol);				\
    uint32_t X(hash)(const X(mat)*A, uint32_t key);			\
    void X(cp)(X(mat) **out0, const X(mat) *in);			\
    X(mat) *X(trans)(const X(mat) *A) CHECK_UNUSED_RESULT;		\
    void X(set)(X(mat) *A, const T val);				\
    void X(show)(const X(mat) *A, const char *format,...) CHECK_ARG(2);	\
    T X(sum)(const X(mat) *A) CHECK_UNUSED_RESULT;			\
    T X(trace)(const X(mat) *A) CHECK_UNUSED_RESULT;			\
    void X(sort)(X(mat) *A, int ascend);				\
    void X(maxmin)(const T *restrict p, long N, R *max, R *min);	\
    R X(max)(const X(mat) *A) CHECK_UNUSED_RESULT;			\
    R X(maxabs)(const X(mat) *A) CHECK_UNUSED_RESULT;			\
    R X(min)(const X(mat) *A) CHECK_UNUSED_RESULT;			\
    R X(sumabs)(const X(mat) *in) CHECK_UNUSED_RESULT;			\
    R X(sumsq)(const X(mat) *in) CHECK_UNUSED_RESULT;			\
    void X(fftshift)(X(mat) *A);					\
    void X(cpcorner2center)(X(mat) *A, const X(mat)*B);			\
    X(cell) *X(cell_cast)(const TwoDim *A) CHECK_UNUSED_RESULT;		\
    X(cell) *X(cellnew2)(const X(cell) *A);				\
    X(cell) *X(cellnew3)(long nx, long ny, long *nnx, long *nny);	\
    X(cell) *X(cellnew_same)(long nx, long ny, long mx, long my);	\
    void X(cellzero)(TwoDim *dc);						\
    void X(cellset)(X(cell)*dc, T alpha);				\
    X(cell) *X(celltrans)(const X(cell) *A);				\
    X(cell) *X(cellref)(const X(cell) *in);				\
    X(cell) *X(celldup)(const X(cell) *in);				\
    void X(cellcp)(void* out0, const TwoDim*in);			\
    X(cell) *X(cellreduce)(const X(cell) *A, int dim);			\
    X(cell) *X(cellcat)(const X(cell) *A, const X(cell) *B, int dim);	\
    X(cell) *X(cellcat_each)(const X(cell) *A, const X(cell) *B, int dim); \
    X(mat) *X(cell2m)(const TwoDim *A);					\
    X(cell)* X(2cellref)(const X(mat) *A, long*dims, long ndim);	\
    void X(2cell)(X(cell) **B, const X(mat) *A, const X(cell) *ref);	\

#endif
