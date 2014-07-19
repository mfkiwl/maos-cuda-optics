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

#define MAT_VERBOSE 0
/*#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>*/
#include "../sys/sys.h"
#include "random.h"
#include "mathdef.h"
#include "defs.h"/*Defines T, X, etc */

/**
   Work horse function that creats the matrix object. if p is
   NULL, memory is allocated. Allocation for X(mat) objects
   should always be created directly or indirectly through
   X(new_do) so that changing data structure is straight
   forward.
*/
static inline X(mat) *X(new_do)(long nx, long ny, T *p, int ref){
    X(mat) *out=calloc(1, sizeof(X(mat)));
    out->id=M_T;
    out->nx=nx;
    out->ny=ny;
    if(ref){/*the data does not belong to us. */
	if(!p){
	    error("When ref is 1, p must not be NULL.\n");
	}
	out->p=p;
    }else{
	if(!p){
	    p=calloc((nx*ny), sizeof(T));
	}
	out->p=p;
	out->nref=calloc(1, sizeof(int));
	out->nref[0]=1;
    }
    return out;
}
/**
   Creat a X(mat) object to reference an already existing
   vector.  free the X(mat) object won't free the existing
   vector.
*/
X(mat) *X(new_ref)(long nx, long ny, T *p){
    return X(new_do)(nx,ny,p,1);
}

/**
   Creat a X(mat) object with already allocated memory
   chunk. the memory is freed when the X(mat) is freed.
*/
X(mat) *X(new_data)(long nx, long ny, T *p){
    return X(new_do)(nx,ny,p,0);
}

/**
   Create a new T matrix object. initialized all to zero.
*/
X(mat) *X(new)(long nx, long ny){
    return X(new_do)(nx,ny,NULL,0);
}
/**
   check the size of matrix if exist. Otherwise create it
*/
void X(init)(X(mat)**A, long nx, long ny){
    if(*A){
	assert((*A)->nx == nx && (*A)->ny==ny);
    }else{
	*A=X(new)(nx, ny);
    }
}
/**
   Free the X(mat), but keep the data.
*/
void X(free_keepdata)(X(mat) *A){
    X(free_do)(A,1);
}
/**
   free a X(mat) object. if keepdata!=0, will not free A->p.
*/
void X(free_do)(X(mat) *A, int keepdata){
    if(!A) return;
    int free_extra=0;
    if(A->nref){
	int nref=atomicadd(A->nref, -1);
	if(!nref){
	    if(!keepdata && A->p){
		if(A->mmap){/*data is mmap'ed. */
		    mmap_unref(A->mmap);
		}else{
		    free_extra=1;
		    free(A->p);
		}
	    }
	    free(A->nref);
	}else if(nref<0){
	    error("The ref is less than 0. something wrong!!!:%d\n",A->nref[0]);
	}
    }else{
	free_extra=1;
    }
    if(free_extra){
#ifndef USE_LONG
	X(fft_free_plan)(A->fft);
#endif
	free(A->header);
    }
    free(A);
}

/**
   creat a X(mat) reference an existing X(mat). Use the reference carefully.
*/
X(mat) *X(ref)(const X(mat) *in){
    if(!in) return NULL;
    X(mat) *out=calloc(1, sizeof(X(mat)));
    memcpy(out,in,sizeof(X(mat)));
    if(!in->nref){
	extern quitfun_t quitfun;
	if(quitfun==&default_quitfun){
	    warning_once("Referencing non-referrable data. This may cause error.\n");
	}
    }else{
	atomicadd(in->nref, 1);
    }
    return out;
}
/**
   create an new X(mat) reference another with different shape.
*/
X(mat) *X(ref_reshape)(X(mat) *in, long nx, long ny){
    X(mat) *out=X(ref)(in);
    if(in->nx*in->ny!=nx*ny){
	error("Must not change number of elements\n");
    }
    out->nx=nx;
    out->ny=ny;
    return out;
}

/**
   creat a new X(mat) referencing columns in existing
   X(mat). reference counted. not used
*/
X(mat)* X(refcols)(X(mat) *in, long icol, long ncol){
    X(mat) *out=calloc(1, sizeof(X(mat)));
    out->nx=in->nx;
    out->ny=ncol;
    out->p=in->p+icol*out->nx;
    return out;
}

/**
   Create a new sub matrix of nx*ny starting from(sx,sy)
*/
X(mat) *X(sub)(const X(mat) *in, long sx, long nx, long sy, long ny){
    if(nx<=0){
	nx=in->nx-sx;
    }
    if(ny<=0){
	ny=in->ny-sy;
    }
    X(mat)*out=X(new)(nx, ny);
    if(sx+nx>in->nx || sy+ny>in->ny){
	error("Invalid parameter range\n");
    }
    PMAT(in, pin);
    PMAT(out, pout);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    pout[iy][ix]=pin[iy+sy][ix+sx];
	}
    }
    return out;
}
/**
   Resize a matrix by adding or removing columns or rows. Data is kept
   whever possible.
*/
void X(resize)(X(mat) *A, long nx, long ny){
    if(!A->nref || A->nref[0]>1){
	error("Resizing a referenced vector\n");
    }
    if(A->nx==nx || A->ny==1){
	A->p=realloc(A->p, sizeof(T)*nx*ny);
	if(nx*ny>A->nx*A->ny){
	    memset(A->p+A->nx*A->ny, 0, (nx*ny-A->nx*A->ny)*sizeof(T));
	}
    }else{/*copy memory to preserve data*/
	T *p=calloc(nx*ny,sizeof(T));
	long minx=A->nx<nx?A->nx:nx;
	long miny=A->ny<ny?A->ny:ny;
	for(long iy=0; iy<miny; iy++){
	    memcpy(p+iy*nx, A->p+iy*A->nx, sizeof(T)*minx);
	}
	free(A->p);
	A->p=p;
    }
    A->nx=nx;
    A->ny=ny;
}

/**
   concatenate two matrixes into 1 along dimension "dim"
*/
X(mat)* X(cat)(const X(mat) *in1, const X(mat) *in2, int dim){
    if(!in2){
	if(in1){
	    return X(dup)(in1);
	}else{
	    return NULL;
	}
    }else if(!in1){
	return X(dup)(in2);
    }

    X(mat) *out=NULL;

    if(dim==1){
	/*along x. */
	if(in1->ny!=in2->ny){
	    error("Mismatch. in1 is (%ld, %ld), in2 is (%ld, %ld)\n", 
		  in1->nx, in1->ny, in2->nx, in2->ny);
	}
	out=X(new)(in1->nx+in2->nx, in1->ny);
	PMAT(out,pout);
	PMAT(in1,pin1);
	PMAT(in2,pin2);
	for(long iy=0; iy<in1->ny; iy++){
	    memcpy(pout[iy],pin1[iy], in1->nx*sizeof(T));
	    memcpy(pout[iy]+in1->nx, pin2[iy], in2->nx*sizeof(T));
	}
    }else if(dim==2){
	/*along y. */
	if(in1->nx!=in2->nx){
	    error("Mismatch. in1 is (%ld, %ld), in2 is (%ld, %ld)\n", 
		  in1->nx, in1->ny, in2->nx, in2->ny);
	}
	out=X(new)(in1->nx, in1->ny+in2->ny);
	memcpy(out->p, in1->p, in1->nx*in1->ny*sizeof(T));
	memcpy(out->p+in1->nx*in1->ny,in2->p,in2->nx*in2->ny*sizeof(T));
    }else{
	error("Invalid dim\n");
    }
    return out;
}


/**
   duplicate a X(mat) array
*/
X(mat) *X(dup)(const X(mat) *in){
    X(mat) *out=NULL;
    X(cp)(&out, in);
    return out;
}

/**
   copy the values from one X(mat) to another.
*/
void X(cp)(X(mat) **out0, const X(mat) *in){
    if(in && in->nx && in->ny){
	if(!*out0) 
	    *out0=X(new)(in->nx, in->ny);
	else{
	    assert(in->nx==(*out0)->nx && in->ny == (*out0)->ny);
	}
	X(mat) *out=*out0;
	if(out->p!=in->p){
	    memcpy(out->p, in->p, in->nx*in->ny*sizeof(T));
	}
    }else{
	X(zero)(*out0);
    }
}

/**
   transpose a X(mat) object
*/
X(mat) *X(trans)(const X(mat) *A){
    if(!A) return NULL;
    X(mat) *B=X(new)(A->ny,A->nx);
    if(A->nx==1 || A->ny==1){
	memcpy(B->p, A->p, A->nx*A->ny*sizeof(T));
    }else{
	PMAT(B,Bp);
	PMAT(A,Ap);
	for(int ix=0; ix<A->nx; ix++){
	    for(int iy=0; iy<A->ny; iy++){
		Bp[ix][iy]=Ap[iy][ix];
	    }
	}
    }
    return B;
}
/**
   set values of each element in a X(mat) to val.
*/
void X(set)(X(mat) *A, const T val){
    if(A){
	for(long i=0; i<A->nx*A->ny; i++){
	    A->p[i]=val;
	}
    }
}

/**
   display a X(mat) matrix.
*/
void X(show)(const X(mat) *A, const char *format, ...){
    if(!A) return;
    format2fn;
    info2("Displaying content of %s:\n",fn);
    PMAT(A,p);
    int colmax=10;
    int iset,i,j;
    int nset=(A->ny+colmax-1)/colmax;
    for(iset=0; iset<nset; iset++){
	int ncol=(iset+1)*colmax;
	if(ncol>A->ny) ncol=A->ny;
	printf("Cols %d to %d\n", iset, ncol-1);
	for(j=0; j<A->nx; j++){
	    for(i=iset*colmax; i<ncol; i++){
		PRINT(p[i][j]);
	    }
	    printf("\n");
	}
    }
}

/**
   create sum of all the elements in A.
*/
T X(sum)(const X(mat) *A){
    T v=0;
    if(A){
	T *restrict p=A->p;
	/**
	   Loops like this will only be vectorized with -ffast-math because
	   different orders of accumulation give different results for floating
	   point numbers.
	*/
	for(int i=0; i<A->nx*A->ny; i++){
	    if(isfinite(p[i])){
		v+=p[i];
	    }
	}
    }
    return v;
}
/**
   trace
*/
T X(trace)(const X(mat)*A){
    T trace=0;
    long n=MIN(A->nx, A->ny);
    for(long i=0; i<n; i++){
	trace+=A->p[i*(1+A->nx)];
    }
    return trace;
}


/**
   shift frequency components by n/2
*/
void X(fftshift)(X(mat) *A){
    long i;
    const long nx=A->nx;
    const long ny=A->ny;
    if((nx&1)==1){
	warning("nx=%ld is not multiple of 2\n", nx);
    }
    const long nx2=nx/2;
    const long ny2=ny/2;
    const long nx2d=nx2*sizeof(T);
    T *tmp=(T*)malloc(nx2d);
    T *data=A->p;
    if(ny==1){
	memcpy(tmp,data,nx2d);
	memcpy(data,data+nx2,nx2d);
	memcpy(data+nx2,tmp,nx2d);
    }else{
	assert((ny&1)==0);
	for(i=0; i<ny2; i++){
	    memcpy(tmp,data+i*nx,nx2d);
	    memcpy(data+i*nx,data+(i+ny2)*nx+nx2, nx2d);
	    memcpy(data+(i+ny2)*nx+nx2,tmp, nx2d);
	    memcpy(tmp,data+i*nx+nx2,nx2d);
	    memcpy(data+i*nx+nx2,data+(i+ny2)*nx, nx2d); 
	    memcpy(data+(i+ny2)*nx,tmp, nx2d);
	}
    }
    
    free(tmp);
}

/**
   reorder B
   and embed/crop into center of A 
   \verbatim
   4 * * 3
   * * * *
   * * * *
   2 * * 1
   \endverbatim
   to
   \verbatim
   1 2 
   3 4
   \endverbatim
*/
void X(cpcorner2center)(X(mat) *A, const X(mat)*B){
    const long nx=A->nx;
    const long ny=A->ny;
    T *Ap=A->p;
    const long ninx=B->nx;
    const long niny=B->ny;
    if(nx>ninx || ny>niny){
	memset(Ap, 0, sizeof(T)*nx*ny);
    }
    const T * Bp=B->p;
    assert((nx&1)==0 && (ny&1)==0 && (ninx&1)==0 && (niny&1)==0);

    const int ny2=MIN(ny,niny)>>1;
    const int nx2=MIN(nx,ninx)>>1;
    const int xskip=nx/2-nx2;
    const int yskip=ny/2-ny2;
    T* Ap0=Ap+yskip*nx+xskip;
    const int nx2d=nx2*sizeof(T);
    for(int i=0; i<ny2; i++){
	memcpy(Ap0+i*nx, Bp+(niny-ny2+i)*ninx+(ninx-nx2),nx2d); 
	memcpy(Ap0+i*nx+nx2, Bp+(niny-ny2+i)*ninx, nx2d); 
	memcpy(Ap0+(i+ny2)*nx, Bp+i*ninx+(ninx-nx2), nx2d); 
	memcpy(Ap0+(i+ny2)*nx+nx2, Bp+i*ninx, nx2d); 
    }
}

/**
   cyclic shift A by nx and ny to B.
   \verbatim
   4   3     1   2 
      
   2   1 to  3   4
   \endverbatim
*/
void X(shift)(X(mat) **B0, const X(mat) *A, int sx, int sy){
    if(!*B0){
	*B0=X(new)(A->nx, A->ny);
    }
    X(mat) *B=*B0;

    const int nx=A->nx; 
    const int ny=A->ny;
    sx=sx%nx; if(sx<0) sx+=nx;
    sy=sy%ny; if(sy<0) sy+=ny;
    if(sx!=0 || sy!=0){
	int dy=ny-sy;
	int dx=nx-sx;
	for(int iy=0; iy<sy; iy++){
	    memcpy(B->p+iy*nx, A->p+(dy+iy)*nx+dx, sx*sizeof(T));/*3 */
	    memcpy(B->p+iy*nx+sx, A->p+(dy+iy)*nx, dx*sizeof(T));/*4 */
	}
	for(int iy=sy; iy<ny; iy++){
	    memcpy(B->p+iy*nx, A->p+(iy-sy)*nx+dx, sx*sizeof(T));/*1 */
	    memcpy(B->p+iy*nx+sx, A->p+(iy-sy)*nx, dx*sizeof(T));/*2 */
	}
    }else{
	memcpy(B->p, A->p, sizeof(T)*nx*ny);
    }
}

/**
   create an new X(cell) similar to A in shape.
   When a cell is empty, it is created with a (0,0) array and cannot be overriden.
*/
X(cell) *X(cellnew2)(const X(cell) *A){
    X(cell) *out=X(cellnew)(A->nx, A->ny);
    long tot=0;
    for(long i=0; i<A->nx*A->ny; i++){
	if(!isempty(A->p[i])){
	    tot+=A->p[i]->nx*A->p[i]->ny;
	}
    }
    out->m=X(new)(tot,1);
    tot=0;
    for(int i=0; i<A->nx*A->ny; i++){
	if(!isempty(A->p[i])){
	    out->p[i]=X(new_ref)(A->p[i]->nx, A->p[i]->ny, out->m->p+tot);
	    tot+=A->p[i]->nx*A->p[i]->ny;
	}else{
	    out->p[i]=X(new)(0,0);//place holder to avoid been overriden.
	}
    }
    return out;
}
/**
   Create an new X(cell) with X(mat) specified. Each block is stored continuously in memory.
 */
X(cell) *X(cellnew3)(long nx, long ny, long *nnx, long *nny){
    X(cell) *out=X(cellnew)(nx,ny);
    long tot=0;
    for(long i=0; i<nx*ny; i++){
	tot+=nnx[i]*(nny?nny[i]:1);
    }
    if(!tot) return NULL;
    out->m=X(new)(tot,1);
    tot=0;
    for(long i=0; i<nx*ny; i++){
	out->p[i]=X(new_ref)(nnx[i], (nny?nny[i]:1), out->m->p+tot);
	tot+=nnx[i]*(nny?nny[i]:1);
    }
    return out;
}

/**
   check the size of matrix if exist. Otherwise create it
*/
void X(cellinit)(X(cell)**A, long nx, long ny){
    if(*A){
	assert((*A)->nx == nx && (*A)->ny==ny);
    }else{
	*A=X(cellnew)(nx, ny);
    }
}
/**
   creat a X(cell) reference an existing X(cell) by referencing the
   elements.
*/
X(cell) *X(cellref)(const X(cell) *in){
    if(!in) return NULL;
    X(cell) *out=X(cellnew)(in->nx, in->ny);
    if(in->m){
	out->m=X(ref)(in->m);
	for(int i=0; i<in->nx*in->ny; i++){
	    out->p[i]=X(new_ref)(in->p[i]->nx, in->p[i]->ny, in->p[i]->p);
	}
    }else{
	for(int i=0; i<in->nx*in->ny; i++){
	    out->p[i]=X(ref)(in->p[i]);
	}
    }
    return out;
}

/**
   duplicate a X(cell) object.
*/
X(cell) *X(celldup)(const X(cell) *in){
    X(cell) *out=NULL;
    X(cellcp)(&out, in);
    return out;
}

/**
   copy the values from one X(cell) to another.
*/
void X(cellcp)(X(cell)** out0, const X(cell) *in){
    if(in){
	if(!*out0){
	    *out0=X(cellnew2)(in);
	}else{
	    assert((*out0)->nx==in->nx && (*out0)->ny==in->ny);
	}
	X(cell)* out=*out0;
	for(int i=0; i<in->nx*in->ny; i++){
	    X(cp)(&out->p[i], in->p[i]);
	}
    }else{/*if input is empty, zero the output. do not free. */
	X(cellzero)(*out0);
    }
}

/**
   setting all elements of a X(cell) to zero.
*/
void X(cellzero)(X(cell) *dc){
    if(dc){
	for(int ix=0; ix<dc->nx*dc->ny; ix++){
	    X(zero)(dc->p[ix]);
	}
    }
}

/**
   setting all elements of a X(cell) to alpha.
*/
void X(cellset)(X(cell)*dc, T val){
    if(dc){
	for(int ix=0; ix<dc->nx*dc->ny; ix++){
	    X(set)(dc->p[ix], val);
	}
    }
}

/**
   transpose a X(cell) object
*/
X(cell) *X(celltrans)(const X(cell) *A){
    if(!A) return NULL;
    X(cell) *B=X(cellnew)(A->ny, A->nx);
    X(mat)* (*Bp)[B->nx]=(void*)B->p;
    X(mat)* (*Ap)[A->nx]=(void*)A->p;
    for(int ix=0; ix<A->nx; ix++){
	for(int iy=0; iy<A->ny; iy++){
	    Bp[ix][iy]=X(trans)(Ap[iy][ix]);
	}
    }
    return B;
}

/**
   Obtain the dimensions.
*/
static void X(celldim)(const X(cell) *A, long *nx, long *ny,
		       long **nxs, long **nys){
    X(mat) *(*Ap)[A->nx] = (X(mat) *(*)[A->nx])A->p;
    *nxs=calloc(A->nx, sizeof(long));
    *nys=calloc(A->ny, sizeof(long));
    *nx=0;
    *ny=0;
    for(long ix=0; ix<A->nx; ix++){
	for(long iy=0; iy<A->ny; iy++){
	    if(!isempty(Ap[iy][ix])){
		*nx+=Ap[iy][ix]->nx;
		(*nxs)[ix]=Ap[iy][ix]->nx;
		break;
	    }
	}
    }
    for(long iy=0; iy<A->ny; iy++){
	for(long ix=0; ix<A->nx; ix++){
	    if(!isempty(Ap[iy][ix])){
		*ny+=Ap[iy][ix]->ny;
		(*nys)[iy]=Ap[iy][ix]->ny;
		break;
	    }
	}
    }
}

/**
   reduce nx*ny cell matrix to 1*ny if dim=1 and nx*1 if dim=2
*/
X(cell) *X(cellreduce)(const X(cell)*A, int dim){
    if(!A) return NULL;
    X(cell)* out=NULL;
    long nx, ny, *nxs, *nys;
    X(celldim)(A, &nx, &ny, &nxs, &nys);
    if(nx==0 || ny==0) return NULL;
    PCELL(A,pA);
    if(dim==1){
	out=X(cellnew)(1, A->ny);
	for(long iy=0; iy<A->ny; iy++){
	    if(nys[iy]==0) continue;
	    out->p[iy]=X(new)(nx,nys[iy]);
	    for(long icol=0; icol<nys[iy]; icol++){
		long kr=0;
		for(long ix=0; ix<A->nx; ix++){
		    if(!isempty(pA[iy][ix])){
			memcpy(out->p[iy]->p+icol*nx+kr,
			       pA[iy][ix]->p+icol*nxs[ix],
			       nxs[ix]*sizeof(T));
		    }			
		    kr+=nxs[ix];
		}
	    }
	}
    }else if(dim==2){
	out=X(cellnew)(A->nx,1);
	for(long ix=0; ix<A->nx; ix++){
	    if(nxs[ix]==0) continue;
	    out->p[ix]=X(new)(nxs[ix],ny);
	    long kr=0;
	    for(long iy=0; iy<A->ny; iy++){
		if(!isempty(pA[iy][ix])){
		    memcpy(out->p[ix]->p+kr*nxs[ix],
			   pA[iy][ix]->p,
			   nxs[ix]*nys[iy]*sizeof(T));
		}
		kr+=nys[iy];
	    }
	}
    }else{
	error("Invalid dim=%d\n",dim);
    }
    free(nxs);
    free(nys);
    return out;
}

/**
   concatenate two cell matrices along dimenstion 'dim'
*/
X(cell) *X(cellcat)(const X(cell) *A, const X(cell) *B, int dim){
    if(!A){
	if(!B){
	    return NULL;
	}else{
	    return X(celldup)(B);
	}
    }else if(!B){
	return X(celldup)(A);
    }

    X(cell) *out=NULL;
    PCELL(A,pA);
    PCELL(B,pB);

    if(dim==1){
	/*along x. */
	if(A->ny!=B->ny){
	    error("Mismatch: A is (%ld, %ld), B is (%ld, %ld)\n",
		  A->nx, A->ny, B->nx, B->ny);
	}
	out=X(cellnew)(A->nx+B->nx, A->ny);
	PCELL(out,pout);
	for(long iy=0; iy<A->ny; iy++){
	    for(long ix=0; ix<A->nx; ix++){
		pout[iy][ix]=X(dup)(pA[iy][ix]);
	    }
	    for(long ix=0; ix<B->nx; ix++){
		pout[iy][ix+A->nx]=X(dup)(pB[iy][ix]);
	    }
	}
    }else if(dim==2){
	/*along y. */
	if(A->nx!=B->nx){
	    error("Mismatch. A is (%ld, %ld), B is (%ld, %ld)\n", 
		  A->nx, A->ny, B->nx, B->ny);
	}
	out=X(cellnew)(A->nx, A->ny+B->ny);
	PCELL(out,pout);
	for(long iy=0; iy<A->ny; iy++){
	    for(long ix=0; ix<A->nx; ix++){
		pout[iy][ix]=X(dup)(pA[iy][ix]);
	    }
	}
	for(long iy=0; iy<B->ny; iy++){
	    for(long ix=0; ix<B->nx; ix++){
		pout[iy+A->ny][ix]=X(dup)(pB[iy][ix]);
	    }
	}
    }else{
	error("Invalid dim\n");
    }
    return out;
}

/**
   concatenate coresponding elements of each X(cell). They must
   have the same shape.
*/
X(cell) *X(cellcat_each)(const X(cell) *A, const X(cell) *B, int dim){
    if(!A){
	if(!B){
	    return NULL;
	}else{
	    return X(celldup)(B);
	}
    }else if(!B){
	return X(celldup)(A);
    }
    if(A->nx!=B->nx || A->ny!=B->ny){
	error("Mismatch: (%ld %ld), (%ld %ld)\n",A->nx, A->ny, B->nx, B->ny);
    }
    X(cell) *out=X(cellnew)(A->nx, A->ny);
    for(long ix=0; ix<A->nx*A->ny; ix++){
	out->p[ix]=X(cat)(A->p[ix], B->p[ix], dim);
    }
    return out;
}

/**
   drop empty rows or columns. (size of *A0 is changed.
*/
void X(celldropempty)(X(cell) **A0, int dim){
    X(cell) *A=*A0;
    if(!A) return;
    PCELL(A,pA);
    if(dim==1){
	/*drop rows */
	int keep[A->nx];
	int ndrop=0;
	for(int ix=0; ix<A->nx; ix++){
	    keep[ix]=0;
	    for(int iy=0; iy<A->ny; iy++){
		if(!isempty(pA[iy][ix])){
		    keep[ix]=1;
		    break;
		}
	    }
	    if(keep[ix]==0){
		ndrop++;
	    }
	}
	if(ndrop){
	    if(ndrop==A->nx){
		X(cellfree)(A);
		*A0=NULL;
	    }else{
		X(cell) *B=calloc(1, sizeof(X(cell)*));
		B->p=calloc((A->nx-ndrop)*A->ny, sizeof(X(mat)*));
		PCELL(B,pB);
		int count=0;
		for(int ix=0; ix<A->nx; ix++){
		    if(keep[ix]){
			if(count!=ix){
			    for(int iy=0; iy<A->ny; iy++){
				pB[iy][count]=pA[iy][ix];
			    }
			}
			count++;
		    }else{
			warning("row %d dropped\n", ix);
		    }
		}
		free(A->p); free(A);
		A=B;
	    }
	}
    }else if(dim==2){
	/*drop cols */
	int count=0;
	for(int iy=0; iy<A->ny; iy++){
	    int keep=0;
	    for(int ix=0; ix<A->nx; ix++){
		if(!isempty(pA[iy][ix])){
		    keep=1;
		    break;
		}
	    }
	    if(keep){
		for(int ix=0; ix<A->nx; ix++){
		    pA[count][ix]=pA[iy][ix];
		}
		count++;
	    }else{
		/*warning("Col %d dropped\n", iy); */
	    }
	}
	A->ny=count;
	if(count==0){
	    X(cellfree)(A);
	    *A0=NULL;
	}else{
	    A->p=realloc(A->p,sizeof(X(mat)*)*A->ny*A->nx);
	}
    }else{
	error("Invalid dim: %d\n",dim);
    }
 
}

/**
   Convert a block matrix to a matrix.
*/
X(mat) *X(cell2m)(const X(cell) *A){
    if(A->nx*A->ny==1){
	return X(dup)(A->p[0]);
    }
    PCELL(A, Ap);
    //X(mat) *(*Ap)[A->nx] = (X(mat) *(*)[A->nx])A->p;
    long nx,ny,*nxs,*nys;
    X(celldim)(A,&nx,&ny,&nxs,&nys);
    X(mat) *out=X(new)(nx,ny);
    long jcol=0;
    for(long iy=0; iy<A->ny; iy++){
	for(long icol=0; icol<nys[iy]; icol++){
	    long kr=0;
	    for(long ix=0; ix<A->nx; ix++){
		if(!isempty(Ap[iy][ix])){
		    memcpy(out->p+((icol+jcol)*nx+kr),
			   Ap[iy][ix]->p+icol*nxs[ix],
			   nxs[ix]*sizeof(T));
		}
		kr+=nxs[ix];
	    }
	}
	jcol+=nys[iy];
    }
    free(nxs);
    free(nys);
    return out;
}

/**
   convert a vector to cell using dimensions specified in dims. Reference the vector
*/
X(cell)* X(2cellref)(const X(mat) *A, long*dims, long ndim){
    long nx=0;
    for(int ix=0; ix<ndim; ix++){
	nx+=dims[ix];
    }
    if(A->ny!=1){
	error("Use d2cell2 instead for non vectors\n");
    }
    if(nx!=A->nx ){
	error("Shape doesn't agree. nx=%ld, nx=%ld\n", nx,A->nx);
    }
    long kr=0;
    X(cell) *B=X(cellnew)(ndim,1);
    B->m=X(ref)(A);
    for(long ix=0; ix<ndim; ix++){
	B->p[ix]=X(new_ref)(dims[ix],1,A->p+kr);
	kr+=dims[ix];
    }
    return B;
}

/**
   make A a cell array using shape information from ref if *B is NULL
*/
void X(2cell)(X(cell) **B, const X(mat) *A, const X(cell) *ref){
    long nx,ny,*nxs,*nys;
    if(*B) ref=*B;/*use B as reference. */
    X(celldim)(ref, &nx, &ny, &nxs, &nys);
    if(nx!=A->nx || ny!=A->ny){
	error("Shape doesn't agree. Reference is %ldx%ld but input is %ldx%ld\n",
	      nx,ny,A->nx,A->ny);
    }
    if(!*B){
	*B=X(cellnew2)(ref);
	PCELL((*B),Bp);
	for(long iy=0; iy<ref->ny; iy++){
	    for(long ix=0; ix<ref->nx; ix++){
		Bp[iy][ix]=X(new)(nxs[ix],nys[iy]);
	    }
	}
    }
    PCELL((*B),Bp);
    long jcol=0;
    for(long iy=0; iy<ref->ny; iy++){
	for(long icol=0; icol<nys[iy]; icol++){
	    long kr=0;
	    for(long ix=0; ix<ref->nx; ix++){
		if(nxs[ix]>0){
		    memcpy(Bp[iy][ix]->p+icol*Bp[iy][ix]->nx,
			   A->p+((icol+jcol)*nx+kr),
			   nxs[ix]*sizeof(T));
		    kr+=nxs[ix];
		}
	    }
	}
	jcol+=nys[iy];
    }
    free(nxs);
    free(nys);
}



/**
   Create a new sub cell matrix of nx*ny starting from(sx,sy)
*/
X(cell) *X(cellsub)(const X(cell) *in, long sx, long nx, long sy, long ny){
    if(nx<=0){
	nx=in->nx-sx;
    }
    if(ny<=0){
	ny=in->ny-sy;
    }
    X(cell)*out=X(cellnew)(nx, ny);
    if(sx+nx>in->nx || sy+ny>in->ny){
	error("Invalid parameter range\n");
    }
    PCELL(in, pin);
    PCELL(out, pout);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    pout[iy][ix]=pin[iy+sy][ix+sx];
	}
    }
    return out;
}
