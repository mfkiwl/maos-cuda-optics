/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "curmat.h"
#include "cucmat.h"
#include "utils.h"
#include "kernel.h"

void curset(curmat *A, float alpha, cudaStream_t stream){
    if(A && A->p){
	set_do<<<DIM(A->nx*A->ny,256),0,stream>>>(A->p, alpha, A->nx*A->ny);
    }
}

void curcp(curmat **out, const curmat *in, cudaStream_t stream){
    if(!in){
	curzero(*out, stream);
    }else{
	if(!*out){
	    *out=curnew(in->nx, in->ny);
	}else{
	    assert((*out)->nx * (*out)->ny==in->nx * in->ny);
	}
	cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(float), cudaMemcpyDeviceToDevice, stream);
    }
}

/**
   add a vector to another, scaled by alpha and beta. all in device memory.
   a=a*alpha+b*beta;
*/
__global__ void add_do(float *restrict a, float alpha, 
		       const float * b, float beta, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=a[i]*alpha+b[i]*beta;
    }
}

__global__ static void scale_do(float *restrict in, int n, float alpha){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	in[i]*=alpha;
    }
}

/**
   out=out*beta+in*alpha;
*/
void curadd(curmat **out, float alpha, curmat *in, float beta, cudaStream_t stream){
    if(!in) return;
    if(!*out || fabsf(alpha)<1e-7){
	curcp(out, in, stream);
	scale_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	    ((*out)->p, in->nx*in->ny, beta);
    }else{
	add_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	    ((*out)->p, alpha, in->p, beta, in->nx*in->ny);
    }
}
__global__ void addcabs2_do(float *restrict a, float alpha, 
			    const fcomplex *restrict b, float beta, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=a[i]*alpha+CABS2(b[i])*beta;
    }
}
/**
   out=out*beta+abs2(in)*alpha;
*/
void curaddcabs2(curmat **out, float alpha, cucmat *in, float beta, cudaStream_t stream){
    if(!*out){
	*out=curnew(in->nx,in->ny);
    }else if(fabsf(alpha)<1e-7){
	curzero(*out, stream);
    }
    addcabs2_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	((*out)->p, alpha, in->p, beta, in->nx*in->ny);
}
void curscale(curmat *in, float alpha, cudaStream_t stream){
    if(!in) return;
    if(fabsf(alpha)<1e-5) {
	curzero(in, stream);
    }else if(fabsf(alpha-1.f)>1e-5){
	int n=in->nx*in->ny;
	scale_do<<<DIM(n,256), 0, stream>>>(in->p, n, alpha); 
    }
}

/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmm(curmat **C, float alpha, const curmat *A, const curmat *B, char trans[2], float beta, cublasHandle_t handle){
    int m,n,k,k2;
    cublasOperation_t transa, transb;
    if(trans[0]=='t'){
	m=A->ny;
	k=A->nx;
	transa=CUBLAS_OP_T;
    }else{
	m=A->nx;
	k=A->ny;
	transa=CUBLAS_OP_N;
    }
    if(trans[1]=='t'){
	n=B->nx;
	k2=B->ny;
	transb=CUBLAS_OP_T;
    }else{
	n=B->ny;
	k2=B->nx;
	transb=CUBLAS_OP_N;
    }
    if(!*C){
	*C=curnew(m,n);
    }else{
	assert((*C)->nx==m && (*C)->ny==n);
    }
    assert(k==k2);
    DO(cublasSgemm(handle, transa, transb, m,n,k,
		       &beta, A->p, A->nx, B->p, B->nx, &alpha, (*C)->p, (*C)->nx));
}
/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmv(curmat **C, float alpha, const curmat *A, const curmat *B, char trans, float beta, cublasHandle_t handle){
    if(!*C){
	*C=curnew(trans=='t'?A->ny:A->nx, 1);
    }else{
	assert((*C)->nx==(trans=='t'?A->ny:A->nx) && (*C)->ny==1);
    }
    cublasSgemv(handle, trans=='t'?CUBLAS_OP_T:CUBLAS_OP_N, A->nx, A->ny, &beta, A->p, A->nx, B->p, 1, &alpha, (*C)->p, 1);
}

cuspcell* cuspcellnew(int nx, int ny){
    cuspcell *out=(cuspcell*)calloc(1, sizeof(cuspcell));
    out->p=(cusp**)calloc(nx*ny, sizeof(void*));
    out->nx=nx;
    out->ny=ny;
    return out;
}

/*
  A=A*beta+B*alpha;
*/
void curcelladd(curcell **A, float beta, const curcell *B, float alpha, cudaStream_t stream){
    if(!B) return;
    if(!*A){
	*A=curcellnew(B);
    }else{
	assert((*A)->nx==B->nx && (*A)->ny==B->ny);
    }
    if((*A)->m && B->m){
	curadd(&(*A)->m, beta, B->m, alpha, stream);
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    curadd(&((*A)->p[i]), beta, B->p[i], alpha,stream);
	}
    }
}


/**
   add a beta to a vector. 
*/
__global__ void adds_do(float *vec, float beta, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	vec[i]+=beta;
    }
}
void curadds(curmat *A, float beta, cudaStream_t stream){
    const int n=A->nx*A->ny;
    adds_do<<<DIM(n, 256), 0, stream>>>(A->p, beta, n);
}
/**
   add a vector to another, scaled by alpha and beta. all in device memory.
   a=a+b*alpha*beta;
*/
__global__ void add2_do(float *restrict a, const float * b, const float *restrict b_sc1, float b_sc2, int n){
    float alpha=*b_sc1*b_sc2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]+=b[i]*alpha;
    }
}


/**
   out=out+in*alpha; beta, alpha lives in device memory.
*/
void curadd2(curmat **out, const curmat *in, float *alpha, float alpha2, cudaStream_t stream){
    if(!*out){
	*out=curnew(in->nx, in->ny);
    }
    add2_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	((*out)->p, in->p, alpha, alpha2, in->nx*in->ny);
}


/**
   A=A*beta+B*alpha; beta, alpha lives in device memory.
*/
void curcelladd2(curcell **A, const curcell *B, float* alpha, float alpha2, cudaStream_t stream){
    if(!B) return;
    if(!*A){
	*A=curcellnew(B);
    }else{
	assert((*A)->nx==B->nx && (*A)->ny==B->ny);
    }
    if((*A)->m && B->m){
	curadd2(&(*A)->m, B->m, alpha, alpha2, stream);
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    curadd2((*A)->p+i, B->p[i], alpha, alpha2,  stream);
	}
    }
}

__global__ static void add3_do(float *restrict a, const float *restrict a_sc, const float a_sc2, const float *restrict b, int n){
    float alpha=*a_sc*a_sc2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=a[i]*alpha+b[i];
    }
}
/**
   out=out*beta+in; beta, alpha lives in device memory.
*/
void curadd3(curmat **out, float *beta, const curmat *in, cudaStream_t stream){
    if(!*out){
	*out=curnew(in->nx, in->ny);
    }
    add3_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	((*out)->p, beta, 1.f, in->p, in->nx*in->ny);
}

/**
   A=A*beta+B*alpha; beta, alpha lives in device memory.
*/
void curcelladd3(curcell **A, float* beta, const curcell *B, cudaStream_t stream){
    if(!B) return;
    if(!*A){
	*A=curcellnew(B);
    }else{
	assert((*A)->nx==B->nx && (*A)->ny==B->ny);
    }
    if((*A)->m && B->m){
	curadd3(&(*A)->m, beta, B->m, stream);
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    curadd3(&((*A)->p[i]), beta, B->p[i],  stream);
	}
    }
}
__global__ void inn_do(float *restrict res, const float *a, const float *b, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sb[threadIdx.x]+=a[i]*b[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    sb[threadIdx.x]+=sb[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	atomicAdd(res, sb[0]);
    }
}
__global__ void sum_do(float *restrict res, const float *a, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sb[threadIdx.x]+=a[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    sb[threadIdx.x]+=sb[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	atomicAdd(res, sb[0]);
    }
}

inline static void inn_wrap(float *res, const float * a, const float * b, 
			    const int n, cudaStream_t stream){
    inn_do<<<DIM(n, DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>(res,a,b,n);
}
inline static void sum_wrap(float *res, const float * a, const int n, cudaStream_t stream){
    sum_do<<<DIM(n, DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>(res,a,n);
}
float curinn(const curmat *a, const curmat *b, cudaStream_t stream){
    float *res;
    cudaMalloc(&res, sizeof(float));
    cudaMemsetAsync(res, 0,sizeof(float), stream);
    inn_wrap(res, a->p, b->p, a->nx*a->ny, stream);
    float out;
    cudaMemcpyAsync(&out, res, sizeof(float), cudaMemcpyDeviceToHost, stream);
    CUDA_SYNC_STREAM;
    return out;
}

/**
   Sum all the elements in an array.
 */
void cursum2(float *restrict res, const curmat *a, cudaStream_t stream){
    cudaMemsetAsync(res, 0,sizeof(float), stream);
    sum_wrap(res, a->p, a->nx*a->ny, stream);
}

__global__ void max_do(float *restrict res, const float *a, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	if(sb[threadIdx.x]<a[i]) sb[threadIdx.x]=a[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    if(sb[threadIdx.x]<sb[threadIdx.x+step]){
		sb[threadIdx.x]=sb[threadIdx.x+step];
	    }
	}
    }
    if(threadIdx.x==0){
	*res=sb[0];
    }
}
/**
   Find the maximum value
*/
float curmax(const curmat *a, cudaStream_t stream){
    float *res;
    cudaMallocHost(&res, sizeof(float));
    int n=a->nx*a->ny;
    max_do<<<DIM(n, DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>> (res, a->p, n);
    CUDA_SYNC_STREAM;
    float out=res[0];
    cudaFreeHost(res);
    return out;
}

/**
   Scale elements
*/
void curcellscale(curcell *A, float alpha, cudaStream_t stream){
    if(!A) return;
    if(A->m){
	curscale(A->m, alpha, stream);
    }else{
	for(int i=0; i<A->nx*A->ny; i++){
	    curscale(A->p[i], alpha, stream);
	}
    }
}