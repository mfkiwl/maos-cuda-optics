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
#ifndef AOS_CUDA_TYPES_H
#define AOS_CUDA_TYPES_H
#include "../math/array.h"
#include "common.h"
#include "kernel.h"
class nonCopyable{
private:
    nonCopyable& operator=(nonCopyable&);
    nonCopyable(const nonCopyable&);
protected:
    nonCopyable(){}
};
/**
   Cpu, Pinned, and GPU template classes are created to allow different memory
   allocation.
*/


//Pinned (page locked) CPU memory
template<typename T>
class Pinned:public Cpu<T>{
public:
    void *operator new[](size_t size){
	void *p=0;
	cudaMallocHost(&p, size);
	memset(p, 0, size);
	return p;
    }
    void operator delete[](void*p){
	cudaFreeHost(p);
    }
};

//Gpu memory
template<typename T>
class Gpu{
    T val;
public:
    void *operator new[](size_t size){
	void *p=0;
	DO(cudaMalloc(&p, size));
	DO(cudaMemset(p, 0, size));
	return p;
    }
    void operator delete[](void*p){
	DO(cudaFree(p));
    }
    static void zero(T *p, size_t size, cudaStream_t stream=(cudaStream_t)-1){
	if(p){
	    if(stream==(cudaStream_t)-1){
		DO(cudaMemset(p, 0, size));
	    }else{
		DO(cudaMemsetAsync(p, 0, size*sizeof(T), stream));
	    }
	}
    }
    static void memcpy(T *pout, T* pin, size_t size){
	cudaMemcpy(pout, pin, size, cudaMemcpyHostToDevice);
    }
};
/**
   A special Cell with pm.
 */
template <typename T>
class CuCell:public Cell<T, Gpu>{
    typedef Cell<T, Gpu> Parent;
public:
    using Parent::operator();
    using Parent::p;
    using Parent::nx;
    using Parent::ny;
    using Parent::M;
    Array<T*,Gpu>pm;/*contains the data pointer of each cell.*/
    Array<T*,Pinned>pm_cpu;/*contains the data pointer of each cell.*/
    
    void p2pm(cudaStream_t stream=(cudaStream_t)-1){
	if(nx && ny){
	    pm_cpu.init(nx,ny);
	    pm.init(nx,ny);
	    for(long i=0; i<nx*ny; i++){
		pm_cpu[i]=p[i]();
	    }
	    if(stream==(cudaStream_t)-1){
		cudaMemcpy(pm(), pm_cpu(), sizeof(T*)*nx*ny,cudaMemcpyHostToDevice);
	    }else{
		cudaMemcpyAsync(pm(), pm_cpu(), sizeof(T*)*nx*ny,cudaMemcpyHostToDevice, stream);
	    }
	}
    }
  
    CuCell(long nxi=0, long nyi=1):Parent(nxi, nyi){
	p2pm();
    }
    CuCell(long nxi, long nyi, long mx, long my, T *pin=NULL):Parent(nxi,nyi, mx, my, pin){
	p2pm();
    }
    //template <typename L>
    CuCell(const long nxi, const long nyi, int *mx, int *my, T *pin=NULL):Parent(nxi,nyi,mx, my, pin){
	p2pm();
    }
    //template <typename L>
    CuCell(const long nxi, const long nyi, long *mx, long *my, T *pin=NULL):Parent(nxi,nyi,mx, my, pin){
	p2pm();
    }
    CuCell& operator=(const CuCell<T>&A){
	Parent::operator=(A);
	p2pm();
	return *this;
    }
    void init(long nxi=0, long nyi=1){
	Parent::init(nxi, nyi);
	pm_cpu.init();
	pm.init();
    }
    void Replace(T *pnew, cudaStream_t stream=(cudaStream_t)-1){
	/*replace the data with a new set. we don't own the pnew.*/
	Parent::Replace(pnew);
	p2pm(stream);
    }
};
typedef class Array<int, Gpu>   cuimat;
typedef class Array<Real, Gpu>   curmat;
typedef class Array<Comp, Gpu>   cucmat;
typedef class CuCell<Real>  curcell;
typedef class CuCell<Comp>  cuccell;
enum TYPE_SP{
    SP_CSC,//compressed sparse column major. 
    SP_CSR,//compressed sparse row major. 
};
class cusp{
    int *p;
    int *i;
    Real *x;
    int nx;
    int ny;
    int nzmax;
    int *nref;
    enum TYPE_SP type;
public:
    enum TYPE_SP Type()const{
	return type;
    }
    long Nx()const{
	return nx;
    }
    long Ny()const{
	return ny;
    }
    long Nzmax()const{
	return nzmax;
    }
    int *Pp(){
	return p;
    }
    const int *Pp() const{
	return p;
    }
    int *Pi(){
	return i;
    }
    const int *Pi() const{
	return i;
    }

    Real *Px(){
	return x;
    }
    const Real *Px() const{
	return x;
    }
    cusp():p(0),i(0),x(0),nx(0),ny(0),nzmax(0),nref(0),type(SP_CSC){}
    cusp(const dsp *in, int tocsr);
    cusp(const cusp&in):p(in.p),i(in.i),x(in.x),nx(in.nx),ny(in.ny),nzmax(in.nzmax),nref(in.nref),type(in.type){
	if(nref) nref[0]++;
    }
    cusp &operator=(const cusp&in){
	if(this!=&in){
	    reset();
	    p=in.p;
	    i=in.i;
	    x=in.x;
	    nx=in.nx;
	    ny=in.ny;
	    nzmax=in.nzmax;
	    nref=in.nref;
	    if(nref) nref[0]++;
	    type=in.type;
	}
	return *this;
    }
    ~cusp(){
	reset();
    }
    void reset(){
	if(nref && !atomicadd(nref, -1)){
	    cudaFree(p);
	    cudaFree(i);
	    cudaFree(x);
	    delete nref;
	}
	p=0; i=0; x=0; nref=0;
    }
 
    void trans();/*covnert to CSR mode by transpose*/
    /*cusp* ref(void){
      if(nref) atomicadd(nref, 1);
      cusp* res=(cusp*)malloc(sizeof(cusp));
      memcpy(res, this, sizeof(*this));
      return res;
      }*/
    operator bool()const{
	return nx && ny;
    }
};
typedef Array<cusp> cuspcell;
typedef Array<curcell> curccell;
typedef Array<curccell> curcccell;
void cp2gpu(curmat &dest, const loc_t *src);
class culoc_t{
private:
    curmat p;/*in device. */
    Real dx;
    Real dy;
public:
    Real2* operator()(){
	return (Real2*)(p());
    }
    const Real2*operator()()const{
	return (Real2*)(p());
    }
    long Nloc()const{
	return p.Ny();
    }
    Real Dx()const{
	return dx;
    }
    Real Dy()const{
	return dy;
    }
    //No need custom copy assignment operator or copy constructor.
    culoc_t(const loc_t *in=0):dx(0),dy(0){
	if(in){
	    dx=in->dx;
	    dy=in->dy;
	    cp2gpu(p, in);
	}
    }
};
class cupts_t:public culoc_t{
    Real dxsa;
    long nxsa;
public:
    Real Dxsa(){
	return dxsa;
    }
    long Nxsa(){
	return nxsa;
    }
    cupts_t(pts_t *in=0):culoc_t((loc_t*)in),dxsa(0),nxsa(0){
	if(!in) return;
	dxsa=in->dx;
	nxsa=in->nx;
    }
};
curmat gpu_dmcubic_cc(Real iac);
/**
   Specifies the grid.
*/
class cugrid_t{
public:
    long  nx, ny;
    Real ox, oy;
    Real dx, dy;
    Real ht;
    Real vx, vy;
    Real iac;
    curmat cubic_cc; /*coefficients for cubic influence function. */
    //use default copy assignment operator and copy constructor

    cugrid_t &operator=(const map_t *in){
	if(in){
	    nx=in->nx;
	    ny=in->ny;
	    ox=in->ox;
	    oy=in->oy;
	    dx=in->dx;
	    dy=in->dy;
	    ht=in->h;
	    vx=in->vx;
	    vy=in->vy;
	    if(fabs(iac-in->iac)>fabs(iac+in->iac)*1e-5){
		iac=in->iac;
		cubic_cc=gpu_dmcubic_cc(in->iac);
	    }
	}
	return *this;
    }
    cugrid_t():nx(0),ny(0),ox(0),oy(0),dx(0),dy(0),ht(0),vx(0),vy(0),iac(0){}
    cugrid_t Scale(Real sc)const{
	cugrid_t tmp(*this);
	tmp.ox*=sc;
	tmp.oy*=sc;
	tmp.dx*=sc;
	tmp.dy*=sc;
	return tmp;
    }
    operator bool(){
	return nx&&ny;
    }
};
class cumap_t:public cugrid_t{
public:
    curmat p;
    /*Init the data, p*/
    cumap_t(){}
    cumap_t(cugrid_t &grid):cugrid_t(grid){
	if(nx>0 && ny>0) {
	    p=curmat(nx,ny);
	}
    }
    cumap_t& operator=(map_t*in){
	if(in){
	    cugrid_t::operator=(in);
	    if(p.Nx()!=nx || p.Ny()!=ny){
		p=curmat(nx, ny);
	    }
	}
	return *this;
    }
    operator const curmat&()const{
	return p;
    }
    operator curmat&(){
	return p;
    }
    Real *operator()(){
	return p();
    }
    const Real *operator()()const{
	return p();
    }
};

typedef Array<cumap_t> cumapcell;
typedef Array<cugrid_t> cugridcell;


template <typename T>
void Zero(Array<T, Gpu> &A, cudaStream_t stream){
    Gpu<T>::zero(A(), A.N(), stream);
    if(sizeof(T)>8){
	warning("call dev zero for %s. This will lead to code error.\n", typeid(T).name());
    }
}

template <typename T>
void Zero(CuCell<T> &A, cudaStream_t stream){
    for(long i=0; i<A.N(); i++){
	Zero(A[i], stream);
    }
}
template <typename T>
//Create a similar cell from this cell.
CuCell<T> New(const CuCell<T> &A){
    if(!A.M()){
	return CuCell<T>(A.Nx(), A.Ny());
    }else{
	long mx[A.N()];
	long my[A.N()];
	for(long i=0; i<A.N(); i++){
	    mx[i]=A[i].Nx();
	    my[i]=A[i].Ny();
	}
	return CuCell<T>(A.Nx(), A.Ny(), mx, my);
    }
}
/*Transpose a matrix in naive way. Faster way is to use shared memory and handle
  a block each time.*/
template <typename T>
__global__ void transpose(T *restrict out, const T *restrict in, int nx, int ny){
    const int stepx=blockDim.x * gridDim.x;
    const int stepy=blockDim.y * gridDim.y;
    const int ix0=threadIdx.x+blockDim.x*blockIdx.x;
    const int iy0=threadIdx.y+blockDim.y*blockIdx.y;
    for(int iy=iy0; iy<ny; iy+=stepy){
	for(int ix=ix0; ix<nx; ix+=stepx){
	    out[iy+ix*ny]=in[ix+iy*nx];
	}
    }
}
template <typename T>
Array<T, Gpu> Transpose(Array<T, Gpu> &A, cudaStream_t stream){
    Array<T, Gpu> B(A.Ny(),A.Nx());;
    transpose<<<dim3(16,16),dim3(16,16),0,stream>>>
	(B(), A(), A.Nx(), A.Ny());
    return B;
}
#endif

