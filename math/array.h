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

#ifndef AOS_MATH_ARRAY_H
#define AOS_MATH_ARRAY_H

#include "../sys/sys.h"

//Standard CPU memory. 
template<typename T>
class Cpu{
    T val;
public:
    void *operator new[](size_t size){
	return ::calloc(size, 1);
    }
    void operator delete[](void*p){
	::free(p);
    }
    static void zero(T *p, size_t size){
	if(p) memset(p, 0, size*sizeof(T));
    }
    static void memcpy(T *pout, const T*pin, size_t size){
	::memcpy(pout, pin, size);
    }
};



/**
   RefP is a reference counting for pointers. 
*/
template <typename T, template<typename> class Dev>
class RefP{
private:
    Dev<T> *p0; //The allocated memory address. p>=p0;
    int *nref;  //The reference counter. Delete p0 only if nref is valid and has value of 1.
protected:
    T *p;  //The memory address for access
private:
    void _init(long n){
	if(n>0){
	    p0=new Dev<T>[n];
	    nref=new int;
	    nref[0]=1;
	}
	p=(T*)p0;
    }
    void _deinit(){
	if(nref && !atomicadd(nref, -1)){
	    delete[] p0;
	    delete nref;
	}
    }
public:
    //Constructors and related
    RefP():p(0),p0(0),nref(0){
    }
    RefP(long n, T *pin=0, int own=1):p0((Dev<T>*)pin),nref(0){
	if(n>0){
	    if(!p0){
		p0=new Dev<T>[n];
		own=1;
	    }
	    if(own){
		nref=new int;
		nref[0]=1;
	    }
	}
	p=(T*)p0;
    }
    //Create a new pointer with offset for p.
    RefP(const RefP& pin, long offset=0):p0(pin.p0),p(pin.p+offset),nref(pin.nref){
	if(nref) atomicadd(nref, 1);
    }
    
    RefP &operator=(const RefP &in){
	if(this!=&in){
	    p=in.p;
	    p0=in.p0;
	    nref=in.nref;
	    if(nref) atomicadd(nref, 1);
	}
	return *this;
    }
    void init(long n=0){
	_deinit();
	if(n>0){
	    p0=new Dev<T>[n];
	    nref=new int;
	    nref[0]=1;
	}else{
	    p0=0;
	    nref=0;
	}
	p=(T*)p0;	
    }
    //Destructors and related
    ~RefP(){
	_deinit();
    }
    
    //Access operators
    //() operator
    T* operator()(){
	return (T*)p;
    }
    const T* operator()() const{
	return (T*)p;
    }
    //conversion operator
    operator T*(){
	return (T*)p;
    }
    operator const T*()const{
	return (T*)p;
    }
    T& operator()(int i){
	return p[i];
    }
    const T& operator()(int i)const{
	return p[i];
    }
    bool operator==(const RefP&in){
	return p==in.p;
    }
    T*operator+(int i){
	return p+i;
    }
    const T* operator+(int i)const {
	return p+i;
    }
};
class TwoDim{
protected:
    long nx;
    long ny;
public:
    TwoDim(long nxi, long nyi):nx(nxi), ny(nyi){}
    long Nx()const{
	return nx;
    }
    long Ny()const{
	return ny;
    }
    long N()const{
	return nx*ny;
    }
    operator bool()const{
	return (nx && ny);
    }
    bool operator==(const TwoDim &in){
	return nx==in.nx && ny==in.ny;
    }
};
/**
   Generic array of basic types and classes.
*/
template <typename T, template<typename> class Dev=Cpu>
class Array:public TwoDim, public RefP<T, Dev>{
    typedef RefP<T, Dev> RefPT;
public:
    string header;
    using RefPT::operator();
    using RefPT::p;
    
    T&operator ()(int ix, int iy){
	return p[ix+nx*iy];
    }
    const T&operator ()(int ix, int iy)const{
	return p[ix+nx*iy];
    }
    bool operator==(const Array&in){
	return RefPT::operator==(in) && TwoDim::operator==(in);
    }
    T *Col(int icol){
	return p+nx*icol;
    }
    const T *Col(int icol)const{
	return p+nx*icol;
    }
 
    void init(long nxi=0, long nyi=1){
	nx=nxi;
	ny=nyi;
	RefPT::init(nxi*nyi);
    }
    //Constructors 
    explicit Array(long nxi=0, long nyi=1, T *pi=NULL, int own=1)
	:TwoDim(nxi, nyi),RefPT(nxi*nyi, pi, own){
    }
    //Create a reference with offset.
    Array(long nxi,long nyi,const RefPT& pi,long offset=0)
	:TwoDim(nxi, nyi),RefPT(pi,offset){
    }
    //Use default destructor

    Array(const Array &in):TwoDim(in),RefPT(in),header(in.header){
    }
    Array &operator=(const Array &in){
	if(this!=&in){
	    RefPT::operator=(in);
	    nx=in.nx;
	    ny=in.ny;
	    header=in.header;
	}
	return *this;
    }
 
    Array Vector(){
	Array tmp=*this;
	tmp.nx=tmp.nx*tmp.ny;
	tmp.ny=1;
	return tmp;
    }
};


/**
   Cell is a special Array that can stores multiple Arrays of data in continuous
   region.
*/

template <typename T, template<typename> class Dev=Cpu >
class Cell:public Array<Array<T,Dev>,Cpu>{
private:
    typedef Array<T,Dev> TMat;
    typedef Array<Array<T,Dev>,Cpu> Parent;
    TMat m; /*optional: contains the continuous data*/
public:
    using Parent::operator();
    using Parent::p;
    using Parent::nx;
    using Parent::ny;
    TMat &M(){
	return m;
    }
    const TMat&M() const{
	return m;
    }
   
    explicit Cell(long nxi=0, long nyi=1):Parent(nxi,nyi){
    }
    Cell(long nxi, long nyi, long mx, long my, T *pin=NULL):Parent(nxi,nyi){
	if(mx && my){
	    m=TMat(mx*my*nxi*nyi,1,pin, 0);
	    for(int i=0; i<nxi*nyi; i++){
		p[i]=TMat(mx, my, m, i*(mx*my));
	    }
	}
    }
    template <typename L>
    Cell(const long nxi, const long nyi, L *mx, L *my, T *pin=NULL):Parent(nxi,nyi){
	long tot=0;
	for(long i=0; i<nxi*nyi; i++){
	    tot+=mx[i]*(my?my[i]:1);
	}
	m=TMat(tot,1,pin,0);
	tot=0;
	for(long i=0; i<nxi*nyi; i++){
	    if(mx[i]){
		p[i]=TMat(mx[i],(my?my[i]:1),m, tot);
		tot+=p[i].N();
	    }
	}
    }
    //Use default copy constructor and copy assignment operator
    void init(long nxi=0, long nyi=1){
	Parent::init(nxi*nyi);
    }
    //Replace memory by a new one. We don't own the memory
    void Replace(T *pnew){
	if(m){
	    m=TMat(m.Nx(), m.Ny(), pnew, 0);//replace m
	}

	for(long i=0; i<nx*ny; i++){
	    p[i]=TMat(p[i].Nx(), p[i].Ny(), pnew, 0);//don't own the data pnew
	    pnew+=p[i].N();
	}
    }
};

/**
   Zeroing the array of fundemental type
*/
template <typename T, template<typename> class Dev >
void Zero(Array<T, Dev> &A){
    Dev<T>::zero(A(), A.N());
    if(sizeof(T)>8){
	warning("call dev zero for %s. This will lead to code error.\n", typeid(T).name());
    }
}
/**
   Zeroing the array of derived type
*/
template <typename T, template<typename> class Dev >
void Zero(Cell<T, Dev> &A){
    for(long i=0; i<A.N(); i++){
	Zero(A[i]);
    }
}
template <typename T, template<typename> class Dev >
//Create a similar cell from this cell.
Cell<T, Dev>New(const Cell<T,Dev> &A){
    if(!A.M()){
	return Cell<T, Dev>(A.Nx(), A.Ny());
    }else{
	long mx[A.N()];
	long my[A.N()];
	for(long i=0; i<A.N(); i++){
	    mx[i]=A[i].Nx();
	    my[i]=A[i].Ny();
	}
	return Cell<T, Dev>(A.Nx(), A.Ny(), mx, my);
    }
}
#endif
