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
#include <typeinfo>
#include <cstddef>
#include <string>
#include "../sys/sys.h"
#include "numtype.h"
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

///Magic number denoting the data type. Specialized for each type.
template <typename T>
struct Magic{
    enum { magic=MCC_ANY};
};
template <>
struct Magic<double>{
    enum { magic=M_DBL};
};
template <>
struct Magic<float>{
    enum { magic=M_FLT};
};
template <>
struct Magic<dcomplex>{
    enum { magic=M_CMP};
};
template <>
struct Magic<fcomplex>{
    enum { magic=M_ZMP};
};
template <>
struct Magic<int>{//always 32 bit int
    enum { magic=M_INT32};
};
template <>
struct Magic<long int>{//always 64 bit int. long is not.
    enum { magic=M_INT64};
};

/**
   RefP is a reference counting for pointers. 

   It can be used to store dynamically allocated memory and use as a regular
   pointer. When passing around, it will do automatically reference accounting
   and release the memory when the last reference is released.

   It can store different kind of dynamic memory:
   1) heap memory
   2) GPU memory (cuda)
   3) pinned memory (paged locked, cuda)
   4) mmaped from file or shared memory.
*/

template <typename T, template<typename> class Dev=Cpu>
class RefP{
private:
    Dev<T> *p0; //The allocated memory address. p>=p0;
    int *nref;  //The reference counter. Delete p0 only if nref is valid and has value of 1.
    long n;
protected:
    T *p;  //The memory address for access
 
public:
    void deinit(){
	if(nref && !atomicadd(nref, -1)){
	    delete[] p0;
	    delete nref;
	}
	p0=0;
	nref=0;
	n=0;
	p=0;
    }

    //Constructors and related
    RefP():p0(0),nref(0),n(0),p(0){
    }
    //Interface to be deprecated /todo cast from T to Dev<T> is dangerous.
    RefP(long ni, T *pin=0, int own=1):p0((Dev<T>*)pin),nref(0),n(ni){
	if(n>0){
	    if(!p0){
		//Need to use new because T may be class types
		p0=new Dev<T>[n];
		own=1;
	    }
	    if(own){
		if(own>1){
		    error("This is not enticipated\n");
		}
		nref=new int;
		nref[0]=1;
	    }
	}
	p=(T*)p0;
    }
    //Create a new pointer with offset for p.
    RefP(const RefP& pin, long offset=0):p0(pin.p0),nref(pin.nref),n(pin.n),p(pin.p+offset){
	if(nref) atomicadd(nref, 1);
    }
 
    RefP &operator=(const RefP &in){
	if(this!=&in){
	    deinit();
	    p0=in.p0;
	    nref=in.nref;
	    n=in.n;
	    p=in.p;
	    if(nref) atomicadd(nref, 1);
	}
	return *this;
    }
    RefP &operator=(T*in){
	deinit();
	p0=0;
	nref=0;
	p=in;
	return *this;
    }
    void init(long ni=0){
	deinit();
	n=ni;
	if(n>0){
	    p0=new Dev<T>[n];
	    nref=new int;
	    nref[0]=1;
	}
	p=(T*)p0;	
    }
    //Destructors and related
    ~RefP(){
	deinit();
    }
    //Replace vector without reference counting.
    void Replace(T* pnew){
	deinit();
	p=pnew;
    }
    //Resize the memory. We use element-wise copy because T may be class time
    void Resize(long ni){
	if(nref && *nref==1){
	    long offset=p-(T*)p0;
	    Dev<T> *p0new=new Dev<T>[ni];
	    long nmin=MIN(n, ni);
	    for(long i=0; i<nmin; i++){
		p0new[i]=p0[i];
	    }
	    delete[] p0;
	    p0=p0new;
	    n=ni;
	    p=(T*)p0+offset;
	}else{
	    error("Canot Resize referenced or not owned data.\n");
	}
    }
    //Access operators
    //() operator
    T* operator()(){
	return (T*)p;
    }
    const T* operator()() const{
	return (T*)p;
    }
    //conversion operator. Automatically satisfying [] operator
    operator T*(){
	return (T*)p;
    }
    operator const T*()const{
	return (T*)p;
    }
    //indexing operator
    T& operator()(int i){
	assert(i>=0 && i<n);
	return p[i];
    }
    const T& operator()(int i)const{
	assert(i>=0 && i<n);
	return p[i];
    }
    //comparison
    bool operator==(const RefP&in){
	return p==in.p;
    }
    //ponter offset
    /*T*operator+(int i){
	return p+i;
    }
    const T* operator+(int i)const {
	return p+i;
	}*/
    long RefCount() const{
	return nref?nref[0]:0;
    }
};

/**
   Base class for arrays and cells. 
 */
class TwoDim{
public://temporary. for backward compatibility before conversion is done. /todo
    uint32_t id;
public://temporary. Change to protected after conversion is done. /todo
    long nx;
    long ny;
public:
    TwoDim(long nxi=0, long nyi=1):id(0),nx(nxi),ny(nyi){}
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
    //make the class polymorphic
    virtual ~TwoDim(){}
};
/**
   Generic array of basic types and classes.
*/
template <typename T, template<typename> class Dev=Cpu>
class Array:public TwoDim{
//, public RefP<T, Dev>{
    typedef RefP<T, Dev> RefPT;
public:
    RefPT p;
public:
    std::string desc;
public: //temporary for backward compatibility before conversion is done. /todo
    mmap_t mmap;/**< not NULL if mmaped.*/ 
    struct fft_t *fft;
    char *header; //temporary for backward compatibility. Convert to desc. /todo
public:
    //using RefPT::operator();
    //using RefPT::p;

    T*operator()(){
	return p();
    }
    const T*operator()()const{
	return p();
    }

    operator T*(){
	return p();
    }
    operator const T*()const {
	return p();
    }


    T&operator ()(long ix, long iy){
	assert(ix>=0 && ix<nx && iy>=0 && iy<ny);
	return p[ix+nx*iy];
    }
    const T&operator ()(long ix, long iy)const{
	assert(ix>=0 && ix<nx && iy>=0 && iy<ny);
	return p[ix+nx*iy];
    }
    bool operator==(const Array&in){
	return RefPT::operator==(in) && TwoDim::operator==(in);
    }
    T *Col(int icol){
	assert(icol>=0 && icol<ny);
	return p+nx*icol;
    }
    const T *Col(int icol)const{
	assert(icol>=0 && icol<ny);
	return p+nx*icol;
    }
 
    void init(long nxi=0, long nyi=1){
	nx=nxi;
	ny=nyi;
	p.init(nxi*nyi);
    }
    virtual void deinit(){
	p.deinit();
	nx=0;
	ny=0;
    }
    virtual ~Array(){
	if(!mmap) free(header);
	//dfft_free_plan(fft);//temporary;/todo
	/*if(std::is_pointer<T>::value){
	    info("T is pointer\n");
	    for(long i=0; i<nx*ny; i++){
		delete (*this)[i];
	    }
	    }*/
    }
    //Constructors 
    explicit Array(long nxi=0, long nyi=1, T *pi=NULL, long own=1)
	:TwoDim(nxi, nyi),p(nxi*nyi, pi, own){
	fft=0;
	header=0;
	id=Magic<T>::magic;
    }
    //Create a reference with offset.
    Array(long nxi,long nyi,const Array& pi,long offset=0)
	:TwoDim(nxi, nyi),p(pi.p,offset){
	fft=0;
	header=0;
	id=Magic<T>::magic;
    }
    //Copy constructor
    Array(const Array &in):TwoDim(in),p(in.p),desc(in.desc),mmap(in.mmap){
	fft=0;
	header=0;
	id=Magic<T>::magic;
    }
    //Use default destructor
    //Copy assignment operator
    Array &operator=(const Array &in){
	if(this!=&in){
	    p=in.p;
	    //RefPT::operator=(in);
	    nx=in.nx;
	    ny=in.ny;
	    desc=in.desc;
	    mmap=in.mmap;
	}
	return *this;
    }
    //Convert matrix into Vector
    Array Vector(){
	return Array(nx*ny,1, *this, 0);
	/*Array tmp=*this;
	tmp.nx=tmp.nx*tmp.ny;
	tmp.ny=1;
	return tmp;*/
    }
    //Resize array in place
    void Resize(long nxi, long nyi){
	Array B(nxi, nyi);
	for(long iy=0; iy<MIN(nyi, Ny()); iy++){
	    for(long ix=0; ix<MIN(nxi, Nx()); ix++){
		B(ix, iy)=(*this)(ix, iy);
	    }
	}
	*this=B;
    }
    //Replace underlining vector
    void Replace(T* pnew){
	p.Replace(pnew);
    }
};

template <typename T, template<class> class Dev>
struct Magic<Array<T, Dev> >{
    enum { magic=MCC_ANY};
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
