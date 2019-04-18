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

#ifndef AOS_MATH_TYPE_H
#define AOS_MATH_TYPE_H
#include "numtype.h"
#include "array.h"
/**
   \file type.h Defines the math data types like dmat, cmat, dcell, ccell,
   dsp, csp data types.

   Don't use ulong for dimensions because subtracting a bigger ulong from a
   smaller ulong overflows.  */

typedef enum CEMBED{
    C_FULL,
    C_ABS2,
    C_REAL,
    C_ABS,
    C_LITERAL
}CEMBED;
static inline uint32_t ID(const TwoDim *A){
    return A->id;
}
static inline int iscell(const TwoDim *A){
    return magic_iscell(ID(A));
}
/*
  We use pointers for reference counter because different array
  may use the same pointer, but with different nx or ny
  partition. */

template <typename T>
class Sparse:public TwoDim{
public:
    RefP<T> x;     /**<The value*/
    RefP<spint> p; /**< column pointers (size n+1) or col indices (size nzmax) when nz!=-1 */ 
    RefP<spint> i; /**< row indices, size nzmax */ 
    //T *restrict x;
    //spint *restrict p ; 
    //spint *restrict i ; 
    char *header;
    long nzmax;
    //int *nref;           /**< reference counting like dmat */ 
public:
    Sparse(const Sparse &in):TwoDim(in),x(in.x),p(in.p),i(in.i),header(0),nzmax(in.nzmax){
    }
    Sparse(long nxi=0, long nyi=1, long nzmaxi=0):TwoDim(nxi, nyi),x(0),p(0),i(0),header(0),nzmax(nzmaxi){
	if(nzmax){
	 p.init(ny+1);
	 i.init(nzmax);
	 x.init(nzmax);
	}
    }
    virtual void deinit(){
	/*if(nref && !atomicadd(nref, -1)){
	    free(p);
	    free(i);
	    free(x);
	    free(nref);
	}
	nref=0;*/
	if(header) free(header);
    }
    virtual ~Sparse(){deinit();}
    Sparse &operator=(const Sparse &in){
	if (this != &in){
	    deinit();
	    x=in.x;
	    p=in.p;
	    i=in.i;
	    nzmax=in.nzmax;
	    if(in.header) header=strdup(in.header);
	}
	return *this;
    }

};

#define DefSparse(T, C)					\
    class C: public Sparse<T>{				\
    public:						\
    C(long nxi=0, long nyi=1, long nzmaxi=0):Sparse<T>(nxi, nyi,nzmaxi){}; \
    }								
    
/*typedef MATARR(double) dmat;//a double matrix object contains 2-d array of double numbers
  typedef MATARR(float) smat;
  typedef MATARR(dcomplex) cmat;
  typedef MATARR(fcomplex) zmat;
  typedef MATARR(long) lmat;
*/

typedef Array<long>   lmat;
typedef Array<double> dmat;
typedef Array<float>  smat;
typedef Array<dcomplex> cmat;
typedef Array<fcomplex> zmat;

DefSparse(double, dsp);
DefSparse(float, ssp);
DefSparse(dcomplex, csp);
DefSparse(fcomplex, zsp);
#undef DefSparse

/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. Can be casted to dmat
*/
class map_t:public Array<double>{
    typedef Array<double> Parent;
    /*The OPD, takes the same form of dmat so can be casted. */
    //ARR(double);
public:
    double dx=0;      /**<Sampling along x*/
    double dy=0;      /**<Sampling along y*/
    double ox=0;      /**<Origin in x*/
    double oy=0;      /**<Origin in y*/
    double h=0;       /**<Heigh conjugation of this surface*/
    double vx=0;      /**Wind velocity. Useful for atmospheric grid*/
    double vy=0;      /**Wind velocity. Useful for atmospheric grid*/
    double iac=0;     /**<Inter-actuator coupling. >0: use cubic influence function*/
public:
    map_t(long nxi, long nyi, double dxi=1./64., double dyi=1./64., double *pi=0)
	:Parent(nxi, nyi, pi, 1), dx(dxi), dy(dyi),ox(-nx/2*dx),oy(-ny/2*dy),h(0),vx(0),vy(0),iac(0){
    }
    map_t(Parent in):Parent(in), dx(1./64), dy(1./64), ox(-nx/2*dx), oy(-ny/2*dy), h(0), vx(0), vy(0), iac(0){}
};

/**
   Map with different x/y sampling. Can be cased to dmat
*/
class rmap_t:public Array<double>{
    typedef Array<double> Parent;
public:
    double dx=0;      /**<Sampling along x (first dimension)*/
    double dy=0;      /**<Sampling along y (second dimension)*/
    double ox=0;      /**<Origin in x*/
    double oy=0;      /**<Origin in y*/
    double txdeg=0;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    double tydeg=0;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    double ftel=0;    /**<Effective focal length of the telescope*/
    double fexit=0;   /**<The distance between the exit pupil and the focus*/
    double fsurf=0;   /**<The distance between the tilted surface (M3) and the focus*/
public:
    rmap_t(long nxi, long nyi, double dxi=1./64., double dyi=1./64., double *pi=0)
	:Parent(nxi, nyi, pi, 1), dx(dxi), dy(dyi),ox(-nxi/2*dxi),oy(-nyi/2*dyi),
	 txdeg(0),tydeg(0),ftel(0),fexit(0),fsurf(0){
    }
};

/**
   Store starting x,y for each col
*/
typedef struct locstatcol_t{
    double xstart; /**<starting x of this column*/
    double ystart; /**<starting y of this column*/
    long   pos;    /**<starting index of this column*/
}locstatcol_t;

/**
   Stores array of locstatcol_t

*/
typedef struct locstat_t{
    locstatcol_t *cols; /**<Information about each column*/
    double dx;          /**<Sampling of the grid along x*/
    double dy;          /**<Sampling of the grid along y*/
    double xmin;        /**<Minimum x*/
    double ymin;        /**<Minimum y*/
    long   ncol;        /**<Number of consecutive columns found*/
    long   nx,ny;       /**<Size for embedding*/
}locstat_t;

/**
   Struct for coordinates like plocs, xloc, aloc etc.
*/
typedef struct loc_t:public TwoDim{
    union{
	double *locx=0;  /**< x coordinates of each point*/
	double *origx;
    };
    union{
	double *locy=0;  /**< y coordinates of each point*/
	double *origy;
    };
    union{
	long nloc=0;   /**< number of points*/
	long nsa;
    };
    union{
	double dx=0;     /**< Sampling along x*/
	double dsa;
	double dsax;
    };
    union{
	double dy=0;     /**< Sampling along y*/
	double dsay;
    };
    double ht=0;     /**< Conjugation height of the loc grid.*/
    double iac=0;    /**<Inter-actuator coupling. >0: use cubic influence function for ray tracing*/
    locstat_t *stat=0;/**<points to column statistics*/
    map_t *map=0;    /**< point to the map used for identifying neihboring points.*/
    int npad=0;      /*padding when create map*/
    int *nref=0;       /**<Reference counting*/
    loc_t(const loc_t&loc):TwoDim(0,0),locx(loc.locx),locy(loc.locy),nloc(loc.nloc),dx(loc.dx),dy(loc.dy),ht(loc.ht),
			   iac(loc.iac),stat(loc.stat),map(loc.map),npad(loc.npad), nref(loc.nref){
	if(nref) (*nref)++;
    }
    loc_t(){
	id=M_LOC64;
    }
}loc_t;
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be used as loc_t.
*/
typedef struct pts_t:public loc_t{
    int nx=0;        /**<number of cols per subaperture*/
    int ny=0;        /**<number of rows per subaperture*/
    double dx=0;     /**<sampling of points in each subaperture*/
    double dy=0;     /**<sampling of points in each subaperture. dy=dx normally required.*/
    pts_t(){
	id=M_LOC64;
    }
    pts_t(const loc_t& in):loc_t(in){}
}pts_t;

//temporary; To be converted to Cell /todo.
//Cannot use template here, otherwise functions using those needs explicit instantiation

//Do not implement destructor for cell
/*class cell:public Array<cell*>{				
public:						
    C(long nxi, long nyi):Array<T*>(nxi, nyi){}	
    };*/


class cell:public Array<TwoDim*>{				
public:						
    TwoDim* m;						
    cell(long nxi, long nyi):Array<TwoDim*>(nxi, nyi),m(0){}	
    virtual void deinit()override { 
	for(long ii=0; ii<nx*ny; ii++){			
	    delete p[ii];
	    p[ii]=0;
	}						
	if(m) delete m;					
    }							
    virtual ~cell(){deinit();}				
};



#define PCell(T, C)					\
    class C:public Array<T*>{				\
    public:						\
    T* m;						\
    C(long nxi, long nyi):Array<T*>(nxi, nyi),m(0){}	\
    virtual void deinit()override {			\
	for(long ii=0; ii<nx*ny; ii++){			\
	    delete p[ii];				\
	}						\
	if(m) delete m;					\
	Array<T*>::deinit();				\
    }							\
    virtual ~C(){deinit();}				\
    }

    //PCell(TwoDim, cell);
PCell(cmat, ccell);
PCell(zmat, zcell);
PCell(dmat, dcell);
PCell(smat, scell);
PCell(lmat, lcell);

PCell(dsp, dspcell);
PCell(ssp, sspcell);
PCell(csp, cspcell);
PCell(zsp, zspcell);

PCell(ccell, cccell);
PCell(zcell, zccell);
PCell(dcell, dccell);
PCell(scell, sccell);
PCell(lcell, iccell);

PCell(cccell, ccccell);
PCell(zccell, zcccell);
PCell(dccell, dcccell);
PCell(sccell, scccell);
PCell(iccell, icccell);

PCell(map_t, mapcell);
PCell(rmap_t, rmapcell);
PCell(loc_t, loccell);

PCell(mapcell, mapccell);
PCell(rmapcell, rmapccell);
PCell(loccell, locccell);


#undef ARR
#undef CELLARR
#undef MATARR

/*A method to simulate operator overloading for indexing arrys*/
#if DEBUG
static inline void assert_1d(long i, long nx, long ny){
    if(i<0 || i>=nx*ny){
	error("%ld is out of range for (%ld,%ld) array\n", i, nx, ny);
    }
}
static inline void assert_2d(long ix, long iy, long nx, long ny){
    if(ix<0 || ix>=nx || iy<0 || iy>=ny){
	error("(%ld,%ld) is out of range for (%ld,%ld) array\n", ix, iy, nx, ny);
    }
}
#define P1(A,i) ((A)->p[assert_1d((i), (A)->nx, (A)->ny),(i)])
#define P2(A,ix,iy) ((A)->p[assert_2d((ix), (iy), (A)->nx, (A)->ny),(ix)+(A)->nx*(iy)])
//#define PP1(A,i) ((A)->p+(assert_1d(i, (A)->nx, (A)->ny),(i)))
//#define PP2(A,ix,iy) ((A)->p+(assert_2d((ix), (iy), (A)->nx, (A)->ny),(ix)+(A)->nx*(iy)))
#else
#define P1(A,i) ((A)->p[(i)])
#define P2(A,ix,iy) ((A)->p[(ix)+(A)->nx*(iy)])
#endif
#define P0(A) ((A)->p)
#define PP0(A) ((A)->p)
#define PP1(A,i) ((A)->p+(i))
#define PP2(A,ix,iy) ((A)->p+(ix)+(A)->nx*(iy))
//#endif
//#define P0(A) _Pragma("#error Invalid use. Use P(A,i) or P(A,ix,iy)\n")
//#define PP0(A) _Pragma("#error Invalid use. Use PP(A,i) or PP(A,ix,iy)\n")
#define P_GET(_0,_1,_2,_3,NAME,...) NAME
#define P(...) P_GET(_0,__VA_ARGS__,P2,P1,P0)(__VA_ARGS__)
#define PP(...) P_GET(_0,__VA_ARGS__,PP2,PP1,PP0)(__VA_ARGS__)
#define PCOL(A,iy) ((A)->p+(iy)*(A)->nx)

//Define indexing using wrapping. See wrap()
//#define P1R(A,i) _Pragma("#error Invalid use. Use PR(A,i,j)")
#define PR(A,ix,iy) P2(A, wrap(ix, A->nx), wrap(iy, A->ny))
//#define PR(...) P_GET(_0,__VA_ARGS__,P2R,P1R,P1R,P1R)(__VA_ARGS__)
//#define PP1R(A,i) _Pragma("#error Invalid use. Use PPR(A,i,j)")
#define PPR(A,ix,iy) PP2(A, wrap(ix, A->nx), wrap(iy, A->ny))
//#define PPR(...) P_GET(_0,__VA_ARGS__,PP2R,PP1R,PP1R,PP1R)(__VA_ARGS__)
#define PCOLR(A,iy) ((A)->p+wrap(iy, A->ny)*(A)->nx)
#endif
