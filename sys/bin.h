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

#ifndef AOS_BIN_H_
#define AOS_BIN_H_
#include <sys/time.h>
#include <sys/mman.h>
#include <limits.h>
#include <stdint.h>
#include <pthread.h>
#include "common.h"
/**
   \file bin.h

   Defines our custom file format .bin that may be gzipped and the basic IO
   functions. All file read/write operators are through functions in this
   file. The routines can also operate on .fits files.
 */
/*The definitions here should not be changed once set for backward/foreward compatibility. */
#define M_CSP64  0x6400  /*sparse complex */
#define M_DSP64  0x6401  /*sparse double*/
#define M_DBL    0x6402  /*double */
#define M_INT64  0x6403  /*int 64 array */
#define M_CMP    0x6404  /*complex array */
#define M_INT32  0x6405  /*int 32 array */
#define M_CSP32  0x6406  /*sparse complex with 32 bit integer */
#define M_DSP32  0x6407  /*sparse double  with 32 bit integer*/
#define M_FLT    0x6408  /*single precision float. */
#define M_ZMP    0x6409  /*single precision complex */
#define M_INT8   0x640A  /*int 8  array */
#define M_INT16  0x640B  /*int 16 array */

#define M_SSP64  0x6430  /*single precision float + int64 */
#define M_SSP32  0x6431  /*single precision float + int32 */
#define M_ZSP64  0x6432  /*single precision complex + int64 */
#define M_ZSP32  0x6433  /*single precision complex + int32 */

/*The individual MC_* and MCC_* have been deprecated. Use MCC_ANY for cell arrays of any type for storage */
#define MCC_ANY  0x6421  /*cell of any thing */
#define M_HEADER 0x6500  /*header. */
#define M_SKIP   0x6600  /*the padding of magic number. */

#define M_MAP64     0x016402/*map_t, compatible with M_DBL*/
#define M_RECTMAP64 0x026402/*map_t, compatible with M_DBL*/
#define M_LOC64     0x036402/*loc_t with double data*/

static inline int magic_iscell(uint32_t magic){
    return (((magic)&0x6410)==0x6410 || ((magic)&0x6420) == 0x6420);
}
#if LONG_MAX==2147483647L //long is 32 bit
#define M_LONG M_INT32
#elif LONG_MAX==9223372036854775807L
#define M_LONG M_INT64 //long is 64 bit
#else
#error "Unknown long size"
#endif
#define M_INT M_INT32
#define M_STR M_INT8

#ifdef DLONG
#define M_SPINT M_INT64
#define M_DSP M_DSP64
#define M_CSP M_CSP64
#define M_SSP M_SSP64
#define M_ZSP M_ZSP64
#else
#define M_SPINT M_INT32
#define M_DSP M_DSP32
#define M_CSP M_CSP32
#define M_SSP M_SSP32
#define M_ZSP M_ZSP32
#endif

#define USE_ZLIB_H 0
#if USE_ZLIB_H
#include <zlib.h> /*zlib.h in ubuntu sucks */
#else
#ifdef __cplusplus
extern "C"{
#endif
typedef void* voidp;
voidp gzopen(const char *path, const char *mod);
voidp gzdopen(int fd, const char *mod);
long  gztell(voidp gzfile);
int   gzclose(voidp gzfile);
int   gzwrite(voidp gzfile, const void* buf, unsigned len);
long  gztell(voidp gzfile);
int   gzread(voidp gzfile, voidp buf, unsigned len);
int   gzseek(voidp file, long offset, int whence);
int   gzrewind(voidp file);
int   gzflush(voidp gzfile, int flush);
const char*gzerror(voidp gzfile, int *error);
#ifdef __cplusplus
}
#endif
#endif
typedef struct file_t file_t;
typedef struct mmap_t mmap_t;
typedef struct {
    uint32_t magic;//this must be the first element because we cast header_t to uint32_t.
    uint64_t nx;
    uint64_t ny;
    char *str;
}header_t;
/*
  Describes the information about mmaped data. Don't unmap a segment of mmaped
  memory, which causes the whole page to be unmapped. Instead, reference count
  the mmaped file and unmap the segment when the nref dropes to 1.
*/
struct mmap_t{
    int fd;   /**<file descriptor. close it if not -1.*/
    void *p;  /**<points to the beginning of mmaped memory for this type of data.*/
    long n;   /**<length of mmaped memory.*/
    long *nref;/**<Number of reference.*/

public:
    /*mmap_t(int fdi, void*pi, long ni):fd(fdi),p(pi),n(ni){
	nref=new long;
	*nref=1;
	}*/
    void _deinit(){
	if(nref){
	    (*nref)--;
	    if((*nref)<=0){
		munmap(p, n);
		if(fd!=-1) close(fd);
		delete nref;
	    }
	}
    }
    mmap_t():fd(0),p(0),n(0),nref(0){}
    void init(int fdi, void*pi, long ni){
	fd=fdi;
	p=pi;
	n=ni;
	nref=new long;
	*nref=1;
    }
    mmap_t(const mmap_t &in):fd(in.fd),p(in.p),n(in.n),nref(in.nref){
	if(nref) (*nref)++;
    }
    mmap_t& operator=(const mmap_t &in){
	if(this!=&in){
	    _deinit();
	    fd=in.fd;
	    p=in.p;
	    n=in.n;
	    nref=in.nref;
	    if(nref) (*nref)++;
	}
	return *this;
    }
    ~mmap_t(){
	_deinit();
    }
    operator bool(){
	return p?1:0;
    }
};

/*
  wrapping for standard and zlib functions that 
  open, close, write and read files, test end of file.
  Make the following function private so that other 
  routine won't call them directly to avoid conflict.
*/
extern int disable_save; ///saving to disk will be disabled when set to nonzero.
/*
  The following functions takes long type integers.
*/
int zfexist(const char *format, ...); 
void zftouch(const char *format, ...);
int zfeof(file_t *fp);
int zfpos(file_t *fp);
int zfseek(file_t *fp, long offset, int whence);
void zfrewind(file_t *fp);
file_t *zfopen_try(const char *fn, const char *mod);
file_t *zfopen(const char *fn, const char *mod);
const char *zfname(file_t *fp);
int zfisfits(file_t *fp);
void zfclose(file_t *fp);
void zflush(file_t *fp);
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp);
int zfread2(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
uint64_t bytes_header(const char *header);
void write_timestamp(file_t *fp);
void write_header(const header_t *header, file_t *fp);
int read_header2(header_t *header, file_t *fp);
void read_header(header_t *header, file_t *fp);
/**Check whether the header refers to a cell. If yes, return NULL. nx, ny are assigned to the dimension.*/
void writearr(const void *fpn, const int isfile, const size_t size, const uint32_t magic,
	      const char *header, const void *p, const uint64_t nx, const uint64_t ny);
void writedbl(const double *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void writeflt(const float *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
int mmap_open(char *fn, int rw);
void mmap_header_rw(char **p0, char **header0, uint32_t magic, long nx, long ny, const char *header);
void mmap_header_ro(char **p0, uint32_t *magic, long *nx, long *ny, char **header0);
#endif
