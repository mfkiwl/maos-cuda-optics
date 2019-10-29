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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "type.h"
#include "mathdef.h"
#include "defs.h"


/**
   Create a new X(mat) matrix object, mmapped from file. The file is truncated if already exists in rw mode.*/
X(mat)* X(new_mmap)(long nx, long ny, const char *header, const char *format, ...){
    if(!nx || !ny) return NULL;
    if(disable_save){
	return X(new)(nx, ny);
    }
    format2fn;
    size_t metasize=3*8+bytes_header(header);//size of meta data. 
    size_t msize=nx*ny*sizeof(T)+metasize;//total size of file/memory
    mem_t *mem=mmap_open(fn, msize, 1);
    char *map=(char*)mem_p(mem);//save value of map
    mmap_write_header(&map, M_T, nx, ny, header);
    X(mat) *out=X(new_do)(nx, ny, (T*)map, mem);
    memset(map, 0, nx*ny*sizeof(T));//Is this necessary?
    if(header) out->header=strdup(header);
    return out;
}
/**
   Create a new X(cell) matrix cell object, mmapped from file. The file is
   truncated if already exists. 
*/
X(cell)* X(cellnew_mmap)(long nx, long ny, long *nnx, long *nny, 
			 const char *header1, const char *header2[],
			 const char *format, ...){
    if(!nx || !ny) return NULL;
    if(disable_save){
	return X(cellnew3)(nx, ny, nnx, nny);
    }
    format2fn;
    long metasize=3*8;
    long msize=metasize+bytes_header(header1);
    for(long ix=0; ix<nx*ny; ix++){
	long mx=(long)nnx>0?nnx[ix]:(-(long)nnx);
	long my=nny?((long)nny>0?nny[ix]:(-(long)nny)):1;
	msize+=metasize+mx*my*sizeof(T)+bytes_header(header2?header2[ix]:NULL);
    }
    mem_t *mem=mmap_open(fn, msize, 1);
    char *map=mem_p(mem);
    mmap_write_header(&map, MCC_ANY, nx, ny, header1);
    X(cell) *out=X(cellnew)(nx,ny);
    if(header1) out->header=strdup(header1);

    for(long ix=0; ix<nx*ny; ix++){
	long mx=(long)nnx>0?nnx[ix]:(-(long)nnx);
	long my=nny?((long)nny>0?nny[ix]:(-(long)nny)):1;
	const char *header2i=header2?header2[ix]:NULL;
	mmap_write_header(&map, M_T, mx, my, header2i);
	out->p[ix]=X(new_do)(mx, my, (T*)map, mem);
	memset(map, 0, mx*my*sizeof(T));//Is this necessary?
	map+=nnx[ix]*my*sizeof(T);
	if(header2i) out->p[ix]->header=strdup(header2i);
    }
    
    return out;
}
/**
   Create a new X(cell) matrix cell object, with identical blocks, mmapped from
   file. be aware that the data is not 8-byte aligned. The file is truncated if
   already exists.  */
X(cell)* X(cellnewsame_mmap)(long nx, long ny, long mx, long my, const char *header,
			     const char *format, ...){
    format2fn;
    return X(cellnew_mmap)(nx, ny, (long*)-mx, (long*)-my, header, NULL, "%s", fn);
}

static X(mat*) X(readdata_mmap)(char **map, mem_t *mem){
    long nx, ny;
    uint32_t magic;
    const char *header;
    mmap_read_header(map, &magic, &nx, &ny, &header);
    if(magic!=M_T){
	error("File has magic %d, we want %d\n", (int)magic, M_T);
    }
    X(mat)* out=X(new_do)(nx, ny, (T*)(*map), mem);
    *map+=nx*ny*sizeof(T);
    if(out){
	out->header=strdup(header);
    }
    return out;
}
/**
   Map the file to memory in read only, shared mode.
*/
X(mat*) X(read_mmap)(const char *format, ...){
    format2fn;
    mem_t *mem=mmap_open(fn, 0, 0);
    char *map=(char*)mem_p(mem);
    X(mat) *out=X(readdata_mmap)(&map, mem);//mem_new(fd, map0,msize));
    return out;
}

/**
   Map the file to memory in read only, shared mode.
*/
X(cell*) X(cellread_mmap)(const char *format, ...){
    format2fn;
    mem_t *mem=mmap_open(fn, 0, 0);
    char *map=mem_p(mem);
    long nx, ny;
    uint32_t magic;
    const char *header;
    mmap_read_header(&map, &magic, &nx, &ny, &header);
    if(!iscell(&magic)){
	error("We want a cell array, File has %x\n", (int)magic);
    }
    X(cell) *out=X(cellnew)(nx, ny);
    if(header) out->header=strdup(header);
    for(long ix=0; ix<nx*ny; ix++){
	out->p[ix]=X(readdata_mmap)(&map, mem);
    }
    return out;
}
