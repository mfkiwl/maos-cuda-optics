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

#include "utils.h"
#include "tomo.h"
#include "recon.h"
#include "cucmat.h"
#define TIMING 0
namespace cuda_recon{
void cufdpcg_t::update(FDPCG_T *fdpcg){
    //copy or update Mb. 
    int nxsave=fdpcg->Mbinv->nx;
    fdpcg->Mbinv->nx=nb;
    cp2gpu(Mb, fdpcg->Mbinv);
    fdpcg->Mbinv->nx=nxsave;
}
cufdpcg_t::cufdpcg_t(FDPCG_T *fdpcg, const curecon_geom *_grid)
    :grid(_grid),fftnc(0),fftips(0),nb(0),bs(0),nby(0),nbz(0),scale(0){
    if(!fdpcg) return;
    scale=fdpcg->scale;
    bs=fdpcg->bs;//linear size of each block.
    nb=fdpcg->permhf->nx/bs;//size of non-redundant block
    update(fdpcg);
    cp2gpu(perm, fdpcg->permhf->p, nb*bs, 1);
    int nps=grid->npsr;
    int count=0;
    int osi=-1;
    int start[nps];
    for(int ips=0; ips<nps; ips++){
	/*group layers with the same os together in a batch fft.*/
	if(osi != grid->xnx[ips]){
	    osi = grid->xnx[ips];
	    start[count]=ips;
	    count++;
	}
    }
    fft=Array<cufftHandle>(count,1);
    ffti=Array<cufftHandle>(count,1);
    fftnc=count;
    fftips=Array<int>(count+1,1);
    for(int ic=0; ic<count; ic++){
	fftips[ic]=start[ic];
    }
    fftips[count]=nps;
    for(int ic=0; ic<count; ic++){
	int ncomp[2];
	/*Notice the reverse in specifying dimensions. THe first element is outmost rank.*/
	ncomp[0]=grid->xny[start[ic]];
	ncomp[1]=grid->xnx[start[ic]];

	int nembed[2];
	nembed[0]=grid->xnx[start[ic]]*grid->xny[start[ic]];
	nembed[1]=grid->xnx[start[ic]];
	DO(cufftPlanMany(&fft[ic], 2, ncomp, 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 FFT_T_R2C, fftips[ic+1]-fftips[ic]));
	DO(cufftPlanMany(&ffti[ic], 2, ncomp, 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 nembed, 1, ncomp[0]*ncomp[1],
			 FFT_T_C2R, fftips[ic+1]-fftips[ic]));
    }
    xhat1=cuccell(grid->npsr, 1, grid->xnx, grid->xny);
    xhat2=cuccell(grid->npsr, 1, grid->xnx, grid->xny);
    {
	/*nby is blockDim.y: number of FD blocks in each threading block*/
	/*nbz is gridDim.x: number of threading blocks to launch*/
	nby=256/bs;//number of blocks in each grid
	nbz=nb/nby;//number of grids to launch.
    }
    /* notice: performance may be improved by using
       R2C FFTs instead of C2C. Need to update perm
       and Mbinv to use R2C.*/
    GPU_FDPCG_T *FDDATA=new GPU_FDPCG_T[nps];
    for(int ips=0; ips<nps; ips++){
	FDDATA[ips].nx=grid->xnx[ips];
	FDDATA[ips].ny=grid->xny[ips];
	if(scale){
	    FDDATA[ips].scale=1.f/sqrtf((Real)(grid->xnx[ips]*grid->xny[ips]));
	}else{
	    FDDATA[ips].scale=1.f;
	}
    }
    fddata=Array<GPU_FDPCG_T,Gpu>(nps, 1);
    cudaMemcpy(fddata(), FDDATA, sizeof(GPU_FDPCG_T)*nps, cudaMemcpyHostToDevice);
    delete [] FDDATA;
}
/*
  2012-08-01: Have tried the following with no improvement:
  Unroll with template.
  Each 32 threads to handle each bs for bs<32.

  Changed 2013-08-13:
  Input and output may overlap for different blocks due to FFT based mirroring. So separate intput/output.
*/

__global__ static void fdpcg_mul_block_sync_half(Comp *xout, const Comp *xin, Comp *Mi, int *restrict perm, int nbb){
    extern __shared__ Comp v[];
    int bs=blockDim.x;//size of each block
    Comp *vin=v+threadIdx.y*2*bs;//stores reordered input
    Comp *vout=vin+bs;//stores output before reorder again
    int nstep=blockDim.y*gridDim.x;
    for(int ib=blockIdx.x*blockDim.y+threadIdx.y; ib<nbb; ib+=nstep){
	const Comp *M=Mi+ib*bs*bs;
	const int ix=threadIdx.x;
	const int pm=perm[ib*bs+ix];//last 2 bit are flags
	const int pm2=abs(pm);
	vin[ix]=xin[pm2];
	if(pm<0){//last bit is on: doing conjugate
	    vin[ix]=Z(cuConj)(vin[ix]);
	}
	vout[ix].x=vout[ix].y=0;
	__syncthreads();//wait for loading of vin
	for(int iy=0; iy<bs; iy++){
	    vout[ix]=Z(cuCfma)(M[ix+iy*bs], vin[iy], vout[ix]);
	}
	if(pm<0){
	    vout[ix]=Z(cuConj)(vout[ix]);
	}
	xout[pm2]=vout[ix];
    }
}
__device__ static inline void do_scale(Real &a, Real b){
    a*=b;
}
__device__ static inline void do_scale(Comp &a, Real b){
    a.x*=b;
    a.y*=b;
}
template<typename T> __global__ static void 
fdpcg_scale(GPU_FDPCG_T *fddata, T *const*xall){
    int ips=blockIdx.z;
    int nx=fddata[ips].nx*fddata[ips].ny;
    const int step=blockDim.x * gridDim.x; 
    Real scale=fddata[ips].scale;
    T *restrict x=xall[ips];
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nx; i+=step){
	do_scale(x[i], scale);
    }
}

/**
   The positive (right half) base frequency (os=1) couples to both positive and
   negative frequencies to layers with os=2. Only apply the block matrix to
   positive frequencies maynot be good.

   If all lays have the same oversampling, either os=1 or os=2, then postive
   frequencies only couple with postive frequencies, so does negative
   frequencies. We can actually skip the multiplication of negative frequencies.

*/
#define DBG_FD 0
void cufdpcg_t::Pre(curcell &xout, const curcell &xin, stream_t &stream){
#if TIMING
    EVENT_INIT(4)
#define RECORD(i) EVENT_TIC(i)
#else
#define RECORD(i)
#endif
    RECORD(0);
    if(!xin.M()){
	error("xin is not continuous");
    }
    if(!xout){
	xout=curcell(grid->npsr, 1, grid->xnx, grid->xny);
    }else if(!(xout).M()){
	error("xout is not continuous");
    }
#if DBG_FD
    cuwrite(xin, "fdg_xin");
#endif
    /*2012-07-11: Use real to complex fft*/
    for(int ic=0; ic<fftnc; ic++){
	int ips=fftips[ic];
	DO(cufftSetStream(fft[ic], stream));
	CUFFTR2C(fft[ic], xin[ips](), xhat1[ips]());
    }
    RECORD(1);
    if(scale){
	fdpcg_scale<<<dim3(9,1,grid->npsr), dim3(256,1),0,stream>>>
	    (fddata(), xhat1.pm());
    }
#if DBG_FD
	cuccellwrite(xhat1, "fdg_fft");
#endif

    fdpcg_mul_block_sync_half<<<nbz, dim3(bs,nby), sizeof(Comp)*bs*2*nby, stream>>>
	(xhat2.M()(), xhat1.M()(), Mb.M()(), perm, nb);
    RECORD(2);
#if DBG_FD
    cuccellwrite(xhat2, "fdg_mul");
#endif
    for(int ic=0; ic<fftnc; ic++){
	int ips=fftips[ic];
	DO(cufftSetStream(ffti[ic], stream));
	CUFFTC2R(ffti[ic], xhat2[ips](), xout[ips]());
    }
    if(scale){
	fdpcg_scale<<<dim3(9,1,grid->npsr),dim3(256,1),0,stream>>>
	    (fddata(), xout.pm());
    }
#if DBG_FD
    cuwrite(xout, "fdg_xout");
#endif
    RECORD(3);
#if TIMING
    EVENT_TOC;
    info("FDPCG: FFT %3.0f MUL %3.0f FFTI %3.0f Total %3.0f\n", 
	  times[1], times[2], times[3], times[0]);
#endif
}
}//namespace
