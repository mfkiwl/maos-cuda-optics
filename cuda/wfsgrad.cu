/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#define TIMING 0

#include "utils.h"
#include "accphi.h"
#include <curand_kernel.h>
#include "wfs.h"
#include "cudata.h"
#include "cucmat.h"
#if !USE_CPP
extern "C"
{
#endif
#include "../maos/sim.h"
#include "../maos/utils.h"
#include "../maos/pywfs.h"
#if !USE_CPP
}
#endif

extern const char *dirskysim;
/*
  Notice that both blocks and threads are partitioning isa
*/
static __global__ void
add_geom_noise_do(Real *restrict g, const Real *restrict nea, 
					 int nsa, curandState *restrict rstat){
    const int id=threadIdx.x + blockIdx.x * blockDim.x;
    curandState lstat=rstat[id];
    const int nstep=blockDim.x * gridDim.x;
    for(int i=id; i<nsa; i+=nstep){
	Real n1=curand_normal(&lstat);
	Real n2=curand_normal(&lstat);
	g[i]+=n1*nea[i];
	g[i+nsa]+=n2*nea[i+nsa]+n1*nea[i+nsa*2];/*cross term. */
    }
    rstat[id]=lstat;
}

/**
   Compute ztilt over a square subaperture.
*/
static __global__ void 
cuztilt_do(Real *restrict g, Real *restrict opd, 
	   const int nsa, const Real dx, const int nx, Real *imcc,
	   const Real (*orig)[2], const Real*restrict amp, Real alpha){
    extern __shared__ Real a0[];
    int idx=threadIdx.x+threadIdx.y*blockDim.x;
    Real *a[3];
    for(int i=0; i<3; i++){
	a[i]=a0+blockDim.x*blockDim.y*i;
	a[i][idx]=0;
    }
    const int isa=blockIdx.x;
    const int skip=isa*nx*nx;
    const Real ox=orig[isa][0];
    const Real oy=orig[isa][1];
    for(int iy=threadIdx.y; iy<nx; iy+=blockDim.y){
	const int skip2=skip+iy*nx;
	const Real y=iy*dx+oy;
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    const int ind=skip2+ix;
	    const Real tmp=amp[ind]*opd[ind];
	    a[0][idx]+=tmp;
	    a[1][idx]+=tmp*(dx*ix+ox);
	    a[2][idx]+=tmp*y;
	}
    }
    for(int step=(blockDim.x*blockDim.y)>>1;step>0;step>>=1){
	__syncthreads();
	if(idx<step){
	    for(int i=0; i<3; i++){
		a[i][idx]+=a[i][idx+step];
	    }
	}
    }
    __syncthreads();
    if(threadIdx.x<2 && threadIdx.y==0){
	Real (*restrict A)[3]=(Real(*)[3])(imcc+isa*9);
	atomicAdd(&g[isa+threadIdx.x*nsa], 
		  alpha*(+a[0][0]*A[0][threadIdx.x+1]
			 +a[1][0]*A[1][threadIdx.x+1]
			 +a[2][0]*A[2][threadIdx.x+1]));
	
    }
}
void cuztilt(Real *restrict g, Real *restrict opd, 
	     const int nsa, const Real dx, const int nx, Real *imcc,
	     const Real (*orig)[2], const Real*restrict amp, Real alpha, cudaStream_t stream){
    const int tx=32;
    cuztilt_do<<<nsa, dim3(tx,tx), tx*tx*3*sizeof(Real), stream>>>
	(g, opd, nsa, dx, nx, imcc, orig, amp, alpha);
}
/**
   Apply matched filter. \todo this implementation relies on shared variable. It
   is probably causing competition.  */
__global__ static void
mtche_do(Real *restrict grad, Real (*mtches)[2], const Real *restrict ints,
	 int sigmatch, const Real *restrict i0sum, Real scale, int pixpsa, int nsa){
    extern __shared__ Real g0[];
    Real *g[3];
    for(int i=0; i<3; i++){
	g[i]=g0+blockDim.x*i;
	g[i][threadIdx.x]=0;
    }
    int isa=blockIdx.x;
    ints+=isa*pixpsa;
    const Real (*const restrict mtche)[2]=mtches+pixpsa*isa;
 
    for (int ipix=threadIdx.x; ipix<pixpsa; ipix+=blockDim.x){
	g[0][threadIdx.x]+=mtche[ipix][0]*ints[ipix];
	g[1][threadIdx.x]+=mtche[ipix][1]*ints[ipix];
	g[2][threadIdx.x]+=ints[ipix];
    }
    for(int step=(blockDim.x)>>1;step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    for(int i=0; i<3; i++){
		g[i][threadIdx.x]+=g[i][threadIdx.x+step];
	    }
	}
    }
    __syncthreads();
    if(threadIdx.x<2){
	if(sigmatch==1){/*normalize gradients according to siglev.*/
	    g[threadIdx.x][0]*=i0sum[isa]/g[2][0];
	}else if(sigmatch==2){
	    g[threadIdx.x][0]*=scale;
	}
	grad[isa+nsa*threadIdx.x]=g[threadIdx.x][0];
    }
}

static void
mtche(Real *restrict grad, Real (*mtches)[2],
      const Real *restrict ints, int sigmatch, const Real *restrict i0sum, Real scale,
      int pixpsa, int nsa, int msa, cudaStream_t stream){
    for(int isa=0; isa<nsa; isa+=msa){
	int ksa=MIN(msa, nsa-isa);
	mtche_do<<<ksa, 16, 16*3*sizeof(Real), stream>>>
	    (grad+isa, mtches+pixpsa*isa, ints+pixpsa*isa, sigmatch, i0sum?i0sum+isa:0, scale, pixpsa, nsa);
    }
}
/**
   Apply tCoG. /todo: replace atomicAdd by reduction.
*/
__global__ static void
tcog_do(Real *grad, const Real *restrict ints, Real siglev, Real *saa,
	int nx, int ny, Real pixthetax, Real pixthetay, int nsa, Real (*cogcoeff)[2], Real *srot){
    __shared__ Real sum[3];
    if(threadIdx.x<3 && threadIdx.y==0) sum[threadIdx.x]=0.f;
    __syncthreads();//is this necessary?
    int isa=blockIdx.x;
    ints+=isa*nx*ny;
    Real thres=cogcoeff[isa][0];
    Real bkgrnd=cogcoeff[isa][1];
    for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    Real im=ints[ix+iy*nx]-bkgrnd;
	    if(im>thres){
		atomicAdd(&sum[0], im);
		atomicAdd(&sum[1], im*ix);
		atomicAdd(&sum[2], im*iy);
	    }
	}
    }
    __syncthreads();
    if(threadIdx.x==0 && threadIdx.y==0){
	if(saa){
	    Real sum2=siglev*saa[isa];
	    if(sum[0]<sum2) sum[0]=sum2;
	}
	if(sum[0]>thres){
	    Real gx=(sum[1]/sum[0]-(nx-1)*0.5)*pixthetax;
	    Real gy=(sum[2]/sum[0]-(ny-1)*0.5)*pixthetay;
	    if(srot){
		Real s,c;
		Z(sincos)(srot[isa], &s, &c);
		Real tmp=gx*c-gy*s;
		gy=gx*s+gy*c;
		gx=tmp;
	    }
	    grad[isa]=gx;
	    grad[isa+nsa]=gy;
	}else{
	    grad[isa]=0;
	    grad[isa+nsa]=0;
	}
    }
}
/**
   Poisson random generator.
*/
__device__ static Real curandp(curandState *rstat, Real xm){
    Real g, t, xmu;
    int x=0, xu;
    if(xm>200){
	x=(int)round(xm+curand_normal(rstat)*sqrt(xm));
    }else{
	while(xm>0){
	    xmu = xm > 12.f ? 12.f : xm;
	    xm -= xmu;
	    g   = __expf(-xmu);
	    xu  = -1;
	    t   = 1.f;
	    while(t>g){
		xu++;
		t *= curand_uniform(rstat);
	    }
	    x += xu;
	}
    }
    return x;
}
/**
   Add noise to pix images.
*/
__global__ static void addnoise_do(Real *restrict ints0, int nsa, int pixpsa, Real bkgrnd, Real bkgrndc, 
				   const Real* restrict bkgrnd2s,
				   const Real* restrict bkgrnd2cs,
				   const Real* restrict qe, Real rne, curandState *rstat){
    const int id=threadIdx.x + blockIdx.x * blockDim.x;
    const int nstep=blockDim.x * gridDim.x;
    curandState lstat=rstat[id];
    for(int isa=id; isa<nsa; isa+=nstep){
	Real *restrict ints=ints0+isa*pixpsa;
	const Real *restrict bkgrnd2=bkgrnd2s?(bkgrnd2s+isa*pixpsa):NULL;
	const Real *restrict bkgrnd2c=bkgrnd2cs?(bkgrnd2cs+isa*pixpsa):NULL;
	for(int ipix=0; ipix<pixpsa; ipix++){
	    Real tot=(ints[ipix]+bkgrnd+(bkgrnd2?bkgrnd2[ipix]:0));
	    Real corr=bkgrnd2c?(bkgrnd2c[ipix]+bkgrndc):bkgrndc;
	    if(qe){
		ints[ipix]=(curandp(&lstat, tot*qe[ipix])+rne*curand_normal(&lstat))/qe[ipix]-corr;
	    }else{
		ints[ipix]=curandp(&lstat, tot)+rne*curand_normal(&lstat)-corr;
	    }
	}
    }
    rstat[id]=lstat;
}
static void
cu_fieldstop(curmat &opd, curmat &amp, int *embed, int nembed, 
	  curmat &fieldstop, Real wvl, cufftHandle fftplan, cudaStream_t stream){
    cucmat wvf(nembed, nembed);
    embed_wvf_do<<<DIM(opd.Nx(), 256), 0, stream>>> (wvf, opd, amp, embed, opd.Nx(), wvl);
    CUFFT(fftplan, wvf, CUFFT_FORWARD);
    cwm_do<<<DIM(wvf.N(), 256),0,stream>>> (wvf(), fieldstop(), wvf.N());
    CUFFT(fftplan, wvf, CUFFT_INVERSE);
    unwrap_phase_do<<<DIM2(wvf.Nx(), wvf.Ny(), 16),0,stream>>> (wvf, opd, embed, opd.Nx(), wvl);
}

__global__ static void
dither_acc_do(Real *restrict *imb, Real *restrict *imx, Real *restrict *imy, 
	      Real *restrict const *pints, Real cd, Real sd, int pixpsa, int nsa){
    for(int isa=blockIdx.x; isa<nsa; isa+=gridDim.x){
	const Real *ints=pints[isa];
	Real *restrict acc_ints=imb[isa];
	Real *restrict acc_intsx=imx[isa];
	Real *restrict acc_intsy=imy[isa];
	for(int ipix=threadIdx.x; ipix<pixpsa; ipix+=blockDim.x){
	    acc_ints[ipix]+=ints[ipix];
	    acc_intsx[ipix]+=ints[ipix]*cd;
	    acc_intsy[ipix]+=ints[ipix]*sd;
	}
    }
}
dither_t::dither_t(int nsa, int pixpsax, int pixpsay):imc(0){
    imb=curcell(nsa,1,pixpsax,pixpsay);
    imx=curcell(nsa,1,pixpsax,pixpsay);
    imy=curcell(nsa,1,pixpsax,pixpsay);
}

/**Accumulate for matched filter updating*/
void dither_t::acc(DITHER_T *dither, curcell &ints, Real cs, Real ss, int npll, cudaStream_t stream){
    const int nsa=ints.N();
    const int pixpsa=ints[0].N();
    dither_acc_do<<<nsa, pixpsa, 0, stream>>>
	(imb.pm, imx.pm, imy.pm, ints.pm, cs, ss, pixpsa, nsa);
    imc++;
    if(imc%npll==0){
	cp2cpu(&dither->imb, imb, stream);
	cp2cpu(&dither->imx, imx, stream);
	cp2cpu(&dither->imy, imy, stream);
	cuzero(imb);
	cuzero(imx);
	cuzero(imy);
    }
}
/**
   Calculate SHWFS gradients
*/
static void shwfs_grad(curmat&gradcalc, const curcell&ints, Array<cuwfs_t>&cuwfs, Array<cupowfs_t>&cupowfs, 
		       const PARMS_T *parms, const POWFS_T *powfs, SIM_T *simu, int iwfs, int ipowfs, stream_t& stream){
    const int nsa=powfs[ipowfs].saloc->nloc;
    CUDA_CHECK_ERROR;
    cuzero(gradcalc, stream);
    const int totpix=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
    switch(parms->powfs[ipowfs].phytype_sim){
    case 0:
	break; //no-op
    case 1://Matched filter
	{
	    Real sigratio=parms->powfs[ipowfs].sigmatch==2?(cuwfs[iwfs].i0sumsum/cursum(ints.M(),stream)):0;
	    mtche(gradcalc, (Real(*)[2])(cuwfs[iwfs].mtche()), ints.M(), 
		  parms->powfs[ipowfs].sigmatch, cuwfs[iwfs].i0sum(), sigratio,
		  totpix, nsa, cuwfs[iwfs].msa, stream);
	}
	break;
    case 2:{//CoG
	Real pixthetax=(Real)parms->powfs[ipowfs].radpixtheta;
	Real pixthetay=(Real)parms->powfs[ipowfs].pixtheta;
	int pixpsax=powfs[ipowfs].pixpsax;
	int pixpsay=powfs[ipowfs].pixpsay;
	Real scale1=0;
	Real *scale2=0;
	switch(parms->powfs[ipowfs].sigmatch){
	case 0://No signal level match. Use sum(i0) as denominator. Linear.
	    scale1=1.f;
	    scale2=cuwfs[iwfs].i0sum();
	    break;
	case 1://Use instantaneous intensity of each sa
	    break;
	case 2://Use averaged instantaneous intensity.
	    scale1=cursum(ints.M(),stream)/powfs[ipowfs].saasum;
	    scale2=cupowfs[ipowfs].saa();
	    break;
	default:
	    error("Invalid sigmatch\n");
	}
	Real *srot=parms->powfs[ipowfs].radpix?cuwfs[iwfs].srot():NULL;
	tcog_do<<<nsa, dim3(pixpsax, pixpsay),0,stream>>>
	    (gradcalc, ints[0], scale1, scale2,
	     pixpsax, pixpsay, pixthetax, pixthetay, nsa, 
	     (Real(*)[2])cuwfs[iwfs].cogcoeff(), srot);
    }
	break;
    default://Use CPU version.
	cp2cpu(&simu->ints->p[iwfs], ints, stream);
	CUDA_SYNC_STREAM;
	shwfs_grad(&simu->gradcl->p[iwfs],simu->ints->p[iwfs]->p,
		   parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
    }
    CUDA_CHECK_ERROR;
	    
}
/**
   Ray tracing and gradient computation for WFS. \todo Expand to do gradients in GPU without transfering
   data back to CPU.
*/
void gpu_wfsgrad_queue(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	//info("thread %ld gpu %d iwfs %d start\n", thread_id(), cudata->igpu, iwfs);
	gpu_set(cuglobal->wfsgpu[iwfs]);
	Array<cupowfs_t> &cupowfs=cudata->powfs;
	Array<cuwfs_t> &cuwfs=cuglobal->wfs;
	const PARMS_T *parms=simu->parms;
	const POWFS_T *powfs=simu->powfs;
	const RECON_T *recon=simu->recon;
	/*output */
	const int CL=parms->sim.closeloop;
	const int isim=simu->wfsisim;
	/*The following are truly constants for this powfs */
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(simu->wfsisim<parms->powfs[ipowfs].step) return;
	const int imoao=parms->powfs[ipowfs].moao;
	const int nsa=powfs[ipowfs].saloc->nloc;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const Real hs=parms->wfs[iwfs].hs;
	const Real hc=parms->wfs[iwfs].hc;
	const int dtrat=parms->powfs[ipowfs].dtrat;
	const int save_gradgeom=parms->save.gradgeom->p[iwfs];
	const int save_opd =parms->save.wfsopd->p[iwfs];
	const int save_ints=parms->save.ints->p[iwfs];
	const int noisy=parms->powfs[ipowfs].noisy;
	/*The following depends on isim */
	const int dtrat_output=((isim+1)%dtrat==0);
	const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
	const int do_pistatout=parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart;
	const int do_geom=!do_phy || save_gradgeom || do_pistatout;
	const Real thetax=parms->wfs[iwfs].thetax;
	const Real thetay=parms->wfs[iwfs].thetay;
	Real2 *loc=cupowfs[ipowfs].loc();
	/*Out to host for now. \todo : keep grad in device when do reconstruction on device. */
	stream_t &stream=cuwfs[iwfs].stream;
	curmat phiout=cuwfs[iwfs].phiout;
	curmat gradacc=cuwfs[iwfs].gradacc;
	curmat gradcalc=cuwfs[iwfs].gradcalc;
	curmat gradref=0;
	ctoc_init(30);
	if(isim%dtrat==0){
	    cuzero(cuwfs[iwfs].ints, stream);
	    cuzero(cuwfs[iwfs].gradacc, stream);
	}
	if(cuwfs[iwfs].opdadd){ /*copy to phiout. */
	    curcp(phiout, cuwfs[iwfs].opdadd, stream);
	}else{
	    cuzero(phiout, stream);
	}
	if(simu->atm && ((!parms->sim.idealwfs && !parms->powfs[ipowfs].lo)
			 || (!parms->sim.wfsalias && parms->powfs[ipowfs].lo))){
	    atm2loc(phiout, cuwfs[iwfs].loc_tel, hs, hc, thetax, thetay, 0, 0, parms->sim.dt, isim, 1, stream);
	}
	if(!parms->powfs[ipowfs].lo && (parms->sim.idealwfs || parms->sim.wfsalias)){
	    Real alpha=parms->sim.idealwfs?1:-1;
	    if(parms->sim.idealwfs==2 || parms->sim.wfsalias==2){
		dmat *opd=dnew(phiout.Nx(), 1);
		add2cpu(&opd, 0, phiout, 1, stream, 0);
		wfs_ideal_atm(simu, opd, iwfs, alpha);
		cp2gpu(phiout, opd);
		dfree(opd);
	    }else{
		dm2loc(phiout, cuwfs[iwfs].loc_dm, cudata->dmproj, parms->ndm,
		       hs, hc, thetax, thetay, 0, 0, alpha, stream);
	    }
	}
	if(simu->telws){
	    Real tt=simu->telws->p[isim];
	    Real angle=simu->winddir?simu->winddir->p[0]:0;
	    curaddptt(phiout, loc, 0, tt*cosf(angle), tt*sinf(angle), stream);
	}
	if(save_opd){
	    zfarr_push(simu->save->wfsopdol[iwfs], simu->wfsisim, phiout, stream);
	}
	if(CL){
	    wait_dmreal(simu, simu->wfsisim);
	    dm2loc(phiout, cuwfs[iwfs].loc_dm, cudata->dmreal, parms->ndm,
		   hs, hc, thetax, thetay, 0, 0, -1, stream);
	    Real ttx=0, tty=0;
	    if(simu->ttmreal){
		ttx+=simu->ttmreal->p[0];
		tty+=simu->ttmreal->p[1];
	    }
	    if(simu->fsmreal && simu->fsmreal->p[iwfs] && !powfs[ipowfs].llt){
		ttx+=simu->fsmreal->p[iwfs]->p[0];
		tty+=simu->fsmreal->p[iwfs]->p[1];
	    }
	    if(ttx || tty){
		curaddptt(phiout, loc, 0, -ttx, -tty, stream);
	    }
	}

	if(parms->tomo.ahst_idealngs==1 && parms->powfs[ipowfs].lo){
	    const real *cleNGSm=simu->cleNGSm->p+isim*recon->ngsmod->nmod;
	    ngsmod2loc(phiout, cupowfs[ipowfs].loc(), recon->ngsmod, cleNGSm, 
		       parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, 
		       -1, stream);
	}
	if(imoao>-1){
	    dm2loc(phiout, cuwfs[iwfs].loc_dm, cudata->dm_wfs[iwfs], 1,
		   INFINITY, hc, 0, 0, 0, 0, -1, stream);
	}
  
	Real focus=(Real)wfsfocusadj(simu, iwfs);
	if(Z(fabs)(focus)>1e-20){
	    const int nloc=cupowfs[ipowfs].loc.Nloc();
	    add_focus_do<<<DIM(nloc, 256), 0, stream>>>(phiout, loc, nloc, focus);
	}

	if(cupowfs[ipowfs].fieldstop){
	    if(parms->powfs[ipowfs].nwvl>1){
		error("Implement broadband case\n");
	    }
	    cu_fieldstop(phiout, cuwfs[iwfs].amp, cupowfs[ipowfs].embed[0], cupowfs[ipowfs].nembed[0], 
			 cupowfs[ipowfs].fieldstop[0], parms->powfs[ipowfs].wvl->p[0], cuwfs[iwfs].plan_fs, stream);
	}
	if(save_opd){
	    zfarr_push(simu->save->wfsopd[iwfs], simu->wfsisim, phiout, stream);
	}
	if(parms->plot.run){
	    const real *realamp=powfs[ipowfs].realamp->p[wfsind]->p;
	    drawopdamp_gpu("wfsopd",powfs[ipowfs].loc,phiout,stream,realamp,NULL,
			   "WFS OPD","x (m)", "y (m)", "WFS %d", iwfs);
	}
	ctoc("opd");
	if(parms->powfs[ipowfs].type==1){
	    CUDA_CHECK_ERROR;
	    pywfs_ints(cuwfs[iwfs].ints[0], phiout, cuwfs[iwfs],parms->wfs[iwfs].sigsim);
	    ctoc("pywfs");
	    CUDA_CHECK_ERROR;
	}else{
	    if(do_geom){
		real ratio;
		if(do_pistatout && dtrat>1){
		    gradref=gradcalc; 
		    cuzero(gradcalc, stream);
		    ratio=1;
		}else{
		    gradref=gradacc;
		    ratio=1.f/(Real)dtrat;
		}

		if(parms->powfs[ipowfs].gtype_sim==1){
		    cuztilt(gradref, phiout, 
			    cupowfs[ipowfs].pts.Nloc(), 
			    cupowfs[ipowfs].pts.Dxsa(), 
			    cupowfs[ipowfs].pts.Nxsa(), cuwfs[iwfs].imcc,
			    cupowfs[ipowfs].pts(), cuwfs[iwfs].amp, ratio, stream);
		}else{
		    cuspmul(gradref, cuwfs[iwfs].GS0, phiout, 1, 'n', ratio, stream);
		}
		if(gradacc()!=gradref()){
		    curadd(gradacc, 1, gradref, 1.f/(Real)dtrat, stream);
		}
	    }   
	    if(parms->powfs[ipowfs].psfout){
		zfarr_push(simu->save->ztiltout[iwfs], simu->wfsisim, gradcalc, stream);
	    }
	    if(do_phy || parms->powfs[ipowfs].psfout || parms->powfs[ipowfs].dither || do_pistatout){/*physical optics */
		CUDA_CHECK_ERROR;
		wfsints(simu, phiout, gradref, iwfs, isim);
		CUDA_CHECK_ERROR;
		ctoc("shwfs");
	    }/*do phy */
	}
	if(dtrat_output){
	    Real rne=0, bkgrnd=0;
	    if(do_phy || parms->powfs[ipowfs].dither){
		/*signal level was already multiplied in ints. */
		curcell &ints=cuwfs[iwfs].ints;
		const int totpix=(powfs[ipowfs].pywfs)?powfs[ipowfs].pywfs->nside:(ints[0].N());//PyWFs and SHWFS
		if(save_ints){
		    zfarr_push(simu->save->intsnf[iwfs], simu->wfsisim, ints, stream);
		}
		if(noisy){
		    if(parms->save.gradnf->p[iwfs]){
			if(parms->powfs[ipowfs].type==1){//PWFS
			    pywfs_grad(gradcalc, cuwfs[iwfs].ints[0], cupowfs[ipowfs].saa, 
				       cuwfs[iwfs].isum, cupowfs[ipowfs].pyoff, powfs[ipowfs].pywfs,stream);
			}else{
			    shwfs_grad(gradcalc, cuwfs[iwfs].ints, cuwfs, cupowfs, parms, powfs, simu, iwfs, ipowfs, stream);
			}
			if(parms->powfs[ipowfs].phytype_sim<3){
			    zfarr_push(simu->save->gradnf[iwfs], isim, gradcalc, stream);
			}else{//CPU version
			    zfarr_push(simu->save->gradnf[iwfs], isim, simu->gradcl->p[iwfs]);
			}
		    }
		    rne=parms->powfs[ipowfs].rne;
		    bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
		    addnoise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
			(ints[0], nsa, totpix, bkgrnd, bkgrnd*parms->powfs[ipowfs].bkgrndc,
			 cuwfs[iwfs].bkgrnd2(), cuwfs[iwfs].bkgrnd2c(), 
			 cuwfs[iwfs].qe, rne, cuwfs[iwfs].custat);
		    ctoc("noise");
		    if(save_ints){
			zfarr_push(simu->save->intsny[iwfs], simu->wfsisim, ints, stream);
		    }
		}
		if(parms->powfs[ipowfs].i0save){
		    curcelladd(cuwfs[iwfs].intsout, 1, ints, 1, stream);
		}
		if(parms->powfs[ipowfs].dither && isim>=parms->powfs[ipowfs].dither_ogskip 
		   && parms->powfs[ipowfs].type==0 && parms->powfs[ipowfs].phytype_sim2==1){
		    real cs, ss;
		    dither_position(&cs, &ss, parms->sim.alfsm, parms->powfs[ipowfs].dtrat,
				    parms->powfs[ipowfs].dither_npoint, isim, simu->dither[iwfs]->deltam);
		    int npll=parms->powfs[ipowfs].dither_pllrat;
		    cuwfs[iwfs].dither.acc(simu->dither[iwfs], ints, cs, ss, npll, stream);
		    ctoc("dither");
		}
	    }
	    if(do_phy){
		if(parms->powfs[ipowfs].type==1){
		    pywfs_grad(gradcalc, cuwfs[iwfs].ints[0], cupowfs[ipowfs].saa, 
			       cuwfs[iwfs].isum, cupowfs[ipowfs].pyoff, powfs[ipowfs].pywfs,stream);
		    //cuwrite(gradcalc, "gradcalc"); exit(0);
		}else{
		    shwfs_grad(gradcalc, cuwfs[iwfs].ints, cuwfs, cupowfs, parms, powfs, simu, iwfs, ipowfs, stream);
		}
		ctoc("grad");
	    }else{
		if(noisy){
		    if(parms->save.gradnf->p[iwfs]){
			zfarr_push(simu->save->gradnf[iwfs], isim, gradacc, stream);	
		    }
		    if(!parms->powfs[ipowfs].usephy){//do not add noise for presimulation to physical optics
			add_geom_noise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
			    (gradacc, cuwfs[iwfs].neasim, nsa,cuwfs[iwfs].custat);
			ctoc("geom_noise");
		    }
		}
	    }

	}/*dtrat_output */
	//info("thread %ld gpu %d iwfs %d queued\n", thread_id(), cudata->igpu, iwfs);
	ctoc_final("wfs %d", iwfs);
	CUDA_CHECK_ERROR;
    }//for iwfs
}

void gpu_wfsgrad_sync(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    gpu_set(cuglobal->wfsgpu[iwfs]);
    Array<cuwfs_t> &cuwfs=cuglobal->wfs;
    stream_t &stream=cuwfs[iwfs].stream;
    const int isim=simu->wfsisim;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const int dtrat_output=((isim+1)%dtrat==0);
    CUDA_SYNC_STREAM;
    if(dtrat_output){
	const int save_gradgeom=parms->save.gradgeom->p[iwfs];
	const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
	curmat &gradacc=cuwfs[iwfs].gradacc;
	curmat &gradcalc=cuwfs[iwfs].gradcalc;
	dmat *gradcl=simu->gradcl->p[iwfs];
	if(do_phy){
	    if(parms->powfs[ipowfs].phytype_sim<3){//3 is handled in cpu.
		cp2cpu(&gradcl, gradcalc, stream);
	    }
	    if(save_gradgeom){//also do geom grad during phy grad sims
		zfarr_push(simu->save->gradgeom[iwfs], simu->wfsisim, gradacc, stream);
	    }
	    if(parms->plot.run && !draw_disabled){// && parms->powfs[ipowfs].lo){
		cp2cpu(&simu->ints->p[iwfs], cuwfs[iwfs].ints, stream);
	    }
	}else{
	    cp2cpu(&gradcl, gradacc, stream);
	}
    }
}
void gpu_save_pistat(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->wfsisim;
    for(int iwfs=0; iwfs<simu->parms->nwfs; iwfs++){
	gpu_set(cuglobal->wfsgpu[iwfs]);
	Array<cuwfs_t> &cuwfs=cuglobal->wfs;
	const int ipowfs=simu->parms->wfs[iwfs].powfs;
	stream_t &stream=cuwfs[iwfs].stream;
	if(parms->powfs[ipowfs].pistatout){
	    int nstep=(isim+1-parms->powfs[ipowfs].pistatstart);
	    if(nstep>0){
		curcell tmp=cuwfs[iwfs].pistatout;
		curcellscale(tmp, 1.f/(Real)nstep, stream);
		if(parms->sim.skysim){
		    cuwrite(tmp, "%s/pistat/pistat_seed%d_sa%d_x%g_y%g.bin",
			    dirskysim,simu->seed,
			    parms->powfs[ipowfs].order,
			    parms->wfs[iwfs].thetax*206265,
			    parms->wfs[iwfs].thetay*206265);
		}else{
		    cuwrite(tmp,"pistat_seed%d_wfs%d.bin", simu->seed,iwfs);
		}
		curcellscale(tmp, nstep, stream);
	    }
	}
	
	if(parms->powfs[ipowfs].i0save){
	    const int dtrat=parms->powfs[ipowfs].dtrat;
	    const int nstep=(isim+1-parms->powfs[ipowfs].phystep)/dtrat;
	    if(nstep>0){
		curcell tmp=cuwfs[iwfs].intsout;
		curcellscale(tmp, 1.f/(Real)nstep, stream);
		cuwrite(tmp,"ints_%d_wfs%d.bin", simu->seed,iwfs);
		curcellscale(tmp, nstep, stream);
	    }
	}
    }
}
