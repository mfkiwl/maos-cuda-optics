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
#include "../math/mathdef.h"
#include "locfft.h"
#define TIMING 0
#if TIMING == 1
#define TIM(A) real tk##A=myclockd()
#else
#define TIM(A)
#endif

/**
   For routines that embed OPDs defined on loc to square(rectangular) array and do fft on it.
*/


locfft_t *locfft_init(loc_t *loc,       /**<[in] The loc*/
		      const dmat *amp,        /**<[in] The amplitude*/
		      const dmat *wvl,        /**<[in] The wavelength*/
		      const lmat *fftsize,    /**<[in] The suggested size for FFT*/
		      const real oversize,  /**<[in] Factor of oversize. 2 fot FFT*/
		      real fieldstop        /**<[in] Size of field stop (radian) if used*/
    ){
    const int nwvl=wvl->nx*wvl->ny;
    locfft_t *locfft=mycalloc(1, locfft_t);
    locfft->embed=lcellnew(nwvl, 1);
    locfft->nembed=lnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	if(iwvl==0 || (fftsize && fftsize->p[iwvl]>0 && fftsize->p[iwvl]!=locfft->nembed->p[0])){
	    locfft->nembed->p[iwvl]=fftsize?fftsize->p[iwvl]:0;
	    locfft->embed->p[iwvl]=loc_create_embed(&locfft->nembed->p[iwvl], loc, oversize, 1);
	}else{
	    locfft->embed->p[iwvl]=lref(locfft->embed->p[0]);
	    locfft->nembed->p[iwvl]=locfft->nembed->p[0];
	}
    }
    locfft->wvl=dref_reshape(wvl, nwvl, 1);
    locfft->amp=amp;
    locfft->loc=loc;
    locfft->ampsum=dsum(amp);
    locfft->ampnorm=dsumsq(amp);
    if(fieldstop){
	locfft->fieldstop=fieldstop;
	locfft->fieldmask=dcellnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    int nembed=locfft->nembed->p[iwvl];
	    locfft->fieldmask->p[iwvl]=dnew(nembed, nembed);
	    real dtheta=wvl->p[iwvl]/(loc->dx*nembed);//sampling of psf
	    real radius=fieldstop/(dtheta*2);
	    dcircle(locfft->fieldmask->p[iwvl], nembed/2+1, nembed/2+1, 1, 1, radius, 1);
	    dfftshift(locfft->fieldmask->p[iwvl]);
	}
    }
    return locfft;
}
/**
   Frees the struct
*/
void locfft_free(locfft_t *locfft){
    if(!locfft) return;
    cellfree(locfft->embed);
    cellfree(locfft->nembed);
    cellfree(locfft->fieldmask);
    cellfree(locfft->wvl);
    free(locfft);
}
/**
   Computes strehl from OPD without doing FFT. The strehl is simply 
   
   \f$s=\sum(A*exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]) \f$
   
   where A is the amplitude map.
*/
static comp strehlcomp(const dmat *iopdevl, const dmat *amp, const real wvl){
    comp i2pi=COMPLEX(0, 2*M_PI/wvl);
    comp strehl=0;
    for(int iloc=0; iloc<iopdevl->nx; iloc++){
	strehl+=amp->p[iloc]*cexp(i2pi*iopdevl->p[iloc]);
    }
    return strehl;
}
/**
   Computes PSF from OPD by FFT. The PSF is computed as

   \f$\textrm{PSF}=\mathcal{F}[A\times exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]]\f$

   The peak value (center) in the computed PSF is normalized by the peak value
   in the differaction limited PSF. In other words, the peak value in the
   computed PSF is the Strehl. Keep this in mind when you compute enclosed
   energy.
   
   Extract center part of psfsize.
*/
void locfft_psf(ccell **psf2s, const locfft_t *locfft, const dmat *opd, const lmat *psfsize, int sum2one){
    long nwvl=locfft->wvl->nx;
    if(!*psf2s){
	*psf2s=ccellnew(nwvl, 1);
    }
    if(opd->nx!=locfft->amp->nx){
	error("The length of opd should be %ld, but is %ld\n", locfft->amp->nx, opd->nx);
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++)
#if _OPENMP>=200805
#pragma omp task
#endif
    {
	if(psfsize && psfsize->p[iwvl]==1){
	    if(!(*psf2s)->p[iwvl]){
		(*psf2s)->p[iwvl]=cnew(1,1);
	    }
	    (*psf2s)->p[iwvl]->p[0]=strehlcomp(opd, locfft->amp, locfft->wvl->p[iwvl]);
	}else{
	    TIM(0);
	    long nembed=locfft->nembed->p[iwvl];
	    long *embed=locfft->embed->p[iwvl]->p;
	    const real *amp=locfft->amp->p;
	    const int ref=!psfsize || psfsize->p[iwvl]==nembed;
	    cmat *psf2=0;
	    if(ref){//Full PSF is returned
		if(!(*psf2s)->p[iwvl]){
		    (*psf2s)->p[iwvl]=cnew(nembed,nembed);
		}
		psf2=(*psf2s)->p[iwvl];
	    }else{//Crop of PSF is returned.
		psf2=cnew(nembed,nembed);
	    }

	    int use1d=0;
#define USE1D_ENABLED 1
#if     USE1D_ENABLED
	    if(psfsize && psfsize->p[iwvl]+200<nembed){/*Want smaller PSF. */
		use1d=1;
	    }
#endif

	    comp i2pi=COMPLEX(0, 2*M_PI/locfft->wvl->p[iwvl]);
	    for(int iloc=0; iloc<opd->nx; iloc++){
		psf2->p[embed[iloc]]=amp[iloc]*cexp(i2pi*opd->p[iloc]);
	    }
	    TIM(1);
	    if(use1d==1){
		cfft2partial(psf2, psfsize->p[iwvl], -1);
	    }else{
		cfft2(psf2,-1);
	    }
	    TIM(2);
	    if(ref){/*just reference */
		cfftshift(psf2);
	    }else{/*create a new array, smaller. */
		if(!(*psf2s)->p[iwvl]){
		    (*psf2s)->p[iwvl]=cnew(psfsize->p[iwvl], psfsize->p[iwvl]);
		}
		ccpcorner2center((*psf2s)->p[iwvl], psf2);
		cfree(psf2); 
	    }
	    TIM(3);
#if TIMING
	    info("locfft_psf(%d:%ldx%ld): exp %.4f, fft %.4f (%.2f GFLOPS), abs2 %.2f.\n", iwvl,nembed, nembed,
		 tk1-tk0, tk2-tk1, 8L*(use1d?psfsize->p[iwvl]:nembed)*nembed*log2(nembed)/(tk2-tk1)*1e-9, tk3-tk2);
#endif		
	}
	real psfnorm;
	if(sum2one){/**PSF sum to 1*/
	    psfnorm=1./(sqrt(locfft->ampnorm)*locfft->nembed->p[iwvl]);
	}else{/**PSF max is strehl*/
	    psfnorm=1./locfft->ampsum;
	}
	if(fabs(psfnorm-1)>1.e-15) {
	    cscale((*psf2s)->p[iwvl], psfnorm);
	}
    }
#if _OPENMP>=200805
#pragma omp taskwait
#endif
}
/**
   Apply a field stop to the OPD.
*/
void locfft_fieldstop(const locfft_t *locfft, dmat *opd, const dmat *wvlwts){
    int nwvl=locfft->wvl->nx;
    if(nwvl>1){
	warning("Not tested for multi-wavelength case yet.\n");
    }
    ccell *wvfs=ccellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	int nembed=locfft->nembed->p[iwvl];
	lmat *embed=locfft->embed->p[iwvl];
	cmat *wvf=cnew(nembed, nembed);
	wvfs->p[iwvl]=wvf;
	//cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
	real wvl=locfft->wvl->p[iwvl];
	comp i2pi=COMPLEX(0,2*M_PI/wvl);
	const real *amp=locfft->amp->p;
	for(int iloc=0; iloc<opd->nx; iloc++){
	    wvf->p[embed->p[iloc]]=amp[iloc]*cexp(i2pi*opd->p[iloc]);
	}
	cfft2(wvf, -1);
	ccwmd(wvf, locfft->fieldmask->p[iwvl], 1);
	cfft2(wvf, 1);
    }
    if(nwvl>1){
	/*Please fix this case. Need to do phase unwrapping first and average OPD
	 * for different wavelength result*/
	error("Not implemented yet\n");
    }
    dmat *opdold=ddup(opd); dzero(opd);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	real wvl=locfft->wvl->p[iwvl];
	real wvlh=wvl*0.5;
	real kki=wvl/(2*M_PI);
	cmat *wvf=wvfs->p[iwvl];
	lmat *embed=locfft->embed->p[iwvl];
	for(int iloc=0; iloc<opd->nx; iloc++){
	    real val=carg(wvf->p[embed->p[iloc]])*kki;
	    if(fabs(val-opdold->p[iloc])>wvlh){//need phase unwrapping
		warning_once("phase unwrapping is needed\n");
		real diff=fmod(val-opdold->p[iloc]+wvlh, wvl);
		if(diff<0) diff+=wvl;
		val=(diff-wvlh)+opdold->p[iloc];
	    }
	    opd->p[iloc]+=wvlwts->p[iwvl]*val;
	}
    }
    ccellfree(wvfs);
}
