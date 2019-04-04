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

#include "common.h"
#include "setup_recon.h"
#include "recon_utils.h"

/**
   Setup ray tracing operator HXF from xloc to aperture ploc along DM fiting directions*/
static dspcell *
setup_fit_HXF(const FIT_T *fit){
    info("Generating HXF");TIC;tic;
    if(!fit->xloc) return 0;
    const int nfit=fit->thetax->nx;
    const int npsr=fit->xloc->nx;
    dspcell *HXF=dspcellnew(nfit, npsr);
#pragma omp parallel for collapse(2)
    for(int ifit=0; ifit<nfit; ifit++){
	for(int ips=0; ips<npsr; ips++){
	    const double hsi=fit->hs->p[ifit];
	    const double ht = fit->xloc->p[ips]->ht-fit->floc->ht;
	    const double scale=1.-ht/hsi;
	    double displace[2];
	    displace[0]=fit->thetax->p[ifit]*ht;
	    displace[1]=fit->thetay->p[ifit]*ht;
	    P(HXF,ifit,ips)=mkh(fit->xloc->p[ips], fit->floc, displace[0], displace[1], scale);
	}
    }
    toc2(" ");
    return HXF;
}

/**
   Setup ray tracing operator HA from aloc to aperture ploc along DM fiting direction*/
static dspcell *
setup_fit_HA(FIT_T *fit){
    const int nfit=fit->thetax->nx;
    const int ndm=fit->aloc->nx;
    dspcell *HA=dspcellnew(nfit, ndm);
    info("Generating HA ");TIC;tic;
#pragma omp parallel for collapse(2)
    for(int ifit=0; ifit<nfit; ifit++){
	for(int idm=0; idm<ndm; idm++){
	    const double hs=fit->hs->p[ifit];
	    const double ht=fit->aloc->p[idm]->ht-fit->floc->ht;
	    const double scale=1.-ht/hs;
	    double displace[2];
	    displace[0]=fit->thetax->p[ifit]*ht;
	    displace[1]=fit->thetay->p[ifit]*ht;
	    loc_t *loc=fit->floc;
	    if(fit->misreg && fit->misreg[ifit+idm*nfit]){
		loc=loctransform(loc, fit->misreg[ifit+idm*nfit]);
	    }
	    P(HA,ifit,idm)=mkh(fit->aloc->p[idm], loc, 
				 displace[0], displace[1], scale);
	    if(loc!=fit->floc){
		locfree(loc);
	    }
	}
    }
    toc2(" ");
    fit->actcpl=genactcpl(HA, fit->W1);
    //cpl accounts for floating actuators, but not stuck actuators.
    act_stuck(fit->aloc, fit->actcpl, fit->actfloat);
    //Do not modify HA by floating actuators, otherwise, HA*actinterp will not work.
    act_stuck(fit->aloc, HA, fit->actstuck);
  
    if(fit->flag.actinterp){
	fit->actinterp=act_extrap(fit->aloc, fit->actcpl, fit->flag.actthres);
    }else if(fit->actfloat){
	warning("There are float actuators, but fit.actinterp is off\n");
    }
    if(fit->actinterp){
	/*
	  DM fitting output a is extrapolated to edge actuators by
	  actinterp*a. The corresponding ray tracing from DM would be
	  HA*actinterp*a. We replace HA by HA*actinterp to take this into
	  account during DM fitting.
	*/
	info("Replacing HA by HA*fit->interp\n");
	
	dspcell *HA2=0;
	dcellmm(&HA2, HA, fit->actinterp, "nn", 1);
	dspcellfree(HA);
	HA=HA2;
    }
   
    return HA;
}
/**
   Setup fitting low rank terms that are in the NULL space of DM fitting
   operator. typically include piston on each DM and tip/tilt on certain
   DMs. Becareful with tip/tilt contraint when using CBS.  */
static void 
setup_fit_lrt(FIT_T *fit){
    const int ndm=fit->aloc->nx;
    fit->NW=dcellnew(ndm,1);
    //double fitscl;     /**<strength of fitting FLM low rank terms (vectors)*/
    double fitscl=1./fit->floc->nloc;
    if(fabs(fitscl)<1.e-15){
	error("fit->fitscl is too small\n");
    }
    int nnw=0;
    if(fit->flag.lrt_piston){
	nnw+=ndm;
    }
    if(fit->flag.lrt_tt){
	nnw+=2*(ndm-1);
    }
    if(nnw==0) return;
    dcell* actcpl=dcelldup(fit->actcpl);
    //include stuck actuator
    act_stuck(fit->aloc, actcpl, fit->actstuck);
    for(int idm=0; idm<ndm; idm++){
	int nloc=fit->aloc->p[idm]->nloc;
	fit->NW->p[idm]=dnew(nloc, nnw);
    }
    int inw=0;/*current column */
    if(fit->flag.lrt_piston){
	info("Adding piston cr to fit matrix\n");
	for(int idm=0; idm<ndm; idm++){
	    int nloc=fit->aloc->p[idm]->nloc;
	    double *p=fit->NW->p[idm]->p+(inw+idm)*nloc;
	    const double *cpl=actcpl->p[idm]->p;
	    for(int iloc=0; iloc<nloc; iloc++){
		if(cpl[iloc]>0.1){ //don't count floating or stuck actuators
		    p[iloc]=fitscl;
		}
	    }
	}
	inw+=ndm;
    }
    if(fit->flag.lrt_tt){
	double factor=0;
	info("Adding TT cr on upper DMs to fit matrix.\n");
	factor=fitscl*2./loc_diam(fit->aloc->p[0]);
	for(int idm=1; idm<ndm; idm++){
	    int nloc=fit->aloc->p[idm]->nloc;
	    double *p=fit->NW->p[idm]->p+(inw+(idm-1)*2)*nloc;
	    double *p2x=p;
	    double *p2y=p+nloc;
	    const double *cpl=actcpl->p[idm]->p;
	    for(int iloc=0; iloc<nloc; iloc++){
		if(cpl[iloc]>0.1){
		    p2x[iloc]=fit->aloc->p[idm]->locx[iloc]*factor;/*x tilt */
		    p2y[iloc]=fit->aloc->p[idm]->locy[iloc]*factor;/*y tilt */
		}
	    }
	}
	inw+=2*(ndm-1);
    }
    if(fit->flag.actslave){
	/*
	  2011-07-19: When doing PSFR study for MVR with SCAO, NGS. Found
	  that slaving is causing mis-measurement of a few edge
	  actuators. First try to remove W1. Or lower the weight. Revert
	  back.
	  1./floc->nloc is on the same order of norm of Ha'*W*Ha. 
	*/
	TIC;tic;
	fit->actslave=slaving(fit->aloc, fit->actcpl,
			      fit->NW, fit->actstuck,
			      fit->actfloat, fit->flag.actthres, 1./fit->floc->nloc);
	toc2("slaving");
    }
    cellfree(actcpl);
}
/**
   Assemble the DM fitting matrix

   The fitting is done by minimizing \f$||H_X x - H_A a||^2_W\f$ where \f$H_X,
   H_A\f$ are ray tracing operator from tomography grid xloc, and deformable
   mirror grid aloc to pupil grid ploc. The norm is weighted using bilinear
   influence functions within the telescope aperture. We have
   
   \f$a=\left[H_A^T(W_0-W_1 W_1^T)H_A\right]^{-1} H_A^T (W_0-W_1) H_X x\f$

   For details see www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803 
*/
static void
setup_fit_matrix(FIT_T *fit){
    const int nfit=fit->thetax->nx;
    const int ndm=fit->aloc->nx;
    if(ndm==0) return;
    
    dspcell* HA=fit->HA;
    dspcell *HAT=dspcelltrans(HA);
    
    info("Before assembling fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    /*Assemble Fit matrix. */
    if(!fit->FR.M && fit->flag.assemble){
	if(fit->HXF){//not idealfit.
	    const int npsr=fit->xloc->nx;
	    info("Building fit->FR\n");
	    fit->FR.M=cellnew(ndm, npsr);
	    dspcell* FRM=(dspcell*)fit->FR.M;
	    dspcell* HXF=fit->HXF;
	    //FRM
	    for(int ips=0; ips<npsr; ips++){
		for(int ifit=0; ifit<nfit; ifit++){
		    if(fabs(fit->wt->p[ifit])<1.e-12) continue;
		    dsp *tmp=dspmulsp(fit->W0, P(HXF,ifit,ips),"nn");
		    for(int idm=0; idm<ndm; idm++){
			dspmulsp2(PP(FRM,idm,ips),P(HAT,idm,ifit), tmp, "nn",
				  fit->wt->p[ifit]);
		    }
		    dspfree(tmp);
		}
	    }
	    fit->FR.V=dcellnew(npsr, 1);
	    dmat **FRV=fit->FR.V->p;  
	    //FRV
	    for(int ips=0; ips<npsr; ips++){
		int nloc=fit->xloc->p[ips]->nloc;
		FRV[ips]=dnew(nloc,nfit);
		for(int ifit=0; ifit<nfit; ifit++){
		    /*notice the sqrt. */
		    if(fabs(fit->wt->p[ifit])<1.e-12) continue;
		    dspmulvec(FRV[ips]->p+ifit*nloc, 
			      P(HXF,ifit,ips), fit->W1->p, 't',
			      sqrt(fit->wt->p[ifit]));
		}
	    }
	    cellfree(fit->HXF);  
	}else{
	    dbg("Avoid building fit->FR.M\n");
	    fit->FR.M=NULL;
	    fit->FR.V=NULL;
	}
	/*Always need FR.U as it is used to do FL.U, FL.V */
	fit->FR.U=dcellnew(ndm, 1);
	dmat **FRU=fit->FR.U->p;
	
	for(int idm=0; idm<ndm; idm++){    
	    int nloc=fit->aloc->p[idm]->nloc;
	    FRU[idm]=dnew(nloc, nfit);
	    for(int ifit=0; ifit<nfit; ifit++){
		/*notice the sqrt. */
		if(fabs(fit->wt->p[ifit])<1.e-12) continue;
		dspmulvec(FRU[idm]->p+ifit*nloc, 
			  P(HA,ifit,idm), fit->W1->p,'t',
			  sqrt(fit->wt->p[ifit]));
	    }
	}
    }

    if(!fit->FL.M){
	info("Building fit->FL\n");
	fit->FL.M=cellnew(ndm, ndm);
	dspcell *FLM=(dspcell*)fit->FL.M;
	for(int idm=0; idm<ndm; idm++){
	    for(int ifit=0; ifit<nfit; ifit++){
		if(fabs(fit->wt->p[ifit])<1.e-12) continue;
		dsp *tmp=dspmulsp(fit->W0, P(HA,ifit,idm),"nn");
		for(int jdm=0; jdm<ndm; jdm++){
		    dspmulsp2(PP(FLM,jdm,idm),P(HAT,jdm,ifit), tmp,"nn",
			      fit->wt->p[ifit]);
		}
		dspfree(tmp);
	    }
	}

	if(fabs(fit->flag.tikcr)>1.e-15){
	    double tikcr=fit->flag.tikcr;
	    /*Estimated from the formula.  1/nloc is due to W0, the other
	      scaling is due to ray tracing between different sampling freq.*/
	    int nact=0;
	    for(int idm=0; idm<ndm; idm++){
		nact+=fit->aloc->p[idm]->nloc;
	    }
	    double maxeig=4./nact;
	    info("Adding tikhonov constraint of %g to FLM\n", tikcr);
	    info("The maximum eigen value is estimated to be around %e\n", maxeig);
	    dcelladdI(fit->FL.M,tikcr*maxeig);
	}

	{/*Low rank terms. */
	    fit->FL.U=dcellcat_each(fit->FR.U, fit->NW, 2);
	    dcell *tmp=NULL;/*negative NW. */
	    dcelladd(&tmp, 1, fit->NW, -1);
	    fit->FL.V=dcellcat_each(fit->FR.U, tmp, 2);
	    dcellfree(tmp);
	}
	if(fit->actslave){
	    dcelladd(&fit->FL.M, 1, fit->actslave, 1);
	}
	/*dspcellsym(fit->FL.M); */
	info("DM Fit number of Low rank terms: %ld in LHS\n", fit->FL.U->p[0]->ny);
    }
    dspcellfree(HAT);
    if(fit->flag.alg==0 || fit->flag.alg==2){
	if(fit->flag.alg==0 && fabs(fit->flag.tikcr)<1.e-14){
	    warning("tickcr=%g is too small, chol may fail.\n", fit->flag.tikcr);
	}
	if(fit->flag.bgs){
	    muv_direct_diag_prep(&(fit->FL),(fit->flag.alg==2)*fit->flag.svdthres);
	}else{
	    muv_direct_prep(&(fit->FL),(fit->flag.alg==2)*fit->flag.svdthres);
	    cellfree(fit->FL.M);
	    dcellfree(fit->FL.U);
	    dcellfree(fit->FL.V);
	}
	info("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }
  
    info("After assemble fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
}
/**
   A generic DM fitting routine.
 */
void setup_fit(FIT_T *fit, int idealfit){
    TIC;tic;
    if(!idealfit && fit->xloc){
	fit->HXF=setup_fit_HXF(fit);
    }
    fit->HA=setup_fit_HA(fit);
    setup_fit_lrt(fit);
    
    /*always assemble fit matrix, faster if many directions */
    if(fit->flag.assemble || fit->flag.alg!=1){
	setup_fit_matrix(fit);
    }
    /*Fall back function method if FR.M is NULL (!HXF<-idealfit) */
    fit->FR.Mfun  = FitR;
    fit->FR.Mdata = fit;
    /*Fall back function method if FL.M is NULL */
    fit->FL.Mfun  = FitL;
    fit->FL.Mdata = fit;
    fit->FL.alg   = fit->flag.alg;
    fit->FL.bgs   = fit->flag.bgs;
    fit->FL.warm  = fit->flag.cgwarm;
    fit->FL.maxit = fit->flag.maxit;
    toc2("Setting up DM Fitting.");
}
void free_fit(FIT_T *fit, int nfit){
    if(!fit) return;
    for(int ifit=0; ifit<nfit; ifit++){
	cellfree(fit[ifit].HXF);
	cellfree(fit[ifit].HA);
	cellfree(fit[ifit].actcpl);
	cellfree(fit[ifit].actinterp);
	cellfree(fit[ifit].actslave);
	cellfree(fit[ifit].NW);
	muv_free(&fit[ifit].FR);
	muv_free(&fit[ifit].FL);
    }
    free(fit);
}
void setup_recon_fit(RECON_T *recon, const PARMS_T *parms){
    FIT_T *fit=mycalloc(1, FIT_T);
    recon->fit=fit;
    fit->thetax=parms->fit.thetax;
    fit->thetay=parms->fit.thetay;
    fit->wt=parms->fit.wt;
    fit->hs=parms->fit.hs;

    fit->xloc=recon->xloc;
    fit->floc=recon->floc;
    fit->aloc=recon->aloc;
    fit->W0=recon->W0;
    fit->W1=recon->W1;

    fit->actfloat=recon->actfloat;
    fit->actstuck=recon->actstuck;
    fit->misreg=parms->recon.misreg_dm2sci;
    memcpy(&fit->flag, &parms->fit, sizeof(FIT_CFG_T));//use parms->fit.
    if(parms->fit.assemble){
	if(parms->load.fit){
	    if(!(zfexist("FRM") && zfexist("FRU") && zfexist("FRV"))){
		error("FRM, FRU, FRV (.bin) not all exist\n");
	    }
	    if(!(zfexist("FLM") && zfexist("FLU") && zfexist("FLV"))){
		error("FLM, FLU, FLV (.bin) not all exist\n");
	    }
	    fit->FR.M=readbin("FRM");
	    fit->FR.U=dcellread("FRU");
	    fit->FR.V=dcellread("FRV");
	    fit->FL.M=readbin("FLM");
	    fit->FL.U=dcellread("FLU");
	    fit->FL.V=dcellread("FLV");
	    
	}
    }
    setup_fit(fit, parms->sim.idealfit);
    
    dcellfree(recon->actcpl);
    recon->actcpl=dcellref(fit->actcpl);
    dcellfree(recon->actinterp);
    recon->actinterp=dspcellref(fit->actinterp);
    
    if(parms->save.setup){
	writebin(fit->HA,"HA");
	writebin(recon->actinterp, "actinterp");
	writebin(recon->actcpl, "actcpl");
    }
    if(parms->save.recon){
	writebin(fit->FR.M,"FRM");
	writebin(fit->FR.V,"FRV");
	writebin(fit->FR.U,"FRU");
	  
       	if(fit->FL.C){
	    chol_convert(fit->FL.C, 1);
	    chol_save(fit->FL.C,"FLC.bin");
	}
	if(fit->FL.MI)
	    writebin(fit->FL.MI,"FLMI");
	if(fit->FL.Up)
	    writebin(fit->FL.Up, "FLUp");
	if(fit->FL.Vp)
	    writebin(fit->FL.Vp, "FLVp");
	if(fit->FL.CB){
	    for(int ib=0; ib<fit->FL.nb; ib++){
		chol_save(fit->FL.CB[ib],"FLCB_%d.bin", ib);
	    }
	}
	if(fit->FL.MIB){
	    writebin(fit->FL.MIB,"FLMIB");
	}
   
    }
}
/**
   Setting fitting parameter for turbulence to WFS lenslet grid.
 */
void setup_powfs_fit(POWFS_T *powfs, const RECON_T *recon, const PARMS_T *parms){
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].lo) continue;
	int nwfs=parms->powfs[ipowfs].nwfs;
	FIT_T *fitall=powfs[ipowfs].fit=mycalloc(nwfs, FIT_T);
	loc_t *wfsloc=mkannloc(parms->aper.d+parms->powfs[ipowfs].dsa*2, 0, parms->powfs[ipowfs].dsa, 0);
	wfsloc->ht=parms->powfs[ipowfs].hc;
	wfsloc->iac=parms->dbg.wfs_iac;//cubic spline better fits the turbulence.
	for(int jwfs=0; jwfs<nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    FIT_T *fit=fitall+jwfs;
	    if(jwfs==0){
		memcpy(&fit->flag, &parms->fit, sizeof(FIT_CFG_T));//use parms->fit.
		fit->flag.alg=0;
		fit->flag.assemble=0;
		fit->notrecon=1; //not for reconstruction
		fit->wt=dnew(1,1); fit->wt->p[0]=1;
		fit->hs=dnew(1,1); fit->hs->p[0]=parms->powfs[ipowfs].hs;
	
		fit->aloc=loccellnew(1,1); fit->aloc->p[0]=locref(wfsloc);
		fit->floc=locref(recon->floc);
		fit->W0=recon->W0;
		fit->W1=recon->W1;

		fit->thetax=dnew(1,1);fit->thetax->p[0]=parms->wfs[iwfs].thetax;
		fit->thetay=dnew(1,1);fit->thetay->p[0]=parms->wfs[iwfs].thetay;
		setup_fit(fit, 1);
	    }else{
		memcpy(fitall+jwfs, fitall, sizeof(FIT_T));
		fit->FR.Mdata=fit;
		fit->FL.Mdata=fit;
		fit->thetax=dnew(1,1);fit->thetax->p[0]=parms->wfs[iwfs].thetax;
		fit->thetay=dnew(1,1);fit->thetay->p[0]=parms->wfs[iwfs].thetay;
	    }
	}
	locfree(wfsloc);
    }
}
