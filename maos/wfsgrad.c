/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "sim.h"
#include "sim_utils.h"
#include "ahst.h"
#include "mtch.h"
#include "save.h"
#include "setup_recon.h"
#include "setup_powfs.h"
#include "pywfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file wfsgrad.c
   contains functions that computes WFS gradients in geometric or physical optics mode.
*/
#define TIMING 0
#if TIMING == 1
#define TIM(A) double tk##A=myclockd()
#else
#define TIM(A)
#endif

static void wfs_ideal_atm(SIM_T *simu, dmat *opd, int iwfs, double alpha){
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const double hs=parms->wfs[iwfs].hs;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    for(int idm=0; idm<parms->ndm; idm++){
	loc_t *loc=powfs[ipowfs].loc_dm?powfs[ipowfs].loc_dm->p[wfsind+idm*parms->nwfs]:powfs[ipowfs].loc;
	double *amp=powfs[ipowfs].realamp->p[wfsind]->p;
	const double ht = parms->dm[idm].ht+parms->dm[idm].vmisreg;
	double dispx=ht*parms->wfs[iwfs].thetax;
	double dispy=ht*parms->wfs[iwfs].thetay;
	double scale=1.-ht/hs;
	if(parms->dm[idm].cubic){
	    prop_grid_cubic(simu->dmprojsq->p[idm],
			    loc, amp, opd->p, 
			    alpha, dispx, dispy, scale, parms->dm[idm].iac, 
			    0, 0);
	}else{
	    prop_grid(simu->dmprojsq->p[idm],
		      loc, amp, opd->p, 
		      alpha, dispx, dispy, scale, 0,
		      0, 0);
	}
    }
}

/**
   computes close loop and pseudo open loop gradidents for both gometric and
   physical optics WFS. Calls wfsints() to accumulate WFS subapertures images in
   physical optics mode.  */

void wfsgrad_iwfs(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const int isim=simu->isim;
    const int iwfs=info->start;
    const PARMS_T *parms=simu->parms;
    const int ipowfs=parms->wfs[iwfs].powfs;
    if(isim<parms->powfs[ipowfs].step) return;
    assert(iwfs<parms->nwfs);
    /*
      simu->gradcl is CL grad output (also for warm-restart of maxapriori
      simu->gradacc is internal, to accumulate geometric grads.
      do not accumulate opd. accumate ints for phy, g for GS
    */
    /*input */
    
    mapcell *atm=simu->atm;
    const RECON_T *recon=simu->recon;
    const POWFS_T *powfs=simu->powfs;
    /*output */
    const int CL=parms->sim.closeloop;
    const int nps=parms->atm.nps;
    const double atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
    const double dt=simu->dt;
    TIM(0);
    /*The following are truly constants for this powfs */
    const int imoao=parms->powfs[ipowfs].moao;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const double hs=parms->wfs[iwfs].hs;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const int save_gradgeom=parms->save.gradgeom->p[iwfs];
    const int save_opd =parms->save.wfsopd->p[iwfs];
    const int save_ints=parms->save.ints->p[iwfs];
    const int noisy=parms->powfs[ipowfs].noisy;
    /*The following depends on isim */
    /*const int dtrat_reset=(isim%dtrat==0); */
    const int dtrat_output=(isim+1-parms->powfs[ipowfs].step)%dtrat==0;
    const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
    const int do_pistatout=parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart;
    const int do_geom=!do_phy || save_gradgeom || do_pistatout;
    const double *realamp=powfs[ipowfs].realamp?powfs[ipowfs].realamp->p[wfsind]->p:0;
    dmat *gradcalc=NULL;
    dmat **gradacc=&simu->gradacc->p[iwfs];
    dmat **gradout=&simu->gradcl->p[iwfs];
    dcell *ints=simu->ints->p[iwfs];
    dmat  *opd=simu->wfsopd->p[iwfs];
    dzero(opd);
    if((isim-parms->powfs[ipowfs].step)%dtrat==0){
	dcellzero(ints);
	dzero(*gradacc);
    }
    if(simu->telws){/*Wind shake */
	double tmp=simu->telws->p[isim];
	double angle=simu->winddir?simu->winddir->p[0]:0;
	double ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
	loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
    }

    /* Add surface error*/
    if(powfs[ipowfs].opdadd && powfs[ipowfs].opdadd->p[wfsind]){
	dadd(&opd,1, powfs[ipowfs].opdadd->p[wfsind],1);
    }

    /* Now begin ray tracing. */
    if(parms->sim.idealwfs && !parms->powfs[ipowfs].lo){
	wfs_ideal_atm(simu, opd, iwfs, 1);
    }else if(atm){
	for(int ips=0; ips<nps; ips++){
	    thread_t *wfs_prop=simu->wfs_prop_atm[iwfs+parms->nwfs*ips];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
	    wfs_propdata->phiout=opd->p;
	    wfs_propdata->displacex1=-atm->p[ips]->vx*dt*isim;
	    wfs_propdata->displacey1=-atm->p[ips]->vy*dt*isim;
	    wfs_propdata->alpha=atmscale;
	    /* have to wait to finish before another phase screen. */
	    CALL_THREAD(wfs_prop, 0);
	}
	/* most expensive 0.10 per LGS for*/
	if(parms->sim.wfsalias){
	    /* Remove subspace of atm projected onto range of DM.*/
	    wfs_ideal_atm(simu, opd, iwfs, -1);
	}
    }
    if(save_opd){
	cellarr_dmat(simu->save->wfsopdol[iwfs], isim, opd);
    }
 
    if(CL){
	for(int idm=0; idm<parms->ndm; idm++){
	    thread_t *wfs_prop=simu->wfs_prop_dm[iwfs+parms->nwfs*idm];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_dm[iwfs+parms->nwfs*idm];
	    wfs_propdata->phiout=opd->p;
	    CALL_THREAD(wfs_prop, 0);
	}/*idm */
	double ptt[3]={0,0,0};
	if(simu->ttmreal){
	    ptt[1]-=simu->ttmreal->p[0];
	    ptt[2]-=simu->ttmreal->p[1];
	}
	//For dithering with downlink instead of uplink FSM
	if(simu->fsmreal && simu->fsmreal->p[iwfs] && !powfs[ipowfs].llt){
	    ptt[1]-=simu->fsmreal->p[iwfs]->p[0];
	    ptt[2]-=simu->fsmreal->p[iwfs]->p[1];
	}
	if(ptt[1] || ptt[2]){
	    loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
	}
    }
    if(parms->powfs[ipowfs].skip && parms->tomo.ahst_idealngs){
	//apply ideal NGS modes to NGS WFS
	ngsmod2science(opd, powfs[ipowfs].loc, recon->ngsmod, 
		       parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay,
		       simu->cleNGSm->p+isim*recon->ngsmod->nmod, -1);
    }
    if(imoao>-1){
	dmat **dmwfs=simu->dm_wfs->p;
	if(dmwfs[iwfs]){
	    /* No need to do mis registration here since the MOAO DM is attached
	       to close to the WFS.*/
	    if(parms->moao[imoao].cubic){
		prop_nongrid_pts_cubic(recon->moao[imoao].aloc->p[0], dmwfs[iwfs]->p,
				       powfs[ipowfs].pts, realamp, opd->p, -1, 0, 0, 1, 
				       parms->moao[imoao].iac, 0, 0);
	    }else{
		prop_nongrid_pts(recon->moao[imoao].aloc->p[0], dmwfs[iwfs]->p,
				 powfs[ipowfs].pts, realamp, opd->p, -1, 0, 0, 1, 
				 0, 0);
	    }
	}
    }
    /* Add defocus to OPD if needed. */
    if(parms->powfs[ipowfs].llt){
	double focus=wfsfocusadj(simu, iwfs);
	if(fabs(focus)>1e-20){
	    loc_add_focus(opd->p, powfs[ipowfs].loc, focus);
	}
    }
    if(parms->powfs[ipowfs].fieldstop>0 && parms->powfs[ipowfs].type==0){
	locfft_fieldstop(powfs[ipowfs].fieldstop, opd, parms->powfs[ipowfs].wvlwts);
    }

    if(save_opd){
	cellarr_dmat(simu->save->wfsopd[iwfs], isim, opd);
    }
    if(parms->plot.run){
	drawopdamp("wfsopd",powfs[ipowfs].loc,opd->p,realamp,NULL,
		   "WFS OPD","x (m)", "y (m)", "WFS %d", iwfs);
    }
    if(do_geom){
	/* Now Geometric Optics gradient calculations. if dtrat==1, we compute
	   gradients directly to gradacc, which is the same as gradcalc. If
	   dtrat>1, we compute gradients to gradcalc, and accumulate to
	   gradacc. gradcalc is used to shift pistat. We DONOT include gradoff
	   adjustment to gradref, but only do it on gradcl. This will make the
	   pistat always peak in center no matter what NCPA is present.
	*/
	if(!do_pistatout || parms->powfs[ipowfs].pistatstc || dtrat==1){
	    gradcalc=dref(*gradacc);
	}//else: calculate first to gradcalc then add to gradacc
	if(parms->powfs[ipowfs].gtype_sim==1){ /*compute ztilt. */
	    pts_ztilt(&gradcalc,powfs[ipowfs].pts,
		      powfs[ipowfs].saimcc->p[powfs[ipowfs].nsaimcc>1?wfsind:0], 
		      realamp, opd->p);
	}else{/*G tilt */
	    dspmm(&gradcalc,adpind(powfs[ipowfs].GS0,wfsind),opd,'n',1);
	}
	if(gradcalc->p!=(*gradacc)->p){
	    dadd(gradacc, 1, gradcalc, 1);
	}
    }

    ccell *psfout=NULL;
    cellarr *psfoutcellarr=NULL;
    cellarr *ztiltoutcellarr=NULL;
    if(parms->powfs[ipowfs].psfout){
	psfout=simu->wfspsfout->p[iwfs];
	psfoutcellarr=simu->save->wfspsfout[iwfs];
	ztiltoutcellarr=simu->save->ztiltout[iwfs];
    }
    TIM(1);
    /* Now begin Physical Optics Intensity calculations */
    if(do_phy || psfout || do_pistatout || parms->powfs[ipowfs].dither){
	dmat *lltopd=NULL;
	if(powfs[ipowfs].llt){//If there is LLT, apply FSM onto LLT
	    if(powfs[ipowfs].llt->ncpa){
		int iotf=powfs[ipowfs].llt->ncpa->nx==1?0:wfsind;
		lltopd=ddup(powfs[ipowfs].llt->ncpa->p[iotf]);
	    }else{
		lltopd=dnew(powfs[ipowfs].llt->pts->nx,
			    powfs[ipowfs].llt->pts->nx);
	    }
	    const long illt=parms->powfs[ipowfs].llt->i->p[wfsind];
	    if(atm){/*LLT OPD */
		for(int ips=0; ips<nps; ips++){
		    const double hl=atm->p[ips]->h;
		    const double scale=1.-hl/hs;
		    const double thetax=parms->wfs[iwfs].thetax-parms->powfs[ipowfs].llt->ox->p[illt]/hs;
		    const double thetay=parms->wfs[iwfs].thetay-parms->powfs[ipowfs].llt->oy->p[illt]/hs;
		    const double displacex=-atm->p[ips]->vx*isim*dt+thetax*hl+parms->powfs[ipowfs].llt->misreg->p[0];
		    const double displacey=-atm->p[ips]->vy*isim*dt+thetay*hl+parms->powfs[ipowfs].llt->misreg->p[1];
		    prop_grid_pts(atm->p[ips],powfs[ipowfs].llt->pts,NULL,
				  lltopd->p,atmscale,displacex,displacey,
				  scale, 1., 0, 0);
		}
	    }
	    double ttx=0, tty=0;//FSM + wind shake induced jitter
	    if((simu->fsmreal && simu->fsmreal->p[iwfs]) ||do_pistatout||parms->sim.fsmideal){
		if(do_pistatout||parms->sim.fsmideal){
		    /* remove tip/tilt completely */
		    dmat *lltg=dnew(2,1);
		    pts_ztilt(&lltg,powfs[ipowfs].llt->pts,
			      powfs[ipowfs].llt->imcc,
			      powfs[ipowfs].llt->amp->p,
			      lltopd->p);
		    simu->fsmreal->p[iwfs]->p[0]=-lltg->p[0];
		    simu->fsmreal->p[iwfs]->p[1]=-lltg->p[1];
		    dfree(lltg);
		}
		
		ttx=simu->fsmreal->p[iwfs]->p[0];
		tty=simu->fsmreal->p[iwfs]->p[1];
	    }
	    if(simu->telws){
		double tmp=simu->telws->p[isim]*parms->powfs[ipowfs].llt->ttrat;
		double angle=simu->winddir?simu->winddir->p[0]:0;
		ttx+=tmp*cos(angle);
		tty+=tmp*sin(angle);
	    }
	    if(ttx !=0 || tty != 0){ /* add tip/tilt to llt opd */
		double ptt[3]={0, ttx, tty};
		loc_add_ptt(lltopd->p, ptt, powfs[ipowfs].llt->loc);
	    }
	    if(save_opd){
		cellarr_dmat(simu->save->wfslltopd[iwfs], isim, lltopd);
	    }
	}
	if(parms->powfs[ipowfs].type==0){
	    WFSINTS_T *intsdata=simu->wfs_intsdata+iwfs;
	    intsdata->ints=ints;
	    intsdata->psfout=psfout;
	    intsdata->pistatout=simu->pistatout->p[iwfs];
	    if(parms->powfs[ipowfs].pistatout==1){
		intsdata->gradref=gradcalc;
	    }
	    intsdata->opd=opd;
	    intsdata->lltopd=lltopd;
	    intsdata->isim=isim;
	    CALL_THREAD(simu->wfs_ints[iwfs], 0);
	    dfree(lltopd);
	    intsdata->opd=0;
	    intsdata->lltopd=0;
	    if(psfout){
		cellarr_ccell(psfoutcellarr, isim, psfout);
		cellarr_dmat(ztiltoutcellarr, isim, *gradacc);
	    }
	}else{//Pywfs
	    pywfs_fft(&simu->ints->p[iwfs]->p[0], powfs[ipowfs].pywfs, opd);
	    dscale(simu->ints->p[iwfs]->p[0], parms->wfs[iwfs].siglevsim);
	}
    }
    TIM(2);
    if(dtrat_output){
	const double rne=parms->powfs[ipowfs].rne;
	const double bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
	if(do_phy){
	    /* In Physical optics mode, do integration and compute
	       gradients. The matched filter are in x/y coordinate even if
	       radpix=1. */
	    if(save_ints){
		cellarr_dcell(simu->save->intsnf[iwfs], isim, ints);
	    }
	    if(noisy){/*add noise */
		const double bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
		dmat **bkgrnd2=NULL;
		dmat **bkgrnd2c=NULL;
		if(powfs[ipowfs].bkgrnd){
		    if(powfs[ipowfs].bkgrnd->ny==1){
			bkgrnd2=powfs[ipowfs].bkgrnd->p;
		    }else{
			bkgrnd2=powfs[ipowfs].bkgrnd->p+nsa*wfsind;
		    }
		}
		if(powfs[ipowfs].bkgrndc){
		    if(powfs[ipowfs].bkgrndc->ny==1){
			bkgrnd2c=powfs[ipowfs].bkgrndc->p;
		    }else{
			bkgrnd2c=powfs[ipowfs].bkgrndc->p+nsa*wfsind;
		    }
		}
		for(int isa=0; isa<ints->nx; isa++){
		    dmat *bkgrnd2i=(bkgrnd2)?bkgrnd2[isa]:NULL;
		    dmat *bkgrnd2ic=(bkgrnd2c)?bkgrnd2c[isa]:NULL;
		    addnoise(ints->p[isa], &simu->wfs_rand[iwfs],
			     bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
		}
		if(save_ints){
		    cellarr_dcell(simu->save->intsny[iwfs], isim, ints);
		}
	    }
	    if(parms->powfs[ipowfs].dither && isim>=parms->powfs[ipowfs].dither_nskip
	       && parms->powfs[ipowfs].type==0 && parms->powfs[ipowfs].phytypesim2==1){
		/*Collect statistics with dithering*/
		DITHER_T *pd=simu->dither[iwfs];
		double cs, ss;
		dither_position(&cs, &ss, parms, ipowfs, isim, pd->deltam);
		//accumulate for matched filter
		dcelladd(&pd->imb, 1, ints, 1.);
		dcelladd(&pd->imx, 1, ints, cs);
		dcelladd(&pd->imy, 1, ints, ss);
	    }
	}
	if(do_phy){
	    if(parms->powfs[ipowfs].type==0){
		calc_phygrads(gradout, ints->p, parms, powfs, iwfs, parms->powfs[ipowfs].phytypesim);
	    }else{
		pywfs_grad(gradout, powfs[ipowfs].pywfs, ints->p[0]);
	    }
	}else{
	    /* geomtric optics accumulation mode. scale and copy results to output. */
	    dcp(gradout,*gradacc);
	    if(dtrat!=1){
		dscale(*gradout,1./dtrat);/*average */
	    }
	    if(noisy && !parms->powfs[ipowfs].usephy){
		const dmat *nea=powfs[ipowfs].neasim->p[wfsind];
		const double *neax=nea->p;
		const double *neay=nea->p+nsa;
		const double *neaxy=nea->p+nsa*2;
		double *restrict ggx=(*gradout)->p;
		double *restrict ggy=(*gradout)->p+nsa;
		for(int isa=0; isa<nsa; isa++){
		    /*Preserve the random sequence. */
		    double n1=randn(&simu->wfs_rand[iwfs]);
		    double n2=randn(&simu->wfs_rand[iwfs]);
		    double errx=neax[isa]*n1;
		    double erry=neay[isa]*n2+neaxy[isa]*n1;/*cross term. */
		    ggx[isa]+=errx;
		    ggy[isa]+=erry;
		}
	    }
	}
	if(save_gradgeom){
	    dmat *gradtmp=NULL;
	    dadd(&gradtmp, 1, *gradacc, 1./dtrat);
	    cellarr_dmat(simu->save->gradgeom[iwfs], isim, gradtmp);/*noise free. */
	    dfree(gradtmp);
	}
    }//dtrat_out
    dfree(gradcalc);
    TIM(3);
#if TIMING==1
    info("wfs %d grad timing: ray %.2f ints %.2f grad %.2f\n",iwfs,tk1-tk0,tk2-tk1,tk3-tk2);
#endif
}
static double calc_dither_amp(dmat *signal, /**<array of data. nmod*nsim */
			      long dtrat,   /**<skip columns due to wfs/sim dt ratio*/
			      long npoint,  /**<number of points during dithering*/
			      int detrend   /**<flag for detrending (remove linear signal)*/
    ){
    const long nmod=signal->nx;
    long nframe=(signal->ny-1)/dtrat+1;//number of wavefront frames
    double slope=0;
    long offset=(nframe/npoint-1)*npoint;//number of WFS frame separations between first and last cycle
    if(detrend && offset){//detrending
	for(long ip=0; ip<npoint; ip++){
	    for(long im=0; im<nmod; im++){
		long i0=ip*dtrat*nmod+im;
		long i1=(ip+offset)*dtrat*nmod+im;
		slope+=signal->p[i1]-signal->p[i0];
	    }
	}
	slope/=(npoint*nmod*offset);
	//info("slope=%g. npoint=%ld, nmod=%ld, nframe=%ld, offset=%ld\n", slope, npoint, nmod, nframe, offset);
    }
    double anglei=M_PI*2/npoint;
    double ipv=0, qdv=0;
    if(nmod==2){//tip and tilt
	for(int iframe=0; iframe<nframe; iframe++){
	    double angle=anglei*iframe;//position of dithering
	    double cs=cos(angle);
	    double ss=sin(angle);
	    double ttx=signal->p[iframe*dtrat*2]-slope*iframe;
	    double tty=signal->p[iframe*dtrat*2+1]-slope*iframe;
	    ipv+=(ttx*cs+tty*ss);
	    qdv+=(ttx*ss-tty*cs);
	}
    }else if(nmod==1){//single mode dithering
	for(int iframe=0; iframe<nframe; iframe++){
	    double angle=anglei*iframe;//position of dithering
	    double cs=cos(angle);
	    double ss=sin(angle);
	    double mod=signal->p[iframe*dtrat]-slope*iframe;
	    ipv+=(mod*cs);
	    qdv+=(mod*ss);
	}
    }
    double a2m=sqrt(ipv*ipv+qdv*qdv)/nframe;
    //info("a2m=%g mas\n", a2m * 206265000);
    //writebin(signal, "signal"); exit(0);
    return a2m;
}

/*Fast steering mirror for each WFS*/
void wfsgrad_fsm(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int isim=simu->isim;
    /*Uplink FSM*/
    const int index=parms->recon.glao?(ipowfs+ipowfs*parms->npowfs):(iwfs+iwfs*parms->nwfs);
    dmat *PTT=recon->PTT?(recon->PTT->p[index]):0;
    if(!PTT){
	error("powfs %d has FSM, but PTT is empty\n", ipowfs);
    }
    /* Compute FSM error. */
    simu->fsmerr=simu->fsmerr_store;
    dmm(&simu->fsmerr->p[iwfs], 0, PTT, simu->gradcl->p[iwfs], "nn", 1);
    //Save data
    simu->fsmerrs->p[iwfs]->p[isim*2  ]=simu->fsmerr->p[iwfs]->p[0];
    simu->fsmerrs->p[iwfs]->p[isim*2+1]=simu->fsmerr->p[iwfs]->p[1];
    simu->fsmcmds->p[iwfs]->p[isim*2  ]=simu->fsmreal->p[iwfs]->p[0];
    simu->fsmcmds->p[iwfs]->p[isim*2+1]=simu->fsmreal->p[iwfs]->p[1];
}

/*Postprocessing for dithering signal extraction*/
void wfsgrad_dither(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int isim=simu->isim;
    if(parms->powfs[ipowfs].dither && isim>=parms->powfs[ipowfs].dither_nskip){
	//Compute estimate of dithering signal per subaperture from WFS gradients
	DITHER_T *pd=simu->dither[iwfs];
	double cs, ss;
	dither_position(&cs, &ss, parms, ipowfs, isim, pd->deltam);
	if(parms->powfs[ipowfs].type==0 && parms->powfs[ipowfs].phytypesim2!=1){
	    //CoG use different gain for each.
	    const int nsa=powfs[ipowfs].saloc->nloc;
	    if(!pd->ggm){
		pd->ggm=dnew(nsa*2,1);
	    }
	    for(int isa=0; isa<nsa; isa++){
		pd->ggm->p[isa]+=cs*simu->gradcl->p[iwfs]->p[isa];
		pd->ggm->p[isa+nsa]+=ss*simu->gradcl->p[iwfs]->p[isa+nsa];
	    }
	}
    }
    
    if(parms->powfs[ipowfs].dither && isim>=parms->powfs[ipowfs].dither_pllskip){
	DITHER_T *pd=simu->dither[iwfs];
	/* 
	   The phase locked loop determines the phase of actual dithering signal
	   from WFS error signal.
	*/
	if(parms->powfs[ipowfs].dither==1){//dither in tip/tilt
	    double err, cd, sd;
	    //Average position of expected spot during integration. Uplink propagation is accounted for in LGS.
	    dither_position(&cd, &sd, parms, ipowfs, isim, pd->deltam);
	    //Correlate with error signal from WFS measurements.
	    err=(-sd*(simu->fsmerr->p[iwfs]->p[0])
		 +cd*(simu->fsmerr->p[iwfs]->p[1]))/(parms->powfs[ipowfs].dither_amp);
	    pd->delta+=parms->powfs[ipowfs].dither_gpll*err;
	}

	/* Determine the dither signal strength from sensor as well as measurements*/
	const int npll=parms->powfs[ipowfs].dither_npll;
	int npllacc=(simu->isim-parms->powfs[ipowfs].dither_pllskip+1)/parms->powfs[ipowfs].dtrat;
	//This only executes when PLL should have output use saved time history.
	if(npllacc>0 && npllacc%npll==0){
	    pd->deltam=pd->delta;
	    const int npoint=parms->powfs[ipowfs].dither_npoint;
	    int detrend=1;//parms->powfs[ipowfs].dither==1?0:1;//detrend not needed in tip/tilt mode
	    
	    int ncol=(npll-1)*parms->powfs[ipowfs].dtrat+1;
	    dmat *tmp=0;
	    tmp=drefcols(simu->fsmerrs->p[iwfs], simu->isim-ncol+1, ncol);
	    pd->a2me=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend);
	    dfree(tmp);
	    tmp=drefcols(simu->fsmcmds->p[iwfs], simu->isim-ncol+1, ncol);
	    pd->a2m=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend);
	    dfree(tmp);
	    
	    if(iwfs==parms->powfs[ipowfs].wfs->p[0]){
		const double anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
		info2("PLL step%d, wfs%d: deltam=%.2f frame, a2m=%.1f mas, a2me=%.1f mas\n",
		      isim, iwfs, pd->deltam/anglei, pd->a2m*206265000, pd->a2me*206265000);
	    }
	    if(simu->resdither){
		int ic=(npllacc-1)/(npll);
		simu->resdither->p[iwfs]->p[ic*2+0]=pd->deltam;
		simu->resdither->p[iwfs]->p[ic*2+1]=pd->a2m;
	    }
	}
    }//PLL loop

    if(parms->powfs[ipowfs].dither){
	DITHER_T *pd=simu->dither[iwfs];
	/* Update drift mode computation. Only useful when wfs t/t is removed*/
	int npll=parms->powfs[ipowfs].dither_npll;
	int npllacc=(simu->isim-parms->powfs[ipowfs].dither_nskip+1)/parms->powfs[ipowfs].dtrat;
	if(npllacc>0 && npllacc % npll==0){
	    double scale1=1./npll;
	    double scale2=scale1*2./(pd->a2m);
	    if(pd->imb){//matched filter
		dcellscale(pd->imb, scale1);
		dmat *ibgrad=0;
		calc_phygrads(&ibgrad, pd->imb->p, parms, powfs, iwfs, 2);
		if(parms->powfs[ipowfs].trs){//tip/tilt drift signal
		    dmat *tt=dnew(2,1);
		    const int index=parms->recon.glao?(ipowfs+ipowfs*parms->npowfs):(iwfs+iwfs*parms->nwfs);
		    dmat *PTT=recon->PTT?(recon->PTT->p[index]):0;
		    dmm(&tt, 0, PTT, ibgrad, "nn", 1);
		    simu->fsmerr->p[iwfs]->p[0]+=tt->p[0];
		    simu->fsmerr->p[iwfs]->p[1]+=tt->p[1];
		    dfree(tt);
		}
		//Smooth trombone movement by provide continuous err.
		if(parms->powfs[ipowfs].llt){
		    dmat *focus=dnew(1,1);
		    dmat *RFlgsg=recon->RFlgsg->p[parms->recon.glao
						  ?(ipowfs+ipowfs*parms->npowfs)
						  :(iwfs+iwfs*parms->nwfs)];
		    dmm(&focus, 0, RFlgsg, ibgrad, "nn", 1);
		    simu->zoomerr->p[iwfs]=focus->p[0];//2014-12-17: removed /npll;
		    dfree(focus);
		}
		dfree(ibgrad);
	    
		dcelladd(&pd->i0, 1, pd->imb, 1);//imb was already scaled
		dcelladd(&pd->gx, 1, pd->imx, scale2);
		dcelladd(&pd->gy, 1, pd->imy, scale2);
		dcellzero(pd->imb);
		dcellzero(pd->imx);
		dcellzero(pd->imy);
	    }else if(pd->ggm){//cog
		dadd(&pd->gg0, 1, pd->ggm, scale2);
		dzero(pd->ggm);
	    }
	}//Drift mode computation
    }

    if(!parms->powfs[ipowfs].trs){
	/*when WFS t/t is used for reconstruction, do now close FSM
	 * loop. Subtract actual dithering signal.*/
	if(parms->powfs[ipowfs].dither && powfs[ipowfs].gain){
	    double amp,cs,ss; 
	    if(0){
		amp=parms->powfs[ipowfs].dither_amp;
		const double anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
		const double angle=anglei*(isim/parms->powfs[ipowfs].dtrat);
		cs=cos(angle);
		ss=sin(angle);

	    }else{
		DITHER_T *pd=simu->dither[iwfs];
		dither_position(&cs, &ss, parms, ipowfs, isim, pd->deltam);
		amp=pd->a2me;
	    }
	    double ptt[2]={-cs*amp, -ss*amp};
	    dmulvec(simu->gradcl->p[iwfs]->p, recon->TT->p[iwfs+iwfs*parms->nwfsr], ptt, 1);
	    //pfsmerrs[isim][0]+=ptt[0];
	    //pfsmerrs[isim][1]+=ptt[1];
	}
	simu->fsmerr->p[iwfs]->p[0]=0;
	simu->fsmerr->p[iwfs]->p[1]=0;
    }
}

/**
   Accomplish Two tasks:
   
   1) High pass filter lgs focus to remove sodium range variation effact.
   2) Average LGS focus measurement to drive the trombone.
   
   We trust the focus measurement of the LGS WFS at high temporal frequency
   where NGS cannot provide due to low frame rate. After the HPF on lgs
   gradients, our system is NO LONGER affected by sodium layer variation.

   if sim.mffocus==1: The HPF is applied to each LGS WFS indepently. This largely
   removes the effect of differential focus. powfs.dfrs is no longer needed. (preferred).

   if sim.mffocus==2: We apply a LPF on the average focus from six LGS WFS, and
   then remove this value from all LGS WFS. The differential focus is still
   present and powfs.dfrs need to be set to 1 to handle it in tomography. This
   is the original focus tracking method.
*/
void wfsgrad_mffocus(SIM_T* simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    
    /*New plate mode focus offset for LGS WFS. Not needed*/
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].llt && parms->sim.ahstfocus==2 
	   && simu->Mint_lo && simu->Mint_lo->mint->p[1]
	   && (simu->isim+1)%parms->powfs[ipowfs].dtrat==0){
	    /*In new ahst mode, the first plate scale mode contains focus for
	      lgs. But it turns out to be not necessary to remove it because the
	      HPF in the LGS path removed the influence of this focus mode. set
	      sim.ahstfocus=2 to enable adjust gradients.*/
	    double scale=simu->recon->ngsmod->scale;
	    double focus=-simu->Mint_lo->mint->p[1]->p[0]->p[2]*(scale-1);
	    dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[ipowfs], focus);
	}
    }

    int hi_output=(!parms->sim.closeloop || (simu->isim+1)%parms->sim.dtrat_hi==0);
    if(hi_output){
	dcellzero(simu->LGSfocus);
	/*residual focus along ngs estimated from LGS measurement.*/
	dcellmm(&simu->LGSfocus, recon->RFlgsg, simu->gradcl,"nn",1);
	dcell *LGSfocus=simu->LGSfocus;
	long nwfsllt=0; 
	double lpfocusm=0;
	double lgsfocusm=0;
	if(simu->isim==parms->sim.start){
	    /*Here we set trombone position according to focus in the first
	     * measurement. And adjust the focus content of this * measurement. */
	    lgsfocusm=0;
	    if(parms->sim.zoomshare){
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		    if(!LGSfocus->p[iwfs]) continue;
		    lgsfocusm+=LGSfocus->p[iwfs]->p[0];
		    nwfsllt++;
		}
		lgsfocusm/=nwfsllt;
	    }
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		if(!LGSfocus->p[iwfs]) continue;
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->sim.zoomshare){
		    simu->zoomint->p[iwfs]=lgsfocusm;
		}else{
		    simu->zoomint->p[iwfs]=LGSfocus->p[iwfs]->p[0]; 
		}
		LGSfocus->p[iwfs]->p[0]-=simu->zoomint->p[iwfs];
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[ipowfs], -simu->zoomint->p[iwfs]);
		info2("wfs %d: Set trombone position to %g.\n", iwfs, simu->zoomint->p[iwfs]);
	    }
	    lgsfocusm=0;
	    nwfsllt=0;
	}
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    if(!LGSfocus->p[iwfs]) continue;
	    int ipowfs=parms->wfs[iwfs].powfs;
	    //In RTC. LPF can be put after using the value to put it off critical path.
	    double lpfocus=parms->sim.lpfocushi;
	    simu->lgsfocuslpf->p[iwfs]=simu->lgsfocuslpf->p[iwfs]*(1-lpfocus)+LGSfocus->p[iwfs]->p[0]*lpfocus;
	    if(parms->sim.mffocus==1){//remove LPF focus from each lgs
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[ipowfs], -simu->lgsfocuslpf->p[iwfs]);
	    }
	    lgsfocusm+=LGSfocus->p[iwfs]->p[0];
	    //Averaged LPF focus
	    lpfocusm+=simu->lgsfocuslpf->p[iwfs]; 
	    nwfsllt++;
	    
	}
	lgsfocusm/=nwfsllt;
	if(parms->sim.mffocus==2){//remove LPF GLOBAL focus from each lgs
	    lpfocusm/=nwfsllt;
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		if(!LGSfocus->p[iwfs]) continue;
		int ipowfs=parms->wfs[iwfs].powfs;
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[ipowfs], -lpfocusm);
	    }
	}
    
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    if(LGSfocus->p[iwfs]){
		if(parms->sim.zoomshare){
		    //all lgs share the same trombone so take the average value.
		    simu->zoomavg->p[iwfs]+=lgsfocusm;
		}else{
		    simu->zoomavg->p[iwfs]+=LGSfocus->p[iwfs]->p[0];
		}
	    }
	}
	/*zoom error is zero order hold even if no output from averager*/
	if((simu->reconisim+1)%parms->sim.zoomdtrat==0){
	    int dtrat=parms->sim.zoomdtrat;
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].llt && (!parms->powfs[ipowfs].dither || parms->powfs[ipowfs].phytypesim!=1)){
		    //For those WFS that dither and run mtch, use focus error from ib instead
		    simu->zoomerr->p[iwfs]=simu->zoomavg->p[iwfs]/dtrat;
		}
	    }
	    dzero(simu->zoomavg);
	}
    }
}

/**
   Every operation here should be in the Simulator not the Controller 
*/
void wfsgrad_post(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const PARMS_T *parms=simu->parms;
#if USE_CUDA
    if(parms->gpu.wfs){
	gpu_wfsgrad_sync(info);
    }
#endif
    //Postprocessing gradients
    const int isim=simu->isim;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const int dtrat=parms->powfs[ipowfs].dtrat;
	const int dtrat_output=((isim+1-parms->powfs[ipowfs].step)%dtrat==0);
	const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
	dmat **gradout=&simu->gradcl->p[iwfs];
	if(dtrat_output){
	    if(simu->powfs[ipowfs].gain && simu->powfs[ipowfs].gain->p[wfsind]){
		dcwm(*gradout, simu->powfs[ipowfs].gain->p[wfsind]);
	    }
	    if(do_phy){
		if(simu->fsmerr_store->p[iwfs]){
		    wfsgrad_fsm(simu, iwfs);
		}
		if(parms->powfs[ipowfs].dither){
		    wfsgrad_dither(simu, iwfs);
		}
	    }
	    //Gradient offset due to mainly NCPA calibration
	    if(simu->powfs[ipowfs].gradoff){
		dadd(gradout, 1, simu->powfs[ipowfs].gradoff->p[wfsind], -1);
	    }
	    if(parms->save.grad->p[iwfs]){
		cellarr_dmat(simu->save->gradcl[iwfs], isim, simu->gradcl->p[iwfs]);
	    }
	    if(parms->plot.run){
		drawopd("Gclx", simu->powfs[ipowfs].saloc, simu->gradcl->p[iwfs]->p, NULL,
			"WFS Closeloop Gradients (x)","x (m)", "y (m)",
			"x %d",  iwfs);
		drawopd("Gcly", simu->powfs[ipowfs].saloc, simu->gradcl->p[iwfs]->p+
			simu->powfs[ipowfs].saloc->nloc, NULL,
			"WFS Closeloop Gradients (y)","x (m)", "y (m)",
			"y %d",  iwfs);
	    }
	}
    }//for iwfs
}

/**
   Dither update: zoom corrector, matched filter, TWFS
*/
static void dither_update(SIM_T *simu){
    POWFS_T *powfs=simu->powfs;
    const PARMS_T *parms=simu->parms;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!parms->powfs[ipowfs].dither) continue;
	if((simu->isim+1)%parms->powfs[ipowfs].dtrat!=0) continue;
	const int nacc=(simu->isim-parms->powfs[ipowfs].dither_nskip+1)/parms->powfs[ipowfs].dtrat;
	const int nwfs=parms->powfs[ipowfs].nwfs;
	const int npll=parms->powfs[ipowfs].dither_npll;
	
	if(parms->sim.zoomshare && parms->powfs[ipowfs].llt //this is LLT
	   && nacc>0 && nacc % npll==0){//There is drift mode computation
	    double sum=0;
	    for(int jwfs=0; jwfs<nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		sum+=simu->zoomerr->p[iwfs];
	    }
	    sum/=(double)nwfs;
	    for(int jwfs=0; jwfs<nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		simu->zoomerr->p[iwfs]=sum;
	    }
	}
	int ngrad=parms->powfs[ipowfs].dither_ngrad;
	if(nacc>0 && nacc % ngrad==0){//This is matched filter or cog update
	    const int nsa=powfs[ipowfs].saloc->nloc;
	    double scale1=(double)parms->powfs[ipowfs].dither_npll/(double)parms->powfs[ipowfs].dither_ngrad;
	    if(parms->powfs[ipowfs].phytypesim2==1 && parms->powfs[ipowfs].type==0){
		info2("Step %d: Update matched filter for powfs %d\n", simu->isim, ipowfs);
		//For matched filter
		if(!powfs[ipowfs].intstat){
		    powfs[ipowfs].intstat=calloc(1, sizeof(INTSTAT_T));
		}
		parms->powfs[ipowfs].radgx=0;
		if(!powfs[ipowfs].intstat->i0 || powfs[ipowfs].intstat->i0->ny!=nwfs){
		    dcellfree(powfs[ipowfs].intstat->i0);
		    dcellfree(powfs[ipowfs].intstat->gx);
		    dcellfree(powfs[ipowfs].intstat->gy);

		    powfs[ipowfs].intstat->i0=dcellnew(nsa, nwfs);
		    powfs[ipowfs].intstat->gx=dcellnew(nsa, nwfs);
		    powfs[ipowfs].intstat->gy=dcellnew(nsa, nwfs);
		}
		if(!powfs[ipowfs].gradoff){
		    powfs[ipowfs].gradoff=dcellnew(nwfs, 1);
		}

		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    //End of accumulation
		    DITHER_T *pd=simu->dither[iwfs];
		    //Scale the output due to accumulation
		    for(int isa=0; isa<nsa; isa++){
			dadd(powfs[ipowfs].intstat->i0->p+isa+jwfs*nsa, 0, pd->i0->p[isa], scale1);
			dadd(powfs[ipowfs].intstat->gx->p+isa+jwfs*nsa, 0, pd->gx->p[isa], scale1);
			dadd(powfs[ipowfs].intstat->gy->p+isa+jwfs*nsa, 0, pd->gy->p[isa], scale1);
		    }
		    dcellzero(pd->i0);
		    dcellzero(pd->gx);
		    dcellzero(pd->gy);
		    //Compute the gradient of i0 using old gradient algorithm and subtract from the gradient offset.
		    dmat *goff=0;
		    calc_phygrads(&goff, powfs[ipowfs].intstat->i0->p+jwfs*nsa, parms, powfs, iwfs, parms->powfs[ipowfs].phytypesim);
		    dadd(&powfs[ipowfs].gradoff->p[jwfs], 1, goff, -1);
		    if(parms->dbg.i0drift){
			dzero(goff);
			//Compute CoG of i0 + goff and drive it toward gradncpa with low gain (0.1)
			calc_phygrads(&goff, powfs[ipowfs].intstat->i0->p+jwfs*nsa, parms, powfs, iwfs, 2);
			dadd(&goff, 1, powfs[ipowfs].gradoff->p[jwfs], 1);
			if(powfs[ipowfs].gradncpa){
			    dadd(&goff, 1, powfs[ipowfs].gradncpa->p[jwfs], -1);
			}
			dadd(&powfs[ipowfs].gradoff->p[jwfs], 1, goff, -0.1);
		    }
		    dfree(goff);
		}
		if(parms->save.dither){
		    writebin(powfs[ipowfs].gradoff, "powfs%d_gradoff_adjusted_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].intstat->i0, "powfs%d_i0_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].intstat->gx, "powfs%d_gx_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].intstat->gy, "powfs%d_gy_%d", ipowfs, simu->isim);
		}
	    }else{
		//For CoG gain
		if(!powfs[ipowfs].gain){
		    powfs[ipowfs].gain=cellnew(nwfs, 1);
		}
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    //End of accumulation
		    DITHER_T *pd=simu->dither[iwfs];
		    const int ng=powfs[ipowfs].saloc->nloc*2;
		    if(!powfs[ipowfs].gain->p[jwfs]){
			powfs[ipowfs].gain->p[jwfs]=dnew(ng,1);
			dset(powfs[ipowfs].gain->p[jwfs], 1);
		    }
		    double mgold=dsum(powfs[ipowfs].gain->p[jwfs])/ng;
		    //gg0 is output/input of dither signal.
		    if(!pd->gg0){//single gain for all subapertures.
			double adj=parms->powfs[ipowfs].dither_gcog*(1-pd->a2me/pd->a2m);
			dadds(powfs[ipowfs].gain->p[jwfs], adj);
		    }else{//separate gain for each gradient
			dscale(pd->gg0, scale1);
			for(long ig=0; ig<ng; ig++){
			    double adj=parms->powfs[ipowfs].dither_gcog*(1.-pd->gg0->p[ig]);
			    powfs[ipowfs].gain->p[jwfs]->p[ig]+=adj;
			}
		    }
		    double mgain=dsum(powfs[ipowfs].gain->p[jwfs])/ng;
		    info2("Step %d: Update CoG gain for wfs %d. Averaged gain adjustment is %g\n", 
			  simu->isim, iwfs, mgain);
			     
		    //adjust WFS measurement dither signal by gain adjustment.
		    pd->a2me*=(mgain/mgold);

		    if(parms->save.dither){
			writebin(powfs[ipowfs].gain, "powfs%d_gain_%d", ipowfs, simu->isim);
		    }
		    dzero(pd->gg0);
		}
	    }
	    if(parms->powfs[ipowfs].phytypesim != parms->powfs[ipowfs].phytypesim2){
		parms->powfs[ipowfs].phytypesim=parms->powfs[ipowfs].phytypesim2;
		parms->powfs[ipowfs].phytype=parms->powfs[ipowfs].phytypesim;
		info2("Step %d: powfs %d changed to %s\n", simu->isim, ipowfs, 
		      parms->powfs[ipowfs].phytypesim==1?"matched filter":"CoG");
	    }
	    if(parms->powfs[ipowfs].phytypesim==1){//Matched filter
		if(parms->powfs[ipowfs].neareconfile || parms->powfs[ipowfs].phyusenea){
		    warning("Disable neareconfile and phyusenea\n");
		    parms->powfs[ipowfs].neareconfile=NULL;
		    parms->powfs[ipowfs].phyusenea=0;
		}
		parms->powfs[ipowfs].phytype=1;//Make sure MF is used for reconstruction.
		genmtch(parms, powfs, ipowfs);
		if(parms->save.dither){
		    writebin(powfs[ipowfs].intstat->mtche, "powfs%d_mtche_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].intstat->i0sum, "powfs%d_i0sum_%d", ipowfs, simu->isim);
		}
#if USE_CUDA
		if(parms->gpu.wfs){
		    info2("Update matched filter in GPU\n");
		    gpu_wfsgrad_update_mtche(parms, powfs);
		}
#endif
	    }
	    if(parms->sim.epdm->p[0]<=0){
		parms->sim.epdm->p[0]=0.5;
		warning("set ephi to 0.5\n");
	    }
	    int UPDATE_TOMO=1;
	    READ_ENV_INT(UPDATE_TOMO,0,1);
	    if(UPDATE_TOMO && parms->recon.alg==0){//no need to update LSR.
		setup_recon(simu->recon, parms, powfs, simu->aper);
#if USE_CUDA
		if(!parms->sim.evlol && (parms->gpu.tomo || parms->gpu.fit)){
		    gpu_update_recon(parms, powfs, simu->recon);
		}
#endif
	    }

	}
	int itpowfs=parms->itpowfs;
	if(itpowfs!=-1){
	    int ntacc=(simu->isim-parms->powfs[itpowfs].step+1);
	    if(ntacc>0 && ntacc%parms->powfs[itpowfs].dtrat==0){
		info2("Step %d: TWFS has output\n", simu->isim);
		dcell *Rmod=0;
		//Build radial mode error using closed loop TWFS measurements from this time step.
		dcellmm(&Rmod, simu->recon->RRtwfs, simu->gradcl, "nn", 1);
		//writebin(simu->gradcl->p[parms->nwfsr-1], "twfs_gcl_%d", simu->isim);
		//writebin(Rmod, "twfs_rmod_%d", simu->isim);
		if(!powfs[ipowfs].gradoff){
		    powfs[ipowfs].gradoff=dcellnew(nwfs, 1);
		}
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    dmm(&powfs[ipowfs].gradoff->p[jwfs], 1, simu->recon->GRall->p[ipowfs], Rmod->p[0], "nn", -1);
		    if(parms->plot.run){
			const int nsa=powfs[ipowfs].saloc->nloc;
			drawopd("Goffx",(loc_t*)powfs[ipowfs].pts, powfs[ipowfs].gradoff->p[jwfs]->p,NULL,
				"WFS Offset (x)","x (m)", "y (m)", "x %d",  iwfs);
			drawopd("Goffy",(loc_t*)powfs[ipowfs].pts, powfs[ipowfs].gradoff->p[jwfs]->p+nsa, NULL,
				"WFS Offset (y)","x (m)", "y (m)", "y %d",  iwfs);
		    }
		}
		if(parms->save.dither){
		    writebin(powfs[ipowfs].gradoff, "powfs%d_gradoff_%d", ipowfs, simu->isim);
		}
		dcellfree(Rmod);
	    }
	}
    }
}
/**
   Calls wfsgrad_iwfs() to computes WFS gradient in parallel.
   It also includes operations on Gradients before tomography.
*/
void wfsgrad(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    if(parms->sim.idealfit || parms->sim.evlol) return;
    // call the task in parallel and wait for them to finish. It may be done in CPU or GPU.
    extern int PARALLEL;
    if(!PARALLEL || parms->tomo.ahst_idealngs || !parms->gpu.wfs){
	CALL_THREAD(simu->wfs_grad_pre, 0);
    }
    CALL_THREAD(simu->wfs_grad_post, 0);
    dither_update(simu);
    if(parms->sim.mffocus){
	//high pass filter lgs focus to remove sodium range variation effect
	wfsgrad_mffocus(simu);
    }
    if(1+simu->isim==parms->sim.end){
#if USE_CUDA
	if(parms->gpu.wfs){
	    gpu_save_gradstat(simu);
	}else
#endif
	    save_gradstat(simu);
    }
    simu->tk_wfs=myclockd()-tk_start;
}
