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
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "skysim_utils.h"
/**
   \file skyc/skysim_utils.c
   Utilities for skysim.c
*/

/**
   Compute Open loop NGS mode wavefront error from mode vectors.  */
real calc_rms(const dmat *mod, const dmat *mcc, int istep0){
    real rms=0;
    for(long istep=istep0; istep<mod->ny; istep++){
	rms+=dwdot(PCOL(mod,istep), mcc, PCOL(mod,istep));
    }
    return rms/(mod->ny-istep0);
}

/**
   convert mod vector to ngs WFS cloc and add to it's opd or complex pupil function.
   Notice the the coordinate for each subaperture is different for TTF.
*/
void ngsmod2wvf(cmat *wvf,            /**<[in/out] complex pupil function*/
		real wvl,           /**<[in] the wavelength*/
		const dmat *modm,     /**<[in] the NGS mode vector*/
		const POWFS_S *powfs,       /**<[in] the powfs configuration*/
		int isa,              /**<[in] index of subaperture*/
		real thetax,        /**<[in] direction of WFS*/
		real thetay,        /**<[in] direction of WFS*/
		const PARMS_S *parms  /**<[in] the parms*/
    ){
    const real *mod=modm->p;
    const comp ik=COMPLEX(0, 2*M_PI/wvl);
    real dx=powfs->dxwvf;
    real ox=powfs->saloc->locx[isa]+dx*0.5;
    real oy=powfs->saloc->locy[isa]+dx*0.5;
    int nx=powfs->nxwvf;
    assert(wvf->nx==wvf->ny);
    comp *p=wvf->p+(wvf->nx-nx)/2*(1+wvf->nx);
    if(modm->nx==2){
	for(int iy=0; iy<nx; iy++){
	    real ym=(oy+iy*dx)*mod[1];
	    for(int ix=0; ix<nx; ix++){
		real x=ox+ix*dx;
		real tmp=x*mod[0]+ym;
		p[ix+wvf->nx*iy]*=cexp(ik*tmp);
	    }
	}
    }else{
	const real hc=parms->maos.hc;
	const real hs=parms->maos.hs;
	const real scale=pow(1.-hc/hs, -2);
	const real scale1=1.-scale;
	real focus;
	if(modm->nx>5){
	    focus=mod[5];
	    if(!parms->maos.ahstfocus){
		focus+=mod[2]*scale1;
	    }
	}else{
	    focus=mod[2]*scale1;
	}
	for(int iy=0; iy<nx; iy++){
	    real y=oy+iy*dx;
	    for(int ix=0; ix<nx; ix++){
		real x=ox+ix*dx;
		real xy=x*y;
		real x2=x*x;
		real y2=y*y;
		real tmp= 
		    +x*mod[0]
		    +y*mod[1]
		    +focus*(x2+y2)
		    +mod[2]*(-2*scale*hc*(thetax*x+thetay*y))
		    +mod[3]*((x2-y2)*scale1 - 2*scale*hc*(thetax*x-thetay*y))
		    +mod[4]*(xy*scale1-scale*hc*(thetay*x+thetax*y));
		p[ix+wvf->nx*iy]*=cexp(ik*tmp);
	    }
	}
    }
}


/**
   Time domain physical simulation.
   
   noisy: 
   - 0: no noise at all; 
   - 1: poisson and read out noise. 
   - 2: only poisson noise.   
*/
dmat *skysim_sim(dmat **mresout, const dmat *mideal, const dmat *mideal_oa, real ngsol, 
		 ASTER_S *aster, const POWFS_S *powfs, 
		 const PARMS_S *parms, int idtratc, int noisy, int phystart){
    int dtratc=0;
    if(!parms->skyc.multirate){
	dtratc=parms->skyc.dtrats->p[idtratc];
    }
    int hasphy;
    if(phystart>-1 && phystart<aster->nstep){
	hasphy=1;
    }else{
	hasphy=0;
    }
    const int nmod=mideal->nx;
    dmat *res=dnew(6,1);/*Results. 1-2: NGS and TT modes., 
			  3-4:On axis NGS and TT modes,
			  4-6: On axis NGS and TT wihtout considering un-orthogonality.*/
    dmat *mreal=NULL;/*modal correction at this step. */
    dmat *merr=dnew(nmod,1);/*modal error */
    dcell *merrm=dcellnew(1,1); merrm->p[0]=dnew(nmod, 1);
    dcell *pmerrm=NULL;
    const int nstep=aster->nstep?aster->nstep:parms->maos.nstep;
    dmat *mres=dnew(nmod,nstep);
    dmat* rnefs=parms->skyc.rnefs;
    dcell *zgradc=dcellnew3(aster->nwfs, 1, aster->ngs, 0);
    dcell *gradout=dcellnew3(aster->nwfs, 1, aster->ngs, 0);
    dmat *gradsave=0;
    if(parms->skyc.dbg){
	gradsave=dnew(aster->tsa*2,nstep);
    }
       
    SERVO_T *st2t=0;
    kalman_t *kalman=0;
    dcell *mpsol=0;
    int multirate=parms->skyc.multirate;
    dmat *moffset=multirate?dnew(nmod,1):0;
    lmat *dtrats=aster->dtrats;//only in multirate case. dtrat of each wfs in the asterism
    if(parms->skyc.servo>0){
	if(multirate){//only supports integrator
	    dmat *gtmp=dnew(1,1); gtmp->p[0]=1;
	    st2t=servo_new(merrm, NULL, 0, parms->maos.dt*dtratc, gtmp);
	    dfree(gtmp);
	}else{
	    //dmat *gtmp=dnew(1,1); gtmp->p[0]=0.5;
	    st2t=servo_new(merrm, NULL, 0, parms->maos.dt*dtratc,aster->gain->p[idtratc]);
	    //dfree(gtmp);
	}
    }else{
	if(multirate){
	    kalman=aster->kalman[0];
	}else{
	    kalman=aster->kalman[idtratc];
	}
    }
    if(kalman){
	kalman_init(kalman);
	mpsol=dcellnew(aster->nwfs, 1); //for psol grad.
    }
    const long nwvl=parms->maos.nwvl;
    dcell **psf=0, **mtche=0, **ints=0, **i0s=0;
    ccell *wvf=0,*wvfc=0, *otf=0;
    dmat *corr=0;
    if(hasphy){
	psf=mycalloc(aster->nwfs,dcell*);
	wvf=ccellnew(aster->nwfs,1);
	wvfc=ccellnew(aster->nwfs,1);
	mtche=mycalloc(aster->nwfs,dcell*);
	ints=mycalloc(aster->nwfs,dcell*);
	otf=ccellnew(aster->nwfs,1);
	i0s=mycalloc(aster->nwfs, dcell*);
	for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
	    const int ipowfs=aster->wfs[iwfs].ipowfs;
	    const long ncomp=parms->maos.ncomp[ipowfs];
	    const long nsa=parms->maos.nsa[ipowfs];
	    wvf->p[iwfs]=cnew(ncomp,ncomp);
	    wvfc->p[iwfs]=NULL;
	    psf[iwfs]=dcellnew(nsa,nwvl);
	    //cfft2plan(wvf->p[iwfs], -1);
	    if(parms->skyc.multirate){
		mtche[iwfs]=aster->wfs[iwfs].pistat->mtche[(int)aster->idtrats->p[iwfs]];
	    }else{
		mtche[iwfs]=aster->wfs[iwfs].pistat->mtche[idtratc];
	    }
	    i0s[iwfs]=aster->wfs[iwfs].pistat->i0s;
	    otf->p[iwfs]=cnew(ncomp,ncomp);
	    //cfft2plan(otf->p[iwfs],-1);
	    //cfft2plan(otf->p[iwfs],1);
	    ints[iwfs]=dcellnew(nsa,1);
	    int pixpsa=parms->skyc.pixpsa[ipowfs];
	    for(long isa=0; isa<nsa; isa++){
		ints[iwfs]->p[isa]=dnew(pixpsa,pixpsa);
	    }
	}
    }
    zfarr *zfmerr=0;
    if(parms->skyc.dbg){
	int dtrati=(multirate?dtrats->p[0]:dtratc);
	zfmerr=zfarr_init(nstep, 1, "%s/skysim_merr_aster%d_dtrat%d", dirsetup, aster->iaster, dtrati);
    }
    for(int irep=0; irep<parms->skyc.navg; irep++){
	if(kalman){
	    kalman_init(kalman);
	}else{
	    servo_reset(st2t);
	}
	dcellzero(zgradc);
	dcellzero(gradout);
	if(ints){
	    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		dcellzero(ints[iwfs]);
	    }
	}
	int plotted=0;
	for(int istep=0; istep<nstep; istep++){
	    memcpy(merr->p, PCOL(mideal,istep), nmod*sizeof(real));
	    dadd(&merr, 1, mreal, -1);/*form NGS mode error; */
	    memcpy(PCOL(mres,istep),merr->p,sizeof(real)*nmod);
	    if(mpsol){//collect averaged modes for PSOL.
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    dadd(&mpsol->p[iwfs], 1, mreal, 1);
		}
	    }
	    pmerrm=0;
	    if(istep>=parms->skyc.evlstart){/*performance evaluation*/
		real res_ngs=dwdot(merr->p,parms->maos.mcc,merr->p);
		if(res_ngs>ngsol*100){
		    dfree(res); res=NULL;
		    break;
		}
		{
		    res->p[0]+=res_ngs;
		    res->p[1]+=dwdot2(merr->p,parms->maos.mcc_tt,merr->p);
		    real dot_oa=dwdot(merr->p, parms->maos.mcc_oa, merr->p);
		    real dot_res_ideal=dwdot(merr->p, parms->maos.mcc_oa, PCOL(mideal,istep));
		    real dot_res_oa=0;
		    for(int imod=0; imod<nmod; imod++){
			dot_res_oa+=merr->p[imod]*P(mideal_oa,imod,istep);
		    }
		    res->p[2]+=dot_oa-2*dot_res_ideal+2*dot_res_oa;
		    res->p[4]+=dot_oa;
		}
		{
		    real dot_oa_tt=dwdot2(merr->p, parms->maos.mcc_oa_tt, merr->p);
		    /*Notice that mcc_oa_tt2 is 2x5 marix. */
		    real dot_res_ideal_tt=dwdot(merr->p, parms->maos.mcc_oa_tt2, PCOL(mideal,istep));
		    real dot_res_oa_tt=0;
		    for(int imod=0; imod<2; imod++){
			dot_res_oa_tt+=merr->p[imod]*P(mideal_oa,imod,istep);
		    }
		    res->p[3]+=dot_oa_tt-2*dot_res_ideal_tt+2*dot_res_oa_tt;
		    res->p[5]+=dot_oa_tt;
		}
	    }//if evl

	    if(istep<phystart || phystart<0){
		/*Ztilt, noise free simulation for acquisition. */
		dmm(&zgradc->m, 1, aster->gm, merr, "nn", 1);/*grad due to residual NGS mode. */
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		    const int ipowfs=aster->wfs[iwfs].ipowfs;
		    const long ng=parms->maos.nsa[ipowfs]*2;
		    for(long ig=0; ig<ng; ig++){
			zgradc->p[iwfs]->p[ig]+=aster->wfs[iwfs].ztiltout->p[istep*ng+ig];
		    }
		}
	
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		    int dtrati=(multirate?dtrats->p[iwfs]:dtratc);
		    if((istep+1) % dtrati==0){
			dadd(&gradout->p[iwfs], 0, zgradc->p[iwfs], 1./dtrati);
			dzero(zgradc->p[iwfs]);
			if(noisy){
			    int idtrati=(multirate?aster->idtrats->p[iwfs]:idtratc);
			    dmat *nea=aster->wfs[iwfs].pistat->sanea->p[idtrati];
			    for(int i=0; i<nea->nx; i++){
				gradout->p[iwfs]->p[i]+=nea->p[i]*randn(&aster->rand);
			    }
			}
			pmerrm=merrm;//record output.
		    }
		}
	    }else{
		/*Accumulate PSF intensities*/
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    const real thetax=aster->wfs[iwfs].thetax;
		    const real thetay=aster->wfs[iwfs].thetay;
		    const int ipowfs=aster->wfs[iwfs].ipowfs;
		    const long nsa=parms->maos.nsa[ipowfs];
		    ccell* wvfout=aster->wfs[iwfs].wvfout[istep];
		    for(long iwvl=0; iwvl<nwvl; iwvl++){
			real wvl=parms->maos.wvl[iwvl];
			for(long isa=0; isa<nsa; isa++){
			    ccp(&wvfc->p[iwfs], P(wvfout,isa,iwvl));
			    /*Apply NGS mode error to PSF. */
			    ngsmod2wvf(wvfc->p[iwfs], wvl, merr, powfs+ipowfs, isa,
				       thetax, thetay, parms);
			    cembedc(wvf->p[iwfs],wvfc->p[iwfs],0,C_FULL);
			    cfft2(wvf->p[iwfs],-1);
			    /*peak in corner. */
			    cabs22d(&psf[iwfs]->p[isa+nsa*iwvl], 1., wvf->p[iwfs], 1.);
			}/*isa */
		    }/*iwvl */
		}/*iwfs */
	
		/*Form detector image from accumulated PSFs*/
		real igrad[2];
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    int dtrati=dtratc, idtrat=idtratc;
		    if(multirate){//multirate
			idtrat=aster->idtrats->p[iwfs];
			dtrati=dtrats->p[iwfs];
		    }
		    if((istep+1) % dtrati == 0){/*has output */
			dcellzero(ints[iwfs]);
			const int ipowfs=aster->wfs[iwfs].ipowfs;
			const long nsa=parms->maos.nsa[ipowfs];
			for(long isa=0; isa<nsa; isa++){
			    for(long iwvl=0; iwvl<nwvl; iwvl++){
				real siglev=aster->wfs[iwfs].siglev->p[iwvl];
				ccpd(&otf->p[iwfs],psf[iwfs]->p[isa+nsa*iwvl]);
				cfft2i(otf->p[iwfs], 1); /*turn to OTF, peak in corner */
				ccwm(otf->p[iwfs], powfs[ipowfs].dtf[iwvl].nominal);
				cfft2(otf->p[iwfs], -1);
				dspmulcreal(ints[iwfs]->p[isa]->p, powfs[ipowfs].dtf[iwvl].si, 
					   otf->p[iwfs]->p, siglev);
			    }
		
			    /*Add noise and apply matched filter. */
#if _OPENMP >= 200805 
#pragma omp critical 
#endif
			    switch(noisy){
			    case 0:/*no noise at all. */
				break;
			    case 1:/*both poisson and read out noise. */
				{
				    real bkgrnd=aster->wfs[iwfs].bkgrnd*dtrati;
				    addnoise(ints[iwfs]->p[isa], &aster->rand, bkgrnd, bkgrnd, 0,0,0,P(rnefs,idtrat,ipowfs), parms->skyc.excess);
				}
				break;
			    case 2:/*there is still poisson noise. */
				addnoise(ints[iwfs]->p[isa], &aster->rand, 0, 0, 0,0,0,0,parms->skyc.excess);
				break;
			    default:
				error("Invalid noisy\n");
			    }
		
			    igrad[0]=0;
			    igrad[1]=0;
			    real pixtheta=parms->skyc.pixtheta[ipowfs];
		
			    switch(parms->skyc.phytype){
			    case 1:
				dmulvec(igrad, mtche[iwfs]->p[isa], ints[iwfs]->p[isa]->p, 1);
				break;
			    case 2:
				dcog(igrad, ints[iwfs]->p[isa], 0, 0, 0, 3*P(rnefs,idtrat,ipowfs), 0);
				igrad[0]*=pixtheta;
				igrad[1]*=pixtheta;
				break;
			    case 3:
				dcorr(&corr, ints[iwfs]->p[isa], i0s[iwfs]->p[isa]);
				dpara3(igrad, corr);
				igrad[0]*=pixtheta;
				igrad[1]*=pixtheta;
				break;
			    default:
				error("Invalid phytype\n");
			    }
			    gradout->p[iwfs]->p[isa]=igrad[0];
			    gradout->p[iwfs]->p[isa+nsa]=igrad[1];
			}/*isa */
			pmerrm=merrm;
			dcellzero(psf[iwfs]);/*reset accumulation.*/
		    }/*if iwfs has output*/
		}/*for wfs*/
	    }/*if phystart */
	    //output to mreal after using it to ensure two cycle delay.
	    if(st2t){//Type I or II control.
		if(st2t->mint->p[0]){//has output.
		    dcp(&mreal, st2t->mint->p[0]->p[0]);
		}
	    }else{//LQG control
		kalman_output(kalman, &mreal, 0, 1);
	    }
	    if(parms->skyc.servo<0){//LQG control
		int indk=0;
		//Form PSOL grads and obtain index to LQG M
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		    int dtrati=(multirate?dtrats->p[iwfs]:dtratc);
		    if((istep+1) % dtrati==0){
			indk|=1<<iwfs;
			dmm(&gradout->p[iwfs], 1, aster->g->p[iwfs], mpsol->p[iwfs], "nn", 1./dtrati);
			dzero(mpsol->p[iwfs]);
		    }
		}
		if(indk){
		    kalman_update(kalman, gradout->m, indk-1);
		}
	    }else{
		if(pmerrm){
		    if(!multirate){//single rate
			dmm(&merrm->p[0], 0, aster->pgm->p[idtratc], gradout->m, "nn", 1);
		    }else{
			int indk=0;
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			    int dtrati=(multirate?dtrats->p[iwfs]:dtratc);
			    if((istep+1) % dtrati==0){
				indk|=1<<iwfs;
			    }
			}
			dzero(merrm->p[0]);
			if(indk==aster->pgm->nx){//slower loop
			    dmm(&moffset, 1, aster->pgm->p[indk-1], gradout->m, "nn", 0.5);
			    for(int imod=0; imod<nmod; imod++){
				if((aster->pgm->p[0] && imod==2)
				   || (aster->pgm->p[1] && imod>1)
				   || (!aster->pgm->p[0] && !aster->pgm->p[1]))
				{//direct output some modes
				    if(!plotted){
					//info("directly output %d\n", imod);
				    }
				    merrm->p[0]->p[imod]=moffset->p[imod];
				    moffset->p[imod]=0;
				}
			    }
			    plotted=1;
			}else{
			    dmm(&merrm->p[0], 1, aster->pgm->p[indk-1], gradout->m, "nn", 0.5);//aster->gain->p[indk-1]->p[0]);
			    dadd(&merrm->p[0], 1, moffset, 1);
			}
			if(zfmerr){
			    zfarr_push(zfmerr, istep, merrm->p[0]);
			}
		    }
		}
		servo_filter(st2t, pmerrm);//do even if merrm is zero. to simulate additional latency
	    }
	    if(parms->skyc.dbg){
		memcpy(PCOL(gradsave, istep), gradout->m->p, sizeof(real)*gradsave->nx);
	    }
	}/*istep; */
    }
    if(parms->skyc.dbg){
	int dtrati=(multirate?dtrats->p[0]:dtratc);
	writebin(gradsave,"%s/skysim_grads_aster%d_dtrat%d",dirsetup, aster->iaster,dtrati);
	writebin(mres,"%s/skysim_mres_aster%d_dtrat%d",dirsetup,aster->iaster,dtrati);
    }
    if(zfmerr){
	zfarr_close(zfmerr);
    }
    dfree(mreal);
    dcellfree(mpsol);
    dfree(merr);
    dcellfree(merrm);
    dcellfree(zgradc);
    dcellfree(gradout);
    dfree(gradsave);
    dfree(corr);
    if(hasphy){
	dcellfreearr(psf, aster->nwfs);
	dcellfreearr(ints, aster->nwfs);
        ccellfree(wvf);
	ccellfree(wvfc);
	ccellfree(otf);
	free(mtche);
	free(i0s);
    }
    servo_free(st2t);
    /*dfree(mres); */
    if(mresout) {
	*mresout=mres;
    }else{
	dfree(mres);
    }
    dscale(res, 1./((nstep-parms->skyc.evlstart)*parms->skyc.navg));
    return res;
}

/**
   Save NGS WFS and other information for later use in MAOS simulations.*/
void skysim_save(const SIM_S *simu, const ASTER_S *aster, const real *ipres, int selaster, int seldtrat, int isky){
    const PARMS_S* parms=simu->parms;
    const int nwvl=parms->maos.nwvl;
    char path[PATH_MAX-100];
    snprintf(path,sizeof(path),"Res%d_%d_maos/sky%d",simu->seed_maos,parms->skyc.seed,isky);
    mymkdir("%s",path);
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	dcell *sepsf=dcelldup(aster[selaster].wfs[iwfs].pistat->psf);
	for(int ic=0; ic<sepsf->nx*sepsf->ny; ic++){
	    dfftshift(sepsf->p[ic]);/*put peak in center. required by MAOS. */
	}
	writebin(sepsf, "%s/pistat_wfs%d",path,iwfs+6);
	dcellfree(sepsf);
	writebin(aster->wfs[iwfs].pistat->sanea->p[seldtrat], "%s/nea_tot_wfs%d",path,iwfs+6);
	writebin(aster[selaster].wfs[iwfs].pistat->sanea->p[seldtrat], 
	       "%s/nea_wfs%d",path,iwfs+6);
	writebin(aster[selaster].wfs[iwfs].pistat->sanea, 
		   "%s/neafull_wfs%d",path,iwfs+6);
    }
    if(parms->skyc.servo>0 && !parms->skyc.multirate){
	writebin(aster[selaster].gain->p[seldtrat], "%s/gain",path);
    }
    writebin(simu->mres->p[isky], "%s/mres",path);
    writebin(simu->psds,"%s/psds",path);
    char fnconf[PATH_MAX];
    snprintf(fnconf,sizeof(fnconf),"%s/base.conf",path);
    FILE *fp=fopen(fnconf,"w");

    fprintf(fp,"sim.seeds=[%d]\n",simu->seed_maos);
    fprintf(fp,"sim.end=%d\n", parms->maos.nstep);
    fprintf(fp,"sim.dt=%g\n", parms->maos.dt);
    fprintf(fp,"sim.zadeg=%g\n", parms->maos.zadeg);
    fprintf(fp,"sim.mffocus=%d\n", parms->maos.mffocus);
    fprintf(fp,"tomo.ahst_focus=%d\n", parms->maos.ahstfocus);
    fprintf(fp,"tomo.ahst_wt=3\n");
    if(parms->skyc.servo>0){
	fprintf(fp,"sim.eplo='gain.bin'\n");
    }
    fprintf(fp,"powfs0_llt.fnrange='%s'\n", parms->maos.fnrange);
    fprintf(fp,"atm.r0z=%.4f\n", parms->maos.r0z);
    fprintf(fp,"atm.size=[128 128]\n");
    if(parms->maos.wddeg){
	fprintf(fp, "atm.wddeg=[");
	for(int ips=0; ips<parms->maos.nwddeg; ips++){
	    fprintf(fp, "%.2f ", parms->maos.wddeg[ips]);
	}
	fprintf(fp, "]\n");
    }
    fprintf(fp,"wfs.thetax=[0 0  -33.287 -20.5725  20.5725 33.287");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp," %.4f", aster[selaster].wfs[iwfs].thetax*206265);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.thetay=[0 35 10.8156 -28.3156 -28.3156 10.8156");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp," %.4f", aster[selaster].wfs[iwfs].thetay*206265);
    }
    fprintf(fp,"]\n");

    fprintf(fp,"wfs.siglev=[900 900 900 900 900 900");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp, " %.2f", aster[selaster].wfs[iwfs].siglevtot);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.wvlwts=[1 1 1 1 1 1");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    fprintf(fp," %.2f ", aster[selaster].wfs[iwfs].siglev->p[iwvl]
		    /aster[selaster].wfs[iwfs].siglevtot);
	}
    }
    fprintf(fp,"]\n");

    if(parms->maos.npowfs<3){
	int nwfs[2]={0,0}; 
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	    nwfs[aster[selaster].wfs[iwfs].ipowfs]++;
	}
	fprintf(fp, "powfs.nwfs=[6");
	for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	    fprintf(fp, " %d", nwfs[ipowfs]);
	}
	fprintf(fp,"]\n");
    }else{
	error("Fill this out please\n");
    }
    dmat* rnefs=parms->skyc.rnefs;
    real rne=P(rnefs,seldtrat,0);
    real bkgrnd=aster[selaster].wfs[0].bkgrnd;

    if(parms->maos.npowfs==1){
	fprintf(fp, "powfs.nwfs=[6 1]\n");
	fprintf(fp, "powfs.piinfile=[\"\" \"pistat\" ]\n");
	fprintf(fp, "powfs.neareconfile=[\"\" \"nea_tot\"]\n");
	fprintf(fp, "powfs.phyusenea=[0 1]\n");
	fprintf(fp, "powfs.dtrat=[1 %d]\n", (int)parms->skyc.dtrats->p[seldtrat]);
	fprintf(fp, "powfs.bkgrnd=[0 %.2f]\n", bkgrnd);
	fprintf(fp, "powfs.rne=[3 %.2f]\n", rne);
	fprintf(fp, "powfs.phystep=[0 %ld]\n", 50+(long)parms->skyc.dtrats->p[seldtrat]*20);
	fprintf(fp, "powfs.noisy=[1 1 ]\n");
	fprintf(fp, "powfs.pixtheta=[0.8 %g]\n", parms->skyc.pixtheta[1]);
	fprintf(fp, "powfs.pixpsa=[6 %d]\n", parms->skyc.pixpsa[0]);
	fprintf(fp, "powfs.ncomp=[64 %d]\n", parms->maos.ncomp[0]);
	fprintf(fp, "powfs.nwvl=[1 %d]\n",nwvl);
	fprintf(fp, "powfs.wvl=[0.589e-6");
	for(int ip=0; ip<1; ip++){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		fprintf(fp, " %.4g", parms->maos.wvl[iwvl]);
	    }
	}
	fprintf(fp,"]\n");
    }else if(parms->maos.npowfs==2){
	fprintf(fp, "powfs.piinfile=[\"\" \"pistat\" \"pistat\"]\n");
	fprintf(fp, "powfs.neareconfile=[\"\" \"nea_tot\" \"nea_tot\"]\n");
	fprintf(fp, "powfs.phyusenea=[0 1 1]\n");
	fprintf(fp, "powfs.dtrat=[1 %d %d]\n", (int)parms->skyc.dtrats->p[seldtrat],
		(int)parms->skyc.dtrats->p[seldtrat]);
	fprintf(fp, "powfs.bkgrnd=[0 %.2f %.2f]\n", bkgrnd, bkgrnd);
	fprintf(fp, "powfs.rne=[3 %.2f %.2f]\n", rne,rne);
	fprintf(fp, "powfs.phystep=[0 %ld %ld]\n", 
		50+(long)parms->skyc.dtrats->p[seldtrat]*20, 
		50+(long)parms->skyc.dtrats->p[seldtrat]*20);
	fprintf(fp, "powfs.noisy=[1 1 1]\n");
	fprintf(fp, "powfs.pixtheta=[0.8 %g %g]\n",
		parms->skyc.pixtheta[0]*206265,
		parms->skyc.pixtheta[1]*206265);
	fprintf(fp, "powfs.pixpsa=[6 %d %d]\n",
		parms->skyc.pixpsa[0], 
		parms->skyc.pixpsa[1]);
	fprintf(fp, "powfs.ncomp=[64 %d %d]\n", 
		parms->maos.ncomp[0], parms->maos.ncomp[1]);
	fprintf(fp, "powfs.nwvl=[1 %d %d]\n",nwvl,nwvl);
	fprintf(fp, "powfs.wvl=[0.589e-6");
	for(int ip=0; ip<2; ip++){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		fprintf(fp, " %.4g", parms->maos.wvl[iwvl]);
	    }
	}
	fprintf(fp,"]\n");
	fprintf(fp, "powfs.wvlwts=[]\n");
    }else{
	error("Fill this out please\n");
    }
 
    fclose(fp);
    snprintf(fnconf,sizeof(fnconf),"%s/skyres.txt",path);
    fp=fopen(fnconf,"w");
    fprintf(fp, "TotAll\tNGS\tTT\n");
    fprintf(fp, "%g\t%g\t%g\n",
	    sqrt(ipres[0])*1e9, sqrt(ipres[1])*1e9, sqrt(ipres[2])*1e9);
    fclose(fp);
}
