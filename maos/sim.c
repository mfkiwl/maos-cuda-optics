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

/*
   Call various functions to do the simulation and evaluate performance.
*/
#include <search.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "common.h"
#include "recon.h"
#include "recon_utils.h"
#include "setup_powfs.h"
#include "setup_recon.h"
#include "sim.h"
#include "sim_utils.h"
#include "fdpcg.h"
#include "save.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
#if HAVE_NUMA_H
#include <numa.h>
#endif
extern int KEEP_MEM;
extern int draw_single;
static real tk_0;
static real tk_1;
static real tk_atm=0;
/**
   Initialize the simulation runtime data struct.
 */
SIM_T *maos_iseed(int iseed){
    if(iseed==0) tk_0=myclockd();
    tk_1=myclockd();
    const PARMS_T *parms=global->parms;
    POWFS_T *powfs=global->powfs;
    APER_T  *aper =global->aper;
    RECON_T *recon=global->recon;
    if(parms->fdlock && parms->fdlock->p[iseed]<0){
	warning("Another MAOS is already running. Skip seed %ld\n", 
		parms->sim.seeds->p[iseed]);
	return 0;
    }
    SIM_T *simu=init_simu(parms,powfs,aper,recon,iseed);
    global->simu=simu;
    if(!simu->pause){
	draw_single=1;//Only draw active frame.
    }else{
	draw_single=0;
    }
    global->iseed=iseed;

    if(parms->atm.frozenflow){
	genatm(simu);/*Generating atmospheric screen(s) that frozen flows.*/
	if(parms->tomo.predict){
	    if(recon->HXWtomo){
		setup_recon_HXW_predict(simu);
	    }
	    if(parms->tomo.precond==1){
		fdpcg_free(recon->fdpcg);
		recon->fdpcg=fdpcg_prepare(parms, recon, powfs, parms->tomo.predict?simu->atm:NULL);
	    }
	}
    }
#if USE_CUDA
    if(parms->gpu.evl || parms->gpu.wfs){
	/*Put here for the initial transfer to avoid messing up timing due to transfering. */
	gpu_atm2gpu(simu->atm, simu->atmscale, parms, iseed, parms->sim.start);/*takes 0.4s for NFIRAOS. */
	if(parms->tomo.predict){
	    gpu_update_recon_cn2(parms, recon);
	}
    }
#endif
    return simu;
}
/**
   Simulation for each time step. 

   Callable from matlab.
*/
void maos_isim(int isim){
    const PARMS_T *parms=global->parms;
    RECON_T *recon=global->recon;
    SIM_T   *simu =global->simu;
    int iseed=global->iseed;
    int simstart=parms->sim.start;
    int simend=parms->sim.end;
    long group=0;
    extern int NO_RECON, NO_WFS, NO_EVL;
    if(isim==simstart+1){//skip slow first step.
	tk_atm=myclockd();
    }
    if(isim+2+parms->sim.dtrat_hi>=simend){
	draw_single=0;
    }
    real ck_0=myclockd();
    simu->wfsisim=isim;
    simu->perfisim=isim;
    simu->status->isim=isim;
    if(!parms->sim.closeloop){
	simu->reconisim=isim;
    }else{//work on gradients from last time step for parallelization.
	simu->reconisim=isim-1;
    }
    sim_update_etf(simu);
    if(!parms->atm.frozenflow){
	//Do not put this one inside parallel 
	genatm(simu);
	/*re-seed the atmosphere in case atm is loaded from shm/file */
	seed_rand(simu->atm_rand, lrand(simu->init_rand));
    }
#if USE_CUDA
    if(parms->gpu.evl || parms->gpu.wfs){
	/*may need to copy another part */
	gpu_atm2gpu(simu->atm, simu->atmscale, parms, iseed, isim);
    }
#endif
    OMPTASK_SINGLE{
	if(parms->sim.dmproj){
	    /* temporarily disable FR.M so that Mfun is used.*/
	    cell *FRM=recon->fit->FR.M; recon->fit->FR.M=NULL; 
	    muv_solve(&simu->dmproj, &recon->fit->FL, &recon->fit->FR, NULL);
	    recon->fit->FR.M=FRM;/*set FR.M back*/
	    save_dmproj(simu);
	    if(!parms->fit.square){
		/* Embed DM commands to a square array for fast ray tracing */
		for(int idm=0; idm<parms->ndm; idm++){
		    loc_embed(simu->dmprojsq->p[idm], recon->aloc->p[idm], simu->dmproj->p[idm]->p);
		}
	    }
#if USE_CUDA
	    if(parms->gpu.evl || parms->gpu.wfs){
		gpu_dmproj2gpu(simu->dmprojsq);
	    }
#endif
	}
	if(PARALLEL==1){
	    simu->tk_0=myclockd();
	    /*
	      We do the big loop in parallel to make better use the
	      CPUs. Notice that the reconstructor is working on grad from
	      last time step so that there is no confliction in data access.

	      Launch perfevl_pre and wfsgrad_pre to allow gpu code to execute in advance. 
	    */

	    if(parms->gpu.evl && !NO_EVL){
		//Queue tasks on GPU, no stream sync is done
		QUEUE_THREAD(&group, simu->perfevl_pre, 0);
	    }
	    if(parms->tomo.ahst_idealngs!=1 && parms->gpu.wfs && !NO_WFS){
		//task for each wfs
		QUEUE_THREAD(&group, simu->wfsgrad_pre, 0);
	    }
	    if(!NO_RECON){
		//don't put this first. It has cpu overhead in computing gradol
		QUEUE(&group, (thread_wrapfun)reconstruct, simu, 1, 0);
	    }
	    if(!NO_EVL){
		if(parms->gpu.evl){
		    //wait for GPU tasks to be queued before calling sync
		    WAIT(group);
		}
		QUEUE(&group, (thread_wrapfun)perfevl, simu, 1, 0);
	    }
	    if(!NO_WFS){
		if(parms->tomo.ahst_idealngs==1 || (parms->gpu.wfs && !parms->gpu.evl)){
		    /*when we want to apply idealngs correction, wfsgrad need to wait for perfevl. */
		    WAIT(group);
		}
		QUEUE(&group, (thread_wrapfun)wfsgrad, simu, 1, 0);
	    }
	    if(!NO_RECON){
		//wait for all tasks to finish before modifying dmreal
		WAIT(group);
		shift_grad(simu);/*before filter() */
		filter_dm(simu);/*updates dmreal, so has to be after prefevl/wfsgrad is done. */
	    }
	    WAIT(group);
	}else{/*do the big loop in serial mode. */
	    if(parms->sim.closeloop){
		if(!NO_EVL) perfevl(simu);/*before wfsgrad so we can apply ideal NGS modes */
		if(!NO_WFS) wfsgrad(simu);/*output grads to gradcl, gradol */
		if(!NO_RECON) {
		    reconstruct(simu);/*uses grads from gradlast cl, gradlast ol. */
		    shift_grad(simu);
		    filter_dm(simu);
		}
	    }else{/*in OL mode,  */
		if(!NO_WFS) wfsgrad(simu);
		if(!NO_RECON) {
		    shift_grad(simu);
		    reconstruct(simu);
		    filter_dm(simu);
		}
		if(!NO_EVL) perfevl(simu);
	    }
	}
	if(simu->tomo_update){//This part causes random CUDA error in Geforce.
	    if(simu->tomo_update==1){//Only update cn2
		setup_recon_tomo_update(simu->recon, simu->parms);
	    }else{//Also update noise.
		setup_recon_update(simu->recon, simu->parms, simu->powfs);
#if USE_CUDA
		if(!parms->sim.evlol && (parms->gpu.tomo || parms->gpu.fit)){
		    gpu_update_recon(parms, simu->recon);
		}
#endif
	    }
	    simu->tomo_update=0;
	}
    }
    real ck_end=myclockd();
    long steps_done=iseed*(simend-simstart)+(isim+1-simstart);
    long steps_rest=parms->sim.nseed*(simend-simstart)-steps_done;
    if(isim==simstart){//first step, rough estimate.
	simu->status->mean=ck_end-ck_0;
	simu->status->rest=simu->status->mean*parms->sim.nseed*(simend-simstart);
    }else{
	simu->status->rest=(long)((ck_end-tk_0-(tk_atm-tk_1)*(iseed+1))/steps_done*steps_rest
				  +(tk_atm-tk_1)*(parms->sim.nseed-iseed-1));
	simu->status->mean=(ck_end-tk_atm)/(real)(isim-simstart);
    }
    simu->status->laps=(long)(ck_end-tk_0);
    simu->status->tot  =ck_end-ck_0;
   
    print_progress(simu);
}

/**
   Closed loop simulation main loop. 

   It calls init_simu() to initialize the simulation struct, and then calls
   maos_isim() for each simulation time step. Arranged this way so that
   maos_isim() can be called from matlab.
*/
void maos_sim(){
    const PARMS_T *parms=global->parms;
    RECON_T *recon=global->recon;
    int simend=parms->sim.end;
    int simstart=parms->sim.start;
   
    dbg("PARALLEL=%d\n", PARALLEL);
    real restot=0; long rescount=0;
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	SIM_T *simu=NULL;
	while(!(simu=maos_iseed(iseed))){
	    iseed++;
	}
	if(recon && recon->cn2est){//temporary. Should put runtime data in simu.
	    cn2est_reset(recon->cn2est);
	}
	    
#ifdef HAVE_NUMA_H
	numa_set_localalloc();
#endif
	if(PARALLEL==2){//event driven synchronization
#pragma omp parallel
#pragma omp sections
	    {
#pragma omp section
		for(int isim=simstart; isim<simend; isim++){
		    simu->perfisim=isim;
		    perfevl(simu);
		    print_progress(simu);
		}
#pragma omp section
		for(int isim=simstart; isim<simend; isim++){
		    simu->wfsisim=isim;
		    wfsgrad(simu);
		    shift_grad(simu);
		}
#pragma omp section
		for(int isim=simstart; isim<simend; isim++){
		    simu->reconisim=isim-1;
		    reconstruct(simu);
		    filter_dm(simu);
		}
	    }
	}else{
	    for(int isim=simstart; isim<simend; isim++){
		maos_isim(isim);
		if(simu->pause>0 && isim%simu->pause==0){
		    info("Press enter to step, c to resume:\n"); 
		    int key;
		    while((key=getchar())!=0x0a){
			if(key=='c'){
			    simu->pause=0;
			}
		    }
		    info("continuing...\n"); 
		}
	    }/*isim */
	}
	{
	    /*Compute average performance*/
	    int isim0;
	    if(parms->sim.closeloop){
		if(parms->sim.end>100){
		    isim0=MAX(50,parms->sim.end*0.2);
		}else{
		    isim0=MIN(20,parms->sim.end*0.5);
		}
	    }else{
		isim0=0;
	    }
	    real sum=0;
	    for(int i=isim0; i<parms->sim.end; i++){
		if(parms->sim.evlol){
		    sum+=simu->ole->p[i*parms->evl.nmod];
		}else{
		    sum+=simu->cle->p[i*parms->evl.nmod];
		}
	    }
	    restot+=sum/(parms->sim.end-isim0);
	    rescount++;
	}
	free_simu(simu);
	global->simu=NULL;
    }/*seed */
    info("Mean: %.2f\n", sqrt(restot/rescount)*1e9);
}
