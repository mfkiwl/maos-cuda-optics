/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "maos.h"
#include "setup_recon.h"
#include "mvm_client.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   Use the various algorithms recon.alg to assemble a final matrix to multiply
   to gradients to get DM commands.
 */
static void 
setup_recon_lsr_mvm(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    info2("Assembling LSR MVM in CPU\n");
    if(recon->LR.Mfun || parms->lsr.alg==1){
	/*
	  First create an identity matrix. then solve each column one by one. 
	*/
	const int ndm=parms->ndm;
	const int nwfs=parms->nwfsr;
	int ntotgrad=0;
	long *ngrad=recon->ngrad;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    ntotgrad+=powfs[ipowfs].saloc->nloc*2;
	}
	recon->MVM=dcellnew(ndm, nwfs);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(!parms->powfs[ipowfs].skip){
		for(int idm=0; idm<ndm; idm++){
		    recon->MVM->p[idm+ndm*iwfs]=dnew(recon->anloc[idm], powfs[ipowfs].saloc->nloc*2);
		}
	    }
	}

	dcell *res=NULL;
	int curg=0, curwfs=0;
	dmat *eye=dnew(ntotgrad, 1);
	dcell *eyec=d2cellref(eye, ngrad, nwfs);
	for(int ig=0; ig<ntotgrad; ig++){
	    if(!detached){
		info2("%6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", ig, ntotgrad);
	    }else if(ig%100==0){
		info2("%6d of %6d\n", ig, ntotgrad);
	    }
	    if(ig) eye->p[ig-1]=0; eye->p[ig]=1;
	    if(!parms->powfs[parms->wfsr[curwfs].powfs].skip){
		dcellzero(res);
		muv_solve(&res, &recon->LL, &recon->LR, eyec);
	    }
	    for(int idm=0; idm<ndm; idm++){
		dmat *to=recon->MVM->p[idm+curwfs*ndm];
		if(to){
		    int nact=to->nx;
		    memcpy(to->p+curg*nact,res->p[idm]->p, nact*sizeof(double));
		}
	    }
	    curg++;
	    if(curg>=ngrad[curwfs]){
		curwfs++;
		curg=0;
	    }
	}
	dcellfree(res);
	dcellfree(eyec);
	dfree(eye);
    }else{
	dcell *LR=NULL;
	spcellfull(&LR, recon->LR.M, 1);
	if(recon->LR.U && recon->LR.V){
	    dcellmm(&LR, recon->LR.U, recon->LR.V, "nt", -1);
	}
	muv_solve(&recon->MVM, &recon->LL, NULL,  LR);
	dcellfree(LR);
    }
}


typedef struct {
    const PARMS_T *parms;
    RECON_T *recon;
    dcell *MVMt;
    long (*curp)[2];
    long ntotact;
}MVR_MVM_T;
static void 
setup_recon_mvr_mvm_iact(thread_t *info){
    MVR_MVM_T *data=info->data;
    const PARMS_T *parms=data->parms;
    RECON_T *recon=data->recon;
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfsr;
    const long ntotact=data->ntotact;
    dcell *FLI=NULL;
    dcell *FRT=NULL;
    dcell *RLT=NULL;
    dcell *RRT=NULL;
    dmat *eye=dnew(ntotact, 1);
    dcell *eyec=d2cellref(eye, recon->anloc, ndm);
    long (*curp)[2]=data->curp;
    dcell *MVMt=data->MVMt;
    int nthread=recon->nthread;
    for(long iact=info->start; iact<info->end; iact++){
	int curdm=curp[iact][0];
	int curact=curp[iact][1];
	if(recon->actcpl && recon->actcpl->p[curdm]->p[curact]<EPS){
	    continue;
	}
	if(info->ithread==0){
	    if(!detached){
		info2("%6ld of %6ld\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", iact*nthread, ntotact);
	    }else if(iact % 100==0){
		info2("%6ld of %6ld\n", iact*nthread, ntotact);
	    }
	}
	//TIC;tic;
	dcellzero(FRT);
	dcellzero(RRT);
	/*Apply F_L*/
	eye->p[iact]=1;
	muv_solve(&FLI, &recon->FL, NULL, eyec);
	eye->p[iact]=0;
	/*Apply F_R'*/
	muv_trans(&FRT, &recon->FR, FLI, 1);
	//toc2("fit");
	/*Apply R_L*/
	dcellzero(RLT);//warm restart.
	muv_solve(&RLT, &recon->RL, NULL, FRT);
	/*Apply R_R'*/
	muv_trans(&RRT, &recon->RR, RLT, 1);
	//toc2("tomo");
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    dmat *to=MVMt->p[iwfs+curdm*nwfs];
	    if(to){
		int ng=to->nx;
		memcpy(to->p+curact*ng, RRT->p[iwfs]->p, ng*sizeof(double));
	    }
	}
	/*{
	    dcellwrite(FLI, "cpu_dmfit_%ld", iact);
	    dcellwrite(FRT, "cpu_opdx_%ld", iact);
	    dcellwrite(RLT, "cpu_opdr_%ld", iact);
	    dcellwrite(RRT, "cpu_grad_%ld", iact);
	    }*/
    }
    dcellfree(FLI);
    dcellfree(FRT);
    dcellfree(RLT);
    dcellfree(RRT);
    dfree(eye);
    dcellfree(eyec);
}


/**
   Use the various algorithms recon.alg to assemble a final matrix to multiply
   to gradients to get DM commands.
 */
static void 
setup_recon_mvr_mvm(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    info2("Assembling MVR MVM in CPU\n");
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfsr;
    long ntotact=0;
    for(int idm=0; idm<ndm; idm++){
	ntotact+=recon->anloc[idm];
    }
    long (*curp)[2]=malloc(ntotact*2*sizeof(long));
    int nact=0;
    for(int idm=0; idm<ndm; idm++){
	for(int iact=0; iact<recon->anloc[idm]; iact++){
	    curp[nact+iact][0]=idm;
	    curp[nact+iact][1]=iact;
	}
	nact+=recon->anloc[idm];
    }
    dcell *MVMt=dcellnew(nwfs, ndm);
    for(int idm=0; idm<ndm; idm++){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(!parms->powfs[ipowfs].skip){
		MVMt->p[iwfs+idm*nwfs]=dnew(powfs[ipowfs].saloc->nloc*2, recon->anloc[idm]);
	    }
	}
    }
    MVR_MVM_T data={parms, recon, MVMt, curp, ntotact};
    int nthread=recon->nthread;
    thread_t info[nthread];
    thread_prep(info, 0, ntotact, nthread, setup_recon_mvr_mvm_iact, &data);
    CALL_THREAD(info, nthread, 1);
    recon->MVM=dcelltrans(MVMt);
    free(curp);
    dcellfree(MVMt);
}

/*assemble matrix to do matrix vector multiply. Split from setup_recon because GPU may be used.*/
void setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    TIC;tic;
    if(parms->recon.mvm){
	if(parms->load.mvm){
	    recon->MVM=dcellread("%s", parms->load.mvm);
	    for(int idm=0; idm<parms->ndm; idm++){
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		    int ipowfs=parms->wfsr[iwfs].powfs;
		    if(!parms->powfs[ipowfs].skip){
			dmat *temp=recon->MVM->p[idm+iwfs*parms->ndm];
			if(temp->nx!=recon->aloc[idm]->nloc 
			   || temp->ny !=powfs[ipowfs].pts->nsa*2){
			    error("MVM[%d, %d] should have dimension (%ld, %ld), but have (%ld, %ld)\n", idm, iwfs, recon->aloc[idm]->nloc , powfs[ipowfs].pts->nsa*2, temp->nx, temp->ny);
			}
		    }
		}
	    }
	}
	if(parms->gpu.tomo && parms->gpu.fit){
#if USE_CUDA
	    gpu_setup_recon_mvm(parms, recon, powfs);
#endif
	}else if(!parms->load.mvm){
	    if(parms->load.mvmi){
		error("Not handled yet.\n");
	    }
	    if(parms->recon.alg==0){
		setup_recon_mvr_mvm(recon, parms, powfs);
	    }else{
		setup_recon_lsr_mvm(recon, parms, powfs);   
	    }
	}
	if(!parms->load.mvm && (parms->save.setup || parms->save.mvm)){
	    dcellwrite(recon->MVM, "MVM.bin");
	}
	if(parms->sim.mvmport){
	    mvm_client_send_m(parms, recon->MVM);
	}
	if(parms->gpu.tomo && parms->gpu.fit){
	    dcellfree(recon->MVM);
	}
    }
    toc2("setup_recon_mvm");
}