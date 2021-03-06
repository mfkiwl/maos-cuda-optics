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


#include <search.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include <fcntl.h>           /* For O_* constants */
#include <errno.h>
#include <getopt.h>
#include "common.h"
#include "mvm_client.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
#include "setup_powfs.h"
#include "genseotf.h"
/*
   A few utility routines
*/

/**
   Plot the loc, together with all beams
*/
void plotloc(const char *fig, const PARMS_T *parms, 
	     loc_t *loc, real ht, const char *format,...){
    format2fn;
    int ncir=parms->evl.nevl + parms->fit.nfit + parms->nwfs;
    if(parms->sim.ncpa_calib){
	ncir+=parms->sim.ncpa_ndir;
    }
    dmat *cir=dnew(4, ncir);
    int count=0;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	real hs=parms->evl.hs->p[ievl];
	P(cir,0,count)=ht*parms->evl.thetax->p[ievl];
	P(cir,1,count)=ht*parms->evl.thetay->p[ievl];
	P(cir,2,count)=parms->aper.d*0.5*(1-ht/hs);
	P(cir,3,count)=0xFF0000;/*rgb color */
	count++;
    }
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	real hs=parms->fit.hs->p[ifit];
	P(cir,0,count)=ht*parms->fit.thetax->p[ifit];
	P(cir,1,count)=ht*parms->fit.thetay->p[ifit];
	P(cir,2,count)=parms->aper.d*0.5*(1-ht/hs);
	P(cir,3,count)=0xFF22DD;/*rgb color */
	count++;
    }
    for(int idir=0; idir<parms->sim.ncpa_ndir; idir++){
	real hs=parms->sim.ncpa_hs->p[idir];
	P(cir,0,count)=ht*parms->sim.ncpa_thetax->p[idir];
	P(cir,1,count)=ht*parms->sim.ncpa_thetay->p[idir];
	P(cir,2,count)=parms->aper.d*0.5*(1-ht/hs);
	P(cir,3,count)=0x22FF00;/*rgb color */
	count++;
    }

    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	real hs=parms->wfs[iwfs].hs;
	int ipowfs=parms->wfs[iwfs].powfs;
	P(cir,0,count)=parms->wfs[iwfs].thetax*ht;
	P(cir,1,count)=parms->wfs[iwfs].thetay*ht;
	P(cir,2,count)=parms->aper.d*0.5*(1.-ht/hs);
	if(isfinite(hs)){//LGS
	    P(cir,3,count)=0xFF8800;
	}else if(!parms->powfs[ipowfs].lo){//Hi NGS
	    P(cir,3,count)=0xFFFF00;
	}else if(parms->powfs[ipowfs].order>1){//TTF
	    P(cir,3,count)=0x0000FF;//TTF
	}else{
	    P(cir,3,count)=0x0000FF;//TT
	}
	count++;
    }
    plot_points(fig, 1, &loc, NULL ,NULL, NULL,NULL, cir, NULL,
		"Coordinate","x (m)","y (m)", "%s",fn);
    dfree(cir);
}
/**
   ploted all the different beam directions as points. */
void plotdir(const char *fig, const PARMS_T *parms, real totfov, const char *format,...){
    format2fn;
    int ncir=1;
    dmat *cir=dnew(4, ncir);
    P(cir,0,0)=0;
    P(cir,1,0)=0;
    P(cir,2,0)=totfov/2;
    P(cir,3,0)=0x000000;/*rgb color */
    int ngroup=2+parms->npowfs;
    ngroup+=1;
    const char *legend[ngroup];
    loccell *locs=(loccell*)cellnew(ngroup, 1);
    int32_t *style=mycalloc(ngroup,int32_t);
    int count=0;
    legend[count]="Evaluation";
    style[count]=(0xFF0000<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->evl.nevl, 0, 0);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	locs->p[count]->locx[ievl]=parms->evl.thetax->p[ievl]*206265;
	locs->p[count]->locy[ievl]=parms->evl.thetay->p[ievl]*206265;
    }
    count++;
    legend[count]="DM Fitting";
    style[count]=(0xFF22DD<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->fit.nfit, 0, 0);
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	locs->p[count]->locx[ifit]=parms->fit.thetax->p[ifit]*206265;
	locs->p[count]->locy[ifit]=parms->fit.thetay->p[ifit]*206265;
    }
    count++;
    legend[count]="NCPA";
    style[count]=(0x22FF00<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->sim.ncpa_ndir, 0, 0);
    for(int ifit=0; ifit<parms->sim.ncpa_ndir; ifit++){
	locs->p[count]->locx[ifit]=parms->sim.ncpa_thetax->p[ifit]*206265;
	locs->p[count]->locy[ifit]=parms->sim.ncpa_thetay->p[ifit]*206265;
    }
    count++;
    const char* const legwfs[]={
	"LGS WFS",
	"NGS WFS",
	"PWFS",
	"TTF WFS",
	"TT WFS",
	"Other WFS",
    };
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	int ilegwfs=6;
	if(parms->powfs[ipowfs].lo){
	    if(parms->powfs[ipowfs].order==1){
		ilegwfs=4;
	    }else{
		ilegwfs=3;
	    }
	}else{
	    if(parms->powfs[ipowfs].trs){
		ilegwfs=0;
	    }else{
		if(parms->powfs[ipowfs].type==1){
		    ilegwfs=2;
		}else{
		    ilegwfs=1;
		}
	    }
	}
	legend[count]=legwfs[ilegwfs];
	locs->p[count]=locnew(parms->powfs[ipowfs].nwfs, 0, 0);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    locs->p[count]->locx[jwfs]=parms->wfs[iwfs].thetax*206265;
	    locs->p[count]->locy[jwfs]=parms->wfs[iwfs].thetay*206265;
	}
	if(isfinite(parms->powfs[ipowfs].hs)){
	    style[count]=(0xFF8800<<8)+(4<<4)+2;
	}else if(!parms->powfs[ipowfs].lo){
	    style[count]=(0xFFFF00<<8)+(4<<4)+1;
	}else if(parms->powfs[ipowfs].order>1){
	    style[count]=(0x0000FF<<8)+(4<<4)+4;
	}else{
	    style[count]=(0x0000FF<<8)+(4<<4)+1;
	}
	count++;
    }
    if(count!=ngroup){
	error("count=%d, ngroup=%d. they should equal.\n", count, ngroup);
    }
    real limit[4];
    limit[0]=limit[2]=-totfov/2;
    limit[1]=limit[3]=totfov/2;
    plot_points(fig, ngroup, locs->p, NULL, style,limit,NULL, cir, legend,
		"Asterism","x (arcsec)", "y (arcsec)", "%s",fn);
    dfree(cir);
    cellfree(locs);
    free(style);
}
/**
   Rename the log files when simulation exits.
*/
void rename_file(int sig){
    draw_final(1);
    if(disable_save) return;
    if(sig==0){
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "run_%s_%ld.log", HOST, (long)getpid());
	remove("run_done.log");
	mysymlink(fn, "run_done.log");
	snprintf(fn, PATH_MAX, "maos_%s_%ld.conf", HOST, (long)getpid());
	remove("maos_done.conf");
	mysymlink(fn, "maos_done.conf");
    }
    if(global && global->parms && global->parms->fdlock && sig!=0){
	char fn[80];
	const PARMS_T *parms=global->parms;
	for(int iseed=global->iseed; iseed<parms->sim.nseed; iseed++){
	    if(parms->fdlock->p[iseed]>0){
		close(parms->fdlock->p[iseed]);
		int seed=parms->sim.seeds->p[iseed];
		snprintf(fn, 80, "Res_%d.lock",seed);
		if(exist(fn)){
		    (void) remove(fn);
		}
	    }
	}
    }
}
/**
   Handles signals.
*/
int maos_signal_handler(int sig){
    info("maos: %s", sys_siglist[sig]);
    rename_file(sig);/*handles signal */
    if(global && global->parms && global->parms->sim.mvmport){
	mvm_client_close();
    }
    scheduler_finish(sig);
    return 0;
}
/**
   Print out usage information.
*/
static void print_usage(void){
    info(
	"Usage: maos [OPTION...] [FILE]...\n"
	"maos is a simulation tool developed to adaptive optics systems\n\n"
	"Examples:\n"
	"maos   # Run the default configuration from default.conf with all CPUs and GPUs\n\n"
	"maos -c scao_ngs.conf sim.seeds=[1 10] -d -o scao_ngs override.conf \n"
	"       # Run a single conjugate natural guide star case, with seed 1 and 10\n"
	"       # detach from the terminal and output results to folder scao_ngs\n"
	"       # and read in overriding parameters stored in override.conf and chol.conf\n"
	"\n"
	"Options: \n"
	"-h, --help        to print out this message\n"
	"-d, --detach      to detach from terminal and run in background\n"
	"-f, --force       force starting simulation without scheduler\n"
	"-n, --nthread=N   Use N threads, default is 1\n"
	"-o, --output=DIR  output the results to DIR.\n"
	"-c, --conf=FILE.conf\n"
	"                  Use FILE.conf as the baseline config instead of nfiraos.conf\n"
	"-p, --path=dir    Add dir to the internal PATH\n"
	"-g, --gpu=i       Use the i'th gpu. 0 for the first. -1 to disable. default: automatic\n"
	"-G, --ngpu=N      Use a total of N gpus.\n"
	"-r, --run=host    Run the job in another host.\n"
	"\n"
	"The following environment variables are supported\n"
	"MAOS_TOMO_SCALE=1e12 Rebalance tomography terms for single precision calculation\n"
	"MAOS_PARALLEL=1      Set to 0 to disable parallel launch\n"
	"MAOS_NO_WFS=0        Set to 1 to disable all WFS calls\n"
	"MAOS_NO_EVL=0        Set to 1 to disable evaluation calls\n"
	"MAOS_NO_RECON=0      Set to 1 to disable reconstruction calls\n"
	"MAOS_KEEP_MEM=0      Set to 1 to keep temporary memory between steps\n"
	"MAOS_MEM_DEBUG=0     Set to 1 to enable malloc/free accounting\n"
	"MAOS_MEM_VERBOSE=0   Set to 1 to print detailed malloc/free info\n"
	"MAOS_LOG_LEVEL=0     Set logging level. -3: error only, -2:  warning, -1  essential info, 0:  useful info, \n"
	"                     1:  debugging info, 2: debugging info, 3: everything.\n"
	);

    exit(0);
}

/**
   Parse command line arguments argc, argv
*/
ARG_T * parse_args(int argc, const char *argv[]){
    ARG_T *arg=mycalloc(1,ARG_T);
    char *host=NULL;
    int nthread=0;
    ARGOPT_T options[]={
	{"help",   'h',M_INT, 0, 1, (void*)print_usage, NULL},
	{"detach", 'd',M_INT, 0, 0, &arg->detach, NULL},
	{"force",  'f',M_INT, 0, 0, &arg->force, NULL},
	{"override",'O',M_INT,0, 0, &arg->over, NULL},
	{"output", 'o',M_STR, 1, 0, &arg->dirout, NULL},
	{"nthread",'n',M_INT, 1, 0, &nthread,NULL},
	{"gpu",    'g',M_INT, 2, 0, &arg->gpus, &arg->ngpu},
	{"ngpu",   'G',M_INT, 1, 0, &arg->ngpu2, NULL},
	{"conf",   'c',M_STR, 1, 0, &arg->conf, NULL},
	{"path",   'p',M_STR, 1, 1, (void*)addpath, NULL},
	{"run",    'r',M_STR, 1, 0, &host, NULL},
	{"server", 'S',M_INT, 0, 0, &arg->server,NULL},
	{NULL,     0,  0,     0, 0, NULL, NULL}
    };
    char *cmds=strnadd(argc-1, argv+1, " ");
    parse_argopt(cmds, options);
  
    if((!host || !strcmp(host, "localhost")) && !arg->detach){//forground running.
	arg->force=1;
    }else{
	if(!arg->dirout){
	    error("detach mode requires specifying -o\n");
	}
	if(getenv("MAOS_DIRECT_LAUNCH")){
	    /*being lanuched by scheduler. We are already detached, so don't daemonize again.*/
	    arg->detach=0;
	    arg->force=0;
	    detached=1;
	}else{
#ifndef MAOS_DISABLE_SCHEDULER
	    /*Detached version. Always launch through scheduler if available.*/
	    if(!host){//launch locally
		if(scheduler_launch_exe("localhost", argc, argv)){
		    warning("Launch locally without scheduler.\n");
		}
	    }else {
		const char *hostend=host+strlen(host);
		//Host maybe coma separated hostnames
		//Host2 is hostname found. Host3 is after coma.
		char *host3=NULL;
		for(char *host2=host; host2 && host2<hostend; host2=host3){
		    char* coma=strchr(host2,',');
		    if(coma){
			*coma='\0';
			host3=coma+1;
		    }else{
			host3=NULL;
		    }
		    if(scheduler_launch_exe(host2, argc, argv)){
			warning("Unable to launch maos at server %s.\n", host2);
		    }
		}
		exit(EXIT_SUCCESS);
	    }
#else
	    arg->force=1;//launch directly when scheduler is disabled.
#endif
	}
    }
    free(host);
    if(nthread<MAXTHREAD && nthread>0){
        NTHREAD=nthread;
    }
    if(arg->ngpu2>0){
	if(!arg->gpus || arg->ngpu==0){
	    arg->ngpu=arg->ngpu2;
	}else{
	    error("-g and -G cannot be specified simultaneously\n");
	}
    }else if(arg->ngpu){//check for -g-1
	for(int ig=arg->ngpu-1; ig>=0; ig--){
	    if(arg->gpus[ig]<0){
		if(ig+1==arg->ngpu){//-g-1 appear last
		    arg->gpus[0]=-1;
		    arg->ngpu=-1;
		}else{
		    //-g-1 is not last. It invalides previous -g's
		    arg->ngpu=arg->ngpu-(1+ig);
		    memcpy(arg->gpus, arg->gpus+ig+1, arg->ngpu*sizeof(int));
		}
		break;
	    }
	}
    }
    arg->confcmd=cmds;
    addpath2(".", 2);
    if(arg->dirout){
	mymkdir("%s",arg->dirout);
	if(chdir(arg->dirout)){
	    error("Unable to chdir to %s\n", arg->dirout);
	}
    }else{
	//warning("Disable saving when no -o is supplied.\n");
	disable_save=1;
    }
    return arg;
}
/**
   Creates header for saving PSFs.
 */
char *evl_header(const PARMS_T *parms, const APER_T *aper, int ievl, int iwvl, int isim){
    char header[400];
    int nembed=aper->embed->nembed->p[iwvl];
    real wvl=parms->evl.wvl->p[iwvl];
    real sumamp2=aper->sumamp2;
    snprintf(header, sizeof(header), 
	     "theta=(%.15g, %.15g) / Field location (arcsec)\n"
	     "r0=%g / Fried parameter (meter)\n"
	     "L0=%g / Outer scale (meter)\n"
	     "wvl=%g / Wavelength (m)\n"
	     "dx=%g / OPD sampling (m)\n"
	     "nfft=(%d,%d) / FFT grid size\n"
	     "dp=%g / PSF sampling (arcsec)\n"
	     "sum=%g / PSF total intensity\n"
	     "dt=%g / Total exposure (seconds)\n", 
	     ievl<0?0:parms->evl.thetax->p[ievl]*206265, ievl<0?0:parms->evl.thetay->p[ievl]*206265,
	     parms->atm.r0, parms->atm.L0->p[0],
	     wvl, parms->evl.dx, nembed, nembed, wvl/(nembed*parms->evl.dx)*206265,
	     sumamp2*nembed*nembed, parms->sim.dt*(isim-parms->evl.psfisim+1));
    return strdup(header);
}
/**
   Apply field stop to OPD using FFT.
 */
void apply_fieldstop(dmat *opd, const dmat *amp, const lmat *embed, long nembed, const dmat *fieldstop, real wvl){
    cmat *wvf=cnew(nembed, nembed);
    //cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
    real kk=2*M_PI/wvl;
    real kki=1./kk;
    real wvlh=wvl*0.5;
    comp i2pi=COMPLEX(0, kk);
    for(int iloc=0; iloc<opd->nx; iloc++){
	wvf->p[embed->p[iloc]]=amp->p[iloc]*cexp(i2pi*opd->p[iloc]);
    }
    cfft2(wvf, -1);
    ccwmd(wvf, fieldstop, 1);
    cfft2(wvf, 1);
    for(int iloc=0; iloc<opd->nx; iloc++){
	real val=carg(wvf->p[embed->p[iloc]])*kki;
	real diff=fmod(val-opd->p[iloc]+wvlh, wvl);
	if(diff<0) diff+=wvl;
	opd->p[iloc]+=diff-wvlh;
    }
    cfree(wvf);
}
/**
   Plot grid points, amplitude maps and NCPA.
 */
void plot_setup(const PARMS_T *parms, const POWFS_T *powfs,
		const APER_T *aper, const RECON_T *recon){
    if(!parms->plot.setup) return;
    plotdir("Aperture",parms,parms->sim.fov*206265,"fov");/*plot wfs/evaluation direction */
    plotloc("Aperture",parms,recon->ploc,0, "ploc");
    plotloc("Aperture",parms,recon->floc,0, "floc");
    for(int idm=0; idm<parms->ndm; idm++){
	real ht=parms->dm[idm].ht;
	plotloc("Aperture", parms, recon->aloc->p[idm], ht, "aloc%d", idm);
    }
    for(int ips=0; ips<recon->npsr; ips++){
	const real ht=recon->ht->p[ips];
	plotloc("Aperture",parms,recon->xloc->p[ips],ht, "xloc%d",ips);
    }
    drawopd("Aperture",aper->locs,aper->amp1->p,NULL,"Aperture Amplitude Map",
	    "x (m)","y (m)","aper");
    
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	drawopd("Aperture", powfs[ipowfs].loc, powfs[ipowfs].amp->p,NULL,
		"WFS Amplitude Map","x (m)","y (m)","powfs %d", ipowfs);
	if(powfs[ipowfs].amp_tel){
	    for(int wfsind=0; wfsind<parms->powfs[ipowfs].nwfs; wfsind++){
		drawopd("Aperture", powfs[ipowfs].loc, powfs[ipowfs].amp_tel->p[wfsind]->p,NULL,
			"WFS Amplitude Map","x (m)","y (m)","powfs %d tel2wfs", ipowfs);
	    }
	}
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    if(powfs[ipowfs].gradncpa){
		drawgrad("Goff",powfs[ipowfs].saloc, powfs[ipowfs].gradncpa->p[jwfs],
			 NULL, parms->plot.grad2opd,
			 "WFS Offset","x (m)", "y (m)", "Goff %d",  iwfs);
	    }
	}
    }
}


/**
   Create WFS amplitude map from coordinate, masked with annular defined by (D,Din). 
*/
dmat *mkamp(loc_t *loc, map_t *ampground, real misregx, real misregy, real D, real Din){
    dmat *amp=dnew(loc->nloc, 1);
    if(ampground){
	prop_grid(ampground, loc, amp->p, 1, misregx, misregy,1,0,0,0);
    }else{
	locannular(amp->p, loc, -misregx, -misregy, D*0.5, Din*0.5, 1);
    }
    return amp;
}
/**
   Test wfs linearity.
 */
void wfslinearity(const PARMS_T *parms, POWFS_T *powfs, const int iwfs){
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int nsa=powfs[ipowfs].saloc->nloc;
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    ccell *potf=intstat->potf->p[intstat->nsepsf>1?wfsind:0];
    cmat *potf2=0;
    ccell *otf=ccellnew(nwvl,1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	otf->p[iwvl]=cnew(potf->p[0]->nx, potf->p[0]->ny);
    }
    real pixthetax=parms->powfs[ipowfs].radpixtheta;
    real pixthetay=parms->powfs[ipowfs].pixtheta;
    dmat **mtche=NULL;
    if(parms->powfs[ipowfs].phytype_sim==1){
	if(powfs[ipowfs].intstat->mtche->ny==1){
	    mtche=powfs[ipowfs].intstat->mtche->p;
	}else{
	    mtche=powfs[ipowfs].intstat->mtche->p+nsa*wfsind;
	}
    }
    int nllt=parms->powfs[ipowfs].llt?parms->powfs[ipowfs].llt->n:0;
    real *srot=NULL;
    cmat ***petf=NULL;
    if(nllt){
	srot=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p;
	petf=mycalloc(nwvl,cmat**);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    petf[iwvl]=PPR(powfs[ipowfs].etfsim[iwvl].etf, 0, wfsind);
	}
    }

    const int nsep=41;
    const real dg=0.1;
    real gx=0, gy=0, dgx=0, dgy=0;
    dmat *ints=dnew(powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay);
    real theta=0, cx=1, sx=0;
    real isep0=-(nsep-1)*0.5;
    dmat *xg=dlinspace(isep0*dg, dg, nsep);
    writebin(xg, "wfslinearity_wfs%d_sep", iwfs);
    dfree(xg);
    dmat *gnfra=0;
    if(srot){
	gnfra=dnew(nsep,nsa*2);
    }
    dmat *gnfxy=dnew(nsep,nsa*2);
    dmat* pgnfxy=gnfxy/*PDMAT*/;
    dmat* pgnfra=gnfra/*PDMAT*/;
    const int ndir=4;
    const char *dirs[]={"x", "y", "r", "a"};
    const char *types[]={"","MF", "CoG", "MAP"};
    if(parms->powfs[ipowfs].mtchcr){
	types[1]="MFC";
    }
    int type=parms->powfs[ipowfs].phytype_sim;
    for(int dir=0; dir<ndir; dir++){
	dzero(gnfra);
	dzero(gnfxy);
	for(int isa=0; isa<nsa; isa++){
	    info_console("isa=%4d\b\b\b\b\b\b\b\b", nsa);
	    if(srot){
		theta=srot[isa];
		cx=cos(theta);
		sx=sin(theta);
	    }
	    switch(dir){
	    case 0://x
		dgx=dg*pixthetax;
		dgy=0;
		break;
	    case 1://y
		dgx=0;
		dgy=dg*pixthetay;
		break;
	    case 2://r
		dgx=cx*dg*pixthetax;
		dgy=sx*dg*pixthetax;
		break;
	    case 3://a
		dgx=-sx*dg*pixthetay;
		dgy=cx*dg*pixthetay;
		break;
	    default:
		error("Invalid dir\n");
	    }
	    for(int isep=0; isep<nsep; isep++){
		gx=dgx*(isep+isep0);
		gy=dgy*(isep+isep0);
		dzero(ints);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    real wvlsig=parms->wfs[iwfs].wvlwts->p[iwvl]
			*parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
		    int idtf=powfs[ipowfs].dtf[iwvl].si->ny>1?wfsind:0;
		    int idtfsa=powfs[ipowfs].dtf[iwvl].si->nx>1?isa:0;
		    dspcell*  psi=powfs[ipowfs].dtf[iwvl].si/*PDSPCELL*/;
		    dsp *sis=P(psi,idtfsa,idtf);
		    real wvl=parms->powfs[ipowfs].wvl->p[iwvl];
		    real dtheta1=powfs[ipowfs].pts->nx*powfs[ipowfs].pts->dx*parms->powfs[ipowfs].embfac/wvl;
		    if(petf){
			ccp(&potf2, potf->p[isa+nsa*iwvl]);
			ccwm(potf2,petf[iwvl][isa]);
		    }else{
			potf2=potf->p[isa+nsa*iwvl];
		    }
		    ctilt2(otf->p[iwvl], potf2, gx*dtheta1, gy*dtheta1, 0);
		    cfft2(otf->p[iwvl], 1);
		    dspmulcreal(ints->p, sis, otf->p[iwvl]->p, wvlsig);
		}
		//ddraw("ints", ints, NULL, NULL, "ints", "x", "y", "ints"); PAUSE;
		real g[3]={0,0,0};
		//notice that all the following gives gradients in x/y coord only.
		switch(type){
		case 0://no-op
		    break;
		case 1:{/*(constraint) Matched filter give gradients along x/y*/
		    dmulvec(g, mtche[isa], ints->p,1.);
		}
		    break;
		case 2:{/*tCoG gives gradients along r/a*/
		    dcog(g,ints,0.,0.,
			 powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2],
			 powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2+1], 0);
		    g[0]*=pixthetax;
		    g[1]*=pixthetay;
		}
		    break;
		case 3:{/*MAP gives gradient along r/a(?)*/
		    g[0]=gx+pixthetax*0.1;
		    g[1]=gy+pixthetay*0.1;
		    g[2]=1;
		    maxapriori(g, ints, parms, powfs, iwfs, isa, 1, 0, 1);
		}
		    break;
		default:
		    error("Invalid");
		}
		if(type==1 || !srot){
		    P(pgnfxy,isep,isa)=g[0]/pixthetax;
		    P(pgnfxy,isep,isa+nsa)=g[1]/pixthetay;
		    
		    if(srot){/*obtain gradients in r/a coord*/
			P(pgnfra,isep,isa)=(g[0]*cx+g[1]*sx)/pixthetax;
			P(pgnfra,isep,isa+nsa)=(-g[0]*sx+g[1]*cx)/pixthetay;
		    }
		}else{
		    P(pgnfra,isep,isa)=g[0]/pixthetax;
		    P(pgnfra,isep,isa+nsa)=g[1]/pixthetay;
		    P(pgnfxy,isep,isa)=(g[0]*cx-g[1]*sx)/pixthetax;
		    P(pgnfxy,isep,isa+nsa)=(g[0]*sx+g[1]*cx)/pixthetay;
		}
	    }
	}/*for isa*/
	writebin(gnfxy, "wfslinearity_xy_wfs%d_%s_%s", iwfs, types[type],dirs[dir]);
	if(srot){
	    writebin(gnfra, "wfslinearity_ra_wfs%d_%s_%s", iwfs, types[type],dirs[dir]);
	}
    }
    dfree(gnfxy);
    dfree(gnfra);
    dfree(ints);
    ccellfree(otf);
    if(petf){
	free(petf);
	cfree(potf2);
    }
}
/**
   Compute spherical aberration
*/
void lgs_wfs_sph_psd(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon, const int iwfs){
    int ipowfs=parms->wfs[iwfs].powfs;
    /*First save matched filter and i0*/
    dcell *mtche=powfs[ipowfs].intstat->mtche;
    dcell *i0=powfs[ipowfs].intstat->i0;
    dmat *i0sum=powfs[ipowfs].intstat->i0sum;
    powfs[ipowfs].intstat->mtche=0;
    powfs[ipowfs].intstat->i0=0;
    powfs[ipowfs].intstat->i0sum=0;
    //writebin(i0, "i0_initial");
    //writebin(mtche, "mtche_initial");
    dmat *opd=zernike(recon->ploc, parms->aper.d, 2, 15, 1);
    //writebin(opd, "ropd");
    dmat *GR=0;
    dspmm(&GR, recon->GP->p[iwfs], opd, "nn", 1);
    dfree(opd);
    dmat *RR=dpinv(GR, recon->saneai->p[iwfs+iwfs*parms->nwfsr]);
    //writebin(GR, "GR");
    //writebin(RR, "RR");
    int nsa=powfs[ipowfs].saloc->nloc;
    dmat *gradmf=dnew(nsa, 2);
    dmat *gradcg=dnew(nsa, 2);
    int ncol=1000;
    int dtrat=parms->powfs[ipowfs].llt->coldtrat;
    dmat *rmodmf=dnew(RR->nx, ncol/dtrat);
    dmat *rmodcg=dnew(RR->nx, ncol/dtrat);
    real scale=1;
    real pixthetax=parms->powfs[ipowfs].radpixtheta;
    real pixthetay=parms->powfs[ipowfs].pixtheta;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    real *srot=(parms->powfs[ipowfs].radpix)?
	powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p:NULL;
    for(int icol=0; icol<1000; icol+=dtrat){
	setup_powfs_etf(powfs, parms, 0, ipowfs, 0, icol);
	gensei(parms, powfs, ipowfs);
	dcell *i0_new=powfs[ipowfs].intstat->i0;
	//writebin(i0_new, "i0_%d", icol);
	for(int isa=0; isa<nsa; isa++){
	    real geach[3]={0,0,1};
	    dmulvec(geach, mtche->p[isa], i0_new->p[isa]->p, 1);
	    if(parms->powfs[ipowfs].sigmatch){
		scale=i0sum->p[isa]/dsum(i0_new->p[isa]);
	    }
	    gradmf->p[isa]=geach[0]*scale;
	    gradmf->p[isa+nsa]=geach[1]*scale;
	    {
		dcog(geach, i0_new->p[isa], 0, 0, 
		     powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2],
		     powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2+1], 0);
		geach[0]*=pixthetax;
		geach[1]*=pixthetay;
		if(srot){
		    real theta=srot[isa];
		    real cx=cos(theta);
		    real sx=sin(theta);
		    real tmp=geach[0]*cx-geach[1]*sx;
		    geach[1]=geach[0]*sx+geach[1]*cx;
		    geach[0]=tmp;
		}
		gradcg->p[isa]=geach[0];
		gradcg->p[isa+nsa]=geach[1];
	    }
	    
	}
	dmulvec(rmodmf->p+icol/dtrat*rmodmf->nx, RR, gradmf->p, 1);
	dmulvec(rmodcg->p+icol/dtrat*rmodcg->nx, RR, gradcg->p, 1);
    }
    dfree(GR);
    dfree(RR);
    writebin(rmodmf, "rmod_mf");
    writebin(rmodcg, "rmod_cg");
    dfree(rmodmf);
    dfree(rmodcg);
    dfree(gradmf);
    dfree(gradcg);
    dcellfree(mtche);
    dcellfree(i0);
    dfree(i0sum);
}
typedef struct {
    const PARMS_T *parms;
    const POWFS_T *powfs;
    const dmat *ints;
    ccell *fotf;
    ccell *otf;//temporary.
    real bkgrnd;
    real rne;
    int noisy;
    int iwfs;
    int isa;
}mapdata_t;
/**
   The function to evaluate the result at x.
*/
static real mapfun(real *x, mapdata_t *info){
    const dmat *ints=info->ints;
    ccell *fotf=info->fotf;
    ccell *otf=info->otf;
    const PARMS_T *parms=info->parms;
    const POWFS_T *powfs=info->powfs;
    int iwfs=info->iwfs;
    int ipowfs=parms->wfs[iwfs].powfs;
    int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    int isa=info->isa;
    int nsa=fotf->nx;
    int nwvl=fotf->ny;
    if(!otf){
	info->otf=ccellnew(nwvl,1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    info->otf->p[iwvl]=cnew(info->fotf->p[0]->nx, info->fotf->p[0]->ny);
	    //cfft2plan(info->otf->p[iwvl], 1);
	    //cfft2plan(info->otf->p[iwvl], -1);
	}
	otf=info->otf;
    }
    dmat *ints2=dnew(ints->nx, ints->ny);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	real wvlsig=parms->wfs[iwfs].wvlwts->p[iwvl]
	    *parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
	dspcell*  psi=powfs[ipowfs].dtf[iwvl].si/*PDSPCELL*/;
	int idtf=powfs[ipowfs].dtf[iwvl].si->ny>1?wfsind:0;
	int idtfsa=powfs[ipowfs].dtf[iwvl].si->nx>1?isa:0;
	dsp *sis=P(psi,idtfsa,idtf);
	real wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	real dtheta1=powfs[ipowfs].pts->nx*powfs[ipowfs].pts->dx*parms->powfs[ipowfs].embfac/wvl;
	ctilt2(info->otf->p[iwvl], info->fotf->p[isa+nsa*iwvl], x[0]*dtheta1, x[1]*dtheta1, 0);
	cfft2(info->otf->p[iwvl], 1);
	dspmulcreal(ints2->p, sis, info->otf->p[iwvl]->p, wvlsig*x[2]);
    }
 
    real sigma=0;
    if(info->noisy){
	real noise=info->rne*info->rne+info->bkgrnd;
	for(int i=0; i<ints->nx*ints->ny; i++){
	    sigma+=pow(ints->p[i]-ints2->p[i],2)/(ints2->p[i]+noise);
	}
    }else{
	for(int i=0; i<ints->nx*ints->ny; i++){
	    sigma+=pow(ints->p[i]-ints2->p[i],2);
	}
    }
    /*dbg("Map fun called with [%g %g] %g, sigma=%g. noisy=%d\n", x[0], x[1], x[2], sigma, info->noisy);*/
    dfree(ints2);
    return sigma;
}
/**
   Implements MAP tracking algorithm. The polar coordinate is implicitly taken care of in mapfun
*/
void maxapriori(real *g, const dmat *ints, const PARMS_T *parms, 
		const POWFS_T *powfs, int iwfs, int isa, int noisy,
		real bkgrnd, real rne){
    int ipowfs=parms->wfs[iwfs].powfs;
    int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    real pixthetax=parms->powfs[ipowfs].radpixtheta;
    real pixthetay=parms->powfs[ipowfs].pixtheta;
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    ccell *fotf=intstat->fotf->p[intstat->nsepsf>1?wfsind:0];
    mapdata_t data={parms, powfs, ints, fotf, NULL, bkgrnd, rne, noisy, iwfs, isa};
    //info("isa %d: %.4e %.4e %.2f", isa, g[0], g[1], g[2]);
    int ncall=dminsearch(g, 3, MIN(pixthetax, pixthetay)*1e-2, 5000, (dminsearch_fun)mapfun, &data);
    ccellfree(data.otf);
    /* convert to native format along x/y or r/a to check for overflow*/
    if(parms->powfs[ipowfs].radpix){
	real theta=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p[isa];
	real cx=cos(theta);
	real sx=sin(theta);
	real tmp=g[0]*cx+g[1]*sx;
	g[1]=-g[0]*sx+g[1]*cx;
	g[0]=tmp;
    }
    real gx=g[0]/pixthetax*2./ints->nx;
    real gy=g[1]/pixthetay*2./ints->ny;
    if(fabs(gx)>0.55||fabs(gy)>0.55){
	warning("sa %4d iter %3d: wrapped: gx=%6.3f, gy=%6.3f ==> ", isa, ncall, gx, gy);
	gx=gx-floor(gx+0.5);
	gy=gy-floor(gy+0.5);
	warning("gx=%6.3f, gy=%6.3f\n", gx, gy);
	g[0]=pixthetax*ints->nx/2*gx;
	g[1]=pixthetay*ints->ny/2*gy;
    }
    //info("==> %.4e %.4e %.2f after %d iter\n", g[0], g[1], g[2], ncall);
    if(parms->powfs[ipowfs].radpix){
	real theta=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p[isa];
	real cx=cos(theta);
	real sx=sin(theta);
	real tmp=g[0]*cx-g[1]*sx;
	g[1]=g[0]*sx+g[1]*cx;
	g[0]=tmp;
    }
}

/**
   Compute the focus adjustment need to apply to OPD of wfs. Used in both CPU and GPU code.
*/
real wfsfocusadj(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    const POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const int isim=simu->wfsisim;
    real focus=0;
    if(parms->powfs[ipowfs].llt){
	if(powfs[ipowfs].focus){//input focus error due to range variation
	    focus+=PR(powfs[ipowfs].focus, isim, wfsind);
	}
	/*
	  The following has been migrated to sim_update_etf to change the hs directly.
	  if(simu->zoomreal && parms->powfs[ipowfs].llt){
	  if(simu->zoompos && simu->zoompos->p[iwfs]){//save
	  simu->zoompos->p[iwfs]->p[isim]=simu->zoomreal->p[iwfs];
	  }
	  focus-=simu->zoomreal->p[iwfs];
	  }*/
    }
    if(simu->telfocusreal){
	focus-=simu->telfocusreal->p[0]->p[0];
    }
    return focus;
}
/**
   Expected averaged position of dithering signal during WFS integration. Called when (isim+1)%dtrat=0
   
   
*/

void dither_position(real *cs, real *ss, int alfsm, int dtrat, int npoint, int isim, real deltam){
    //adjust for delay due to propagation, and computation delay if delay is not 2 frame.
    const int adjust=alfsm+1-dtrat;
    //adjust to get delay at beginning of integration
    const int adjust2=dtrat-1;
    const real anglei=(2*M_PI/npoint);
    const real angle=((isim-adjust-adjust2)/dtrat)*anglei+deltam;
    const real angle2=((isim-adjust)/dtrat)*anglei+deltam;
    const real delay=(real)adjust/dtrat;
    const real beta=1+delay+floor(-delay);
    const real scale=1./(beta*beta+(1-beta)*(1-beta));
    //use average of two places during accumulation and scale
    *cs=(beta*cos(angle)+(1-beta)*cos(angle2))*scale;
    *ss=(beta*sin(angle)+(1-beta)*sin(angle2))*scale;
}
/**
   Find peak, then using parabolic fit on 3x3 window around it.
*/
/*
  void parabolic_peak_fit(real *grad, dmat *corr){
  real valmax=0;
  int jy=0, jx=0;
  //Find Peak location (jx, jy)
  for(int iy=1; iy<corr->ny-1; iy++){
  for(int ix=1; ix<corr->nx-1; ix++){
  if(P(corr,ix,iy)>valmax){
  jy=iy; jx=ix;
  valmax=P(corr,ix,iy);
  }
  }
  }
  //Calculate 1d sum of 3 row/columns.
  real vx[3], vy[3];
  for(long iy=0; iy<3; iy++){
  vy[iy]=0; vx[iy]=0;
  for(long ix=0; ix<3; ix++){
  vy[iy]+=P(corr, ix+jx-1, iy+jy-1);
  vx[iy]+=P(corr, iy+jx-1, ix+jy-1);
  }
  }
  //Parabolic fit.
  real px[2], py[2];
  px[0]=(vx[0]+vx[2])*0.5-vx[1];
  py[0]=(vy[0]+vy[2])*0.5-vy[1];
  px[1]=(vx[2]-vx[0])*0.5;
  py[1]=(vy[2]-vy[0])*0.5;
  //Center
  grad[0]=px[0]==0?0:(-px[1]/(2*px[0])+jx-(corr->nx-1)*0.5);
  grad[1]=py[0]==0?0:(-py[1]/(2*py[0])+jy-(corr->ny-1)*0.5);
  }
*/
/**
   Fit 3 points around the peak of a 1-d profile.
*/
real parabolic_peak_1d(dmat *corr){
    real valmax=0;
    int jx=0;
    for(int ix=1; ix<corr->nx-1; ix++){
	if(P(corr,ix)>valmax){
	    jx=ix;
	    valmax=P(corr,ix);
	}
    }
    //Parabolic fit.
    real px[2];
    px[0]=(P(corr, jx+1)+P(corr, jx-1))*0.5-P(corr, jx);
    px[1]=(P(corr, jx+1)-P(corr, jx-1))*0.5;
    return px[0]==0?0:(-px[1]/(2*px[0])+jx-(corr->nx-1)*0.5);
}
/**
   First sum along 1 dimension, then fit 3 points around the peak. More robust than the old method.
*/
void parabolic_peak_sum(real *grad, dmat *corr, int nbox){
    const long nx=corr->nx;
    const long ny=corr->ny;
    if(nbox<=0 || nbox>nx) nbox=nx;
    if(nbox<=0 || nbox>ny) nbox=ny;
    dmat *corrx=dnew(nbox,1);
    dmat *corry=dnew(nbox,1);
    const long offx=(nx-nbox)/2;
    const long offy=(ny-nbox)/2;
    for(int iy=0; iy<nbox; iy++){
	for(int ix=0; ix<nbox; ix++){
	    P(corrx, ix)+=P(corr, ix+offx, iy+offy);
	    P(corry, iy)+=P(corr, ix+offx, iy+offy);
	}
    }
    grad[0]=parabolic_peak_1d(corrx);
    grad[1]=parabolic_peak_1d(corry);
    if (0){
	info("grad=%g, %g\n", grad[0], grad[1]);
	writebin(corr, "corr");
	exit(0);
    }
}
/**
   Calculate gradients using current specified algorithm
*/
void shwfs_grad(dmat **pgrad, dmat *ints[], const PARMS_T *parms, const POWFS_T *powfs, const int iwfs, const int phytype){
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const real rne=parms->powfs[ipowfs].rne;
    const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    dmat **mtche=NULL;
    real *i0sum=NULL;
    real i0sumg=0;
    if(phytype==1){
	mtche=PPR(powfs[ipowfs].intstat->mtche, 0, wfsind);
    }
    if(powfs[ipowfs].intstat && powfs[ipowfs].intstat->i0sum){
	i0sum=PPR(powfs[ipowfs].intstat->i0sum, 0, wfsind);
	i0sumg=PR(powfs[ipowfs].intstat->i0sumsum, wfsind, 0);
    }
 
    const real *srot=(parms->powfs[ipowfs].radpix)?PR(powfs[ipowfs].srot, wfsind, 0)->p:NULL;
    real pixthetax=parms->powfs[ipowfs].radpixtheta;
    real pixthetay=parms->powfs[ipowfs].pixtheta;
    /*output directly to simu->gradcl. replace */
    if(!*pgrad){
	*pgrad=dnew(nsa*2, 1);
    }
    real *pgradx=(*pgrad)->p;
    real *pgrady=pgradx+nsa;
    real i1sum=0;
    dmat *corr=0;
    /**
       note on powfs.sigmatch:
       0: Does not normalize measurement by intensity (linear)
       1: normalize measurement per subaperture by intensity 
       2: normalize measurement globally for all subapertures
     */
    if(parms->powfs[ipowfs].sigmatch==2){
	for(int isa=0; isa<nsa; isa++){
	    i1sum+=dsum(ints[isa]);
	}
    }
    //sigtot: expected total signal level for full subapertures
    real sigtot=NAN;
    switch(parms->powfs[ipowfs].sigmatch){
    case 0:
	sigtot=parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
	break;
    case 2:
	sigtot=(i1sum/powfs[ipowfs].saasum);
	break;
    }

    for(int isa=0; isa<nsa; isa++){
	real geach[3]={0,0,1};
	switch(phytype){
	case 1:{//matched filter
	    real scale=1.;
	    switch(parms->powfs[ipowfs].sigmatch){
	    case 0://no normalization
		break;
	    case 1://match per subaperture
		scale=i0sum[isa]/dsum(ints[isa]);
		break;
	    case 2://match globally.
		scale=i0sumg/i1sum;
		break;
	    }
	    dmulvec(geach, mtche[isa],ints[isa]->p,scale);
	}
	    break;
	case 2:{//CoG
	    real sumi=0;
	    switch(parms->powfs[ipowfs].sigmatch){
	    case 0://normalization use model intensity (linear model)
		if(i0sum){
		    sumi=i0sum[isa];
		}else{
		    sumi=sigtot*powfs[ipowfs].saa->p[isa];
		}
		break;
	    case 1://normalization use current intensity (usual method, non-linear)
		break;
	    case 2://normalized use scaled current intensity (non-linear)
		sumi=sigtot*powfs[ipowfs].saa->p[isa];
		break;
	    }
	    dcog(geach,ints[isa],0.,0.,
		 P(powfs[ipowfs].cogcoeff->p[wfsind],0,isa),
		 P(powfs[ipowfs].cogcoeff->p[wfsind],1,isa), sumi);
	    geach[0]*=pixthetax;
	    geach[1]*=pixthetay;
	  
	}
	    break;
	case 3:{//MAP: (to be removed. This algorithm is not very useful.)
	    geach[0]=pgradx[isa];//warm restart
	    geach[1]=pgrady[isa];
	    maxapriori(geach, ints[isa], parms, powfs, iwfs, isa, 1, bkgrnd, rne);
	}
	    break;
	case 4:{//Correlation. Peak first.
	    dcorr(&corr, ints[isa], powfs[ipowfs].intstat->i0->p[isa]);
	    dpara3(geach, corr);
	    if((fabs(geach[0])+fabs(geach[1]))>powfs[ipowfs].pixpsax){
		warning_once("wfs %d, subaperture %d has incorrect gradient measurement\n", iwfs, isa);
		geach[0]=0;
		geach[1]=0;
	    }
	    geach[0]*=pixthetax;
	    geach[1]*=pixthetay;
	}
	    break;
	case 5:{//Correlation. Sum first (to be removed. Not working well)
	    dcorr(&corr, ints[isa], powfs[ipowfs].intstat->i0->p[isa]);
	    parabolic_peak_sum(geach, corr, 5);
	    geach[0]*=pixthetax;
	    geach[1]*=pixthetay;
	}
	    break;
	default:
	    error("Invalid");
	}
	if(phytype>1 && srot){
	    real theta=srot[isa];
	    real cx=cos(theta);
	    real sx=sin(theta);
	    real tmp=geach[0]*cx-geach[1]*sx;
	    geach[1]=geach[0]*sx+geach[1]*cx;
	    geach[0]=tmp;
	}
	pgradx[isa]=geach[0];
	pgrady[isa]=geach[1];
    }//for isa
    dfree(corr);
}
/**
   Read cell array from file specified by whole name or prefix.
*/
dcell *dcellread_prefix(const char *file, const PARMS_T *parms, int ipowfs){
    dcell *nea=0;
    int iwfs0=parms->powfs[ipowfs].wfs->p[0];
    if(!file){
	return nea;
    }else if(zfexist("%s_powfs%d.bin", file, ipowfs)){
	//info("using %s_powfs%d.bin\n", file, ipowfs);
	nea=dcellread("%s_powfs%d.bin", file, ipowfs);
    }else if(zfexist("%s_wfs%d.bin", file, iwfs0)){
	nea=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<nea->nx; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    //info("using %s_wfs%d.bin\n", file, iwfs);
	    nea->p[jwfs]=dread("%s_wfs%d.bin", file, iwfs);
	}
    }else if(zfexist(file)){
	//info("using %s\n", file);
	nea=dcellread("%s", file);
    }else{
	error("%s_powfs%d.bin or %s_wfs%d.bin not found\n", file, ipowfs, file, iwfs0);
    }
    return nea;
}
/**
   Wait for dmreal to be available in event driven simulation.
*/
void wait_dmreal(SIM_T *simu, int isim){
    if(PARALLEL==2){
	//if(simu->dmreal_isim!=-1){
	while(simu->dmreal_isim!=isim){
	    //if(simu->dmreal_isim+1!=isim){
	    //	error("dmreal_isim=%d, isim=%d\n", simu->dmreal_isim, isim);
	    //}
	    //dbg("waiting: dmreal_isim is %d need %d...\n", simu->dmreal_isim, isim);
	    pthread_cond_wait(&simu->dmreal_condr, &simu->dmreal_mutex);
	}
	//dbg("ready: dmreal_isim is %d need %d\n", simu->dmreal_isim, isim);
	atomicadd(&simu->dmreal_count,1);
	//}
	pthread_mutex_unlock(&simu->dmreal_mutex);
	pthread_cond_signal(&simu->dmreal_condw);
    }
}
/**
   Concatenate and plot subaperture images.
 */
void draw_ints(const dcell *ints, const loc_t *saloc, int iwfs){
    dmat *ints2=0;
    if(ints->nx==1){//T
	ints2=dref(ints->p[0]);
    }else if(ints->nx==4){//TTF
	dcell *ints3=dcellref(ints);
	cellreshape(ints3, 2, 2);
	ints2=dcell2m(ints3);
	dcellfree(ints3);
    }else{
	dcell *ints3=0;
	loc_embed_cell(&ints3, saloc, ints);
	ints2=dcell2m(ints3);
	dcellfree(ints3);
    }
    ddraw("Ints", ints2, NULL, NULL, "WFS Subaperture Images",
	  "x", "y", "wfs %d", iwfs);
    dfree(ints2);
}
