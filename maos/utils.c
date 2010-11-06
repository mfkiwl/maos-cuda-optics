/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <math.h>
#include <search.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <signal.h>
#include <fcntl.h>           /* For O_* constants */
#include <errno.h>
#include <getopt.h>
#include "maos.h"
char *dirsetup=NULL;
char *dirskysim=NULL;
/**
   \file maos/utils.c
   A few utility routines
 */
/**
   add photon and read out noise.  pcalib part of bkgrnd is calibrated
out. pcalib2 part of bkgrnd2 is calibrated out.  */
void addnoise(dmat *A,              /**<The pixel intensity array*/
	      rand_t* rstat,   /**<The random stream*/
	      const double bkgrnd,  /**<background in PDEs per pixel per frame*/
	      const double pcalib,  /**<Fraction of bkgrnd that is calibrated out*/
	      const double *bkgrnd2,/**<background in PDEs of each pixel per frame.*/
	      const double pcalib2, /**<Fraction of bkgrnd taht is calibrated out*/
	      const double rne      /**<Read out noise per pixel per read*/
	      ){
 
    if(bkgrnd2){
	for(int ix=0; ix<A->nx*A->ny; ix++){
	    A->p[ix]=randp(rstat,A->p[ix]+bkgrnd+bkgrnd2[ix])
		-bkgrnd*pcalib-bkgrnd2[ix]*pcalib2
		+rne*randn(rstat);
	}
    }else{
	for(int ix=0; ix<A->nx*A->ny; ix++){
	    A->p[ix]=randp(rstat,A->p[ix]+bkgrnd)
			   -bkgrnd*pcalib+rne*randn(rstat);
	}
    }
}
/**
   Calls create_metapupil is simplified interface by returning a map_t object.
 */
map_t * create_metapupil_wrap(const PARMS_T *parms, double ht,double dx,
				double offset,double guard, long nin,
				T_TYPE type, int pad,int square){
    map_t * amp=calloc(sizeof(map_t),1);
    create_metapupil(parms,ht,dx,offset,&(amp->nx),&(amp->ny),
		     &(amp->ox),&(amp->oy),&(amp->p),guard, nin, type,pad,square);
    amp->dx=dx;
    amp->h=ht;
    return amp;
}
/**
   create a metapupil map, with size nx*ny, origin at (ox,oy), sampling of dx,
   height of ht, that can cover all the WFS and science beams.  

   offset: distance in pixel from the point closest to the origin to origin (right
   side).  
   0: there is a point on the origin.  
   1/2: the closest point to origin is 1/2 pixel.  
   
   pad: 1: round nx, ny to power of 2.  */

void create_metapupil(const PARMS_T *parms, double ht,double dx,
		      double offset, long* nxout, long* nyout,
		      double *oxout, double *oyout, double**map, 
		      double guard, long nin, T_TYPE type,int pad,int square){
    double R=parms->aper.d/2;
    double maxx=0,maxy=0;
    double sx,sy;//temporary variables
    int i;
    int use_wfs_hi=1;
    int use_wfs_lo=1;
    int use_evl=1;
    int use_fit=1;
    const char *typestr[]={
	"PLOC","ALOC","XLOC","ATM"
    };
    if(parms->tomo.alg==0 && parms->tomo.split==1){
	switch(type){
	case T_PLOC:
	    break;
	case T_ALOC:
	    break;
	case T_XLOC:
	    /*use_evl=0;
	    use_fit=1;
	    if(parms->tomo.split==1){
		use_wfs_lo=0;
		}*/
	    break;
	case T_ATM:
	    break;
	default:
	    warning("Unknown type\n");
	}
    }
    //find minimum map size to cover all the beams
    for(i=0; i<parms->nwfs; i++){
	int ipowfs=parms->wfs[i].powfs;
	if((parms->powfs[ipowfs].lo && !use_wfs_lo) 
	   ||(!parms->powfs[ipowfs].lo && !use_wfs_hi)){
	    warning("Skip wfs %d when making grid for type %s\n",i,typestr[type]);
	    continue;
	}
	sx=fabs(parms->wfs[i].thetax*ht)
	    +(1.-ht/parms->powfs[ipowfs].hs)*R+guard;
	sy=fabs(parms->wfs[i].thetay*ht)
	    +(1.-ht/parms->powfs[ipowfs].hs)*R+guard;
	if(sx>maxx) maxx=sx;
	if(sy>maxy) maxy=sy;
    }
    if(use_evl){
	for(i=0; i<parms->evl.nevl; i++){
	    sx=fabs(parms->evl.thetax[i]*ht)+R+guard;
	    sy=fabs(parms->evl.thetay[i]*ht)+R+guard;
	    if(sx>maxx) maxx=sx;
	    if(sy>maxy) maxy=sy;
	}
    }else{
	warning("Skip evl when making grid for type %s\n", typestr[type]);
    }
    if(use_fit){
	for(i=0; i<parms->fit.nfit; i++){
	    sx=fabs(parms->fit.thetax[i]*ht)+R+guard;
	    sy=fabs(parms->fit.thetay[i]*ht)+R+guard;
	    if(sx>maxx) maxx=sx;
	    if(sy>maxy) maxy=sy;
	}
    }else{
	warning("Skip fit when making grid for type %s\n", typestr[type]);
    }
    //normalized grid size +3 is extra guard band
    maxx=ceil(maxx/dx)+2;
    maxy=ceil(maxy/dx)+2;
    long nx,ny;
    nx=iceil(maxx+offset+1)*2;
    ny=iceil(maxy+offset+1)*2;
    if(pad){//pad to power of 2
	nx=1<<iceil(log2((double)nx));
	ny=1<<iceil(log2((double)ny));
    }
    //Make it square
    nx=(nx<ny)?ny:nx;
    ny=nx;
    if(nin>1){
	if(nin<nx){
	    warning("nin=%ld is too small\n",nin);
	}
	//warning("overriding nx=%ld, ny=%ld by supplied nin=%ld\n",nx,ny,nin);
	nx=nin;
	ny=nin;
    }
    double ox,oy;
    ox=((nx/2)-offset);
    oy=((ny/2)-offset);

    if(nxout)
	*nxout=nx;
    if(nyout)
	*nyout=ny;
    if(oxout)
	*oxout=-ox*dx;
    if(oyout)
	*oyout=-oy*dx;

    if(map && square){
	dmat *dmap=dnew(nx,ny);
	*map=dmap->p;
	dset(dmap,1);
	dfree_keepdata(dmap);
    }else if(map){
	//guard+=0.707;//necessary with dcircle_deprecated. not with current dcircle
	dmat *dmap=dnew(nx,ny);
	*map=dmap->p;
	double Rn=R/dx;
	double Rg=guard/dx;
	for(i=0; i<parms->nwfs; i++){
	    int ipowfs=parms->wfs[i].powfs;
	    if((parms->powfs[ipowfs].lo && !use_wfs_lo) 
	       ||(!parms->powfs[ipowfs].lo && !use_wfs_hi)){
		//warning("Skip wfs %d when making grid for type %s\n",i,typestr[type]);
		continue;
	    }
	    sx=ox+(parms->wfs[i].thetax*ht)/dx;
	    sy=oy+(parms->wfs[i].thetay*ht)/dx;
	    double RR=Rn*(1.-ht/parms->powfs[ipowfs].hs)+Rg;
	    dcircle_symbolic(dmap,sx,sy,RR);
	}

	if(use_evl){
	    for(i=0; i<parms->evl.nevl; i++){
		sx=ox+(parms->evl.thetax[i]*ht)/dx;
		sy=oy+(parms->evl.thetay[i]*ht)/dx;
		double RR=Rn+Rg;
		dcircle_symbolic(dmap,sx,sy,RR);
	    }
	}else{
	    warning("Skip evl when making grid for type %s\n", typestr[type]);
	}
	if(use_fit){
	    for(i=0; i<parms->fit.nfit; i++){
		sx=ox+(parms->fit.thetax[i]*ht)/dx;
		sy=oy+(parms->fit.thetay[i]*ht)/dx;
		double RR=Rn+Rg;
		dcircle_symbolic(dmap,sx,sy,RR);
	    }
	}else{
	    warning("Skip fit when making grid for type %s\n", typestr[type]);
	}

	for(i=0; i<nx*ny; i++){
	    dmap->p[i]=(dmap->p[i])>1.e-15?1:0;
	}
	dfree_keepdata(dmap);
	//free(dmap->nref);
	//free(dmap);//don't free data
    }
}
/**
   Plot the projection of all the different FoVs on a certain grid. used in setup_recon
 */
void plotloc(char *fig, const PARMS_T *parms, 
	     loc_t *loc, double ht, char *format,...){
    format2fn;
    int ncir=parms->evl.nevl + parms->nwfs;
    double (*cir)[4];
    cir=(double(*)[4])calloc(ncir*4,sizeof(double));
    int count=0;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	cir[count][0]=ht*parms->evl.thetax[ievl];
	cir[count][1]=ht*parms->evl.thetay[ievl];
	cir[count][2]=parms->aper.d*0.5;
	cir[count][3]=900;//rgb color
	count++;
    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	double hs=parms->powfs[parms->wfs[iwfs].powfs].hs;
	cir[count][0]=parms->wfs[iwfs].thetax*ht;
	cir[count][1]=parms->wfs[iwfs].thetay*ht;
	cir[count][2]=parms->aper.d*0.5*(1.-ht/hs);
	if(isinf(hs)){
	    cir[count][3]=290;//rgb color
	}else{
	    cir[count][3]=992;
	}
	count++;
    }
    plot_coord(fig, loc->nloc, loc->locx, loc->locy,NULL, NULL,ncir, cir, 
	       "Coordinate","x (m)","y (m)", "%s",fn);
    free(cir);
}
/**
   ploted all the different beam directions as points. used in setup_parms */
void plotdir(char *fig, const PARMS_T *parms, double totfov, char *format,...){
    format2fn;
    int ncir=1;
    double (*cir)[4];
    cir=(double(*)[4])calloc(ncir*4,sizeof(double));
    long npts=parms->evl.nevl + parms->nwfs + parms->fit.nfit;
    cir[0][0]=0;
    cir[0][1]=0;
    cir[0][2]=totfov/2;
    cir[0][3]=000;//rgb color
    double *ptsx=calloc(npts,sizeof(double));
    double *ptsy=calloc(npts,sizeof(double));
    long *style=calloc(npts, sizeof(long));
    int nevl=parms->evl.nevl;
    int ind=0;
    for(int ievl=0; ievl<nevl; ievl++){
	ptsx[ind]=parms->evl.thetax[ievl]*206265;
	ptsy[ind]=parms->evl.thetay[ievl]*206265;
	style[ind]=(900<<8)+(4<<4)+3;
	ind++;
    }
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	ptsx[ind]=parms->fit.thetax[ifit]*206265;
	ptsy[ind]=parms->fit.thetay[ifit]*206265;
	style[ind]=(918<<8)+(4<<4)+3;
	ind++;
    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	ptsx[ind]=parms->wfs[iwfs].thetax*206265;
	ptsy[ind]=parms->wfs[iwfs].thetay*206265;
	int ipowfs=parms->wfs[iwfs].powfs;
	if(!isinf(parms->powfs[ipowfs].hs)){
	    style[ind]=(940<<8)+(4<<4)+2;
	}else if(!parms->powfs[ipowfs].lo){
	    style[ind]=(990<<8)+(4<<4)+1;
	}else if(parms->powfs[ipowfs].order>1){
	    style[ind]=(9<<8)+(4<<4)+4;
	}else{
	    style[ind]=(9<<8)+(4<<4)+1;
	}
	ind++;
    }
    double limit[4];
    limit[0]=limit[2]=-totfov/2;
    limit[1]=limit[3]=totfov/2;
    plot_coord(fig, npts, ptsx, ptsy, style,limit,ncir,cir,
	       "Asterism","x (arcsec)", "y (arcsec)", "%s",fn);
    free(cir);
    free(ptsx);
    free(ptsy);
    free(style);
}
/**
   Rename the log files when simulation exits.
 */
void rename_file(int sig){
    char fnnew[256];
    char fnold[256];
    char suffix[16];
    int pid=getpid();
    switch(sig){
    case 0:
	sprintf(suffix,"done");
	break;
    case SIGBUS:
    case SIGILL:
    case SIGSEGV:
    case SIGABRT:
	sprintf(suffix,"err");
	break;
    case SIGKILL:
    case SIGINT: //Ctrl-C
    case SIGTERM:
    case SIGQUIT: //Ctrl-'\'
	sprintf(suffix,"killed");
	break;
    default:
	sprintf(suffix,"unknown");
    }
    snprintf(fnnew,256,"kill_%d",pid);
    if(exist(fnnew)) remove(fnnew);
    
    snprintf(fnold,256,"run_%d.log",pid);
    snprintf(fnnew,256,"run_%d.%s", pid,suffix);
    rename(fnold,fnnew);
    mysymlink(fnnew, "run_recent.log");
}
/**
   Handles signals. We don't want to exit the simulation when SIGPIPE happens
   (when we are writing to closed sockets)
 */
void maos_signal_handler(int sig){
    if(sig==SIGPIPE){
	warning3("Program received signal SIGPIPE, broken pipe.\n");
	return;
    }
    disable_signal_handler;
    rename_file(sig);//handles signal
    if(sig!=0){
	info2("Caught signal %d\n",sig);
	if(sig == SIGSEGV){
	    print_backtrace(0);
	}
	scheduler_finish(1);
	exit(sig);
    }
}
/**
   Print out usage information.
 */
static void print_usage(void){
    info2(
"Usage: maos [OPTION...] [FILE]...\n"
"maos is a simulation tool developed to simulate scao/mcao systems\n\n"
"Examples:\n"
"maos   # Run the default configuration of NFIRAOS: nfiaros.conf as the baseline.\n"
"maos -c scao_ngs.conf -s 2 -n 2 -d -o scao_ngs override.conf chol.conf\n"
"       # Run a single conjugate natural guide star case, with seed 2, 2 threads\n"
"       # detach from the terminal and output results to folder scao_ngs\n"
"       # and read in overriding parameters stored in override.conf and chol.conf\n"
"\n"
"Options: \n"
"-h, --help        to print out this message\n"
"-d, --detach      to detach from terminal and run in background\n"
"-f, --force       force starting simulation without scheduler\n"
"-n N, --nthread=N Use N threads, default is 1\n"
"-s N, --seed=N    Use seed N instead of the numbers specified in .conf files\n"
"-o DIR, --output=DIR\n"
"                  output the results to DIR.\n"
"-c FILE.conf, --conf=FILE.conf\n"
"                  Use FILE.conf as the baseline config instead of nfiraos.conf\n"
	  );
}
/**
   Parse command line arguments argc, argv
 */
ARG_T * parse_args(int argc, char **argv){
    ARG_T *arg=calloc(1, sizeof(ARG_T));
    arg->nthread=1;//default is 1 thread.
    static struct option long_options[]={
	{"help",0,0,'h'},
	{"detach",0,0,'d'},
	{"force",0,0,'f'},
	{"output",1,0,'o'},
	{"nthread",1,0,'n'},
	{"conf",1,0,'c'},
	{"seed",1,0,'s'},
	{"path",1,0,'p'},
	{NULL,0,0,0}
    };
    while(1){
	int option_index = 0;
	int c = getopt_long(argc, argv, "hdfo:n:c:s:p:",
                        long_options, &option_index);
	if(c==-1) break;
	switch(c){
	case 'h':
	    print_usage();
	    exit(0);
	    break;
	case 'd':
	    arg->detach=1;
	    break;
	case 'f':
	    arg->force=1; 
	    break;
	case 'o':
	    if(arg->dirout){
		error("Duplicate argument for dirout\n");
	    }
	    arg->dirout=strdup(optarg);
	    break;
	case 'n':
	    arg->nthread=strtol(optarg,NULL,10);
	    if(arg->nthread<=0){
		warning("illigal nthread. set to 1.\n");
		arg->nthread=1;
	    }else if(arg->nthread>NCPU){
		warning("nthread is larger than number of cpus, reset to %d\n",
			NCPU);
		arg->nthread=NCPU;
	    }
	    break;
	case 'c':
	    if(arg->conf){
		error("Multiple -c switches found\n");
	    }
	    arg->conf=strdup(optarg);
	    break;
	case 's':{
	    arg->nseed++;
	    arg->seeds=realloc(arg->seeds, sizeof(int)*arg->nseed);
	    arg->seeds[arg->nseed-1]=strtol(optarg,NULL,10);
	    //info2("Command line supplied seed %d\n",arg->seeds[arg->nseed-1]);
	}
	    break;
	case 'p':{
	    addpath(optarg);
	}
	case '?':
	    warning("Unregonized option. exit.\n");exit(1);
	    break;
	default:
	    printf("?? getopt returned 0%o ??\n", c);exit(1);
	    break;
	}
    }
    if(!arg->detach){//foreground task will start immediately.
	arg->force=1;
    }
    arg->iconf=optind;
    arg->argc=argc;
    arg->argv=argv;
    if(!arg->dirout){
	arg->dirout=strtime();
    }
    
    if(!arg->conf){
	/*If -c is not specifid in path, will use nfiraos.conf*/
	arg->conf=strdup("nfiraos.conf");
    }
    info2("Main config file is %s\n",arg->conf);
    //Setup PATH and result directory
    char *config_path=find_config("maos");
    char *bin_path=stradd(config_path, "/bin", NULL);
    addpath(config_path);
    addpath(bin_path);
    free(bin_path);
    free(config_path);
    addpath(".");
    mymkdir("%s",arg->dirout);
    if(chdir(arg->dirout)){
	error("Unable to chdir to %s\n", arg->dirout);
    }
    info2("Output folder is '%s' %d threads\n",arg->dirout, arg->nthread);

    return arg;
}
/**
   Computes strehl from OPD without doing FFT. The strehl is simply 

   \f$s=\sum(A*exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]) \f$

   where A is the amplitude map.
 */
ccell *strehlcomp(const dmat *iopdevl, const double *amp,const int nwvl, const double *wvl){
    ccell *psf2s=ccellnew(nwvl,1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	dcomplex i2pi=I*2*M_PI/wvl[iwvl];
	dcomplex strehl=0;
	for(int iloc=0; iloc<iopdevl->nx; iloc++){
	    strehl+=amp[iloc]*cexp(i2pi*iopdevl->p[iloc]);
	}
	psf2s->p[iwvl]=cnew(1,1);
	psf2s->p[iwvl]->p[0]=strehl;
    }
    return psf2s;
}
/**
   Computes PSF from OPD by doing FFT. The PSF is computed as

   \f$\textrm{PSF}=\mathcal{F}[A\times exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]]\f$

   The peak value (center) in the computed PSF is normalized by the peak value
   in the differaction limited PSF. In other words, the peak value in the
   computed PSF is the Strehl. Keep this in mind when you compute enclosed
   energy.
 */
ccell *psfcomp(const dmat *iopdevl, const double *amp,
	       int **embeds, const int *nembeds, const int psfsize,
	       const int nwvl, const double *wvl){
    if(psfsize==1){
	return strehlcomp(iopdevl, amp, nwvl, wvl);
    }
    ccell *psf2s=ccellnew(nwvl,1);

    for(int iwvl=0; iwvl<nwvl; iwvl++){
	int nembed=nembeds[iwvl];
	int *embed=embeds[iwvl];
	//warning("Output PSF for direction %d\n",ievl);
	cmat *psf2=cnew(nembed,nembed);
	int ncomp;
	int use1d;
	int use1d_enable=0;
	if(nembed<psfsize){
	    ncomp=nembed;
	    if(use1d_enable){
		use1d=1;
		cfft2partialplan(psf2, ncomp, -1);
	    }else{
		use1d=0;
		cfft2plan(psf2, -1);
	    }
	}else{
	    ncomp=psfsize;
	    use1d=0;
	    cfft2plan(psf2, -1);
	}
	//PCMAT(psf2, ppsf2);
	//double psfnorm=1./(sqrt(aper->sumamp2)*aper->nembed);//This makes sum(psf)=1;
	double psfnorm=1; /*since sum(aper->amp)=1, this makes psf normalized by
			    differaction limited PSF. max(psf) is strehl.*/
	dcomplex i2pi=I*2*M_PI/wvl[iwvl];
	psf2s->p[iwvl]=cnew(ncomp,ncomp);

	czero(psf2);
	for(int iloc=0; iloc<iopdevl->nx; iloc++){
	    psf2->p[embed[iloc]]=amp[iloc]*cexp(i2pi*iopdevl->p[iloc]);
	    /*
	    int ix=embed[iloc][0];
	    int iy=embed[iloc][1];
	    ppsf2[iy][ix]=amp[iloc]*cexp(i2pi*iopdevl->p[iloc]);
	    */
	}
	if(use1d==1){
	    cfft2partial(psf2,ncomp,-1);
	}else{
	    cfft2(psf2,-1);
	}
	ccpcorner2center(psf2s->p[iwvl], psf2);
	
	if(fabs(psfnorm-1)>1.e-15)
	    cscale(psf2s->p[iwvl], psfnorm);
	cfree(psf2);
    }
    return psf2s;
}
/**
   Simple embed and accumulation
 */
void embed_in(double *out, const double *in, long nin, long *embed){
    for(long i=0; i<nin; i++){
	out[embed[i]]+=in[i];
    }
}
/**
   Simple embed and accumulation
 */
void embed_out(const double *out, double *in, long nin, long *embed){
    for(long i=0; i<nin; i++){
	in[i]+=out[embed[i]];
    }
}

/**
   Try to lock file for each seed. If file is already locked by other process,
   will not run that seed tox not step over other process's feet.
*/
int lock_seeds(PARMS_T *parms){
    char fn[80];
    int to_run=0;
    parms->fdlock=calloc(parms->sim.nseed, sizeof(int));
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	snprintf(fn, 80, "Res_%d.lock",parms->sim.seeds[iseed]);
	parms->fdlock[iseed]=lock_file(fn, 0);
	if(parms->fdlock[iseed]<0){
	    warning("Another MAOS is already running with seed %d. Skip\n",
		    parms->sim.seeds[iseed]);
	}else{
	    to_run++;
	}
    }
    return to_run;
}
/**
   Estimate anisoplanatic angle theta0 from Fried parameter r0, layer height and
   weights.  */
double calc_aniso(double r0, int nht, double *ht, double *wt){
    double wh=0;
    for(int iht=0; iht<nht; iht++){
	wh+=pow(ht[iht],5./3.)*wt[iht];
    }
    return 0.3144*r0*pow(wh,-3./5.);
}
/**
   Estimate generalized aniso angle theta2 from Fried parameter r0, and layer
   height and weights, and deformable mirror conjugation heights hc1 hc2 of the
   ground and altitude DMs. */
double calc_aniso2(double r0, int nht, double *ht, double *wt, double hc1, double hc2){
    double wh=0;
    double hh=pow(hc2-hc1,5./3.);
    for(int iht=0; iht<nht; iht++){
	double t1=0.5*pow(fabs(ht[iht]-hc1),5./3.)+0.5*pow(fabs(ht[iht]-hc2),5./3.);
	double t2=-0.25*hh-0.25/hh*pow(pow(fabs(ht[iht]-hc1),5./3.)-pow(fabs(ht[iht]-hc2),5./3.),2);
	wh+=wt[iht]*(t1+t2);
    }
    return 0.3144*r0*pow(wh,-3./5.);
}