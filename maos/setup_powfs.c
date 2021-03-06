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

#include "common.h"
#include "setup_powfs.h"
#include "genseotf.h"
#include "pywfs.h"
#include "setup_recon.h"
#include "recon_utils.h"
#define MOVES(p,i,j) p[i]=p[j]
#define MOVED(p,n,i,j) memcpy(p+n*i, p+n*j, sizeof(real)*n)
#define MOVEPTS(pts,count,isa)			\
    if(pts){					\
	MOVES(pts->origx, count,isa);		\
	MOVES(pts->origy, count,isa);		\
    }
#define MOVELOC(loc,count,isa)			\
    if(loc){					\
	MOVES(loc->locx,count,isa);		\
	MOVES(loc->locy,count,isa);		\
    }
#define MOVEDLOC(loc,nxsa,count,isa)		\
    if(loc){					\
	MOVED(loc->locx,nxsa,count,isa);	\
	MOVED(loc->locy,nxsa,count,isa);	\
    }

/**
   Free the powfs geometric parameters
*/
static void
free_powfs_geom(POWFS_T *powfs, int ipowfs){
    ptsfree(powfs[ipowfs].pts);
    dfree(powfs[ipowfs].saa);
    locfree(powfs[ipowfs].saloc);
    locfree(powfs[ipowfs].loc);
    dfree(powfs[ipowfs].amp);
    cellfree(powfs[ipowfs].loc_tel);
    dcellfree(powfs[ipowfs].amp_tel);
    dcellfree(powfs[ipowfs].saa_tel);
    cellfree(powfs[ipowfs].loc_dm);
    dcellfree(powfs[ipowfs].realamp);
    dcellfree(powfs[ipowfs].realsaa);
    dfree(powfs[ipowfs].sumamp);
    dfree(powfs[ipowfs].sumamp2);
    locfft_free(powfs[ipowfs].fieldstop);
}
/**
   Convert amplitude map of subapertures to normalized subaperture area. Used only in setup_powfs_geom*/
static dmat *wfsamp2saa(dmat *wfsamp, long nxsa){
    real areanormfactor=1./(real)nxsa;
    long nsa=wfsamp->nx/nxsa;
    dmat* saa=dnew(nsa, 1);
    for(int isa=0; isa<nsa; isa++){
	real *amp=wfsamp->p+isa*nxsa;
	real area=0;
	for(int iamp=0; iamp<nxsa; iamp++){
	    area+=amp[iamp]*amp[iamp];
	}
	saa->p[isa]=area*areanormfactor;
    }
    return saa;
}
/**
   Creates WFS pupil mask.
 */
void wfspupmask(const PARMS_T *parms, loc_t *loc, dmat *amp, int iwfs){
    long nloc=loc->nloc;
    dmat *ampmask=dnew(nloc, 1);
    real ht=parms->atm.hmax*0.7;
    for(int jwfs=0; jwfs<parms->nwfs; jwfs++){
	int jpowfs=parms->wfs[jwfs].powfs;
	if(parms->powfs[jpowfs].lo) continue;
	real hs=parms->wfs[iwfs].hs;
	real r=parms->aper.d*0.5*(1.-ht/hs)/(1.-ht/hs);
	real sx=(parms->wfs[jwfs].thetax-parms->wfs[iwfs].thetax)*ht;
	real sy=(parms->wfs[jwfs].thetay-parms->wfs[iwfs].thetay)*ht;
	loccircle(ampmask->p, loc, sx,sy, r, 1);
    }
    for(int i=0; i<nloc; i++){
	if(ampmask->p[i]<0.5) amp->p[i]=0;
    }
    dfree(ampmask);
}
static void
sa_reduce(POWFS_T *powfs, int ipowfs, real thresarea){
    dmat *saa=NULL;//create a temporary sa area array
    if(powfs[ipowfs].saa_tel){
	//We expact subaperture to union of all wfs
	long nsa=powfs[ipowfs].saa_tel->p[0]->nx;
	saa=dnew(nsa,1);
	for(long iwfs=0; iwfs<powfs[ipowfs].saa_tel->nx; iwfs++){
	    for(long isa=0; isa<nsa; isa++){
		if(P(saa,isa)<powfs[ipowfs].saa_tel->p[iwfs]->p[isa]){
		    P(saa,isa)=powfs[ipowfs].saa_tel->p[iwfs]->p[isa];
		}
	    }
	}
    }else{
	saa=ddup(powfs[ipowfs].saa);
    }
    if(dmax(saa)>1.01){
	warning("The sa area maxes to %g, which should be leq 1 (misregistration can cause this).\n", 
		dmax(saa));
    }
   
    if(powfs[ipowfs].pts->nsa>4){
	loc_t *ptsloc=(loc_t*)powfs[ipowfs].pts;
	loc_create_map(ptsloc);
	real dx1=1./ptsloc->dx;
	real dy1=1./ptsloc->dy;
	int changed=0;
	const real thresarea2=0.1;//secondary threshold to enable being surrounded subaperturs.
	do{//validata subapertures that have enough neighbors. Disable isolated subapertures
	    changed=0;
	    for(long isa=0; isa<ptsloc->nloc; isa++){
		long ix=round((ptsloc->locx[isa]-ptsloc->map->ox)*dx1);
		long iy=round((ptsloc->locy[isa]-ptsloc->map->oy)*dy1);
		int nedge=0;
		int ncorner=0;
		int nself=0;
		if(saa->p[isa]>=thresarea){
		    nself=1;
		}
		for(int jx=-1; jx<2; jx++){
		    for(int jy=-1; jy<2; jy++){
			if (jx==0 && jy==0){
			    continue;
			}
			long jsa=loc_map_get(ptsloc->map, ix+jx, iy+jy);
			if(jsa && saa->p[jsa-1]>=thresarea){
			    if(abs(jx+jy)==1){//edge
				nedge++;
			    }else{//corner
				ncorner++;
			    }
			}
		    }
		}
		if(nself){//disable isolated valid subaperture
		    if(ncorner+nedge<=2){
			saa->p[isa]=0;
			changed++;
		    }
		}else if(0){//enable isolated in-valid subaperture
		    if(nedge+ncorner>=6 && saa->p[isa]>=thresarea2){
			saa->p[isa]=1;
			changed++;
		    }
		}
	    }
	}while(changed);
	loc_free_map(ptsloc);
    }
    
    int count=0;
    const int nxsa=powfs[ipowfs].pts->nx * powfs[ipowfs].pts->nx;
    for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	if(saa->p[isa]>=thresarea){
	    /*Area is above threshold, keep.  Shift pts, ptsm, loc, locm, amp,
	      ampm, saloc area is already normalized that maxes to 1. The MOVE*
	      are defined in the beginining of this file.*/
	    if(count!=isa){
		MOVEPTS(powfs[ipowfs].pts,  count, isa);
		MOVES(powfs[ipowfs].saa->p, count, isa);
		MOVELOC(powfs[ipowfs].saloc, count, isa);
		MOVEDLOC(powfs[ipowfs].loc, nxsa, count, isa);
		MOVED(powfs[ipowfs].amp->p, nxsa, count, isa);
		if(powfs[ipowfs].saa_tel){
		    for(int jwfs=0; jwfs<powfs[ipowfs].saa_tel->nx; jwfs++){
			MOVES(powfs[ipowfs].saa_tel->p[jwfs]->p, count, isa);
			MOVEDLOC(powfs[ipowfs].loc_tel->p[jwfs], nxsa, count, isa);
			MOVED(powfs[ipowfs].amp_tel->p[jwfs]->p,nxsa, count, isa);
		    }
		}
	    }
	    count++;
	}
    }
    if(count==0){
	error("there are no subapertures above threshold.\n");
    }
    dfree(saa);

    ptsresize(powfs[ipowfs].pts, count);
    locresize(powfs[ipowfs].saloc, count);
    locresize(powfs[ipowfs].loc, count*nxsa);
    dresize(powfs[ipowfs].saa, count, 1);
    dresize(powfs[ipowfs].amp, count*nxsa, 1);
    powfs[ipowfs].saasum=dsum(powfs[ipowfs].saa);
    if(powfs[ipowfs].loc_tel){
	for(int jwfs=0; jwfs<powfs[ipowfs].saa_tel->nx; jwfs++){
	    locresize(powfs[ipowfs].loc_tel->p[jwfs], count*nxsa);
	    dresize(powfs[ipowfs].amp_tel->p[jwfs], count*nxsa, 1);
	    dresize(powfs[ipowfs].saa_tel->p[jwfs], count, 1);
	}
    }	
}
/**
   setting up subaperture geometry.
   
   - powfs->pts: the location of the lower left grid coordinate in this subaperture >=saloc.
   
   - powfs->saloc: the lower left corner of the subaperture
   
   - powfs->loc: the x,y coordinates of all the points, ordered acoording to subaperture.
   
   2011-02-11: In the matched filter and geometric gradient computation, the
   amplitude should use the real (misregistered) amplitude map. The geometric
   optics operator should just the normal loc but misregistered amplitude map
   (on locm).  The subaperture selector should, however, use assumed amp
   (non-misregistered).
   
*/
static void 
setup_shwfs_geom(POWFS_T *powfs, const PARMS_T *parms, 
		 APER_T *aper, int ipowfs){
    free_powfs_geom(powfs,ipowfs);
    int nwfsp=parms->powfs[ipowfs].nwfs;
    /*order of the system. 60 for TMT */
    const int order  = parms->powfs[ipowfs].order;
    /*Subaperture lateral length. 0.5 for TMT. */
    const real dsa = parms->powfs[ipowfs].dsa;
    /*The subaperture area that are considered as full. */
    real areafulli;
    if(order>2){
	areafulli=1.;
    }else{
	areafulli=4./M_PI;
    }

    /*Number of OPD pixels in 1 direction in each
      subaperture. Make it even to do fft.*/
    int nx = 2*(int)round(0.5*dsa/parms->powfs[ipowfs].dx);
    const real dx=dsa/nx;/*adjust dx. */
    const real dxoffset=dx*0.5;//Always keep points inside subaperture for simulation.
    if(fabs(parms->powfs[ipowfs].dx-dx)>EPS)
	info2("Adjusting dx from %g to %g\n", parms->powfs[ipowfs].dx,dx);
    if(fabs(dsa - nx * dx)>EPS){
	warning("nx=%d,dsa=%f,dx=%f not agree\n", nx, dsa, dx);
    }
    info2("There are %d points in each subaperture of %gm.\n", nx, dsa);
    const int nxsa=nx*nx;/*Total Number of OPD points. */
    if(parms->powfs[ipowfs].saloc){
	powfs[ipowfs].saloc=locread("%s", parms->powfs[ipowfs].saloc);
	if(fabs(powfs[ipowfs].saloc->dx-dsa)>1e-6*fabs(dsa)){
	    error("loaded saloc has dx=%g, while powfs.dsa=%g\n",
		  powfs[ipowfs].saloc->dx, dsa);
	}
    }else{
	/*The coordinate of the subaperture (lower left coordinate) */
	powfs[ipowfs].saloc=locnew(order*order, dsa, dsa);
	int count = 0;
	/*Offset of the coordinate of the center most subaperture from the center. */
	real offset;
	if(order & 1){/*odd */
	    offset = -0.5;
	}else{
	    offset = 0.0;
	}
	/*r2max: Maximum distance^2 from the center to keep a subaperture */
	real r2max=pow(order*0.5, 2);
	real r2min=dsa<parms->aper.din?pow(parms->aper.din/dsa/2,2):-1;
	/*the lower left *grid* coordinate of the subaperture */

	/*Collect all the subapertures that are within the allowed radius*/
	for(int j=-order/2; j<=(order-1)/2; j++){
	    for(int i=-order/2; i<=(order-1)/2; i++){
		//Normalized coordinate in uniq of sa size
		real xc=((real)i+offset);
		real yc=((real)j+offset);
		//Radius of four corners.
		real r1=pow(xc+1,2)+pow(yc+1,2);
		real r2=pow(xc+1,2)+pow(yc,2);
		real r3=pow(xc,2)+pow(yc+1,2);
		real r4=pow(xc,2)+pow(yc,2);
		if(r1<r2max || r2<r2max || r3<r2max || r4<r2max||
		   r1>r2min || r2>r2min || r3>r2min || r4>r2min){
		    powfs[ipowfs].saloc->locx[count]=xc*dsa;
		    powfs[ipowfs].saloc->locy[count]=yc*dsa;
		    count++;
		}
	    }
	}
	powfs[ipowfs].saloc->nloc=count;
    }
    /*convert saloc to pts*/
    powfs[ipowfs].pts=ptsnew(powfs[ipowfs].saloc->nloc, dsa, dsa, nx, dx, dx);
    for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	powfs[ipowfs].pts->origx[isa]=powfs[ipowfs].saloc->locx[isa]+dxoffset;
	powfs[ipowfs].pts->origy[isa]=powfs[ipowfs].saloc->locy[isa]+dxoffset;
    }
    /*Calculate the amplitude for each subaperture OPD point by ray tracing from
      the pupil amplitude map. Pupil distortion is accounted for.*/
    powfs[ipowfs].loc=pts2loc(powfs[ipowfs].pts);
    /*The assumed amp. */
    if(parms->powfs[ipowfs].amp){
	powfs[ipowfs].amp=dread("%s", parms->powfs[ipowfs].amp);
	if(powfs[ipowfs].amp->nx!=powfs[ipowfs].loc->nloc){
	    error("%s is in wrong format. Need %ld, has %ld.\n", 
		  parms->powfs[ipowfs].amp, powfs[ipowfs].loc->nloc, powfs[ipowfs].amp->nx);
	}
    }else{
	powfs[ipowfs].amp=mkamp(powfs[ipowfs].loc, aper->ampground, 
				-parms->misreg.pupil->p[0],-parms->misreg.pupil->p[1], 
				parms->aper.d, parms->aper.din);
    }
    /*The threashold for normalized area (by areafulli) to keep subaperture. */
    real thresarea=parms->powfs[ipowfs].saat;
    dmat *ampi=NULL;
    if(parms->powfs[ipowfs].safill2d<1){
	/*subaperture amplitude map to simulate lenslet fill factor*/
	int nedge=1;
	while ((nx-2*nedge)*(nx-2*nedge)>nx*nx*parms->powfs[ipowfs].safill2d){
	    nedge++;
	}
	ampi=dnew(nx,nx);
	real alpha=(nx*nx*parms->powfs[ipowfs].safill2d-(nx-2*nedge)*(nx-2*nedge))
	    /((nx-2*(nedge-1))*(nx-2*(nedge-1))-(nx-2*nedge)*(nx-2*nedge));
	real tot=0;
	for(int iy=nedge-1; iy<nx-nedge+1; iy++){
	    for(int ix=nedge-1; ix<nx-nedge+1; ix++){
		P(ampi,ix,iy)=1;
		if(ix==nedge-1 || ix==nx-nedge || iy==nedge-1 || iy==nx-nedge){
		    P(ampi,ix,iy)=alpha;
		}
		tot+=P(ampi,ix,iy);
	    }
	}
	if(parms->save.setup){
	    writebin(ampi,"powfs%d_ampi", ipowfs);
	}
	for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	    for(int i=0; i<nx*nx; i++){
		powfs[ipowfs].amp->p[nx*nx*isa+i]*=ampi->p[i];
	    }
	}
	/*do not multiply to siglev. Already handled automatically*/
	thresarea*=parms->powfs[ipowfs].safill2d;
    }

    powfs[ipowfs].saa=wfsamp2saa(powfs[ipowfs].amp, nxsa);

   
    setup_powfs_misreg_tel(powfs, parms, aper, ipowfs);
    /*Go over all the subapertures, calculate the normalized
      subaperture illumination area and remove all that are below
      the are threshold*/

    /*
      About physical optics imaging: The subaperture area (saa) has been
      normalized against a square subaperture to scale the physical optics
      i0. After FFT, the PSF sum to 1 for full square subapertures, but sum to
      PI/4 for full circular TT or full quadrant circular TTF subapertures. The
      areascale variable (1 for square subaperture and 4/M_PI for (quad)circular
      sa) is therefore used to normalize the PSF so that it sums to 1 for either
      full square or full (quadrant) circular subaperture. Used in wfsints.

      The subaperture area (saa) is subsequently scaled so that the max is 1 for
      full square or (quadrant) circular subaperture. saa is then used in
      setup_recon to scale geometric sanea.
    */

    powfs[ipowfs].areascale=areafulli;
    if(fabs(areafulli-1)>EPS){
	dscale(powfs[ipowfs].saa, areafulli);
	dcellscale(powfs[ipowfs].saa_tel, areafulli);
    }
    if(!parms->powfs[ipowfs].saloc && !parms->powfs[ipowfs].amp){
	sa_reduce(powfs, ipowfs, thresarea);
    }
    info("There are %ld valid subaperture.\n", powfs[ipowfs].saloc->nloc);
    setup_powfs_misreg_dm(powfs, parms, aper, ipowfs);
    powfs[ipowfs].realsaa=dcellnew(nwfsp, 1);
    for(int jwfs=0; jwfs<nwfsp; jwfs++){
	if(powfs[ipowfs].loc_tel){
	    powfs[ipowfs].realsaa->p[jwfs]=dref(powfs[ipowfs].saa_tel->p[jwfs]);
	}else{
	    powfs[ipowfs].realsaa->p[jwfs]=dref(powfs[ipowfs].saa);
	}
    }
    powfs[ipowfs].realamp=dcellnew(nwfsp, 1);
    powfs[ipowfs].sumamp=dnew(nwfsp, 1);
    powfs[ipowfs].sumamp2=dnew(nwfsp, 1);
    for(int jwfs=0; jwfs<nwfsp; jwfs++){
	dmat *realamp;
	if(powfs[ipowfs].loc_tel){
	    realamp=powfs[ipowfs].amp_tel->p[jwfs];
	}else{
	    realamp=powfs[ipowfs].amp;
	}
	real sumamp2=0;
	real sumamp=0;
	for(long i=0; i<powfs[ipowfs].loc->nloc; i++){
	    sumamp2+=realamp->p[i]*realamp->p[i];
	    sumamp+=realamp->p[i];
	}
	powfs[ipowfs].sumamp2->p[jwfs]=sumamp2;
	powfs[ipowfs].sumamp->p[jwfs]=sumamp;
	powfs[ipowfs].realamp->p[jwfs]=dref(realamp);
    }
    if(parms->powfs[ipowfs].fieldstop){
	warning("powfs%d: generating field stop \n", ipowfs);
	if(parms->powfs[ipowfs].nwvl>1){
	    error("Not implemented yet. need to do phase unwrap in wfsgrad.\n");
	}
	powfs[ipowfs].fieldstop=locfft_init(powfs[ipowfs].loc, powfs[ipowfs].amp, 
					    parms->powfs[ipowfs].wvl, NULL, 2,  
					    parms->powfs[ipowfs].fieldstop);
    }
    dfree(ampi);
 
    if(parms->save.setup){
	locwrite((loc_t*)powfs[ipowfs].pts,"powfs%d_pts",ipowfs);
	locwrite(powfs[ipowfs].saloc,"powfs%d_saloc",ipowfs); 
	writebin(powfs[ipowfs].saa,"powfs%d_saa",ipowfs);
	locwrite(powfs[ipowfs].loc,"powfs%d_loc",ipowfs);
	writebin(powfs[ipowfs].amp,"powfs%d_amp",ipowfs);
	if(powfs[ipowfs].loc_tel){
	    writebin(powfs[ipowfs].saa_tel,"powfs%d_saa_tel",ipowfs);
	    writebin(powfs[ipowfs].amp_tel,"powfs%d_amp_tel",ipowfs);
	    writebin(powfs[ipowfs].loc_tel,"powfs%d_loc_tel",ipowfs);
	}
	if(powfs[ipowfs].loc_dm){
	    for(int idm=0; idm<parms->ndm; idm++){
		for(int jwfs=0; jwfs<nwfsp; jwfs++){
		    writebin(powfs[ipowfs].loc_dm,"powfs%d_loc_dm", ipowfs);
		}
	    }
	}
    }
}
/**
   Setup telescope to WFS pupil misregistration.
   
   It is tricky to implement the misregistration of the deformable mirror and the WFS.
   
   Misregistration of the deformable mirror.  
   
   Denote the grid of actuators as saloc, and with misregistration, it becomes
   salocm. 
   
   - We should use saloc for DM fitting and pseudo open loop gradient
   computations because we are not aware of the misregistration. 

   - We should use salocm for raytracing to WFS and performance evaluation
   because this is what happens in the system, no matter whether we know the
   misregistration or not.

   Misregistration of the WFS.
   
   Denote the grid of the WFS entry pupil as loc, and with misregistration, it
   becomes locm. Project the telescope amplitude map onto loc, we get amp, which
   is what we know. Project the telescope amplitude map onto locm, we get ampm,
   which is what happens but we don't know.

   - We should use loc and amp to compute the reconstructor.

   - We should use locm and ampm for tracing to WFS and performance evaluation.
   
   - We should use loc and ampm for matched filter, averaging gradient, or
   zernike tip/tilt operator computation, since this is the process that happens
   in the system, and the matched filter is updated through dithering to the
   true amplitude map. We use loc because it is our model within the WFS.
   
   - We should use ptsm->area (real area) to select subapertures for use in
   reconstruction, because in real system, the subapertures are selected by
   their illumination level, and is affected by the real amplitude.


 */
void 
setup_powfs_misreg_tel(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    int nwfsp=parms->powfs[ipowfs].nwfs;
    if(parms->misreg.tel2wfs){
	TIC;tic;
	/*
	  Misregistration/distortion from Telescope pupil to WFS pupil. The
	  amplitude map after misregistration/distortion will be used for
	  wavefront sensing.
	  They are not used for wavefront reconstruction
	*/
	int isset=0;
	powfs[ipowfs].loc_tel=(loccell*)cellnew(nwfsp, 1);
	powfs[ipowfs].amp_tel=dcellnew(nwfsp, 1);
	if(parms->powfs[ipowfs].type==0){
	    powfs[ipowfs].saa_tel=dcellnew(nwfsp, 1);
	}
#pragma omp parallel for shared(isset)
	for(int jwfs=0; jwfs<nwfsp; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    if(parms->misreg.tel2wfs[iwfs]){
		isset=1;
		powfs[ipowfs].loc_tel->p[jwfs]
		    =loctransform(powfs[ipowfs].loc, parms->misreg.tel2wfs[iwfs]);
		powfs[ipowfs].amp_tel->p[jwfs]
		    =mkamp(powfs[ipowfs].loc_tel->p[jwfs], aper->ampground,
			   -parms->misreg.pupil->p[0], -parms->misreg.pupil->p[1], 
			   parms->aper.d, parms->aper.din);
		if(parms->powfs[ipowfs].type==0){
		    const int nxsa=powfs[ipowfs].pts->nx * powfs[ipowfs].pts->nx;
		    powfs[ipowfs].saa_tel->p[jwfs]=wfsamp2saa(powfs[ipowfs].amp_tel->p[jwfs], nxsa);
		}
	    }
	}
	if(!isset){
	    cellfree(powfs[ipowfs].loc_tel);
	    dcellfree(powfs[ipowfs].saa_tel); 
	    dcellfree(powfs[ipowfs].amp_tel);
	}else{
	    toc("misreg.tel2wfs");
	}
    }/*if misreg */ 
}
/**
   setup DM to WFS misregistration.
*/
void
setup_powfs_misreg_dm(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    (void)aper;
    int nwfsp=parms->powfs[ipowfs].nwfs;
    if(parms->misreg.dm2wfs){
	TIC;tic;
	/*
	  Misregistration/distortion from DM to WFS pupil. Not about telescope
	  pupil. The distorted grid are used for ray tracing from DM to WFS.
	  They are not used for wavefront reconstruction
	*/
	powfs[ipowfs].loc_dm=(loccell*)cellnew(nwfsp*parms->ndm, 1);
	int isset=0;
#pragma omp parallel for collapse(2) shared(isset)
	for(int idm=0; idm<parms->ndm; idm++){
	    for(int jwfs=0; jwfs<nwfsp; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		if(parms->misreg.dm2wfs[iwfs+idm*parms->nwfs]){
		    powfs[ipowfs].loc_dm->p[jwfs+nwfsp*idm]
			=loctransform(powfs[ipowfs].loc, parms->misreg.dm2wfs[iwfs+idm*parms->nwfs]);
		    isset=1;
		}
	    }
	}
	if(!isset){
	    cellfree(powfs[ipowfs].loc_dm);
	}else{
	    toc("misreg.dm2wfs");
	}
    }
}
/**
   Creating geometric wavefront gradient operator GS0 and ZA0 from WFS OPD to
   subaperture grads. Use loc, and ampm. Do not locm. locm is only used to do
   ray tracing.  This is simulation, not RTC, so use real system distortion,
   misregistration, etc. */
static void 
setup_shwfs_grad(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(parms->powfs[ipowfs].gtype_recon==0 ||parms->powfs[ipowfs].gtype_sim==0){
	dspcellfree(powfs[ipowfs].GS0);
	/*Setting up every gradient tilt (g-ztilt) */
	if(parms->load.GS0){
	    powfs[ipowfs].GS0=dspcellread("powfs%d_GS0",ipowfs);
	    if(powfs[ipowfs].amp_tel){
		assert(powfs[ipowfs].GS0->nx==parms->powfs[ipowfs].nwfs);
	    }else{
		assert(powfs[ipowfs].GS0->nx==1);
	    }
	}else{
	    /*This mkg takes about 5 seconds. */
	    if(powfs[ipowfs].amp_tel){
		powfs[ipowfs].GS0=dspcellnew(parms->powfs[ipowfs].nwfs, 1);
	    }else{
		powfs[ipowfs].GS0=dspcellnew(1, 1);
	    }
	    for(int iwfs=0; iwfs<powfs[ipowfs].GS0->nx; iwfs++){
		powfs[ipowfs].GS0->p[iwfs]=mkg(powfs[ipowfs].loc, 
					       powfs[ipowfs].loc,
					       powfs[ipowfs].realamp->p[iwfs],
					       powfs[ipowfs].saloc,
					       1, 0, 0, 1);
	    }
	    if(parms->save.setup && powfs[ipowfs].GS0){
		writebin(powfs[ipowfs].GS0,"powfs%d_GS0",ipowfs);
	    }
	}
    }
    if(parms->powfs[ipowfs].gtype_recon==1 ||parms->powfs[ipowfs].gtype_sim==1){
	/*setting up zernike best fit (ztilt) inv(M'*W*M). good for NGS. */
	if(parms->powfs[ipowfs].order>4){ 
	    warning("Ztilt for high order powfs %d is not good\n", ipowfs);
	}
	powfs[ipowfs].nsaimcc=MAX(1,(powfs[ipowfs].loc_tel?parms->powfs[ipowfs].nwfs:1));
	int nsaimcc=powfs[ipowfs].nsaimcc;
	cellfree(powfs[ipowfs].saimcc);
	powfs[ipowfs].saimcc=(dccell*)cellnew(nsaimcc, 1);
	for(int imcc=0; imcc<nsaimcc; imcc++){
	    dcell *mcc=pts_mcc_ptt(powfs[ipowfs].pts, powfs[ipowfs].realamp->p[imcc]->p);
	    powfs[ipowfs].saimcc->p[imcc]=dcellinvspd_each(mcc);
	    dcellfree(mcc);
	}
    }
}
void setup_powfs_neasim(const PARMS_T *parms, POWFS_T *powfs){
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	const long nsa=powfs[ipowfs].saloc->nloc;
	dcell *nea=0;
	//if(parms->powfs[ipowfs].neaphy || parms->powfs[ipowfs].phystep>-1){
	if(powfs[ipowfs].sanea){
	    info2("Use sanea to derive neasim\n");
	    nea=dcelldup(powfs[ipowfs].sanea);
	    for(int ii=0; ii<nea->nx; ii++){
		nea_chol(&nea->p[ii], nea->p[ii]);
	    }
	}else{
	    //neasimfile/neareconfile is saved by skyc in rad, not in rad^2.
	    char *file=parms->powfs[ipowfs].neasimfile; 
	    if(!file && parms->powfs[ipowfs].neasim==-1){
		file=parms->powfs[ipowfs].neareconfile;
	    }
	    if(file){
		nea=dcellread_prefix(file, parms, ipowfs);
	    }
	}
	if(nea){
	    for(int ii=0; ii<nea->nx*nea->ny; ii++){
		check_nea(nea->p[ii],nsa);
	    }
	}else{
	    int nnea=powfs[ipowfs].loc_tel?parms->powfs[ipowfs].nwfs:1;
	    nea=dcellnew(nnea, 1);
	    for(int jwfs=0; jwfs<nnea; jwfs++){
		real nea_rad;
		if(parms->powfs[ipowfs].neasim<0){
		    nea_rad=parms->powfs[ipowfs].nearecon;//in mas
		}else{
		    nea_rad=parms->powfs[ipowfs].neasim;//in mas
		}
		nea_rad=nea_rad/206265000./sqrt(parms->powfs[ipowfs].dtrat);//in rad
		real *saa=powfs[ipowfs].realsaa?powfs[ipowfs].realsaa->p[jwfs]->p:powfs[ipowfs].saa->p;
		dmat *nea_each=nea->p[jwfs]=dnew(nsa, 3);
		for(int isa=0; isa<nsa; isa++){
		    P(nea_each,isa,0)=P(nea_each,isa,1)=nea_rad/sqrt(saa[isa]);
		}
	    }
	}
	powfs[ipowfs].neasim=nea; nea=0;
	if(parms->save.setup){
	    writebin(powfs[ipowfs].neasim,"powfs%d_neasim", ipowfs);
	}
    }
}
/**
   Prepare the parameters for physical optics setup.

   \todo Make LGS number of pixels depend on subapertue elongation in radial
   coordinate CCD mode. Make notfx, notfy dependent on number of pixels for radial ccd.

*/
static void 
setup_powfs_prep_phy(POWFS_T *powfs,const PARMS_T *parms,int ipowfs){
    const real pixthetax=parms->powfs[ipowfs].radpixtheta;
    const real pixthetay=parms->powfs[ipowfs].pixtheta;
    const int pixpsay=parms->powfs[ipowfs].pixpsa;
    const int radpix=parms->powfs[ipowfs].radpix;
    const real dsa=powfs[ipowfs].pts->dsa;
    const int nsa=powfs[ipowfs].saloc->nloc;
    if(parms->powfs[ipowfs].llt){
	const int nllt=parms->powfs[ipowfs].llt->n;
	real rsa2, rsa2max=0;
	dcellfree(powfs[ipowfs].srot);
	dcellfree(powfs[ipowfs].srsa);
	dcellfree(powfs[ipowfs].sprint);

	powfs[ipowfs].srot=dcellnew(nllt,1);
	powfs[ipowfs].srsa=dcellnew(nllt,1);
	powfs[ipowfs].srsamax=dnew(nllt,1);
	powfs[ipowfs].sprint=dcellnew(nllt,1);

	for(int illt=0;illt<nllt;illt++){
	    /*adjusted llt center because saloc->locx/y is corner */
	    real ox2=parms->powfs[ipowfs].llt->ox->p[illt]-dsa*0.5;
	    real oy2=parms->powfs[ipowfs].llt->oy->p[illt]-dsa*0.5;
	    powfs[ipowfs].srot->p[illt]=dnew(nsa,1);
	    powfs[ipowfs].srsa->p[illt]=dnew(nsa,1);

	    for(int isa=0; isa<nsa; isa++){
		real ddx=(powfs[ipowfs].saloc->locx[isa]-ox2);
		real ddy=(powfs[ipowfs].saloc->locy[isa]-oy2);
		rsa2=pow(ddx,2)+pow(ddy,2);
		powfs[ipowfs].srot->p[illt]->p[isa]=atan2(ddy,ddx);
		powfs[ipowfs].srsa->p[illt]->p[isa]=sqrt(rsa2);
		if(rsa2>rsa2max) rsa2max=rsa2;
	    }
	    powfs[ipowfs].srsamax->p[illt]=sqrt(rsa2max);
	    real dprint=MAX(parms->aper.d*0.1, dsa);
	    int pnsa=(int)ceil(sqrt(rsa2max)/dprint)+1;
	    
	    real prot[pnsa];
	    powfs[ipowfs].sprint->p[illt]=dnew(pnsa,1);
	    real *pp=powfs[ipowfs].sprint->p[illt]->p;
	    for(int ind=0; ind<pnsa; ind++){
		prot[ind]=INFINITY;
		pp[ind]=-1;
	    }
	    real desrot=0;
	    if(fabs(parms->powfs[ipowfs].llt->ox->p[illt])<dsa 
	       && fabs(parms->powfs[ipowfs].llt->oy->p[illt])<dsa){
		desrot=0;
	    }else{
		real ddx=(0-parms->powfs[ipowfs].llt->ox->p[illt]);
		real ddy=(0-parms->powfs[ipowfs].llt->oy->p[illt]);
		desrot=atan2(ddy,ddx);
	    }
	    for(int isa=0; isa<nsa; isa++){
		int ind=(int)round(powfs[ipowfs].srsa->p[illt]->p[isa]/dprint);
		real irot=fabs(powfs[ipowfs].srot->p[illt]->p[isa]-desrot);
		if(ind>=pnsa){
		    error("ind=%d>=pnsa=%d\n", ind, pnsa);
		}
		if(irot<prot[ind]){
		    prot[ind]=irot;
		    pp[ind]=(real)isa;
		}
	    }
	}
	if(parms->save.setup){
	    writebin(powfs[ipowfs].sprint,"powfs%d_sprint",ipowfs);
	    writebin(powfs[ipowfs].srot,"powfs%d_srot",ipowfs);
	    writebin(powfs[ipowfs].srsa,"powfs%d_srsa",ipowfs);
	}
    }
    int pixpsax=pixpsay;
    if(radpix<0){/*determine radpix adaptively */
	/*use rsa2max to determine radpix */
	error("Not emplemented yet\n");
    }
    pixpsax=(radpix)?radpix:pixpsay;
    
    /*record # of detector pixels. */
    powfs[ipowfs].pixpsax=pixpsax;
    powfs[ipowfs].pixpsay=pixpsay;

    {
	/*to avoid aliasing, fft warping. usually 2. */
	const real embfac=parms->powfs[ipowfs].embfac;
 	real safov=MAX(pixpsax*pixthetax, pixpsay*pixthetay);
	real wvlmin=parms->powfs[ipowfs].wvl->p[0];
	real dtheta=wvlmin/(dsa*embfac);/*PSF sampling. */

	/*size required to form detector image. */
	int notf;
	const char *okind;
	if(parms->powfs[ipowfs].notf){
	    notf=(parms->powfs[ipowfs].notf+1)/2*2;
	    okind="input";
	}else{
	    notf=MAX(powfs[ipowfs].pts->nx*embfac, ceil(safov/dtheta));
	    /*
	    if(notf>8 && notf<16){
		notf=16;
	    }else if(notf>16 && notf<32){
		notf=32;
	    }else if(notf>64 && notf<67){
		notf=64;
	    }else if(notf<128 && notf>120){
		notf=128;
		}*/
	    notf=nextfftsize(notf);
	    okind="calculated";
	}
	powfs[ipowfs].notfx=notf;
	powfs[ipowfs].notfy=notf;
	info2("notf is %dx%d (%s)\n", powfs[ipowfs].notfx, powfs[ipowfs].notfy, okind);
	if(safov > dtheta*notf){
	    warning("Subaperture PSF size (%g\") is smaller than detector FoV (%g\").\n", 
		    dtheta*notf*206265, safov*206265);
	}
    }/*notf */


    if(parms->powfs[ipowfs].bkgrndfn){
	char *fn=parms->powfs[ipowfs].bkgrndfn;
	info2("Loading sky background/rayleigh backscatter from %s\n",fn);
	dcellfree(powfs[ipowfs].bkgrnd);
	powfs[ipowfs].bkgrnd=dcellread("%s",fn);
    }
    if(parms->powfs[ipowfs].bkgrndfnc){
	char *fn=parms->powfs[ipowfs].bkgrndfnc;
	info2("Loading sky background/rayleigh backscatter correction from %s\n",fn);
	dcellfree(powfs[ipowfs].bkgrndc);
	powfs[ipowfs].bkgrndc=dcellread("%s",fn);
    }
    if(parms->powfs[ipowfs].bkgrndfn || parms->powfs[ipowfs].bkgrndfnc){
	real bkscale = parms->sim.dt*800*parms->powfs[ipowfs].dtrat;
	if(fabs(bkscale-1)>1.e-20){
	    dcellscale(powfs[ipowfs].bkgrnd, bkscale);
	    dcellscale(powfs[ipowfs].bkgrndc, bkscale);
	    warning("Scaling bkgrnd by %g", bkscale);
	}
	if(parms->powfs[ipowfs].bkgrndfn && 
	   (powfs[ipowfs].bkgrnd->nx!=powfs[ipowfs].saloc->nloc
	    ||(powfs[ipowfs].bkgrnd->ny!=1
	       && powfs[ipowfs].bkgrnd->ny!=parms->powfs[ipowfs].nwfs))){
	    error("powfs%d: bkgrnd is of dimension %ld x %ld, "
		  "but should be %ld x 1 or %d\n",
		  ipowfs, powfs[ipowfs].bkgrnd->nx, powfs[ipowfs].bkgrnd->ny,
		  powfs[ipowfs].saloc->nloc, parms->powfs[ipowfs].nwfs);
	}
	if(parms->powfs[ipowfs].bkgrndfnc && 
	   (powfs[ipowfs].bkgrndc->nx!=powfs[ipowfs].saloc->nloc
	    ||(powfs[ipowfs].bkgrndc->ny!=1
	       && powfs[ipowfs].bkgrndc->ny!=parms->powfs[ipowfs].nwfs))){
	    error("powfs%d: bkgrndc is of dimension %ld x %ld, "
		  "but should be %ld x 1 or %d\n",
		  ipowfs, powfs[ipowfs].bkgrndc->nx, powfs[ipowfs].bkgrndc->ny,
		  powfs[ipowfs].saloc->nloc, parms->powfs[ipowfs].nwfs);
	}
    }
}

/**
   Setting up Detector transfer function used to translate PSFs from FFTs to
   detector pixel images. Integration over pixel size and leaking between pixels
   are considered.
   
   We are going to sample the PSF \f$\phi(\theta)\f$ onto detectors. For a detector
   pixel \f$i\f$ at \f$\theta_{i}\f$, oriented along angle \f$\alpha\f$, we have
   \f[
   I(\theta_{i})=\int\textrm{d}\theta\phi(\theta)S(\theta-\theta_{i})=\int\textrm{d}\theta\phi(\theta-\theta_{i})S(\theta)\f]
   where \f$S(\theta)\f$ is the pixel weighting function.  Rewrite
   the equation in Fourier domain, 

   \f{eqnarray*}{ I(\theta_{i}) & = &
   \int\textrm{d}\theta\int\textrm{d} u\hat{\phi}(u)e^{2\pi i(\theta-\theta_{i})u}S(\theta)\\ & =
   & \int\textrm{d} u\hat{\phi}(u)\hat{S}(-u)e^{2\pi i\theta_{i}u}\\ & = &
   \mathcal{F}^{-1}[\hat{\phi}(u)\hat{S}(-u)](\theta_{i})\f} 

   For a polar coordinate CCD with pixel edges alonging , we have

   The detector pixels can be modeled as a square, convolved with a Gaussian to
   model the pixel blurring (charge leaking)
   \f[
   S(\theta)=S^{\prime}(\theta_{ra})=\sqcap(\theta_{r}/\Delta)\sqcap(\theta_{a}/\Delta)*\frac{1}{2\pi
   \theta_b^2}\exp(-\frac{\theta_r^2+\theta_u^2}{2\theta_b^2})
   \f], where \f$R\f$ is the rotation matrix and 
   \f$\theta_{ra}=R^{T}\theta\f$ is the coordinate in the local radial-azimuthal
   coordinate system for a Polar coordinate CCD, \f$\theta_b\f$ is a measure of the blurring of the pixel. For a pixel aligned along \f$x/y\f$ directions, \f$R\f$ is simply identity.

   The Fourier transform of the detector pixel function \f$S(\theta)\f$ is then

   \f{eqnarray*}{ 
   \hat{S}(u) & = & \int\textrm{d}\theta S(\theta)e^{-2\pi i\theta^{T}u}\\ 
   & = & \int\textrm{d}\theta
   S^{\prime}(R^{T}\theta)e^{-2\pi i\theta^{T}u}\\ 
   & = &\int\textrm{d}(R^{T}\theta)S^{\prime}(R^{T}\theta)e^{-2\pi i(R^{T}\theta)^{T}(R^{T}u)}\\
   & = & \hat{S^{\prime}}(R^{T}u)\f}

   The Fourier transform of \f$\hat{S^{\prime}}\f$ is 
   \f[
   \hat{S^{\prime}}(u_{ra})=\Delta^{2}\textrm{sinc}(u_{r}\Delta)\textrm{sinc}(u_{\theta}\Delta)
   \exp[-2\pi^2\theta_b^2(u_r^2+u_a^2)]
   \f]

   Finally, we have \f[
   I(\theta_{i})=\mathcal{F}^{-1}[\hat{\phi}(u)\Delta^{2}\textrm{sinc}(\Delta(u_{x}\cos\theta+u_{y}\sin\theta))\textrm{sinc}(\Delta(-u_{x}\sin\theta+u_{y}\cos\theta))\exp[-2\pi^2\theta_b^2(u_r^2+u_a^2)](\theta_{i})\f]
   the last \f$\theta_{i}\f$ means evaluate the function at \f$\theta_{i}\f$ by
   interpolation. When there is leaking between actuators. We can multiply the
   sinc functions with a Gaussian blurring function.
*/
/*
   \section WFS gradient pixel offset
   
   Use powfs.pixoffx and powfs.pixoffy to set WFs gradient pix offset
   
   - If abs(pixoffx) is less than 1. All subapertures have uniform pixel offset along x/r: pixoffx and along y/a: pixoffy (pixel).
   
   - If pixoffx==1 : There is a rotational offset with maximum value of pixoffy (pixel) at the edge.

   - If pixoffx==2: There is a global offset along x and y of pixoffy (pixel).

*/
static void 
setup_powfs_dtf(POWFS_T *powfs,const PARMS_T *parms,int ipowfs){
    dmat *pixoffx=0;
    dmat *pixoffy=0;
    if(parms->powfs[ipowfs].pixoffx||parms->powfs[ipowfs].pixoffy){
	const int nsa=powfs[ipowfs].pts->nsa;
	if(fabs(parms->powfs[ipowfs].pixoffx)<1){
	    info2("powfs%d: uniform pixel offset\n", ipowfs);
	    //both pixoff within 1 denotes constant offset in unit of pixel.
	    pixoffx=dnew(nsa,1); dset(pixoffx, parms->powfs[ipowfs].pixoffx);
	    pixoffy=dnew(nsa,1); dset(pixoffy, parms->powfs[ipowfs].pixoffy);
	}else{
	    const int nwfs=1; //parms->powfs[ipowfs].nwfs;
	    pixoffx=dnew(nsa, nwfs);
	    pixoffy=dnew(nsa, nwfs);
	    if(parms->powfs[ipowfs].pixoffx==1){
		info2("powfs%d: CCD rotation wrt lenslet @ %g pixel.\n", ipowfs, parms->powfs[ipowfs].pixoffy);
		//pixoffx==1 denotes clocking effect with pixoffy denotes maximum amount in unit of pixel CCW.
		real pixtheta=0.5*(parms->powfs[ipowfs].pixtheta+parms->powfs[ipowfs].radpixtheta);
		const real dsah=powfs[ipowfs].pts->dsa*0.5;
		const real angle=parms->powfs[ipowfs].pixoffy*pixtheta/(parms->aper.d*0.5-dsah);
		const real ct=cos(angle);
		const real st=sin(angle);
	
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    for(int isa=0; isa<nsa; isa++){
			real ox=powfs[ipowfs].saloc->locx[isa]+dsah;
			real oy=powfs[ipowfs].saloc->locy[isa]+dsah;
			real dx=ox*ct-oy*st-ox;
			real dy=ox*st+oy*ct-oy;
			if(parms->powfs[ipowfs].radpix){//radial coordinate.
			    dy=sqrt(dx*dx+dy*dy);
			    dx=0;
			}
			P(pixoffx, isa, jwfs)=dx/pixtheta;
			P(pixoffy, isa, jwfs)=dy/pixtheta;
		    }
		}
	
	    }else if(parms->powfs[ipowfs].pixoffx==2){
		info2("powfs%d: CCD global shift wrt lenslet @ %g pixel.\n", ipowfs, parms->powfs[ipowfs].pixoffy);
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    for(int isa=0; isa<nsa; isa++){
			real gx=parms->powfs[ipowfs].pixoffy;
			real gy=0;
			if(parms->powfs[ipowfs].radpix){
			    real angle=PR(PR(powfs[ipowfs].srot, jwfs, 1), isa, 1);
			    real ct=cos(angle);
			    real st=sin(angle);
			    real gx2=gx*ct+gy*st;//XY->RA: CW
			    gy=-gx*st+gy*ct;
			    gx=gx2;
			}
			P(pixoffx, isa, jwfs)=gx;
			P(pixoffy, isa, jwfs)=gy;
		    }
		}
	    }else{
		error("Invalid input in pixoffx\n");
	    }
	}
	if(parms->save.setup>1){
	    writebin(pixoffx, "powfs%d_pixoffx", ipowfs);
	    writebin(pixoffy, "powfs%d_pixoffy", ipowfs);
	}
    }

    powfs[ipowfs].pixoffx=pixoffx;
    powfs[ipowfs].pixoffy=pixoffy;
    
    powfs[ipowfs].dtf=mkdtf(parms->powfs[ipowfs].wvl, 
			    powfs[ipowfs].pts->dsa,
			    parms->powfs[ipowfs].embfac,
			    powfs[ipowfs].notfx,
			    powfs[ipowfs].notfy,
			    powfs[ipowfs].pixpsax,
			    powfs[ipowfs].pixpsay,
			    parms->powfs[ipowfs].radpixtheta,
			    parms->powfs[ipowfs].pixtheta,
			    powfs[ipowfs].pixoffx,
			    powfs[ipowfs].pixoffy,
			    parms->powfs[ipowfs].pixblur,
			    powfs[ipowfs].srot,
			    parms->powfs[ipowfs].radpix);

    int nwvl=parms->powfs[ipowfs].nwvl;
    powfs[ipowfs].dtheta=dnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	powfs[ipowfs].dtheta->p[iwvl]=powfs[ipowfs].dtf[iwvl].dtheta;
	if(parms->save.setup>1){
	    writebin(powfs[ipowfs].dtf[iwvl].nominal,
		     "powfs%d_dtf%d_nominal",ipowfs,iwvl);
	    writebin(powfs[ipowfs].dtf[iwvl].si,
		     "powfs%d_dtf%d_si",ipowfs,iwvl);
	}
    }
}
/**
   setup the range to sodium layer as an additional parameter.
*/
static void setup_powfs_focus(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(!parms->powfs[ipowfs].llt || !parms->powfs[ipowfs].llt->fnrange) return;
    char *fnrange=parms->powfs[ipowfs].llt->fnrange;
    warning("loading sodium range variation from %s\n", fnrange);
    if(powfs[ipowfs].focus) dfree(powfs[ipowfs].focus);
    powfs[ipowfs].focus=dread("%s",fnrange);
    if(powfs[ipowfs].focus->ny!=1 &&powfs[ipowfs].focus->ny!=parms->powfs[ipowfs].nwfs){
	error("fnrange has wrong format. Must be column vectors of 1 or %d columns\n",
	      parms->powfs[ipowfs].nwfs);
    }
    /*(D/h)^2/(16*sqrt(3)) convert from range to WFE in m but here we want
      focus mode, so just do 1/(2*h^2).*/
    /*1./cos() is for zenith angle adjustment of the range.*/
    real range2focus=0.5*pow(1./parms->powfs[ipowfs].hs,2)*(1./cos(parms->sim.za));
    dscale(powfs[ipowfs].focus, range2focus);
    if(parms->save.setup){
	writebin(powfs[ipowfs].focus, "powfs%d_focus", ipowfs);
    }
}

/**
   Load and smooth out sodium profile and adjust for telescope altitude and
   zenith angle. We preserve the sum of the sodium profile, which represents a
   scaling factor of the signal level.
   
*/
static void setup_powfs_sodium(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    char *fnprof=parms->powfs[ipowfs].llt->fnprof;
    dcell *Nains=dcellread("%s",fnprof);
    int nprof=Nains->nx*Nains->ny;
    if(nprof!=1 && nprof!=parms->powfs[ipowfs].nwfs){
	error("The sodium profile input %s is in wrong fromat\n", fnprof);
    }
    powfs[ipowfs].sodium=dcellnew(nprof, 1);
    double nhtel=-parms->sim.htel;
    double secza=1./cos(parms->sim.za);
    for(int i=0; i<Nains->nx*Nains->ny; i++){
	dmat *Nain=Nains->p[i];
	if(Nain->ny<2 || Nain->nx!=Nains->p[0]->nx){
	    error("The sodium profile input %s is in wrong format\n", fnprof);
	}
	if(parms->dbg.na_smooth){/*resampling the sodium profile by binning. */
	    /*Make new sampling: */
	    const real rsamax=dmax(powfs[ipowfs].srsamax);
	    const real dthetamin=dmin(powfs[ipowfs].dtheta);
	    /*minimum sampling required. */
	    const real dxnew=pow(parms->powfs[ipowfs].hs,2)/rsamax*dthetamin;
	    powfs[ipowfs].sodium->p[i]=smooth(Nain, dxnew);
	}else{
	    info2("Not smoothing sodium profile\n");
	    powfs[ipowfs].sodium->p[i]=dref(Nain);
	}
	for(int ih=0; ih<powfs[ipowfs].sodium->p[i]->nx; ih++){
	    powfs[ipowfs].sodium->p[i]->p[ih]=(powfs[ipowfs].sodium->p[i]->p[ih]+nhtel)*secza;
	}
    }
    dcellfree(Nains);
    if(parms->save.setup){
	writebin(powfs[ipowfs].sodium,"powfs%d_sodium",ipowfs);
    }
}
typedef struct{
    DTF_T *dtfs;  /**<The dtfs*/
    real hs;    /**<Guide star focus range*/
    dcell *sodium;/**<The sodium profile. First column is coordinate.*/
    int icol;     /**<Which sodium profile to use*/
    dcell *srot;  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
    dcell *srsa;  /**<Subaperture to LLT distance*/
    int no_interp;/**<Use direct sum instead of interpolation + FFT. Slower */
    int free;     /**<Free this array after using?*/
}mketf_t;
ETF_T* mketf_wrap(mketf_t *data){
    ETF_T*result=mketf(data->dtfs, data->hs, data->sodium, data->icol,
		       data->srot, data->srsa, data->no_interp);
    if(data->free) free(data);
    return result;
}
/**
   Compute Elongation Transfer function.
   - mode=0: for preparation.
   - mode=1: for simulation.
*/
void setup_powfs_etf(POWFS_T *powfs, const PARMS_T *parms, double deltah, int ipowfs, int mode, int icol){
    if(!parms->powfs[ipowfs].llt) return;
    mketf_t etfdata={powfs[ipowfs].dtf, 
		     parms->powfs[ipowfs].hs+deltah,
		     powfs[ipowfs].sodium,
		     icol,
		     powfs[ipowfs].srot,
		     powfs[ipowfs].srsa,
		     !parms->dbg.na_interp, 0};
    if(mode==0){/*preparation. */
	if(powfs[ipowfs].etfprep && powfs[ipowfs].etfsim!=powfs[ipowfs].etfprep){
	    etf_free(powfs[ipowfs].etfprep);
	}
	powfs[ipowfs].etfprep=mketf_wrap(&etfdata);
    }else{/*simulation*/
	if(mode==1){/*first pair for interpolation*/
	    if(powfs[ipowfs].etfsim==powfs[ipowfs].etfprep){
		powfs[ipowfs].etfsim=0;
	    }
	    etf_free(powfs[ipowfs].etfsim); powfs[ipowfs].etfsim=0;
	    if(powfs[ipowfs].etfsim2 && icol==powfs[ipowfs].etfsim2->icol && deltah==0){//reuse etfsim2 as etfsim
		powfs[ipowfs].etfsim=powfs[ipowfs].etfsim2;
		powfs[ipowfs].etfsim2=0;
	    }else{
		powfs[ipowfs].etfsim=mketf_wrap(&etfdata);
	    }
	}else if(mode==2){/*second pair for interpolation*/
	    static pthread_t etfthread=0;
	    ETF_T *etfasync=0;
	    etf_free(powfs[ipowfs].etfsim2); powfs[ipowfs].etfsim2=0;
	    if(etfthread){//preparation already running in a thread
		pthread_join(etfthread, (void**)(void*)&etfasync);
		if(etfasync->icol!=icol){
		    warning("Async prepared etfsim2 (%d) is not correct (%d)\n",
			    etfasync->icol, icol);
		}else{
		    powfs[ipowfs].etfsim2=etfasync;
		}
		etfthread=0;
	    }
	    if(!powfs[ipowfs].etfsim2){//no async available
		powfs[ipowfs].etfsim2=mketf_wrap(&etfdata);
	    }
	    if(0){//Disable incase deltah is not zero.
		//asynchronously preparing for next update.
		//Copy data to heap so that they don't disappear during thread execution
		mketf_t *etfdata2=mycalloc(1,mketf_t);//freed by mketf_wrap.
		memcpy(etfdata2, &etfdata, sizeof(mketf_t));
		etfdata2->icol++;
		etfdata2->free=1;
		if(pthread_create(&etfthread, NULL, (void*(*)(void *))mketf_wrap, etfdata2)){
		    warning_once("Thread creation failed\n");
		    free(etfdata2);
		    etfthread=0;
		}
	    }
	}else{
	    error("Invalid mode=%d\n", mode);
	}
    }
    print_mem("After setup_powfs_etf");
}

/**
   setting up uplink pts/amp lotf
*/
static void 
setup_powfs_llt(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(!parms->powfs[ipowfs].llt) return;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    real wvl0=parms->powfs[ipowfs].wvl->p[0];
    LLT_T *llt=powfs[ipowfs].llt=mycalloc(1,LLT_T);
    const LLT_CFG_T *lltcfg=parms->powfs[ipowfs].llt;

    real lltd=lltcfg->d;
    real dx;
    if(1){//avoid extending LLT psf
	int notf=MAX(powfs[ipowfs].notfx, powfs[ipowfs].notfy);
	dx=parms->powfs[ipowfs].embfac*powfs[ipowfs].pts->dsa/notf;
    }else{
	dx=powfs[ipowfs].pts->dx;
    }
    if(lltd>powfs[ipowfs].pts->dsa){
	error("Please check the implementation for this case.\n");
    }
    real lltdsa=MAX(lltd, powfs[ipowfs].pts->dsa);
    int nx=round(lltdsa/dx); lltdsa=dx*nx;
    pts_t *lpts=llt->pts=ptsnew(1, lltdsa, lltdsa, nx, dx, dx);
    llt->pts->origx[0]=(dx-lltdsa)*0.5;
    llt->pts->origy[0]=(dx-lltdsa)*0.5;
    
    real sumamp2=0;
    llt->amp=dnew(nx*nx, 1);
    if(lltcfg->fnamp){
	map_t *lltamp=mapread("%s", lltcfg->fnamp);
	prop_grid_pts(lltamp, llt->pts, llt->amp->p, 1, 0, 0, 1, 1, 0, 0);
	sumamp2=dinn(llt->amp, llt->amp);
	mapfree(lltamp);
    }else{
	real *amps=llt->amp->p;
	const real l2max =pow(lltd*0.5,2);
	/*the waist is defined as the radius where amplitude
	  drop to 1/e or intensity to 1/e^2.*/
	real r2waist=pow(lltd*0.5*parms->powfs[ipowfs].llt->widthp,2);
	//misreg is treated as pupil centering error. It clips the laser
	const real misregx=parms->powfs[ipowfs].llt->misreg->p[0];
	const real misregy=parms->powfs[ipowfs].llt->misreg->p[1];
	const real ox=llt->pts->origx[0];
	const real oy=llt->pts->origy[0];
	for(int iy=0; iy<nx; iy++){
	    real yy=iy*dx+oy;
	    for(int ix=0; ix<nx; ix++){
		real xx=ix*dx+ox;
		real r2=xx*xx+yy*yy;
		real r2m=pow(xx-misregx,2)+pow(yy-misregy,2);
		real amp=exp(-r2m/r2waist);
		sumamp2+=pow(amp, 2);
		if(r2<=l2max){
		    amps[iy*nx+ix]=amp;
		}
	    }
	}
    }
    /*normalized so that max(otf)=1 for unclipped beam; */
    sumamp2=1./(sqrt(sumamp2));
    dscale(llt->amp, sumamp2);
    if(lltcfg->fnsurf){
	mapcell *ncpa=mapcellread("%s",lltcfg->fnsurf);
	int nlotf=ncpa->nx*ncpa->ny;
	assert(nlotf==1 || nlotf==parms->powfs[ipowfs].nwfs);
	llt->ncpa=dcellnew(nlotf, 1);
	for(int ilotf=0; ilotf<nlotf; ilotf++){
	    llt->ncpa->p[ilotf]=dnew(nx,nx);
	    prop_grid_pts(ncpa->p[ilotf], llt->pts, llt->ncpa->p[ilotf]->p, 1, 0, 0, 1, 0, 0, 0);
	}
	cellfree(ncpa);
    }
    /*find negative values in llt->amp and transfer then to surface error.*/
    for(int i=0; i<llt->amp->nx*llt->amp->ny; i++){
	if(llt->amp->p[i]<0){
	    if(!llt->ncpa){
		llt->ncpa=dcellnew(1,1);
		llt->ncpa->p[0]=dnew(nx,nx);
	    }
	    for(int ic=0; ic<llt->ncpa->nx; ic++){
		if(nwvl>1) error("Please implement\n");
		llt->ncpa->p[ic]->p[i]+=wvl0/2;
	    }
	    llt->amp->p[i]=-llt->amp->p[i];
	}
    }
  
    llt->loc=mksqloc(nx,nx,dx,dx, lpts->origx[0], lpts->origy[0]);
    llt->mcc =pts_mcc_ptt(llt->pts, llt->amp->p);
    llt->imcc =dcellinvspd_each(llt->mcc);
    if(lltcfg->focus){
	if(!llt->ncpa){
	    llt->ncpa=dcellnew(1, 1);
	    llt->ncpa->p[0]=dnew(nx,nx);
	}
	real nm2rad=1e-9*2*sqrt(3)*pow(lltcfg->d,-2);
	dmat *focus=dnew(nx, nx);
	loc_add_focus(focus->p, llt->loc, lltcfg->focus*nm2rad);
	real var=0, piston=0;
	long count=0;
	for(long i=0; i<llt->loc->nloc; i++){
	    if(llt->amp->p[i]>0){
		count++;
		var+=focus->p[i]*focus->p[i];
		piston+=focus->p[i];
	    }
	}
	var/=count;
	piston/=count;
	dadds(focus, -piston);
	var=sqrt(var-piston*piston);
	for(int ic=0; ic<llt->ncpa->nx; ic++){
	    dadd(&llt->ncpa->p[ic], 1, focus, lltcfg->focus*1e-9/var);
	}
	cellfree(focus);
	info2("Adding focus %g nm to LLT ncpa\n", lltcfg->focus);
    }
    /*Remove tip/tilt/focus from ncpa */
    if(llt->ncpa && parms->powfs[ipowfs].llt->ttfr){
	int nmod=parms->powfs[ipowfs].llt->ttfr==2?4:3;
	dmat *pttfa=zernike(llt->loc, parms->powfs[ipowfs].llt->d, 0, 2, 0);
	dmat *pttf=drefcols(pttfa, 0, nmod);
	dmat *proj=dpinv(pttf, llt->amp);
	dmat *res=dnew(nmod, 1);
	for(int ilotf=0; ilotf<llt->ncpa->nx*llt->ncpa->ny; ilotf++){
	    dzero(res);
	    dmulvec(res->p, proj, llt->ncpa->p[ilotf]->p, 1);
	    dmulvec(llt->ncpa->p[ilotf]->p, pttf, res->p, -1);
	}
	if(parms->save.setup){
	    writebin(pttf, "powfs%d_llt_pttf", ipowfs);
	    writebin(proj, "powfs%d_llt_proj", ipowfs);
	}
	dfree(pttfa);
	dfree(pttf);
	dfree(proj);
	dfree(res);
    }
    if(parms->save.setup){
	locwrite(llt->loc,"powfs%d_llt_loc",ipowfs);
	writebin(llt->amp,"powfs%d_llt_amp",ipowfs);
	writebin(llt->imcc,"powfs%d_llt_imcc",ipowfs);
	if(llt->ncpa){
	    writebin(llt->ncpa,"powfs%d_llt_ncpa", ipowfs);
	}
    }
    if(parms->save.setup>1){
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    if(powfs[ipowfs].etfprep[iwvl].etf){
		writebin(powfs[ipowfs].etfprep[iwvl].etf, 
			 "powfs%d_etfprep%d",ipowfs,iwvl);
	    }
	}
	if(powfs[ipowfs].etfsim && powfs[ipowfs].etfsim != powfs[ipowfs].etfprep){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(powfs[ipowfs].etfsim[iwvl].etf){
		    writebin(powfs[ipowfs].etfsim[iwvl].etf, 
			       "powfs%d_etfsim%d",ipowfs,iwvl);
		}
	    }
	}
    }
}
/*compute cog NEA using Monte Carlo realizations of noise*/
static void cog_nea(real *nea, dmat *ints, real cogthres, real cogoff, int ntry, 
		    rand_t *rstat, real bkgrnd, real bkgrndc, dmat *bkgrnd2i, dmat *bkgrnd2ic, real rne
		    ){
    dmat *ints2=dnew(ints->nx, ints->ny);
    real gnf[2]={0,0};
    real gny[2]={0,0};
    dcog(gnf, ints, 0, 0, cogthres, cogoff, 0);
    seed_rand(rstat, 1);/*reset the seed each time.*/
    nea[0]=0; nea[1]=0; nea[2]=0; nea[3]=0;
    for(int i=0; i<ntry; i++){
	dcp(&ints2, ints);
	addnoise(ints2, rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, 0, rne, 1);
	dcog(gny, ints2, 0, 0, cogthres, cogoff, 0);
	real errx=gny[0]-gnf[0];
	real erry=gny[1]-gnf[1];
	nea[0]+=errx*errx;
	nea[1]+=errx*erry;
	nea[3]+=erry*erry;
    }
    dfree(ints2);
    real stry=1./ntry;
    nea[0]=nea[0]*stry;
    nea[3]=nea[3]*stry;
    nea[1]=nea[1]*stry;
    nea[2]=nea[1];
}
typedef struct {
    dmat *ints;
    real bkgrnd;
    real bkgrndc;
    dmat *bkgrnd2i;
    dmat *bkgrnd2ic;
    real rne;
    rand_t *rstat;
    int ntry;
}cogdata_t;

/**
   Setup CoG gradient offset for simulation and NEA for reconstruction.
*/
static void 
setup_powfs_cog(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    TIC;tic;
    const int nwfs=parms->powfs[ipowfs].nwfs;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const int ntry=500;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const real pixthetax=parms->powfs[ipowfs].radpixtheta;
    const real pixthetay=parms->powfs[ipowfs].pixtheta;
    const real rne=parms->powfs[ipowfs].rne;
    const real bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
    const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    int do_nea=0;
    rand_t rstat;
    if(!intstat || !intstat->i0){
	if(!parms->powfs[ipowfs].phyusenea){
	    error("powfs[%d].i0 is not available, please enable phyusenea.\n", ipowfs);
	}
    }
    if(parms->powfs[ipowfs].phytype_recon==2 && parms->powfs[ipowfs].skip!=3 && !parms->powfs[ipowfs].phyusenea){
	/*need nea in reconstruction*/
	do_nea=1;
	powfs[ipowfs].sanea=dcellnew(intstat->i0->ny, 1);
	seed_rand(&rstat, 1);
    }
    powfs[ipowfs].cogcoeff=dcellnew(nwfs,1);
    
    for(int jwfs=0; jwfs<nwfs; jwfs++){
	//int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	if(jwfs==0 || (intstat && intstat->i0->ny>1)){
	    real *srot=NULL;
	    if(parms->powfs[ipowfs].radpix){
		srot=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?jwfs:0]->p;
	    }
	    dmat **bkgrnd2=NULL;
	    dmat **bkgrnd2c=NULL;
	    dmat *psanea=0;
	    if(do_nea){
		psanea=powfs[ipowfs].sanea->p[jwfs]=dnew(nsa,3);
		warning("Compute NEA in CoG using Monte Carlo simulation\n");
	    }
	    if(powfs[ipowfs].bkgrnd){
		if(powfs[ipowfs].bkgrnd->ny==1){
		    bkgrnd2=powfs[ipowfs].bkgrnd->p;
		}else{
		    bkgrnd2=powfs[ipowfs].bkgrnd->p+nsa*jwfs;
		}
	    }
	    if(powfs[ipowfs].bkgrndc){
		if(powfs[ipowfs].bkgrndc->ny==1){
		    bkgrnd2c=powfs[ipowfs].bkgrndc->p;
		}else{
		    bkgrnd2c=powfs[ipowfs].bkgrndc->p+nsa*jwfs;
		}
	    }

	    powfs[ipowfs].cogcoeff->p[jwfs]=dnew(2,nsa);
	    dmat*  cogcoeff=powfs[ipowfs].cogcoeff->p[jwfs]/*PDMAT*/;

	    for(int isa=0; isa<nsa; isa++){
		dmat *bkgrnd2i=bkgrnd2?bkgrnd2[isa]:NULL;
		dmat *bkgrnd2ic=bkgrnd2c?bkgrnd2c[isa]:NULL;
	
		P(cogcoeff,0,isa)=parms->powfs[ipowfs].cogthres;
		P(cogcoeff,1,isa)=parms->powfs[ipowfs].cogoff;
		  
		if(do_nea){
		    dmat *nea=dnew(2,2);
		    dmat *ints=intstat->i0->p[isa+jwfs*nsa];/*equivalent noise*/
		    cog_nea(nea->p, ints, P(cogcoeff,0,isa), P(cogcoeff,1,isa), ntry, &rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
		    nea->p[0]=nea->p[0]*pixthetax*pixthetax;
		    nea->p[3]=nea->p[3]*pixthetay*pixthetay;
		    nea->p[1]=nea->p[1]*pixthetax*pixthetay;
		    nea->p[2]=nea->p[1];
		    
		    if(srot){
			dmat *nea2=0;
			drotvecnn(&nea2, nea, srot[isa]);
			dfree(nea); nea=nea2; nea2=0;
		    }
		    P(psanea,isa,0)=nea->p[0];
		    P(psanea,isa,1)=nea->p[3];
		    P(psanea,isa,2)=nea->p[1];
		    dfree(nea);
		}
	    }
	}else{
	    powfs[ipowfs].cogcoeff->p[jwfs]=dref(powfs[ipowfs].cogcoeff->p[0]);
	}
    }//for jwfs
    if(parms->save.setup){
	writebin(powfs[ipowfs].cogcoeff,"powfs%d_cogcoeff", ipowfs);
    }
    toc("setup_powfs_cog");
}


/**
   Setup the (matched filter or CoG) pixel processing parameters for physical optics wfs.
*/
static void 
setup_powfs_phygrad(POWFS_T *powfs,const PARMS_T *parms, int ipowfs){
    long nsa=powfs[ipowfs].saloc->nloc;
    if(powfs[ipowfs].intstat){
	error("Should only be called once\n");
    }
    if(parms->powfs[ipowfs].phytype_recon==1 || parms->powfs[ipowfs].phytype_sim==1 || !parms->powfs[ipowfs].phyusenea
       || (powfs[ipowfs].opdbias && parms->powfs[ipowfs].ncpa_method==2)
       ||parms->powfs[ipowfs].phytype_sim==4
	){
	INTSTAT_T *intstat=powfs[ipowfs].intstat=mycalloc(1,INTSTAT_T);
	if(parms->powfs[ipowfs].i0load){
	    info2("Loading i0, gx, gy\n");
	    if(zfexist("%s/powfs%d_i0",parms->powfs[ipowfs].i0load, ipowfs)){
		intstat->i0=dcellread("%s/powfs%d_i0",parms->powfs[ipowfs].i0load, ipowfs);
	    }else{
		//convert runtime saved ints.
		intstat->i0=dcellnew(nsa, parms->powfs[ipowfs].nwfs);
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    dcell *i0=dcellread("%s/ints_1_wfs%d", parms->powfs[ipowfs].i0load, iwfs);
		    for(int isa=0; isa<nsa; isa++){
			P(intstat->i0, isa, jwfs)=dref(P(i0, isa));
		    }
		    dcellfree(i0);
		}
	    }
	    if((parms->powfs[ipowfs].phytype_recon==1||parms->powfs[ipowfs].phytype_sim==1)
		&& !parms->powfs[ipowfs].mtchfft){
		intstat->gx=dcellread("%s/powfs%d_gx",parms->powfs[ipowfs].i0load, ipowfs);
		intstat->gy=dcellread("%s/powfs%d_gy",parms->powfs[ipowfs].i0load, ipowfs);
	    }
	    if(intstat->i0->header){
		real dt=search_header_num(intstat->i0->header,"dt");
		if(isfinite(dt)){
		    real ratio=parms->sim.dt*parms->powfs[ipowfs].dtrat/dt;
		    warning("Scale loaded i0 by %g\n", ratio);
		    dcellscale(intstat->i0, ratio);
		    dcellscale(intstat->gx, ratio);
		    dcellscale(intstat->gy, ratio);
		}else{
		    warning("Loaded i0 header does not have dt\n");
		}
	    }else{
		warning("Loaded i0 does not have header\n");
	    }
	}else{
	    if(parms->powfs[ipowfs].piinfile){
		/*load psf. 1 for each wavefront sensor. */
		info2("Using 1 sepsf for each wfs when loading sepsf\n");
		intstat->nsepsf=parms->powfs[ipowfs].nwfs;
		intstat->sepsf=dccellnew(parms->powfs[ipowfs].nwfs, 1);
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    dcell *sepsf=intstat->sepsf->p[jwfs]=dcellread("%s_wfs%d",parms->powfs[ipowfs].piinfile,iwfs);
		    real pmax=dmax(sepsf->p[0]);
		    if(sepsf->p[0]->p[0]>pmax*0.5){
			error("wfs %d:  psf must have peak at center, not corner.\n", iwfs);
		    }
		    if(sepsf->nx!=nsa){
			error("piinfile doesn't match\n");
		    }
		}
	    }else if(parms->powfs[ipowfs].lo && parms->powfs[ipowfs].order<=2){
		error("Please specify piinfile for lo order phy wfs\n");
	    }else{
		/*Gen short exposure OTFs due to atmosphere*/
		genseotf(parms, powfs, ipowfs);
		if(parms->powfs[ipowfs].llt){
		    genselotf(parms, powfs, ipowfs);
		}
		/*Generating short exposure psfs for both uplink and downlink turbulence effect. */
		gensepsf(parms,powfs,ipowfs);
		if(parms->save.setup>1 && intstat){
		    writebin(intstat->sepsf->p[0],"powfs%d_sepsf",ipowfs);
		}
		/*Free short exposure otf. */
		ccellfree(intstat->lotf);
		cellfree(intstat->otf);
	    }
	    /*generate short exposure i0,gx,gy from psf. */
	    gensei(parms,powfs,ipowfs);
	}
	if(parms->save.setup){
	    writebin(intstat->i0,"powfs%d_i0",ipowfs);
	    writebin(intstat->gx,"powfs%d_gx",ipowfs);
	    writebin(intstat->gy,"powfs%d_gy",ipowfs);
	}
    }
    /*Generating Matched filter */
    if(parms->powfs[ipowfs].phytype_recon==1 || parms->powfs[ipowfs].phytype_sim==1){
	genmtch(parms,powfs,ipowfs);
	if(parms->save.setup){
	    writebin(powfs[ipowfs].intstat->mtche,"powfs%d_mtche",ipowfs);
	}
    }
    if(parms->powfs[ipowfs].phytype_recon==2 || parms->powfs[ipowfs].phytype_sim==2 || parms->powfs[ipowfs].dither){
	setup_powfs_cog(parms, powfs, ipowfs);
    }
    if(parms->save.setup){
	writebin(powfs[ipowfs].sanea,"powfs%d_sanea",ipowfs);
    }
}
/*
  Setup gradient offset for calibration. opdadd is the wavefront aberration in
  WFS due to optics, without DM correction. opdbias is the wavefront aberration
  in WFS after DM system flat is applied.
*/
void setup_powfs_calib(const PARMS_T *parms, POWFS_T *powfs){
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	//if opdadd is null, but dm_ncpa is not, there will still be opdbias.
	if(powfs[ipowfs].opdbias){
	    //Always compute gradncpa. May not use it in CMF non updated case.
	    if(!powfs[ipowfs].gradncpa){
		powfs[ipowfs].gradncpa=dcellnew(parms->powfs[ipowfs].nwfs,1);
	    }else{
		warning("gradncpa already exists\n");
		dcellzero(powfs[ipowfs].gradncpa);
	    }
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		if(powfs[ipowfs].opdbias->p[jwfs]){
		    real *realamp=powfs[ipowfs].realamp->p[jwfs]->p;
		    if(parms->powfs[ipowfs].type==1){//pywfs
			dmat *ints=0;
			pywfs_fft(&ints, powfs[ipowfs].pywfs, powfs[ipowfs].opdbias->p[jwfs]);
			//writebin(powfs[ipowfs].opdbias->p[jwfs], "opdbias\n");exit(0);
			pywfs_grad(&powfs[ipowfs].gradncpa->p[jwfs], powfs[ipowfs].pywfs, ints);
			dfree(ints);
		    }else if(parms->powfs[ipowfs].gtype_sim==1){//Ztilt
			pts_ztilt(&powfs[ipowfs].gradncpa->p[jwfs], powfs[ipowfs].pts,
				  powfs[ipowfs].saimcc->p[powfs[ipowfs].nsaimcc>1?jwfs:0], 
				  realamp, powfs[ipowfs].opdbias->p[jwfs]->p);
		    }else{//Gtilt
			if(parms->powfs[ipowfs].ncpa_method==1){//GS0*opd
			    dspmm(&powfs[ipowfs].gradncpa->p[jwfs],PR(powfs[ipowfs].GS0, jwfs, 0),
				  powfs[ipowfs].opdbias->p[jwfs],"nn",1);
			}else if(parms->powfs[ipowfs].ncpa_method==2){//CoG(i0)
			    if(!powfs[ipowfs].gradncpa){
				powfs[ipowfs].gradncpa=dcellnew(parms->powfs[ipowfs].nwfs,1);
			    }
			    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
			    shwfs_grad(&powfs[ipowfs].gradncpa->p[jwfs],
					  PCOLR(powfs[ipowfs].intstat->i0, jwfs),
					  parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
			    
			}
		    }
		}
	    }
	}
	if(powfs[ipowfs].pixoffx){
	    if(!powfs[ipowfs].gradncpa){
		powfs[ipowfs].gradncpa=dcellnew(parms->powfs[ipowfs].nwfs,1);
	    }
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int nsa=powfs[ipowfs].pts->nsa;
		if(!powfs[ipowfs].gradncpa->p[jwfs]){
		    powfs[ipowfs].gradncpa->p[jwfs]=dnew(nsa*2, 1);
		}
		
		real pixthetax=parms->powfs[ipowfs].pixtheta;
		real pixthetay=parms->powfs[ipowfs].radpixtheta;
		for(int isa=0; isa<nsa; isa++){
		    real gx=pixthetax*PR(powfs[ipowfs].pixoffx, isa, jwfs);
		    real gy=pixthetay*PR(powfs[ipowfs].pixoffy, isa, jwfs);
		    if(parms->powfs[ipowfs].radpix){
			real angle=PR(PR(powfs[ipowfs].srot, jwfs, 1), isa, 1);
			real ct=cos(angle);
			real st=sin(angle);
			real gx2=gx*ct-gy*st;//RA->XY; CCW
			gy=gx*st+gy*ct;
			gx=gx2;
		    }
		    P(powfs[ipowfs].gradncpa->p[jwfs], isa)+=gx;
		    P(powfs[ipowfs].gradncpa->p[jwfs], isa+nsa)+=gy;
		}
	    }
	}
	if(parms->save.setup){
	    writebin(powfs[ipowfs].gradncpa, "powfs%d_gradncpa", ipowfs);
	}
    }
}

/**
   Setup the powfs struct based on parms and aper. Everything about wfs are
   setup here.  \callgraph */
POWFS_T * setup_powfs_init(const PARMS_T *parms, APER_T *aper){
    TIC;tic;
    POWFS_T *powfs=mycalloc(parms->npowfs,POWFS_T);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].type==0){
	    info("\n%sSetting up powfs %d geom%s\n\n", GREEN,ipowfs,BLACK);
	    setup_shwfs_geom(powfs,parms,aper,ipowfs);
	    setup_shwfs_grad(powfs,parms,ipowfs);
	}else if(parms->powfs[ipowfs].type==1){
	    info("\n%sSetting up powfs %d in Pyramid mode%s\n\n", GREEN,ipowfs,BLACK);
	    pywfs_setup(powfs, parms, aper, ipowfs);
	}else{
	    error("powfs %d: invalid wfstype=%d\n", ipowfs, parms->powfs[ipowfs].type);
	}
    }
    toc("setup_powfs_init");
    return powfs;
}
/**
   Setup physical optics parameters for SHWFS, such as DTF, ETF, LLT, pixel processing.
*/
void setup_powfs_phy(const PARMS_T *parms, POWFS_T *powfs){
    TIC;tic;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs && parms->powfs[ipowfs].type==0 
	   &&(parms->powfs[ipowfs].usephy
	      ||parms->powfs[ipowfs].psfout
	      ||parms->powfs[ipowfs].pistatout
	      ||parms->powfs[ipowfs].neaphy)){
	    info("\n%sSetting up powfs %d PO WFS%s\n\n", GREEN, ipowfs, BLACK);
	    /*We have physical optics. setup necessary struct */
	    setup_powfs_prep_phy(powfs,parms,ipowfs);
	    setup_powfs_dtf(powfs,parms,ipowfs);
	    if(parms->powfs[ipowfs].llt){
		/*prepare Laser launch telescope. */
		setup_powfs_sodium(powfs,parms,ipowfs);/*read sodium profile and smooth it */
		setup_powfs_etf(powfs,parms,0,ipowfs,0,parms->powfs[ipowfs].llt->colprep);/*etf for prep */
		if(!parms->powfs[ipowfs].llt->coldtrat){/*const etf for sim */
		    if(parms->powfs[ipowfs].llt->colprep==parms->powfs[ipowfs].llt->colsim){
			powfs[ipowfs].etfsim=powfs[ipowfs].etfprep;
		    }else{
			setup_powfs_etf(powfs,parms,0,ipowfs,1,parms->powfs[ipowfs].llt->colsim);
		    }
		}
		setup_powfs_llt(powfs,parms,ipowfs);
	    }
	
	    if(parms->powfs[ipowfs].llt){
		/*If there is LLT, setup the extra focus term if needed. */
		setup_powfs_focus(powfs,parms,ipowfs);
	    }
	    if(parms->powfs[ipowfs].usephy || parms->powfs[ipowfs].neaphy){
		setup_powfs_phygrad(powfs,parms,ipowfs);
	    }
	}
    }/*ipowfs */
    toc("setup_powfs_phy");
}
/**
   free unused parameters before simulation starts
*/
void free_powfs_unused(const PARMS_T *parms, POWFS_T *powfs){
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(powfs[ipowfs].intstat){
	    cellfree(powfs[ipowfs].intstat->sepsf);
	    cellfree(powfs[ipowfs].intstat->gx);
	    cellfree(powfs[ipowfs].intstat->gy);
	}
	if(!parms->powfs[ipowfs].needGS0 && powfs[ipowfs].GS0){
	    dspcellfree(powfs[ipowfs].GS0);
	    powfs[ipowfs].GS0=NULL;
	}
    }
}

static void free_powfs_shwfs(POWFS_T *powfs, int ipowfs){
    dtf_free(powfs[ipowfs].dtf);
    dspcellfree(powfs[ipowfs].GS0);
    dcellfree(powfs[ipowfs].neasim);
    dcellfree(powfs[ipowfs].sanea);
    cellfree(powfs[ipowfs].cogcoeff);
    if(powfs[ipowfs].intstat){
	INTSTAT_T *intstat=powfs[ipowfs].intstat;
	cellfree(intstat->fotf);
	cellfree(intstat->potf);
	cellfree(intstat->mtche);
	cellfree(powfs[ipowfs].intstat->i0);
	dfree(intstat->i0sum);
	dfree(intstat->i0sumsum);
	free(intstat);
	powfs[ipowfs].intstat=NULL;
    }
    {
	dcellfree(powfs[ipowfs].srot);
	dcellfree(powfs[ipowfs].srsa);
	dfree(powfs[ipowfs].srsamax);
	dcellfree(powfs[ipowfs].sprint);
	dfree(powfs[ipowfs].pixoffx);
	dfree(powfs[ipowfs].pixoffy);
    }

    cellfree(powfs[ipowfs].saimcc);
    if(powfs[ipowfs].llt){
	ptsfree(powfs[ipowfs].llt->pts);
	dfree(powfs[ipowfs].llt->amp);
	locfree(powfs[ipowfs].llt->loc);
	dcellfree(powfs[ipowfs].llt->mcc);
	dcellfree(powfs[ipowfs].llt->imcc);
	dcellfree(powfs[ipowfs].llt->ncpa);
	free(powfs[ipowfs].llt);
    }
    dcellfree(powfs[ipowfs].sodium);
    if(powfs[ipowfs].etfprep!=powfs[ipowfs].etfsim){
	etf_free(powfs[ipowfs].etfprep);
    }
    etf_free(powfs[ipowfs].etfsim);
    etf_free(powfs[ipowfs].etfsim2);
    dcellfree(powfs[ipowfs].opdadd);
    dcellfree(powfs[ipowfs].opdbias);
    dcellfree(powfs[ipowfs].gradncpa);
    dfree(powfs[ipowfs].dtheta);
}
/**
   Free all parameters of powfs at the end of simulation.
*/
void free_powfs(const PARMS_T *parms, POWFS_T *powfs){
    free_powfs_unused(parms, powfs);
    free_powfs_fit(powfs,parms);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	free_powfs_geom(powfs,ipowfs);
	free_powfs_shwfs(powfs, ipowfs);
	pywfs_free(powfs[ipowfs].pywfs);
    }
    free(powfs);
}
