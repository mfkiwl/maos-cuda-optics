/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

/**
   \file setup_surf.c
   Setup additional NCPA surface.
*/
#include "maos.h"
#include "setup_surf.h"
#include "recon_utils.h"
/**
   Setup tilted surface (M3) by ray tracing from the tilted surface to WFS and
   Science grid.

   \todo Merge setup_tsurf, setup_surf, and setup_powfs_ncpa. Be careful about
   which NCPA is compensated by WFS offsets. setup_tsurf and setup_surf act on
   surfaces intersect with most WFS and science, while setup_powfs_ncpa act on
   individual surfaces for each WFS. Maybe it is good idea to keep them
   separate.


   2011-09-07: 
   Relocate simu->surfwfs to powfs.opdadd In order to perserve surface for different seeds
*/
/**
   Propagate tilt surface from star at hs, along direction (thetax, thetay), to loc.
*/
static void tsurf2loc(rectmap_t **tsurf, int ntsurf, dmat *opd, loc_t *locin, double thetax, double thetay, double hs, double rot){
    loc_t *locuse=NULL;
    if(fabs(rot)>1e-10){
	locuse=locdup(locin);
	locrot(locuse, rot);
    }else{
	locuse=locin;
    }
    for(int itsurf=0; itsurf<ntsurf; itsurf++){
	const double alx=tsurf[itsurf]->txdeg/180*M_PI;
	const double aly=tsurf[itsurf]->tydeg/180*M_PI;
	const double ftel=tsurf[itsurf]->ftel;
	const double fexit=tsurf[itsurf]->fexit;
	const double fsurf=tsurf[itsurf]->fsurf;
	const double mag=fexit/ftel;
	const double scalex=-mag;
	const double scaley=mag;
	const double scaleopd=-2;
	const double het=fexit-fsurf;/*distance between exit pupil and M3. */
	rectmap_t *mapsurf=tsurf[itsurf];

	double d_img_focus=1./(1./ftel-1./hs)-ftel;
	/*info2("iwfs%d: d_img_focus=%g\n",iwfs,d_img_focus); */
	double d_img_exit=fexit+d_img_focus;
		
	/*2010-04-02: do not put - sign */
	double bx=thetax*(d_img_focus+ftel)/d_img_exit;
	double by=thetay*(d_img_focus+ftel)/d_img_exit;
	proj_rect_grid(mapsurf,alx,aly,locuse,scalex,scaley, NULL,opd->p,scaleopd, d_img_exit, het, bx, by);
    }
    if(locuse!=locin){
	locfree(locuse);
    }
}

static void 
setup_surf_tilt(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs){
    info("Setting up tilt surface (M3)\n");
    rectmap_t **tsurf=calloc(parms->ntsurf, sizeof(rectmap_t*));
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	char *fn=parms->tsurf[itsurf];
	info("Loading tilt surface from %s\n", fn);
	tsurf[itsurf]=rectmapread("%s",fn); 
    }
    const double rot=-parms->aper.rotdeg/180.*M_PI;

    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	tsurf2loc(tsurf, parms->ntsurf, aper->opdadd->p[ievl], aper->locs, 
		  parms->evl.thetax[ievl], parms->evl.thetay[ievl], parms->evl.hs[ievl], rot);
    }

    for(int iwfs=0; iwfs<parms->nwfs && powfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	loc_t *locwfsin;
	    
	if(powfs[ipowfs].nlocm){
	    error("We don't handle this case yet. Think carefully when to apply shift.\n");
	    int ilocm=powfs[ipowfs].nlocm>1?parms->powfs[ipowfs].wfsind[iwfs]:0;
	    locwfsin=powfs[ipowfs].locm[ilocm];
	}else{
	    locwfsin=powfs[ipowfs].loc;
	}
	tsurf2loc(tsurf, parms->ntsurf, powfs[ipowfs].opdadd->p[wfsind], locwfsin, 
		  parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, parms->powfs[ipowfs].hs, rot);
    }
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	rectmapfree(tsurf[itsurf]);
    }
    free(tsurf);
}

/**
   Setup surface perpendicular to the beam by ray tracing from the surface to
   WFS and Science grid
*/
static void 
setup_surf_perp(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs, RECON_T *recon){
    info2("Setting up surface OPD (M1/M2/M3)\n");
    if(fabs(parms->aper.misreg[0])>EPS || fabs(parms->aper.misreg[1])>EPS){
	warning("Please adjust telescope surface ox, oy to account for misregistration. Not doing "
		"in maos because some surfaces may belong to instrument.\n");
    }
    loc_t *locevl;
    const double rot=-parms->aper.rotdeg/180.*M_PI;
    const int nevl=parms->evl.nevl;
    const int nwfs=parms->nwfs;
    int *evlcover=malloc(nevl*sizeof(int));
    int *wfscover=malloc(nwfs*sizeof(int));
    int opdxcover;
    for(int isurf=0; isurf<parms->nsurf; isurf++){
	char *fn=parms->surf[isurf];
	if(!fn) continue;
	info("Loading surface OPD from %s\n", fn);
	map_t *surf=mapread("%s",fn);
	dwrite((dmat*)surf, "surf_%d", isurf);
	const char *strname=search_header(surf->header, "SURFNAME");
	const char *strevl=search_header(surf->header, "SURFEVL");
	const char *strwfs=search_header(surf->header, "SURFWFS");
	const char *stropdx=search_header(surf->header, "SURFOPDX");
	int do_rot=0;
	if(strname && !strcmp(strname, "M1")){
	    warning("Rotate loc for M1\n");
	    do_rot=(fabs(rot)>1.e-10);
	}
	if(do_rot){
	    locevl=locdup(aper->locs);
	    locrot(locevl,rot);
	}else{
	    locevl=aper->locs;
	}
	if(!strevl){
	    warning2("surf[%d] does not contain SURFEVL\n", isurf);
	    for(int ievl=0; ievl<nevl; ievl++){
		evlcover[ievl]=1;
	    }
	}else{
	    readstr_intarr_nmax(&evlcover, nevl, strevl);
	}
	if(!strwfs){
	    warning2("surf[%d] does not contain SURFWFS\n", isurf);
	    for(int iwfs=0;iwfs<nwfs; iwfs++){
		wfscover[iwfs]=1;
	    }
	}else{
	    readstr_intarr_nmax(&wfscover, nwfs, strwfs);
	}
	if(!stropdx){
	    opdxcover=1;
	}else{
	    opdxcover=(int)readstr_num(stropdx, NULL);
	}
	double hl=surf->h;
	for(int ievl=0; ievl<nevl; ievl++){
	    if(!evlcover[ievl]){
		warning2("Skip evl %d for surface %s\n", ievl, parms->surf[isurf]);
		continue;
	    }
	    const double displacex=parms->evl.thetax[ievl]*hl;
	    const double displacey=parms->evl.thetay[ievl]*hl;
	    
	    if(do_rot){
		prop_grid(surf, locevl, NULL, aper->opdadd->p[ievl]->p, 
			  1, displacex, displacey, 1, 0, 0, 0);
	    }else{
		prop_grid_stat(surf, aper->locs->stat, 
			       aper->opdadd->p[ievl]->p, 
			       1, displacex, displacey, 1, 0, 0, 0);
	    }
	}
	for(int iwfs=0; iwfs<parms->nwfs && powfs; iwfs++){
	    if(!wfscover[iwfs]){
		warning2("Skip wfs %d for surface %s\n", iwfs, parms->surf[isurf]);
		continue;
	    }
	    const int ipowfs=parms->wfs[iwfs].powfs;
	    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	    const double hs=parms->powfs[ipowfs].hs;
	    const double scale=1.-hl/hs;
	    const double displacex=parms->wfs[iwfs].thetax*hl+powfs[ipowfs].misreg[wfsind][0];
	    const double displacey=parms->wfs[iwfs].thetay*hl+powfs[ipowfs].misreg[wfsind][1];

	    loc_t *locwfs, *locwfsin;
	    if(powfs[ipowfs].locm){
		int ilocm=powfs[ipowfs].nlocm>1?parms->powfs[ipowfs].wfsind[iwfs]:0;
		locwfsin=powfs[ipowfs].locm[ilocm];
	    }else{
		locwfsin=powfs[ipowfs].loc;
	    }
	    if(do_rot){
		locwfs=locdup(locwfsin);
		locrot(locwfs,rot);
	    }else{
		locwfs=locwfsin;
	    }
	    prop_grid(surf, locwfs, NULL, powfs[ipowfs].opdadd->p[wfsind]->p, 
		      1, displacex, displacey, scale, 1., 0, 0); 
	    if(do_rot){
		locfree(locwfs);
	    }
	}
	if(parms->sim.idealfit && recon && opdxcover){
	    if(!recon->opdxadd){
		recon->opdxadd=dcellnew(parms->atmr.nps, 1);
	    }
	    double distmin=INFINITY;
	    int jpsr=-1;
	    /*Select the layer that is closed to the surface. */
	    for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
		double dist=fabs(parms->atmr.ht[ipsr]-hl);
		if(dist < distmin){
		    jpsr=ipsr;
		    distmin=dist;
		}
	    }
	    loc_t *xloc;
	    if(do_rot){
		xloc = locdup(recon->xloc[jpsr]);
		locrot(xloc, rot);
	    }else{
		xloc = recon->xloc[jpsr];
	    }
	    loc_t *surfloc=mksqloc_map(surf);
	    dsp *H=mkhb(xloc, surfloc, NULL, 0, 0, 1, 0, 0);
	    double scale=pow(surf->dx/xloc->dx,2);
	    if(!recon->opdxadd->p[jpsr]){
		recon->opdxadd->p[jpsr]=dnew(xloc->nloc, 1);
	    }
	    spmulvec(recon->opdxadd->p[jpsr]->p, H, surf->p, scale);
	    if(do_rot) locfree(xloc);
	    locfree(surfloc);
	}
	mapfree(surf);
    }
    if(locevl!=aper->locs){
	locfree(locevl);
    }
    free(evlcover);
    free(wfscover);
}

/** We trace rays from Science focal plan OPD to ploc along evaluation
    directions (type=1) or on axis only (type=2).*/
static void FitR_NCPA(dcell **xout, RECON_T *recon, APER_T *aper){
    const PARMS_T *parms=recon->parms;
    dcell *xp=dcellnew(parms->evl.nevl, 1);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	xp->p[ievl]=dnew(recon->floc->nloc,1);
	prop_nongrid(aper->locs, aper->opdadd->p[ievl]->p,
		     recon->floc, NULL, xp->p[ievl]->p, 1, 0, 0, 1, 0, 0);
    }
    applyW(xp, recon->W0, recon->W1, parms->evl.wt);
    sptcellmulmat_thread(xout, recon->HAevl, xp, 1);
    dcellfree(xp);
}
void FitL_NCPA(dcell **xout, const void *A, 
	       const dcell *xin, const double alpha){
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    dcell *xp=NULL;
    spcellmulmat_thread(&xp, recon->HAevl, xin, 1.);
    applyW(xp, recon->W0, recon->W1, parms->evl.wt);
    sptcellmulmat_thread(xout, recon->HAevl, xp, alpha);
    dcellfree(xp);xp=NULL;
    dcellmm(&xp,recon->fitNW, xin, "tn", 1);
    dcellmm(xout,recon->fitNW, xp, "nn", alpha);
    dcellfree(xp);
    if(recon->actslave){
	spcellmulmat(xout, recon->actslave, xin, 1);
    }
}
static void setup_recon_HAevl(RECON_T *recon, const PARMS_T *parms){
    const int nevl=parms->evl.nevl;
    const int ndm=parms->ndm;
    recon->HAevl=spcellnew(nevl, ndm);
    PDSPCELL(recon->HAevl,HA);
    info2("Generating HA ");TIC;tic;
    for(int ievl=0; ievl<nevl; ievl++){
	double hs=parms->evl.hs[ievl];
	for(int idm=0; idm<ndm; idm++){
	    const double ht=parms->dm[idm].ht;
	    const double scale=1.-ht/hs;
	    double displace[2];
	    displace[0]=parms->evl.thetax[ievl]*ht;
	    displace[1]=parms->evl.thetay[ievl]*ht;
	    HA[idm][ievl]=mkh(recon->aloc[idm], recon->floc, NULL,
			      displace[0], displace[1], 
			      scale,parms->dm[idm].cubic,parms->dm[idm].iac);
	}
    }
    toc2(" ");
    if(parms->save.setup){
	spcellwrite(recon->HA,"%s/HA",dirsetup);
    }
}
#include "mtch.h"
#include "genseotf.h"
void setup_surf(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs, RECON_T *recon){
    if(parms->nsurf<=0 && parms->ntsurf<=0){
	info2("No surfaces to setup\n");
	return;
    }
    if(aper->opdadd){
	error("Who sets aper->opdadd?\n");
    }
    aper->opdadd=dcellnew(parms->evl.nevl,1);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	aper->opdadd->p[ievl]=dnew(aper->locs->nloc, 1);
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	powfs[ipowfs].opdadd=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    powfs[ipowfs].opdadd->p[jwfs]=dnew(powfs[ipowfs].npts, 1);
	}
    }

    if(parms->ntsurf>0){
	setup_surf_tilt(parms, aper, powfs);
    }
    if(parms->nsurf>0){
	setup_surf_perp(parms, aper, powfs, recon);
    }
    if(parms->sim.ncpa_calib){//calibrate NCPA
	info("calibrating NCPA\n");
	setup_recon_HAevl(recon, parms);
	dcell *rhs=NULL;
	dcell *dmn=NULL;
	FitR_NCPA(&rhs, recon, aper);
	int maxit=40;
	pcg(&dmn, FitL_NCPA, recon, NULL, NULL, rhs, 1, maxit);

	dcellwrite(dmn, "dm_ncpa");
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    dcellcp(&powfs[ipowfs].opdbias, powfs[ipowfs].opdadd);
	    for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
		int iwfs0=parms->powfs[ipowfs].wfs[iwfs];
		double hs=parms->powfs[ipowfs].hs;
		double thetax=parms->wfs[iwfs0].thetax;
		double thetay=parms->wfs[iwfs0].thetay;
		for(int idm=0; idm<parms->ndm; idm++){
		    double ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
		    double scale=1-ht/hs;
		    double dispx=ht*thetax;
		    double dispy=ht*thetay;
		    if(parms->dm[idm].cubic){
			prop_nongrid_pts_cubic(recon->aloc[idm], dmn->p[idm]->p, 
					       powfs[ipowfs].pts, NULL, powfs[ipowfs].opdbias->p[iwfs]->p, 
					       -1, dispx, dispy, 1, parms->dm[idm].iac, 0, 0);
		    }else{
			prop_nongrid_pts(recon->aloc[idm], dmn->p[idm]->p, 
					 powfs[ipowfs].pts, NULL, powfs[ipowfs].opdbias->p[iwfs]->p, 
					 -1, dispx, dispy, 1, 0, 0);
		    }
		}
	    }

	    if(parms->powfs[ipowfs].ncpa_method==1){
		if(!powfs[ipowfs].gradoff){
		    powfs[ipowfs].gradoff=dcellnew(parms->powfs[ipowfs].nwfs,1);
		}
		for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
		    double *realamp=powfs[ipowfs].realamp[iwfs];
		    if(parms->powfs[ipowfs].gtype_sim==1){
			pts_ztilt(&powfs[ipowfs].gradoff->p[iwfs], powfs[ipowfs].pts,
				  powfs[ipowfs].nsaimcc>1?powfs[ipowfs].saimcc[iwfs]:powfs[ipowfs].saimcc[0], 
				  realamp, powfs[ipowfs].opdbias->p[iwfs]->p);
		    }else{
			spmulmat(&powfs[ipowfs].gradoff->p[iwfs],adpind(powfs[ipowfs].GS0, iwfs),
				 powfs[ipowfs].opdbias->p[iwfs],1);
		    }
		}
	    
	    }else if(parms->powfs[ipowfs].ncpa_method==2){
		genseotf(parms,powfs,ipowfs);
		gensepsf(parms,powfs,ipowfs);
		gensei(parms,powfs,ipowfs);
		genmtch(parms,powfs,ipowfs);
		if(parms->save.setup){
		    dcellwrite(powfs[ipowfs].intstat->i0,"%s/powfs%d_i0_2",dirsetup,ipowfs);
		    dcellwrite(powfs[ipowfs].intstat->gx,"%s/powfs%d_gx_2",dirsetup,ipowfs);
		    dcellwrite(powfs[ipowfs].intstat->gy,"%s/powfs%d_gy_2",dirsetup,ipowfs);
		}
	    }
	}
	/* to do 
	   dm flat
	   matched filter
	   genseotf() reentrant. ok
	   genselotf() reentrant. ok
	   gensepsf() reentrant. ok
	   gensei() reentrant. ok.
	   genmtch() reentrant. ok
	*/
    }
    if(parms->save.setup){
	dcellwrite(aper->opdadd, "%s/surfevl.bin", dirsetup);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    dcellwrite(powfs[ipowfs].opdadd, "%s/surfpowfs_%d.bin", dirsetup, ipowfs);
	}
	if(recon->opdxadd) dcellwrite(recon->opdxadd, "%s/surfopdx", dirsetup);
    }
}
