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
#include "mkdtf.h"
#include "mkh.h"
/**
   mkdtf() computes the parameters used to sample PSFs onto detectors in the
   Fourier domain. It incorporates size of the detector pixels, and charge
   diffusion. For polar coordinate detectors, it also rotates the pixels and
   pixel coordinates so that they are aranged along the radial and azimuthal
   direction.
 */
DTF_T *mkdtf(const dmat *wvls, /**<List of wavelength*/
	     real dxsa,        /**<Subaperture size*/
	     real embfac,      /**<Embedding factor (2)*/
	     long notfx,       /**<FFT size along x*/
	     long notfy,       /**<FFT size along y*/
	     long pixpsax,     /**<Number of pixels along x(r)*/
	     long pixpsay,     /**<Number of pixels along y(a)*/
	     real pixthetax,   /**<Pixel size along x (r)*/
	     real pixthetay,   /**<Pixel size along y (a)*/
	     const dmat* pixoffx,  /**<offset of image center from center of pixel array, along x or radial*/
	     const dmat* pixoffy,  /**<offset of image center from center of pixel array, along y or azimuthal*/
	     real pixblur,     /**<Pixel blur sigma(fraction of pixel)*/
	     const dcell *srot,/**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     int radpix        /**<1: Pixels are along radial/azimuthal direction*/
    ){
    int nwvl=wvls->nx*wvls->ny;
    DTF_T *dtfs=mycalloc(nwvl,DTF_T);
    dtfs->nwvl=nwvl;
    const real blurx=pixblur*pixthetax;
    const real blury=pixblur*pixthetay;
    const real e0x=-2*M_PI*M_PI*blurx*blurx;//blurring factors
    const real e0y=-2*M_PI*M_PI*blury*blury;
    const int do_blur=fabs(blurx)>EPS && fabs(blury)>EPS;
    const long notfx2=notfx>>1;
    const long notfy2=notfy>>1;
    const real pxo=-(pixpsax*0.5-0.5)*pixthetax;
    const real pyo=-(pixpsay*0.5-0.5)*pixthetay;
    real pxo2=pxo;
    real pyo2=pyo;
    if(pixoffx || pixoffy){
	int nwfs=1; int nsa=1;
	if(srot){
	    nwfs=srot->nx;
	    nsa=srot->p[0]->nx;
	}
	if(!pixoffx || !pixoffy){
	    error("pixoffx and pixoffy must simultaneously be supplied\n");
	}
	if(pixoffx->nx!=1 && pixoffx->nx!=nsa && pixoffx->ny!=1 && pixoffx->ny!=nwfs){
	    error("pixoffx has wrong format\n");
	}
	if(pixoffy->nx!=pixoffx->nx || pixoffy->ny!=pixoffx->ny){
	    error("pixoffy must have the same format as pixoffx\n");
	}
    }

    /*both nominal and si depends on wavelength.*/
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const real wvl=wvls->p[iwvl];
	const real dtheta=wvl/(dxsa*embfac);/*PSF sampling. */
	const real dux=1./(dtheta*notfx);
	const real duy=1./(dtheta*notfy);
	const real dux2=dux*dux;
	const real duy2=duy*duy;
	const real pdtheta=pixthetax*pixthetay/(dtheta*dtheta);//scaling factor due to binning into detectors
	const real duxp=dux*pixthetax;
	const real duyp=duy*pixthetay;
	/*For LGS WFS (srot!=NULL), there is elongation effect. For radial pixel
	  CCD without using rotating psf/otf method, we need to create
	  nominal/si for each subaperture. The PSF/OTF and hence DTF are on x/y
	  coordinate, so pixels and pixel coordinates need to rotated to r/a
	  direction. */
	const int multi_dtf=(srot && radpix);//DTF is along x/y coord, pixel is r/a coord.
	const int ndtf=(multi_dtf || (pixoffx && pixoffx->nx>1))?srot->p[0]->nx:1;
	const int multi_wfs=((srot && srot->nx>1 && radpix) || (pixoffx && pixoffx->ny>1));
	const int nwfs=multi_wfs?MAX(srot->nx, pixoffx->ny):1;

	dtfs[iwvl].dtheta=dtheta;
	dtfs[iwvl].dxsa=dxsa;
	dtfs[iwvl].radpix=radpix;
	dtfs[iwvl].wvl=wvl;
	dtfs[iwvl].notfx=notfx;
	dtfs[iwvl].notfy=notfy;
	dtfs[iwvl].nominal=ccellnew_same(ndtf,nwfs,notfx,notfy);
	dtfs[iwvl].si=dspcellnew(ndtf,nwfs);
	ccell*  nominals=dtfs[iwvl].nominal;
	dspcell*  sis=dtfs[iwvl].si;
	cmat *nominal=cnew(notfx,notfy);
	//Coordinate of PSF points
	loc_t *loc_psf=mksqloc(notfx,notfy,dtheta,dtheta,-notfx2*dtheta, -notfy2*dtheta);
	real theta=0;
	real ct=cos(theta);
	real st=sin(theta);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    for(int isa=0; isa<ndtf; isa++){
		if(multi_dtf){
		    theta=PR(srot, iwfs, 1)->p[isa];
		    ct=cos(theta);
		    st=sin(theta);
		}
		
		for(int iy=0; iy<notfy; iy++){
		    int jy=iy-notfy2;
		    for(int ix=0; ix<notfx; ix++){
			int jx=ix-notfx2;
			real ir=ct*jx+st*jy;
			real ia=-st*jx+ct*jy;
			//Pixel function
			P(nominal,ix,iy)=sinc(ir*duxp)*sinc(ia*duyp)*pdtheta;
			if(do_blur){//Charge diffusion.
			    P(nominal,ix,iy)*=exp(e0x*(ir*ir*dux2)+e0y*(ia*ia*duy2));
			}
		    }
		}
		/*put peak in corner. */
		cfftshift(nominal);
		/*Embed a fftshift operation in nominal so that we avoid fftshift when forming i0.*/
		cfft2(nominal,-1);
		cfftshift(nominal);
		cfft2(nominal,1);
		//cancel FFT scaling effect.
		cscale(nominal,1./(real)(nominal->nx*nominal->ny));
		ccp(PP(nominals,isa,iwfs), nominal);
		//Coordinate of pixels
		if(pixoffx){
		    real dx=PR(pixoffx, isa, iwfs)*pixthetax;
		    real dy=PR(pixoffy, isa, iwfs)*pixthetay;
		    pxo2=pxo-dx;
		    pyo2=pyo-dy;
		}
		loc_t *loc_ccd=mksqloc(pixpsax,pixpsay, pixthetax,pixthetay,pxo2,pyo2);
		locrot(loc_ccd, theta);
		P(sis,isa,iwfs)=mkh(loc_psf,loc_ccd,0,0,1);
		locfree(loc_ccd);
	    }/*isa */
	}/*iwfs */
	cfree(nominal);
	locfree(loc_psf);

	/*Create an excessive high frequency in nominal so that
	  we don't have to do fftshift later.*/
	dtfs[iwvl].Ux=cnew(notfx,1);
	dtfs[iwvl].Uy=cnew(notfy,1);
	comp *Ux=dtfs[iwvl].Ux->p;
	comp *Uy=dtfs[iwvl].Uy->p;

	/*The following is used in genseotf to compute shifted i0.*/
	for(int ix=0; ix<notfx; ix++){
	    int jx=ix<notfx2?ix:(ix-notfx);
	    Ux[ix]=COMPLEX(0, -2.*M_PI*jx*dux);
	}
	for(int iy=0; iy<notfy; iy++){
	    int jy=iy<notfy2?iy:(iy-notfy);
	    Uy[iy]=COMPLEX(0, -2.*M_PI*jy*duy);
	}

    }/*iwvl */
    return dtfs;
}

/**
   Wrap the index for dataset with total of n frames for continuity. The actual data becomes
   0, 1, 2, ..., n-2, n-1, n-2, ..., 0, 1
*/
static inline int wrap_seq(long index, long n){
    long m=n*2-1;
    index=index%m;
    if(index<0) index+=m;
    if(index>=n) index=m-1-index;
    return index;
}

ETF_T *mketf(DTF_T *dtfs,  /**<The dtfs*/
	     real hs,      /**<LGS WFS focus altitude*/
	     const dcell *sodium,/**<The sodium profile. In each cell First column is coordinate.*/
	     int icol,     /**<Which sodium profile to use*/
	     const dcell *srot,  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     const dcell *srsa,  /**<Subaperture to LLT distance*/
	     int no_interp /**<Use direct sum instead of interpolation + FFT. Slower */
    ){
    int nwvl=dtfs[0].nwvl;
    ETF_T *etfs=mycalloc(nwvl,ETF_T);
    etfs->nwvl=nwvl;
    /*setup elongation along radial direction. don't care azimuthal. */
    if(!srot) error("srot is required");
    const int nllt=MAX(srot->nx, sodium->nx);
    const int nsa=srot->p[0]->nx;
    dmat *sodium0=sodium->p[0];
    const int ncol=sodium0->ny-1;
    const int nhp=sodium0->nx; 
    const real* px=sodium0->p;

    icol=wrap_seq(icol, ncol);
    //info("Na using column %d.\n",icol);

    //adjusting sodium height for the zenith angle;
    real hpmin=0, dhp1=0;

    if(!no_interp){/**Requires linear spacing*/
	hpmin=px[0];
	dhp1=1./(px[1]-px[0]);
	/*assume linear spacing. check the assumption valid */
	real diff;
	if((diff=fabs((px[nhp-1]-px[0])/((nhp-1)*(px[1]-px[0]))-1.))>1.e-5){
	    error("sodium profile is not evenly spaced: %g\n", diff);
	}
    }//else: Any spacing is ok.

    /*the sum of pp determines the scaling of the pixel intensity. */
    real *psrot[nllt];
    real *pna[nllt];
    real i0scale[nllt];
    for(int illt=0; illt<nllt; illt++){
	psrot[illt]=srot->p[srot->nx>1?illt:0]->p;
	pna[illt]=PR(sodium, illt, 0)->p+nhp*(1+icol);
	i0scale[illt]=dvecsum(pna[illt], nhp);
	if(fabs(i0scale[illt]-1)>0.01){
	    warning("Siglev is scaled by %g by sodium profile\n", i0scale[illt]);
	}
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const real dtheta=dtfs[iwvl].dtheta;
	const long notfx=dtfs[iwvl].notfx;
	const long notfy=dtfs[iwvl].notfy;
	const long notfx2=notfx>>1;
	const long notfy2=notfy>>1;
	const real dux=1./(dtheta*notfx);
	const real duy=1./(dtheta*notfy);
	ccell *petf=0;

	/*
	  PSF/OTF/ETF is defined in x/y coordinate. 
	  2010-01-04: Fuse dtf nominal into etf for this case.
	*/
	petf=etfs[iwvl].etf=ccellnew(nsa,nllt);
	
	etfs[iwvl].icol=icol;
	if(no_interp){ /* No interpolation, no fft. intensity scaling is automatically taken care of */
	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    real rsa=srsa->p[illt]->p[isa];
		   
		    const real theta=psrot[illt][isa];
		    const real ct=cos(theta);
		    const real st=sin(theta);
		    P(petf,isa,illt)=cnew(notfx,notfy);
		    cmat *etf2d=P(petf,isa,illt);
#pragma omp parallel for default(shared)
		    for(int icompy=0; icompy<notfy; icompy++){
			const real ky=duy*(icompy>=notfy2?(icompy-notfy):icompy);
			for(int icompx=0; icompx<notfx; icompx++){
			    const real kx=dux*(icompx>=notfx2?(icompx-notfx):icompx);
			    const real kr=(ct*kx+st*ky);/*along radial*/
			    for(int ih=0; ih<nhp; ih++){
				const real tmp=(-2*M_PI*(kr*(rsa/px[ih]-rsa/hs)));
				P(etf2d,icompx,icompy)+=COMPLEX(pna[illt][ih]*cos(tmp), pna[illt][ih]*sin(tmp));
			    }
			}
		    }
		    
		}//isa
	    }//illt
	}else{
	    /*
	      The ETF is computed as DFT
	      ETF(k_i)=\sum_j exp(-2*\pi*I*k_i*{\theta}_j)P({\theta}_j)
	      Where \theta_j is the radial coord of pixel j, and
	      P({\theta}_j}=sodium(h_j) with h_j=rsa*hs/(rsa+hs*{\theta}_j)
	      where hs is GS distance and rsa is subaperture to LLT distance
	      We replace the above interpolation (and FFT) by doing summation directly
	      ETF(k_i)=\sum_j exp(-2*\pi*I*k_i*(rsa/h_j-rsa/hs))P(h_j)
	    */
	    const int npad=2;/*zero padding to reduce aliasing */
	    const int nover=2;/*enough size for rotation*/
	    const int netf=notfx*nover*npad;
	    const real dtetf=dtheta/nover;
	    const real dusc=(netf*dtetf)/(dtheta*notfx);
	    cmat *etf=cnew(netf,1);
	    real *thetas=mycalloc(netf,real);
	    const int netf2=netf>>1;
	    /*Only interpolating the center part. the rest is padding. */
	    const int etf0=netf2-(int)round(notfx2*(dtheta/dtetf));
	    const int etf1=etf0+(int)round(notfx*(dtheta/dtetf));
	    if(etf0<0) error("Invalid configuration\n");
	    for(int it=etf0; it<etf1; it++){
		thetas[it]=(it-netf2)*dtetf;
	    }
	    //info("etf0=%d, etf1=%d, netf=%d\n", etf0, etf1, netf);
	    //cfft2plan(etf, -1);

	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    /*1d ETF along radius. */
		    real rsa=srsa->p[illt]->p[isa];
		    real etf2sum=0;
		    czero(etf);
		    for(int icomp=etf0; icomp<etf1; icomp++){
			/*peak in center */
			const real itheta=thetas[icomp];
			/*non linear mapping. changed from - to + on 2014-05-07.*/
			const real ih=hs*rsa/(rsa+itheta*hs);
			/*interpolating to get Na profile strenght. */
			/*this is bilinear interpolation. not good. need to do averaging */
			const real iih=(ih-hpmin)*dhp1;
			const int iihf=ifloor(iih);
			const real iihw=iih-iihf;
			if(iihf<0 || iihf>nhp-2){
			    etf->p[icomp]=0.;
			}else{
			    real tmp=pna[illt][iihf]*(1.-iihw)+pna[illt][iihf+1]*iihw;
			    /*neglected rsa1 due to renormalization. */
			    etf->p[icomp]=tmp;
			    etf2sum+=tmp;
			}
		    }
		    if(fabs(etf2sum)>1.e-20){
			/*2010-11-09:

			  We used to normalize the etf before fft so that after fft
			  it max to 1. The strength of original profile doesn't
			  matter.
		    
			  Changed: We normalize the etf by sum(profile), so we can model
			  the variation of the intensity and meteor trails.
		      
			*/
			cscale(etf,i0scale[illt]/etf2sum);

			//check for truncation
			double ratio_edge=0.5*creal(etf->p[etf0]+etf->p[etf1-1])/creal(etf->p[netf2]);
			if(ratio_edge>0.1){
			    writebin(etf, "etf_%d", isa);
			    error("sa %d: sodium profile is cropped when computing etf. Increase powfs.pixpsa or powfs.notf.\n", isa);
			}
			cfftshift(etf);/*put peak in corner; */
			cfft2(etf, -1);
		
			/*Rotate the ETF. */
			/*need to put peak in center. */
			cfftshift(etf);
			real theta=psrot[illt][isa];
			real ct=cos(theta);
			real st=sin(theta);
			P(petf,isa,illt)=cnew(notfx,notfy);
			cmat *etf2d=P(petf,isa,illt);
#pragma omp parallel for default(shared)
			for(int icompy=0; icompy<notfy; icompy++){
			    real iy=(icompy-notfy2);
			    for(int icompx=0; icompx<notfx; icompx++){
				real ix=(icompx-notfx2);
				real ir=(dusc*(ct*ix+st*iy))+netf2;/*index in etf */
				int iir=ifloor(ir);
				ir=ir-iir;
				if(iir>=0 && iir<netf-1){
				    /*bilinear interpolation. */
				    P(etf2d,icompx,icompy)=etf->p[iir]*(1.-ir)
					+etf->p[iir+1]*ir;
				}/*else{P(etf2d,icompx,icompy)=0;}*/
			    }
			}
			cfftshift(P(petf,isa,illt));/*peak in corner; */
		    }else{
			warning_once("Wrong focus!\n");
			P(petf,isa,illt)=cnew(notfx,notfy);
			cset(P(petf,isa,illt),1);
		    }
		}//for isa
	    }//for illt.
	    cfree(etf);
	    free(thetas);
	}//if na_interp
	//fuse nominal to etf to avoid multiply again.
	ccell*  pnominal=dtfs[iwvl].nominal;
	int mnominal=0;/*multiply with isa to get index into pnominal. */
	if(dtfs[iwvl].nominal->nx>1){
	    mnominal=1;//polar ccd
	}
	int mllt=0;
	if(dtfs[iwvl].nominal->ny>1){
	    mllt=1;
	}
	for(int illt=0; illt<nllt; illt++){
	    for(int isa=0; isa<nsa; isa++){
		ccwm(P(petf,isa,illt), P(pnominal,isa*mnominal,illt*mllt));
	    }
	}
	dtfs[iwvl].fused=1;
	
    }//for iwvl
    return etfs;
}

void dtf_free_do(DTF_T *dtfs){
    if(!dtfs) return;
    for(int iwvl=0;iwvl<dtfs->nwvl;iwvl++){
	ccellfree(dtfs[iwvl].nominal);
	dspcellfree(dtfs[iwvl].si);
	cfree(dtfs[iwvl].Ux);
	cfree(dtfs[iwvl].Uy);
    }
    free(dtfs);
}

void etf_free_do(ETF_T *etfs){
    if(!etfs) return;
    for(int iwvl=0;iwvl<etfs->nwvl;iwvl++){
	ccellfree(etfs[iwvl].etf);
    }
    free(etfs);
}
/**
   profile is nxm array with m>=2. First column is coordinate and the remaining
   is function of coordinate. This routine smoothes the profile on new
   coordinate with spacing of dxnew.
 */
dmat* smooth(const dmat *prof, real dxnew){
    const long nxin=prof->nx;
    const real x0in=prof->p[0];
    const real dxin=(prof->p[nxin-1]-x0in)/(nxin-1);
    dmat *out;
    if(dxnew > dxin * 2){
	const long nxnew=ceil((prof->p[nxin-1]-x0in)/dxnew);
	loc_t *loc_in=mk1dloc_vec(prof->p, nxin);
	loc_t *loc_out=mk1dloc(x0in, dxnew, nxnew);
	dsp *ht=mkhb(loc_out, loc_in, 0, 0, 1);
	out=dnew(nxnew, prof->ny);
	memcpy(out->p, loc_out->locx, sizeof(real)*nxnew);
#pragma omp parallel for default(shared)
	for(long icol=1; icol<prof->ny; icol++){
	    /*input profile */
	    const real *pin=PCOL(prof,icol);
	    /*output profile */
	    real *pout=PCOL(out, icol);
	    /*preserve sum of input profile */
	    real Nasum=dvecsum(pin, nxin);
	    dspmulvec(pout, ht, pin, 'n', 1);
	    dnormalize_sumabs(pout, nxnew, Nasum);
	}
	dspfree(ht);
	locfree(loc_in);
	locfree(loc_out);
    }else{
	out=dref(prof);
    }
    return out;
}
