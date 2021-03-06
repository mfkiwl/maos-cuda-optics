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
#ifndef AOS_LIB_DTF_H
#define AOS_LIB_DTF_H
#include "../math/mathdef.h"
/**
   \file mkdtf.h
   Routine to generate detector transfer function and elongation transfer function due to sodium layer.
*/

/**
   contains the data associated with a detector transfer function for a
   subaperture. The PSF is computed as
   \f$\textrm{PSF}=\frac{1}{N^2\sum(\textrm{amp}^2)}|\textrm{fftshift}\mathcal{F}[A
   \exp(-\frac{2\pi}{\lambda}\textrm{opd})]|^2\f$.  The subaperture image is
   computed as
   \f$I=\textrm{si}*\mathcal{F}^{-1}[\mathcal{F}[\textrm{PSF}\times\textrm{nominal}]]\f$
*/
typedef struct DTF_T{
    ccell *nominal;      /**<The FFT of the pixel functions*/
    dspcell *si;         /**<The pixel selection*/
    real wvl;          /**<Wavelength*/
    real dtheta;       /**<Sampling of PSF*/
    cmat *Ux;            /**<Special frequency vector along x*/
    cmat *Uy;            /**<Special frequency vector along y*/
    real dxsa;         /**<Subaperture size*/
    long notfx;         /**<FFT size along x*/
    long notfy;         /**<FFT size along y*/
    int radpix;          /**<1: Pixels are along radial/azimuthal direction*/
    int fused;           /**<Whether the DTF has been fused to ETF*/
    int nwvl;            /**<Number of DTF_T*/
}DTF_T;

typedef struct ETF_T{
    ccell *etf;          /**<Store the 2D ETF when*/
    int icol;            /**<Store the column index*/
    int nwvl;            /**<Number of DTF_T*/
}ETF_T;

DTF_T *mkdtf(const dmat *wvls, /**<List of wavelength*/
	     real dxsa,/**<Subaperture size*/
	     real embfac,/**<Embedding factor (2)*/
	     long notfx,/**<FFT size along x*/
	     long notfy,/**<FFT size along y*/
	     long pixpsax,/**<Number of pixels along x(r)*/
	     long pixpsay,/**<Number of pixels along y(a)*/
	     real pixthetax,/**<Pixel size along x (r)*/
	     real pixthetay,/**<Pixel size along y (a)*/
	     const dmat* pixoffx,  /**<offset of image center from center of detector*/
	     const dmat* pixoffy,  /**<offset of image center from center of detector*/
	     real pixblur,  /**<Pixel blur sigma(fraction of pixel)*/
	     const dcell *srot, /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     int radpix         /**<1: Pixels are along radial/azimuthal direction*/
    );
ETF_T *mketf(DTF_T *dtfs,  /**<The dtfs*/
	     real hs,      /**<Guide star focus range*/
	     const dcell *sodium, /**<The sodium profile. First column is coordinate.*/
	     int icol,     /**<Which sodium profile to use*/
	     const dcell *srot,  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     const dcell *srsa,  /**<Subaperture to LLT distance*/
	     int no_interp /**<Use direct sum instead of interpolation + FFT. Slower */
    );
void dtf_free_do(DTF_T *dtfs);
void etf_free_do(ETF_T *etfs);
/**frees DTF_T */
#define dtf_free(dtfs) ({dtf_free_do(dtfs); dtfs=NULL;})
/**frees ETF_T */
#define etf_free(etfs) ({etf_free_do(etfs); etfs=NULL;})
dmat* smooth(const dmat *profile, real dxnew);
#endif
