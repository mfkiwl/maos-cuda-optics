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

typedef struct{
    const loc_t *loc; /*reference to the grid*/
    const dmat *amp;  /*reference to the amplitude map*/
    const dmat *wvl;  /*reference to the wavelength*/
    imat  *nembed;/**<size of embedding array (square)*/
    icell *embed; /**<embedding index*/
    double ampsum;/**<sum(amp)*/
    double ampnorm;/**<sum(amp.*amp)*/
    dcell *fieldmask;/**<Masking the PSF in fourier domain*/
    double fieldstop;
}locfft_t;

locfft_t *locfft_init(const loc_t *loc, const dmat *amp, const imat *fftsize,
		      const dmat *wvl, double fieldstop); 
void locfft_free(locfft_t *locfft);
ccell* locfft_psf(locfft_t *locfft, dmat *opd, imat *psfsize);
void locfft_fieldstop(locfft_t *locfft, dmat *opd, dmat *wvlwts);