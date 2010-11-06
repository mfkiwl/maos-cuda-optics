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

#ifndef SKYC_PHOTON_H
#define SKYC_PHOTON_H

void photon_flux(double *Np, double *Nptot, double *Nbtot, double *QCSNR, double *QCNEA,
		 int nwvl, double* wvls, double *mags, 
		 double dxsa, int iscircle, double pixtheta, 
		 double dt, double za, double *strehl, double imperrnm,
		 double *thruput, double *qe, double rne);
extern const double Z_J;
extern const double Z_H;

#endif