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
/**
   \file setup_powfs.h
   
   Setting up WFS geometry. like the subaperture location, subaperture grid
   points, physical optics detection transfer function, LGS elongation transfer
   function, etc.

   \todo isolate DTF, ETF, routines, make then generic interface and relocate to the lib folder.

   Do not use sparse interpolation to replace ray tracing for fine sampled
   destination grid, especially cubic splines. The interpolation marix takes too
   much space.  

   TODO: This routine and POWFS_T should only contain information about the
   simulation, not about any model used during reconstruction (RTC) to avoid
   leaking information from the "real world (simulation)" to our knowledge (RTC).
*/

#ifndef AOS_POWFS_H
#define AOS_POWFS_H
#include "common.h"
POWFS_T * setup_powfs_init(const PARMS_T *parms, APER_T *aper);
void setup_powfs_misreg_tel(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs);
void setup_powfs_misreg_dm(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs);
void setup_powfs_phy(const PARMS_T *parms,  POWFS_T *powfs);
void setup_powfs_neasim(const PARMS_T *parms,  POWFS_T *powfs);
void setup_powfs_calib(const PARMS_T *parms, POWFS_T *powfs);
void free_powfs_unused(const PARMS_T *parms, POWFS_T *powfs);
void free_powfs(const PARMS_T *parms, POWFS_T *powfs);
//void test_powfs(const PARMS_T *parms, POWFS_T *powfs);
void setup_powfs_etf(POWFS_T *powfs, const PARMS_T *parms, double deltah,
		     int ipowfs, int mode, int icol);
void wfspupmask(const PARMS_T *parms, loc_t *loc, dmat *amp, int iwfs);

#endif
