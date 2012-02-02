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
#ifndef AOS_CUDA_CUCMAT_H
#define AOS_CUDA_CUCMAT_H

#include "utils.h"
#include "types.h"

#define cucnew  cunew<fcomplex> 
#define cucref  curef<fcomplex>
#define cucfree cufree<fcomplex>
#define cuczero cuzero<fcomplex>
#define cuccellnew  cucellnew<fcomplex>
#define cuccellfree cucellfree<fcomplex>
#define cuccellzero cucellzero<fcomplex>
#define cucwrite     cuwrite<fcomplex, (uint32_t)M_ZMP>
#define cuccellwrite cucellwrite<fcomplex, (uint32_t)M_ZMP>

#endif