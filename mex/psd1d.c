/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
  Wrap of the function psd1d to mex routines.
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_DATA,
	P_LSEG,
	P_TOT,
    };
    enum{
	PL_PSD,
	PL_TOT,
    }; 
    if(nrhs < P_TOT){
	mexErrMsgTxt("Usage: psd=psd1d(data, lseg)");
    }
    
    dmat *data=mx2d(prhs[P_DATA]);
    long lseg=(long)mxGetScalar(prhs[P_LSEG]);
    dmat *psd;
    if(nrhs > P_TOT){
	double dt=(double)mxGetScalar(prhs[P_TOT]);
	psd=psd1dt(data, lseg, dt);
    }else{
	psd=psd1d(data, lseg);
    }
    dfree(data);
    plhs[PL_PSD]=d2mx(psd);
    dfree(psd);
    dfree(data);
}