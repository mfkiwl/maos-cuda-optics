/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "../lib/aos.h"
/*
int main(){
    dspcell *FLM=dspcellread("FLM.bin");
    dsp *F=FLM->p[0];
    spchol *Fchol=chol_factorize(F);
 
    rand_t rrand;
    seed_rand(&rrand,1);
    dmat *y=dnew(F->m,10);
    drandu(y,1,&rrand);
    dmat *x=NULL;
    chol_solve(&x,Fchol,y);
    writebin(y,"y");
    writebin(x,"x");
    chol_save(Fchol,"Chol");
    chol_free(Fchol);
    dfree(x);
    dfree(y);
    dspcellfree(FLM);
}
*/
TIC;
int main(int argc, char* argv[]){
    /*dsp *RLMc1=dspread("RLMc_old.bin"); */
    if(argc!=2){
	error("Need 1 argument\n");
    }
    dspcell *RLM=dspcellread("%s",argv[1]);
    dsp *RLMc=dspcell2sp(RLM);
    tic;info("chol ...");
    spchol *R1=chol_factorize(RLMc);
    toc("done");
    rand_t rstat;
    seed_rand(&rstat,1);
    dmat *y=dnew(RLMc->nx, 1);
    drandn(y, 1, &rstat);
    dmat *x=NULL, *x2=NULL, *x3=NULL;
    chol_convert(R1, 1);
    tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_solve_upper(&x3, R1, y);
    toc("upper");tic;
    chol_solve_upper(&x3, R1, y);
    toc("upper");tic;
    chol_solve_lower(&x2, R1,y);
    toc("lower");tic;
    chol_solve_lower(&x2, R1,y);
    toc("lower");tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    writebin(y,"y");
    writebin(x,"x");
    writebin(x2,"x2");
    writebin(x3,"x3");
    chol_free(R1);
    dspfree(RLMc);
    dspcellfree(RLM);
}
