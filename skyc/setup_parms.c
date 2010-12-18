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

#include "skyc.h"
#include "parms.h"
#include "nafocus.h"
#undef malloc
#undef calloc
#undef free
#undef realloc
static void setup_config(const ARG_S*arg){
    /**
       This routine setups the internally mained list of PATHS in
       which all the configurations files are looked for.

       Environment variable AOS_CONFIG_PATH will be added if availab.e.
     */
    open_config(arg->conf,NULL,0);
    /*Read in any additional overriding .conf files. They should
      contain keys already existed in nfiraos.conf*/

    if(arg->iconf<arg->argc){
	char fntmp[PATH_MAX];
	snprintf(fntmp,PATH_MAX,"/tmp/skyc_%ld.conf",(long)getpid());
	FILE *fptmp=fopen(fntmp,"w");
	int inline_conf=0;
	for(int iconf=arg->iconf; iconf<arg->argc; iconf++){
	    char *fno=arg->argv[iconf];
	    //info("Embed (%s)\n",fno);
	    if(index(fno,'=')){
		inline_conf++;
		fprintf(fptmp,"%s\n",fno);
	    }else if(check_suffix(fno,".conf")){
		open_config(fno,NULL,1);/*1 means protected. will not be overriden by
				     base .conf's, but can be overriden by user
				     supplied options.*/
	    }else{
		error("Invalid command line option: %s\n",fno);
	    }
	}
	fclose(fptmp);
	if(inline_conf>0){
	    open_config(fntmp,NULL,1);
	}
	if(remove(fntmp)){
	    perror("remove");
	    warning("Unable to remove file %s\n",fntmp);
	}
    }else{
	warning2("No override file is used\n");
    }
}
#define READ_INT(A) parms->A = readcfg_int(#A) //read a key with int value.
#define READ_DBL(A) parms->A = readcfg_dbl(#A) //read a key with double value
#define READ_STR(A) parms->A = readcfg_str(#A) //read a key with string value.
static int readcfg_dblarr_n(int n, double **p, const char*format,...){
    /**
       Read in double array of n elements.
     */
    format2key;
    int m=readcfg_dblarr(p,"%s",key);
    if(m!=n){
	error("Wrong # of elements for 'key' %s.\n"
	      "Required %d, got %d\n",key,n,m);
    }
    return 0;
}
static int readcfg_intarr_n(int n,int **p, const char*format,...){
    /**
       Read in int array of n elements.
     */
    format2key;
    int m=readcfg_intarr(p,"%s",key);
    if(m!=n){
	error("Wrong # of elements for 'key' %s.\n"
	      "Required %d, got %d\n",key,n,m);
    }
    return 0;
}
static int readcfg_strarr_n(int n,char ***p, const char*format,...){
    /**
       Read in int array of n elements.
     */
    format2key;
    int m=readcfg_strarr(p,"%s",key);
    if(m!=n){
	error("Wrong # of elements for 'key' %s.\n"
	      "Required %d, got %d\n",key,n,m);
    }
    return 0;
}
static void setup_parms_skyc(PARMS_S *parms){
    READ_INT(skyc.dbg);
    READ_INT(skyc.dbgsky);
    READ_INT(skyc.dbgaster);
    READ_INT(skyc.verbose);
    READ_INT(skyc.save);
    READ_INT(skyc.start);
    READ_INT(skyc.nsky);
    READ_DBL(skyc.lat);
    READ_DBL(skyc.lon);
    READ_DBL(skyc.catscl);
    READ_DBL(skyc.patfov);
    READ_DBL(skyc.keepout);
    READ_INT(skyc.ngsalign);
    READ_INT(skyc.npowfs);
    READ_INT(skyc.psd_scale);
    READ_INT(skyc.noisy);
    READ_INT(skyc.demote);
    READ_INT(skyc.ttfbrightest);
    READ_INT(skyc.bspstrehl);
    READ_INT(skyc.maxaster);
    READ_INT(skyc.maxstar);
    readcfg_intarr_n(parms->skyc.npowfs,&parms->skyc.nwfsmax,"skyc.nwfsmax");
    readcfg_dblarr_n(parms->skyc.npowfs,&parms->skyc.pixratio,"skyc.pixratio");
    readcfg_dblarr_n(parms->skyc.npowfs,&parms->skyc.pixblur,"skyc.pixblur");
    readcfg_intarr_n(parms->skyc.npowfs,&parms->skyc.pixpsa,"skyc.pixpsa");
    readcfg_dblarr_n(parms->skyc.npowfs,&parms->skyc.pixoffx,"skyc.pixoffx");
    readcfg_dblarr_n(parms->skyc.npowfs,&parms->skyc.pixoffy,"skyc.pixoffy");
    READ_INT(skyc.limitnstep);
    READ_DBL(skyc.intgain);
    READ_DBL(skyc.rne);
    READ_INT(skyc.noisefull);
    READ_DBL(skyc.imperrnm);
    READ_INT(skyc.mtchcr);
    READ_INT(skyc.mtch);
    parms->skyc.nwvl=readcfg_dblarr(&parms->skyc.qe,"skyc.qe");
    readcfg_dblarr_n(parms->skyc.nwvl, &parms->skyc.telthruput, "skyc.telthruput");
    parms->skyc.ndtrat=readcfg_intarr(&parms->skyc.dtrats,"skyc.dtrats");
    parms->skyc.nseed=readcfg_intarr(&parms->skyc.seeds, "skyc.seeds");
    READ_INT(skyc.servo);
    READ_INT(skyc.gsplit);
    READ_INT(skyc.evlstart);
    READ_INT(skyc.phystart);
    READ_INT(skyc.gradnea);
    char *temp, *temp1;
    temp=readcfg_str("skyc.psd_ws"); temp1=find_file(temp);
    parms->skyc.psd_ws=dread("%s",temp1); free(temp); free(temp1);
 
    READ_DBL(skyc.zc_f);
    READ_DBL(skyc.zc_zeta);
    READ_DBL(skyc.na_alpha);
    READ_DBL(skyc.na_beta);

    READ_STR(skyc.stars);
}
/**
   Setup infromation output from maos run.
*/
static void setup_parms_maos(PARMS_S *parms){
    READ_DBL(maos.r0z);
    READ_DBL(maos.dt);
    READ_DBL(maos.zadeg); parms->maos.za=parms->maos.zadeg*M_PI/180.;
    READ_DBL(maos.hc);
    READ_DBL(maos.hs);
    READ_DBL(maos.D);
    READ_INT(maos.nmod);
    parms->maos.nseed=readcfg_intarr(&parms->maos.seeds, "maos.seeds");

    parms->maos.nwvl=readcfg_dblarr(&parms->maos.wvl, "maos.wvl");
    READ_INT(maos.npowfs);
    readcfg_intarr_n(parms->maos.npowfs,&parms->maos.msa,"maos.msa");
    readcfg_intarr_n(parms->maos.npowfs,&parms->maos.nsa,"maos.nsa");

    readcfg_intarr_n(parms->maos.npowfs,&parms->maos.ncomp,"maos.ncomp");
    readcfg_intarr_n(parms->maos.npowfs,&parms->maos.embfac,"maos.embfac");
    readcfg_dblarr_n(parms->maos.npowfs,&parms->maos.dxsa,"maos.dxsa");

    readcfg_strarr_n(parms->maos.npowfs, &parms->maos.fnwfsloc, "maos.fnwfsloc");
    readcfg_strarr_n(parms->maos.npowfs, &parms->maos.fnwfsamp, "maos.fnwfsamp");
    readcfg_strarr_n(parms->maos.npowfs, &parms->maos.fnsaloc, "maos.fnsaloc");
    char *temp, *temp1;
    temp=readcfg_str("maos.fnmcc"); temp1=find_file(temp);
    parms->maos.mcc=dread("%s",temp1); free(temp); free(temp1);
    parms->maos.mcc_tt=dsub(parms->maos.mcc,0,2,0,2);

    temp=readcfg_str("maos.fnmcc_oa"); temp1=find_file(temp);
    parms->maos.mcc_oa=dread("%s",temp1); free(temp); free(temp1);
    parms->maos.mcc_oa_tt=dsub(parms->maos.mcc_oa, 0, 2, 0, 2);
    parms->maos.mcc_oa_tt2=dsub(parms->maos.mcc_oa,0, 2, 0, 5);
    
    READ_STR(maos.fnmideal);
    READ_STR(maos.fnmidealp);
    READ_INT(maos.evlindoa);
    READ_DBL(maos.ngsgrid);
    READ_INT(maos.nstep);
    if(readcfg_peek("maos.wddeg")){
	parms->maos.nwddeg=readcfg_dblarr(&parms->maos.wddeg,"maos.wddeg");
    }else{
	parms->maos.nwddeg=0;
	parms->maos.wddeg=NULL;
    }
}
PARMS_S *setup_parms(const ARG_S *arg){
    setup_config(arg);
   
    PARMS_S *parms=calloc(1, sizeof(PARMS_S));
    parms->skyc.nthread=arg->nthread;
    setup_parms_skyc(parms);
    setup_parms_maos(parms);
    
    info2("There are %d MAOS PSF seeds: ",parms->maos.nseed);
    for(int i=0; i<parms->maos.nseed; i++){
	info2(" %d", parms->maos.seeds[i]);
    }
    info2("\n");
    info2("There are %d asterism seeds:", parms->skyc.nseed);
    for(int i=0; i<parms->skyc.nseed; i++){
	info2(" %d", parms->skyc.seeds[i]);
    }
    info2("\n");
    const double wvlmin=parms->maos.wvl[0];
    parms->skyc.pixtheta=calloc(parms->skyc.npowfs, sizeof(double));
    for(int ipowfs=0; ipowfs<parms->skyc.npowfs; ipowfs++){
	const double dxsa=parms->maos.dxsa[ipowfs];
	parms->skyc.pixtheta[ipowfs]=parms->skyc.pixratio[ipowfs]*wvlmin/dxsa;
	info2("powfs %d, pixtheta=%g mas\n", ipowfs, parms->skyc.pixtheta[ipowfs]*206265000);
    }
    parms->skyc.fss=calloc(parms->skyc.ndtrat, sizeof(double));
    parms->skyc.rnefs=dnew(parms->skyc.ndtrat, parms->maos.npowfs);
    if(parms->skyc.rne<0){
	PDMAT(parms->skyc.rnefs,rnefs);
	info2("Using frame rate dependent read out noise:\n");
	for(long idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
	    int dtrat=parms->skyc.dtrats[idtrat];
	    double fs=1./(parms->maos.dt*dtrat);
	    parms->skyc.fss[idtrat]=fs;
	    info2("%5.1f Hz: ", fs);
	    for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
		int N=parms->skyc.pixpsa[ipowfs];
		double raw_frame_time=((0.00604*N+0.01033)*N+0.00528)*1e-3;
		double K=1./(raw_frame_time*fs);
		rnefs[ipowfs][idtrat]=sqrt(pow(10.6*pow(K,-0.45),2) + 2.7*2.7 + 0.0034*K);
		info2("%5.1f ",rnefs[ipowfs][idtrat]);
	    }
	    info2("\n");
	}
    }else{
	info2("Using constant read out noise of %g\n", parms->skyc.rne);
	dset(parms->skyc.rnefs, parms->skyc.rne);
    }
    parms->skyc.resfocus=dnew(parms->skyc.ndtrat, 1);
    for(long idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
	int dtrat=parms->skyc.dtrats[idtrat];
	double fs=1./(parms->maos.dt*dtrat);

	parms->skyc.resfocus->p[idtrat]=
	    pow(nafocus_residual(fs, parms->maos.dt, parms->skyc.zc_f, parms->skyc.zc_zeta,
				 parms->maos.D, parms->maos.hs, 
				 parms->skyc.na_alpha, parms->skyc.na_beta),2);
    }
 
    if(parms->skyc.npowfs != parms->maos.npowfs){
	error("skyc.npowfs should match maos.npowfs\n");
    }
    if(parms->skyc.ngsalign){
	warning2("NGS are aligned to grid spaced by %g\"\n", parms->maos.ngsgrid);
    }
    if(0){
	char temp[80]; 
	char *temp1;
	snprintf(temp,80, "PSD/PSD_NGS_r0z_%.4f_za%g.bin.gz",parms->maos.r0z, parms->maos.zadeg);
	temp1=find_file(temp); parms->skyc.psd_ngs=dread("%s",temp1); 
	info2("Loading PSD of NGS modes from %s\n", temp1);free(temp1);

	snprintf(temp,80, "PSD/PSD_TT_r0z_%.4f_za%g.bin.gz",parms->maos.r0z, parms->maos.zadeg);
	temp1=find_file(temp); parms->skyc.psd_tt=dread("%s",temp1); 
	info2("Loading PSD of TT modes from %s\n", temp1);free(temp1);

	snprintf(temp,80, "PSD/PSD_PS_r0z_%.4f_za%g.bin.gz",parms->maos.r0z, parms->maos.zadeg);
	temp1=find_file(temp); parms->skyc.psd_ps=dread("%s",temp1); 
	info2("Loading PSD of PS modes from %s\n", temp1);free(temp1);
    }else{
	char temp[80]; 
	char *temp1;
	snprintf(temp,80, "PSD/PSD_NGS.bin.gz");
	temp1=find_file(temp); parms->skyc.psd_ngs=dread("%s",temp1); 
	info2("Loading PSD of NGS modes from %s\n", temp1);free(temp1);
	snprintf(temp,80, "PSD/PSD_TT.bin.gz");
	temp1=find_file(temp); parms->skyc.psd_tt=dread("%s",temp1); 
	info2("Loading PSD of TT modes from %s\n", temp1);free(temp1);
	snprintf(temp,80, "PSD/PSD_PS.bin.gz");
	temp1=find_file(temp); parms->skyc.psd_ps=dread("%s",temp1); 
	info2("Loading PSD of PS modes from %s\n", temp1);free(temp1);
    }
    char fnconf[PATH_MAX];
    snprintf(fnconf, PATH_MAX, "skyc_%ld.conf", (long)getpid());
    close_config("%s",fnconf);
    mysymlink(fnconf, "skyc_recent.conf");
    if(parms->skyc.gradnea){
	info2("Variance of the gradients in stored PSF is added to NEA\n");
    }
    info2("Maximum number of asterisms in each star field is %d\n", parms->skyc.maxaster);
    
    if(arg->detach || parms->skyc.nthread>1){
	parms->skyc.verbose=0;
    }else if(parms->skyc.verbose<=0){
	parms->skyc.verbose=1;
    }
    parms->skyc.nwfstot=0;
    for(int ipowfs=0; ipowfs<parms->skyc.npowfs; ipowfs++){
	parms->skyc.nwfstot+=parms->skyc.nwfsmax[ipowfs];
    }
    return parms;
}
/**
   Free the data in parms.
 */
void free_parms(PARMS_S *parms){
    dfree(parms->skyc.psd_ngs);
    dfree(parms->skyc.psd_tt);
    dfree(parms->skyc.psd_ws);
    dfree(parms->skyc.psd_ps);
    dfree(parms->maos.mcc);
    dfree(parms->maos.mcc_oa);
    dfree(parms->maos.mcc_tt);
    dfree(parms->maos.mcc_oa_tt);
    dfree(parms->maos.mcc_oa_tt2);
    dfree(parms->skyc.rnefs);
    dfree(parms->skyc.resfocus);
}
