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
#include "common.h"
#include "sim_utils.h"
#include "maos.h"
#include "version.h"

//Return maos data pointer by name. Currently the numbers are filled manually. The next step is to automatically generate them from maos/types.h

static void *find_var(const char *name){
    struct VAR_MAP{
	char *name;
	void *var;
    };
#define VAR_GET(STRUCT,FIELD) {#STRUCT "." #FIELD, global->STRUCT->FIELD}
#define VAR_GET_2(STRUCT,FIELD,NUM) {#STRUCT "[" #NUM "]." #FIELD, global->STRUCT[NUM].FIELD}
    struct VAR_MAP simu_map[]={
	VAR_GET(simu,atm),
	VAR_GET(simu,wfsopd),
	VAR_GET(simu,ints),
	VAR_GET(simu,gradcl),
	VAR_GET(simu,gradacc),
	VAR_GET(simu,gradlastcl),
	VAR_GET(simu,gradlastol),
	VAR_GET(simu,cn2res),
	VAR_GET(simu,gradoff),
	VAR_GET(simu,gradscale),
	VAR_GET(simu,opdr),
	VAR_GET(simu,cgres),
	VAR_GET(simu,clem),
	VAR_GET(simu,clemp),
	VAR_GET(simu,corrNGSm),
	VAR_GET(simu,cleNGSm),
	VAR_GET(simu,dmpsol),
	VAR_GET(simu,dmcmd),
	VAR_GET(simu,dmreal),
	VAR_GET(simu,dmadd),
	VAR_GET(simu,dmfit),
	VAR_GET(simu,dmerr),
	VAR_GET(simu,dmerr_store),
	VAR_GET(simu,Merr_lo),
	VAR_GET(simu,Merr_lo2),
	VAR_GET(simu,Mngs),
	VAR_GET(simu,fsmerr),
	VAR_GET(simu,fsmreal),
	VAR_GET(simu,LGSfocus),
	VAR_GET(simu,zoomerr),
	VAR_GET(simu,zoomint),
	VAR_GET(simu,zoomreal),
	VAR_GET(simu,evlopd),
	{NULL, NULL}//mark the end
    };
    
    struct VAR_MAP aper_map[]={
	VAR_GET(aper,locs),
	VAR_GET(aper,locs_dm),
	VAR_GET(aper,amp),
	VAR_GET(aper,amp1),
	VAR_GET(aper,ampground),
	VAR_GET(aper,mod),
	VAR_GET(aper,mcc),
	VAR_GET(aper,imcc),
	VAR_GET(aper,opdadd),
	VAR_GET(aper,opdfloc),
	{NULL, NULL}
    };

    struct VAR_MAP recon_map[]={
	VAR_GET(recon,ht),
	VAR_GET(recon,wt),
	VAR_GET(recon,os),
	VAR_GET(recon,dx),
	VAR_GET(recon,saloc),
	VAR_GET(recon,ploc),
	VAR_GET(recon,pmap),
	VAR_GET(recon,ploc_tel),
	VAR_GET(recon,xloc),
	VAR_GET(recon,xmap),
	VAR_GET(recon,xcmap),
	VAR_GET(recon,xmcc),
	VAR_GET(recon,fmap),
	VAR_GET(recon,floc),
	VAR_GET(recon,aloc),
	VAR_GET(recon,amap),
	VAR_GET(recon,acmap),
	VAR_GET(recon,actfloat),
	VAR_GET(recon,actstuck),
	VAR_GET(recon,amod),
	VAR_GET(recon,W0),
	VAR_GET(recon,W1),
	VAR_GET(recon,L2),
	VAR_GET(recon,GP),
	VAR_GET(recon,HXW),
	VAR_GET(recon,HXWtomo),
	VAR_GET(recon,GX),
	VAR_GET(recon,GXtomo),
	VAR_GET(recon,GXlo),
	VAR_GET(recon,GXL),
	VAR_GET(recon,GA),
	VAR_GET(recon,GAlo),
	VAR_GET(recon,GAhi),
	VAR_GET(recon,GM),
	VAR_GET(recon,GMhi),
	VAR_GET(recon,HA_ncpa),
	{NULL, NULL}
    };
    struct VAR_MAP powfs_map[]={
	VAR_GET_2(powfs,saloc,0),
	VAR_GET_2(powfs,pts,0),
	VAR_GET_2(powfs,saa,0),
	VAR_GET_2(powfs,loc,0),
	VAR_GET_2(powfs,amp,0),
	VAR_GET_2(powfs,srot,0),
	VAR_GET_2(powfs,srsa,0),

	VAR_GET_2(powfs,neasim,0),
	VAR_GET_2(powfs,opdadd,0),
	{NULL,NULL}
    };
    struct MAP_MAP{
	char *name;
	struct VAR_MAP *map;
    };
#define MAP_ALL(NAME) {#NAME "_", NAME ## _map}
    struct MAP_MAP map_map[]={
	MAP_ALL(simu),
	MAP_ALL(aper),
	MAP_ALL(recon),
	MAP_ALL(powfs),
	{NULL, NULL}
    };
    if(global){
	if(global->simu){
	    if(name){
		for(int j=0; map_map[j].name;  j++){
		    char *div=strchr(name, '.');
		    if(!div) continue;
		    char *div2=strchr(name,'[');
		    if(div2 && div2<div){
			div=div2;
		    }
		    if(!strncmp(map_map[j].name, name, div-name)){//first find the correct map
			struct VAR_MAP *var_map=map_map[j].map;
			for(int i=0; var_map[i].name ; i++){
			    if(!strcmp(var_map[i].name, name)){//then find the correct variable
				info("%s is found at %p\n", name, var_map[i].var);
				return var_map[i].var;
			    }
			}
			info("%s not found.\n", name);
		    }
		}
	    }
	    {
		static cell* dummy=NULL;
		if(!dummy){
		    dummy=cellnew(0,0);
		    const char *msg0="Available variables are:\n";
		    long count=strlen(msg0);
		    for(int j=0; map_map[j].name;  j++){
			struct VAR_MAP *var_map=map_map[j].map;
			for(int i=0; var_map[i].name ; i++){
			    count+=strlen(var_map[i].name)+1;
			}
		    }
		    dummy->header=mycalloc(count, char);
		    strcat(dummy->header, msg0);
		    for(int j=0; map_map[j].name;  j++){
			struct VAR_MAP *var_map=map_map[j].map;
			for(int i=0; var_map[i].name ; i++){
			    strcat(dummy->header, var_map[i].name);
			    strcat(dummy->header, "\n");
			}
		    }
		}
		return dummy;
	    }
	}else{
	    warning("global->simu is NULL\n");
	}
    }else{
	warning("global is NULL\n");
    }
    return NULL;
}
//Listen to requests coming from other client.
static void *maos_var(void* psock){
    thread_block_signal();
    int cmd[2];
    int sock=(int)(long)psock;
    while(!streadintarr(sock, cmd, 2)){
	info("maos_var: cmd=%d, %d\n", cmd[0], cmd[1]);
	switch(cmd[0]){
	case MAOS_VAR:
	    {
		if(cmd[1]==1 || cmd[1]==2){//get or put
		    char *name=NULL;
		    streadstr(sock, &name);
		    cell *var=(cell*)find_var(name);
		    info("maos_var: request[%d] %s %p\n", cmd[1], name, var);
		    {
			if(cmd[1]==1){//client to get
			    writesock(var, sock);
			}else if(cmd[1]==2){//client to put
			    cell *newdata=readsock(sock);
			    dcellcp(&var, newdata);
			}else{
			    warning("maos_var: unknown operation %d\n", cmd[1]);
			}
		    }
		}
	    }
	    break;
	case MAOS_PAUSE:
	    {
		if(global && global->simu){
		    global->simu->pause=cmd[1];
		    putchar('\a');
		}
	    }
	    break;
	}//switch
    }//while
    info("maos_var: client closed\n");
    close(sock);
    return NULL;
}

//Listen to commands coming from scheduler
static void* maos_listener(void *psock){
    int sock=(int)(long)psock;
    thread_block_signal();
    int cmd[2];
    while(!streadintarr(sock, cmd, 2)){
	switch(cmd[0]){
	case MAOS_DRAW://Starting draw to received fd.
	    {
		int fd;
		if(streadfd(sock, &fd)){
		    warning("unable to read fd from %d\n", sock);
		    continue;
		}else{
		    info("got fd=%d\n", fd);
		}
		draw_add(fd);
		if(global){
		    PARMS_T *parms=(PARMS_T*)global->parms;//cast away constness
		    parms->plot.setup=1;
		    parms->plot.run=1;
		    if(global->setupdone){//setup is already finished. request plot setup.
			plot_setup(global->parms, global->powfs, global->aper, global->recon);
		    }
		}
	    }break;
	case MAOS_VAR:
	    {
		int fd;
		if(streadfd(sock, &fd)){
		    warning("unable to read fd from %d\n", sock);
		    continue;
		}else{
		    info("got fd=%d\n", fd);
		}
		thread_new(maos_var, (void*)(long)fd);
		
	    }break;

	default:
	    warning("unknown cmd %d\n", cmd[0]);
	    break;
	}
    }
    return NULL;
}
void maos_version(void){
    info2("SRC: %s v%s %s\n", SRCDIR, PACKAGE_VERSION, GIT_VERSION);
    info2("BUILD: %s by %s on %s %s", BUILDDIR, COMPILER, __DATE__, __TIME__);
#if CPU_SINGLE 
    info2(" CPU(single)");
#else
    info2(" CPU(double)");
#endif
#if USE_CUDA
#if CUDA_DOUBLE
    info2(" +CUDA(double)");
#else
    info2(" +CUDA(single)");
#endif
#else
    info2(" -CUDA");
#endif
#ifdef __OPTIMIZE__
    info2(" +optimization.\n");
#else
    info2(" -optimization\n");
#endif
    info2("Launched at %s in %s with PID %ld.\n",myasctime(),HOST, (long)getpid());
#if HAS_LWS
    extern uint16_t PORT;
    info2("The web based job monitor can be accessed at http://localhost:%d\n", 1+PORT);
#endif
}

/**
   This is the standard entrance routine to the program.  It first calls
   setup_parms() to setup the simulation parameters and check for possible
   errors. It then waits for starting signal from the scheduler if in batch
   mode. Finally it hands the control to maos() to start the actual simulation.

   Call maos with overriding *.conf files or embed the overriding parameters in
   the command line to override the default parameters, e.g.

   <p><code>maos base.conf save.setup=1 'powfs.phystep=[0 100 100]'</code><p>

   Any duplicate parameters will override the pervious specified value. The
   configure file default.conf will be loaded as the master .conf unless a -c
   switch is used with another .conf file. For scao simulations, call maos with
   -c switch and the right base .conf file.
   
   <p><code>maos -c scao_ngs.conf override.conf</code><p>

   for scao NGS simulations 

   <p><code>maos -c scao_lgs.conf override.conf</code><p>

   for scao LGS simulations.  With -c switch, default.conf will not be read,
   instead scao_ngs.conf or scao_lgs.conf are read as the master config file.
   Do not specify any parameter that are not understood by the code, otherwise
   maos will complain and exit to prevent accidental mistakes.
       
   Generally you link the maos executable into a folder that is in your PATH
   evironment or into the folder where you run simulations.

   Other optional parameters:
   \verbatim
   -d          do detach from console and not exit when logged out
   -s 2 -s 4   set seeds to [2 4]
   -n 4        launch 4 threads.
   -f          To disable job scheduler and force proceed
   \endverbatim
   In detached mode, drawing is automatically disabled.
   \callgraph
*/
int main(int argc, const char *argv[]){
    char *scmd=argv2str(argc,argv," ");
    ARG_T* arg=parse_args(argc,argv);/*does chdir */
    if(arg->detach){
	daemonize();
    }else{
	redirect();
    }
    /*Launch the scheduler if it is not running and report about our process */
    int ngpu;
#if USE_CUDA
    ngpu=arg->ngpu;
    if(!ngpu) ngpu=0xFFFFFF; 
    else if(ngpu<0) ngpu=0;
#else
    ngpu=0;
#endif
    scheduler_start(scmd,NTHREAD,ngpu,!arg->force);
    info("%s\n", scmd);
    info("Output folder is '%s'. %d threads\n",arg->dirout, NTHREAD);
    free(arg->dirout); arg->dirout = NULL;

    maos_version();
    /*setting up parameters before asking scheduler to check for any errors. */
    PARMS_T *parms=setup_parms(arg->conf, arg->confcmd, arg->over);
    free(arg->conf); arg->conf=NULL;
    if(arg->confcmd){
	remove(arg->confcmd); free(arg->confcmd); arg->confcmd=NULL;
    }
    if(parms->sim.nseed>0){
	/*register signal handler */
	register_signal_handler(maos_signal_handler);
 
	if(!arg->force){
	    /*
	      Ask job scheduler for permission to proceed. If no CPUs are
	      available, will block until ones are available.
	      if arg->force==1, will run immediately.
	    */
	    info("Waiting start signal from the scheduler ...\n");
	    int count=0;
	    while(scheduler_wait()&& count<60){
		/*Failed to wait. fall back to own checking.*/
		warning_time("failed to get reply from scheduler. retry\n");
		sleep(10);
		count++;
		scheduler_start(scmd,NTHREAD,ngpu,!arg->force);
	    }
	    if(count>=60){
		warning_time("fall back to own checker\n");
		wait_cpu(NTHREAD);
	    }
	}
	if(scheduler_listen(maos_listener)){
	    info("Failed to start maos_listener\n");
	}
	setup_parms_gpu(parms, arg->gpus, arg->ngpu);
    
 
	/* do not use prallel single in maos(). It causes blas to run single threaded
	 * during preparation. Selective enable parallel for certain setup functions
	 * that doesn't use blas*/

	info("\n*** Preparation started at %s in %s. ***\n\n",myasctime(),HOST);
	maos_setup(parms);
	if(parms->sim.end>parms->sim.start){
	    info("\n*** Simulation  started at %s in %s. ***\n\n",myasctime(),HOST);
	    maos_sim();
	}
	rename_file(0);
    }
    
    maos_reset();
    free_parms(parms);
    info("\n*** Simulation finished at %s in %s. ***\n\n",myasctime(),HOST);
    free(scmd);
    free(arg->gpus);
    free(arg);
    scheduler_finish(0);
    return 0;
}
