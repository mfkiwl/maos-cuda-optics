#include "../lib/aos.h"
/*
static void test_w1(){
    loc_t *loc=mksqloc2(200,200,1./64.);
    double *amp=calloc(loc->nloc, sizeof(double));
    loccircle(amp,loc,0,0, 1,1);
    drawopd("test_loc",loc,amp,"loc");
    normalize(amp,loc->nloc,1);
    dsp *W0=NULL;
    dmat *W1=NULL;
    mkw_amp(loc,amp,&W0,&W1);
    info("Sum W1=%g",dsum(W1));
    writedbl(amp,loc->nloc,1,"amp");
    dwrite(W1,"W1");
    spwrite(W0,"W0");
    locwrite(loc,"loc");
}
static void test_wploc(){
    loc_t *loc=locread("ploc.bin");
    double *amp=NULL; long nx,ny;
    readdbl(&amp,&nx,&ny,"pamp.bin");
    double ampmax=maxdbl(amp,nx);
    double ampmax1=1./ampmax;
    for(int i=0; i<nx; i++){
	amp[i]=round(amp[i]*ampmax1)*ampmax;
    }
    normalize(amp,nx,1);
    drawopd("test_loc",loc,amp,"loc");
    normalize(amp,loc->nloc,1);
    dsp *W0=NULL;
    dmat *W1=NULL;
    mkw_amp(loc,amp,&W0,&W1);
    info("Sum W1=%g\n",dsum(W1));
    writedbl(amp,loc->nloc,1,"amp");
    dwrite(W1,"pW1");
    spwrite(W0,"pW0");
}
static void test_wcir(){
    loc_t *loc=mksqloc2(100,100,1./4.);
    dsp *W0=NULL;
    dmat *W1=NULL;
    mkw_circular(loc,0,0,10,&W0,&W1);
    locwrite(loc,"loc");
    spwrite(W0,"W0");
    dwrite(W1,"W1");
}
static void test_int_reg(){
    int n=100;
    double dx=1./n;
    double cy=0;
    double cx=0;
    double cr=0.3;
    double cr2=cr*cr;
    int area=0;
    for(int i=0; i<n; i++){
	double y=dx*(i+0.5);
	double r2y=pow(y-cy,2);
	for(int j=0; j<n; j++){
	    double x=dx*(j+0.5);
	    if(pow(x-cx,2)+r2y < cr2){
		area++;
	    }
	}
    }
    double ar=(double)area/(double)(n*n);
    info("Area=%g, Err=%g\n",ar,ar-M_PI*cr2/4);
}
static void test_int_rand(){
    int n=10000;
    double cy=0;
    double cx=0;
    double cr=0.3;
    double cr2=cr*cr;
    int area=0;
    rand_t stat;
    seed_rand(&stat,1);
    for(int i=0; i<n; i++){
	double x=drand48();
	double y=drand48();
	if(pow(x-cx,2)+pow(y-cy,2) <= cr2){
	    area++;
	}
    } 
    double ar=(double)area/(double)(n);
    info("Area=%g, Err=%g\n",ar,ar-M_PI*cr2/4);
}
*/
/*
static void test_sqlocrot(void){
    int nn=64;
    double dx=1./64.;
    loc_t *loc=mksqlocrot(nn,nn,dx,-nn/2*dx,-nn/2*dx,0);
    locwrite(loc,"loc_0deg.bin");
    loc_t *loc2=mksqlocrot(nn,nn,dx,-nn/2*dx,-nn/2*dx,M_PI/10);
    locwrite(loc2,"loc_18deg.bin");
    }*/
/*
void test_loc_reduce_sp(void){
    int nloc;
    loc_t **xloc=locarrread(&nloc,"xloc.bin");
    spcell *G0=spcellread("G0.bin");
    loc_reduce_sp(xloc[0],G0->p[0],2,1);
    locarrwrite(xloc,nloc,"xloc2");
    spcellwrite(G0,"G02.bin");
    locarrfree(xloc,nloc);
    spcellfree(G0);
    G0=spcellread("G0.bin");
    xloc=locarrread(&nloc,"xloc.bin");
    dsp *G0t=sptrans(G0->p[0]);
    loc_reduce_sp(xloc[0],G0t,1,1);
    dsp *G03=sptrans(G0t);
    spwrite(G03,"G03.bin");
    locarrwrite(xloc,nloc,"xloc3");
    }

static void test_loc_reduce_spcell(void){
    int nloc;
    loc_t **xloc=locarrread(&nloc,"xloc.bin");
    spcell *G0=spcellread("G0.bin");
    spcell *G0t=spcelltrans(G0);
    loc_reduce_spcell(xloc[0],G0t,1,1);
    spcell *G04=spcelltrans(G0t);
    spcellwrite(G04,"G05.bin");
    locarrwrite(xloc,nloc,"xloc5");
    }
static void test_zernike(){
    loc_t *loc=locread("loc_circle.bin");
    dmat *mod=zernike(loc,1,10);
    locwrite(loc,"loc");
    dwrite(mod,"mod");
    dgramschmidt(mod,NULL);
    dwrite(mod,"mod_norm");
    }
static void test_embed(){
    loc_t *loc=mkcirloc(10, 1./64);
    locwrite(loc, "loccir");
    dmat *opd=dnew(645,645);
    dembed_locstat(&opd, loc, NULL);
    dwrite(opd, "locopd");
    }*/
int main(){
    /*  
    test_zernike();
    test_loc_reduce_spcell();
    test_sqlocrot();
	test_w1();
    test_wploc();
    test_wcir();
    test_int_reg();
    test_int_rand();
    test_embed();
    */
    for(long i=1; i<1000;i++){
	info("i=%ld, %ld\n", i, nextpow2(i));
    }
}
