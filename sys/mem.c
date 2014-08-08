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

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <search.h>
#include <execinfo.h>
#include <sys/stat.h>
#include <unistd.h>

int exit_success=0;
#include "mem.h"
#include "thread.h"
#include "scheduler_client.h"

/*
  Record allocated memory and their size.  Warn if some memory is not freed.
  Notice that the pointer return by tsearch and the pointer passwd to the
  function called by twalk is the addess of the tree node, which points to the
  key. so **key is required.

  We use the backtrace to see the chain of calling when memory is
  allocated. This is better than using BASEFILE, __LINE__ which does not tell
  where the call from eventually from ans is not useful.
  */

#include "misc.h"
#undef malloc
#undef calloc
#undef free
#undef realloc
void *(*CALLOC)(size_t nmemb, size_t size)=0;
void *(*MALLOC)(size_t size)=0;
void *(*REALLOC)(void*p0, size_t size)=0;
void  (*FREE)(void *p)=0;
void *(*calloc_default)(size_t, size_t)=calloc;
void *(*malloc_default)(size_t)=malloc;
void *(*realloc_default)(void *, size_t)=realloc;
void  (*free_default)(void *)=free;

static int MEM_VERBOSE=0;
static int MEM_DEBUG=0;
PNEW(mutex_mem);
static void *MROOT=NULL;
static long long memcnt=0;
static void *MSTATROOT=NULL;
/*max depth in backtrace */
#define DT 16
typedef struct{
    void *p;
    void *func[DT];
    int nfunc;
    size_t size;
}T_MEMKEY;
typedef struct{
    void *p;
    void **func;
    int nfunc;
    size_t size;
}T_STATKEY;
static int stat_cmp(const void *a, const void *b){
    const T_STATKEY *pa=a;
    const T_STATKEY *pb=b;
    int nfunc=MIN(pa->nfunc, pb->nfunc);
    for(int i=0; i<nfunc; i++){
	if(pa->func[i]<pb->func[i]){
	    return -1;
	}else if(pa->func[i]>pb->func[i]){
	    return 1;
	}
    }
    if(pa->nfunc<pb->nfunc){
	return -1;
    }else if(pa->nfunc>pb->nfunc){
	return 1;
    }else{
	return 0;
    }
}
static void stat_usage(const void *key, VISIT which, int level){
    const T_MEMKEY* key2=*((const T_MEMKEY**)key);
    (void) level;
    if(which==leaf || which==postorder){
	/*printf("%p: size %8zu B", key2->p, (key2->size));
	  print_backtrace_symbol(key2->func, key2->nfunc-2);*/
	void **found;
	T_STATKEY key3;
	key3.func=(void**)key2->func;
	key3.nfunc=key2->nfunc-2;
	if(!(found=tfind(&key3, &MSTATROOT, stat_cmp))){
	    T_STATKEY *keynew=calloc_default(1, sizeof(T_STATKEY));
	    keynew->p=key2->p;
	    keynew->func=malloc_default(key3.nfunc*sizeof(void*));
	    memcpy(keynew->func, key3.func, sizeof(void*)*key3.nfunc);
	    keynew->nfunc=key3.nfunc;
	    keynew->size=key2->size;
	    if(!tsearch(keynew, &MSTATROOT, stat_cmp)){
		error("Error inserting to tree\n");
	    }
	}else{
	    T_STATKEY *keynew=*found;
	    keynew->size+=key2->size;
	}
    }
}
static void print_usage(const void *key, VISIT which, int level){
    const T_STATKEY *key2=*((const T_STATKEY**)key);
    (void) level;
    if(which==leaf || which==postorder){
	info2("size %4zu B@%p", (key2->size), key2->p);
	print_backtrace_symbol(key2->func, key2->nfunc);
    }
}
typedef struct T_DEINIT{/*contains either fun or data that need to be freed. */
    void (*fun)(void);
    void *data;
    struct T_DEINIT *next;
}T_DEINIT;
T_DEINIT *DEINIT=NULL;

static int key_cmp(const void *a, const void *b){
    void *p1,*p2;
    if(!a || !b) return 1; 
    p1=((T_MEMKEY*)a)->p;
    p2=((T_MEMKEY*)b)->p;
    if(p1<p2) return -1;
    else if(p1>p2) return 1;
    else return 0;
}
static void memkey_add(void *p,size_t size){
    if(!p){
	if(size){
	    error("memory allocation for %zu failed\n", size);
	}
	return;
    }
    if(MEM_VERBOSE){
	info("%p malloced with %zu bytes\n",p, size);
    }
    LOCK(mutex_mem);
    T_MEMKEY *key=calloc_default(1,sizeof(T_MEMKEY));
    key->p=p;
    key->size=size;
    key->nfunc=backtrace(key->func,DT);
    if(!tsearch(key, &MROOT, key_cmp)){
	error("Error inserting to tree\n");
    }
    memcnt++;
    UNLOCK(mutex_mem);
}
/*static void memkey_update(void*p, size_t size){
    if(!p) return;
    LOCK(mutex_mem);
    T_MEMKEY key;
    key.p=p;
    void **found;
    if(!(found=tfind(&key, &MROOT, key_cmp))){
	error("Record not found for memkey_update %p\n",p);
    }
    T_MEMKEY*key0=*found;
    key0->size=size;
    key0->nfunc=backtrace(key0->func,DT);
    meminfo("%p realloced with %lu bytes\n", p, size);
    UNLOCK(mutex_mem);
}*/
static void memkey_del(void*p){
    if(!p) return;
    LOCK(mutex_mem);
    void **found=0;
    T_MEMKEY key;
    key.p=p;
    if((found=tfind(&key, &MROOT, key_cmp))){/*found. */
	T_MEMKEY* key1=*found;/*the address of allocated T_MEMKEY. */
	if(!tdelete(&key, &MROOT, key_cmp)){/*return parent. */
	    error("Error deleting old record\n");
	}
	free(key1);
	memcnt--;
    }
    UNLOCK(mutex_mem);
}
static void *calloc_dbg(size_t nmemb, size_t size){
    void *p=calloc_default(nmemb,size);
    memkey_add(p,size);
    return p;
}
static void *malloc_dbg(size_t size){
    void *p=malloc_default(size);
    memkey_add(p,size);
    return p;
}
static void *realloc_dbg(void*p0, size_t size){
    if(!p0) return malloc_dbg(size);
    memkey_del(p0);
    void *p=realloc_default(p0,size);
    memkey_add(p,size);
    return p;
}
static void free_dbg(void *p){
    if(!p) return;
    memkey_del(p);
    free_default(p);
}

/**
   Register a function or data to call or free upon exit
*/
void register_deinit(void (*fun)(void), void *data){
    T_DEINIT *node=calloc(1, sizeof(T_DEINIT));
    node->fun=fun;
    node->data=data;
    node->next=DEINIT;
    DEINIT=node;
}
static __attribute__((constructor)) void init(){
    READ_ENV_INT(MEM_DEBUG, 0, 1);
    READ_ENV_INT(MEM_VERBOSE, 0, 1);
    if(!CALLOC){
	if(MEM_DEBUG){
	    CALLOC=calloc_dbg;
	    MALLOC=malloc_dbg;
	    REALLOC=realloc_dbg;
	    FREE=free_dbg;
	}else{
	    CALLOC=calloc_default;
	    MALLOC=malloc_default;
	    REALLOC=realloc_default;
	    FREE=free_default;
	}
    }
    void init_process(void);
    void init_scheduler(void);
    init_process();
    init_scheduler();
    if(not_nan(NAN)){
	error("NAN check failed\n");
    }
}
/**
   Register routines to be called with mem.c is unloading (deinit).
 */
static __attribute__((destructor)) void deinit(){
    void freepath();
    void thread_pool_destroy();
    freepath();
    thread_pool_destroy();
    for(T_DEINIT *p1=DEINIT;p1;p1=DEINIT){
	DEINIT=p1->next;
	if(p1->fun) p1->fun();
	if(p1->data) FREE(p1->data);
	FREE(p1);
    }
    if(MALLOC==malloc_dbg){
	if(exit_success){
	    if(MROOT){
		warning("%lld allocated memory not freed!!!\n",memcnt);
		twalk(MROOT,stat_usage);
		twalk(MSTATROOT, print_usage);
	    }else{
		info("All allocated memory are freed.\n");
		if(memcnt>0){
		    warning("But memory count is still none zero: %lld\n",memcnt);
		}
	    }
	}else{
	    info("exit_success=%d\n", exit_success);
	}
    }
}
