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

/*Taken from suitesparse package from UFL. The operations performs ok with
 * non-ascending ordering in each column and may produce results such way. Sort
 * the columns in ascending order before exporting to matlab.*/


#include <limits.h>
#include "../sys/sys.h"
#include "type.h"
/*Obtain from SuiteSparse ss_multiply. Gather together */
typedef spint SS_INT;
#define SS_MAX(a,b) (((a) > (b)) ? (a) : (b))

#define SS_ENTRY T
#define cs X(sp)

#define SS_FLIP(i) (-(i)-2)
#define SS_UNFLIP(i) (((i) < 0) ? SS_FLIP(i) : (i))
#define SS_MARKED(w,j) (w [j] < 0)
#define SS_MARK(w,j) { w [j] = SS_FLIP (w [j]) ; }
#define SS_CSC(A) (A)

/**
 free a sparse matrix */
static cs *ss_spfree (cs *A)
{
    delete A;
    /*if(A){
	free (A->p) ;
	free (A->i) ;
	free (A->x) ;
	}*/
    return NULL;
}
/**
 wrapper for malloc */
static void *ss_malloc (SS_INT n, size_t size)
{
    return (malloc (SS_MAX (n,1) * size)) ;
}

/**
 wrapper for calloc */
static void *ss_calloc (SS_INT n, size_t size)
{
    return (calloc (SS_MAX (n,1), size)) ;
}


/**
 change the max # of entries sparse matrix */
/*static SS_INT ss_sprealloc (cs *A, SS_INT nzmax)
{
    SS_INT ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (SS_CSC (A)) ? (A->p [A->ny]) : A->nz ;
    A->i = (SS_INT*)ss_realloc (A->i, nzmax, sizeof (SS_INT), &oki) ;
    if (A->x) A->x = (SS_ENTRY*)ss_realloc (A->x, nzmax, sizeof (SS_ENTRY), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}*/

/**
 free workspace and return a sparse matrix result */
static cs *ss_done (cs *C, void *w, void *x, SS_INT ok)
{
    free (w) ;                       /* free workspace */
    free (x) ;
    return (ok ? C : ss_spfree (C)) ;   /* return result if OK, else free it */
}

/**
 x = x + beta * A(:,j), where x is a dense vector and A(:,j) is dsp */
static SS_INT ss_scatter (const cs *A, SS_INT j, SS_ENTRY beta, SS_INT *w, SS_ENTRY *x, SS_INT mark,
    cs *C, SS_INT nz)
{
    SS_INT i, p, *Ci;
    const  SS_INT *Ap, *Ai ;
    const SS_ENTRY *Ax ;
    if (!SS_CSC (A) || !w || !SS_CSC (C)) return (-1) ;     /* check inputs */
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}

/**
 C = A*B */
cs* X(ss_multiply) (const cs *A, const cs *B)
{
    SS_INT p, j, nz = 0, anz, *Cp, *Ci, m, n, bnz, *w;
    const SS_INT *Bp, *Bi;
    SS_ENTRY *x, *Cx;
    const SS_ENTRY *Bx ;
    cs *C ;
    if (!SS_CSC (A) || !SS_CSC (B)) {
	error("A or B is not in CSC format\n");
	return (NULL) ;      /* check inputs */
    }
    if (A->ny != B->nx) {
	error("Matrix mismatch\n");
	return (NULL) ;
    }
    m = A->nx ; anz = A->p [A->ny] ;
    n = B->ny ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
    w = (SS_INT*)ss_calloc (m, sizeof (SS_INT)) ;                    /* get workspace */
    x = (SS_ENTRY*)ss_malloc (m, sizeof (SS_ENTRY)) ; /* get workspace */
    //C = ss_spalloc (m, n, anz + bnz, values) ;        /* allocate result */
    C = X(spnew) (m ,n, anz+bnz);
    if (!C || !w || (!x)){
	error("Out of memory\n");
	return (ss_done (C, w, x, 0)) ;
    }
    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !X(spsetnzmax) (C, 2*(C->nzmax)+m))
        {
	    error("Out of memory. \n");
            return (ss_done (C, w, x, 0)) ;             /* out of memory */
        } 
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = ss_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    X(spsetnzmax) (C, 0) ;               /* remove extra space from C */
    return (ss_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/**
 C = alpha*A + beta*B */
cs* X(ss_add) (const cs *A, const cs *B, SS_ENTRY alpha, SS_ENTRY beta)
{
    SS_INT p, j, nz = 0, anz, *Cp, *Ci, m, n, bnz, *w ;
    const SS_INT *Bp;
    SS_ENTRY *x, *Cx ;
    cs *C ;
    if (!SS_CSC (A) || !SS_CSC (B)){
	error("A or B is not in CSC format\n");
	return (NULL) ;         /* check inputs */
    }
    if (A->nx != B->nx || A->ny != B->ny) {
	error("Matrix mismatch\n");
	return (NULL) ;
    }
    m = A->nx ; anz = A->p [A->ny] ;
    n = B->ny ; Bp = B->p ; ; bnz = Bp [n] ;
    w = (SS_INT*)ss_calloc (m, sizeof (SS_INT)) ;                       /* get workspace */
    x = (SS_ENTRY*)ss_malloc (m, sizeof (SS_ENTRY)) ;    /* get workspace */
    SS_INT cnz=anz + bnz;
    if(cnz>m*n){
	cnz=m*n;
    }
    C = X(spnew)(m, n, cnz) ;           /* allocate result*/
    if (!C || !w || (!x)) {
	error("Out of memory\n");
	return (ss_done (C, w, x, 0)) ;
    }
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)
    {
        Cp [j] = nz ;                   /* column j of C starts here */
        nz = ss_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        nz = ss_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    X(spsetnzmax) (C, 0) ;               /* remove extra space from C */
    return (ss_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/**
 drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
static SS_INT ss_fkeep (cs *A, SS_INT (*fkeep) (SS_INT, SS_INT, SS_ENTRY, void *), void *other)
{
    SS_INT j, p, nz = 0, n, *Ap, *Ai ;
    SS_ENTRY *Ax ;
    if (!SS_CSC (A) || !fkeep) return (-1) ;    /* check inputs */
    n = A->ny ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        Ap [j] = nz ;                       /* record new location of col j */
        for ( ; p < Ap [j+1] ; p++)
        {
            if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other))
            {
                if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
                Ai [nz++] = Ai [p] ;
            }
        }
    }
    Ap [n] = nz ;                           /* finalize A */
    X(spsetnzmax) (A, 0) ;                   /* remove extra space from A */
    return (nz) ;
}
/**
   return true of is not zero.
*/
static SS_INT ss_nonzero (SS_INT i, SS_INT j, SS_ENTRY aij, void *other)
{
    (void) i;
    (void) j;
    (void) other;
    return fabs(aij)>1e-50 ;
}
/**
   drop zeros in the sparse matrix.
 */
SS_INT X(ss_dropzeros) (cs *A)
{
    return (ss_fkeep (A, &ss_nonzero, NULL)) ;  /* keep all nonzero entries */
}
/**
   whether value is below threshold.
 */
static SS_INT ss_tol (SS_INT i, SS_INT j, SS_ENTRY aij, void *tol)
{
    (void) i;
    (void) j;
    return (fabs (aij) > *((double *) tol)) ;
}
/**
   drop values below threashold of tol.
*/ 
SS_INT X(ss_droptol) (cs *A, double tol)
{
    return (ss_fkeep (A, &ss_tol, &tol)) ;    /* keep all large entries */
}
static double ss_cumsum (SS_INT *p, SS_INT *c, SS_INT n)
{
    SS_INT i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid SS_INT overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

cs* X(ss_transpose) (const cs *A, SS_INT values)
{
    SS_INT p, q, j, *Cp, *Ci, n, m, *w ;
    SS_ENTRY *Cx;
    const SS_INT *Ap, *Ai;
    const SS_ENTRY *Ax;
    cs *C ;
    if (!SS_CSC (A)) return (NULL) ;    /* check inputs */
    m = A->nx ; n = A->ny ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = X(spnew) (n, m, Ap [n]) ;       /* allocate result */
    w = (SS_INT*)ss_calloc (m, sizeof (SS_INT)) ;                      /* get workspace */
    if (!C || !w) return (ss_done (C, w, NULL, 0)) ;       /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    ss_cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = (values > 0) ? CONJ (Ax [p]) : Ax [p] ;
        }
    }
    return (ss_done (C, w, NULL, 1)) ;  /* success; free w and return C */
}
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
