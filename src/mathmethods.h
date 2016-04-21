#ifndef MATHMETHODS_H
#define MATHMETHODS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct 
{
	double* ve;
	unsigned int dim;
} VEC;

typedef struct
{
	int* ve;
	unsigned int dim;
} IVEC;

typedef struct
{
	double* storage;
	double** me;
	unsigned int m,n;
} MAT;

//Initialialization
VEC v_get(int size);
IVEC iv_get(int size);
MAT m_get(int m,int n);

//Resize
void v_resize(VEC *v, int size);

//Release
void v_free(VEC *v);
void iv_free(IVEC *v);
void m_free(MAT *A);

//Copy
void v_copy(VEC src, VEC dest);
void v_copy_n(VEC src, VEC dest, unsigned int n);
void m_copy(MAT src, MAT dest);

void v_add(VEC u,VEC v,VEC w,unsigned int start);
void m_add(MAT A,MAT B,MAT C);
void mv_mlt(MAT A,VEC u,VEC v);
void mTv_mlt(MAT A,VEC u,VEC v);
void mm_mlt(MAT A,MAT B,MAT C);
void mTm_mlt(MAT A,MAT B,MAT C);
void mmT_mlt(MAT A,MAT B,MAT C);
void v_sub(VEC u,VEC v,VEC w,unsigned int start);
void sm_mlt(double alpha,MAT A,MAT B,unsigned int startrow,unsigned int startcol);
void dipaa(double alpha,MAT A,MAT B,unsigned int startrow,unsigned int startcol);
void diag_mm_mlt(MAT A,MAT B,MAT C);
void diag_mm_mltT(MAT A,MAT B,MAT C);
void diag_m_add(MAT A,MAT B,MAT C);
void diag_invert(MAT A,MAT B);
void diag_mv_mlt(MAT A,VEC u,VEC v);
void mTdiag_mlt(MAT A,VEC B,MAT C);
void mdiag_mlt(MAT A,VEC B,MAT C);
void diagm_mlt(VEC A,MAT B,MAT C);

void VDVT_mlt(MAT V,VEC D,MAT A);
void VDinvVT_mlt(MAT V,VEC D,MAT A,VEC temp);
void VsqrtDinvVT_mlt(MAT V,VEC D,MAT A,VEC temp);
void VsqrtDVT_mlt(MAT V,VEC D,MAT A,VEC temp);
double vTAv(MAT A,VEC v);

VEC lagrange_weights(short unsigned int s, VEC c);
//double lagrange_bary(double theta,short unsigned int s,VEC c,VEC* Q,VEC w);
void lagrange_bary(double theta,VEC c,VEC* Z,VEC w,VEC sum);
double lagrange(double theta,short unsigned int s,VEC c,VEC* Q);
double norm_inf(VEC v,VEC w,unsigned int start);
double norm_inf_u(VEC v,VEC w,unsigned int start,unsigned int t);
double vector_norminf(VEC v,unsigned int start);
void Print_Vector(VEC v);
void Print_Matrix(MAT A);
void Print_VectorC(VEC v);
void Print_MatrixC(MAT A);

#if !defined(_MSC_VER)
#define min(val1,val2)	(((val1) > (val2)) ? (val2) : (val1))
#define max(val1,val2)	(((val1) < (val2)) ? (val2) : (val1))
#endif

double sq(double val);
void daxpy(double alpha, VEC x, VEC y, unsigned int begin);
void daxpy_u(double alpha, VEC x, VEC y, unsigned int begin, unsigned int end);
void sv_mlt(double val, VEC v, unsigned int begin);
void sv_mlt_u(double val, VEC v, unsigned int begin, unsigned int end);


#endif //MATHMETHODS_H