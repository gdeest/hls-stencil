#include <queue>

typedef float F;

typedef struct {
        F vals[2];
} pack_t;

typedef std::queue<pack_t> fifo_t;


struct data_array {
	F *B0;
	F *B1;
	F *B2;
	F *Baux;
};

typedef struct data_array data_array;


//#define DUMP_INPUT
/* #define DUMP_OUTPUT */
//#define GOLDEN

//#define DUMP_REG_IN
//#define DUMP_REG_OUT

/* PROTECTED REGION ID(parameters) DISABLED START */
#define T 100
#define N 200
/* PROTECTED REGION END */

/* PROTECTED REGION ID(tilesizes) DISABLED START */
#define S0 8
#define S1 8
#define S2 8
/* PROTECTED REGION END */

#define TILE_TO_FULL(t0, t1, t2, p0, p1, p2) (t0*S0+p0, t1*S1+p1, t2*S2+p2)
#define SKEW(p0, p1, p2) ((p0), (p0) + (p1), (p0) + (p2))
#define SKEW_INV(p0, p1, p2) ((p0), -(p0) + (p1), -(p0) + (p2))

//Skewed problem sizes
#define N0 (T)
#define N1 (T) + (N)
#define N2 (T) + (N)

//Number of tiles
#define NT0 (((T) + S0 - 1) / S0)        // floord(N0, S0)
#define NT1 (((T) + (N) + S1 - 1) / S1)  // floord(N1, S1)
#define NT2 (((T) + (N) + S2 - 1) / S2)  // floord(N2, S2)
#define NT (NT0)*(NT1)*(NT2)

//Halo
#define HALO_0 1
#define HALO_1 2
#define HALO_2 2

//Tile faces

// Projection along dim 0 (t)
// 0 extends 1
#define FACE_P0_0 HALO_0
#define FACE_P0_1 S1+HALO_1
#define FACE_P0_2 S2
#define FACE_P0 (FACE_P0_0)*(FACE_P0_1)*(FACE_P0_2)

// 1 extends 2
#define FACE_P1_0 S0
#define FACE_P1_1 HALO_1
#define FACE_P1_2 S2+HALO_2
#define FACE_P1 (FACE_P1_0)*(FACE_P1_1)*(FACE_P1_2)

// 2 extends 0
#define FACE_P2_0 S0+HALO_0
#define FACE_P2_1 S1
#define FACE_P2_2 HALO_2
#define FACE_P2 (FACE_P2_0)*(FACE_P2_1)*(FACE_P2_2)


#define LIN2D(i0, i1, n0, n1) (i0)*(n1)+(i1)
#define LIN3D(i0, i1, i2, n0, n1, n2) (i0)*(n1*n2)+(i1)*(n2)+(i2)

//Memory Access
#define B_P0_SIZE (HALO_0)*(N1+HALO_1*(NT0+1))*(N2)
#define B_P1_SIZE (N0)*(HALO_1)*(N2+HALO_2*(NT1+1))
#define B_P2_SIZE (N0+HALO_0*(NT2+1))*(N1)*(HALO_2)
#define B_AUX_SIZE (HALO_0*HALO_1*HALO_2)*(NT0+NT1)*(NT0+NT2)


//most generic form for P0 access is the following
//#define B_P0_acc(A, tt, ti, tj, t, i, j) A[(tj)*(HALO_0*(N1+HALO_1)*S2)+(ti)*(HALO_0*S1*S2)+HALO_0*HALO_1*S2+LIN3D(i,j,(t)%HALO_0,S1,S2,HALO_0)]

//

#define B_P0_offset(tt, ti, tj, t, i, j) (tj)*(HALO_0*(N1+HALO_1)*S2)+(ti)*(HALO_0*S1*S2)+(NT0+1-(tt))*HALO_0*HALO_1*S2+LIN2D(i,j,S1,S2)
#define B_P0_acc(A, tt, ti, tj, t, i, j) A[B_P0_offset(tt, ti, tj, t, i, j)]
#define B_P0_accFull(A, t, i, j) B_P0_acc(A, (t)/S0, (i)/S1, (j)/S2, (t)%S0, (i)%S1, (j)%S2)

#define B_P1_offset(tt, ti, tj, t, i, j) (tt)*(S0*HALO_1*(N2+HALO_2))+(tj)*(S0*HALO_1*S2)+(NT1+1-(ti))*S0*HALO_1*HALO_2+LIN3D(j,t,(i)%2,S2,S0,2)
#define B_P1_acc(A, tt, ti, tj, t, i, j) A[B_P1_offset(tt, ti, tj, t, i, j)]
#define B_P1_accFull(A, t, i, j) B_P1_acc(A, (t)/S0, (i)/S1, (j)/S2, (t)%S0, (i)%S1, (j)%S2)

#define B_P2_offset(tt, ti, tj, t, i, j) ((ti)*((N0+HALO_0)*S1*HALO_2)+(tt)*(S0*S1*HALO_2)+(NT2+1-(tj))*HALO_0*S1*HALO_2+LIN3D(t,i,(j)%2,S0,S1,2))
#define B_P2_acc(A, tt, ti, tj, t, i, j) A[B_P2_offset(tt, ti, tj, t, i, j)]
#define B_P2_accFull(A, t, i, j) B_P2_acc(A, (t)/S0, (i)/S1, (j)/S2, (t)%S0, (i)%S1, (j)%S2)



#define B_P2_offsetFull(t, i, j) B_P2_offset((t)/S0, (i)/S1, (j)/S2, (t)%S0, (i)%S1, (j)%S2)

#define B_AUX_offset(tt, ti, tj, x) (LIN3D((tt)-(ti)+NT1, (tt)-(tj)+NT2, (x), (NT0)+(NT1), (NT0)+(NT2), (HALO_0*HALO_1*HALO_2)))
#define B_AUX_acc(A, tt, ti, tj, x) A[B_AUX_offset(tt, ti, tj, x)]


// Select alloc/free functions based on whether we're compiling with SDSoC or not.
#if defined(__SDSCC__) || defined(__SDSVHLS__)
#include "sds_lib.h"
#define ALLOC sds_alloc
#define FREE sds_free
#else
#define ALLOC malloc
#define FREE  free
#endif

#define SEED 4784


