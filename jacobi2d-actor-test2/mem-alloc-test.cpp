#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "jacobi2d.h"

#include "common-util.cpp"

using namespace std;

void load_face(data_array *da, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], int t0, int t1, int t2) {
	int p0, p1, p2;
	int x;
	F *offset;

	//load
	//P0
	offset = &B_P0_acc(da->B0, t0, t1, t2, 0, -HALO_1, 0);
	x = 0;
	for (p1=0; p1<S1+HALO_1; p1++)
		for (p2=HALO_2; p2<S2+HALO_2; p2++)
			for (p0=0; p0<HALO_0; p0++) {
				reg[p0][p1][p2] = offset[x];
				//printf("regP0[%d][%d][%d] = %.02f\n", p0, p1, p2, offset[x]);
				x++;
			}
	//P1
	offset = &B_P1_acc(da->B1, t0, t1, t2, 0, 0, -HALO_2);
//printf("offsetP1 = %d\n", B_P1_offset(t0, t1, t2, 0, 0, -HALO_2));
	x = 0;
	for (p2=0; p2<S2+HALO_2; p2++)
		for (p0=HALO_0; p0<S0+HALO_0; p0++)
			for (p1=0; p1<HALO_1; p1++) {
				reg[p0][p1][p2] = offset[x];
				//printf("regP1[%d][%d][%d] = %.02f\n", p0, p1, p2, offset[x]);
				x++;
			}

	//P2
	offset = &B_P2_acc(da->B2, t0, t1, t2, -HALO_0, 0, 0);
	x = 0;
	for (p0=0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2++) {
				reg[p0][p1][p2] = offset[x];
				//printf("regP2[%d][%d][%d] = %.02f\n", p0, p1, p2, offset[x]);
				x++;
			}

	//aux
	offset = &B_AUX_acc(da->Baux, t0, t1, t2, 0);
	x = 0;
	for (p0=0; p0<HALO_0; p0++)
		for (p1=0; p1<HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2++) {
				reg[p0][p1][p2] = offset[x];
				//printf("regAUX[%d][%d][%d] = %.02f @ %d\n", p0, p1, p2, offset[x], B_AUX_offset(t0, t1, t2, x));
				x++;
			}
}

void store_face(data_array *da, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], int t0, int t1, int t2) {
	int p0, p1, p2, x;
	F *offset;

	//store
	//P0
	offset = &B_P0_acc(da->B0, t0+1, t1, t2, 0, 0, 0);
	x = 0;
	for (p1=HALO_1; p1<S1+HALO_1; p1++)
		for (p2=HALO_2; p2<S2+HALO_2; p2++)
			for (p0=S0; p0<S0+HALO_0; p0++) {
				offset[x] = reg[p0][p1][p2];
//printf("offsetP0[%d] = %.02f @ %d\n", x, offset[x], B_P0_offset(t0+1,t1,t2,0,0,0)+ x);
				x++;
			}
	//P1
	offset = &B_P1_acc(da->B1, t0, t1+1, t2, 0, 0, -HALO_2);
//printf("offsetP1 = %d\n", B_P1_offset(t0, t1+1, t2, 0, 0, -HALO_2));
	x = 0;
	for (p2=0; p2<S2+HALO_2; p2++)
		for (p0=HALO_0; p0<S0+HALO_0; p0++)
			for (p1=S1; p1<S1+HALO_1; p1++) {
				offset[x] = reg[p0][p1][p2];
//printf("offsetP1[%d] = %.02f @ %d\n", x, offset[x], B_P1_offset(t0, t1+1, t2, 0, 0, -HALO_2)+x);
				x++;
			}

	//P2
	offset = &B_P2_acc(da->B2, t0, t1, t2+1, -HALO_0, 0, 0);
	x = 0;
	for (p0=0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=S2; p2<S2+HALO_2; p2++) {
				offset[x] = reg[p0][p1][p2];
//printf("offsetP2[%d] = %.02f\n", x, offset[x]);
				x++;
			}

	//aux
	offset = &B_AUX_acc(da->Baux, t0, t1, t2, 0);
	x = 0;
	for (p0=0; p0<HALO_0; p0++)
		for (p1=0; p1<HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2++) {
				offset[x] = reg[p0+S0][p1+S1][p2+S2];
//printf("offsetAUX[%d] = %.02f @ %d\n", x, offset[x], B_AUX_offset(t0, t1, t2, x));
				x++;
			}


}

void compute_tile(data_array *da, int t0, int t1, int t2) {
	int p0, p1, p2;
	int x;
	F *offset;

	F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2];

	//CMP for intra-tile; useful only for handling partial tiles
	//F out0[HALO_0][S1+HALO_1][S2];
	//F out1[S0][HALO_1][S2+HALO_2];
	//F out2[S0+HALO_0][S1][HALO_2];

	load_face(da, reg, t0, t1, t2);

#ifdef DUMP_REG_IN
	printf("tile (%d %d %d) input \n", t0, t1, t2);
	dump_reg(reg);
#endif

	//for (p0=HALO_0; p0<S0+HALO_0; p0++)
	//	for (p1=HALO_1; p1<S1+HALO_1; p1++)
	//		for (p2=HALO_2; p2<S2+HALO_2; p2++) {
	for (p0=1; p0<S0+HALO_0; p0++)
		for (p1=1; p1<S1+HALO_1; p1++)
			for (p2=1; p2<S2+HALO_2; p2++) {
				//guard for partial tiles
				int st = t0*S0+p0;
				int si = t1*S1+p1-HALO_1;
				int sj = t2*S2+p2-HALO_2;
				int t = st;
				int i = si - st;
				int j = sj - st;

				if (1<=t && t<T && 1<=i && i<N-1 && 1<=j && j<N-1 && HALO_0<p0 && HALO_1<p1 && HALO_2<p2) {
//printf("%d %d %d -> %d %d %d -> %d %d %d\n", p0, p1, p2, st,si,sj,t,i,j);
					reg[p0][p1][p2] = compute_point(reg[p0-1][p1-1][p2], reg[p0-1][p1][p2-1], reg[p0-1][p1-1][p2-1], reg[p0-1][p1-1][p2-2], reg[p0-1][p1-2][p2-1]);
				} else {
					reg[p0][p1][p2] = reg[p0-1][p1-1][p2-1];
				}
			}

#ifdef DUMP_REG_OUT
	printf("tile (%d %d %d) output \n", t0, t1, t2);
	dump_reg(reg);
#endif

	store_face(da, reg, t0, t1, t2);
}

void compute(data_array *da) {
	int t0, t1, t2;
	for (t0=0; t0<(N0+S0-1)/S0; t0++)
		for (t1=0; t1<(N1+S1-1)/S1; t1++)
			for (t2=0; t2<(N2+S2-1)/S2; t2++) {
				compute_tile(da, t0, t1, t2);
			}
	
}


int main() {

	data_array data;
	allocate_data_array(&data);
	init_data_array(&data, SEED);
	compute(&data);
// #ifdef DUMP_OUTPUT
	dump_output(&data);
// #endif
#ifdef GOLDEN
	F* B = (F*) alloc_array_F(2*N*N);
	golden((F[2][N][N])B, SEED);
#endif
	
	return 0;	
}



