#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <queue>

#include "jacobi2d.h"


using namespace std;

typedef struct {
        F vals[2];
} pack_t;

typedef queue<pack_t> fifo_t;


struct data_array {
	F *B0;
	F *B1;
	F *B2;
	F *Baux;
};

typedef struct data_array data_array;

F* alloc_array_F(int size) {
	printf("Allocating array of size: %d\n", size*sizeof(F));
	F *arr = (F*)ALLOC(size*sizeof(F));
	
	if (!arr) {
		printf("Error allocating array.\n");
		exit(-1);
	}
	
	return arr;
}

void allocate_data_array(data_array *da) {
	da->B0 =  alloc_array_F(B_P0_SIZE);
	da->B1 =  alloc_array_F(B_P1_SIZE);
	da->B2 =  alloc_array_F(B_P2_SIZE);
	da->Baux =  alloc_array_F(B_AUX_SIZE);
}


void dump_input(data_array *da) {
	int t, i, j;
	t = T-1;

	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			printf ("%.02f ", B_P0_accFull(da->B0, t, i, j));
		}
		printf("\n");
	}
}

void dump_output(data_array *da) {
	int t, i, j;
	t = T-1;

	printf("out(cmp)\n");

	for (i=T; i<T+N; i++) {
		for (j=T; j<T+N; j++) {
			printf ("%.02f ", B_P0_accFull(da->B0, t+S0, i, j));
		}
		printf("\n");
	}
}

void init_data_array(data_array *da, unsigned int seed) {
	int i,j;
	srand(seed);
	F count = 1;

	//P0 
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			B_P0_accFull(da->B0, 0, i, j) = count;//(1.0*rand()/RAND_MAX)*300.0;
			count++;
		}
	}

	//P2 : might have some issue due to % of -1
	/*for (int i=0; i<N; i++)
		for (int j=0; j<N; j++) {
			B_P2_accFull(da->B2, -1, i, j+S2) = B_P0_accFull(da->B0, 0, i, j+S2-HALO_2);
		}

	for (int ti=0; ti<N/S1+1; ti++)
		for (int tj=1; tj<N/S2+1; tj++) {
			int x =0;
			for (int i=ti*S1; i<(ti+1)*S1; i++)
				for (int j=tj*S2+S2-HALO_2; j<(tj+1)*S2; j++) {
	B_P2_accFull(da->B2, 0, i, j) = B_P0_accFull(da->B0, 0, i, j);
				x++;
			}
		}
*/
	//aux
	for (int ti=1; ti<N/S1+1; ti++)
		for (int tj=1; tj<N/S2+1; tj++) {
			int x =0;
			for (int i=(ti-1)*S1+S1-HALO_1; i<(ti)*S1; i++)
				for (int j=(tj-1)*S2+S2-HALO_2; j<(tj)*S2; j++) {
	B_AUX_acc(da->Baux, 0, ti, tj, x) = B_P0_accFull(da->B0, 0, i, j);
					x++;
				}
	}

/*
	for (int i=0; i < B_P2_SIZE; i++) {
		printf("%.02f ", da->B2[i]);
	}
	printf("\n");
	for (int i=0; i < B_P0_SIZE; i++) {
		printf("%.02f ", da->B0[i]);
	}
	printf("\n");

	for (int i=0; i < B_P1_SIZE; i++) {
		printf("%.02f ", da->B1[i]);
	}
	printf("\n");

	for (int i=0; i < B_P2_SIZE; i++) {
		printf("%.02f ", da->B2[i]);
	}
	printf("\n");
*/

// #ifdef DUMP_INPUT
	printf("input (cmp)");
	dump_input(da);
// #endif
}

void golden(F B[2][N][N], unsigned int seed) {
	int t,i,j;
	srand(seed);
	int count = 1;

	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			B[0][i][j] = count;//(1.0*rand()/RAND_MAX)*300.0;
			B[1][i][j] = count;//(1.0*rand()/RAND_MAX)*300.0;
			count++;
		}
	}


	for (int t=0; t<T; t++) {
		for (int i=1; i<N-1; i++) {
			for (int j=1; j<N-1; j++) {
				B[(t+1)%2][i][j] = (B[t%2][i][j] + B[t%2][i+1][j] + B[t%2][i-1][j] + B[t%2][i][j+1] + B[t%2][i][j-1])*0.2;
			}
		}
	}

	
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			fprintf(stderr, "%.2lf ", B[1][i][j]);
			//if (j%20==0) fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}
}


F compute_point(F a, F b, F c, F d, F e) {

	F res =  (a + b + c + d + e)*0.2;
	//printf("%.02f = %.02f %.02f %.02f %.02f %.02f\n", res, a, b, c, d, e);
	return res;
}

void dump_reg(F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2]) {

	for (int i0 = 0; i0<S0+HALO_0; i0++)  {
		for (int i1 = 0; i1<S1+HALO_1; i1++) {
			for (int i2 = 0; i2<S2+HALO_2; i2++) {
				printf("%.02f ", reg[i0][i1][i2]);
			}
			printf("\n");
		}
		printf("\n");
	}

}
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
#ifdef DUMP_OUTPUT
	dump_output(&data);
#endif
#ifdef GOLDEN
	F* B = (F*) alloc_array_F(2*N*N);
	golden((F[2][N][N])B, SEED);
#endif
	
	return 0;	
}



