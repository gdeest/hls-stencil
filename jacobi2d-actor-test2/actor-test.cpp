#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "jacobi2d.h"

#include "common-util.cpp"

using namespace std;

/***

load_B_XX is reads from the main memory to the FIFO for each projection
(B is because the variable being updated was named B; you may have multiple arrays updated in a stencil)

The current code is for unroll factor = 2. You can generate the same for other
unroll factors using the perl script. However, for testing with this code,
calls to FIFO (std:queue) must be manually updated, to repalce write by push,
and read by front + pop.

For some of the faces, you need to store the packs in a temp. storage to feed
the fifo in lex order.

***/

void load_B_Pt(F *base, fifo_t &inputs) {
	int t, i, j;
	int count = 0;

	pack_t input[HALO_0];
	for (i=0; i<S1+HALO_1; i++) {
		for (j=HALO_2; j<S2+HALO_2; j++) {
			for (t=0; t<HALO_0; t++) {
				int index = (j + 2 - HALO_2) % 2;
				input[t-0].vals[index] = base[count];
				count++;
				if (index == 1) inputs.push(input[t]); 
			}
		}
	}

}
void load_B_Pi(F *base, fifo_t &inputs) {
	int t, i, j;
	int count = 0;
	

	int jSize = (S2+HALO_2+2-1) / 2; //ceil
	pack_t input[S0][HALO_1][jSize];
	int jCount = 0;
	for (j=0; j<S2+HALO_2; j++) {
		int index = (j + 2 - HALO_2) % 2;
		for (t=HALO_0; t<S0+HALO_0; t++) {
			for (i=0; i<HALO_1; i++) {
				input[t-HALO_0][i-0][jCount].vals[index] = base[count];
				count++;
			}
		}
		if (index == 1) jCount++;
	}


	//must load the entire face to populate the fifo in lex order

	for (t=HALO_0; t<S0+HALO_0; t++) {
		for (i=0; i<HALO_1; i++) {
			for (j=0; j<jSize; j++) {
				inputs.push(input[t-HALO_0][i-0][j]); 
			}
		}
	}

}
void load_B_Pj(F *base, fifo_t &inputs) {
	int t, i, j;
	int count = 0;
	

	pack_t input;
	for (t=0; t<S0+HALO_0; t++) {
		for (i=HALO_1; i<S1+HALO_1; i++) {
			for (j=0; j<HALO_2; j++) {
				int index = (j + 2 - HALO_2) % 2;
				input.vals[index] = base[count];
				count++;
				if (index == 1) inputs.push(input); 
			}
		}
	}

}
void load_B_aux(F *base, fifo_t &inputs) {
	int t, i, j;
	int count = 0;
	
	pack_t input;
	for (t=0; t<HALO_0; t++) {
		for (i=0; i<HALO_1; i++) {
			for (j=0; j<HALO_2; j++) {
				int index = (j + 2 - HALO_2) % 2;
				input.vals[index] = base[count];
				count++;
				if (index == 1) inputs.push(input); 
			}
		}
	}

}

/**
 This function is mostly a wrapper for the above functions. The only thing it
does is to compute the base address for each face.  The base address is the
tile to compute with a negative offset in the corresponding halo dimension to
read the halo regions.
**/
void load_face_fifos(data_array *da, int t0, int t1, int t2, fifo_t &inputPt, fifo_t &inputPi, fifo_t &inputPj, fifo_t &inputAux) {
	F *offset;

	//load
	//P0
	offset = &B_P0_acc(da->B0, t0, t1, t2, 0, -HALO_1, 0);
	load_B_Pt(offset, inputPt);

	//P1
	offset = &B_P1_acc(da->B1, t0, t1, t2, 0, 0, -HALO_2);
	load_B_Pi(offset, inputPi);

	//P2
	offset = &B_P2_acc(da->B2, t0, t1, t2, -HALO_0, 0, 0);
	load_B_Pj(offset, inputPj);

	//aux
	offset = &B_AUX_acc(da->Baux, t0, t1, t2, 0);
	load_B_aux(offset, inputAux);
}


/**
 This function performs the juggling of the FIFOs for each face to form a single FIFO.
 It is functionally correct, but probably not adequate for HLS yet.
**/
void load_face_fifo(data_array *da, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], fifo_t &inputs, int t0, int t1, int t2) {

	int t, i, j;

	fifo_t inputPt, inputPi, inputPj, inputAux;
	load_face_fifos(da, t0, t1, t2, inputPt, inputPi, inputPj, inputAux);

	int offset = (HALO_2>2)?(HALO_2 - 2):(2 - HALO_2);


	for (t=0; t<S0+HALO_0; t++)
		for (i=0; i<S1+HALO_1; i++)
			for (j=offset; j<S2+HALO_2; j+=2) {
				if (t == 0) {
					if (i < HALO_1 && j < HALO_2) {
						inputs.push(inputAux.front());
						inputAux.pop();
					}
					if (i >= HALO_2 && j < HALO_2) {
						inputs.push(inputPj.front());
						inputPj.pop();
					}
					if (j >= HALO_2) {
						inputs.push(inputPt.front());
						inputPt.pop();
					}
				} else {
					if (i < HALO_1) {
						inputs.push(inputPi.front());
						inputPi.pop();
					}
					if (i >= HALO_2 && j < HALO_2) {
						inputs.push(inputPj.front());
						inputPj.pop();
					}
				}

			}
}

/**
 This is a version that uses the face FIFOs to populate the SA registers. This
is used to initialize the tile inputs for partial tiles, which does not use
FIFO for computing.
**/
void load_face_regs(data_array *da, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], int t0, int t1, int t2) {
	int p0, p1, p2;

	fifo_t inputPt, inputPi, inputPj, inputAux;
	load_face_fifos(da, t0, t1, t2, inputPt, inputPi, inputPj, inputAux);

	//P0
	for (p0=0; p0<HALO_0; p0++)
		for (p1=0; p1<S1+HALO_1; p1++)
			for (p2=HALO_2; p2<S2+HALO_2; p2+=2) {
				pack_t pack = inputPt.front(); inputPt.pop();
				reg[p0][p1][p2] = pack.vals[0];
				reg[p0][p1][p2+1] = pack.vals[1];
			}

	//P1
	for (p0=HALO_0; p0<S0+HALO_0; p0++)
		for (p1=0; p1<HALO_1; p1++)
			for (p2=0; p2<S2+HALO_2; p2+=2) {
				pack_t pack = inputPi.front(); inputPi.pop();
				reg[p0][p1][p2] = pack.vals[0];
				reg[p0][p1][p2+1] = pack.vals[1];
			}

	//P2
	for (p0=0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2+=2) {
				pack_t pack = inputPj.front(); inputPj.pop();
				reg[p0][p1][p2] = pack.vals[0];
				reg[p0][p1][p2+1] = pack.vals[1];
			}

	//aux
	for (p0=0; p0<HALO_0; p0++)
		for (p1=0; p1<HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2+=2) {
				pack_t pack = inputAux.front(); inputAux.pop();
				reg[p0][p1][p2] = pack.vals[0];
				reg[p0][p1][p2+1] = pack.vals[1];
			}
}

/**
 This function is responsible for writing the faces in the SA register array to main memory.
Currently, the writes do not use FIFOs.
**/
void store_face_regs(data_array *da, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], int t0, int t1, int t2) {
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

void store_face_fifo(data_array *da, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], fifo_t &outputPt, fifo_t &outputPi, fifo_t &outputPj, fifo_t &outputAux, int t0, int t1, int t2) {
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
				x++;
			}
	//P1
	offset = &B_P1_acc(da->B1, t0, t1+1, t2, 0, 0, 0);
	x = 0;
	for (p2=HALO_2; p2<S2+HALO_2; p2++)
		for (p0=HALO_0; p0<S0+HALO_0; p0++)
			for (p1=S1; p1<S1+HALO_1; p1++) {
				offset[x] = reg[p0][p1][p2];
				x++;
			}

	//P2
	offset = &B_P2_acc(da->B2, t0, t1, t2+1, 0, 0, 0);
	x = 0;
	for (p0=HALO_0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=S2; p2<S2+HALO_2; p2+=2) {

				pack_t pack = outputPj.front();
				outputPj.pop();
				offset[x] = pack.vals[0];
				offset[x+1] = pack.vals[1];
				x+=2;

				//offset[x] = reg[p0][p1][p2];
				//x++;
			}

	//aux
	offset = &B_AUX_acc(da->Baux, t0, t1, t2, 0);
	x = 0;
	for (p0=0; p0<HALO_0; p0++)
		for (p1=0; p1<HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2+=2) {
				pack_t pack = outputAux.front();
				outputAux.pop();
				offset[x] = pack.vals[0];
				offset[x+1] = pack.vals[1];
				x+=2;
			}


}



pack_t compute_pack(pack_t win[3][2]) {
        pack_t result;

        F d[8];

        d[0] = win[2][1].vals[0];
        d[1] = win[2][0].vals[1];
        d[2] = win[1][1].vals[1];
        d[3] = win[1][1].vals[0];
        d[4] = win[1][0].vals[1];
        d[5] = win[1][0].vals[0];
        d[6] = win[0][1].vals[0];
        d[7] = win[0][0].vals[1];

        result.vals[0] = compute_point(d[1], d[3], d[4], d[5], d[7]);
        result.vals[1] = compute_point(d[0], d[2], d[3], d[4], d[6]);

        return result;
}

/**
This code was copied from the HLS version. The tile sizes were renamed to match my convention.
The output values are also written to the SA register array for correctly performing the writes.
**/
void compute_full_tile(data_array *da, fifo_t &inputs, fifo_t &outputs, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2]) {


#define BUFF(tt,xx0,xx1) buff[(tt)%2][xx0][xx1]

        pack_t buff[2][S1][S2/2];

 	pack_t lines[2][(S2/2)+1];

	pack_t win[3][2];

        for (int tt=0; tt<S0; tt++) {
    		for (int xx0=-2; xx0<S1; xx0++) {
			for (int xx1=-1; xx1<S2/2; xx1++) {
        bool outside = xx0<0 || xx1<0;
        bool read_fifo = outside || tt==0;
        bool output_result = (tt==S0-1 || xx0>=S1-2 || xx1==S2/2-1);

        pack_t input;

        {   
          // Shift window down
          win[0][0] = win[0][1];
          win[1][0] = win[1][1];
          win[2][0] = win[2][1];

          /* Read top window row */
          // Retrieve top left and middle pack_t from line buffers
          win[0][1] = lines[0][xx1+1];
          win[1][1] = lines[1][xx1+1];

          // Retrieve top right pack_t from either previous timestep or input FIFO.
          if (read_fifo) {
            	input = inputs.front();
		inputs.pop();

		//needed for testing only
  		reg[tt][xx0+HALO_1][xx1*2+HALO_2] = input.vals[0];
  		reg[tt][xx0+HALO_1][xx1*2+HALO_2+1] = input.vals[1];
          } else {
          	input = BUFF(tt-1,xx0,xx1);
          }   
          win[2][1] = input;

          // Shift lines right and store input.
          lines[0][xx1+1] = lines[1][xx1+1];
          lines[1][xx1+1] = input;
        }

        if (!outside) {
          pack_t result = compute_pack(win);

          BUFF(tt, xx0, xx1) = result;

	  reg[tt+HALO_0][xx0+HALO_1][xx1*2+HALO_2] = result.vals[0];
	  reg[tt+HALO_0][xx0+HALO_1][xx1*2+HALO_2+1] = result.vals[1];

          if (output_result)
	    outputs.push(result);
        }
      }
    }
        }
//printf("\n");


	/*
	int p0, p1, p2;
	for (p0=HALO_0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=HALO_2; p2<S2+HALO_2; p2++) {
				reg[p0][p1][p2] = compute_point(reg[p0-1][p1-1][p2], reg[p0-1][p1][p2-1], reg[p0-1][p1-1][p2-1], reg[p0-1][p1-1][p2-2], reg[p0-1][p1-2][p2-1]);
			}*/
}

/**
 This function takes care of the partial tiles. It is inefficient but works for
testing.  When there is no computation, the values are forwarded so that the
boundary values are always available.
*/
void compute_partial_tile(data_array *da, F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], int t0, int t1, int t2) {
	int p0, p1, p2;
	
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

				if (1<=t && t<T && 1<=i && i<N-1 && 1<=j && j<N-1 && HALO_0<=p0 && HALO_1<=p1 && HALO_2<=p2) {
					reg[p0][p1][p2] = compute_point(reg[p0-1][p1-1][p2], reg[p0-1][p1][p2-1], reg[p0-1][p1-1][p2-1], reg[p0-1][p1-1][p2-2], reg[p0-1][p1-2][p2-1]);
				} else {
					reg[p0][p1][p2] = reg[p0-1][p1-1][p2-1];
				}
			}
}

void output_shuffle(fifo_t &outputs, fifo_t &outputPt, fifo_t &outputPi, fifo_t &outputPj, fifo_t &outputAux) {

	int t, i, j;

int count = 0;
	for (t=0; t<S0; t++)
		for (i=0; i<S1; i++)
			for (j=0; j<S2; j+=2) {
				bool i1border = (i  >= S1-HALO_1);
				bool i2border = (j + 2-1 >= S2-HALO_2);
				if (t == S0-1) {
					pack_t output = outputs.front();
					outputs.pop();


					outputPt.push(output);

					if (i1border && i2border) {
						outputAux.push(output);
					}

					if (i1border) {
						outputPi.push(output);
					}
					if (i2border) {
						outputPj.push(output);
					}
				} else {
					if ((i1border) || i2border) {
						pack_t output = outputs.front();
						outputs.pop();

						if (i1border) {
							outputPi.push(output);
						}
						if (i2border) {
							outputPj.push(output);

						}
					}
				}

			}

}

void compute_tile(data_array *da, int t0, int t1, int t2) {
	int p0, p1, p2;
	int x;
	F *offset;

	F reg[S0+HALO_0][S1+HALO_1][S2+HALO_2];
	bool isFullTile;

	fifo_t inputs, outputs;


	//This part executes the same test as in the partial tile code to see if the tile is full or partial. Inefficient but works
	isFullTile = true;
	for (p0=HALO_0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=HALO_2; p2<S2+HALO_2; p2++) {
				//guard for partial tiles
				int st = t0*S0+p0;
				int si = t1*S1+p1-HALO_1;
				int sj = t2*S2+p2-HALO_2;
				int t = st;
				int i = si - st;
				int j = sj - st;
				isFullTile &= (1<=t && t<T && 1<=i && i<N-1 && 1<=j && j<N-1);
			}

	if (isFullTile) {
		load_face_fifo(da, reg, inputs, t0, t1, t2);
		compute_full_tile(da, inputs, outputs, reg);
		fifo_t outputPt, outputPi, outputPj, outputAux;
		output_shuffle(outputs, outputPt, outputPi, outputPj, outputAux);
		store_face_fifo(da, reg, outputPt, outputPi, outputPj, outputAux, t0, t1, t2);
	} else {
		load_face_regs(da, reg, t0, t1, t2);
		compute_partial_tile(da, reg, t0, t1, t2);
		store_face_regs(da, reg, t0, t1, t2);
	}


/*
#ifdef DUMP_REG_IN
	printf("tile (%d %d %d) input \n", t0, t1, t2);
	dump_reg(reg);
#endif
#ifdef DUMP_REG_OUT
	printf("tile (%d %d %d) output \n", t0, t1, t2);
	dump_reg(reg);
#endif
*/
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



