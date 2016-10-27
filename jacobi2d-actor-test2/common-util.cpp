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
		}*/



	for (int ti=0; ti<N/S1+1; ti++)
		for (int tj=1; tj<N/S2+1; tj++) {
			F* offset = &B_P2_acc(da->B2, 0, ti, tj+1, -HALO_0, 0, 0);
			int x = 0;
			for (int i=ti*S1; i<(ti+1)*S1; i++)
				for (int j=tj*S2+S2-HALO_2; j<(tj+1)*S2; j++) {
					offset[x] = B_P0_accFull(da->B0, 0, i, j);
					x++;
				}
		}

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

#ifdef DUMP_INPUT
	printf("input (cmp)");
	dump_input(da);
#endif
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


