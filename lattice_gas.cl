//old opencl in a machine
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//optimal blocks size for GTX630 
#define BSIZE 10


////////////////////////////////////////////////////////////////////////
/*
	blocking update
*/
////////////////////////////////////////////////////////////////////////
__kernel void lattice_gas_0( 
					__global int* vec1,
					__global int* vec2,
					write_only image2d_t output,
					__const unsigned int N
					)
{

	const int gid = get_global_id(0);
	int i,j,ii,jj;

	int maxblock = (N-2)/(BSIZE-2);

	//there are just full blocks now
	//calculated subblock starting pos:
	ii=1 + (gid/maxblock)*(BSIZE-2);
	jj=1 + (gid%maxblock)*(BSIZE-2);

	int row1, offs2;

	/*
		copy state to local
	*/
	int A[BSIZE*BSIZE];
	int B[BSIZE*BSIZE];
	for (i = 0; i<BSIZE; i++){
		row1 = i*BSIZE;
		offs2 = (ii + i -1)*N + jj - 1;
		for (j = 0; j<BSIZE; j++){
			A[row1 + j] = vec1[offs2 + j];
		}
	}

	
	for (i = 0; i<BSIZE; i++){
		row1 = i*BSIZE;
		for (j = 0; j<BSIZE; j++){

			int state=A[row1+j];

			/* 
				collision detection
			*/
			if (state == 10){
				state = 5;
			}
			else if (state == 5){
				state = 10;
			}


			/*
				reflect particles at edges
			*/
			int cell_pos=(ii+i-1)*N+jj+j-1;
			//top
			if (cell_pos/N == 1  && ((state & 8) == 8)){
				state = state | 2;
				state = state & 7;
			}
			//bottom
			if (cell_pos/N == N-2 && ((state & 2) == 2)){
				state = state | 8;
				state = state & 13;
			}
			//left
			if (cell_pos%N == 1 && ((state & 1) == 1)){
				state = state | 4;
				state = state & 14;
			}
			//right
				if ( cell_pos%N == N-2 && ((state & 4) == 4)){
				state = state | 1;
				state = state & 11;
			}

			A[row1+j]=state;
		}
	}
	
	/*
		Transport
			different for loop because it is
			not done for the edges of the local block
	*/
	for (i = 1; i<BSIZE-1; i++){
		for (j = 1; j<BSIZE-1; j++){
			int snew = 0;
			if ((A[(i-1)*BSIZE+j] & 2) == 2){ // from top going down
				snew = snew | 2;
			}
			if ((A[i*BSIZE+j-1] & 4) == 4){ // from left going right
				snew = snew | 4;
			}
			if ((A[i*BSIZE+j+1] & 1) == 1){ // from right going left
				snew = snew | 1;
			}
			if ((A[(i+1)*BSIZE+j] & 8) == 8){ // from bottom going up
				snew = snew | 8;
			}
			B[i*BSIZE+j] = snew;
		}
	}

	/* 
		copy back to global
	*/
	for (i = 1; i<BSIZE-1; i++){
		row1 = i*BSIZE;
		offs2 = (ii + i - 1 )*N + jj - 1;
		for (j = 1; j<BSIZE-1; j++){
			vec2[offs2 + j] = B[row1 + j];
		}
	}

	/*
		write to image output
	*/
	//lookup table for particle number
	__const uint num_of_parts[16]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
	//__const int num_of_parts[16]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
	int iim;
	for (i = 1; i<BSIZE-1; i++){
		row1 = i*BSIZE;
		iim = ii + i - 1  ;
		for (j = 1; j<BSIZE-1; j++){
			uint lum=num_of_parts[B[row1 + j]];
			float lumf=0.24 * num_of_parts[B[row1 + j]];
			//write_imageui (output, (int2)(iim, jj+j-1), (uint4)(lum,lum,lum,lum) );
			//write_imageui (output, (int2)(iim, jj+j-1), num_of_parts[B[row1 + j]]);
			//write_imagef (output, (int2)(iim, jj+j-1), 2* num_of_parts[B[row1 + j]]);
			write_imagef (output, (int2)(iim, jj+j-1),(float4)(1-lumf,1-lumf,1-lumf,1));
		}
	}

	return;
}