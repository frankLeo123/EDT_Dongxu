#include <device_functions.h>

#include "ma.h"

// Parameters for CUDA kernel executions
#define BLOCKX		16
#define BLOCKY		16
#define BLOCKSIZE	64
#define TILE_DIM	32
#define BLOCK_ROWS	8


/****** Global Variables *******/
short2 **maTextures;       // Two textures used to compute 2D Voronoi Diagram

int maMemSize;             // Size (in bytes) of a texture
int maTexSize;             // Texture size (squared texture)

texture<short2> maTexColor; 
texture<short2> maTexLinks;
texture<short2> maTexNext;


/********* Kernels ********/
#include "maKernel.h"

///////////////////////////////////////////////////////////////////////////
//
// Initialize necessary memory for 2D Voronoi Diagram computation
// - textureSize: The size of the Discrete Voronoi Diagram (width = height)
//
///////////////////////////////////////////////////////////////////////////
void maInitialization(int textureSize)
{
	maTexSize = textureSize; 
	maMemSize = maTexSize * maTexSize * sizeof(short2); 

	maTextures = (short2 **) malloc(3 * sizeof(short2 *)); 

	// Allocate 2 textures
	cudaMalloc((void **) &maTextures[0], maMemSize); 
	cudaMalloc((void **) &maTextures[1], maMemSize); 
	cudaMalloc((void **) &maTextures[2], maMemSize); 
}

///////////////////////////////////////////////////////////////////////////
//
// Deallocate all allocated memory
//
///////////////////////////////////////////////////////////////////////////
void maDeInitialization()
{
	cudaFree(maTextures[0]); 
	cudaFree(maTextures[1]); 
	cudaFree(maTextures[2]); 
	free(maTextures); 
}


// Copy input to GPU 
void maInitializeInput(short *input)
{
	cudaMemcpy(maTextures[0], input, maMemSize, cudaMemcpyHostToDevice); 
}

// In-place transpose a squared texture. 
// Block orders are modified to optimize memory access. 
// Point coordinates are also swapped. 
void maTranspose(short2 *texture)
{
	dim3 block(TILE_DIM, BLOCK_ROWS); 
	dim3 grid(maTexSize / TILE_DIM, maTexSize / TILE_DIM); 

	cudaBindTexture(0, maTexColor, texture); 
	maKernelTranspose<<< grid, block >>>(texture, maTexSize); 
	cudaUnbindTexture(maTexColor); 
}

// Phase 1 of PBA. m1 must divides texture size
void maPhase1(int m1) 
{
	dim3 block = dim3(BLOCKSIZE);   
	dim3 grid = dim3((maTexSize / block.x)*m1); 

	// Flood vertically in their own bands
	cudaBindTexture(0, maTexColor, maTextures[0]); 
	maKernelFloodDown<<< grid, block >>>(maTextures[1], maTexSize, maTexSize / m1); 
	cudaUnbindTexture(maTexColor); 

	cudaBindTexture(0, maTexColor, maTextures[1]); 
	maKernelFloodUp<<< grid, block >>>(maTextures[1], maTexSize, maTexSize / m1); 

	// Passing information between bands
	grid = dim3(maTexSize / block.x, m1); 
	maKernelPropagateInterband<<< grid, block >>>(maTextures[0], maTexSize, maTexSize / m1); 

	cudaBindTexture(0, maTexLinks, maTextures[0]); 
	maKernelUpdateVertical<<< grid, block >>>(maTextures[1], maTexSize, m1, maTexSize / m1); 
	cudaUnbindTexture(maTexLinks); 
	cudaUnbindTexture(maTexColor); 
}

short xx1[4000000],xx2[4000000],xx3[4000000];

int todo(int x,int y, int size)
{
	return y*size + x;
}

void checkCal()
{
	for (int i=0;i<maTexSize;i++)
	{
		int pre = MARKER;
		for (int j=0;j<maTexSize;j++)
		{
			int x=i,y=j;
			int id = todo(x,y,maTexSize);
			xx2[id*2]=pre;
			if (xx1[id*2]!=MARKER) pre=y;
		}
		int bck = MARKER;
		for (int j=maTexSize-1;j>=0;j--)
		{
			int x=i,y=j;
			int id = todo(x,y,maTexSize);
			xx2[id*2+1]=bck;
			if (xx1[id*2]!=MARKER) bck=y;
		}
	}

}

int check()
{
	int ans=0;
	for (int i=0;i<maTexSize;i++)
	{
		for (int j=0;j<maTexSize;j++)
		{
			int x = i,y=j;
			int id =  todo(x,y,maTexSize);
			if (xx1[id*2]!=xx2[id*2]) 
				ans++;
		}
	}
	return ans;
}

int check1()
{
	int ans=0;
	for (int i=0;i<maTexSize;i++)
	{
		for (int k=0;k<16;k++)
		{
			int pre = MARKER;
			for (int j=k*(maTexSize/16);j<(k+1)*(maTexSize/16);j++)
			{
				int x=i,y=j;
				int id = todo(x,y,maTexSize);
				if (xx3[id*2]!=pre) 
					ans++;
				if (xx1[id*2]!=MARKER) pre=y;
			}
		}
		
	}
	return ans;
}

int check2()
{
	int ans=0;
	for (int i=0;i<maTexSize;i++)
	{
		for (int k=0;k<16;k++)
		{
			int x=i,ty=k*(maTexSize/16);
			int tid = todo(x,ty,maTexSize);
			if (xx3[tid*2]!=xx2[tid*2]) 
				ans++;
		}
	}
	return ans;
}

void maPhase2pre(short2 * phase1Res , int m2pre)
{
	int num;
	cudaMemcpy(xx1, phase1Res, maMemSize, cudaMemcpyDeviceToHost);
	checkCal();
	cudaMemcpy(maTextures[2],xx2,maMemSize,cudaMemcpyHostToDevice);
/*
	dim3 block = dim3(BLOCKSIZE);
	dim3 grid = dim3(maTexSize/block.x,m2pre);

	cudaBindTexture(0,maTexColor,phase1Res);
	_B3_maKernelLinkDownScan<<<grid,block>>>(maTextures[2],maTexSize,maTexSize/m2pre);
	cudaUnbindTexture(maTexColor);
//	cudaMemcpy(xx3, maTextures[2], maMemSize, cudaMemcpyDeviceToHost); 
//	num = check1();
//	printf("\nnum:%d\n",num);
	cudaBindTexture(0,maTexLinks,maTextures[1]);
	cudaBindTexture(0,maTexColor,maTextures[2]);
	_B4_maKernelLinkDownEndPoint<<<grid,block>>>(maTextures[0],maTexSize,maTexSize/m2pre);
	cudaUnbindTexture(maTexLinks);
//	cudaMemcpy(xx3, maTextures[0], maMemSize, cudaMemcpyDeviceToHost); 
//	num = check2();
//	printf("\nnum:%d\n",num);
	cudaBindTexture(0,maTexLinks,maTextures[0]);
	_B5_maKernelLinkDownInner<<<grid,block>>>(maTextures[2],maTexSize,maTexSize/m2pre);
	cudaUnbindTexture(maTexLinks);
	cudaUnbindTexture(maTexColor);

	cudaMemcpy(xx1, maTextures[2], maMemSize, cudaMemcpyDeviceToHost); 
	num = check();
	printf("\nnum:%d\n",num);
*/
}

// Phase 2 of PBA. m2 must divides texture size
void maPhase2(int m2) 
{
	dim3 block = dim3(BLOCKSIZE);
	dim3 grid = dim3(maTexSize/block.x,m2);

	cudaBindTexture(0,maTexNext,maTextures[2]);
	cudaBindTexture(0,maTexColor,maTextures[1]);

	_B1_maKernelCalEndPoint<<<grid,block>>>(maTextures[0],maTexSize,maTexSize/m2);

	cudaBindTexture(0,maTexLinks,maTextures[0]);

	_B2_maKernelCalInnerBand<<<grid,block>>>(maTextures[0],maTexSize,maTexSize/m2);

	cudaUnbindTexture(maTexColor);
	cudaUnbindTexture(maTexLinks);
	cudaUnbindTexture(maTexNext);
}


void maCompute(int floodBand, int linkBand, int maurerBand,short* output)
{
	maPhase1(floodBand); 

	maTranspose(maTextures[1]); 

	maPhase2pre(maTextures[1] , linkBand);

	maPhase2(maurerBand); 

	maTranspose(maTextures[0]); 
	
}

// Compute 2D Voronoi diagram
// Input: a 2D texture. Each pixel is represented as two "short" integer. 
//    For each site at (x, y), the pixel at coordinate (x, y) should contain 
//    the pair (x, y). Pixels that are not sites should contain the pair (MARKER, MARKER)
// See original paper for the effect of the three parameters: 
//    phase1Band, phase2Band, phase3Band
// Parameters must divide textureSize
void ma(short *input, short *output, int floodBand, int linkBand , int maurerBand) 
{
	// Initialization
	maInitializeInput(input); 

	// Computation
	maCompute(floodBand, linkBand , maurerBand,output); 

	// Copy back the result
	cudaMemcpy(output, maTextures[0], maMemSize, cudaMemcpyDeviceToHost); 

}