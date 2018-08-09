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
int *dindex,*dreindex;
int *hindex,*hreindex;
int *pnt,*npnt,*pppl,*icpl;
int maMemSize;             // Size (in bytes) of a texture
int maTexSize;             // Texture size (squared texture)

texture<short2> maTexColor; 
texture<short2> maTexLinks;
texture<short2> maTexP1;
texture<short2> maTexIcp;
texture<short2> maTexPpp;
texture<int> maTexIndex;
texture<int> maTexReindex;
texture<int> maTexPnt;
texture<int> maTexNpnt;
texture<int> maTexPppl;
texture<int> maTexIcpl;

/********* Kernels ********/
#include "maKernel.h"

///////////////////////////////////////////////////////////////////////////
//
// Initialize necessary memory for 2D Voronoi Diagram computation
// - textureSize: The size of the Discrete Voronoi Diagram (width = height)
//
///////////////////////////////////////////////////////////////////////////

void makePattern(int step)
{
	int mstep=step,mstep2;
	int num=0;
	for (int i=0;i<maTexSize;i+=step) 
	{
		hindex[i]=num;
		hreindex[num]=i;
		num++;
	}
	while(mstep>=2)
	{
		mstep2 = mstep/2;
		for (int i=mstep2;i<maTexSize;i+=mstep)
		{
			hindex[i]=num;
			hreindex[num]=i;
			num++;
		}
		mstep = mstep2;
	}
	cudaMemcpy(dindex,hindex,maTexSize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dreindex,hreindex,maTexSize*sizeof(int),cudaMemcpyHostToDevice);
}

void maInitialization(int textureSize)
{
	maTexSize = textureSize; 
	maMemSize = maTexSize *( maTexSize ) * sizeof(short2) ; 

	maTextures = (short2 **) malloc(4 * sizeof(short2 *)); 

	// Allocate 2 textures
	cudaMalloc((void **) &maTextures[0], maMemSize); 
	cudaMalloc((void **) &maTextures[1], maMemSize); 
	cudaMalloc((void **) &maTextures[2], maMemSize); 
	cudaMalloc((void **) &maTextures[3], maMemSize); 
	cudaMalloc((void **) &dindex, textureSize*sizeof(int)); 
	cudaMalloc((void **) &dreindex, textureSize*sizeof(int));
	cudaMalloc((void **) &pnt, textureSize*sizeof(int)); 
	cudaMalloc((void **) &npnt, textureSize*sizeof(int)); 
	cudaMalloc((void **) &pppl, textureSize*sizeof(int)); 
	cudaMalloc((void **) &icpl, textureSize*sizeof(int)); 
	hindex = (int *) malloc(maTexSize*sizeof(int));
	hreindex = (int *) malloc(maTexSize*sizeof(int));
	cudaMemset(icpl,-1,textureSize*sizeof(int));
	cudaMemset(maTextures[2],127,maMemSize);
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
	cudaFree(maTextures[3]); 
	free(maTextures); 
	cudaFree(dindex);
	cudaFree(dreindex);
	cudaFree(pnt);
	cudaFree(npnt);
	cudaFree(pppl);
	cudaFree(icpl);
}


// Copy input to GPU 
void maInitializeInput(short *input)
{
	cudaMemcpy(maTextures[0], input, maMemSize, cudaMemcpyHostToDevice); 
}

// In-place transpose a squared texture. 
// Block orders are modified to optimize memory access. 
// Point coordinates are also swapped. 
/*void maTranspose(short2 *texture)
{
	dim3 block(TILE_DIM, BLOCK_ROWS); 
	dim3 grid(maTexSize / TILE_DIM, maTexSize / TILE_DIM); 

	cudaBindTexture(0, maTexColor, texture); 
	maKernelTranspose<<< grid, block >>>(texture, maTexSize); 
	cudaUnbindTexture(maTexColor); 
}*/
void maTranspose(short2 *texture)
{
	short *tmp1 = (short *)malloc(maMemSize);
	short *tmp2 = (short *)malloc(maMemSize);
	cudaMemcpy(tmp1,texture,maMemSize,cudaMemcpyDeviceToHost);
	for (int i=0;i<maTexSize;i++)
	{
		for (int j=0;j<maTexSize;j++)
		{
			int idx1 = i + j*maTexSize;
			int idx2 = j + i*maTexSize;
			tmp2[idx2*2] = tmp1[idx1*2+1];
			tmp2[idx2*2+1] = tmp1[idx1*2];
		}
	}
	cudaMemcpy(texture,tmp2,maMemSize,cudaMemcpyHostToDevice);
}

// Phase 1 of PBA. m1 must divides texture size
void maPhase1() 
{
	dim3 block = dim3(BLOCKSIZE);   
	dim3 grid = dim3((maTexSize / block.x)); 

	// Flood vertically in their own bands
	cudaBindTexture(0, maTexColor, maTextures[0]); 
	maKernelFloodDown<<< grid, block >>>(maTextures[1], maTexSize); 
	cudaUnbindTexture(maTexColor); 

	cudaBindTexture(0, maTexColor, maTextures[1]); 
	maKernelFloodUp<<< grid, block >>>(maTextures[1], maTexSize); 
	cudaUnbindTexture(maTexColor); 
}

// Phase 2 of PBA. m2 must divides texture size
void maPhase2(int STEP) 
{	
	int step = STEP;
	makePattern(step);
	cudaBindTexture(0,maTexIndex,dindex);
	cudaBindTexture(0,maTexReindex,dreindex);

	cudaBindTexture(0,maTexColor,maTextures[1]);
	dim3 block,grid;
	block = dim3(BLOCKSIZE);
	grid = dim3(((maTexSize/STEP)/block.x) + 1 );
	maKernelColorInit1<<<grid,block>>>(maTextures[2],icpl,maTexSize,0,step);
	cudaUnbindTexture(maTexColor);
	/*cudaBindTexture(0,maTexIcp,maTextures[2]);
	cudaBindTexture(0,maTexIcpl,icpl);
	maKernelColorInit2<<<grid,block>>>(maTextures[0],maTexSize,0,step);
	cudaUnbindTexture(maTexIcpl);
	cudaUnbindTexture(maTexIcp);*/

	while (step>=2)
	{
		int step2 = step/2;
		int taskNum = (maTexSize - step2)/step;
		taskNum++;

		block = dim3(BLOCKSIZE);
		grid = dim3((taskNum/BLOCKSIZE)  );
		cudaBindTexture(0,maTexIcpl,icpl);
		cudaBindTexture(0,maTexP1,maTextures[1]);
		cudaBindTexture(0,maTexIcp,maTextures[2]);
		maKernelColor15<<<grid,block>>>(maTextures[3],pppl,maTexSize,step2,step);
		cudaUnbindTexture(maTexIcp);
		cudaUnbindTexture(maTexIcpl);
		cudaUnbindTexture(maTexP1);


		//cudaMemcpy(xxx, maTextures[5], maMemSize, cudaMemcpyDeviceToHost); 

		block = dim3(BLOCKSIZE);
		grid = dim3((taskNum/BLOCKSIZE)   );
		cudaBindTexture(0,maTexPppl,pppl);
		cudaBindTexture(0,maTexPpp,maTextures[3]);
		maKernelColor2<<<grid,block>>>(maTextures[2],icpl,maTexSize,step2,step);
		cudaUnbindTexture(maTexPppl);
		cudaUnbindTexture(maTexPpp);




		block = dim3(BLOCKSIZE);
		grid = dim3((maTexSize/BLOCKSIZE)  );
		cudaBindTexture(0,maTexIcp,maTextures[2]);
		cudaBindTexture(0,maTexIcpl,icpl);
		maKernelColor3<<<grid,block>>>(maTextures[0],maTexSize);
		cudaUnbindTexture(maTexIcp);
		cudaUnbindTexture(maTexIcpl);
		step = step2;
	}

	cudaUnbindTexture(maTexIndex);
	cudaUnbindTexture(maTexReindex);

	/**
	int step = STEP;
	cudaBindTexture(0,maTexColor,maTextures[1]);
	cudaBindTexture(0,maTexIndex,index);
	dim3 block,grid;
	block = dim3(BLOCKSIZE);
	grid = dim3(((maTexSize/STEP)/block.x)+1);
	maKernelColorInit<<<grid,block>>>(maTextures[0],maTextures[2],maTexSize,0,step);

	while(step>=2)
	{
		int step2 = step/2;
		int taskNum = (maTexSize-step2)/step;
		taskNum+=1;
		block = dim3(BLOCKSIZE);
		grid = dim3((taskNum/BLOCKSIZE)+1);
		cudaBindTexture(0,maTexLinks,maTextures[0]);
		maKernelColorLine<<<grid,block>>>(maTextures[0],maTextures[2],maTexSize,step2,step);
		cudaUnbindTexture(maTexLinks);
		step=step2;
	}

	cudaUnbindTexture(maTexColor);
	cudaUnbindTexture(maTexIndex);
	/**/
	/**
	cudaBindTexture(0,maTexColor,maTextures[1]);
	dim3 block,grid;
	block = dim3(BLOCKSIZE);
	grid = dim3(maTexSize/block.x);
	maKernelTest1<<<grid,block>>>(maTextures[2],index,maTexSize);
	cudaBindTexture(0,maTexLinks,maTextures[2]);
	cudaBindTexture(0,maTexIndex,index);
	maKernelTest2<<<grid,block>>>(maTextures[0],maTexSize);
	cudaUnbindTexture(maTexColor);
	cudaUnbindTexture(maTexLinks);
	cudaUnbindTexture(maTexIndex);/**/
}

void maCompute(int STEP,short* output)
{
	//cudaMemcpy(xx, maTextures[0], maMemSize, cudaMemcpyDeviceToHost);
	maPhase1(); 
//cudaMemcpy(xx, maTextures[1], maMemSize, cudaMemcpyDeviceToHost); 
	maTranspose(maTextures[1]); 
//	cudaMemcpy(xx, maTextures[1], maMemSize, cudaMemcpyDeviceToHost); 
	maPhase2(STEP); 

	maTranspose(maTextures[0]); 
	
}

void rerange(short *outData)
{
	short *tmp = (short *)malloc(maMemSize);
	for (int i=0;i<maTexSize;i++)
	{
		int k = hreindex[i];
		for (int j=0;j<maTexSize;j++)
		{
			int idtmp = j + k*maTexSize;
			int id = j+ i*maTexSize;
			tmp[idtmp*2] = outData[id*2];
			tmp[idtmp*2+1] = outData[id*2+1];
		}
	}
	memcpy(outData,tmp,maMemSize);
}

// Compute 2D Voronoi diagram
// Input: a 2D texture. Each pixel is represented as two "short" integer. 
//    For each site at (x, y), the pixel at coordinate (x, y) should contain 
//    the pair (x, y). Pixels that are not sites should contain the pair (MARKER, MARKER)
// See original paper for the effect of the three parameters: 
//    phase1Band, phase2Band, phase3Band
// Parameters must divide textureSize
void ma(short *input, short *output, int STEP) 
{
	// Initialization
	maInitializeInput(input); 
	makePattern(STEP);
	// Computation
	maCompute(STEP,output); 

	// Copy back the result
	cudaMemcpy(output, maTextures[0], maMemSize, cudaMemcpyDeviceToHost); 
	//rerange(output);
	
}