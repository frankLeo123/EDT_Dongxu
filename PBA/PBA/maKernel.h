#ifndef __MAKERNEL__
#define __MAKERNEL__


#define TOID(x, y, size)    (__mul24((y), (size)) + (x))


/************************************************************************/
/* phase 1                                                              */
/************************************************************************/

// Transpose a square matrix
__global__ void maKernelTranspose(short2 *data, int size)
{
	__shared__ short2 block1[TILE_DIM][TILE_DIM + 1];
	__shared__ short2 block2[TILE_DIM][TILE_DIM + 1];

	int blockIdx_y = blockIdx.x;
	int blockIdx_x = blockIdx.x+blockIdx.y;

	if (blockIdx_x >= gridDim.x)
		return ; 

	int blkX, blkY, x, y, id1, id2; 
	short2 pixel; 

	blkX = __mul24(blockIdx_x, TILE_DIM); 
	blkY = __mul24(blockIdx_y, TILE_DIM); 

	x = blkX + threadIdx.x;
	y = blkY + threadIdx.y;
	id1 = __mul24(y, size) + x;

	x = blkY + threadIdx.x;
	y = blkX + threadIdx.y;
	id2 = __mul24(y, size) + x;

	// read the matrix tile into shared memory
	for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
		block1[threadIdx.y + i][threadIdx.x] = tex1Dfetch(maTexColor, id1 + __mul24(i, size));
		block2[threadIdx.y + i][threadIdx.x] = tex1Dfetch(maTexColor, id2 + __mul24(i, size));
	}

	__syncthreads();

	// write the transposed matrix tile to global memory
	for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
		pixel = block1[threadIdx.x][threadIdx.y + i];
		data[id2 + __mul24(i, size)] = make_short2(pixel.y, pixel.x); 
		pixel = block2[threadIdx.x][threadIdx.y + i];
		data[id1 + __mul24(i, size)] = make_short2(pixel.y, pixel.x); 
	}
}

__global__ void maKernelFloodDown(short2 *output, int size, int bandSize) 
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x; 
	int ty = blockIdx.y * bandSize; 
	int id = TOID(tx, ty, size); 

	short2 pixel1, pixel2; 

	pixel1 = make_short2(MARKER, MARKER); 

	for (int i = 0; i < bandSize; i++, id += size) {
		pixel2 = tex1Dfetch(maTexColor, id); 

		if (pixel2.x != MARKER) 
			pixel1 = pixel2; 

		output[id] = pixel1; 
	}
}

__global__ void maKernelFloodUp(short2 *output, int size, int bandSize) 
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x; 
	int ty = (blockIdx.y+1) * bandSize - 1; 
	int id = TOID(tx, ty, size); 

	short2 pixel1, pixel2; 
	int dist1, dist2; 

	pixel1 = make_short2(MARKER, MARKER); 

	for (int i = 0; i < bandSize; i++, id -= size) {
		dist1 = abs(pixel1.y - ty + i); 

		pixel2 = tex1Dfetch(maTexColor, id); 
		dist2 = abs(pixel2.y - ty + i); 

		if (dist2 < dist1) 
			pixel1 = pixel2; 

		output[id] = pixel1; 
	}
}

__global__ void maKernelPropagateInterband(short2 *output, int size, int bandSize) 
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x; 
	int inc = __mul24(bandSize, size); 
	int ny, nid, nDist; 
	short2 pixel; 

	// Top row, look backward
	int ty = __mul24(blockIdx.y, bandSize); 
	int topId = TOID(tx, ty, size); 
	int bottomId = TOID(tx, ty + bandSize - 1, size); 

	pixel = tex1Dfetch(maTexColor, topId); 
	int myDist = abs(pixel.y - ty); 

	for (nid = bottomId - inc; nid >= 0; nid -= inc) {
		pixel = tex1Dfetch(maTexColor, nid); 

		if (pixel.x != MARKER) { 
			nDist = abs(pixel.y - ty); 

			if (nDist < myDist) 
				output[topId] = pixel; 

			break;	
		}
	}

	// Last row, look downward
	ty = ty + bandSize - 1; 
	pixel = tex1Dfetch(maTexColor, bottomId); 
	myDist = abs(pixel.y - ty); 

	for (ny = ty + 1, nid = topId + inc; ny < size; ny += bandSize, nid += inc) {
		pixel = tex1Dfetch(maTexColor, nid); 

		if (pixel.x != MARKER) { 
			nDist = abs(pixel.y - ty); 

			if (nDist < myDist) 
				output[bottomId] = pixel; 

			break; 
		}
	}
}

__global__ void maKernelUpdateVertical(short2 *output, int size, int band, int bandSize) 
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x; 
	int ty = blockIdx.y * bandSize; 

	short2 top = tex1Dfetch(maTexLinks, TOID(tx, ty, size)); 
	short2 bottom = tex1Dfetch(maTexLinks, TOID(tx, ty + bandSize - 1, size)); 
	short2 pixel; 

	int dist, myDist; 

	int id = TOID(tx, ty, size); 

	for (int i = 0; i < bandSize; i++, id += size) {
		pixel = tex1Dfetch(maTexColor, id); 
		myDist = abs(pixel.y - (ty + i)); 

		dist = abs(top.y - (ty + i)); 
		if (dist < myDist) { myDist = dist; pixel = top; }

		dist = abs(bottom.y - (ty + i)); 
		if (dist < myDist) pixel = bottom; 

		output[id] = pixel; 
	}
}

/************************************************************************/
/* phase 2 precalculated                                                */
/************************************************************************/


__global__ void  _B3_maKernelLinkDownScan(short2 *out,int size,int bandSize)
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * bandSize;
	int tid = TOID(tx,ty,size);
	int bx = tx;
	int by = ty + bandSize - 1;
	int bid = TOID(bx,by,size);	

	short pre = MARKER;

	for (int i=tid,k=ty;i<=bid;i+=size,k++)
	{
		short2 p = tex1Dfetch(maTexColor,i);
		out[i] =make_short2( pre , MARKER);
		if (p.x!=MARKER) pre = k;
	}
	//out[bid].x = pre;
}

__global__ void  _B4_maKernelLinkDownEndPoint(short2 *out,int size,int bandSize)
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * bandSize;
	int tid = TOID(tx,ty,size);
	int bx = tx;
	int by = ty + bandSize - 1;
	int bid = TOID(bx,by,size);

	int inc = size*bandSize;
	int prebid = bid - inc;
	int preby = by - bandSize;
	
	short pre=MARKER;
	for (int k=preby,i=prebid;k>=0;k-=bandSize,i-=inc)
	{
		short2 p = tex1Dfetch(maTexColor,i);
		short2 q = tex1Dfetch(maTexLinks,i);
		if (p.x!=MARKER)
		{
			pre = p.x;
			if (q.x!=MARKER) pre = k;
			break;
		}
	}
	out[tid] = make_short2(pre,MARKER); 
}


__global__ void  _B5_maKernelLinkDownInner(short2 *out,int size,int bandSize)
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * bandSize;
	int tid = TOID(tx,ty,size);
	int bx = tx;
	int by = ty + bandSize - 1;
	int bid = TOID(bx,by,size);

	short2 p = tex1Dfetch(maTexLinks,tid);
	short pre = p.x;

	for (int i=tid,k=ty;i<=bid;i+=size,k++)
	{
		short2 n = tex1Dfetch(maTexColor,i);
		if (n.x != MARKER) break;
		out[i] = make_short2(pre,MARKER);
	}
}

__global__ void  _B6_maKernelLinkUpScan(short2 *out,int size,int bandSize)
{

}

__global__ void  _B7_maKernelLinkUpEndPoint(short2 *out,int size,int bandSize)
{

}
__global__ void  _B8_maKernelLinkUpInner(short2 *out,int size,int bandSize)
{

}

/************************************************************************/
/* phase 2                                                              */
/************************************************************************/

/**
__device__ short2 getNearestPoint(int x,int bx,int ex,int h,int size)
{
	short2 ans;
	int minDis = INF;
	int rad;
	int k=0;
	if (x<bx) k = bx-x;
	if (x>ex) k = x-ex;
	int kk=k*k;
	int lb = bx,rb = ex;
	while (1)
	{
		int f = 0;
		if (x+k<=rb)
		{
			f=1;
			short2 np = tex1Dfetch(maTexColor,TOID(h,x+k,size));
			int dx = np.x - h; 
			int dis1 = dx*dx + kk;
			if (dis1<minDis)
			{
				minDis = dis1;
				ans = np;
				rad = (int)sqrt(float(dis1));
				rb = x+rad<ex?x+rad:ex;
				lb = x-rad>bx?x-rad:bx;
			}
		}
		if (x-k>=lb)
		{
			f=1;
			short2 np = tex1Dfetch(maTexColor,TOID(h,x-k,size));
			int dx = np.x - h;
			int dis2 = dx*dx + kk;
			if (dis2<minDis)
			{
				minDis = dis2;
				ans = np;
				rad = (int)sqrt(float(dis2));
				rb = x+rad<ex?x+rad:ex;
				lb = x-rad>bx?x-rad:bx;
			}
		}
		kk=kk+2*k+1;
		k++;
		if (f==0) break;
	}
	return ans;
}
/**/

/**/
__device__ short2 getNearestPoint(int x,int bx,int ex,int h,int size)
{
	short2 ans;
	int minDis = INF;
	int rad;
	int k=0;
	if (x<bx) k = bx-x;
	if (x>ex) k = x-ex;
	//int kk=k*k;
	int lb = bx,rb = ex;
	int up=x-k,down=x+k;
	while (1)
	{
		int f = 0;
		if (down<=rb&&down!=MARKER)
		{
			f=1;
			short2 np = tex1Dfetch(maTexColor,TOID(h,down,size));
			short2 lp = np;
			if (np.x==MARKER)
			{
				lp = tex1Dfetch(maTexNext,TOID(h,down,size));
				np = tex1Dfetch(maTexColor,TOID(h,lp.y,size));
				down = lp.y;
			}
			else down++;
			int dx = np.x - h; 
			int dy = np.y - x;
			int dis1 = dx*dx + dy*dy;
			if (dis1<minDis)
			{
				minDis = dis1;
				ans = np;
				rad = (int)sqrt(float(dis1));
				rb = x+rad<ex?x+rad:ex;
				lb = x-rad>bx?x-rad:bx;
			}
			
		}
		if (up>=lb)
		{
			f=1;
			short2 np = tex1Dfetch(maTexColor,TOID(h,up,size));
			short2 lp = np;	
			if (np.x==MARKER)
			{
				lp = tex1Dfetch(maTexNext,TOID(h,up,size));
				np = tex1Dfetch(maTexColor,TOID(h,lp.x,size));
				up = lp.x;
			}
			else up--;
			int dx = np.x - h;
			int dy = np.y - x;
			int dis2 = dx*dx + dy*dy;
			if (dis2<minDis)
			{
				minDis = dis2;
				ans = np;
				rad = (int)sqrt(float(dis2));
				rb = x+rad<ex?x+rad:ex;
				lb = x-rad>bx?x-rad:bx;
			}
			
		}
		//kk=kk+2*k+1;
		//k++;
		if (f==0) break;
	}
	return ans;
}
/**/

__global__ void _B1_maKernelCalEndPoint(short2 *out,int size, int bandSize)
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * bandSize;
	int tid = TOID(tx,ty,size);
	int bx = tx;
	int by = ty + bandSize - 1;
	int bid = TOID(bx,by,size);	
	out[tid] = getNearestPoint(ty,0,size-1,tx,size);
	out[bid] = getNearestPoint(by,0,size-1,bx,size);
	/*for (int i=ty,k=tid;i<=by;i++,k+=size)
	{
		out[k] = getNearestPoint(i,0,size-1,tx,size);
	}*/
}

__global__ void _B2_maKernelCalInnerBand(short2 *out,int size,int bandSize)
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * bandSize;
	int tid = TOID(tx,ty,size);
	int bx = tx;
	int by = ty + bandSize - 1;
	int bid = TOID(bx,by,size);

	short2 itp = tex1Dfetch(maTexLinks,tid);
	short2 ibp = tex1Dfetch(maTexLinks,bid);
	
	if (itp.x==ibp.x&&itp.y==ibp.y)
	{
		for (int i = tid + size ; i<bid ; i+=size) out[i]=itp;
		return ;
	}

	int sck[15][3];
	short2 sckp[15][2];
	int st = 0 ;
	sck[st][0] = ty;
	sck[st][1] = by;
	sck[st][2] = 0;
	sckp[st][0] = itp;
	sckp[st][1] = ibp;
	st++;
	
	while (st>=1)
	{
		
		int t = sck[st-1][0];
		short2 tp = sckp[st-1][0];
		int b = sck[st-1][1];
		short2 bp = sckp[st-1][1];
		int ip = sck[st-1][2];
		
		if (ip==0)
		{
			int midU = (t+b)/2;
			if (t<midU)
			{
				short2 midUp = getNearestPoint(midU,tp.y,bp.y,tx,size);
				if (midUp.x==tp.x&&midUp.y==tp.y)
				{
					for (int i=t+1;i<=midU;i++)
						out[TOID(tx,i,size)] = midUp;
					sck[st-1][2]++;
				}
				else
				{
					out[TOID(tx,midU,size)] = midUp;
					sck[st][0]=t;
					sckp[st][0]=tp;
					sck[st][1]=midU;
					sckp[st][1]=midUp;
					sck[st][2]=0;
					sck[st-1][2]++;
					st++;
					continue;
				}
			}
		}
		
		if (ip==1||ip==0)
		{
			int midD = ((t+b)/2) + 1;
			if (midD<b)
			{
				short2 midDp = getNearestPoint(midD,tp.y,bp.y,tx,size); 
				
				if (midDp.x==bp.x&&midDp.y==bp.y)
				{
					for (int i = midD ; i<=b-1;i++)
						out[TOID(tx,i,size)] = midDp;
					sck[st-1][2]++;
				}
				else
				{
					out[TOID(tx,midD,size)] = midDp;
					sck[st][0]=midD;
					sckp[st][0]=midDp;
					sck[st][1]=b;
					sckp[st][1]=bp;
					sck[st][2]=0;
					sck[st-1][2]++;
					st++;
					continue;
				}
			}
		}
		st--;
	}
	
}


#endif