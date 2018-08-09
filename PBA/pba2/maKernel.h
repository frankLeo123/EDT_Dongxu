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

__global__ void maKernelFloodDown(short2 *output, int size) 
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x; 
	if (tx>=size) return ;
	int ty = 0; 
	int id = TOID(tx, ty, size); 

	short2 pixel1, pixel2; 

	pixel1 = make_short2(MARKER, MARKER); 

	for (int i = 0; i < size; i++, id += size) {
		pixel2 = tex1Dfetch(maTexColor, id); 

		if (pixel2.x != MARKER) 
			pixel1 = pixel2; 

		output[id] = pixel1; 
	}
}

__global__ void maKernelFloodUp(short2 *output, int size) 
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x; 
	if (tx>=size) return ;
	int ty = size - 1; 
	int id = TOID(tx, ty, size); 

	short2 pixel1, pixel2; 
	int dist1, dist2; 

	pixel1 = make_short2(MARKER, MARKER); 

	for (int i = 0; i < size; i++, id -= size) {
		dist1 = abs(pixel1.y - ty + i); 

		pixel2 = tex1Dfetch(maTexColor, id); 
		dist2 = abs(pixel2.y - ty + i); 

		if (dist2 < dist1) 
			pixel1 = pixel2; 

		output[id] = pixel1; 
	}
}

__device__ float insertpnt(float x1,float y1,float x2, float y2,float x0)
{
	float xM = (x1+x2)/2.0f;
	float yM = (y1+y2)/2.0f;
	float nx = x2 - x1;
	float ny = y2 - y1;
	if (ny==0) return 1000000.0;
	else 
		return yM+nx*(xM-x0)/ny;
}

__device__ int remove(short2 lst2,short2 lst1,short2 now,int h)
{
	float i1 = insertpnt(lst2.x,lst2.y,lst1.x,lst1.y,h);
	float i2 = insertpnt(lst1.x,lst1.y,now.x,now.y,h);
	return (i1>i2);
}

__device__ void maColor(short2 *out,short2 *g,int size,int h,int b,int by,int e,int ey)
{
	//short2 g[200];
	int l=-1;
	for (int i=by;i<=ey;i++)
	{
		int id = h+__mul24(i,size);
		short2 dm_row = tex1Dfetch(maTexColor,id);
		if (dm_row.x!=MARKER)
		{
			while (l>=1&&remove(g[h+__mul24(l-1,size)],g[h+__mul24(l,size)],dm_row,h))
			{
				--l;
			}
			l++;g[h+__mul24(l,size)] = dm_row;
		}
	}

	int ns = l;
	if (ns==-1) return ;
	int tmpx,tmpy,tmp1,tmp2;
	short2 tmpp,tmpg;
	l =0 ;
	for (int i=b;i<=e;i++)
	{
		tmpg = g[h+__mul24(l,size)];
		tmpy = tmpg.y - i; 
		tmpx = tmpg.x - h;
		tmp1 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
		tmpp = tmpg;
		while(1)
		{
			if (l>=ns) break;
			tmpg = g[h+__mul24(l+1,size)];
			tmpy = tmpg.y - i;
			tmpx = tmpg.x - h;
			tmp2 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
			if (tmp1<=tmp2) break;
			++l;
			tmp1 = tmp2;
			tmpp = tmpg;
		}
		out[h+__mul24(i,size)]=tmpp;
	}
}

__global__ void maKernelColorInit(short2 *out,short2 *g,int size,int kb,int ks)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int idx = kb + tx*ks;
	if (idx>=size) return ;
	maColor(out,g,size,idx,0,0,size-1,size-1);
}

__global__ void maKernelColorLine(short2 *out,short2 *g,int size,int kb,int ks)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int idx = kb + tx*ks;
	if (idx>=size) return ;
	int b = MARKER,by=MARKER;
	short prey;
	for (int i=0;i<size;i++)
	{
		short2 lf = tex1Dfetch(maTexLinks,TOID(idx-kb,i,size));
		short2 rg = tex1Dfetch(maTexLinks,TOID(idx+kb,i,size));
		if (lf.x==rg.x&&lf.y==rg.y)
		{
			out[idx + __mul24(i,size)] = lf;
			if (b!=MARKER)
			{
				maColor(out,g,size,idx,b,by,i-1,lf.y);
				b=by=MARKER;
			}
		}
		else
		{
			if (b==MARKER)
			{
				b=i;
				if (i==0) by=0;
				else by = prey;
			}
		}
		prey=lf.y;
	}
	if (b!=MARKER) maColor(out,g,size,idx,b,by,size-1,size-1);
}

__global__ void maKernelTest1(short2 *g,int *lx,int size)
{
	int h = blockDim.x * blockIdx.x + threadIdx.x;
	int l=-1;
	for (int i=0;i<size;i++)
	{
		int id = TOID(h,i,size);
		short2 dm_row = tex1Dfetch(maTexColor,id);
		if (dm_row.x!=MARKER)
		{
			while (l>=1&&remove(g[TOID(h,l-1,size)],g[TOID(h,l,size)],dm_row,h))
			{
				--l;
			}
			l++;g[TOID(h,l,size)] = dm_row;
		}
	}
	lx[h]=l;
}

__global__ void maKernelTest2(short2 *out ,int size)
{
	int h = blockDim.x * blockIdx.x + threadIdx.x;

	int ns = tex1Dfetch(maTexIndex,h);
	if (ns==-1) return ;
	int tmpx,tmpy,tmp1,tmp2;
	short2 tmpp,tmpg;
	int l =0 ;
	for (int i=0;i<size;i++)
	{
		tmpg = tex1Dfetch(maTexLinks,h+__mul24(l,size));
		//tmpg = links[h+__mul24(l,size)];
		tmpy = tmpg.y - i; 
		tmpx = tmpg.x - h;
		tmp1 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
		tmpp = tmpg;
		while(1)
		{
			if (l>=ns) break;
			tmpg = tex1Dfetch(maTexLinks,h+__mul24(l+1,size));	
			//tmpg = links[h+__mul24(l+1,size)];
			tmpy = tmpg.y - i;
			tmpx = tmpg.x - h;
			tmp2 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
			if (tmp1<=tmp2) break;
			++l;
			tmp1 = tmp2;
			tmpp = tmpg;
		}
		out[TOID(h,i,size)]=tmpp;
	}
}


__global__ void maKernelColorInit1(short2 *g,int *lx,int size,int kb,int ks)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int h = tx*ks+kb;
	if (h>=size) return ;
	int h_index = tex1Dfetch(maTexIndex,h);
	//h_index = h;
	int l=-1;
	short2 lst2=make_short2(MARKER,MARKER),lst1=make_short2(MARKER,MARKER);
	for (int i=0;i<size;i++)
	{
		int id = TOID(h,i,size);
		short2 dm_row = tex1Dfetch(maTexColor,id);
		if (dm_row.x!=MARKER)
		{
			//while (l>=1&&remove(g[TOID(h_index,l-1,size)],g[TOID(h_index,l,size)],dm_row,h))
			while (l>=1&&remove(lst2,lst1,dm_row,h))
			{
				--l;
				lst1 = lst2;
				if (i>=1) lst2 = g[TOID(h_index,l-1,size)];
			}
			l++;
			g[TOID(h_index,l,size)] = dm_row;
			lst2 = lst1;lst1 = dm_row;
		}
	}
	lx[h_index]=l;
}

__global__ void maKernelColorInit2(short2 *out ,int size,int kb,int ks)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int h = tx*ks+kb;
	if (h>=size) return ;
	int h_index = tex1Dfetch(maTexIndex,h);
	//h_index = h;
	int ns = tex1Dfetch(maTexNpnt,h_index);
	if (ns==-1) return ;
	int tmpx,tmpy,tmp1,tmp2;
	short2 tmpp,tmpg;
	int l =0 ;
	for (int i=0;i<size;i++)
	{
		tmpg = tex1Dfetch(maTexColor,h_index+__mul24(l,size));
		//tmpg = links[h+__mul24(l,size)];
		tmpy = tmpg.y - i; 
		tmpx = tmpg.x - h;
		tmp1 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
		tmpp = tmpg;
		while(1)
		{
			if (l>=ns) break;
			tmpg = tex1Dfetch(maTexColor,h_index+__mul24(l+1,size));	
			//tmpg = links[h+__mul24(l+1,size)];
			tmpy = tmpg.y - i;
			tmpx = tmpg.x - h;
			tmp2 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
			if (tmp1<=tmp2) break;
			++l;
			tmp1 = tmp2;
			tmpp = tmpg;
		}
		out[TOID(h_index,i,size)]=tmpp;
	}
}
__global__ void maKernelColor1(short2 *out,
							   short2 *pntOut,int *pntl,
							   short2 *npntOut,int *npntl,
							   int size,int kb,int ks)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int idx = kb + tx*ks;
	if (idx>=size) return ;
	int idx_index  = tex1Dfetch(maTexIndex,idx);
	int b = -1;
	int tmpPntl = 0;
	for (int i=0;i<size;i++)
	{
		int lfidx = tex1Dfetch(maTexIndex,idx-kb);
		int rgidx = tex1Dfetch(maTexIndex,idx+kb);
		short2 lf = tex1Dfetch(maTexLinks,lfidx + __mul24(i,size));
		short2 rg = tex1Dfetch(maTexLinks,rgidx + __mul24(i,size));
		if (lf.x==rg.x&&lf.y==rg.y)
		{
			out[idx_index + __mul24(i,size)] = lf;
			if (b!=-1)
			{
				pntOut[idx_index + __mul24(tmpPntl,size)] = make_short2(b,i-1);
				tmpPntl++;
				b=-1;
			}
		}
		else
		{
			if (b==-1)
			{
				b=i;
			}
		}
	}
	if (b!=-1)
	{
		pntOut[idx_index + __mul24(tmpPntl,size)] = make_short2(b,size-1);
		tmpPntl++;
	}
	pntl[idx_index] = tmpPntl;
}

__global__ void maKernelColor2(short2 *out,int *npntl,int size,int kb,int ks)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int h = tx*ks+kb;
	if (h>=size) return;
	int h_index = tex1Dfetch(maTexIndex,h);
	int ol = tex1Dfetch(maTexNpnt,h_index);
	int l=-1;
	short2 lst2=make_short2(MARKER,MARKER),lst1=make_short2(MARKER,MARKER);
	for (int i=0;i<size;i++)
	{
		int id = h+__mul24(i,size);
		short2 dm_row = tex1Dfetch(maTexColor,id);
		if (dm_row.x!=MARKER)
		{
			while (l>=1&&remove(lst2,lst1,dm_row,h))
			{
				--l;
				lst1 = lst2;
				if (i>=1) lst2 = out[h_index + __mul24(l-1,size)];
			}
			l++;
			out[h_index + __mul24(l,size)] = dm_row;
			lst2 = lst1;lst1 = dm_row;
		}
	}
	npntl[h_index]=l;
}

__global__ void maKernelColor3(short2 *out,int size,int kb,int ks)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int h = tx*ks+kb;
	if (h>=size) return ;
	int h_index = tex1Dfetch(maTexIndex,h);
	int ns = tex1Dfetch(maTexNpnt,h_index);
	if (ns==-1) return ;
	int tmpx,tmpy,tmp1,tmp2;
	short2 tmpp,tmpg,now;
	int l =0 ;
	int pntNum = tex1Dfetch(maTexPnt,h_index);
	for (int k=0;k<pntNum;k++)
	{
		//int k=0;
		short2 tmp = tex1Dfetch(maTexLinks,h_index + __mul24(k,size));
		//if (pntNum>0) tmp = tex1Dfetch(maTexLinks,h_index + __mul24(k,size));
		int b = tmp.x;
		int e = tmp.y;
		//int i=b;
		for (int i=b;i<=e;i++)
		//while (k<pntNum)
		{
			tmpg = tex1Dfetch(maTexColor,h_index + __mul24(l,size));
			tmpy = tmpg.y - i; 
			tmpx = tmpg.x - h;
			tmp1 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
			tmpp = tmpg;
			while(1)
			{
				if (l>=ns) break;
				tmpg = tex1Dfetch(maTexColor,h_index + __mul24(l+1,size));	
				tmpy = tmpg.y - i;
				tmpx = tmpg.x - h;
				tmp2 = __mul24(tmpx,tmpx) + __mul24(tmpy,tmpy);
				if (tmp1<=tmp2) break;
				++l;
				tmp1 = tmp2;
				tmpp = tmpg;
			}
			out[TOID(h_index,i,size)]=tmpp;
			/*
			i++;
			if (i==e)
			{
				k++;
				if (k==pntNum) break;
				tmp = tex1Dfetch(maTexLinks,h_index + __mul24(k,size));
				b = tmp.x;
				e = tmp.y;
				i=b;
			}*/
		}
	}
}

#endif