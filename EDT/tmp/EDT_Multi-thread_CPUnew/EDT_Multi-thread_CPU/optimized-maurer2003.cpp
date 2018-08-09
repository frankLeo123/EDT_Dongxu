/*----FUNCTION-----------------------------------------------------------------
 * 
 *  Description:  Calvin Maurer's EDT (Euclidean DistanceTransform): 
 *
 * REFERENCE
 *  Paper: Calvin Maurer et. al. PAMI feb. 2003 
 *
 *  Implementation "inspired" in the LTI lib:
 *    http://ltilib.sourceforge.net
 *  
 *  Distances are calculated only where pixels are equal to 1
 *
 * - PARAMETER ----------------------------------------------------------------
 *      Mode   Type             Name            Description
 * ----------------------------------------------------------------------------
 *    in-out:  ImgPUInt32         im            binary image on input;
 *                                              grayscale on output, with
 *                                              squared distance transform.
 *----------------------------------------------------------------------------*/

#include "optimized-maurer2003.h"
#include "common.h"


int getYb(Pnt* icp,int icpl,int bound)
{
	int b=0,e=icpl;
	while (b<e)
	{
		int mid=(b+e)/2;
		if (icp[mid].y<bound) b=mid+1;
		else e=mid;
	}
	return b;
}
//y1<y2
inline float opt_interpnt(int x1,int y1,int x2,int y2,int x0)
{
	float xM = (float)(x1 + x2)/2.0f;
	float yM = (float)(y1 + y2)/2.0f;
	float nx = x2 - x1;
	float ny = y2 - y1;
	return yM + nx*(xM-x0)/ny;
}

#ifdef ST

int opt_edt_1d_vertical(DistMap *dm)
{
	int   rows=dm->rows, r,
		  cols=dm->cols, c;
	int b;
	for (c=0; c < cols; c++) {
		Pnt pre = RC(dm,0,c);
		for (r=1; r<rows; r++)
		{
			if (RC(dm,r,c).x==MARKER) RC(dm,r,c) = pre;
			else pre = RC(dm,r,c);
		}
		for (r=rows-2; r >= 0; r--) 
		{
			if (RC(dm,r,c).x==MARKER) RC(dm,r,c) = RC(dm,r+1,c);
			else
			{
				Pnt now;
				now.x = r;now.y=c;
				int dis1 = getDist(&now,&RC(dm,r,c));
				int dis2 = getDist(&now,&RC(dm,r+1,c));
				if (dis2<dis1) RC(dm,r,c) = RC(dm,r+1,c);
			}
		}
	}
	return 0;
}

#else

struct Phrase1
{
	DistMap *dm;
	int c;
};

DWORD WINAPI opt_edt_1d_vertical_kernel(LPVOID lpParam)
{
	Phrase1 *tmp = (Phrase1 *)lpParam;
	DistMap *dm = tmp->dm;
	int c = tmp->c;
	int rows = dm->rows;
	Pnt pre = RC(dm,0,c);
	int r;
	for (r=1; r<rows; r++)
	{
		if (RC(dm,r,c).x==MARKER) RC(dm,r,c) = pre;
		else pre = RC(dm,r,c);
	}
	for (r=rows-2; r >= 0; r--) 
	{
		if (RC(dm,r,c).x==MARKER) RC(dm,r,c) = RC(dm,r+1,c);
		else
		{
			Pnt now;
			now.x = r;now.y=c;
			int dis1 = getDist(&now,&RC(dm,r,c));
			int dis2 = getDist(&now,&RC(dm,r+1,c));
			if (dis2<dis1) RC(dm,r,c) = RC(dm,r+1,c);
		}
	}
	return 0;
}

int opt_edt_1d_vertical(DistMap *dm)
{
	int   rows=dm->rows,
		  cols=dm->cols, c;
	Phrase1 tmp[THREADNUM];
	HANDLE handle[THREADNUM];
	int num=0;
	for (c=0; c < cols; c++) {
		
		tmp[num].dm = dm;
		tmp[num].c = c;
		handle[num]=CreateThread(NULL,0,opt_edt_1d_vertical_kernel,&tmp[num],0,NULL);
		num++;
		if (num==THREADNUM||c==(cols-1))
		{
			num=0;
			WaitForMultipleObjects(THREADNUM,handle,true,INFINITE);
		}
	}
	for (int i=0;i<THREADNUM;i++) CloseHandle(handle[i]);
	return 0;
}
#endif


#ifdef ST


int opt_maurer_voronoi_edt_2D_init(DistMap *dm,DistMap *dmo,Pnt **icp,int *icpl,float *yyy,Pnt *g,int step)
{
	int r = dm->rows;
	int c = dm->cols;
	int i,l,ynum,ns,yl,fi;
	Pnt tmpp;
	for (int k=0;k<c;k+=step)
	{
		ynum=-1;
		l = -1;
		Pnt *dm_row=dm->data+k*c;
		Pnt *dmo_row = dmo->data+k*c;
		for (i=0; i <= c-1; ++i){
			if ((fi = dm_row[i].x) != MARKER) {
				Pnt now;
				now.x = fi;
				now.y = i;
				float y1=-1,y2=-1;
				if (l>=1)
				{
					y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
					y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
				}
				while (l>=1&& y1>y2)
				{
					--l;
					--ynum;
					if (l>=1)
					{
						y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
						y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
					}
				}
				++l; g[l].x = fi; g[l].y = i;
				if (l>0)
				{
					if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
					else yyy[++ynum]=y2;
				}
			}
		}
		for (i=0;i<=l;i++) icp[k][i]=g[i];
		icpl[k]=l;
		if ((ns=l) == -1) return 1;
		l = yl=0;
		for (i=0; i <= c-1; ++i) {
			tmpp = g[l];
			while (i>yyy[yl])
			{
				if (l>=ns) break;
				yl++;
				l++;
				tmpp=g[l];
			}

			dmo_row[i] = tmpp;
		}
	}
	return 0;
}

int cmp(const void *a,const void *b)
{
	Pnt *c = (Pnt *)a;
	Pnt *d = (Pnt *)b;
	return   (c->y) - (d->y)  ;
}

int opt_edt_maurer_2D_from_1D(DistMap *dm,DistMap *dmo,Pnt **icp,int *icpl,int step)
{	
	int j, ns,fi,e,ey,yl,ynum,l,pl;
	int c = dm->cols;
	int r = dm->rows;
	Pnt *g = (Pnt*)malloc(sizeof(Pnt)*c);
	Pnt *ppp = (Pnt*)malloc(sizeof(Pnt)*c);
	float *yyy= (float*)malloc(sizeof(float)*c);
	int *xxx= (int*)malloc(sizeof(int)*c);
	
	Pnt *dm_row,*dmo_row,tmpp;
	opt_maurer_voronoi_edt_2D_init(dm,dmo,icp,icpl,yyy,g,step);
	while (step>=2)
	{
		int step2 = step/2;
		for (int k=step2;k<r;k+=step)
		{

			//for (j=0;j<c;j++) RC(dmo,k,j).x = RC(dmo,k,j).y = MARKER;
			for (j=0;j<c;j++) xxx[j]=0;		
			int ul = k-step2;
			int dl = k+step2;
			int uyb = 0;
			int dyb = 0;
			int ticpl=-1;
			while (uyb<=icpl[ul]||dyb<=icpl[dl])
			{
				if (icp[ul][uyb].y<icp[dl][dyb].y)
				{
					uyb++;
				}
				else if (icp[ul][uyb].y>icp[dl][dyb].y)
				{
					dyb++;
				}
				else
				{
					if (icp[ul][uyb].x==icp[dl][dyb].x)
					{
						float yul1,yul2,ydl1,ydl2;
						if (uyb==0) yul1 = 0;
						else yul1 = opt_interpnt(icp[ul][uyb].x,icp[ul][uyb].y,icp[ul][uyb-1].x,icp[ul][uyb-1].y,ul);
						if (uyb>=icpl[ul]) yul2 = c-1;
						else yul2 = opt_interpnt(icp[ul][uyb].x,icp[ul][uyb].y,icp[ul][uyb+1].x,icp[ul][uyb+1].y,ul);
						if (dyb==0) ydl1 = 0;
						else ydl1 = opt_interpnt(icp[dl][dyb].x,icp[dl][dyb].y,icp[dl][dyb-1].x,icp[dl][dyb-1].y,dl);
						if (dyb>=icpl[dl]) ydl2 = c-1;
						else ydl2 = opt_interpnt(icp[dl][dyb].x,icp[dl][dyb].y,icp[dl][dyb+1].x,icp[dl][dyb+1].y,dl);
						
						float yk1 = (yul1 + ydl1)/2.0;
						float yk2 = (yul2 + ydl2)/2.0;
						int iyk1;
						int iyk2 = int(min(yk2,c-1));
						if (yk1>EPS) 
						{
							iyk1 = int(yk1+ EPS) + 1;
						}
						else iyk1 = 0;
						for (j=iyk1;j<=iyk2;j++)
						{
							RC(dmo,k,j)=icp[ul][uyb];
						}
						xxx[icp[ul][uyb].y]=1;
					}
					uyb++;dyb++;
				}
			}

			dm_row=dm->data+k*c;
			dmo_row = dmo->data+k*c;
 			int b=MARKER,by=MARKER,e,ey;
			for (int i=0;i<c;i++)
			{
				if (RC(dmo,k,i).x!=MARKER)
				{
					if (b!=MARKER)
					{
						e=i-1;
						ey=RC(dmo,k,i).y;
						
						ynum = l = pl = -1;
						for (j=by;j<=ey;j++) if (dm_row[j].x!=MARKER) ppp[++pl]=dm_row[j];
						for (j=0; j <= pl; ++j)
						{
							Pnt now;
							now = ppp[j];
							float y1=-1,y2=-1;
							if (l>=1)
							{
								y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
								y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
							}
							while (l>=1&& y1>y2)
							{
								--l;
								--ynum;
								if (l>=1)
								{
									y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
									y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
								}
							}
							++l; g[l]=now;
							if (l>0)
							{
								if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
								else yyy[++ynum]=y2;
							}

						}
						for (j=0;j<=l;j++) xxx[g[j].y]=1;
						if ((ns=l) == -1) return 1;
						l = yl=0;
						for (j=b; j <= e; ++j) {
							tmpp = g[l];
							while (j>yyy[yl])
							{
								
								if (l>=ns) break;
								yl++;
								l++;
								tmpp=g[l];
							}
							dmo_row[j] = tmpp;
						}

						b=MARKER;by=MARKER;
					}
				}
				else
				{
					if (b==MARKER)
					{
						b=i;
						if (i==0) by=0;
						else by=RC(dmo,k,i-1).y;
					}
				}
			}
			if (b!=MARKER)
			{
				dm_row=dm->data+k*c;
				dmo_row = dmo->data+k*c;
				e=c-1;ey=c-1;
				ynum = l = pl = -1;
				for (j=by;j<=ey;j++) if (dm_row[j].x!=MARKER) ppp[++pl]=dm_row[j];
				for (j=0; j <= pl; ++j)
				{
					Pnt now;
					now = ppp[j];
					float y1=-1,y2=-1;
					if (l>=1)
					{
						y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
						y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
					}
					while (l>=1&& y1>y2)
					{
						--l;
						--ynum;
						if (l>=1)
						{
							y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
							y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
						}
					}
					++l; g[l]=now;
					if (l>0)
					{
						if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
						else yyy[++ynum]=y2;
					}

				}
				for (j=0;j<=l;j++) xxx[g[j].y]=1;
				if ((ns=l) == -1) return 1;
				l = yl=0;
				for (j=b; j <= e; ++j) {
					tmpp = g[l];
					while (j>yyy[yl])
					{
						if (l>=ns) break;
						yl++;
						l++;
						tmpp=g[l];
					}
					dmo_row[j] = tmpp;
				}
			}
			ticpl=-1;
			for (j=0;j<c;j++) if (xxx[j]==1) icp[k][++ticpl]=dm_row[j];
			icpl[k]=ticpl;
			
		}
		step=step2;
	}
	free(yyy);
	free(ppp);
	free(g); 
	return 1;
}

/*
int opt_edt_maurer_2D_from_1D(DistMap *dm,DistMap *dmo,Pnt **icp,int *icpl,int step)
{	
	int j, ns,fi,e,ey,yl,ynum,l,pl;
	int c = dm->cols;
	int r = dm->rows;
	Pnt *g = (Pnt*)malloc(sizeof(Pnt)*c);
	Pnt *ppp = (Pnt*)malloc(sizeof(Pnt)*c);
	float *yyy= (float*)malloc(sizeof(float)*c);
	Pnt *dm_row,*dmo_row,tmpp;
	opt_maurer_voronoi_edt_2D_init(dm,dmo,icp,icpl,yyy,g,step);
	while (step>=2)
	{
		int step2 = step/2;
		for (int k=step2;k<r;k+=step)
		{

			int b=MARKER,by=MARKER;
			for (int i=0;i<c;i++)
			{
				
 				Pnt up = RC(dmo,k-step2,i);
				Pnt down = RC(dmo,k+step2,i);
				if (up.x==down.x&&up.y==down.y) 
				{
					RC(dmo,k,i)=up;
					if (b!=MARKER)
					{
						dm_row=dm->data+k*c;
						dmo_row = dmo->data+k*c;
						e=i-1,ey=up.y;
						ynum=l=pl=-1;
						if (by<=b&&e<=ey)
						{
							int ul = k-step2;
							int dl = k+step2;
							int uyb = getYb(icp[ul],icpl[ul],by);
							int dyb = getYb(icp[dl],icpl[dl],by);
							while (icp[ul][uyb].y<b||icp[dl][dyb].y<b)
							{
								if (icp[ul][uyb].y<icp[dl][dyb].y)
								{
									ppp[++pl]=dm_row[icp[ul][uyb].y];
									uyb++;
								}
								else if (icp[ul][uyb].y>icp[dl][dyb].y)
								{
									ppp[++pl]=dm_row[icp[dl][dyb].y];
									dyb++;
								}
								else
								{
									ppp[++pl]=dm_row[icp[dl][dyb].y];
									uyb++;dyb++;
								}
							}
							for (j=b;j<=e;j++)
							{
								if (dm_row[j].x!=MARKER) ppp[++pl]=dm_row[j];
							}
							uyb = getYb(icp[ul],icpl[ul],e+1);
							dyb = getYb(icp[dl],icpl[dl],e+1);
							while (icp[ul][uyb].y<ey||icp[dl][dyb].y<ey)
							{
								if (icp[ul][uyb].y<icp[dl][dyb].y)
								{
									ppp[++pl]=dm_row[icp[ul][uyb].y];
									uyb++;
								}
								else if (icp[ul][uyb].y>icp[dl][dyb].y)
								{
									ppp[++pl]=dm_row[icp[dl][dyb].y];
									dyb++;
								}
								else
								{
									ppp[++pl]=dm_row[icp[dl][dyb].y];
									uyb++;dyb++;
								}
							}
							if (e<ey&&dm_row[ey].x!=MARKER) ppp[++pl]=dm_row[ey];
						}
						else
						{
							for (j=by;j<=ey;j++) if (dm_row[j].x!=MARKER) ppp[++pl]=dm_row[j];
						}
						
						for (j=0; j <= pl; ++j)
						{
							Pnt now;
							now.x = ppp[j].x;
							now.y = ppp[j].y;
							float y1=-1,y2=-1;
							if (l>=1)
							{
								y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
								y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
							}
							while (l>=1&& y1>y2)
							{
								--l;
								--ynum;
								if (l>=1)
								{
									y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
									y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
								}
							}
							++l; g[l]=now;
							if (l>0)
							{
								if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
								else yyy[++ynum]=y2;
							}
							
						}
						int jn = 0;
						for(j=0;j<=l;j++)
						{
							if (icpl[k]!=-1&&g[j].x==icp[k][icpl[k]+jn].x&&g[j].y==icp[k][icpl[k]+jn].y) continue;
							icp[k][icpl[k]+1+jn]=g[j];
							jn++;
						}
						icpl[k]+=jn;
						if ((ns=l) == -1) return 1;
						l = yl =0;
						for (j=b; j <= e; ++j) {
							tmpp = g[l];
							while (j>yyy[yl])
							{
								if (l>=ns) break;
								yl++;
								l++;
								tmpp=g[l];
							}

							dmo_row[j] = tmpp;
						}

						b=by=MARKER;
					}
				}
				else
				{
					if (b==MARKER)
					{
						b=i;
						if (i==0) by=0;
						else by=RC(dmo,k,i-1).y;
					}
				}
			}
			if (b!=MARKER) 
			{
				dm_row=dm->data+k*c;
				dmo_row = dmo->data+k*c;
				e=c-1;ey=c-1;
				ynum = l = pl = -1;
				for (j=by;j<=ey;j++) if (dm_row[j].x!=MARKER) ppp[++pl]=dm_row[j];
				for (j=0; j <= pl; ++j)
				{
					Pnt now;
					now = ppp[j];
					float y1=-1,y2=-1;
					if (l>=1)
					{
						y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
						y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
					}
					while (l>=1&& y1>y2)
					{
						--l;
						--ynum;
						if (l>=1)
						{
							y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
							y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
						}
					}
					++l; g[l]=now;
					if (l>0)
					{
						if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
						else yyy[++ynum]=y2;
					}
					
				}
				int jn = 0;
				for(j=0;j<=l;j++) 
				{
					if (icpl[k]!=-1&&g[j].x==icp[k][icpl[k]+jn].x&&g[j].y==icp[k][icpl[k]+jn].y) continue;
					icp[k][icpl[k]+1+jn]=g[j];
					jn++;
				}
				icpl[k]+=jn;
				if ((ns=l) == -1) return 1;
				l = yl=0;
				for (j=b; j <= e; ++j) {
					tmpp = g[l];
					while (j>yyy[yl])
					{
						if (l>=ns) break;
						yl++;
						l++;
						tmpp=g[l];
					}
					dmo_row[j] = tmpp;
				}
			}
		}
		step=step2;
	}
	free(yyy);
	free(ppp);
	free(g); 
	return 1;
}
*/


#else

DistMap *dm,*dmo;
Pnt **icp;
int *icpl;

struct Phrase2
{
	float *yyy;
	Pnt *g;
	int step,k;
};


DWORD WINAPI opt_maurer_voronoi_edt_2D_init_kernel(LPVOID lpParam)
{
	Phrase2 *tmp = (Phrase2 *)lpParam;
	float *yyy = tmp->yyy;
	Pnt *g = tmp->g;
	int step = tmp->step;
	int k = tmp->k;
	//opt_maurer_voronoi_edt_2D(dm,dmo_row,im_row,g,h,b,by,e,ey);

	int r = dm->rows,c = dm->rows,fi,ns,yl,i;
	Pnt tmpp;
	int ynum=-1;
	int l = -1;
	Pnt *dm_row=dm->data+k*c;
	Pnt *dmo_row = dmo->data+k*c;
	for (i=0; i <= c-1; ++i){
		if ((fi = dm_row[i].x) != MARKER) {
			Pnt now;
			now.x = fi;
			now.y = i;
			float y1=-1,y2=-1;
			if (l>=1)
			{
				y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
				y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
			}
			while (l>=1&& y1>y2)
			{
				--l;
				--ynum;
				if (l>=1)
				{
					y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
					y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
				}
			}
			++l; g[l].x = fi; g[l].y = i;
			if (l>0)
			{
				if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
				else yyy[++ynum]=y2;
			}
		}
	}
	for (i=0;i<=l;i++) icp[k][i]=g[i];
	icpl[k]=l;
	if ((ns=l) == -1) return 1;
	l = yl=0;
	for (i=0; i <= c-1; ++i) {
		tmpp = g[l];
		while (i>yyy[yl])
		{
			if (l>=ns) break;
			yl++;
			l++;
			tmpp=g[l];
		}

		dmo_row[i] = tmpp;
	}

	return 0;
}


int opt_maurer_voronoi_edt_2D_init(DistMap *dm,DistMap *dmo,Pnt **icp,int *icpl,float *yyy,Pnt *g,int step)
{
	int r = dm->rows;
	int c = dm->cols;
	HANDLE hdl[THREADNUM];
	Phrase2 tmp[THREADNUM];
	int num=0;
	for (int k=0;k<c;k+=step)
	{
		tmp[num].g = g + num*c;
		tmp[num].yyy = yyy + num*c;
		tmp[num].k = k;
		tmp[num].step = step;

		hdl[num]=CreateThread(NULL,0,opt_maurer_voronoi_edt_2D_init_kernel,&tmp[num],0,NULL);
		num++;
		if (num==THREADNUM)
		{
			num=0;
			WaitForMultipleObjects(THREADNUM,hdl,true,INFINITE);
		}
	}
	WaitForMultipleObjects(num,hdl,true,INFINITE);
	for (int i=0;i<THREADNUM;i++) CloseHandle(hdl[i]);
	return 0;
}

struct Phrase3
{
	float *yyy;
	int *xxx;
	int k,step2;
	Pnt *g;
};

DWORD WINAPI opt_edt_maurer_2D_from_1D_kernel(LPVOID lpParam)
{
	Phrase3 *tmp = (Phrase3 *)lpParam;
	int *xxx = tmp->xxx;
	Pnt *g = tmp->g;
	float *yyy = tmp->yyy;
	int k = tmp->k;
	int step2 = tmp->step2;


	int c = dm->cols;
	int r = dm->rows;

	int ynum,l,pl,j,ns,yl;
	Pnt *dm_row,*dmo_row,tmpp;

	for (j=0;j<c;j++) xxx[j]=0;		
	int ul = k-step2;
	int dl = k+step2;
	int uyb = 0;
	int dyb = 0;
	int ticpl=-1;
	while (uyb<=icpl[ul]||dyb<=icpl[dl])
	{
		if (icp[ul][uyb].y<icp[dl][dyb].y)
		{
			uyb++;
		}
		else if (icp[ul][uyb].y>icp[dl][dyb].y)
		{
			dyb++;
		}
		else
		{
			if (icp[ul][uyb].x==icp[dl][dyb].x)
			{
				float yul1,yul2,ydl1,ydl2;
				if (uyb==0) yul1 = 0;
				else yul1 = opt_interpnt(icp[ul][uyb].x,icp[ul][uyb].y,icp[ul][uyb-1].x,icp[ul][uyb-1].y,ul);
				if (uyb>=icpl[ul]) yul2 = c-1;
				else yul2 = opt_interpnt(icp[ul][uyb].x,icp[ul][uyb].y,icp[ul][uyb+1].x,icp[ul][uyb+1].y,ul);
				if (dyb==0) ydl1 = 0;
				else ydl1 = opt_interpnt(icp[dl][dyb].x,icp[dl][dyb].y,icp[dl][dyb-1].x,icp[dl][dyb-1].y,dl);
				if (dyb>=icpl[dl]) ydl2 = c-1;
				else ydl2 = opt_interpnt(icp[dl][dyb].x,icp[dl][dyb].y,icp[dl][dyb+1].x,icp[dl][dyb+1].y,dl);

				float yk1 = (yul1 + ydl1)/2.0;
				float yk2 = (yul2 + ydl2)/2.0;
				int iyk1;
				int iyk2 = int(min(yk2,c-1));
				if (yk1>EPS) 
				{
					iyk1 = int(yk1+ EPS) + 1;
				}
				else iyk1 = 0;
				for (j=iyk1;j<=iyk2;j++)
				{
					RC(dmo,k,j)=icp[ul][uyb];
				}
				xxx[icp[ul][uyb].y]=1;
			}
			uyb++;dyb++;
		}
	}

	dm_row=dm->data+k*c;
	dmo_row = dmo->data+k*c;
	int b=MARKER,by=MARKER,e,ey;
	for (int i=0;i<c;i++)
	{
		if (RC(dmo,k,i).x!=MARKER)
		{
			if (b!=MARKER)
			{
				e=i-1;
				ey=RC(dmo,k,i).y;

				ynum = l = pl = -1;
				//for (j=by;j<=ey;j++) if (dm_row[j].x!=MARKER) ppp[++pl]=dm_row[j];
				for (j=by; j <= ey; ++j)
				{
					if (dm_row[j].x!=MARKER)
					{
						Pnt now;
						now = dm_row[j];
						float y1=-1,y2=-1;
						if (l>=1)
						{
							y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
							y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
						}
						while (l>=1&& y1>y2)
						{
							--l;
							--ynum;
							if (l>=1)
							{
								y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
								y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
							}
						}
						++l; g[l]=now;
						if (l>0)
						{
							if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
							else yyy[++ynum]=y2;
						}
					}
				}
				for (j=0;j<=l;j++) xxx[g[j].y]=1;
				if ((ns=l) == -1) return 1;
				l = yl=0;
				for (j=b; j <= e; ++j) {
					tmpp = g[l];
					while (j>yyy[yl])
					{

						if (l>=ns) break;
						yl++;
						l++;
						tmpp=g[l];
					}
					dmo_row[j] = tmpp;
				}

				b=MARKER;by=MARKER;
			}
		}
		else
		{
			if (b==MARKER)
			{
				b=i;
				if (i==0) by=0;
				else by=RC(dmo,k,i-1).y;
			}
		}
	}
	if (b!=MARKER)
	{
		dm_row=dm->data+k*c;
		dmo_row = dmo->data+k*c;
		e=c-1;ey=c-1;
		ynum = l = pl = -1;
		//for (j=by;j<=ey;j++) if (dm_row[j].x!=MARKER) ppp[++pl]=dm_row[j];
		for (j=by; j <= ey; ++j)
		{
			if (dm_row[j].x!=MARKER)
			{
				Pnt now;
				now = dm_row[j];
				float y1=-1,y2=-1;
				if (l>=1)
				{
					y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
					y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k);
				}
				while (l>=1&& y1>y2)
				{
					--l;
					--ynum;
					if (l>=1)
					{
						y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,k);
						y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,k); 
					}
				}
				++l; g[l]=now;
				if (l>0)
				{
					if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,k);
					else yyy[++ynum]=y2;
				}
			}

		}
		for (j=0;j<=l;j++) xxx[g[j].y]=1;
		if ((ns=l) == -1) return 1;
		l = yl=0;
		for (j=b; j <= e; ++j) {
			tmpp = g[l];
			while (j>yyy[yl])
			{
				if (l>=ns) break;
				yl++;
				l++;
				tmpp=g[l];
			}
			dmo_row[j] = tmpp;
		}
	}
	ticpl=-1;
	for (j=0;j<c;j++) if (xxx[j]==1) icp[k][++ticpl]=dm_row[j];
	icpl[k]=ticpl;


	return 0;
}


int opt_edt_maurer_2D_from_1D(DistMap *_dm,DistMap *_dmo,Pnt **_icp,int *_icpl,int step)
{
	dm = _dm;
	dmo = _dmo;
	icp = _icp;
	icpl = _icpl;

	int c = dm->cols;
	int r = dm->rows;

	Pnt *g = (Pnt*)malloc(sizeof(Pnt)*c*THREADNUM);
	int *xxx = (int*)malloc(sizeof(int)*c*THREADNUM);
	float *yyy = (float*)malloc(sizeof(float)*c*THREADNUM);

	opt_maurer_voronoi_edt_2D_init(dm,dmo,icp,icpl,yyy,g,step);

	HANDLE hdl[THREADNUM];
	Phrase3 tmp[THREADNUM];

	while (step>=2)
	{
		int step2 = step/2;
		int num=0;
		for (int k=step2;k<r;k+=step)
		{
			tmp[num].g = g + num*c;
			tmp[num].yyy = yyy + num*c;
			tmp[num].xxx = xxx + num*c;
			tmp[num].k = k;
			tmp[num].step2 = step2;

			hdl[num] = CreateThread(NULL,0,opt_edt_maurer_2D_from_1D_kernel,&tmp[num],0,NULL);
			num++;
			if (num==THREADNUM)
			{
				num=0;
				WaitForMultipleObjects(THREADNUM,hdl,true,INFINITE);
			}
		}
		WaitForMultipleObjects(num,hdl,true,INFINITE);
		step=step2;
	}
	free(g);
	free(xxx);
	free(yyy);
	return 1;
}

#endif


int opt_edt_maurer2003(DistMap * dm,DistMap *dmo,Pnt **icp,int *icpl,int step)
{
	opt_edt_1d_vertical(dm);
	opt_edt_maurer_2D_from_1D(dm,dmo,icp,icpl,step);
	return 0;
}

