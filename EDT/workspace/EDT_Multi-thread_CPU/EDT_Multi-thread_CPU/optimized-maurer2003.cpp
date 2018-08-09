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
	int   rows=dm->rows, r,
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

int opt_edt_maurer_2D_from_1D(DistMap *dm,DistMap *dmo,int step);

int opt_edt_maurer2003(DistMap * dm,DistMap *dmo,int step)
{
	char *fname="edt_maurer2003";
	int i,j,r,c;
	int stat;
	int infty;

	//assert(im->isbinary);
	r = dm->rows; c = dm->cols;
	infty = MAXINT - r*r - c*c;
	// Vertical columnwise EDT
	stat = opt_edt_1d_vertical(dm);
	// Lotufo's 1D EDT is equivalent to Maurer's D1.
	// There is a remark in section 5 of Maurer's paper that says
	// D1 can be calculated using a routine like this.
	//CHECK_RET_STATUS(false);
	stat = opt_edt_maurer_2D_from_1D(dm,dmo,step);    
	//CHECK_RET_STATUS(false);

	/*for (i=0; i < r; i++)
	{
		for (j=0;j<c;j++)
		{
			Pnt tmp = RC(dm,i,j);
			int tmpx = tmp.x - i;
			int tmpy = tmp.y - j;
			RC(im,i,j) = tmpx*tmpx + tmpy*tmpy;
		}
	}*/
	//im->isbinary = false;
	return 0;
}

#ifdef ST

int opt_remove_edt(Pnt *mlst2,Pnt *mlst1,Pnt *now,int h,float *y1,float *y2);
float opt_interpnt(int x1,int y1,int x2,int y2,int x0);
int opt_maurer_voronoi_edt_2D(DistMap *dm,Pnt *dmo_row,Pnt *im_row, Pnt *g, int h,int b,int by,int e,int ey);

Pnt ppp[5000];
float yyy[5000];
Pnt icp[5000][5000];
int icpl[5000];


int opt_maurer_voronoi_edt_2D_init(DistMap *dm,DistMap *dmo,Pnt *g,int step)
{
	int r = dm->rows;
	int c = dm->cols;
	int i,l,ynum,ns,yl;
	int infty = MAXINT - r*r - c*c, fi;
	Pnt tmpp;
	for (int k=0;k<c;k+=step)
	{
		//opt_maurer_voronoi_edt_2D(dm,dmo->data+k*c,dm->data+k*c,g,k,0,0,c-1,c-1);

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

int getYb(int h,int bound)
{
	int b=0,e=icpl[h];
	while (b<e)
	{
		int mid=(b+e)/2;
		if (icp[h][mid].y<bound) b=mid+1;
		else e=mid;
	}
	return b;
	/*for (int i=0;i<=icpl[h];i++)
	{
		if (icp[h][i].y>=bound) return i;
	}
	return -1;*/
}

int opt_edt_maurer_2D_from_1D(DistMap *dm,DistMap *dmo,int step)
{	
	int j, ns,fi,e,ey,yl,ynum,l,pl;
	int c = dm->cols;
	int r = dm->rows;
	for (int i=0;i<r;i++) 
	{
		for (j=0;j<c;j++) icp[i][j].x=icp[i][j].y=-MARKER;
		icpl[i]=-1;
	}
	Pnt *g = (Pnt*)malloc(sizeof(Pnt)*c);
	Pnt *dm_row,*dmo_row,tmpp;
	opt_maurer_voronoi_edt_2D_init(dm,dmo,g,step);
	while (step>=2)
	{
		int step2 = step/2;
		for (int k=step2;k<r;k+=step)
		{
			int b=MARKER,by=MARKER;
			for (int i=0;i<c;i++)
			{
				//printf("%d %d\n",k,i);
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
						int ul = k-step2;
						int dl = k+step2;
						int uyb = getYb(ul,by);
						int dyb = getYb(dl,by);
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
						if (dm_row[ey].x!=MARKER) ppp[++pl]=dm_row[ey];
						
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
				int ul = k-step2;
				int dl = k+step2;
				int uyb = getYb(ul,by);
				int dyb = getYb(dl,by);
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
				if (dm_row[ey].x!=MARKER) ppp[++pl]=dm_row[ey];

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

	free(g); 
	return 1;
}

#else

struct Phrase2
{
	DistMap *dm;
	int k,h,b,by,e,ey;
	Pnt *dmo_row,*im_row,*g;
};

int opt_maurer_voronoi_edt_2D(DistMap *dm,Pnt *dmo_row,Pnt *im_row, Pnt *g, int h,int b,int by,int e,int ey);


DWORD WINAPI opt_maurer_voronoi_edt_2D_init_kernel(LPVOID lpParam)
{
	Phrase2 *tmp = (Phrase2 *)lpParam;
	DistMap *dm = tmp->dm;
	Pnt *dmo_row = tmp->dmo_row;
	Pnt *im_row = tmp->im_row;
	Pnt *g = tmp->g;
	int h = tmp->h;
	int b = tmp->b;
	int by = tmp->by;
	int e = tmp->e;
	int ey = tmp->ey;
	opt_maurer_voronoi_edt_2D(dm,dmo_row,im_row,g,h,b,by,e,ey);
	return 0;
}


int opt_maurer_voronoi_edt_2D_init(DistMap *dm,DistMap *dmo,Pnt *g,int step)
{
	int r = dm->rows;
	int c = dm->cols;
	HANDLE hdl[THREADNUM];
	Phrase2 tmp[THREADNUM];
	int num=0;
	for (int k=0;k<c;k+=step)
	{
		tmp[num].dm = dm;
		tmp[num].dmo_row = dmo->data+ k*c;
		tmp[num].im_row = dm->data + k*c;
		tmp[num].g = g + num*c;
		tmp[num].h = k;
		tmp[num].b = 0;
		tmp[num].by = 0;
		tmp[num].e = c-1;
		tmp[num].ey = c-1;

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
	DistMap *dm,*dmo;
	int k,step2;
	Pnt *g;
};

DWORD WINAPI opt_edt_maurer_2D_from_1D_kernel(LPVOID lpParam)
{
	Phrase3 *tmp = (Phrase3 *)lpParam;
	DistMap *dm = tmp->dm;
	DistMap *dmo = tmp->dmo;
	Pnt *g = tmp->g;
	int k = tmp->k;
	int step2 = tmp->step2;
	int c = dm->cols;
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
				opt_maurer_voronoi_edt_2D(dm,dmo->data+k*c,dm->data + k*c,g,k,b,by,i-1,up.y);
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
	if (b!=MARKER) opt_maurer_voronoi_edt_2D(dm,dmo->data+k*c,dm->data+k*c,g,k,b,by,c-1,c-1);
	return 0;
}


int opt_edt_maurer_2D_from_1D(DistMap *dm,DistMap *dmo,int step)
{
	int stat;
	int i1; // same naming as in the paper
	Pnt *g;
	char *fname="edt_maurer_2D_from_1D";

	int c = dm->cols;
	int r = dm->rows;

	g = (Pnt*)malloc(sizeof(Pnt)*c*THREADNUM);

	opt_maurer_voronoi_edt_2D_init(dm,dmo,g,step);

/**/
	HANDLE hdl[THREADNUM];
	Phrase3 tmp[THREADNUM];

	Pnt *dm_row;
	dm_row = dm->data;
	while (step>=2)
	{
		int step2 = step/2;
		int num=0;
		for (int k=step2;k<r;k+=step)
		{
			tmp[num].dm = dm;
			tmp[num].dmo = dmo;
			tmp[num].g = g + num*c;
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
	}/**/
	/*for (i1=0; i1 < nrows; ++i1, dm_row += ncols ) {
		stat = opt_maurer_voronoi_edt_2D(dm, dm_row,   g,i1);
		//CHECK_RET_STATUS(false);
	}*/
	free(g); //free(h);

	return 1;
}

#endif


int opt_maurer_voronoi_edt_2D(DistMap *dm,Pnt *dmo_row ,Pnt *dm_row, Pnt *g,int h,int b,int by,int e,int ey)
{
	int l, i, ns, tmpx,tmpy, tmp1, tmp2, r,c;
	int infty, fi;
	Pnt tmpp;

	r = dm->rows; c=dm->cols;
	infty = MAXINT - r*r - c*c;
	int ynum=-1;
	l = -1;
	//printf("%d %d %d %d\n",b,e,by,ey);
	/*for (i=by;i<b;i++)
	{
		if (l>=0&&g[l].x==dmo_row[i].x&&g[l].y==dmo_row[i].y) continue;
		++l;g[l]=dmo_row[i];
		if (l>0)
		{
			ynum++;
			yyy[ynum]=opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,h);
		}
	}*/
	for (i=by; i <= ey; ++i){
		if ((fi = dm_row[i].x) != MARKER) {
			Pnt now;
			now.x = fi;
			now.y = i;
			float y1=-1,y2=-1;
			if (l>=1)
			{
				y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,h);
				y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,h);
			}
			//while ( l >= 1 && opt_remove_edt(&g[l-1], &g[l], &now , h,&y1,&y2) )
			while (l>=1&& y1>y2)
			{
				--l;
				--ynum;
				if (l>=1)
				{
					y1 = opt_interpnt(g[l-1].x,g[l-1].y,g[l].x,g[l].y,h);
					y2 = opt_interpnt(g[l].x,g[l].y,now.x,now.y,h); 
				}
			}
			++l; g[l].x = fi; g[l].y = i;
			if (l>0)
			{
				if (l==1) yyy[++ynum] = opt_interpnt(g[0].x,g[0].y,g[1].x,g[1].y,h);
				else yyy[++ynum]=y2;
			}
		}
	}

	
	if ((ns=l) == -1) return 1;
	l = 0;
	int yl=0;
	for (i=b; i <= e; ++i) {
		/*
		tmpy = g[l].y - i;
		tmpx = g[l].x-h;
		tmp1 = tmpx*tmpx + tmpy*tmpy;
		tmpp = g[l];
		while(1) {
			if (l >= ns) break;

			tmpy = g[l+1].y - i;
			tmpx = g[l+1].x - h;
			tmp2 = tmpx*tmpx + tmpy*tmpy; 
			if (tmp1 <= tmp2) break;
			++l;
			tmp1 = tmp2;
			tmpp = g[l];
		}*/
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

	return 1;
}




inline float opt_interpnt(int x1,int y1,int x2,int y2,int x0)
{
	float xM = (float)(x1 + x2)/2.0f;
	float yM = (float)(y1 + y2)/2.0f;
	float nx = x2 - x1;
	float ny = y2 - y1;
	return yM + nx*(xM-x0)/ny;
}


inline int opt_remove_edt(Pnt *mlst2,Pnt *mlst1,Pnt *now,int h,float *y1,float *y2)
{
	/*
	int d2 = mlst2->x - h;
	int d1 = mlst1->x - h;
	int dn = now->x - h;
	if (d2>=0 && d1>=0 && dn>=0 && d1<=d2 && d1<=dn) return 0;
	if (d2<=0 && d1<=0 && dn<=0 && d1>=d2 && d1>=dn) return 0;
	*/
	/*float x1=mlst2->x,y11=mlst2->y,x2=mlst1->x,y22=mlst1->y,x0=h;
	float xM = (x1 + x2)/2.0f;
	float yM = (y11 + y22)/2.0f;
	float nx = x2 - x1;
	float ny = y22 - y11;
	float i1 =  yM + nx*(xM-x0)/ny;


	x1 = mlst1->x;
	y11 = mlst1->y;
	x2 = now->x;
	y22=now->y;
	xM = (float)(x1 + x2)/2.0f;
	yM = (float)(y11 + y22)/2.0f;
	nx = x2 - x1;
	ny = y22 - y11;
	float i2 = yM + nx*(xM-x0)/ny;
	*/

	float i1 = opt_interpnt(mlst2->x,mlst2->y,mlst1->x,mlst1->y,h);
	float i2 = opt_interpnt(mlst1->x,mlst1->y,now->x,now->y,h);
	*y1 = i1;
	*y2 = i2;
	return (i1>i2);
}

