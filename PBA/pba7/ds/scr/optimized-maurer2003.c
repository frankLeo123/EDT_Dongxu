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
	/* NOTE: Lotufo's implementation (obtained by requesting him) of this first 
	 * part  is much less readable. Although pointers could be used more 
	 * efficiently, this first part is much faster than the 2nd part and is not 
	 * worth optimizing.  So I kept it readable, close to the paper's pseudocode.  
	 */

	return 0;
}

int opt_edt_maurer_2D_from_1D(DistMap *dm,DistMap *dmo);

int opt_edt_maurer2003(DistMap * dm,DistMap *dmo)
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
	stat = opt_edt_maurer_2D_from_1D(dm,dmo);    
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

int opt_maurer_voronoi_edt_2D(DistMap *dm,Pnt *dmo_row,Pnt *im_row, Pnt *g, int h,int b,int by,int e,int ey);

int opt_maurer_voronoi_edt_2D_init(DistMap *dm,DistMap *dmo,Pnt *g,int step)
{
	int r = dm->rows;
	int c = dm->cols;
	for (int k=0;k<c;k+=step)
	{
		opt_maurer_voronoi_edt_2D(dm,dmo->data+k*c,dm->data+k*c,g,k,0,0,c-1,c-1);
	}
}


int opt_edt_maurer_2D_from_1D(DistMap *dm,DistMap *dmo)
{
	int stat;
	int i1; // same naming as in the paper
	Pnt *g;
	char *fname="edt_maurer_2D_from_1D";

	int c = dm->cols;
	int r = dm->rows;
	int step = STEP;

	g = (Pnt*)malloc(sizeof(Pnt)*c);

	opt_maurer_voronoi_edt_2D_init(dm,dmo,g,step);

/**/
	Pnt *dm_row;
	dm_row = dm->data;
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
		}
		step=step2;
	}/**/
	/*for (i1=0; i1 < nrows; ++i1, dm_row += ncols ) {
		stat = opt_maurer_voronoi_edt_2D(dm, dm_row,   g,i1);
		//CHECK_RET_STATUS(false);
	}*/
	free(g); //free(h);

	return 1;
}

int opt_remove_edt(Pnt *lst2,Pnt *lst1,Pnt *now,int h);

int opt_maurer_voronoi_edt_2D(DistMap *dm,Pnt *dmo_row ,Pnt *dm_row, Pnt *g,int h,int b,int by,int e,int ey)
{
	int l, i, ns, tmpx,tmpy, tmp1, tmp2, r,c;
	int infty, fi;
	Pnt tmpp;

	r = dm->rows; c=dm->cols;
	infty = MAXINT - r*r - c*c;

	l = -1;
	for (i=by; i <= ey; ++i){
		if ((fi = dm_row[i].x) != MARKER) {
			Pnt now;
			now.x = fi;
			now.y = i;
			while ( l >= 1 && opt_remove_edt(&g[l-1], &g[l], &now , h) )
			{
				--l;
			}
			++l; g[l].x = fi; g[l].y = i;
		}
	}

	// The following are lines 15-25 of the article
	if ((ns=l) == -1) return 1;
	l = 0;
	for (i=b; i <= e; ++i) {
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
		}

		dmo_row[i] = tmpp;
	}

	return 1;
}


float opt_interpnt(int x1,int y1,int x2,int y2,int x0);

inline float opt_interpnt(int x1,int y1,int x2,int y2,int x0)
{
	float xM = (float)(x1 + x2)/2.0f;
	float yM = (float)(y1 + y2)/2.0f;
	float nx = x2 - x1;
	float ny = y2 - y1;
	return yM + nx*(xM-x0)/ny;
}


inline int opt_remove_edt(Pnt *lst2,Pnt *lst1,Pnt *now,int h)
{
	/*
	int d2 = lst2->x - h;
	int d1 = lst1->x - h;
	int dn = now->x - h;
	if (d2>=0 && d1>=0 && dn>=0 && d1<=d2 && d1<=dn) return 0;
	if (d2<=0 && d1<=0 && dn<=0 && d1>=d2 && d1>=dn) return 0;
	*/
	float i1 = opt_interpnt(lst2->x,lst2->y,lst1->x,lst1->y,h);
	float i2 = opt_interpnt(lst1->x,lst1->y,now->x,now->y,h);
	return (i1>i2);
}


