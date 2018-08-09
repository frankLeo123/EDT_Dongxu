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

#include "maurer2003.h"
#include "common.h"

int edt_1d_vertical(DistMap *dm)
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

int edt_maurer_2D_from_1D(DistMap *im);

int edt_maurer2003(DistMap * dm)
{
	char *fname="edt_maurer2003";
	int i,j,r,c;
	int stat;
	int infty;

	//assert(im->isbinary);
	r = dm->rows; c = dm->cols;
	infty = MAXINT - r*r - c*c;
	// Vertical columnwise EDT
	stat = edt_1d_vertical(dm);
	// Lotufo's 1D EDT is equivalent to Maurer's D1.
	// There is a remark in section 5 of Maurer's paper that says
	// D1 can be calculated using a routine like this.
	//CHECK_RET_STATUS(false);
	stat = edt_maurer_2D_from_1D(dm);    
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

int maurer_voronoi_edt_2D(DistMap *dm, Pnt *im_row, Pnt *g, int h);


int edt_maurer_2D_from_1D(DistMap *dm)
{
	int stat;
	int i1; // same naming as in the paper
	Pnt *g;
	char *fname="edt_maurer_2D_from_1D";

	int ncols = dm->cols;
	int nrows = dm->rows;

	// Call voronoi_edt_2D for every row.
	// OBS: g and h are internal to maurer_voronoi_edt_2d and are
	// pre-allocated here for efficiency.
	g = (Pnt*)malloc(sizeof(Pnt)*ncols);
	/*animal_malloc_array(int, ncols);
	  if (!g) {
	  animal_err_flush_trace();                                
	  animal_err_register(fname, ANIMAL_ERROR_MALLOC_FAILED,"");  
	  return false;                                          
	  } */                                                       
	//h = (int*)malloc(sizeof(int)*ncols);
	/*animal_malloc_array(int, ncols);
	  if (!h) {
	  animal_err_flush_trace();                                
	  animal_err_register(fname, ANIMAL_ERROR_MALLOC_FAILED,"");  
	  return false;                                          
	  } */                                                       


	Pnt *dm_row;
	dm_row = dm->data;
	for (i1=0; i1 < nrows; ++i1, dm_row += ncols ) {
		stat = maurer_voronoi_edt_2D(dm, dm_row,  /* internal: */ g,i1);
		//CHECK_RET_STATUS(false);
	}
	free(g); //free(h);

	return 1;
}

int remove_edt(Pnt *lst2,Pnt *lst1,Pnt *now,int h);

int maurer_voronoi_edt_2D(DistMap *dm, Pnt *dm_row, Pnt *g,int h)
{
	int l, i, ns, tmpx,tmpy, tmp1, tmp2, r,c;
	int infty, fi;
	Pnt tmpp;

	r = dm->rows; c=dm->cols;
	infty = MAXINT - r*r - c*c;

	l = -1;
	for (i=0; i < c; ++i){
		if ((fi = dm_row[i].x) != MARKER) {
			Pnt now;
			now.x = fi;
			now.y = i;
			while ( l >= 1 && remove_edt(&g[l-1], &g[l], &now , h) )
			{
				--l;
			}
			++l; g[l].x = fi; g[l].y = i;
		}
	}

	// The following are lines 15-25 of the article
	if ((ns=l) == -1) return 1;
	l = 0;
	for (i=0; i < c; ++i) {
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

		dm_row[i] = tmpp;
	}

	return 1;
}


float interpnt(int x1,int y1,int x2,int y2,int x0);

inline float interpnt(int x1,int y1,int x2,int y2,int x0)
{
	float xM = (float)(x1 + x2)/2.0f;
	float yM = (float)(y1 + y2)/2.0f;
	float nx = x2 - x1;
	float ny = y2 - y1;
	return yM + nx*(xM-x0)/ny;
}


inline int remove_edt(Pnt *lst2,Pnt *lst1,Pnt *now,int h)
{
	float i1 = interpnt(lst2->x,lst2->y,lst1->x,lst1->y,h);
	float i2 = interpnt(lst1->x,lst1->y,now->x,now->y,h);
	return (i1>i2);
}


