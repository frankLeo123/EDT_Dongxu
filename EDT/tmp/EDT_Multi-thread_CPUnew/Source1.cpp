inline float interpnt(int x1,int y1,int x2,int y2,int x0)
{
	float xM = (float)(x1 + x2)/2.0f;
	float yM = (float)(y1 + y2)/2.0f;
	float nx = x2 - x1;
	float ny = y2 - y1;
	return yM + nx*(xM-x0)/ny;
}


inline int remove_edt(Pnt *mlst2,Pnt *mlst1,Pnt *now,int h)
{
	float i1 = interpnt(mlst2->x,mlst2->y,mlst1->x,mlst1->y,h);
	float i2 = interpnt(mlst1->x,mlst1->y,now->x,now->y,h);
	return (i1>i2);
}

int edt_maurer2003(DistMap * dm)
{
	int ncols = dm->cols;
	int nrows = dm->rows;
	for (int c=0; c < ncols; c++) {
		Pnt pre = RC(dm,0,c);
		for (int r=1; r<nrows; r++)
		{
			if (RC(dm,r,c).x==MARKER) RC(dm,r,c) = pre;
			else pre = RC(dm,r,c);
		}
		for (int r=nrows-2; r >= 0; r--) 
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
	
	Pnt *g = (Pnt*)malloc(sizeof(Pnt)*ncols);
	Pnt *dm_row = dm->data;
	for (int i1=0; i1 < nrows; ++i1, dm_row += ncols ) {
		int l = -1;
		for (int i=0; i < ncols; ++i){
			if (dm_row[i].x != MARKER) {
				Pnt now;
				now.x = dm_row[i].x;
				now.y = i;
				while ( l >= 1 && remove_edt(&g[l-1], &g[l], &now , h) )
				{
					--l;
				}
				++l; g[l].x = dm_row[i].x; g[l].y = i;
			}
		}
		int ns = l
		if (ns == -1) continue;
		l = 0;
		for (int i=0; i < ncols; ++i) {
			int tmpy = g[l].y - i;
			int tmpx = g[l].x-h;
			int tmp1 = tmpx*tmpx + tmpy*tmpy;
			Pnt tmpp = g[l];
			while(1) {
				if (l >= ns) break;

				tmpy = g[l+1].y - i;
				tmpx = g[l+1].x - h;
				int tmp2 = tmpx*tmpx + tmpy*tmpy; 
				if (tmp1 <= tmp2) break;
				++l;
				tmp1 = tmp2;
				tmpp = g[l];
			}

			dm_row[i] = tmpp;
		}
	}
	free(g);
}





