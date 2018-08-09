#include "common.h"

Img* CreateImage(int rows,int cols)
{
	Img * tmp = (Img*)malloc(sizeof(Img));
	tmp->rows = rows;
	tmp->cols = cols;
	tmp->data = (int*)malloc(rows*cols*sizeof(int));
	return tmp;
}

DistMap* CreateDistMap(int rows,int cols)
{
	DistMap *tmp = (DistMap*)malloc(sizeof(DistMap));
	tmp->rows = rows;
	tmp->cols = cols;
	tmp->data = (Pnt*)malloc(rows*cols*sizeof(Pnt));
	return tmp;
}

int getDist(Pnt *a,Pnt *b)
{
	int xa = a->x;
	int ya = a->y;
	int xb = b->x;
	int yb = b->y;;
	return (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
}

