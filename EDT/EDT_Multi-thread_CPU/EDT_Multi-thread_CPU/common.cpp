#include "common.h"

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

