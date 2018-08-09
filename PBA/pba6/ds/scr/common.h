#ifndef COMMON_H
#define COMMON_H 1

#define FG 1
#define MAXINT 0x7fffffff
#define MARKER -30000
#define RC(im,i,j) ((im)->data[(j) + (i)*(im->cols)])

#include <stdlib.h>
#include <stdio.h>
extern int STEP;
int STEP;
typedef struct 
{
	int rows,cols;
	int * data;
} Img;
Img * CreateImage(int rows,int cols);

typedef struct
{
	short x,y;
}Pnt;

typedef struct
{
	int rows,cols;
	Pnt *data;
} DistMap;

DistMap *CreateDistMap(int rows,int cols);
int getDist(Pnt *a,Pnt *b);

#endif
