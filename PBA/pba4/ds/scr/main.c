#include "common.h"
#include "maurer2003.h"
#include "optimized-maurer2003.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void out(FILE *fp,DistMap *dm)
{
	int r=dm->rows;
	int c=dm->cols;
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c-1;j++)
		{
			Pnt tmp = RC(dm,i,j);
			int tmpx = tmp.x-i;
			int tmpy = tmp.y-j;
			int dis = tmpx*tmpx + tmpy*tmpy ; 
			fprintf(fp,"%d ",dis);
		}
		Pnt tmp = RC(dm,i,c-1);
		int tmpx = tmp.x-i;
		int tmpy = tmp.y-c+1;
		int dis = tmpx*tmpx + tmpy*tmpy ; 
		fprintf(fp,"%d\n",dis);
	}
	/*for (int i = 0 ;i < im->rows;++i)
	{
		for (int j = 0;j<im->cols-1;j++) fprintf(fp,"%d ",RC(im,i,j));
		fprintf(fp,"%d\n",RC(im,i,im->cols-1));
	}*/
}



int *genorateData(int r,int c,float perc)
{
	int *in  = (int*)malloc(sizeof(int)*r*c);
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			in[i*c+j]=1;
		}
	}
	srand(time(0));
	int num = r*c*perc;
	//printf("num:%d\n",num);
	for (int i=0;i<num;i++)
	{
		int x = rand()%r;
		int y = rand()%c;
		while (in[x*c+y]==0)
		{
			x = rand()%r;
			y = rand()%c;
		}
		in[x*c+y]=0;
	}
	return in;
}

void check(DistMap *dm1,DistMap *dm2)
{
	int r = dm1->rows;
	int c = dm1->cols;
	int num=0;
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			Pnt a = RC(dm1,i,j);
			Pnt b = RC(dm2,i,j);
			if (a.x!=b.x||a.y!=b.y) num++;
		}
	}
	float perc = (float)num/(float)(r*c);
	printf("%d %f\n",num,perc);
}

int main(int argc,char* argv[])
{
	if (argc<=1) return 0;
	/*char inPath[] = "../data/in%s.txt";
	char moutPath[] = "../data/out%s_maurer2003.txt";
	char omoutPath[] = "../data/out%s_optimized-maurer2003.txt";
	sprintf(inPath,inPath,argv[1]);
	sprintf(moutPath,moutPath,argv[1]);
	sprintf(omoutPath,omoutPath,argv[1]);
	printf("%s\n%s\n%s\n",inPath,moutPath,omoutPath);
	freopen(inPath,"r",stdin);
	FILE *mfp = fopen(moutPath,"w");
	FILE *omfp = fopen(omoutPath,"w");*/
	int r,c;
	float p;
	r = atoi(argv[1]);
	c = atoi(argv[2]);
	p = atof(argv[3]);
	STEP = atoi(argv[4]);
	//printf("%d %d %f\n",r,c,p);
	int *in = genorateData(r,c,p);
	//scanf("%d%d",&r,&c);
	//Img * mImg = CreateImage(r,c);
	//Img * omImg = CreateImage(r,c);
	DistMap * mDistMap = CreateDistMap(r,c);
	DistMap * omDistMap = CreateDistMap(r,c);
	DistMap *omOut = CreateDistMap(r,c);
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++) 
		{
			int tmp;
			//scanf("%d",&tmp);
			//RC(omImg,i,j) = RC(mImg,i,j);
			//RC(mImg,i,j) = tmp;
			tmp = in[i*c+j];
			if (tmp == FG)
			{
				RC(omDistMap,i,j).x = MARKER;
				RC(omDistMap,i,j).y = MARKER;
				RC(mDistMap,i,j).x = MARKER;
				RC(mDistMap,i,j).y = MARKER;
			}
			else
			{
				RC(omDistMap,i,j).x = i;
				RC(omDistMap,i,j).y = j;
				RC(mDistMap,i,j).x = i;
				RC(mDistMap,i,j).y = j;
			}

		}
	}
	time_t s = clock();
	time_t e ;
	edt_maurer2003(mDistMap);
	e = clock();
	double cst = (double)(e-s)/CLOCKS_PER_SEC * 1000.0;
	printf("edt_maurer2003:%lf\n",cst);
	//out(mfp,mDistMap);
	s=clock();
	opt_edt_maurer2003(omDistMap,omOut);
	e = clock();
	cst = (double)(e-s)/CLOCKS_PER_SEC * 1000.0;
	printf("opt_edt_maurer2003:%lf\n",cst);
	//out(omfp,omOut);
	check(mDistMap,omOut);
	return 0;
}
