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



DistMap *genorateData1(int r,int c,int pn)
{
	DistMap *in  = (DistMap*)malloc(sizeof(DistMap));
	in->rows = r;in->cols = c;
	in->data = (Pnt*)malloc(r*c*sizeof(Pnt));
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			in->data[i*c+j].x=MARKER;
			in->data[i*c+j].y=MARKER;
		}
	}
	srand(time(0));
	int num = pn;
	//printf("num:%d\n",num);
	for (int i=0;i<num;i++)
	{
		int x = rand()%r;
		int y = rand()%c;
		while (in->data[x*c+y].x!=MARKER)
		{
			x = rand()%r;
			y = rand()%c;
		}
		in->data[x*c+y].x=x;
		in->data[x*c+y].y=y;
	}
	return in;
}

DistMap *genorateData2(int r,int c,DistMap *in)
{
	int or=in->rows,oc=in->cols;
	DistMap *in1  = (DistMap*)malloc(sizeof(DistMap));
	in1->rows = r;in1->cols = c;
	in1->data = (Pnt*)malloc(r*c*sizeof(Pnt));
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			in1->data[i*c+j].x=MARKER;
			in1->data[i*c+j].y=MARKER;
		}
	}
	double px = (double)(r)/in->rows;
	double py = (double)(c)/in->cols;
	for (int i=0;i<or;i++)
	{
		for (int j=0;j<oc;j++)
		{
			if (in->data[i*oc+j].x!=MARKER)
			{
				double nx = in->data[i*oc+j].x;
				double ny = in->data[i*oc+j].y;
			   	int x = nx*px;
				int y = ny*py;
				in1->data[x*c+y].x = x;
				in1->data[x*c+y].y = y;
			}
		}
	}
	return in1;
}

bool gcheck3(DistMap *in,int x,int y1,int y2)
{
	int c=in->cols;
	for (int i=y1;i<=y2;i++)
	{
		if ((in->data[x*c + i]).x!=MARKER) return false;
	}
	return true;
}

DistMap *genorateData3(int r,int c,int pn)
{
	DistMap *in  = (DistMap*)malloc(sizeof(DistMap));
	in->rows = r;in->cols = c;
	in->data = (Pnt*)malloc(r*c*sizeof(Pnt));
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			in->data[i*c+j].x=MARKER;
			in->data[i*c+j].y=MARKER;
		}
	}
	srand(time(0));
	int num = pn;
	//printf("num:%d\n",num);
	for (int i=0;i<num;i++)
	{
		int x = rand()%r;
		int y1 = rand()%c;
		int y2 = (rand()%(c-y1))+y1;
		while (gcheck3(in,x,y1,y2)==false)
		{
			x = rand()%r;
			y1 = rand()%c;
			y2 = y2 = (rand()%(c-y1))+y1;
		}
		for (int j = y1;j<=y2;j++) 
		{
			in->data[x*c+j].x = x;
			in->data[x*c+j].y = j;
		}
	}
	return in;
}


bool gcheck4(DistMap *in,int x,int y,int d)
{
	int c=in->cols;
	for (int i=0;i<d;i++)
	{
		if ((in->data[(x+i)*c + y + i]).x!=MARKER) return false;
	}
	return true;
}

DistMap *genorateData4(int r,int c,int pn)
{
	DistMap *in  = (DistMap*)malloc(sizeof(DistMap));
	in->rows = r;in->cols = c;
	in->data = (Pnt*)malloc(r*c*sizeof(Pnt));
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			in->data[i*c+j].x=MARKER;
			in->data[i*c+j].y=MARKER;
		}
	}
	srand(time(0));
	int num = pn;
	//printf("num:%d\n",num);
	for (int i=0;i<num;i++)
	{
		int x = rand()%r;
		int y = rand()%c;
		int tmp = min(r-x,c-y);
		int d = rand()%tmp;
		while (gcheck4(in,x,y,d)==false)
		{
			x = rand()%r;
			y = rand()%c;
			tmp = min(r-x,c-y);
			d =  rand()%tmp;
		}
		for (int j = 0;j<d;j++) 
		{
			in->data[(x+j)*c+y+j].x = x+j;
			in->data[(x+j)*c+y+j].y = y+j;
		}
	}
	return in;
}

bool gcheck5(DistMap *in,int x1,int x2,int y)
{
	int c=in->cols;
	for (int i=x1;i<=x2;i++)
	{
		if ((in->data[i*c + y]).x!=MARKER) return false;
	}
	return true;
}

DistMap *genorateData5(int r,int c,int pn)
{
	DistMap *in  = (DistMap*)malloc(sizeof(DistMap));
	in->rows = r;in->cols = c;
	in->data = (Pnt*)malloc(r*c*sizeof(Pnt));
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			in->data[i*c+j].x=MARKER;
			in->data[i*c+j].y=MARKER;
		}
	}
	srand(time(0));
	int num = pn;
	//printf("num:%d\n",num);
	for (int i=0;i<num;i++)
	{
		int y = rand()%c;;
		int x1 = rand()%r;
		int x2 = (rand()%(r-x1))+x1;
		while (gcheck5(in,x1,x2,y)==false)
		{
			y = rand()%c;
			x1 = rand()%r;
			x2 = (rand()%(r-x1))+x1;
		}
		for (int j = x1;j<=x2;j++) 
		{
			in->data[j*c+y].x = j;
			in->data[j*c+y].y = y;
		}
	}
	return in;
}

DistMap *genorateData6(int r,int c,int pn)
{
	DistMap *in  = (DistMap*)malloc(sizeof(DistMap));
	in->rows = r;in->cols = c;
	in->data = (Pnt*)malloc(r*c*sizeof(Pnt));
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<c;j++)
		{
			in->data[i*c+j].x=MARKER;
			in->data[i*c+j].y=MARKER;
		}
	}
	srand(time(0));
	int num1 = pn*0.5;
	for (int i=0;i<num1;i++)
	{
		int x = rand()%r;
		int y = rand()%c;
		while (in->data[x*c+y].x!=MARKER&&x<(r/2)&&y<(c/2))
		{
			x = rand()%r;
			y = rand()%c;
		}
		in->data[x*c+y].x=x;
		in->data[x*c+y].y=y;
	}
	int num2 = pn*0.1;
	for (int i=0;i<num2;i++)
	{
		int x = rand()%r;
		int y = rand()%c;
		while (in->data[x*c+y].x!=MARKER&&x>(r/2)&&y<(c/2))
		{
			x = rand()%r;
			y = rand()%c;
		}
		in->data[x*c+y].x=x;
		in->data[x*c+y].y=y;
	}
	int num3 = pn*0.2;
	for (int i=0;i<num3;i++)
	{
		int x = rand()%r;
		int y = rand()%c;
		while (in->data[x*c+y].x!=MARKER&&x<(r/2)&&y>(c/2))
		{
			x = rand()%r;
			y = rand()%c;
		}
		in->data[x*c+y].x=x;
		in->data[x*c+y].y=y;
	}
	int num4 = pn*0.2;
	for (int i=0;i<num4;i++)
	{
		int x = rand()%r;
		int y = rand()%c;
		while (in->data[x*c+y].x!=MARKER&&x>(r/2)&&y>(c/2))
		{
			x = rand()%r;
			y = rand()%c;
		}
		in->data[x*c+y].x=x;
		in->data[x*c+y].y=y;
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
			if (a.x!=b.x||a.y!=b.y) 
			{
				int disa= (i-a.x)*(i-a.x) + (j-a.y)*(j-a.y);
				int disb= (i-b.x)*(i-b.x) + (j-b.y)*(j-b.y);
				if (disa!=disb) 
				{
					
					num++;
					//printf("%d %d\n",i,j);
				}
			}
		}
	}
	float perc = (float)num/(float)(r*c);
	printf("%d %f\n",num,perc);
}

double runCase(DistMap *in,int step,double *tm1,double *tm2)
{
	int r=in->rows,c=in->cols;
	DistMap * mDistMap = CreateDistMap(r,c);
	DistMap * omDistMap = CreateDistMap(r,c);
	DistMap *omOut = CreateDistMap(r,c);
	Pnt **icp = (Pnt **)malloc(r*sizeof(Pnt*));
	for (int i=0;i<r;i++) icp[i]=(Pnt*)malloc(c*sizeof(Pnt));
	int *icpl = (int *)malloc(r*sizeof(int));
	for (int i=0;i<r;i++)
	{
		icpl[i]=-1;
		for (int j=0;j<c;j++) 
		{
			Pnt tmp;
			tmp = in->data[i*c+j];
			if (tmp.x == MARKER)
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
			icp[i][j].x=icp[i][j].y=-MARKER;
			RC(omOut,i,j).x = RC(omOut,i,j).y = MARKER;
		}
	}
	time_t s = clock();
	time_t e ;
	edt_maurer2003(mDistMap);
	e = clock();
	*tm1 = (double)(e-s)/CLOCKS_PER_SEC * 1000.0;




	s=clock();
	opt_edt_maurer2003(omDistMap,omOut,icp,icpl,step);
	e = clock();
	*tm2 = (double)(e-s)/CLOCKS_PER_SEC * 1000.0;
	
	check(mDistMap,omOut);
	for (int i=0;i<r;i++) free(icp[i]);
	free(icp);
	free(icpl);
	free(mDistMap->data);
	free(omDistMap->data);
	free(omOut->data);
	free(mDistMap);
	free(omDistMap);
	free(omOut);
	return 0;
}

void runAndSatistics(DistMap *in,int step,int caseNum,double *atm1,double *atm2,double *aperc)
{
	(*aperc ) =-1000.0;
	for (int i=0;i<caseNum;i++)
	{
		double tm1,tm2,perc;
		runCase(in,step,&tm1,&tm2);
		perc = (tm1-tm2)/tm2;
		if (perc>(*aperc))
		{
			(*atm1) = tm1;
			(*atm2) = tm2;
			(*aperc) = perc;
		}
	}
}

void saveInput(DistMap *in)
{
	FILE *fp;
	fp=fopen("in.txt","w");
	fprintf(fp,"%d %d\n",in->rows,in->cols);
	for (int i=0;i<in->rows;i++)
	{
		for (int j=0;j<(in->cols)-1;j++)
		{
			fprintf(fp,"%d %d ",in->data[ i*in->cols + j].x,in->data[ i*in->cols + j].y);
		}
		fprintf(fp,"%d %d\n",in->data[ i*in->cols + in->cols -1].x,in->data[ i*in->cols + in->cols -1].y);
	}
	fclose(fp);
}

void getInput(DistMap *in)
{
	FILE *fp;
	fp=fopen("T:\\tmp\\EDT_Multi-thread_CPUnew\\EDT_Multi-thread_CPU\\in.txt","r");
	fscanf(fp,"%d %d",&(in->rows),&(in->cols));
	for (int i=0;i<in->rows;i++)
	{
		for (int j=0;j<(in->cols);j++)
		{
			fscanf(fp,"%d%d",&(in->data[ i*in->cols + j].x),&(in->data[ i*in->cols + j].y));
		}
	}
	fclose(fp);
}

/*
 *./main pn caseNum size1 step1 size2 step2 .....
 *
 *
 */

int main(int argc,char* argv[])
{
	if (argc<5) 
	{
		printf("./main pn caseNum size1 step1 size2 step2..\n");
		return 0;
	}
	int pn = atoi(argv[1]);
	int caseNum = atoi(argv[2]);
	DistMap *in;
	for (int i=3;i<argc;i+=2)
	{
		int size = atoi(argv[i]);
		int step = atoi(argv[i+1]);
		if (i==3) in = genorateData1(size,size,pn);
		else 
		{
			DistMap *in1 = genorateData2(size,size,in);
			free(in);
			in = in1;
		}
		//in = genorateData5(size,size,pn);
		double atm1,atm2,aperc;
		//saveInput(in);
		//getInput(in);
		in = genorateData6(size,size,pn);
		runAndSatistics(in,step,caseNum,&atm1,&atm2,&aperc);
		printf("%d*%d image with %d points : %.0lfms %.0lfms(new) %.2lf%%\n",size,size,step,atm1,atm2,aperc*100.0);
	}
	return 0;
}



/*
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXN 1030
#define MAXM 1030
#define MAXVAL 2100000000

#define MAX_THREADS 4
#define QUEMAX 60000
#define MAXALLNUM 600000
struct SedQP
{
	int x,y;
	SedQP(){}
	SedQP(int x,int y):x(x),y(y){}
};


struct rp
{
	int x,y,fx,fy;
	rp *next;
	rp(){next = NULL;}
	rp(int x,int y):x(x),y(y){next = NULL;}
};

struct List
{
	rp *head,*tail;
	List(){ this->head = this->tail = NULL ;}
	void clear()
	{
		for (rp *i = this->head ; i != NULL ; ) 
		{
			rp * tmp = i;
			i = i->next;
			free(tmp);
		}
		this->head = this->tail = NULL ;
	}
	void push(rp *node)
	{
		if (this->head==NULL) this->head = this->tail =  node;
		else 
		{
			this->tail->next = node;
			this->tail = node;
		}
		this->tail->next = NULL;
	}
	bool empty()
	{
		if (this->head == this->tail && this->head == NULL) return true;
		else return false;
	}
};


int N,M;
int ori[MAXN][MAXM];
int res[MAXN][MAXM];

rp* ique[QUEMAX];
int iqh,iqr;

rp allrp[MAX_THREADS][MAXALLNUM];
int allnum[MAX_THREADS] = {0};

List que[MAX_THREADS][QUEMAX];
int queMax[MAX_THREADS];

int dx[8]={0,0,-1,1,-1,-1,1,1};
int dy[8]={-1,1,0,0,-1,1,-1,1};
int num[MAX_THREADS] = {0} ;

CRITICAL_SECTION csRes[MAXN][MAXM];
HANDLE synLocker;

int sed[MAXN][MAXM];
int sedDis2Level[MAXN][MAXM];
int sedLevel2Dis[MAXN*MAXM];
SedQP sedQueue[MAXN*MAXM];
int sedQueueNum;

int cmp(const void *a,const void *b)
{
	SedQP *c = (SedQP *)a;
	SedQP *d = (SedQP *)b;
	int numc = c->x * c->x + c->y * c->y;
	int numd = d->x * d->x + d->y * d->y;
	return numc-numd;
}


struct FindInitKernelParam_t
{
	int i;
};

FindInitKernelParam_t fkParam[MAX_THREADS];

void findInitKernel(int i)
{
	for (int j=0;j<=M+1;j++)
	{
		res[i][j]=MAXVAL;
        if (ori[i][j]==0)
        {
			bool flag = false;
			for (int k=0;k<8;k++)
			{
				int nx = i+ dx[k];
				int ny = j+dy[k];
				if (ori[nx][ny]==1)
				{
					flag =true;
					break;
				}
			}
			if (flag)
			{
				rp * nr = &allrp[0][allnum[0]++];
				nr->x = nr->fx = i;
				nr->y = nr->fy = j;
				nr->next = NULL;
				ique[iqr] = nr;
				iqr++;
			}
        }
	}
	return ;
}

void findInit()
{
	iqh=iqr;
	for (int i=0;i<=N+1;i++) findInitKernel(i);
}



DWORD WINAPI EDTKernel(LPVOID lpParam)
{
	WaitForSingleObject(synLocker,INFINITE);
	FindInitKernelParam_t *tmp = (FindInitKernelParam_t*)lpParam;
	int ti = tmp->i;
	int d = -1;
	while (1)
	{
		//if (d!=-1) que[d].clear();
		d++;
		if (que[ti][d].empty()) continue;
		if (d>queMax[ti]) break;
		for (rp *i = que[ti][d].head ; i != NULL ; i = i->next)
		{
			int x = i->x;
			int y = i->y;
			int fx = i->fx;
			int fy = i->fy;
			for (int j=0;j<4;j++)
			{
				int nx = x + dx[j];
				int ny = y + dy[j];
				if (nx<0||nx>N+1||ny<0||ny>M+1) continue;
				int rx = nx-fx,ry = ny-fy;
				if (rx<0) rx = -rx;
				if (ry<0) ry = -ry;
				int l = sedDis2Level[rx][ry];
				if (l<d) continue;
				int dis = sed[rx][ry];
				//EnterCriticalSection(&csRes[nx][ny]);
				if (dis<res[nx][ny])
				{
					res[nx][ny]=dis;
					//LeaveCriticalSection(&csRes[nx][ny]);
					rp *nxt = &allrp[ti][allnum[ti]++];
					num[ti]++;
					nxt->x = nx;
					nxt->y = ny;
					nxt->fx = fx;
					nxt->fy = fy;
					nxt->next = NULL;
					que[ti][l].push(nxt);
					if (l>queMax[ti]) queMax[ti] = l;
				}
				//else LeaveCriticalSection(&csRes[nx][ny]);
			}
		}
	}
	return 0;
}

void EDT()
{
	HANDLE hThreads[MAX_THREADS];
	for (int i=0;i<MAX_THREADS;i++)
	{
		fkParam[i].i = i;
		hThreads[i] = CreateThread(NULL,0,EDTKernel,&fkParam[i],0,NULL);
		if (hThreads[i]==NULL)
		{
			printf("Thread Create Error At %d\n",i);
		}
	}
	ReleaseSemaphore(synLocker,MAX_THREADS,NULL);
	WaitForMultipleObjects(MAX_THREADS,hThreads,true,INFINITE);
	for (int i=0;i<MAX_THREADS;i++) CloseHandle(hThreads[i]);
}

void out()
{
     for (int i=1;i<=N;i++)
     {
         for (int j=1;j<M;j++)
         {
			 if (ori[i][j]==0) printf("0 ");
             else printf("%d ",res[i][j]);
         }
         if (ori[i][M]==0) printf("0\n");
         else printf("%d\n",res[i][M]);
     }
     
}

void Init()
{
	sedQueueNum = 0;
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
		{
			sed[i][j] = i*i+j*j;
			sedQueue[sedQueueNum++] = SedQP(i,j);
		}
	qsort(sedQueue,sedQueueNum,sizeof(SedQP),cmp);
	int idx = 0;
	int idxv = sedQueue[0].x*sedQueue[0].x + sedQueue[0].y*sedQueue[0].y;
	for (int i=0;i<sedQueueNum;i++)
	{
		int nowv = sedQueue[i].x*sedQueue[i].x + sedQueue[i].y*sedQueue[i].y;
		if (nowv != idxv)
		{
			idxv = nowv;
			idx++;
		}
		sedDis2Level[sedQueue[i].x][sedQueue[i].y] = idx;
	}
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
		{
			sedLevel2Dis[sedDis2Level[i][j]] = sed[i][j];
		}
	memset(que,0,sizeof(que));
	memset(queMax,0,sizeof(queMax));

}

int main()
{
    freopen("T:\\workspace\\Data\\in22.txt","r",stdin);

    scanf("%d %d",&N,&M);
    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
            scanf("%d",&ori[i+1][j+1]);
	synLocker = CreateSemaphore(NULL,0,MAX_THREADS,NULL);
	for (int i=0;i<N+2;i++)
        for (int j=0;j<M+2;j++)
		{
           	if (!InitializeCriticalSectionAndSpinCount(&csRes[i][j],0x00000010))
			{
				printf("InitializeCriticalSectionAndSpinCount Error!\n");
				return 0;
			}
		}
	Init();
	clock_t t1 = clock();
    findInit();
	
	int idx = 0;
	for (int i=0;i<iqr;i++)
	{
		//que[idx][qr[idx]++]=ique[i];
		que[idx][0].push(ique[i]);
		idx++;
		idx%=MAX_THREADS;
	}
	
    EDT();
    clock_t t2 = clock();
    printf("%d\n",t2-t1);
	int sum =0 ;
	int sm = -1;
	for (int i=0;i<MAX_THREADS;i++) {sum+=num[i];printf("%d\n",num[i]);sm = sm<queMax[i]?queMax[i]:sm;}
	printf("%d\n",sum);
	printf("QUEMAX:%d\n",sm);
	sm = -1;
	for (int i=0;i<MAX_THREADS;i++) {sm = sm<allnum[i]?allnum[i]:sm;}
	printf("MAXALLNUM:%d\n",sm);

	freopen("T:\\workspace\\Data\\out22_MP_CPU.txt","w",stdout);
    out();
	//printf("%d\n",qr);
	CloseHandle(synLocker);
	for (int i=0;i<N+2;i++)
        for (int j=0;j<M+2;j++)
			DeleteCriticalSection(&csRes[i][j]);
    return 0;
}*/
/*#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXN 2000+2
#define MAXM 2000+2
#define MAXVAL 2100000000
#define QUEMAX 3000000
#define MAX_THREADS 8

int sed[MAXN][MAXM];

struct qp
{
    int x,y,fx,fy;
    qp(){}
    qp(int x,int y,int fx,int fy):x(x),y(y),fx(fx),fy(fy){}
};

int N,M;
int ori[MAXN][MAXM];
int res[MAXN][MAXM];

qp que[MAX_THREADS][QUEMAX];
int qh[MAX_THREADS],qr[MAX_THREADS];
qp ique[QUEMAX];
int iqh,iqr;

int dx[8]={0,0,-1,1,-1,-1,1,1};
int dy[8]={-1,1,0,0,-1,1,-1,1};
int num = 0 ;
CRITICAL_SECTION csRes[MAXN][MAXM];


struct FindInitKernelParam_t
{
	int i;
};

FindInitKernelParam_t fkParam[MAX_THREADS];

void findInitKernel(int i)
{
	for (int j=0;j<=M+1;j++)
	{
		res[i][j]=MAXVAL;
        if (ori[i][j]==0)
        {
			bool flag = false;
			for (int k=0;k<8;k++)
			{
				int nx = i+ dx[k];
				int ny = j+dy[k];
				if (ori[nx][ny]==1)
				{
					flag =true;
					break;
				}
			}
			if (flag)
			{
				qp tmp = qp(i,j,i,j);
				ique[iqr] = tmp;
				iqr++;
			}
        }
	}
	return ;
}

void findInit()
{
	memset(qh,0,sizeof(qh));
	memset(qr,0,sizeof(qr));
	iqh=iqr;
	for (int i=0;i<=N+1;i++) findInitKernel(i);
}



DWORD WINAPI EDTKernel(LPVOID lpParam)
{
	FindInitKernelParam_t *tmp = (FindInitKernelParam_t*)lpParam;
	int i = tmp->i;
	while(1)
	{
		
		qp now;
		if (qh[i]!=qr[i])
		{
			now = que[i][qh[i]++];
			if (qh[i]>=QUEMAX) qh[i]=0;
		}
		else break;
		
		bool flag = false;
		for (int k=0;k<4;k++)
		{
			int nx = now.x+dx[k];
			int ny = now.y+dy[k];

			if (nx<0||nx>N+1||ny<0||ny>M+1) continue;
			qp nxt = qp(nx,ny,now.fx,now.fy);
			int t1,t2;
			if (nxt.x>nxt.fx) t1=nxt.x-nxt.fx;
			else t1=nxt.fx-nxt.x;
			if (nxt.y>nxt.fy) t2 = nxt.y-nxt.fy;
			else t2=nxt.fy-nxt.y;
			//int dis = nxt.getSED();
			int dis = sed[t1][t2];
			//EnterCriticalSection(&csRes[nx][ny]);
			if (dis<res[nx][ny])
			{
				res[nx][ny]=dis;
				//LeaveCriticalSection(&csRes[nx][ny]);
				que[i][qr[i]++]=nxt;
				if (qr[i]>=QUEMAX) qr[i]=0;
				flag=true;
				num++;
			}
			//else LeaveCriticalSection(&csRes[nx][ny]);
			
		}
	}
	return 0;
}

void EDT()
{
	HANDLE hThreads[MAX_THREADS];
	for (int i=0;i<MAX_THREADS;i++)
	{
		fkParam[i].i = i;
		hThreads[i] = CreateThread(NULL,0,EDTKernel,&fkParam[i],0,NULL);
		if (hThreads[i]==NULL)
		{
			printf("Thread Create Error At %d\n",i);
		}
	}
	WaitForMultipleObjects(MAX_THREADS,hThreads,true,INFINITE);
	for (int i=0;i<MAX_THREADS;i++) CloseHandle(hThreads[i]);
}

void out()
{
     for (int i=1;i<=N;i++)
     {
         for (int j=1;j<M;j++)
         {
			 if (ori[i][j]==0) printf("0 ");
             else printf("%d ",res[i][j]);
         }
         if (ori[i][M]==0) printf("0\n");
         else printf("%d\n",res[i][M]);
     }
     
}

void Init()
{
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			sed[i][j] = i*i+j*j;

}

int main()
{
    freopen("T:\\workspace\\Data\\in8.txt","r",stdin);

    scanf("%d %d",&N,&M);
    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
            scanf("%d",&ori[i+1][j+1]);
	for (int i=0;i<N+2;i++)
        for (int j=0;j<M+2;j++)
		{
           	if (!InitializeCriticalSectionAndSpinCount(&csRes[i][j],0x00000010))
			{
				printf("InitializeCriticalSectionAndSpinCount Error!\n");
				return 0;
			}
		}
	Init();
	
    findInit();
	
	int idx = 0;
	for (int i=0;i<iqr;i++)
	{
		que[idx][qr[idx]++]=ique[i];
		idx++;
		idx%=MAX_THREADS;
	}
	clock_t t1 = clock();
    EDT();
    clock_t t2 = clock();
    printf("%d\n",t2-t1);
	printf("%d\n",num);
	freopen("T:\\workspace\\Data\\out7_MP_CPU.txt","w",stdout);
    out();
	//printf("%d\n",qr);

	for (int i=0;i<N+2;i++)
        for (int j=0;j<M+2;j++)
			DeleteCriticalSection(&csRes[i][j]);
    return 0;
}*/