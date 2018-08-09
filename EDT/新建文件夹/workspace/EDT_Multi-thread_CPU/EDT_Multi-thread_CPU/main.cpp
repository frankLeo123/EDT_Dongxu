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
}
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