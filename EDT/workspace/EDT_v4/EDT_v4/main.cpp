#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define MAXN 1030
#define MAXM 1030
#define MAXVAL 2100000000
#define QUEMAX 30000
#define MAXALLNUM 1100000

int num=0;

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

rp allrp[MAXALLNUM];
int allnum = 0;

int sed[MAXN][MAXM];
int sedDis2Level[MAXN][MAXM];
int sedLevel2Dis[MAXN*MAXM];
SedQP sedQueue[MAXN*MAXM];
int sedQueueNum;

int N,M;
int ori[MAXN][MAXM];
int res[MAXN][MAXM];

List que[QUEMAX];
int queMax;

int dx[8]={0,0,-1,1,-1,-1,1,1};
int dy[8]={-1,1,0,0,-1,1,-1,1};

int cmp(const void *a,const void *b)
{
	SedQP *c = (SedQP *)a;
	SedQP *d = (SedQP *)b;
	int numc = c->x * c->x + c->y * c->y;
	int numd = d->x * d->x + d->y * d->y;
	return numc-numd;
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
	memset(que,0,sizeof(0));
	queMax = 0;
}

void findInit()
{
	for (int i=0;i<=N+1;i++)
		for (int j=0;j<=M+1;j++)
		{
			if (ori[i][j]==0)
			{
				res[i][j]=0;
				bool f = false;
				for (int k=0;k<8;k++)
				{
					int nx = i+dx[k];
					int ny = j+dy[k];
					if (nx<0||nx>N+1||ny<0||ny>M+1) continue;
					if (ori[nx][ny]==1) 
					{
						f=true;
						break;
					}
				}
				if (f==true)
				{
					rp *tmp  = (rp*)malloc(sizeof(rp));
					tmp->x = tmp->fx = i;
					tmp->y = tmp->fy = j;
					tmp->next = NULL;
					que[0].push(tmp);
				}
			}
			else res[i][j]=MAXVAL;
		}
}
void out()
{
	for (int i=1;i<=N;i++)
	{
		for (int j=1;j<M;j++)
		{
			printf("%d ",res[i][j]);
		}
		printf("%d\n",res[i][M]);
	}

}

void EDT()
{
	int d = -1;
	while (1)
	{
		//if (d!=-1) que[d].clear();
		d++;
		if (que[d].empty()) continue;
		if (d>queMax) break;
		for (rp *i = que[d].head ; i != NULL ; i = i->next)
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
				if (dis<res[nx][ny])
				{
					res[nx][ny]=dis;
					//rp *nxt = (rp*)malloc(sizeof(rp));
					rp *nxt = &allrp[allnum++];
					num++;
					nxt->x = nx;
					nxt->y = ny;
					nxt->fx = fx;
					nxt->fy = fy;
					nxt->next = NULL;
					que[l].push(nxt);
					if (l>queMax) queMax = l;
				}
			}
		}
	}
	//printf("%d\n",num);
}





int main()
{
	int *a,*b,c;
	freopen("T:\\workspace\\Data\\in22.txt","r",stdin);
	
	scanf("%d %d",&N,&M);
	memset(ori,0,sizeof(ori));
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			scanf("%d",&ori[i+1][j+1]);
	 
	Init();
	clock_t t1 = clock();
	findInit();
	// while(1) printf("%d\n",qr);

	EDT();
	clock_t t2 = clock();
	printf("%d\n",t2-t1);
	printf("%d\n",queMax);
	printf("%d\n",num);
	printf("%d\n",allnum);
	freopen("T:\\workspace\\Data\\out22_4.txt","w",stdout);
	out();
	return 0;
}
