#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define MAXN 5010
#define MAXM 5010
#define MAXVAL 2100000000
#define QUEMAX 10000000

struct qp
{
	int x,y,fx,fy;
	qp(){}
	qp(int x,int y,int fx,int fy):x(x),y(y),fx(fx),fy(fy){}
};

struct SedQP
{
	int x,y;
	SedQP(){}
	SedQP(int x,int y):x(x),y(y){}
};

int sed[MAXN][MAXM];
int sedLevel[MAXN][MAXM];
SedQP sedQueue[MAXN*MAXM];
int sedQueueNum;

int N,M;
int ori[MAXN][MAXM];
int res[MAXN][MAXM];



qp queData[2][QUEMAX];
int qh1,qr1;
int qh,qr;

int dx[8]={0,0,-1,1,-1,-1,1,1};
int dy[8]={-1,1,0,0,-1,1,-1,1};
int num = 0 ;

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
		sedLevel[sedQueue[i].x][sedQueue[i].y] = idx;
	}
}
int oneNum =0 ;
void findInit(qp *que)
{
	qh=qr=0;
	qh1=qr1=0;
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
					if (nx<0||nx>N+1||j<0||j>M+1) continue;
					if (ori[nx][ny]==1) 
					{
						f=true;
						break;
					}
				}
				if (f==true)
				{
					qp tmp = qp(i,j,i,j);
					que[qr++] = tmp;
				}
			}
			else res[i][j]=MAXVAL,oneNum++;
		}
}
int upNum =0;
void EDT(qp *que,qp *que1)
{
	int tarLevel = 0;
	while(qh!=qr)
	{
		int tqr = qr;
		tarLevel++;
		
		for (int i=qh;i!=qr;i++)
		{
			if (i>=QUEMAX) i=0;
			qp now = que[i];
			bool flag =false;
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
				int dis = sed[t1][t2];
				int nxtLevel = sedLevel[t1][t2];
				if (nxtLevel <= tarLevel && dis<res[nx][ny])
				{
					upNum++;
					res[nx][ny]=dis;
					que[tqr++]=nxt;
					if (tqr>=QUEMAX) tqr=0;
					num++;
				}
				else if (nxtLevel > tarLevel)
				{
					flag =true;
				}
			}
			if (flag)
			{
				que[tqr++]=now;
				if (tqr>=QUEMAX) tqr=0;
				num++;
			}
		}
		qh=qr;
		qr=tqr;
	}
	printf("%d %d\n",oneNum,upNum);
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



int main()
{
	freopen("T:\\workspace\\Data\\in8.txt","r",stdin);
	//freopen("T:\\workspace\\Data\\out22_2.txt","w",stdout);
	scanf("%d %d",&N,&M);
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			scanf("%d",&ori[i+1][j+1]);

	Init();
	clock_t t1 = clock();
	qp *que=queData[0],*que1=queData[1];
	findInit(que);
	// while(1) printf("%d\n",qr);

	EDT(que,que1);
	clock_t t2 = clock();
	printf("%d\n",t2-t1);
	printf("%d\n",num);
	//out();
	return 0;
}
