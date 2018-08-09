#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define MAXN 5000+2
#define MAXM 5000+2
#define MAXVAL 2100000000
#define QUEMAX 30000000

int sed[MAXN][MAXM];

struct qp
{
	int x,y,fx,fy,level;
	qp(){}
	qp(int x,int y,int fx,int fy,int level):x(x),y(y),fx(fx),fy(fy),level(level){}
	int getSED()
	{
		//return (x-fx)*(x-fx)+(y-fy)*(y-fy);
		int t1,t2;
		if (x>fx) t1=x-fx;
		else t1=fx-x;
		if (y>fy) t2 = y-fy;
		else t2=fy-y;
		return sed[t1][t2];
	}
};

int N,M;
int ori[MAXN][MAXM];
int res[MAXN][MAXM];



qp queData[2][QUEMAX];
int qh1,qr1;
int qh,qr;

int dx[8]={0,0,-1,1,-1,-1,1,1};
int dy[8]={-1,1,0,0,-1,1,-1,1};
int num = 0 ,maxLevel=0;

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
					qp tmp = qp(i,j,i,j,0);
					que[qr++] = tmp;
				}
			}
			else res[i][j]=MAXVAL;
		}
}

void EDT(qp *que,qp *que1)
{

	while(qh!=qr)
	{
		qp now = que[qh++];
		if (qh>=QUEMAX) qh=0;
		bool flag = false;
		for (int k=0;k<4;k++)
		{
			int nx = now.x+dx[k];
			int ny = now.y+dy[k];

			if (nx<0||nx>N+1||ny<0||ny>M+1) continue;
			qp nxt = qp(nx,ny,now.fx,now.fy,now.level+1);
			int t1,t2;
			if (nxt.x>nxt.fx) t1=nxt.x-nxt.fx;
			else t1=nxt.fx-nxt.x;
			if (nxt.y>nxt.fy) t2 = nxt.y-nxt.fy;
			else t2=nxt.fy-nxt.y;
			//int dis = nxt.getSED();
			int dis = sed[t1][t2];
			if (dis<res[nx][ny])
			{
				res[nx][ny]=dis;
				que[qr++]=nxt;
				flag=true;
				num++;
				if (nxt.level>maxLevel) maxLevel = nxt.level;
				if (qr>=QUEMAX) qr=0;
			}
		}
		if (flag==false)
		{
			for (int k=4;k<8;k++)
			{
				int nx = now.x+dx[k];
				int ny = now.y+dy[k];
				if (nx<0||nx>N+1||ny<0||ny>M+1) continue;
				qp nxt = qp(nx,ny,now.fx,now.fy,now.level+1);
				int t1,t2;
				if (nxt.x>nxt.fx) t1=nxt.x-nxt.fx;
				else t1=nxt.fx-nxt.x;
				if (nxt.y>nxt.fy) t2 = nxt.y-nxt.fy;
				else t2=nxt.fy-nxt.y;
				//int dis = nxt.getSED();
				int dis = sed[t1][t2];
				if (dis<res[nx][ny])
				{
					res[nx][ny]=dis;
					que1[qr1++]=nxt;
					num++;
					if (nxt.level>maxLevel) maxLevel = nxt.level;
					if (qr>=QUEMAX) qr=0;
				}
			}
		}
	}
	if (qh1!=qr1)
	{
		qh=qh1;
		qr=qr1;
		qh1=qr1=0;
		EDT(que1,que);
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

void Init()
{
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			sed[i][j] = i*i+j*j;
}

int main()
{
	freopen("T:\\workspace\\Data\\in7.txt","r",stdin);
	//freopen("T:\\workspace\\Data\\out18_3.txt","w",stdout);
	scanf("%d %d",&N,&M);
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			scanf("%d",&ori[i+1][j+1]);

	Init();
	clock_t t1 = clock();
	qp *que=queData[0],*que1=queData[1];
	findInit(que);
	// while(1) printf("%d\n",qr);
    printf("%d\n",qr);
	EDT(que,que1);
	clock_t t2 = clock();
	printf("%d\n",t2-t1);
	printf("%d\n",num);
	printf("%d\n",maxLevel);
	//out();
	while(1){printf("");}
	return 0;
}
