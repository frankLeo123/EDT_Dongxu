#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define MAXN 1000+2
#define MAXM 1000+2
#define MAXVAL 2100000000
#define QUEMAX 3000000

struct qp
{
    int x,y,fx,fy,level;
    qp(){}
    qp(int x,int y,int fx,int fy,int level):x(x),y(y),fx(fx),fy(fy),level(level){}
    int getSED()
    {
        return (x-fx)*(x-fx)+(y-fy)*(y-fy);
    }
};

int N,M;
int ori[MAXN][MAXM];
int res[MAXN][MAXM];

qp que[QUEMAX];
int qh,qr;

int dx[8]={0,0,-1,1,-1,-1,1,1};
int dy[8]={-1,1,0,0,-1,1,-1,1};


void findInit()
{
     qh=qr=0;
     for (int i=0;i<=N+1;i++)
         for (int j=0;j<=M+1;j++)
         {
             res[i][j]=MAXVAL;
             if (ori[i][j]==0)
             {
                qp tmp = qp(i,j,i,j,0);
                que[qr++] = tmp;
             }
         }
}
int num=0;
void EDT()
{
     
     while(qh!=qr)
     {
         qp now = que[qh++];
         for (int k=0;k<4;k++)
         {
             int nx = now.x+dx[k];
             int ny = now.y+dy[k];
             if (nx<0||nx>N+1||ny<0||ny>M+1) continue;
             qp nxt = qp(nx,ny,now.fx,now.fy,now.level+1);
             int dis = nxt.getSED();
             if (dis<res[nx][ny])
             {
                 res[nx][ny]=dis;
                 que[qr++]=nxt;
                 num++;
             }
         }
     }
}

void out()
{
     for (int i=1;i<=N;i++)
     {
         for (int j=1;j<M;j++)
         {
             printf("%4d ",res[i][j]);
         }
         printf("%4d\n",res[i][M]);
     }
     
}

int main()
{
    freopen("Data//in7.txt","r",stdin);
    //freopen("out7_2.txt","w",stdout);
    scanf("%d %d",&N,&M);
    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
            scanf("%d",&ori[i+1][j+1]);
    clock_t t1 = clock();
    findInit();
    printf("%d\n",qr);
    EDT();
    clock_t t2 = clock();
    printf("%d\n",t2-t1);
    printf("%d",num);
    //out();
    while(1){printf("");}
    return 0;
}
