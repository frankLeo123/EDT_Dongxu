#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <time.h>

#define MAXN 5000+2
#define MAXM 5000+2
#define FOREGROUND_PIXEL 1
#define MAXVAL 2100000000

int N,M;
int oriIma[MAXN][MAXM];
int resIma[MAXN][MAXM];
int Rx[MAXN][MAXM],Ry[MAXN][MAXM];
int qx[9]={0,0,-1,-1,-1,0,1,1,1};
int qy[9]={0,-1,-1,0,1,1,1,0,-1};
int Gx[9]={0,1,1,0,1,1,1,0,1};
int Gy[9]={0,0,1,1,1,0,1,1,1};

int getHfunc(int num,int x,int y)
{
	if (num==1||num==5) return 2*Rx[x][y]+1;
	if (num==3||num==7) return 2*Ry[x][y]+1;
	if (num%2==0) return 2*(Rx[x][y]+Ry[x][y]+1);
}

void fastSEDT()
{
	memset(Rx,0,sizeof(Rx));
	memset(Ry,0,sizeof(Ry));
	for (int i=1;i<=N;i++)
		for (int j=1;j<=M;j++)
		{
			if (oriIma[i][j]==FOREGROUND_PIXEL)
			{
				resIma[i][j]=MAXVAL;
				for (int k=1;k<=4;k++)
				{
					int nx = i + qx[k];
					int ny = j + qy[k];
					if (resIma[i][j]>resIma[nx][ny]+getHfunc(k,nx,ny))
					{
						resIma[i][j] = resIma[nx][ny]+getHfunc(k,nx,ny);
						Rx[i][j] = Rx[nx][ny] + Gx[k];
						Ry[i][j] = Ry[nx][ny] + Gy[k];
					}
				}
			}
		}

	for (int i=N;i>=1;i--)
		for (int j=M;j>=1;j--)
		{
			if (oriIma[i][j]==FOREGROUND_PIXEL)
			{
				for (int k=5;k<=8;k++)
				{
					int nx = i + qx[k];
					int ny = j + qy[k];
					if (resIma[i][j]>resIma[nx][ny]+getHfunc(k,nx,ny))
					{
						resIma[i][j] = resIma[nx][ny]+getHfunc(k,nx,ny);
						Rx[i][j] = Rx[nx][ny] + Gx[k];
						Ry[i][j] = Ry[nx][ny] + Gy[k];
					}
				}
			}
		}
	return ;
}

int main()
{
	freopen("T:\\workspace\\Data\\in17.txt","r",stdin);
//	freopen("T:\\workspace\\Data\\out19Frank.txt","w",stdout);
	scanf("%d%d",&N,&M);
	for (int i=1;i<=N;i++)
		for (int j=1;j<=M;j++)
			scanf("%d",&oriIma[i][j]);
	clock_t t1 = clock();
	fastSEDT();
	clock_t t2 = clock();
	printf("%d\n",t2-t1);
	while(1){printf("");}
/*	for (int i=1;i<=N;i++)
	{
		for (int j = 1;j<M;j++) printf("%d ",resIma[i][j]);
		printf("%d\n",resIma[i][M]);
	}*/
	return 0;
}
