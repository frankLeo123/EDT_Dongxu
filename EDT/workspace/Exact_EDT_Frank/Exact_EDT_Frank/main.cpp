#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <time.h>

#define MAXN 1000+2
#define MAXM 1000+2
#define FOREGROUND_PIXEL 1
#define MAXVAL 2100000000
#define DEEP 5

struct rp
{
	int x,y;
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
	rp * pop()
	{
		rp * ans = NULL;
		if (this->head == this->tail && this->head != NULL)
		{
			ans = this->head;
			this->head = this->tail = NULL;
		}
		else if (this->head != this->tail && this->head != NULL)
		{
			ans = this->head;
			this->head = this->head->next;
		}
		return ans;
	}

	void append(List *list)
	{
		if (this->head == NULL)
		{
			this->head = list->head;
			this->tail = list->tail;
		}
		else
		{
			this->tail->next = list->head;
		}
	}
};

int N,M;
int oriIma[MAXN][MAXM];
int resIma[MAXN][MAXM];
List R[MAXN][MAXM];
int qx[9]={0,0,-1,-1,-1,0,1,1,1};
int qy[9]={0,-1,-1,0,1,1,1,0,-1};
int Gx[9]={0,1,1,0,1,1,1,0,1};
int Gy[9]={0,0,1,1,1,0,1,1,1};

int getHfunc(int num,int x,int y,List *hList)
{
	int ans = MAXVAL;
	List *now = &R[x][y];
	for (rp *i = now->head ; i != NULL ;i = i->next)
	{
		int tmp;
		if (num==1||num==5) tmp = 2*(i->x)+1;
		if (num==3||num==7) tmp = 2*(i->y)+1;
		if (num%2==0) tmp = 2*((i->x)+(i->y)+1);
		if (tmp<ans) 
		{
				ans = tmp;
				hList->clear();
				rp *tmp = (rp*)malloc(sizeof(rp));
				tmp->x = i->x;
				tmp->y = i->y;
				tmp->next = NULL;
				hList->push(tmp);
		}
		else if (tmp == ans)
		{
			rp *tmp = (rp*)malloc(sizeof(rp));
			tmp->x = i->x;
			tmp->y = i->y;
			tmp->next = NULL;
			hList->push(tmp);
		}
	}
	return ans;
}

void fastSEDT()
{
	memset(R,0,sizeof(R));
	for (int i=0;i<=N+1;i++)
		for (int j=0;j<=M+1;j++)
		{
			if (oriIma[i][j]==FOREGROUND_PIXEL)
			{
				resIma[i][j]=MAXVAL;
				for (int k=1;k<=4;k++)
				{
					int nx = i + qx[k];
					int ny = j + qy[k];
					List *idx = (List *)malloc(sizeof(List));
					idx->head = idx->tail = NULL;
					int H = getHfunc(k,nx,ny,idx);
					if (resIma[i][j] > resIma[nx][ny] + H )
					{
						resIma[i][j] = resIma[nx][ny] + H;
						R[i][j].clear();
						for (rp *it = idx->head; it != NULL ; it=it->next)
						{
							rp *tmp = (rp*)malloc(sizeof(rp));
							tmp->x = it->x + Gx[k];
							tmp->y = it->y + Gy[k];
							tmp->next = NULL;
							R[i][j].push(tmp);
						}
					}
					else if (resIma[i][j] == ( resIma[nx][ny] + H ) )
					{

						for (rp *it = idx->head; it != NULL ; it=it->next)
						{
							rp *tmp = (rp*)malloc(sizeof(rp));
							tmp->x = it->x + Gx[k];
							tmp->y = it->y + Gy[k];
							tmp->next = NULL;
							bool flag =false;
							for (rp *jt = R[i][j].head ; jt != NULL ; jt = jt->next)
							{
								if (tmp->x == jt->x && tmp->y == jt->y) {flag=true;break;}
							}
							if (flag) continue;
							R[i][j].push(tmp);
						}
					}
					idx->clear();
					free(idx);
				}
			}
			else
			{
				rp *tmp = (rp*)malloc(sizeof(rp));
				tmp->x = 0;
				tmp->y = 0;
				tmp->next = NULL;
				R[i][j].push(tmp);
				resIma[i][j]=0;
			}
		}



	for (int i=N;i>=1;i--)
		for (int j=M;j>=1;j--)
		{
			if (oriIma[i][j]==FOREGROUND_PIXEL)
			{
				if (i==32&&j==480)
				{
					printf("");
				}
				for (int k=5;k<=8;k++)
				{
					int nx = i + qx[k];
					int ny = j + qy[k];
					List *idx = (List *)malloc(sizeof(List));
					idx->head = idx->tail = NULL;
					int H = getHfunc(k,nx,ny,idx);
					if (resIma[i][j]>resIma[nx][ny]+H)
					{
						resIma[i][j] = resIma[nx][ny] + H;
						R[i][j].clear();
						for (rp *it = idx->head; it != NULL ; it=it->next)
						{
							rp *tmp = (rp*)malloc(sizeof(rp));
							tmp->x = it->x + Gx[k];
							tmp->y = it->y + Gy[k];
							tmp->next = NULL;
							R[i][j].push(tmp);
						}
					}
					else if (resIma[i][j] == ( resIma[nx][ny]+H ))
					{
						for (rp *it = idx->head; it != NULL ; it=it->next)
						{
							rp *tmp = (rp*)malloc(sizeof(rp));
							tmp->x = it->x + Gx[k];
							tmp->y = it->y + Gy[k];
							tmp->next = NULL;
							bool flag =false;
							for (rp *jt = R[i][j].head ; jt != NULL ; jt = jt->next)
							{
								if (tmp->x == jt->x && tmp->y == jt->y) {flag=true;break;}
							}
							if (flag) continue;
							R[i][j].push(tmp);
						}
					}
					idx->clear();
					free(idx);
				}
			}
		}
	return ;
}

int main()
{
    freopen("T:\\workspace\\Data\\in21.txt","r",stdin);
	freopen("T:\\workspace\\Data\\out21Frank_v2.txt","w",stdout);
	scanf("%d%d",&N,&M);
	for (int i=1;i<=N;i++)
		for (int j=1;j<=M;j++)
			scanf("%d",&oriIma[i][j]);
	clock_t t1 = clock();
	fastSEDT();
	clock_t t2 = clock();
	//printf("%d\n",t2-t1);
	for (int i=1;i<=N;i++)
	{
		for (int j = 1;j<M;j++) printf("%d ",resIma[i][j]);
		printf("%d\n",resIma[i][M]);
	}
	return 0;
}
