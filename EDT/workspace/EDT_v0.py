class Point:
    def __init__(self,x,y,fx,fy,l):
        self.x = x
        self.y = y
        self.fx = fx
        self.fy = fy
        self.level = l
    def getSED(self):
        return (self.x-self.fx)*(self.x-self.fx)+(self.y-self.fy)*(self.y-self.fy)
        
def findInit():
    for i in range(0,N+2):
        for j in range(0,M+2):
            if ori[i][j] == 0:
                que.append(Point(i,j,i,j,0))
                res[i][j]=0

def EDT():
    dx = [0,1,0,-1]
    dy = [-1,0,1,0]
    while len(que)!=0:
        now = que.pop(0)
        nowLevel = now.level
        x=now.x
        y=now.y
        fx=now.fx
        fy=now.fy
        for i in range(0,4):
            nx = x+dx[i]
            ny = y+dy[i]
            if nx<0 or nx>N+1 or ny<0 or ny>M+1:continue
            np = Point(nx,ny,fx,fy,nowLevel+1)
            dis = np.getSED()
            if dis<res[nx][ny]:
                res[nx][ny]=dis
                que.append(np)

inf = open("in-large.txt","r")
outf = open("out12.txt","w");
N,M  = map(int, inf.readline().split())
print(N,M)

que = []
MAXVAL = 2000000000
ori = [[0 for x in range(0,N+2)] for y in range(0,M+2)]
res = [[MAXVAL for x in range(0,N+2)] for y in range(0,M+2)]

for i in range(1,N+1):
    tmp = list(map(int,inf.readline().split()))
    for j in range(1,M+1):
        ori[i][j]=tmp[j-1]

findInit()
EDT()

for i in range(1,N+1):
    for j in range(1,M):
        outf.write("%d " % res[i][j])
    outf.write("%d\n" % res[i][M])

inf.close()
outf.close()

