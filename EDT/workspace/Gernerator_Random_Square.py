import random
import math
name = "Data//in22.txt"
N = 1024
M = 1024
thete = math.pi /6.0 ;
Num = 10
MaxSize = 50

f = open(name,"w")
f.write(str(N)+" "+str(M)+"\n")

res = [[0 for x in range (0,N)] for y in range(0,M)]
random.seed()


for i in range(0,N):
    for j in range(0,M):
        res[i][j] = 1

for k in range(0,Num):
    size = int(random.uniform(MaxSize/4,MaxSize))
    ox = int(random.uniform(size,N))
    oy = int(random.uniform(size,M))
    bsize = int(size/2)
    for i in range(-bsize,bsize):
        for j in range(-bsize,bsize):
            nx = int(i*math.cos(thete)-j*math.sin(thete))
            ny = int(i*math.sin(thete)+j*math.cos(thete))
            res[nx+ox][ny+oy]=0
            nx = int(i*math.cos(thete)-j*math.sin(thete)+0.5)
            ny = int(i*math.sin(thete)+j*math.cos(thete)+0.5)
            res[nx+ox][ny+oy]=0

for i in range(0,N):
    for j in range(0,M):
        f.write(str(res[i][j])+" ")
    f.write("\n")
f.close()
