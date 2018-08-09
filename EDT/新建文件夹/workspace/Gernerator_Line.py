import random
name = "Data//in21.txt"
N = 1024
M = 1024
A = -1
B = 3
C = 2048

f = open(name,"w")
f.write(str(N)+" "+str(M)+"\n")

res = [[0 for x in range (0,N)] for y in range(0,M)]
random.seed()

for i in range(0,N):
    for j in range(0,M):
        res[i][j] = 1
for i in range(0,N):
    for j in range(0,M):
        if A*i+B*j>=C-2 and A*i+B*j<=C+2:
            res[i][j] = 0

for i in range(0,N):
    for j in range(0,M):
        f.write(str(res[i][j])+" ")
    f.write("\n")
f.close()
