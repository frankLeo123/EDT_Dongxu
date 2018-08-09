import random
name = "Data//in20.txt"
N = 1024
M = 1024
R = 400

f = open(name,"w")
f.write(str(N)+" "+str(M)+"\n")

res = [[0 for x in range (0,N)] for y in range(0,M)]
random.seed()

for i in range(0,N):
    for j in range(0,M):
        dis = (i - N/2)*(i - N/2) + (j - M/2)*(j - M/2)
        if dis<=R*R : res[i][j] = 1

for i in range(0,N):
    for j in range(0,M):
        f.write(str(res[i][j])+" ")
    f.write("\n")
f.close()
