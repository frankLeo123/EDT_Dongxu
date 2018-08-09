import random
name = "in19.txt"
N = 5000
M = 5000


f = open(name,"w")
f.write(str(N)+" "+str(M)+"\n")

res = [[0 for x in range (0,N)] for y in range(0,M)]
random.seed()
num = int(random.uniform(5,1000));
#print(num)
for case in range(0,num):
    print(case/num)
    x1 = int(random.uniform(0,N))
    y1 = int(random.uniform(0,M))
    x2 = int(random.uniform(x1,N))
    y2 = int(random.uniform(y1,M))
    #print(x1,x2,y1,y2)
    for i in range(x1,min(x2+1,N)):
        for j in range(y1,min(y2+1,M)):
            res[i][j]=1
for i in range(0,N):
    for j in range(0,M):
        f.write(str(res[i][j])+" ")
    f.write("\n")
f.close()
