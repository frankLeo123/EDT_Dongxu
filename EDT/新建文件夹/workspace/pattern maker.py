N = 21
M = 21
res = [[0 for x in range(0,N)] for y in range(0,M)]
cx = int(N/2)
cy = int(M/2)
for i in range(0,N):
    for j in range(0,M-1):
        res[i][j] = (i-cx)*(i-cx)+(j-cy)*(j-cy)
        print("%3d " % res[i][j],end="")
    res[i][M-1] = (i-cx)*(i-cx)+(M-1-cy)*(M-1-cy)
    print(str(res[i][M-1])+"\n")
