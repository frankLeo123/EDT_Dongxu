case = "22"
outS = "Data//out"+case+"mat.txt"
out = ["Data//out"+case+"_4.txt",
       "Data//out"+case+"_MP_CPU.txt",
       "Data//out"+case+"_MP_CPU_UNSYN.txt",
       "Data//out"+case+"Frank.txt",
       "Data//out"+case+"Frank_v2.txt",
       "Data//out"+case+"_3.txt",
       "Data//out"+case+"_2.txt",
       ]
fs = open(outS,"r")
N,M = map(int,fs.readline().split())
outsArr = [[0 for x in range(0,N)] for y in range(0,M)]

for i in range(0,N):
    tmp = list(map(int,fs.readline().split()))
    for j in range(0,M):
        outsArr[i][j]=tmp[j]
        
for i in range(0,len(out)):
    try:
        f = open(out[i],"r")
    except IOError:
        print(out[i]+"   Can't be open!")
        continue
    try:
        outArr = [[0 for x in range(0,N)] for y in range(0,M)]
        for j in range(0,N):
            tmp = list(map(int,f.readline().split()))
            for k in range(0,M):
                outArr[j][k]=tmp[k]
        num =0 ;
        for j in range(0,N):
            for k in range(0,M):
                if outArr[j][k]!=outsArr[j][k]:
                    num+=1
        res=float(float(num)/float(N*M))
        print(out[i]+"   "+str(res)+" "+str(num))
        f.close()
    except IndexError:
        print(out[i]+"   Error occur!")
        continue
fs.close()
