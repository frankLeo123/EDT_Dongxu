fid=fopen('T:\\workspace\\Data\\in7.txt','r');
outf=fopen('T:\\workspace\\Data\\out7mat.txt','w');
[data,count]=fscanf(fid,'%d');
N=data(1);
M=data(2);
A=zeros(N+2,M+2);
for i=1:N+2
    for j=1:M+2
        if j==M+2 || i==N+2 || i==1|| M==1
            tmp = 0;
        else 
            tmp = data(2+(i-2)*M+j-1);
        end
        
        if tmp==0
            A(i,j)=1;
        else
            A(i,j)=0;
        end
    end
end
[D,L]=bwdist(A);
fprintf(outf,'%d %d\r\n',N,M);
for i=2:N+1
    for j=2:M
        D(i,j)=round(D(i,j)*D(i,j));
        fprintf(outf,'%d ',D(i,j));
    end
    fprintf(outf,'%d\r\n',D(i,M+1));
end