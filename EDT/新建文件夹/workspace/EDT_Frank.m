fid=fopen('Data//in7.txt','r');
[data,count]=fscanf(fid,'%d');
N=data(1);
M=data(2);
A=zeros(N+2,M+2);
Rx = zeros(N+2,M+2);
Ry = zeros(N+2,M+2);
DD = zeros(N+2,M+2);
qx=[0,-1,-1,-1,0,1,1,1];
qy=[-1,-1,0,1,1,1,0,-1];
Gx=[1,1,0,1,1,1,0,1];
Gy=[0,1,1,1,0,1,1,1];
res = zeros(N+2,M+2);
resf = zeros(N+2,M+2);
for i=1:N
    for j=1:M 
        A(i,j) = data(2+(i-1)*M+j);
    end
end

for i=1:N
    for j=1:M
        if A(i,j)==1
            res(i,j)=2100000000;
            for k=1:4
                nx = i + qx(k);
                ny = j + qy(k);
                if nx==0 || nx >N || ny==0||ny>M
                    rx = 0;
                    ry = 0;
                    rest = 0;
                else
                    rx = Rx(nx,ny);
                    ry = Ry(nx,ny);
                    rest = res(nx,ny);
                end
                if k==1 || k==5
                    H = 2*rx+1;
                end
                if k==3 || k==7
                    H = 2*ry+1;
                end
                if mod(k,2)==0
                    H = 2*(rx+ry+1);
                end
                if res(i,j)>rest+H
                    res(i,j)=rest+H;
                    resf(i,j)=res(i,j);
                    DD(i,j)=k;
                    Rx(i,j)=rx+Gx(k);
                    Ry(i,j)=ry+Gy(k);
                end
            end
        end
    end
end

for i=N:-1:1
    for j=M:-1:1
        if A(i,j)==1
            for k=5:8
                nx = i + qx(k);
                ny = j + qy(k);
                if nx==0 || nx >N || ny==0||ny>M
                    rx = 0;
                    ry = 0;
                    rest = 0;
                else
                    rx = Rx(nx,ny);
                    ry = Ry(nx,ny);
                    rest = res(nx,ny);
                end
                if k==1 || k==5
                    H = 2*rx+1;
                end
                if k==3 || k==7
                    H = 2*ry+1;
                end
                if mod(k,2)==0
                    H = 2*(rx+ry+1);
                end
                if res(i,j)>rest+H
                    res(i,j)=rest+H;
                    Rx(i,j)=rx+Gx(k);
                    Ry(i,j)=ry+Gy(k);
                end
            end
        end
    end
end