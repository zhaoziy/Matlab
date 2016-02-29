clc;
clear;

n=6;
x=1;
y=0;
eps=10e-8;
count=1;

for i=1:(n^2-1)
x=[x ; 1];
end
a_t=2;
b_t=-1;
for i=1:(n-1)
a_t=[a_t 2];
end
for i=1:(n-2)
b_t=[b_t -1];
end
T=diag(a_t)+diag(b_t,1)+diag(b_t,-1);

syms e f;
c_t=T+2*eye(n);
d_t=-1*eye(n);
g_t=e;
h_t=f;
for i=1:(n-1)
g_t=[g_t e];
end
for i=1:(n-2)
h_t=[h_t f];
end
A=diag(g_t)+diag(h_t,1)+diag(h_t,-1);
A=subs(A,{e,f},{c_t,d_t});
A=double(A)
b=A*x;

n=n^2;
for i=1:n
x(i)=10;
end

U=-1*triu(A,1);
L=-1*tril(A,-1);
D=triu(A,0)-triu(A,1);

BG=inv(D-L)*U;
f=inv(D-L)*b;

flag=norm(x,2)
while(flag>eps)
    temp=x;
    x=BG*x+f
    flag=norm((x-temp),2)
    fprintf('所用迭代次数:%d',count);
    count=count+1;
end
x'