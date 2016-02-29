clc;
clear;

n=6;
x=1;
y=0;

for i=1:(n^2-1)
x=[x ; 1];
y=[y ; 0];
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
b=A*x

n=n^2;
L=zeros(n);

for j=1:n
    if(j==1)
        k=1;
        L(j,j)=sqrt(A(j,j)-L(j,k)*L(j,k));
        for i=(j+1):n
            L(i,j)=(A(i,j)-L(i,k)*L(j,k))/L(j,j);
        end
    else
        temp1=0;
        temp2=0;
        for k=1:(j-1)
            temp1=temp1+L(j,k)*L(j,k);
            temp2=temp2+L(i,k)*L(j,k);
        end
        L(j,j)=sqrt(A(j,j)-temp1);
        for i=(j+1):n
            L(i,j)=(A(i,j)-temp2)/L(j,j);
        end
    end
end

y=inv(L)*b;
x=inv(L')*y;
x'