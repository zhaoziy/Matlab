function x=fun1(A,b)
clc;
clear;
n=7;
i=1;
j=1;
while(i<n+1)
    for b=1:n
     A(i,j)=i^(j-1);
     j=j+1; 
    end
    i=i+1;
    j=j-n;
end
A
b=zeros(n,1);
i=1;
while(i<=n)
    if(i==1)
        b(i,1)=n;
        i=i+1;
    else
        b(i,1)=(i^n-1)/(i-1);
        i=i+1;
    end
end
b
x=zeros(n,1);
y=zeros(n,1);
temprow=zeros(n,1);
tempconstant=0;
Pvector=zeros(n,1);
for col=1:n-1
    [max_element,index]=max(abs(A(col:n,col)));
    temprow=A(col,:);
    A(col,:)=A(index+col-1,:);
    A(index+col-1,:)=temprow;
    tempconstant=b(col);
    b(col)=b(index+col-1);
    b(index+col-1)=tempconstant;
    if A(col,col)==0
        disp('A is singular.no unique solution');
        return
    end
    for row=col+1:n
        mult=A(row,col)/A(col,col);
        A(row,col)=mult;
        A(row,col+1:n)=A(row,col+1:n)-mult*A(col,col+1:n);
    end
end
y(1)=b(1);
for k=2:n
    y(k)=b(k)-A(k,1:k-1)*y(1:k-1);
end
x(n)=y(n)/A(n,n);
for k=n-1:-1:1
    x(k)=(y(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
