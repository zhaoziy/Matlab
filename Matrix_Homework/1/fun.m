function[x]=fun(f)
%牛顿法
clc;
clear;
syms x;
f=2*x^3-5*x^2-19*x+42;
df=diff(f,x);
x0=50;
epsilon=0.00001;
count=1;
fvalue=subs(f,x0);
dfvalue=subs(df,x0);
fvalue=double(fvalue);
dfvalue=double(dfvalue);
x1=x0-fvalue/dfvalue;
reerr=abs((x1-x0)/x1);
p(count)=reerr;
fprintf('%f , %f\n',x1,reerr)
count=count+1;
while(reerr>epsilon)
    x0=x1;
    fvalue=subs(f,x0);
    dfvalue=subs(df,x0);
    fvalue=double(fvalue);
    dfvalue=double(dfvalue);
    x1=x0-fvalue/dfvalue;
    reerr=2*abs(x1-x0)/x1;
    p(count)=reerr;
    fprintf('%f , %f\n',x1,reerr)
    count=count+1;
end
b=1:1:count-1;
plot(b,p,'-+')
fprintf('所用迭代次数:%d',count-1);