function s=simpson()
clc;
clear;
syms x
eps=0.5*10e-8;
f=2*x/(x^2-3);
a=2;
b=3;
n=10;
h=(b-a)/n;  %步长
sum_temp1=0;
sum_temp2=0;
for i=(a+h/2):h:(b-h/2)
    sum_temp1=sum_temp1+subs(f,i);
end
sum_temp1=double(sum_temp1);
for i=(a+h):h:(b-h)
    sum_temp2=sum_temp2+subs(f,i);
end
sum_temp2=double(sum_temp2);

fa=subs(f,a);
fb=subs(f,b);
fa=double(fa);
fb=double(fb);
Sim=(b-a)/(6*n)*(fa+4*sum_temp1+2*sum_temp2+fb)
e=abs(Sim-log(6));

while(e>eps)
    n=n+10;
    h=(b-a)/n;  %步长
    sum_temp1=0;
    sum_temp2=0;
    for i=(a+h/2):h:(b-h/2)
        sum_temp1=sum_temp1+subs(f,i);
    end
    sum_temp1=double(sum_temp1);
    for i=(a+h):h:(b-h)
        sum_temp2=sum_temp2+subs(f,i);
    end
    sum_temp2=double(sum_temp2);
    fa=subs(f,a);
    fb=subs(f,b);
    fa=double(fa);
    fb=double(fb);
    Sim=(b-a)/(6*n)*(fa+4*sum_temp1+2*sum_temp2+fb)
    e=abs(Sim-log(6))
end
Sim
n