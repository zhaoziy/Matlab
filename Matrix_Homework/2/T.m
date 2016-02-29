function s=T()
clc;
clear;
syms x
eps=0.5*10e-8;
f=2*x/(x^2-3);
a=2;
b=3;
n=3000;

h=(b-a)/n;  %步长
sum_temp=0;
for i=(a+h):h:(b-h)
    sum_temp=sum_temp+subs(f,i);
end
sum_temp=double(sum_temp);
fa=subs(f,a);
fb=subs(f,b);
fa=double(fa);
fb=double(fb);
Tn=(b-a)/(2*n)*(fa+2*sum_temp+fb)
e=abs(Tn-log(6));

while(e>eps)
    n=n+500;
    h=(b-a)/n;  %步长
    sum_temp=0;
    for i=(a+h):h:(b-h)
        sum_temp=sum_temp+subs(f,i);
    end
    sum_temp=double(sum_temp);
    fa=subs(f,a);
    fb=subs(f,b);
    fa=double(fa);
    fb=double(fb);
    Tn=(b-a)/(2*n)*(fa+2*sum_temp+fb)
    e=abs(Tn-log(6))
end
Tn
n