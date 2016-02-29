clc;
clear;
syms t u
U=(1/10)*exp(-2*t)+t*exp(-2*t);
f=exp(-2*t)-2*u;
h=0.05;
n=40;
a=0;b=h*n;
T=zeros(1,n+1);
Y=zeros(1,n+1);
T=a:h:b;
Y(1)=1/10;
for k=1:n
   K1=h*subs(f,{t,u},{T(k),Y(k)});
   K2=h*subs(f,{t,u},{T(k)+h/2,Y(k)+K1/2});
   K3=h*subs(f,{t,u},{T(k)+h/2,Y(k)+K2/2});
   K4=h*subs(f,{t,u},{T(k)+h,Y(k)+K3});
   Y(k+1)=Y(k)+(K1+2*K2+2*K3+K4)/6;
end
for m=1:(n+1)
    P(m)=subs(U,T(m));
end
P=double(P);
R=[T' Y'];
Y'
plot(T,Y,'--r^')
hold on
plot(T,P,'-b*')
hold off



    


