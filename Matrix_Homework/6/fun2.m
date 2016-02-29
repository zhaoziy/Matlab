function x=fun2(f)
syms x
% f=1/(1+25*x^2);
f=atan(x);
% f=x/(1+x^4);
n=3;
for i=1:n
    X(1)=-1;
    X(i+1)=-1+2*i/n;
    y(1)=subs(f,-1);
    y(i+1)=subs(f,X(i));
end
X;
%y;
p=0;
for k=1:n
    index=setdiff(1:n,k);
    p=p+y(k)*prod((x-X(index))./(X(k)-X(index)));
end
p
Y=-1:0.1:1;
%u=subs(p,Y)
subplot(1,2,1);
plot(Y,subs(p,Y),'--^')
subplot(1,2,2);
plot(Y,subs(f,Y),'--*')
